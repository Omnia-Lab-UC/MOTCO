"""Shared pooled preprocessing for semi-synthetic integration and diagnostics."""

from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Any, Mapping, Sequence, TypeVar, cast

import numpy as np
import pandas as pd

from motco.simulations.generator import logit
from motco.simulations.semisynthetic import OmicsLayer, PopulationTrajectories, SemiSyntheticTrajectoryDataset

OMIC_LAYERS: tuple[OmicsLayer, ...] = ("methylation", "expression", "proteomics")
_SCALE_TOLERANCE = 1e-10

#: A block name: any string, or the synthetic path's :data:`OmicsLayer` literal.
BlockName = TypeVar("BlockName", bound=str)

#: Per-block unit transforms applied before pooled standardization.
UNIT_TRANSFORMS: tuple[str, ...] = ("identity", "beta_to_mvalue")

#: The synthetic datasets' unit contract: methylation betas are analysed as M-values.
_SYNTHETIC_TRANSFORMS: Mapping[str, str] = {"methylation": "beta_to_mvalue"}

#: Label-blind per-block weightings of the pooled-standardized blocks.
#: ``equal_variance`` multiplies a block by 1/√p_b (unit total variance);
#: ``mfa`` by 1/√λ₁, λ₁ the leading eigenvalue of the standardized block's
#: covariance (unit variance of its leading principal component).
BLOCK_WEIGHTINGS: tuple[str, ...] = ("equal_variance", "mfa")


def validate_block_weighting(weighting: Any) -> str | None:
    """Return a requested block weighting kind, or ``None`` when absent.

    An unknown kind raises :class:`ValueError` naming it and the allowed kinds.
    """

    if weighting is None:
        return None
    if not isinstance(weighting, str) or weighting not in BLOCK_WEIGHTINGS:
        raise ValueError(f"Unknown block weighting {weighting!r}; allowed: {list(BLOCK_WEIGHTINGS)} or None.")
    return weighting


def selected_layers(integration_params: Mapping[str, Any] | None) -> tuple[OmicsLayer, ...]:
    """Resolve ``integration_params["layers"]`` to canonical omic order.

    Absent (or ``None``) means every layer. A present selection must be a
    non-empty, duplicate-free subset of :data:`OMIC_LAYERS`.
    """

    raw = None if integration_params is None else integration_params.get("layers")
    if raw is None:
        return OMIC_LAYERS
    if isinstance(raw, str | bytes) or not isinstance(raw, Sequence):
        raise ValueError(f"integration_params['layers'] must be a list of layer names from {list(OMIC_LAYERS)}.")
    names = [str(name) for name in raw]
    if not names:
        raise ValueError(f"integration_params['layers'] must not be empty; allowed layers: {list(OMIC_LAYERS)}.")
    unknown = sorted(set(names) - set(OMIC_LAYERS))
    if unknown:
        raise ValueError(
            f"integration_params['layers'] names unknown layer(s) {unknown}; allowed layers: {list(OMIC_LAYERS)}."
        )
    if len(set(names)) != len(names):
        raise ValueError(
            f"integration_params['layers'] repeats a layer: {names}; allowed layers: {list(OMIC_LAYERS)}."
        )
    return tuple(cast(OmicsLayer, layer) for layer in OMIC_LAYERS if layer in names)


@dataclass(frozen=True)
class BlockScaler:
    """Fitted pooled location and scale for one aligned omic block.

    ``weight`` multiplies every standardized feature of the block (the block
    weighting; 1.0, exact in floating point, when no weighting was requested).
    """

    feature_names: tuple[str, ...]
    mean: np.ndarray
    scale: np.ndarray
    weight: float = 1.0


@dataclass(frozen=True)
class FittedOmicsPreprocessor:
    """Pooled per-feature scalers in canonical omic order.

    ``scalers`` holds exactly the layers the preprocessor was fitted on;
    :attr:`layers` lists them in canonical order, and every transform covers
    those layers only. ``block_weighting`` names the kind the scalers' weights
    were fitted with (``None``: every weight is 1).
    """

    scalers: Mapping[OmicsLayer, BlockScaler]
    methylation_units: str = "mvalue"
    block_weighting: str | None = None

    @property
    def layers(self) -> tuple[OmicsLayer, ...]:
        return tuple(layer for layer in OMIC_LAYERS if layer in self.scalers)

    def require_layers(self, layers: Sequence[OmicsLayer]) -> None:
        """Raise unless the preprocessor was fitted on exactly ``layers``."""

        if tuple(layers) != self.layers:
            raise ValueError(
                f"Preprocessor was fitted on layers {list(self.layers)}, not the requested {list(layers)}."
            )

    def require_block_weighting(self, weighting: str | None) -> None:
        """Raise unless the preprocessor was fitted with block weighting ``weighting``."""

        if weighting != self.block_weighting:
            raise ValueError(
                f"Preprocessor was fitted with block_weighting={self.block_weighting!r}, not the requested "
                f"{weighting!r}."
            )

    def transform_dataset(self, dataset: SemiSyntheticTrajectoryDataset) -> dict[OmicsLayer, pd.DataFrame]:
        """Transform observed blocks using the fitted feature contract."""

        blocks = {layer: getattr(dataset, layer) for layer in self.layers}
        transformed = transform_blocks(blocks, self.scalers, _SYNTHETIC_TRANSFORMS)
        return {layer: frame.set_axis(frame.index.astype(str), axis=0) for layer, frame in transformed.items()}

    def transform_population(self, population: PopulationTrajectories) -> dict[OmicsLayer, pd.DataFrame]:
        """Transform analytic means with the same observed-fitted scalers.

        Population means are already on the integration scale (M-values for
        methylation), so no unit transform is applied.
        """

        blocks = {layer: population.layers[layer] for layer in self.layers}
        return transform_blocks(blocks, self.scalers, {})


def fit_omics_preprocessor(
    dataset: SemiSyntheticTrajectoryDataset,
    layers: Sequence[OmicsLayer] = OMIC_LAYERS,
    weighting: str | None = None,
) -> FittedOmicsPreprocessor:
    """Fit pooled block-wise scalers on the observed dataset's ``layers``.

    ``weighting`` (one of :data:`BLOCK_WEIGHTINGS`, or ``None``) fits each
    block's weight with the scalers (:func:`fit_block_scalers`).
    """

    weighting = validate_block_weighting(weighting)
    blocks = {layer: getattr(dataset, layer) for layer in layers}
    return FittedOmicsPreprocessor(
        scalers=fit_block_scalers(blocks, _SYNTHETIC_TRANSFORMS, weighting=weighting),
        block_weighting=weighting,
    )


def fit_block_scalers(
    blocks: Mapping[BlockName, pd.DataFrame],
    transforms: Mapping[str, str],
    *,
    weighting: str | None = None,
) -> dict[BlockName, BlockScaler]:
    """Fit pooled per-feature location and scale for each named block.

    ``transforms`` maps a block name to its unit transform (one of
    :data:`UNIT_TRANSFORMS`); a block it does not name is ``identity``. Scalers
    are fitted on the transformed values, and a feature whose pooled standard
    deviation is (numerically) zero keeps unit scale, so it is centered only.

    ``weighting`` (one of :data:`BLOCK_WEIGHTINGS`; ``None``, the default,
    leaves every weight at 1) fits one label-blind scalar per block on its
    pooled-standardized values ``Z_b`` (n × p_b): ``equal_variance`` is
    1/√p_b, constant features included, so a block of non-constant features
    has unit total variance; ``mfa`` is √n / σ₁(Z_b), so the block's leading
    principal component has unit variance. An unknown kind is rejected before
    any block is transformed, and an ``mfa`` block whose features are all
    constant is an error naming it.
    """

    weighting = validate_block_weighting(weighting)
    scalers: dict[BlockName, BlockScaler] = {}
    for name, matrix in blocks.items():
        values = apply_unit_transform(name, matrix.astype(float).to_numpy(dtype=float), transforms.get(name))
        mean = values.mean(axis=0)
        scale = values.std(axis=0)
        scale[scale < _SCALE_TOLERANCE] = 1.0
        scalers[name] = BlockScaler(
            feature_names=tuple(matrix.columns.astype(str)),
            mean=mean,
            scale=scale,
            weight=1.0 if weighting is None else _block_weight(name, (values - mean) / scale, weighting),
        )
    return scalers


def _block_weight(name: str, standardized: np.ndarray, weighting: str) -> float:
    """The weight of one pooled-standardized block (n × p_b) under ``weighting``."""

    n_samples, n_features = standardized.shape
    if weighting == "equal_variance":
        return 1.0 / math.sqrt(n_features)
    leading = float(np.linalg.svd(standardized, compute_uv=False)[0]) if standardized.size else 0.0
    leading_sd = leading / math.sqrt(n_samples) if n_samples else 0.0
    if leading_sd < _SCALE_TOLERANCE:
        raise ValueError(
            f"Block {name!r} has no variance, so its mfa weight (1 / the leading principal component's "
            "standard deviation) is undefined."
        )
    return 1.0 / leading_sd


def transform_blocks(
    blocks: Mapping[BlockName, pd.DataFrame],
    scalers: Mapping[BlockName, BlockScaler],
    transforms: Mapping[str, str],
) -> dict[BlockName, pd.DataFrame]:
    """Standardize each named block with its fitted scaler, times its weight.

    Every block must carry exactly the fitted feature order; rows keep the
    block's own index. A weight of 1 (no weighting) multiplies exactly, so the
    unweighted output is unchanged.
    """

    transformed: dict[BlockName, pd.DataFrame] = {}
    for name, matrix in blocks.items():
        scaler = scalers[name]
        _validate_features(name, matrix.columns, scaler.feature_names)
        values = apply_unit_transform(name, matrix.astype(float).to_numpy(dtype=float), transforms.get(name))
        transformed[name] = pd.DataFrame(
            (values - scaler.mean) / scaler.scale * scaler.weight,
            index=matrix.index,
            columns=scaler.feature_names,
        )
    return transformed


def apply_unit_transform(name: str, values: np.ndarray, transform: str | None) -> np.ndarray:
    """Map one block's values through its declared unit transform.

    ``identity`` (or ``None``) returns the values unchanged. ``beta_to_mvalue``
    is the generator's clipped logit of methylation beta values. The value range
    is not checked here — the synthetic path has always clipped — so callers
    taking real betas validate the open interval ``(0, 1)`` themselves.
    """

    if transform is None or transform == "identity":
        return values
    if transform == "beta_to_mvalue":
        return logit(values)
    raise ValueError(
        f"Block {name!r} declares unknown unit transform {transform!r}; allowed: {list(UNIT_TRANSFORMS)}."
    )


def concatenate_blocks(
    blocks: Mapping[BlockName, pd.DataFrame],
    order: Sequence[BlockName] | None = None,
) -> pd.DataFrame:
    """Concatenate aligned standardized blocks with collision-proof names.

    Blocks are taken in ``order`` when given; otherwise in canonical omic order,
    concatenating only the canonical layers present in ``blocks``. Every joint
    column is named ``f"{block}__{feature}"``.
    """

    names: list[Any] = [layer for layer in OMIC_LAYERS if layer in blocks] if order is None else list(order)
    frames = []
    for name in names:
        matrix = blocks[name]
        renamed = matrix.copy()
        renamed.columns = [f"{name}__{column}" for column in matrix.columns.astype(str)]
        frames.append(renamed)
    return pd.concat(frames, axis=1)


def _validate_features(layer: str, columns: pd.Index, expected: tuple[str, ...]) -> None:
    actual = tuple(columns.astype(str))
    if actual != expected:
        raise ValueError(f"{layer} feature order does not match the fitted preprocessing artifact.")
