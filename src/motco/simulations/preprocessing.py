"""Shared pooled preprocessing for semi-synthetic integration and diagnostics."""

from __future__ import annotations

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
    """Fitted pooled location and scale for one aligned omic block."""

    feature_names: tuple[str, ...]
    mean: np.ndarray
    scale: np.ndarray


@dataclass(frozen=True)
class FittedOmicsPreprocessor:
    """Pooled per-feature scalers in canonical omic order.

    ``scalers`` holds exactly the layers the preprocessor was fitted on;
    :attr:`layers` lists them in canonical order, and every transform covers
    those layers only.
    """

    scalers: Mapping[OmicsLayer, BlockScaler]
    methylation_units: str = "mvalue"

    @property
    def layers(self) -> tuple[OmicsLayer, ...]:
        return tuple(layer for layer in OMIC_LAYERS if layer in self.scalers)

    def require_layers(self, layers: Sequence[OmicsLayer]) -> None:
        """Raise unless the preprocessor was fitted on exactly ``layers``."""

        if tuple(layers) != self.layers:
            raise ValueError(
                f"Preprocessor was fitted on layers {list(self.layers)}, not the requested {list(layers)}."
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
) -> FittedOmicsPreprocessor:
    """Fit pooled block-wise scalers on the observed dataset's ``layers``."""

    blocks = {layer: getattr(dataset, layer) for layer in layers}
    return FittedOmicsPreprocessor(scalers=fit_block_scalers(blocks, _SYNTHETIC_TRANSFORMS))


def fit_block_scalers(
    blocks: Mapping[BlockName, pd.DataFrame],
    transforms: Mapping[str, str],
) -> dict[BlockName, BlockScaler]:
    """Fit pooled per-feature location and scale for each named block.

    ``transforms`` maps a block name to its unit transform (one of
    :data:`UNIT_TRANSFORMS`); a block it does not name is ``identity``. Scalers
    are fitted on the transformed values, and a feature whose pooled standard
    deviation is (numerically) zero keeps unit scale, so it is centered only.
    """

    scalers: dict[BlockName, BlockScaler] = {}
    for name, matrix in blocks.items():
        values = apply_unit_transform(name, matrix.astype(float).to_numpy(dtype=float), transforms.get(name))
        scale = values.std(axis=0)
        scale[scale < _SCALE_TOLERANCE] = 1.0
        scalers[name] = BlockScaler(
            feature_names=tuple(matrix.columns.astype(str)),
            mean=values.mean(axis=0),
            scale=scale,
        )
    return scalers


def transform_blocks(
    blocks: Mapping[BlockName, pd.DataFrame],
    scalers: Mapping[BlockName, BlockScaler],
    transforms: Mapping[str, str],
) -> dict[BlockName, pd.DataFrame]:
    """Standardize each named block with its fitted scaler.

    Every block must carry exactly the fitted feature order; rows keep the
    block's own index.
    """

    transformed: dict[BlockName, pd.DataFrame] = {}
    for name, matrix in blocks.items():
        scaler = scalers[name]
        _validate_features(name, matrix.columns, scaler.feature_names)
        values = apply_unit_transform(name, matrix.astype(float).to_numpy(dtype=float), transforms.get(name))
        transformed[name] = pd.DataFrame(
            (values - scaler.mean) / scaler.scale,
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
