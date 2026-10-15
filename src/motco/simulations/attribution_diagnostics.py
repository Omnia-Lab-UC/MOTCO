"""Compact, JSON-safe orientation-attribution diagnostics for simulation studies.

This module is the *simulation adapter* over :mod:`motco.stats.attribution`. It
never changes the attribution contract: it hands the frozen fitted PLS model and
the exact pooled standardized feature matrix used for trajectory measurement to
:func:`~motco.stats.attribution.analyze_orientation_attribution`, and then
reduces the full in-memory result to a bounded record suitable for JSONL
persistence.

What is *not* persisted is as load-bearing as what is: no fitted estimator, no
full standardized matrix, no bootstrap matrix, and no unrestricted feature
table. Aggregate stability and truth-recovery metrics are computed from the
complete in-memory result and only their summaries are written, while feature
identifiers are truncated to the configured ``top_k`` so record size does not
scale with the feature count.

Decompositions
--------------
The bounded record carries both of attribution's decompositions, distinguished by
the ``transition_id`` field: the adjacent stage transitions (``s -> s+1``), and
the principal-orientation block under the reserved identifier from
``result.config.principal_component_id``. The latter decomposes the tested
``angle`` estimand; its truth-recovery row is scored against the union of
per-transition drivers, since it summarizes the whole configuration.

Units
-----
Feature effects are reported in pooled standardized units and, where the fitted
preprocessor exposes a positive per-feature scale, in original units too. The
methylation block's original units are **M-values**, not beta values — the
record labels this explicitly so a reader cannot mistake one for the other.
"""

from __future__ import annotations

from dataclasses import dataclass
from time import perf_counter
from typing import Any, Mapping, Sequence

import numpy as np
import pandas as pd

from motco.simulations.preprocessing import OMIC_LAYERS
from motco.simulations.semisynthetic import OmicsLayer, SemiSyntheticTrajectoryDataset
from motco.stats.attribution import (
    AttributionError,
    OrientationAttributionResult,
    analyze_orientation_attribution,
)

#: Version of the persisted attribution diagnostic record contract.
#: 1 = per-adjacent-transition decompositions only;
#: 2 = the principal-orientation block is reported beside them, so ``top_features``,
#: ``truth_recovery``, and ``stability`` each carry one extra identifier — records
#: written under the two contracts must never be mixed within a shard.
ATTRIBUTION_SCHEMA_VERSION = 2

#: Original-unit basis per omic layer, recorded so M-values are never read as betas.
UNIT_BASIS: dict[str, str] = {
    "methylation": "mvalue",
    "expression": "expression",
    "proteomics": "abundance",
}

_COMPONENTS: tuple[str, ...] = ("observed", "pls_captured", "residual")
_TRUTH_TOLERANCE = 1e-12


class AttributionDiagnosticError(ValueError):
    """Raised when attribution diagnostics cannot be computed for a replicate."""


@dataclass(frozen=True)
class TruthDrivers:
    """Generator-truth driver features aligned to the joint feature order."""

    feature_names: tuple[str, ...]
    global_mask: np.ndarray
    per_transition: dict[str, np.ndarray]
    available: bool
    definition: str

    def names(self, transition_id: str | None = None) -> list[str]:
        mask = self.global_mask if transition_id is None else self.per_transition[transition_id]
        return [self.feature_names[index] for index in np.flatnonzero(mask)]


def derive_truth_driver_features(
    dataset: SemiSyntheticTrajectoryDataset,
    feature_names: Sequence[str],
    stage_levels: Sequence[str],
    layers: Sequence[OmicsLayer] = OMIC_LAYERS,
) -> TruthDrivers:
    """Aligned truth drivers from changed group-stage differential patterns.

    A feature is a truth driver for the transition ``s -> s+1`` when the two
    groups' *population* mean changes over that transition differ, i.e. when

    ``delta_A * (v_A[s+1] - v_A[s]) != delta_B * (v_B[s+1] - v_B[s])``

    for that feature's per-omic effect size ``delta`` and per-stage differential
    indicator ``v``. Because expression and protein indicators are re-derived
    from methylation through the CpG→gene→protein incidence maps, propagated
    cascade drivers satisfy this definition on their own layers — a real
    downstream driver is therefore never counted as a false positive. The
    definition also covers effect-size-only constructions such as
    ``magnitude_kind='all'``, where the indicators are identical between groups
    but the per-omic delta is not.

    The global driver set is the union over transitions. ``layers`` names the
    measured blocks, in canonical order, that make up the joint feature order.
    """

    truth = dataset.truth
    indicators = truth.get("indicators")
    group_labels = truth.get("group_labels")
    deltas = truth.get("deltas")
    names = tuple(str(name) for name in feature_names)
    n_transitions = max(len(stage_levels) - 1, 0)
    transition_ids = [
        f"{stage_levels[index]}->{stage_levels[index + 1]}" for index in range(n_transitions)
    ]
    if not indicators or not group_labels or not deltas:
        empty = np.zeros(len(names), dtype=bool)
        return TruthDrivers(
            feature_names=names,
            global_mask=empty,
            per_transition={transition_id: empty.copy() for transition_id in transition_ids},
            available=False,
            definition="unavailable: dataset truth carries no differential indicators",
        )

    label_a, label_b = (str(label) for label in group_labels)
    per_transition: dict[str, np.ndarray] = {}
    blocks: list[list[np.ndarray]] = [[] for _ in transition_ids]
    for layer in layers:
        layer_index = OMIC_LAYERS.index(layer)
        ind_a = np.asarray(indicators[label_a][layer], dtype=float)
        ind_b = np.asarray(indicators[label_b][layer], dtype=float)
        delta_a = float(deltas[label_a][layer_index])
        delta_b = float(deltas[label_b][layer_index])
        for index in range(n_transitions):
            step_a = delta_a * (ind_a[:, index + 1] - ind_a[:, index])
            step_b = delta_b * (ind_b[:, index + 1] - ind_b[:, index])
            blocks[index].append(np.abs(step_a - step_b) > _TRUTH_TOLERANCE)

    for index, transition_id in enumerate(transition_ids):
        mask = np.concatenate(blocks[index]) if blocks[index] else np.zeros(0, dtype=bool)
        if mask.shape[0] != len(names):
            raise AttributionDiagnosticError(
                f"Truth indicators ({mask.shape[0]}) do not align to the joint feature order ({len(names)})."
            )
        per_transition[transition_id] = mask

    global_mask = (
        np.logical_or.reduce(list(per_transition.values()))
        if per_transition
        else np.zeros(len(names), dtype=bool)
    )
    return TruthDrivers(
        feature_names=names,
        global_mask=global_mask,
        per_transition=per_transition,
        available=True,
        definition=(
            "features whose group-stage differential mean change differs between groups, "
            "including CpG->gene->protein propagated effects"
        ),
    )


def compute_attribution_diagnostics(
    dataset: SemiSyntheticTrajectoryDataset | None = None,
    *,
    model: Any,
    features: pd.DataFrame,
    original_scales: np.ndarray | None,
    settings: Any,
    group_col: str,
    stage_col: str,
    selected_components: int | None,
    feature_order_signature: str | None,
    methylation_units: str = "mvalue",
    layers: Sequence[OmicsLayer] = OMIC_LAYERS,
    metadata: pd.DataFrame | None = None,
    groups: Sequence[str] | None = None,
    stages: Sequence[str] | None = None,
    truth: TruthDrivers | None = None,
    original_basis: Mapping[str, str] | None = None,
) -> dict[str, Any]:
    """Run frozen-model attribution and reduce it to a compact JSON-safe record.

    See :func:`run_attribution_diagnostics` for the arguments; this returns only
    the compact record.
    """

    record, _ = run_attribution_diagnostics(
        dataset,
        model=model,
        features=features,
        original_scales=original_scales,
        settings=settings,
        group_col=group_col,
        stage_col=stage_col,
        selected_components=selected_components,
        feature_order_signature=feature_order_signature,
        methylation_units=methylation_units,
        layers=layers,
        metadata=metadata,
        groups=groups,
        stages=stages,
        truth=truth,
        original_basis=original_basis,
    )
    return record


def run_attribution_diagnostics(
    dataset: SemiSyntheticTrajectoryDataset | None = None,
    *,
    model: Any,
    features: pd.DataFrame,
    original_scales: np.ndarray | None,
    settings: Any,
    group_col: str,
    stage_col: str,
    selected_components: int | None,
    feature_order_signature: str | None,
    methylation_units: str = "mvalue",
    layers: Sequence[OmicsLayer] = OMIC_LAYERS,
    metadata: pd.DataFrame | None = None,
    groups: Sequence[str] | None = None,
    stages: Sequence[str] | None = None,
    truth: TruthDrivers | None = None,
    original_basis: Mapping[str, str] | None = None,
) -> tuple[dict[str, Any], OrientationAttributionResult]:
    """Run frozen-model attribution; return the compact record and the full result.

    With a synthetic ``dataset``, the metadata and generator-truth drivers are
    read from it (``layers`` names its measured blocks). Without one, the caller
    supplies row-aligned ``metadata`` and, optionally, explicit ``groups`` and
    ``stages`` orders and a ``truth`` (an unavailable :class:`TruthDrivers` for
    real data). ``original_basis`` overrides the per-block original-unit labels
    (default :data:`UNIT_BASIS`).
    """

    if metadata is None:
        if dataset is None:
            raise AttributionDiagnosticError("Attribution needs either a dataset or row-aligned metadata.")
        metadata = dataset.metadata
    if truth is None and dataset is None:
        raise AttributionDiagnosticError("Attribution without a dataset needs an explicit truth record.")
    start = perf_counter()
    metadata = metadata.reset_index(drop=True)
    aligned_features = features.reset_index(drop=True)
    try:
        result = analyze_orientation_attribution(
            aligned_features,
            metadata,
            model,
            group_col=group_col,
            stage_col=stage_col,
            groups=groups,
            stages=stages,
            original_scales=(
                None if original_scales is None else np.asarray(original_scales, dtype=float).tolist()
            ),
            bootstrap_replicates=int(settings.bootstrap_replicates),
            bootstrap_seed=settings.bootstrap_seed,
            top_k=int(settings.top_k),
            zero_tolerance=float(settings.zero_tolerance),
        )
    except AttributionError as exc:
        raise AttributionDiagnosticError(f"Orientation attribution failed: {exc}") from exc
    attribution_seconds = perf_counter() - start

    if truth is None:
        assert dataset is not None
        truth = derive_truth_driver_features(
            dataset,
            result.feature_names,
            list(result.config.stages),
            layers,
        )
    transitions = _transition_records(result)
    top_features = _top_feature_records(result, int(settings.top_k))
    recovery = _truth_recovery_records(result, truth, int(settings.top_k))
    stability = _stability_records(result)

    record = {
        "schema_version": ATTRIBUTION_SCHEMA_VERSION,
        "status": "computed",
        "settings": {
            "enabled": True,
            "bootstrap_replicates": int(settings.bootstrap_replicates),
            "bootstrap_seed": settings.bootstrap_seed,
            "top_k": int(settings.top_k),
            "zero_tolerance": float(settings.zero_tolerance),
        },
        "model": {
            "selected_components": (
                int(selected_components) if selected_components is not None else int(result.config.model_components)
            ),
            "model_components": int(result.config.model_components),
            "feature_order_signature": feature_order_signature,
            "n_features": int(len(result.feature_names)),
            "groups": list(result.config.groups),
            "stages": list(result.config.stages),
        },
        "units": {
            "standardized": "pooled_standardized_input",
            "original_available": bool(result.config.original_units_available),
            "original_basis": dict(UNIT_BASIS if original_basis is None else original_basis),
            "methylation_units": methylation_units,
        },
        "truth": {
            "available": truth.available,
            "definition": truth.definition,
            "n_drivers": int(truth.global_mask.sum()),
            "n_drivers_by_transition": {
                transition_id: int(mask.sum()) for transition_id, mask in truth.per_transition.items()
            },
        },
        "transitions": transitions,
        "top_features": top_features,
        "truth_recovery": recovery,
        "stability": stability,
        "runtime": {
            "attribution_seconds": float(attribution_seconds),
            "bootstrap_replicates": int(settings.bootstrap_replicates),
        },
    }
    return record, result


def unavailable_record(reason: str, *, status: str = "not_requested") -> dict[str, Any]:
    """A JSON-safe marker for replicates that carry no attribution result.

    ``status`` distinguishes a cell that was never selected (``not_requested``)
    from one that was eligible but failed (``failed``), so an ineligible record
    is never read as a missing diagnostic.
    """

    return {
        "schema_version": ATTRIBUTION_SCHEMA_VERSION,
        "status": status,
        "reason": reason,
    }


def _transition_records(result: OrientationAttributionResult) -> list[dict[str, Any]]:
    """Transition identity, path lengths, and observed/captured/residual metrics."""

    groups = list(result.config.groups)
    records: list[dict[str, Any]] = []
    for transition in result.transitions:
        observed = transition.observed
        captured = transition.pls_captured
        entry: dict[str, Any] = {
            "transition_id": transition.transition_id,
            "from_stage": transition.from_stage,
            "to_stage": transition.to_stage,
            "groups": groups,
            "components": {},
            "retention": _retention(observed, captured),
        }
        for name in _COMPONENTS:
            component = getattr(transition, name)
            entry["components"][name] = {
                "path_length_group_1": _finite(component.path_lengths[0]),
                "path_length_group_2": _finite(component.path_lengths[1]),
                "contrast_available": bool(component.contrast_available),
                "contrast_norm": (
                    _finite(float(np.linalg.norm(component.directional_contrast)))
                    if component.directional_contrast is not None
                    else None
                ),
            }
        records.append(entry)
    return records


def _retention(observed: Any, captured: Any) -> dict[str, Any]:
    """Observed-versus-PLS-captured retention for one transition.

    Unavailable values stay ``None``; a transition whose group path length is at
    or below the zero tolerance has no direction, and treating that as zero
    retention would fabricate a measurement.
    """

    observed_contrast = observed.directional_contrast
    captured_contrast = captured.directional_contrast
    if observed_contrast is None or captured_contrast is None:
        return {"available": False, "cosine": None, "norm_ratio": None, "residual_fraction": None}
    observed_norm = float(np.linalg.norm(observed_contrast))
    captured_norm = float(np.linalg.norm(captured_contrast))
    if observed_norm <= 0.0:
        return {"available": False, "cosine": None, "norm_ratio": None, "residual_fraction": None}
    cosine = (
        float(np.dot(observed_contrast, captured_contrast) / (observed_norm * captured_norm))
        if captured_norm > 0.0
        else None
    )
    residual_norm = float(np.linalg.norm(observed_contrast - captured_contrast))
    return {
        "available": True,
        "cosine": _finite(cosine) if cosine is not None else None,
        "norm_ratio": _finite(captured_norm / observed_norm),
        "residual_fraction": _finite(residual_norm / observed_norm),
    }


def _top_feature_records(result: OrientationAttributionResult, top_k: int) -> list[dict[str, Any]]:
    """Bounded signed top-k feature records per transition and component."""

    effects = result.feature_effects
    if effects.empty:
        return []
    records: list[dict[str, Any]] = []
    for (transition_id, component), block in effects.groupby(["transition_id", "component"], sort=False):
        ranked = block.assign(_magnitude=block["effect_standardized"].abs())
        ranked = ranked[np.isfinite(ranked["_magnitude"].to_numpy(dtype=float))]
        ranked = ranked.sort_values(["_magnitude", "feature"], ascending=[False, True]).head(top_k)
        for rank, row in enumerate(ranked.itertuples(index=False), start=1):
            records.append(
                {
                    "transition_id": str(transition_id),
                    "component": str(component),
                    "rank": rank,
                    "feature": str(row.feature),
                    "sign": int(np.sign(float(row.effect_standardized))),
                    "effect_standardized": _finite(float(row.effect_standardized)),
                    "effect_original": _finite(row.effect_original),
                    "sign_stability": _finite(row.sign_stability),
                    "top_k_selection_frequency": _finite(row.top_k_selection_frequency),
                }
            )
    return records


def _truth_recovery_records(
    result: OrientationAttributionResult,
    truth: TruthDrivers,
    top_k: int,
) -> list[dict[str, Any]]:
    """Top-k precision, recall, and selected counts against generator truth."""

    effects = result.feature_effects
    if effects.empty:
        return []
    principal_id = result.config.principal_component_id
    records: list[dict[str, Any]] = []
    for (transition_id, component), block in effects.groupby(["transition_id", "component"], sort=False):
        transition_id = str(transition_id)
        # The principal-orientation contrast summarizes the whole configuration,
        # so its truth set is the union over transitions, not any single step.
        mask = truth.global_mask if transition_id == principal_id else truth.per_transition.get(transition_id)
        driver_names = (
            set() if mask is None else {truth.feature_names[index] for index in np.flatnonzero(mask)}
        )
        ranked = block.assign(_magnitude=block["effect_standardized"].abs())
        ranked = ranked[np.isfinite(ranked["_magnitude"].to_numpy(dtype=float))]
        ranked = ranked.sort_values(["_magnitude", "feature"], ascending=[False, True]).head(top_k)
        selected = [str(name) for name in ranked["feature"].tolist()]
        hits = sum(1 for name in selected if name in driver_names)
        available = truth.available and mask is not None and bool(selected)
        records.append(
            {
                "transition_id": transition_id,
                "component": str(component),
                "available": available,
                "selected_count": len(selected),
                "truth_count": len(driver_names),
                "hit_count": hits,
                "precision": (hits / len(selected)) if available else None,
                "recall": (hits / len(driver_names)) if available and driver_names else None,
            }
        )
    return records


def _stability_records(result: OrientationAttributionResult) -> list[dict[str, Any]]:
    """Bootstrap sign and top-k selection stability, aggregated over features."""

    records: list[dict[str, Any]] = []
    for record in result.bootstrap_records:
        records.append(
            {
                "transition_id": record.transition_id,
                "component": record.component,
                "requested_replicates": int(record.requested_replicates),
                "valid_replicates": int(record.valid_replicates),
                "top_k": int(record.top_k),
                "seed": record.seed,
                "mean_sign_stability": _nanmean(record.sign_stability),
                "mean_rank_stability": _nanmean(record.rank_stability),
                "mean_top_k_selection_frequency": _nanmean(record.top_k_selection_frequency),
            }
        )
    return records


def _nanmean(values: np.ndarray) -> float | None:
    array = np.asarray(values, dtype=float)
    if array.size == 0 or not np.isfinite(array).any():
        return None
    return float(np.nanmean(array))


def _finite(value: Any) -> float | None:
    if value is None:
        return None
    numeric = float(value)
    return numeric if np.isfinite(numeric) else None


def flatten_attribution_diagnostics(
    diagnostics: Mapping[str, Any],
) -> dict[str, Any]:
    """Flatten one attribution record's scalar fields for tabular summaries."""

    flat: dict[str, Any] = {
        "attribution_schema_version": diagnostics.get("schema_version"),
        "attribution_status": diagnostics.get("status"),
    }
    settings = diagnostics.get("settings", {})
    flat["attribution_bootstrap_replicates"] = settings.get("bootstrap_replicates")
    flat["attribution_top_k"] = settings.get("top_k")
    model = diagnostics.get("model", {})
    flat["attribution_selected_components"] = model.get("selected_components")
    flat["attribution_feature_order_signature"] = model.get("feature_order_signature")
    runtime = diagnostics.get("runtime", {})
    flat["attribution_seconds"] = runtime.get("attribution_seconds")
    return flat


__all__ = [
    "ATTRIBUTION_SCHEMA_VERSION",
    "UNIT_BASIS",
    "AttributionDiagnosticError",
    "TruthDrivers",
    "compute_attribution_diagnostics",
    "derive_truth_driver_features",
    "run_attribution_diagnostics",
    "flatten_attribution_diagnostics",
    "unavailable_record",
]
