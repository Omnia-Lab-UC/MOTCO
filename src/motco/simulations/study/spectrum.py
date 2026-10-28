"""Per-cell latent configuration-spectrum summaries and eigengap stratification.

The relative eigengap ``(l1 - l2) / sum(l)`` of the centered stage-mean
configuration in the latent space is the observable the geometry audit
(``docs/reports/geometry-audit-2026-09-01.md``) found governs whether an
orientation is resolvable: a near-isotropic configuration has no dominant axis
for PC1 to lock onto, its ``angle`` permutation null is wide, and the replicate
cannot reject however large its observed angle.

The continuity-resolved table extends the same idea along the generator's
baseline-continuity axis (``semisynthetic.SemiSyntheticTrajectoryParams.
baseline_continuity``): a baseline whose stage programs persist gives the
configuration a direction to differ in, and the table pairs the resulting
orientation power with the eigengap that is supposed to explain it — so a power
difference along the axis is read off the recorded geometry rather than off the
knob.

Every table here reads **recorded** values only — they never regenerate a
dataset or recompute a spectrum from raw data — and they qualify reporting
rather than changing any statistic or decision.
"""

from __future__ import annotations

import math
from collections.abc import Iterable, Mapping, Sequence
from typing import Any

import numpy as np
import pandas as pd

from motco.simulations.grid import SimulationReplicateResult
from motco.stats.trajectory import pooled_relative_eigengap

#: Modes whose power is stratified by eigengap. Orientation is the mode the
#: eigengap governs; the others move statistics the eigengap does not gate.
STRATIFIED_MODES: tuple[str, ...] = ("orientation",)

#: Phases carrying power cells (as opposed to the Type I controls).
POWER_PHASES: tuple[str, ...] = ("power_primary", "power_ofat", "power_design")

#: Phases whose cells sit at a design point: the primary grid (baseline column,
#: stamped with its coordinates when a design grid is declared) and the
#: design-point grids. OFAT cells vary one factor off the baseline and are not
#: design points.
DESIGN_PHASES: tuple[str, ...] = ("power_primary", "power_design")

#: Cell-metadata key carrying a cell's design-grid coordinates.
DESIGN_POINT_KEY = "design_point"

#: Number of within-cell strata. Terciles: enough to show a monotone trend at
#: 100 replicates per cell without splitting the counts past usefulness.
DEFAULT_N_STRATA = 3

_QUANTILES: tuple[tuple[str, float], ...] = (
    ("q25", 0.25),
    ("median", 0.50),
    ("q75", 0.75),
)


def pooled_eigengap(record: SimulationReplicateResult) -> float | None:
    """Recorded pooled relative eigengap for one replicate, or ``None``.

    ``None`` covers both a record written before the field existed and a
    recorded configuration whose eigengap is undefined; the two are separated by
    :func:`has_spectrum`.
    """

    return pooled_relative_eigengap(record.config_spectrum)


def group_eigengaps(record: SimulationReplicateResult) -> dict[str, float | None]:
    """Recorded per-group relative eigengaps for one replicate, keyed by group."""

    groups = (record.config_spectrum or {}).get("groups") or {}
    return {str(name): _eigengap(entry) for name, entry in dict(groups).items()}


def has_spectrum(record: SimulationReplicateResult) -> bool:
    """Whether the record carries a spectrum block at all."""

    return bool(record.config_spectrum)


def _eigengap(entry: Any) -> float | None:
    if not isinstance(entry, Mapping):
        return None
    value = entry.get("relative_eigengap")
    if value is None:
        return None
    value = float(value)
    return value if math.isfinite(value) else None


def summarize_config_spectrum(
    records: Sequence[SimulationReplicateResult],
) -> pd.DataFrame:
    """Per-cell summary of the recorded pooled and per-group eigengaps.

    One row per (cell, configuration), where configuration is ``pooled`` or a
    group label. ``n_recorded`` counts replicates carrying a spectrum block and
    ``n_available`` those whose eigengap is defined, so records predating the
    field and degenerate configurations stay separately visible.
    """

    buckets: dict[tuple[str, str], list[float | None]] = {}
    identity: dict[tuple[str, str], dict[str, Any]] = {}
    recorded: dict[tuple[str, str], int] = {}

    for record in records:
        if record.status != "completed":
            continue
        meta = dict(record.cell_metadata)
        mode = _resolve_mode(meta, record.phase)
        effect_size = meta.get("effect_size")
        entries: list[tuple[str, float | None]] = [("pooled", pooled_eigengap(record))]
        entries.extend(sorted(group_eigengaps(record).items()))
        if not has_spectrum(record):
            # Still counted, so a cell of legacy records reports zeros rather
            # than disappearing from the table.
            entries = [("pooled", None)]
        for configuration, value in entries:
            key = (record.cell_id, configuration)
            buckets.setdefault(key, []).append(value)
            recorded[key] = recorded.get(key, 0) + int(has_spectrum(record))
            identity.setdefault(
                key,
                {
                    "trajectory_mode": mode,
                    "effect_size": None if effect_size is None else float(effect_size),
                    "cell_id": record.cell_id,
                    "phase": record.phase,
                    "configuration": configuration,
                },
            )

    rows: list[dict[str, Any]] = []
    for key, values in buckets.items():
        available = [value for value in values if value is not None]
        row = dict(identity[key])
        row.update(
            {
                "n_replicates": len(values),
                "n_recorded": recorded[key],
                "n_available": len(available),
                "mean_eigengap": _mean(available),
                "sd_eigengap": _sd(available),
            }
        )
        for name, q in _QUANTILES:
            row[f"{name}_eigengap"] = (
                float(np.quantile(available, q)) if available else None
            )
        row["min_eigengap"] = min(available) if available else None
        row["max_eigengap"] = max(available) if available else None
        rows.append(row)

    frame = pd.DataFrame(
        rows,
        columns=[
            "trajectory_mode",
            "effect_size",
            "cell_id",
            "phase",
            "configuration",
            "n_replicates",
            "n_recorded",
            "n_available",
            "mean_eigengap",
            "sd_eigengap",
            "q25_eigengap",
            "median_eigengap",
            "q75_eigengap",
            "min_eigengap",
            "max_eigengap",
        ],
    )
    if frame.empty:
        return frame
    return frame.sort_values(
        ["trajectory_mode", "effect_size", "cell_id", "configuration"]
    ).reset_index(drop=True)


def stratify_power_by_eigengap(
    records: Sequence[SimulationReplicateResult],
    *,
    alpha: float = 0.05,
    statistics: Sequence[str] = ("delta", "angle", "shape"),
    modes: Sequence[str] = STRATIFIED_MODES,
    n_strata: int = DEFAULT_N_STRATA,
) -> pd.DataFrame:
    """Rejection rates within eigengap strata, for the stratified power cells.

    Strata are within-cell quantile bins of the **recorded** pooled eigengap, so
    every replicate in a cell is compared against its own cell's distribution.
    A cell whose records carry no spectrum yields one row per statistic marked
    ``status='unavailable'`` rather than silently vanishing — the covariate is
    missing, which is a different fact from a degenerate geometry.
    """

    if not (0 < alpha < 1):
        raise ValueError("alpha must be between 0 and 1.")
    if n_strata < 2:
        raise ValueError("n_strata must be at least 2.")

    rows: list[dict[str, Any]] = []
    for cell_id, group in _stratified_cells(records, modes).items():
        meta = dict(group[0].cell_metadata)
        effect_size = meta.get("effect_size")
        identity = {
            "trajectory_mode": _resolve_mode(meta, group[0].phase),
            "effect_size": None if effect_size is None else float(effect_size),
            "cell_id": cell_id,
            "phase": group[0].phase,
        }
        gaps = [(record, pooled_eigengap(record)) for record in group]
        usable = [(record, gap) for record, gap in gaps if gap is not None]
        if not usable:
            for statistic in statistics:
                rows.append(
                    {
                        **identity,
                        "statistic": statistic,
                        "stratum": None,
                        "n_strata": n_strata,
                        "eigengap_low": None,
                        "eigengap_high": None,
                        "mean_eigengap": None,
                        "n_replicates": len(group),
                        "n_available": 0,
                        "n_rejected": 0,
                        "rejection_rate": None,
                        "monte_carlo_se": None,
                        "status": "unavailable",
                    }
                )
            continue

        assignments = _quantile_strata([gap for _, gap in usable], n_strata)
        for stratum in range(n_strata):
            members = [
                (record, gap)
                for (record, gap), index in zip(usable, assignments)
                if index == stratum
            ]
            stratum_gaps = [gap for _, gap in members]
            for statistic in statistics:
                p_values = [
                    float(p)
                    for record, _ in members
                    if (p := record.p_values.get(statistic)) is not None
                    and math.isfinite(float(p))
                ]
                rejected = sum(p < alpha for p in p_values)
                rate = rejected / len(p_values) if p_values else None
                rows.append(
                    {
                        **identity,
                        "statistic": statistic,
                        "stratum": stratum,
                        "n_strata": n_strata,
                        "eigengap_low": min(stratum_gaps) if stratum_gaps else None,
                        "eigengap_high": max(stratum_gaps) if stratum_gaps else None,
                        "mean_eigengap": _mean(stratum_gaps),
                        "n_replicates": len(members),
                        "n_available": len(p_values),
                        "n_rejected": rejected,
                        "rejection_rate": rate,
                        "monte_carlo_se": (
                            math.sqrt(rate * (1.0 - rate) / len(p_values))
                            if rate is not None
                            else None
                        ),
                        "status": "ok",
                    }
                )

    frame = pd.DataFrame(
        rows,
        columns=[
            "trajectory_mode",
            "effect_size",
            "cell_id",
            "phase",
            "statistic",
            "stratum",
            "n_strata",
            "eigengap_low",
            "eigengap_high",
            "mean_eigengap",
            "n_replicates",
            "n_available",
            "n_rejected",
            "rejection_rate",
            "monte_carlo_se",
            "status",
        ],
    )
    if frame.empty:
        return frame
    return frame.sort_values(
        ["trajectory_mode", "effect_size", "cell_id", "statistic", "stratum"],
        na_position="last",
    ).reset_index(drop=True)


#: Truth-metadata key recording the baseline stage-program continuity a record
#: was generated under. Absent from records written before the axis existed.
CONTINUITY_KEY = "baseline_continuity"

#: Terciles of the recorded eigengap, reported per continuity value.
_CONTINUITY_QUANTILES: tuple[tuple[str, float], ...] = (
    ("q33", 1.0 / 3.0),
    ("median", 0.50),
    ("q67", 2.0 / 3.0),
)


def record_continuity(record: SimulationReplicateResult) -> float | None:
    """Baseline continuity the record was generated under, or ``None``.

    ``None`` covers records written before the continuity axis existed; they are
    a different fact from a recorded ρ of 0 and are excluded from the
    continuity-resolved view rather than folded into the zero bin.
    """

    value = (record.truth_metadata or {}).get(CONTINUITY_KEY)
    if value is None:
        return None
    value = float(value)
    return value if math.isfinite(value) else None


def resolve_orientation_by_continuity(
    records: Sequence[SimulationReplicateResult],
    *,
    alpha: float = 0.05,
    statistics: Sequence[str] = ("delta", "angle", "shape"),
    modes: Sequence[str] = STRATIFIED_MODES,
) -> pd.DataFrame:
    """Orientation operating characteristics resolved by baseline continuity.

    One row per (continuity value, mode, effect size, statistic), pooling the
    cells that share those coordinates. Alongside the rejection rate each row
    carries the recorded pooled-eigengap distribution and the dispersion of the
    per-replicate ``angle`` permutation-null width, because the *eigengap* — not
    the knob — is the observable expected to carry a continuity-conditioned
    orientation claim to real data: power that rises along ρ should be traceable
    to a configuration that acquired a dominant axis.

    Returns an **empty frame** unless the record set spans at least two distinct
    continuity values, so a study that holds the axis fixed reports exactly what
    it reported before the axis existed.
    """

    if not (0 < alpha < 1):
        raise ValueError("alpha must be between 0 and 1.")

    wanted = set(modes)
    # Records at a design point are additionally keyed on every design
    # coordinate other than continuity, so a grid crossing ρ with (say) sample
    # size never pools three operating points into one row. Without a design
    # grid the extra key is empty and the frame is exactly the pre-grid one.
    buckets: dict[
        tuple[float, tuple[tuple[str, Any], ...], str, float | None],
        list[SimulationReplicateResult],
    ] = {}
    extra_axes: list[str] = []
    for record in records:
        if record.status != "completed" or record.phase not in POWER_PHASES:
            continue
        continuity = record_continuity(record)
        if continuity is None:
            continue
        meta = dict(record.cell_metadata)
        mode = _resolve_mode(meta, record.phase)
        if mode not in wanted:
            continue
        effect = meta.get("effect_size")
        point = record_design_point(record) or {}
        others = tuple(
            (axis, value) for axis, value in point.items() if axis != CONTINUITY_AXIS
        )
        for axis, _ in others:
            if axis not in extra_axes:
                extra_axes.append(axis)
        key = (continuity, others, mode, None if effect is None else float(effect))
        buckets.setdefault(key, []).append(record)

    columns = list(_CONTINUITY_COLUMNS[:1]) + extra_axes + list(_CONTINUITY_COLUMNS[1:])
    if len({key[0] for key in buckets}) < 2:
        return pd.DataFrame(columns=columns)

    rows: list[dict[str, Any]] = []
    for (continuity, others, mode, effect), group in sorted(
        buckets.items(), key=_continuity_sort_key
    ):
        gaps = [gap for gap in (pooled_eigengap(record) for record in group) if gap is not None]
        widths = [
            width
            for width in (_angle_null_width(record) for record in group)
            if width is not None
        ]
        identity = {
            "baseline_continuity": continuity,
            **{axis: None for axis in extra_axes},
            **dict(others),
            "trajectory_mode": mode,
            "effect_size": effect,
            "n_cells": len({record.cell_id for record in group}),
            "n_replicates": len(group),
        }
        geometry = {
            "n_eigengap_available": len(gaps),
            "mean_eigengap": _mean(gaps),
            **{
                f"{name}_eigengap": (float(np.quantile(gaps, q)) if gaps else None)
                for name, q in _CONTINUITY_QUANTILES
            },
            "n_angle_null_available": len(widths),
            "median_angle_null_q95": (
                float(np.quantile(widths, 0.5)) if widths else None
            ),
            "iqr_angle_null_q95": (
                float(np.quantile(widths, 0.75) - np.quantile(widths, 0.25))
                if widths
                else None
            ),
            "sd_angle_null_q95": _sd(widths),
        }
        for statistic in statistics:
            p_values = [
                float(p)
                for record in group
                if (p := record.p_values.get(statistic)) is not None
                and math.isfinite(float(p))
            ]
            rejected = sum(p < alpha for p in p_values)
            rate = rejected / len(p_values) if p_values else None
            rows.append(
                {
                    **identity,
                    "statistic": statistic,
                    "n_available": len(p_values),
                    "n_rejected": rejected,
                    "rejection_rate": rate,
                    "monte_carlo_se": (
                        math.sqrt(rate * (1.0 - rate) / len(p_values))
                        if rate is not None
                        else None
                    ),
                    **geometry,
                }
            )

    return pd.DataFrame(rows, columns=columns)


#: Design-grid axis name of the baseline-continuity generator parameter.
CONTINUITY_AXIS = f"generator.{CONTINUITY_KEY}"


_CONTINUITY_COLUMNS: tuple[str, ...] = (
    "baseline_continuity",
    "trajectory_mode",
    "effect_size",
    "statistic",
    "n_cells",
    "n_replicates",
    "n_available",
    "n_rejected",
    "rejection_rate",
    "monte_carlo_se",
    "n_eigengap_available",
    "mean_eigengap",
    "q33_eigengap",
    "median_eigengap",
    "q67_eigengap",
    "n_angle_null_available",
    "median_angle_null_q95",
    "iqr_angle_null_q95",
    "sd_angle_null_q95",
)


def _continuity_sort_key(
    item: tuple[tuple[float, tuple[tuple[str, Any], ...], str, float | None], Any],
) -> tuple[Any, ...]:
    continuity, others, mode, effect = item[0]
    return (continuity, _sortable(others), mode, math.inf if effect is None else effect)


def _sortable(values: Any) -> Any:
    """Make design coordinates orderable across mixed value types."""

    if isinstance(values, tuple):
        # Tagged like the scalars, so a sequence-valued coordinate (a size table)
        # and a scalar or ``None`` on the same axis still compare.
        return (2, tuple(_sortable(value) for value in values))
    if isinstance(values, bool | int | float):
        return (0, float(values))
    return (1, str(values))


def record_design_point(record: SimulationReplicateResult) -> dict[str, Any] | None:
    """Design-grid coordinates the record's cell sits at, or ``None``.

    Only cells of :data:`DESIGN_PHASES` are design points; OFAT cells vary one
    factor off the baseline and return ``None`` even if a grid was declared.
    """

    if record.phase not in DESIGN_PHASES:
        return None
    point = (record.cell_metadata or {}).get(DESIGN_POINT_KEY)
    if not isinstance(point, Mapping):
        return None
    return {axis: _freeze(value) for axis, value in point.items()}


def _freeze(value: Any) -> Any:
    """Sequence-valued coordinates as tuples, so a design point can key a dict.

    JSON persistence turns a tuple-valued axis (``generator.group_stage_sizes``,
    ``evaluation.integration_params.layers``) into nested lists.
    """

    if isinstance(value, list | tuple):
        return tuple(_freeze(item) for item in value)
    return value


#: Design-grid axis name of the retained PLS rank (nested evaluation parameter).
#: ``null`` on this axis is the cross-validated reference column.
RANK_AXIS = "evaluation.integration_params.forced_components"

#: Component-selection modes an integration can record.
COMPONENT_SELECTION_MODES: tuple[str, ...] = ("cv", "forced")


def record_component_selection(record: SimulationReplicateResult) -> str | None:
    """Recorded component-selection mode (``cv``/``forced``), or ``None`` when absent.

    Records written before the marker existed carry no field; they are a
    different fact from a recorded ``cv`` and are reported as missing rather
    than assumed.
    """

    value = (record.integration_metadata or {}).get("component_selection")
    return None if value is None else str(value)


def record_forced_rank(record: SimulationReplicateResult) -> int | None:
    """Forced retained rank declared at the record's design point, or ``None``.

    Reads the design coordinates rather than the integration metadata, so a
    forced rank is recognized only when the study *declared* it as a design
    axis — the condition under which a forced-rank record is a legitimate
    study measurement rather than a stray diagnostic.
    """

    point = record_design_point(record) or {}
    value = point.get(RANK_AXIS)
    return None if value is None else int(value)


def _selected_components(record: SimulationReplicateResult) -> int | None:
    """Recorded selected latent dimensionality, or ``None`` when unavailable."""

    integration = record.integration_metadata or {}
    selected = integration.get("selected_lv")
    if selected is None:
        selected = (integration.get("integration_params") or {}).get("selected_lv")
    return None if selected is None else int(selected)


def resolve_operating_by_design_point(
    records: Sequence[SimulationReplicateResult],
    *,
    alpha: float = 0.05,
    statistics: Sequence[str] = ("delta", "angle", "shape"),
) -> pd.DataFrame:
    """Operating characteristics resolved on every design-grid coordinate.

    One row per (design point, mode, effect size, statistic), the baseline
    column and each point's zero-effect anchor included (mode ``none`` at effect
    ``0.0``). Beside the rejection rate each row carries the recorded
    pooled-eigengap distribution, the dispersion of the per-replicate ``angle``
    null width, and the distribution of the selected latent dimensionality —
    the covariates the design-point decision is read against.

    Returns an **empty frame** unless at least one completed ``power_design``
    record is present, so a study without a design grid reports nothing new.
    """

    if not (0 < alpha < 1):
        raise ValueError("alpha must be between 0 and 1.")

    buckets, axes, saw_design_phase = _design_point_buckets(records)
    columns = axes + list(_DESIGN_POINT_COLUMNS)
    if not saw_design_phase:
        return pd.DataFrame(columns=columns)

    rows: list[dict[str, Any]] = []
    for (coords, mode, effect), group in _sorted_design_buckets(buckets):
        gaps = [gap for gap in (pooled_eigengap(record) for record in group) if gap is not None]
        widths = [
            width for width in (_angle_null_width(record) for record in group) if width is not None
        ]
        selected = [
            value for value in (_selected_components(record) for record in group) if value is not None
        ]
        identity = _design_identity(axes, coords, mode, effect, group)
        geometry = {
            "n_eigengap_available": len(gaps),
            "mean_eigengap": _mean(gaps),
            **{
                f"{name}_eigengap": (float(np.quantile(gaps, q)) if gaps else None)
                for name, q in _CONTINUITY_QUANTILES
            },
            "n_angle_null_available": len(widths),
            "median_angle_null_q95": (float(np.quantile(widths, 0.5)) if widths else None),
            "iqr_angle_null_q95": (
                float(np.quantile(widths, 0.75) - np.quantile(widths, 0.25)) if widths else None
            ),
            "sd_angle_null_q95": _sd(widths),
            "n_selected_lv_available": len(selected),
            "median_selected_lv": (float(np.median(selected)) if selected else None),
            "min_selected_lv": (min(selected) if selected else None),
            "max_selected_lv": (max(selected) if selected else None),
            "component_selection": _component_selection_of(group),
        }
        for statistic in statistics:
            p_values = [
                float(p)
                for record in group
                if (p := record.p_values.get(statistic)) is not None and math.isfinite(float(p))
            ]
            rejected = sum(p < alpha for p in p_values)
            rate = rejected / len(p_values) if p_values else None
            rows.append(
                {
                    **identity,
                    "statistic": statistic,
                    "n_available": len(p_values),
                    "n_rejected": rejected,
                    "rejection_rate": rate,
                    "monte_carlo_se": (
                        math.sqrt(rate * (1.0 - rate) / len(p_values)) if rate is not None else None
                    ),
                    **geometry,
                }
            )
    return pd.DataFrame(rows, columns=columns)


#: Bucket key of the design-point tables: the design coordinates, the mode and
#: the effect size.
_DesignKey = tuple[tuple[tuple[str, Any], ...], str, float | None]


def _design_point_buckets(
    records: Sequence[SimulationReplicateResult],
) -> tuple[dict[_DesignKey, list[SimulationReplicateResult]], list[str], bool]:
    """Completed design-point records bucketed by (design point, mode, effect).

    Shared by every design-point table so their rows line up: the axes are
    listed in order of first appearance, and the third value says whether any
    ``power_design`` record was seen — without one the study declared no design
    grid and the tables are empty.
    """

    buckets: dict[_DesignKey, list[SimulationReplicateResult]] = {}
    axes: list[str] = []
    saw_design_phase = False
    for record in records:
        if record.status != "completed":
            continue
        point = record_design_point(record)
        if point is None:
            continue
        saw_design_phase = saw_design_phase or record.phase == "power_design"
        for axis in point:
            if axis not in axes:
                axes.append(axis)
        meta = dict(record.cell_metadata)
        mode = _resolve_mode(meta, record.phase)
        effect = meta.get("effect_size")
        key = (tuple(point.items()), mode, None if effect is None else float(effect))
        buckets.setdefault(key, []).append(record)
    return buckets, axes, saw_design_phase


def _sorted_design_buckets(
    buckets: Mapping[_DesignKey, list[SimulationReplicateResult]],
) -> list[tuple[_DesignKey, list[SimulationReplicateResult]]]:
    return sorted(
        buckets.items(),
        key=lambda item: (_sortable(item[0][0]), item[0][1], math.inf if item[0][2] is None else item[0][2]),
    )


def _design_identity(
    axes: Sequence[str],
    coords: tuple[tuple[str, Any], ...],
    mode: str,
    effect: float | None,
    group: Sequence[SimulationReplicateResult],
) -> dict[str, Any]:
    """Identity columns shared by the design-point tables."""

    phases = sorted({record.phase for record in group})
    return {
        **{axis: None for axis in axes},
        **dict(coords),
        "phase": phases[0] if len(phases) == 1 else "+".join(phases),
        "is_baseline": all(record.phase == "power_primary" for record in group),
        "trajectory_mode": mode,
        "effect_size": effect,
        "n_cells": len({record.cell_id for record in group}),
        "n_replicates": len(group),
    }


def record_cv_auroc(record: SimulationReplicateResult) -> float | None:
    """Recorded held-out stage AUROC of the PLS double CV, or ``None``.

    ``None`` covers a record without the field (``concat``, ``pca``, records
    written before it existed), a forced rank (which runs no CV and records
    ``None``) and a non-finite value; all are reported as unavailable rather
    than as a low AUROC.
    """

    value = (record.integration_metadata or {}).get("cv_mean_auroc")
    if value is None:
        return None
    value = float(value)
    return value if math.isfinite(value) else None


def resolve_stage_signal_by_design_point(
    records: Sequence[SimulationReplicateResult],
) -> pd.DataFrame:
    """Stage separability resolved on every design-grid coordinate.

    One row per (design point, mode, effect size), bucketed exactly as
    :func:`resolve_operating_by_design_point` so the rows line up with the
    operating table, the baseline column and every anchor included. Each row
    carries the number of records with a recorded held-out CV stage AUROC and
    that AUROC's mean and terciles — the stage-signal number that compares
    directly with a real analysis's recorded ``cv_mean_auroc``. A column whose
    records carry no CV AUROC (a forced rank, ``concat``, ``pca``) reports zero
    available and empty summaries.

    Returns an **empty frame** unless at least one completed ``power_design``
    record is present, so a study without a design grid reports nothing new.
    """

    buckets, axes, saw_design_phase = _design_point_buckets(records)
    columns = axes + list(_STAGE_SIGNAL_COLUMNS)
    if not saw_design_phase:
        return pd.DataFrame(columns=columns)

    rows: list[dict[str, Any]] = []
    for (coords, mode, effect), group in _sorted_design_buckets(buckets):
        aurocs = [value for value in (record_cv_auroc(record) for record in group) if value is not None]
        row = _design_identity(axes, coords, mode, effect, group)
        row.update(
            {
                "n_auroc_available": len(aurocs),
                "mean_cv_auroc": _mean(aurocs),
                **{
                    f"{name}_cv_auroc": (float(np.quantile(aurocs, q)) if aurocs else None)
                    for name, q in _CONTINUITY_QUANTILES
                },
            }
        )
        rows.append(row)
    return pd.DataFrame(rows, columns=columns)


_STAGE_SIGNAL_COLUMNS: tuple[str, ...] = (
    "phase",
    "is_baseline",
    "trajectory_mode",
    "effect_size",
    "n_cells",
    "n_replicates",
    "n_auroc_available",
    "mean_cv_auroc",
    "q33_cv_auroc",
    "median_cv_auroc",
    "q67_cv_auroc",
)


def _component_selection_of(group: Sequence[SimulationReplicateResult]) -> str | None:
    """``cv``/``forced`` when every record agrees, ``mixed`` otherwise, ``None`` if unrecorded.

    ``mixed`` cannot arise from a well-formed study — one design point has one
    evaluation configuration — so its appearance flags a configuration error
    rather than being averaged away.
    """

    modes = {mode for mode in (record_component_selection(record) for record in group) if mode is not None}
    if not modes:
        return None
    if len(modes) == 1:
        return next(iter(modes))
    return "mixed"


_DESIGN_POINT_COLUMNS: tuple[str, ...] = (
    "phase",
    "is_baseline",
    "trajectory_mode",
    "effect_size",
    "statistic",
    "n_cells",
    "n_replicates",
    "n_available",
    "n_rejected",
    "rejection_rate",
    "monte_carlo_se",
    "n_eigengap_available",
    "mean_eigengap",
    "q33_eigengap",
    "median_eigengap",
    "q67_eigengap",
    "n_angle_null_available",
    "median_angle_null_q95",
    "iqr_angle_null_q95",
    "sd_angle_null_q95",
    "n_selected_lv_available",
    "median_selected_lv",
    "min_selected_lv",
    "max_selected_lv",
    "component_selection",
)


def _angle_null_width(record: SimulationReplicateResult) -> float | None:
    """Recorded 95th-percentile width of this replicate's ``angle`` null."""

    entry = (record.null_summary or {}).get("angle") or {}
    value = entry.get("q95")
    if value is None:
        return None
    value = float(value)
    return value if math.isfinite(value) else None


def _stratified_cells(
    records: Iterable[SimulationReplicateResult],
    modes: Sequence[str],
) -> dict[str, list[SimulationReplicateResult]]:
    wanted = set(modes)
    groups: dict[str, list[SimulationReplicateResult]] = {}
    for record in records:
        if record.status != "completed" or record.phase not in POWER_PHASES:
            continue
        meta = dict(record.cell_metadata)
        if _resolve_mode(meta, record.phase) not in wanted:
            continue
        groups.setdefault(record.cell_id, []).append(record)
    return dict(sorted(groups.items()))


def _quantile_strata(values: Sequence[float], n_strata: int) -> list[int]:
    """Assign each value to a within-cell quantile stratum.

    Ranks are used rather than raw cut points so ties cannot pile every
    replicate into one bin, and the strata stay as close to equal-sized as the
    data allow.
    """

    order = np.argsort(np.asarray(values, dtype=float), kind="stable")
    ranks = np.empty(len(values), dtype=int)
    ranks[order] = np.arange(len(values))
    return [int(rank * n_strata // max(len(values), 1)) for rank in ranks]


def _resolve_mode(meta: Mapping[str, Any], phase: str) -> str:
    mode = meta.get("trajectory_mode")
    if mode is not None:
        return str(mode)
    return "none" if phase.startswith("type_i") else "unknown"


def _mean(values: Sequence[float]) -> float | None:
    return float(np.mean(values)) if values else None


def _sd(values: Sequence[float]) -> float | None:
    return float(np.std(values, ddof=1)) if len(values) > 1 else None


__all__ = [
    "COMPONENT_SELECTION_MODES",
    "CONTINUITY_AXIS",
    "CONTINUITY_KEY",
    "DEFAULT_N_STRATA",
    "DESIGN_PHASES",
    "DESIGN_POINT_KEY",
    "POWER_PHASES",
    "RANK_AXIS",
    "STRATIFIED_MODES",
    "group_eigengaps",
    "has_spectrum",
    "pooled_eigengap",
    "record_component_selection",
    "record_continuity",
    "record_cv_auroc",
    "record_design_point",
    "record_forced_rank",
    "resolve_operating_by_design_point",
    "resolve_orientation_by_continuity",
    "resolve_stage_signal_by_design_point",
    "stratify_power_by_eigengap",
    "summarize_config_spectrum",
]
