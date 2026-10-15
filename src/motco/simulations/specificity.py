"""Descriptive specificity instrumentation for the feature-surgery modes.

Characterizes how MOTCO *responds* to each realistic trajectory difference —
which statistics move and how much they cross-talk — rather than gating on a
clean diagonal. For each mode it generates several replicates, evaluates them
through the MOTCO trajectory test with RRPP, and reports the per-statistic
rejection rate plus a group-vs-stage projection diagnostic (how much of the
injected group signal lands in the disease/stage-discriminant subspace).

This is a *descriptive* tool, not a pass/fail gate: cross-talk and even
non-detection are findings, not failures. The heavier cluster-run study
produces the definitive specificity matrix and power curves.

The magnitude mode's off-target ``angle``/``shape`` response is **not** the
methylation ``rev.logit`` nonlinearity, as this module previously stated: the
per-omic population geometry of a scaled-delta trajectory is exactly
size-only (``angle`` and ``shape`` are 0 to machine precision at
``population_native``). It is block asymmetry — ``magnitude_kind='all'``
scales ``delta_methyl`` alone while the measurement space standardizes and
concatenates all three omic blocks, so the pooled trajectory rotates purely
because one component grew and the others did not. See
:func:`decompose_block_response`, which reads that decomposition out of
recorded geometry.
"""

from __future__ import annotations

from collections.abc import Sequence
from dataclasses import dataclass
from pathlib import Path
from typing import TYPE_CHECKING, Any

import numpy as np
import pandas as pd

from motco.simulations.evaluation import (
    SimulationEvaluationParams,
    evaluate_semisynthetic_trajectory,
    integrate_semisynthetic_dataset,
)
from motco.simulations.reference import IntersimReference, load_reference
from motco.simulations.semisynthetic import (
    SemiSyntheticTrajectoryParams,
    TrajectoryMode,
    generate_semisynthetic_trajectory,
)

if TYPE_CHECKING:  # imported lazily at runtime to keep the module import cycle-free
    from motco.simulations.grid import SimulationReplicateResult

STATISTICS: tuple[str, ...] = ("delta", "angle", "shape")

#: The statistic each non-null mode is designed to move.
TARGET_STATISTIC: dict[str, str] = {
    "magnitude": "delta",
    "orientation": "angle",
    "shape": "shape",
}


#: Shape is degenerate with a single trajectory step, so the 2-stage isolation
#: pass only exercises the size/orientation modes plus the negative controls.
SHAPE_FREE_MODES: tuple[str, ...] = ("none", "translation", "magnitude", "orientation")


@dataclass(frozen=True)
class ModeSpecificity:
    """Per-mode rejection rates and group-vs-stage projection diagnostic."""

    mode: str
    rejection_rates: dict[str, float]
    group_in_stage_fraction: float
    n_replicates: int
    integration_method: str = "concat"


@dataclass(frozen=True)
class ShapeNullDiagnostic:
    """Observed Procrustes distance vs its RRPP permutation null, per mode.

    Splits the ``shape`` rejection into its two ingredients so the saturation
    can be attributed: ``observed_mean`` is the average observed group-vs-group
    Procrustes distance, and the ``null_*`` summaries describe the spread of the
    permutation null. An anti-conservative test shows a collapsing null (small
    ``null_spread_mean``) rather than an extreme observed distance.
    """

    mode: str
    integration_method: str
    standardize: bool
    observed_mean: float
    null_q025_mean: float
    null_median_mean: float
    null_q975_mean: float
    null_spread_mean: float
    rejection_rate: float
    n_replicates: int


def _group_in_stage_fraction(dataset, params: SimulationEvaluationParams) -> float:
    """Fraction of the group mean-difference energy lying in the stage subspace.

    Projects the concatenated group mean-difference vector onto the span of the
    stage LS-mean directions; ~1 means the injected group signal is disease
    (stage) relevant, ~0 means it is orthogonal to it.
    """

    latent = integrate_semisynthetic_dataset(dataset, params)
    Y = latent.matrix.to_numpy(dtype=float)
    meta = dataset.metadata
    groups = meta[params.group_col].to_numpy()
    stages = meta[params.stage_col].to_numpy()
    g_labels = sorted(set(groups))
    group_diff = Y[groups == g_labels[1]].mean(0) - Y[groups == g_labels[0]].mean(0)
    norm = np.linalg.norm(group_diff)
    if norm < 1e-12:
        return 0.0

    # Stage subspace: centered per-stage means (disease-relevant directions).
    stage_means = np.vstack([Y[stages == s].mean(0) for s in sorted(set(stages))])
    stage_means = stage_means - stage_means.mean(0, keepdims=True)
    basis, _ = np.linalg.qr(stage_means.T)
    projected = basis @ (basis.T @ group_diff)
    return float(np.linalg.norm(projected) / norm)


def evaluate_mode_specificity(
    mode: TrajectoryMode,
    *,
    n_replicates: int = 10,
    n_samples: int = 180,
    n_stages: int = 4,
    effect_size: float = 1.0,
    p_dmp: float = 0.2,
    shape_kind: str = "relocate",
    magnitude_kind: str = "all",
    permutations: int = 99,
    alpha: float = 0.05,
    n_jobs: int | None = -1,
    base_seed: int = 0,
    reference: IntersimReference | None = None,
    integration_method: str = "concat",
    integration_params: dict[str, object] | None = None,
) -> ModeSpecificity:
    """Run replicates for one mode and report per-statistic rejection rates.

    ``integration_method`` selects the latent space the trajectory is measured
    in — the ``concat`` baseline (default), ``snf``, or the ``pls`` production
    latent space — and ``integration_params`` is forwarded to it (e.g. the PLS
    cross-validation knobs). The same selection drives both the RRPP
    rejection-rate evaluation and the group-in-stage projection.
    """

    ref = reference if reference is not None else load_reference()
    int_params = dict(integration_params or {})
    eval_params = SimulationEvaluationParams(
        permutations=permutations,
        n_jobs=n_jobs,
        integration_method=integration_method,  # type: ignore[arg-type]
        integration_params=int_params,
        orientation_estimator="pc1",  # historical diagnostic: PC1 angle
    )
    rejections = {stat: 0 for stat in STATISTICS}
    available = {stat: 0 for stat in STATISTICS}
    gis_values: list[float] = []

    for rep in range(n_replicates):
        params = SemiSyntheticTrajectoryParams(
            seed=base_seed + rep,
            trajectory_mode=mode,
            n_samples=n_samples,
            n_stages=n_stages,
            group_effect_size=effect_size,
            p_dmp=p_dmp,
            shape_kind=shape_kind,  # type: ignore[arg-type]
            magnitude_kind=magnitude_kind,  # type: ignore[arg-type]
        )
        dataset = generate_semisynthetic_trajectory(params, reference=ref)
        result = evaluate_semisynthetic_trajectory(
            dataset,
            SimulationEvaluationParams(
                permutations=permutations,
                n_jobs=n_jobs,
                seed=base_seed + rep,
                integration_method=integration_method,  # type: ignore[arg-type]
                integration_params=int_params,
                orientation_estimator="pc1",  # historical diagnostic: PC1 angle
            ),
        )
        for stat in STATISTICS:
            p = result.p_values.get(stat)
            if p is not None and np.isfinite(p):
                available[stat] += 1
                if p < alpha:
                    rejections[stat] += 1
        gis_values.append(_group_in_stage_fraction(dataset, eval_params))

    rates = {
        stat: (rejections[stat] / available[stat] if available[stat] else float("nan"))
        for stat in STATISTICS
    }
    return ModeSpecificity(
        mode=mode,
        rejection_rates=rates,
        group_in_stage_fraction=float(np.mean(gis_values)),
        n_replicates=n_replicates,
        integration_method=integration_method,
    )


def characterize_two_stage(
    *,
    modes: tuple[str, ...] = SHAPE_FREE_MODES,
    n_replicates: int = 10,
    n_samples: int = 180,
    effect_size: float = 1.0,
    p_dmp: float = 0.2,
    permutations: int = 99,
    alpha: float = 0.05,
    n_jobs: int | None = -1,
    base_seed: int = 0,
    reference: IntersimReference | None = None,
    integration_method: str = "concat",
    integration_params: dict[str, object] | None = None,
) -> dict[str, ModeSpecificity]:
    """Run the shape-free (``n_stages=2``) isolation pass for each mode.

    With two stages the trajectory is a single step, so Procrustes ``shape`` is
    degenerate (the evaluation reports it as ``nan`` and it never counts). This
    isolates ``magnitude``→``delta`` and ``orientation``→``angle`` with shape out
    of the picture — the cleanest test of whether the constructions are specific
    or whether the 3/4-stage cross-talk was shape contaminating them.

    ``integration_method``/``integration_params`` select the latent space, as in
    :func:`evaluate_mode_specificity`.
    """

    ref = reference if reference is not None else load_reference()
    return {
        mode: evaluate_mode_specificity(
            mode,  # type: ignore[arg-type]
            n_replicates=n_replicates,
            n_samples=n_samples,
            n_stages=2,
            effect_size=effect_size,
            p_dmp=p_dmp,
            permutations=permutations,
            alpha=alpha,
            n_jobs=n_jobs,
            base_seed=base_seed,
            reference=ref,
            integration_method=integration_method,
            integration_params=integration_params,
        )
        for mode in modes
    }


def evaluate_shape_null(
    mode: TrajectoryMode,
    *,
    integration_method: str = "concat",
    standardize: bool = True,
    n_replicates: int = 10,
    n_samples: int = 180,
    n_stages: int = 4,
    effect_size: float = 1.0,
    p_dmp: float = 0.2,
    shape_kind: str = "relocate",
    magnitude_kind: str = "all",
    permutations: int = 99,
    alpha: float = 0.05,
    n_jobs: int | None = -1,
    base_seed: int = 0,
    reference: IntersimReference | None = None,
) -> ShapeNullDiagnostic:
    """Split the ``shape`` rejection into observed distance vs permutation null.

    For each replicate, records the observed group-vs-group Procrustes distance
    and the quantiles/spread of its RRPP permutation null (via
    ``include_null_distributions``). ``integration_method``/``standardize`` select
    raw concat (``standardize=False``), concat-standardize (the default), or SNF
    — the probe for whether per-feature standardization is what breaks
    Procrustes scale-invariance.
    """

    ref = reference if reference is not None else load_reference()
    integration_params: dict[str, object] = (
        {} if integration_method == "snf" else {"standardize": standardize}
    )

    observed: list[float] = []
    q025: list[float] = []
    medians: list[float] = []
    q975: list[float] = []
    spreads: list[float] = []
    rejections = 0
    available = 0

    for rep in range(n_replicates):
        params = SemiSyntheticTrajectoryParams(
            seed=base_seed + rep,
            trajectory_mode=mode,
            n_samples=n_samples,
            n_stages=n_stages,
            group_effect_size=effect_size,
            p_dmp=p_dmp,
            shape_kind=shape_kind,  # type: ignore[arg-type]
            magnitude_kind=magnitude_kind,  # type: ignore[arg-type]
        )
        dataset = generate_semisynthetic_trajectory(params, reference=ref)
        result = evaluate_semisynthetic_trajectory(
            dataset,
            SimulationEvaluationParams(
                integration_method=integration_method,  # type: ignore[arg-type]
                integration_params=integration_params,
                permutations=permutations,
                n_jobs=n_jobs,
                seed=base_seed + rep,
                include_null_distributions=True,
                orientation_estimator="pc1",  # historical diagnostic: PC1 angle
            ),
        )
        obs = result.pair_statistics.get("shape")
        null = (result.null_distributions or {}).get("shape")
        if obs is None or not np.isfinite(obs) or not null:
            continue
        null_arr = np.asarray(null, dtype=float)
        observed.append(float(obs))
        q025.append(float(np.quantile(null_arr, 0.025)))
        medians.append(float(np.median(null_arr)))
        q975.append(float(np.quantile(null_arr, 0.975)))
        spreads.append(float(null_arr.std()))
        available += 1
        p = result.p_values.get("shape")
        if p is not None and np.isfinite(p) and p < alpha:
            rejections += 1

    def _mean(values: list[float]) -> float:
        return float(np.mean(values)) if values else float("nan")

    return ShapeNullDiagnostic(
        mode=mode,
        integration_method=integration_method,
        standardize=standardize,
        observed_mean=_mean(observed),
        null_q025_mean=_mean(q025),
        null_median_mean=_mean(medians),
        null_q975_mean=_mean(q975),
        null_spread_mean=_mean(spreads),
        rejection_rate=(rejections / available if available else float("nan")),
        n_replicates=available,
    )


def target_leads(report: ModeSpecificity) -> bool:
    """Descriptive flag: does the mode's target statistic lead the response?

    Informational only — **not** a pass/fail gate. For non-null modes it reports
    whether the target statistic has the highest rejection rate and clears 0.5;
    for the negative controls (``none``/``translation``) it reports whether no
    statistic rejects strongly. Cross-talk is expected; a ``False`` here is a
    finding to characterize, not a failure to fix.
    """

    rates = report.rejection_rates
    target = TARGET_STATISTIC.get(report.mode)
    if target is None:  # none / translation -> negative controls
        return max(rates.values()) <= 0.5
    others = [v for k, v in rates.items() if k != target]
    return rates[target] >= max(others) and rates[target] >= 0.5


#: Scope name for the concatenated (all-omic) measurement, as recorded in
#: ``realized_geometry``. Every other scope is one omic block.
JOINT_SCOPE = "joint"


def _anchor_cell_ids(records: Sequence[SimulationReplicateResult]) -> set[str]:
    """Cell ids serving as the shared zero-effect reference.

    Mirrors the rule used by the study's localization reader: the explicitly
    flagged anchor, plus any ``type_i_*`` cell carrying no trajectory mode.
    """

    return {
        record.cell_id
        for record in records
        if (record.cell_metadata or {}).get("zero_effect_anchor")
        or (record.phase.startswith("type_i_") and not (record.cell_metadata or {}).get("trajectory_mode"))
    }


def decompose_block_response(
    records: Sequence[SimulationReplicateResult] | str | Path,
    *,
    mode: str,
    statistics: Sequence[str] = STATISTICS,
) -> pd.DataFrame:
    """Pooled-versus-per-block realized geometry for one trajectory mode.

    Answers "does this mode's response live in a single omic block, or only in
    the concatenation of the blocks?" — the question that separates a construction
    that is impure *within* an omic from one that is impure only *across* omics.

    ``records`` is either a sequence of replicate results or a path to a merged
    JSONL file. The returned frame carries one row per (checkpoint, scope,
    statistic, effect_size) for ``mode``, with the shared zero-effect anchor's
    value at the same checkpoint/scope/statistic alongside, so an effect-driven
    response is distinguishable from a sampling noise floor. ``is_joint`` marks
    the concatenated scope.

    A checkpoint that records some scopes but not others (``pls_latent`` is
    joint-only) yields rows only for the scopes present: a missing scope is
    absent, never zero-filled, because 0 is a meaningful value here.
    """

    from motco.simulations.grid import read_replicate_results
    from motco.simulations.study.phase4 import summarize_realized_geometry

    resolved: Sequence[SimulationReplicateResult] = (
        read_replicate_results(Path(records))
        if isinstance(records, (str, Path))
        else records
    )
    if not resolved:
        raise ValueError("No replicate records supplied.")

    geometry = summarize_realized_geometry(resolved)
    wanted = [s for s in statistics if s in set(geometry["statistic"])]
    if not wanted:
        raise ValueError(f"None of {tuple(statistics)} present in the recorded geometry.")

    anchors = _anchor_cell_ids(resolved)
    anchor_geometry = geometry[geometry["cell_id"].isin(anchors)]
    anchor_value = {
        (row.checkpoint, row.scope, row.statistic): row.median
        for row in anchor_geometry.itertuples(index=False)
    }

    block = geometry[
        (geometry["trajectory_mode"] == mode)
        & (~geometry["cell_id"].isin(anchors))
        & (geometry["statistic"].isin(wanted))
    ]
    if block.empty:
        raise ValueError(f"No non-anchor records for trajectory mode {mode!r}.")

    rows: list[dict[str, Any]] = []
    for row in block.itertuples(index=False):
        rows.append(
            {
                "trajectory_mode": mode,
                "effect_size": row.effect_size,
                "checkpoint": row.checkpoint,
                "measurement_space": row.measurement_space,
                "scope": row.scope,
                "is_joint": row.scope == JOINT_SCOPE,
                "statistic": row.statistic,
                "median": row.median,
                "mean": row.mean,
                "sd": row.sd,
                "n_available": row.n_available,
                "anchor_median": anchor_value.get((row.checkpoint, row.scope, row.statistic)),
            }
        )

    frame = pd.DataFrame(rows)
    frame["excess_over_anchor"] = frame["median"] - frame["anchor_median"]
    return frame.sort_values(
        ["statistic", "checkpoint", "is_joint", "scope", "effect_size"],
        ignore_index=True,
    )


#: A per-block response at or below this magnitude counts as flat. The
#: analytically-zero cases arrive as floating-point dust from the Procrustes and
#: eigen routines (``shape`` at ``population_native`` lands near 1e-17), so an
#: exact ``== 0`` test would miss precisely the constructions that *are* pure.
#: Well below any real response: recorded ``shape`` values run ~1e-2.
BLOCK_FLAT_ATOL = 1e-12


def summarize_block_localization(
    frame: pd.DataFrame,
    *,
    flat_atol: float = BLOCK_FLAT_ATOL,
) -> pd.DataFrame:
    """Per (checkpoint, statistic): does the response require the concatenation?

    Collapses :func:`decompose_block_response` to one row per checkpoint and
    statistic, reporting the largest per-block response and the joint response
    at the top effect size. ``joint_only`` is True when every individual block
    is flat (within ``flat_atol``) while the joint scope moves — the signature of
    a response created by concatenating blocks that were scaled unequally.

    A checkpoint recording no joint scope (``population_native``) cannot be
    joint-only: there is no concatenated measurement to compare against, so the
    flag stays False rather than being inferred.
    """

    if frame.empty:
        return pd.DataFrame(
            columns=[
                "checkpoint",
                "statistic",
                "top_effect_size",
                "max_block_response",
                "joint_response",
                "joint_only",
            ]
        )

    rows: list[dict[str, Any]] = []
    for (checkpoint, statistic), group in frame.groupby(["checkpoint", "statistic"], sort=True):
        effects = [e for e in group["effect_size"].dropna().unique()]
        top = max(effects) if effects else None
        at_top = group[group["effect_size"] == top] if top is not None else group
        blocks = at_top[~at_top["is_joint"]]
        joint = at_top[at_top["is_joint"]]
        # Per-block response is measured against the anchor so a constant
        # sampling floor does not read as a construction response.
        block_excess = blocks["excess_over_anchor"].abs()
        max_block = float(block_excess.max()) if not block_excess.empty else None
        joint_excess = joint["excess_over_anchor"].abs()
        joint_value = float(joint_excess.max()) if not joint_excess.empty else None
        joint_only = (
            max_block is not None
            and joint_value is not None
            and max_block <= flat_atol
            and joint_value > flat_atol
        )
        rows.append(
            {
                "checkpoint": checkpoint,
                "statistic": statistic,
                "top_effect_size": top,
                "max_block_response": max_block,
                "joint_response": joint_value,
                "joint_only": joint_only,
            }
        )
    return pd.DataFrame(rows)


@dataclass(frozen=True)
class UniformDeltaComparison:
    """Population-standardized geometry of one magnitude construction.

    ``construction`` is the ``magnitude_kind`` value that produced the row:
    ``"all"`` (methylation's delta only, the production default) or ``"joint"``
    (every omic's delta scaled together).
    ``joint_*`` is the concatenated measurement — the space the trajectory
    statistics are actually computed in — and ``max_block_*`` is the largest
    single-omic value, which is 0 to machine precision for a construction that is
    size-pure within each block.
    """

    construction: str
    effect_size: float
    joint_delta: float
    joint_angle: float | None
    joint_shape: float | None
    max_block_angle: float
    max_block_shape: float


def compare_uniform_delta_construction(
    *,
    effect_sizes: Sequence[float] = (0.0, 0.25, 0.5, 0.75, 1.0),
    n_samples: int = 1200,
    n_stages: int = 4,
    p_dmp: float = 0.1,
    baseline_continuity: float = 0.0,
    seed: int = 2,
    magnitude_kinds: Sequence[str] = ("all", "joint"),
    reference: IntersimReference | None = None,
) -> list[UniformDeltaComparison]:
    """Is a size-pure magnitude change realizable after per-block standardization?

    Measures the *population* trajectory geometry at the standardized checkpoint
    for each magnitude construction in ``magnitude_kinds`` (by default the
    production ``all`` and the joint-delta ``joint``), selected through the
    public ``magnitude_kind`` generator parameter. No sampling, no RRPP and no
    PLS fit are involved: the question is whether the concatenated trajectory
    rotates at all, which the analytic population means answer directly.

    An effect size of 0 is the anchor — both constructions reduce to the identity
    there, so its joint angle and shape are the floating-point floor against which
    the nonzero effects are read.
    """

    from motco.simulations.diagnostics import geometry_from_means
    from motco.simulations.preprocessing import (
        OMIC_LAYERS,
        concatenate_blocks,
        fit_omics_preprocessor,
    )

    ref = reference if reference is not None else load_reference()
    out: list[UniformDeltaComparison] = []
    for effect in effect_sizes:
        for construction in magnitude_kinds:
            params = SemiSyntheticTrajectoryParams(
                seed=seed,
                trajectory_mode="magnitude",
                n_samples=n_samples,
                n_stages=n_stages,
                group_effect_size=float(effect),
                p_dmp=p_dmp,
                baseline_continuity=baseline_continuity,
                magnitude_kind=construction,  # type: ignore[arg-type]
            )
            dataset = generate_semisynthetic_trajectory(params, reference=ref)
            population = dataset.population_trajectories
            if population is None:  # pragma: no cover - generator always builds these
                raise ValueError("Generator did not expose population trajectories.")
            groups = sorted(dataset.metadata["group"].astype(str).unique().tolist())
            stages = sorted(dataset.metadata["stage"].astype(str).unique().tolist())
            standardized = fit_omics_preprocessor(dataset).transform_population(population)
            joint = geometry_from_means(concatenate_blocks(standardized), groups, stages)
            blocks = [
                geometry_from_means(standardized[layer], groups, stages)
                for layer in OMIC_LAYERS
            ]
            out.append(
                UniformDeltaComparison(
                    construction=construction,
                    effect_size=float(effect),
                    joint_delta=joint.delta,
                    joint_angle=joint.angle,
                    joint_shape=joint.shape,
                    max_block_angle=max(abs(b.angle or 0.0) for b in blocks),
                    max_block_shape=max(abs(b.shape or 0.0) for b in blocks),
                )
            )
    return out
