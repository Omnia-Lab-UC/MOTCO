"""Semi-synthetic trajectory datasets from the numpy generator.

Trajectory differences between two groups are defined as **feature-set surgery
on the methylation differential indicators only** — gene-expression and protein
differential indicators are always *re-derived* from the (group-specific)
methylation indicators through the cached CpG→gene→protein incidence maps. This
honours the biological cascade (methylation drives expression drives protein)
and keeps the datasets realistic rather than tailored to MOTCO: we manipulate
only the original methylation features and never the latent space.

Group A inherits a random baseline trajectory whose per-stage methylation
indicators are drawn along a declared *continuity axis*
(``baseline_continuity``, ρ ∈ [0, 1)): each CpG's differential status follows a
stationary first-order Markov chain across stages, so the per-stage Bernoulli
(``p_dmp``) marginal is preserved at every ρ while the cross-stage correlation
is ρ^|t−s|. ρ = 0 is the independent, deliberately *not* forced continuous
endpoint — the isotropic stress test, byte-identical to the pre-axis generator —
and larger ρ makes stage means *trend*: pairwise distances grow with stage
separation, giving the configuration a dominant PC1 (a direction to differ in)
instead of a near-regular simplex with a near-zero eigengap. Group B is a
deterministic transform of A's **methylation** indicators:

- ``none``        -- identical indicators (null).
- ``translation`` -- A's stage-changing sites unchanged, plus an extra set ``U``
  of methylation sites (disjoint from the stage-changing sites) differential at
  *every* B stage and at none of A's. A constant group offset → moves only the
  (untested) group main effect, not size/orientation/shape.
- ``magnitude``   -- scaled effect size. ``magnitude_kind='joint'`` (the
  default) keeps the indicators and scales **every** omic's δ by the same
  ``1 + e`` → B's native-space trajectory is exactly ``1 + e`` times A's in
  every block, a size change that is pure in the joint standardized space as
  well as within each block. ``magnitude_kind='all'`` keeps the indicators and
  scales the global methylation δ alone (``δ_methyl_B = (1 + e)·δ_methyl``) →
  uniformly enlarges every methylation step but rotates the concatenated
  trajectory; it was the default before the Phase 5 magnitude re-measurement
  and is pinned explicitly by the historical study configs.
  ``magnitude_kind='extremes'`` leaves δ and scales A's methylation indicators
  at the first and last stages only → a size change localized to the
  endpoints. ``extremes`` (like ``shape_kind='magnitude'`` below) remains a
  methylation-only scaling; the joint default does not extend to it.
- ``orientation`` -- relocate a fraction ``e`` of the stage-changing sites to
  different CpGs, the **same relocation at every stage** → the per-stage pattern
  runs along different feature axes (a rotation).
- ``shape``       -- perturb a **single interior stage** (≥3 stages): relocate a
  fraction ``e`` of that stage's sites (``shape_kind='relocate'``) or scale that
  stage's methylation effect (``shape_kind='magnitude'``) → bends one vertex.

``group_effect_size`` (``e``) is the unified knob (``e = 0`` is null for every
mode). Cross-talk between statistics (e.g. magnitude bending shape via the
methylation ``rev.logit`` nonlinearity) is *expected and reported*, not
engineered away — how well MOTCO separates the modes is an open question the
study characterizes. Generation runs entirely on the numpy generator and cached
reference data — no R at runtime.

**Pool-limited surgeries and censoring.** ``orientation``, ``translation``, and
``shape`` with ``shape_kind='relocate'`` each draw their surgery from a finite
destination/candidate pool, so a large ``e`` can request more sites than the
pool holds. ``surgery_censoring`` makes that explicit: ``"error"`` (the default)
raises rather than realize a partial surgery, ``"clamp"`` realizes as much as
the pool allows. The pools depend on the per-replicate baseline draw, so whether
the limit binds is a property of the dataset, not of the config alone — which is
why the default fails loudly and truth metadata always records the nominal size,
the realized size, and a ``censored`` flag. Silent clamping made distinct
requested effects produce identical datasets and was reported as independent
evidence; see ``docs/reports/geometry-audit-2026-09-01.md`` finding F2.

**Stage-span construction (opt-in).** ``construction='stage_span'`` replaces
the surgeries above with a measurement-defined sensitivity construction. It
leaves both groups' indicators and effects at the baseline and instead edits
group B's exact null population stage means inside the span of B's own stage
configuration. The edit is made over ``span_layers`` (default: all three
blocks) in analysed units (methylation M-values), each feature divided by its
null population pooled SD — the within-cell variance (the reference covariance
diagonal) plus the between-cell variance of the null means at the realized cell
proportions, the population analogue of the pooled sample scaler. With ``z_s``
B's standardized stage means, centroid ``c``, path length ``L``, mean step
``L̄ = L / (S − 1)`` and ``u`` the unit net displacement:

- ``magnitude``   -- ``c + (1 + e)(z_s − c)``: size ×(1 + e), as the joint surgery.
- ``orientation`` -- a rigid rotation about ``c`` by ``θ = 90°·e`` (``e ≤ 2``)
  in the plane of ``u`` and ``v``, the leading principal direction of the
  configuration once its ``u`` components are removed: net-displacement
  ``angle`` = θ, ``delta`` = ``shape`` = 0.
- ``shape``       -- the stage the ``shape`` surgery perturbs moves ``e·L̄`` off
  the chord between its neighbours (inside the span), then the configuration is
  rescaled to path length ``L`` and recentred on ``c``: ``delta`` = ``angle`` = 0.
- ``translation`` -- every stage moves ``e·L̄`` along ``u``: a constant offset
  along the progression, the hardest in-span negative control.

The plant is mapped back to analysed units and added to group B's cell means
before sampling (``generate_omics(..., cell_shift_*=)``); blocks outside
``span_layers`` are untouched. The construction draws no random numbers, so a
stage-span dataset is the ``none`` dataset at the same seed plus the plant. It
is circular by design — it plants where the stage-supervised measurement looks
— and answers how sensitive the test is to a difference it can see, not how
visible the biologically grounded surgeries are. ``orientation`` and ``shape``
need a second in-span direction: fewer than three stages, or a configuration
whose second singular value is at most ``1e-8`` times its first, is an error.
Truth records the construction, the span layers, the effect parameters, the
intended population geometry in the null-standardized coordinates and the
plant's norm per stage; the population trajectories include the plant.
``construction='surgery'`` normalizes to ``None``, the default, so every
surgery dataset, cell id and signature is unchanged.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any, Literal, Mapping, Sequence

import numpy as np
import pandas as pd

from motco.simulations.generator import (
    GeneratedOmics,
    derive_coupled_indicators,
    generate_omics,
    markov_indicators,
    omic_population_means,
)
from motco.simulations.reference import IntersimReference, load_reference
from motco.stats.trajectory import estimate_difference

OmicsLayer = Literal["methylation", "expression", "proteomics"]
TrajectoryMode = Literal["none", "translation", "magnitude", "orientation", "shape"]
ShapeKind = Literal["relocate", "magnitude"]
MagnitudeKind = Literal["all", "extremes", "joint"]
SurgeryCensoring = Literal["error", "clamp"]
Construction = Literal["surgery", "stage_span"]

_OMICS_LAYERS: tuple[OmicsLayer, ...] = ("methylation", "expression", "proteomics")
_MODES = frozenset({"none", "translation", "magnitude", "orientation", "shape"})
_SHAPE_KINDS = frozenset({"relocate", "magnitude"})
_MAGNITUDE_KINDS = frozenset({"all", "extremes", "joint"})
_SURGERY_CENSORINGS = frozenset({"error", "clamp"})
_CONSTRUCTIONS = frozenset({"surgery", "stage_span"})

#: Largest stage-span orientation effect: ``e = 2`` rotates the configuration by 180°.
STAGE_SPAN_MAX_ORIENTATION_EFFECT = 2.0

#: A stage-span configuration whose second singular value is at most this
#: fraction of its first has no second in-span direction.
STAGE_SPAN_RANK_TOLERANCE = 1e-8

#: A pooled population SD below this keeps unit scale, as the sample scaler does.
_SCALE_TOLERANCE = 1e-10


class SemiSyntheticTrajectoryError(ValueError):
    """Raised when semi-synthetic trajectory generation parameters are invalid."""


@dataclass(frozen=True)
class SemiSyntheticTrajectoryParams:
    """Parameters for semi-synthetic trajectory generation.

    ``group_ratio`` is the proportion assigned to the first label in
    ``group_labels`` within every stage. ``group_effect_size`` is the unified
    effect knob described in the module docstring (``0`` is null for all modes).
    ``p_dmp`` is the per-stage probability that a methylation feature is
    differential (InterSIM's ``p.DMP``); expression/protein indicators are
    derived from it via the cross-omic maps. ``baseline_continuity`` (ρ) is the
    baseline stage-program continuity described in the module docstring: 0 draws
    each stage's program independently, higher values make a CpG's differential
    status persist across adjacent stages with correlation ρ while holding every
    stage's marginal at ``p_dmp``. ``delta_*`` are the per-omic
    mean-shift sizes (InterSIM's ``delta.*``). ``shape_kind`` selects the
    single-interior-stage perturbation used by ``shape``. ``magnitude_kind``
    selects whether ``magnitude`` scales every omic's effect size together
    (``joint``, the default), group B's methylation effect alone at *all*
    stages (``all``, a uniform methylation δ scale), or the methylation effect
    only at the *extreme* stages (``extremes``, first and last) leaving interior
    stages at the baseline effect.
    ``group_stage_sizes`` is the explicit alternative to ``n_samples`` /
    ``stage_sample_prop`` / ``group_ratio``: one tuple of per-stage counts per
    group, rows ordered as ``group_labels`` and columns as stages, realized
    exactly. It may not be combined with a non-default proportional setting.
    ``None`` (the default) keeps proportional sizing, byte-identical to the
    generator before the table existed.
    ``surgery_censoring`` is the policy for pool-limited surgeries described in
    the module docstring: ``"error"`` (default) refuses to realize a partial
    surgery, ``"clamp"`` realizes the largest surgery the pool allows. Use
    ``"clamp"`` only to reproduce a historical run or for exploratory work where
    a partial surgery is acceptable — under it the requested effect no longer
    labels the realized construction.
    ``construction`` selects how the effect modes are realized: ``None`` (the
    default; ``"surgery"`` normalizes to it) is the feature-set surgery,
    ``"stage_span"`` the opt-in stage-span construction described in the module
    docstring. ``span_layers`` names the blocks the stage span is computed over
    (``None`` = all three, canonical order) and is accepted only under
    ``stage_span``. Both are omitted from cell ids and signatures while
    ``None``.
    """

    seed: int
    trajectory_mode: TrajectoryMode = "none"
    n_samples: int = 600
    n_stages: int = 3
    group_effect_size: float = 0.0
    group_ratio: float = 0.5
    group_labels: tuple[str, str] = ("A", "B")
    p_dmp: float = 0.2
    baseline_continuity: float = 0.0
    delta_methyl: float = 2.0
    delta_expr: float = 2.0
    delta_protein: float = 2.0
    shape_kind: ShapeKind = "relocate"
    magnitude_kind: MagnitudeKind = "joint"
    surgery_censoring: SurgeryCensoring = "error"
    stage_sample_prop: tuple[float, ...] | None = None
    group_stage_sizes: tuple[tuple[int, ...], tuple[int, ...]] | None = None
    construction: Construction | None = None
    span_layers: tuple[OmicsLayer, ...] | None = None

    def __post_init__(self) -> None:
        if self.construction == "surgery":
            object.__setattr__(self, "construction", None)
        layers = self.span_layers
        if layers is not None and not isinstance(layers, str | bytes) and isinstance(layers, Sequence):
            names = tuple(str(layer) for layer in layers)
            if set(names) <= set(_OMICS_LAYERS) and len(set(names)) == len(names):
                # Canonical order, so equal span sets hash equally.
                names = tuple(layer for layer in _OMICS_LAYERS if layer in names)
            object.__setattr__(self, "span_layers", names)


#: Proportional sizing fields that conflict with an explicit ``group_stage_sizes``.
_PROPORTIONAL_SIZING_FIELDS: tuple[str, ...] = ("n_samples", "stage_sample_prop", "group_ratio")


@dataclass(frozen=True)
class PopulationTrajectories:
    """Exact group-stage means in integration units, keyed by omic layer.

    Each frame uses a ``(group, stage)`` MultiIndex and canonical feature
    columns. Methylation values are M-values rather than stored B values.
    """

    layers: dict[OmicsLayer, pd.DataFrame]


@dataclass(frozen=True)
class SemiSyntheticTrajectoryDataset:
    """MOTCO-ready semi-synthetic trajectory dataset."""

    methylation: pd.DataFrame
    expression: pd.DataFrame
    proteomics: pd.DataFrame
    metadata: pd.DataFrame
    truth: dict[str, Any] = field(default_factory=dict)
    population_trajectories: PopulationTrajectories | None = None


@dataclass(frozen=True)
class _GroupIndicators:
    methyl: np.ndarray
    expr: np.ndarray
    protein: np.ndarray


def generate_semisynthetic_trajectory(
    params: SemiSyntheticTrajectoryParams,
    *,
    reference: IntersimReference | None = None,
) -> SemiSyntheticTrajectoryDataset:
    """Generate a semi-synthetic trajectory dataset using the numpy generator."""

    _validate_params(params)
    ref = reference if reference is not None else load_reference()
    rng = np.random.default_rng(params.seed)

    if params.group_stage_sizes is not None:
        group_a_sizes, group_b_sizes = (list(row) for row in params.group_stage_sizes)
    else:
        stage_sizes = _stage_sizes(params)
        group_a_sizes, group_b_sizes = _group_stage_sizes(stage_sizes, params)

    methyl_a = _baseline_methyl(rng, ref, params)
    deltas_a = (params.delta_methyl, params.delta_expr, params.delta_protein)
    transform_meta: dict[str, Any]
    if params.construction == "stage_span":
        # Draws nothing: B keeps A's indicators and effects, so the stream is the
        # ``none`` stream and B's plant is added to its cell means below.
        methyl_b, deltas_b, transform_meta = methyl_a.copy(), deltas_a, {}
    else:
        methyl_b, deltas_b, transform_meta = _transform_group_b(rng, ref, params, methyl_a)

    indicators_a = _derive_group(methyl_a, ref)
    indicators_b = _derive_group(methyl_b, ref)

    plant: _StageSpanPlant | None = None
    if params.construction == "stage_span":
        plant = _stage_span_plant(ref, params, indicators_a, deltas_a, group_a_sizes, group_b_sizes)
        transform_meta = plant.transform
    shifts_b = plant.shifts if plant is not None else None

    gen_a = _generate_group(rng, ref, indicators_a, deltas_a, group_a_sizes)
    gen_b = _generate_group(rng, ref, indicators_b, deltas_b, group_b_sizes, shifts=shifts_b)
    population_trajectories = _build_population_trajectories(
        ref, params, indicators_a, indicators_b, deltas_a, deltas_b, shifts_b=shifts_b
    )

    methylation, expression, proteomics, metadata = _assemble(
        ref, params, gen_a, gen_b, group_a_sizes, group_b_sizes
    )
    truth = _build_truth(params, indicators_a, indicators_b, deltas_a, deltas_b, transform_meta)
    if params.group_stage_sizes is not None:
        truth["group_stage_sizes"] = {
            params.group_labels[0]: list(group_a_sizes),
            params.group_labels[1]: list(group_b_sizes),
        }
    if plant is not None:
        truth.update(plant.truth)

    return SemiSyntheticTrajectoryDataset(
        methylation=methylation,
        expression=expression,
        proteomics=proteomics,
        metadata=metadata,
        truth=truth,
        population_trajectories=population_trajectories,
    )


# --------------------------------------------------------------------------- #
# Validation
# --------------------------------------------------------------------------- #


def _validate_params(params: SemiSyntheticTrajectoryParams) -> None:
    if params.trajectory_mode not in _MODES:
        raise SemiSyntheticTrajectoryError(f"Unknown trajectory_mode: {params.trajectory_mode}")
    if params.shape_kind not in _SHAPE_KINDS:
        raise SemiSyntheticTrajectoryError(f"Unknown shape_kind: {params.shape_kind}")
    if params.magnitude_kind not in _MAGNITUDE_KINDS:
        raise SemiSyntheticTrajectoryError(f"Unknown magnitude_kind: {params.magnitude_kind}")
    if params.surgery_censoring not in _SURGERY_CENSORINGS:
        raise SemiSyntheticTrajectoryError(
            f"Unknown surgery_censoring: {params.surgery_censoring}; "
            f"expected one of {sorted(_SURGERY_CENSORINGS)}."
        )
    if len(params.group_labels) != 2 or params.group_labels[0] == params.group_labels[1]:
        raise SemiSyntheticTrajectoryError("group_labels must contain two distinct labels.")
    if not (0 < params.group_ratio < 1):
        raise SemiSyntheticTrajectoryError("group_ratio must be between 0 and 1.")
    if params.group_effect_size < 0:
        raise SemiSyntheticTrajectoryError("group_effect_size must be non-negative.")
    if not (0 <= params.p_dmp <= 1):
        raise SemiSyntheticTrajectoryError("p_dmp must be between 0 and 1.")
    if not (0 <= params.baseline_continuity < 1):
        raise SemiSyntheticTrajectoryError(
            "baseline_continuity must be in [0, 1); "
            f"got {params.baseline_continuity}. (1 is a degenerate constant "
            "stage program and is excluded.)"
        )
    if params.n_stages < 2:
        raise SemiSyntheticTrajectoryError("n_stages must be at least 2.")
    validate_construction(params)
    if params.trajectory_mode == "shape" and params.n_stages < 3:
        raise SemiSyntheticTrajectoryError("trajectory_mode='shape' requires at least three stages.")
    for name in ("delta_methyl", "delta_expr", "delta_protein"):
        if getattr(params, name) < 0:
            raise SemiSyntheticTrajectoryError(f"{name} must be non-negative.")
    if params.stage_sample_prop is not None:
        if len(params.stage_sample_prop) != params.n_stages:
            raise SemiSyntheticTrajectoryError("stage_sample_prop must have one entry per stage.")
        if abs(sum(params.stage_sample_prop) - 1.0) > 1e-6:
            raise SemiSyntheticTrajectoryError("stage_sample_prop must sum to 1.")
    if params.group_stage_sizes is not None:
        _validate_group_stage_sizes(params)


def validate_construction(params: SemiSyntheticTrajectoryParams) -> None:
    """Reject an unknown construction, misplaced or malformed span layers, or an infeasible stage-span edit."""

    construction = params.construction
    if construction is not None and construction not in _CONSTRUCTIONS:
        raise SemiSyntheticTrajectoryError(
            f"Unknown construction: {construction!r}; expected one of {sorted(_CONSTRUCTIONS)}."
        )
    if construction != "stage_span":
        if params.span_layers is not None:
            raise SemiSyntheticTrajectoryError(
                f"span_layers {_display_layers(params.span_layers)} is set but construction is not "
                "'stage_span'; span_layers names the blocks of the stage-span construction only."
            )
        return
    resolve_span_layers(params)
    mode = params.trajectory_mode
    if mode in ("orientation", "shape") and params.n_stages < 3:
        raise SemiSyntheticTrajectoryError(
            f"construction='stage_span' with trajectory_mode={mode!r} needs at least three stages: "
            f"the stage span of a {params.n_stages}-stage configuration has no second direction."
        )
    if mode == "orientation" and params.group_effect_size > STAGE_SPAN_MAX_ORIENTATION_EFFECT:
        raise SemiSyntheticTrajectoryError(
            f"construction='stage_span' orientation rotates by 90°·e, so group_effect_size must be at most "
            f"{STAGE_SPAN_MAX_ORIENTATION_EFFECT:g} (180°); got {params.group_effect_size}."
        )


def resolve_span_layers(params: SemiSyntheticTrajectoryParams) -> tuple[OmicsLayer, ...]:
    """The blocks a stage-span construction is computed over, in canonical order.

    ``None`` means all three blocks. A declared set must be non-empty, known
    and duplicate-free; anything else is an error naming the value.
    """

    raw = params.span_layers
    if raw is None:
        return _OMICS_LAYERS
    if isinstance(raw, str | bytes) or not isinstance(raw, Sequence):
        raise SemiSyntheticTrajectoryError(
            f"span_layers must be a list of layer names from {list(_OMICS_LAYERS)}; got {raw!r}."
        )
    names = [str(layer) for layer in raw]
    if not names:
        raise SemiSyntheticTrajectoryError(
            f"span_layers must not be empty; allowed layers: {list(_OMICS_LAYERS)}."
        )
    unknown = sorted(set(names) - set(_OMICS_LAYERS))
    if unknown:
        raise SemiSyntheticTrajectoryError(
            f"span_layers names unknown layer(s) {unknown}; allowed layers: {list(_OMICS_LAYERS)}."
        )
    if len(set(names)) != len(names):
        raise SemiSyntheticTrajectoryError(f"span_layers repeats a layer: {names}.")
    return tuple(layer for layer in _OMICS_LAYERS if layer in names)


def _display_layers(layers: Any) -> Any:
    return list(layers) if isinstance(layers, tuple | list) else layers


def _validate_group_stage_sizes(params: SemiSyntheticTrajectoryParams) -> None:
    defaults = SemiSyntheticTrajectoryParams(seed=params.seed)
    conflicts = [
        name for name in _PROPORTIONAL_SIZING_FIELDS if getattr(params, name) != getattr(defaults, name)
    ]
    if conflicts:
        raise SemiSyntheticTrajectoryError(
            "group_stage_sizes fixes every group-stage cell size and cannot be combined with "
            f"non-default proportional sizing: {', '.join(conflicts)}."
        )
    table = params.group_stage_sizes
    assert table is not None
    if len(table) != 2:
        raise SemiSyntheticTrajectoryError(
            f"group_stage_sizes must have exactly two rows (one per group label); got {len(table)}."
        )
    for label, row in zip(params.group_labels, table, strict=True):
        if len(row) != params.n_stages:
            raise SemiSyntheticTrajectoryError(
                f"group_stage_sizes row for group {label!r} must have one entry per stage "
                f"(n_stages={params.n_stages}); got {len(row)}."
            )
        if any(isinstance(n, bool) or not isinstance(n, int | np.integer) or n < 1 for n in row):
            raise SemiSyntheticTrajectoryError(
                f"group_stage_sizes cells must be integers >= 1; group {label!r} has {list(row)}."
            )


# --------------------------------------------------------------------------- #
# Sizing
# --------------------------------------------------------------------------- #


def _stage_sizes(params: SemiSyntheticTrajectoryParams) -> list[int]:
    k = params.n_stages
    prop = params.stage_sample_prop or tuple([1.0 / k] * k)
    sizes = [int(round(params.n_samples * prop[i])) for i in range(k - 1)]
    sizes.append(params.n_samples - sum(sizes))
    if any(n < 2 for n in sizes):
        raise SemiSyntheticTrajectoryError(
            "Each stage needs at least 2 samples (one per group); increase n_samples."
        )
    return sizes


def _group_stage_sizes(
    stage_sizes: list[int], params: SemiSyntheticTrajectoryParams
) -> tuple[list[int], list[int]]:
    a_sizes, b_sizes = [], []
    for n in stage_sizes:
        a = int(round(n * params.group_ratio))
        a = min(max(a, 1), n - 1)
        a_sizes.append(a)
        b_sizes.append(n - a)
    return a_sizes, b_sizes


# --------------------------------------------------------------------------- #
# Methylation baseline and the group-B transform (methylation only)
# --------------------------------------------------------------------------- #


def _baseline_methyl(
    rng: np.random.Generator,
    ref: IntersimReference,
    params: SemiSyntheticTrajectoryParams,
) -> np.ndarray:
    """Group A's per-stage methylation differential indicators (binary).

    Drawn as a stationary Markov chain along the stage axis with persistence
    ``params.baseline_continuity``; at the default 0 this is exactly the
    independent per-stage draw (see :func:`~motco.simulations.generator.markov_indicators`).
    """

    return markov_indicators(
        rng, ref.n_cpg, params.n_stages, params.p_dmp, params.baseline_continuity
    )


def _derive_group(methyl: np.ndarray, ref: IntersimReference) -> _GroupIndicators:
    """Derive expression and protein indicators from methylation via the cascade."""

    expr, protein = derive_coupled_indicators((methyl > 0).astype(float), ref)
    return _GroupIndicators(methyl=methyl, expr=expr, protein=protein)


def _transform_group_b(
    rng: np.random.Generator,
    ref: IntersimReference,
    params: SemiSyntheticTrajectoryParams,
    methyl_a: np.ndarray,
) -> tuple[np.ndarray, tuple[float, float, float], dict[str, Any]]:
    """Return group B's methylation indicators, per-omic deltas, and truth notes."""

    deltas_a = (params.delta_methyl, params.delta_expr, params.delta_protein)
    e = params.group_effect_size
    mode = params.trajectory_mode

    if mode == "none" or e == 0:
        return methyl_a.copy(), deltas_a, {}

    if mode == "translation":
        return _translation_methyl(rng, ref, params, methyl_a)

    if mode == "magnitude":
        return _magnitude_methyl(methyl_a, e, params)

    if mode == "orientation":
        return _orientation_methyl(rng, methyl_a, e, params)

    if mode == "shape":
        return _shape_methyl(rng, methyl_a, e, params)

    raise SemiSyntheticTrajectoryError(f"Unknown trajectory_mode: {mode}")


def _translation_methyl(
    rng: np.random.Generator,
    ref: IntersimReference,
    params: SemiSyntheticTrajectoryParams,
    methyl_a: np.ndarray,
) -> tuple[np.ndarray, tuple[float, float, float], dict[str, Any]]:
    """Add an extra constant set U of DMPs on an independent gene program.

    For a clean location-only control, U's CpGs must regulate genes that the
    stage trajectory does *not* touch — otherwise the CpG→gene OR-derivation
    would let U saturate (flatten) stage-varying genes in group B, deforming the
    derived trajectory. So U is drawn from stage-inactive CpGs whose mapped gene
    is also absent from the stage program.
    """

    methyl_b = methyl_a.copy()
    stage_active = methyl_a.sum(1) > 0
    used_genes = ref.incidence_cpg_gene[stage_active].sum(0) > 0  # genes in the stage program
    cpg_gene = ref.incidence_cpg_gene.argmax(1)  # each CpG's mapped gene
    fresh = (~stage_active) & (~used_genes[cpg_gene])  # stage-inactive CpGs on fresh genes
    candidates = np.where(fresh)[0]
    requested = params.p_dmp * ref.n_cpg
    nominal = int(round(params.group_effect_size * requested))
    n_extra = _resolve_surgery_size(
        policy=params.surgery_censoring,
        surgery="translation",
        nominal=nominal,
        pool=len(candidates),
        saturating_effect=len(candidates) / requested if requested else 0.0,
    )
    if n_extra >= 1:
        u = rng.choice(candidates, size=n_extra, replace=False)
        methyl_b[u, :] = 1.0  # differential at every B stage, none of A's
    deltas = (params.delta_methyl, params.delta_expr, params.delta_protein)
    return (
        methyl_b,
        deltas,
        {
            "translation_set_size": int(n_extra),
            "translation_nominal": nominal,
            "censored": bool(n_extra < nominal),
        },
    )


def _magnitude_methyl(
    methyl_a: np.ndarray,
    e: float,
    params: SemiSyntheticTrajectoryParams,
) -> tuple[np.ndarray, tuple[float, float, float], dict[str, Any]]:
    """Scale group B's effect size: methylation δ, endpoint indicators, or every δ.

    ``magnitude_kind='all'`` scales the global methylation δ, uniformly
    enlarging every methylation step (the original behavior). ``'extremes'``
    instead leaves δ unchanged and scales A's methylation *indicators* at the
    first and last stages only — a localized size change at the endpoints that
    leaves interior vertices at the baseline effect, probing whether confining
    the scale reduces the shape co-movement seen with the all-stage variant.
    ``'joint'`` (default) scales every omic's δ by the same ``1 + e``: ``'all'`` is
    size-pure within each omic block but rotates the *concatenated* trajectory
    (one block grew while the others did not), whereas ``'joint'`` is size-pure
    in the joint standardized space as well. None of the branches consumes
    randomness, so the RNG stream is identical across kinds.
    """

    if params.magnitude_kind == "joint":
        scale = 1.0 + e
        scaled = (
            float(scale * params.delta_methyl),
            float(scale * params.delta_expr),
            float(scale * params.delta_protein),
        )
        return (
            methyl_a.copy(),
            scaled,
            {
                "magnitude_kind": "joint",
                "delta_methyl_scale": scale,
                "delta_expr_scale": scale,
                "delta_protein_scale": scale,
            },
        )

    if params.magnitude_kind == "extremes":
        methyl_b = methyl_a.astype(float).copy()
        endpoints = (0, params.n_stages - 1)
        for stage in endpoints:
            methyl_b[:, stage] = methyl_a[:, stage] * (1.0 + e)
        deltas = (params.delta_methyl, params.delta_expr, params.delta_protein)
        return methyl_b, deltas, {"magnitude_kind": "extremes", "magnitude_scale": 1.0 + e}

    scaled = (float((1.0 + e) * params.delta_methyl), params.delta_expr, params.delta_protein)
    return methyl_a.copy(), scaled, {"magnitude_kind": "all", "delta_methyl_scale": 1.0 + e}


def _resolve_surgery_size(
    *,
    policy: SurgeryCensoring,
    surgery: str,
    nominal: int,
    pool: int,
    saturating_effect: float,
) -> int:
    """Resolve a requested surgery size against its pool under the censoring policy.

    Reads pool sizes only and consumes no randomness, so the RNG call sequence
    under ``"clamp"`` is identical to the pre-policy generator's at every seed.
    """

    if nominal <= pool:
        return nominal
    if policy == "clamp":
        return pool
    raise SemiSyntheticTrajectoryError(
        f"trajectory_mode={surgery!r} requested a surgery of {nominal} site(s) but only {pool} "
        f"are available in the destination pool, so the requested effect cannot be realized in "
        f"full: this draw saturates at group_effect_size={saturating_effect:.4f}. Lower "
        "group_effect_size (or p_dmp), or set surgery_censoring='clamp' to accept a partial "
        "surgery whose realized size no longer matches the requested effect."
    )


def _relocate_rows(
    rng: np.random.Generator,
    src_pool: np.ndarray,
    dst_pool: np.ndarray,
    fraction: float,
    *,
    policy: SurgeryCensoring,
    surgery: str,
) -> tuple[np.ndarray, np.ndarray, int]:
    """Pick ``fraction`` of ``src_pool`` and an equal number from ``dst_pool``.

    Returns the source rows, the destination rows, and the *nominal* (requested)
    size; the realized size is ``src.size``.
    """

    nominal = int(round(min(max(fraction, 0.0), 1.0) * len(src_pool)))
    k = _resolve_surgery_size(
        policy=policy,
        surgery=surgery,
        nominal=nominal,
        pool=len(dst_pool),
        saturating_effect=len(dst_pool) / len(src_pool) if len(src_pool) else 0.0,
    )
    if k < 1:
        return np.empty(0, dtype=int), np.empty(0, dtype=int), nominal
    src = rng.choice(src_pool, size=k, replace=False)
    dst = rng.choice(dst_pool, size=k, replace=False)
    return src, dst, nominal


def _orientation_methyl(
    rng: np.random.Generator,
    methyl_a: np.ndarray,
    e: float,
    params: SemiSyntheticTrajectoryParams,
) -> tuple[np.ndarray, tuple[float, float, float], dict[str, Any]]:
    """Relocate a fraction of stage-changing sites, the same relocation at every stage."""

    active = np.where(methyl_a.sum(1) > 0)[0]
    inactive = np.where(methyl_a.sum(1) == 0)[0]
    src, dst, nominal = _relocate_rows(
        rng, active, inactive, e, policy=params.surgery_censoring, surgery="orientation"
    )
    methyl_b = methyl_a.copy()
    if src.size:
        methyl_b[dst, :] = methyl_a[src, :]  # move whole rows → same relocation per stage
        methyl_b[src, :] = 0.0
    deltas = (params.delta_methyl, params.delta_expr, params.delta_protein)
    return (
        methyl_b,
        deltas,
        {
            "orientation_relocated": int(src.size),
            "orientation_nominal": nominal,
            "censored": bool(src.size < nominal),
        },
    )


def _shape_methyl(
    rng: np.random.Generator,
    methyl_a: np.ndarray,
    e: float,
    params: SemiSyntheticTrajectoryParams,
) -> tuple[np.ndarray, tuple[float, float, float], dict[str, Any]]:
    """Perturb a single interior stage: relocate its sites or scale its effect."""

    methyl_b = methyl_a.astype(float).copy()
    interior = list(range(1, params.n_stages - 1))
    stage = interior[len(interior) // 2]  # the single interior stage to perturb
    meta: dict[str, Any] = {"shape_stage": stage, "shape_kind": params.shape_kind}

    if params.shape_kind == "magnitude":
        methyl_b[:, stage] = methyl_a[:, stage] * (1.0 + e)
        meta["shape_scale"] = 1.0 + e
    else:  # relocate this stage's sites to globally-inactive CpGs → bends the vertex
        active_here = np.where(methyl_a[:, stage] > 0)[0]
        global_inactive = np.where(methyl_a.sum(1) == 0)[0]
        src, dst, nominal = _relocate_rows(
            rng, active_here, global_inactive, e, policy=params.surgery_censoring, surgery="shape"
        )
        if src.size:
            methyl_b[dst, stage] = methyl_a[src, stage]
            methyl_b[src, stage] = 0.0
        meta["shape_relocated"] = int(src.size)
        meta["shape_nominal"] = nominal
        meta["censored"] = bool(src.size < nominal)
    deltas = (params.delta_methyl, params.delta_expr, params.delta_protein)
    return methyl_b, deltas, meta


# --------------------------------------------------------------------------- #
# Generation + assembly
# --------------------------------------------------------------------------- #


def _generate_group(
    rng: np.random.Generator,
    ref: IntersimReference,
    indicators: _GroupIndicators,
    deltas: tuple[float, float, float],
    cell_sizes: list[int],
    *,
    shifts: Mapping[OmicsLayer, np.ndarray] | None = None,
) -> GeneratedOmics:
    return generate_omics(
        cell_sizes=cell_sizes,
        indicators_methyl=indicators.methyl,
        indicators_expr=indicators.expr,
        indicators_protein=indicators.protein,
        delta_methyl=deltas[0],
        delta_expr=deltas[1],
        delta_protein=deltas[2],
        rng=rng,
        reference=ref,
        **_cell_shift_kwargs(shifts),
    )


def _cell_shift_kwargs(shifts: Mapping[OmicsLayer, np.ndarray] | None) -> dict[str, np.ndarray | None]:
    shifts = shifts or {}
    return {
        "cell_shift_methyl": shifts.get("methylation"),
        "cell_shift_expr": shifts.get("expression"),
        "cell_shift_protein": shifts.get("proteomics"),
    }


def _assemble(
    ref: IntersimReference,
    params: SemiSyntheticTrajectoryParams,
    gen_a: GeneratedOmics,
    gen_b: GeneratedOmics,
    a_sizes: list[int],
    b_sizes: list[int],
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    label_a, label_b = params.group_labels
    groups = [label_a] * sum(a_sizes) + [label_b] * sum(b_sizes)
    stages = np.concatenate([gen_a.cell_ids, gen_b.cell_ids]).astype(int)
    sample_ids = [f"sample{i}" for i in range(len(groups))]

    metadata = pd.DataFrame(
        {
            "sample_id": sample_ids,
            "group": groups,
            "stage": stages,
            "cluster": stages,
        },
        index=sample_ids,
    )

    methylation = pd.DataFrame(
        np.vstack([gen_a.methylation, gen_b.methylation]),
        index=sample_ids,
        columns=ref.cpg_names,
    )
    expression = pd.DataFrame(
        np.vstack([gen_a.expression, gen_b.expression]),
        index=sample_ids,
        columns=ref.gene_names,
    )
    proteomics = pd.DataFrame(
        np.vstack([gen_a.proteomics, gen_b.proteomics]),
        index=sample_ids,
        columns=ref.protein_names,
    )
    return methylation, expression, proteomics, metadata


def _build_truth(
    params: SemiSyntheticTrajectoryParams,
    indicators_a: _GroupIndicators,
    indicators_b: _GroupIndicators,
    deltas_a: tuple[float, float, float],
    deltas_b: tuple[float, float, float],
    transform_meta: dict[str, Any],
) -> dict[str, Any]:
    def counts(ind: _GroupIndicators) -> dict[str, list[int]]:
        return {
            "methylation": (ind.methyl != 0).sum(0).astype(int).tolist(),
            "expression": (ind.expr != 0).sum(0).astype(int).tolist(),
            "proteomics": (ind.protein != 0).sum(0).astype(int).tolist(),
        }

    label_a, label_b = params.group_labels
    return {
        "trajectory_mode": params.trajectory_mode,
        "group_effect_size": params.group_effect_size,
        "group_labels": [label_a, label_b],
        "group_ratio": params.group_ratio,
        "n_stages": params.n_stages,
        "p_dmp": params.p_dmp,
        "baseline_continuity": params.baseline_continuity,
        "shape_kind": params.shape_kind,
        "magnitude_kind": params.magnitude_kind,
        "surgery_censoring": params.surgery_censoring,
        "seed": params.seed,
        "stage_assumption": "clusters-as-stages",
        "deltas": {label_a: list(deltas_a), label_b: list(deltas_b)},
        "indicator_counts": {label_a: counts(indicators_a), label_b: counts(indicators_b)},
        "indicators": {
            label_a: {
                "methylation": indicators_a.methyl,
                "expression": indicators_a.expr,
                "proteomics": indicators_a.protein,
            },
            label_b: {
                "methylation": indicators_b.methyl,
                "expression": indicators_b.expr,
                "proteomics": indicators_b.protein,
            },
        },
        "transform": transform_meta,
    }


def _build_population_trajectories(
    ref: IntersimReference,
    params: SemiSyntheticTrajectoryParams,
    indicators_a: _GroupIndicators,
    indicators_b: _GroupIndicators,
    deltas_a: tuple[float, float, float],
    deltas_b: tuple[float, float, float],
    *,
    shifts_b: Mapping[OmicsLayer, np.ndarray] | None = None,
) -> PopulationTrajectories:
    """Build exact population trajectories using the generator's mean contract.

    ``shifts_b`` is group B's stage-span plant, if any, so its population means
    are the ones its samples were drawn from.
    """

    def means(
        indicators: _GroupIndicators,
        deltas: tuple[float, float, float],
        shifts: Mapping[OmicsLayer, np.ndarray] | None = None,
    ) -> dict[str, np.ndarray]:
        return omic_population_means(
            indicators_methyl=indicators.methyl,
            indicators_expr=indicators.expr,
            indicators_protein=indicators.protein,
            delta_methyl=deltas[0],
            delta_expr=deltas[1],
            delta_protein=deltas[2],
            reference=ref,
            **_cell_shift_kwargs(shifts),
        )

    by_group = {
        params.group_labels[0]: means(indicators_a, deltas_a),
        params.group_labels[1]: means(indicators_b, deltas_b, shifts_b),
    }
    feature_names = {
        "methylation": ref.cpg_names,
        "expression": ref.gene_names,
        "proteomics": ref.protein_names,
    }
    index = pd.MultiIndex.from_product(
        [params.group_labels, [str(stage) for stage in range(params.n_stages)]],
        names=["group", "stage"],
    )
    layers: dict[OmicsLayer, pd.DataFrame] = {}
    for layer in _OMICS_LAYERS:
        values = np.vstack([by_group[group][layer] for group in params.group_labels])
        layers[layer] = pd.DataFrame(values, index=index, columns=feature_names[layer])
    return PopulationTrajectories(layers=layers)


# --------------------------------------------------------------------------- #
# Stage-span construction
# --------------------------------------------------------------------------- #


@dataclass(frozen=True)
class _StageSpanPlant:
    """Group B's stage-span plant in analysed units, with its truth notes.

    ``shifts`` maps each span layer to a ``(n_features, n_stages)`` mean shift
    (methylation in M-values) and is empty for a null plant, so a null cell
    samples exactly as ``none`` does.
    """

    shifts: dict[OmicsLayer, np.ndarray]
    transform: dict[str, Any]
    truth: dict[str, Any]


def null_pooled_scale(
    null_means: Mapping[str, np.ndarray],
    stage_weights: Sequence[float] | np.ndarray,
    *,
    layers: Sequence[str],
    reference: IntersimReference,
) -> np.ndarray:
    """Null population pooled SD per feature, concatenated over ``layers``.

    ``null_means[layer]`` holds the null stage means (``(n_stages, n_features)``,
    analysed units), shared by both groups, and ``stage_weights`` each stage's
    share of all samples. The variance is the within-cell variance (the
    reference covariance diagonal) plus the between-cell variance of the null
    means at those weights — the population analogue of the pooled sample
    scaler (``ddof = 0``). A feature whose SD is (numerically) zero keeps unit
    scale, as the sample scaler does.
    """

    weights = np.asarray(stage_weights, dtype=float)
    weights = weights / weights.sum()
    within = {
        "methylation": np.diag(reference.cov_M),
        "expression": np.diag(reference.cov_expr),
        "proteomics": np.diag(reference.cov_protein),
    }
    scales = []
    for layer in layers:
        means = np.asarray(null_means[layer], dtype=float)
        pooled_mean = weights @ means
        between = weights @ (means - pooled_mean) ** 2
        scale = np.sqrt(within[layer] + between)
        scale[scale < _SCALE_TOLERANCE] = 1.0
        scales.append(scale)
    return np.concatenate(scales)


def _stage_span_plant(
    ref: IntersimReference,
    params: SemiSyntheticTrajectoryParams,
    indicators: _GroupIndicators,
    deltas: tuple[float, float, float],
    a_sizes: Sequence[int],
    b_sizes: Sequence[int],
) -> _StageSpanPlant:
    """Edit B's standardized null stage means and map the edit back to analysed units."""

    layers = resolve_span_layers(params)
    null_means = omic_population_means(
        indicators_methyl=indicators.methyl,
        indicators_expr=indicators.expr,
        indicators_protein=indicators.protein,
        delta_methyl=deltas[0],
        delta_expr=deltas[1],
        delta_protein=deltas[2],
        reference=ref,
    )
    stage_weights = np.asarray(a_sizes, dtype=float) + np.asarray(b_sizes, dtype=float)
    scale = null_pooled_scale(null_means, stage_weights, layers=layers, reference=ref)
    configuration = np.hstack([null_means[layer] for layer in layers])
    z = (configuration - stage_weights @ configuration / stage_weights.sum()) / scale

    null = params.trajectory_mode == "none" or params.group_effect_size == 0
    effect: dict[str, Any]
    if null:
        edited, effect = z, {}
    else:
        edited, effect = edit_stage_configuration(z, params.trajectory_mode, float(params.group_effect_size))
    plant = edited - z

    shifts: dict[OmicsLayer, np.ndarray] = {}
    if not null:
        offset = 0
        plant_units = plant * scale
        for layer in layers:
            width = null_means[layer].shape[1]
            shifts[layer] = plant_units[:, offset : offset + width].T.copy()
            offset += width

    label_a, label_b = params.group_labels
    geometry = _intended_geometry(z, edited, (label_a, label_b))
    transform = {} if null else {"construction": "stage_span", **effect}
    truth = {
        "construction": "stage_span",
        "span_layers": list(layers),
        "stage_span": {
            "units": "null_standardized",
            "path_length": _path_length(z),
            "mean_step": _path_length(z) / (z.shape[0] - 1),
            "intended_geometry": geometry,
            "plant_norm": [float(norm) for norm in np.linalg.norm(plant, axis=1)],
        },
    }
    return _StageSpanPlant(shifts=shifts, transform=transform, truth=truth)


def edit_stage_configuration(
    z: np.ndarray,
    mode: str,
    effect: float,
) -> tuple[np.ndarray, dict[str, Any]]:
    """Apply one stage-span edit to a stage configuration ``z`` (stages × features).

    Returns the edited configuration and the effect parameters, in the units of
    ``z``. See the module docstring for the four edits; ``orientation`` and
    ``shape`` raise when the configuration has no second in-span direction.
    """

    z = np.asarray(z, dtype=float)
    n_stages = z.shape[0]
    centroid = z.mean(axis=0)
    centered = z - centroid
    path_length = _path_length(z)
    mean_step = path_length / (n_stages - 1)

    if mode == "magnitude":
        scale = 1.0 + effect
        return centroid + scale * centered, {"scale": scale}

    net = z[-1] - z[0]
    net_norm = float(np.linalg.norm(net))
    if net_norm <= STAGE_SPAN_RANK_TOLERANCE * max(path_length, np.finfo(float).tiny):
        raise SemiSyntheticTrajectoryError(
            f"construction='stage_span' {mode} needs a net displacement, but the stage configuration "
            "returns to its first stage."
        )
    u = net / net_norm

    if mode == "translation":
        offset = effect * mean_step
        return z + offset * u, {"offset_length": offset}

    if mode not in ("orientation", "shape"):
        raise SemiSyntheticTrajectoryError(f"Unknown stage-span trajectory_mode: {mode}")
    if n_stages < 3:
        raise SemiSyntheticTrajectoryError(
            f"construction='stage_span' {mode} needs at least three stages: the stage span of a "
            f"{n_stages}-stage configuration has no second direction."
        )
    singular = np.linalg.svd(centered, compute_uv=False)
    if singular.size < 2 or singular[1] <= STAGE_SPAN_RANK_TOLERANCE * singular[0]:
        raise SemiSyntheticTrajectoryError(
            f"construction='stage_span' {mode}: the stage span has no second direction (second singular "
            f"value {float(singular[1]) if singular.size > 1 else 0.0:.3g} against a first of "
            f"{float(singular[0]):.3g})."
        )

    if mode == "orientation":
        v = _leading_direction(centered - np.outer(centered @ u, u))
        theta = np.deg2rad(90.0 * effect)
        along_u = centered @ u
        along_v = centered @ v
        rotated = (
            centered
            + (np.cos(theta) - 1.0) * (np.outer(along_u, u) + np.outer(along_v, v))
            + np.sin(theta) * (np.outer(along_u, v) - np.outer(along_v, u))
        )
        return centroid + rotated, {"rotation_degrees": float(90.0 * effect)}

    interior = list(range(1, n_stages - 1))
    stage = interior[len(interior) // 2]  # the stage the shape surgery perturbs
    chord = z[stage + 1] - z[stage - 1]
    chord_norm = float(np.linalg.norm(chord))
    off_chord = z[stage] - z[stage - 1]
    if chord_norm > 0:
        chord_unit = chord / chord_norm
        off_chord = off_chord - (off_chord @ chord_unit) * chord_unit
    else:
        chord_unit = np.zeros_like(chord)
    if float(np.linalg.norm(off_chord)) > STAGE_SPAN_RANK_TOLERANCE * mean_step:
        direction = off_chord / np.linalg.norm(off_chord)
    else:
        # The stage sits on its chord: bend along the configuration's leading
        # in-span direction perpendicular to the chord.
        direction = _leading_direction(centered - np.outer(centered @ chord_unit, chord_unit))
    displacement = effect * mean_step
    bent = z.copy()
    bent[stage] = bent[stage] + displacement * direction
    rescale = path_length / _path_length(bent)
    edited = centroid + rescale * (bent - bent.mean(axis=0))
    return edited, {"shape_stage": stage, "displacement": float(displacement), "rescale": float(rescale)}


def _leading_direction(matrix: np.ndarray) -> np.ndarray:
    """Unit leading right singular vector, signed so its largest stage score is positive."""

    left, singular, right = np.linalg.svd(matrix, full_matrices=False)
    if singular.size == 0 or singular[0] <= 0:
        raise SemiSyntheticTrajectoryError("construction='stage_span': the stage span has no second direction.")
    direction = right[0]
    scores = left[:, 0]
    if scores[int(np.argmax(np.abs(scores)))] < 0:
        direction = -direction
    return direction


def _path_length(configuration: np.ndarray) -> float:
    return float(np.linalg.norm(np.diff(configuration, axis=0), axis=1).sum())


def _intended_geometry(z: np.ndarray, edited: np.ndarray, labels: tuple[str, str]) -> dict[str, Any]:
    """Population ``delta`` / net-displacement ``angle`` / ``shape`` of B's edit against A."""

    n_stages = z.shape[0]
    ordered = np.vstack([z, edited])
    identity = np.eye(2 * n_stages)
    contrast = [list(range(n_stages)), list(range(n_stages, 2 * n_stages))]
    deltas, angles, shapes = estimate_difference(
        ordered, identity, identity, contrast, orientation_estimator="net_displacement"
    )
    angle = float(angles[0, 1])
    return {
        "orientation_estimator": "net_displacement",
        "path_lengths": {labels[0]: _path_length(z), labels[1]: _path_length(edited)},
        "delta": float(deltas[0, 1]),
        "angle": angle if np.isfinite(angle) else None,
        "shape": float(shapes[0, 1]) if n_stages >= 3 else None,
    }


# --------------------------------------------------------------------------- #
# Convenience helpers consumed downstream
# --------------------------------------------------------------------------- #


#: Trajectory modes whose surgery size is limited by a destination/candidate pool.
POOL_LIMITED_MODES: frozenset[str] = frozenset({"orientation", "translation", "shape"})

#: Guard band, in standard deviations of the destination pool, applied by
#: :func:`expected_surgery_headroom`. The pools are random per replicate, so a
#: config that fits only *in expectation* would still fail loudly at runtime on
#: unlucky draws; three SDs is conservative and cheap.
HEADROOM_GUARD_SIGMAS: float = 3.0


def expected_stage_active_fraction(p_dmp: float, n_stages: int, continuity: float = 0.0) -> float:
    """Probability that a CpG is differential in at least one stage.

    Under the baseline's stationary Markov chain the CpG stays non-differential
    across every stage with probability ``(1 − p)·(1 − p·(1 − ρ))^(n − 1)``
    (start non-differential, then take ``n − 1`` non-entering transitions), so
    the stage-active union is its complement. At ρ = 0 this is the independence
    union ``1 − (1 − p)^n``.
    """

    if n_stages <= 0:
        return 0.0
    stay_inactive = 1.0 - p_dmp * (1.0 - continuity)
    return 1.0 - (1.0 - p_dmp) * stay_inactive ** (n_stages - 1)


@dataclass(frozen=True)
class SurgeryHeadroom:
    """Expected pool headroom for one cell's pool-limited surgery.

    All sizes are *expectations* over the baseline indicator draw, in CpG counts.
    ``available`` is the expected pool less the guard band, and
    ``saturating_effect`` is the ``group_effect_size`` at which the requested
    surgery would exactly consume it. ``fits`` is the enumeration-time verdict.
    """

    trajectory_mode: str
    group_effect_size: float
    nominal: float
    pool: float
    guard_band: float
    available: float
    saturating_effect: float
    fits: bool


def expected_surgery_headroom(
    params: SemiSyntheticTrajectoryParams,
    *,
    reference: IntersimReference | None = None,
) -> SurgeryHeadroom | None:
    """Expected destination-pool headroom for ``params``'s pool-limited surgery.

    Returns ``None`` for modes that perform no pool-limited surgery (``none``,
    ``magnitude``, and ``shape`` with ``shape_kind='magnitude'``), for a
    zero effect, where the transform short-circuits to group A's baseline, and
    for the stage-span construction, which draws from no feature pool.

    Orientation and shape use the expected stage-active fraction under the
    baseline's Markov continuity ρ,
    ``a = 1 − (1 − p_dmp)·(1 − p_dmp·(1 − ρ))^(n_stages − 1)`` — the complement
    of never visiting the differential state in ``n_stages`` steps — which
    reduces to the independence union ``1 − (1 − p_dmp)^n_stages`` at ρ = 0.
    Higher continuity makes stage programs overlap, shrinking the active union
    and *growing* the destination pool, so headroom improves along the axis.
    Translation's candidate pool depends on the CpG→gene incidence — a CpG
    qualifies only if it is stage-inactive *and* every CpG mapped to its gene is
    too — so it is computed from the cached reference maps rather than
    approximated.
    """

    _validate_params(params)
    mode = params.trajectory_mode
    effect = float(params.group_effect_size)
    if params.construction == "stage_span" or mode not in POOL_LIMITED_MODES or effect == 0.0:
        return None
    if mode == "shape" and params.shape_kind != "relocate":
        return None

    ref = reference if reference is not None else load_reference()
    n = float(ref.n_cpg)
    active_fraction = expected_stage_active_fraction(
        params.p_dmp, params.n_stages, params.baseline_continuity
    )

    if mode == "translation":
        pool = _expected_translation_pool(ref, active_fraction)
        source = params.p_dmp * n
    else:
        pool = (1.0 - active_fraction) * n
        # orientation draws from every stage-active CpG; shape from one stage's.
        source = active_fraction * n if mode == "orientation" else params.p_dmp * n

    band = HEADROOM_GUARD_SIGMAS * _pool_sd(pool, n)
    available = max(pool - band, 0.0)
    nominal = effect * source
    saturating = available / source if source > 0 else 0.0
    return SurgeryHeadroom(
        trajectory_mode=mode,
        group_effect_size=effect,
        nominal=nominal,
        pool=pool,
        guard_band=band,
        available=available,
        saturating_effect=saturating,
        fits=nominal <= available,
    )


def _pool_sd(pool: float, n_cpg: float) -> float:
    """Binomial SD of a pool of expected size ``pool`` drawn over ``n_cpg`` CpGs."""

    if n_cpg <= 0:
        return 0.0
    p = min(max(pool / n_cpg, 0.0), 1.0)
    return float(np.sqrt(n_cpg * p * (1.0 - p)))


def _expected_translation_pool(ref: IntersimReference, active_fraction: float) -> float:
    """Expected size of translation's fresh-CpG candidate pool.

    A CpG qualifies when it is stage-inactive *and* its mapped gene carries no
    stage-active CpG — i.e. when every CpG incident to that gene is inactive.
    Stage-activity is independent across CpGs, so the probability for CpG ``i``
    is ``(1 − a)`` raised to the size of that CpG's gene neighbourhood
    (including ``i`` itself, which the generator's ``argmax`` mapping may leave
    outside the incidence column for a CpG mapped to no gene).

    The *mean* this returns is exact. Note that the candidate indicators are
    positively correlated through shared genes, so the binomial guard band
    :func:`expected_surgery_headroom` puts around it is a slight under-estimate
    of this pool's spread — a translation cell sitting exactly on the headroom
    boundary can still censor on an unlucky draw, which the runtime policy then
    reports.
    """

    incidence = ref.incidence_cpg_gene
    mapped_gene = incidence.argmax(1)
    gene_sizes = incidence.sum(0)
    rows = np.arange(incidence.shape[0])
    # +1 wherever the CpG is not itself incident to the gene argmax selected.
    neighbourhood = gene_sizes[mapped_gene] + (1.0 - incidence[rows, mapped_gene])
    return float(np.sum((1.0 - active_fraction) ** neighbourhood))


def affected_omics_layers() -> tuple[OmicsLayer, ...]:
    """Canonical omic-layer order."""

    return _OMICS_LAYERS


def list_constructions() -> Sequence[Construction]:
    """All supported constructions (``surgery`` is the default)."""

    return ("surgery", "stage_span")


def list_trajectory_modes() -> Sequence[TrajectoryMode]:
    """All supported trajectory modes."""

    return ("none", "translation", "magnitude", "orientation", "shape")
