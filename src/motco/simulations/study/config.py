"""Declarative study configuration schema and loader."""

from __future__ import annotations

import dataclasses
import itertools
import json
from collections.abc import Mapping, Sequence
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

from motco.simulations.evaluation import AttributionDiagnosticSettings, SimulationEvaluationParams
from motco.simulations.grid import SimulationGridError, _split_axis
from motco.simulations.semisynthetic import (
    _MAGNITUDE_KINDS,
    SemiSyntheticTrajectoryParams,
)

_TRAJECTORY_MODES = {"none", "translation", "magnitude", "orientation", "shape"}
_AXIS_NAMESPACES = {"generator", "evaluation"}
_TARGET_KINDS = {"type_i_control", "power_monotonicity", "specificity"}
_TARGET_BOUNDS = {"two_sided", "upper"}
_GATE_ROLES = {"mandatory_power", "mandatory_control", "descriptive"}
_STATISTICS = {"delta", "angle", "shape"}
_ATTRIBUTION_COMPONENTS = ("observed", "pls_captured", "residual")
_N_JOBS_OVERRIDE_POLICIES = ("forbid", "warn")
_DRIVER_AGREEMENT_DESCRIPTIVE = "descriptive"
_ROOT_KEYS = {
    "generator",
    "evaluation",
    "trajectory_modes",
    "effect_sizes",
    "axes",
    "design_grid",
    "n_replicates",
    "base_seed",
    "alpha",
    "acceptance",
    "attribution",
    "matched_seeds",
    "report_contract",
    "metadata",
}


class StudyConfigError(ValueError):
    """Raised when a study configuration is invalid."""


@dataclass(frozen=True)
class TypeIControlTarget:
    """Type I control target evaluated on null cells."""

    alpha: float
    se_tolerance: float = 2.0
    name: str = "type_i_control"
    bound: str = "two_sided"
    kind: str = field(default="type_i_control", init=False)


@dataclass(frozen=True)
class PowerMonotonicityTarget:
    """Monotone power and floor at the largest effect size."""

    trajectory_mode: str
    statistic: str
    min_power_at_top: float
    name: str = "power_monotonicity"
    kind: str = field(default="power_monotonicity", init=False)


@dataclass(frozen=True)
class SpecificityTarget:
    """Off-diagonal specificity (rate ≈ alpha within tolerance)."""

    trajectory_mode: str
    statistic: str
    alpha: float
    se_tolerance: float = 2.0
    name: str = "specificity"
    bound: str = "two_sided"
    kind: str = field(default="specificity", init=False)


@dataclass(frozen=True)
class AttributionSelector:
    """Which enumerated cells receive orientation-attribution diagnostics.

    The selector is resolved into per-cell evaluation parameters during
    enumeration, so eligibility and the effective bootstrap settings are part of
    the cell's parameter signature rather than something inferred after a run.
    Selection is never conditioned on an observed p-value.
    """

    enabled: bool = False
    trajectory_modes: tuple[str, ...] = ("orientation",)
    phases: tuple[str, ...] = ("power_primary",)
    nonzero_effects_only: bool = True
    effect_sizes: tuple[float, ...] | None = None
    bootstrap_replicates: int = 100
    bootstrap_seed: int | None = 0
    top_k: int = 20
    zero_tolerance: float = 1e-12

    def __post_init__(self) -> None:
        if not self.enabled:
            return
        unknown = sorted(set(self.trajectory_modes) - _TRAJECTORY_MODES)
        if unknown:
            raise StudyConfigError(f"attribution.trajectory_modes has unknown mode(s): {unknown}.")
        if not self.trajectory_modes:
            raise StudyConfigError("attribution.trajectory_modes must contain at least one mode.")
        if not self.phases:
            raise StudyConfigError("attribution.phases must contain at least one phase.")
        if self.effect_sizes is not None:
            if not self.effect_sizes:
                raise StudyConfigError("attribution.effect_sizes must be null or a non-empty sequence.")
            for value in self.effect_sizes:
                if float(value) < 0:
                    raise StudyConfigError(f"attribution.effect_sizes must be non-negative; got {value}.")
        if self.bootstrap_replicates < 0:
            raise StudyConfigError("attribution.bootstrap_replicates must be non-negative.")
        if self.top_k < 1:
            raise StudyConfigError("attribution.top_k must be at least 1.")
        if not (self.zero_tolerance > 0):
            raise StudyConfigError("attribution.zero_tolerance must be positive.")

    def settings(self) -> AttributionDiagnosticSettings:
        """Effective evaluation-level settings implied by this selector."""

        return AttributionDiagnosticSettings(
            enabled=True,
            bootstrap_replicates=int(self.bootstrap_replicates),
            bootstrap_seed=None if self.bootstrap_seed is None else int(self.bootstrap_seed),
            top_k=int(self.top_k),
            zero_tolerance=float(self.zero_tolerance),
        )

    def selects(self, *, phase: str, trajectory_mode: str | None, effect_size: float | None) -> bool:
        """Whether a cell with this identity is eligible for attribution."""

        if not self.enabled:
            return False
        if phase not in self.phases:
            return False
        if trajectory_mode is None or trajectory_mode not in self.trajectory_modes:
            return False
        if effect_size is None:
            return False
        value = float(effect_size)
        if self.nonzero_effects_only and value == 0.0:
            return False
        if self.effect_sizes is not None:
            return any(abs(value - float(candidate)) < 1e-12 for candidate in self.effect_sizes)
        return True


@dataclass(frozen=True)
class MatchedSeedPolicy:
    """Opt-in matched generator seeds shared across primary cells.

    When enabled, every primary power cell draws its generator seed from a
    shared ``(seed family, replicate index)`` key instead of its own cell id, so
    requested-effect comparisons are paired at the generated-reference level.
    Negative-control and OFAT cells keep their own seed families so they stay
    independent. The policy version enters the parameter signature.
    """

    enabled: bool = False
    version: int = 1
    primary_family: str = "primary"
    shared_zero_effect_anchor: bool = True

    def __post_init__(self) -> None:
        if not self.enabled:
            return
        if self.version < 1:
            raise StudyConfigError("matched_seeds.version must be at least 1.")
        if not str(self.primary_family).strip():
            raise StudyConfigError("matched_seeds.primary_family must be a non-empty string.")


@dataclass(frozen=True)
class GateRule:
    """One predeclared mode/statistic gate role."""

    trajectory_mode: str
    statistic: str
    role: str
    min_power_at_top: float | None = None

    def __post_init__(self) -> None:
        if self.trajectory_mode not in _TRAJECTORY_MODES:
            raise StudyConfigError(f"gate rule trajectory_mode {self.trajectory_mode!r} is unknown.")
        if self.statistic not in _STATISTICS:
            raise StudyConfigError(f"gate rule statistic {self.statistic!r} is unknown.")
        if self.role not in _GATE_ROLES:
            raise StudyConfigError(f"gate rule role {self.role!r} is unknown; expected one of {sorted(_GATE_ROLES)}.")

    @property
    def name(self) -> str:
        return f"{self.role}[{self.trajectory_mode},{self.statistic}]"


@dataclass(frozen=True)
class Phase4GateConfig:
    """Every parameter the Phase 4 gate consumes.

    Summary and report code reads all thresholds from here and hard-codes none
    of its own, so a gate can be re-specified by editing the committed config.
    """

    enabled: bool = False
    alpha: float = 0.05
    control_se_tolerance: float = 2.0
    monotonicity_se_tolerance: float = 2.0
    min_power_at_top: float = 0.80
    confirmation_se_threshold: float = 1.0
    max_marginal_exceedances: int = 1
    control_modes: tuple[str, ...] = ("none", "translation")
    rules: tuple[GateRule, ...] = ()
    require_complete_records: bool = True

    def __post_init__(self) -> None:
        if not self.enabled:
            return
        if not (0 < self.alpha < 1):
            raise StudyConfigError("acceptance.gate.alpha must be between 0 and 1.")
        for name in ("control_se_tolerance", "monotonicity_se_tolerance", "confirmation_se_threshold"):
            value = float(getattr(self, name))
            if value < 0:
                raise StudyConfigError(f"acceptance.gate.{name} must be non-negative.")
        if not (0 <= self.min_power_at_top <= 1):
            raise StudyConfigError("acceptance.gate.min_power_at_top must be between 0 and 1.")
        if self.max_marginal_exceedances < 0:
            raise StudyConfigError("acceptance.gate.max_marginal_exceedances must be non-negative.")
        unknown = sorted(set(self.control_modes) - _TRAJECTORY_MODES)
        if unknown:
            raise StudyConfigError(f"acceptance.gate.control_modes has unknown mode(s): {unknown}.")
        if not self.control_modes:
            raise StudyConfigError("acceptance.gate.control_modes must contain at least one mode.")
        if not self.rules:
            raise StudyConfigError("acceptance.gate.rules must declare every mandatory and descriptive pair.")
        if not any(rule.role == "mandatory_power" for rule in self.rules):
            raise StudyConfigError("acceptance.gate.rules must declare at least one 'mandatory_power' rule.")
        seen: set[tuple[str, str]] = set()
        for rule in self.rules:
            key = (rule.trajectory_mode, rule.statistic)
            if key in seen:
                raise StudyConfigError(f"acceptance.gate.rules has a duplicate entry for {key}.")
            seen.add(key)

    def power_floor(self, rule: GateRule) -> float:
        """Minimum power at the top effect for ``rule``."""

        return float(self.min_power_at_top if rule.min_power_at_top is None else rule.min_power_at_top)


@dataclass(frozen=True)
class DesignGrid:
    """Design grid: crossed axis values, or an explicit list of design points.

    Unlike ``StudyConfig.axes`` (one factor at a time off the baseline), the
    design grid enumerates the full power grid — one zero-effect anchor plus
    every (mode × nonzero effect) — at each *design point*. In the crossed form
    (``axes``) every combination of the declared values is a point, and every
    axis MUST list the baseline value. In the explicit form (``points``) the
    points are listed directly, every point names the same axes, and the
    baseline point MUST be one of them; ``axes`` is then derived from the
    points (axes in the first point's order, values in order of first
    appearance) so every reader of ``axes`` works unchanged. Either way the
    baseline column is served by the primary grid and never re-emitted.
    """

    axes: Mapping[str, tuple[Any, ...]] = field(default_factory=dict)
    points: tuple[Mapping[str, Any], ...] = ()

    def __post_init__(self) -> None:
        if self.points:
            derived = _axes_from_points(self.points)
            if self.axes and _to_jsonable(dict(self.axes)) != _to_jsonable(derived):
                raise StudyConfigError(
                    "design_grid declares both 'axes' and 'points'; use exactly one form."
                )
            object.__setattr__(self, "axes", derived)
            return
        for axis, values in self.axes.items():
            _validate_axis_namespace(axis, label="design_grid.axes")
            if not values:
                raise StudyConfigError(f"design_grid.axes[{axis!r}] must list at least one value.")

    @property
    def enabled(self) -> bool:
        return bool(self.axes)

    @property
    def explicit(self) -> bool:
        return bool(self.points)


def _axes_from_points(points: Sequence[Mapping[str, Any]]) -> dict[str, tuple[Any, ...]]:
    """Validate explicit design points and derive the per-axis value lists."""

    first = points[0]
    if not isinstance(first, Mapping) or not first:
        raise StudyConfigError("design_grid.points[0] must be a non-empty mapping of axis to value.")
    for axis in first:
        _validate_axis_namespace(axis, label="design_grid.points")
    expected = set(first)
    for index, point in enumerate(points):
        if not isinstance(point, Mapping):
            raise StudyConfigError(f"design_grid.points[{index}] must be a mapping of axis to value.")
        missing = sorted(expected - set(point))
        extra = sorted(set(point) - expected)
        if missing or extra:
            detail = "; ".join(
                part
                for part in (
                    f"missing axis(es) {missing}" if missing else "",
                    f"unexpected axis(es) {extra}" if extra else "",
                )
                if part
            )
            raise StudyConfigError(
                f"design_grid.points[{index}] {_to_jsonable(dict(point))} must name exactly the axes of "
                f"design_grid.points[0]: {detail}."
            )
        for other_index in range(index):
            other = points[other_index]
            if all(_same_axis_value(point[axis], other[axis]) for axis in first):
                raise StudyConfigError(
                    f"design_grid.points[{index}] {_display_point(point)} duplicates "
                    f"design_grid.points[{other_index}].{_estimator_hint(point)}"
                )
    axes: dict[str, tuple[Any, ...]] = {}
    for axis in first:
        values: list[Any] = []
        for point in points:
            if not any(_same_axis_value(point[axis], seen) for seen in values):
                values.append(point[axis])
        axes[axis] = tuple(values)
    return axes


@dataclass(frozen=True)
class DesignPointDecisionRule:
    """Predeclared rule choosing the study design point from a design grid.

    Per design point the target statistic's rejection rate at the largest
    enumerated effect is compared with ``min_power_at_top``; a point *meets* the
    floor when ``rate − confirmation_se_threshold·SE ≥ floor``, is *marginal*
    when only the point estimate clears it, and *fails* otherwise. The chosen
    point is the first confirmed one in ``prefer`` order (each named axis
    ascending, then the remaining design axes in declared order). The rule is
    advisory: it never feeds the Phase 4 gate or the acceptance targets.
    """

    trajectory_mode: str
    statistic: str
    min_power_at_top: float
    prefer: tuple[str, ...]
    confirmation_se_threshold: float = 1.0
    name: str = "design_point"
    kind: str = field(default="design_point", init=False)

    def __post_init__(self) -> None:
        if self.trajectory_mode not in _TRAJECTORY_MODES:
            raise StudyConfigError(
                f"acceptance.design_point.trajectory_mode {self.trajectory_mode!r} is unknown."
            )
        if self.statistic not in _STATISTICS:
            raise StudyConfigError(f"acceptance.design_point.statistic {self.statistic!r} is unknown.")
        if not (0 <= self.min_power_at_top <= 1):
            raise StudyConfigError("acceptance.design_point.min_power_at_top must be between 0 and 1.")
        if self.confirmation_se_threshold < 0:
            raise StudyConfigError(
                "acceptance.design_point.confirmation_se_threshold must be non-negative."
            )
        if not self.prefer:
            raise StudyConfigError("acceptance.design_point.prefer must name at least one design axis.")
        for axis in self.prefer:
            _validate_axis_namespace(axis, label="acceptance.design_point.prefer")


@dataclass(frozen=True)
class StatisticPair:
    """One (trajectory mode, statistic) pair named by a decision rule."""

    trajectory_mode: str
    statistic: str

    def __post_init__(self) -> None:
        if self.trajectory_mode not in _TRAJECTORY_MODES:
            raise StudyConfigError(f"trajectory_mode {self.trajectory_mode!r} is unknown.")
        if self.statistic not in _STATISTICS:
            raise StudyConfigError(f"statistic {self.statistic!r} is unknown.")

    @property
    def label(self) -> str:
        return f"{self.trajectory_mode}/{self.statistic}"


@dataclass(frozen=True)
class TypeIBound:
    """Anchor Type I bound: ``rate ≤ alpha + se_tolerance · sqrt(alpha(1−alpha)/n)``."""

    alpha: float = 0.05
    se_tolerance: float = 2.0

    def __post_init__(self) -> None:
        if not (0 < self.alpha < 1):
            raise StudyConfigError("type_i_bound.alpha must be between 0 and 1.")
        if self.se_tolerance < 0:
            raise StudyConfigError("type_i_bound.se_tolerance must be non-negative.")


#: JSON/`null` baseline value of the retained-rank axis: the key is absent and
#: the evaluator runs its stage-supervised double cross-validation.
RANK_AXIS_BASELINE = None


@dataclass(frozen=True)
class RankDecisionRule:
    """Predeclared, advisory retained-rank rule evaluated against the CV column.

    ``axis`` names a design-grid axis over the nested evaluation parameter
    ``evaluation.integration_params.forced_components`` whose values include
    ``null`` (the cross-validated reference column). A forced-rank column
    *qualifies* when (a) its ``target`` power at the top effect beats the
    reference by more than ``gain_se_multiplier`` pooled MC standard errors,
    (b) its zero-effect anchor is within ``type_i_bound`` on every statistic,
    and (c) no ``protected`` power at the top effect falls below the reference
    by more than ``loss_se_multiplier`` pooled SEs. The verdict is ``keep_cv``
    when no column qualifies, else ``adopt_fixed_rank`` with the smallest
    qualifying rank (parsimony is predeclared). It never feeds the Phase 4 gate
    or the acceptance targets.
    """

    axis: str
    target: StatisticPair
    protected: tuple[StatisticPair, ...]
    type_i_bound: TypeIBound = field(default_factory=TypeIBound)
    gain_se_multiplier: float = 2.0
    loss_se_multiplier: float = 2.0
    name: str = "rank_decision"
    kind: str = field(default="rank_decision", init=False)

    def __post_init__(self) -> None:
        _validate_axis_namespace(self.axis, label="acceptance.rank_decision.axis")
        _, path = _split_axis(self.axis)
        if len(path) != 2:
            raise StudyConfigError(
                f"acceptance.rank_decision.axis {self.axis!r} must name a nested evaluation "
                "integration parameter (evaluation.integration_params.<key>)."
            )
        if not self.protected:
            raise StudyConfigError("acceptance.rank_decision.protected must name at least one pair.")
        if self.target in self.protected:
            raise StudyConfigError(
                f"acceptance.rank_decision.target {self.target.label} cannot also be protected."
            )
        for name in ("gain_se_multiplier", "loss_se_multiplier"):
            if float(getattr(self, name)) < 0:
                raise StudyConfigError(f"acceptance.rank_decision.{name} must be non-negative.")


#: Version of the anchor Type I decision rule; written to ``anchor_decision.json``.
ANCHOR_DECISION_RULE_VERSION = 1


@dataclass(frozen=True)
class AnchorDecisionRule:
    """Predeclared, advisory Type I rule over every design point's zero-effect anchor.

    For each design point (the baseline column included) and each named
    statistic, the anchor's rejection rate over its ``R`` available replicates
    is compared with ``alpha + se_multiplier · sqrt(alpha(1 − alpha)/R)``: at or
    below it the statistic *holds*, above it *exceeds*. A point's verdict is
    ``size_holds`` when every named statistic holds. It never feeds the Phase 4
    gate or the acceptance targets.
    """

    alpha: float = 0.05
    se_multiplier: float = 2.0
    statistics: tuple[str, ...] = ("delta", "angle", "shape")
    name: str = "anchor_decision"
    kind: str = field(default="anchor_decision", init=False)

    def __post_init__(self) -> None:
        if not (0 < self.alpha < 1):
            raise StudyConfigError("acceptance.anchor_decision.alpha must be between 0 and 1.")
        if not (self.se_multiplier > 0):
            raise StudyConfigError("acceptance.anchor_decision.se_multiplier must be positive.")
        if not self.statistics:
            raise StudyConfigError("acceptance.anchor_decision.statistics must name at least one statistic.")
        unknown = sorted(set(self.statistics) - _STATISTICS)
        if unknown:
            raise StudyConfigError(
                f"acceptance.anchor_decision.statistics has unknown statistic(s): {unknown}; "
                f"expected a subset of {sorted(_STATISTICS)}."
            )
        if len(set(self.statistics)) != len(self.statistics):
            raise StudyConfigError("acceptance.anchor_decision.statistics lists a statistic twice.")


#: Version of the column decision rule; written to ``column_decision.json``.
COLUMN_DECISION_RULE_VERSION = 1

#: Target criterion kinds of the column decision rule.
COLUMN_DECISION_CRITERIA: tuple[str, ...] = ("gain", "no_loss")


@dataclass(frozen=True)
class ColumnDecisionRule:
    """Predeclared, advisory rule judging one candidate design column against the baseline.

    ``candidate`` names design-axis coordinates (a subset of the design axes)
    that resolve to exactly one non-baseline design point. At the largest
    enumerated effect of each named mode, the candidate column qualifies when
    (a) its ``target`` rate beats the baseline's by more than
    ``gain_se_multiplier`` pooled MC SEs (``criterion="gain"``) or falls below
    it by at most ``loss_se_multiplier`` pooled SEs (``criterion="no_loss"``),
    (b) its zero-effect anchor is within ``type_i_bound`` on every statistic,
    (c) every mandatory ``controls`` pair is within the same one-sided bound,
    and (d) no ``protected`` pair falls below the baseline by more than
    ``loss_se_multiplier`` pooled SEs. It never feeds the Phase 4 gate or the
    acceptance targets.
    """

    candidate: Mapping[str, Any]
    target: StatisticPair
    criterion: str
    controls: tuple[StatisticPair, ...]
    protected: tuple[StatisticPair, ...] = ()
    type_i_bound: TypeIBound = field(default_factory=TypeIBound)
    gain_se_multiplier: float = 2.0
    loss_se_multiplier: float = 2.0
    name: str = "column_decision"
    kind: str = field(default="column_decision", init=False)

    def __post_init__(self) -> None:
        if not self.candidate:
            raise StudyConfigError("acceptance.column_decision.candidate must name at least one design axis.")
        for axis in self.candidate:
            _validate_axis_namespace(axis, label="acceptance.column_decision.candidate")
        if self.criterion not in COLUMN_DECISION_CRITERIA:
            raise StudyConfigError(
                f"acceptance.column_decision.criterion {self.criterion!r} is unknown; "
                f"expected one of {list(COLUMN_DECISION_CRITERIA)}."
            )
        if not self.controls:
            raise StudyConfigError("acceptance.column_decision.controls must name at least one pair.")
        for role, pairs in (("controls", self.controls), ("protected", self.protected)):
            if self.target in pairs:
                raise StudyConfigError(
                    f"acceptance.column_decision.target {self.target.label} cannot also be in {role}."
                )
            if len(set(pairs)) != len(pairs):
                raise StudyConfigError(f"acceptance.column_decision.{role} lists a pair twice.")
        for name in ("gain_se_multiplier", "loss_se_multiplier"):
            if float(getattr(self, name)) < 0:
                raise StudyConfigError(f"acceptance.column_decision.{name} must be non-negative.")

    @property
    def pairs(self) -> tuple[StatisticPair, ...]:
        """Target, controls and protected pairs, in that order."""

        return (self.target, *self.controls, *self.protected)


@dataclass(frozen=True)
class ReportContract:
    """Declared reporting and execution rules the findings report must honor.

    The block is optional; when absent, reporting and the shard runner behave
    exactly as before it existed. When declared it is *enforced*, not merely
    recorded: ``driver_component`` selects the one attribution component that
    driver tables and the attribution figure present (``phase4_attribution.csv``
    keeps every component); ``cross_replicate_driver_agreement`` is fixed to
    ``"descriptive"`` because the study re-draws the driver set every replicate,
    so cross-replicate top-k Jaccard and sign agreement describe sampling
    variation and can never support a stability claim (``"claim"`` is not a
    legal value); ``n_jobs_override`` says whether ``run_study_shard.py`` refuses
    (``"forbid"``) or warns about (``"warn"``, the historical behavior) a
    ``--n-jobs`` that differs from ``evaluation.n_jobs``. The report echoes the
    resolved contract to ``report/report_contract.json``.
    """

    driver_component: str
    cross_replicate_driver_agreement: str = _DRIVER_AGREEMENT_DESCRIPTIVE
    n_jobs_override: str = "warn"

    def __post_init__(self) -> None:
        if self.driver_component not in _ATTRIBUTION_COMPONENTS:
            raise StudyConfigError(
                f"report_contract.driver_component {self.driver_component!r} is unknown; "
                f"expected one of {list(_ATTRIBUTION_COMPONENTS)}."
            )
        if self.cross_replicate_driver_agreement != _DRIVER_AGREEMENT_DESCRIPTIVE:
            raise StudyConfigError(
                "report_contract.cross_replicate_driver_agreement must be "
                f"{_DRIVER_AGREEMENT_DESCRIPTIVE!r}; got {self.cross_replicate_driver_agreement!r}. "
                "Cross-replicate driver agreement is never claimable: the driver set is "
                "re-drawn every replicate."
            )
        if self.n_jobs_override not in _N_JOBS_OVERRIDE_POLICIES:
            raise StudyConfigError(
                f"report_contract.n_jobs_override {self.n_jobs_override!r} is unknown; "
                f"expected one of {list(_N_JOBS_OVERRIDE_POLICIES)}."
            )

    @property
    def forbids_n_jobs_override(self) -> bool:
        return self.n_jobs_override == "forbid"


@dataclass(frozen=True)
class AcceptanceTargets:
    """Collection of pre-specified acceptance targets."""

    type_i: tuple[TypeIControlTarget, ...] = ()
    power: tuple[PowerMonotonicityTarget, ...] = ()
    specificity: tuple[SpecificityTarget, ...] = ()
    gate: Phase4GateConfig = field(default_factory=Phase4GateConfig)
    design_point: DesignPointDecisionRule | None = None
    rank_decision: RankDecisionRule | None = None
    anchor_decision: AnchorDecisionRule | None = None
    column_decision: ColumnDecisionRule | None = None


@dataclass(frozen=True)
class StudyConfig:
    """Declarative study definition."""

    generator: SemiSyntheticTrajectoryParams
    evaluation: SimulationEvaluationParams
    trajectory_modes: tuple[str, ...]
    effect_sizes: tuple[float, ...]
    axes: Mapping[str, tuple[Any, ...]] = field(default_factory=dict)
    design_grid: DesignGrid = field(default_factory=DesignGrid)
    n_replicates: int = 1
    base_seed: int = 0
    alpha: float = 0.05
    acceptance: AcceptanceTargets = field(default_factory=AcceptanceTargets)
    attribution: AttributionSelector = field(default_factory=AttributionSelector)
    matched_seeds: MatchedSeedPolicy = field(default_factory=MatchedSeedPolicy)
    report_contract: ReportContract | None = None
    metadata: Mapping[str, Any] = field(default_factory=dict)

    def __post_init__(self) -> None:
        if self.n_replicates < 0:
            raise StudyConfigError("n_replicates must be non-negative.")
        if not self.trajectory_modes:
            raise StudyConfigError("trajectory_modes must contain at least one mode.")
        unknown = sorted(set(self.trajectory_modes) - _TRAJECTORY_MODES)
        if unknown:
            raise StudyConfigError(f"Unknown trajectory mode(s): {unknown}.")
        if not self.effect_sizes:
            raise StudyConfigError("effect_sizes must contain at least one value.")
        for value in self.effect_sizes:
            if float(value) < 0:
                raise StudyConfigError(f"effect_sizes must be non-negative; got {value}.")
        for axis in self.axes:
            _validate_axis_namespace(axis)
        self._validate_design_grid()
        self._validate_rank_decision()
        self._validate_column_decision()
        if not (0 < self.alpha < 1):
            raise StudyConfigError("alpha must be between 0 and 1.")
        if self.attribution.enabled:
            if self.evaluation.integration_method != "pls":
                raise StudyConfigError(
                    "attribution diagnostics require evaluation.integration_method='pls'; "
                    f"got {self.evaluation.integration_method!r}."
                )
            unselectable = sorted(set(self.attribution.trajectory_modes) - set(self.trajectory_modes))
            if unselectable:
                raise StudyConfigError(
                    f"attribution.trajectory_modes select mode(s) absent from trajectory_modes: {unselectable}."
                )
            if self.attribution.effect_sizes is not None:
                declared = [float(value) for value in self.effect_sizes]
                missing = [
                    float(value)
                    for value in self.attribution.effect_sizes
                    if not any(abs(float(value) - candidate) < 1e-12 for candidate in declared)
                ]
                if missing:
                    raise StudyConfigError(
                        f"attribution.effect_sizes select effect(s) absent from effect_sizes: {missing}."
                    )


    def _validate_design_grid(self) -> None:
        grid = self.design_grid
        for axis, values in grid.axes.items():
            if axis in self.axes:
                raise StudyConfigError(
                    f"design_grid.axes[{axis!r}] is also declared under axes; an axis is either "
                    "crossed in the design grid or varied one factor at a time, not both."
                )
            if grid.explicit:
                continue
            baseline = self.axis_baseline_value(axis)
            if not any(_same_axis_value(value, baseline) for value in values):
                raise StudyConfigError(
                    f"design_grid.axes[{axis!r}] must include the baseline value "
                    f"{_display_axis_value(axis, baseline)!r} so the primary power cells are one of the "
                    "design points."
                )
        if grid.explicit and not any(self.is_baseline_design_point(point) for point in grid.points):
            raise StudyConfigError(
                "design_grid.points must include the baseline point "
                f"{_display_point(self.baseline_design_point())} so the primary power cells are one of "
                f"the design points.{_estimator_hint(self.baseline_design_point())}"
            )
        rule = self.acceptance.design_point
        if rule is None:
            return
        if not grid.enabled:
            raise StudyConfigError("acceptance.design_point requires a non-empty design_grid.")
        unknown = [axis for axis in rule.prefer if axis not in grid.axes]
        if unknown:
            raise StudyConfigError(
                f"acceptance.design_point.prefer names axis(es) absent from design_grid.axes: {unknown}."
            )
        if rule.trajectory_mode not in self.trajectory_modes:
            raise StudyConfigError(
                f"acceptance.design_point.trajectory_mode {rule.trajectory_mode!r} is absent from "
                "trajectory_modes."
            )

    def _validate_rank_decision(self) -> None:
        rule = self.acceptance.rank_decision
        if rule is None:
            return
        grid = self.design_grid
        if rule.axis not in grid.axes:
            raise StudyConfigError(
                f"acceptance.rank_decision.axis {rule.axis!r} is not declared under design_grid.axes."
            )
        if not any(_same_axis_value(value, RANK_AXIS_BASELINE) for value in grid.axes[rule.axis]):
            raise StudyConfigError(
                f"design_grid.axes[{rule.axis!r}] must include null (the cross-validated reference "
                "column) for acceptance.rank_decision to compare against."
            )
        for pair in (rule.target, *rule.protected):
            if pair.trajectory_mode not in self.trajectory_modes:
                raise StudyConfigError(
                    f"acceptance.rank_decision names mode {pair.trajectory_mode!r} ({pair.label}), "
                    "which is absent from trajectory_modes."
                )

    def _validate_column_decision(self) -> None:
        rule = self.acceptance.column_decision
        if rule is None:
            return
        grid = self.design_grid
        if not grid.enabled:
            raise StudyConfigError("acceptance.column_decision requires a non-empty design_grid.")
        candidate = _to_jsonable(dict(rule.candidate))
        unknown = [axis for axis in rule.candidate if axis not in grid.axes]
        if unknown:
            raise StudyConfigError(
                f"acceptance.column_decision.candidate {candidate} names axis(es) absent from "
                f"design_grid.axes: {unknown}."
            )
        for pair in rule.pairs:
            if pair.trajectory_mode not in self.trajectory_modes:
                raise StudyConfigError(
                    f"acceptance.column_decision names mode {pair.trajectory_mode!r} ({pair.label}), "
                    "which is absent from trajectory_modes."
                )
        matches = self.matching_design_points(rule.candidate)
        if not matches:
            raise StudyConfigError(
                f"acceptance.column_decision.candidate {candidate} matches no design point."
            )
        if len(matches) > 1:
            raise StudyConfigError(
                f"acceptance.column_decision.candidate {candidate} matches {len(matches)} design points "
                f"{_to_jsonable(matches)}; name enough axes to select exactly one."
            )
        if self.is_baseline_design_point(matches[0]):
            raise StudyConfigError(
                f"acceptance.column_decision.candidate {candidate} matches the baseline point; the "
                "candidate must be a non-baseline design point."
            )

    def matching_design_points(self, coordinates: Mapping[str, Any]) -> list[dict[str, Any]]:
        """Design points whose value on every axis named in ``coordinates`` equals it."""

        return [
            point
            for point in self.design_points()
            if all(axis in point and _same_axis_value(point[axis], value) for axis, value in coordinates.items())
        ]

    def axis_baseline_value(self, axis: str) -> Any:
        """Baseline value of a namespaced ``generator.``/``evaluation.`` axis.

        For the nested form ``evaluation.integration_params.<key>`` this is the
        mapping's current value, or ``None`` when the key is absent — which is
        how the cross-validated column of a retained-rank ladder is addressed.
        """

        namespace, path = _split_axis_or_raise(axis)
        if len(path) == 2:
            return dict(self.evaluation.integration_params).get(path[1])
        (field_name,) = path
        target = self.generator if namespace == "generator" else self.evaluation
        if not hasattr(target, field_name):
            raise StudyConfigError(f"axis {axis!r} names an unknown {namespace} field {field_name!r}.")
        return getattr(target, field_name)

    def design_points(self) -> list[dict[str, Any]]:
        """Every design point as ``{axis: value}``, baseline included, in declared order."""

        if self.design_grid.explicit:
            return [dict(point) for point in self.design_grid.points]
        axes = list(self.design_grid.axes.items())
        if not axes:
            return []
        points: list[dict[str, Any]] = []
        for combination in itertools.product(*(values for _, values in axes)):
            points.append({axis: value for (axis, _), value in zip(axes, combination)})
        return points

    def baseline_design_point(self) -> dict[str, Any]:
        """The design point occupied by the primary power cells."""

        return {axis: self.axis_baseline_value(axis) for axis in self.design_grid.axes}

    def is_baseline_design_point(self, point: Mapping[str, Any]) -> bool:
        baseline = self.baseline_design_point()
        return all(_same_axis_value(point[axis], baseline[axis]) for axis in baseline)


def _same_axis_value(left: Any, right: Any) -> bool:
    return _to_jsonable(left) == _to_jsonable(right)


def load_study_config(path: str | Path) -> StudyConfig:
    """Load a study configuration from a YAML or JSON file."""

    path = Path(path)
    if not path.exists():
        raise StudyConfigError(f"Study configuration not found: {path}.")
    text = path.read_text(encoding="utf-8")
    data = _parse_text(text, path.suffix.lower())
    if not isinstance(data, Mapping):
        raise StudyConfigError(f"Study configuration root must be a mapping, got {type(data).__name__}.")
    return _build_config(data)


def dump_study_config(config: StudyConfig, path: str | Path) -> None:
    """Write a study configuration to JSON (always JSON for portability)."""

    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    payload = _config_to_dict(config)
    path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")


def _parse_text(text: str, suffix: str) -> Any:
    if suffix in {".yaml", ".yml"}:
        try:
            import yaml  # type: ignore[import-untyped,unused-ignore]
        except ImportError as exc:  # pragma: no cover - exercised only without PyYAML
            raise StudyConfigError(
                "PyYAML is required to load YAML study configurations. "
                "Install pyyaml or rewrite the config as JSON."
            ) from exc
        return yaml.safe_load(text)
    if suffix == ".json":
        return json.loads(text)
    # try YAML first, fall back to JSON
    try:
        import yaml  # type: ignore[import-untyped,unused-ignore]

        return yaml.safe_load(text)
    except ImportError:
        return json.loads(text)


def _build_config(data: Mapping[str, Any]) -> StudyConfig:
    required = {"generator", "evaluation", "trajectory_modes", "effect_sizes"}
    missing = sorted(required - set(data))
    if missing:
        raise StudyConfigError(f"Study configuration is missing required field(s): {missing}.")
    unknown = sorted(set(data) - _ROOT_KEYS)
    if unknown:
        raise StudyConfigError(
            f"Study configuration has unknown root key(s): {unknown}; "
            f"expected only {sorted(_ROOT_KEYS)}."
        )

    generator = _build_generator(data["generator"])
    evaluation = _build_evaluation(data.get("evaluation") or {})
    trajectory_modes = tuple(str(mode) for mode in data["trajectory_modes"])
    effect_sizes = tuple(float(value) for value in data["effect_sizes"])
    axes = _build_axes(data.get("axes") or {})
    design_grid = _build_design_grid(data.get("design_grid") or {})
    n_replicates = int(data.get("n_replicates", 1))
    base_seed = int(data.get("base_seed", 0))
    alpha = float(data.get("alpha", 0.05))
    acceptance = _build_acceptance(data.get("acceptance") or {})
    attribution = _build_attribution(data.get("attribution") or {})
    matched_seeds = _build_matched_seeds(data.get("matched_seeds") or {})
    report_contract = _build_report_contract(data.get("report_contract"))
    metadata = dict(data.get("metadata") or {})
    return StudyConfig(
        generator=generator,
        evaluation=evaluation,
        trajectory_modes=trajectory_modes,
        effect_sizes=effect_sizes,
        axes=axes,
        design_grid=design_grid,
        n_replicates=n_replicates,
        base_seed=base_seed,
        alpha=alpha,
        acceptance=acceptance,
        attribution=attribution,
        matched_seeds=matched_seeds,
        report_contract=report_contract,
        metadata=metadata,
    )


def _build_generator(raw: Mapping[str, Any]) -> SemiSyntheticTrajectoryParams:
    if "seed" not in raw:
        raise StudyConfigError("generator.seed is required.")
    field_names = {f.name for f in dataclasses.fields(SemiSyntheticTrajectoryParams)}
    unknown = sorted(set(raw) - field_names)
    if unknown:
        raise StudyConfigError(f"generator has unknown field(s): {unknown}.")
    kwargs: dict[str, Any] = {}
    for f in dataclasses.fields(SemiSyntheticTrajectoryParams):
        if f.name not in raw:
            continue
        value = raw[f.name]
        if f.name == "group_labels" and value is not None:
            value = tuple(str(v) for v in value)
        if f.name == "stage_sample_prop" and value is not None:
            value = tuple(float(v) for v in value)
        if f.name == "group_stage_sizes":
            value = _normalize_group_stage_sizes(value)
        kwargs[f.name] = value
    # Reject an unknown construction variant at load time rather than letting
    # every work unit fail at generation. The selectable set is the generator's
    # own (`all`, `extremes`, `joint`); a profile that omits the key resolves to
    # the generator default (`all`), so historical configs are unchanged.
    magnitude_kind = kwargs.get("magnitude_kind")
    if magnitude_kind is not None and magnitude_kind not in _MAGNITUDE_KINDS:
        raise StudyConfigError(
            f"generator.magnitude_kind {magnitude_kind!r} is unknown; "
            f"expected one of {sorted(_MAGNITUDE_KINDS)}."
        )
    return SemiSyntheticTrajectoryParams(**kwargs)


def _build_evaluation(raw: Mapping[str, Any]) -> SimulationEvaluationParams:
    field_names = {f.name for f in dataclasses.fields(SimulationEvaluationParams)}
    unknown = sorted(set(raw) - field_names)
    if unknown:
        raise StudyConfigError(f"evaluation has unknown field(s): {unknown}.")
    kwargs: dict[str, Any] = {}
    for f in dataclasses.fields(SimulationEvaluationParams):
        if f.name not in raw:
            continue
        value = raw[f.name]
        if f.name == "attribution" and value is not None and not isinstance(value, AttributionDiagnosticSettings):
            value = _build_attribution_settings(value)
        if f.name == "orientation_estimator":
            value = _normalize_orientation_estimator(value)
        kwargs[f.name] = value
    return SimulationEvaluationParams(**kwargs)


def _build_attribution_settings(raw: Mapping[str, Any]) -> AttributionDiagnosticSettings:
    field_names = {f.name for f in dataclasses.fields(AttributionDiagnosticSettings)}
    unknown = sorted(set(raw) - field_names)
    if unknown:
        raise StudyConfigError(f"evaluation.attribution has unknown field(s): {unknown}.")
    return AttributionDiagnosticSettings(**{key: raw[key] for key in raw})


def _normalize_group_stage_sizes(value: Any) -> Any:
    """Nested JSON lists → a tuple of int tuples, so equal tables hash equally.

    Malformed shapes pass through unchanged for the generator to reject with
    its own message.
    """

    if value is None or isinstance(value, str | bytes) or not isinstance(value, Sequence):
        return value
    rows: list[Any] = []
    for row in value:
        if isinstance(row, str | bytes) or not isinstance(row, Sequence):
            return value
        rows.append(tuple(row))
    return tuple(rows)


_ORIENTATION_ESTIMATOR_DEFAULT = next(
    f.default for f in dataclasses.fields(SimulationEvaluationParams) if f.name == "orientation_estimator"
)


def _normalize_orientation_estimator(value: Any) -> Any:
    """Config spelling of ``evaluation.orientation_estimator`` → evaluator value.

    ``null`` names the evaluator's default, as it does for the nested
    ``evaluation.integration_params.<key>`` form; the default is read from the
    dataclass so a later default change needs no edit here. ``"pc1"`` names PC1,
    whose evaluator value is the hash-omitted ``None``. ``null`` and ``"pc1"``
    are therefore distinct axis values.
    """

    if value is None:
        return _ORIENTATION_ESTIMATOR_DEFAULT
    return None if value == "pc1" else value


#: Axis values that need the same normalization as their baseline field.
_AXIS_VALUE_NORMALIZERS = {
    "generator.group_stage_sizes": _normalize_group_stage_sizes,
    "evaluation.orientation_estimator": _normalize_orientation_estimator,
}


def _display_axis_value(axis: str, value: Any) -> Any:
    """A normalized axis value in config spelling, for error messages.

    The PC1 estimator is held as ``None`` after normalization; show it as
    ``"pc1"`` so a message never reads as ``null`` (the evaluator's default).
    """

    if axis == "evaluation.orientation_estimator" and value is None:
        return "pc1"
    return _to_jsonable(value)


def _display_point(point: Mapping[str, Any]) -> dict[str, Any]:
    return {axis: _display_axis_value(axis, value) for axis, value in point.items()}


def _estimator_hint(point: Mapping[str, Any]) -> str:
    if "evaluation.orientation_estimator" not in point:
        return ""
    return (
        " For evaluation.orientation_estimator, null names the evaluator's default "
        f"({_ORIENTATION_ESTIMATOR_DEFAULT!r}) and 'pc1' names PC1; the baseline is the baseline "
        "evaluation's estimator."
    )


def _normalize_axis_values(axis: str, values: Sequence[Any]) -> tuple[Any, ...]:
    normalize = _AXIS_VALUE_NORMALIZERS.get(axis)
    return tuple(values) if normalize is None else tuple(normalize(value) for value in values)


def _build_axes(raw: Mapping[str, Any]) -> Mapping[str, tuple[Any, ...]]:
    axes: dict[str, tuple[Any, ...]] = {}
    for axis, values in raw.items():
        _validate_axis_namespace(axis)
        if not isinstance(values, Sequence) or isinstance(values, str | bytes):
            raise StudyConfigError(f"axis {axis!r} values must be a sequence.")
        axes[axis] = _normalize_axis_values(axis, values)
    return axes


def _build_design_grid(raw: Mapping[str, Any]) -> DesignGrid:
    if not raw:
        return DesignGrid()
    unknown = sorted(set(raw) - {"axes", "points"})
    if unknown:
        raise StudyConfigError(f"design_grid has unknown field(s): {unknown}.")
    if "points" in raw:
        if "axes" in raw:
            raise StudyConfigError("design_grid declares both 'axes' and 'points'; use exactly one form.")
        raw_points = raw["points"]
        if isinstance(raw_points, Mapping | str | bytes) or not isinstance(raw_points, Sequence):
            raise StudyConfigError("design_grid.points must be a sequence of mappings.")
        if not raw_points:
            raise StudyConfigError("design_grid.points must list at least one point.")
        points: list[dict[str, Any]] = []
        for index, point in enumerate(raw_points):
            if not isinstance(point, Mapping):
                raise StudyConfigError(f"design_grid.points[{index}] must be a mapping of axis to value.")
            for axis in point:
                _validate_axis_namespace(axis, label=f"design_grid.points[{index}]")
            points.append({axis: _normalize_axis_values(axis, [value])[0] for axis, value in point.items()})
        return DesignGrid(points=tuple(points))
    axes: dict[str, tuple[Any, ...]] = {}
    for axis, values in (raw.get("axes") or {}).items():
        _validate_axis_namespace(axis, label="design_grid.axes")
        if not isinstance(values, Sequence) or isinstance(values, str | bytes):
            raise StudyConfigError(f"design_grid.axes[{axis!r}] values must be a sequence.")
        axes[axis] = _normalize_axis_values(axis, values)
    return DesignGrid(axes=axes)


def _build_design_point_rule(raw: Mapping[str, Any] | None) -> DesignPointDecisionRule | None:
    if not raw:
        return None
    required = {"trajectory_mode", "statistic", "min_power_at_top", "prefer"}
    missing = sorted(required - set(raw))
    if missing:
        raise StudyConfigError(f"acceptance.design_point is missing required field(s): {missing}.")
    unknown = sorted(set(raw) - (required | {"confirmation_se_threshold", "name", "kind"}))
    if unknown:
        raise StudyConfigError(f"acceptance.design_point has unknown field(s): {unknown}.")
    prefer = raw["prefer"]
    if isinstance(prefer, str) or not isinstance(prefer, Sequence):
        raise StudyConfigError("acceptance.design_point.prefer must be a sequence of axis names.")
    return DesignPointDecisionRule(
        trajectory_mode=str(raw["trajectory_mode"]),
        statistic=str(raw["statistic"]),
        min_power_at_top=float(raw["min_power_at_top"]),
        prefer=tuple(str(axis) for axis in prefer),
        confirmation_se_threshold=float(raw.get("confirmation_se_threshold", 1.0)),
        name=str(raw.get("name", "design_point")),
    )


def _build_statistic_pair(raw: Any, *, label: str) -> StatisticPair:
    if not isinstance(raw, Mapping):
        raise StudyConfigError(f"{label} must be a mapping with 'trajectory_mode' and 'statistic'.")
    missing = sorted({"trajectory_mode", "statistic"} - set(raw))
    if missing:
        raise StudyConfigError(f"{label} is missing required field(s): {missing}.")
    unknown = sorted(set(raw) - {"trajectory_mode", "statistic"})
    if unknown:
        raise StudyConfigError(f"{label} has unknown field(s): {unknown}.")
    return StatisticPair(trajectory_mode=str(raw["trajectory_mode"]), statistic=str(raw["statistic"]))


def _build_rank_decision_rule(raw: Mapping[str, Any] | None) -> RankDecisionRule | None:
    if not raw:
        return None
    required = {"axis", "target", "protected"}
    missing = sorted(required - set(raw))
    if missing:
        raise StudyConfigError(f"acceptance.rank_decision is missing required field(s): {missing}.")
    known = required | {"type_i_bound", "gain_se_multiplier", "loss_se_multiplier", "name", "kind"}
    unknown = sorted(set(raw) - known)
    if unknown:
        raise StudyConfigError(f"acceptance.rank_decision has unknown field(s): {unknown}.")
    protected_raw = raw["protected"]
    if isinstance(protected_raw, Mapping) or not isinstance(protected_raw, Sequence):
        raise StudyConfigError("acceptance.rank_decision.protected must be a sequence of pairs.")
    bound_raw = raw.get("type_i_bound") or {}
    if not isinstance(bound_raw, Mapping):
        raise StudyConfigError("acceptance.rank_decision.type_i_bound must be a mapping.")
    unknown_bound = sorted(set(bound_raw) - {"alpha", "se_tolerance"})
    if unknown_bound:
        raise StudyConfigError(f"acceptance.rank_decision.type_i_bound has unknown field(s): {unknown_bound}.")
    return RankDecisionRule(
        axis=str(raw["axis"]),
        target=_build_statistic_pair(raw["target"], label="acceptance.rank_decision.target"),
        protected=tuple(
            _build_statistic_pair(entry, label="acceptance.rank_decision.protected[]")
            for entry in protected_raw
        ),
        type_i_bound=TypeIBound(
            alpha=float(bound_raw.get("alpha", 0.05)),
            se_tolerance=float(bound_raw.get("se_tolerance", 2.0)),
        ),
        gain_se_multiplier=float(raw.get("gain_se_multiplier", 2.0)),
        loss_se_multiplier=float(raw.get("loss_se_multiplier", 2.0)),
        name=str(raw.get("name", "rank_decision")),
    )


def _build_anchor_decision_rule(raw: Mapping[str, Any] | None) -> AnchorDecisionRule | None:
    if raw is None:
        return None
    if not isinstance(raw, Mapping):
        raise StudyConfigError("acceptance.anchor_decision must be a mapping.")
    unknown = sorted(set(raw) - {"alpha", "se_multiplier", "statistics", "name", "kind"})
    if unknown:
        raise StudyConfigError(f"acceptance.anchor_decision has unknown field(s): {unknown}.")
    statistics = raw.get("statistics", ("delta", "angle", "shape"))
    if isinstance(statistics, str) or not isinstance(statistics, Sequence):
        raise StudyConfigError("acceptance.anchor_decision.statistics must be a sequence of statistics.")
    return AnchorDecisionRule(
        alpha=float(raw.get("alpha", 0.05)),
        se_multiplier=float(raw.get("se_multiplier", 2.0)),
        statistics=tuple(str(statistic) for statistic in statistics),
        name=str(raw.get("name", "anchor_decision")),
    )


def _build_pair_list(raw: Any, *, label: str) -> tuple[StatisticPair, ...]:
    if isinstance(raw, Mapping | str | bytes) or not isinstance(raw, Sequence):
        raise StudyConfigError(f"{label} must be a sequence of pairs.")
    return tuple(_build_statistic_pair(entry, label=f"{label}[]") for entry in raw)


def _build_column_decision_rule(raw: Any) -> ColumnDecisionRule | None:
    if raw is None:
        return None
    if not isinstance(raw, Mapping):
        raise StudyConfigError("acceptance.column_decision must be a mapping.")
    required = {"candidate", "target", "criterion", "controls"}
    missing = sorted(required - set(raw))
    if missing:
        raise StudyConfigError(f"acceptance.column_decision is missing required field(s): {missing}.")
    known = required | {"protected", "type_i_bound", "gain_se_multiplier", "loss_se_multiplier", "name", "kind"}
    unknown = sorted(set(raw) - known)
    if unknown:
        raise StudyConfigError(f"acceptance.column_decision has unknown field(s): {unknown}.")
    candidate_raw = raw["candidate"]
    if not isinstance(candidate_raw, Mapping):
        raise StudyConfigError("acceptance.column_decision.candidate must be a mapping of axis to value.")
    for axis in candidate_raw:
        _validate_axis_namespace(axis, label="acceptance.column_decision.candidate")
    bound_raw = raw.get("type_i_bound") or {}
    if not isinstance(bound_raw, Mapping):
        raise StudyConfigError("acceptance.column_decision.type_i_bound must be a mapping.")
    unknown_bound = sorted(set(bound_raw) - {"alpha", "se_tolerance"})
    if unknown_bound:
        raise StudyConfigError(
            f"acceptance.column_decision.type_i_bound has unknown field(s): {unknown_bound}."
        )
    return ColumnDecisionRule(
        candidate={axis: _normalize_axis_values(axis, [value])[0] for axis, value in candidate_raw.items()},
        target=_build_statistic_pair(raw["target"], label="acceptance.column_decision.target"),
        criterion=str(raw["criterion"]),
        controls=_build_pair_list(raw["controls"], label="acceptance.column_decision.controls"),
        protected=_build_pair_list(raw.get("protected", []), label="acceptance.column_decision.protected"),
        type_i_bound=TypeIBound(
            alpha=float(bound_raw.get("alpha", 0.05)),
            se_tolerance=float(bound_raw.get("se_tolerance", 2.0)),
        ),
        gain_se_multiplier=float(raw.get("gain_se_multiplier", 2.0)),
        loss_se_multiplier=float(raw.get("loss_se_multiplier", 2.0)),
        name=str(raw.get("name", "column_decision")),
    )


def _build_acceptance(raw: Mapping[str, Any]) -> AcceptanceTargets:
    unknown = sorted(
        set(raw)
        - {
            "type_i",
            "power",
            "specificity",
            "gate",
            "design_point",
            "rank_decision",
            "anchor_decision",
            "column_decision",
        }
    )
    if unknown:
        raise StudyConfigError(f"acceptance has unknown block(s): {unknown}.")
    type_i = tuple(
        _build_type_i_target(entry, index) for index, entry in enumerate(raw.get("type_i", []) or [])
    )
    power = tuple(_build_power_target(entry) for entry in raw.get("power", []) or [])
    specificity = tuple(
        _build_specificity_target(entry, index)
        for index, entry in enumerate(raw.get("specificity", []) or [])
    )
    gate = _build_gate(raw.get("gate") or {})
    design_point = _build_design_point_rule(raw.get("design_point"))
    rank_decision = _build_rank_decision_rule(raw.get("rank_decision"))
    anchor_decision = _build_anchor_decision_rule(raw.get("anchor_decision"))
    column_decision = _build_column_decision_rule(raw.get("column_decision"))
    return AcceptanceTargets(
        type_i=type_i,
        power=power,
        specificity=specificity,
        gate=gate,
        design_point=design_point,
        rank_decision=rank_decision,
        anchor_decision=anchor_decision,
        column_decision=column_decision,
    )


def _build_gate(raw: Mapping[str, Any]) -> Phase4GateConfig:
    if not raw:
        return Phase4GateConfig()
    known = {
        "enabled",
        "alpha",
        "control_se_tolerance",
        "monotonicity_se_tolerance",
        "min_power_at_top",
        "confirmation_se_threshold",
        "max_marginal_exceedances",
        "control_modes",
        "rules",
        "require_complete_records",
    }
    unknown = sorted(set(raw) - known)
    if unknown:
        raise StudyConfigError(f"acceptance.gate has unknown field(s): {unknown}.")
    enabled = bool(raw.get("enabled", True))
    if enabled:
        required = {
            "alpha",
            "control_se_tolerance",
            "monotonicity_se_tolerance",
            "min_power_at_top",
            "confirmation_se_threshold",
            "rules",
        }
        missing = sorted(required - set(raw))
        if missing:
            raise StudyConfigError(f"acceptance.gate is missing required parameter(s): {missing}.")
    rules = tuple(_build_gate_rule(entry) for entry in raw.get("rules", []) or [])
    control_modes = (
        tuple(str(mode) for mode in raw["control_modes"])
        if "control_modes" in raw
        else ("none", "translation")
    )
    return Phase4GateConfig(
        enabled=enabled,
        alpha=float(raw.get("alpha", 0.05)),
        control_se_tolerance=float(raw.get("control_se_tolerance", 2.0)),
        monotonicity_se_tolerance=float(raw.get("monotonicity_se_tolerance", 2.0)),
        min_power_at_top=float(raw.get("min_power_at_top", 0.80)),
        confirmation_se_threshold=float(raw.get("confirmation_se_threshold", 1.0)),
        max_marginal_exceedances=int(raw.get("max_marginal_exceedances", 1)),
        control_modes=control_modes,
        rules=rules,
        require_complete_records=bool(raw.get("require_complete_records", True)),
    )


def _build_gate_rule(raw: Mapping[str, Any]) -> GateRule:
    required = {"trajectory_mode", "statistic", "role"}
    missing = sorted(required - set(raw))
    if missing:
        raise StudyConfigError(f"acceptance.gate.rules entry missing field(s): {missing}.")
    unknown = sorted(set(raw) - (required | {"min_power_at_top"}))
    if unknown:
        raise StudyConfigError(f"acceptance.gate.rules entry has unknown field(s): {unknown}.")
    floor = raw.get("min_power_at_top")
    return GateRule(
        trajectory_mode=str(raw["trajectory_mode"]),
        statistic=str(raw["statistic"]),
        role=str(raw["role"]),
        min_power_at_top=None if floor is None else float(floor),
    )


def _build_attribution(raw: Mapping[str, Any]) -> AttributionSelector:
    if not raw:
        return AttributionSelector()
    field_names = {f.name for f in dataclasses.fields(AttributionSelector)}
    unknown = sorted(set(raw) - field_names)
    if unknown:
        raise StudyConfigError(f"attribution has unknown field(s): {unknown}.")
    effect_sizes = raw.get("effect_sizes")
    return AttributionSelector(
        enabled=bool(raw.get("enabled", True)),
        trajectory_modes=tuple(str(mode) for mode in raw.get("trajectory_modes", ("orientation",))),
        phases=tuple(str(phase) for phase in raw.get("phases", ("power_primary",))),
        nonzero_effects_only=bool(raw.get("nonzero_effects_only", True)),
        effect_sizes=None if effect_sizes is None else tuple(float(value) for value in effect_sizes),
        bootstrap_replicates=int(raw.get("bootstrap_replicates", 100)),
        bootstrap_seed=None if raw.get("bootstrap_seed", 0) is None else int(raw.get("bootstrap_seed", 0)),
        top_k=int(raw.get("top_k", 20)),
        zero_tolerance=float(raw.get("zero_tolerance", 1e-12)),
    )


def _build_matched_seeds(raw: Mapping[str, Any]) -> MatchedSeedPolicy:
    if not raw:
        return MatchedSeedPolicy()
    field_names = {f.name for f in dataclasses.fields(MatchedSeedPolicy)}
    unknown = sorted(set(raw) - field_names)
    if unknown:
        raise StudyConfigError(f"matched_seeds has unknown field(s): {unknown}.")
    return MatchedSeedPolicy(
        enabled=bool(raw.get("enabled", True)),
        version=int(raw.get("version", 1)),
        primary_family=str(raw.get("primary_family", "primary")),
        shared_zero_effect_anchor=bool(raw.get("shared_zero_effect_anchor", True)),
    )


def _build_report_contract(raw: Any) -> ReportContract | None:
    if raw is None:
        return None
    if not isinstance(raw, Mapping):
        raise StudyConfigError(f"report_contract must be a mapping, got {type(raw).__name__}.")
    field_names = {f.name for f in dataclasses.fields(ReportContract)}
    unknown = sorted(set(raw) - field_names)
    if unknown:
        raise StudyConfigError(f"report_contract has unknown field(s): {unknown}.")
    if "driver_component" not in raw:
        raise StudyConfigError("report_contract.driver_component is required when the block is declared.")
    return ReportContract(
        driver_component=str(raw["driver_component"]),
        cross_replicate_driver_agreement=str(
            raw.get("cross_replicate_driver_agreement", _DRIVER_AGREEMENT_DESCRIPTIVE)
        ),
        n_jobs_override=str(raw.get("n_jobs_override", "warn")),
    )


def _parse_target_bound(raw: Mapping[str, Any], label: str) -> str:
    """Read a control target's bound form, rejecting anything but the two forms."""

    bound = raw.get("bound", "two_sided")
    if not isinstance(bound, str) or bound not in _TARGET_BOUNDS:
        raise StudyConfigError(
            f"{label}.bound must be one of {sorted(_TARGET_BOUNDS)}; got {bound!r}."
        )
    return str(bound)


def _build_type_i_target(raw: Mapping[str, Any], index: int = 0) -> TypeIControlTarget:
    if "alpha" not in raw:
        raise StudyConfigError("acceptance.type_i entries require 'alpha'.")
    return TypeIControlTarget(
        alpha=float(raw["alpha"]),
        se_tolerance=float(raw.get("se_tolerance", 2.0)),
        name=str(raw.get("name", "type_i_control")),
        bound=_parse_target_bound(raw, f"acceptance.type_i[{index}]"),
    )


def _build_power_target(raw: Mapping[str, Any]) -> PowerMonotonicityTarget:
    required = {"trajectory_mode", "statistic", "min_power_at_top"}
    missing = sorted(required - set(raw))
    if missing:
        raise StudyConfigError(f"acceptance.power entry missing field(s): {missing}.")
    mode = str(raw["trajectory_mode"])
    if mode not in _TRAJECTORY_MODES:
        raise StudyConfigError(f"acceptance.power trajectory_mode {mode!r} is unknown.")
    return PowerMonotonicityTarget(
        trajectory_mode=mode,
        statistic=str(raw["statistic"]),
        min_power_at_top=float(raw["min_power_at_top"]),
        name=str(raw.get("name", f"power[{mode},{raw['statistic']}]")),
    )


def _build_specificity_target(raw: Mapping[str, Any], index: int = 0) -> SpecificityTarget:
    required = {"trajectory_mode", "statistic", "alpha"}
    missing = sorted(required - set(raw))
    if missing:
        raise StudyConfigError(f"acceptance.specificity entry missing field(s): {missing}.")
    mode = str(raw["trajectory_mode"])
    if mode not in _TRAJECTORY_MODES:
        raise StudyConfigError(f"acceptance.specificity trajectory_mode {mode!r} is unknown.")
    return SpecificityTarget(
        trajectory_mode=mode,
        statistic=str(raw["statistic"]),
        alpha=float(raw["alpha"]),
        se_tolerance=float(raw.get("se_tolerance", 2.0)),
        name=str(raw.get("name", f"specificity[{mode},{raw['statistic']}]")),
        bound=_parse_target_bound(raw, f"acceptance.specificity[{index}]"),
    )


def _validate_axis_namespace(axis: str, *, label: str = "axis") -> None:
    """Accept ``generator.<field>``, ``evaluation.<field>``, or the one nested form.

    The nested form is ``evaluation.integration_params.<key>``; any other
    multi-segment axis is rejected by name so a typo cannot address a key the
    evaluator never reads.
    """

    if "." not in axis:
        raise StudyConfigError(
            f"{label} {axis!r} must use a namespace prefix: 'generator.' or 'evaluation.'."
        )
    namespace, _, field_name = axis.partition(".")
    if namespace not in _AXIS_NAMESPACES:
        raise StudyConfigError(f"{label} {axis!r} has unsupported namespace {namespace!r}.")
    if not field_name:
        raise StudyConfigError(f"{label} {axis!r} is missing a field name.")
    try:
        _split_axis(axis)
    except SimulationGridError as exc:
        raise StudyConfigError(f"{label} {axis!r}: {exc}") from exc


def _split_axis_or_raise(axis: str) -> tuple[str, tuple[str, ...]]:
    try:
        return _split_axis(axis)
    except SimulationGridError as exc:
        raise StudyConfigError(str(exc)) from exc


def _config_to_dict(config: StudyConfig) -> dict[str, Any]:
    payload = _config_to_dict_base(config)
    # Written only when declared, so a dumped config without the rule is
    # unchanged.
    if config.acceptance.anchor_decision is not None:
        payload["acceptance"]["anchor_decision"] = _dataclass_dict(config.acceptance.anchor_decision)
    if config.acceptance.column_decision is not None:
        payload["acceptance"]["column_decision"] = _dataclass_dict(config.acceptance.column_decision)
    return payload


def _config_to_dict_base(config: StudyConfig) -> dict[str, Any]:
    return {
        # An unset size table is omitted, as in the parameter signature, so a
        # dumped config that never used it is unchanged.
        "generator": {
            key: value
            for key, value in _dataclass_dict(config.generator).items()
            if not (key == "group_stage_sizes" and value is None)
        },
        # PC1 is written ``"pc1"``: a config ``null`` names the evaluator's default.
        "evaluation": {
            key: _display_axis_value(f"evaluation.{key}", value)
            for key, value in _dataclass_dict(config.evaluation).items()
        },
        "trajectory_modes": list(config.trajectory_modes),
        "effect_sizes": list(config.effect_sizes),
        "axes": {axis: [_display_axis_value(axis, v) for v in values] for axis, values in config.axes.items()},
        # The explicit form round-trips as the listed points; ``axes`` is derived.
        "design_grid": (
            {"points": [_display_point(point) for point in config.design_grid.points]}
            if config.design_grid.explicit
            else {
                "axes": {
                    axis: [_display_axis_value(axis, v) for v in values]
                    for axis, values in config.design_grid.axes.items()
                }
            }
        ),
        "n_replicates": config.n_replicates,
        "base_seed": config.base_seed,
        "alpha": config.alpha,
        "acceptance": {
            "type_i": [_dataclass_dict(t) for t in config.acceptance.type_i],
            "power": [_dataclass_dict(t) for t in config.acceptance.power],
            "specificity": [_dataclass_dict(t) for t in config.acceptance.specificity],
            "gate": _dataclass_dict(config.acceptance.gate),
            "design_point": (
                None
                if config.acceptance.design_point is None
                else _dataclass_dict(config.acceptance.design_point)
            ),
            "rank_decision": (
                None
                if config.acceptance.rank_decision is None
                else _dataclass_dict(config.acceptance.rank_decision)
            ),
        },
        "attribution": _dataclass_dict(config.attribution),
        "matched_seeds": _dataclass_dict(config.matched_seeds),
        "report_contract": (
            None if config.report_contract is None else _dataclass_dict(config.report_contract)
        ),
        "metadata": dict(config.metadata),
    }


def _dataclass_dict(value: Any) -> dict[str, Any]:
    return {f.name: _to_jsonable(getattr(value, f.name)) for f in dataclasses.fields(value)}


def _to_jsonable(value: Any) -> Any:
    if dataclasses.is_dataclass(value) and not isinstance(value, type):
        return _dataclass_dict(value)
    if isinstance(value, Mapping):
        return {str(k): _to_jsonable(v) for k, v in value.items()}
    if isinstance(value, tuple | list):
        return [_to_jsonable(v) for v in value]
    if isinstance(value, Path):
        return str(value)
    return value


__all__ = [
    "ANCHOR_DECISION_RULE_VERSION",
    "COLUMN_DECISION_CRITERIA",
    "COLUMN_DECISION_RULE_VERSION",
    "AcceptanceTargets",
    "AnchorDecisionRule",
    "AttributionSelector",
    "ColumnDecisionRule",
    "DesignGrid",
    "DesignPointDecisionRule",
    "GateRule",
    "MatchedSeedPolicy",
    "Phase4GateConfig",
    "PowerMonotonicityTarget",
    "RANK_AXIS_BASELINE",
    "RankDecisionRule",
    "ReportContract",
    "SpecificityTarget",
    "StatisticPair",
    "StudyConfig",
    "StudyConfigError",
    "TypeIBound",
    "TypeIControlTarget",
    "dump_study_config",
    "load_study_config",
]
