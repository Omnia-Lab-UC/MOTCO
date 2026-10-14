"""Evaluate study summaries against pre-specified acceptance targets."""

from __future__ import annotations

import json
import math
from collections.abc import Mapping, Sequence
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Any

import pandas as pd

from motco.simulations.grid import (
    SimulationReplicateResult,
    SimulationSummaryResult,
    _record_trajectory_mode,
)
from motco.simulations.study.config import (
    ANCHOR_DECISION_RULE_VERSION,
    AcceptanceTargets,
    AnchorDecisionRule,
    DesignPointDecisionRule,
    PowerMonotonicityTarget,
    RankDecisionRule,
    SpecificityTarget,
    TypeIControlTarget,
)
from motco.simulations.study.phase4 import type_i_inflation_bound
from motco.simulations.study.spectrum import record_design_point


@dataclass(frozen=True)
class TargetEvaluation:
    """Result of evaluating one acceptance target against summaries."""

    target_name: str
    target_kind: str
    met: bool | None
    rationale: str
    observations: dict[str, Any]


def evaluate_targets(
    acceptance: AcceptanceTargets,
    summaries: Sequence[SimulationSummaryResult],
    records: Sequence[SimulationReplicateResult],
) -> list[TargetEvaluation]:
    """Evaluate every acceptance target against the supplied summaries."""

    cell_meta = _cell_metadata_index(records)
    evaluations: list[TargetEvaluation] = []
    for type_i_target in acceptance.type_i:
        evaluations.extend(_evaluate_type_i(type_i_target, summaries, cell_meta))
    for power_target in acceptance.power:
        evaluations.append(_evaluate_power(power_target, summaries, cell_meta))
    for specificity_target in acceptance.specificity:
        evaluations.append(_evaluate_specificity(specificity_target, summaries, cell_meta))
    return evaluations


def write_target_report(
    evaluations: Sequence[TargetEvaluation],
    out_dir: Path,
) -> dict[str, Path]:
    """Write acceptance target evaluations as CSV + JSON."""

    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    csv_path = out_dir / "acceptance_report.csv"
    json_path = out_dir / "acceptance_report.json"
    rows = []
    for ev in evaluations:
        row = {
            "target_name": ev.target_name,
            "target_kind": ev.target_kind,
            "met": ev.met,
            "rationale": ev.rationale,
        }
        for key, value in ev.observations.items():
            row[f"obs_{key}"] = value
        rows.append(row)
    pd.DataFrame(rows).to_csv(csv_path, index=False)
    json_path.write_text(
        json.dumps([asdict(ev) for ev in evaluations], indent=2, default=_json_default) + "\n",
        encoding="utf-8",
    )
    return {"csv": csv_path, "json": json_path}


def _evaluate_type_i(
    target: TypeIControlTarget,
    summaries: Sequence[SimulationSummaryResult],
    cell_meta: dict[str, dict],
) -> list[TargetEvaluation]:
    """Type I target: per null cell × statistic, |rate - alpha| ≤ k · SE.

    Under ``bound="upper"`` the rule is the gate's one-sided control bound
    instead (see :func:`_judge_upper`).
    """

    out: list[TargetEvaluation] = []
    for summary in summaries:
        if not summary.phase.startswith("type_i_"):
            continue
        meta = cell_meta.get(summary.cell_id, {})
        observations = {
            "cell_id": summary.cell_id,
            "phase": summary.phase,
            "trajectory_mode": meta.get("trajectory_mode"),
            "statistic": summary.statistic,
            "alpha_target": target.alpha,
            "se_tolerance": target.se_tolerance,
            "rejection_rate": summary.rejection_rate,
            "monte_carlo_se": summary.monte_carlo_se,
            "available_replicates": summary.available_replicates,
        }
        if summary.rejection_rate is None or summary.monte_carlo_se is None:
            out.append(
                TargetEvaluation(
                    target_name=f"{target.name}[{summary.cell_id},{summary.statistic}]",
                    target_kind=target.kind,
                    met=None,
                    rationale="No available replicates",
                    observations=observations,
                )
            )
            continue
        if target.bound == "upper":
            met, rationale = _judge_upper(
                target, summary.rejection_rate, summary.available_replicates, observations
            )
        else:
            bound = target.se_tolerance * summary.monte_carlo_se
            deviation = abs(summary.rejection_rate - target.alpha)
            met = deviation <= bound
            rationale = (
                f"|rate - alpha| = {deviation:.4f} {'<=' if met else '>'} "
                f"{target.se_tolerance:.2f}·SE ({bound:.4f})"
            )
            observations["deviation"] = deviation
            observations["bound"] = bound
        out.append(
            TargetEvaluation(
                target_name=f"{target.name}[{summary.cell_id},{summary.statistic}]",
                target_kind=target.kind,
                met=met,
                rationale=rationale,
                observations=observations,
            )
        )
    return out


def _judge_upper(
    target: TypeIControlTarget | SpecificityTarget,
    rate: float,
    n_available: int,
    observations: dict[str, Any],
) -> tuple[bool, str]:
    """One-sided control rule shared with the Phase 4 gate.

    Met when ``rate <= alpha + k · sqrt(alpha(1 − alpha)/R)``, R the available
    replicates. Records the bound kind, the alpha-based SE and the *signed*
    deviation (negative is conservative) on ``observations``.
    """

    alpha_se = math.sqrt(target.alpha * (1.0 - target.alpha) / n_available)
    bound = type_i_inflation_bound(target.alpha, target.se_tolerance, n_available)
    met = rate <= bound
    observations["bound_kind"] = "upper"
    observations["alpha_se"] = alpha_se
    observations["deviation"] = rate - target.alpha
    observations["bound"] = bound
    rationale = (
        f"rate = {rate:.4f} {'<=' if met else '>'} "
        f"alpha + {target.se_tolerance:.2f}·SE_alpha ({bound:.4f})"
    )
    return met, rationale


def _evaluate_power(
    target: PowerMonotonicityTarget,
    summaries: Sequence[SimulationSummaryResult],
    cell_meta: dict[str, dict],
) -> TargetEvaluation:
    """Power target: rate non-decreasing in effect_size and ≥ floor at top effect."""

    points: list[tuple[float, float, float | None]] = []
    for summary in summaries:
        meta = cell_meta.get(summary.cell_id, {})
        if meta.get("varied_axis") is not None:
            continue
        if meta.get("trajectory_mode") != target.trajectory_mode:
            continue
        if summary.statistic != target.statistic:
            continue
        effect_size = meta.get("effect_size")
        if effect_size is None or summary.rejection_rate is None:
            continue
        points.append((float(effect_size), float(summary.rejection_rate), summary.monte_carlo_se))
    points.sort(key=lambda x: x[0])
    observations: dict[str, Any] = {
        "trajectory_mode": target.trajectory_mode,
        "statistic": target.statistic,
        "min_power_at_top": target.min_power_at_top,
        "points": [{"effect_size": e, "rate": r, "se": s} for e, r, s in points],
    }
    if not points:
        return TargetEvaluation(
            target_name=target.name,
            target_kind=target.kind,
            met=None,
            rationale="No matching power summaries",
            observations=observations,
        )
    rates = [r for _, r, _ in points]
    monotone = all(rates[i] <= rates[i + 1] for i in range(len(rates) - 1))
    top_rate = rates[-1]
    meets_floor = top_rate >= target.min_power_at_top
    met = monotone and meets_floor
    rationale = (
        f"monotone={monotone}; top_rate={top_rate:.3f} "
        f"{'>=' if meets_floor else '<'} floor={target.min_power_at_top:.3f}"
    )
    observations["monotone"] = monotone
    observations["top_rate"] = top_rate
    return TargetEvaluation(
        target_name=target.name,
        target_kind=target.kind,
        met=met,
        rationale=rationale,
        observations=observations,
    )


def _evaluate_specificity(
    target: SpecificityTarget,
    summaries: Sequence[SimulationSummaryResult],
    cell_meta: dict[str, dict],
) -> TargetEvaluation:
    """Specificity target: off-diagonal rate ≈ alpha within k · SE.

    Under ``bound="upper"`` the rule is the gate's one-sided control bound
    instead (see :func:`_judge_upper`). Uses the largest effect-size cell for the given (mode, statistic) pair.
    """

    candidates: list[tuple[float, SimulationSummaryResult]] = []
    for summary in summaries:
        meta = cell_meta.get(summary.cell_id, {})
        if meta.get("varied_axis") is not None:
            continue
        if meta.get("trajectory_mode") != target.trajectory_mode:
            continue
        if summary.statistic != target.statistic:
            continue
        effect_size = meta.get("effect_size")
        if effect_size is None:
            continue
        candidates.append((float(effect_size), summary))
    candidates.sort(key=lambda x: x[0])
    observations: dict[str, Any] = {
        "trajectory_mode": target.trajectory_mode,
        "statistic": target.statistic,
        "alpha_target": target.alpha,
        "se_tolerance": target.se_tolerance,
    }
    if not candidates:
        return TargetEvaluation(
            target_name=target.name,
            target_kind=target.kind,
            met=None,
            rationale="No matching specificity summaries",
            observations=observations,
        )
    effect_size, summary = candidates[-1]
    observations.update(
        {
            "effect_size": effect_size,
            "cell_id": summary.cell_id,
            "rejection_rate": summary.rejection_rate,
            "monte_carlo_se": summary.monte_carlo_se,
            "available_replicates": summary.available_replicates,
        }
    )
    if summary.rejection_rate is None or summary.monte_carlo_se is None:
        return TargetEvaluation(
            target_name=target.name,
            target_kind=target.kind,
            met=None,
            rationale="No available replicates at top effect size",
            observations=observations,
        )
    if target.bound == "upper":
        met, rationale = _judge_upper(
            target, summary.rejection_rate, summary.available_replicates, observations
        )
        return TargetEvaluation(
            target_name=target.name,
            target_kind=target.kind,
            met=met,
            rationale=rationale,
            observations=observations,
        )
    bound = target.se_tolerance * summary.monte_carlo_se
    deviation = abs(summary.rejection_rate - target.alpha)
    met = deviation <= bound
    rationale = (
        f"|rate - alpha| = {deviation:.4f} {'<=' if met else '>'} "
        f"{target.se_tolerance:.2f}·SE ({bound:.4f})"
    )
    observations["deviation"] = deviation
    observations["bound"] = bound
    return TargetEvaluation(
        target_name=target.name,
        target_kind=target.kind,
        met=met,
        rationale=rationale,
        observations=observations,
    )


# ── Design-point decision ─────────────────────────────────────────────────────

#: Per-column classifications, in the order they are reported.
DESIGN_POINT_STATUSES: tuple[str, ...] = ("meets", "marginal", "fails", "unavailable")


@dataclass(frozen=True)
class DesignPointColumn:
    """One design point's standing against the predeclared rule."""

    design_point: dict[str, Any]
    is_baseline: bool
    top_effect: float | None
    n_available: int
    rejection_rate: float | None
    monte_carlo_se: float | None
    lower_bound: float | None
    status: str
    anchor_rates: dict[str, float | None]
    anchor_available: int


@dataclass(frozen=True)
class DesignPointDecision:
    """Outcome of the design-point rule over every enumerated design point.

    ``verdict`` is ``"chosen"`` when a confirmed column exists (``chosen`` names
    it), ``"revise_claim"`` when none is confirmed, and ``"no_design_grid"``
    when the records carry no design points at all. Advisory only: it feeds
    neither the Phase 4 gate nor the acceptance-target report.
    """

    rule: dict[str, Any]
    alpha: float
    verdict: str
    chosen: dict[str, Any] | None
    rationale: str
    columns: tuple[DesignPointColumn, ...]

    def to_frame(self) -> pd.DataFrame:
        rows = []
        for column in self.columns:
            row: dict[str, Any] = {**column.design_point}
            row.update(
                {
                    "is_baseline": column.is_baseline,
                    "top_effect": column.top_effect,
                    "n_available": column.n_available,
                    "rejection_rate": column.rejection_rate,
                    "monte_carlo_se": column.monte_carlo_se,
                    "lower_bound": column.lower_bound,
                    "status": column.status,
                    "chosen": self.chosen is not None and column.design_point == self.chosen,
                    "anchor_available": column.anchor_available,
                }
            )
            for statistic, rate in column.anchor_rates.items():
                row[f"anchor_{statistic}_rate"] = rate
            rows.append(row)
        return pd.DataFrame(rows)


def evaluate_design_point_decision(
    records: Sequence[SimulationReplicateResult],
    rule: DesignPointDecisionRule,
    *,
    alpha: float = 0.05,
    statistics: Sequence[str] = ("delta", "angle", "shape"),
) -> DesignPointDecision:
    """Classify every design point against ``rule`` and pick the first confirmed one.

    For each design point the target statistic's rejection rate is taken at the
    **largest effect enumerated in that column** for the target mode. A column
    *meets* the floor when ``rate − k·SE ≥ floor`` (``k`` the rule's
    confirmation threshold), is *marginal* when ``rate ≥ floor`` without
    confirmation, *fails* below the floor, and is *unavailable* when no target
    record exists. Columns are ordered by the rule's ``prefer`` axes (ascending,
    in the declared order) and then by the remaining design axes; the chosen
    point is the first *meets* column in that order. Every threshold comes from
    ``rule``; nothing is hard-coded here.
    """

    if not (0 < alpha < 1):
        raise ValueError("alpha must be between 0 and 1.")

    by_point: dict[tuple[tuple[str, Any], ...], list[SimulationReplicateResult]] = {}
    axes: list[str] = []
    for record in records:
        if record.status != "completed":
            continue
        point = record_design_point(record)
        if point is None:
            continue
        for axis in point:
            if axis not in axes:
                axes.append(axis)
        by_point.setdefault(tuple(point.items()), []).append(record)

    rule_payload = {
        "trajectory_mode": rule.trajectory_mode,
        "statistic": rule.statistic,
        "min_power_at_top": rule.min_power_at_top,
        "confirmation_se_threshold": rule.confirmation_se_threshold,
        "prefer": list(rule.prefer),
    }
    if not by_point or not any(r.phase == "power_design" for group in by_point.values() for r in group):
        return DesignPointDecision(
            rule=rule_payload,
            alpha=alpha,
            verdict="no_design_grid",
            chosen=None,
            rationale="No design-point records are present; the rule has nothing to evaluate.",
            columns=(),
        )

    order_axes = [axis for axis in rule.prefer if axis in axes] + [
        axis for axis in axes if axis not in rule.prefer
    ]
    columns: list[DesignPointColumn] = []
    for coords, group in by_point.items():
        point = dict(coords)
        target = [
            r
            for r in group
            if str(r.cell_metadata.get("trajectory_mode")) == rule.trajectory_mode
            and r.cell_metadata.get("effect_size") is not None
            and float(r.cell_metadata["effect_size"]) != 0.0
        ]
        anchors = [
            r
            for r in group
            if r.cell_metadata.get("effect_size") is not None
            and float(r.cell_metadata["effect_size"]) == 0.0
        ]
        anchor_rates: dict[str, float | None] = {}
        anchor_available = 0
        for statistic in statistics:
            p_values = _finite_p_values(anchors, statistic)
            anchor_rates[statistic] = (
                sum(p < alpha for p in p_values) / len(p_values) if p_values else None
            )
            anchor_available = max(anchor_available, len(p_values))

        if not target:
            columns.append(
                DesignPointColumn(
                    design_point=point,
                    is_baseline=all(r.phase == "power_primary" for r in group),
                    top_effect=None,
                    n_available=0,
                    rejection_rate=None,
                    monte_carlo_se=None,
                    lower_bound=None,
                    status="unavailable",
                    anchor_rates=anchor_rates,
                    anchor_available=anchor_available,
                )
            )
            continue

        top_effect = max(float(r.cell_metadata["effect_size"]) for r in target)
        at_top = [r for r in target if float(r.cell_metadata["effect_size"]) == top_effect]
        p_values = _finite_p_values(at_top, rule.statistic)
        if not p_values:
            rate = se = lower = None
            status = "unavailable"
        else:
            rate = sum(p < alpha for p in p_values) / len(p_values)
            se = math.sqrt(rate * (1.0 - rate) / len(p_values))
            lower = rate - rule.confirmation_se_threshold * se
            if lower >= rule.min_power_at_top:
                status = "meets"
            elif rate >= rule.min_power_at_top:
                status = "marginal"
            else:
                status = "fails"
        columns.append(
            DesignPointColumn(
                design_point=point,
                is_baseline=all(r.phase == "power_primary" for r in group),
                top_effect=top_effect,
                n_available=len(p_values),
                rejection_rate=rate,
                monte_carlo_se=se,
                lower_bound=lower,
                status=status,
                anchor_rates=anchor_rates,
                anchor_available=anchor_available,
            )
        )

    columns.sort(key=lambda column: tuple(_orderable(column.design_point.get(axis)) for axis in order_axes))
    chosen = next((column for column in columns if column.status == "meets"), None)
    if chosen is not None:
        verdict = "chosen"
        rationale = (
            f"First design point in preference order {order_axes} whose {rule.trajectory_mode}/"
            f"{rule.statistic} power at its top effect ({chosen.top_effect:g}) clears the "
            f"{rule.min_power_at_top:.2f} floor with {rule.confirmation_se_threshold:g}·SE "
            f"confirmation: rate {chosen.rejection_rate:.3f}, lower bound {chosen.lower_bound:.3f}."
        )
    else:
        verdict = "revise_claim"
        counts = {status: sum(c.status == status for c in columns) for status in DESIGN_POINT_STATUSES}
        rationale = (
            f"No design point clears the {rule.min_power_at_top:.2f} floor for {rule.trajectory_mode}/"
            f"{rule.statistic} with {rule.confirmation_se_threshold:g}·SE confirmation "
            f"({counts['marginal']} marginal, {counts['fails']} failing, {counts['unavailable']} "
            "unavailable). Per the readiness worklist, revise the method or the claim — not the "
            "Monte Carlo size."
        )
    return DesignPointDecision(
        rule=rule_payload,
        alpha=alpha,
        verdict=verdict,
        chosen=None if chosen is None else dict(chosen.design_point),
        rationale=rationale,
        columns=tuple(columns),
    )


def write_design_point_decision(decision: DesignPointDecision, out_dir: Path) -> dict[str, Path]:
    """Write the design-point decision as JSON (verdict + columns) and CSV (columns)."""

    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    json_path = out_dir / "design_point_decision.json"
    csv_path = out_dir / "design_point_decision.csv"
    payload = {
        "verdict": decision.verdict,
        "chosen": decision.chosen,
        "rationale": decision.rationale,
        "alpha": decision.alpha,
        "rule": decision.rule,
        "columns": [asdict(column) for column in decision.columns],
    }
    json_path.write_text(json.dumps(payload, indent=2, default=_json_default) + "\n", encoding="utf-8")
    decision.to_frame().to_csv(csv_path, index=False)
    return {"design_point_decision": json_path, "design_point_decision_csv": csv_path}


# ── Retained-rank decision ────────────────────────────────────────────────────

#: Names of the three qualifying criteria, in the order they are reported.
RANK_CRITERIA: tuple[str, ...] = ("gain", "anchor", "protected")


@dataclass(frozen=True)
class RateEstimate:
    """A rejection rate with its Monte Carlo standard error and denominator."""

    rate: float | None
    monte_carlo_se: float | None
    n_available: int

    @property
    def available(self) -> bool:
        return self.rate is not None and self.monte_carlo_se is not None


@dataclass(frozen=True)
class RankCriterion:
    """One criterion's standing for one forced-rank column.

    ``passed`` is ``None`` when a quantity it needs is unavailable; a column
    with any unavailable criterion does not qualify. ``failing`` names the
    pairs (protected) or statistics (anchor) that broke the criterion, and
    ``observations`` carries every rate, SE, and threshold the verdict used.
    """

    name: str
    passed: bool | None
    detail: str
    failing: tuple[str, ...]
    observations: dict[str, Any]


@dataclass(frozen=True)
class RankColumn:
    """One retained-rank column's standing against the predeclared rule."""

    rank: int | None
    is_reference: bool
    top_effect: float | None
    target: RateEstimate
    anchor: dict[str, RateEstimate]
    protected: dict[str, RateEstimate]
    criteria: tuple[RankCriterion, ...]
    qualifies: bool | None

    def criterion(self, name: str) -> RankCriterion | None:
        return next((criterion for criterion in self.criteria if criterion.name == name), None)


@dataclass(frozen=True)
class RankDecision:
    """Outcome of the retained-rank rule over every column of the rank axis.

    ``verdict`` is ``"adopt_fixed_rank"`` (``chosen_rank`` names the smallest
    qualifying rank), ``"keep_cv"`` when no column qualifies,
    ``"no_reference"`` when the records carry no cross-validated (``null``)
    column, and ``"no_design_grid"`` when no design-point record is present.
    Advisory only: it feeds neither the Phase 4 gate nor the acceptance targets.
    """

    rule: dict[str, Any]
    alpha: float
    verdict: str
    chosen_rank: int | None
    rationale: str
    reference: RankColumn | None
    columns: tuple[RankColumn, ...]

    def to_frame(self) -> pd.DataFrame:
        rows = []
        for column in self.columns:
            row: dict[str, Any] = {
                "rank": column.rank,
                "component_selection": "cv" if column.is_reference else "forced",
                "is_reference": column.is_reference,
                "top_effect": column.top_effect,
                "target_rate": column.target.rate,
                "target_se": column.target.monte_carlo_se,
                "target_available": column.target.n_available,
            }
            for statistic, estimate in column.anchor.items():
                row[f"anchor_{statistic}_rate"] = estimate.rate
                row[f"anchor_{statistic}_available"] = estimate.n_available
            for label, estimate in column.protected.items():
                key = label.replace("/", "_")
                row[f"protected_{key}_rate"] = estimate.rate
                row[f"protected_{key}_se"] = estimate.monte_carlo_se
            for criterion in column.criteria:
                row[f"{criterion.name}_passed"] = criterion.passed
                row[f"{criterion.name}_failing"] = ", ".join(criterion.failing) or None
                for key, value in criterion.observations.items():
                    if isinstance(value, bool | int | float) or value is None:
                        row[f"{criterion.name}_{key}"] = value
            row["qualifies"] = column.qualifies
            row["chosen"] = self.chosen_rank is not None and column.rank == self.chosen_rank
            rows.append(row)
        return pd.DataFrame(rows)


def evaluate_rank_decision(
    records: Sequence[SimulationReplicateResult],
    rule: RankDecisionRule,
    *,
    alpha: float = 0.05,
    statistics: Sequence[str] = ("delta", "angle", "shape"),
) -> RankDecision:
    """Evaluate the predeclared retained-rank rule against the cross-validated column.

    Columns are the values of ``rule.axis`` across the design points; the
    reference is the column whose value is ``null``. For every forced column
    ``k`` and the reference ``cv``, at the largest effect enumerated for each
    mode:

    - **gain**: ``rate_k − rate_cv > gain_se_multiplier · sqrt(SE_k² + SE_cv²)``
      on the target pair;
    - **anchor**: for every statistic, the column's zero-effect rate is at most
      ``alpha_b + se_tolerance · sqrt(alpha_b(1 − alpha_b)/n)`` (``alpha_b`` the
      rule's Type I bound), the same bound the Type I control target uses;
    - **protected**: for every protected pair,
      ``rate_cv − rate_k ≤ loss_se_multiplier · sqrt(SE_k² + SE_cv²)``.

    Columns are paired on identical datasets (same matched-seed family, same
    generator parameters), so the pooled independent SE is conservative — the
    right direction for a rule that would change the production configuration.
    Every threshold is read from ``rule``; nothing is hard-coded here.
    """

    if not (0 < alpha < 1):
        raise ValueError("alpha must be between 0 and 1.")

    rule_payload = {
        "axis": rule.axis,
        "target": {"trajectory_mode": rule.target.trajectory_mode, "statistic": rule.target.statistic},
        "protected": [
            {"trajectory_mode": pair.trajectory_mode, "statistic": pair.statistic} for pair in rule.protected
        ],
        "type_i_bound": {"alpha": rule.type_i_bound.alpha, "se_tolerance": rule.type_i_bound.se_tolerance},
        "gain_se_multiplier": rule.gain_se_multiplier,
        "loss_se_multiplier": rule.loss_se_multiplier,
    }

    by_rank: dict[int | None, list[SimulationReplicateResult]] = {}
    saw_design = False
    for record in records:
        if record.status != "completed":
            continue
        point = record_design_point(record)
        if point is None or rule.axis not in point:
            continue
        saw_design = saw_design or record.phase == "power_design"
        value = point[rule.axis]
        by_rank.setdefault(None if value is None else int(value), []).append(record)

    if not saw_design:
        return RankDecision(
            rule=rule_payload,
            alpha=alpha,
            verdict="no_design_grid",
            chosen_rank=None,
            rationale="No design-point records on the rank axis are present; the rule has nothing to evaluate.",
            reference=None,
            columns=(),
        )

    measured = {
        rank: _measure_rank_column(group, rule, alpha=alpha, statistics=statistics)
        for rank, group in by_rank.items()
    }
    reference_measure = measured.get(None)
    forced_ranks = sorted(k for k in measured if k is not None)
    if reference_measure is None:
        return RankDecision(
            rule=rule_payload,
            alpha=alpha,
            verdict="no_reference",
            chosen_rank=None,
            rationale=(
                f"No cross-validated reference column (axis {rule.axis!r} = null) is present, so no "
                "forced rank can be compared."
            ),
            reference=None,
            columns=tuple(
                _rank_column(rank, measured[rank], None, rule, statistics=statistics) for rank in forced_ranks
            ),
        )

    reference = _rank_column(None, reference_measure, None, rule, statistics=statistics)
    columns: list[RankColumn] = [reference]
    for rank in forced_ranks:
        columns.append(_rank_column(rank, measured[rank], reference_measure, rule, statistics=statistics))
    qualifying = [column.rank for column in columns if column.qualifies and column.rank is not None]
    if qualifying:
        chosen: int | None = min(qualifying)
        verdict = "adopt_fixed_rank"
        rationale = (
            f"Rank(s) {qualifying} clear all three criteria against the cross-validated column on "
            f"{rule.target.label} at the top effect; the smallest, {chosen}, is adopted (parsimony "
            "predeclared). Advisory: the readiness worklist records the committed rule."
        )
    else:
        chosen = None
        verdict = "keep_cv"
        failing = {
            column.rank: [criterion.name for criterion in column.criteria if criterion.passed is not True]
            for column in columns
            if column.rank is not None
        }
        rationale = (
            f"No forced rank clears every criterion against the cross-validated column on "
            f"{rule.target.label} (failing criteria per rank: {failing}); stage-supervised double "
            "cross-validation remains the committed rank rule."
        )
    return RankDecision(
        rule=rule_payload,
        alpha=alpha,
        verdict=verdict,
        chosen_rank=chosen,
        rationale=rationale,
        reference=reference,
        columns=tuple(columns),
    )


@dataclass(frozen=True)
class _RankMeasure:
    top_effect: float | None
    target: RateEstimate
    anchor: dict[str, RateEstimate]
    protected: dict[str, RateEstimate]


def _measure_rank_column(
    group: Sequence[SimulationReplicateResult],
    rule: RankDecisionRule,
    *,
    alpha: float,
    statistics: Sequence[str],
) -> _RankMeasure:
    anchors = [
        r
        for r in group
        if r.cell_metadata.get("effect_size") is not None and float(r.cell_metadata["effect_size"]) == 0.0
    ]
    anchor = {statistic: _rate_estimate(anchors, statistic, alpha) for statistic in statistics}
    target_records, top_effect = _top_effect_records(group, rule.target.trajectory_mode)
    target = _rate_estimate(target_records, rule.target.statistic, alpha)
    protected: dict[str, RateEstimate] = {}
    for pair in rule.protected:
        pair_records, _ = _top_effect_records(group, pair.trajectory_mode)
        protected[pair.label] = _rate_estimate(pair_records, pair.statistic, alpha)
    return _RankMeasure(top_effect=top_effect, target=target, anchor=anchor, protected=protected)


def _top_effect_records(
    group: Sequence[SimulationReplicateResult], mode: str
) -> tuple[list[SimulationReplicateResult], float | None]:
    candidates = [
        r
        for r in group
        if str(r.cell_metadata.get("trajectory_mode")) == mode
        and r.cell_metadata.get("effect_size") is not None
        and float(r.cell_metadata["effect_size"]) != 0.0
    ]
    if not candidates:
        return [], None
    top = max(float(r.cell_metadata["effect_size"]) for r in candidates)
    return [r for r in candidates if float(r.cell_metadata["effect_size"]) == top], top


def _rate_estimate(records: Sequence[SimulationReplicateResult], statistic: str, alpha: float) -> RateEstimate:
    p_values = _finite_p_values(records, statistic)
    if not p_values:
        return RateEstimate(rate=None, monte_carlo_se=None, n_available=0)
    rate = sum(p < alpha for p in p_values) / len(p_values)
    return RateEstimate(
        rate=rate,
        monte_carlo_se=math.sqrt(rate * (1.0 - rate) / len(p_values)),
        n_available=len(p_values),
    )


def _pooled_se(left: RateEstimate, right: RateEstimate) -> float:
    return math.sqrt(float(left.monte_carlo_se or 0.0) ** 2 + float(right.monte_carlo_se or 0.0) ** 2)


def _rank_column(
    rank: int | None,
    measure: _RankMeasure,
    reference: _RankMeasure | None,
    rule: RankDecisionRule,
    *,
    statistics: Sequence[str],
) -> RankColumn:
    criteria: tuple[RankCriterion, ...] = ()
    qualifies: bool | None = None
    if rank is not None and reference is not None:
        criteria = (
            _gain_criterion(measure, reference, rule),
            _anchor_criterion(measure, rule, statistics=statistics),
            _protected_criterion(measure, reference, rule),
        )
        if any(criterion.passed is None for criterion in criteria):
            qualifies = None
        else:
            qualifies = all(bool(criterion.passed) for criterion in criteria)
    return RankColumn(
        rank=rank,
        is_reference=rank is None,
        top_effect=measure.top_effect,
        target=measure.target,
        anchor=dict(measure.anchor),
        protected=dict(measure.protected),
        criteria=criteria,
        qualifies=qualifies,
    )


def _gain_criterion(measure: _RankMeasure, reference: _RankMeasure, rule: RankDecisionRule) -> RankCriterion:
    observations: dict[str, Any] = {
        "target": rule.target.label,
        "rate": measure.target.rate,
        "se": measure.target.monte_carlo_se,
        "reference_rate": reference.target.rate,
        "reference_se": reference.target.monte_carlo_se,
        "se_multiplier": rule.gain_se_multiplier,
    }
    if measure.target.rate is None or reference.target.rate is None:
        return RankCriterion(
            name="gain",
            passed=None,
            detail=f"{rule.target.label} rate unavailable in this column or the reference.",
            failing=(),
            observations=observations,
        )
    difference = measure.target.rate - reference.target.rate
    pooled = _pooled_se(measure.target, reference.target)
    threshold = rule.gain_se_multiplier * pooled
    passed = difference > threshold
    observations.update({"difference": difference, "pooled_se": pooled, "threshold": threshold})
    return RankCriterion(
        name="gain",
        passed=passed,
        detail=(
            f"{rule.target.label}: rate − reference = {difference:+.3f} "
            f"{'>' if passed else '<='} {rule.gain_se_multiplier:g}·pooled SE ({threshold:.3f})"
        ),
        failing=() if passed else (rule.target.label,),
        observations=observations,
    )


def _anchor_criterion(
    measure: _RankMeasure, rule: RankDecisionRule, *, statistics: Sequence[str]
) -> RankCriterion:
    bound_alpha = rule.type_i_bound.alpha
    tolerance = rule.type_i_bound.se_tolerance
    observations: dict[str, Any] = {"alpha": bound_alpha, "se_tolerance": tolerance}
    failing: list[str] = []
    unavailable: list[str] = []
    for statistic in statistics:
        estimate = measure.anchor.get(statistic)
        if estimate is None or estimate.rate is None:
            unavailable.append(statistic)
            observations[f"{statistic}_rate"] = None
            observations[f"{statistic}_bound"] = None
            continue
        bound = bound_alpha + tolerance * math.sqrt(bound_alpha * (1.0 - bound_alpha) / estimate.n_available)
        observations[f"{statistic}_rate"] = estimate.rate
        observations[f"{statistic}_bound"] = bound
        observations[f"{statistic}_n"] = estimate.n_available
        if estimate.rate > bound:
            failing.append(statistic)
    if unavailable:
        return RankCriterion(
            name="anchor",
            passed=None,
            detail=f"Anchor rate unavailable for {unavailable}.",
            failing=(),
            observations=observations,
        )
    passed = not failing
    return RankCriterion(
        name="anchor",
        passed=passed,
        detail=(
            "Every statistic's zero-effect rate is within the Type I bound."
            if passed
            else f"Anchor inflated on {failing}: rate exceeds alpha + {tolerance:g}·SE."
        ),
        failing=tuple(failing),
        observations=observations,
    )


def _protected_criterion(
    measure: _RankMeasure, reference: _RankMeasure, rule: RankDecisionRule
) -> RankCriterion:
    observations: dict[str, Any] = {"se_multiplier": rule.loss_se_multiplier}
    failing: list[str] = []
    unavailable: list[str] = []
    for pair in rule.protected:
        estimate = measure.protected.get(pair.label)
        baseline = reference.protected.get(pair.label)
        key = pair.label.replace("/", "_")
        if estimate is None or baseline is None or estimate.rate is None or baseline.rate is None:
            unavailable.append(pair.label)
            observations[f"{key}_loss"] = None
            observations[f"{key}_threshold"] = None
            continue
        loss = baseline.rate - estimate.rate
        threshold = rule.loss_se_multiplier * _pooled_se(estimate, baseline)
        observations[f"{key}_rate"] = estimate.rate
        observations[f"{key}_reference_rate"] = baseline.rate
        observations[f"{key}_loss"] = loss
        observations[f"{key}_threshold"] = threshold
        if loss > threshold:
            failing.append(pair.label)
    if unavailable:
        return RankCriterion(
            name="protected",
            passed=None,
            detail=f"Protected rate unavailable for {unavailable}.",
            failing=(),
            observations=observations,
        )
    passed = not failing
    return RankCriterion(
        name="protected",
        passed=passed,
        detail=(
            "No protected power falls below the reference by more than the allowed loss."
            if passed
            else f"Protected loss on {failing}: reference − rate exceeds {rule.loss_se_multiplier:g}·pooled SE."
        ),
        failing=tuple(failing),
        observations=observations,
    )


def write_rank_decision(decision: RankDecision, out_dir: Path) -> dict[str, Path]:
    """Write the retained-rank decision as JSON (verdict + columns) and CSV (columns)."""

    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    json_path = out_dir / "rank_decision.json"
    csv_path = out_dir / "rank_decision.csv"
    payload = {
        "verdict": decision.verdict,
        "chosen_rank": decision.chosen_rank,
        "rationale": decision.rationale,
        "alpha": decision.alpha,
        "rule": decision.rule,
        "pairing": (
            "Columns share the primary matched-seed family and generator parameters: at every "
            "replicate index they evaluate the same generated dataset under different retained "
            "ranks, so differences are measurement differences, not sampling differences."
        ),
        "reference": None if decision.reference is None else asdict(decision.reference),
        "columns": [asdict(column) for column in decision.columns],
    }
    json_path.write_text(json.dumps(payload, indent=2, default=_json_default) + "\n", encoding="utf-8")
    decision.to_frame().to_csv(csv_path, index=False)
    return {"rank_decision": json_path, "rank_decision_csv": csv_path}


#: Classification of one (design point, statistic) anchor rate.
ANCHOR_CLASSIFICATIONS: tuple[str, ...] = ("holds", "exceeds", "unavailable")


@dataclass(frozen=True)
class AnchorStatistic:
    """One statistic's anchor rejection rate against the Type I bound."""

    statistic: str
    rate: float | None
    monte_carlo_se: float | None
    n_available: int
    bound: float | None
    classification: str


@dataclass(frozen=True)
class ControlRate:
    """A non-anchor negative control at a design point, reported descriptively."""

    phase: str
    trajectory_mode: str
    effect_size: float | None
    rates: dict[str, RateEstimate]


@dataclass(frozen=True)
class AnchorPoint:
    """One design point's anchor standing and its other negative controls.

    ``verdict`` is ``size_holds`` when every named statistic holds, ``excess``
    when any exceeds, and ``unavailable`` when none exceeds but some statistic
    has no available replicate (a missing anchor is never read as either).
    """

    design_point: dict[str, Any]
    is_baseline: bool
    statistics: tuple[AnchorStatistic, ...]
    verdict: str
    other_controls: tuple[ControlRate, ...]


@dataclass(frozen=True)
class AnchorDecision:
    """Outcome of the anchor Type I rule over every design point.

    ``verdict`` summarizes the points: ``size_holds`` when every point holds,
    ``excess`` when any point exceeds, ``unavailable`` otherwise, and
    ``no_anchor`` when the records carry no zero-effect anchor. Advisory only.
    """

    rule: dict[str, Any]
    alpha: float
    verdict: str
    rationale: str
    points: tuple[AnchorPoint, ...]

    def to_frame(self) -> pd.DataFrame:
        axes: list[str] = []
        for point in self.points:
            for axis in point.design_point:
                if axis not in axes:
                    axes.append(axis)
        rows: list[dict[str, Any]] = []
        for point in self.points:
            coords = {axis: point.design_point.get(axis) for axis in axes}
            for entry in point.statistics:
                rows.append(
                    {
                        **coords,
                        "is_baseline": point.is_baseline,
                        "role": "anchor",
                        "phase": None,
                        "trajectory_mode": "none",
                        "effect_size": 0.0,
                        "statistic": entry.statistic,
                        "rate": entry.rate,
                        "monte_carlo_se": entry.monte_carlo_se,
                        "n_available": entry.n_available,
                        "bound": entry.bound,
                        "classification": entry.classification,
                        "point_verdict": point.verdict,
                    }
                )
            for control in point.other_controls:
                for statistic, estimate in control.rates.items():
                    rows.append(
                        {
                            **coords,
                            "is_baseline": point.is_baseline,
                            "role": "control",
                            "phase": control.phase,
                            "trajectory_mode": control.trajectory_mode,
                            "effect_size": control.effect_size,
                            "statistic": statistic,
                            "rate": estimate.rate,
                            "monte_carlo_se": estimate.monte_carlo_se,
                            "n_available": estimate.n_available,
                            "bound": None,
                            "classification": None,
                            "point_verdict": point.verdict,
                        }
                    )
        return pd.DataFrame(rows)


def anchor_bound(alpha: float, se_multiplier: float, n_available: int) -> float:
    """``alpha + se_multiplier · sqrt(alpha(1 − alpha)/n_available)``."""

    return alpha + se_multiplier * math.sqrt(alpha * (1.0 - alpha) / n_available)


def evaluate_anchor_decision(
    records: Sequence[SimulationReplicateResult],
    rule: AnchorDecisionRule,
    *,
    alpha: float = 0.05,
) -> AnchorDecision:
    """Evaluate the predeclared anchor Type I rule at every design point.

    The anchor at a point is its zero-effect cell (mode ``none`` at effect 0) in
    the ``power_primary``/``power_design`` phases; a study without a design grid
    has one point, ``{}``, served by the primary anchor. ``alpha`` is the level
    a p-value is rejected at; the bound uses ``rule.alpha`` and the anchor's
    *available* replicate count, so a failed replicate widens the bound rather
    than silently tightening the test. Every other cell at the point — and the
    ``type_i_baseline`` cells, which evaluate in the baseline column — is listed
    under ``other_controls`` with its rates and no classification.
    """

    if not (0 < alpha < 1):
        raise ValueError("alpha must be between 0 and 1.")
    rule_payload = {
        "version": ANCHOR_DECISION_RULE_VERSION,
        "name": rule.name,
        "alpha": rule.alpha,
        "se_multiplier": rule.se_multiplier,
        "statistics": list(rule.statistics),
        "bound": "alpha + se_multiplier * sqrt(alpha * (1 - alpha) / R), R = available anchor replicates",
    }

    def key_of(point: Mapping[str, Any]) -> str:
        return json.dumps(point, sort_keys=True, default=str)

    points: dict[str, dict[str, Any]] = {}
    groups: dict[str, dict[tuple[str, str, float | None], list[SimulationReplicateResult]]] = {}
    baseline_key: str | None = None
    unplaced: list[SimulationReplicateResult] = []
    for record in records:
        if record.status != "completed":
            continue
        point = record_design_point(record)
        if point is None:
            if record.phase == "power_primary":
                point = {}
            else:
                if record.phase == "type_i_baseline":
                    unplaced.append(record)
                continue
        key = key_of(point)
        points.setdefault(key, point)
        if record.phase == "power_primary":
            baseline_key = key
        meta = record.cell_metadata
        effect = meta.get("effect_size")
        group_key = (record.phase, _record_trajectory_mode(record), None if effect is None else float(effect))
        groups.setdefault(key, {}).setdefault(group_key, []).append(record)
    if baseline_key is not None:
        for record in unplaced:
            effect = record.cell_metadata.get("effect_size")
            group_key = (record.phase, _record_trajectory_mode(record), None if effect is None else float(effect))
            groups[baseline_key].setdefault(group_key, []).append(record)

    entries: list[AnchorPoint] = []
    for key in sorted(points, key=lambda k: (k != baseline_key, k)):
        anchors: list[SimulationReplicateResult] = []
        controls: list[ControlRate] = []
        for (phase, mode, effect), group in sorted(groups[key].items(), key=lambda item: _group_order(item[0])):
            if phase in ("power_primary", "power_design") and mode == "none" and effect == 0.0:
                anchors.extend(group)
                continue
            controls.append(
                ControlRate(
                    phase=phase,
                    trajectory_mode=mode,
                    effect_size=effect,
                    rates={statistic: _rate_estimate(group, statistic, alpha) for statistic in rule.statistics},
                )
            )
        if not anchors:
            continue
        statistics: list[AnchorStatistic] = []
        for statistic in rule.statistics:
            estimate = _rate_estimate(anchors, statistic, alpha)
            if estimate.rate is None:
                statistics.append(AnchorStatistic(statistic, None, None, 0, None, "unavailable"))
                continue
            bound = anchor_bound(rule.alpha, rule.se_multiplier, estimate.n_available)
            statistics.append(
                AnchorStatistic(
                    statistic=statistic,
                    rate=estimate.rate,
                    monte_carlo_se=estimate.monte_carlo_se,
                    n_available=estimate.n_available,
                    bound=bound,
                    classification="holds" if estimate.rate <= bound else "exceeds",
                )
            )
        classes = {entry.classification for entry in statistics}
        verdict = "excess" if "exceeds" in classes else ("unavailable" if "unavailable" in classes else "size_holds")
        entries.append(
            AnchorPoint(
                design_point=dict(points[key]),
                is_baseline=key == baseline_key,
                statistics=tuple(statistics),
                verdict=verdict,
                other_controls=tuple(controls),
            )
        )

    if not entries:
        return AnchorDecision(
            rule=rule_payload,
            alpha=alpha,
            verdict="no_anchor",
            rationale="No zero-effect anchor records are present; the rule has nothing to evaluate.",
            points=(),
        )
    verdicts = [entry.verdict for entry in entries]
    if "excess" in verdicts:
        verdict = "excess"
        exceeding = [
            f"{entry.design_point} ({', '.join(_exceeding_statistics(entry))})"
            for entry in entries
            if entry.verdict == "excess"
        ]
        rationale = (
            f"Anchor rate above the bound at {len(exceeding)} of {len(entries)} point(s): "
            f"{'; '.join(exceeding)}."
        )
    elif "unavailable" in verdicts:
        verdict = "unavailable"
        rationale = "No point exceeds the bound, but some anchor statistic has no available replicate."
    else:
        verdict = "size_holds"
        rationale = f"Every named statistic holds the bound at all {len(entries)} point(s)."
    return AnchorDecision(rule=rule_payload, alpha=alpha, verdict=verdict, rationale=rationale, points=tuple(entries))


def _group_order(group: tuple[str, str, float | None]) -> tuple[str, str, float]:
    phase, mode, effect = group
    return (phase, mode, math.inf if effect is None else effect)


def _exceeding_statistics(point: AnchorPoint) -> list[str]:
    return [entry.statistic for entry in point.statistics if entry.classification == "exceeds"]


def write_anchor_decision(decision: AnchorDecision, out_dir: Path) -> dict[str, Path]:
    """Write the anchor decision as JSON (verdicts + every point) and CSV (one row per rate)."""

    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    json_path = out_dir / "anchor_decision.json"
    csv_path = out_dir / "anchor_decision.csv"
    payload = {
        "verdict": decision.verdict,
        "rationale": decision.rationale,
        "alpha": decision.alpha,
        "rule": decision.rule,
        "points": [asdict(point) for point in decision.points],
    }
    json_path.write_text(json.dumps(payload, indent=2, default=_json_default) + "\n", encoding="utf-8")
    decision.to_frame().to_csv(csv_path, index=False)
    return {"anchor_decision": json_path, "anchor_decision_csv": csv_path}


def _finite_p_values(records: Sequence[SimulationReplicateResult], statistic: str) -> list[float]:
    return [
        float(p)
        for record in records
        if (p := record.p_values.get(statistic)) is not None and math.isfinite(float(p))
    ]


def _orderable(value: Any) -> tuple[int, Any]:
    if isinstance(value, bool | int | float):
        return (0, float(value))
    return (1, str(value))


def _cell_metadata_index(records: Sequence[SimulationReplicateResult]) -> dict[str, dict]:
    index: dict[str, dict] = {}
    for record in records:
        if record.cell_id in index:
            continue
        index[record.cell_id] = dict(record.cell_metadata)
    return index


def _json_default(value: Any) -> Any:
    if isinstance(value, float) and (math.isnan(value) or math.isinf(value)):
        return None
    raise TypeError(f"Object of type {type(value).__name__} is not JSON serializable")


__all__ = [
    "ANCHOR_CLASSIFICATIONS",
    "AnchorDecision",
    "AnchorPoint",
    "AnchorStatistic",
    "ControlRate",
    "DESIGN_POINT_STATUSES",
    "RANK_CRITERIA",
    "DesignPointColumn",
    "DesignPointDecision",
    "RankColumn",
    "RankCriterion",
    "RankDecision",
    "RateEstimate",
    "TargetEvaluation",
    "anchor_bound",
    "evaluate_anchor_decision",
    "evaluate_design_point_decision",
    "evaluate_rank_decision",
    "evaluate_targets",
    "write_anchor_decision",
    "write_design_point_decision",
    "write_rank_decision",
    "write_target_report",
]
