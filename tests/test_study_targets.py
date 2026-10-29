from __future__ import annotations

import json
import math
from pathlib import Path

import pandas as pd
import pytest

from motco.simulations import SimulationReplicateResult, SimulationSummaryResult
from motco.simulations.study.config import (
    AcceptanceTargets,
    PowerMonotonicityTarget,
    SpecificityTarget,
    TypeIControlTarget,
)
from motco.simulations.study.phase4 import type_i_inflation_bound
from motco.simulations.study.targets import evaluate_targets, write_target_report


def _summary(
    cell_id: str,
    phase: str,
    statistic: str,
    rate: float,
    se: float,
    available: int = 100,
) -> SimulationSummaryResult:
    return SimulationSummaryResult(
        cell_id=cell_id,
        phase=phase,
        statistic=statistic,
        alpha=0.05,
        completed_replicates=available,
        available_replicates=available,
        rejected_replicates=int(rate * available),
        rejection_rate=rate,
        monte_carlo_se=se,
        unavailable_replicates=0,
    )


def _record(
    cell_id: str,
    phase: str,
    mode: str | None,
    effect_size: float | None,
    varied_axis=None,
) -> SimulationReplicateResult:
    metadata = {"varied_axis": varied_axis}
    if mode is not None:
        metadata["trajectory_mode"] = mode
    if effect_size is not None:
        metadata["effect_size"] = effect_size
    return SimulationReplicateResult(
        cell_id=cell_id,
        phase=phase,
        replicate_index=0,
        replicate_seed=0,
        generator_seed=0,
        evaluation_seed=0,
        parameter_signature="sig",
        status="completed",
        cell_metadata=metadata,
    )


def test_type_i_control_target_met_and_not_met() -> None:
    summaries = [
        _summary("null-met", "type_i_baseline", "delta", 0.06, 0.02),
        _summary("null-fail", "type_i_baseline", "delta", 0.20, 0.02),
    ]
    records = [
        _record("null-met", "type_i_baseline", mode=None, effect_size=None),
        _record("null-fail", "type_i_baseline", mode=None, effect_size=None),
    ]
    acceptance = AcceptanceTargets(type_i=(TypeIControlTarget(alpha=0.05, se_tolerance=2.0),))
    evaluations = evaluate_targets(acceptance, summaries, records)
    by_cell = {ev.observations["cell_id"]: ev for ev in evaluations}
    assert by_cell["null-met"].met is True
    assert by_cell["null-fail"].met is False


def test_power_monotonicity_target_met() -> None:
    summaries = [
        _summary("magn-0.1", "power_primary", "delta", 0.20, 0.04),
        _summary("magn-0.5", "power_primary", "delta", 0.85, 0.03),
    ]
    records = [
        _record("magn-0.1", "power_primary", mode="magnitude", effect_size=0.1),
        _record("magn-0.5", "power_primary", mode="magnitude", effect_size=0.5),
    ]
    acceptance = AcceptanceTargets(
        power=(PowerMonotonicityTarget(trajectory_mode="magnitude", statistic="delta", min_power_at_top=0.80),)
    )
    [evaluation] = evaluate_targets(acceptance, summaries, records)
    assert evaluation.met is True
    assert evaluation.observations["monotone"] is True
    assert evaluation.observations["top_rate"] == 0.85


def test_power_monotonicity_target_not_met_due_to_floor() -> None:
    summaries = [
        _summary("magn-0.1", "power_primary", "delta", 0.20, 0.04),
        _summary("magn-0.5", "power_primary", "delta", 0.55, 0.03),
    ]
    records = [
        _record("magn-0.1", "power_primary", mode="magnitude", effect_size=0.1),
        _record("magn-0.5", "power_primary", mode="magnitude", effect_size=0.5),
    ]
    acceptance = AcceptanceTargets(
        power=(PowerMonotonicityTarget(trajectory_mode="magnitude", statistic="delta", min_power_at_top=0.80),)
    )
    [evaluation] = evaluate_targets(acceptance, summaries, records)
    assert evaluation.met is False


def test_power_monotonicity_target_not_met_due_to_non_monotone() -> None:
    summaries = [
        _summary("magn-0.1", "power_primary", "delta", 0.95, 0.02),
        _summary("magn-0.5", "power_primary", "delta", 0.40, 0.03),
    ]
    records = [
        _record("magn-0.1", "power_primary", mode="magnitude", effect_size=0.1),
        _record("magn-0.5", "power_primary", mode="magnitude", effect_size=0.5),
    ]
    acceptance = AcceptanceTargets(
        power=(PowerMonotonicityTarget(trajectory_mode="magnitude", statistic="delta", min_power_at_top=0.30),)
    )
    [evaluation] = evaluate_targets(acceptance, summaries, records)
    assert evaluation.met is False
    assert evaluation.observations["monotone"] is False


def test_specificity_target_uses_largest_effect_size_and_met() -> None:
    summaries = [
        _summary("translation-0.1", "power_primary", "angle", 0.06, 0.02),
        _summary("translation-0.5", "power_primary", "angle", 0.07, 0.02),
    ]
    records = [
        _record("translation-0.1", "power_primary", mode="translation", effect_size=0.1),
        _record("translation-0.5", "power_primary", mode="translation", effect_size=0.5),
    ]
    acceptance = AcceptanceTargets(
        specificity=(SpecificityTarget(trajectory_mode="translation", statistic="angle", alpha=0.05, se_tolerance=2.0),)
    )
    [evaluation] = evaluate_targets(acceptance, summaries, records)
    assert evaluation.met is True
    assert evaluation.observations["effect_size"] == 0.5


def test_specificity_target_not_met_when_far_from_alpha() -> None:
    summaries = [_summary("translation-0.5", "power_primary", "angle", 0.50, 0.02)]
    records = [_record("translation-0.5", "power_primary", mode="translation", effect_size=0.5)]
    acceptance = AcceptanceTargets(
        specificity=(SpecificityTarget(trajectory_mode="translation", statistic="angle", alpha=0.05, se_tolerance=2.0),)
    )
    [evaluation] = evaluate_targets(acceptance, summaries, records)
    assert evaluation.met is False


def test_write_target_report_emits_csv_and_json(tmp_path: Path) -> None:
    summaries = [_summary("null-met", "type_i_baseline", "delta", 0.05, 0.02)]
    records = [_record("null-met", "type_i_baseline", mode=None, effect_size=None)]
    acceptance = AcceptanceTargets(type_i=(TypeIControlTarget(alpha=0.05),))
    evaluations = evaluate_targets(acceptance, summaries, records)
    paths = write_target_report(evaluations, tmp_path)
    assert paths["csv"].exists()
    assert paths["json"].exists()
    payload = json.loads(paths["json"].read_text(encoding="utf-8"))
    assert payload[0]["target_kind"] == "type_i_control"


# --- one-sided (upper) bound form -------------------------------------------------

_PRE_CHANGE_TYPE_I_KEYS = [
    "cell_id",
    "phase",
    "trajectory_mode",
    "statistic",
    "alpha_target",
    "se_tolerance",
    "rejection_rate",
    "monte_carlo_se",
    "available_replicates",
    "deviation",
    "bound",
]


def _empirical_se(rate: float, n: int) -> float:
    return math.sqrt(rate * (1.0 - rate) / n)


def _type_i_eval(rate: float, bound: str, n: int = 500):
    summaries = [_summary("null", "type_i_baseline", "shape", rate, _empirical_se(rate, n), available=n)]
    records = [_record("null", "type_i_baseline", mode=None, effect_size=None)]
    target = TypeIControlTarget(alpha=0.05, se_tolerance=2.0, bound=bound)
    [evaluation] = evaluate_targets(AcceptanceTargets(type_i=(target,)), summaries, records)
    return evaluation


def _specificity_eval(rate: float, bound: str, n: int = 500):
    summaries = [_summary("magnitude-1.0", "power_primary", "shape", rate, _empirical_se(rate, n), available=n)]
    records = [_record("magnitude-1.0", "power_primary", mode="magnitude", effect_size=1.0)]
    target = SpecificityTarget(
        trajectory_mode="magnitude", statistic="shape", alpha=0.05, se_tolerance=2.0, bound=bound
    )
    [evaluation] = evaluate_targets(AcceptanceTargets(specificity=(target,)), summaries, records)
    return evaluation


@pytest.mark.parametrize("evaluate", [_type_i_eval, _specificity_eval])
def test_conservative_control_is_met_under_upper_bound(evaluate) -> None:
    evaluation = evaluate(0.000, "upper")
    assert evaluation.met is True
    obs = evaluation.observations
    assert obs["bound"] == pytest.approx(0.0695, abs=5e-5)
    assert obs["deviation"] == pytest.approx(-0.050)
    assert obs["bound_kind"] == "upper"
    assert obs["alpha_se"] == pytest.approx(math.sqrt(0.05 * 0.95 / 500))
    assert obs["bound"] == pytest.approx(type_i_inflation_bound(0.05, 2.0, 500))
    assert evaluation.rationale == "rate = 0.0000 <= alpha + 2.00·SE_alpha (0.0695)"


@pytest.mark.parametrize("evaluate", [_type_i_eval, _specificity_eval])
def test_zero_rejection_control_is_not_met_under_two_sided_bound(evaluate) -> None:
    evaluation = evaluate(0.000, "two_sided")
    assert evaluation.met is False
    assert evaluation.observations["bound"] == 0.0


@pytest.mark.parametrize("evaluate", [_type_i_eval, _specificity_eval])
def test_inflated_control_is_not_met_under_upper_bound(evaluate) -> None:
    evaluation = evaluate(0.080, "upper")
    assert evaluation.met is False
    assert evaluation.observations["deviation"] == pytest.approx(0.030)
    assert evaluation.rationale == "rate = 0.0800 > alpha + 2.00·SE_alpha (0.0695)"


@pytest.mark.parametrize("evaluate", [_type_i_eval, _specificity_eval])
def test_conservative_nonzero_control_is_met_under_upper_bound(evaluate) -> None:
    evaluation = evaluate(0.024, "upper")
    assert evaluation.met is True
    assert evaluation.observations["deviation"] == pytest.approx(-0.026)


def test_default_form_observations_are_unchanged() -> None:
    evaluation = _type_i_eval(0.024, "two_sided")
    assert list(evaluation.observations) == _PRE_CHANGE_TYPE_I_KEYS
    assert evaluation.rationale.startswith("|rate - alpha| = 0.0260 > 2.00·SE (")
    spec = _specificity_eval(0.024, "two_sided")
    assert "bound_kind" not in spec.observations
    assert "alpha_se" not in spec.observations
    assert spec.observations["deviation"] == pytest.approx(0.026)


def test_mixed_report_adds_bound_kind_columns_blank_on_two_sided_rows(tmp_path: Path) -> None:
    summaries = [
        _summary("null", "type_i_baseline", "shape", 0.0, 0.0, available=500),
        _summary("magnitude-1.0", "power_primary", "shape", 0.0, 0.0, available=500),
    ]
    records = [
        _record("null", "type_i_baseline", mode=None, effect_size=None),
        _record("magnitude-1.0", "power_primary", mode="magnitude", effect_size=1.0),
    ]
    acceptance = AcceptanceTargets(
        type_i=(TypeIControlTarget(alpha=0.05),),
        specificity=(SpecificityTarget(trajectory_mode="magnitude", statistic="shape", alpha=0.05, bound="upper"),),
    )
    paths = write_target_report(evaluate_targets(acceptance, summaries, records), tmp_path)
    frame = pd.read_csv(paths["csv"])
    assert {"obs_bound_kind", "obs_alpha_se"} <= set(frame.columns)
    two_sided = frame[frame["target_kind"] == "type_i_control"].iloc[0]
    upper = frame[frame["target_kind"] == "specificity"].iloc[0]
    assert pd.isna(two_sided["obs_bound_kind"]) and pd.isna(two_sided["obs_alpha_se"])
    assert upper["obs_bound_kind"] == "upper"
    payload = json.loads(paths["json"].read_text(encoding="utf-8"))
    assert "bound_kind" not in payload[0]["observations"]
    assert payload[1]["observations"]["bound_kind"] == "upper"


def test_two_sided_only_report_has_no_bound_kind_column(tmp_path: Path) -> None:
    summaries = [_summary("null", "type_i_baseline", "delta", 0.05, 0.02)]
    records = [_record("null", "type_i_baseline", mode=None, effect_size=None)]
    acceptance = AcceptanceTargets(type_i=(TypeIControlTarget(alpha=0.05),))
    paths = write_target_report(evaluate_targets(acceptance, summaries, records), tmp_path)
    columns = set(pd.read_csv(paths["csv"]).columns)
    assert "obs_bound_kind" not in columns and "obs_alpha_se" not in columns


# --- column-scoped Type I entries -------------------------------------------------

from motco.simulations.study.config import AnchorDecisionRule  # noqa: E402
from motco.simulations.study.enumerate import DESIGN_AXIS_MARKER, DESIGN_PHASE, DESIGN_POINT_KEY  # noqa: E402
from motco.simulations.study.summary import summarize_study  # noqa: E402
from motco.simulations.study.targets import evaluate_anchor_decision  # noqa: E402

N = "generator.n_samples"
LAYERS = "evaluation.integration_params.layers"
TWO = ["methylation", "expression"]
COLUMNS = [{N: n, LAYERS: layers} for n in (80, 200, 500) for layers in (TWO, None)]
BASELINE = {N: 200, LAYERS: TWO}
R = 500


def _column_record(
    point: dict, mode: str, effect: float, index: int, rejected: dict[str, bool]
) -> SimulationReplicateResult:
    phase = "power_primary" if point == BASELINE else DESIGN_PHASE
    blocks = 3 if point[LAYERS] is None else len(point[LAYERS])
    return SimulationReplicateResult(
        cell_id=f"{phase}-{mode}-{effect}-{point[N]}-{blocks}",
        phase=phase,
        replicate_index=index,
        replicate_seed=index,
        generator_seed=index,
        evaluation_seed=3,
        parameter_signature="sig",
        status="completed",
        p_values={s: (0.01 if rejected[s] else 0.5) for s in ("delta", "angle", "shape")},
        cell_metadata={
            "trajectory_mode": mode,
            "effect_size": effect,
            DESIGN_POINT_KEY: dict(point),
            "varied_axis": DESIGN_AXIS_MARKER if phase == DESIGN_PHASE else None,
        },
    )


def _type_i_baseline_record(mode: str, effect: float, index: int) -> SimulationReplicateResult:
    return SimulationReplicateResult(
        cell_id=f"type_i_baseline-{mode}",
        phase="type_i_baseline",
        replicate_index=index,
        replicate_seed=index,
        generator_seed=10_000 + index,
        evaluation_seed=3,
        parameter_signature="sig",
        status="completed",
        p_values={"delta": 0.5, "angle": 0.01, "shape": 0.5},
        cell_metadata={"trajectory_mode": mode, "effect_size": effect, "varied_axis": None},
    )


def _rejections(point: dict, mode: str, effect: float) -> dict[str, int]:
    """Rejections out of R: anchors near alpha (one column over the bound), translation angle rising with n."""

    blocks = 3 if point[LAYERS] is None else len(point[LAYERS])
    if mode == "none":
        angle = 40 if (point[N], blocks) == (500, 3) else 25  # 0.080 exceeds 0.0695; 0.050 holds
        return {"delta": 20, "angle": angle, "shape": 30}
    if mode == "translation":
        return {"delta": 20, "angle": {0.5: 15, 1.0: 30}[effect] + (5 if point[N] == 500 else 0), "shape": 25}
    return {"delta": int(R * effect), "angle": int(R * effect / 2), "shape": int(R * effect / 4)}


def _six_column_records(*, drop_translation_at: dict | None = None) -> list[SimulationReplicateResult]:
    records: list[SimulationReplicateResult] = []
    for point in COLUMNS:
        cells = [("none", 0.0), ("orientation", 0.5), ("orientation", 1.0), ("translation", 0.5), ("translation", 1.0)]
        for mode, effect in cells:
            if mode == "translation" and point == drop_translation_at:
                continue
            counts = _rejections(point, mode, effect)
            for i in range(R):
                records.append(_column_record(point, mode, effect, i, {s: i < counts[s] for s in counts}))
    for i in range(R):
        records.append(_type_i_baseline_record("none", 0.0, i))
        records.append(_type_i_baseline_record("translation", 1.0, i))
    return records


def _evaluate(target: TypeIControlTarget, records: list[SimulationReplicateResult]) -> list:
    summaries = summarize_study(records, alpha=0.05)
    return evaluate_targets(AcceptanceTargets(type_i=(target,)), summaries, records)


def _column_key(point: dict) -> tuple:
    return (point[N], None if point[LAYERS] is None else tuple(point[LAYERS]))


def test_column_scoped_translation_angle_entry_is_evaluated_at_every_column() -> None:
    records = _six_column_records()
    target = TypeIControlTarget(
        alpha=0.05, se_tolerance=2.0, bound="upper", trajectory_mode="translation", statistic="angle"
    )
    evaluations = _evaluate(target, records)
    assert len(evaluations) == 6
    assert all(ev.target_kind == "type_i_control" for ev in evaluations)
    assert all(ev.observations["statistic"] == "angle" for ev in evaluations)
    assert all(ev.observations["trajectory_mode"] == "translation" for ev in evaluations)
    assert all(ev.observations["effect_size"] == 1.0 for ev in evaluations)  # the largest effect, never 0.5
    assert [ev.observations["is_baseline"] for ev in evaluations] == [True] + [False] * 5
    assert evaluations[0].observations["design_point"] == BASELINE
    assert evaluations[0].observations["phase"] == "power_primary"
    assert {_column_key(ev.observations["design_point"]) for ev in evaluations} == {_column_key(p) for p in COLUMNS}
    for ev in evaluations:
        obs = ev.observations
        assert obs["available_replicates"] == R
        assert obs["bound"] == pytest.approx(type_i_inflation_bound(0.05, 2.0, R))
        assert obs["bound"] == pytest.approx(0.0695, abs=5e-5)
        assert obs["bound_kind"] == "upper"
        expected_rate = _rejections(obs["design_point"], "translation", 1.0)["angle"] / R
        assert obs["rejection_rate"] == pytest.approx(expected_rate)
        assert ev.met is (expected_rate <= obs["bound"])
        assert ev.target_name.startswith("type_i_control[translation@1.0, angle, {")
        assert json.loads(ev.target_name[len("type_i_control[translation@1.0, angle, "):-1]) == obs["design_point"]
    # n = 500 columns: 35/500 = 0.070 > 0.0695; the other four hold at 0.060.
    assert sum(ev.met is False for ev in evaluations) == 2
    assert {ev.observations["design_point"][N] for ev in evaluations if ev.met is False} == {500}


def test_column_scoped_none_entry_without_statistic_covers_every_statistic_and_not_the_type_i_phase() -> None:
    records = _six_column_records()
    target = TypeIControlTarget(alpha=0.05, se_tolerance=2.0, bound="upper", trajectory_mode="none")
    evaluations = _evaluate(target, records)
    assert len(evaluations) == 6 * 3
    assert {ev.observations["phase"] for ev in evaluations} == {"power_primary", DESIGN_PHASE}
    assert all(ev.observations["effect_size"] == 0.0 for ev in evaluations)
    per_column = {}
    for ev in evaluations:
        per_column.setdefault(_column_key(ev.observations["design_point"]), []).append(ev.observations["statistic"])
    assert all(stats == ["delta", "angle", "shape"] for stats in per_column.values())
    assert len(per_column) == 6
    # The Type I baseline anchor (angle rejects every replicate) is never read.
    assert all(ev.observations["rejection_rate"] < 0.1 for ev in evaluations)

    # Two-sided form on the same entry: the phase-scoped observations without bound_kind.
    two_sided = _evaluate(TypeIControlTarget(alpha=0.05, trajectory_mode="none"), records)
    assert len(two_sided) == 18
    assert all("bound_kind" not in ev.observations for ev in two_sided)
    assert all(ev.observations["deviation"] >= 0 for ev in two_sided)


def test_column_scoped_anchor_evaluations_agree_with_the_anchor_decision() -> None:
    records = _six_column_records()
    rule = AnchorDecisionRule(alpha=0.05, se_multiplier=2.0, statistics=("delta", "angle", "shape"))
    decision = evaluate_anchor_decision(records, rule, alpha=0.05)
    target = TypeIControlTarget(alpha=0.05, se_tolerance=2.0, bound="upper", trajectory_mode="none")
    evaluations = _evaluate(target, records)
    anchor = {
        (_column_key(point.design_point), stat.statistic): (
            point.is_baseline,
            stat.rate,
            stat.n_available,
            stat.bound,
            stat.classification,
        )
        for point in decision.points
        for stat in point.statistics
    }
    assert len(anchor) == 18
    for ev in evaluations:
        obs = ev.observations
        key = (_column_key(obs["design_point"]), obs["statistic"])
        is_baseline, rate, n_available, bound, classification = anchor[key]
        assert obs["is_baseline"] is is_baseline
        assert obs["rejection_rate"] == pytest.approx(rate)
        assert obs["available_replicates"] == n_available
        assert obs["bound"] == pytest.approx(bound)
        assert ev.met is (classification == "holds")
    # The construction puts one anchor over the bound; both readers see it at the same column.
    assert [ev.observations["design_point"] for ev in evaluations if ev.met is False] == [{N: 500, LAYERS: None}]
    assert decision.verdict == "excess"


def test_column_without_the_control_is_reported_unavailable() -> None:
    missing = {N: 80, LAYERS: None}
    records = _six_column_records(drop_translation_at=missing)
    target = TypeIControlTarget(alpha=0.05, bound="upper", trajectory_mode="translation", statistic="angle")
    evaluations = _evaluate(target, records)
    assert len(evaluations) == 6
    [unavailable] = [ev for ev in evaluations if ev.met is None]
    assert unavailable.observations["design_point"] == missing
    assert unavailable.observations["cell_id"] is None
    assert unavailable.observations["available_replicates"] == 0
    assert "No available replicates" in unavailable.rationale
    assert unavailable.target_name == f"type_i_control[translation, angle, {json.dumps(missing, sort_keys=True)}]"
    assert sum(ev.met is not None for ev in evaluations) == 5


def test_column_scoped_report_rows_carry_the_column_columns(tmp_path: Path) -> None:
    records = _six_column_records()
    acceptance = AcceptanceTargets(
        type_i=(
            TypeIControlTarget(alpha=0.05),
            TypeIControlTarget(alpha=0.05, bound="upper", trajectory_mode="translation", statistic="angle"),
        )
    )
    paths = write_target_report(evaluate_targets(acceptance, summarize_study(records, alpha=0.05), records), tmp_path)
    frame = pd.read_csv(paths["csv"])
    assert {"obs_design_point", "obs_is_baseline"} <= set(frame.columns)
    phase_scoped = frame[frame["obs_phase"] == "type_i_baseline"]
    assert len(phase_scoped) == 6  # two Type I baseline cells x three statistics, as before
    assert phase_scoped["obs_design_point"].isna().all()
    column_scoped = frame[frame["obs_trajectory_mode"] == "translation"]
    column_scoped = column_scoped[column_scoped["obs_phase"] != "type_i_baseline"]
    assert len(column_scoped) == 6
    assert column_scoped["obs_is_baseline"].astype(bool).sum() == 1
    payload = json.loads(paths["json"].read_text(encoding="utf-8"))
    column_rows = [row for row in payload if row["observations"].get("design_point") is not None]
    assert len(column_rows) == 6
    assert column_rows[0]["observations"]["design_point"] == BASELINE
