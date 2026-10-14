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
