"""Column decision rule: config parsing and validation, criteria, verdicts, files.

Records are built by hand with known p-values so every rate is constructed:
``n_reject`` of ``n`` replicates reject a statistic, the rest do not.
"""

from __future__ import annotations

import json
from pathlib import Path

import pandas as pd
import pytest

from motco.simulations.grid import SimulationReplicateResult
from motco.simulations.study import (
    ColumnDecisionRule,
    StatisticPair,
    StudyConfigError,
    dump_study_config,
    evaluate_column_decision,
    load_study_config,
    write_column_decision,
)
from motco.simulations.study.enumerate import DESIGN_PHASE, DESIGN_POINT_KEY

METHOD = "evaluation.integration_method"
EST = "evaluation.orientation_estimator"
# Loaded (normalized) baseline point: the PC1 estimator is held as ``None``.
BASE = {METHOD: "pls", EST: None}
# Config spelling of the same point: PC1 is ``"pc1"``; ``null`` names the default.
RAW_BASE = {METHOD: "pls", EST: "pc1"}
ND = {METHOD: "pls", EST: "net_displacement"}
CONCAT = {METHOD: "concat", EST: "net_displacement"}
MODES = ("magnitude", "orientation", "shape", "translation")


def _rule_raw(**overrides) -> dict:
    rule = {
        "candidate": {EST: "net_displacement"},
        "target": {"trajectory_mode": "orientation", "statistic": "angle"},
        "criterion": "gain",
        "controls": [
            {"trajectory_mode": "translation", "statistic": "angle"},
            {"trajectory_mode": "magnitude", "statistic": "angle"},
            {"trajectory_mode": "shape", "statistic": "angle"},
        ],
        "protected": [
            {"trajectory_mode": "magnitude", "statistic": "delta"},
            {"trajectory_mode": "shape", "statistic": "shape"},
        ],
        "type_i_bound": {"alpha": 0.05, "se_tolerance": 2.0},
        "gain_se_multiplier": 2.0,
        "loss_se_multiplier": 2.0,
    }
    rule.update(overrides)
    return rule


def _raw(rule: dict | None, *, points: list[dict] | None = None, modes=MODES) -> dict:
    acceptance = {} if rule is None else {"column_decision": rule}
    return {
        "generator": {"seed": 2, "trajectory_mode": "magnitude", "n_samples": 60, "p_dmp": 0.1},
        "evaluation": {"integration_method": "pls", "permutations": 0, "seed": 3, "orientation_estimator": "pc1"},
        "trajectory_modes": list(modes),
        "effect_sizes": [0.0, 1.0],
        "design_grid": {"points": points if points is not None else [RAW_BASE, ND]},
        "acceptance": acceptance,
    }


def _load(tmp_path: Path, raw: dict):
    path = tmp_path / "config.json"
    path.write_text(json.dumps(raw), encoding="utf-8")
    return load_study_config(path)


# ── Config ────────────────────────────────────────────────────────────────────


def test_rule_loads_and_round_trips(tmp_path: Path) -> None:
    config = _load(tmp_path, _raw(_rule_raw()))
    rule = config.acceptance.column_decision
    assert isinstance(rule, ColumnDecisionRule)
    assert rule.candidate == {EST: "net_displacement"}
    assert rule.target == StatisticPair("orientation", "angle")
    assert [pair.label for pair in rule.controls] == ["translation/angle", "magnitude/angle", "shape/angle"]
    assert rule.criterion == "gain"
    dumped = tmp_path / "dumped.json"
    dump_study_config(config, dumped)
    payload = json.loads(dumped.read_text(encoding="utf-8"))
    assert payload["acceptance"]["column_decision"]["criterion"] == "gain"
    assert load_study_config(dumped).acceptance.column_decision == rule


def test_rule_absent_is_not_dumped(tmp_path: Path) -> None:
    config = _load(tmp_path, _raw(None))
    assert config.acceptance.column_decision is None
    dumped = tmp_path / "dumped.json"
    dump_study_config(config, dumped)
    assert "column_decision" not in json.loads(dumped.read_text(encoding="utf-8"))["acceptance"]


def test_protected_is_optional(tmp_path: Path) -> None:
    raw = _rule_raw()
    raw.pop("protected")
    assert _load(tmp_path, _raw(raw)).acceptance.column_decision.protected == ()


def test_full_coordinates_resolve_after_normalization(tmp_path: Path) -> None:
    # "pc1" names the PC1 baseline column, so this is rejected...
    with pytest.raises(StudyConfigError, match="matches the baseline point"):
        _load(tmp_path, _raw(_rule_raw(candidate={METHOD: "pls", EST: "pc1"})))
    # ...while the full net-displacement coordinates resolve to the one candidate,
    for estimator in ("net_displacement", None):  # null names the default, net displacement
        config = _load(tmp_path, _raw(_rule_raw(candidate={METHOD: "pls", EST: estimator})))
        assert config.matching_design_points(config.acceptance.column_decision.candidate) == [ND]


@pytest.mark.parametrize(
    ("overrides", "message"),
    [
        ({"bogus": 1}, "unknown field"),
        ({"criterion": "superiority"}, "criterion 'superiority' is unknown"),
        ({"controls": []}, "controls must name at least one pair"),
        (
            {"controls": [{"trajectory_mode": "orientation", "statistic": "angle"}]},
            "cannot also be in controls",
        ),
        (
            {"protected": [{"trajectory_mode": "orientation", "statistic": "angle"}]},
            "cannot also be in protected",
        ),
        ({"type_i_bound": {"alpha": 0.05, "k": 2}}, "type_i_bound has unknown field"),
        ({"gain_se_multiplier": -1.0}, "gain_se_multiplier must be non-negative"),
        ({"candidate": {"evaluation.integration_params.forced_components": 4}}, "absent from design_grid.axes"),
        ({"candidate": {EST: "pc2"}}, "matches no design point"),
        ({"candidate": {EST: "pc1"}}, "matches the baseline point"),
    ],
)
def test_rule_rejections(tmp_path: Path, overrides: dict, message: str) -> None:
    with pytest.raises(StudyConfigError, match=message):
        _load(tmp_path, _raw(_rule_raw(**overrides)))


def test_missing_required_field_is_rejected(tmp_path: Path) -> None:
    raw = _rule_raw()
    raw.pop("criterion")
    with pytest.raises(StudyConfigError, match="missing required field"):
        _load(tmp_path, _raw(raw))


def test_mode_absent_from_trajectory_modes_is_rejected(tmp_path: Path) -> None:
    with pytest.raises(StudyConfigError, match="'magnitude'.*absent from trajectory_modes"):
        _load(tmp_path, _raw(_rule_raw(), modes=("orientation", "shape", "translation")))


def test_candidate_matching_several_points_is_rejected(tmp_path: Path) -> None:
    with pytest.raises(StudyConfigError, match="matches 2 design points"):
        _load(tmp_path, _raw(_rule_raw(), points=[RAW_BASE, ND, CONCAT]))


def test_rule_requires_a_design_grid(tmp_path: Path) -> None:
    raw = _raw(_rule_raw())
    raw.pop("design_grid")
    with pytest.raises(StudyConfigError, match="requires a non-empty design_grid"):
        _load(tmp_path, raw)


# ── Evaluation ────────────────────────────────────────────────────────────────


def _records(
    point: dict | None,
    phase: str,
    mode: str,
    effect: float,
    *,
    n: int,
    rejects: dict[str, int],
) -> list[SimulationReplicateResult]:
    meta: dict = {"effect_size": effect, "trajectory_mode": mode}
    if point is not None:
        meta[DESIGN_POINT_KEY] = dict(point)
    cell = f"{phase}-{mode}-{effect}-{json.dumps(point)}"
    return [
        SimulationReplicateResult(
            cell_id=cell,
            phase=phase,
            replicate_index=index,
            replicate_seed=index,
            generator_seed=index,
            evaluation_seed=3,
            parameter_signature="sig",
            status="completed",
            p_values={s: (0.01 if index < rejects.get(s, 0) else 0.5) for s in ("delta", "angle", "shape")},
            cell_metadata=meta,
        )
        for index in range(n)
    ]


def _column(point: dict, phase: str, *, n: int = 500, overrides: dict | None = None) -> list:
    """A column at the top effect 1.0 (plus a lower effect that must not be read) and its anchor."""

    rejects = {
        ("none", 0.0): {"delta": 25, "angle": 25, "shape": 25},
        ("orientation", 1.0): {"angle": 150},
        ("orientation", 0.5): {"angle": 500},
        ("translation", 1.0): {"delta": 500, "angle": 30},
        ("magnitude", 1.0): {"delta": 500, "angle": 25},
        ("shape", 1.0): {"shape": 400, "angle": 25},
    }
    rejects.update(overrides or {})
    out: list = []
    for (mode, effect), counts in rejects.items():
        out += _records(point, phase, mode, effect, n=n, rejects=counts)
    return out


def _rule(**overrides) -> ColumnDecisionRule:
    fields = {
        "candidate": {EST: "net_displacement"},
        "target": StatisticPair("orientation", "angle"),
        "criterion": "gain",
        "controls": (
            StatisticPair("translation", "angle"),
            StatisticPair("magnitude", "angle"),
            StatisticPair("shape", "angle"),
        ),
        "protected": (StatisticPair("magnitude", "delta"), StatisticPair("shape", "shape")),
    }
    fields.update(overrides)
    return ColumnDecisionRule(**fields)


def _criterion(decision, name):
    return decision.column("candidate").criteria[[c.name for c in decision.column("candidate").criteria].index(name)]


def test_gain_qualifies() -> None:
    records = _column(BASE, "power_primary") + _column(
        ND, DESIGN_PHASE, overrides={("orientation", 1.0): {"angle": 335}}
    )
    decision = evaluate_column_decision(records, _rule())
    target = _criterion(decision, "target")
    assert decision.column("baseline").target.rate == pytest.approx(0.300)
    assert target.observations["rate"] == pytest.approx(0.670)
    assert target.observations["difference"] == pytest.approx(0.370)
    assert target.observations["pooled_se"] == pytest.approx(0.0294, abs=5e-5)
    assert target.passed is True
    assert decision.verdict == "candidate_qualifies"
    assert decision.column("candidate").qualifies is True
    assert decision.failing == {}


def test_gain_below_threshold_fails_the_target() -> None:
    records = _column(BASE, "power_primary") + _column(
        ND, DESIGN_PHASE, overrides={("orientation", 1.0): {"angle": 170}}
    )
    decision = evaluate_column_decision(records, _rule())
    assert _criterion(decision, "target").passed is False
    assert decision.verdict == "candidate_fails"
    assert decision.failing == {"target": ["orientation/angle"]}


def test_control_above_the_bound_fails_whatever_the_rest() -> None:
    records = _column(BASE, "power_primary") + _column(
        ND,
        DESIGN_PHASE,
        overrides={("orientation", 1.0): {"angle": 335}, ("translation", 1.0): {"delta": 500, "angle": 40}},
    )
    decision = evaluate_column_decision(records, _rule())
    controls = _criterion(decision, "controls")
    assert controls.passed is False
    assert controls.failing == ("translation/angle",)
    assert controls.observations["translation_angle_rate"] == pytest.approx(0.080)
    assert controls.observations["translation_angle_bound"] == pytest.approx(0.0695, abs=5e-5)
    assert _criterion(decision, "target").passed is True
    assert decision.verdict == "candidate_fails"
    assert decision.failing == {"controls": ["translation/angle"]}


def test_baseline_leak_does_not_judge_the_candidate() -> None:
    # PC1's shape → angle leak (0.89 at Phase 5) is the baseline's; only the candidate is judged.
    records = _column(BASE, "power_primary", overrides={("shape", 1.0): {"shape": 400, "angle": 445}}) + _column(
        ND, DESIGN_PHASE, overrides={("orientation", 1.0): {"angle": 335}}
    )
    decision = evaluate_column_decision(records, _rule())
    assert decision.verdict == "candidate_qualifies"
    assert decision.column("baseline").controls["shape/angle"].rate == pytest.approx(0.89)


def test_no_loss_qualifies() -> None:
    records = _column(BASE, "power_primary", overrides={("orientation", 1.0): {"angle": 430}}) + _column(
        ND, DESIGN_PHASE, overrides={("orientation", 1.0): {"angle": 420}}
    )
    decision = evaluate_column_decision(records, _rule(criterion="no_loss"))
    target = _criterion(decision, "target")
    assert target.observations["difference"] == pytest.approx(0.020)
    assert target.observations["difference"] <= target.observations["threshold"]
    assert target.passed is True
    assert decision.verdict == "candidate_qualifies"


def test_no_loss_fails_on_a_large_loss() -> None:
    records = _column(BASE, "power_primary", overrides={("orientation", 1.0): {"angle": 430}}) + _column(
        ND, DESIGN_PHASE, overrides={("orientation", 1.0): {"angle": 350}}
    )
    decision = evaluate_column_decision(records, _rule(criterion="no_loss"))
    assert _criterion(decision, "target").passed is False
    assert decision.verdict == "candidate_fails"


def test_anchor_and_protected_failures_are_named() -> None:
    records = _column(BASE, "power_primary") + _column(
        ND,
        DESIGN_PHASE,
        overrides={
            ("orientation", 1.0): {"angle": 335},
            ("none", 0.0): {"delta": 25, "angle": 45, "shape": 25},
            ("shape", 1.0): {"shape": 300, "angle": 25},
        },
    )
    decision = evaluate_column_decision(records, _rule())
    assert decision.failing == {"anchor": ["angle"], "protected": ["shape/shape"]}
    assert decision.verdict == "candidate_fails"


def test_unavailable_control_yields_unavailable() -> None:
    candidate = _column(ND, DESIGN_PHASE, overrides={("orientation", 1.0): {"angle": 335}})
    candidate = [r for r in candidate if r.cell_metadata["trajectory_mode"] != "magnitude"]
    decision = evaluate_column_decision(_column(BASE, "power_primary") + candidate, _rule())
    controls = _criterion(decision, "controls")
    assert controls.passed is None
    assert decision.verdict == "unavailable"
    assert decision.column("candidate").qualifies is None


def test_bound_uses_the_available_replicates() -> None:
    candidate = _column(ND, DESIGN_PHASE, overrides={("orientation", 1.0): {"angle": 335}})
    for record in candidate:
        if record.cell_metadata["trajectory_mode"] == "translation" and record.replicate_index >= 400:
            record.p_values["angle"] = float("nan")
    decision = evaluate_column_decision(_column(BASE, "power_primary") + candidate, _rule())
    controls = _criterion(decision, "controls")
    assert controls.observations["translation_angle_n"] == 400
    assert controls.observations["translation_angle_bound"] == pytest.approx(0.05 + 2 * (0.05 * 0.95 / 400) ** 0.5)


def test_no_candidate_and_no_design_grid() -> None:
    assert evaluate_column_decision(_column(BASE, "power_primary"), _rule()).verdict == "no_design_grid"
    records = _column(BASE, "power_primary") + _column(CONCAT, DESIGN_PHASE)
    decision = evaluate_column_decision(records, _rule(candidate={METHOD: "pls", EST: "net_displacement"}))
    assert decision.verdict == "no_candidate"
    assert [column.role for column in decision.columns] == ["baseline", "other"]


def test_other_columns_are_measured_without_a_verdict(tmp_path: Path) -> None:
    records = (
        _column(BASE, "power_primary")
        + _column(ND, DESIGN_PHASE, overrides={("orientation", 1.0): {"angle": 335}})
        + _column(CONCAT, DESIGN_PHASE, overrides={("orientation", 1.0): {"angle": 500}})
    )
    decision = evaluate_column_decision(records, _rule(candidate={METHOD: "pls", EST: "net_displacement"}))
    assert [column.role for column in decision.columns] == ["baseline", "candidate", "other"]
    other = decision.column("other")
    assert other.criteria == () and other.qualifies is None
    assert other.target.rate == pytest.approx(1.0)

    paths = write_column_decision(decision, tmp_path)
    payload = json.loads(paths["column_decision"].read_text(encoding="utf-8"))
    assert payload["verdict"] == "candidate_qualifies"
    assert payload["advisory"] is True
    assert payload["rule"]["candidate"] == {METHOD: "pls", EST: "net_displacement"}
    assert [column["role"] for column in payload["columns"]] == ["baseline", "candidate", "other"]
    frame = pd.read_csv(paths["column_decision_csv"])
    assert list(frame["role"]) == ["baseline", "candidate", "other"]
    assert frame.loc[1, "criterion_controls_translation_angle_bound"] == pytest.approx(0.0695, abs=5e-5)
    assert frame.loc[1, "target_rate"] == pytest.approx(0.670)
    assert frame.loc[1, "criterion_target_difference"] == pytest.approx(0.370)
    assert frame.loc[1, "qualifies"] and pd.isna(frame.loc[0, "qualifies"])
