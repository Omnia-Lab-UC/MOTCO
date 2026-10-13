"""Anchor Type I decision rule: config parsing, bound, classification, files.

Records are built by hand with known p-values so every rate is constructed:
``n_reject`` of ``n`` anchor replicates reject a statistic, the rest do not.
"""

from __future__ import annotations

import json
from pathlib import Path

import pandas as pd
import pytest

from motco.simulations.grid import SimulationReplicateResult
from motco.simulations.study import (
    AnchorDecisionRule,
    StudyConfigError,
    dump_study_config,
    evaluate_anchor_decision,
    load_study_config,
    write_anchor_decision,
)
from motco.simulations.study.enumerate import DESIGN_PHASE, DESIGN_POINT_KEY
from motco.simulations.study.targets import anchor_bound

EST = "evaluation.orientation_estimator"
BASE = {EST: None}
ND = {EST: "net_displacement"}


def _raw(acceptance: dict) -> dict:
    return {
        "generator": {"seed": 2, "trajectory_mode": "magnitude", "n_samples": 60, "p_dmp": 0.1},
        "evaluation": {"integration_method": "concat", "permutations": 0, "seed": 3},
        "trajectory_modes": ["translation"],
        "effect_sizes": [0.0, 1.0],
        "acceptance": acceptance,
    }


def _load(tmp_path: Path, raw: dict):
    path = tmp_path / "config.json"
    path.write_text(json.dumps(raw), encoding="utf-8")
    return load_study_config(path)


# ── Config ────────────────────────────────────────────────────────────────────


def test_rule_loads_and_round_trips(tmp_path: Path) -> None:
    rule = {"alpha": 0.05, "se_multiplier": 2.0, "statistics": ["delta", "angle", "shape"]}
    config = _load(tmp_path, _raw({"anchor_decision": rule}))
    assert config.acceptance.anchor_decision == AnchorDecisionRule()
    dumped = tmp_path / "dumped.json"
    dump_study_config(config, dumped)
    payload = json.loads(dumped.read_text(encoding="utf-8"))
    assert payload["acceptance"]["anchor_decision"]["statistics"] == ["delta", "angle", "shape"]
    assert load_study_config(dumped).acceptance.anchor_decision == config.acceptance.anchor_decision


def test_rule_absent_is_not_dumped(tmp_path: Path) -> None:
    config = _load(tmp_path, _raw({}))
    assert config.acceptance.anchor_decision is None
    dumped = tmp_path / "dumped.json"
    dump_study_config(config, dumped)
    assert "anchor_decision" not in json.loads(dumped.read_text(encoding="utf-8"))["acceptance"]


@pytest.mark.parametrize(
    ("rule", "message"),
    [
        ({"alpha": 0.05, "bogus": 1}, "unknown field"),
        ({"statistics": ["angle", "eigengap"]}, "unknown statistic"),
        ({"statistics": []}, "at least one statistic"),
        ({"se_multiplier": 0.0}, "se_multiplier must be positive"),
        ({"se_multiplier": -1.0}, "se_multiplier must be positive"),
        ({"alpha": 1.5}, "alpha must be between 0 and 1"),
    ],
)
def test_rule_rejections(tmp_path: Path, rule: dict, message: str) -> None:
    with pytest.raises(StudyConfigError, match=message):
        _load(tmp_path, _raw({"anchor_decision": rule}))


# ── Evaluation ────────────────────────────────────────────────────────────────


def _records(
    point: dict | None,
    phase: str,
    mode: str,
    effect: float,
    *,
    n: int,
    rejects: dict[str, int],
    status: str = "completed",
) -> list[SimulationReplicateResult]:
    meta: dict = {"effect_size": effect}
    if mode != "none" or phase != "type_i_baseline":
        meta["trajectory_mode"] = mode
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
            status=status,  # type: ignore[arg-type]
            p_values={s: (0.01 if index < rejects.get(s, 0) else 0.5) for s in ("delta", "angle", "shape")},
            cell_metadata=meta,
        )
        for index in range(n)
    ]


def test_bound_matches_the_probe_value() -> None:
    assert anchor_bound(0.05, 2.0, 500) == pytest.approx(0.0695, abs=5e-5)
    assert anchor_bound(0.05, 2.0, 100) == pytest.approx(0.0936, abs=5e-5)


def test_rate_at_the_bound_holds_and_above_it_exceeds() -> None:
    records = (
        _records(BASE, "power_primary", "none", 0.0, n=500, rejects={"angle": 34, "delta": 25, "shape": 20})
        + _records(ND, DESIGN_PHASE, "none", 0.0, n=500, rejects={"angle": 40, "delta": 25, "shape": 20})
    )
    decision = evaluate_anchor_decision(records, AnchorDecisionRule())
    baseline, nd = decision.points
    assert baseline.is_baseline and baseline.design_point == BASE
    angle = next(s for s in baseline.statistics if s.statistic == "angle")
    assert angle.rate == pytest.approx(0.068)
    assert angle.bound == pytest.approx(0.0695, abs=5e-5)
    assert angle.classification == "holds"
    assert baseline.verdict == "size_holds"
    nd_angle = next(s for s in nd.statistics if s.statistic == "angle")
    assert nd_angle.rate == pytest.approx(0.080)
    assert nd_angle.classification == "exceeds"
    assert all(s.classification == "holds" for s in nd.statistics if s.statistic != "angle")
    assert nd.verdict == "excess"
    assert decision.verdict == "excess"


def test_bound_uses_the_available_replicate_count() -> None:
    records = _records(BASE, "power_primary", "none", 0.0, n=480, rejects={"angle": 33}) + _records(
        BASE, "power_primary", "none", 0.0, n=20, rejects={}, status="failed"
    )
    (point,) = evaluate_anchor_decision(records, AnchorDecisionRule(statistics=("angle",))).points
    (angle,) = point.statistics
    assert angle.n_available == 480
    assert angle.bound == pytest.approx(anchor_bound(0.05, 2.0, 480))
    assert angle.rate == pytest.approx(33 / 480)


def test_point_without_available_replicates_is_unavailable() -> None:
    records = _records(BASE, "power_primary", "none", 0.0, n=50, rejects={})
    for record in records:
        record.p_values["angle"] = None
    decision = evaluate_anchor_decision(records, AnchorDecisionRule())
    (point,) = decision.points
    angle = next(s for s in point.statistics if s.statistic == "angle")
    assert angle.classification == "unavailable" and angle.bound is None and angle.n_available == 0
    assert point.verdict == "unavailable"
    assert decision.verdict == "unavailable"


def test_no_anchor_records() -> None:
    records = _records(BASE, "power_primary", "translation", 1.0, n=10, rejects={})
    decision = evaluate_anchor_decision(records, AnchorDecisionRule())
    assert decision.verdict == "no_anchor" and decision.points == ()


def test_other_controls_are_descriptive_and_type_i_cells_join_the_baseline() -> None:
    records = (
        _records(BASE, "power_primary", "none", 0.0, n=100, rejects={})
        + _records(BASE, "power_primary", "translation", 1.0, n=100, rejects={"delta": 90, "angle": 30})
        + _records(None, "type_i_baseline", "none", 0.0, n=100, rejects={"angle": 4})
        + _records(None, "type_i_baseline", "translation", 1.0, n=100, rejects={"delta": 88})
        + _records(ND, DESIGN_PHASE, "none", 0.0, n=100, rejects={})
        + _records(ND, DESIGN_PHASE, "translation", 1.0, n=100, rejects={"delta": 91})
    )
    decision = evaluate_anchor_decision(records, AnchorDecisionRule())
    baseline, nd = decision.points
    # The translation control's 0.30 angle rate is context only: the verdict reads the anchor.
    assert baseline.verdict == "size_holds"
    controls = {(c.phase, c.trajectory_mode): c for c in baseline.other_controls}
    assert set(controls) == {
        ("power_primary", "translation"),
        ("type_i_baseline", "none"),
        ("type_i_baseline", "translation"),
    }
    assert controls[("power_primary", "translation")].rates["angle"].rate == pytest.approx(0.30)
    assert controls[("type_i_baseline", "none")].rates["angle"].rate == pytest.approx(0.04)
    assert [(c.phase, c.trajectory_mode) for c in nd.other_controls] == [(DESIGN_PHASE, "translation")]


def test_written_files(tmp_path: Path) -> None:
    records = _records(BASE, "power_primary", "none", 0.0, n=500, rejects={"angle": 34}) + _records(
        BASE, "power_primary", "translation", 1.0, n=500, rejects={}
    )
    paths = write_anchor_decision(evaluate_anchor_decision(records, AnchorDecisionRule()), tmp_path)
    payload = json.loads(paths["anchor_decision"].read_text(encoding="utf-8"))
    assert payload["verdict"] == "size_holds"
    assert payload["rule"]["version"] == 1
    assert payload["points"][0]["design_point"] == {EST: None}
    frame = pd.read_csv(paths["anchor_decision_csv"])
    anchors = frame[frame["role"] == "anchor"]
    assert list(anchors["statistic"]) == ["delta", "angle", "shape"]
    assert list(anchors["classification"]) == ["holds", "holds", "holds"]
    assert set(frame.loc[frame["role"] == "control", "trajectory_mode"]) == {"translation"}
