"""Estimator selection over the two net-displacement power column decisions."""

from __future__ import annotations

import hashlib
import json
import sys
from pathlib import Path

import pytest


def _script():
    sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
    try:
        return __import__("net_displacement_power_follow_up")
    finally:
        sys.path.pop(0)


follow_up = _script()
CANDIDATE = dict(follow_up.CANDIDATE)


def _decision(verdict: str, *, criterion: str = "gain", failing: dict | None = None) -> dict:
    passed = {"target": True, "anchor": True, "controls": True, "protected": True}
    for name in failing or {}:
        passed[name] = False
    if verdict == "unavailable":
        passed["controls"] = None
    return {
        "verdict": verdict,
        "failing": failing or {},
        "rule": {"candidate": CANDIDATE, "criterion": criterion},
        "columns": [
            {"role": "baseline", "criteria": []},
            {
                "role": "candidate",
                "criteria": [
                    {
                        "name": name,
                        "passed": value,
                        "observations": (
                            {"target": "orientation/angle", "rate": 0.67, "baseline_rate": 0.30, "difference": 0.37}
                            if name == "target"
                            else {}
                        ),
                    }
                    for name, value in passed.items()
                ],
            },
        ],
    }


def test_candidate_qualifies_at_both_designs() -> None:
    result = follow_up.select_estimator(
        {
            "phase6": _decision("candidate_qualifies"),
            "phase5": _decision("candidate_qualifies", criterion="no_loss"),
        }
    )
    assert result["verdict"] == "adopt_net_displacement"
    assert result["rule_version"] == 1 and result["advisory"] is True
    assert result["follow_up_change"] is not None
    assert result["designs"]["phase6"]["target"]["difference"] == pytest.approx(0.37)
    assert result["designs"]["phase5"]["criterion"] == "no_loss"


@pytest.mark.parametrize("phase5_verdict", ["candidate_fails", "unavailable"])
def test_candidate_fails_at_one_design(phase5_verdict: str) -> None:
    failing = {"controls": ["translation/angle"]} if phase5_verdict == "candidate_fails" else None
    result = follow_up.select_estimator(
        {
            "phase6": _decision("candidate_qualifies"),
            "phase5": _decision(phase5_verdict, criterion="no_loss", failing=failing),
        }
    )
    assert result["verdict"] == "keep_pc1"
    assert result["follow_up_change"] is None
    phase5 = result["designs"]["phase5"]
    assert phase5["verdict"] == phase5_verdict
    if failing:
        assert phase5["failing"] == {"controls": ["translation/angle"]}
        assert "translation/angle" in result["rationale"]
    else:
        assert phase5["unavailable"] == ["controls"]
    assert result["designs"]["phase6"]["failing"] == {}


def test_failing_criteria_are_named_at_both_designs() -> None:
    result = follow_up.select_estimator(
        {
            "phase6": _decision("candidate_fails", failing={"target": ["orientation/angle"]}),
            "phase5": _decision("candidate_fails", criterion="no_loss", failing={"anchor": ["angle"]}),
        }
    )
    assert result["verdict"] == "keep_pc1"
    assert "phase6 candidate_fails (target: orientation/angle)" in result["rationale"]
    assert "phase5 candidate_fails (anchor: angle)" in result["rationale"]


@pytest.mark.parametrize(
    ("decision", "message"),
    [
        ({"columns": []}, "not a column decision"),
        ({**_decision("candidate_qualifies"), "verdict": "no_candidate"}, "verdict is 'no_candidate'"),
        (
            {**_decision("candidate_qualifies"), "rule": {"candidate": {"evaluation.integration_method": "concat"}}},
            "expected",
        ),
        ({**_decision("candidate_qualifies"), "columns": []}, "exactly one candidate column"),
    ],
)
def test_malformed_input_is_a_clear_error(decision: dict, message: str) -> None:
    with pytest.raises(follow_up.ComparisonError, match=message):
        follow_up.select_estimator({"phase6": decision, "phase5": _decision("candidate_qualifies")})


def test_missing_design_is_an_error() -> None:
    with pytest.raises(follow_up.ComparisonError, match="missing for design"):
        follow_up.select_estimator({"phase6": _decision("candidate_qualifies")})


def test_cli_writes_inputs_with_hashes(tmp_path: Path) -> None:
    paths = {}
    for design in ("phase6", "phase5"):
        path = tmp_path / f"{design}.json"
        path.write_text(json.dumps(_decision("candidate_qualifies")), encoding="utf-8")
        paths[design] = path
    out = tmp_path / "comparison_decision.json"
    assert follow_up.main(["--phase6", str(paths["phase6"]), "--phase5", str(paths["phase5"]), "--out", str(out)]) == 0
    payload = json.loads(out.read_text(encoding="utf-8"))
    assert payload["verdict"] == "adopt_net_displacement"
    for design, path in paths.items():
        assert payload["inputs"][design] == {
            "path": str(path),
            "sha256": hashlib.sha256(path.read_bytes()).hexdigest(),
        }


def test_cli_reports_unreadable_input(tmp_path: Path) -> None:
    bad = tmp_path / "bad.json"
    bad.write_text("{not json", encoding="utf-8")
    good = tmp_path / "good.json"
    good.write_text(json.dumps(_decision("candidate_qualifies")), encoding="utf-8")
    assert follow_up.main(["--phase6", str(bad), "--phase5", str(good), "--out", str(tmp_path / "o.json")]) == 1
