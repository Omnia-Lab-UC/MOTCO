"""Follow-up selection over the two net-displacement anchor decisions."""

from __future__ import annotations

import hashlib
import json
import sys
from pathlib import Path

import pytest


def _script():
    sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
    try:
        return __import__("net_displacement_follow_up")
    finally:
        sys.path.pop(0)


follow_up = _script()
M, F, N, E = follow_up.METHOD, follow_up.FORCED, follow_up.NCOMP, follow_up.ESTIMATOR
COLUMNS = {
    "baseline": {M: "pls", F: None, N: None, E: None},
    "lead": {M: "pls", F: None, N: None, E: "net_displacement"},
    "forced": {M: "pls", F: 4, N: None, E: "net_displacement"},
    "pca": {M: "pca", F: None, N: 4, E: "net_displacement"},
    "concat": {M: "concat", F: None, N: None, E: "net_displacement"},
}


def _decision(verdicts: dict[str, str]) -> dict:
    return {
        "verdict": "excess" if "excess" in verdicts.values() else "size_holds",
        "points": [
            {
                "design_point": COLUMNS[name],
                "verdict": verdict,
                "statistics": [
                    {"statistic": "angle", "rate": 0.06, "monte_carlo_se": 0.01, "n_available": 500,
                     "bound": 0.0695, "classification": "holds" if verdict == "size_holds" else "exceeds"}
                ],
            }
            for name, verdict in verdicts.items()
        ],
    }


HOLDS = {name: "size_holds" for name in COLUMNS}


def test_lead_holds_at_both_designs() -> None:
    result = follow_up.select_follow_up({"phase6": _decision(HOLDS), "phase5": _decision(HOLDS)})
    assert result["verdict"] == "paper_grade_power_comparison"
    assert result["rule_version"] == 1
    assert {d["localization"] for d in result["designs"].values()} == {"none"}
    assert len(result["designs"]["phase6"]["net_displacement_columns"]) == 4


def test_lead_exceeds_at_one_design_with_localization() -> None:
    some = {**HOLDS, "lead": "excess", "pca": "excess"}
    every = {**HOLDS, "lead": "excess", "forced": "excess", "pca": "excess", "concat": "excess"}
    result = follow_up.select_follow_up({"phase6": _decision(some), "phase5": _decision(HOLDS)})
    assert result["verdict"] == "diagnose_anchor_excess"
    assert result["designs"]["phase6"]["localization"] == "some_columns"
    assert result["designs"]["phase5"]["localization"] == "none"
    result = follow_up.select_follow_up({"phase6": _decision(HOLDS), "phase5": _decision(every)})
    assert result["verdict"] == "diagnose_anchor_excess"
    assert result["designs"]["phase5"]["localization"] == "every_column"


def test_missing_lead_column_is_a_clear_error() -> None:
    no_lead = {name: verdict for name, verdict in HOLDS.items() if name != "lead"}
    with pytest.raises(follow_up.FollowUpError, match="phase5: expected exactly one lead column"):
        follow_up.select_follow_up({"phase6": _decision(HOLDS), "phase5": _decision(no_lead)})


def test_unavailable_lead_is_a_clear_error() -> None:
    with pytest.raises(follow_up.FollowUpError, match="unavailable anchor is a failed run"):
        follow_up.select_follow_up({"phase6": _decision({**HOLDS, "lead": "unavailable"}), "phase5": _decision(HOLDS)})


def test_cli_writes_every_declared_field(tmp_path: Path) -> None:
    paths = {}
    for design in ("phase6", "phase5"):
        paths[design] = tmp_path / f"{design}.json"
        paths[design].write_text(json.dumps(_decision(HOLDS)), encoding="utf-8")
    out = tmp_path / "follow_up_decision.json"
    assert follow_up.main(["--phase6", str(paths["phase6"]), "--phase5", str(paths["phase5"]), "--out", str(out)]) == 0
    payload = json.loads(out.read_text(encoding="utf-8"))
    assert set(payload) == {"rule_version", "lead_column", "verdict", "rationale", "designs", "advisory", "inputs"}
    for design, path in paths.items():
        assert payload["inputs"][design] == {
            "path": str(path),
            "sha256": hashlib.sha256(path.read_bytes()).hexdigest(),
        }
