"""``scripts/rank_floor_check.py``: reproducibility, arms, resume, merge and the decision rule."""

from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import pytest


def _script():
    sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
    try:
        return __import__("rank_floor_check")
    finally:
        sys.path.pop(0)


check = _script()
SMALL_CV = {"cv1_splits": 2, "cv2_splits": 3, "n_repeats": 2, "max_components": 5, "random_state": 1203}
#: One effect per mode keeps the fast test at five cells.
SMALL_EFFECTS = {"magnitude": [0.25], "orientation": [1.0], "shape": [1.0]}


def _settings(out_dir: Path) -> dict:
    args = check.build_parser().parse_args(["run", "--reps", "2", "--perms", "9", "--out-dir", str(out_dir)])
    settings = check.build_settings(args)
    settings["cv"] = dict(SMALL_CV)
    settings["effects"] = {mode: list(values) for mode, values in SMALL_EFFECTS.items()}
    return settings


def _records(out_dir: Path) -> list[dict]:
    lines = check.shard_path(out_dir, 0, 1).read_text(encoding="utf-8").splitlines()
    records = sorted((json.loads(line) for line in lines), key=lambda record: record["key"])
    for record in records:
        record.pop("runtime_seconds")
    return records


@pytest.fixture(scope="module")
def runs(tmp_path_factory):
    outputs = []
    for name in ("first", "second"):
        out_dir = tmp_path_factory.mktemp(name)
        settings = _settings(out_dir)
        check.ensure_settings(out_dir, settings)
        assert check.run_shard(out_dir, settings, 0, 1) == 10  # 5 cells x 2 replicates
        outputs.append((out_dir, settings))
    return outputs


def test_units_are_reproducible(runs) -> None:
    (first, _), (second, _) = runs
    assert _records(first) == _records(second)


def test_arms_share_the_data_and_differ_in_the_floor(runs) -> None:
    out_dir, settings = runs[0]
    records = _records(out_dir)
    assert {r["mode"] for r in records} == {"anchor", "translation", "magnitude", "orientation", "shape"}
    for record in records:
        assert set(record["arms"]) == {"floor_1", "floor_2"}
        one, two = record["arms"]["floor_1"], record["arms"]["floor_2"]
        assert one["min_components"] == 1 and two["min_components"] == 2
        assert two["selected_lv"] >= 2 and two["angle_one_dimensional"] is False
        assert one["selected_lv"] >= 1 and one["angle_one_dimensional"] == (one["selected_lv"] == 1)
        for entry in (one, two):
            assert entry["error"] is None
            assert entry["component_selection"] == "cv"
            assert set(entry["p_values"]) == {"delta", "angle", "shape"}
            assert set(entry["null_summary"]) == {"delta", "angle", "shape"}
            assert 0.0 <= entry["cv_mean_auroc"] <= 1.0
        # Both arms were analysed on one built unit; the build is deterministic in the seeds.
        again = check.build_unit(record["mode"], record["effect"], record["replicate"], settings)
        assert again["seeds"] == record["seeds"]
        rebuilt = check.build_unit(record["mode"], record["effect"], record["replicate"], settings)
        for name in settings["layers"]:
            np.testing.assert_array_equal(again["blocks"][name].to_numpy(), rebuilt["blocks"][name].to_numpy())
    # Matched seeds: every cell of a replicate shares the generator and RRPP seeds.
    by_replicate: dict[int, set[str]] = {}
    for record in records:
        by_replicate.setdefault(record["replicate"], set()).add(json.dumps(record["seeds"], sort_keys=True))
    assert all(len(seeds) == 1 for seeds in by_replicate.values())


def test_resume_skips_recorded_units(runs) -> None:
    out_dir, settings = runs[0]
    before = _records(out_dir)
    assert check.run_shard(out_dir, settings, 0, 1) == 0
    assert _records(out_dir) == before


def test_merge_writes_the_decision(runs) -> None:
    out_dir, settings = runs[0]
    summary = check.merge(out_dir, settings)
    assert summary["units"] == 10 and summary["missing"] == []
    decision = json.loads((out_dir / "rank_floor_decision.json").read_text())
    assert decision["decision"] in {"adopt_floor", "diagnose_floor_excess"}
    assert decision["judged_arm"] == "floor_2"
    assert set(decision["controls"]) == {"anchor", "translation"}
    assert set(decision["verdicts"]["anchor"]) == {"delta", "angle", "shape"}
    assert decision["calibration"]["replicates"] == 2
    assert decision["construction"]["strengths"] == {
        "delta_methyl": 2.0, "delta_expr": 2.0, "delta_protein": 2.0, "p_dmp": 0.1
    }
    assert decision["failures"] == {"floor_1": 0, "floor_2": 0}
    assert len(decision["power"]) == 3 * 3  # three power cells x three statistics
    for name in ("pvalues", "replicates", "rates"):
        assert (out_dir / f"rank_floor_{name}.csv").exists()
    rates = (out_dir / "rank_floor_rates.csv").read_text(encoding="utf-8").splitlines()
    assert rates[0].split(",")[:3] == ["mode", "effect", "statistic"]
    assert len(rates) == 1 + 5 * 3


def test_changed_settings_are_refused(runs) -> None:
    out_dir, settings = runs[0]
    with pytest.raises(SystemExit, match="different settings"):
        check.ensure_settings(out_dir, {**settings, "permutations": 19})


def test_headroom_report_covers_the_pool_limited_cells(runs) -> None:
    _, settings = runs[0]
    report = check.headroom_report(settings)
    assert set(report) == {'["translation", 1.0]', '["orientation", 1.0]', '["shape", 1.0]'}
    assert all(entry["fits"] for entry in report.values())


# ── Decision rule ────────────────────────────────────────────────────────────


def _uniform(n: int, rejections: int) -> list[float]:
    return [0.01] * rejections + [0.5] * (n - rejections)


def _rows(rates: dict[tuple[str, float], dict[str, dict[str, int]]], n: int = 500) -> list[dict]:
    """``{(mode, effect): {arm: {statistic: rejections}}}`` → p-value rows."""

    rows = []
    for (mode, effect), per_arm in rates.items():
        for arm, per_statistic in per_arm.items():
            for statistic, rejections in per_statistic.items():
                for replicate, p_value in enumerate(_uniform(n, rejections)):
                    rows.append({"mode": mode, "effect": effect, "replicate": replicate, "arm": arm,
                                 "statistic": statistic, "p_value": p_value})
    return rows


HOLD = {"delta": 25, "angle": 30, "shape": 34}  # 0.068 <= 0.0695 at R = 500
EXCESS_ANGLE = {"delta": 25, "angle": 36, "shape": 34}  # 0.072 > 0.0695
INFLATED = {"delta": 200, "angle": 200, "shape": 200}
ANCHOR, TRANSLATION = ("anchor", 0.0), ("translation", 1.0)


def test_floor_holds_at_both_controls_adopts_the_floor() -> None:
    decision = check.rank_floor_decision(
        _rows({ANCHOR: {"floor_1": HOLD, "floor_2": HOLD}, TRANSLATION: {"floor_1": HOLD, "floor_2": HOLD}}),
        alpha=0.05, k=2.0,
    )
    assert decision["controls"]["anchor"]["floor_2"]["shape"]["bound"] == pytest.approx(0.0695, abs=1e-4)
    assert decision["verdicts"] == {"anchor": dict.fromkeys(HOLD, "holds"), "translation": dict.fromkeys(HOLD, "holds")}
    assert decision["exceeding"] == []
    assert decision["decision"] == "adopt_floor"


def test_anchor_angle_excess_is_diagnosed() -> None:
    decision = check.rank_floor_decision(
        _rows({ANCHOR: {"floor_1": HOLD, "floor_2": EXCESS_ANGLE}, TRANSLATION: {"floor_1": HOLD, "floor_2": HOLD}}),
        alpha=0.05, k=2.0,
    )
    assert decision["verdicts"]["anchor"]["angle"] == "exceeds"
    assert decision["exceeding"] == [{"mode": "anchor", "statistic": "angle", "rate": pytest.approx(0.072),
                                      "bound": pytest.approx(0.0695, abs=1e-4)}]
    assert decision["decision"] == "diagnose_floor_excess"


def test_translation_excess_is_diagnosed() -> None:
    decision = check.rank_floor_decision(
        _rows({ANCHOR: {"floor_1": HOLD, "floor_2": HOLD}, TRANSLATION: {"floor_1": HOLD, "floor_2": EXCESS_ANGLE}}),
        alpha=0.05, k=2.0,
    )
    assert decision["verdicts"]["anchor"] == dict.fromkeys(HOLD, "holds")
    assert decision["verdicts"]["translation"]["angle"] == "exceeds"
    assert decision["decision"] == "diagnose_floor_excess"


def test_floor_one_arm_does_not_enter_the_decision() -> None:
    held = check.rank_floor_decision(
        _rows({ANCHOR: {"floor_1": HOLD, "floor_2": HOLD}, TRANSLATION: {"floor_1": HOLD, "floor_2": HOLD}}),
        alpha=0.05, k=2.0,
    )
    inflated = check.rank_floor_decision(
        _rows({ANCHOR: {"floor_1": INFLATED, "floor_2": HOLD}, TRANSLATION: {"floor_1": INFLATED, "floor_2": HOLD}}),
        alpha=0.05, k=2.0,
    )
    assert inflated["controls"]["anchor"]["floor_1"]["angle"]["rate"] == pytest.approx(0.4)
    assert inflated["controls"]["anchor"]["floor_1"]["angle"]["verdict"] == "descriptive"
    assert "floor_1" not in inflated["verdicts"]["anchor"]
    assert inflated["decision"] == held["decision"] == "adopt_floor"


def test_missing_p_values_reduce_r_and_an_absent_control_is_unavailable() -> None:
    rows = _rows({ANCHOR: {"floor_1": HOLD, "floor_2": HOLD}})
    rows += [{"mode": "anchor", "effect": 0.0, "replicate": 0, "arm": "floor_2", "statistic": "delta", "p_value": None}]
    decision = check.rank_floor_decision(rows, alpha=0.05, k=2.0)
    assert decision["controls"]["anchor"]["floor_2"]["delta"]["replicates"] == 500
    assert set(decision["verdicts"]["translation"].values()) == {"unavailable"}
    assert decision["decision"] == "diagnose_floor_excess"


def test_calibration_fraction_is_read_from_the_unfloored_anchor_ranks() -> None:
    rows = _rows({ANCHOR: {"floor_1": HOLD, "floor_2": HOLD}, TRANSLATION: {"floor_1": HOLD, "floor_2": HOLD}})
    ranks = [{"mode": "anchor", "effect": 0.0, "replicate": i, "arm": "floor_1", "selected_lv": 1 if i < 3 else 2}
             for i in range(4)]
    ranks += [{"mode": "anchor", "effect": 0.0, "replicate": i, "arm": "floor_2", "selected_lv": 2} for i in range(4)]
    ranks += [{"mode": "shape", "effect": 1.0, "replicate": 0, "arm": "floor_1", "selected_lv": 1}]
    met = check.rank_floor_decision(rows, alpha=0.05, k=2.0, ranks=ranks)
    assert met["calibration"]["rank1_fraction"] == pytest.approx(0.75)
    assert met["calibration"]["rank_distribution"] == {"1": 3, "2": 1}
    assert met["calibration"]["met"] is True
    unmet = check.rank_floor_decision(rows, alpha=0.05, k=2.0, ranks=ranks, calibration_threshold=0.9)
    assert unmet["calibration"]["met"] is False
    assert unmet["decision"] == met["decision"] == "adopt_floor"


def test_rate_rows_report_both_arms_and_their_difference() -> None:
    rows = _rows(
        {ANCHOR: {"floor_1": HOLD, "floor_2": HOLD}, ("orientation", 1.0): {"floor_1": HOLD, "floor_2": INFLATED}}
    )
    ranks = [{"mode": "orientation", "effect": 1.0, "replicate": i, "arm": arm, "selected_lv": lv}
             for i, (arm, lv) in enumerate([("floor_1", 1), ("floor_1", 2), ("floor_2", 2), ("floor_2", 3)])]
    table = {(row["mode"], row["statistic"]): row for row in check.rate_rows(rows, alpha=0.05, k=2.0, ranks=ranks)}
    power = table[("orientation", "angle")]
    assert power["verdict"] == "descriptive" and power["bound"] is None
    assert power["rate_floor_1"] == pytest.approx(0.06) and power["rate_floor_2"] == pytest.approx(0.4)
    assert power["difference"] == pytest.approx(0.34)
    assert power["pooled_se"] == pytest.approx((0.06 * 0.94 / 500 + 0.4 * 0.6 / 500) ** 0.5)
    assert json.loads(power["rank_distribution_floor_1"]) == {"1": 1, "2": 1}
    assert power["rank1_fraction_floor_1"] == pytest.approx(0.5) and power["rank1_fraction_floor_2"] == 0.0
    anchor = table[("anchor", "angle")]
    assert anchor["verdict"] == "holds" and anchor["bound"] == pytest.approx(0.0695, abs=1e-4)
