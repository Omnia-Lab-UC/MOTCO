"""Tests for ``scripts/orientation_measurement_probe.py``: seeds, resume, rate table, advancement rule."""

from __future__ import annotations

import importlib.util
import json
import math
import sys
from pathlib import Path

import pytest

_SCRIPT = Path(__file__).resolve().parents[1] / "scripts" / "orientation_measurement_probe.py"
_spec = importlib.util.spec_from_file_location("orientation_measurement_probe", _SCRIPT)
assert _spec is not None and _spec.loader is not None
probe = importlib.util.module_from_spec(_spec)
sys.modules[_spec.name] = probe  # dataclasses resolve their module through sys.modules
_spec.loader.exec_module(probe)

SMALL = ["--reps", "3", "--perms", "19", "--designs", "phase6"]


# ── Seeds and arms on a real unit ────────────────────────────────────────────


def test_every_arm_of_a_unit_shares_its_seeds() -> None:
    settings = probe.build_settings(probe.build_parser().parse_args(["--out-dir", "x", *SMALL]))
    records = probe.evaluate_unit(probe.Unit("phase6", "orientation_1.0", 0), settings)
    assert [r["arm"] for r in records] == list(probe.ARMS)
    assert all(r["status"] == "completed" for r in records)
    assert {(r["generator_seed"], r["rrpp_seed"]) for r in records} == {
        (probe.generator_seed(900, "phase6", 0), probe.rrpp_seed(900, "phase6", "orientation_1.0", 0))
    }
    # Generator seeds are matched across cells; RRPP seeds are per unit.
    assert probe.generator_seed(900, "phase6", 0) == probe.generator_seed(900, "phase6", 0)
    assert probe.rrpp_seed(900, "phase6", "anchor", 0) != probe.rrpp_seed(900, "phase6", "orientation_1.0", 0)
    by_arm = {r["arm"]: r for r in records}
    # Estimators share the space fit: everything but angle agrees.
    for space in probe.SPACES:
        pc1, net = by_arm[f"{space}__pc1"], by_arm[f"{space}__net_displacement"]
        assert pc1["retained_dimensionality"] == net["retained_dimensionality"]
        assert pc1["pair_statistics"]["delta"] == net["pair_statistics"]["delta"]
        assert pc1["pair_statistics"]["shape"] == net["pair_statistics"]["shape"]
    assert by_arm["pca__pc1"]["retained_dimensionality"] == 4
    assert by_arm["pls_forced__pc1"]["retained_dimensionality"] == 4
    assert by_arm["concat__pc1"]["realized_angle"]["latent"] is None
    assert by_arm["pca__pc1"]["realized_angle"]["latent"] is not None


# ── Resume, settings header, shards ──────────────────────────────────────────


def _fake_unit(unit, settings):
    """Deterministic stand-in for a unit: records derived from its seeds only."""

    gen = probe.generator_seed(settings["base_seed"], unit.design, unit.replicate)
    perm = probe.rrpp_seed(settings["base_seed"], unit.design, unit.cell, unit.replicate)
    records = []
    for index, arm in enumerate(probe.arm_list(settings)):
        space, estimator = arm.split("__")
        p = ((perm + 31 * index) % 100) / 100
        records.append(
            {
                "design": unit.design,
                "cell": unit.cell,
                "replicate": unit.replicate,
                "arm": arm,
                "space": space,
                "estimator": estimator,
                "generator_seed": gen,
                "rrpp_seed": perm,
                "status": "completed",
                "retained_dimensionality": 3,
                "p_values": {"delta": p, "angle": 1 - p, "shape": p / 2},
                "realized_angle": {"population": 90.0, "observed_standardized": 80.0, "latent": 30.0},
            }
        )
    return records


def _outputs(out_dir: Path) -> dict[str, str]:
    names = ("merged.jsonl", "rate_table.csv", "checkpoint_table.csv", "advancement_decision.json", "summary.md")
    return {name: (out_dir / name).read_text(encoding="utf-8") for name in names}


def test_interrupted_then_resumed_run_equals_uninterrupted(tmp_path: Path, monkeypatch) -> None:
    monkeypatch.setattr(probe, "evaluate_unit", _fake_unit)
    straight = tmp_path / "straight"
    assert probe.main(["--out-dir", str(straight), *SMALL]) == 0

    resumed = tmp_path / "resumed"
    assert probe.main(["--out-dir", str(resumed), *SMALL, "--limit", "4"]) == 0
    shard = probe.shard_path(resumed, 0, 1)
    # Simulate an interruption mid-unit: a partial unit and a torn line.
    partial = _fake_unit(probe.Unit("phase6", "shape_1.0", 2), probe.load_settings(resumed))[:3]
    with shard.open("a", encoding="utf-8") as handle:
        handle.write("".join(json.dumps(r) + "\n" for r in partial) + '{"design": "pha')
    calls = []
    monkeypatch.setattr(probe, "evaluate_unit", lambda unit, settings: calls.append(unit) or _fake_unit(unit, settings))
    assert probe.main(["--out-dir", str(resumed), *SMALL]) == 0
    assert len(calls) == 15 - 4  # only the missing units ran
    assert _outputs(resumed) == _outputs(straight)


def test_resume_with_changed_settings_is_refused(tmp_path: Path, monkeypatch) -> None:
    monkeypatch.setattr(probe, "evaluate_unit", _fake_unit)
    assert probe.main(["--out-dir", str(tmp_path), *SMALL, "--limit", "1"]) == 0
    with pytest.raises(SystemExit, match=r"differing: permutations, replicates"):
        probe.main(["--out-dir", str(tmp_path), "--reps", "4", "--perms", "99", "--designs", "phase6"])


def test_shards_partition_the_units_and_merge(tmp_path: Path, monkeypatch) -> None:
    monkeypatch.setattr(probe, "evaluate_unit", _fake_unit)
    for index in range(3):
        assert probe.main(["--out-dir", str(tmp_path), *SMALL, "--shard", f"{index}/3"]) == 0
    assert not (tmp_path / "merged.jsonl").exists()
    assert probe.main(["--out-dir", str(tmp_path), "--merge"]) == 0
    merged = [json.loads(line) for line in (tmp_path / "merged.jsonl").read_text().splitlines()]
    assert len(merged) == 15 * len(probe.ARMS)
    assert json.loads((tmp_path / "advancement_decision.json").read_text())["complete"] is True


def test_smoke_summary_has_every_section(tmp_path: Path, monkeypatch) -> None:
    monkeypatch.setattr(probe, "evaluate_unit", _fake_unit)
    assert probe.main(["--out-dir", str(tmp_path), *SMALL]) == 0
    summary = (tmp_path / "summary.md").read_text(encoding="utf-8")
    for heading in (
        "## Settings",
        "## Phase 5 deviation",
        "## Rejection rates",
        "### phase6",
        "## Checkpoint table",
        "## Advancement decision",
    ):
        assert heading in summary
    assert "joint magnitude construction" in summary


# ── Rate table and advancement rule ──────────────────────────────────────────


def test_rate_table_counts_strict_rejections() -> None:
    records = [
        {"design": "d", "cell": "c", "arm": "a", "status": "completed", "p_values": {"angle": p}}
        for p in (0.01, 0.05, 0.2, None)
    ]
    (row,) = probe.rate_table(records, alpha=0.05)
    assert (row["count"], row["rejections"]) == (3, 1)
    assert row["se"] == pytest.approx(math.sqrt((1 / 3) * (2 / 3) / 3))


R = 100


def _row(design, cell, arm, statistic, rate):
    return {
        "design": design,
        "cell": cell,
        "arm": arm,
        "statistic": statistic,
        "rate": rate,
        "se": math.sqrt(rate * (1 - rate) / R),
        "count": R,
        "rejections": round(rate * R),
    }


def _table(arm_rates: dict[str, dict[tuple[str, str, str], float]]) -> list[dict]:
    """A complete rate table: every arm at the reference's rates unless overridden."""

    base = {}
    for design in ("phase6", "phase5"):
        for statistic in probe.STATISTICS:
            base[(design, "anchor", statistic)] = 0.04
        base[(design, "orientation_1.0", "angle")] = 0.30
        base[(design, "magnitude_1.0", "delta")] = 0.80
        base[(design, "magnitude_1.0", "angle")] = 0.10
        base[(design, "shape_1.0", "shape")] = 0.50
        base[(design, "shape_1.0", "angle")] = 0.10
    rows = []
    for arm in ("pls_cv__pc1", *arm_rates):
        rates = {**base, **arm_rates.get(arm, {})}
        rows += [_row(design, cell, arm, statistic, rate) for (design, cell, statistic), rate in rates.items()]
    return rows


GAIN = {("phase6", "orientation_1.0", "angle"): 0.60}


def test_an_arm_that_gains_and_keeps_everything_advances() -> None:
    decision = probe.advancement_decision(_table({"pca__net_displacement": GAIN}), alpha=0.05)
    assert decision["advancing"] == ["pca__net_displacement"]
    assert decision["arms"]["pls_cv__pc1"]["role"] == "reference"


def test_anchor_inflation_blocks_advancement() -> None:
    rates = {**GAIN, ("phase5", "anchor", "angle"): 0.10}
    decision = probe.advancement_decision(_table({"pca__net_displacement": rates}), alpha=0.05)
    entry = decision["arms"]["pca__net_displacement"]
    assert entry["criteria"]["anchor"]["pass"] is False
    assert entry["criteria"]["gain"]["pass"] is True
    assert entry["advances"] is False
    bound = probe.anchor_bound(0.05, R)
    assert bound == pytest.approx(0.05 + 2 * math.sqrt(0.05 * 0.95 / R))
    failing = [c for c in entry["criteria"]["anchor"]["checks"] if not c["pass"]]
    assert [(c["design"], c["statistic"]) for c in failing] == [("phase5", "angle")]
    assert failing[0]["threshold"] == pytest.approx(bound)


def test_specificity_failure_blocks_advancement() -> None:
    rates = {**GAIN, ("phase6", "magnitude_1.0", "angle"): 0.40}
    decision = probe.advancement_decision(_table({"pls_within_group__pc1": rates}), alpha=0.05)
    entry = decision["arms"]["pls_within_group__pc1"]
    assert entry["criteria"]["specificity"]["pass"] is False
    assert entry["advances"] is False
    assert "specificity" in entry["reason"]


def test_no_arm_advancing_records_every_reason() -> None:
    decision = probe.advancement_decision(
        _table({"pca__pc1": {}, "pls_forced__pc1": {}, "concat__net_displacement": GAIN}), alpha=0.05
    )
    assert decision["advancing"] == []
    for arm in ("pca__pc1", "pls_forced__pc1"):
        assert decision["arms"][arm]["reason"] == "failed: gain"
        assert set(decision["arms"][arm]["criteria"]) == {"anchor", "gain", "no_loss", "protected", "specificity"}


def test_concat_never_advances() -> None:
    decision = probe.advancement_decision(_table({"concat__net_displacement": GAIN}), alpha=0.05)
    entry = decision["arms"]["concat__net_displacement"]
    assert all(criterion["pass"] for criterion in entry["criteria"].values())
    assert entry["role"] == "ceiling_reference"
    assert entry["advances"] is False
    assert decision["advancing"] == []


def test_two_pooled_se_thresholds() -> None:
    reference = 0.30
    margin = 2 * math.sqrt(2 * reference * (1 - reference) / R)  # equal rates: pooled SE of the reference

    def gain_passes(rate: float) -> bool:
        table = _table({"pca__pc1": {("phase6", "orientation_1.0", "angle"): rate}})
        check = probe.advancement_decision(table, alpha=0.05)["arms"]["pca__pc1"]["criteria"]["gain"]["checks"][0]
        return check["pass"]

    # The margin is computed from both arms' SEs, so probe either side of it.
    assert not gain_passes(reference + 0.9 * margin)
    assert gain_passes(reference + 1.3 * margin)

    def loss_passes(rate: float) -> bool:
        table = _table({"pca__pc1": {("phase5", "orientation_1.0", "angle"): rate}})
        return probe.advancement_decision(table, alpha=0.05)["arms"]["pca__pc1"]["criteria"]["no_loss"]["pass"]

    assert loss_passes(reference - 0.9 * margin)
    assert not loss_passes(reference - 1.3 * margin)


def test_missing_rates_fail_their_criterion() -> None:
    table = [r for r in _table({"pca__pc1": GAIN}) if not (r["arm"] == "pca__pc1" and r["design"] == "phase5")]
    entry = probe.advancement_decision(table, alpha=0.05)["arms"]["pca__pc1"]
    assert entry["criteria"]["no_loss"]["pass"] is False
    assert entry["criteria"]["no_loss"]["checks"][0]["missing"] is True
