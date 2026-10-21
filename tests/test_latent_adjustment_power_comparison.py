"""``scripts/latent_adjustment_power_comparison.py``: matched units, rates with ranks, and the decision rule."""

from __future__ import annotations

import csv
import json
import sys
from pathlib import Path

import numpy as np
import pytest


def _script():
    sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
    try:
        return __import__("latent_adjustment_power_comparison")
    finally:
        sys.path.pop(0)


power = _script()
SMALL_CV = {"cv1_splits": 2, "cv2_splits": 3, "n_repeats": 2, "max_components": 5, "random_state": 1203}


def _settings(out_dir: Path) -> dict:
    parser = power.build_parser()
    args = parser.parse_args(
        [
            "run", "--reps", "1", "--perms", "9", "--out-dir", str(out_dir),
            "--magnitude-effects", "0.1", "--orientation-effects", "1.0", "--shape-effects", "1.0",
        ]
    )
    args.calibration_threshold = 0.8
    settings = power.build_settings(args)
    settings["cv"] = dict(SMALL_CV)
    return settings


def test_grid_and_rule_are_in_the_settings(tmp_path: Path) -> None:
    settings = _settings(tmp_path)
    assert settings["columns"]["nuisance"] == {"s_z": 40.0, "s_c": 2.0}
    assert settings["columns"]["clean"] == {"s_z": 0.0, "s_c": 0.0}
    assert power.cells(settings) == [
        ("anchor", 0.0), ("translation", 1.0), ("magnitude", 0.1), ("orientation", 1.0), ("shape", 1.0)
    ]
    assert len(power.units(settings)) == 1 * 2 * 5
    assert settings["rule"]["gain_targets"] == [["orientation", "angle"]]
    assert settings["rule"]["gain_se_multiplier"] == settings["rule"]["loss_se_multiplier"] == 2.0
    assert "calibration" not in settings


def test_columns_share_the_data_and_differ_only_in_the_injection(tmp_path: Path) -> None:
    settings = _settings(tmp_path)
    nuisance = power.build_unit("nuisance", "orientation", 1.0, 0, settings)
    clean = power.build_unit("clean", "orientation", 1.0, 0, settings)
    assert nuisance["seeds"] == clean["seeds"]
    assert nuisance["metadata"].equals(clean["metadata"])
    for name in settings["layers"]:
        difference = nuisance["blocks"][name].to_numpy() - clean["blocks"][name].to_numpy()
        assert np.abs(difference).max() > 0
    # Cells share the replicate's seeds too (matched seeds across effects).
    anchor = power.build_unit("nuisance", "anchor", 0.0, 0, settings)
    assert anchor["seeds"] == nuisance["seeds"]
    assert anchor["metadata"]["z"].equals(nuisance["metadata"]["z"])


def test_arms_of_a_unit_differ_only_in_the_setting(tmp_path: Path, monkeypatch) -> None:
    import motco.analysis as analysis

    calls = []
    original = analysis.analyze_trajectories

    def recording(dataset, params):
        calls.append((dataset, params))
        return original(dataset, params)

    monkeypatch.setattr(analysis, "analyze_trajectories", recording)
    settings = _settings(tmp_path)
    record = power.evaluate_unit("nuisance", "orientation", 1.0, 0, settings)
    assert [params.covariate_latent_adjustment for _, params in calls] == [None, "stage", "group_stage"]
    datasets = {id(dataset) for dataset, _ in calls}
    assert len(datasets) == 1
    from dataclasses import replace

    stripped = {replace(params, covariate_latent_adjustment=None) for _, params in calls}
    assert len(stripped) == 1  # same CV settings, permutations and RRPP seed
    assert set(record["arms"]) == {"baseline", "stage", "group_stage"}
    assert record["arms"]["baseline"]["latent_adjustment"] is None


def test_run_and_merge_write_rates_with_rank_columns(tmp_path: Path) -> None:
    settings = _settings(tmp_path)
    power.ensure_settings(tmp_path, settings)
    assert power.run_shard(tmp_path, settings, 0, 1, limit=2) == 2
    assert power.run_shard(tmp_path, settings, 0, 1, limit=2) == 2  # the next two pending units
    summary = power.merge(tmp_path, settings)
    assert summary["units"] == 4 and summary["expected"] == 10
    with (tmp_path / "latent_adjustment_power_rates.csv").open() as handle:
        rows = list(csv.DictReader(handle))
    assert {"rank_distribution", "rank1_fraction"} <= set(rows[0])
    for row in rows:
        counts = dict(part.split(":") for part in row["rank_distribution"].split(";") if part)
        assert sum(int(v) for v in counts.values()) in {0, 1}
    decision = json.loads((tmp_path / "latent_adjustment_power_decision.json").read_text())
    assert set(decision["verdicts"]) == {"stage", "group_stage"}
    assert decision["verdicts"]["stage"] == "unavailable"  # the grid is incomplete


# ── Decision rule ────────────────────────────────────────────────────────────


EFFECTS = {"magnitude": [0.1, 0.25], "orientation": [0.5, 1.0], "shape": [0.5, 1.0]}
CONTROLS = {"anchor": 0.0, "translation": 1.0}
R = 500


def _rows(overrides: dict[tuple, float] | None = None, *, default: float = 0.3, control: float = 0.04,
          drop_arm: str | None = None) -> list[dict]:
    """R replicates per (column, mode, effect, arm, statistic) at the given rejection rates."""

    overrides = overrides or {}
    rows = []
    cells = [(mode, effect) for mode, effects in EFFECTS.items() for effect in effects] + list(CONTROLS.items())
    for column in ("nuisance", "clean"):
        for mode, effect in cells:
            for arm in ("baseline", "stage", "group_stage"):
                if arm == drop_arm:
                    continue
                for statistic in ("delta", "angle", "shape"):
                    rate = overrides.get(
                        (column, mode, effect, arm, statistic), control if mode in CONTROLS else default
                    )
                    rejections = round(rate * R)
                    for replicate in range(R):
                        rows.append({
                            "column": column, "mode": mode, "effect": effect, "replicate": replicate, "arm": arm,
                            "statistic": statistic, "p_value": 0.01 if replicate < rejections else 0.5,
                            "selected_lv": 1 if replicate < 100 else 2,
                        })
    return rows


def _decide(rows: list[dict]) -> dict:
    return power.latent_adjustment_power_decision(rows, alpha=0.05, k=2.0, effects=EFFECTS, controls=CONTROLS)


def test_gain_with_the_nuisance_axis_and_no_loss_without_it_qualifies() -> None:
    gain = {("nuisance", "orientation", 1.0, arm, "angle"): 0.45 for arm in ("stage", "group_stage")}
    gain[("nuisance", "orientation", 1.0, "baseline", "angle")] = 0.30
    decision = _decide(_rows(gain))
    assert decision["verdicts"] == {"stage": "candidate_qualifies", "group_stage": "candidate_qualifies"}
    read = decision["arms"]["stage"]["reads"]["gain"][0]
    assert read["difference"] == pytest.approx(0.15)
    assert read["threshold"] == pytest.approx(2.0 * np.sqrt(0.45 * 0.55 / R + 0.30 * 0.70 / R))
    assert "baseline" not in decision["verdicts"]


def test_no_gain_fails_the_candidate() -> None:
    decision = _decide(_rows())
    assert decision["verdicts"]["stage"] == "candidate_fails"
    assert list(decision["arms"]["stage"]["failing"]) == ["gain"]


def test_loss_without_the_nuisance_axis_fails_the_candidate() -> None:
    overrides = {("nuisance", "orientation", 1.0, arm, "angle"): 0.45 for arm in ("stage", "group_stage")}
    overrides[("nuisance", "orientation", 1.0, "baseline", "angle")] = 0.30
    overrides[("clean", "magnitude", 0.1, "baseline", "delta")] = 0.50
    overrides[("clean", "magnitude", 0.1, "group_stage", "delta")] = 0.40
    overrides[("clean", "magnitude", 0.1, "stage", "delta")] = 0.48
    decision = _decide(_rows(overrides))
    assert decision["verdicts"] == {"stage": "candidate_qualifies", "group_stage": "candidate_fails"}
    failing = decision["arms"]["group_stage"]["failing"]
    assert list(failing) == ["no_loss"]
    assert failing["no_loss"][0]["mode"] == "magnitude" and failing["no_loss"][0]["statistic"] == "delta"


def test_control_excess_fails_the_candidate() -> None:
    overrides = {("nuisance", "orientation", 1.0, arm, "angle"): 0.45 for arm in ("stage", "group_stage")}
    overrides[("nuisance", "orientation", 1.0, "baseline", "angle")] = 0.30
    overrides[("clean", "translation", 1.0, "stage", "angle")] = 0.09
    decision = _decide(_rows(overrides))
    assert decision["verdicts"] == {"stage": "candidate_fails", "group_stage": "candidate_qualifies"}
    control = decision["arms"]["stage"]["failing"]["controls"][0]
    assert control["bound"] == pytest.approx(0.0695, abs=1e-4)


def test_baseline_control_excess_does_not_judge_the_baseline() -> None:
    overrides = {("nuisance", "orientation", 1.0, arm, "angle"): 0.45 for arm in ("stage", "group_stage")}
    overrides[("nuisance", "orientation", 1.0, "baseline", "angle")] = 0.30
    overrides[("nuisance", "anchor", 0.0, "baseline", "shape")] = 0.30
    decision = _decide(_rows(overrides))
    assert decision["verdicts"] == {"stage": "candidate_qualifies", "group_stage": "candidate_qualifies"}


def test_missing_arm_is_unavailable() -> None:
    decision = _decide(_rows(drop_arm="group_stage"))
    assert decision["verdicts"]["group_stage"] == "unavailable"


def test_rate_table_carries_rank_distribution_and_rank1_fraction() -> None:
    table = power.rate_table(_rows(), alpha=0.05)
    row = next(r for r in table if r["column"] == "clean" and r["mode"] == "shape" and r["arm"] == "stage")
    assert row["rank_distribution"] == "1:100;2:400"
    assert row["rank1_fraction"] == pytest.approx(0.2)
    assert row["replicates"] == R and row["rate"] == pytest.approx(0.3)
    assert row["se"] == pytest.approx(np.sqrt(0.3 * 0.7 / R))
