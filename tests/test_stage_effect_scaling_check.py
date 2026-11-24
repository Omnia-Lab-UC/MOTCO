"""``scripts/stage_effect_scaling_check.py``: designs, seeds, calibration, the unit, merge, reproduction and decision.

Units run at fast CV knobs and few permutations on the smallest design; the
decision runs on synthetic tables over the real settings. No real data are read.
"""

from __future__ import annotations

import json
import math
import os
import shutil
import subprocess
import sys
from pathlib import Path
from unittest import mock

import pytest


def _script(name: str):
    sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
    try:
        return __import__(name)
    finally:
        sys.path.pop(0)


ses = _script("stage_effect_scaling_check")
ubc = ses.ubc
lsc = ses.lsc
ssc = ses.ssc

REPO = Path(__file__).resolve().parents[1]
FAST_CV = {"cv1_splits": 3, "cv2_splits": 4, "n_repeats": 1, "max_components": 6, "random_state": 1203}
PERMS = 39
NEW_SATURATED = ["n320_delta1", "n1280_delta0.5", "n800_delta1", "n800_delta0.25"]


def _steps(medians: dict[float, float | None]) -> list[dict]:
    return [{"delta": delta, "replicates": 40, "completed": 40, "auroc_median": auroc, "auroc_q25": auroc,
             "auroc_q75": auroc, "rank_median": 5.0, "rank_q10": 4.0, "rank_q90": 7.0}
            for delta, auroc in medians.items()]


def _write_calibration(out_dir: Path, medians: dict[float, float | None] | None = None, **extra) -> Path:
    settings = ses.calibration_settings()
    steps = _steps(medians or {0.03: 0.55, 0.04: 0.59, 0.05: 0.63, 0.06: 0.66, 0.07: 0.70, 0.08: 0.74, 0.1: 0.8})
    payload = {"settings": settings, "hash": lsc.settings_hash(settings), "steps": steps,
               "adopted": ses.adopt_calibration(steps), "failures": 0, "core_hours": 1.0, **extra}
    path = Path(out_dir) / ses.CALIBRATION_FILE
    path.parent.mkdir(parents=True, exist_ok=True)
    lsc._write_json(path, payload)
    return path


@pytest.fixture(scope="module")
def committed():
    return ssc.load_parent(ses.PARENT_DIR)["settings"]


@pytest.fixture(scope="module")
def saturated():
    return ses.saturated_settings()


@pytest.fixture()
def band(tmp_path):
    _write_calibration(tmp_path)
    return ses.band_settings(tmp_path)


def _seeds(settings, *, rrpp: bool = False) -> set[int]:
    return {ubc.unit_fields(unit, settings)["rrpp_seed" if rrpp else "generator_seed"]
            for unit in ubc.profile_units(settings)}


# ── Designs, unit counts and seeds ───────────────────────────────────────────


def test_unit_counts_per_profile(saturated, band) -> None:
    assert len(ses.profile_units(saturated)) == 16_000
    assert len(ses.profile_units(band)) == 8_400
    sat_signal = ubc.signal_cells(ses.restricted(saturated["unit_settings"], NEW_SATURATED))
    band_signal = ubc.signal_cells(band["unit_settings"])
    assert len(sat_signal) + len(band_signal) == ses.MAP_CELLS_EXPECTED == 71
    for settings in (saturated, band):
        units = settings["unit_settings"]
        assert units["permutations"] == 999 and units["replicates"] == {"control": 500, "other": 200}
        assert units["recorded_spaces"] == ["pca_matched"] and units["units"]["base_seed"] == 2043
        assert len({ubc.unit_key(unit) for unit in ses.profile_units(settings)}) == len(ses.profile_units(settings))


def test_unit_settings_are_the_parents_with_designs_and_cells_replaced() -> None:
    mine = ses.unit_settings({"n200_delta0.5": (200, 0.5)}, ubc.STAGE_SPAN_CELLS, replicates=200,
                             control_replicates=500, permutations=999)
    assert mine == ubc.stage_span_settings(designs=["n200_delta0.5"])


def test_the_rerun_design_is_the_parents(saturated, committed) -> None:
    units = saturated["unit_settings"]
    assert units["designs"][0] == ses.RERUN_DESIGN == "n200_delta0.5"
    assert units["units"]["designs"][ses.RERUN_DESIGN] == committed["units"]["designs"][ses.RERUN_DESIGN]
    assert units["units"]["cells"] == committed["units"]["cells"]
    for key in ("base_seed", "pls_cv_params", "min_components", "orientation_estimator"):
        assert units["units"][key] == committed["units"][key]
    for key in ("permutations", "replicates", "rules", "arms", "recorded_spaces", "population_null_tolerance"):
        assert units[key] == committed[key]
    rerun = [unit for unit in ses.profile_units(saturated) if unit.design == ses.RERUN_DESIGN]
    parent = [unit for unit in ubc.profile_units(committed) if unit.design == ses.RERUN_DESIGN]
    assert rerun == parent and len(rerun) == 3_200
    for unit in (rerun[0], rerun[777], rerun[1500], rerun[-1]):
        assert ubc.unit_fields(unit, units) == ubc.unit_fields(unit, committed)
        assert ubc.stage_span_generator_params(unit, units) == ubc.stage_span_generator_params(unit, committed)


def test_new_designs_differ_from_their_match_only_in_stage_effect_and_size(saturated, committed) -> None:
    mine = saturated["unit_settings"]["units"]["designs"]
    theirs = committed["units"]["designs"]
    allowed = {"n_samples", *ubc.STAGE_EFFECT_KEYS}
    assert saturated["pairs"] == [list(pair) for pair in ses.SATURATED_PAIRS]
    for weaker, match in ses.SATURATED_PAIRS:
        new, old = mine[weaker]["generator"], theirs[match]["generator"]
        assert set(new) == set(old)
        assert {key for key in new if new[key] != old[key]} == allowed
        assert mine[weaker]["layers"] is None and "span_layers" not in new
        assert len({new[key] for key in ubc.STAGE_EFFECT_KEYS}) == 1
        assert ses.design_delta(saturated["unit_settings"], weaker) < ses.design_delta(committed, match)
        assert math.isclose(ses.delta_sqrt_n(saturated["unit_settings"], weaker), ses.delta_sqrt_n(committed, match))


def test_a_new_design_measures_the_three_blocks(saturated) -> None:
    units = saturated["unit_settings"]
    unit = lsc.Unit("n320_delta1", "anchor", 0)
    dataset, params = ubc.unit_inputs(unit, units, ubc.unit_fields(unit, units))
    assert sum(getattr(dataset, layer).shape[1] for layer in ("methylation", "expression", "proteomics")) == 658
    assert dataset.metadata.shape[0] == 320
    assert "layers" not in params.integration_params
    assert dataset.truth["span_layers"] == ["methylation", "expression", "proteomics"]


def test_new_datasets_are_fresh_and_cells_share_their_null_dataset(saturated, committed) -> None:
    new = ses.restricted(saturated["unit_settings"], NEW_SATURATED)
    parent_seeds = _seeds(committed)
    assert not _seeds(new) & parent_seeds
    assert not _seeds(new, rrpp=True) & _seeds(committed, rrpp=True)
    assert _seeds(ses.restricted(saturated["unit_settings"], [ses.RERUN_DESIGN])) <= parent_seeds
    units = saturated["unit_settings"]
    for design in NEW_SATURATED:
        seeds = {ubc.unit_fields(lsc.Unit(design, cell, 3), units)["generator_seed"]
                 for cell in ubc._cell_names(units)}
        assert len(seeds) == 1


# ── Calibration ──────────────────────────────────────────────────────────────


def test_calibration_adopts_the_value_nearest_the_target() -> None:
    adopted = ses.adopt_calibration(_steps({0.03: 0.55, 0.05: 0.62, 0.06: 0.64, 0.1: 0.8, 0.2: None}))
    assert (adopted["delta"], adopted["auroc_median"], adopted["in_band"]) == (0.06, 0.64, True)
    assert adopted["distance"] == pytest.approx(0.004)


def test_a_tie_takes_the_smaller_delta() -> None:
    adopted = ses.adopt_calibration(_steps({0.07: 0.646, 0.05: 0.626, 0.1: 0.9}))
    assert adopted["delta"] == 0.05


def test_a_nearest_value_outside_the_band_is_refused_with_the_ladder(tmp_path) -> None:
    adopted = ses.adopt_calibration(_steps({0.03: 0.55, 0.04: 0.58, 0.1: 0.72}))
    assert (adopted["delta"], adopted["in_band"]) == (0.04, False)
    _write_calibration(tmp_path, {0.03: 0.55, 0.04: 0.58, 0.1: 0.72})
    with pytest.raises(SystemExit, match=r"outside the band 0.6–0.67.*δ = 0.03: median AUROC 0.550; "
                                         r"δ = 0.04: median AUROC 0.580; δ = 0.1: median AUROC 0.720"):
        ses.band_settings(tmp_path)


def test_calibration_settings_are_the_comparisons_target_band_and_stream() -> None:
    settings = ses.calibration_settings()
    assert settings["ladder"] == [0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.1]
    assert (settings["replicates"], settings["n_samples"]) == (40, 320)
    assert (settings["target"], settings["band"]) == (lsc.CALIBRATION_TARGET, list(lsc.CALIBRATION_BAND))
    assert (settings["target"], settings["band"]) == (0.636, [0.60, 0.67])
    assert settings["anchor"] == list(ubc.STAGE_SPAN_CELLS[0])


def test_a_calibration_unit_fits_production_on_three_blocks_without_permutations() -> None:
    import motco.simulations.evaluation as evaluation

    settings = {**ses.calibration_settings(), "pls_cv_params": dict(FAST_CV)}
    real = evaluation.integrate_semisynthetic_dataset
    seen = []

    def spy(dataset, params, *args, **kwargs):
        seen.append((dataset, params))
        return real(dataset, params, *args, **kwargs)

    with mock.patch.object(evaluation, "integrate_semisynthetic_dataset", spy):
        entry = ses.calibration_unit(0.05, 3, settings)
    assert entry["status"] == "completed", entry.get("error_message")
    assert entry["generator_seed"] == lsc.calibration_seed(2043, 3)
    assert 0.0 <= entry["auroc"] <= 1.0 and entry["rank"] >= 2
    (dataset, params), = seen
    assert params.integration_method == "pls" and "layers" not in params.integration_params
    assert params.permutations == 0 and params.min_components == 2
    assert params.orientation_estimator == "net_displacement"
    assert dataset.metadata.shape[0] == 320
    assert all(deltas == [0.05, 0.05, 0.05] for deltas in dataset.truth["deltas"].values())


def test_calibrate_writes_both_files_and_refuses_to_overwrite(tmp_path, monkeypatch) -> None:
    auroc = {0.03: 0.55, 0.04: 0.6, 0.05: 0.64, 0.06: 0.7, 0.07: 0.75, 0.08: 0.8, 0.1: 0.9}

    def fake(delta, replicate, settings):
        return {"delta": delta, "replicate": replicate, "generator_seed": replicate, "status": "completed",
                "auroc": auroc[delta] + (0.001 if replicate % 2 else -0.001), "rank": 3, "seconds": 1.0}

    monkeypatch.setattr(ses, "calibration_unit", fake)
    payload = ses.calibrate(tmp_path, ses.calibration_settings(replicates=4))
    assert payload["adopted"]["delta"] == 0.05 and payload["adopted"]["in_band"]
    assert len(payload["steps"]) == 7 and payload["steps"][0]["rank_q10"] == 3
    assert ses.load_calibration(tmp_path / ses.CALIBRATION_FILE)["adopted"] == payload["adopted"]
    rows = lsc._read_csv(tmp_path / "calibration" / "calibration_units.csv")
    assert len(rows) == 28 and list(rows[0]) == list(ses.CALIBRATION_UNIT_COLUMNS)
    with pytest.raises(SystemExit, match=r"calibration.json exists"):
        ses.calibrate(tmp_path, ses.calibration_settings(replicates=4))


# ── The band profile ─────────────────────────────────────────────────────────


def test_band_designs_cells_and_matched_pair(band) -> None:
    units = band["unit_settings"]
    assert units["designs"] == ["band_n320", "band_n1280_half", "band_n1280"]
    assert [ses.design_delta(units, d) for d in units["designs"]] == [0.05, 0.025, 0.05]
    assert [ses.design_size(units, d) for d in units["designs"]] == [320, 1280, 1280]
    assert math.isclose(ses.delta_sqrt_n(units, "band_n320"), ses.delta_sqrt_n(units, "band_n1280_half"))
    assert band["pairs"] == [["band_n1280_half", "band_n320"]]
    cells = units["units"]["cells"]
    assert len(cells) == 11 and all(cell[5] == "planted" for cell in cells)
    assert all(cell[3] in (None, 0.0) for cell in cells)
    assert {(cell[1], cell[2]) for cell in cells} == {
        ("none", 0.0), ("translation", 1.0), ("magnitude", 0.5), ("magnitude", 1.0), ("magnitude", 2.0),
        ("orientation", 0.25), ("orientation", 0.5), ("orientation", 1.0),
        ("shape", 0.25), ("shape", 0.5), ("shape", 1.0)}
    assert band["calibration"]["delta"] == 0.05 and band["calibration"]["file"] == "calibration/calibration.json"


def test_band_seeds_are_fresh(band, saturated, committed) -> None:
    ours = _seeds(band["unit_settings"])
    assert not ours & _seeds(saturated["unit_settings"])
    assert not ours & _seeds(committed)
    assert not ours & {lsc.calibration_seed(2043, r) for r in range(40)}


def test_a_missing_or_changed_calibration_is_refused(tmp_path, monkeypatch) -> None:
    with pytest.raises(SystemExit, match=r"calibration.json does not exist; run `calibrate` first"):
        ses.band_settings(tmp_path)
    _write_calibration(tmp_path)
    settings = ses.band_settings(tmp_path, replicates=1, control_replicates=1, permutations=9)
    monkeypatch.setattr(ses, "evaluate_unit", lambda unit, units: _fake_record(unit, units))
    ses.run_shard(tmp_path, settings, 0, 1, limit=2)
    _write_calibration(tmp_path, core_hours=2.0)
    with pytest.raises(SystemExit, match=r"calibration.json has changed since the band settings were written"):
        ses.run_shard(tmp_path, settings, 0, 1)
    with pytest.raises(SystemExit, match=r"calibration.json has changed"):
        ses.merge(tmp_path, "band")
    (tmp_path / ses.CALIBRATION_FILE).unlink()
    with pytest.raises(SystemExit, match=r"calibration.json does not exist; the band profile is bound"):
        ses.merge(tmp_path, "band")
    _write_calibration(tmp_path)
    corrupt = json.loads((tmp_path / ses.CALIBRATION_FILE).read_text())
    corrupt["settings"]["target"] = 0.7
    (tmp_path / ses.CALIBRATION_FILE).write_text(json.dumps(corrupt))
    with pytest.raises(SystemExit, match=r"is corrupt"):
        ses.band_settings(tmp_path)


# ── Resume, settings and parent binding (fast fake units) ────────────────────


def _fake_record(unit, units, *, rejections=None, status="completed", rank=2):
    fields = ubc.unit_fields(unit, units)
    base = {**fields, "settings_hash": lsc.settings_hash(units), "status": status,
            "runtime_seconds": {"unit": 1.0, "fit": 0.5}}
    if status != "completed":
        return {**base, "error_type": "RuntimeError", "error_message": "patched"}
    rejections = rejections or {}
    arms = {arm: {"p_value": 0.01 if rejections.get(arm) else 0.5, "rejects": bool(rejections.get(arm, False))}
            for arm in ubc.ARMS}
    return {
        **base,
        "blocks": ["methylation", "expression", "proteomics"],
        "n_features": 658,
        "n_population_null": 100,
        "production": {"rank": rank, "cv_mean_auroc": 0.9,
                       "span_capture": {"group_difference": 0.5, "pooled_stage": 0.9, "random_share": 0.01}},
        "scan": {"hits": 0, "hits_population_null": 0},
        "fdr": {"discoveries": 0, "discoveries_population_null": 0},
        "spaces": {"pca_matched": {"dimension": rank}},
        "arms": arms,
    }


def _small_saturated(**overrides):
    overrides.setdefault("designs", [ses.RERUN_DESIGN, "n320_delta1"])
    return ses.saturated_settings(replicates=2, control_replicates=2, permutations=PERMS, **overrides)


@pytest.fixture()
def fake_units(monkeypatch):
    monkeypatch.setattr(ses, "evaluate_unit", lambda unit, units: _fake_record(unit, units))


def test_resume_skips_complete_units(tmp_path, monkeypatch) -> None:
    settings = _small_saturated()
    calls = []

    def fake(unit, units):
        calls.append(ubc.unit_key(unit))
        return _fake_record(unit, units)

    monkeypatch.setattr(ses, "evaluate_unit", fake)
    assert ses.run_shard(tmp_path, settings, 0, 1, limit=5) == 5
    assert ses.run_shard(tmp_path, settings, 0, 1) == 52 - 5
    assert len(calls) == 52 and len(set(calls)) == 52
    assert ses.run_shard(tmp_path, settings, 0, 1) == 0
    summary = ses.merge(tmp_path, "saturated")
    assert summary["complete_units"] == 52 and summary["missing_units"] == 0
    assert len((tmp_path / "saturated" / "merged.jsonl").read_text().splitlines()) == 52


def test_shards_are_round_robin(tmp_path, fake_units) -> None:
    settings = _small_saturated()
    assert [ses.run_shard(tmp_path, settings, index, 3) for index in range(3)] == [18, 17, 17]
    assert ses.merge(tmp_path, "saturated")["complete_units"] == 52


def test_a_resume_under_other_settings_is_refused_naming_the_keys(tmp_path, fake_units) -> None:
    ses.run_shard(tmp_path, _small_saturated(), 0, 1, limit=1)
    with pytest.raises(SystemExit, match=r"differing: unit_settings.permutations"):
        ses.run_shard(tmp_path, ses.saturated_settings(designs=[ses.RERUN_DESIGN, "n320_delta1"], replicates=2,
                                                       control_replicates=2, permutations=19), 0, 1)
    with pytest.raises(SystemExit, match=r"differing: unit_settings.replicates.other"):
        ses.run_shard(tmp_path, ses.saturated_settings(designs=[ses.RERUN_DESIGN, "n320_delta1"], replicates=3,
                                                       control_replicates=2, permutations=PERMS), 0, 1)


def test_merge_refuses_records_written_under_other_settings(tmp_path, monkeypatch) -> None:
    monkeypatch.setattr(ses, "evaluate_unit",
                        lambda unit, units: {**_fake_record(unit, units), "settings_hash": "other"})
    ses.run_shard(tmp_path, _small_saturated(), 0, 1, limit=2)
    with pytest.raises(SystemExit, match=r"written under other settings"):
        ses.merge(tmp_path, "saturated")


def _parent_copy(tmp_path: Path) -> Path:
    target = tmp_path / "parent" / "stage_span"
    target.mkdir(parents=True)
    for name in ses.PARENT_FILE_NAMES:
        shutil.copy(ses.PARENT_DIR / name, target / name)
    return target


@pytest.mark.parametrize("name", ["rate_table.csv", "cell_table.csv"])
def test_a_changed_parent_table_is_refused_by_name(tmp_path, fake_units, name) -> None:
    parent = _parent_copy(tmp_path)
    settings = _small_saturated(parent_dir=parent)
    assert set(settings["parent_files"]) == {str(parent / n) for n in ses.PARENT_FILE_NAMES}
    out = tmp_path / "run"
    ses.run_shard(out, settings, 0, 1, limit=3)
    with (parent / name).open("a") as handle:
        handle.write("\n")
    with pytest.raises(SystemExit, match=rf"stage_span/{name} has changed since the settings were written"):
        ses.run_shard(out, settings, 0, 1)
    with pytest.raises(SystemExit, match=rf"{name} has changed"):
        ses.merge(out, "saturated")


def test_a_changed_parent_settings_file_is_refused(tmp_path, fake_units) -> None:
    parent = _parent_copy(tmp_path)
    settings = _small_saturated(parent_dir=parent)
    stored = json.loads((parent / "settings.json").read_text())
    stored["settings"]["alpha"] = 0.1
    stored["hash"] = lsc.settings_hash(stored["settings"])
    (parent / "settings.json").write_text(json.dumps(stored))
    with pytest.raises(SystemExit, match=r"stage_span/settings.json has changed since the settings were written"):
        ses.run_shard(tmp_path / "run", settings, 0, 1, limit=0)
    stored["hash"] = "0" * 64
    (parent / "settings.json").write_text(json.dumps(stored))
    with pytest.raises(SystemExit, match=r"settings.json fails its hash"):
        _small_saturated(parent_dir=parent)


def test_the_repository_parent_is_bound(saturated) -> None:
    assert saturated["parent_files"] == {
        f"results/univariate-baseline-2026-11-23/stage_span/{name}": ssc.sha256_file(ses.PARENT_DIR / name)
        for name in ses.PARENT_FILE_NAMES
    }
    assert saturated["parent_dir"] == "results/univariate-baseline-2026-11-23/stage_span"


# ── The unit (real units, fast knobs) ────────────────────────────────────────


UNIT = lsc.Unit(ses.RERUN_DESIGN, "shape_e0.5", 1)


def _fast(settings):
    settings["unit_settings"]["units"]["pls_cv_params"] = dict(FAST_CV)
    return settings


@pytest.fixture(scope="module")
def unit_records():
    import motco.simulations.evaluation as evaluation

    mine = _fast(ses.saturated_settings(designs=[ses.RERUN_DESIGN, "n320_delta1"], replicates=2,
                                        control_replicates=2, permutations=PERMS))
    parent = ubc.stage_span_settings(designs=[ses.RERUN_DESIGN], replicates=2, control_replicates=2,
                                     permutations=PERMS)
    parent["units"]["pls_cv_params"] = dict(FAST_CV)
    real = evaluation.plsda_doubleCV
    calls = []

    def counting(*args, **kwargs):
        calls.append(1)
        return real(*args, **kwargs)

    with mock.patch.object(evaluation, "plsda_doubleCV", counting):
        record = ses.evaluate_unit(UNIT, mine["unit_settings"])
    return record, ubc.evaluate_unit(UNIT, parent), len(calls), mine


def test_a_unit_record_equals_the_parent_drivers(unit_records) -> None:
    record, parent, _, mine = unit_records
    assert record["status"] == "completed", record.get("error_message")
    assert record["settings_hash"] == lsc.settings_hash(mine["unit_settings"])

    def strip(r):
        return {key: value for key, value in r.items() if key not in ("settings_hash", "runtime_seconds")}

    assert strip(record) == strip(parent)
    assert record["generator_seed"] == lsc.generator_seed(2043, UNIT.design, UNIT.replicate)
    assert set(record["spaces"]) == {"pca_matched"}


def test_one_cross_validated_fit_per_unit(unit_records) -> None:
    assert unit_records[2] == 1


def test_a_failed_unit_is_recorded_and_the_shard_continues(tmp_path) -> None:
    settings = _fast(_small_saturated())

    def broken(*args, **kwargs):
        raise RuntimeError("patched scan failure")

    with mock.patch("motco.stats.interaction.interaction_feature_scan", broken):
        assert ses.run_shard(tmp_path, settings, 0, 26, limit=2) == 2
    records = [json.loads(line) for line in next((tmp_path / "saturated" / "records").glob("*.jsonl")).open()]
    assert len(records) == 2
    for record in records:
        assert record["status"] == "failed"
        assert (record["error_type"], record["error_message"]) == ("RuntimeError", "patched scan failure")


# ── Merge and reproduction ───────────────────────────────────────────────────


def _write_shard(out_dir: Path, settings, records, name="shard_000_of_001.jsonl") -> None:
    directory = ses.profile_dir(out_dir, settings["profile"])
    lsc.ensure_settings(directory, settings)
    (directory / "records").mkdir(parents=True, exist_ok=True)
    with (directory / "records" / name).open("a", encoding="utf-8") as handle:
        for record in records:
            handle.write(lsc._line(record))


def test_merge_writes_one_row_per_design_cell_and_arm(tmp_path) -> None:
    settings = _small_saturated()
    units = settings["unit_settings"]
    records = []
    for unit in ses.profile_units(settings):
        if (unit.design, unit.cell, unit.replicate) == ("n320_delta1", "magnitude_e0.1", 1):
            records.append(_fake_record(unit, units, status="failed"))
        elif unit.cell == "magnitude_e0.1":
            records.append(_fake_record(unit, units, rejections={"production": True, "fdr": True}, rank=4))
        else:
            records.append(_fake_record(unit, units, rank=2 + unit.replicate))
    _write_shard(tmp_path, settings, records)
    _write_shard(tmp_path, settings, records[:3], name="shard_dup.jsonl")
    summary = ses.merge(tmp_path, "saturated")
    directory = tmp_path / "saturated"
    assert summary["failed_units"] == 1 and summary["complete_units"] == 52
    accounting = json.loads((directory / "accounting.json").read_text())
    assert accounting["failures_by_cell"] == {"n320_delta1/magnitude_e0.1": 1}
    assert accounting["duplicate_records"] == 3 and accounting["profile"] == "saturated"
    assert accounting["settings_hash"] == lsc.settings_hash(settings)
    rates = lsc._read_csv(directory / "rate_table.csv")
    assert list(rates[0]) == list(ubc.RATE_COLUMNS) and len(rates) == 2 * 13 * len(ubc.ARMS)
    assert len({(r["design"], r["cell"], r["arm"]) for r in rates}) == len(rates)
    by = {(r["design"], r["cell"], r["arm"]): r for r in rates}
    failed = by["n320_delta1", "magnitude_e0.1", "production"]
    assert (failed["count"], failed["rejections"], failed["failures"]) == ("1", "1", "1")
    assert len(lsc._read_csv(directory / "paired_table.csv")) == 2 * 13 * len(ubc.BASELINE_ARMS)
    cells = lsc._read_csv(directory / "cell_table.csv")
    assert len(cells) == 26 and list(cells[0]) == list(ubc.CELL_COLUMNS)
    ranks = lsc._read_csv(directory / "rank_table.csv")
    assert list(ranks[0]) == list(ses.RANK_COLUMNS)
    stratum = {(r["rank"], r["arm"]): r for r in ranks if (r["design"], r["cell"]) == ("n320_delta1", "anchor")}
    assert {key[0] for key in stratum} == {"2", "3"}
    assert (stratum["2", "production"]["count"], stratum["3", "production"]["count"]) == ("1", "1")
    magnitude = [r for r in ranks if (r["design"], r["cell"], r["arm"]) == ("n320_delta1", "magnitude_e0.1", "fdr")]
    assert [(r["rank"], r["count"], r["rejections"]) for r in magnitude] == [("4", "1", "1")]
    reproduction = json.loads((directory / "reproduction.json").read_text())
    assert reproduction["passed"] is False and reproduction["rows_compared"] == 13 * len(ubc.ARMS)


def test_the_band_merge_writes_no_reproduction(tmp_path, fake_units) -> None:
    _write_calibration(tmp_path)
    settings = ses.band_settings(tmp_path, replicates=1, control_replicates=1, permutations=9)
    ses.run_shard(tmp_path, settings, 0, 1)
    summary = ses.merge(tmp_path, "band")
    assert summary["complete_units"] == 33 and summary["reproduction_passed"] is None
    assert not (tmp_path / "band" / "reproduction.json").exists()


@pytest.fixture(scope="module")
def committed_rates():
    return lsc._read_csv(ses.PARENT_DIR / "rate_table.csv")


def test_the_committed_table_resolves_every_compared_row(committed_rates) -> None:
    rows = {(r["cell"], r["arm"]) for r in committed_rates if r["design"] == ses.RERUN_DESIGN}
    assert rows == {(cell[0], arm) for cell in ubc.STAGE_SPAN_CELLS for arm in ubc.ARMS}


def test_the_reproduction_passes_on_the_committed_counts(saturated, committed_rates) -> None:
    mine = [{**r, "count": int(r["count"]), "rejections": int(r["rejections"])}
            for r in committed_rates if r["design"] == ses.RERUN_DESIGN]
    result = ses.rerun_reproduction(mine, committed_rates, saturated)
    assert result["passed"] is True and result["rows_compared"] == 234
    assert result["differing"] == [] and result["unresolved_rows"] == []


def test_the_reproduction_fails_and_lists_every_differing_row(saturated, committed_rates) -> None:
    mine = [dict(r) for r in committed_rates if r["design"] == ses.RERUN_DESIGN]
    target = next(r for r in mine if (r["cell"], r["arm"]) == ("shape_e0.25", "production"))
    wanted = int(target["rejections"])
    target["rejections"] = str(wanted + 1)
    dropped = mine.pop(next(i for i, r in enumerate(mine) if (r["cell"], r["arm"]) == ("anchor", "fdr")))
    result = ses.rerun_reproduction(mine, committed_rates, saturated)
    assert result["passed"] is False and result["rows_compared"] == 233
    assert result["differing"] == [{"cell": "shape_e0.25", "arm": "production",
                                    "committed": {"count": 200, "rejections": wanted},
                                    "reproduced": {"count": 200, "rejections": wanted + 1}}]
    assert result["unresolved_rows"] == [{"cell": dropped["cell"], "arm": "fdr", "missing_from": "this run"}]
    absent = ses.rerun_reproduction([], committed_rates, ses.saturated_settings(designs=["n320_delta1"]))
    assert absent["passed"] is False and "not in this run" in absent["reason"]


# ── The decision on synthetic tables ─────────────────────────────────────────


def _rows(units, *, signal=100):
    rates, paired = [], []
    for design in units["designs"]:
        for cell in ubc._cell_names(units):
            control = cell in ubc.CONTROL_CELLS
            count = 500 if control else 200
            for arm in ubc.ARMS:
                rates.append({"design": design, "cell": cell, "arm": arm, "count": str(count),
                              "rejections": str(25 if control else signal)})
            for arm in ubc.BASELINE_ARMS:
                paired.append({"design": design, "cell": cell, "arm": arm, "n_paired": str(count),
                               "rate_production": "0.5", "rate_arm": "0.5", "b": "0", "c": "0"})
    return rates, paired


def _block(settings, **extra):
    rates, paired = _rows(settings["unit_settings"])
    return {"settings": settings, "accounting": {"missing_units": 0, "failed_units": 0}, "rates": rates,
            "paired": paired, "cells": [], "ranks": [], **extra}


@pytest.fixture()
def blocks(saturated, band, committed):
    committed_rates, _ = _rows(committed)
    return {
        "saturated": _block(saturated, reproduction={"passed": True, "rows_compared": 234, "differing": []}),
        "band": _block(band),
        "committed": {"settings": committed, "rates": committed_rates, "cells": []},
        "calibration": {"steps": _steps({0.05: 0.63}), "adopted": {"delta": 0.05}, "failures": 0},
    }


def _decide(blocks):
    return ses.decision(blocks["saturated"], blocks["band"], blocks["committed"], blocks["calibration"])


def _set(block, design, cell, arm, rejections, count=None):
    for row in block["rates"]:
        if (row["design"], row["cell"], row["arm"]) == (design, cell, arm):
            row["rejections"] = str(rejections)
            if count is not None:
                row["count"] = str(count)
            return
    raise AssertionError((design, cell, arm))


def _move(blocks, weaker, cell, *, weak, match, arm="production"):
    """Set one saturated pair-cell; the other design matched to the same committed design follows the match."""

    reference = dict(ses.SATURATED_PAIRS)[weaker]
    _set(blocks["saturated"], weaker, cell, arm, weak)
    _set(blocks["committed"], reference, cell, arm, match)
    for other, shared in ses.SATURATED_PAIRS:
        if shared == reference and other != weaker:
            _set(blocks["saturated"], other, cell, arm, match)


def _set_pair(block, design, cell, arm, *, b=0, c=0, n=200):
    for row in block["paired"]:
        if (row["design"], row["cell"], row["arm"]) == (design, cell, arm):
            row.update({"n_paired": str(n), "b": str(b), "c": str(c),
                        "rate_production": str((100 + c) / n), "rate_arm": str((100 + b) / n)})
            return
    raise AssertionError((design, cell, arm))


def _entry(result, weaker, cell):
    return next(e for e in result["criteria"]["pair_cells"] if (e["weaker"], e["cell"]) == (weaker, cell))


def test_the_pairs_are_the_fifty_three_pair_cells(blocks) -> None:
    result = _decide(blocks)
    cells = result["criteria"]["pair_cells"]
    assert result["m_pair_cells"] == len(cells) == ses.PAIR_CELLS_EXPECTED == 53
    assert result["label_level"] == pytest.approx(0.025 / 53)
    by_pair: dict[str, int] = {}
    for entry in cells:
        by_pair[entry["pair"]] = by_pair.get(entry["pair"], 0) + 1
    assert by_pair == {"n320_delta1~n80": 11, "n1280_delta0.5~n80": 11, "n800_delta1~n200": 11,
                       "n800_delta0.25~n200_delta0.5": 11, "band_n1280_half~band_n320": 9}
    assert not {entry["cell"] for entry in cells} & set(ubc.CONTROL_CELLS)
    assert {entry["family"] for entry in cells} == {"magnitude", "orientation", "shape", "out_of_span"}


def test_no_label_reads_scaling_unrefuted(blocks) -> None:
    result = _decide(blocks)
    assert result["verdict"] == "scaling_unrefuted"
    assert "not excluded" in result["consequence"]
    assert all(entry["label"] == "unresolved" for entry in result["criteria"]["pair_cells"])
    assert result["named"] == {"weak_costs": [], "weak_gains": []}


def test_the_hand_computed_exact_test_in_each_direction() -> None:
    def lower_tail(a, n1, k, n2):
        total = n1 + n2
        return sum(math.comb(k, x) * math.comb(total - k, n1 - x) for x in range(0, a + 1)) / math.comb(total, n1)

    def upper_tail(a, n1, k, n2):
        total = n1 + n2
        return sum(math.comb(k, x) * math.comb(total - k, n1 - x)
                   for x in range(a, min(k, n1) + 1)) / math.comb(total, n1)

    p = ses.fisher_p_values((20, 6), (20, 14))
    assert p["p_costs"] == pytest.approx(lower_tail(6, 20, 20, 20), rel=1e-9)
    assert p["p_gains"] == pytest.approx(upper_tail(6, 20, 20, 20), rel=1e-9)
    costs = ses.exact_label({"count": 200, "rejections": 82}, {"count": 200, "rejections": 164}, level=0.025 / 53)
    assert costs["label"] == "weak_costs" and costs["difference"] == pytest.approx(-0.41)
    assert costs["p_costs"] == pytest.approx(lower_tail(82, 200, 246, 200), rel=1e-6)
    gains = ses.exact_label({"count": 200, "rejections": 164}, {"count": 200, "rejections": 82}, level=0.025 / 53)
    assert gains["label"] == "weak_gains" and gains["p_gains"] == pytest.approx(costs["p_costs"])


def test_a_significant_difference_below_the_margin_is_unresolved() -> None:
    reading = ses.exact_label({"count": 20_000, "rejections": 9_000}, {"count": 20_000, "rejections": 10_000},
                              level=0.025 / 53)
    assert reading["p_costs"] < 1e-20 and reading["difference"] == pytest.approx(-0.05)
    assert reading["label"] == "unresolved"
    edge = ses.exact_label({"count": 200, "rejections": 80}, {"count": 200, "rejections": 100}, level=1.0)
    assert edge["label"] == "weak_costs"


def test_a_cost_in_a_saturated_pair_reads_weak_stage_costs(blocks) -> None:
    _move(blocks, "n1280_delta0.5", "magnitude_e0.25", weak=82, match=164)
    result = _decide(blocks)
    entry = _entry(result, "n1280_delta0.5", "magnitude_e0.25")
    assert (entry["label"], entry["control_label"], entry["counted"]) == ("weak_costs", "unresolved", True)
    assert entry["production"]["match"] == {"count": 200, "rejections": 164, "rate": 0.82}
    assert result["verdict"] == "weak_stage_costs"
    assert result["named"]["weak_costs"] == [["n1280_delta0.5~n80", "magnitude_e0.25"]]


def test_the_saturated_references_are_the_committed_counts(blocks) -> None:
    # The match's rate in this run's saturated table is never read: only the committed table is.
    _set(blocks["saturated"], "n1280_delta0.5", "magnitude_e0.25", "production", 82)
    for row in blocks["saturated"]["rates"]:
        if row["design"] == ses.RERUN_DESIGN:
            row["rejections"] = "200"
    entry = _entry(_decide(blocks), "n800_delta0.25", "shape_e0.5")
    assert entry["production"]["match"]["rejections"] == 100


def test_a_label_the_control_shares_does_not_count(blocks) -> None:
    for arm in ("production", "fdr"):
        _set(blocks["saturated"], "n800_delta1", "shape_e0.1", arm, 82)
        _set(blocks["committed"], "n200", "shape_e0.1", arm, 164)
    result = _decide(blocks)
    entry = _entry(result, "n800_delta1", "shape_e0.1")
    assert (entry["label"], entry["control_label"], entry["counted"]) == ("weak_costs", "weak_costs", False)
    assert result["verdict"] == "scaling_unrefuted"
    assert result["criteria"]["label_cells"]["shared_with_control"] == [["n800_delta1~n200", "shape_e0.1"]]


def test_a_label_the_control_opposes_counts(blocks) -> None:
    _set(blocks["saturated"], "n800_delta1", "shape_e0.1", "production", 82)
    _set(blocks["committed"], "n200", "shape_e0.1", "production", 164)
    _set(blocks["saturated"], "n800_delta1", "shape_e0.1", "fdr", 164)
    _set(blocks["committed"], "n200", "shape_e0.1", "fdr", 82)
    result = _decide(blocks)
    entry = _entry(result, "n800_delta1", "shape_e0.1")
    assert (entry["label"], entry["control_label"], entry["counted"]) == ("weak_costs", "weak_gains", True)
    assert result["verdict"] == "weak_stage_costs"


def test_gains_and_mixed(blocks) -> None:
    _set(blocks["band"], "band_n1280_half", "orient_e0.5_phi0", "production", 164)
    _set(blocks["band"], "band_n320", "orient_e0.5_phi0", "production", 82)
    result = _decide(blocks)
    assert result["verdict"] == "weak_stage_gains"
    assert result["named"]["weak_gains"] == [["band_n1280_half~band_n320", "orient_e0.5_phi0"]]
    _move(blocks, "n320_delta1", "orient_e0.5_phi90", weak=82, match=164)
    result = _decide(blocks)
    assert result["verdict"] == "mixed"
    assert result["named"]["weak_costs"] == [["n320_delta1~n80", "orient_e0.5_phi90"]]


def test_an_unavailable_reading_is_unresolved_and_named(blocks) -> None:
    _set(blocks["saturated"], "n320_delta1", "shape_e0.25", "production", 0, count=0)
    result = _decide(blocks)
    entry = _entry(result, "n320_delta1", "shape_e0.25")
    assert entry["label"] == "unresolved" and entry["production"]["unavailable"] is True
    assert result["criteria"]["label_cells"]["unavailable"] == [["n320_delta1~n80", "shape_e0.25"]]


@pytest.mark.parametrize("family, design", [("saturated", "n800_delta0.25"), ("band", "band_n1280")])
def test_production_size_excess_in_either_family(blocks, family, design) -> None:
    _set(blocks[family], design, "anchor", "production", 100)
    _move(blocks, "n1280_delta0.5", "magnitude_e0.25", weak=82, match=164)
    result = _decide(blocks)
    assert result["verdict"] == "production_size_excess"
    failing = result["failing"]["production_size"]
    assert list(failing) == [family] and failing[family][0]["design"] == design
    assert result["criteria"]["size"]["production"]["families"]["saturated"]["m"] == 8
    assert result["criteria"]["size"]["production"]["families"]["band"]["m"] == 6
    assert result["named"]["weak_costs"]


def test_the_rerun_design_is_outside_the_size_family(blocks) -> None:
    _set(blocks["saturated"], ses.RERUN_DESIGN, "anchor", "production", 100)
    result = _decide(blocks)
    assert result["verdict"] == "scaling_unrefuted"
    observed = result["criteria"]["size"]["production"]["families"]["saturated"]["observations"]
    assert {row["design"] for row in observed} == set(NEW_SATURATED)


def test_an_inadmissible_per_feature_arm_leaves_the_map(blocks) -> None:
    _set(blocks["band"], "band_n320", "translation", "fdr", 120)
    result = _decide(blocks)
    assert result["verdict"] == "scaling_unrefuted"
    assert list(result["inadmissible_per_feature_arms"]) == ["fdr"]
    assert list(result["inadmissible_per_feature_arms"]["fdr"]) == ["band"]
    assert result["map"]["admissible_arms"] == ["total"]
    assert all(set(entry["against"]) <= {"total"} for entry in result["map"]["labels"])
    _set(blocks["saturated"], "n800_delta1", "anchor", "total", 120)
    result = _decide(blocks)
    assert result["map"]["admissible_arms"] == []
    assert "no admissible per-feature arm" in result["map"]["reason"]
    assert all(entry["label"] == "unresolved" for entry in result["map"]["labels"])


def test_the_gate_still_records_size_labels_grades_and_map(blocks) -> None:
    blocks["saturated"]["reproduction"] = {
        "passed": False, "rows_compared": 234,
        "differing": [{"cell": "anchor", "arm": "production", "committed": {"count": 500, "rejections": 25},
                       "reproduced": {"count": 500, "rejections": 26}}]}
    _move(blocks, "n1280_delta0.5", "magnitude_e0.25", weak=82, match=164)
    result = _decide(blocks)
    assert result["verdict"] == "reproduction_failed"
    assert result["failing"]["reproduction"]["differing"][0]["cell"] == "anchor"
    assert result["reproduction"] == {"passed": False, "rows_compared": 234}
    assert result["criteria"]["size"]["production"]["admissible"] is True
    assert result["named"]["weak_costs"] == [["n1280_delta0.5~n80", "magnitude_e0.25"]]
    assert len(result["grades"]) == 11 and result["map"]["m_cells"] == 71
    assert set(result["consequences"]) == set(ses.VERDICTS)


# ── Grades per design ────────────────────────────────────────────────────────


def _detect(block, design, *, off_target=10):
    units = block["settings"].get("unit_settings", block["settings"])
    for cell in ubc._cell_names(units):
        meta = ubc._cell_meta({"settings": units}, design, cell)
        mode = meta["trajectory_mode"]
        if mode not in ubc.TARGETS or (meta["span_alignment"] or 0) > 0:
            continue
        _set(block, design, cell, "production", 180)
        for statistic in ubc.OFF_TARGET[mode]:
            _set(block, design, cell, statistic, off_target)


def _grade(result, design):
    return next(entry for entry in result["grades"] if entry["design"] == design)


def test_every_design_is_graded_on_its_own(blocks) -> None:
    result = _decide(blocks)
    assert [(e["design"], e["source"]) for e in result["grades"]] == [
        *((d, "committed") for d in ("n80", "n200", "n500", "n200_delta0.5")),
        *((d, "this run, saturated profile") for d in NEW_SATURATED),
        *((d, "this run, band profile") for d in ("band_n320", "band_n1280_half", "band_n1280")),
    ]
    assert {e["grade"] for e in result["grades"]} == {"no_detected_cell"}
    assert _grade(result, "n800_delta1")["considered_cells"] == 9
    assert _grade(result, "band_n320")["considered_cells"] == 9


def test_names_kind_and_detects_only(blocks) -> None:
    _detect(blocks["saturated"], "n800_delta1")
    _detect(blocks["band"], "band_n1280", off_target=60)
    result = _decide(blocks)
    assert _grade(result, "n800_delta1")["grade"] == "names_kind"
    assert len(_grade(result, "n800_delta1")["detected_cells"]) == 9
    assert _grade(result, "n800_delta1")["family"]["m"] == 18
    detects = _grade(result, "band_n1280")
    assert detects["grade"] == "detects_only" and len(detects["failing"]) == 18
    assert _grade(result, "n320_delta1")["grade"] == "no_detected_cell"


def test_a_committed_and_a_new_design_with_the_same_counts_get_the_same_grade(blocks) -> None:
    for block, design in ((blocks["committed"], "n80"), (blocks["saturated"], "n320_delta1")):
        _detect(block, design, off_target=30)
    result = _decide(blocks)
    old, new = _grade(result, "n80"), _grade(result, "n320_delta1")
    assert old["grade"] == new["grade"]
    assert [row["status"] for row in old["family"]["observations"]] == \
        [row["status"] for row in new["family"]["observations"]]


def test_out_of_span_cells_never_enter_the_grade(blocks) -> None:
    for cell in ("orient_e0.5_phi90", "orient_e0.5_phi90_k20"):
        _set(blocks["saturated"], "n800_delta1", cell, "production", 200)
        for statistic in ("delta", "shape"):
            _set(blocks["saturated"], "n800_delta1", cell, statistic, 200)
    result = _decide(blocks)
    assert _grade(result, "n800_delta1")["grade"] == "no_detected_cell"


def test_off_target_rates_never_change_the_verdict(blocks) -> None:
    _move(blocks, "n1280_delta0.5", "magnitude_e0.25", weak=82, match=164)
    _detect(blocks["saturated"], "n800_delta1")
    _detect(blocks["committed"], "n200")
    before = _decide(blocks)
    _detect(blocks["saturated"], "n800_delta1", off_target=150)
    after = _decide(blocks)
    assert _grade(before, "n800_delta1")["grade"] != _grade(after, "n800_delta1")["grade"]
    assert before["verdict"] == after["verdict"] == "weak_stage_costs"
    assert before["criteria"]["pair_cells"] == after["criteria"]["pair_cells"]


# ── The map and the descriptive tables ───────────────────────────────────────


def test_the_map_has_one_row_per_family_and_new_design(blocks) -> None:
    _set_pair(blocks["saturated"], "n800_delta1", "magnitude_e0.25", "fdr", c=60)
    _set_pair(blocks["saturated"], "n800_delta1", "magnitude_e0.25", "total", c=60)
    _set_pair(blocks["band"], "band_n320", "shape_e1", "total", b=60)
    result = _decide(blocks)
    reading = result["map"]
    assert reading["m_cells"] == 71 and len(reading["labels"]) == 71
    assert reading["admissible_arms"] == ["fdr", "total"]
    assert reading["levels"]["production_leads"] == pytest.approx(0.025 / 71)
    assert {(row["profile"], row["design"]) for row in reading["table"]} == {
        *(("saturated", d) for d in NEW_SATURATED), *(("band", d) for d in ses.BAND_DESIGNS)}
    assert len(reading["table"]) == 4 * 4 + 3 * 3
    assert sum(row[label] for row in reading["table"] for label in ubc.LABELS) == 71
    label = {(e["design"], e["cell"]): e["label"] for e in reading["labels"]}
    assert label["n800_delta1", "magnitude_e0.25"] == "production_leads"
    assert label["band_n320", "shape_e1"] == "baseline_leads"
    assert {e["profile"] for e in reading["labels"]} == {"saturated", "band"}
    assert not any(e["design"] == ses.RERUN_DESIGN for e in reading["labels"])


def test_the_descriptive_tables_have_the_rows_they_name(blocks) -> None:
    for block, design in ((blocks["saturated"], "n320_delta1"), (blocks["band"], "band_n1280")):
        block["cells"] = [{"design": design, "cell": "anchor", "auroc_median": "0.64", "rank_q10": "9",
                           "rank_median": "12", "rank_q90": "19", "capture_pooled_stage_median": "0.6",
                           "capture_group_difference_median": "", "completed": "500"}]
    blocks["band"]["ranks"] = [
        {"design": "band_n1280_half", "cell": "magnitude_e1", "rank": "19", "arm": "production", "count": "150",
         "rejections": "30"},
        {"design": "band_n1280_half", "cell": "magnitude_e1", "rank": "18", "arm": "production", "count": "50",
         "rejections": "12"},
    ]
    _set(blocks["band"], "band_n1280_half", "magnitude_e1", "production", 42)
    _set(blocks["band"], "band_n320", "magnitude_e1", "production", 120)
    result = _decide(blocks)
    tables = result["descriptive"]
    index = {(row["design"], row["source"]): row for row in tables["index"]}
    assert len(tables["index"]) == 4 + 5 + 3
    row = index["band_n1280", "this run, band profile"]
    assert (row["auroc_median"], row["rank_median"], row["n_samples"], row["delta"]) == (0.64, 12.0, 1280, 0.05)
    assert row["delta_sqrt_n"] == pytest.approx(0.05 * math.sqrt(1280))
    assert index["n80", "committed"]["delta_sqrt_n"] == pytest.approx(2 * math.sqrt(80))
    assert len(tables["cell_rates"]) == 5 * 13 + 3 * 11
    assert set(tables["cell_rates"][0]["rates"]) == set(ubc.ARMS)
    assert len(tables["stretch"]) == 4 * 3 + 4 * 3 + 3 * 3
    stretch = next(s for s in tables["stretch"] if (s["design"], s["cell"]) == ("band_n1280", "magnitude_e2"))
    assert stretch["e_delta_sqrt_n"] == pytest.approx(2 * 0.05 * math.sqrt(1280))
    assert set(stretch["rates"]) == set(ses.STRETCH_ARMS)
    assert result["named"]["weak_costs"] == [["band_n1280_half~band_n320", "magnitude_e1"]]
    strata, = tables["rank_strata"]
    assert [entry["rank"] for entry in strata["ranks"]] == [18, 19]
    assert strata["ranks"][1]["arms"]["production"]["rate"] == pytest.approx(0.2)
    assert tables["calibration"]["adopted"] == {"delta": 0.05}


def test_recorded_arms_never_change_the_verdict_grades_or_map(blocks) -> None:
    _move(blocks, "n1280_delta0.5", "magnitude_e0.25", weak=82, match=164)
    _detect(blocks["saturated"], "n800_delta1")
    _detect(blocks["committed"], "n200")
    before = _decide(blocks)
    recorded = {"scan", "production_minp", "span", "complement",
                *(arm for arm in ubc.ARMS if arm.startswith(("pca_matched", "concat")))}
    for block in (blocks["saturated"], blocks["band"], blocks["committed"]):
        for row in block["rates"]:
            if row["arm"] in recorded:
                row["rejections"] = str(int(row["count"]))
    for block in (blocks["saturated"], blocks["band"]):
        for row in block["paired"]:
            if row["arm"] == "scan":
                row.update({"b": "90", "c": "0", "rate_arm": "0.95"})
    after = _decide(blocks)
    assert after["verdict"] == before["verdict"] == "weak_stage_costs"
    assert after["criteria"] == before["criteria"]
    assert [(g["design"], g["grade"]) for g in after["grades"]] == [(g["design"], g["grade"]) for g in before["grades"]]
    assert after["map"]["labels"] == before["map"]["labels"]


# ── End to end on fake units ─────────────────────────────────────────────────


def test_calibrate_run_merge_and_decide_end_to_end(tmp_path, monkeypatch) -> None:
    def fake_calibration(delta, replicate, settings):
        return {"delta": delta, "replicate": replicate, "generator_seed": 0, "status": "completed",
                "auroc": 0.5 + 2.5 * delta, "rank": 10, "seconds": 0.1}

    monkeypatch.setattr(ses, "calibration_unit", fake_calibration)
    monkeypatch.setattr(ses, "evaluate_unit", lambda unit, units: _fake_record(unit, units))
    assert ses.main(["calibrate", "--out-dir", str(tmp_path), "--reps", "2"]) == 0
    assert ses.load_calibration(tmp_path / ses.CALIBRATION_FILE)["adopted"]["delta"] == 0.05
    for profile in ses.PROFILES:
        assert ses.main(["run", "--profile", profile, "--out-dir", str(tmp_path), "--reps", "1",
                         "--control-reps", "1", "--perms", "9"]) == 0
        assert ses.main(["merge", "--profile", profile, "--out-dir", str(tmp_path)]) == 0
    assert ses.main(["decide", "--out-dir", str(tmp_path)]) == 0
    result = json.loads((tmp_path / ses.DECISION_NAME).read_text())
    assert result["verdict"] == "reproduction_failed" and result["complete"] is True
    assert result["m_pair_cells"] == 53 and result["calibration"]["delta"] == 0.05
    assert result["descriptive"]["calibration"]["adopted"]["delta"] == 0.05


# ── The sbatch array ─────────────────────────────────────────────────────────


SBATCH = REPO / "scripts" / "stage_effect_scaling_array.sbatch"


def test_the_sbatch_array_parses() -> None:
    subprocess.run(["bash", "-n", str(SBATCH)], check=True)


@pytest.mark.parametrize("profile", ses.PROFILES)
def test_the_sbatch_dry_run_prints_the_run_command(tmp_path, profile) -> None:
    env = {**os.environ, "SES_PROFILE": profile, "SES_OUT": str(tmp_path / "run"), "N_SHARDS": "200",
           "SLURM_ARRAY_TASK_ID": "7", "SES_DRY_RUN": "1", "SES_ARGS": "--limit 1", "SLURM_SUBMIT_DIR": str(REPO)}
    result = subprocess.run(["bash", str(SBATCH)], cwd=tmp_path, env=env, capture_output=True, text=True, check=True)
    command = result.stdout.strip()
    assert command.startswith(f"PYTHONPATH={REPO}/src python {REPO}/scripts/stage_effect_scaling_check.py run")
    assert f"--profile {profile}" in command and "--shard 7/200" in command and "--limit 1" in command
    assert f"--out-dir {tmp_path / 'run'}" in command
    for bad, message in (({"SES_PROFILE": "surgery"}, "SES_PROFILE must be"),
                         ({"SES_OUT": "relative/run"}, "SES_OUT must be an absolute path")):
        failed = subprocess.run(["bash", str(SBATCH)], cwd=tmp_path, env={**env, **bad}, capture_output=True,
                                text=True)
        assert failed.returncode == 2 and message in failed.stderr
