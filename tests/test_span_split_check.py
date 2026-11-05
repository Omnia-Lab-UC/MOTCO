"""``scripts/span_split_check.py``: units, the simulation and plasmode units, resume, merge, reproduction, decision.

The simulation unit runs on the comparison's ``seaad`` design at fast CV knobs;
the plasmode profile runs on the plasmode tests' synthetic stand-in export. No
real data are read.
"""

from __future__ import annotations

import json
import shutil
import subprocess
import sys
from pathlib import Path

import numpy as np
import pytest


def _script(name: str):
    sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
    try:
        return __import__(name)
    finally:
        sys.path.pop(0)


ssc = _script("span_split_check")
lsc = ssc.lsc
pc = ssc.pc

FAST_CV = {"cv1_splits": 3, "cv2_splits": 4, "n_repeats": 1, "max_components": 6, "random_state": 1203}
REPO = Path(__file__).resolve().parents[1]


# ── Simulation settings and units ────────────────────────────────────────────


def test_simulation_units_and_grade() -> None:
    settings = ssc.simulation_settings()
    units = ssc.simulation_units(settings)
    assert len(units) == 3 * (2 * 500 + 12 * 100)
    for design in lsc.DESIGNS:
        for cell, *_ in lsc.CELLS:
            replicates = sorted(u.replicate for u in units if (u.design, u.cell) == (design, cell))
            assert replicates == list(range(500 if cell in ("anchor", "translation") else 100))
    assert settings["permutations"] == 999 and settings["reproduction_permutations"] == 199
    parent = json.loads((ssc.COMPARISON_DIR / "settings.json").read_text())["settings"]
    assert settings["comparison"]["weak_delta"] == parent["weak_delta"]
    assert settings["comparison"]["weak_delta"]["value"] == 0.11
    assert settings["comparison"]["arms"] == ["pls_cv"]
    assert settings["rejection"] == "p < alpha"


def _parent_copy(tmp_path: Path) -> Path:
    target = tmp_path / "comparison"
    (target / "calibration").mkdir(parents=True)
    for name in ("settings.json", "rate_table.csv", "calibration/calibration.json"):
        shutil.copy(ssc.COMPARISON_DIR / name, target / name)
    return target


def test_a_calibration_whose_hash_differs_is_refused(tmp_path) -> None:
    parent = _parent_copy(tmp_path)
    calibration = parent / "calibration" / "calibration.json"
    calibration.write_text(calibration.read_text() + " ")
    with pytest.raises(SystemExit, match="not the comparison's calibration_sha256 2e440c50"):
        ssc.simulation_settings(parent_dir=parent)


def test_a_parent_whose_settings_fail_their_hash_is_refused(tmp_path) -> None:
    parent = _parent_copy(tmp_path)
    stored = json.loads((parent / "settings.json").read_text())
    stored["settings"]["base_seed"] = 1
    (parent / "settings.json").write_text(json.dumps(stored))
    with pytest.raises(SystemExit, match="fails its hash"):
        ssc.simulation_settings(parent_dir=parent)


def test_shared_units_are_the_comparisons_datasets() -> None:
    settings = ssc.simulation_settings(designs=["seaad"])
    unit = lsc.Unit("seaad", "orient_e0.25_phi30", 7)
    mine = lsc.unit_dataset(unit, settings["comparison"])
    parent = json.loads((ssc.COMPARISON_DIR / "settings.json").read_text())["settings"]
    theirs = lsc.unit_dataset(unit, parent)
    for name in ("methylation", "expression", "proteomics", "metadata"):
        assert getattr(mine, name).equals(getattr(theirs, name))


def _small_settings(**overrides):
    settings = ssc.simulation_settings(designs=["seaad"], replicates=2, control_replicates=2, permutations=39)
    settings["comparison"]["pls_cv_params"] = dict(FAST_CV)
    settings["reproduction_permutations"] = 19
    settings.update(overrides)
    return settings


UNIT = lsc.Unit("seaad", "orient_e0.5_phi90", 1)


@pytest.fixture(scope="module")
def simulation_unit():
    settings = _small_settings()
    return settings, ssc.evaluate_simulation_unit(UNIT, settings)


def test_simulation_record_carries_every_field(simulation_unit) -> None:
    settings, record = simulation_unit
    assert record["status"] == "completed", record.get("error_message")
    assert (record["design"], record["cell"], record["replicate"]) == ("seaad", "orient_e0.5_phi90", 1)
    assert record["generator_seed"] == lsc.generator_seed(2035, "seaad", 1)
    assert record["rrpp_seed"] == lsc.rrpp_seed(2035, "seaad", "orient_e0.5_phi90", 1)
    production = record["production"]
    assert {"rank", "cv_mean_auroc", "pair_statistics", "p_values", "p_values_reproduction", "span_capture",
            "population_latent"} <= set(production)
    assert set(production["p_values"]) == set(production["p_values_reproduction"]) == {"delta", "angle", "shape"}
    for part in ssc.PARTS:
        entry = record["split"]["parts"][part]
        assert {"dimension", "ss_interaction", "rss", "f", "p_value_f", "p_value_ss", "p_value_ss_df"} <= set(entry)
    assert record["split"]["parts"]["span"]["dimension"] == production["rank"]
    assert record["split"]["parts"]["total"]["dimension"] == 367 + 131
    assert 0.0 <= record["split"]["span_share"] <= 1.0


def test_production_arm_equals_the_comparisons_unit(simulation_unit) -> None:
    settings, record = simulation_unit
    comparison = {**settings["comparison"], "permutations": 19, "arms": ["pls_cv"]}
    theirs = lsc.evaluate_unit(UNIT, comparison)[0]
    production = record["production"]
    assert production["rank"] == theirs["retained_dimensionality"]
    assert production["pair_statistics"] == theirs["pair_statistics"]
    assert production["span_capture"] == theirs["span_capture"]
    assert production["population_latent"] == theirs["population_latent"]
    # The first 19 of the unit's 39 draws are the comparison's 19-draw sequence.
    assert production["p_values_reproduction"] == theirs["p_values"]


def test_split_record_equals_the_primitive_on_the_units_fit(simulation_unit) -> None:
    from motco.simulations.evaluation import build_simulation_trajectory_design, integrate_semisynthetic_dataset
    from motco.stats.interaction import interaction_split_test

    settings, record = simulation_unit
    dataset = lsc.unit_dataset(UNIT, settings["comparison"])
    params = ssc.production_params(settings, seed=record["rrpp_seed"], n_stages=3)
    latent = integrate_semisynthetic_dataset(dataset, params)
    model = latent.artifacts.model
    features = (latent.artifacts.features.to_numpy(dtype=float) - model._x_mean) / model._x_std
    basis, _ = np.linalg.qr(np.asarray(model.x_rotations_, dtype=float))
    design = build_simulation_trajectory_design(dataset.metadata, group_col="group", stage_col="stage")
    for correction, field in ((None, "p_value_ss"), ("df", "p_value_ss_df")):
        expected = interaction_split_test(
            features, design.model_full, design.model_reduced, basis, permutations=39, seed=record["rrpp_seed"],
            residual_correction=correction,
        )
        for part in ssc.PARTS:
            entry = record["split"]["parts"][part]
            assert entry[field] == expected.parts[part].p_value_ss
            if correction is None:
                assert entry["ss_interaction"] == expected.parts[part].ss_interaction
                assert entry["rss"] == expected.parts[part].rss
                assert entry["f"] == expected.parts[part].f
                assert entry["p_value_f"] == expected.parts[part].p_value_f


def test_a_failing_unit_is_recorded_without_aborting_the_shard(tmp_path, monkeypatch) -> None:
    settings = _small_settings()

    def broken(unit, settings):
        raise RuntimeError(f"no data for {unit.cell}")

    monkeypatch.setattr(lsc, "unit_dataset", broken)
    evaluated = ssc.run_shard(tmp_path, settings, 0, 1, limit=2)
    assert evaluated == 2
    lines = (tmp_path / "simulation" / "records" / "shard_000_of_001.jsonl").read_text().splitlines()
    records = [json.loads(line) for line in lines]
    assert [r["status"] for r in records] == ["failed", "failed"]
    assert records[0]["error_type"] == "RuntimeError" and "no data" in records[0]["error_message"]


def _fake_simulation_record(unit, settings, *, p=0.5):
    return {
        "profile": "simulation",
        "design": unit.design,
        "cell": unit.cell,
        "trajectory_mode": lsc.cell_spec(settings["comparison"], unit.cell)["mode"],
        "effect_size": 0.0,
        "span_alignment": None,
        "span_support": None,
        "replicate": unit.replicate,
        "settings_hash": lsc.settings_hash(settings),
        "status": "completed",
        "production": {
            "rank": 2,
            "cv_mean_auroc": 0.9,
            "p_values": {"delta": p, "angle": p, "shape": p},
            "p_values_reproduction": {"delta": p, "angle": p, "shape": p},
            "span_capture": {"group_difference": 0.5, "pooled_stage": 0.9},
        },
        "split": {
            "parts": {part: {"p_value_f": p, "p_value_ss": p, "p_value_ss_df": p} for part in ssc.PARTS},
            "span_share": 0.25,
        },
        "runtime_seconds": {"unit": 1.0, "fit": 0.5, "evaluation": 0.3, "split": 0.1},
    }


def test_resume_skips_completed_units(tmp_path, monkeypatch) -> None:
    settings = _small_settings()
    calls: list = []

    def fake(unit, settings):
        calls.append(unit.key)
        return _fake_simulation_record(unit, settings)

    monkeypatch.setattr(ssc, "evaluate_simulation_unit", fake)
    assert ssc.run_shard(tmp_path, settings, 0, 2, limit=3) == 3
    with (tmp_path / "simulation" / "records" / "shard_000_of_002.jsonl").open("a") as handle:
        handle.write('{"torn": ')  # an interrupted write
    first = list(calls)
    total = len(ssc.simulation_units(settings))
    assert ssc.run_shard(tmp_path, settings, 0, 2) == total // 2 - 3
    assert not set(first) & set(calls[3:])
    ssc.run_shard(tmp_path, settings, 1, 2)
    summary = ssc.merge(tmp_path, "simulation")
    assert summary["complete_units"] == total and summary["missing_units"] == 0
    merged = (tmp_path / "simulation" / "merged.jsonl").read_text().splitlines()
    assert len(merged) == total


def test_settings_mismatch_is_refused_naming_the_keys(tmp_path, monkeypatch) -> None:
    monkeypatch.setattr(ssc, "evaluate_simulation_unit", _fake_simulation_record)
    ssc.run_shard(tmp_path, _small_settings(), 0, 1, limit=1)
    with pytest.raises(SystemExit, match="differing: .*permutations"):
        ssc.run_shard(tmp_path, _small_settings(permutations=49), 0, 1, limit=1)


# ── Simulation merge ─────────────────────────────────────────────────────────


def test_simulation_tables_count_failures_by_cell() -> None:
    settings = _small_settings()
    units = ssc.simulation_units(settings)
    records = [_fake_simulation_record(u, settings, p=0.01 if u.replicate == 0 else 0.5) for u in units]
    records[0] = {**records[0], "status": "failed", "production": None, "split": None}
    tables = ssc.simulation_tables(records, settings)
    assert set(tables["split"][0]) == set(ssc.SIMULATION_SPLIT_COLUMNS)
    assert len(tables["split"]) == 14 * 3 * 3
    first = [r for r in tables["split"] if (r["cell"], r["part"], r["statistic"]) == (records[0]["cell"], "total", "f")]
    assert first[0]["failures"] == 1 and first[0]["count"] == 1
    other = [r for r in tables["split"] if (r["cell"], r["part"], r["statistic"]) == ("translation", "span", "ss")]
    assert (other[0]["count"], other[0]["rejections"], other[0]["rate"]) == (2, 1, 0.5)
    assert other[0]["span_share_median"] == 0.25
    production = [r for r in tables["production"] if r["cell"] == "translation" and r["statistic"] == "angle"]
    assert production[0]["rejections"] == 1
    capture = [r for r in tables["capture"] if r["cell"] == "translation"][0]
    assert capture["group_difference_median"] == 0.5 and capture["rank_median"] == 2


def _committed_rows(records, settings, *, bump=None):
    pseudo = [{"design": r["design"], "cell": r["cell"], "arm": "pls_cv", "status": r["status"],
               "p_values": r["production"]["p_values_reproduction"]} for r in records]
    rows = [{**row, "count": str(row["count"]), "rejections": str(row["rejections"])}
            for row in lsc.rate_table(pseudo, alpha=0.05)]
    if bump is not None:
        rows[bump]["rejections"] = str(int(rows[bump]["rejections"]) + 1)
    return rows + [{"design": "seaad", "cell": "anchor", "arm": "concat", "statistic": "angle", "count": "9",
                    "rejections": "9"}]


def test_simulation_reproduction_passes_and_fails() -> None:
    settings = _small_settings()
    records = [_fake_simulation_record(u, settings, p=0.01 if u.replicate else 0.5)
               for u in ssc.simulation_units(settings)]
    passed = ssc.simulation_reproduction(records, settings, _committed_rows(records, settings))
    assert passed["passed"] and passed["rows_compared"] == 14 * 3 and not passed["differing"]
    failed = ssc.simulation_reproduction(records, settings, _committed_rows(records, settings, bump=4))
    assert not failed["passed"]
    assert len(failed["differing"]) == 1
    assert failed["differing"][0]["committed"]["rejections"] == failed["differing"][0]["reproduced"]["rejections"] + 1


def test_reproduction_uses_only_the_shared_replicates() -> None:
    settings = _small_settings()
    settings["parent"] = {**settings["parent"], "shared_replicates": 1}
    records = [_fake_simulation_record(u, settings) for u in ssc.simulation_units(settings)]
    result = ssc.simulation_reproduction(records, settings, _committed_rows(records[:0], settings))
    assert not result["passed"]  # no committed pls_cv row exists for an empty run
    shared = [r for r in records if r["replicate"] < 1]
    assert ssc.simulation_reproduction(records, settings, _committed_rows(shared, settings))["passed"]


# ── Plasmode profile ─────────────────────────────────────────────────────────


def _plasmode_tests():
    sys.path.insert(0, str(REPO))
    try:
        from tests import test_plasmode_check as module
    finally:
        sys.path.pop(0)
    return module


plasmode_tests = _plasmode_tests()


@pytest.fixture(scope="module")
def stand_in(tmp_path_factory):
    directory = plasmode_tests.write_export(tmp_path_factory.mktemp("export"))
    export = pc.load_export(directory)
    settings = {**pc.base_settings(export), "cv": dict(plasmode_tests.SMALL_CV), "attribution": None}
    built = pc.build_unit(0.0, 0, {**settings, "stream": "fidelity"}, export)
    _, result = pc.analyse(built["blocks"], built["metadata"], export, pc.analysis_params(settings, None,
                                                                                         permutations=0))
    latent = result.record["latent"]
    fidelity = {"selected_lv": int(latent["selected_lv"]), "cv_mean_auroc": round(float(latent["cv_mean_auroc"]), 6)}
    return directory, export, fidelity


def _write_parent(directory: Path, settings: dict, summary_rows: list[dict]) -> Path:
    directory.mkdir(parents=True, exist_ok=True)
    (directory / "settings.json").write_text(json.dumps({"hash": lsc.settings_hash(settings), "settings": settings}))
    columns = ("effect", "cell", "arm", "block", "component", "subset", "quantity", "n", "mean", "se", "role")
    lsc._write_csv(directory / "plasmode_summary.csv", summary_rows, columns)
    return directory


@pytest.fixture(scope="module")
def plasmode_parent(stand_in, tmp_path_factory):
    _, export, fidelity = stand_in
    parent = pc.build_settings(export, fidelity, replicates=2, permutations=9, bootstrap=5, e_star=3.0,
                               descriptive=2.0)
    parent["cv"] = dict(plasmode_tests.SMALL_CV)
    return parent, _write_parent(tmp_path_factory.mktemp("plasmode"), parent, [])


def _plasmode_settings(stand_in, plasmode_parent, **kwargs):
    _, export, fidelity = stand_in
    _, parent_dir = plasmode_parent
    kwargs.setdefault("replicates", 2)
    kwargs.setdefault("permutations", 9)
    return ssc.plasmode_settings(export, fidelity, parent_dir=parent_dir, **kwargs)


def test_plasmode_units_are_replicate_major(stand_in, plasmode_parent) -> None:
    settings = _plasmode_settings(stand_in, plasmode_parent)
    assert ssc.plasmode_units(settings) == [(0.0, 0), (2.0, 0), (3.0, 0), (0.0, 1), (2.0, 1), (3.0, 1)]
    assert settings["rejection"] == "p <= alpha" and settings["span_split_test"] is True
    assert settings["attribution"] is False
    assert ssc.plasmode_units(_plasmode_settings(stand_in, plasmode_parent, replicates=500))[-1] == (3.0, 499)


def test_an_export_whose_hashes_differ_is_refused(stand_in, plasmode_parent, tmp_path) -> None:
    directory, _, fidelity = stand_in
    copy = shutil.copytree(directory, tmp_path / "export")
    frame = (copy / "atac.csv").read_text().splitlines()
    (copy / "atac.csv").write_text("\n".join(frame[:-1] + [frame[-1].replace(",", ",1", 1)]) + "\n")
    other = pc.load_export(copy)
    with pytest.raises(SystemExit, match=r"differ from .*settings.json for: atac.csv; refusing"):
        ssc.plasmode_settings(other, fidelity, parent_dir=plasmode_parent[1], replicates=2, permutations=9)


def test_a_cell_the_plasmode_did_not_run_is_refused(stand_in, plasmode_parent) -> None:
    with pytest.raises(SystemExit, match=r"Cell\(s\) \[1.5\]"):
        _plasmode_settings(stand_in, plasmode_parent, cells=(0.0, 1.5))


@pytest.fixture(scope="module")
def plasmode_record(stand_in, plasmode_parent):
    _, export, _ = stand_in
    settings = _plasmode_settings(stand_in, plasmode_parent)
    return settings, ssc.evaluate_plasmode_unit(3.0, 1, settings, export)


def test_plasmode_record_holds_seeds_and_statistics_only(stand_in, plasmode_record) -> None:
    _, export, _ = stand_in
    _, record = plasmode_record
    assert record["status"] == "completed", record.get("error_message")
    text = json.dumps(record)
    for donor in export.metadata[export.sample_col]:
        assert donor not in text
    assert "Female" not in text and "Male" not in text
    assert set(record) == {"profile", "key", "effect", "cell", "replicate", "stream", "seeds", "settings_hash",
                           "status", "label_redraws", "analysis", "span_split", "ss_uncorrected", "runtime_seconds"}
    assert record["seeds"] == pc.unit_seeds(2034, "paper", 1)


def test_plasmode_split_fields_equal_the_analysis_block(stand_in, plasmode_record) -> None:
    _, export, _ = stand_in
    settings, record = plasmode_record
    built = pc.build_unit(3.0, 1, settings["plasmode"], export)
    _, result = pc.analyse(built["blocks"], built["metadata"], export, ssc.plasmode_params(settings, built["seeds"]))
    assert record["span_split"] == json.loads(json.dumps(result.record["span_split"]))
    assert record["span_split"]["residual_correction"] == "df"
    assert set(record["ss_uncorrected"]) == set(ssc.PARTS)


def test_plasmode_trajectory_p_values_equal_the_plasmodes_own_unit(stand_in, plasmode_parent, plasmode_record):
    _, export, _ = stand_in
    parent, _ = plasmode_parent
    _, record = plasmode_record
    theirs = pc.evaluate_unit(3.0, 1, parent, export)
    assert record["analysis"]["p_values"] == theirs["analysis"]["p_values"]
    assert record["analysis"]["statistics"] == theirs["analysis"]["statistics"]


def test_the_anchor_unit_passes_the_latent_check(stand_in, plasmode_parent) -> None:
    _, export, _ = stand_in
    settings = _plasmode_settings(stand_in, plasmode_parent)
    record = ssc.evaluate_plasmode_unit(0.0, 0, settings, export)
    assert record["status"] == "completed", record.get("error_message")
    wrong = {**settings, "fidelity": {"selected_lv": 19, "cv_mean_auroc": 0.5}}
    failed = ssc.evaluate_plasmode_unit(0.0, 0, wrong, export)
    assert failed["status"] == "failed" and "anchor_latent" in failed["error_message"]


def _fake_plasmode_record(effect, replicate, settings, p=0.5):
    return {
        "profile": "plasmode",
        "key": pc.unit_key(effect, replicate),
        "effect": effect,
        "cell": pc.cell_label(effect),
        "replicate": replicate,
        "settings_hash": lsc.settings_hash(settings),
        "status": "completed",
        "analysis": {"p_values": {"delta": p, "angle": 0.05, "shape": p}},
        "span_split": {"parts": {part: {"p_value_f": p, "p_value_ss": 0.01} for part in ssc.PARTS},
                       "span_share": 0.1},
        "ss_uncorrected": {part: 0.001 for part in ssc.PARTS},
        "runtime_seconds": 2.0,
    }


def test_plasmode_tables_and_both_reproduction_outcomes(stand_in, plasmode_parent) -> None:
    settings = _plasmode_settings(stand_in, plasmode_parent)
    records = [_fake_plasmode_record(e, r, settings) for e, r in ssc.plasmode_units(settings)]
    tables = ssc.plasmode_tables(records, settings)
    assert set(tables["split"][0]) == set(ssc.PLASMODE_SPLIT_COLUMNS)
    by = {(r["cell"], r["part"], r["statistic"]): r for r in tables["split"]}
    assert by[("anchor", "total", "f")]["rate"] == 0.0
    assert by[("anchor", "total", "ss_df")]["rate"] == 1.0  # the analysis's SS p-value is the df one
    assert by[("e=3", "complement", "ss")]["rate"] == 1.0
    angle = [r for r in tables["trajectory"] if r["statistic"] == "angle"]
    assert all(r["rate"] == 1.0 for r in angle)  # p <= alpha: 0.05 rejects
    committed = [
        {"effect": str(e), "cell": pc.cell_label(e), "arm": "analysis", "subset": "all", "quantity": f"reject_{s}",
         "n": "2", "mean": "1.0" if s == "angle" else "0.0"}
        for e in settings["cells"] for s in ssc.STATISTICS
    ]
    assert ssc.plasmode_reproduction(records, settings, committed)["passed"]
    committed[1] = {**committed[1], "mean": "0.5"}
    failed = ssc.plasmode_reproduction(records, settings, committed)
    assert not failed["passed"]
    assert [(row["cell"], row["statistic"]) for row in failed["differing"]] == [("anchor", "angle")]


def test_plasmode_merge_writes_its_tables(stand_in, plasmode_parent, tmp_path, monkeypatch) -> None:
    _, export, _ = stand_in
    settings = _plasmode_settings(stand_in, plasmode_parent)
    monkeypatch.setattr(ssc, "evaluate_plasmode_unit", lambda e, r, s, x: _fake_plasmode_record(e, r, s))
    assert ssc.run_shard(tmp_path, settings, 0, 1, export=export) == 6
    summary = ssc.merge(tmp_path, "plasmode")
    assert summary["complete_units"] == 6 and summary["failed_units"] == 0
    assert summary["reproduction_passed"] is False  # the stand-in parent has no committed summary rows
    for name in ("split_table.csv", "trajectory_table.csv", "accounting.json", "reproduction.json", "merged.jsonl"):
        assert (tmp_path / "plasmode" / name).exists()


# ── Decision ─────────────────────────────────────────────────────────────────


def _row(rate: float, count: int, **cell) -> dict:
    rejections = round(rate * count)
    rate = rejections / count
    return {**cell, "rate": rate, "se": (rate * (1 - rate) / count) ** 0.5, "count": count, "rejections": rejections}


def _tables(*, size=0.05, complement=0.9, angle=0.1, plasmode_size=0.05, overrides=None):
    overrides = overrides or {}
    cells = [cell[0] for cell in lsc.CELLS]
    sim_settings = {"profile": "simulation", "permutations": 999, "rejection": "p < alpha",
                    "comparison": {"designs": {d: {} for d in lsc.DESIGNS}, "cells": [[c] for c in cells]}}
    split, production, capture = [], [], []
    for design in lsc.DESIGNS:
        for cell in cells:
            count = 500 if cell in ssc.CONTROL_CELLS else 100
            for part in ssc.PARTS:
                for statistic in ssc.SPLIT_STATISTICS:
                    rate = size if cell in ssc.CONTROL_CELLS else (complement if part == "complement" else 0.5)
                    rate = overrides.get((design, cell, part, statistic), rate)
                    split.append(_row(rate, count, design=design, cell=cell, part=part, statistic=statistic))
            for statistic in ssc.STATISTICS:
                rate = angle if statistic == "angle" else 0.05
                production.append(_row(rate, count, design=design, cell=cell, statistic=statistic))
            capture.append({"design": design, "cell": cell, "span_share_median": 0.3, "group_difference_median": 0.1,
                            "pooled_stage_median": 0.9, "rank_median": 2})
    pm_settings = {"profile": "plasmode", "permutations": 999, "rejection": "p <= alpha"}
    pm_split = [
        _row(overrides.get(("plasmode", label, part, statistic), plasmode_size if label == "anchor" else 0.2), 500,
             cell=label, part=part, statistic=statistic)
        for label in ("anchor", "e=2", "e=3") for part in ssc.PARTS for statistic in ssc.SPLIT_STATISTICS
    ]
    trajectory = [
        _row(0.05, 500, cell=label, statistic=s) for label in ("anchor", "e=2", "e=3") for s in ssc.STATISTICS
    ]
    reproduced = {"passed": True, "rows_compared": 42, "differing": [], "source": "x"}
    accounting = {"missing_units": 0, "failed_units": 0}
    simulation = {"settings": sim_settings, "split": split, "production": production, "capture": capture,
                  "reproduction": dict(reproduced), "accounting": accounting}
    plasmode = {"settings": pm_settings, "split": pm_split, "trajectory": trajectory,
                "reproduction": dict(reproduced), "accounting": accounting}
    return simulation, plasmode


def test_validated() -> None:
    result = ssc.decision(*_tables())
    assert result["verdict"] == "split_validated"
    assert result["criteria"]["C1"]["m"] == 18 and result["criteria"]["C2"]["m"] == 3
    assert result["criteria"]["C3"]["pass"]
    assert [check["design"] for check in result["criteria"]["C3"]["checks"]] == ["n200", "seaad"]
    assert "measurement-visibility" in result["consequence"]
    descriptive = result["descriptive"]
    assert set(descriptive) == {"localization", "share_against_capture", "span_against_production",
                                "complement_against_comparison", "plasmode_power", "ss_diagnostic"}
    assert all("source" in row and "grade" in row for row in descriptive["localization"])
    assert [row["cell"] for row in descriptive["plasmode_power"]] == ["e=2", "e=3"]


def test_simulation_size_excess_names_the_rates() -> None:
    overrides = {("seaad_weak", "anchor", "complement", "f"): 0.12, ("n200", "translation", "total", "f"): 0.10}
    result = ssc.decision(*_tables(overrides=overrides))
    assert result["verdict"] == "size_excess"
    failing = {(row["design"], row["cell"], row["part"]) for row in result["failing"]["C1"]}
    assert failing == {("seaad_weak", "anchor", "complement"), ("n200", "translation", "total")}


def test_size_excess_takes_precedence_over_a_failing_gain() -> None:
    result = ssc.decision(*_tables(complement=0.1, overrides={("plasmode", "anchor", "span", "f"): 0.09}))
    assert result["verdict"] == "size_excess"
    assert [row["part"] for row in result["failing"]["C2"]] == ["span"]
    assert result["failing"]["C3"] == ["n200", "seaad"]


def test_gain_at_one_design_only_fails_naming_it() -> None:
    result = ssc.decision(*_tables(overrides={("seaad", "orient_e0.5_phi90", "complement", "f"): 0.15}))
    assert result["verdict"] == "no_out_of_span_gain"
    assert result["failing"] == {"C3": ["seaad"]}
    check = result["criteria"]["C3"]["checks"][0]
    assert check["pass"] and check["difference"] == pytest.approx(0.8)


def test_reproduction_failure_gates_every_criterion() -> None:
    simulation, plasmode = _tables()
    simulation["reproduction"] = {"passed": False, "differing": [{"design": "n200"}], "rows_compared": 42}
    result = ssc.decision(simulation, plasmode)
    assert result["verdict"] == "reproduction_failed"
    assert "criteria" not in result and "descriptive" not in result
    assert result["failing"]["reproduction"][0]["differing"] == [{"design": "n200"}]


def test_descriptive_rates_do_not_change_the_verdict() -> None:
    first = ssc.decision(*_tables())
    second = ssc.decision(*_tables(overrides={("n200", "shape_e0.5", "span", "f"): 0.0,
                                              ("plasmode", "e=3", "complement", "f"): 1.0}))
    assert first["verdict"] == second["verdict"] == "split_validated"
    assert first["criteria"] == second["criteria"]
    assert first["descriptive"] != second["descriptive"]


def test_missing_gain_rate_fails_c3() -> None:
    simulation, plasmode = _tables()
    simulation["production"] = [r for r in simulation["production"] if r["design"] != "n200"]
    result = ssc.decision(simulation, plasmode)
    assert result["verdict"] == "no_out_of_span_gain"
    assert result["criteria"]["C3"]["checks"][0]["missing"] is True


# ── sbatch ───────────────────────────────────────────────────────────────────


def test_sbatch_template_parses() -> None:
    subprocess.run(["bash", "-n", str(REPO / "scripts" / "span_split_check_array.sbatch")], check=True)
