"""``scripts/interaction_scan_check.py``: settings binding, both units, resume, merge, reproduction and decision.

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


isc = _script("interaction_scan_check")
ssc = isc.ssc
lsc = isc.lsc
pc = isc.pc

FAST_CV = {"cv1_splits": 3, "cv2_splits": 4, "n_repeats": 1, "max_components": 6, "random_state": 1203}
REPO = Path(__file__).resolve().parents[1]


# ── Simulation settings and units ────────────────────────────────────────────


def test_simulation_units_equal_the_parents() -> None:
    settings = isc.simulation_settings()
    units = isc.simulation_units(settings)
    committed = json.loads((isc.SPAN_SPLIT_DIR / "simulation" / "settings.json").read_text())["settings"]
    assert units == ssc.simulation_units(committed)
    assert len(units) == 3 * (2 * 500 + 12 * 100)
    for design in lsc.DESIGNS:
        for cell, *_ in lsc.CELLS:
            replicates = sorted(u.replicate for u in units if (u.design, u.cell) == (design, cell))
            assert replicates == list(range(500 if cell in ("anchor", "translation") else 100))
    assert settings["span_split"] == committed
    assert settings["permutations"] == 999
    assert settings["rejection"] == {"scan": "p <= alpha", "parent": "p < alpha"}
    assert set(settings["parent"]["files"]) == {"settings.json", "split_table.csv", "production_table.csv"}
    assert settings["parent"]["files"]["settings.json"] == ssc.sha256_file(
        isc.SPAN_SPLIT_DIR / "simulation" / "settings.json"
    )


def _parent_copy(tmp_path: Path, profile: str = "simulation") -> Path:
    target = tmp_path / "span-split" / profile
    target.mkdir(parents=True)
    for name in isc.PARENT_FILES[profile]:
        shutil.copy(isc.SPAN_SPLIT_DIR / profile / name, target / name)
    return target


def test_a_changed_parent_settings_file_is_refused_by_name(tmp_path) -> None:
    parent = _parent_copy(tmp_path)
    stored = json.loads((parent / "settings.json").read_text())
    stored["settings"]["alpha"] = 0.1
    stored["hash"] = lsc.settings_hash(stored["settings"])
    (parent / "settings.json").write_text(json.dumps(stored))
    with pytest.raises(SystemExit, match=r"span-split/simulation/settings.json differs .*differing: alpha\."):
        isc.simulation_settings(parent_dir=parent)


def test_a_parent_file_whose_hash_differs_is_refused_by_name(tmp_path) -> None:
    parent = _parent_copy(tmp_path)
    settings = isc.simulation_settings(parent_dir=parent, designs=["seaad"], replicates=1, control_replicates=1,
                                       permutations=9)
    isc.check_parent_files(settings)
    with (parent / "split_table.csv").open("a") as handle:
        handle.write("\n")
    with pytest.raises(SystemExit, match=r"split_table.csv has changed since the settings were written"):
        isc.check_parent_files(settings)
    with pytest.raises(SystemExit, match=r"split_table.csv has changed"):
        isc.run_shard(tmp_path / "run", settings, 0, 1, limit=0)


def _small_settings(**overrides):
    settings = isc.simulation_settings(designs=["seaad"], replicates=2, control_replicates=2, permutations=39)
    settings["span_split"]["comparison"]["pls_cv_params"] = dict(FAST_CV)
    settings["span_split"]["reproduction_permutations"] = 19
    settings.update(overrides)
    return settings


UNIT = lsc.Unit("seaad", "orient_e0.5_phi90_k20", 1)


@pytest.fixture(scope="module")
def simulation_unit():
    settings = _small_settings()
    return settings, isc.evaluate_simulation_unit(UNIT, settings)


def test_simulation_record_carries_every_field(simulation_unit) -> None:
    _, record = simulation_unit
    assert record["status"] == "completed", record.get("error_message")
    assert (record["design"], record["cell"], record["replicate"]) == ("seaad", "orient_e0.5_phi90_k20", 1)
    assert record["rrpp_seed"] == lsc.rrpp_seed(2035, "seaad", "orient_e0.5_phi90_k20", 1)
    scan = record["scan"]
    assert {"n_features", "n_degenerate", "max_f", "max_block", "p_value_global", "hits", "blocks",
            "n_population_null", "hits_population_null"} <= set(scan)
    assert scan["n_features"] == 367 + 131
    assert set(scan["blocks"]) == {"methylation", "expression"}
    assert sum(entry["null_max_share"] for entry in scan["blocks"].values()) == pytest.approx(1.0)
    assert sum(entry["hits"] for entry in scan["blocks"].values()) == scan["hits"]
    assert set(record["runtime_seconds"]) == {"fit", "evaluation", "split", "scan", "unit"}


def test_production_and_split_fields_equal_the_parents_unit(simulation_unit) -> None:
    settings, record = simulation_unit
    theirs = ssc.evaluate_simulation_unit(UNIT, settings["span_split"])
    assert record["production"] == json.loads(json.dumps(lsc._jsonable(theirs["production"])))
    assert record["split"] == theirs["split"]
    for key in ("design", "cell", "replicate", "generator_seed", "rrpp_seed", "trajectory_mode"):
        assert record[key] == theirs[key]


def test_scan_fields_equal_the_primitive_on_the_split_inputs(simulation_unit) -> None:
    from motco.simulations.evaluation import integrate_semisynthetic_dataset
    from motco.stats.interaction import interaction_feature_scan

    settings, record = simulation_unit
    dataset = lsc.unit_dataset(UNIT, settings["span_split"]["comparison"])
    params = ssc.production_params(settings["span_split"], seed=record["rrpp_seed"], n_stages=3)
    latent = integrate_semisynthetic_dataset(dataset, params)
    features, _, full, reduced = ssc.split_inputs(dataset, latent, params)
    scan = interaction_feature_scan(features, full, reduced, permutations=39, seed=record["rrpp_seed"])
    columns = list(latent.artifacts.features.columns)
    assert record["scan"]["max_f"] == scan.max_f
    assert record["scan"]["max_block"] == columns[scan.max_column].split("__")[0]
    assert record["scan"]["p_value_global"] == scan.p_value_global
    assert record["scan"]["hits"] == len(scan.hits())


def test_population_null_count_matches_a_direct_computation(simulation_unit) -> None:
    settings, record = simulation_unit
    dataset = lsc.unit_dataset(UNIT, settings["span_split"]["comparison"])
    population = dataset.population_trajectories
    columns = [f"methylation__{c}" for c in dataset.methylation.columns] + [
        f"expression__{c}" for c in dataset.expression.columns
    ]
    expected = 0
    for column in columns:
        layer, feature = column.split("__", 1)
        values = population.layers[layer][feature]
        difference = [values.loc[("B", stage)] - values.loc[("A", stage)] for stage in ("0", "1", "2")]
        expected += (max(difference) - min(difference)) <= 1e-9 * float(np.abs(values).max())
    assert record["scan"]["n_population_null"] == expected
    assert 0 < expected < len(columns)  # a sparse out-of-span plant leaves most out-of-span features null
    assert 0 <= record["scan"]["hits_population_null"] <= record["scan"]["hits"]


def test_population_null_covers_every_feature_under_translation_and_the_anchor() -> None:
    settings = _small_settings()
    for cell in ("anchor", "translation"):
        dataset = lsc.unit_dataset(lsc.Unit("seaad", cell, 0), settings["span_split"]["comparison"])
        columns = [f"methylation__{c}" for c in dataset.methylation.columns]
        mask = isc.population_null_mask(dataset, columns, "group", 1e-9)
        assert mask is not None and mask.all()


def test_a_failing_unit_is_recorded_without_aborting_the_shard(tmp_path, monkeypatch) -> None:
    settings = _small_settings()

    def broken(unit, settings):
        raise RuntimeError(f"no data for {unit.cell}")

    monkeypatch.setattr(lsc, "unit_dataset", broken)
    assert isc.run_shard(tmp_path, settings, 0, 1, limit=2) == 2
    lines = (tmp_path / "simulation" / "records" / "shard_000_of_001.jsonl").read_text().splitlines()
    records = [json.loads(line) for line in lines]
    assert [r["status"] for r in records] == ["failed", "failed"]
    assert records[0]["error_type"] == "RuntimeError"


def _fake_simulation_record(unit, settings, *, p=0.5, hits=0, null_hits=0, n_null=10):
    blocks = {"methylation": {"n_features": 367, "n_degenerate": 0, "hits": hits, "null_max_share": 0.7},
              "expression": {"n_features": 131, "n_degenerate": 0, "hits": 0, "null_max_share": 0.3}}
    return {
        "profile": "simulation",
        "design": unit.design,
        "cell": unit.cell,
        "trajectory_mode": lsc.cell_spec(settings["span_split"]["comparison"], unit.cell)["mode"],
        "effect_size": 0.0,
        "span_alignment": None,
        "span_support": None,
        "replicate": unit.replicate,
        "settings_hash": lsc.settings_hash(settings),
        "status": "completed",
        "production": {"rank": 2, "p_values": {"delta": p, "angle": p, "shape": p}},
        "split": {"parts": {part: {"p_value_f": p} for part in isc.PARTS}, "span_share": 0.25},
        "scan": {"n_features": 498, "n_degenerate": 0, "max_f": 9.0, "max_block": "methylation",
                 "p_value_global": p, "hits": hits, "blocks": blocks, "n_population_null": n_null,
                 "hits_population_null": null_hits},
        "runtime_seconds": {"unit": 2.0, "fit": 1.0, "evaluation": 0.5, "split": 0.2, "scan": 0.1},
    }


def test_resume_skips_completed_units(tmp_path, monkeypatch) -> None:
    settings = _small_settings()
    calls: list = []

    def fake(unit, settings):
        calls.append(unit.key)
        return _fake_simulation_record(unit, settings)

    monkeypatch.setattr(isc, "evaluate_simulation_unit", fake)
    assert isc.run_shard(tmp_path, settings, 0, 2, limit=3) == 3
    with (tmp_path / "simulation" / "records" / "shard_000_of_002.jsonl").open("a") as handle:
        handle.write('{"torn": ')  # an interrupted write
    first = list(calls)
    total = len(isc.simulation_units(settings))
    assert isc.run_shard(tmp_path, settings, 0, 2) == total // 2 - 3
    assert not set(first) & set(calls[3:])
    isc.run_shard(tmp_path, settings, 1, 2)
    summary = isc.merge(tmp_path, "simulation")
    assert summary["complete_units"] == total and summary["missing_units"] == 0
    assert len((tmp_path / "simulation" / "merged.jsonl").read_text().splitlines()) == total
    for name in ("scan_table.csv", "block_balance.csv", "comparator_table.csv", "false_hit_table.csv",
                 "accounting.json", "reproduction.json"):
        assert (tmp_path / "simulation" / name).exists()
    accounting = json.loads((tmp_path / "simulation" / "accounting.json").read_text())
    assert accounting["scan_share_of_unit_seconds"] == pytest.approx(0.05)


def test_settings_mismatch_is_refused_naming_the_keys(tmp_path, monkeypatch) -> None:
    monkeypatch.setattr(isc, "evaluate_simulation_unit", _fake_simulation_record)
    isc.run_shard(tmp_path, _small_settings(), 0, 1, limit=1)
    with pytest.raises(SystemExit, match="differing: .*permutations"):
        isc.run_shard(tmp_path, _small_settings(permutations=49), 0, 1, limit=1)


# ── Simulation merge ─────────────────────────────────────────────────────────


def _simulation_records(settings):
    records = []
    for unit in isc.simulation_units(settings):
        rejects = unit.replicate == 0
        records.append(_fake_simulation_record(unit, settings, p=0.01 if rejects else 0.5, hits=3 if rejects else 0,
                                               null_hits=1 if rejects else 0))
    return records


def test_simulation_tables_count_failures_by_cell() -> None:
    settings = _small_settings()
    records = _simulation_records(settings)
    records[0] = {**records[0], "status": "failed", "production": None, "split": None, "scan": None}
    scan = isc.scan_rows(records, settings)
    assert set(scan[0]) == set(isc.SCAN_COLUMNS["simulation"]) and len(scan) == 14
    anchor = [r for r in scan if r["cell"] == records[0]["cell"]][0]
    assert (anchor["failures"], anchor["count"], anchor["rejections"]) == (1, 1, 0)
    translation = [r for r in scan if r["cell"] == "translation"][0]
    assert (translation["count"], translation["rejections"], translation["rate"]) == (2, 1, 0.5)
    assert translation["hits_median"] == 1.5
    blocks = isc.block_rows(records, settings)
    row = [r for r in blocks if (r["cell"], r["block"]) == ("translation", "methylation")][0]
    assert row["null_max_share_mean"] == pytest.approx(0.7) and row["feature_share"] == pytest.approx(367 / 498)
    false_hits = [r for r in isc.false_hit_rows(records, settings) if r["cell"] == "translation"][0]
    assert (false_hits["count"], false_hits["rejections"]) == (2, 1)
    comparators = isc.comparator_rows(records, settings)
    assert {r["test"] for r in comparators} == set(isc.COMPARATOR_TESTS)
    rules = {r["test"]: r["rule"] for r in comparators}
    assert rules["scan"] == "p <= alpha" and rules["angle"] == "p < alpha"


def _committed(records, settings, *, bump=None):
    comparators = isc.comparator_rows(records, settings)
    split = [{"design": r["design"], "cell": r["cell"], "part": r["test"].split("_", 1)[1], "statistic": "f",
              "count": str(r["count"]), "rejections": str(r["rejections"])}
             for r in comparators if r["test"].startswith("split_")]
    split += [{**row, "statistic": "ss", "rejections": "999"} for row in split[:2]]  # never compared
    trajectory = [{"design": r["design"], "cell": r["cell"], "statistic": r["test"], "count": str(r["count"]),
                   "rejections": str(r["rejections"])} for r in comparators if r["test"] in isc.STATISTICS]
    if bump is not None:
        trajectory[bump]["rejections"] = str(int(trajectory[bump]["rejections"]) + 1)
    return comparators, split, trajectory


def test_simulation_reproduction_passes_and_fails() -> None:
    settings = _small_settings()
    records = _simulation_records(settings)
    comparators, split, trajectory = _committed(records, settings)
    passed = isc.reproduction(comparators, settings, split, trajectory)
    assert passed["passed"] and passed["rows_compared"] == 14 * 6 and not passed["differing"]
    comparators, split, trajectory = _committed(records, settings, bump=4)
    failed = isc.reproduction(comparators, settings, split, trajectory)
    assert not failed["passed"] and len(failed["differing"]) == 1
    entry = failed["differing"][0]
    assert entry["committed"]["rejections"] == entry["reproduced"]["rejections"] + 1
    assert entry["cell"][0] == "seaad" and entry["test"] in isc.STATISTICS


def test_reproduction_reads_only_the_runs_designs() -> None:
    settings = _small_settings()
    records = _simulation_records(settings)
    comparators, split, trajectory = _committed(records, settings)
    other = [{**row, "design": "n200", "rejections": "77"} for row in trajectory]
    assert isc.reproduction(comparators, settings, split, trajectory + other)["passed"]


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
    _, result = pc.analyse(built["blocks"], built["metadata"], export,
                           pc.analysis_params(settings, None, permutations=0))
    latent = result.record["latent"]
    fidelity = {"selected_lv": int(latent["selected_lv"]), "cv_mean_auroc": round(float(latent["cv_mean_auroc"]), 6)}
    return directory, export, fidelity


def _write(directory: Path, settings: dict, tables: dict[str, tuple[list[dict], tuple]]) -> Path:
    directory.mkdir(parents=True, exist_ok=True)
    (directory / "settings.json").write_text(json.dumps({"hash": lsc.settings_hash(settings), "settings": settings}))
    for name, (rows, columns) in tables.items():
        lsc._write_csv(directory / name, rows, columns)
    return directory


@pytest.fixture(scope="module")
def parents(stand_in, tmp_path_factory):
    """A stand-in plasmode run and a stand-in split-test plasmode run at 2 × 9."""

    _, export, fidelity = stand_in
    plasmode = pc.build_settings(export, fidelity, replicates=2, permutations=9, bootstrap=5, e_star=3.0,
                                 descriptive=2.0)
    plasmode["cv"] = dict(plasmode_tests.SMALL_CV)
    summary_columns = ("effect", "cell", "arm", "block", "component", "subset", "quantity", "n", "mean", "se", "role")
    summary = [{"effect": 3.0, "cell": "e=3", "arm": "adjusted", "block": "principal", "component": "residual",
                "subset": "all", "quantity": "planted_precision", "n": 2, "mean": 0.4, "se": 0.1, "role": ""}]
    plasmode_dir = _write(tmp_path_factory.mktemp("plasmode"), plasmode,
                          {"plasmode_summary.csv": (summary, summary_columns)})
    span_split = ssc.plasmode_settings(export, fidelity, parent_dir=plasmode_dir, replicates=2, permutations=9)
    parent_dir = _write(
        tmp_path_factory.mktemp("span-split") / "plasmode",
        span_split,
        {"split_table.csv": ([], ssc.PLASMODE_SPLIT_COLUMNS), "trajectory_table.csv": ([], ssc.TRAJECTORY_COLUMNS)},
    )
    return plasmode, plasmode_dir, parent_dir


def _plasmode_settings(stand_in, parents, **kwargs):
    _, export, fidelity = stand_in
    _, plasmode_dir, parent_dir = parents
    kwargs.setdefault("replicates", 2)
    kwargs.setdefault("permutations", 9)
    return isc.plasmode_settings(export, fidelity, parent_dir=parent_dir, plasmode_dir=plasmode_dir, **kwargs)


def test_plasmode_units_and_settings(stand_in, parents) -> None:
    settings = _plasmode_settings(stand_in, parents)
    assert isc.plasmode_units(settings) == [(0.0, 0), (2.0, 0), (3.0, 0), (0.0, 1), (2.0, 1), (3.0, 1)]
    assert settings["rejection"] == {"scan": "p <= alpha", "parent": "p <= alpha"}
    assert settings["interaction_scan"] is True and settings["span_split"]["attribution"] is False
    assert set(settings["parent"]["files"]) == {"settings.json", "split_table.csv", "trajectory_table.csv"}
    assert isc.plasmode_units(_plasmode_settings(stand_in, parents, replicates=500))[-1] == (3.0, 499)


def test_an_export_whose_hashes_differ_is_refused(stand_in, parents, tmp_path) -> None:
    directory, _, fidelity = stand_in
    _, plasmode_dir, parent_dir = parents
    copy = shutil.copytree(directory, tmp_path / "export")
    frame = (copy / "atac.csv").read_text().splitlines()
    (copy / "atac.csv").write_text("\n".join(frame[:-1] + [frame[-1].replace(",", ",1", 1)]) + "\n")
    other = pc.load_export(copy)
    with pytest.raises(SystemExit, match=r"differ from .*settings.json for: atac.csv; refusing"):
        isc.plasmode_settings(other, fidelity, parent_dir=parent_dir, plasmode_dir=plasmode_dir, replicates=2,
                              permutations=9)


def test_the_committed_plasmode_parent_must_match(stand_in, parents) -> None:
    _, export, fidelity = stand_in
    _, plasmode_dir, _ = parents
    # The real split-test run was made on the real export, not the stand-in.
    with pytest.raises(SystemExit, match=r"plasmode/settings.json differs .*differing: fidelity.cv_mean_auroc"):
        isc.plasmode_settings(export, fidelity, plasmode_dir=plasmode_dir)


@pytest.fixture(scope="module")
def plasmode_record(stand_in, parents):
    _, export, _ = stand_in
    settings = _plasmode_settings(stand_in, parents)
    return settings, isc.evaluate_plasmode_unit(3.0, 1, settings, export)


def test_plasmode_record_holds_seeds_statistics_and_counts_only(stand_in, plasmode_record) -> None:
    _, export, _ = stand_in
    _, record = plasmode_record
    assert record["status"] == "completed", record.get("error_message")
    text = json.dumps(record)
    for donor in export.metadata[export.sample_col]:
        assert donor not in text
    assert "Female" not in text and "Male" not in text
    for name, frame in export.blocks.items():
        for feature in frame.columns:
            assert f'"{feature}"' not in text
    assert set(record) == {"profile", "key", "effect", "cell", "replicate", "stream", "seeds", "settings_hash",
                           "status", "label_redraws", "analysis", "span_split", "scan", "planted", "runtime_seconds"}
    assert record["planted"]["planted"] == 20 and record["planted"]["top_k"] == 20
    assert record["planted"]["hits"] == record["scan"]["hits"]
    assert record["planted"]["hits_in_planted"] + record["planted"]["hits_outside_planted"] == record["scan"]["hits"]


def test_plasmode_scan_fields_equal_the_analysis_block(stand_in, plasmode_record) -> None:
    _, export, _ = stand_in
    settings, record = plasmode_record
    built = pc.build_unit(3.0, 1, settings["span_split"]["plasmode"], export)
    _, result = pc.analyse(built["blocks"], built["metadata"], export, isc.plasmode_params(settings, built["seeds"]))
    block = json.loads(json.dumps(result.record["interaction_scan"]))
    for key in ("n_features", "n_degenerate", "df_interaction", "df_residual", "max_f", "max_block",
                "p_value_global", "blocks", "permutations", "seed", "residual_correction", "draws_match_rrpp"):
        assert record["scan"][key] == block[key]
    assert record["scan"]["hits"] == len(block["hits"])
    assert record["scan"]["residual_correction"] == "df"
    expected = isc.planted_counts(result.interaction_scan, built["plant"])
    assert record["planted"] == json.loads(json.dumps(expected))


def test_plasmode_trajectory_and_split_equal_the_parents_unit(stand_in, plasmode_record) -> None:
    _, export, _ = stand_in
    settings, record = plasmode_record
    theirs = ssc.evaluate_plasmode_unit(3.0, 1, settings["span_split"], export)
    assert record["analysis"] == theirs["analysis"]
    assert record["span_split"] == theirs["span_split"]


def test_the_anchor_unit_passes_the_latent_check(stand_in, parents) -> None:
    _, export, _ = stand_in
    settings = _plasmode_settings(stand_in, parents)
    assert isc.evaluate_plasmode_unit(0.0, 0, settings, export)["status"] == "completed"
    wrong = json.loads(json.dumps(settings))
    wrong["span_split"]["fidelity"] = {"selected_lv": 19, "cv_mean_auroc": 0.5}
    failed = isc.evaluate_plasmode_unit(0.0, 0, wrong, export)
    assert failed["status"] == "failed" and "anchor_latent" in failed["error_message"]


def test_planted_counts_on_a_hand_built_table() -> None:
    import pandas as pd

    plant = {"rna": {"features": ["g1", "g2"]}, "atac": {"features": ["p0"]}}
    table = pd.DataFrame(
        {
            "block": ["rna", "rna", "rna", "atac", "atac"],
            "feature": ["g0", "g1", "g2", "p0", "p1"],
            "f": [9.0, 8.0, 1.0, np.nan, 7.0],
            "p_value_adjusted": [0.01, 0.02, 0.6, np.nan, 0.04],
            "degenerate": [False, False, False, True, False],
        }
    )
    counts = isc.planted_counts(table, plant, top_k=2)
    assert (counts["hits"], counts["hits_in_planted"], counts["hits_outside_planted"]) == (3, 1, 2)
    assert counts["top_k_in_planted"] == 1  # the top 2 by F are g0 and g1
    assert counts["blocks"]["atac"] == {"hits": 1, "hits_in_planted": 0, "hits_outside_planted": 1,
                                        "top_k_in_planted": 0, "planted": 1}


def _fake_plasmode_record(effect, replicate, settings, *, p=0.5, hits=0, inside=0, top=0):
    return {
        "profile": "plasmode",
        "key": pc.unit_key(effect, replicate),
        "effect": effect,
        "cell": pc.cell_label(effect),
        "replicate": replicate,
        "settings_hash": lsc.settings_hash(settings),
        "status": "completed",
        "analysis": {"p_values": {"delta": 0.5, "angle": 0.05, "shape": 0.5}},
        "span_split": {"parts": {part: {"p_value_f": 0.5} for part in isc.PARTS}},
        "scan": {"n_features": 100, "max_f": 12.0, "p_value_global": p, "hits": hits,
                 "blocks": {"rna": {"n_features": 60, "n_degenerate": 0, "hits": hits, "null_max_share": 0.6},
                            "atac": {"n_features": 40, "n_degenerate": 0, "hits": 0, "null_max_share": 0.4}}},
        "planted": {"top_k": 20, "planted": 20, "hits": hits, "hits_in_planted": inside,
                    "hits_outside_planted": hits - inside, "top_k_in_planted": top, "blocks": {}},
        "runtime_seconds": 80.0,
    }


def test_recovery_metrics_on_hand_built_counts(stand_in, parents) -> None:
    settings = _plasmode_settings(stand_in, parents)
    records = [
        _fake_plasmode_record(3.0, 0, settings, p=0.001, hits=5, inside=4, top=10),
        _fake_plasmode_record(3.0, 1, settings, p=0.2, hits=0, inside=0, top=6),
        _fake_plasmode_record(2.0, 0, settings, p=0.01, hits=2, inside=2, top=4),
        {**_fake_plasmode_record(2.0, 1, settings), "status": "failed"},
    ]
    rows = {row["cell"]: row for row in isc.recovery_rows(records, settings)}
    assert set(rows) == {"e=2", "e=3"}
    e3 = rows["e=3"]
    assert (e3["false_hit_count"], e3["false_hit_rejections"], e3["false_hit_rate"]) == (2, 1, 0.5)
    assert e3["recall_mean"] == pytest.approx((4 / 20 + 0) / 2)
    assert (e3["hit_precision_mean"], e3["hit_precision_n"]) == (pytest.approx(0.8), 1)
    assert e3["top_k_precision_mean"] == pytest.approx((10 / 20 + 6 / 20) / 2)
    assert (rows["e=2"]["failures"], rows["e=2"]["false_hit_rejections"]) == (1, 0)
    scan = {row["cell"]: row for row in isc.scan_rows(records, settings)}
    assert (scan["e=3"]["rejections"], scan["e=3"]["count"]) == (1, 2)
    assert scan["anchor"]["count"] == 0


def test_plasmode_reproduction_both_outcomes(stand_in, parents) -> None:
    settings = _plasmode_settings(stand_in, parents)
    records = [_fake_plasmode_record(e, r, settings) for e, r in isc.plasmode_units(settings)]
    comparators = isc.comparator_rows(records, settings)
    angle = {r["cell"]: r for r in comparators if r["test"] == "angle"}
    assert all(row["rate"] == 1.0 for row in angle.values())  # p <= alpha: 0.05 rejects
    split = [{"cell": r["cell"], "part": r["test"].split("_", 1)[1], "statistic": "f", "count": r["count"],
              "rejections": r["rejections"]} for r in comparators if r["test"].startswith("split_")]
    trajectory = [{"cell": r["cell"], "statistic": r["test"], "count": r["count"], "rejections": r["rejections"]}
                  for r in comparators if r["test"] in isc.STATISTICS]
    passed = isc.reproduction(comparators, settings, split, trajectory)
    assert passed["passed"] and passed["rows_compared"] == 3 * 6
    trajectory[1] = {**trajectory[1], "rejections": 1}
    failed = isc.reproduction(comparators, settings, split, trajectory)
    assert not failed["passed"]
    assert [(row["cell"], row["test"]) for row in failed["differing"]] == [(["anchor"], "angle")]


def test_plasmode_merge_writes_its_tables(stand_in, parents, tmp_path, monkeypatch) -> None:
    _, export, _ = stand_in
    settings = _plasmode_settings(stand_in, parents)
    monkeypatch.setattr(isc, "evaluate_plasmode_unit", lambda e, r, s, x: _fake_plasmode_record(e, r, s))
    assert isc.run_shard(tmp_path, settings, 0, 1, export=export) == 6
    summary = isc.merge(tmp_path, "plasmode")
    assert summary["complete_units"] == 6 and summary["failed_units"] == 0
    assert summary["reproduction_passed"] is False  # the stand-in parent has no committed rows
    for name in ("scan_table.csv", "block_balance.csv", "comparator_table.csv", "recovery_table.csv",
                 "accounting.json", "reproduction.json", "merged.jsonl"):
        assert (tmp_path / "plasmode" / name).exists()
    with pytest.raises(SystemExit, match="refusing to append"):
        isc.run_shard(tmp_path, settings, 0, 1, export=None)


# ── Decision ─────────────────────────────────────────────────────────────────


def _row(rate: float, count: int, **cell) -> dict:
    rejections = round(rate * count)
    rate = rejections / count
    return {**cell, "rate": rate, "se": (rate * (1 - rate) / count) ** 0.5, "count": count, "rejections": rejections}


def _tables(*, size=0.05, plasmode_size=0.05, false_hits=0.04, scan_e3=0.6, existing=0.05, overrides=None):
    overrides = overrides or {}
    cells = [cell[0] for cell in lsc.CELLS]
    sim_settings = {"profile": "simulation", "permutations": 999,
                    "rejection": {"scan": "p <= alpha", "parent": "p < alpha"},
                    "span_split": {"comparison": {"designs": {d: {} for d in lsc.DESIGNS}}}}
    scan, comparators, blocks, false = [], [], [], []
    for design in lsc.DESIGNS:
        for cell in cells:
            count = 500 if cell in isc.CONTROL_CELLS else 100
            rate = overrides.get(("simulation", design, cell), size if cell in isc.CONTROL_CELLS else 0.3)
            scan.append(_row(rate, count, design=design, cell=cell))
            for test in isc.COMPARATOR_TESTS:
                comparators.append(_row(rate if test == "scan" else 0.2, count, design=design, cell=cell, test=test))
            for block, share in (("methylation", 0.74), ("expression", 0.26)):
                blocks.append({"design": design, "cell": cell, "block": block, "n": count, "feature_share": share,
                               "null_max_share_mean": share, "null_max_share_se": 0.01})
            false.append({**_row(0.03, count, design=design, cell=cell), "n_population_null_median": 50})
    pm_settings = {"profile": "plasmode", "permutations": 999,
                   "rejection": {"scan": "p <= alpha", "parent": "p <= alpha"}}
    labels = ("anchor", "e=2", "e=3")
    pm_scan = [_row(overrides.get(("plasmode", label), plasmode_size if label == "anchor" else scan_e3), 500,
                    cell=label) for label in labels]
    pm_comparators = []
    for label in labels:
        for test in isc.COMPARATOR_TESTS:
            if test == "scan":
                rate = overrides.get(("plasmode", label), plasmode_size if label == "anchor" else scan_e3)
            else:
                rate = overrides.get(("plasmode", label, test), existing)
            pm_comparators.append(_row(rate, 500, cell=label, test=test))
    recovery = []
    for label in ("e=2", "e=3"):
        rate = overrides.get(("false_hits", label), false_hits)
        base = _row(rate, 500)
        recovery.append({"cell": label, "false_hit_rate": base["rate"], "false_hit_se": base["se"],
                         "false_hit_count": 500, "false_hit_rejections": base["rejections"],
                         "recall_mean": 0.3, "hit_precision_mean": 0.95, "top_k_precision_mean": 0.5})
    pm_blocks = [{"cell": label, "block": block, "n": 500, "feature_share": 0.5, "null_max_share_mean": share,
                  "null_max_share_se": 0.02} for label in labels for block, share in (("rna", 0.45), ("atac", 0.55))]
    reproduced = {"passed": True, "rows_compared": 84, "differing": [], "rule": "p < alpha", "source_dir": "x"}
    accounting = {"missing_units": 0, "failed_units": 0}
    simulation = {"settings": sim_settings, "scan": scan, "comparators": comparators, "blocks": blocks,
                  "false_hits": false, "reproduction": dict(reproduced), "accounting": accounting}
    plasmode = {"settings": pm_settings, "scan": pm_scan, "comparators": pm_comparators, "blocks": pm_blocks,
                "recovery": recovery, "reproduction": dict(reproduced), "accounting": accounting}
    return simulation, plasmode


ATTRIBUTION = [{"cell": "e=3", "arm": "adjusted", "component": "residual", "mean": "0.4323", "se": "0.0104",
                "n": "500"}]


def test_validated() -> None:
    result = isc.decision(*_tables(), attribution_precision=ATTRIBUTION)
    assert result["verdict"] == "scan_validated"
    criteria = result["criteria"]
    assert (criteria["S1"]["m"], criteria["S2"]["m"], criteria["S3"]["m"]) == (6, 1, 2)
    assert criteria["G"]["pass"] and criteria["G"]["comparator"] == "delta"
    assert criteria["G"]["difference"] == pytest.approx(0.55)
    assert "measurement-visibility" in result["consequence"] and "Phase 6" in result["consequence"]
    descriptive = result["descriptive"]
    assert set(descriptive) == {"plasmode_power_and_recovery", "block_balance", "simulation_localization",
                                "simulation_false_hits"}
    power = descriptive["plasmode_power_and_recovery"]
    assert [row["cell"] for row in power] == ["e=2", "e=3"]
    assert power[1]["attribution_planted_precision"][0]["planted_precision"] == 0.4323
    assert power[1]["attribution_planted_precision"][0]["source"] == "plasmode 2026-11-02"
    for table in descriptive.values():
        assert all("source" in row and "grade" in row for row in table)


def test_simulation_size_excess_names_the_cells() -> None:
    result = isc.decision(*_tables(overrides={("simulation", "seaad_weak", "anchor"): 0.12,
                                              ("simulation", "n200", "translation"): 0.11}))
    assert result["verdict"] == "size_excess"
    assert {(row["design"], row["cell"]) for row in result["failing"]["S1"]} == {
        ("seaad_weak", "anchor"), ("n200", "translation")
    }


def test_size_excess_takes_precedence_over_failing_familywise_control() -> None:
    result = isc.decision(*_tables(overrides={("plasmode", "anchor"): 0.09, ("false_hits", "e=3"): 0.2}))
    assert result["verdict"] == "size_excess"
    assert result["failing"]["S2"][0]["cell"] == "anchor" and result["failing"]["S2"][0]["rate"] == 0.09
    assert result["failing"]["S3"][0]["cell"] == "e=3"


def test_familywise_excess_takes_precedence_over_a_failing_gain() -> None:
    result = isc.decision(*_tables(scan_e3=0.06, overrides={("false_hits", "e=3"): 0.2}))
    assert result["verdict"] == "fwer_excess"
    assert result["failing"]["S3"][0]["status"] == "above_ceiling"
    assert result["failing"]["S3"][0]["rate"] == 0.2
    assert "G" in result["failing"]


def test_fwer_excess_when_gain_holds() -> None:
    result = isc.decision(*_tables(overrides={("false_hits", "e=2"): 0.15}))
    assert result["verdict"] == "fwer_excess"
    assert set(result["failing"]) == {"S3"}


def test_no_sparse_gain_names_the_comparator() -> None:
    result = isc.decision(*_tables(scan_e3=0.12, overrides={("plasmode", "e=3", "split_complement"): 0.10}))
    assert result["verdict"] == "no_sparse_gain"
    assert result["failing"]["G"]["comparator"] == "split_complement"
    assert result["failing"]["G"]["difference"] == pytest.approx(0.02)


def test_missing_comparator_fails_g() -> None:
    simulation, plasmode = _tables()
    plasmode["comparators"] = [r for r in plasmode["comparators"] if not (r["cell"] == "e=3" and r["test"] == "shape")]
    result = isc.decision(simulation, plasmode)
    assert result["verdict"] == "no_sparse_gain"
    assert result["criteria"]["G"]["missing"] == ["shape"]


def test_reproduction_failure_gates_every_criterion() -> None:
    simulation, plasmode = _tables()
    plasmode["reproduction"] = {"passed": False, "differing": [{"cell": ["e=3"], "test": "angle"}],
                                "rows_compared": 18}
    result = isc.decision(simulation, plasmode)
    assert result["verdict"] == "reproduction_failed"
    assert "criteria" not in result and "descriptive" not in result
    assert result["failing"]["reproduction"] == [
        {"profile": "plasmode", "differing": [{"cell": ["e=3"], "test": "angle"}]}
    ]


def test_descriptive_rates_do_not_change_the_verdict() -> None:
    first = isc.decision(*_tables())
    simulation, plasmode = _tables(overrides={("simulation", "n200", "shape_e0.5"): 0.0})
    plasmode["recovery"][1]["recall_mean"] = 0.99
    second = isc.decision(simulation, plasmode)
    assert first["verdict"] == second["verdict"] == "scan_validated"
    assert first["criteria"] == second["criteria"]
    assert first["descriptive"] != second["descriptive"]


# ── sbatch ───────────────────────────────────────────────────────────────────


def test_sbatch_template_parses() -> None:
    subprocess.run(["bash", "-n", str(REPO / "scripts" / "interaction_scan_check_array.sbatch")], check=True)
