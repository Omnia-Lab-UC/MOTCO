"""``scripts/interaction_scan_blocks_check.py``: settings binding, the unequal subset, both arms, merge and decision.

The simulation units run on the comparison's ``seaad`` design; the plasmode
profiles run on the plasmode tests' synthetic stand-in export. No real data are
read.
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


isb = _script("interaction_scan_blocks_check")
isc = isb.isc
ssc = isb.ssc
lsc = isb.lsc
pc = isb.pc

REPO = Path(__file__).resolve().parents[1]
FAST_CV = {"cv1_splits": 3, "cv2_splits": 4, "n_repeats": 1, "max_components": 6, "random_state": 1203}


# ── Simulation settings ──────────────────────────────────────────────────────


def test_simulation_settings_rebuild_and_bind_the_parent() -> None:
    settings = isb.simulation_settings()
    committed = json.loads((isb.PARENT_DIR / "simulation" / "settings.json").read_text())["settings"]
    assert settings["interaction_scan"] == committed
    units = isb.simulation_units(settings)
    assert units == isc.simulation_units(committed) and len(units) == 6600
    assert settings["arms"] == ["joint", "by_block"] and settings["rejection"] == "p <= alpha"
    assert settings["scan"]["residual_correction"] is None
    files = settings["parent"]["files"]
    assert set(files) == {"settings.json", "scan_table.csv", "block_balance.csv", "false_hit_table.csv"}
    for name, digest in files.items():
        assert digest == ssc.sha256_file(isb.PARENT_DIR / "simulation" / name)
    isb.check_parent_files(settings)


def _parent_copy(tmp_path: Path, profile: str = "simulation") -> Path:
    target = tmp_path / "interaction-scan" / profile
    target.mkdir(parents=True)
    for name in isb.PARENT_FILES[profile]:
        shutil.copy(isb.PARENT_DIR / profile / name, target / name)
    return target


def test_a_changed_parent_settings_file_is_refused_by_name(tmp_path) -> None:
    parent = _parent_copy(tmp_path)
    stored = json.loads((parent / "settings.json").read_text())
    stored["settings"]["alpha"] = 0.1
    stored["hash"] = lsc.settings_hash(stored["settings"])
    (parent / "settings.json").write_text(json.dumps(stored))
    with pytest.raises(SystemExit, match=r"interaction-scan/simulation/settings.json differs .*differing: alpha\."):
        isb.simulation_settings(parent_dir=parent)


def test_a_changed_parent_table_is_refused_by_name(tmp_path) -> None:
    parent = _parent_copy(tmp_path)
    settings = isb.simulation_settings(parent_dir=parent, designs=["seaad"], replicates=1, control_replicates=1,
                                       permutations=9)
    isb.check_parent_files(settings)
    with (parent / "block_balance.csv").open("a") as handle:
        handle.write("\n")
    with pytest.raises(SystemExit, match=r"block_balance.csv has changed since the settings were written"):
        isb.check_parent_files(settings)
    with pytest.raises(SystemExit, match=r"block_balance.csv has changed"):
        isb.run_shard(tmp_path / "run", settings, 0, 1, limit=0)


# ── Simulation units ─────────────────────────────────────────────────────────


def _small_settings(**overrides):
    settings = isb.simulation_settings(designs=["seaad"], replicates=2, control_replicates=2, permutations=39)
    settings.update(overrides)
    return settings


UNIT = lsc.Unit("seaad", "orient_e0.5_phi90_k20", 1)


@pytest.fixture(scope="module")
def simulation_unit():
    settings = _small_settings()
    return settings, isb.evaluate_simulation_unit(UNIT, settings)


def test_simulation_record_holds_both_arms(simulation_unit) -> None:
    _, record = simulation_unit
    assert record["status"] == "completed", record.get("error_message")
    assert (record["design"], record["cell"], record["replicate"]) == ("seaad", "orient_e0.5_phi90_k20", 1)
    assert record["rrpp_seed"] == lsc.rrpp_seed(2035, "seaad", "orient_e0.5_phi90_k20", 1)
    scan = record["scan"]
    assert scan["n_features"] == 367 + 131 and set(scan["blocks"]) == {"methylation", "expression"}
    assert sum(entry["null_max_share"] for entry in scan["blocks"].values()) == pytest.approx(1.0)
    assert sum(entry["n_population_null"] for entry in scan["blocks"].values()) == scan["n_population_null"]
    assert set(record["arms"]) == {"joint", "by_block"}
    for arm in record["arms"].values():
        assert sum(entry["hits"] for entry in arm["blocks"].values()) == arm["hits"]
        assert sum(entry["hits_population_null"] for entry in arm["blocks"].values()) == arm["hits_population_null"]
    by_block = record["arms"]["by_block"]
    assert by_block["n_blocks_tested"] == 2
    assert by_block["p_value_global"] == min(1.0, 2 * min(e["p_value"] for e in by_block["blocks"].values()))
    assert record["block_maxima"]["n_blocks_tested"] == 2
    assert -1.0 <= record["block_maxima"]["spearman"] <= 1.0
    assert 0 < record["block_maxima"]["minp_p_value_global"] <= 1.0
    assert set(record["runtime_seconds"]) == {"inputs", "scans", "unit"}


def test_simulation_record_holds_no_feature_name(simulation_unit) -> None:
    settings, record = simulation_unit
    dataset = lsc.unit_dataset(UNIT, settings["interaction_scan"]["span_split"]["comparison"])
    text = json.dumps(record)
    for frame in (dataset.methylation, dataset.expression):
        for feature in frame.columns:
            assert f'"{feature}"' not in text


def _inputs(settings, unit=UNIT):
    span_split = settings["interaction_scan"]["span_split"]
    dataset = lsc.unit_dataset(unit, span_split["comparison"])
    seed = lsc.rrpp_seed(2035, unit.design, unit.cell, unit.replicate)
    params = ssc.production_params(span_split, seed=seed, n_stages=3)
    return dataset, params, seed, isb.simulation_inputs(dataset, params)


def test_both_arms_share_every_f(simulation_unit) -> None:
    settings, _ = simulation_unit
    _, _, seed, (features, full, reduced) = _inputs(settings)
    labels = isb.block_labels(features.columns)
    scans = isb.run_arms(features.to_numpy(dtype=float), full, reduced, labels, permutations=19, seed=seed,
                         residual_correction=None)
    np.testing.assert_array_equal(scans["joint"].f, scans["by_block"].f)
    np.testing.assert_array_equal(scans["joint"].null_max_f, scans["by_block"].null_max_f)
    assert scans["joint"].combination == "joint" and scans["by_block"].combination == "block_bonferroni"


def test_simulation_inputs_give_the_pls_input_coordinates_f(simulation_unit) -> None:
    from motco.simulations.evaluation import integrate_semisynthetic_dataset
    from motco.stats.interaction import interaction_feature_scan

    settings, _ = simulation_unit
    fast = json.loads(json.dumps(settings))
    fast["interaction_scan"]["span_split"]["comparison"]["pls_cv_params"] = dict(FAST_CV)
    dataset, params, _, (features, full, reduced) = _inputs(fast)
    latent = integrate_semisynthetic_dataset(dataset, params)
    pls_features, _, pls_full, pls_reduced = ssc.split_inputs(dataset, latent, params)
    assert list(features.columns) == list(latent.artifacts.features.columns)
    np.testing.assert_array_equal(full, pls_full)
    np.testing.assert_array_equal(reduced, pls_reduced)
    mine = interaction_feature_scan(features.to_numpy(dtype=float), full, reduced)
    theirs = interaction_feature_scan(pls_features, full, reduced)
    np.testing.assert_allclose(mine.f, theirs.f, rtol=1e-10)


def test_the_joint_arm_equals_the_parents_unit(simulation_unit) -> None:
    settings, record = simulation_unit
    fast = json.loads(json.dumps(settings["interaction_scan"]))
    fast["span_split"]["comparison"]["pls_cv_params"] = dict(FAST_CV)
    fast["span_split"]["reproduction_permutations"] = 19
    theirs = isc.evaluate_simulation_unit(UNIT, fast)
    assert theirs["status"] == "completed", theirs.get("error_message")
    view = isb.parent_view(record)["scan"]
    for key in ("n_features", "p_value_global", "hits", "n_population_null", "hits_population_null"):
        assert view[key] == theirs["scan"][key], key
    assert view["max_f"] == pytest.approx(theirs["scan"]["max_f"], rel=1e-10)
    for name, entry in theirs["scan"]["blocks"].items():
        assert view["blocks"][name]["hits"] == entry["hits"]
        assert view["blocks"][name]["null_max_share"] == entry["null_max_share"]


def test_a_failing_unit_is_recorded_without_aborting_the_shard(tmp_path, monkeypatch) -> None:
    settings = _small_settings()

    def broken(unit, settings):
        raise RuntimeError(f"no data for {unit.cell}")

    monkeypatch.setattr(lsc, "unit_dataset", broken)
    assert isb.run_shard(tmp_path, settings, 0, 1, limit=2) == 2
    lines = (tmp_path / "simulation" / "records" / "shard_000_of_001.jsonl").read_text().splitlines()
    records = [json.loads(line) for line in lines]
    assert [r["status"] for r in records] == ["failed", "failed"]
    assert records[0]["error_type"] == "RuntimeError" and "no data for" in records[0]["error_message"]


# ── The block maxima ─────────────────────────────────────────────────────────


def _scan(seed=0, permutations=59):
    from motco.stats.interaction import interaction_feature_scan

    rng = np.random.default_rng(seed)
    groups = np.repeat([0, 1], 30)
    stages = np.tile(np.repeat([0, 1, 2], 10), 2)
    full = np.column_stack([np.ones(60), groups, stages == 1, stages == 2, groups * (stages == 1),
                            groups * (stages == 2)]).astype(float)
    Y = rng.normal(size=(60, 16))
    Y[:, 2] += 1.2 * groups * stages
    labels = ["a"] * 6 + ["b"] * 10
    return interaction_feature_scan(Y, full, full[:, :4], permutations=permutations, seed=seed, blocks=labels,
                                    return_null=True)


def test_minp_matches_a_direct_computation() -> None:
    from scipy.stats import spearmanr

    scan = _scan()
    maxima = isb.block_maxima(scan)
    null = scan.null_block_max_f
    B = null.shape[0]
    observed = [(1 + np.sum(null[:, k] >= block.max_f)) / (1 + B) for k, block in enumerate(scan.blocks)]
    reference = [np.concatenate([[block.max_f], null[:, k]]) for k, block in enumerate(scan.blocks)]
    draws = [min(np.sum(reference[k] >= null[b, k]) / (1 + B) for k in range(2)) for b in range(B)]
    expected = (1 + sum(1 for value in draws if value <= min(observed))) / (1 + B)
    assert maxima["minp_p_value_global"] == pytest.approx(expected)
    assert maxima["spearman"] == pytest.approx(spearmanr(null[:, 0], null[:, 1]).statistic)
    assert maxima["minp_p_value_global"] <= scan.p_value_global + 1e-12  # minP never loses to Bonferroni here


def test_block_maxima_without_permutations_are_empty() -> None:
    scan = _scan(permutations=0)
    assert isb.block_maxima(scan) == {"n_blocks_tested": 2, "spearman": None, "minp_p_value_global": None}


# ── Plasmode profiles (stand-in export) ──────────────────────────────────────


def _plasmode_tests():
    sys.path.insert(0, str(REPO))
    try:
        from tests import test_plasmode_check as module
    finally:
        sys.path.pop(0)
    return module


plasmode_tests = _plasmode_tests()
STAND_IN_SUBSET = 12


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
    """Stand-ins for the plasmode run, the split-test plasmode run and the interaction-scan plasmode run, 2 × 9."""

    _, export, fidelity = stand_in
    plasmode = pc.build_settings(export, fidelity, replicates=2, permutations=9, bootstrap=5, e_star=3.0,
                                 descriptive=2.0)
    plasmode["cv"] = dict(plasmode_tests.SMALL_CV)
    summary_columns = ("effect", "cell", "arm", "block", "component", "subset", "quantity", "n", "mean", "se", "role")
    plasmode_dir = _write(
        tmp_path_factory.mktemp("plasmode"), plasmode, {"plasmode_summary.csv": ([], summary_columns)}
    )
    span_split = ssc.plasmode_settings(export, fidelity, parent_dir=plasmode_dir, replicates=2, permutations=9)
    split_dir = _write(
        tmp_path_factory.mktemp("span-split") / "plasmode",
        span_split,
        {"split_table.csv": ([], ssc.PLASMODE_SPLIT_COLUMNS), "trajectory_table.csv": ([], ssc.TRAJECTORY_COLUMNS)},
    )
    scan = isc.plasmode_settings(export, fidelity, parent_dir=split_dir, plasmode_dir=plasmode_dir, replicates=2,
                                 permutations=9)
    scan_dir = _write(
        tmp_path_factory.mktemp("interaction-scan") / "plasmode",
        scan,
        {
            "scan_table.csv": ([], isc.SCAN_COLUMNS["plasmode"]),
            "block_balance.csv": ([], isc.BLOCK_COLUMNS["plasmode"]),
            "recovery_table.csv": ([], isc.RECOVERY_COLUMNS),
        },
    )
    return plasmode_dir, split_dir, scan_dir


def _plasmode_settings(stand_in, parents, profile="plasmode", **kwargs):
    _, export, _ = stand_in
    plasmode_dir, split_dir, scan_dir = parents
    kwargs.setdefault("replicates", 2)
    kwargs.setdefault("permutations", 9)
    kwargs.setdefault("unequal_size", STAND_IN_SUBSET)
    return isb.plasmode_settings(export, profile=profile, parent_dir=scan_dir, span_split_dir=split_dir,
                                 plasmode_dir=plasmode_dir, **kwargs)


def test_plasmode_settings_rebuild_and_bind(stand_in, parents) -> None:
    settings = _plasmode_settings(stand_in, parents)
    assert isb.plasmode_units(settings) == [(0.0, 0), (2.0, 0), (3.0, 0), (0.0, 1), (2.0, 1), (3.0, 1)]
    assert settings["scan"]["residual_correction"] == "df" and "subset" not in settings
    assert set(settings["parent"]["files"]) == {"settings.json", "scan_table.csv", "block_balance.csv",
                                                "recovery_table.csv"}
    isb.check_parent_files(settings)
    unequal = _plasmode_settings(stand_in, parents, profile="plasmode_unequal")
    assert set(unequal["parent"]["files"]) == {"settings.json"}
    assert unequal["interaction_scan"] == settings["interaction_scan"]


def test_an_export_whose_hashes_differ_is_refused_by_name(stand_in, parents, tmp_path) -> None:
    directory, _, _ = stand_in
    plasmode_dir, split_dir, scan_dir = parents
    copy = shutil.copytree(directory, tmp_path / "export")
    frame = (copy / "atac.csv").read_text().splitlines()
    (copy / "atac.csv").write_text("\n".join(frame[:-1] + [frame[-1].replace(",", ",1", 1)]) + "\n")
    with pytest.raises(SystemExit, match=r"differ from .*settings.json for: atac.csv; refusing"):
        isb.plasmode_settings(pc.load_export(copy), parent_dir=scan_dir, span_split_dir=split_dir,
                              plasmode_dir=plasmode_dir, replicates=2, permutations=9)


def test_the_committed_plasmode_parent_must_match(stand_in) -> None:
    _, export, _ = stand_in
    # The real interaction-scan run was made on the real export, not the stand-in.
    with pytest.raises(SystemExit, match=r"differ from .*settings.json for: .*; refusing another export"):
        isb.plasmode_settings(export)


def test_the_unequal_subset_is_fixed_and_recorded(stand_in, parents) -> None:
    _, export, _ = stand_in
    first = _plasmode_settings(stand_in, parents, profile="plasmode_unequal")["subset"]
    second = _plasmode_settings(stand_in, parents, profile="plasmode_unequal")["subset"]
    assert first == second
    assert (first["block"], first["size"], first["of"]) == ("rna", STAND_IN_SUBSET, 30)
    assert first["seed"] == isb.unequal_seed(2034) != isb.unequal_seed(2035)
    reduced = isb.reduced_export(export, first)
    retained = list(reduced.blocks["rna"].columns)
    assert len(retained) == STAND_IN_SUBSET
    assert retained == [name for name in export.blocks["rna"].columns if name in set(retained)]  # export order
    assert isb._features_sha256(retained) == first["features_sha256"]
    pd_atac = reduced.blocks["atac"]
    assert pd_atac is export.blocks["atac"] and export.blocks["rna"].shape[1] == 30
    with pytest.raises(SystemExit, match="reduced feature list differs"):
        isb.reduced_export(export, {**first, "features_sha256": "0" * 64})


def test_planted_features_come_from_the_reduced_block(stand_in, parents) -> None:
    _, export, _ = stand_in
    settings = _plasmode_settings(stand_in, parents, profile="plasmode_unequal")
    reduced = isb.unit_export(export, settings)
    retained = set(reduced.blocks["rna"].columns)
    parent = settings["interaction_scan"]["span_split"]["plasmode"]
    for replicate in range(4):
        built = pc.build_unit(2.0, replicate, parent, reduced)
        assert len(built["plant"]["rna"]["features"]) == 10
        assert set(built["plant"]["rna"]["features"]) <= retained


@pytest.fixture(scope="module")
def plasmode_record(stand_in, parents):
    _, export, _ = stand_in
    settings = _plasmode_settings(stand_in, parents)
    return settings, isb.evaluate_plasmode_unit(3.0, 1, settings, export)


def test_plasmode_record_holds_seeds_statistics_and_counts_only(stand_in, plasmode_record) -> None:
    _, export, _ = stand_in
    _, record = plasmode_record
    assert record["status"] == "completed", record.get("error_message")
    text = json.dumps(record)
    for donor in export.metadata[export.sample_col]:
        assert donor not in text
    assert "Female" not in text and "Male" not in text
    for frame in export.blocks.values():
        for feature in frame.columns:
            assert f'"{feature}"' not in text
    assert set(record) == {"profile", "key", "effect", "cell", "replicate", "stream", "seeds", "settings_hash",
                           "status", "label_redraws", "scan", "arms", "block_maxima", "runtime_seconds"}
    planted = record["scan"]["planted"]
    assert planted["planted"] == 20 and planted["top_k"] == 20
    for arm in record["arms"].values():
        counts = arm["planted"]
        assert counts["hits"] == arm["hits"]
        assert counts["hits_in_planted"] + counts["hits_outside_planted"] == arm["hits"]
        assert sum(entry["hits"] for entry in counts["blocks"].values()) == arm["hits"]


def test_plasmode_scan_inputs_are_the_analysis(stand_in, plasmode_record) -> None:
    _, export, _ = stand_in
    settings, record = plasmode_record
    built = pc.build_unit(3.0, 1, settings["interaction_scan"]["span_split"]["plasmode"], export)
    _, result = pc.analyse(built["blocks"], built["metadata"], export,
                           isc.plasmode_params(settings["interaction_scan"], built["seeds"]))
    block = result.record["interaction_scan"]
    joint = record["arms"]["joint"]
    assert joint["p_value_global"] == block["p_value_global"]
    assert joint["hits"] == len(block["hits"])
    assert record["scan"]["max_f"] == block["max_f"] and record["scan"]["max_block"] == block["max_block"]
    for name, entry in block["blocks"].items():
        assert joint["blocks"][name]["hits"] == entry["hits"]
        assert record["scan"]["blocks"][name]["null_max_share"] == entry["null_max_share"]
    expected = isc.planted_counts(result.interaction_scan, built["plant"])
    view = isb.parent_view(record)["planted"]
    assert {key: view[key] for key in expected if key != "blocks"} == {
        key: value for key, value in expected.items() if key != "blocks"
    }


def test_the_unequal_unit_scans_the_reduced_block(stand_in, parents) -> None:
    _, export, _ = stand_in
    settings = _plasmode_settings(stand_in, parents, profile="plasmode_unequal")
    record = isb.evaluate_plasmode_unit(2.0, 0, settings, isb.unit_export(export, settings))
    assert record["status"] == "completed", record.get("error_message")
    assert record["profile"] == "plasmode_unequal"
    assert record["scan"]["blocks"]["rna"]["n_features"] == STAND_IN_SUBSET
    assert record["scan"]["blocks"]["atac"]["n_features"] == 30
    assert record["scan"]["planted"]["blocks"]["rna"]["planted"] == 10


def test_a_failing_plasmode_unit_is_recorded(stand_in, parents, monkeypatch) -> None:
    _, export, _ = stand_in
    settings = _plasmode_settings(stand_in, parents)

    def broken(*args, **kwargs):
        raise ValueError("no plant")

    monkeypatch.setattr(pc, "build_unit", broken)
    record = isb.evaluate_plasmode_unit(2.0, 0, settings, export)
    assert (record["status"], record["error_type"], record["error_message"]) == ("failed", "ValueError", "no plant")
    assert record["key"] == pc.unit_key(2.0, 0) and "scan" not in record


def test_a_plasmode_run_refuses_another_export(stand_in, parents, tmp_path) -> None:
    settings = _plasmode_settings(stand_in, parents)
    with pytest.raises(SystemExit, match="refusing to append"):
        isb.run_shard(tmp_path, settings, 0, 1, export=None)


# ── Run, resume and settings binding ─────────────────────────────────────────


def _arm(p: float, hits: int, blocks: dict, **extra) -> dict:
    return {"p_value_global": p, "hits": hits, "blocks": blocks, **extra}


def _fake_simulation_record(unit, settings, *, p_joint=0.5, p_block=0.5, hits=0, null_hits=0):
    base = isb._base_simulation(unit, settings)
    scan = {
        "n_features": 498, "n_degenerate": 0, "df_interaction": 2, "df_residual": 194, "max_f": 9.0,
        "max_block": "methylation", "n_population_null": 10,
        "blocks": {"methylation": {"n_features": 367, "n_degenerate": 0, "max_f": 9.0, "null_max_share": 0.7,
                                   "n_population_null": 8},
                   "expression": {"n_features": 131, "n_degenerate": 0, "max_f": 4.0, "null_max_share": 0.3,
                                  "n_population_null": 2}},
    }

    def blocks(p=None):
        out = {"methylation": {"hits": hits, "hits_population_null": null_hits},
               "expression": {"hits": 0, "hits_population_null": 0}}
        if p is not None:
            for entry in out.values():
                entry["p_value"] = p
        return out

    return {
        **base,
        "status": "completed",
        "scan": scan,
        "arms": {"joint": _arm(p_joint, hits, blocks(), hits_population_null=null_hits),
                 "by_block": _arm(p_block, hits, blocks(p_block / 2), hits_population_null=null_hits,
                                  n_blocks_tested=2)},
        "block_maxima": {"n_blocks_tested": 2, "spearman": 0.1, "minp_p_value_global": p_block},
        "runtime_seconds": {"unit": 1.0, "inputs": 0.6, "scans": 0.4},
    }


def test_resume_skips_completed_units(tmp_path, monkeypatch) -> None:
    settings = _small_settings()
    calls: list = []

    def fake(unit, settings):
        calls.append(unit.key)
        return _fake_simulation_record(unit, settings)

    monkeypatch.setattr(isb, "evaluate_simulation_unit", fake)
    assert isb.run_shard(tmp_path, settings, 0, 2, limit=3) == 3
    with (tmp_path / "simulation" / "records" / "shard_000_of_002.jsonl").open("a") as handle:
        handle.write('{"torn": ')  # an interrupted write
    first = list(calls)
    total = len(isb.simulation_units(settings))
    assert isb.run_shard(tmp_path, settings, 0, 2) == total // 2 - 3
    assert not set(first) & set(calls[3:])
    isb.run_shard(tmp_path, settings, 1, 2)
    summary = isb.merge(tmp_path, "simulation")
    assert summary["complete_units"] == total and summary["missing_units"] == 0
    assert summary["reproduction_passed"] is False  # 2 replicates are not the parent's 500
    assert len((tmp_path / "simulation" / "merged.jsonl").read_text().splitlines()) == total
    for name in (*isb.table_columns("simulation"), "accounting.json", "reproduction.json"):
        assert (tmp_path / "simulation" / name).exists(), name
    accounting = json.loads((tmp_path / "simulation" / "accounting.json").read_text())
    assert accounting["stage_seconds_median"] == {"inputs": 0.6, "scans": 0.4}


def test_settings_mismatch_is_refused_naming_the_keys(tmp_path, monkeypatch) -> None:
    monkeypatch.setattr(isb, "evaluate_simulation_unit", _fake_simulation_record)
    isb.run_shard(tmp_path, _small_settings(), 0, 1, limit=1)
    with pytest.raises(SystemExit, match="differing: .*permutations"):
        isb.run_shard(tmp_path, _small_settings(permutations=49), 0, 1, limit=1)


# ── Merge tables and reproduction ────────────────────────────────────────────


def _simulation_records(settings):
    records = []
    for unit in isb.simulation_units(settings):
        rejects = unit.replicate == 0
        records.append(_fake_simulation_record(unit, settings, p_joint=0.01 if rejects else 0.5,
                                               p_block=0.5 if rejects else 0.5, hits=3 if rejects else 0,
                                               null_hits=1 if rejects else 0))
    return records


def test_simulation_tables_count_failures_and_discordance() -> None:
    settings = _small_settings()
    records = _simulation_records(settings)
    records[0] = {**records[0], "status": "failed"}
    scan = {(r["cell"], r["arm"]): r for r in isb.scan_rows(records, settings)}
    assert len(scan) == 14 * 2
    anchor = scan[(records[0]["cell"], "joint")]
    assert (anchor["failures"], anchor["count"], anchor["rejections"]) == (1, 1, 0)
    translation = scan[("translation", "joint")]
    assert (translation["count"], translation["rejections"], translation["hits_median"]) == (2, 1, 1.5)
    assert scan[("translation", "by_block")]["rejections"] == 0
    paired = {r["cell"]: r for r in isb.paired_rows(records, settings)}
    assert (paired["translation"]["joint_only"], paired["translation"]["by_block_only"]) == (1, 0)
    assert paired["translation"]["p_loss"] == 0.5 and paired["translation"]["p_gain"] == 1.0
    blocks = {(r["cell"], r["arm"], r["block"]): r for r in isb.block_rows(records, settings)}
    row = blocks[("translation", "by_block", "methylation")]
    assert row["null_max_share_mean"] == pytest.approx(0.7) and row["feature_share"] == pytest.approx(367 / 498)
    assert row["hit_unit_rate"] == 0.5 and row["population_null_unit_rate"] == 0.5
    false_hits = {(r["cell"], r["arm"]): r for r in isb.false_hit_rows(records, settings)}
    assert (false_hits[("translation", "joint")]["count"], false_hits[("translation", "joint")]["rejections"]) == (2, 1)
    minp = {r["cell"]: r for r in isb.minp_rows(records, settings)}
    assert (minp["translation"]["joint_rate"], minp["translation"]["bonferroni_rate"]) == (0.5, 0.0)
    assert minp["translation"]["spearman_median"] == pytest.approx(0.1)


def _committed_tables(records, settings, tmp_path: Path) -> dict[str, list[dict]]:
    """The parent's tables of these records, written and read back as the committed CSVs are."""

    tables = isb._parent_tables(records, settings)
    columns = {"scan_table.csv": isc.SCAN_COLUMNS["simulation"], "block_balance.csv": isc.BLOCK_COLUMNS["simulation"],
               "false_hit_table.csv": isc.FALSE_HIT_COLUMNS}
    out = {}
    for name, rows in tables.items():
        lsc._write_csv(tmp_path / name, rows, columns[name])
        out[name] = lsc._read_csv(tmp_path / name)
    return out


def test_reproduction_passes_and_fails_naming_the_cell(tmp_path) -> None:
    settings = _small_settings()
    records = _simulation_records(settings)
    committed = _committed_tables(records, settings, tmp_path)
    passed = isb.reproduction(records, settings, committed)
    assert passed["passed"] and not passed["differing"]
    assert passed["rows_compared"] == 14 + 14 * 2 + 14
    bumped = json.loads(json.dumps(committed))
    bumped["scan_table.csv"][3]["rejections"] = str(int(bumped["scan_table.csv"][3]["rejections"]) + 1)
    failed = isb.reproduction(records, settings, bumped)
    assert not failed["passed"] and len(failed["differing"]) == 1
    entry = failed["differing"][0]
    assert entry["table"] == "scan_table.csv" and entry["cell"] == ["seaad", bumped["scan_table.csv"][3]["cell"]]
    assert int(entry["committed"]["rejections"]) == int(entry["reproduced"]["rejections"]) + 1
    other = json.loads(json.dumps(committed))
    other["scan_table.csv"] += [{**row, "design": "n200", "rejections": "77"} for row in committed["scan_table.csv"]]
    assert isb.reproduction(records, settings, other)["passed"]  # only the run's designs are read


def test_the_joint_arm_reproduces_through_the_parents_functions() -> None:
    settings = _small_settings()
    records = _simulation_records(settings)
    views = [isb.parent_view(record) for record in records]
    rows = {r["cell"]: r for r in isc.scan_rows(views, settings["interaction_scan"])}
    assert (rows["translation"]["count"], rows["translation"]["rejections"]) == (2, 1)
    false_hits = {r["cell"]: r for r in isc.false_hit_rows(views, settings["interaction_scan"])}
    assert false_hits["translation"]["rejections"] == 1


def _fake_plasmode_record(effect, replicate, settings, *, p_joint=0.5, p_block=0.5, rna=(0, 0), atac=(0, 0)):
    """``rna``/``atac``: (joint, by_block) planted hits in the block; no hit outside the planted set."""

    def planted(arm_index):
        blocks = {name: {"hits": pair[arm_index], "hits_in_planted": pair[arm_index], "hits_outside_planted": 0}
                  for name, pair in (("rna", rna), ("atac", atac))}
        total = sum(entry["hits"] for entry in blocks.values())
        return {"hits": total, "hits_in_planted": total, "hits_outside_planted": 0, "blocks": blocks}

    arms = {}
    for index, (arm, p) in enumerate((("joint", p_joint), ("by_block", p_block))):
        counts = planted(index)
        arms[arm] = {"p_value_global": p, "hits": counts["hits"],
                     "blocks": {name: {"hits": entry["hits"]} for name, entry in counts["blocks"].items()},
                     "planted": counts}
    return {
        "profile": settings["profile"],
        "key": pc.unit_key(effect, replicate),
        "effect": effect,
        "cell": pc.cell_label(effect),
        "replicate": replicate,
        "settings_hash": lsc.settings_hash(settings),
        "status": "completed",
        "scan": {"n_features": 60, "n_degenerate": 0, "max_f": 12.0, "max_block": "rna",
                 "blocks": {"rna": {"n_features": 30, "n_degenerate": 0, "max_f": 12.0, "null_max_share": 0.4},
                            "atac": {"n_features": 30, "n_degenerate": 0, "max_f": 8.0, "null_max_share": 0.6}},
                 "planted": {"top_k": 20, "planted": 20, "top_k_in_planted": 5,
                             "blocks": {name: {"top_k_in_planted": 2, "planted": 10} for name in ("rna", "atac")}}},
        "arms": arms,
        "block_maxima": {"n_blocks_tested": 2, "spearman": 0.2, "minp_p_value_global": p_block},
        "runtime_seconds": {"unit": 2.0},
    }


def test_plasmode_tables_and_allocation(stand_in, parents) -> None:
    settings = _plasmode_settings(stand_in, parents, profile="plasmode_unequal")
    records = [
        _fake_plasmode_record(2.0, 0, settings, p_joint=0.01, p_block=0.01, rna=(1, 3), atac=(2, 1)),
        _fake_plasmode_record(2.0, 1, settings, p_joint=0.01, p_block=0.2, rna=(0, 0), atac=(2, 2)),
        _fake_plasmode_record(3.0, 0, settings, rna=(2, 4)),
        _fake_plasmode_record(3.0, 1, settings, rna=(2, 5)),
        _fake_plasmode_record(0.0, 0, settings),
    ]
    allocation = {(r["cell"], r["block"]): r for r in isb.allocation_rows(records, settings)}
    assert set(allocation) == {("e=2", "rna"), ("e=2", "atac"), ("e=3", "rna"), ("e=3", "atac")}
    e2 = allocation[("e=2", "rna")]
    assert (e2["n"], e2["gains"], e2["losses"], e2["ties"], e2["p_sign"]) == (2, 1, 0, 1, 0.5)
    assert allocation[("e=3", "rna")]["p_sign"] == 0.25 and allocation[("e=2", "atac")]["losses"] == 1
    blocks = {(r["cell"], r["arm"], r["block"]): r for r in isb.block_rows(records, settings)}
    rna = blocks[("e=3", "by_block", "rna")]
    assert rna["recall_mean"] == pytest.approx((4 + 5) / 2 / 10) and rna["hit_precision_mean"] == 1.0
    assert blocks[("anchor", "joint", "rna")].get("recall_mean") is None
    false_hits = {(r["cell"], r["arm"]): r for r in isb.false_hit_rows(records, settings)}
    assert set(false_hits) == {("e=2", "joint"), ("e=2", "by_block"), ("e=3", "joint"), ("e=3", "by_block")}
    paired = {r["cell"]: r for r in isb.paired_rows(records, settings)}
    assert (paired["e=2"]["joint_only"], paired["e=2"]["by_block_only"], paired["e=2"]["both"]) == (1, 0, 1)


def test_plasmode_merge_writes_its_tables(stand_in, parents, tmp_path, monkeypatch) -> None:
    _, export, _ = stand_in
    for profile in ("plasmode", "plasmode_unequal"):
        settings = _plasmode_settings(stand_in, parents, profile=profile)
        monkeypatch.setattr(isb, "evaluate_plasmode_unit", lambda e, r, s, x: _fake_plasmode_record(e, r, s))
        assert isb.run_shard(tmp_path, settings, 0, 1, export=export) == 6
        summary = isb.merge(tmp_path, profile)
        assert summary["complete_units"] == 6 and summary["failed_units"] == 0
        for name in (*isb.table_columns(profile), "accounting.json", "merged.jsonl"):
            assert (tmp_path / profile / name).exists(), name
        assert (tmp_path / profile / "reproduction.json").exists() == (profile == "plasmode")
        if profile == "plasmode":
            assert summary["reproduction_passed"] is False  # the stand-in parent has no committed rows
        else:
            assert "reproduction_passed" not in summary


def test_decide_reads_the_three_merged_profiles(stand_in, parents, tmp_path, monkeypatch) -> None:
    _, export, _ = stand_in
    monkeypatch.setattr(isb, "evaluate_simulation_unit", _fake_simulation_record)
    monkeypatch.setattr(isb, "evaluate_plasmode_unit", lambda e, r, s, x: _fake_plasmode_record(e, r, s))
    isb.run_shard(tmp_path, _small_settings(), 0, 1)
    isb.merge(tmp_path, "simulation")
    for profile in ("plasmode", "plasmode_unequal"):
        isb.run_shard(tmp_path, _plasmode_settings(stand_in, parents, profile=profile), 0, 1, export=export)
        isb.merge(tmp_path, profile)
    result = isb.decide(tmp_path)
    assert result["verdict"] == "reproduction_failed"  # 2 replicates do not reproduce the parent
    assert json.loads((tmp_path / isb.DECISION_NAME).read_text())["verdict"] == "reproduction_failed"
    assert set(result["settings_hashes"]) == set(isb.PROFILES)
    (tmp_path / "plasmode_unequal" / "accounting.json").unlink()
    with pytest.raises(SystemExit, match="run `merge --profile plasmode_unequal` first"):
        isb.decide(tmp_path)


# ── Decision ─────────────────────────────────────────────────────────────────


DESIGNS = ("n200", "seaad", "seaad_weak")


def _rate_row(rejections: int, count: int = 500, **cell) -> dict:
    return {**cell, "count": str(count), "rejections": str(rejections), "rate": str(rejections / count)}


def _tables(*, size=25, plasmode_size=25, outside=21, loss=(10, 10), gain=(30, 5), reproduction=True, overrides=None):
    """Three profiles' merged tables; every criterion holds by default."""

    overrides = overrides or {}

    def settings(profile):
        out = {"profile": profile, "permutations": 999}
        if profile == "simulation":
            out["interaction_scan"] = {"span_split": {"comparison": {"designs": list(DESIGNS)}}}
        if profile == "plasmode_unequal":
            out["subset"] = {"block": "rna", "size": 500}
        return out

    accounting = {"missing_units": 0, "failed_units": 0}
    simulation = {
        "settings": settings("simulation"),
        "accounting": accounting,
        "reproduction": {"passed": reproduction, "rows_compared": 10, "rule": "p <= alpha", "differing": []
                         if reproduction else [{"table": "scan_table.csv", "cell": ["n200", "anchor"]}]},
        "scan_table.csv": [
            _rate_row(overrides.get(("simulation", design, cell, arm), size), design=design, cell=cell, arm=arm)
            for design in DESIGNS
            for cell in ("anchor", "translation", "magnitude")
            for arm in ("joint", "by_block")
        ],
        "false_hit_table.csv": [_rate_row(5, design="n200", cell="anchor", arm="by_block")],
        "block_table.csv": [{"design": "n200", "cell": "anchor", "arm": "by_block", "block": "methylation",
                             "n": "500", "feature_share": "0.737", "null_max_share_mean": "0.5",
                             "hit_unit_rate": "0.03", "hits_mean": "0.05"}],
        "minp_table.csv": [{"design": "n200", "cell": "anchor", "count": "500", "joint_rate": "0.05",
                            "bonferroni_rate": "0.05", "minp_rate": "0.051", "spearman_median": "0.0"}],
    }

    def plasmode_block(profile):
        scan = [_rate_row(overrides.get((profile, cell, arm), plasmode_size if cell == "anchor" else 400), cell=cell,
                          arm=arm, effect=str(effect))
                for cell, effect in (("anchor", 0), ("e=2", 2), ("e=3", 3)) for arm in ("joint", "by_block")]
        false_hits = [_rate_row(overrides.get((profile, "outside", cell), outside), cell=cell, arm=arm)
                      for cell in ("e=2", "e=3") for arm in ("joint", "by_block")]
        joint_only, by_block_only = overrides.get((profile, "loss"), loss)
        paired = [{"cell": cell, "effect": str(effect), "n": "500", "joint_only": str(joint_only),
                   "by_block_only": str(by_block_only),
                   "p_loss": str(isb._upper_binomial(joint_only, joint_only + by_block_only))}
                  for cell, effect in (("anchor", 0), ("e=2", 2), ("e=3", 3))]
        block = {
            "settings": settings(profile),
            "accounting": accounting,
            "scan_table.csv": scan,
            "false_hit_table.csv": false_hits,
            "paired_table.csv": paired,
            "block_table.csv": [{"cell": "e=3", "effect": "3", "arm": "by_block", "block": "rna", "n": "500",
                                 "recall_mean": "0.6", "hit_precision_mean": "0.99", "outside_unit_rate": "0.02"}],
            "minp_table.csv": [{"cell": "anchor", "effect": "0", "count": "500", "joint_rate": "0.05",
                                "bonferroni_rate": "0.05", "minp_rate": "0.05", "spearman_median": "0.1"}],
        }
        if profile == "plasmode":
            block["reproduction"] = {"passed": True, "rows_compared": 12, "rule": "p <= alpha", "differing": []}
        else:
            gains, losses = overrides.get(("unequal", "gain"), gain)
            block["allocation_table.csv"] = [
                {"cell": cell, "effect": str(effect), "block": "rna", "n": "500", "gains": str(gains),
                 "losses": str(losses), "ties": str(500 - gains - losses),
                 "p_sign": str(isb._upper_binomial(gains, gains + losses))}
                for cell, effect in (("e=2", 2), ("e=3", 3))
            ]
        return block

    return simulation, plasmode_block("plasmode"), plasmode_block("plasmode_unequal")


def test_validated() -> None:
    result = isb.decision(*_tables())
    assert result["verdict"] == "by_block_validated"
    assert result["failing"] == {}
    assert set(result["criteria"]) == {"S1", "S2", "S3", "P", "G"}
    assert result["criteria"]["S1"]["m"] == 6 and result["criteria"]["S2"]["m"] == 2
    assert result["criteria"]["S3"]["m"] == 4 and len(result["criteria"]["P"]["observations"]) == 4
    assert "joint scan stays the default" in result["consequence"]
    assert "validated for exports whose blocks differ in size or noise" in result["consequence"]


def test_simulation_size_excess_names_the_cell() -> None:
    result = isb.decision(*_tables(overrides={("simulation", "seaad", "translation", "by_block"): 60}))
    assert result["verdict"] == "size_excess"
    assert [(row["design"], row["cell"]) for row in result["failing"]["S1"]] == [("seaad", "translation")]


def test_the_joint_arm_is_not_judged_for_size() -> None:
    result = isb.decision(*_tables(overrides={("simulation", "seaad", "translation", "joint"): 60}))
    assert result["verdict"] == "by_block_validated"


def test_real_noise_size_excess() -> None:
    result = isb.decision(*_tables(overrides={("plasmode_unequal", "anchor", "by_block"): 60}))
    assert result["verdict"] == "size_excess"
    assert [row["profile"] for row in result["failing"]["S2"]] == ["plasmode_unequal"]


def test_familywise_excess_names_the_cell() -> None:
    result = isb.decision(*_tables(overrides={("plasmode_unequal", "outside", "e=3"): 60}))
    assert result["verdict"] == "fwer_excess"
    assert [(row["profile"], row["cell"]) for row in result["failing"]["S3"]] == [("plasmode_unequal", "e=3")]


def test_power_loss_names_the_cell_and_both_counts() -> None:
    result = isb.decision(*_tables(overrides={("plasmode", "loss"): (30, 5)}))
    assert result["verdict"] == "power_loss"
    failing = result["failing"]["P"]
    assert {(row["profile"], row["cell"]) for row in failing} == {("plasmode", "e=2"), ("plasmode", "e=3")}
    assert all((row["joint_only"], row["by_block_only"]) == (30, 5) and row["p_loss"] < 0.05 / 4 for row in failing)


def test_no_allocation_gain() -> None:
    result = isb.decision(*_tables(overrides={("unequal", "gain"): (12, 10)}))
    assert result["verdict"] == "no_allocation_gain"
    assert all(not row["pass"] for row in result["failing"]["G"])


def test_allocation_gain_at_one_effect_suffices() -> None:
    simulation, plasmode, unequal = _tables()
    unequal["allocation_table.csv"][1] = {**unequal["allocation_table.csv"][1], "gains": "3", "losses": "3",
                                          "p_sign": str(isb._upper_binomial(3, 6))}
    assert isb.decision(simulation, plasmode, unequal)["verdict"] == "by_block_validated"


def test_precedence() -> None:
    both = {("plasmode", "anchor", "by_block"): 60, ("plasmode", "loss"): (30, 5)}
    assert isb.decision(*_tables(overrides=both))["verdict"] == "size_excess"
    fwer_and_power = {("plasmode", "outside", "e=2"): 60, ("plasmode", "loss"): (30, 5), ("unequal", "gain"): (5, 5)}
    assert isb.decision(*_tables(overrides=fwer_and_power))["verdict"] == "fwer_excess"
    power_and_gain = {("plasmode", "loss"): (30, 5), ("unequal", "gain"): (5, 5)}
    assert isb.decision(*_tables(overrides=power_and_gain))["verdict"] == "power_loss"


def test_reproduction_failure_gates_every_criterion() -> None:
    result = isb.decision(*_tables(reproduction=False, overrides={("plasmode", "loss"): (30, 5)}))
    assert result["verdict"] == "reproduction_failed"
    assert "criteria" not in result and "descriptive" not in result
    assert result["failing"]["reproduction"][0]["profile"] == "simulation"


def test_unavailable_cells_fail_their_criterion() -> None:
    simulation, plasmode, unequal = _tables()
    unequal["allocation_table.csv"] = []
    assert isb.decision(simulation, plasmode, unequal)["verdict"] == "no_allocation_gain"
    simulation, plasmode, unequal = _tables()
    plasmode["paired_table.csv"] = []
    assert isb.decision(simulation, plasmode, unequal)["verdict"] == "power_loss"


def test_descriptive_rates_do_not_change_the_verdict() -> None:
    first = isb.decision(*_tables())
    simulation, plasmode, unequal = _tables(overrides={("simulation", "n200", "magnitude", "by_block"): 400})
    simulation["minp_table.csv"][0]["minp_rate"] = "0.9"
    plasmode["block_table.csv"][0]["recall_mean"] = "0.01"
    second = isb.decision(simulation, plasmode, unequal, parent_recovery=[{"cell": "e=3", "recall_mean": "0.525"}])
    assert first["verdict"] == second["verdict"] == "by_block_validated"
    assert first["criteria"] == second["criteria"]
    assert first["descriptive"] != second["descriptive"]
    recall = second["descriptive"]["recall_by_block"][0]
    assert recall["parent"]["recall"] == 0.525 and recall["parent"]["source"] == isb.PARENT_SOURCE


def test_exact_binomial_tail() -> None:
    assert isb._upper_binomial(0, 0) == 1.0 and isb._upper_binomial(0, 5) == 1.0
    assert isb._upper_binomial(5, 5) == pytest.approx(1 / 32)
    assert isb._upper_binomial(30, 35) == pytest.approx(sum(
        __import__("math").comb(35, k) for k in range(30, 36)) / 2**35)


# ── sbatch ───────────────────────────────────────────────────────────────────


def test_sbatch_template_parses_and_dry_runs(tmp_path) -> None:
    script = REPO / "scripts" / "interaction_scan_blocks_check_array.sbatch"
    subprocess.run(["bash", "-n", str(script)], check=True)
    env = {"PATH": "/usr/bin:/bin:/run/current-system/sw/bin", "ISB_PROFILE": "plasmode_unequal",
           "ISB_OUT": str(tmp_path / "run"), "ISB_EXPORT": "/abs/cohort_30", "N_SHARDS": "50",
           "SLURM_ARRAY_TASK_ID": "7", "ISB_DRY_RUN": "1", "ISB_ARGS": "--limit 1"}
    shown = subprocess.run(["bash", str(script)], check=True, capture_output=True, text=True, env=env,
                           cwd=tmp_path).stdout
    assert "interaction_scan_blocks_check.py run --profile plasmode_unequal --export-dir /abs/cohort_30" in shown
    assert "--shard 7/50 --limit 1" in shown
    env["ISB_PROFILE"] = "other"
    refused = subprocess.run(["bash", str(script)], capture_output=True, text=True, env=env, cwd=tmp_path)
    assert refused.returncode == 2 and "ISB_PROFILE must be" in refused.stderr
