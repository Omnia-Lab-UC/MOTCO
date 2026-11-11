"""``scripts/production_latent_space_check.py``: settings, units, arms, the span plant, merge, reproduction, decision.

The simulation units run on the comparison's ``seaad`` design at fast CV knobs;
the plasmode profile runs on the plasmode tests' synthetic stand-in export. No
real data are read.
"""

from __future__ import annotations

import dataclasses
import json
import math
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


plsc = _script("production_latent_space_check")
lsc = plsc.lsc
pc = plsc.pc

REPO = Path(__file__).resolve().parents[1]
FAST_CV = {"cv1_splits": 3, "cv2_splits": 4, "n_repeats": 1, "max_components": 6, "random_state": 1203}


def _plasmode_tests():
    sys.path.insert(0, str(REPO))
    try:
        from tests import test_plasmode_check as module
    finally:
        sys.path.pop(0)
    return module


plasmode_tests = _plasmode_tests()


# ── Simulation settings, units and arms ──────────────────────────────────────


@pytest.fixture(scope="module")
def full_simulation_settings():
    return plsc.simulation_settings()


def test_simulation_units_and_settings(full_simulation_settings) -> None:
    settings = full_simulation_settings
    units = plsc.simulation_units(settings)
    assert len(units) == 3_000
    cells = [cell[0] for cell in settings["comparison"]["cells"]]
    assert cells == [*plsc.SHARED_CELLS, "magnitude_e0.25", "orient_e1.0_phi0", "shape_e1.0"]
    for design in lsc.DESIGNS:
        for cell in cells:
            assert sorted(u.replicate for u in units if (u.design, u.cell) == (design, cell)) == list(range(100))
    assert settings["permutations"] == 199 and settings["rejection"] == "p < alpha"
    assert set(settings["parents"]) == {"comparison_settings", "comparison_calibration", "comparison_rate_table"}
    assert settings["parents"]["comparison_rate_table"]["path"] == (
        "results/latent-space-comparison-2026-11-04/rate_table.csv"
    )
    assert settings["comparison"]["weak_delta"]["value"] == 0.11


def test_the_four_arms_differ_only_in_weighting_and_sparsity(full_simulation_settings) -> None:
    from motco.stats.sparse_pls import DEFAULT_KEEP_FRACTIONS

    settings = full_simulation_settings
    parent = json.loads((plsc.COMPARISON_DIR / "settings.json").read_text())["settings"]
    reference = lsc.arm_params("pls_cv", parent, seed=11, n_stages=3)
    production = plsc.simulation_arm_params("pls_cv", settings, seed=11, n_stages=3)
    assert production == reference
    for arm, (weighted, sparse) in plsc.ARM_FACTORS.items():
        params = plsc.simulation_arm_params(arm, settings, seed=11, n_stages=3)
        assert dataclasses.replace(params, integration_params={}) == dataclasses.replace(
            production, integration_params={}
        )
        extra = {key: value for key, value in params.integration_params.items()
                 if key not in production.integration_params}
        assert {key: value for key, value in params.integration_params.items() if key not in extra} == dict(
            production.integration_params
        )
        expected = {}
        if weighted:
            expected["block_weighting"] = "equal_variance"
        if sparse:
            expected["sparsity"] = {"scope": "block", "keep_fractions": list(DEFAULT_KEEP_FRACTIONS)}
        assert extra == expected
        assert params.span_capture is True and params.min_components == 2
        assert params.orientation_estimator == "net_displacement"


@pytest.mark.parametrize(
    "unit", [lsc.Unit("seaad", "orient_e0.25_phi0", 7), lsc.Unit("seaad_weak", "shape_e0.5", 42),
             lsc.Unit("n200", "translation", 3)]
)
def test_shared_cells_are_the_comparisons_units(full_simulation_settings, unit) -> None:
    settings = full_simulation_settings
    parent = json.loads((plsc.COMPARISON_DIR / "settings.json").read_text())["settings"]
    mine, theirs = lsc.unit_dataset(unit, settings["comparison"]), lsc.unit_dataset(unit, parent)
    for name in ("methylation", "expression", "proteomics", "metadata"):
        assert getattr(mine, name).equals(getattr(theirs, name))
    assert lsc.cell_spec(settings["comparison"], unit.cell) == lsc.cell_spec(parent, unit.cell)
    seed = lsc.generator_seed(2035, unit.design, unit.replicate)
    assert lsc._generator_params(settings["comparison"]["designs"][unit.design]["generator"], seed=seed) == (
        lsc._generator_params(parent["designs"][unit.design]["generator"], seed=seed)
    )


def test_added_cells_have_their_own_rrpp_seeds(full_simulation_settings) -> None:
    base = full_simulation_settings["comparison"]["base_seed"]
    for design in lsc.DESIGNS:
        shared = {lsc.rrpp_seed(base, design, cell, 5) for cell in plsc.SHARED_CELLS}
        added = {lsc.rrpp_seed(base, design, cell[0], 5) for cell in plsc.ADDED_CELLS}
        assert len(added) == 3 and not shared & added


def _parent_copy(tmp_path: Path) -> Path:
    target = tmp_path / "comparison"
    (target / "calibration").mkdir(parents=True)
    for name in ("settings.json", "rate_table.csv", "calibration/calibration.json"):
        shutil.copy(plsc.COMPARISON_DIR / name, target / name)
    return target


def test_a_changed_parent_is_refused_by_name(tmp_path) -> None:
    parent = _parent_copy(tmp_path)
    settings = plsc.simulation_settings(parent_dir=parent, replicates=2, designs=["seaad"])
    plsc.verify_parents(settings)
    (parent / "rate_table.csv").write_text((parent / "rate_table.csv").read_text() + "\n")
    with pytest.raises(SystemExit, match=r"Parent comparison_rate_table .*rate_table.csv\) changed"):
        plsc.run_shard(tmp_path / "run", settings, 0, 1, limit=0)
    assert not (tmp_path / "run").exists()


def test_a_calibration_the_comparison_did_not_use_is_refused(tmp_path) -> None:
    parent = _parent_copy(tmp_path)
    calibration = parent / "calibration" / "calibration.json"
    calibration.write_text(calibration.read_text() + " ")
    with pytest.raises(SystemExit, match="not the comparison's calibration_sha256"):
        plsc.simulation_settings(parent_dir=parent)


def test_a_resume_under_other_settings_is_refused_naming_the_keys(tmp_path) -> None:
    settings = plsc.simulation_settings(replicates=2, designs=["seaad"])
    plsc.run_shard(tmp_path, settings, 0, 1, limit=0)
    other = plsc.simulation_settings(replicates=3, designs=["seaad"], ladder=(1.0, 0.2))
    with pytest.raises(SystemExit, match=r"differing: .*arms\.sparse_block\.sparsity\.keep_fractions.*replicates"):
        plsc.run_shard(tmp_path, other, 0, 1, limit=0)


def _fake_records(unit, settings, *, p=0.5, failing_arm=None):
    profile = settings["profile"]
    if profile == "simulation":
        mode, effect = plsc.simulation_cell_specs()[unit.cell]
    else:
        effect = plsc.plasmode_effect(settings, unit.cell)
        mode = "none" if unit.cell == plsc.ANCHOR_CELL else plsc.PLANT_MODE
    records = []
    for arm in plsc.ARMS:
        base = {
            "profile": profile, "design": unit.design, "cell": unit.cell, "trajectory_mode": mode,
            "effect_size": effect, "replicate": unit.replicate, "settings_hash": lsc.settings_hash(settings),
            "arm": arm, "arm_seconds": 2.0, "generate_seconds": 1.0, "unit_runtime_seconds": 9.0,
        }
        if arm == failing_arm:
            records.append({**base, "status": "failed", "error_type": "ValueError", "error_message": "x"})
            continue
        record = {
            **base, "status": "completed", "rank": 3, "keep_fraction": 0.5 if "sparse" in arm else None,
            "support": [{"methylation": 10, "expression": 4}, {"methylation": 6, "expression": 2}]
            if "sparse" in arm else None,
            "cv_mean_auroc": 0.7, "pair_statistics": {"delta": 0.1, "angle": 5.0, "shape": 0.01},
            "p_values": {"delta": p, "angle": p, "shape": p},
            "span_capture": {"pooled_stage": 0.9, "group_difference": 0.8},
            "population_latent": {"delta": 0.2, "angle": 4.0, "shape": 0.1},
        }
        if profile == "plasmode" and arm == "pls_cv" and unit.cell == plsc.ANCHOR_CELL:
            record["fidelity"] = {"passed": True, "rank": 5, "cv_mean_auroc_rounded": 0.635556}
        records.append(record)
    return records


def test_resume_runs_only_the_remaining_units(tmp_path, monkeypatch) -> None:
    settings = plsc.simulation_settings(replicates=2, designs=["seaad"])
    calls = []

    def fake(unit, settings):
        calls.append(unit.key)
        return _fake_records(unit, settings)

    monkeypatch.setattr(plsc, "evaluate_simulation_unit", fake)
    assert plsc.run_shard(tmp_path, settings, 0, 2, limit=3) == 3
    path = lsc.shard_path(tmp_path / "simulation", 0, 2)
    with path.open("a") as handle:
        handle.write('{"design": "seaad", "cell": "anchor", "rep')  # a torn line
    first = list(calls)
    assert plsc.run_shard(tmp_path, settings, 0, 2) == 10 - 3
    assert not set(first) & set(calls[3:])
    assert plsc.run_shard(tmp_path, settings, 1, 2) == 10
    assert plsc.run_shard(tmp_path, settings, 0, 2) == 0
    summary = plsc.merge(tmp_path, "simulation")
    assert summary["complete_units"] == 20 and summary["missing_units"] == 0
    merged = [json.loads(line) for line in (tmp_path / "simulation" / "merged.jsonl").read_text().splitlines()]
    keys = [(plsc.record_key(record), record["arm"]) for record in merged]
    assert len(keys) == len(set(keys)) == 20 * 4


# ── One simulation unit ──────────────────────────────────────────────────────


def _small_settings(**overrides):
    settings = plsc.simulation_settings(designs=["seaad"], replicates=2, permutations=19)
    settings["comparison"]["pls_cv_params"] = dict(FAST_CV)
    settings.update(overrides)
    return settings


UNIT = lsc.Unit("seaad", "shape_e0.5", 1)


@pytest.fixture(scope="module")
def simulation_unit():
    settings = _small_settings()
    return settings, plsc.evaluate_simulation_unit(UNIT, settings)


def test_simulation_records_carry_every_field(simulation_unit) -> None:
    settings, records = simulation_unit
    assert [record["arm"] for record in records] == list(plsc.ARMS)
    for record in records:
        assert record["status"] == "completed", record.get("error_message")
        assert record["generator_seed"] == lsc.generator_seed(2035, "seaad", 1)
        assert record["rrpp_seed"] == lsc.rrpp_seed(2035, "seaad", "shape_e0.5", 1)
        assert {"rank", "keep_fraction", "support", "cv_mean_auroc", "pair_statistics", "p_values", "span_capture",
                "population_latent", "arm_seconds", "generate_seconds", "unit_runtime_seconds"} <= set(record)
        assert set(record["p_values"]) == {"delta", "angle", "shape"}
        sparse = "sparse" in record["arm"]
        assert (record["keep_fraction"] is not None) == sparse
        if sparse:
            assert len(record["support"]) == record["rank"]
            assert set(record["support"][0]) == {"methylation", "expression"}
        assert record["block_weighting"] == ("equal_variance" if "weighted" in record["arm"] else None)


def test_production_arm_equals_the_comparisons_unit(simulation_unit) -> None:
    settings, records = simulation_unit
    comparison = {**settings["comparison"], "arms": ["pls_cv"]}
    theirs = lsc.evaluate_unit(UNIT, comparison)[0]
    mine = records[0]
    assert mine["rank"] == theirs["retained_dimensionality"]
    for key in ("pair_statistics", "p_values", "span_capture", "population_latent", "cv_mean_auroc"):
        assert mine[key] == theirs[key], key


def test_each_candidate_equals_its_own_evaluation_from_scratch(simulation_unit) -> None:
    from motco.simulations.evaluation import evaluate_semisynthetic_trajectory

    settings, records = simulation_unit
    dataset = lsc.unit_dataset(UNIT, settings["comparison"])
    seed = lsc.rrpp_seed(2035, "seaad", "shape_e0.5", 1)
    for record in records[1:]:
        params = plsc.simulation_arm_params(record["arm"], settings, seed=seed, n_stages=3)
        expected = plsc.simulation_arm_record(evaluate_semisynthetic_trajectory(dataset, params))
        for key, value in expected.items():
            assert record[key] == value, (record["arm"], key)


def test_a_failing_arm_is_recorded_without_aborting(monkeypatch) -> None:
    settings = _small_settings()
    original = plsc.simulation_arm_params

    def failing(arm, settings, **kwargs):
        if arm == "weighted":
            raise ValueError("planted failure")
        return original(arm, settings, **kwargs)

    monkeypatch.setattr(plsc, "simulation_arm_params", failing)
    records = plsc.evaluate_simulation_unit(lsc.Unit("seaad", "anchor", 0), settings)
    status = {record["arm"]: record["status"] for record in records}
    assert status == {"pls_cv": "completed", "weighted": "failed", "sparse_block": "completed",
                      "weighted_sparse_block": "completed"}
    assert records[1]["error_message"] == "planted failure"


# ── Plasmode profile ─────────────────────────────────────────────────────────


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
    parent = pc.build_settings(export, fidelity, replicates=2, permutations=9, bootstrap=5, e_star=3.0,
                               descriptive=2.0)
    parent["cv"] = dict(plasmode_tests.SMALL_CV)
    parent_dir = tmp_path_factory.mktemp("plasmode")
    (parent_dir / "settings.json").write_text(json.dumps({"hash": lsc.settings_hash(parent), "settings": parent}))
    return directory, export, parent_dir


def _plasmode_settings(stand_in, **kwargs):
    _, export, parent_dir = stand_in
    kwargs.setdefault("replicates", 2)
    kwargs.setdefault("permutations", 9)
    return plsc.plasmode_settings(export, parent_dir=parent_dir, **kwargs)


def test_plasmode_units_and_settings(stand_in) -> None:
    settings = _plasmode_settings(stand_in, replicates=100)
    units = plsc.plasmode_units(settings)
    assert len(units) == 400
    assert [unit.key for unit in units[:5]] == [
        ("plasmode", "anchor", 0), ("plasmode", "span_e0.25", 0), ("plasmode", "span_e0.5", 0),
        ("plasmode", "span_e1", 0), ("plasmode", "anchor", 1),
    ]
    assert settings["rejection"] == "p <= alpha" and settings["attribution"] is False
    assert set(settings["parents"]) == {"plasmode_settings"}
    lever = plsc.plasmode_units(_plasmode_settings(stand_in, replicates=100, plant_replicates=50))
    assert len(lever) == 100 + 3 * 50


def test_an_export_whose_hashes_differ_is_refused(stand_in, tmp_path) -> None:
    directory, _, parent_dir = stand_in
    copy = shutil.copytree(directory, tmp_path / "export")
    lines = (copy / "rna.csv").read_text().splitlines()
    (copy / "rna.csv").write_text("\n".join(lines[:-1] + [lines[-1].replace(",", ",1", 1)]) + "\n")
    other = pc.load_export(copy)
    with pytest.raises(SystemExit, match=r"differ from the settings' for: rna.csv; refusing"):
        plsc.plasmode_settings(other, parent_dir=parent_dir, replicates=2)
    settings = _plasmode_settings(stand_in)
    with pytest.raises(SystemExit, match=r"for: rna.csv; refusing"):
        plsc.run_shard(tmp_path / "run", settings, 0, 1, export=other, limit=0)


def test_the_plasmode_arms_differ_only_in_weighting_and_sparsity(stand_in) -> None:
    settings = _plasmode_settings(stand_in)
    seeds = pc.unit_seeds(2034, "paper", 3)
    production = plsc.plasmode_arm_params("pls_cv", settings, seeds)
    parent = {**settings["plasmode"], "attribution": None}
    assert production == pc.analysis_params(parent, seeds, permutations=9)
    assert production.residual_correction == "df" and production.min_components == 2
    assert production.orientation_estimator == "net_displacement" and not production.attribution.enabled
    for arm, (weighted, sparse) in plsc.ARM_FACTORS.items():
        params = plsc.plasmode_arm_params(arm, settings, seeds)
        assert dataclasses.replace(params, block_weighting=None, sparsity=None) == production
        assert params.block_weighting == ("equal_variance" if weighted else None)
        assert (params.sparsity is not None) == sparse
        if sparse:
            assert params.sparsity.scope == "block"


def test_the_plasmode_labels_are_the_plasmodes_draws(stand_in) -> None:
    _, export, _ = stand_in
    settings = _plasmode_settings(stand_in)
    for replicate in (0, 1):
        mine = plsc.build_plasmode_unit(lsc.Unit("plasmode", "span_e0.5", replicate), settings, export)
        theirs = pc.build_unit(0.0, replicate, settings["plasmode"], export)
        assert mine["metadata"].equals(theirs["metadata"])
        assert mine["seeds"] == theirs["seeds"]


def _centred_stage_means(values, stages, rows, stage_order):
    means = np.vstack([values[rows & (stages == stage)].mean(axis=0) for stage in stage_order])
    return means - means.mean(axis=0)


@pytest.mark.parametrize("effect", [0.25, 1.0])
def test_span_plant_geometry(stand_in, effect) -> None:
    from motco.simulations.preprocessing import fit_block_scalers

    _, export, _ = stand_in
    labels, _ = pc.draw_labels(export.real_labels, export.stages, export.stage_order, 2034, "paper", 4)
    planted = plsc.span_plant(export.blocks, export.stages, export.stage_order, labels, effect, group="Male")
    stages = export.stages.astype(str)
    male = labels.astype(str) == "Male"
    everyone = np.ones_like(male)
    scalers = fit_block_scalers(dict(export.blocks), {})
    for name, frame in export.blocks.items():
        scaler = scalers[name]
        before = (frame.to_numpy(dtype=float) - scaler.mean) / scaler.scale
        after = (planted[name].to_numpy(dtype=float) - scaler.mean) / scaler.scale
        expected = _centred_stage_means(before, stages, male, export.stage_order) + effect * _centred_stage_means(
            before, stages, everyone, export.stage_order
        )
        np.testing.assert_allclose(_centred_stage_means(after, stages, male, export.stage_order), expected,
                                   rtol=0.0, atol=1e-8)
        assert np.array_equal(planted[name].to_numpy()[~male], frame.to_numpy()[~male])
        # In the block's units the shift is e · (pooled stage mean − mean of the stage means), scaler-free.
        raw = frame.to_numpy(dtype=float)
        pooled = _centred_stage_means(raw, stages, everyone, export.stage_order)
        shift = planted[name].to_numpy(dtype=float) - raw
        for index, stage in enumerate(export.stage_order):
            rows = male & (stages == stage)
            np.testing.assert_allclose(shift[rows], np.broadcast_to(effect * pooled[index], shift[rows].shape),
                                       rtol=0.0, atol=1e-8)
        assert list(planted[name].index) == list(frame.index)
        assert list(planted[name].columns) == list(frame.columns)


def test_the_anchor_is_the_exports_blocks(stand_in) -> None:
    _, export, _ = stand_in
    labels, _ = pc.draw_labels(export.real_labels, export.stages, export.stage_order, 2034, "paper", 0)
    blocks = plsc.span_plant(export.blocks, export.stages, export.stage_order, labels, 0.0, group="Male")
    assert all(blocks[name] is export.blocks[name] for name in export.blocks)
    built = plsc.build_plasmode_unit(lsc.Unit("plasmode", "anchor", 0), _plasmode_settings(stand_in), export)
    assert all(built["blocks"][name] is export.blocks[name] for name in export.blocks)


@pytest.fixture(scope="module")
def plasmode_anchor(stand_in):
    _, export, _ = stand_in
    settings = _plasmode_settings(stand_in)
    return settings, plsc.evaluate_plasmode_unit(lsc.Unit("plasmode", "anchor", 1), settings, export)


def test_plasmode_records_hold_seeds_and_statistics_only(stand_in, plasmode_anchor) -> None:
    _, export, _ = stand_in
    _, records = plasmode_anchor
    assert [record["arm"] for record in records] == list(plsc.ARMS)
    text = json.dumps(records)
    for donor in export.metadata[export.sample_col]:
        assert donor not in text
    assert "Female" not in text and "Male" not in text
    for record in records:
        assert record["status"] == "completed", record.get("error_message")
        assert record["seeds"] == pc.unit_seeds(2034, "paper", 1)
        assert record["rrpp_seed"] == record["seeds"]["rrpp"]
        assert not any(isinstance(value, list) and len(value) == len(export.metadata) for value in record.values())
    assert records[0]["fidelity"]["passed"] is True
    assert all("fidelity" not in record for record in records[1:])
    assert records[1]["block_weights"] == {name: 1.0 / math.sqrt(30) for name in pc.BLOCKS}
    assert set(records[2]["support"][0]) == set(pc.BLOCKS)


def test_the_plasmode_production_arm_equals_a_plain_analysis(stand_in, plasmode_anchor) -> None:
    _, export, _ = stand_in
    settings, records = plasmode_anchor
    built = plsc.build_plasmode_unit(lsc.Unit("plasmode", "anchor", 1), settings, export)
    for record in (records[0], records[3]):
        params = plsc.plasmode_arm_params(record["arm"], settings, built["seeds"])
        _, result = pc.analyse(built["blocks"], built["metadata"], export, params)
        expected = json.loads(json.dumps(lsc._jsonable(plsc.plasmode_arm_record(result))))
        assert {key: record[key] for key in expected} == expected


def test_a_fidelity_mismatch_is_recorded_not_raised(stand_in) -> None:
    _, export, _ = stand_in
    settings = _plasmode_settings(stand_in)
    wrong = {**settings, "fidelity": {"selected_lv": 19, "cv_mean_auroc": 0.5}}
    records = plsc.evaluate_plasmode_unit(lsc.Unit("plasmode", "anchor", 0), wrong, export)
    assert records[0]["status"] == "completed"
    assert records[0]["fidelity"]["passed"] is False and "anchor_latent" in records[0]["fidelity"]["message"]
    planted = plsc.evaluate_plasmode_unit(lsc.Unit("plasmode", "span_e0.25", 0), wrong, export)
    assert all("fidelity" not in record for record in planted)


# ── Merge and reproduction ───────────────────────────────────────────────────


def _records(settings, *, p=0.5, overrides=None):
    overrides = overrides or {}
    out = []
    for unit in plsc.enumerate_units(settings):
        out.extend(overrides.get(unit.key, lambda u, s: _fake_records(u, s, p=p))(unit, settings))
    return out


def test_merge_tables_count_failures_and_pairs() -> None:
    settings = plsc.simulation_settings(replicates=3, designs=["seaad"])

    def unit_records(unit, settings):
        records = _fake_records(unit, settings, failing_arm="weighted" if unit.replicate == 2 else None)
        for record in records:
            if record["status"] == "completed":
                reject = (record["arm"] == "pls_cv") == (unit.replicate == 0)
                record["p_values"] = {"delta": 0.01 if reject else 0.5, "angle": 0.05, "shape": None}
        return records

    rows = plsc.rate_table(_records(settings, overrides={u.key: unit_records for u in plsc.simulation_units(settings)}),
                           settings)
    assert set(rows[0]) == set(plsc.RATE_COLUMNS)
    by = {(r["cell"], r["arm"], r["statistic"]): r for r in rows}
    assert by[("anchor", "pls_cv", "delta")]["rejections"] == 1 and by[("anchor", "pls_cv", "delta")]["count"] == 3
    assert by[("anchor", "weighted", "delta")]["count"] == 2 and by[("anchor", "weighted", "delta")]["failures"] == 1
    assert by[("anchor", "weighted", "delta")]["rejections"] == 1
    assert by[("anchor", "pls_cv", "angle")]["rejections"] == 0  # p < alpha: 0.05 does not reject
    assert by[("anchor", "pls_cv", "shape")]["count"] == 0 and by[("anchor", "pls_cv", "shape")]["rate"] is None

    records = _records(settings, overrides={u.key: unit_records for u in plsc.simulation_units(settings)})
    paired = {(r["cell"], r["arm"], r["statistic"]): r for r in plsc.paired_table(records, settings)}
    row = paired[("anchor", "weighted", "delta")]
    assert set(row) == set(plsc.PAIRED_COLUMNS)
    assert (row["n_paired"], row["b"], row["c"]) == (2, 1, 1)  # replicate 2 has no weighted reading
    row = paired[("anchor", "sparse_block", "delta")]
    assert (row["n_paired"], row["b"], row["c"]) == (3, 1, 2)
    assert row["p_loss"] == pytest.approx(7 / 8) and row["p_gain"] == pytest.approx(0.5)
    none = paired[("anchor", "sparse_block", "angle")]
    assert (none["b"], none["c"], none["p_loss"], none["p_gain"], none["p_two_sided"]) == (0, 0, 1.0, 1.0, 1.0)

    cells = {(r["cell"], r["arm"]): r for r in plsc.cell_table(records, settings)}
    assert set(cells[("anchor", "pls_cv")]) == set(plsc.cell_columns(plsc.SIMULATION_BLOCKS))
    sparse = cells[("magnitude", "sparse_block")]
    assert sparse["keep_fraction_median"] == 0.5 and sparse["support_methylation_median"] == 8.0
    assert sparse["target_statistic"] == "delta" and sparse["population_target_median"] == 0.2
    assert cells[("magnitude", "pls_cv")]["support_methylation_median"] is None
    assert cells[("anchor", "weighted")]["n"] == 2

    summary = plsc.accounting(records, settings, expected=30, duplicates=0)
    assert summary["complete_units"] == 30 and summary["failed_records"] == 10
    assert summary["failures_by_arm_and_cell"]["weighted"]["seaad/anchor"] == 1
    assert summary["recorded_core_hours"] == pytest.approx(30 * 9.0 / 3600)
    assert summary["arm_seconds_median"]["pls_cv"] == 2.0


def test_duplicates_across_shards_are_counted(tmp_path, monkeypatch) -> None:
    settings = plsc.simulation_settings(replicates=2, designs=["seaad"])
    monkeypatch.setattr(plsc, "evaluate_simulation_unit", lambda unit, settings: _fake_records(unit, settings))
    plsc.run_shard(tmp_path, settings, 0, 1)
    plsc.run_shard(tmp_path, settings, 0, 2)
    summary = plsc.merge(tmp_path, "simulation")
    assert summary["complete_units"] == 20 and summary["duplicate_records"] == 10 * 4


def _committed(rates, *, bump=None):
    rows = []
    for row in rates:
        if row["arm"] != "pls_cv" or row["cell"] not in plsc.SHARED_CELLS:
            continue
        entry = {key: row[key] for key in ("design", "cell", "statistic", "count", "rejections", "failures")}
        entry["arm"] = "pls_cv"
        if bump == (row["design"], row["cell"], row["statistic"]):
            entry["rejections"] = int(entry["rejections"]) + 1
        rows.append({key: str(value) for key, value in entry.items()})
    return rows


def test_simulation_reproduction_both_outcomes() -> None:
    settings = plsc.simulation_settings(replicates=2)
    rates = plsc.rate_table(_records(settings, p=0.01), settings)
    passed = plsc.simulation_reproduction(rates, settings, _committed(rates))
    assert passed["passed"] and passed["rows_compared"] == 63 and passed["rows_matching"] == 63
    failed = plsc.simulation_reproduction(rates, settings, _committed(rates, bump=("n200", "translation", "angle")))
    assert not failed["passed"] and failed["rows_matching"] == 62
    assert [(r["design"], r["cell"], r["statistic"]) for r in failed["differing"]] == [("n200", "translation", "angle")]
    assert failed["differing"][0]["committed"]["rejections"] == failed["differing"][0]["reproduced"]["rejections"] + 1


def test_plasmode_reproduction_outcomes(stand_in) -> None:
    settings = _plasmode_settings(stand_in, replicates=3)
    records = _records(settings)
    passed = plsc.plasmode_reproduction(records, settings)
    assert passed["passed"] and passed["units_checked"] == 3 and passed["units_passed"] == 3
    anchor = [r for r in records if r["cell"] == "anchor" and r["arm"] == "pls_cv"]
    anchor[1]["fidelity"] = {"passed": False, "rank": 4, "cv_mean_auroc_rounded": 0.6}
    failed = plsc.plasmode_reproduction(records, settings)
    assert not failed["passed"] and failed["differing"] == [
        {"unit": plsc.unit_key("plasmode", "anchor", 1), "rank": 4, "cv_mean_auroc_rounded": 0.6}
    ]
    assert plsc.accounting(records, settings, expected=12, duplicates=0)["fidelity_failures"] == 1
    anchor[1]["fidelity"]["passed"] = True
    anchor[2].update(status="failed")
    unavailable = plsc.plasmode_reproduction(records, settings)
    assert not unavailable["passed"] and unavailable["unavailable"] == [plsc.unit_key("plasmode", "anchor", 2)]


def test_merge_and_decide_end_to_end(stand_in, tmp_path, monkeypatch) -> None:
    _, export, _ = stand_in
    simulation = plsc.simulation_settings(replicates=2)
    plasmode = _plasmode_settings(stand_in)
    monkeypatch.setattr(plsc, "evaluate_simulation_unit", lambda unit, settings: _fake_records(unit, settings))
    monkeypatch.setattr(plsc, "evaluate_plasmode_unit", lambda unit, settings, export: _fake_records(unit, settings))
    assert plsc.run_shard(tmp_path, simulation, 0, 1) == 60
    assert plsc.run_shard(tmp_path, plasmode, 0, 1, export=export) == 8
    assert plsc.merge(tmp_path, "simulation")["reproduction_passed"] is False  # fake rates, not the committed ones
    assert plsc.merge(tmp_path, "plasmode")["reproduction_passed"] is True
    for profile in plsc.PROFILES:
        for name in ("rate_table.csv", "paired_table.csv", "cell_table.csv", "accounting.json", "reproduction.json",
                     "merged.jsonl"):
            assert (tmp_path / profile / name).exists()
    result = plsc.decide(tmp_path)
    assert result["verdict"] == "reproduction_failed"
    assert set(result["arms"]) == set(plsc.CANDIDATE_ARMS)
    assert all(set(entry["criteria"]) == {"S1", "S2", "S3", "G"} for entry in result["arms"].values())
    assert (tmp_path / plsc.DECISION_NAME).exists()


# ── Decision ─────────────────────────────────────────────────────────────────


def _tables(*, control=5, off_target=5, large_off_target=17, count=100, rate_overrides=None, paired_overrides=None):
    """Synthetic merged tables for every cell the plan reads; overrides keyed (design, cell, arm, statistic)."""

    rate_overrides, paired_overrides = rate_overrides or {}, paired_overrides or {}
    specs = plsc.simulation_cell_specs()
    cells = [(design, cell, *specs[cell]) for design in lsc.DESIGNS for cell in specs]
    cells += [(plsc.PLASMODE_DESIGN, name, "none" if name == "anchor" else "magnitude", effect)
              for name, effect in plsc.plasmode_cells()]
    tables = {"simulation": {"rates": [], "paired": [], "cells": []}, "plasmode": {"rates": [], "paired": [],
                                                                                   "cells": []}}
    for design, cell, mode, effect in cells:
        profile = "plasmode" if design == plsc.PLASMODE_DESIGN else "simulation"
        for statistic in plsc.STATISTICS:
            target = plsc.TARGETS.get(mode) == statistic
            for arm in plsc.ARMS:
                if target:
                    rejections = 50
                elif mode in plsc.TARGETS and effect > 0.5:
                    rejections = large_off_target
                else:
                    rejections = control if mode in ("none", "translation") else off_target
                rejections = rate_overrides.get((design, cell, arm, statistic), rejections)
                tables[profile]["rates"].append({
                    "design": design, "cell": cell, "trajectory_mode": mode, "effect_size": effect, "arm": arm,
                    "statistic": statistic, "count": count, "rejections": rejections, "rate": rejections / count,
                    "failures": 0,
                })
                if arm == "pls_cv":
                    continue
                b, c = 0, 0
                row = {"design": design, "cell": cell, "trajectory_mode": mode, "effect_size": effect,
                       "statistic": statistic, "arm": arm, "n_paired": count, "rate_reference": 0.5,
                       "rate_candidate": 0.5, "difference": 0.0, "b": b, "c": c,
                       **plsc.paired_p_values(b, c)}
                row.update(paired_overrides.get((design, cell, arm, statistic), {}))
                tables[profile]["paired"].append(row)
            for arm in plsc.ARMS:
                if statistic == "delta":
                    tables[profile]["cells"].append({"design": design, "cell": cell, "arm": arm, "n": count,
                                                     "rank_median": 3})
    return tables


PASSED = {"simulation": {"passed": True}, "plasmode": {"passed": True}}
GAIN = {"b": 2, "c": 16, **plsc.paired_p_values(2, 16)}


def _gain_for(*arms):
    return {("seaad_weak", "orient_e1.0_phi0", arm, "angle"): GAIN for arm in arms}


def _decide(tables, reproduction=PASSED):
    return plsc.decision(tables["simulation"], tables["plasmode"], reproduction=reproduction)


def test_rule_constants() -> None:
    assert plsc.paired_p_values(0, 0) == {"p_loss": 1.0, "p_gain": 1.0, "p_two_sided": 1.0}
    assert plsc.loss_allowed(16) == 2
    assert GAIN["p_gain"] <= 0.05 / 11


def test_proposal_names_the_preferred_qualifying_arm() -> None:
    tables = _tables(paired_overrides=_gain_for("weighted", "weighted_sparse_block"))
    result = _decide(tables)
    assert result["verdict"] == "propose_paper_grade" and result["named_arm"] == "weighted"
    assert result["qualifying_arms"] == ["weighted", "weighted_sparse_block"]
    assert result["arms"]["sparse_block"]["verdict"] == "no_gain"
    for arm in plsc.CANDIDATE_ARMS:
        assert set(result["arms"][arm]["criteria"]) == {"S1", "S2", "S3", "G"}
    gain = result["arms"]["weighted"]["criteria"]["G"]
    assert gain["holds"] and [(g["cell"], g["b"], g["c"]) for g in gain["gains"]] == [("orient_e1.0_phi0", 2, 16)]
    assert gain["pooled"]["c"] - gain["pooled"]["b"] == 14


def test_without_a_gain_production_is_kept() -> None:
    result = _decide(_tables())
    assert result["verdict"] == "keep_production" and result["named_arm"] is None
    assert {arm: entry["verdict"] for arm, entry in result["arms"].items()} == dict.fromkeys(plsc.CANDIDATE_ARMS,
                                                                                            "no_gain")
    observations = result["arms"]["weighted"]["criteria"]["G"]["observations"]
    assert len(observations) == 11 and all(row["status"] == "no_gain" for row in observations)


def test_gain_at_exactly_the_bonferroni_threshold() -> None:
    threshold = 0.05 / 11
    exact = {("plasmode", "span_e0.5", "weighted", "delta"): {"p_gain": threshold, "b": 1, "c": 12}}
    assert _decide(_tables(paired_overrides=exact))["arms"]["weighted"]["criteria"]["G"]["holds"]
    above = {("plasmode", "span_e0.5", "weighted", "delta"): {"p_gain": threshold * 1.0001, "b": 1, "c": 12}}
    assert not _decide(_tables(paired_overrides=above))["arms"]["weighted"]["criteria"]["G"]["holds"]


def test_size_excess_precedes_a_failing_s2_and_names_the_readings() -> None:
    losses = {("n200", cell, "weighted", stat): {"p_loss": 0.001, "b": 12, "c": 0}
              for cell, stat in (("magnitude", "delta"), ("shape_e0.25", "shape"), ("orient_e0.5_phi0", "angle"))}
    tables = _tables(rate_overrides={("seaad", "anchor", "weighted", "angle"): 25},
                     paired_overrides={**losses, **_gain_for("weighted")})
    result = _decide(tables)
    weighted = result["arms"]["weighted"]
    assert weighted["verdict"] == "size_excess"
    assert weighted["failing_criteria"] == ["S1", "S2"]
    failing = weighted["failing_readings"]["S1"]
    assert [(r["family"], r["design"], r["cell"], r["statistic"], r["status"]) for r in failing] == [
        ("simulation", "seaad", "anchor", "angle", "above_ceiling")
    ]
    assert {(r["cell"], r["status"]) for r in weighted["failing_readings"]["S2"]} == {
        ("magnitude", "below_ceiling"), ("shape_e0.25", "below_ceiling"), ("orient_e0.5_phi0", "below_ceiling")
    }


def test_one_chance_exceedance_is_allowed_and_the_ceiling_is_not() -> None:
    chance = _decide(_tables(rate_overrides={("n200", "translation", "weighted", "delta"): 11},
                             paired_overrides=_gain_for("weighted")))
    family = chance["arms"]["weighted"]["criteria"]["S1"]["families"]["simulation"]
    assert family["holds"] and family["observed_exceedances"] == 1 and family["allowed_exceedances"] >= 1
    assert chance["arms"]["weighted"]["verdict"] == "qualifies"
    plasmode = _decide(_tables(rate_overrides={("plasmode", "anchor", "weighted", "shape"): 30},
                               paired_overrides=_gain_for("weighted")))
    assert plasmode["arms"]["weighted"]["verdict"] == "size_excess"
    assert not plasmode["arms"]["weighted"]["criteria"]["S1"]["families"]["plasmode"]["holds"]


def test_power_loss_counts_against_the_allowed_number() -> None:
    def losses(n, p):
        cells = [("n200", "magnitude", "delta"), ("seaad", "shape_e0.5", "shape"),
                 ("n200", "orient_e1.0_phi0", "angle")]
        return {(d, c, "sparse_block", s): {"p_loss": p, "b": 9, "c": 1} for d, c, s in cells[:n]}

    two = _decide(_tables(paired_overrides={**losses(2, 0.01), **_gain_for("sparse_block")}))
    assert two["arms"]["sparse_block"]["verdict"] == "qualifies"
    assert two["arms"]["sparse_block"]["criteria"]["S2"]["observed_losses"] == 2
    three = _decide(_tables(paired_overrides={**losses(3, 0.01), **_gain_for("sparse_block")}))
    assert three["arms"]["sparse_block"]["verdict"] == "power_loss"
    assert three["arms"]["sparse_block"]["criteria"]["S2"]["m"] == 16


def test_specificity_reads_e_up_to_half_and_reports_the_large_plants() -> None:
    result = _decide(_tables(large_off_target=30, paired_overrides=_gain_for("weighted")))
    assert result["arms"]["weighted"]["verdict"] == "qualifies"
    assert result["arms"]["weighted"]["criteria"]["S3"]["m"] == 24
    large = result["descriptive"]["large_plant_off_target"]
    assert {(row["cell"], row["statistic"]) for row in large} == {
        ("orient_e1.0_phi0", "delta"), ("orient_e1.0_phi0", "shape"), ("shape_e1.0", "delta"), ("shape_e1.0", "angle")
    }
    assert len(large) == 2 * 4 * 4 and all(row["rate"] == 0.3 for row in large)
    excess = _decide(_tables(rate_overrides={("seaad", "shape_e0.5", "weighted", "delta"): 30},
                             paired_overrides=_gain_for("weighted")))
    assert excess["arms"]["weighted"]["verdict"] == "specificity_excess"


def test_the_gate_still_records_every_criterion() -> None:
    tables = _tables(paired_overrides=_gain_for("weighted"))
    result = _decide(tables, reproduction={"simulation": {"passed": False}, "plasmode": {"passed": True}})
    assert result["verdict"] == "reproduction_failed" and result["named_arm"] is None
    assert result["gate"]["failures"] == [{"profile": "simulation", "reason": "reproduction failed"}]
    assert result["arms"]["weighted"]["verdict"] == "qualifies"
    assert result["qualifying_arms"] == ["weighted"]


def test_an_unavailable_reading_fails_its_criterion() -> None:
    tables = _tables(paired_overrides=_gain_for("weighted"))
    tables["simulation"]["paired"] = [
        row for row in tables["simulation"]["paired"]
        if not (row["design"] == "n200" and row["cell"] == "shape_e1.0" and row["arm"] == "weighted")
    ]
    result = _decide(tables)
    assert result["arms"]["weighted"]["verdict"] == "power_loss"
    assert result["arms"]["weighted"]["criteria"]["S2"]["n_unavailable"] == 1


# ── sbatch ───────────────────────────────────────────────────────────────────


def test_sbatch_template_parses_and_dry_runs(tmp_path) -> None:
    script = REPO / "scripts" / "production_latent_space_array.sbatch"
    subprocess.run(["bash", "-n", str(script)], check=True)
    env = {"PATH": "/usr/bin:/bin", "PLS_DRY_RUN": "1", "PLS_PROFILE": "plasmode", "PLS_OUT": str(tmp_path / "run"),
           "PLS_EXPORT": str(tmp_path / "export"), "N_SHARDS": "100", "SLURM_ARRAY_TASK_ID": "3",
           "PLS_ARGS": "--plant-reps 50"}
    shell = shutil.which("bash") or "bash"
    out = subprocess.run([shell, str(script)], env=env, cwd=REPO, check=True, capture_output=True, text=True).stdout
    assert out.split() == ["python", "scripts/production_latent_space_check.py", "run", "--profile", "plasmode",
                           "--export-dir", str(tmp_path / "export"), "--out-dir", str(tmp_path / "run"), "--shard",
                           "3/100", "--plant-reps", "50"]
    assert not (tmp_path / "run").exists()
