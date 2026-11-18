"""``scripts/adni_gain_specificity_check.py``: plants, the simulation profile, checks, merge and the decision.

The plants profile runs on the anchor test's synthetic stand-in export (no real
data are read), the simulation profile on a small registered table, both at
small CV knobs and few permutations. The decision tests run on synthetic
tables at the committed replicates.
"""

from __future__ import annotations

import copy
import dataclasses
import json
import math
import shlex
import shutil
import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from tests.test_adni_anchor_check import (
    ADNI_TABLE,
    SMALL_CV,
    STAND_IN_FEATURES,
    _block,
    _copy_export,
    _rewrite,
    _snapshot,
    _stand_in_contract,
    _walk,
    write_export,
)


def _script(name: str):
    sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
    try:
        return __import__(name)
    finally:
        sys.path.pop(0)


gain = _script("adni_gain_specificity_check")
aac, adc, lsc, plsc, psc, grade = gain.aac, gain.adc, gain.lsc, gain.plsc, gain.psc, gain.grade
d10 = _script("covariate_type_i_check")

REPO = Path(__file__).resolve().parents[1]
PERMUTATIONS = 9
#: Half the stand-in's 176 participants, so the held-out reference is as large as the sample.
STAND_IN_POINT = 88
POINT = f"n{STAND_IN_POINT}"
SIM_TABLE = "gain_test"
SIM_CV = {"cv1_splits": 2, "cv2_splits": 3, "n_repeats": 1, "max_components": 5, "random_state": 1203}
adc.TABLES.setdefault(SIM_TABLE, {
    "group_stage_sizes": [[12, 12, 12], [12, 12, 12]],
    "phase_shares": [[0.6, 0.6, 0.4], [0.55, 0.5, 0.5]],
    "source": {"table": "tests/test_adni_gain_specificity_check.py"},
    "sha256": {},
    "reference": None,
})
SIM_POINT = gain.simulation_point(SIM_TABLE)


# ── Fixtures ─────────────────────────────────────────────────────────────────


@pytest.fixture(scope="module")
def export_dir(tmp_path_factory) -> Path:
    return write_export(tmp_path_factory.mktemp("export") / "motco_input")


@pytest.fixture(scope="module")
def fidelity(export_dir) -> dict:
    """The stand-in's own cohort fit (label-free), as ADNI's 8 and 0.599514."""

    cohort = aac.load_cohort(export_dir, _stand_in_contract(fidelity_rank=2, fidelity_auroc=0.5))
    settings = gain.plants_settings(cohort, control_replicates=2, plant_replicates=2, permutations=PERMUTATIONS,
                                    point_sizes=(STAND_IN_POINT,), cv=SMALL_CV)
    built = aac.build_unit(cohort, settings, gain.COHORT_POINT, 0)
    dataset = aac.make_dataset(built["blocks"], built["metadata"], cohort.contract)
    latent = aac._analyse(dataset, aac.analysis_params(settings, 1, permutations=0)).record["latent"]
    return {"rank": int(latent["selected_lv"]), "auroc": round(float(latent["cv_mean_auroc"]), 6)}


@pytest.fixture(scope="module")
def contract(fidelity) -> aac.CohortContract:
    return _stand_in_contract(fidelity_rank=fidelity["rank"], fidelity_auroc=fidelity["auroc"])


@pytest.fixture(scope="module")
def cohort(export_dir, contract):
    return aac.load_cohort(export_dir, contract)


@pytest.fixture(scope="module")
def design_dir(tmp_path_factory) -> Path:
    """A parent-layout ADNI-design run on the small table: its settings and its anchor p-values for replicates 0, 1."""

    directory = tmp_path_factory.mktemp("design")
    args = adc.build_parser().parse_args(["run", "--out-dir", str(directory), "--table", SIM_TABLE, "--reps", "2",
                                          "--perms", str(PERMUTATIONS)])
    settings = adc.build_settings(args)
    settings["cv"] = dict(SIM_CV)
    adc.ensure_settings(directory, settings)
    records = {}
    for replicate in (0, 1):
        record = adc.evaluate_unit(adc.ANCHOR, 0.0, replicate, settings)
        records[record["key"]] = record
    d10._write_csv(directory / "adni_design_pvalues.csv", adc.p_value_rows(records), adc.PVALUE_COLUMNS)
    return directory


def _plants(cohort, design, **changes) -> dict:
    kwargs = {"control_replicates": 2, "plant_replicates": 2, "permutations": PERMUTATIONS,
              "point_sizes": (STAND_IN_POINT,), "cv": SMALL_CV, "anchor_dir": gain.ANCHOR_DIR, "design_dir": design}
    kwargs.update(changes)
    return gain.plants_settings(cohort, **kwargs)


def _simulation(design, **changes) -> dict:
    kwargs = {"control_replicates": 2, "plant_replicates": 2, "permutations": PERMUTATIONS, "table": SIM_TABLE,
              "cv": SIM_CV, "anchor_dir": gain.ANCHOR_DIR, "design_dir": design}
    kwargs.update(changes)
    return gain.simulation_settings(**kwargs)


@pytest.fixture(scope="module")
def settings(cohort, design_dir) -> dict:
    return _plants(cohort, design_dir)


@pytest.fixture(scope="module")
def sim_settings(design_dir) -> dict:
    return _simulation(design_dir)


@pytest.fixture(scope="module")
def records(cohort, settings) -> dict:
    """Stand-in units of the plants profile, keyed by cell and point."""

    wanted = [(POINT, "anchor"), (gain.COHORT_POINT, "anchor"), (POINT, "sample_magnitude_e0.5"),
              (POINT, "heldout_orientation_e1"), (POINT, "sample_out_of_span_e0.5"),
              (gain.COHORT_POINT, "sample_shape_e1"), (POINT, "sample_translation_e1")]
    return {(point, cell): gain.evaluate_plants_unit(gain.Unit(point, cell, 0), settings, cohort)
            for point, cell in wanted}


@pytest.fixture(scope="module")
def prepared(tmp_path_factory, cohort, settings) -> Path:
    """An output root with a passed plants pre-flight."""

    out_dir = tmp_path_factory.mktemp("prepared")
    result = gain.plants_preflight(out_dir, settings, cohort)
    assert result["passed"], result
    return out_dir


@pytest.fixture(scope="module")
def sim_prepared(tmp_path_factory, sim_settings) -> Path:
    out_dir = tmp_path_factory.mktemp("sim_prepared")
    result = gain.simulation_preflight(out_dir, sim_settings)
    assert result["passed"], result
    return out_dir


def _counting_cv(monkeypatch) -> list:
    import motco.simulations.evaluation as evaluation_module

    calls: list = []
    original = evaluation_module.plsda_doubleCV

    def counting(*args, **kwargs):
        calls.append(1)
        return original(*args, **kwargs)

    monkeypatch.setattr(evaluation_module, "plsda_doubleCV", counting)
    return calls


def _scalers_and_built(cohort, settings, point=POINT, replicate=0):
    built = aac.build_unit(cohort, settings, point, replicate)
    return built, gain.sample_scalers(built["blocks"], cohort.contract)


def _paper_plants() -> dict:
    points = aac.design_points(ADNI_TABLE, (750,))
    return {"profile": "plants", "points": points, "cells": gain.plants_cells(points),
            "replicates": {"control": 500, "plant": 200}, "permutations": 999, "arms": list(gain.ARMS),
            "arm_specs": gain.arm_specs()}


def _paper_simulation() -> dict:
    point = gain.simulation_point("adni")
    return {"profile": "simulation", "points": [{"name": point, "participants": 1501, "table": ADNI_TABLE,
                                                 "subsampled": False}],
            "cells": gain.simulation_cells(point), "replicates": {"control": 500, "plant": 200},
            "permutations": 999, "arms": list(gain.ARMS), "arm_specs": gain.arm_specs()}


# ── 1.1 Resume, settings and parent refusals ────────────────────────────────


def test_settings_bind_the_parents_and_the_export(settings, sim_settings) -> None:
    assert set(settings["parents"]) == set(sim_settings["parents"]) == {
        "anchor_rates", "anchor_decision", "design_settings", "design_pvalues"}
    assert settings["parents"]["anchor_rates"]["path"] == "results/adni-real-noise-anchor-2026-11-17/rate_table.csv"
    assert len(settings["parents"]["anchor_decision"]["sha256"]) == 64
    assert sorted(settings["export"]["files"]) == ["covariates.csv", "manifest.json", "metadata.csv", "nmr.csv",
                                                   "p180.csv"]
    assert "export" not in sim_settings
    assert settings["base_seed"] == aac.BASE_SEED == 2039 and settings["own_base_seed"] == gain.BASE_SEED == 2040
    assert sim_settings["base_seed"] == 2040 and sim_settings["construction"]["base_seed"] == 2040


def test_resume_skips_complete_units(tmp_path, cohort, settings, prepared) -> None:
    out_dir = tmp_path / "out"
    shutil.copytree(prepared, out_dir)
    assert gain.run_shard(out_dir, settings, 1, 20, cohort=cohort, limit=1) == 1
    assert gain.run_shard(out_dir, settings, 1, 20, cohort=cohort) == 2
    assert gain.run_shard(out_dir, settings, 1, 20, cohort=cohort) == 0
    path = lsc.shard_path(out_dir / "plants", 1, 20)
    keys = [gain.record_key(json.loads(line)) for line in path.read_text().splitlines()]
    expected = [gain.unit_key(*unit) for position, unit in enumerate(gain.enumerate_units(settings))
                if position % 20 == 1]
    assert sorted(keys) == sorted(expected) and len(expected) == 3


def test_other_settings_are_refused_naming_the_keys(tmp_path, cohort, design_dir, prepared) -> None:
    out_dir = tmp_path / "out"
    shutil.copytree(prepared, out_dir)
    other = _plants(cohort, design_dir, permutations=PERMUTATIONS + 10)
    with pytest.raises(SystemExit, match="differing: permutations"):
        gain.run_shard(out_dir, other, 0, 1, cohort=cohort)
    with pytest.raises(SystemExit, match="differing: permutations"):
        gain.plants_preflight(out_dir, other, cohort)
    with pytest.raises(SystemExit, match=r"differing: replicates\.plant"):
        gain.run_shard(out_dir, _plants(cohort, design_dir, plant_replicates=3), 0, 1, cohort=cohort)


@pytest.mark.parametrize("role", sorted(gain.PARENT_FILES))
def test_a_changed_parent_is_refused_by_name(tmp_path, cohort, design_dir, export_dir, role) -> None:
    anchor, design = tmp_path / "anchor", tmp_path / "design"
    shutil.copytree(gain.ANCHOR_DIR, anchor)
    shutil.copytree(design_dir, design)
    settings = _plants(cohort, design_dir, anchor_dir=anchor, design_dir=design)
    sim = _simulation(design_dir, anchor_dir=anchor, design_dir=design)
    out_dir = tmp_path / "out"
    lsc.ensure_settings(out_dir / "plants", settings)
    lsc.ensure_settings(out_dir / "simulation", sim)
    source, name = gain.PARENT_FILES[role]
    with ((anchor if source == "anchor" else design) / name).open("a") as handle:
        handle.write("\n")
    calls = (lambda: gain.run_shard(out_dir, settings, 0, 1, cohort=cohort),
             lambda: gain.run_shard(out_dir, sim, 0, 1),
             lambda: gain.merge(out_dir, "plants", export_dir),
             lambda: gain.merge(out_dir, "simulation"),
             lambda: gain.decide(out_dir, export_dir))
    for call in calls:
        with pytest.raises(SystemExit, match=f"Parent {role} .* changed"):
            call()


def test_another_export_is_refused_by_file(tmp_path, export_dir, contract, design_dir, settings, prepared) -> None:
    out_dir = tmp_path / "out"
    shutil.copytree(prepared, out_dir)
    other = _copy_export(export_dir, tmp_path / "other")
    frame = _block(other, "nmr")
    frame.iloc[0, 0] += 0.5
    _rewrite(other, "nmr.csv", frame, index=True)
    moved = aac.load_cohort(other, contract)
    with pytest.raises(SystemExit, match=r"export\.files\.manifest\.json, export\.files\.nmr\.csv"):
        gain.run_shard(out_dir, _plants(moved, design_dir), 0, 1, cohort=moved)
    with pytest.raises(SystemExit, match="differ from the settings' for: manifest.json, nmr.csv"):
        gain.run_shard(out_dir, settings, 0, 1, cohort=moved)
    with pytest.raises(SystemExit, match="differ from the settings' for: manifest.json, nmr.csv"):
        gain.merge(out_dir, "plants", other)
    with pytest.raises(SystemExit, match="--export-dir is required"):
        gain.merge(out_dir, "plants")


def test_a_contract_violation_is_named_before_any_unit(monkeypatch, tmp_path, export_dir, contract) -> None:
    monkeypatch.setitem(aac.CONTRACTS, "adni", contract)
    broken = _copy_export(export_dir, tmp_path / "broken")
    (broken / "nmr.csv").unlink()
    out_dir = tmp_path / "out"
    for command in ("preflight", "run"):
        with pytest.raises(SystemExit, match="nmr.csv"):
            gain.main([command, "--profile", "plants", "--export-dir", str(broken), "--out-dir", str(out_dir)])
    with pytest.raises(SystemExit, match="--export-dir is required"):
        gain.main(["run", "--profile", "plants", "--out-dir", str(out_dir)])
    assert not out_dir.exists()


# ── 1.2 Points, references, cells, seeds ─────────────────────────────────────


def test_adni_point_tables_and_the_heldout_reference() -> None:
    points = aac.design_points(ADNI_TABLE, gain.POINT_SIZES)
    assert [point["name"] for point in points] == ["n750", "cohort"]
    assert points[0]["table"] == [[118, 155, 63], [111, 224, 79]]
    heldout = (np.asarray(ADNI_TABLE) - np.asarray(points[0]["table"])).tolist()
    assert heldout == [[119, 155, 62], [112, 224, 79]] and sum(map(sum, heldout)) == 751
    assert points[1]["table"] == ADNI_TABLE and not points[1]["subsampled"]


def test_the_reference_is_disjoint_from_the_sample(cohort, settings) -> None:
    built = aac.build_unit(cohort, settings, POINT, 0)
    rows = gain.heldout_rows(cohort, built["rows"])
    assert np.intersect1d(rows, built["rows"]).size == 0
    assert rows.size + built["rows"].size == len(cohort.metadata) and rows.size == len(cohort.metadata) - 88
    ids = set(cohort.metadata["rid"].iloc[rows]) | set(cohort.metadata["rid"].iloc[built["rows"]])
    assert ids == set(cohort.metadata["rid"])


def test_the_cohort_point_has_everyone_and_no_heldout_cell(cohort, settings) -> None:
    built = aac.build_unit(cohort, settings, gain.COHORT_POINT, 0)
    assert np.array_equal(built["rows"], np.arange(len(cohort.metadata)))
    assert gain.heldout_rows(cohort, built["rows"]).size == 0
    cells = gain.point_cells(settings, gain.COHORT_POINT)
    assert all(cell["reference"] != gain.HELDOUT for cell in cells)
    assert all(cell["reference"] in (None, gain.SAMPLE) for cell in cells)
    with pytest.raises(ValueError, match="held-out reference is empty"):
        gain.reference_configuration(cohort, built, gain.sample_scalers(built["blocks"], cohort.contract), "heldout")


def test_a_cohort_anchor_unit_is_the_anchor_drivers_unit(cohort, settings) -> None:
    anchor_settings = aac.build_settings(cohort, replicates=2, permutations=PERMUTATIONS, point_sizes=(40, 80),
                                         fixed_ranks=(2, 4), cv=SMALL_CV)
    for replicate in (0, 1):
        ours = aac.build_unit(cohort, settings, gain.COHORT_POINT, replicate)
        theirs = aac.build_unit(cohort, anchor_settings, aac.COHORT_POINT, replicate)
        assert ours["seeds"] == theirs["seeds"]
        assert ours["metadata"]["sex"].tolist() == theirs["metadata"]["sex"].tolist()
        assert ours["label_redraws"] == theirs["label_redraws"]


def test_cells_of_one_replicate_share_the_sample_labels_and_seed(records, cohort, settings) -> None:
    first = records[(POINT, "anchor")]
    for (point, _), record in records.items():
        if point == POINT:
            assert record["rrpp_seed"] == first["rrpp_seed"]
            assert record["seeds"]["subsample"] == first["seeds"]["subsample"]
            assert record["seeds"]["labels"] == first["seeds"]["labels"]
            assert record["cell_table"] == first["cell_table"]
    one, two = aac.build_unit(cohort, settings, POINT, 1), aac.build_unit(cohort, settings, POINT, 1)
    assert np.array_equal(one["rows"], two["rows"]) and one["metadata"].equals(two["metadata"])


def test_unit_counts_at_the_committed_replicates() -> None:
    plants = _paper_plants()
    units = gain.enumerate_units(plants)
    assert sum(unit.point == "n750" for unit in units) == 4200
    assert sum(unit.point == "cohort" for unit in units) == 3000
    assert len(units) == 7200 and len({gain.unit_key(*unit) for unit in units}) == 7200
    assert sum(unit.cell == "anchor" for unit in units) == 1000
    assert sum(unit.cell.startswith("heldout_") for unit in units) == 1200
    assert all(unit.point == "n750" for unit in units if unit.cell.startswith("heldout_"))
    assert len(gain.enumerate_units(_paper_simulation())) == 2600
    assert len(gain.simulation_cells("x")) == 10 and len(plants["cells"]) == 18


def test_cell_names_and_targets() -> None:
    cells = {cell["name"]: cell for cell in _paper_plants()["cells"]}
    assert cells["sample_translation_e1"]["kind"] == "control" and cells["sample_translation_e1"]["target"] is None
    assert cells["heldout_magnitude_e2"]["target"] == "delta" and cells["heldout_magnitude_e2"]["points"] == ["n750"]
    assert cells["sample_orientation_e0.25"]["target"] == "angle"
    assert cells["sample_shape_e0.5"]["target"] == "shape"
    assert cells["sample_out_of_span_e0.5"]["kind"] == "out_of_span"
    assert sorted(name for name in cells if name.startswith("heldout_")) == [
        "heldout_magnitude_e0.5", "heldout_magnitude_e1", "heldout_magnitude_e2", "heldout_orientation_e0.5",
        "heldout_orientation_e1", "heldout_shape_e1"]


def test_own_seeds_are_the_drivers_and_differ_from_the_anchors(settings) -> None:
    import hashlib

    direction = gain.own_seed(settings, POINT, 3, "direction")
    payload = json.dumps(["adni-gain-specificity", 2040, POINT, 3, "direction"], separators=(",", ":"))
    assert direction == int(hashlib.sha256(payload.encode()).hexdigest()[:8], 16) & 0x7FFFFFFF
    assert direction != gain.own_seed(settings, POINT, 3, "reference_split")
    assert direction not in aac.unit_seeds(settings, POINT, 3).values()
    with pytest.raises(ValueError, match="Unknown purpose"):
        gain.own_seed(settings, POINT, 3, "labels")


# ── 1.3 The covariate-adjusted configuration ────────────────────────────────


def _direct_configuration(features: np.ndarray, metadata: pd.DataFrame, stage_order) -> np.ndarray:
    stages = metadata["stage"].astype(str).to_numpy()
    indicators = np.column_stack([(stages == stage).astype(float) for stage in stage_order])
    columns = [metadata["age"].to_numpy(dtype=float), metadata["education"].to_numpy(dtype=float)]
    for name, levels in (("apoe4", ("1", "2")), ("phase", ("ADNI2GO",)), ("statin", ("1",))):
        values = metadata[name].astype(str).to_numpy()
        columns += [(values == level).astype(float) for level in levels if (values == level).any()]
    covariates = np.column_stack(columns)
    design = np.hstack([indicators, covariates - covariates.mean(axis=0)])
    return np.linalg.lstsq(design, features, rcond=None)[0][: len(stage_order)]


def test_the_configuration_equals_a_direct_least_squares_fit(cohort, settings) -> None:
    built, scalers = _scalers_and_built(cohort, settings)
    stage_order = cohort.contract.stage_order
    assert sorted(stage_order) != list(stage_order)  # CN, MCI, AD do not sort in declared order
    features = gain.standardized(built["blocks"], scalers, cohort.contract)
    ours = gain.reference_configuration(cohort, built, scalers, "sample")
    np.testing.assert_allclose(ours, _direct_configuration(features, built["metadata"], stage_order), atol=1e-10)
    rows = gain.heldout_rows(cohort, built["rows"])
    heldout = cohort.metadata.iloc[rows].reset_index(drop=True)
    held_features = gain.standardized({name: frame.iloc[rows] for name, frame in cohort.blocks.items()}, scalers,
                                      cohort.contract)
    np.testing.assert_allclose(gain.reference_configuration(cohort, built, scalers, "heldout"),
                               _direct_configuration(held_features, heldout, stage_order), atol=1e-10)


def test_the_features_are_the_analysis_own(cohort, settings) -> None:
    from motco import analysis

    built, scalers = _scalers_and_built(cohort, settings)
    dataset = aac.make_dataset(built["blocks"], built["metadata"], cohort.contract)
    inputs = analysis._trajectory_inputs(dataset, aac.analysis_params(settings, 1))
    np.testing.assert_array_equal(gain.standardized(built["blocks"], scalers, cohort.contract),
                                  inputs.features.to_numpy(dtype=float))


def test_the_centred_configuration_does_not_depend_on_the_reference_row(cohort, settings) -> None:
    built, scalers = _scalers_and_built(cohort, settings)
    contract = cohort.contract
    features = gain.standardized(built["blocks"], scalers, contract)
    stages = gain.stage_positions(built["metadata"], contract)
    covariates = gain.covariate_columns(built["metadata"], contract)
    mean = gain.adjusted_configuration(features, stages, covariates, 3)
    other = gain.adjusted_configuration(features, stages, covariates, 3, reference_row=covariates[0] + 2.0)
    assert not np.allclose(mean, other)
    np.testing.assert_allclose(mean - mean.mean(axis=0), other - other.mean(axis=0), atol=1e-10)


def test_a_heldout_configuration_ignores_the_order_of_the_analysed_rows(cohort, settings) -> None:
    built, scalers = _scalers_and_built(cohort, settings)
    order = np.random.default_rng(3).permutation(len(built["rows"]))
    shuffled = {**built, "rows": built["rows"][order],
                "blocks": {name: frame.iloc[order] for name, frame in built["blocks"].items()},
                "metadata": built["metadata"].iloc[order].reset_index(drop=True)}
    again = gain.sample_scalers(shuffled["blocks"], cohort.contract)
    np.testing.assert_allclose(gain.reference_configuration(cohort, built, scalers, "heldout"),
                               gain.reference_configuration(cohort, shuffled, again, "heldout"), atol=1e-12)


def test_a_rank_deficient_fit_fails_the_unit(monkeypatch, cohort, settings) -> None:
    features = np.random.default_rng(0).normal(size=(20, 4))
    stages = np.repeat([0, 1, 2], [7, 7, 6])
    duplicated = np.column_stack([np.arange(20.0), np.arange(20.0)])
    with pytest.raises(ValueError, match="rank-deficient"):
        gain.adjusted_configuration(features, stages, duplicated, 3)
    with pytest.raises(ValueError, match="rank-deficient"):
        gain.adjusted_configuration(features, np.repeat([0, 1], 10), np.zeros((20, 0)), 3)
    monkeypatch.setattr(gain, "covariate_columns",
                        lambda metadata, contract: np.column_stack([metadata["age"], metadata["age"]]))
    record = gain.evaluate_plants_unit(gain.Unit(POINT, "sample_magnitude_e1", 0), settings, cohort)
    assert record["status"] == "failed" and "rank-deficient" in record["error_message"]


# ── 1.4 The plant ────────────────────────────────────────────────────────────


def _plant(cohort, settings, cell_name: str, point=POINT, replicate=0):
    built = aac.build_unit(cohort, settings, point, replicate)
    return built, gain.build_plant(cohort, settings, built, gain.cell_spec(settings, cell_name))


def test_only_the_second_group_moves_by_one_vector_per_stage(cohort, settings) -> None:
    built, planted = _plant(cohort, settings, "sample_orientation_e0.5")
    contract = cohort.contract
    labels = built["metadata"]["sex"].to_numpy()
    stages = gain.stage_positions(built["metadata"], contract)
    start = 0
    for name in contract.block_names:
        before, after = built["blocks"][name], planted["blocks"][name]
        assert after.index.equals(before.index) and after.columns.equals(before.columns)
        moved = after.to_numpy() - before.to_numpy()
        width = before.shape[1]
        assert np.all(moved[labels == "Female"] == 0.0)
        for stage in range(3):
            rows = moved[(labels == "Male") & (stages == stage)]
            np.testing.assert_allclose(rows, np.tile(planted["shift_raw"][stage, start:start + width],
                                                     (rows.shape[0], 1)), atol=1e-12)
        start += width
    assert built["metadata"].equals(aac.build_unit(cohort, settings, POINT, 0)["metadata"])
    assert planted["shift_raw"] is not None and np.any(planted["shift_raw"])


def test_each_edit_is_the_declared_one(cohort, settings) -> None:
    built, scalers = _scalers_and_built(cohort, settings)
    z = gain.reference_configuration(cohort, built, scalers, "sample")
    centred = z - z.mean(axis=0)
    shift, _ = gain.plant_shift(z, "magnitude", 0.5)
    np.testing.assert_allclose(shift, 0.5 * centred, atol=1e-12)
    shift, _ = gain.plant_shift(z, "orientation", 0.5)
    edited = z + shift
    gram = lambda m: (m - m.mean(axis=0)) @ (m - m.mean(axis=0)).T  # noqa: E731
    np.testing.assert_allclose(gram(edited), gram(z), atol=1e-10)  # rigid
    net, turned = z[-1] - z[0], edited[-1] - edited[0]
    cosine = net @ turned / (np.linalg.norm(net) * np.linalg.norm(turned))
    assert math.isclose(math.degrees(math.acos(min(1.0, cosine))), 45.0, abs_tol=1e-8)
    path = lambda m: float(np.linalg.norm(np.diff(m, axis=0), axis=1).sum())  # noqa: E731
    shift, details = gain.plant_shift(z, "shape", 1.0)
    assert math.isclose(path(z + shift), path(z), rel_tol=1e-12) and details["shape_stage"] == 1
    shift, _ = gain.plant_shift(z, "translation", 1.0)
    np.testing.assert_allclose(shift - shift[0], 0.0, atol=1e-12)
    assert math.isclose(float(np.linalg.norm(shift[0])), path(z) / 2, rel_tol=1e-12)


def test_a_zero_effect_returns_the_frames_untouched(cohort, settings) -> None:
    built = aac.build_unit(cohort, settings, POINT, 0)
    zero = np.zeros((3, sum(STAND_IN_FEATURES.values())))
    out = gain.apply_plant(built["blocks"], built["metadata"], cohort.contract, zero)
    assert all(out[name] is built["blocks"][name] for name in cohort.contract.block_names)
    built2, scalers = _scalers_and_built(cohort, settings)
    z = gain.reference_configuration(cohort, built2, scalers, "sample")
    for mode in ("magnitude", "orientation", "shape", "translation"):
        shift, _ = gain.plant_shift(z, mode, 0.0)
        assert not np.any(np.abs(shift) > 1e-12), mode


def test_the_planted_groups_adjusted_means_move_by_the_shift(cohort, settings) -> None:
    contract = cohort.contract
    for cell in ("sample_shape_e1", "heldout_magnitude_e1", "sample_out_of_span_e0.25"):
        built, planted = _plant(cohort, settings, cell)
        scalers = gain.sample_scalers(built["blocks"], contract)
        male = (built["metadata"]["sex"] == "Male").to_numpy()
        metadata = built["metadata"].loc[male].reset_index(drop=True)
        stages = gain.stage_positions(metadata, contract)
        covariates = gain.covariate_columns(metadata, contract)
        before = gain.adjusted_configuration(gain.standardized(built["blocks"], scalers, contract)[male], stages,
                                             covariates, 3)
        after = gain.adjusted_configuration(gain.standardized(planted["blocks"], scalers, contract)[male], stages,
                                            covariates, 3)
        np.testing.assert_allclose(after - before, planted["shift_std"], atol=1e-10)


def test_the_out_of_span_direction(cohort, settings) -> None:
    built, scalers = _scalers_and_built(cohort, settings)
    z = gain.reference_configuration(cohort, built, scalers, "sample")
    seed = gain.own_seed(settings, POINT, 0, "direction")
    vector = gain.out_of_span_vector(z, seed)
    assert math.isclose(float(np.linalg.norm(vector)), 1.0, rel_tol=1e-12)
    assert float(np.max(np.abs((z - z.mean(axis=0)) @ vector))) <= 1e-10
    assert np.array_equal(vector, gain.out_of_span_vector(z, seed))
    _, quarter = _plant(cohort, settings, "sample_out_of_span_e0.25")
    _, half = _plant(cohort, settings, "sample_out_of_span_e0.5")
    assert quarter["direction_seed"] == half["direction_seed"] == seed
    for planted in (quarter, half):
        assert planted["edit"]["span_alignment"] == 90.0
        np.testing.assert_allclose(planted["configurations"]["sample"], z, atol=1e-12)


# ── 2.1 Arms on one fit ──────────────────────────────────────────────────────


@pytest.mark.parametrize("point", [POINT, "cohort"])
def test_a_planted_unit_makes_one_cross_validation(monkeypatch, cohort, settings, point) -> None:
    calls = _counting_cv(monkeypatch)
    record = gain.evaluate_plants_unit(gain.Unit(point, "sample_magnitude_e1", 1), settings, cohort)
    assert record["status"] == "completed", record.get("inconsistent_checks") or record.get("error_message")
    assert len(calls) == 1 and list(record["arms"]) == list(gain.ARMS)
    rank = record["production"]["rank"]
    assert record["production"]["fit"] == "cross_validated"
    for arm, entry in record["arms"].items():
        assert entry["status"] == "completed" and entry["rank"] == rank, arm
    assert all(entry["fit"] == "forced_at_cv_rank" for arm, entry in record["arms"].items() if arm != "pls_cv")


def test_a_cohort_anchor_unit_makes_no_cross_validation(monkeypatch, cohort, settings) -> None:
    calls = _counting_cv(monkeypatch)
    record = gain.evaluate_plants_unit(gain.Unit("cohort", "anchor", 1), settings, cohort)
    assert record["status"] == "completed" and calls == []
    rank = int(settings["cohort_rank"])
    assert record["production"] == {**record["production"], "rank": rank, "fit": "forced_at_cohort_rank"}
    assert all(entry["rank"] == rank and entry["fit"] == "forced_at_cohort_rank" for entry in record["arms"].values())


def test_arm_parameters_differ_only_in_the_named_settings(settings) -> None:
    production = dataclasses.asdict(gain.arm_params("pls_cv", settings, 11, rank=None))
    assert production["forced_components"] is None and production["residual_correction"] == "df"
    expected = {
        "predictive_orthonormal": {"latent_coordinates": "predictive_orthonormal"},
        "pls_cv__within_stage": {"permutation_strata": "stage"},
        "predictive_orthonormal__within_stage": {"latent_coordinates": "predictive_orthonormal",
                                                 "permutation_strata": "stage"},
    }
    for arm, changes in expected.items():
        params = dataclasses.asdict(gain.arm_params(arm, settings, 11, rank=3))
        differing = {key: value for key, value in params.items() if production[key] != value}
        assert differing == {"forced_components": 3, **changes}, arm
    assert production["seed"] == 11 and production["permutations"] == PERMUTATIONS
    assert production["min_components"] == 2 and production["orientation_estimator"] == "net_displacement"


def test_a_failing_arm_leaves_the_others(monkeypatch, cohort, settings) -> None:
    original = gain._analyse

    def failing(dataset, params):
        if params.latent_coordinates == "predictive_orthonormal" and params.permutation_strata == "stage":
            raise RuntimeError("planted failure")
        return original(dataset, params)

    monkeypatch.setattr(gain, "_analyse", failing)
    record = gain.evaluate_plants_unit(gain.Unit(POINT, "sample_magnitude_e1", 0), settings, cohort)
    statuses = {arm: entry["status"] for arm, entry in record["arms"].items()}
    assert statuses == {"pls_cv": "completed", "predictive_orthonormal": "completed",
                        "pls_cv__within_stage": "completed", "predictive_orthonormal__within_stage": "failed"}
    assert record["arms"]["predictive_orthonormal__within_stage"]["error_message"] == "planted failure"


def test_a_failing_unit_is_recorded_and_the_shard_goes_on(monkeypatch, tmp_path, cohort, settings, prepared) -> None:
    out_dir = tmp_path / "out"
    shutil.copytree(prepared, out_dir)
    original = aac.build_unit

    def failing(cohort_, settings_, point, replicate):
        if point == POINT and replicate == 0:
            raise RuntimeError("planted unit failure")
        return original(cohort_, settings_, point, replicate)

    monkeypatch.setattr(aac, "build_unit", failing)
    assert gain.run_shard(out_dir, settings, 0, 10, cohort=cohort, limit=4) == 4
    statuses = {gain.record_key(record): record["status"]
                for record in gain.read_records([lsc.shard_path(out_dir / "plants", 0, 10)])}
    assert statuses[f"{POINT}/anchor/0"] == "failed"
    assert any(status == "completed" for key, status in statuses.items() if not key.startswith(f"{POINT}/"))


# ── 2.2 Checks ───────────────────────────────────────────────────────────────


def test_every_check_passes_on_stand_in_units(records) -> None:
    for key, record in records.items():
        assert record["status"] == "completed", (key, record.get("error_message"), record.get("inconsistent_checks"))
        for name in gain.CHECKS["plants"]:
            check = record["consistency"][name]
            assert check["passed"], (key, name, check)
    assert records[(POINT, "anchor")]["consistency"]["plant"]["applicable"] is False
    assert records[(POINT, "sample_magnitude_e0.5")]["consistency"]["plant"]["applicable"] is True


def test_a_perturbed_plant_model_or_table_is_flagged(cohort, settings) -> None:
    built, planted = _plant(cohort, settings, "sample_magnitude_e1")
    contract = cohort.contract
    args = (built["blocks"], planted["blocks"], built["metadata"], contract)
    assert gain.check_plant(*args, planted["shift_raw"], planted["scale"])["passed"]
    wrong = planted["shift_raw"].copy()
    wrong[1, 0] += 1e-6 * planted["scale"][0]
    flagged = gain.check_plant(*args, wrong, planted["scale"])
    assert not flagged["passed"] and set(flagged["failing"]) == {"second_group"}
    leaked = {name: frame.copy() for name, frame in planted["blocks"].items()}
    female = np.flatnonzero((built["metadata"]["sex"] == "Female").to_numpy())[0]
    leaked["p180"].iloc[female, 0] += 1.0
    flagged = gain.check_plant(built["blocks"], leaked, built["metadata"], contract, planted["shift_raw"],
                               planted["scale"])
    assert not flagged["passed"] and "first_group" in flagged["failing"]
    dataset = aac.make_dataset(planted["blocks"], built["metadata"], contract)
    result = aac._analyse(dataset, aac.analysis_params(settings, 5, forced=3))
    model = result.latent_artifacts.model
    perturbed = copy.deepcopy(model)
    perturbed.x_rotations_ = perturbed.x_rotations_ + 1e-6
    flagged = aac.check_same_fit({"a": model, "b": perturbed}, ["a", "b"])
    assert not flagged["passed"] and set(flagged["failing_arms"]) == {"b"}
    assert not aac.check_cell_table([[1, 2], [3, 4]], {"labels": [[1, 2], [4, 3]]})["passed"]


def _synthetic(settings, point, cell, replicate, *, p=0.5, rank=2, overrides=None, status="completed") -> dict:
    overrides = overrides or {}
    arms = {arm: {"status": "completed", "rank": rank, "measured_dimension": 2,
                  "pair_statistics": {"delta": 1.0, "angle": 10.0, "shape": 0.1},
                  "p_values": {statistic: overrides.get(arm, {}).get(statistic, p) for statistic in gain.STATISTICS},
                  "angle_null_q95": 30.0, "seconds": 0.1} for arm in gain.ARMS}
    spec = gain.cell_spec(settings, cell)
    return {
        "profile": settings["profile"], "point": point, "cell": cell, "replicate": replicate, "status": status,
        "inconsistent_checks": [] if status == "completed" else ["plant"],
        "consistency": {"plant": {"deviation": 1.0, "failing": {"second_group": 1.0}}},
        **{key: spec[key] for key in ("kind", "reference", "mode", "effect", "target")},
        "settings_hash": lsc.settings_hash(settings), "subsample_redraws": 0, "label_redraws": 1,
        "production": {"rank": rank, "cv_mean_auroc": 0.6}, "arms": arms, "runtime_seconds": {"unit": 2.0},
        "plant": None if cell == "anchor" else {"norm_per_stage": [1.0, 0.5, 1.5], "max_abs_entry": 0.2,
                                                "scale_drift": 0.01 * (replicate + 1), "edit": {}},
        "capture": {"plant": {"span": 0.9, "predictive": 0.4}, "reference": {"span": 0.8, "predictive": 0.3}},
        "reliability": {"seed": 1, "correlation": 0.2 + 0.1 * replicate,
                        "full_reference": 2 * (0.2 + 0.1 * replicate) / (1.2 + 0.1 * replicate),
                        "eigenvalues": [2.0, 0.1], "eigenvalue_ratio": 0.05}
        if cell == "anchor" and point == POINT else None,
    }


def test_an_inconsistent_unit_leaves_every_rate(settings) -> None:
    records = [_synthetic(settings, POINT, "sample_magnitude_e1", replicate, p=0.01) for replicate in range(3)]
    records[1]["status"] = "inconsistent"
    rows = {(row["arm"], row["statistic"]): row for row in gain.rate_table(records, settings)
            if row["point"] == POINT and row["cell"] == "sample_magnitude_e1"}
    assert rows[("pls_cv", "delta")]["count"] == 2 and rows[("pls_cv", "delta")]["inconsistent"] == 1
    paired = [row for row in gain.paired_table(records, settings) if row["cell"] == "sample_magnitude_e1"]
    assert all(row["n_paired"] == 2 for row in paired if row["point"] == POINT)


# ── 2.3 Readings ─────────────────────────────────────────────────────────────


def test_a_capture_share_is_a_direct_projection(cohort, settings, records) -> None:
    built, planted = _plant(cohort, settings, "sample_magnitude_e1")
    dataset = aac.make_dataset(planted["blocks"], built["metadata"], cohort.contract)
    model = aac._analyse(dataset, aac.analysis_params(settings, 5, forced=3)).latent_artifacts.model
    vectors = np.random.default_rng(1).normal(size=(3, model.x_rotations_.shape[0]))
    basis, _ = np.linalg.qr(model.x_rotations_)
    direct = float(np.sum((vectors @ basis) ** 2) / np.sum(vectors**2))
    assert math.isclose(gain.capture_shares(vectors, model)["span"], direct, rel_tol=1e-12)
    inside = np.vstack([basis[:, 0], basis[:, 1] - basis[:, 2]])
    assert math.isclose(gain.capture_shares(inside, model)["span"], 1.0, rel_tol=1e-12)
    outside = vectors - (vectors @ basis) @ basis.T
    assert gain.capture_shares(outside, model)["span"] <= 1e-12
    assert gain.capture_shares(np.zeros((3, basis.shape[0])), model) == {"span": None, "predictive": None}
    from motco.simulations.diagnostics import LinearLatentMap

    predictive = LinearLatentMap.from_pls_predictive(model, "predictive_orthonormal").basis
    direct = float(np.sum((vectors @ predictive) ** 2) / np.sum(vectors**2))
    assert math.isclose(gain.capture_shares(vectors, model)["predictive"], direct, rel_tol=1e-12)
    capture = records[(POINT, "sample_magnitude_e0.5")]["capture"]
    assert set(capture) == {"plant", "reference"}
    assert math.isclose(capture["plant"]["span"], capture["reference"]["span"], rel_tol=1e-9)  # a stretch
    assert set(records[(POINT, "anchor")]["capture"]) == {"sample", "heldout"}
    assert set(records[("cohort", "anchor")]["capture"]) == {"sample"}


def test_the_reliability_is_reproducible_from_its_seed(cohort, settings, records) -> None:
    record = records[(POINT, "anchor")]
    reading = record["reliability"]
    assert reading["seed"] == record["seeds"]["reference_split"] == gain.own_seed(settings, POINT, 0,
                                                                                    "reference_split")
    assert set(reading) == {"seed", "correlation", "full_reference", "eigenvalues", "eigenvalue_ratio"}
    r = reading["correlation"]
    assert reading["full_reference"] == pytest.approx(2 * r / (1 + r), rel=1e-12)
    built = aac.build_unit(cohort, settings, POINT, 0)
    again = gain.heldout_reliability(cohort, built, reading["seed"])
    assert again["seed"] == reading["seed"]
    assert again["correlation"] == pytest.approx(reading["correlation"], rel=1e-12)
    assert again["eigenvalues"] == pytest.approx(reading["eigenvalues"], rel=1e-12)
    assert -1.0 <= reading["correlation"] <= 1.0 and reading["eigenvalues"][0] >= reading["eigenvalues"][1]
    first, second = gain.split_within_stage(np.repeat([0, 1, 2], [5, 6, 7]), 4, 3)
    assert np.intersect1d(first, second).size == 0 and first.size + second.size == 18 and first.size == 8
    assert all(record.get("reliability") is None for key, record in records.items() if key != (POINT, "anchor"))


def test_the_reliability_eigenvalues_are_the_stage_space_cross_products(cohort) -> None:
    metadata = cohort.metadata
    stages = gain.stage_positions(metadata, cohort.contract)
    truth = np.random.default_rng(2).normal(size=(3, 7)) * np.array([[3.0], [0.5], [-2.0]])
    exact = gain.reliability(truth[stages], stages, metadata, cohort.contract, 11)
    centred = truth - truth.mean(axis=0)
    expected = np.sort(np.linalg.eigvalsh(centred @ centred.T))[::-1][:2]
    assert exact["correlation"] == pytest.approx(1.0, rel=1e-12) and exact["full_reference"] == pytest.approx(1.0)
    np.testing.assert_allclose(exact["eigenvalues"], expected, rtol=1e-10)
    noisy = truth[stages] + np.random.default_rng(3).normal(scale=5.0, size=(stages.size, 7))
    reading = gain.reliability(noisy, stages, metadata, cohort.contract, 11)
    assert -1.0 < reading["correlation"] < 1.0 and reading["eigenvalues"][0] >= reading["eigenvalues"][1]


# ── 2.4 Records, docstring ───────────────────────────────────────────────────


def test_a_record_holds_no_participant_data(records, cohort) -> None:
    ids = cohort.metadata["rid"].tolist()
    features = sum(STAND_IN_FEATURES.values())
    for key, record in records.items():
        text = json.dumps(record)
        assert not any(identifier in text for identifier in ids), key
        n = record["participants"]
        assert _walk(record, n) == [] and _walk(record, len(ids) - n if n < len(ids) else n) == []
        assert _walk(record, features) == [], key
        assert not {"rows", "labels", "shift", "shift_raw", "shift_std", "configuration"} & set(record)


def test_the_docstrings_commands_run_against_the_stand_in(monkeypatch, tmp_path, export_dir, contract,
                                                          design_dir) -> None:
    monkeypatch.setitem(aac.CONTRACTS, "adni", contract)
    monkeypatch.setattr(gain, "POINT_SIZES", (STAND_IN_POINT,))
    monkeypatch.setattr(gain, "CV", SMALL_CV)
    monkeypatch.setattr(gain, "PERMUTATIONS", PERMUTATIONS)
    monkeypatch.setattr(gain, "DESIGN_DIR", design_dir)
    monkeypatch.setattr(gain, "SIMULATION_TABLE", SIM_TABLE)
    monkeypatch.setattr(gain, "SIMULATION_CV", SIM_CV)
    monkeypatch.setattr(gain, "SAMPLE_EFFECTS", {"magnitude": (1.0,)})
    monkeypatch.setattr(gain, "HELDOUT_EFFECTS", {"magnitude": (1.0,)})
    monkeypatch.setattr(gain, "OUT_OF_SPAN_EFFECTS", (0.5,))
    monkeypatch.setattr(gain, "SIMULATION_EFFECTS", {"magnitude": (1.0,)})
    commands = [line.strip() for line in gain.__doc__.splitlines()
                if line.strip().startswith("python scripts/adni_gain_specificity_check.py")]
    assert [shlex.split(command)[2] for command in commands] == ["preflight", "run", "merge", "preflight", "run",
                                                                 "merge", "decide"]
    out_dir = tmp_path / "smoke"
    before = _snapshot(export_dir)
    for command in commands:
        argv = shlex.split(command.replace("$EXPORT", str(export_dir)).replace("$OUT", str(out_dir)))[2:]
        assert gain.main(argv) == 0, command
    decision = json.loads((out_dir / gain.DECISION_NAME).read_text())
    assert decision["complete"]["plants"]["complete_units"] == 4 * (5 + 4)
    assert decision["complete"]["simulation"]["complete_units"] == 4 * 3
    assert decision["complete"]["plants"]["inconsistent_units"] == 0
    assert decision["complete"]["simulation"]["inconsistent_units"] == 0
    assert decision["verdict"] == "gate_failed" and decision["indicative"]["indicative"]
    assert decision["gate"]["failures"] == [
        f"the cohort anchor does not reproduce the anchor run's counts ({len(gain.ARMS) * 3} differing)"]
    assert decision["preflight"]["plants"]["passed"] and decision["preflight"]["simulation"]["passed"]
    assert _snapshot(export_dir) == before
    for path in [*(out_dir / "plants").glob("*.csv"), *(out_dir / "plants").glob("*.json"), out_dir /
                 gain.DECISION_NAME]:
        text = path.read_text()
        assert not any(identifier in text for identifier in aac.load_cohort(export_dir, contract).metadata["rid"])


# ── 3.1–3.2 The simulation profile ───────────────────────────────────────────


def test_simulation_settings_carry_the_construction() -> None:
    settings = gain.simulation_settings()
    construction = settings["construction"]
    assert construction["check"] == "adni_design" and construction["base_seed"] == 2040
    assert construction["generator"]["group_stage_sizes"] == ADNI_TABLE
    assert {construction["generator"][key] for key in ("delta_methyl", "delta_expr", "delta_protein")} == {0.02}
    assert {key: construction["injection"][key] for key in ("s_z", "s_p")} == {"s_z": 0.1, "s_p": 1.75}
    assert {key: construction["generator"][key] for key in ("construction", "span_layers", "span_scale")} == {
        "construction": "stage_span", "span_layers": ["methylation", "expression"], "span_scale": "planted"}
    assert settings["analysis"]["residual_correction"] == "df" and settings["analysis"]["cv"] == aac.CV
    assert settings["points"][0]["name"] == "table_adni" and settings["points"][0]["participants"] == 1501
    assert settings["replicates"] == {"control": 500, "plant": 200} and settings["permutations"] == 999


def test_the_anchor_is_the_parents_anchor_at_a_seed(sim_settings) -> None:
    args = adc.build_parser().parse_args(["run", "--out-dir", "unused", "--table", SIM_TABLE, "--reps", "2",
                                          "--perms", str(PERMUTATIONS), "--base-seed", "2040"])
    parent = adc.build_settings(args)
    ours = adc.build_unit(adc.ANCHOR, 0.0, 1, sim_settings["construction"])
    theirs = adc.build_unit(adc.ANCHOR, 0.0, 1, parent)
    for name in ("methylation", "expression"):
        pd.testing.assert_frame_equal(ours["blocks"][name], theirs["blocks"][name], check_exact=True)
    pd.testing.assert_frame_equal(ours["metadata"], theirs["metadata"], check_exact=True)
    assert ours["seeds"] == theirs["seeds"]


def test_cells_of_a_replicate_share_everything_but_the_plant(sim_settings) -> None:
    construction = sim_settings["construction"]
    anchor = adc.build_unit(adc.ANCHOR, 0.0, 1, construction)
    rotation = adc.build_unit("orientation", 1.0, 1, construction)
    stretch = adc.build_unit("magnitude", 1.0, 1, construction)
    for built in (rotation, stretch):
        assert built["seeds"] == anchor["seeds"]
        pd.testing.assert_frame_equal(built["metadata"], anchor["metadata"], check_exact=True)
    geometry = rotation["dataset"].truth["stage_span"]["intended_geometry"]
    assert math.isclose(geometry["angle"], 90.0, abs_tol=1e-8)
    assert abs(geometry["delta"]) < 1e-10 and abs(geometry["shape"]) < 1e-10
    geometry = stretch["dataset"].truth["stage_span"]["intended_geometry"]
    # A zero angle carries arccos round-off of order 1e-6 degrees.
    assert geometry["delta"] > 0 and abs(geometry["angle"]) < 1e-4 and abs(geometry["shape"]) < 1e-10
    assert rotation["dataset"].truth["span_scale"] == "planted"


def test_a_simulated_unit_makes_one_cross_validation(monkeypatch, sim_settings) -> None:
    calls = _counting_cv(monkeypatch)
    record = gain.evaluate_simulation_unit(gain.Unit(SIM_POINT, "population_orientation_e1", 0), sim_settings)
    assert record["status"] == "completed", record.get("inconsistent_checks") or record.get("error_message")
    assert len(calls) == 1
    rank = record["production"]["rank"]
    assert all(entry["rank"] == rank and entry["status"] == "completed" for entry in record["arms"].values())
    assert record["consistency"]["same_fit"]["passed"] and record["consistency"]["rank"]["passed"]
    assert record["consistency"]["construction"]["passed"]
    assert math.isclose(record["truth"]["intended_geometry"]["angle"], 90.0, abs_tol=1e-8)
    assert record["cell_table"] == [[12, 12, 12], [12, 12, 12]]


def test_a_dataset_of_another_construction_is_flagged(sim_settings) -> None:
    cell = gain.cell_spec(sim_settings, "population_shape_e1")
    truth = {"construction": "stage_span", "span_scale": "planted", "trajectory_mode": "shape",
             "group_effect_size": 1.0}
    assert gain.check_construction(truth, cell)["passed"]
    for change in ({"construction": "feature_surgery"}, {"span_scale": None}, {"trajectory_mode": "orientation"},
                   {"group_effect_size": 2.0}):
        flagged = gain.check_construction({**truth, **change}, cell)
        assert not flagged["passed"] and flagged["differing"] == list(change)
    anchor = gain.cell_spec(sim_settings, "anchor")
    assert gain.check_construction({**truth, "trajectory_mode": "magnitude", "group_effect_size": 0.0},
                                   anchor)["passed"]


# ── 3.3 sbatch ───────────────────────────────────────────────────────────────


def test_the_sbatch_template_dry_runs_both_profiles(tmp_path) -> None:
    script = REPO / "scripts" / "adni_gain_specificity_array.sbatch"
    shell = shutil.which("bash") or "bash"
    subprocess.run([shell, "-n", str(script)], check=True)
    env = {"PATH": "/usr/bin:/bin", "GAIN_DRY_RUN": "1", "GAIN_PROFILE": "simulation",
           "GAIN_OUT": str(tmp_path / "run"), "GAIN_REPO": str(REPO), "N_SHARDS": "100", "SLURM_ARRAY_TASK_ID": "3",
           "GAIN_ARGS": "--limit 1"}
    driver = str(REPO / "scripts" / "adni_gain_specificity_check.py")

    def dry(**changes: str) -> subprocess.CompletedProcess[str]:
        return subprocess.run([shell, str(script)], env={**env, **changes}, cwd=tmp_path, capture_output=True,
                              text=True)

    words = dry().stdout.split()
    assert words[0] == f"PYTHONPATH={REPO / 'src'}"
    assert words[1:] == ["python", driver, "run", "--profile", "simulation", "--out-dir", str(tmp_path / "run"),
                         "--shard", "3/100", "--limit", "1"]
    export = str(tmp_path / "export")
    plants = dry(GAIN_PROFILE="plants", GAIN_EXPORT=export, N_SHARDS="1", SLURM_ARRAY_TASK_ID="0",
                 GAIN_ARGS="--workers 48")
    assert plants.returncode == 0
    assert plants.stdout.split()[1:] == ["python", driver, "run", "--profile", "plants", "--export-dir", export,
                                         "--out-dir", str(tmp_path / "run"), "--shard", "0/1", "--workers", "48"]
    assert not (tmp_path / "run").exists()
    missing = dry(GAIN_PROFILE="plants")
    assert missing.returncode == 2 and "GAIN_EXPORT" in missing.stderr
    stray = dry(GAIN_EXPORT=export)
    assert stray.returncode == 2 and "plants profile" in stray.stderr
    for changes in ({"GAIN_OUT": "run"}, {"GAIN_PROFILE": "plants", "GAIN_EXPORT": "export"}):
        relative = dry(**changes)
        assert relative.returncode == 2 and "absolute" in relative.stderr
    unknown = dry(GAIN_PROFILE="cohort")
    assert unknown.returncode == 2 and "simulation or plants" in unknown.stderr


# ── 4.1–4.2 Pre-flights ──────────────────────────────────────────────────────


def test_the_plants_preflight_passes_on_the_stand_in(prepared, settings) -> None:
    record = json.loads((prepared / "plants" / gain.PREFLIGHT_NAME).read_text())
    assert record["passed"] and record["settings_hash"] == lsc.settings_hash(settings)
    assert [unit["unit"] for unit in record["units"]] == ["cohort/0", "cohort/1"]
    for unit in record["units"]:
        assert unit["fidelity"]["passed"] and unit["forced_rank"] == settings["cohort_rank"]


@pytest.mark.parametrize("change", ["rank", "auroc"])
def test_a_plants_preflight_on_another_fit_fails_and_the_run_refuses(tmp_path, export_dir, contract, design_dir,
                                                                     change) -> None:
    wrong = dataclasses.replace(contract, **(
        {"fidelity_rank": contract.fidelity_rank + 1} if change == "rank"
        else {"fidelity_auroc": round(contract.fidelity_auroc + 0.01, 6)}))
    cohort = aac.load_cohort(export_dir, wrong)
    settings = _plants(cohort, design_dir)
    result = gain.plants_preflight(tmp_path, settings, cohort)
    assert not result["passed"] and not all(unit["fidelity"]["passed"] for unit in result["units"])
    with pytest.raises(SystemExit, match="pre-flight failed"):
        gain.run_shard(tmp_path, settings, 0, 1, cohort=cohort)


def test_a_forced_call_that_differs_fails_the_plants_preflight(monkeypatch, tmp_path, cohort, settings) -> None:
    original = aac._analyse

    def drifting(dataset, params):
        result = original(dataset, params)
        if params.forced_components is not None:
            result.record["p_values"]["angle"] = 0.5
        return result

    monkeypatch.setattr(aac, "_analyse", drifting)
    result = gain.plants_preflight(tmp_path, settings, cohort)
    assert not result["passed"] and result["units"][0]["p_values_differing"] == ["angle"]


def test_the_run_waits_on_each_preflight(tmp_path, cohort, settings, sim_settings) -> None:
    with pytest.raises(SystemExit, match="run `preflight --profile plants` first"):
        gain.run_shard(tmp_path, settings, 0, 1, cohort=cohort)
    with pytest.raises(SystemExit, match="run `preflight --profile simulation` first"):
        gain.run_shard(tmp_path, sim_settings, 0, 1)


def test_the_simulation_preflight_passes(sim_prepared, design_dir) -> None:
    record = json.loads((sim_prepared / "simulation" / gain.PREFLIGHT_NAME).read_text())
    assert record["passed"] and [unit["replicate"] for unit in record["units"]] == [0, 1]
    committed = lsc._read_csv(design_dir / "adni_design_pvalues.csv")
    for unit in record["units"]:
        assert unit["differing"] == []
        assert unit["committed"] == {row["statistic"]: row["p_value"] for row in committed
                                     if int(row["replicate"]) == unit["replicate"]}


def test_a_differing_committed_pvalue_fails_the_simulation_preflight(tmp_path, design_dir) -> None:
    design = tmp_path / "design"
    shutil.copytree(design_dir, design)
    rows = lsc._read_csv(design / "adni_design_pvalues.csv")
    original = rows[1]["p_value"]
    rows[1]["p_value"] = "0.123"
    lsc._write_csv(design / "adni_design_pvalues.csv", rows, adc.PVALUE_COLUMNS)
    settings = _simulation(design_dir, design_dir=design)
    result = gain.simulation_preflight(tmp_path / "out", settings)
    assert not result["passed"]
    differing = result["units"][0]["differing"]
    assert differing == [{"statistic": rows[1]["statistic"], "committed": "0.123", "observed": float(original)}]
    with pytest.raises(SystemExit, match="pre-flight failed"):
        gain.run_shard(tmp_path / "out", settings, 0, 1)


# ── 4.3 The anchor reproduction ──────────────────────────────────────────────


def _committed_cohort_rows() -> list[dict]:
    return [row for row in lsc._read_csv(gain.ANCHOR_DIR / "rate_table.csv") if row["point"] == "cohort"]


def test_the_anchor_reproduction() -> None:
    committed = _committed_cohort_rows()
    ours = [{"point": "cohort", "cell": "anchor", "arm": row["arm"], "statistic": row["statistic"],
             "count": row["count"], "rejections": row["rejections"]} for row in committed if row["arm"] in gain.ARMS]
    assert {(row["arm"], row["statistic"], row["rejections"]) for row in ours if row["arm"] == "pls_cv"} == {
        ("pls_cv", "delta", "25"), ("pls_cv", "angle", "32"), ("pls_cv", "shape", "28")}
    passed = gain.anchor_reproduction(ours, committed)
    assert passed["passed"] and len(passed["rows"]) == 12
    ours[3]["rejections"] = str(int(ours[3]["rejections"]) + 1)
    failed = gain.anchor_reproduction(ours, committed)
    assert not failed["passed"] and len(failed["differing"]) == 1
    entry = failed["differing"][0]
    assert entry["arm"] == ours[3]["arm"] and entry["statistic"] == ours[3]["statistic"]
    assert entry["observed"]["rejections"] == entry["committed"]["rejections"] + 1
    gated = gain.gate({"plants": {"passed": True}, "simulation": {"passed": True}},
                      {"plants": {"inconsistent_units": 0}, "simulation": {"inconsistent_units": 0}}, failed)
    assert not gated["passed"] and "does not reproduce" in gated["failures"][0]


# ── 5.1 Merge on synthetic records ───────────────────────────────────────────


def test_the_rate_table(settings) -> None:
    records = [
        _synthetic(settings, POINT, "heldout_magnitude_e1", 0, overrides={"pls_cv": {"delta": 0.05}}),
        _synthetic(settings, POINT, "heldout_magnitude_e1", 1, overrides={"pls_cv": {"delta": 0.051}}),
        {"profile": "plants", "point": POINT, "cell": "heldout_magnitude_e1", "replicate": 2, "status": "failed"},
    ]
    rows = gain.rate_table(records, settings)
    cell = {(row["arm"], row["statistic"]): row for row in rows if row["cell"] == "heldout_magnitude_e1"}
    assert cell[("pls_cv", "delta")]["rejections"] == 1 and cell[("pls_cv", "delta")]["count"] == 2
    assert cell[("pls_cv", "delta")]["failures"] == 1 and cell[("pls_cv", "delta")]["role"] == "target"
    assert cell[("pls_cv", "angle")]["role"] == "off_target" and cell[("pls_cv", "angle")]["reference"] == "heldout"
    expected = sum(len(gain.point_cells(settings, point)) for point in gain.point_names(settings))
    assert len(rows) == expected * 4 * 3
    roles = {row["role"] for row in rows if row["cell"] in ("anchor", "sample_translation_e1")}
    assert roles == {"control"}


def test_the_paired_table_names_every_reference(settings) -> None:
    records = [_synthetic(settings, POINT, "sample_orientation_e1", replicate,
                          overrides={"predictive_orthonormal": {"angle": 0.01}} if replicate < 3 else {})
               for replicate in range(5)]
    rows = [row for row in gain.paired_table(records, settings)
            if row["cell"] == "sample_orientation_e1" and row["point"] == POINT]
    assert {(row["arm"], row["reference_arm"], row["comparison"]) for row in rows} == set(gain.PAIRS)
    angle = {row["arm"]: row for row in rows if row["statistic"] == "angle" and row["reference_arm"] == "pls_cv"}
    entry = angle["predictive_orthonormal"]
    assert (entry["b"], entry["c"], entry["n_paired"]) == (0, 3, 5) and entry["role"] == "target"
    assert entry["p_gain"] == pytest.approx(0.125) and entry["p_loss"] == 1.0
    draws = {row["arm"]: row for row in rows if row["statistic"] == "angle" and row["comparison"] == "draws"}
    assert draws["predictive_orthonormal__within_stage"]["reference_arm"] == "predictive_orthonormal"
    assert draws["predictive_orthonormal__within_stage"]["b"] == 3


def test_no_discordant_units_give_paired_pvalues_of_one(settings) -> None:
    records = [_synthetic(settings, POINT, "sample_shape_e1", replicate, p=0.01) for replicate in range(4)]
    rows = [row for row in gain.paired_table(records, settings) if row["cell"] == "sample_shape_e1"
            and row["point"] == POINT]
    assert rows and all((row["b"], row["c"]) == (0, 0) for row in rows)
    assert all(row["p_loss"] == row["p_gain"] == row["p_two_sided"] == 1.0 for row in rows)


def test_the_cell_table(settings) -> None:
    records = [_synthetic(settings, POINT, "anchor", replicate, rank=2 + replicate) for replicate in range(3)]
    records += [_synthetic(settings, POINT, "sample_magnitude_e1", replicate) for replicate in range(2)]
    rows = gain.cell_rows(records, settings)

    def value(cell, quantity, **keys):
        found = [row["value"] for row in rows if row["cell"] == cell and row["quantity"] == quantity and
                 row["point"] == POINT and all(row.get(key) == item for key, item in keys.items())]
        assert len(found) == 1, (cell, quantity, keys)
        return found[0]

    assert value("anchor", "cv_rank_units", key="3") == 1 and value("anchor", "cv_rank_median") == 3
    assert value("anchor", "reliability_median") == pytest.approx(0.3)
    assert value("anchor", "reliability_full_reference_median") == pytest.approx(0.6 / 1.3)
    assert value("anchor", "reliability_eigenvalue_ratio_median") == pytest.approx(0.05)
    assert value("anchor", "rank_class_count", arm="pls_cv", statistic="delta", key="k_minus_1") == 1
    assert value("anchor", "rank_class_count", arm="pls_cv", statistic="delta", key="above_k_minus_1") == 2
    assert value("sample_magnitude_e1", "plant_norm_median", key="MCI") == 0.5
    assert value("sample_magnitude_e1", "scale_drift_max") == pytest.approx(0.02)
    assert value("sample_magnitude_e1", "capture_median", key="plant:predictive") == pytest.approx(0.4)
    assert value("sample_magnitude_e1", "capture_median", key="reference:span") == pytest.approx(0.8)
    assert value("sample_magnitude_e1", "cv_auroc_median") == pytest.approx(0.6)
    accounting = gain.accounting(records, settings, duplicates=1)
    by_cell = {(entry["point"], entry["cell"]): entry for entry in accounting["by_cell"]}
    assert by_cell[(POINT, "anchor")]["complete"] == 3 and by_cell[(POINT, "anchor")]["expected"] == 2
    assert accounting["duplicate_records"] == 1 and accounting["label_redraws"] == 5


def test_merge_and_decide_on_the_stand_in(tmp_path, cohort, settings, sim_settings, export_dir) -> None:
    out_dir = tmp_path / "out"
    small = {**settings, "cells": [cell for cell in settings["cells"] if cell["name"] in (
        "anchor", "sample_translation_e1", "heldout_magnitude_e1", "sample_out_of_span_e0.5")]}
    gain.plants_preflight(out_dir, small, cohort)
    gain.run_shard(out_dir, small, 0, 1, cohort=cohort, workers=2)
    merged = gain.merge(out_dir, "plants", export_dir)
    assert merged["complete_units"] == merged["expected_units"] == 2 * (4 + 3)
    assert merged["inconsistent_units"] == 0 and merged["reproduction_passed"] is False
    sim_small = {**sim_settings, "cells": sim_settings["cells"][:3]}
    gain.simulation_preflight(out_dir, sim_small)
    gain.run_shard(out_dir, sim_small, 0, 1)
    assert gain.merge(out_dir, "simulation")["complete_units"] == 6
    result = gain.decide(out_dir, export_dir)
    assert result["verdict"] == "gate_failed" and result["production_statement"] is None
    assert result["indicative"]["indicative"]
    assert set(result["criteria"]) == {"S1_production", "S1_candidate", "S2", "S3", "G", "recorded_gains"}
    assert (out_dir / "plants" / "reproduction.json").exists()
    for name in ("rate_table.csv", "paired_table.csv", "cell_table.csv", "accounting.json", "reproduction.json"):
        text = (out_dir / "plants" / name).read_text()
        assert not any(identifier in text for identifier in cohort.metadata["rid"])


# ── 6.1–6.4 The decision on synthetic tables ─────────────────────────────────


def _tables(plants, simulation, *, rates=None, paired=None) -> dict:
    """Every reading of both profiles: controls at 25 of 500, plants at 10 of 200, no discordant pair."""

    rates, paired = rates or {}, paired or {}
    out = {}
    for settings in (plants, simulation):
        rate_rows, paired_rows = [], []
        for point in gain.point_names(settings):
            for cell in gain.point_cells(settings, point):
                n = gain.cell_replicates(settings, cell)
                for statistic in gain.STATISTICS:
                    for arm in gain.ARMS:
                        hits = rates.get((point, cell["name"], arm, statistic), round(0.05 * n))
                        rate_rows.append({"profile": settings["profile"], "point": point, "cell": cell["name"],
                                          "kind": cell["kind"], "arm": arm, "statistic": statistic, "count": n,
                                          "rejections": hits})
                    for arm, reference, comparison in gain.PAIRS:
                        b, c = paired.get((point, cell["name"], statistic, arm), (0, 0))
                        paired_rows.append({"profile": settings["profile"], "point": point, "cell": cell["name"],
                                            "kind": cell["kind"], "statistic": statistic, "arm": arm,
                                            "reference_arm": reference, "comparison": comparison, "n_paired": n,
                                            "rate_reference": 0.05, "rate_arm": 0.05, "difference": (c - b) / n,
                                            "b": b, "c": c, **plsc.paired_p_values(b, c)})
        out[settings["profile"]] = {"rates": rate_rows, "paired": paired_rows, "cells": []}
    return out


def _decide(*, rates=None, paired=None, gate_ok=True, plants=None, simulation=None) -> dict:
    plants, simulation = plants or _paper_plants(), simulation or _paper_simulation()
    tables = _tables(plants, simulation, rates=rates, paired=paired)
    summaries = {settings["profile"]: {"inconsistent_units": 0, "by_cell": [
        {"point": point, "cell": cell["name"], "usable": gain.cell_replicates(settings, cell)}
        for point in gain.point_names(settings) for cell in gain.point_cells(settings, point)]}
        for settings in (plants, simulation)}
    preflights = {"plants": {"passed": gate_ok}, "simulation": {"passed": True}}
    return gain.decision({"plants": plants, "simulation": simulation}, tables, summaries, preflights,
                         {"passed": True, "differing": []})


#: One simulated target gain at the family threshold, so G holds by default in the verdict tests.
SIM_GAIN = {("table_adni", "population_orientation_e1", "angle", "predictive_orthonormal"): (0, 12)}


def test_the_plan_has_the_designs_families() -> None:
    plan = gain.criteria_plan(_paper_plants(), _paper_simulation())
    sizes = {name: len(readings) for name, readings in plan.items()}
    assert sizes == {"S1_plants": 9, "S1_simulation": 6, "S2": 30, "S3_simulation": 16, "S3_in_sample": 32, "G": 8,
                     "gain_heldout": 6, "gain_in_sample": 16, "heldout_off_target": 12, "out_of_span": 12}
    assert {(point, cell) for _, point, cell, _ in plan["S1_plants"]} == {
        ("n750", "anchor"), ("n750", "sample_translation_e1"), ("cohort", "sample_translation_e1")}
    assert not any(cell.startswith("sample_out_of_span") for readings in plan.values() for _, _, cell, _ in readings
                   if readings is not plan["out_of_span"])


def test_the_size_rule_at_500() -> None:
    assert grade.exceedance_threshold(500) == 34 and grade.exceedance_threshold(200) == 16
    assert grade.ceiling_count(500, 9) == 38 and grade.ceiling_count(500, 6) == 37
    assert grade.ceiling_count(200, 32) == 20


@pytest.mark.parametrize("family, readings, holds", [
    ("plants", {("n750", "anchor", "delta"): 35}, True),
    ("plants", {("n750", "anchor", "delta"): 35, ("cohort", "sample_translation_e1", "shape"): 35}, False),
    ("plants", {("cohort", "sample_translation_e1", "angle"): 39}, False),
    ("plants", {("cohort", "sample_translation_e1", "angle"): 38}, True),
    ("simulation", {("table_adni", "anchor", "shape"): 38}, False),
    ("simulation", {("table_adni", "anchor", "shape"): 37}, True),
    ("simulation", {("table_adni", "anchor", "shape"): 35, ("table_adni", "population_translation_e1", "delta"): 35},
     False),
])
def test_s1_families(family, readings, holds) -> None:
    for arm in ("pls_cv", "predictive_orthonormal"):
        rates = {(point, cell, arm, statistic): hits for (point, cell, statistic), hits in readings.items()}
        decided = _decide(rates=rates, paired=SIM_GAIN)
        key = "S1_production" if arm == "pls_cv" else "S1_candidate"
        judged = decided["criteria"][key][family]
        assert judged["holds"] is holds, (arm, judged["failing_readings"])
        assert decided["criteria"][key]["holds"] is holds
        assert decided["verdict"] == ("reopened" if holds else
                                      "production_excess" if arm == "pls_cv" else "size_excess")


def test_an_unavailable_size_reading_fails_its_family() -> None:
    plants = _paper_plants()
    tables = _tables(plants, _paper_simulation(), paired=SIM_GAIN)
    tables["plants"]["rates"] = [row for row in tables["plants"]["rates"]
                                 if not (row["cell"] == "anchor" and row["point"] == "n750"
                                         and row["arm"] == "predictive_orthonormal" and row["statistic"] == "angle")]
    plan = gain.criteria_plan(plants, _paper_simulation())
    judged = gain.judge_size(gain.index_rates(tables["plants"]["rates"]), plan["S1_plants"], "predictive_orthonormal")
    assert not judged["holds"] and judged["failing_readings"][0]["status"] == "unavailable"


def _losses(cells, statistic_of, value) -> dict:
    return {(point, cell, statistic_of(cell), "predictive_orthonormal"): value for point, cell in cells}


S2_CELLS = [("table_adni", "population_magnitude_e1"), ("n750", "heldout_orientation_e1"),
            ("cohort", "sample_shape_e1")]


def _target(cell: str) -> str:
    return {"magnitude": "delta", "orientation": "angle", "shape": "shape"}[cell.split("_")[1]]


@pytest.mark.parametrize("cells, pair, holds", [
    (S2_CELLS[:2], (6, 0), True),
    (S2_CELLS, (6, 0), False),
    (S2_CELLS[:1], (10, 0), False),
])
def test_s2_protected_power(cells, pair, holds) -> None:
    paired = {**SIM_GAIN, **_losses(cells, _target, pair)}
    decided = _decide(paired=paired)
    judged = decided["criteria"]["S2"]
    assert judged["m"] == 30 and judged["allowed"] == 2 and judged["holds"] is holds
    assert decided["verdict"] == ("reopened" if holds else "power_loss")
    if not holds:
        assert {row["cell"] for row in decided["named_readings"]} == {cell for _, cell in cells}


def _off(cell: str) -> str:
    return "angle" if _target(cell) == "delta" else "delta"


@pytest.mark.parametrize("family, cells, pair, holds", [
    ("in_sample", [("n750", "sample_magnitude_e1"), ("cohort", "sample_orientation_e0.5")], (0, 6), True),
    ("in_sample", [("n750", "sample_magnitude_e1"), ("cohort", "sample_orientation_e0.5"),
                   ("n750", "sample_shape_e1")], (0, 6), False),
    ("in_sample", [("cohort", "sample_shape_e0.5")], (0, 11), False),
    ("simulation", [("table_adni", "population_magnitude_e1"), ("table_adni", "population_shape_e2")], (0, 6),
     True),
    ("simulation", [("table_adni", "population_magnitude_e1"), ("table_adni", "population_shape_e2"),
                    ("table_adni", "population_orientation_e0.5")], (0, 6), False),
    ("simulation", [("table_adni", "population_shape_e1")], (0, 10), False),
])
def test_s3_specificity(family, cells, pair, holds) -> None:
    paired = {**SIM_GAIN, **_losses(cells, _off, pair)}
    decided = _decide(paired=paired)
    judged = decided["criteria"]["S3"][family]
    assert judged["m"] == (32 if family == "in_sample" else 16) and judged["allowed"] == 2
    assert judged["holds"] is holds and decided["criteria"]["S3"]["holds"] is holds
    assert decided["verdict"] == ("reopened" if holds else "specificity_excess")


def test_heldout_off_target_and_out_of_span_readings_change_no_criterion() -> None:
    base = _decide(paired=SIM_GAIN)
    noisy = {("n750", "heldout_magnitude_e1", "angle", "predictive_orthonormal"): (0, 30),
             ("n750", "heldout_orientation_e1", "delta", "predictive_orthonormal"): (0, 30),
             ("cohort", "sample_out_of_span_e0.5", "shape", "predictive_orthonormal"): (0, 40),
             ("n750", "sample_out_of_span_e0.25", "angle", "predictive_orthonormal"): (40, 0)}
    rates = {("cohort", "sample_out_of_span_e0.5", "pls_cv", "delta"): 150,
             ("n750", "heldout_magnitude_e2", "pls_cv", "shape"): 120}
    changed = _decide(paired={**SIM_GAIN, **noisy}, rates=rates)
    assert changed["verdict"] == base["verdict"] == "reopened"
    for name in ("S1_production", "S1_candidate", "S2", "S3", "G"):
        assert changed["criteria"][name] == base["criteria"][name], name
    assert changed["production_statement"]["off_target"] == base["production_statement"]["off_target"]


def test_a_heldout_gain_without_a_simulated_one_fails_g_and_is_recorded() -> None:
    decided = _decide(paired={("n750", "heldout_magnitude_e1", "delta", "predictive_orthonormal"): (0, 12)})
    assert decided["verdict"] == "no_gain" and not decided["criteria"]["G"]["holds"]
    assert decided["criteria"]["G"]["m"] == 8 and decided["criteria"]["G"]["threshold"] == pytest.approx(0.05 / 8)
    recorded = decided["criteria"]["recorded_gains"]
    assert recorded["heldout"]["holds"] and [row["cell"] for row in recorded["heldout"]["gains"]] == [
        "heldout_magnitude_e1"]
    assert recorded["heldout"]["m"] == 6 and recorded["in_sample"]["m"] == 16 and not recorded["in_sample"]["holds"]
    assert "not a qualification" in decided["consequence"]


def test_the_seven_verdicts_and_their_consequences() -> None:
    loss = _losses(S2_CELLS[:1], _target, (10, 0))
    excess = _losses([("n750", "sample_shape_e1")], _off, (0, 11))
    cases = {
        "gate_failed": _decide(gate_ok=False, paired=SIM_GAIN),
        "production_excess": _decide(rates={("n750", "anchor", "pls_cv", "delta"): 40}, paired=SIM_GAIN),
        "size_excess": _decide(rates={("cohort", "sample_translation_e1", "predictive_orthonormal", "angle"): 40},
                               paired=SIM_GAIN),
        "power_loss": _decide(paired={**SIM_GAIN, **loss}),
        "specificity_excess": _decide(paired={**SIM_GAIN, **excess}),
        "no_gain": _decide(),
        "reopened": _decide(paired=SIM_GAIN),
    }
    assert set(cases) == set(gain.VERDICTS)
    for verdict, decided in cases.items():
        assert decided["verdict"] == verdict
        assert decided["consequence"] == f"{gain.CONSEQUENCES[verdict]} {gain.STANDING}"
        assert decided["consequence"].endswith("the out-of-span cells are reported beside the verdict and were not "
                                               "judged.")
        assert set(decided["licenses"]) == {"L1", "L2", "L3", "L4", "L5"}
        assert set(decided["scoping_rule"]) == {"S1", "S2", "S3", "S4", "S5"}
        assert decided["rule_version"] == 1 and not decided["indicative"]["indicative"]
    assert cases["gate_failed"]["criteria"]["G"]["holds"] and cases["gate_failed"]["failing_criteria"] == []
    assert "opt-in diagnostic" in cases["reopened"]["consequence"]
    for verdict in ("size_excess", "power_loss", "specificity_excess"):
        assert "not met for ADNI" in cases[verdict]["consequence"]
    assert cases["production_excess"]["named_readings"][0]["cell"] == "anchor"


def test_the_first_failing_criterion_names_the_verdict_and_all_are_listed() -> None:
    excess = _losses([("n750", "sample_shape_e1")], _off, (0, 11))
    decided = _decide(paired=excess)  # S3 and G fail
    assert decided["verdict"] == "specificity_excess" and decided["failing_criteria"] == ["S3", "G"]
    assert [row["cell"] for row in decided["named_readings"]] == ["sample_shape_e1"]
    gated = _decide(gate_ok=False, paired=excess)
    assert gated["verdict"] == "gate_failed" and gated["failing_criteria"] == ["S3", "G"]
    assert gated["criteria"]["S3"]["in_sample"]["observed"] == 1
    assert gated["production_statement"] is None


def test_a_decision_below_paper_grade_is_indicative() -> None:
    plants = {**_paper_plants(), "replicates": {"control": 4, "plant": 4}, "permutations": 999}
    decided = _decide(plants=plants, paired=SIM_GAIN)
    assert decided["indicative"]["indicative"]
    assert any("replicates per control" in reason for reason in decided["indicative"]["reasons"])
    simulation = {**_paper_simulation(), "permutations": 99}
    assert any("99 permutations" in reason for reason in
               _decide(simulation=simulation, paired=SIM_GAIN)["indicative"]["reasons"])


def test_informative_from_on_ladders_that_reach_the_floor_or_not() -> None:
    rates = {("n750", "heldout_magnitude_e0.5", "pls_cv", "delta"): 100,
             ("n750", "heldout_magnitude_e1", "pls_cv", "delta"): 165,
             ("n750", "heldout_magnitude_e2", "pls_cv", "delta"): 200,
             ("cohort", "sample_orientation_e0.25", "pls_cv", "angle"): 20,
             ("cohort", "sample_orientation_e0.5", "pls_cv", "angle"): 60,
             ("cohort", "sample_orientation_e1", "pls_cv", "angle"): 40}
    statement = _decide(rates=rates, paired=SIM_GAIN)["production_statement"]
    power = {(entry["point"], entry["reference"], entry["mode"]): entry for entry in statement["power"]}
    assert len(power) == 3 + 3 + 3
    held = power[("n750", "heldout", "magnitude")]
    assert held["informative_from"] == 1.0 and held["informative_rate"] == pytest.approx(0.825)
    assert "rejects a held-out stretch of e = 1 in 165 of 200" in held["statement"]
    assert "bounds the detectable effect from above" in held["statement"]
    rotation = power[("cohort", "sample", "orientation")]
    assert rotation["informative_from"] is None and rotation["max_rate"] == pytest.approx(0.3)
    assert rotation["effect_at_max"] == 0.5 and "stays below the 0.80 floor" in rotation["statement"]
    assert "bounds it from below" in rotation["statement"]


@pytest.mark.parametrize("hits, verdict", [
    ({("n750", "sample_magnitude_e1", "angle"): 17, ("cohort", "sample_shape_e1", "delta"): 17}, "specific"),
    ({("n750", "sample_magnitude_e1", "angle"): 17, ("cohort", "sample_shape_e1", "delta"): 17,
      ("cohort", "sample_orientation_e1", "shape"): 17}, "leaks"),
    ({("cohort", "sample_orientation_e1", "shape"): 21}, "leaks"),
    ({("cohort", "sample_orientation_e1", "shape"): 20}, "specific"),
])
def test_productions_off_target_family(hits, verdict) -> None:
    rates = {(point, cell, "pls_cv", statistic): value for (point, cell, statistic), value in hits.items()}
    decided = _decide(rates=rates, paired=SIM_GAIN)
    off_target = decided["production_statement"]["off_target"]
    assert off_target["family"]["m"] == 32 and off_target["family"]["allowed_exceedances"] == 2
    assert off_target["verdict"] == verdict and decided["verdict"] == "reopened"
    if verdict == "leaks":
        assert {(row["point"], row["cell"], row["statistic"]) for row in off_target["named_readings"]} <= set(hits)


def test_descriptive_readings_never_change_the_verdict_or_the_statement() -> None:
    plants, simulation = _paper_plants(), _paper_simulation()
    base_tables = _tables(plants, simulation, paired=SIM_GAIN)
    summaries = {"plants": {"inconsistent_units": 0, "by_cell": []},
                 "simulation": {"inconsistent_units": 0, "by_cell": []}}
    preflights = {"plants": {"passed": True}, "simulation": {"passed": True}}
    reproduction = {"passed": True, "differing": []}
    base = gain.decision({"plants": plants, "simulation": simulation}, base_tables, summaries, preflights,
                         reproduction)
    tables = copy.deepcopy(base_tables)
    tables["plants"]["cells"] = [
        {"profile": "plants", "point": "n750", "cell": "anchor", "quantity": "reliability_median", "arm": "",
         "statistic": "", "key": "", "value": 0.9},
        {"profile": "plants", "point": "n750", "cell": "sample_out_of_span_e0.5", "quantity": "capture_median",
         "arm": "", "statistic": "", "key": "plant:span", "value": 0.99},
    ]
    for row in tables["plants"]["rates"]:
        if row["arm"].endswith("__within_stage"):
            row["rejections"] = 150
    for row in tables["plants"]["paired"]:
        if row["comparison"] == "draws":
            row.update({"b": 50, "c": 0, **plsc.paired_p_values(50, 0)})
    changed = gain.decision({"plants": plants, "simulation": simulation}, tables, summaries, preflights,
                            reproduction)
    assert changed["verdict"] == base["verdict"] == "reopened"
    assert changed["production_statement"] == base["production_statement"]
    assert changed["criteria"] == base["criteria"]
    assert not changed["descriptive"]["within_stage"]["size"]["pls_cv__within_stage"]["plants"]["holds"]
    assert changed["descriptive"]["heldout_off_target"][0]["reliability_median"] == 0.9
    bracket = changed["descriptive"]["bracket"]
    assert {(row["mode"], row["effect"]) for row in bracket} == {("magnitude", 0.5), ("magnitude", 1.0),
                                                                 ("orientation", 0.5), ("orientation", 1.0),
                                                                 ("shape", 1.0)}
