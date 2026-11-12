"""``scripts/predictive_subspace_check.py``: settings, units, seeds, parents, arms, checks, merge, decision.

The simulation units run on the driver's ``seaad`` designs at fast CV knobs; the
plasmode profile runs on the plasmode tests' synthetic stand-in export. No real
data are read.
"""

from __future__ import annotations

import dataclasses
import json
import math
import shutil
import sys
from pathlib import Path

import numpy as np
import pytest
from scipy.stats import binom


def _script(name: str):
    sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
    try:
        return __import__(name)
    finally:
        sys.path.pop(0)


psc = _script("predictive_subspace_check")
lsc, ssc, plsc, pc = psc.lsc, psc.ssc, psc.plsc, psc.pc

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


def _fast(settings: dict, **overrides) -> dict:
    """The settings at fast CV knobs (and any other unit-setting overrides)."""

    settings = json.loads(json.dumps(settings))
    units = psc._harness_units(settings)
    units["pls_cv_params"] = dict(FAST_CV)
    units.update(overrides)
    return settings


# ── 3.1 / 3.2 Settings and units ─────────────────────────────────────────────


@pytest.fixture(scope="module")
def full_simulation():
    return psc.simulation_settings()


@pytest.fixture(scope="module")
def full_weak():
    return psc.weak_size_settings()


def test_unit_counts(full_simulation, full_weak, stand_in) -> None:
    assert len(psc.simulation_units(full_simulation)) == 7_800
    assert len(psc.weak_size_units(full_weak)) == 1_000
    assert len(psc.plasmode_units(_plasmode_settings(stand_in, replicates=500, plant_replicates=100))) == 800


def test_simulation_designs_cells_and_grade(full_simulation) -> None:
    settings = full_simulation
    assert list(settings["design_cells"]) == list(psc.SIMULATION_DESIGNS)
    standard = [cell[0] for cell in psc.CELLS]
    assert standard == [
        "anchor", "translation", "magnitude_e0.05", "magnitude_e0.1", "magnitude_e0.25", "orient_e0.1_phi0",
        "orient_e0.25_phi0", "orient_e0.5_phi0", "shape_e0.1", "shape_e0.25", "shape_e0.5",
    ]
    for design in psc.STAGE_SPAN_DESIGNS:
        assert settings["design_cells"][design] == standard
    assert settings["design_cells"][psc.SURGERY_DESIGN] == ["surgery_orient_e1"]
    assert settings["replicates"] == 100 and settings["permutations"] == 199
    assert settings["units"]["base_seed"] == 2036 and settings["rejection"] == "p < alpha"
    assert settings["arms"] == list(psc.ARMS)
    assert settings["stage_effects"] == {
        "n200": 2.0, "seaad": 2.0, "n200_delta1": 1.0, "n200_delta0.5": 0.5, "n200_delta0.25": 0.25,
        "seaad_delta0.5": 0.5, "seaad_weak": 0.11, "seaad_surgery": 2.0,
    }


def test_designs_differ_from_the_comparisons_only_in_effect_scale_and_construction(full_simulation) -> None:
    parent = json.loads((psc.COMPARISON_DIR / "settings.json").read_text())["settings"]
    for design, (source, delta) in psc.DESIGN_SOURCES.items():
        ours = full_simulation["units"]["designs"][design]["generator"]
        theirs = parent["designs"][source]["generator"]
        differing = {key for key in set(ours) | set(theirs) if ours.get(key) != theirs.get(key)}
        if design == psc.SURGERY_DESIGN:
            assert differing == {"construction", "span_layers"}
            assert "span_scale" not in ours
            continue
        expected = {"span_scale"} | (set(psc.STAGE_EFFECT_KEYS) if delta not in (None, 2.0) else set())
        assert differing == expected, design
        assert ours["span_scale"] == "planted"
        if delta not in (None, 2.0):
            assert all(ours[key] == delta for key in psc.STAGE_EFFECT_KEYS)


def test_the_weak_delta_is_the_committed_calibrations(full_simulation) -> None:
    calibration = json.loads((psc.COMPARISON_DIR / "calibration" / "calibration.json").read_text())
    weak = full_simulation["units"]["weak_delta"]
    assert weak["value"] == calibration["adopted"]["delta"] == 0.11
    generator = full_simulation["units"]["designs"]["seaad_weak"]["generator"]
    assert all(generator[key] == 0.11 for key in psc.STAGE_EFFECT_KEYS)
    assert full_simulation["parents"]["comparison_calibration"]["sha256"] == weak["calibration_sha256"]


def test_a_fresh_seed_family(full_simulation) -> None:
    ours = {
        lsc.generator_seed(2036, design, replicate)
        for design in psc.SIMULATION_DESIGNS
        for replicate in range(100)
    }
    spent = {
        lsc.generator_seed(base, design, replicate)
        for base in (2035, 4242)
        for design in (*lsc.DESIGNS, *psc.SIMULATION_DESIGNS)
        for replicate in range(500)
    }
    assert len(ours) == len(psc.SIMULATION_DESIGNS) * 100
    assert not ours & spent
    for design in ("n200", "seaad", "seaad_weak"):
        for cell in ("anchor", "translation"):
            assert lsc.rrpp_seed(2036, design, cell, 3) != lsc.rrpp_seed(2035, design, cell, 3)


def test_cells_of_one_design_and_replicate_share_their_null_dataset(full_simulation) -> None:
    units = full_simulation["units"]
    anchor = lsc.unit_dataset(lsc.Unit("seaad", "anchor", 4), units)
    rotated = lsc.unit_dataset(lsc.Unit("seaad", "orient_e0.25_phi0", 4), units)
    group_a = anchor.metadata["group"].astype(str) == "A"
    assert anchor.metadata.equals(rotated.metadata)
    for layer in ("methylation", "expression"):
        rows = group_a.to_numpy()
        np.testing.assert_array_equal(getattr(anchor, layer).to_numpy()[rows], getattr(rotated, layer).to_numpy()[rows])
    assert rotated.truth["span_scale"] == "planted"
    surgery = lsc.unit_dataset(lsc.Unit(psc.SURGERY_DESIGN, "surgery_orient_e1", 0), units)
    assert surgery.truth.get("construction") != "stage_span"
    assert surgery.truth["trajectory_mode"] == "orientation"


def test_weak_size_units_are_the_span_split_checks(full_weak) -> None:
    split = json.loads((psc.SPAN_SPLIT_DIR / "settings.json").read_text())["settings"]["comparison"]
    units = psc.weak_size_units(full_weak)
    assert {(unit.design, unit.cell) for unit in units} == {("seaad_weak", "anchor"), ("seaad_weak", "translation")}
    assert full_weak["permutations"] == 999 and full_weak["units"]["comparison"]["permutations"] == 999
    for unit in (units[0], units[517], units[-1]):
        ours = psc._harness_units(full_weak)
        assert lsc._generator_params(
            ours["designs"][unit.design]["generator"],
            seed=lsc.generator_seed(ours["base_seed"], unit.design, unit.replicate),
        ) == lsc._generator_params(
            split["designs"][unit.design]["generator"],
            seed=lsc.generator_seed(split["base_seed"], unit.design, unit.replicate),
        )
        assert lsc.rrpp_seed(ours["base_seed"], unit.design, unit.cell, unit.replicate) == lsc.rrpp_seed(
            split["base_seed"], unit.design, unit.cell, unit.replicate
        )
        assert lsc.cell_spec(ours, unit.cell) == lsc.cell_spec(split, unit.cell)
    assert set(full_weak["parents"]) == {"span_split_settings", "span_split_production", "comparison_calibration"}


def _comparison_copy(tmp_path: Path) -> Path:
    target = tmp_path / "comparison"
    (target / "calibration").mkdir(parents=True)
    for name in ("settings.json", "rate_table.csv", "calibration/calibration.json"):
        shutil.copy(psc.COMPARISON_DIR / name, target / name)
    return target


def test_a_changed_parent_is_refused_by_name(tmp_path) -> None:
    parent = _comparison_copy(tmp_path)
    settings = psc.simulation_settings(parent_dir=parent, replicates=1, designs=["seaad"])
    plsc.verify_parents(settings)
    (parent / "settings.json").write_text((parent / "settings.json").read_text() + "\n")
    with pytest.raises(SystemExit, match=r"Parent comparison_settings .* changed"):
        plsc.verify_parents(settings)
    with pytest.raises(SystemExit, match="Parent comparison_settings"):
        psc.run_shard(tmp_path / "run", settings, 0, 1, limit=0)


def test_an_unknown_design_is_refused() -> None:
    with pytest.raises(SystemExit, match=r"Unknown design\(s\) \['n80'\]"):
        psc.simulation_settings(designs=["n80"])


# ── 3.4 One unit ─────────────────────────────────────────────────────────────


@pytest.fixture(scope="module")
def small_settings():
    return _fast(psc.simulation_settings(replicates=2, permutations=19, designs=["seaad", "seaad_weak"]))


@pytest.fixture(scope="module")
def unit_record(small_settings):
    return psc.evaluate_harness_unit(lsc.Unit("seaad", "orient_e0.25_phi0", 1), small_settings)


def test_the_record_carries_every_field(unit_record, small_settings) -> None:
    record = unit_record
    assert record["status"] == "completed", record.get("inconsistent_checks") or record.get("error_message")
    for key in ("profile", "design", "cell", "trajectory_mode", "effect_size", "replicate", "generator_seed",
                "rrpp_seed", "settings_hash", "consistency", "production", "arms", "runtime_seconds"):
        assert key in record, key
    assert record["generator_seed"] == lsc.generator_seed(2036, "seaad", 1)
    assert record["rrpp_seed"] == lsc.rrpp_seed(2036, "seaad", "orient_e0.25_phi0", 1)
    assert list(record["arms"]) == list(psc.ARMS)
    for arm, entry in record["arms"].items():
        assert entry["status"] == "completed", (arm, entry.get("error_message"))
        for key in ("rank", "measured_dimension", "cv_mean_auroc", "pair_statistics", "p_values", "angle_null_q95",
                    "span_capture", "population_latent"):
            assert key in entry, (arm, key)
    rank = record["production"]["rank"]
    arms = record["arms"]
    for arm in ("pls_cv", *psc.CANDIDATE_ARMS):
        assert arms[arm]["rank"] == rank and arms[arm]["cv_mean_auroc"] == record["production"]["cv_mean_auroc"]
    assert arms["pls_cv"]["measured_dimension"] == rank
    assert all(arms[arm]["measured_dimension"] == 2 for arm in psc.CANDIDATE_ARMS)
    assert arms["forced_k_minus_1"]["rank"] == 2 and arms["forced_k_minus_1"]["measured_dimension"] == 2
    assert arms["predictive_rank_plus_one"]["rank"] == rank + 1
    assert arms["predictive_rank_plus_one"]["measured_dimension"] == 2
    assert record["production"]["loadings_rank"] == min(rank, 2)
    assert 0.0 < record["production"]["stage_r2"] <= 1.0
    for name in psc.CONSISTENCY_CHECKS:
        assert record["consistency"][name]["passed"], name
    assert record["consistency"]["one_subspace"]["applicable"] is True


def test_the_arms_equal_plain_evaluations(unit_record, small_settings) -> None:
    from motco.simulations.evaluation import evaluate_semisynthetic_trajectory

    units = small_settings["units"]
    unit = lsc.Unit("seaad", "orient_e0.25_phi0", 1)
    dataset = lsc.unit_dataset(unit, units)
    params = lsc.arm_params("pls_cv", units, seed=unit_record["rrpp_seed"], n_stages=3)
    plain = psc.harness_arm_record(evaluate_semisynthetic_trajectory(dataset, params), 0.0)
    for value in ("pls_cv", *psc.CANDIDATE_ARMS):
        expected = plain if value == "pls_cv" else psc.harness_arm_record(
            evaluate_semisynthetic_trajectory(dataset, dataclasses.replace(params, latent_coordinates=value)), 0.0
        )
        got = json.loads(json.dumps(lsc._jsonable(
            {key: item for key, item in unit_record["arms"][value].items() if key != "seconds"}
        )))
        want = json.loads(json.dumps(lsc._jsonable({key: item for key, item in expected.items() if key != "seconds"})))
        assert got == want, value


def test_one_cross_validation_and_two_forced_fits(small_settings, monkeypatch) -> None:
    import motco.simulations.evaluation as evaluation_module

    calls = {"cv": 0, "fit": 0}
    cv, fit = evaluation_module.plsda_doubleCV, evaluation_module.fit_plsda_model

    def counting_cv(*args, **kwargs):
        calls["cv"] += 1
        return cv(*args, **kwargs)

    def counting_fit(*args, **kwargs):
        calls["fit"] += 1
        return fit(*args, **kwargs)

    monkeypatch.setattr(evaluation_module, "plsda_doubleCV", counting_cv)
    monkeypatch.setattr(evaluation_module, "fit_plsda_model", counting_fit)
    record = psc.evaluate_harness_unit(lsc.Unit("seaad", "anchor", 0), small_settings)
    assert record["status"] == "completed"
    assert calls == {"cv": 1, "fit": 3}


def test_a_patched_basis_is_flagged_and_left_out_of_the_rates(small_settings, monkeypatch) -> None:
    import motco.simulations.evaluation as evaluation_module

    def wrong(model):
        # Two copies of one direction: a rank-one "subspace" that cannot carry the two-dimensional fitted labels.
        first = np.zeros(np.asarray(model.x_rotations_).shape[1])
        first[0] = 1.0
        return np.column_stack([first, first])

    monkeypatch.setattr(evaluation_module, "predictive_directions", wrong)
    record = psc.evaluate_harness_unit(lsc.Unit("seaad", "translation", 0), small_settings)
    assert record["status"] == "inconsistent"
    assert "fitted_labels" in record["inconsistent_checks"]
    assert record["consistency"]["fitted_labels"]["deviation"] > psc.CONSISTENCY_TOLERANCE
    rows = psc.rate_table([record], small_settings)
    assert all(row["count"] == 0 and row["inconsistent"] == 1 for row in rows)


def test_a_failed_unit_is_recorded(small_settings, monkeypatch) -> None:
    def broken(*args, **kwargs):
        raise RuntimeError("generator broke")

    monkeypatch.setattr(lsc, "unit_dataset", broken)
    record = psc.evaluate_harness_unit(lsc.Unit("seaad", "anchor", 0), small_settings)
    assert record["status"] == "failed" and record["error_message"] == "generator broke"


def test_consistency_checks_name_each_failure() -> None:
    class Model:
        x_rotations_ = np.eye(4)[:, :3]
        y_loadings_ = np.array([[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [-1.0, -1.0, 0.0]])

    scores = np.random.default_rng(1).normal(size=(12, 3))
    good = scores @ np.eye(3)[:, :2]
    arms = {
        "pls_cv": {"status": "completed", "rank": 3, "cv_mean_auroc": 0.9, "p_values": {"delta": 0.1}},
        "predictive": {"status": "completed", "rank": 3, "cv_mean_auroc": 0.9, "p_values": {"delta": 0.1},
                       "span_capture": {"pooled_stage": 0.5, "group_difference": 0.2, "rank": 2}},
        "predictive_fitted": {"status": "completed", "rank": 4, "cv_mean_auroc": 0.9, "p_values": {"delta": 0.1},
                              "span_capture": {"pooled_stage": 0.6, "group_difference": 0.2, "rank": 2}},
    }
    checks = psc.consistency_checks(
        model=Model(), scores=scores, n_stages=3, candidates={"predictive": good, "predictive_fitted": scores[:, :1]},
        arms=arms, with_capture=True,
    )
    assert checks["rank"]["passed"]
    assert not checks["fitted_labels"]["passed"]
    assert not checks["same_fit"]["passed"] and checks["same_fit"]["deviation"] == math.inf
    assert not checks["one_subspace"]["passed"] and checks["one_subspace"]["deviation"] == pytest.approx(0.1)
    assert checks["identity_at_low_rank"]["applicable"] is False

    class Shifted(Model):
        x_rotations_ = Model.x_rotations_ + 1e-6

    same_rank = {**arms, "predictive_fitted": {**arms["predictive_fitted"], "rank": 3}}
    refit = psc.consistency_checks(
        model=Model(), scores=scores, n_stages=3, candidates={"predictive": good}, arms=same_rank, with_capture=False,
        models={"predictive": Shifted()},
    )
    assert not refit["same_fit"]["passed"] and refit["same_fit"]["deviation"] == pytest.approx(1e-6)
    shared = psc.consistency_checks(
        model=Model(), scores=scores, n_stages=3, candidates={"predictive": good}, arms=same_rank, with_capture=False,
    )
    assert shared["same_fit"]["passed"] and shared["same_fit"]["deviation"] == 0.0


# ── Resume, refusal and merge on a small run ─────────────────────────────────


def test_resume_merge_and_settings_refusal(tmp_path, small_settings) -> None:
    settings = _fast(psc.simulation_settings(replicates=1, permutations=9, designs=["seaad"]))
    out = tmp_path / "run"
    assert psc.run_shard(out, settings, 0, 1, limit=2) == 2
    assert psc.run_shard(out, settings, 0, 1, limit=1) == 1
    remaining = len(psc.simulation_units(settings)) - 3
    assert psc.run_shard(out, settings, 0, 1) == remaining
    assert psc.run_shard(out, settings, 0, 1) == 0
    summary = psc.merge(out, "simulation")
    assert summary["complete_units"] == 11 and summary["missing_units"] == 0
    assert summary["duplicate_records"] == 0 and summary["reproduction_passed"] is True
    directory = out / "simulation"
    for name in ("rate_table.csv", "paired_table.csv", "cell_table.csv", "accounting.json", "reproduction.json",
                 "merged.jsonl"):
        assert (directory / name).exists(), name
    assert len((directory / "merged.jsonl").read_text().splitlines()) == 11
    other = _fast(psc.simulation_settings(replicates=2, permutations=9, designs=["seaad"]))
    with pytest.raises(SystemExit, match="differing: .*replicates"):
        psc.run_shard(out, other, 0, 1)


# ── 3.3 The plasmode profile ─────────────────────────────────────────────────


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
    shutil.copy(psc.PLASMODE_DIR / "plasmode_summary.csv", parent_dir / "plasmode_summary.csv")
    production = plsc.plasmode_settings(export, parent_dir=parent_dir, replicates=2, permutations=9)
    production_dir = tmp_path_factory.mktemp("production")
    (production_dir / "settings.json").write_text(
        json.dumps({"hash": lsc.settings_hash(production), "settings": production})
    )
    shutil.copy(psc.PRODUCTION_PLASMODE_DIR / "rate_table.csv", production_dir / "rate_table.csv")
    return directory, export, parent_dir, production_dir, production


def _plasmode_settings(stand_in, **kwargs):
    _, export, parent_dir, production_dir, _ = stand_in
    kwargs.setdefault("replicates", 2)
    kwargs.setdefault("plant_replicates", 1)
    kwargs.setdefault("permutations", 9)
    return psc.plasmode_settings(export, parent_dir=parent_dir, production_dir=production_dir, **kwargs)


def test_plasmode_settings_and_units(stand_in) -> None:
    settings = _plasmode_settings(stand_in, replicates=500, plant_replicates=100)
    units = psc.plasmode_units(settings)
    assert sum(unit.cell == "anchor" for unit in units) == 500
    assert sum(unit.cell != "anchor" for unit in units) == 300
    assert {unit.cell for unit in units} == {"anchor", "span_e0.25", "span_e0.5", "span_e1"}
    assert settings["rejection"] == "p <= alpha" and settings["arms"] == list(psc.ARMS[:-1])
    assert set(settings["parents"]) == {
        "plasmode_settings", "plasmode_summary", "production_plasmode_settings", "production_plasmode_rates",
    }
    assert settings["predictive_arm_fit"] == psc.PLASMODE_PREDICTIVE_FIT


def test_a_changed_export_is_refused(stand_in, tmp_path) -> None:
    directory, export, *_ = stand_in
    settings = _plasmode_settings(stand_in)
    altered = dataclasses.replace(export, files={**export.files, "rna.csv": "0" * 64})
    with pytest.raises(SystemExit, match="rna.csv"):
        psc.run_shard(tmp_path / "run", settings, 0, 1, export=altered, limit=0)


def test_a_production_selection_that_differs_is_refused(stand_in, tmp_path) -> None:
    _, export, parent_dir, _, production = stand_in
    changed = {**production, "cells": production["cells"][:2]}
    directory = tmp_path / "production"
    directory.mkdir()
    (directory / "settings.json").write_text(json.dumps({"hash": lsc.settings_hash(changed), "settings": changed}))
    shutil.copy(psc.PRODUCTION_PLASMODE_DIR / "rate_table.csv", directory / "rate_table.csv")
    with pytest.raises(SystemExit, match="differ from this profile's in: cells"):
        psc.plasmode_settings(export, parent_dir=parent_dir, production_dir=directory)


@pytest.mark.parametrize("key", [("plasmode", "anchor", 1), ("plasmode", "span_e0.5", 0)])
def test_plasmode_units_are_the_production_selections(stand_in, key) -> None:
    _, export, _, _, production = stand_in
    settings = _plasmode_settings(stand_in)
    unit = lsc.Unit(*key)
    ours, theirs = plsc.build_plasmode_unit(unit, settings, export), plsc.build_plasmode_unit(unit, production, export)
    assert ours["metadata"].equals(theirs["metadata"])
    for name in pc.BLOCKS:
        assert ours["blocks"][name].equals(theirs["blocks"][name])
    assert ours["seeds"] == theirs["seeds"]
    params = psc.plasmode_production_params(settings, ours["seeds"])
    assert params == plsc.plasmode_arm_params("pls_cv", production, theirs["seeds"])
    assert params.latent_coordinates is None and params.forced_components is None
    predictive = psc.plasmode_arm_params("predictive", params, n_stages=3, production_rank=5)
    assert predictive == dataclasses.replace(params, latent_coordinates="predictive", forced_components=5)
    assert psc.plasmode_arm_params("forced_k_minus_1", params, n_stages=3) == dataclasses.replace(
        params, forced_components=2
    )
    with pytest.raises(ValueError, match="needs the unit's pls_cv rank"):
        psc.plasmode_arm_params("predictive", params, n_stages=3)


@pytest.fixture(scope="module")
def plasmode_record(stand_in):
    _, export, *_ = stand_in
    return psc.evaluate_plasmode_unit(lsc.Unit("plasmode", "anchor", 1), _plasmode_settings(stand_in), export)


def test_plasmode_record_holds_seeds_and_statistics_only(stand_in, plasmode_record) -> None:
    _, export, *_ = stand_in
    record = plasmode_record
    assert record["status"] == "completed", record.get("inconsistent_checks") or record.get("error_message")
    assert list(record["arms"]) == list(psc.ARMS[:-1])
    assert record["seeds"] == pc.unit_seeds(2034, "paper", 1) and record["rrpp_seed"] == record["seeds"]["rrpp"]
    text = json.dumps(record)
    for donor in export.metadata[export.sample_col]:
        assert donor not in text
    assert "Female" not in text and "Male" not in text

    def no_vectors(value):
        if isinstance(value, dict):
            return all(no_vectors(item) for item in value.values())
        if isinstance(value, list):
            return len(value) != len(export.metadata) and all(no_vectors(item) for item in value)
        return True

    assert no_vectors(record)
    assert record["arms"]["pls_cv"]["fidelity"]["passed"] is True
    rank = record["production"]["rank"]
    for arm in psc.CANDIDATE_ARMS:
        assert record["arms"][arm]["rank"] == rank
        assert record["arms"][arm]["measured_dimension"] == min(rank, 2)
        assert record["arms"][arm]["fit"] == "forced_at_pls_cv_rank"
        assert record["arms"][arm]["cv_mean_auroc"] == record["production"]["cv_mean_auroc"]
    assert record["consistency"]["same_fit"]["deviation"] == 0.0
    assert record["arms"]["forced_k_minus_1"]["rank"] == 2
    assert record["consistency"]["same_fit"]["passed"] and record["consistency"]["fitted_labels"]["passed"]
    assert record["consistency"]["one_subspace"]["applicable"] is False


def test_the_forced_predictive_arms_equal_a_cross_validated_analysis(stand_in, plasmode_record) -> None:
    _, export, *_ = stand_in
    settings = _plasmode_settings(stand_in)
    built = plsc.build_plasmode_unit(lsc.Unit("plasmode", "anchor", 1), settings, export)
    production = psc.plasmode_production_params(settings, built["seeds"])
    for arm in psc.CANDIDATE_ARMS:
        _, result = pc.analyse(built["blocks"], built["metadata"], export,
                               dataclasses.replace(production, latent_coordinates=arm))
        expected = json.loads(json.dumps(lsc._jsonable(psc.analysis_arm_record(result, 0.0))))
        ours = plasmode_record["arms"][arm]
        for key in ("rank", "cv_mean_auroc", "measured_dimension", "pair_statistics", "p_values", "angle_null_q95"):
            assert ours[key] == expected[key], (arm, key)


def test_a_plasmode_unit_runs_one_cross_validation(stand_in, monkeypatch) -> None:
    import motco.simulations.evaluation as evaluation_module

    _, export, *_ = stand_in
    calls = []
    original = evaluation_module.plsda_doubleCV

    def counting(*args, **kwargs):
        calls.append(1)
        return original(*args, **kwargs)

    monkeypatch.setattr(evaluation_module, "plsda_doubleCV", counting)
    record = psc.evaluate_plasmode_unit(lsc.Unit("plasmode", "span_e0.5", 0), _plasmode_settings(stand_in), export)
    assert record["status"] == "completed" and len(calls) == 1


def test_the_plasmode_production_arm_equals_the_production_selections(stand_in, plasmode_record) -> None:
    _, export, _, _, production = stand_in
    unit = lsc.Unit("plasmode", "anchor", 1)
    built = plsc.build_plasmode_unit(unit, production, export)
    params = plsc.plasmode_arm_params("pls_cv", production, built["seeds"])
    _, result = pc.analyse(built["blocks"], built["metadata"], export, params)
    expected = plsc.plasmode_arm_record(result)
    ours = plasmode_record["arms"]["pls_cv"]
    for key in ("rank", "cv_mean_auroc", "pair_statistics", "p_values"):
        assert ours[key] == json.loads(json.dumps(lsc._jsonable(expected[key]))), key


# ── 4.1 Merge on synthetic records ───────────────────────────────────────────


def _synthetic_record(unit, settings, *, status="completed", p=0.5, overrides=None):
    base = {
        "profile": settings["profile"], "design": unit.design, "cell": unit.cell,
        "trajectory_mode": "none", "effect_size": 0.0, "replicate": unit.replicate,
        "settings_hash": lsc.settings_hash(settings), "runtime_seconds": {"unit": 10.0, "fit": 4.0},
    }
    if status == "failed":
        return {**base, "status": "failed", "error_type": "RuntimeError", "error_message": "x"}
    arms = {}
    for arm in settings["arms"]:
        arms[arm] = {
            "status": "completed", "rank": 3, "measured_dimension": 2 if arm in psc.CANDIDATE_ARMS else 3,
            "cv_mean_auroc": 0.99, "pair_statistics": {"delta": 0.5, "angle": 10.0, "shape": 0.1},
            "p_values": {"delta": p, "angle": p, "shape": None}, "angle_null_q95": 30.0,
            "span_capture": {"pooled_stage": 0.9, "group_difference": 0.5, "rank": 3},
            "population_latent": {"delta": 0.4, "angle": 12.0, "shape": 0.2}, "seconds": 1.0,
        }
    arms.update(overrides or {})
    checks = {name: {"applicable": True, "passed": True, "deviation": 0.0, "tolerance": 1e-8}
              for name in psc.CONSISTENCY_CHECKS}
    failing = []
    if status == "inconsistent":
        checks["fitted_labels"] = {"applicable": True, "passed": False, "deviation": 0.3, "tolerance": 1e-8}
        failing = ["fitted_labels"]
    return {**base, "status": status, "inconsistent_checks": failing, "consistency": checks,
            "production": {"rank": 3, "cv_mean_auroc": 0.99, "loadings_rank": 2, "stage_r2": 0.9}, "arms": arms}


def test_merge_tables_count_failures_inconsistent_units_and_pairs(small_settings) -> None:
    settings = small_settings
    units = [lsc.Unit("seaad", "anchor", r) for r in range(5)]
    records = [
        _synthetic_record(units[0], settings, p=0.01),
        _synthetic_record(units[1], settings, p=0.5, overrides={"predictive": {
            **_synthetic_record(units[1], settings)["arms"]["predictive"], "p_values": {"delta": 0.01, "angle": 0.5}
        }}),
        _synthetic_record(units[2], settings, status="failed"),
        _synthetic_record(units[3], settings, status="inconsistent", p=0.01),
        _synthetic_record(units[4], settings, overrides={"forced_k_minus_1": {"status": "failed"}}),
    ]
    rows = {(r["arm"], r["statistic"]): r for r in psc.rate_table(records, settings)}
    assert set(next(iter(rows.values()))) == set(psc.RATE_COLUMNS)
    assert rows[("pls_cv", "delta")]["count"] == 3 and rows[("pls_cv", "delta")]["rejections"] == 1
    assert rows[("pls_cv", "delta")]["failures"] == 1 and rows[("pls_cv", "delta")]["inconsistent"] == 1
    assert rows[("predictive", "delta")]["rejections"] == 2
    assert rows[("forced_k_minus_1", "delta")]["count"] == 2 and rows[("forced_k_minus_1", "delta")]["failures"] == 2
    assert rows[("pls_cv", "shape")]["count"] == 0 and rows[("pls_cv", "shape")]["rate"] is None

    paired = {(r["arm"], r["reference_arm"], r["statistic"]): r for r in psc.paired_table(records, settings)}
    assert set(next(iter(paired.values()))) == set(psc.PAIRED_COLUMNS)
    row = paired[("predictive", "pls_cv", "delta")]
    assert (row["n_paired"], row["b"], row["c"]) == (3, 0, 1)
    assert row["p_gain"] == pytest.approx(0.5) and row["p_loss"] == 1.0
    assert paired[("forced_k_minus_1", "pls_cv", "delta")]["n_paired"] == 2
    assert ("forced_k_minus_1", "predictive", "delta") in paired
    assert ("predictive_rank_plus_one", "predictive", "delta") in paired
    none = paired[("predictive_orthonormal", "pls_cv", "angle")]
    assert (none["b"], none["c"], none["p_loss"], none["p_gain"], none["p_two_sided"]) == (0, 0, 1.0, 1.0, 1.0)

    cells = {r["arm"]: r for r in psc.cell_table(records, settings)}
    assert set(cells["pls_cv"]) == set(psc.CELL_COLUMNS)
    assert cells["pls_cv"]["n"] == 3 and cells["predictive"]["measured_dimension_median"] == 2.0
    assert cells["pls_cv"]["observed_angle_median"] == 10.0 and cells["pls_cv"]["population_angle_median"] == 12.0
    assert cells["pls_cv"]["stage_r2_median"] == 0.9 and cells["pls_cv"]["angle_null_q95_median"] == 30.0

    summary = psc.accounting(records, settings, expected=6, duplicates=1)
    assert summary["failed_units"] == 1 and summary["failed_units_by_cell"] == {"seaad/anchor": 1}
    assert summary["inconsistent_units"] == 1 and summary["inconsistent_units_by_cell"] == {"seaad/anchor": 1}
    assert summary["inconsistent_unit_checks"][0]["checks"] == {"fitted_labels": 0.3}
    assert summary["arm_failures_by_cell"] == {"forced_k_minus_1": {"seaad/anchor": 1}}
    assert summary["missing_units"] == 1 and summary["duplicate_records"] == 1


def test_merge_keeps_one_record_per_unit_and_counts_duplicates(tmp_path) -> None:
    settings = _fast(psc.simulation_settings(replicates=1, permutations=9, designs=["seaad"]))
    directory = tmp_path / "simulation"
    lsc.ensure_settings(directory, settings)
    units = psc.simulation_units(settings)
    shard = lsc.shard_path(directory, 0, 2)
    shard.parent.mkdir(parents=True)
    lines = [lsc._line(_synthetic_record(unit, settings)) for unit in units]
    shard.write_text("".join(lines) + lines[0] + "{torn")
    summary = psc.merge(tmp_path, "simulation")
    assert summary["complete_units"] == len(units) and summary["duplicate_records"] == 1
    assert len((directory / "merged.jsonl").read_text().splitlines()) == len(units)
    other = _synthetic_record(units[0], {**settings, "replicates": 9})
    shard.write_text(lsc._line(other))
    with pytest.raises(SystemExit, match="written under other settings"):
        psc.merge(tmp_path, "simulation")


# ── 4.2 Reproduction ─────────────────────────────────────────────────────────


def _counts(row) -> dict:
    return {"count": int(row["count"]), "rejections": int(row["rejections"])}


def _rates_from(committed: dict, *, bump=None) -> list[dict]:
    rows = []
    for (design, cell, statistic), values in committed.items():
        rejections = values["rejections"] + (1 if bump == (design, cell, statistic) else 0)
        rows.append({"design": design, "cell": cell, "arm": "pls_cv", "statistic": statistic,
                     "count": values["count"], "rejections": rejections})
    return rows


def test_weak_size_reproduction_both_outcomes(full_weak) -> None:
    committed_rows = lsc._read_csv(psc.SPAN_SPLIT_DIR / "production_table.csv")
    committed = {
        (row["design"], row["cell"], row["statistic"]): _counts(row)
        for row in committed_rows if row["design"] == "seaad_weak" and row["cell"] in ("anchor", "translation")
    }
    assert len(committed) == 6 and all(values["count"] == 500 for values in committed.values())
    passed = psc.weak_size_reproduction(_rates_from(committed), full_weak, committed_rows)
    check = passed["checks"]["span_split_production"]
    assert passed["passed"] and check["rows_compared"] == 6 and check["rows_matching"] == 6
    failed = psc.weak_size_reproduction(
        _rates_from(committed, bump=("seaad_weak", "translation", "shape")), full_weak, committed_rows
    )
    differing = failed["checks"]["span_split_production"]["differing"]
    assert not failed["passed"] and len(differing) == 1
    assert differing[0]["key"] == ["seaad_weak", "translation", "shape"]
    assert differing[0]["reproduced"]["rejections"] == differing[0]["committed"]["rejections"] + 1


def _anchor_records(n: int, *, passed: bool = True) -> list[dict]:
    return [
        {"design": "plasmode", "cell": "anchor", "replicate": r, "status": "completed",
         "arms": {"pls_cv": {"status": "completed", "fidelity": {"passed": passed, "rank": 5,
                                                                   "cv_mean_auroc_rounded": 0.635556}}}}
        for r in range(n)
    ]


def test_plasmode_reproduction_both_outcomes(stand_in) -> None:
    settings = _plasmode_settings(stand_in)
    summary_rows = lsc._read_csv(psc.PLASMODE_DIR / "plasmode_summary.csv")
    production_rows = lsc._read_csv(psc.PRODUCTION_PLASMODE_DIR / "rate_table.csv")
    anchor = psc.plasmode_summary_rejections(summary_rows)
    assert anchor == {
        ("plasmode", "anchor", "delta"): {"count": 500, "rejections": 21},
        ("plasmode", "anchor", "angle"): {"count": 500, "rejections": 24},
        ("plasmode", "anchor", "shape"): {"count": 500, "rejections": 13},
    }
    plant = {
        (row["design"], row["cell"], row["statistic"]): _counts(row)
        for row in production_rows if row["arm"] == "pls_cv" and row["cell"] != "anchor"
    }
    assert len(plant) == 9
    rates = _rates_from({**anchor, **plant})
    passed = psc.plasmode_reproduction(_anchor_records(500), rates, settings, summary_rows, production_rows)
    assert passed["passed"]
    assert passed["checks"]["plasmode_anchor"]["rows_compared"] == 3
    assert passed["checks"]["production_plant"]["rows_compared"] == 9
    assert passed["checks"]["fidelity"]["units_passed"] == 500
    bumped = _rates_from({**anchor, **plant}, bump=("plasmode", "span_e0.5", "delta"))
    failed = psc.plasmode_reproduction(_anchor_records(500), bumped, settings, summary_rows, production_rows)
    assert not failed["passed"] and not failed["checks"]["production_plant"]["passed"]
    assert failed["checks"]["plasmode_anchor"]["passed"] and failed["checks"]["fidelity"]["passed"]
    infidel = psc.plasmode_reproduction(_anchor_records(3, passed=False), rates, settings, summary_rows,
                                        production_rows)
    assert not infidel["checks"]["fidelity"]["passed"] and len(infidel["checks"]["fidelity"]["differing"]) == 3
    short = psc.plasmode_reproduction(_anchor_records(2), _rates_from({**anchor, **plant},
                                      bump=("plasmode", "anchor", "angle")), settings, summary_rows, production_rows)
    assert not short["checks"]["plasmode_anchor"]["passed"]


# ── 5.1 / 5.2 Decision on synthetic tables ───────────────────────────────────


def _profile_cells(profile: str) -> list[tuple[str, str]]:
    if profile == "simulation":
        return [(design, cell) for design in psc.SIMULATION_DESIGNS for cell in psc.design_cells(design)]
    if profile == "weak_size":
        return [("seaad_weak", cell) for cell in psc.CONTROL_CELLS]
    return [("plasmode", "anchor"), *(("plasmode", cell) for cell in psc.PLANT_CELLS)]


def _default_count(profile: str, cell: str) -> int:
    if profile == "weak_size" or (profile == "plasmode" and cell == "anchor"):
        return 500
    return 100


def synthetic_tables(*, rate=None, pairs=None) -> dict:
    """Tables for every reading. ``rate(profile, design, cell, arm, statistic)`` gives (count, rejections) or None;
    ``pairs(profile, design, cell, statistic, arm, reference)`` gives (n, b, c) or None (the defaults)."""

    tables = {}
    for profile in psc.PROFILES:
        rows, paired, cells = [], [], []
        for design, cell in _profile_cells(profile):
            mode = psc._cell_mode(cell)
            fields = {"design": design, "cell": cell, "trajectory_mode": mode, "effect_size": 0.0}
            for arm in psc.PROFILE_ARMS[profile]:
                cells.append({**fields, "arm": arm, "n": 100, "rank_median": 3.0, "rank_q10": 2.0, "rank_q90": 4.0,
                              "cv_mean_auroc_median": 0.99, "pooled_stage_median": 0.9,
                              "group_difference_median": 0.5, "stage_r2_median": 0.95,
                              "observed_delta_median": 0.4, "observed_angle_median": 20.0,
                              "population_delta_median": 0.5, "population_angle_median": 22.0})
                for statistic in psc.STATISTICS:
                    count = _default_count(profile, cell)
                    value = rate(profile, design, cell, arm, statistic) if rate else None
                    count, hits = value if value is not None else (count, int(0.04 * count))
                    rows.append({**fields, "arm": arm, "statistic": statistic, "count": count, "rejections": hits,
                                 "rate": hits / count if count else None})
            for statistic in psc.STATISTICS:
                for arm, reference in psc.PAIRS:
                    if arm not in psc.PROFILE_ARMS[profile]:
                        continue
                    value = pairs(profile, design, cell, statistic, arm, reference) if pairs else None
                    n, b, c = value if value is not None else (100, 0, 0)
                    paired.append({**fields, "statistic": statistic, "arm": arm, "reference_arm": reference,
                                   "n_paired": n, "rate_reference": 0.5, "rate_arm": 0.5 + (c - b) / max(n, 1),
                                   "difference": (c - b) / max(n, 1), "b": b, "c": c,
                                   **plsc.paired_p_values(b, c)})
        tables[profile] = {"rates": rows, "paired": paired, "cells": cells}
    return tables


GAIN_CELL = ("simulation", "n200", "orient_e0.25_phi0", "angle")
PASSING_REPRODUCTION = {profile: {"passed": True, "checks": {}} for profile in psc.PROFILES}
NO_INCONSISTENT = {profile: {"inconsistent_units": 0} for profile in psc.PROFILES}


def _gains(arms=psc.CANDIDATE_ARMS, cell=GAIN_CELL):
    def pairs(profile, design, cell_name, statistic, arm, reference):
        if (profile, design, cell_name, statistic) == cell and arm in arms and reference == "pls_cv":
            return (100, 0, 20)
        return None

    return pairs


def _decide(tables, reproduction=PASSING_REPRODUCTION, summaries=NO_INCONSISTENT):
    return psc.decision(tables, reproduction=reproduction, accounting_summaries=summaries)


def test_plan_sizes_and_allowed_counts() -> None:
    plan = psc.criteria_plan()
    sizes = {name: len(readings) for name, readings in plan.items()}
    assert sizes == {"S1_simulation": 42, "S1_weak_size": 6, "S1_plasmode": 3, "S2": 57, "S3_strong": 36,
                     "S3_band": 72, "G": 18, "G_band": 36, "G_weak": 9, "G_plasmode": 3}
    assert psc.allowed_count(57) == int(binom.ppf(0.95, 57, 0.025)) == 4
    assert psc.allowed_count(72) == int(binom.ppf(0.95, 72, 0.025)) == 4


def test_every_arm_qualifies_and_the_nearest_is_named() -> None:
    result = _decide(synthetic_tables(pairs=_gains()))
    assert result["verdict"] == "propose_paper_grade" and result["rule_version"] == 1
    assert result["qualifying_arms"] == list(psc.CANDIDATE_ARMS)
    assert result["named_arm"] == "predictive" and "preferred as nearer production" in result["preference_reason"]
    assert result["arms"]["predictive"]["criteria"]["G"]["gains"][0]["cell"] == "orient_e0.25_phi0"


def test_the_arm_with_fewer_off_target_readings_above_the_bound_is_named() -> None:
    def rate(profile, design, cell, arm, statistic):
        if arm == "predictive" and (profile, design, cell, statistic) == ("simulation", "n200", "shape_e0.1", "delta"):
            return (100, 10)  # above the bound, within the count rule
        return None

    result = _decide(synthetic_tables(rate=rate, pairs=_gains()))
    assert result["arms"]["predictive"]["verdict"] == "qualifies"
    assert result["arms"]["predictive"]["criteria"]["S3"]["off_target_above_bound"] == 1
    assert result["named_arm"] == "predictive_orthonormal"


def test_size_excess_names_the_readings_and_takes_precedence() -> None:
    def rate(profile, design, cell, arm, statistic):
        if arm == "predictive" and (profile, cell, statistic) == ("weak_size", "anchor", "shape"):
            return (500, 60)
        return None

    def pairs(profile, design, cell, statistic, arm, reference):
        if arm == "predictive" and reference == "pls_cv" and (design, cell) == ("seaad", "magnitude_e0.25"):
            return (100, 30, 0)
        return _gains()(profile, design, cell, statistic, arm, reference)

    result = _decide(synthetic_tables(rate=rate, pairs=pairs))
    arm = result["arms"]["predictive"]
    assert arm["verdict"] == "size_excess" and arm["failing_criteria"] == ["S1", "S2"]
    assert arm["failing_readings"]["S1"][0]["family"] == "weak_size"
    assert arm["failing_readings"]["S1"][0]["statistic"] == "shape"
    assert arm["failing_readings"]["S2"][0]["cell"] == "magnitude_e0.25"
    assert result["named_arm"] == "predictive_orthonormal"


def test_power_loss() -> None:
    def pairs(profile, design, cell, statistic, arm, reference):
        if arm == "predictive_fitted" and reference == "pls_cv" and (profile, cell) == ("plasmode", "span_e0.5"):
            return (100, 25, 0)
        return _gains()(profile, design, cell, statistic, arm, reference)

    arm = _decide(synthetic_tables(pairs=pairs))["arms"]["predictive_fitted"]
    assert arm["verdict"] == "power_loss"
    assert arm["failing_readings"]["S2"][0]["status"] == "below_ceiling"
    assert arm["failing_readings"]["S2"][0]["profile"] == "plasmode"


@pytest.mark.parametrize("family", ["strong", "band"])
def test_specificity_fails_in_each_family_separately(family) -> None:
    def rate(profile, design, cell, arm, statistic):
        if family == "strong" and arm == "predictive" and (design, cell, statistic) == ("seaad", "orient_e0.5_phi0",
                                                                                         "shape"):
            return (100, 30)
        return None

    def pairs(profile, design, cell, statistic, arm, reference):
        if family == "band" and arm == "predictive" and reference == "pls_cv" and (design, cell, statistic) == (
            "n200_delta0.25", "orient_e0.5_phi0", "shape"
        ):
            return (100, 0, 30)
        return _gains()(profile, design, cell, statistic, arm, reference)

    result = _decide(synthetic_tables(rate=rate, pairs=pairs))
    s3 = result["arms"]["predictive"]["criteria"]["S3"]
    assert result["arms"]["predictive"]["verdict"] == "specificity_excess"
    assert s3["families"][family]["holds"] is False
    assert s3["families"]["band" if family == "strong" else "strong"]["holds"] is True
    readings = result["arms"]["predictive"]["failing_readings"]["S3"]
    assert readings[0]["family"] == family and readings[0]["cell"] == "orient_e0.5_phi0"
    if family == "band":
        assert readings[0]["production_rate"] == pytest.approx(0.04)


def test_a_gain_only_outside_delta_two_fails_g() -> None:
    band_gain = _gains(cell=("simulation", "n200_delta0.25", "orient_e0.25_phi0", "angle"))
    result = _decide(synthetic_tables(pairs=band_gain))
    arm = result["arms"]["predictive"]
    assert arm["verdict"] == "no_gain" and result["verdict"] == "keep_production" and result["named_arm"] is None
    assert arm["recorded_gains"]["band"]["holds"] is True
    assert arm["recorded_gains"]["weak"]["holds"] is False


def test_an_unavailable_reading_fails_its_criterion() -> None:
    def rate(profile, design, cell, arm, statistic):
        if arm == "predictive" and (profile, cell, statistic) == ("plasmode", "anchor", "angle"):
            return (0, 0)
        return None

    def pairs(profile, design, cell, statistic, arm, reference):
        if arm == "predictive" and reference == "pls_cv" and (design, cell, statistic) == ("n200", "shape_e0.5",
                                                                                            "shape"):
            return (0, 0, 0)
        return _gains()(profile, design, cell, statistic, arm, reference)

    arm = _decide(synthetic_tables(rate=rate, pairs=pairs))["arms"]["predictive"]
    assert arm["failing_criteria"] == ["S1", "S2", "G"]
    assert arm["criteria"]["S1"]["families"]["plasmode"]["n_unavailable"] == 1
    assert arm["criteria"]["G"]["n_unavailable"] == 1


def test_the_gate_still_records_every_criterion() -> None:
    reproduction = {**PASSING_REPRODUCTION, "plasmode": {"passed": False, "checks": {"fidelity": {"passed": False}}}}
    result = _decide(synthetic_tables(pairs=_gains()), reproduction=reproduction)
    assert result["verdict"] == "reproduction_failed" and result["named_arm"] is None
    assert result["gate"]["failures"] == [{"profile": "plasmode", "reason": "fidelity failed"}]
    assert all(set(arm["criteria"]) == {"S1", "S2", "S3", "G"} for arm in result["arms"].values())
    assert result["arms"]["predictive"]["verdict"] == "qualifies"
    inconsistent = _decide(synthetic_tables(pairs=_gains()),
                           summaries={**NO_INCONSISTENT, "simulation": {"inconsistent_units": 2}})
    assert inconsistent["verdict"] == "reproduction_failed"
    assert inconsistent["gate"]["failures"] == [{"profile": "simulation", "reason": "2 inconsistent unit(s)"}]


def test_descriptive_tables_have_one_row_per_cell_and_arm() -> None:
    descriptive = _decide(synthetic_tables(pairs=_gains()))["descriptive"]
    assert len(descriptive["control_arm"]) == 7 * 9 + 3
    assert {"rate_pls_cv", "rate_predictive", "rate_forced_k_minus_1", "against_pls_cv",
            "against_predictive"} <= set(descriptive["control_arm"][0])
    assert len(descriptive["rank_plus_one"]) == 7 * 9 * 3
    assert len(descriptive["rank_sensitivity"]) == 6 * 3
    assert all(row["statistic"] == "shape" for row in descriptive["rank_sensitivity"])
    assert len(descriptive["surgery_cell"]) == len(psc.ARMS) * 3
    assert len(descriptive["weak_targets"]) == 9
    assert len(descriptive["band_off_target_rates"]) == 72
    assert len(descriptive["orientation_estimation"]) == 7 * 3 * len(psc.ARMS) + len(psc.ARMS)
    designs = {(row["profile"], row["design"]) for row in descriptive["rank_and_auroc"]}
    assert designs == {
        *(("simulation", design) for design in psc.SIMULATION_DESIGNS), ("weak_size", "seaad_weak"),
        ("plasmode", "plasmode"),
    }
    assert all(set(row) >= {"design", "cell", "arm", "observed_angle_median", "population_angle_median"}
               for row in descriptive["orientation_estimation"])
