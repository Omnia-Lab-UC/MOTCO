"""``scripts/predictive_size_excess_check.py``: settings, units, arms, the two checks, merge, reproduction, decision.

Every unit runs on the plasmode tests' synthetic stand-in export, with stand-in
copies of the parents built from it. No real data are read.
"""

from __future__ import annotations

import dataclasses
import json
import shutil
import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest
from scipy.stats import binom


def _script(name: str):
    sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
    try:
        return __import__(name)
    finally:
        sys.path.pop(0)


pse = _script("predictive_size_excess_check")
psc, plsc, pc, lsc, grade = pse.psc, pse.plsc, pse.pc, pse.lsc, pse.grade

REPO = Path(__file__).resolve().parents[1]
PERMUTATIONS = 9


def _plasmode_tests():
    sys.path.insert(0, str(REPO))
    try:
        from tests import test_plasmode_check as module
    finally:
        sys.path.pop(0)
    return module


plasmode_tests = _plasmode_tests()


def _write_settings(directory: Path, settings: dict) -> None:
    directory.mkdir(parents=True, exist_ok=True)
    (directory / "settings.json").write_text(json.dumps({"hash": lsc.settings_hash(settings), "settings": settings}))


# ── Stand-in export and parents ──────────────────────────────────────────────


@pytest.fixture(scope="module")
def stand_in(tmp_path_factory):
    """The stand-in export, its plasmode and production-selection parents, and a stand-in pilot (500/100)."""

    directory = plasmode_tests.write_export(tmp_path_factory.mktemp("export"))
    export = pc.load_export(directory)
    base = {**pc.base_settings(export), "cv": dict(plasmode_tests.SMALL_CV), "attribution": None}
    built = pc.build_unit(0.0, 0, {**base, "stream": "fidelity"}, export)
    _, result = pc.analyse(built["blocks"], built["metadata"], export, pc.analysis_params(base, None, permutations=0))
    latent = result.record["latent"]
    fidelity = {"selected_lv": int(latent["selected_lv"]), "cv_mean_auroc": round(float(latent["cv_mean_auroc"]), 6)}
    parent = pc.build_settings(export, fidelity, replicates=2, permutations=PERMUTATIONS, bootstrap=5, e_star=3.0,
                               descriptive=2.0)
    parent["cv"] = dict(plasmode_tests.SMALL_CV)
    parent_dir = tmp_path_factory.mktemp("plasmode")
    _write_settings(parent_dir, parent)
    shutil.copy(psc.PLASMODE_DIR / "plasmode_summary.csv", parent_dir / "plasmode_summary.csv")
    production = plsc.plasmode_settings(export, parent_dir=parent_dir, replicates=2, permutations=PERMUTATIONS)
    production_dir = tmp_path_factory.mktemp("production")
    _write_settings(production_dir, production)
    shutil.copy(psc.PRODUCTION_PLASMODE_DIR / "rate_table.csv", production_dir / "rate_table.csv")
    pilot = psc.plasmode_settings(export, parent_dir=parent_dir, production_dir=production_dir, replicates=500,
                                  plant_replicates=100, permutations=PERMUTATIONS)
    pilot_dir = tmp_path_factory.mktemp("pilot")
    _write_settings(pilot_dir, pilot)
    shutil.copy(pse.PILOT_PLASMODE_DIR / "rate_table.csv", pilot_dir / "rate_table.csv")
    return {"directory": directory, "export": export, "parent_dir": parent_dir, "production_dir": production_dir,
            "pilot_dir": pilot_dir, "pilot": pilot}


def _settings(stand_in, **kwargs) -> dict:
    kwargs.setdefault("replicates", 2)
    kwargs.setdefault("plant_replicates", 1)
    kwargs.setdefault("permutations", PERMUTATIONS)
    kwargs.setdefault("pilot_dir", stand_in["pilot_dir"])
    return pse.build_settings(stand_in["export"], parent_dir=stand_in["parent_dir"],
                              production_dir=stand_in["production_dir"], **kwargs)


@pytest.fixture(scope="module")
def settings(stand_in):
    return _settings(stand_in)


@pytest.fixture(scope="module")
def anchor_record(stand_in, settings):
    return pse.evaluate_unit(lsc.Unit("plasmode", "anchor", 1), settings, stand_in["export"])


@pytest.fixture(scope="module")
def plant_record(stand_in, settings):
    return pse.evaluate_unit(lsc.Unit("plasmode", "span_e0.5", 0), settings, stand_in["export"])


def _preflighted(tmp_path: Path, stand_in, settings: dict) -> Path:
    out = tmp_path / "run"
    result = pse.preflight(out, settings, stand_in["export"])
    assert result["passed"], result
    return out


# ── 1.1 Settings, parents, resume ────────────────────────────────────────────


def test_settings_bind_the_pilots_files_and_its_parents(settings) -> None:
    assert set(settings["parents"]) == {
        "plasmode_settings", "plasmode_summary", "production_plasmode_settings", "production_plasmode_rates",
        "pilot_plasmode_settings", "pilot_plasmode_rates",
    }
    assert settings["arms"] == list(pse.ARMS) and len(pse.ARMS) == 14
    assert settings["decisive_arms"] == list(pse.DECISIVE_ARMS)
    assert settings["anchor_rank"] == settings["pilot"]["fidelity"]["selected_lv"]
    assert settings["rejection"] == "p <= alpha" and settings["residual_correction"] == "df"


def test_the_default_parents_are_the_committed_pilots() -> None:
    committed = json.loads((pse.PILOT_PLASMODE_DIR / "settings.json").read_text())["settings"]
    pse.check_pilot(committed, committed, source=pse.PILOT_PLASMODE_DIR)
    assert committed["replicates"] == {"anchor": 500, "plant": 100} and committed["permutations"] == 999
    assert committed["plasmode"]["base_seed"] == 2034 and committed["plasmode"]["stream"] == "paper"
    assert len(psc.plasmode_units(committed)) == 800
    assert (pse.PILOT_PLASMODE_DIR / "rate_table.csv").exists()


def test_other_units_than_the_pilots_are_refused(stand_in) -> None:
    with pytest.raises(SystemExit, match="replicates.anchor"):
        _settings(stand_in, replicates=501)
    with pytest.raises(SystemExit, match="permutations"):
        _settings(stand_in, permutations=PERMUTATIONS + 1)


def test_resume_skips_complete_units_and_other_settings_are_refused(tmp_path, stand_in, settings) -> None:
    export = stand_in["export"]
    out = _preflighted(tmp_path, stand_in, settings)
    shard = pse.enumerate_units(settings)[0::2]
    assert len(shard) == 3
    assert pse.run_shard(out, settings, 0, 2, export=export, limit=1) == 1
    assert pse.run_shard(out, settings, 0, 2, export=export) == 2
    assert pse.run_shard(out, settings, 0, 2, export=export) == 0
    records = pse.read_records([lsc.shard_path(out, 0, 2)])
    assert sorted(pse.record_key(r) for r in records) == sorted(pse.unit_key(*unit.key) for unit in shard)
    other = _settings(stand_in, replicates=1)
    with pytest.raises(SystemExit, match="pilot.replicates.anchor"):
        pse.run_shard(out, other, 0, 2, export=export)


def test_the_run_waits_on_a_passed_preflight(tmp_path, stand_in, settings) -> None:
    with pytest.raises(SystemExit, match="run `preflight` first"):
        pse.run_shard(tmp_path / "bare", settings, 0, 1, export=stand_in["export"], limit=0)
    out = _preflighted(tmp_path, stand_in, settings)
    record = json.loads((out / pse.PREFLIGHT_NAME).read_text())
    record["passed"] = False
    (out / pse.PREFLIGHT_NAME).write_text(json.dumps(record))
    with pytest.raises(SystemExit, match="pre-flight failed"):
        pse.run_shard(out, settings, 0, 1, export=stand_in["export"], limit=0)


def test_the_preflight_compares_the_cv_with_the_forced_fit(tmp_path, stand_in, settings) -> None:
    out = _preflighted(tmp_path, stand_in, settings)
    record = json.loads((out / pse.PREFLIGHT_NAME).read_text())
    assert [unit["unit"] for unit in record["units"]] == [pse.unit_key("plasmode", "anchor", r) for r in (0, 1)]
    for unit in record["units"]:
        assert unit["fidelity"]["passed"] and unit["forced_rank"] == settings["anchor_rank"]
        assert unit["statistics_differing"] == [] and unit["p_values_differing"] == []
        assert unit["cross_validated"] == unit["forced"]
    assert pse.load_preflight(out, settings)["passed"] is True
    assert pse.load_preflight(out, _settings(stand_in, replicates=1)) is None


def test_a_changed_parent_is_refused_by_name(tmp_path, stand_in) -> None:
    pilot_dir = tmp_path / "pilot"
    shutil.copytree(stand_in["pilot_dir"], pilot_dir)
    settings = _settings(stand_in, pilot_dir=pilot_dir)
    out = _preflighted(tmp_path, stand_in, settings)
    pse.run_shard(out, settings, 0, 1, export=stand_in["export"], limit=1)
    pse.merge(out)
    with (pilot_dir / "rate_table.csv").open("a") as handle:
        handle.write("\n")
    for call in (
        lambda: pse.run_shard(out, settings, 0, 1, export=stand_in["export"], limit=0),
        lambda: pse.merge(out),
        lambda: pse.decide(out),
        lambda: pse.preflight(out, settings, stand_in["export"]),
    ):
        with pytest.raises(SystemExit, match="pilot_plasmode_rates"):
            call()


# ── 1.2 Units and seeds ──────────────────────────────────────────────────────


def test_units_and_seeds_are_the_pilots(stand_in) -> None:
    full = _settings(stand_in, replicates=500, plant_replicates=100)
    units = pse.enumerate_units(full)
    assert len(units) == 800
    assert units == psc.plasmode_units(stand_in["pilot"])
    assert sum(unit.cell == "anchor" for unit in units) == 500
    assert {unit.cell for unit in units} == {"anchor", "span_e0.25", "span_e0.5", "span_e1"}
    parent = full["pilot"]["plasmode"]
    assert parent == stand_in["pilot"]["plasmode"]
    export = stand_in["export"]
    for unit in (units[0], units[1], units[5], units[-1]):
        ours = plsc.build_plasmode_unit(unit, full["pilot"], export)
        theirs = plsc.build_plasmode_unit(unit, stand_in["pilot"], export)
        assert ours["seeds"] == theirs["seeds"] == pc.unit_seeds(parent["base_seed"], parent["stream"], unit.replicate)
        assert ours["metadata"].equals(theirs["metadata"])
        for name in pc.BLOCKS:
            assert ours["blocks"][name].equals(theirs["blocks"][name])


def test_a_changed_export_is_refused_before_any_unit(tmp_path, stand_in, settings) -> None:
    out = _preflighted(tmp_path, stand_in, settings)
    altered = dataclasses.replace(stand_in["export"], files={**stand_in["export"].files, "rna.csv": "0" * 64})
    with pytest.raises(SystemExit, match="rna.csv"):
        pse.run_shard(out, settings, 0, 1, export=altered)
    assert not (out / "records").exists()
    with pytest.raises(SystemExit, match="rna.csv"):
        _settings({**stand_in, "export": altered})


# ── 1.3 Reference arms ───────────────────────────────────────────────────────


def _counting_cv(monkeypatch) -> list:
    import motco.simulations.evaluation as evaluation_module

    calls: list = []
    original = evaluation_module.plsda_doubleCV

    def counting(*args, **kwargs):
        calls.append(1)
        return original(*args, **kwargs)

    monkeypatch.setattr(evaluation_module, "plsda_doubleCV", counting)
    return calls


def test_an_anchor_unit_makes_no_cross_validation(stand_in, settings, monkeypatch) -> None:
    calls = _counting_cv(monkeypatch)
    record = pse.evaluate_unit(lsc.Unit("plasmode", "anchor", 0), settings, stand_in["export"])
    assert record["status"] == "completed" and calls == []
    assert record["rank"] == settings["anchor_rank"] and record["component_selection"] == "forced"
    for arm, entry in record["arms"].items():
        assert entry["rank"] == settings["anchor_rank"], arm


def test_a_plant_unit_makes_exactly_one_cross_validation(stand_in, settings, monkeypatch) -> None:
    calls = _counting_cv(monkeypatch)
    record = pse.evaluate_unit(lsc.Unit("plasmode", "span_e1", 0), settings, stand_in["export"])
    assert record["status"] == "completed" and len(calls) == 1
    assert record["component_selection"] == "cv"
    assert {entry["rank"] for entry in record["arms"].values()} == {record["rank"]}


@pytest.mark.parametrize("key", [("plasmode", "anchor", 1), ("plasmode", "span_e0.5", 0)])
def test_the_reference_arms_equal_the_pilots_records(stand_in, settings, anchor_record, plant_record, key) -> None:
    ours = anchor_record if key[1] == "anchor" else plant_record
    theirs = psc.evaluate_plasmode_unit(lsc.Unit(*key), settings["pilot"], stand_in["export"])
    assert theirs["status"] == "completed"
    for arm in pse.REFERENCE_ARMS:
        assert ours["arms"][arm]["pair_statistics"] == theirs["arms"][arm]["pair_statistics"], arm
        assert ours["arms"][arm]["p_values"] == theirs["arms"][arm]["p_values"], arm
        assert ours["arms"][arm]["rank"] == theirs["arms"][arm]["rank"], arm
    assert ours["seeds"] == theirs["seeds"] and ours["rrpp_seed"] == theirs["rrpp_seed"]


def test_a_record_holds_seeds_and_statistics_only(stand_in, anchor_record, plant_record) -> None:
    export = stand_in["export"]
    n = len(export.metadata)

    def no_vectors(value) -> bool:
        if isinstance(value, dict):
            return all(no_vectors(item) for item in value.values())
        if isinstance(value, list):
            return len(value) != n and all(no_vectors(item) for item in value)
        return True

    for record in (anchor_record, plant_record):
        text = json.dumps(record)
        for donor in export.metadata[export.sample_col]:
            assert donor not in text
        assert "Female" not in text and "Male" not in text
        assert no_vectors(record)
        assert record["seeds"] == pc.unit_seeds(2034, "paper", record["replicate"])
        assert record["fold_seed"] == pse.fold_seed(record["replicate"]) == record["folds"]["seed"]
        assert set(record["folds"]) == {"seed", "count", "training_stage_counts"}


def test_a_failed_unit_is_recorded_and_the_shard_goes_on(stand_in, settings, monkeypatch) -> None:
    def broken(*args, **kwargs):
        raise RuntimeError("boom")

    monkeypatch.setattr(pse.plsc, "build_plasmode_unit", broken)
    record = pse.evaluate_unit(lsc.Unit("plasmode", "anchor", 0), settings, stand_in["export"])
    assert record["status"] == "failed" and record["error_message"] == "boom"


# ── 2.1 The inference path and its identity ──────────────────────────────────


@pytest.mark.parametrize("which", ["anchor", "plant"])
def test_path_identity_holds_on_both_cell_kinds(anchor_record, plant_record, which) -> None:
    record = anchor_record if which == "anchor" else plant_record
    assert record["status"] == "completed", record["inconsistent_checks"]
    identity = record["consistency"]["path_identity"]
    assert identity["passed"] and identity["deviation"] <= pse.PATH_IDENTITY_TOLERANCE
    assert identity["p_values_equal"] and set(identity["arms"]) == set(pse.REFERENCE_ARMS)
    for check in pse.CONSISTENCY_CHECKS:
        assert record["consistency"][check]["passed"], check


def test_a_perturbed_latent_matrix_is_flagged_and_leaves_every_rate(stand_in, settings, anchor_record,
                                                                     monkeypatch) -> None:
    original = pse.reference_latent

    def perturbed(result):
        latent = original(result)
        if result.record["measurement"].get("latent_coordinates", {}).get("value") == "predictive":
            # Row-varying: a constant shift would be absorbed by the intercept and change no statistic.
            return latent + 1e-3 * np.random.default_rng(0).normal(size=latent.shape)
        return latent

    monkeypatch.setattr(pse, "reference_latent", perturbed)
    record = pse.evaluate_unit(lsc.Unit("plasmode", "anchor", 0), settings, stand_in["export"])
    assert record["status"] == "inconsistent"
    assert "path_identity" in record["inconsistent_checks"]
    assert record["consistency"]["path_identity"]["failing_arms"] == ["predictive"]
    summary = pse.accounting([record, anchor_record], settings, expected=2, duplicates=0)
    assert summary["inconsistent_units"] == 1
    assert "predictive" in summary["inconsistent_unit_checks"][0]["checks"]["path_identity"]["arms"]
    rates = pse.rate_table([record, anchor_record], settings)
    assert {row["count"] for row in rates} == {1} and {row["inconsistent"] for row in rates} == {1}


# ── 2.2 Within-stage draws ───────────────────────────────────────────────────


def _design(n_per_cell: int = 6, seed: int = 3):
    from motco.simulations.evaluation import build_simulation_trajectory_design

    rng = np.random.default_rng(seed)
    working = pd.DataFrame({
        "group": np.repeat(["0", "1"], 3 * n_per_cell),
        "stage": np.tile(np.repeat(["0", "1", "2"], n_per_cell), 2),
    }).sample(frac=1.0, random_state=seed).reset_index(drop=True)
    design = build_simulation_trajectory_design(working, group_col="group", stage_col="stage")
    Y = rng.normal(size=(len(working), 4)) * np.where(working["stage"].to_numpy() == "2", 3.0, 1.0)[:, None]
    return working, design, Y


def test_every_drawn_index_stays_in_its_rows_stage() -> None:
    rng = np.random.default_rng(0)
    strata = rng.choice(["0", "1", "2"], size=50)
    draws = pse.within_stage_indices(strata, ["0", "1", "2"], 300, seed=11)
    assert draws.shape == (300, 50)
    assert np.all(strata[draws] == strata[None, :])
    for row in draws:
        assert sorted(row.tolist()) == list(range(50))
    assert len({row.tobytes() for row in draws}) > 250
    with pytest.raises(ValueError, match="exactly one"):
        pse.within_stage_indices(strata, ["0", "1"], 3, seed=1)


@pytest.mark.parametrize("correction", [None, "df"])
def test_one_stratum_is_rrpp_exactly(correction) -> None:
    from motco.stats.permutation import RRPP

    working, design, Y = _design()
    indices = pse.within_stage_indices(np.full(len(working), "all"), ["all"], 25, seed=97)
    ours = pse.restricted_null(Y, design.model_full, design.model_reduced, design.ls_means, design.contrast, indices,
                               orientation_estimator="net_displacement", residual_correction=correction)
    theirs = RRPP(Y, design.model_full, design.model_reduced, design.ls_means, design.contrast, permutations=25,
                  n_jobs=1, progress=False, seed=97, orientation_estimator="net_displacement",
                  residual_correction=correction)
    for left, right in zip(ours, theirs):
        assert len(left) == len(right) == 25
        for a, b in zip(left, right):
            assert np.array_equal(np.asarray(a), np.asarray(b), equal_nan=True)


def test_within_stage_residuals_sum_to_zero_per_stage_and_keep_the_df_scale() -> None:
    working, design, Y = _design()
    stages = working["stage"].to_numpy()
    _, residuals = pse.reduced_fit(Y, design.model_reduced, None)
    _, corrected = pse.reduced_fit(Y, design.model_reduced, "df")
    n, rank = design.model_reduced.shape[0], np.linalg.matrix_rank(np.asarray(design.model_reduced, dtype=float))
    np.testing.assert_allclose(corrected, residuals * np.sqrt(n / (n - rank)), rtol=1e-12)
    for idx in pse.within_stage_indices(stages, ["0", "1", "2"], 20, seed=5):
        permuted = corrected[idx]
        for stage in ("0", "1", "2"):
            np.testing.assert_allclose(permuted[stages == stage].sum(axis=0), 0.0, atol=1e-10)


def test_the_within_stage_arms_of_a_unit_share_one_sequence(stand_in, settings, monkeypatch) -> None:
    seen: list[np.ndarray] = []
    original = pse.restricted_null

    def recording(*args, **kwargs):
        seen.append(np.array(args[5], copy=True))
        return original(*args, **kwargs)

    monkeypatch.setattr(pse, "restricted_null", recording)
    record = pse.evaluate_unit(lsc.Unit("plasmode", "anchor", 0), settings, stand_in["export"])
    assert record["status"] == "completed"
    assert len(seen) == 5  # four reference coordinate systems and the out-of-fold coordinates
    assert all(np.array_equal(seen[0], other) for other in seen[1:])
    assert seen[0].shape == (PERMUTATIONS, len(stand_in["export"].metadata))


# ── 2.3 Out-of-fold coordinates ──────────────────────────────────────────────


@pytest.fixture(scope="module")
def unit_inputs(stand_in, settings):
    """The anchor unit's features, stage codes and rank, and its in-sample ``predictive_fitted`` matrix."""

    export = stand_in["export"]
    unit = lsc.Unit("plasmode", "anchor", 0)
    built = plsc.build_plasmode_unit(unit, settings["pilot"], export)
    production = psc.plasmode_production_params(settings["pilot"], built["seeds"])
    params = dataclasses.replace(production, forced_components=settings["anchor_rank"],
                                 latent_coordinates="predictive_fitted")
    dataset, result = pc.analyse(built["blocks"], built["metadata"], export, params)
    inputs = pse.analysis._trajectory_inputs(dataset, production)
    return {
        "features": result.latent_artifacts.features.to_numpy(dtype=float),
        "strata": inputs.working["stage"].astype(str).to_numpy(),
        "rank": settings["anchor_rank"],
        "fitted": result.latent_scores.to_numpy(dtype=float),
        "scale": bool(result.latent_artifacts.model.scale),
    }


def test_every_donor_is_held_out_once_and_never_trains_its_own_row(unit_inputs, monkeypatch) -> None:
    strata = unit_inputs["strata"]
    assignment = pse.stratified_folds(strata, ["0", "1", "2"], pse.fold_seed(0))
    assert sorted(set(assignment.tolist())) == [0, 1, 2, 3]
    for fold in range(4):
        training = strata[assignment != fold]
        assert sorted(set(training.tolist())) == ["0", "1", "2"]
    trained: list[np.ndarray] = []
    original = pse.fit_plsda_model

    def recording(X, y, **kwargs):
        trained.append(np.array(X, copy=True))
        return original(X, y, **kwargs)

    monkeypatch.setattr(pse, "fit_plsda_model", recording)
    features = unit_inputs["features"]
    built = pse.out_of_fold_fitted(features, strata, assignment, rank=unit_inputs["rank"], scale=True)
    assert len(trained) == 4
    for fold, X_train in enumerate(trained):
        np.testing.assert_array_equal(X_train, features[assignment != fold])
        held = {row.tobytes() for row in features[assignment == fold]}
        assert not held & {row.tobytes() for row in X_train}
    assert np.isfinite(built["matrix"]).all() and built["matrix"].shape == (len(strata), 3)
    assert built["training_stage_counts"] == [[15, 15, 30]] * 4


def test_the_all_donor_limit_is_predictive_fitted(unit_inputs) -> None:
    rows = np.ones(len(unit_inputs["strata"]), dtype=bool)
    matrix, _ = pse.fitted_indicators(unit_inputs["features"], unit_inputs["strata"], rows, rows,
                                      rank=unit_inputs["rank"], scale=unit_inputs["scale"])
    np.testing.assert_allclose(matrix, unit_inputs["fitted"], rtol=0.0, atol=1e-8)


def test_folds_are_reproducible_per_replicate(unit_inputs) -> None:
    strata, order = unit_inputs["strata"], ["0", "1", "2"]
    first = pse.stratified_folds(strata, order, pse.fold_seed(3))
    again = pse.stratified_folds(strata, order, pse.fold_seed(3))
    other = pse.stratified_folds(strata, order, pse.fold_seed(4))
    assert np.array_equal(first, again) and not np.array_equal(first, other)
    kwargs = {"rank": unit_inputs["rank"], "scale": True}
    left = pse.out_of_fold_fitted(unit_inputs["features"], strata, first, **kwargs)["matrix"]
    right = pse.out_of_fold_fitted(unit_inputs["features"], strata, again, **kwargs)["matrix"]
    assert np.array_equal(left, right)
    assert pse.fold_seed(3) != pse.fold_seed(4)
    assert pse.fold_seed(3) == pse._seed(2037, 3, "folds")


# ── 2.4 The fourteen arms ────────────────────────────────────────────────────


def test_the_record_holds_every_arm_and_reading(anchor_record) -> None:
    record = anchor_record
    assert list(record["arms"]) == list(pse.ARMS)
    for arm, entry in record["arms"].items():
        assert entry["status"] == "completed", (arm, entry.get("error_message"))
        assert set(entry["pair_statistics"]) == set(pse.STATISTICS)
        assert set(entry["p_values"]) == set(pse.STATISTICS)
        assert entry["angle_null_q95"] is not None
        assert {key: entry[key] for key in ("coordinates", "draws", "covariates", "role")} == pse.ARM_SPECS[arm]
    assert set(record["coordinate_systems"]) == set(pse.COORDINATE_SYSTEMS)
    for system in record["coordinate_systems"].values():
        assert set(system["dispersion"]) == {"stage_0", "stage_1", "stage_2", "max_min_ratio"}
    assert set(record["stage_fit"]) == {"in_sample_r2", "out_of_fold_share"}
    assert record["stage_fit"]["out_of_fold_share"] < record["stage_fit"]["in_sample_r2"]
    assert set(record["consistency"]) == set(pse.CONSISTENCY_CHECKS)


def test_the_covariate_free_arms_measure_the_reference_latent_matrices(anchor_record, plant_record) -> None:
    for record in (anchor_record, plant_record):
        assert record["consistency"]["same_fit_covariate_free"]["passed"]
        for arm in pse.COVARIATE_FREE_ARMS:
            assert record["arms"][arm]["latent_deviation"] == 0.0, arm


def test_a_failing_arm_leaves_the_other_thirteen(stand_in, settings, monkeypatch) -> None:
    original = pse.analyse_without_covariates

    def failing(blocks, metadata, export, params):
        if params.latent_coordinates == "predictive":
            raise RuntimeError("arm boom")
        return original(blocks, metadata, export, params)

    monkeypatch.setattr(pse, "analyse_without_covariates", failing)
    record = pse.evaluate_unit(lsc.Unit("plasmode", "anchor", 0), settings, stand_in["export"])
    failed = [arm for arm, entry in record["arms"].items() if entry["status"] != "completed"]
    assert failed == ["predictive__no_covariates"]
    assert record["arms"]["predictive__no_covariates"]["error_message"] == "arm boom"
    assert record["status"] == "completed" and len(record["arms"]) == 14


# ── 2.5 sbatch ───────────────────────────────────────────────────────────────


def test_sbatch_template_parses_and_dry_runs(tmp_path) -> None:
    script = REPO / "scripts" / "predictive_size_excess_array.sbatch"
    shell = shutil.which("bash") or "bash"
    subprocess.run([shell, "-n", str(script)], check=True)
    env = {"PATH": "/usr/bin:/bin", "PSE_DRY_RUN": "1", "PSE_OUT": str(tmp_path / "run"),
           "PSE_EXPORT": str(tmp_path / "export"), "PSE_REPO": str(REPO), "N_SHARDS": "100",
           "SLURM_ARRAY_TASK_ID": "3", "PSE_ARGS": "--limit 1"}
    out = subprocess.run([shell, str(script)], env=env, cwd=tmp_path, check=True, capture_output=True,
                         text=True).stdout
    words = out.split()
    assert words[0] == f"PYTHONPATH={REPO / 'src'}"
    assert words[1:] == ["python", str(REPO / "scripts" / "predictive_size_excess_check.py"), "run", "--export-dir",
                         str(tmp_path / "export"), "--out-dir", str(tmp_path / "run"), "--shard", "3/100", "--limit",
                         "1"]
    assert not (tmp_path / "run").exists()
    relative = subprocess.run([shell, str(script)], env={**env, "PSE_OUT": "run"}, cwd=tmp_path,
                              capture_output=True, text=True)
    assert relative.returncode == 2 and "absolute" in relative.stderr


# ── 3.1 Merge on synthetic records ───────────────────────────────────────────


def _synthetic(cell: str, replicate: int, p_values: dict, *, status: str = "completed") -> dict:
    arms = {}
    for arm in pse.ARMS:
        p = p_values.get(arm, 0.5)
        arms[arm] = {"status": "completed", **pse.ARM_SPECS[arm], "rank": 5,
                     "pair_statistics": dict.fromkeys(pse.STATISTICS, 1.0),
                     "p_values": dict.fromkeys(pse.STATISTICS, p) if not isinstance(p, dict) else p,
                     "angle_null_q95": 30.0 + replicate, "seconds": 1.0}
    return {
        "design": "plasmode", "cell": cell, "trajectory_mode": "none" if cell == "anchor" else "magnitude",
        "effect_size": 0.0 if cell == "anchor" else 0.5, "replicate": replicate, "status": status,
        "inconsistent_checks": [], "consistency": {}, "rank": 5,
        "stage_fit": {"in_sample_r2": 0.84, "out_of_fold_share": 0.2 + 0.01 * replicate},
        "coordinate_systems": {
            system: {"dispersion": {"stage_0": 0.5, "stage_1": 1.0, "stage_2": 1.5, "max_min_ratio": 3.0}}
            for system in pse.COORDINATE_SYSTEMS
        },
        "arms": arms, "runtime_seconds": {"unit": 10.0, "production": 5.0},
    }


def test_the_counterparts_of_the_ten_non_reference_arms() -> None:
    counterparts = {(arm, ref) for arm, ref, relation in pse.pairs() if relation == "counterpart"}
    assert counterparts == {
        ("predictive__within_stage", "predictive"),
        ("predictive_orthonormal__within_stage", "predictive_orthonormal"),
        ("predictive_fitted__within_stage", "predictive_fitted"),
        ("pls_cv__within_stage", "pls_cv"),
        ("fitted_out_of_fold", "predictive_fitted"),
        ("fitted_out_of_fold__within_stage", "predictive_fitted__within_stage"),
        ("pls_cv__no_covariates", "pls_cv"),
        ("predictive__no_covariates", "predictive"),
        ("predictive_orthonormal__no_covariates", "predictive_orthonormal"),
        ("predictive_fitted__no_covariates", "predictive_fitted"),
    }
    production = {arm for arm, ref, relation in pse.pairs() if relation == "production"}
    assert production == set(pse.ARMS) - {"pls_cv", "pls_cv__within_stage", "pls_cv__no_covariates"}
    assert all(ref == "pls_cv" for _, ref, relation in pse.pairs() if relation == "production")


def test_merge_tables_on_synthetic_records(settings) -> None:
    records = [
        _synthetic("anchor", 0, {"predictive": 0.01, "predictive__within_stage": 0.5}),
        _synthetic("anchor", 1, {"predictive": 0.01, "predictive__within_stage": 0.01}),
        _synthetic("anchor", 2, {"predictive": 0.5, "predictive__within_stage": 0.04}),
        _synthetic("anchor", 3, {}, status="inconsistent"),
        {"design": "plasmode", "cell": "anchor", "trajectory_mode": "none", "effect_size": 0.0, "replicate": 4,
         "status": "failed"},
    ]
    rates = {(r["arm"], r["statistic"]): r for r in pse.rate_table(records, settings)}
    assert rates[("predictive", "delta")]["count"] == 3 and rates[("predictive", "delta")]["rejections"] == 2
    assert rates[("predictive__within_stage", "delta")]["rejections"] == 2
    assert rates[("pls_cv", "delta")]["failures"] == 1 and rates[("pls_cv", "delta")]["inconsistent"] == 1
    paired = {(r["arm"], r["reference_arm"], r["statistic"]): r for r in pse.paired_table(records, settings)}
    row = paired[("predictive__within_stage", "predictive", "delta")]
    assert (row["n_paired"], row["b"], row["c"], row["relation"]) == (3, 1, 1, "counterpart")
    assert row["p_fewer"] == pytest.approx(float(binom.sf(0, 2, 0.5)))
    against = paired[("predictive", "pls_cv", "delta")]
    assert (against["b"], against["c"], against["relation"]) == (0, 2, "production")
    assert against["p_more"] == pytest.approx(0.25) and against["p_fewer"] == 1.0
    quiet = paired[("pls_cv__no_covariates", "pls_cv", "angle")]
    assert (quiet["b"], quiet["c"]) == (0, 0)
    assert quiet["p_fewer"] == quiet["p_more"] == quiet["p_two_sided"] == 1.0
    cells = {r["arm"]: r for r in pse.cell_table(records, settings)}
    assert cells["fitted_out_of_fold"]["coordinates"] == "out_of_fold" and cells["fitted_out_of_fold"]["n"] == 3
    assert cells["pls_cv"]["dispersion_stage_2_median"] == 1.5 and cells["pls_cv"]["dispersion_ratio_median"] == 3.0
    assert cells["pls_cv"]["stage_r2_median"] == 0.84
    assert cells["pls_cv"]["out_of_fold_share_median"] == pytest.approx(0.21)
    assert cells["pls_cv"]["angle_null_q95_median"] == 31.0 and cells["pls_cv"]["rank_median"] == 5


def test_merge_writes_every_table_and_counts_duplicates(tmp_path, stand_in, settings) -> None:
    out = _preflighted(tmp_path, stand_in, settings)
    pse.run_shard(out, settings, 0, 1, export=stand_in["export"], limit=2)
    path = lsc.shard_path(out, 0, 1)
    lines = path.read_text().splitlines()
    path.write_text("\n".join([*lines, lines[0]]) + "\n")
    summary = pse.merge(out)
    assert summary["complete_units"] == 2 and summary["duplicate_records"] == 1
    assert summary["missing_units"] == len(pse.enumerate_units(settings)) - 2
    for name in ("rate_table.csv", "paired_table.csv", "cell_table.csv", "accounting.json", "reproduction.json",
                 "merged.jsonl"):
        assert (out / name).exists(), name
    accounting = json.loads((out / "accounting.json").read_text())
    assert accounting["settings_hash"] == lsc.settings_hash(settings)
    assert set(accounting["arm_seconds_median"]) == set(pse.ARMS)


# ── 3.2 Reproduction ─────────────────────────────────────────────────────────


def _committed_rows() -> list[dict]:
    return lsc._read_csv(pse.PILOT_PLASMODE_DIR / "rate_table.csv")


def test_the_committed_table_resolves_all_48_rows() -> None:
    counts = pse._committed_counts(_committed_rows())
    assert len(counts) == 48
    assert counts[("anchor", "predictive", "delta")] == {"count": 500, "rejections": 68}
    assert counts[("anchor", "pls_cv", "shape")] == {"count": 500, "rejections": 13}


def test_reproduction_passes_on_the_committed_counts_and_lists_every_difference(settings) -> None:
    committed = _committed_rows()
    passed = pse.reproduction(committed, settings, committed)
    assert passed["passed"] and passed["rows_compared"] == 48 and passed["rows_matching"] == 48
    bumped = [dict(row) for row in committed]
    for row in bumped:
        if (row["cell"], row["arm"], row["statistic"]) == ("span_e0.5", "predictive_fitted", "angle"):
            row["rejections"] = str(int(row["rejections"]) + 1)
    failed = pse.reproduction(bumped, settings, committed)
    assert not failed["passed"] and failed["rows_matching"] == 47
    [row] = failed["differing"]
    assert (row["cell"], row["arm"], row["statistic"]) == ("span_e0.5", "predictive_fitted", "angle")
    assert row["reproduced"]["rejections"] == row["committed"]["rejections"] + 1
    missing = pse.reproduction([r for r in committed if r["arm"] != "pls_cv"], settings, committed)
    assert not missing["passed"] and len(missing["differing"]) == 12


# ── 4.1 / 4.2 Decision on synthetic tables ───────────────────────────────────

RESTORED = (25, 25, 25)
EXCESS = (66, 10, 38)
REDUCING = (30, 2)  # b, c against the counterpart: fewer rejections, p well below 0.025
NEUTRAL = (5, 5)
PASSING = {"passed": True, "differing": []}
CLEAN = {"inconsistent_units": 0}
PREFLIGHT = {"passed": True}


def synthetic_tables(grades: dict[str, str], *, descriptive_noise: bool = False) -> dict:
    """Rate and paired rows realising a grade for each decisive arm; every other reading is neutral."""

    rates, paired, cells = [], [], []
    for cell in ("anchor", *pse.PLANT_CELLS):
        n = 500 if cell == "anchor" else 100
        for arm in pse.ARMS:
            counts = dict(zip(pse.STATISTICS, RESTORED if cell == "anchor" else (50, 5, 5)))
            if cell == "anchor" and grades.get(arm, "restored") != "restored":
                counts = dict(zip(pse.STATISTICS, EXCESS))
            if descriptive_noise and arm not in pse.DECISIVE_ARMS and arm not in pse.REFERENCE_ARMS:
                counts = {statistic: n - 1 for statistic in pse.STATISTICS}
            if descriptive_noise and cell != "anchor":
                counts = {statistic: 3 for statistic in pse.STATISTICS}
            for statistic in pse.STATISTICS:
                rates.append({"design": "plasmode", "cell": cell, "arm": arm, "statistic": statistic,
                              "count": str(n), "rejections": str(counts[statistic]),
                              "rate": str(counts[statistic] / n)})
            cells.append({"cell": cell, "arm": arm, **pse.ARM_SPECS[arm], "angle_null_q95_median": "30",
                          "dispersion_stage_0_median": "0.5", "dispersion_ratio_median": "3"})
        for arm, reference, _ in pse.pairs():
            for statistic in pse.STATISTICS:
                b, c = NEUTRAL
                if cell == "anchor" and grades.get(arm) == "reduced" and reference == pse.COUNTERPARTS.get(arm):
                    b, c = REDUCING if statistic == "delta" else NEUTRAL
                if descriptive_noise and arm not in pse.DECISIVE_ARMS:
                    b, c = (0, 40)
                paired.append({"cell": cell, "statistic": statistic, "arm": arm, "reference_arm": reference,
                               "n_paired": "500", "b": str(b), "c": str(c), "rate_reference": "0.1",
                               "rate_arm": "0.1"})
    return {"rates": rates, "paired": paired, "cells": cells}


def _decide(grades: dict[str, str], *, reproduced=PASSING, summary=CLEAN, checked=PREFLIGHT, noise=False) -> dict:
    return pse.decision(synthetic_tables(grades, descriptive_noise=noise), reproduced=reproduced, summary=summary,
                        checked=checked)


W_ARMS = pse.WITHIN_STAGE_ARMS
O_ARM = pse.OUT_OF_FOLD_ARM


def test_the_synthetic_counts_sit_on_the_intended_side_of_the_size_rule() -> None:
    def judged(counts):
        return grade.judge_rates([{"n_available": 500, "n_rejected": k} for k in counts])["holds"]

    assert judged(RESTORED) and not judged(EXCESS)
    assert plsc.paired_p_values(*REDUCING)["p_loss"] <= 0.025 < plsc.paired_p_values(*NEUTRAL)["p_loss"]


@pytest.mark.parametrize(
    ("grades", "verdict"),
    [
        ({}, "both_restore"),
        ({O_ARM: "unchanged"}, "null_restores"),
        ({**dict.fromkeys(W_ARMS, "unchanged")}, "out_of_fold_restores"),
        ({**dict.fromkeys(W_ARMS, "unchanged"), O_ARM: "unchanged"}, "unexplained"),
        ({**dict.fromkeys(W_ARMS, "reduced"), O_ARM: "unchanged"}, "partly_explained"),
        ({**dict.fromkeys(W_ARMS, "unchanged"), O_ARM: "reduced"}, "partly_explained"),
    ],
)
def test_each_verdict(grades, verdict) -> None:
    result = _decide(grades)
    assert result["verdict"] == verdict
    assert result["consequence"] == pse.CONSEQUENCES[verdict]
    for arm in pse.DECISIVE_ARMS:
        assert result["arms"][arm]["grade"] == grades.get(arm, "restored"), arm
    if verdict == "partly_explained":
        assert result["reduced_checks"] == (["W"] if grades.get(O_ARM) == "unchanged" else ["O"])


def test_the_reduced_grade_names_the_statistic_and_the_counts() -> None:
    result = _decide({O_ARM: "reduced", **dict.fromkeys(W_ARMS, "unchanged")})
    graded = result["arms"][O_ARM]
    assert graded["grade"] == "reduced" and graded["reduced_statistics"] == ["delta"]
    assert graded["counterpart"] == "predictive_fitted"
    assert (graded["paired"]["delta"]["b"], graded["paired"]["delta"]["c"]) == REDUCING
    assert "b = 30, c = 2" in graded["reason"]


def test_the_restored_grade_records_counts_bound_and_ceiling() -> None:
    graded = _decide({})["arms"][W_ARMS[0]]
    rows = graded["size_rule"]["observations"]
    assert [row["n_rejected"] for row in rows] == list(RESTORED)
    assert all("bound" in row and "ceiling" in row for row in rows)


def test_one_metric_is_not_enough() -> None:
    result = _decide({W_ARMS[0]: "restored", W_ARMS[1]: "restored", W_ARMS[2]: "unchanged", O_ARM: "unchanged"})
    assert result["checks"]["W"]["grade"] == "unchanged" and result["verdict"] == "unexplained"
    reduced = _decide({W_ARMS[2]: "reduced", O_ARM: "unchanged"})
    assert reduced["checks"]["W"]["grade"] == "reduced" and reduced["verdict"] == "partly_explained"
    assert pse.weakest(["restored", "reduced", "restored"]) == "reduced"
    assert pse.weakest(["restored", "unchanged", "reduced"]) == "unchanged"


def test_an_unavailable_reading_is_unchanged_and_named() -> None:
    tables = synthetic_tables({})
    tables["rates"] = [row for row in tables["rates"]
                       if not (row["arm"] == O_ARM and row["cell"] == "anchor" and row["statistic"] == "angle")]
    result = pse.decision(tables, reproduced=PASSING, summary=CLEAN, checked=PREFLIGHT)
    assert result["arms"][O_ARM]["grade"] == "unchanged" and result["arms"][O_ARM]["unavailable"] == ["angle"]
    assert "angle" in result["arms"][O_ARM]["reason"] and result["verdict"] == "null_restores"


@pytest.mark.parametrize(
    ("kwargs", "reason"),
    [
        ({"reproduced": {"passed": False, "differing": [{}]}}, "reproduction failed"),
        ({"summary": {"inconsistent_units": 2}}, "2 inconsistent unit(s)"),
        ({"checked": None}, "pre-flight record is missing"),
        ({"checked": {"passed": False}}, "pre-flight did not pass"),
    ],
)
def test_the_gate_still_records_the_grades(kwargs, reason) -> None:
    result = _decide({O_ARM: "unchanged"}, **kwargs)
    assert result["verdict"] == "reproduction_failed" and not result["gate"]["passed"]
    assert any(reason in failure for failure in result["gate"]["failures"])
    assert result["checks"]["W"]["grade"] == "restored" and result["checks"]["O"]["grade"] == "unchanged"
    assert result["consequence"] == pse.CONSEQUENCES["reproduction_failed"]


def test_the_consequences_are_the_designs() -> None:
    assert set(pse.CONSEQUENCES) == set(pse.VERDICTS)
    assert "within-stage draws as an opt-in option" in pse.CONSEQUENCES["null_restores"]
    assert "out-of-fold reading confirms" in pse.CONSEQUENCES["both_restore"]
    assert "cross-fitted measurement" in pse.CONSEQUENCES["out_of_fold_restores"]
    assert "open item 2's cohorts" in pse.CONSEQUENCES["partly_explained"]
    assert "prediction for ADNI is withdrawn" in pse.CONSEQUENCES["unexplained"]


@pytest.mark.parametrize(
    "grades", [{}, {O_ARM: "unchanged"}, {**dict.fromkeys(W_ARMS, "reduced"), O_ARM: "unchanged"}]
)
def test_descriptive_readings_never_change_the_verdict(grades) -> None:
    plain, noisy = _decide(grades), _decide(grades, noise=True)
    assert noisy["verdict"] == plain["verdict"]
    assert {arm: entry["grade"] for arm, entry in noisy["arms"].items()} == {
        arm: entry["grade"] for arm, entry in plain["arms"].items()
    }
    assert noisy["descriptive"] != plain["descriptive"]


def test_a_restored_arm_states_its_plant_delta_beside_production() -> None:
    result = _decide({O_ARM: "unchanged"})
    stated = {row["arm"]: row for row in result["descriptive"]["restored_arms_plant_delta"]}
    assert set(stated) == set(W_ARMS)
    row = stated[W_ARMS[0]]
    assert row["statement_cell"] == "span_e0.5" and set(row["delta_rates"]) == set(pse.PLANT_CELLS)
    assert row["delta_rates"]["span_e0.5"] == {"rate_arm": 0.5, "rate_pls_cv": 0.5, "rate_counterpart": 0.5}
    assert row["delta_above_production"] is False
    descriptive = result["descriptive"]
    assert [entry["arm"] for entry in descriptive["anchor_descriptive_arms"]] == list(pse.DESCRIPTIVE_ARMS[:2])
    assert len(descriptive["covariate_free"]) == 4 * 4 * 3
    assert {(row["cell"], row["coordinates"]) for row in descriptive["dispersion"]} == {
        (cell, system) for cell in ("anchor", *pse.PLANT_CELLS) for system in pse.COORDINATE_SYSTEMS
    }
    assert len(descriptive["null_widths"]) == 4 * 14


def test_decide_writes_the_decision_file(tmp_path, stand_in, settings) -> None:
    out = _preflighted(tmp_path, stand_in, settings)
    pse.run_shard(out, settings, 0, 1, export=stand_in["export"])
    pse.merge(out)
    result = pse.decide(out)
    written = json.loads((out / pse.DECISION_NAME).read_text())
    assert written["verdict"] == result["verdict"] == "reproduction_failed"  # 2 units cannot reproduce 500
    assert written["rule_version"] == 1 and written["preflight"]["passed"] is True
    assert written["complete"]["complete_units"] == len(pse.enumerate_units(settings))
    assert set(written["inputs"]) == {"settings.json", "rate_table.csv", "paired_table.csv", "cell_table.csv",
                                      "accounting.json", "reproduction.json", "preflight.json"}
    assert set(written["arms"]) == set(pse.DECISIVE_ARMS)
