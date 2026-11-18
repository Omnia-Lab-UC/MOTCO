"""``scripts/adni_latent_plasmode_check.py``: units, parts, the estimator, spaces, checks, the calibration,
pre-flights, merge, reproductions and the decision.

Units run on the anchor test's synthetic stand-in export (no real data are
read) at small CV knobs and few splits. The estimator's properties are checked
on Gaussian configurations around a known truth; the merge and the decision on
hand-built records and tables.
"""

from __future__ import annotations

import copy
import dataclasses
import json
import math
import re
import shlex
import shutil
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


alp = _script("adni_latent_plasmode_check")
aac, gain, lsc, plsc = alp.aac, alp.gain, alp.lsc, alp.plsc

REPO = Path(__file__).resolve().parents[1]
STAND_IN_POINTS = (40, 88)
FIXED = (2, 4)
SPLITS = 4
N_STAGES = 3


# ── Fixtures ─────────────────────────────────────────────────────────────────


@pytest.fixture(scope="module")
def export_dir(tmp_path_factory) -> Path:
    return write_export(tmp_path_factory.mktemp("export") / "motco_input")


@pytest.fixture(scope="module")
def contract(export_dir) -> aac.CohortContract:
    """The stand-in's own cohort fit (label-free), as ADNI's 8 and 0.599514."""

    cohort = aac.load_cohort(export_dir, _stand_in_contract(fidelity_rank=2, fidelity_auroc=0.5))
    settings = alp.cohort_settings(cohort, replicates=2, splits=SPLITS, point_sizes=STAND_IN_POINTS,
                                   fixed_ranks=FIXED, cv=SMALL_CV)
    built = aac.build_unit(cohort, settings, aac.COHORT_POINT, 0)
    dataset = aac.make_dataset(built["blocks"], built["metadata"], cohort.contract)
    latent = alp.fit_analysis(dataset, settings, 1, None).record["latent"]
    return _stand_in_contract(fidelity_rank=int(latent["selected_lv"]),
                              fidelity_auroc=round(float(latent["cv_mean_auroc"]), 6))


@pytest.fixture(scope="module")
def cohort(export_dir, contract):
    return aac.load_cohort(export_dir, contract)


def _cohort_settings(cohort, **changes) -> dict:
    kwargs = {"replicates": 2, "splits": SPLITS, "point_sizes": STAND_IN_POINTS, "fixed_ranks": FIXED, "cv": SMALL_CV}
    kwargs.update(changes)
    return alp.cohort_settings(cohort, **kwargs)


def _calibration_settings(cohort, **changes) -> dict:
    kwargs = {"cohorts": 2, "units": 2, "splits": SPLITS, "point_sizes": STAND_IN_POINTS, "fixed_ranks": FIXED,
              "cv": SMALL_CV}
    kwargs.update(changes)
    return alp.calibration_settings(cohort, **kwargs)


@pytest.fixture(scope="module")
def settings(cohort) -> dict:
    return _cohort_settings(cohort)


@pytest.fixture(scope="module")
def cal_settings(cohort) -> dict:
    return _calibration_settings(cohort)


@pytest.fixture(scope="module")
def records(cohort, settings) -> dict:
    return {point: alp.evaluate_cohort_unit(alp.Unit(point, 0), settings, cohort) for point in ("n40", "n88")}


@pytest.fixture(scope="module")
def cal_records(cohort, cal_settings) -> dict:
    units = [alp.Unit("n40", 0, 1.0, 0), alp.Unit("n88", 1, 0.0, 1), alp.Unit("n88", 0, 0.67, 1)]
    return {alp.unit_key(unit): alp.evaluate_calibration_unit(unit, cal_settings, cohort) for unit in units}


@pytest.fixture(scope="module")
def prepared(tmp_path_factory, cohort, settings, cal_settings) -> Path:
    """An output root with both pre-flights passed."""

    out_dir = tmp_path_factory.mktemp("prepared")
    assert alp.preflight(out_dir, settings, cohort)["passed"]
    assert alp.preflight(out_dir, cal_settings, cohort)["passed"]
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


def _fits(cohort, settings, point="n40", replicate=0):
    built = aac.build_unit(cohort, settings, point, replicate)
    dataset = aac.make_dataset(built["blocks"], built["metadata"], cohort.contract)
    production = alp.fit_analysis(dataset, settings, 1, None)
    forced = {rank: alp.fit_analysis(dataset, settings, 1, rank) for rank in settings["fixed_ranks"]}
    return built, production, forced


def _orthonormal(rng, p: int, d: int) -> np.ndarray:
    return np.linalg.qr(rng.normal(size=(p, d)))[0]


# ── 1.1 Settings, resume and refusals ───────────────────────────────────────


def test_settings_bind_the_parents_and_the_export(settings, cal_settings) -> None:
    for current in (settings, cal_settings):
        assert set(current["parents"]) == {"anchor_points", "gain_cells"}
        assert current["parents"]["anchor_points"]["path"] == \
            "results/adni-real-noise-anchor-2026-11-17/point_table.csv"
        assert current["parents"]["gain_cells"]["path"] == \
            "results/adni-gain-specificity-2026-11-17/plants/cell_table.csv"
        assert sorted(current["export"]["files"]) == ["covariates.csv", "manifest.json", "metadata.csv", "nmr.csv",
                                                      "p180.csv"]
        assert current["own_base_seed"] == alp.BASE_SEED == 2041
        assert current["permutations"] == 0 and current["analysis"]["permutations"] == 0
        assert list(current["spaces"]) == ["pls_cv", "predictive", "pls_rank2", "pls_rank4", "pca"]
    assert settings["profile"] == "cohort" and settings["base_seed"] == aac.BASE_SEED == 2039
    assert cal_settings["profile"] == "calibration" and "base_seed" not in cal_settings
    assert cal_settings["strengths"] == [0.0, 0.67, 1.0]


def test_resume_skips_complete_units(tmp_path, cohort, settings, prepared) -> None:
    out_dir = tmp_path / "out"
    shutil.copytree(prepared, out_dir)
    assert alp.run_shard(out_dir, settings, cohort, 0, 2, limit=1) == 1
    assert alp.run_shard(out_dir, settings, cohort, 0, 2) == 1
    assert alp.run_shard(out_dir, settings, cohort, 0, 2) == 0
    path = lsc.shard_path(out_dir / "cohort", 0, 2)
    keys = [alp.record_key(json.loads(line)) for line in path.read_text().splitlines()]
    expected = [alp.unit_key(unit) for position, unit in enumerate(alp.enumerate_units(settings)) if position % 2 == 0]
    assert sorted(keys) == sorted(expected) == sorted(["n40/0", "n40/1"])


def test_other_settings_are_refused_naming_the_keys(tmp_path, cohort, prepared) -> None:
    out_dir = tmp_path / "out"
    shutil.copytree(prepared, out_dir)
    other = _cohort_settings(cohort, replicates=3)
    with pytest.raises(SystemExit, match="differing: replicates"):
        alp.run_shard(out_dir, other, cohort, 0, 1)
    with pytest.raises(SystemExit, match="differing: replicates"):
        alp.preflight(out_dir, other, cohort)
    with pytest.raises(SystemExit, match="differing: cohorts"):
        alp.run_shard(out_dir, _calibration_settings(cohort, cohorts=3), cohort, 0, 1)
    with pytest.raises(SystemExit, match="differing: splits"):
        alp.run_shard(out_dir, _calibration_settings(cohort, splits=SPLITS + 1), cohort, 0, 1)


@pytest.mark.parametrize("role", sorted(alp.PARENT_FILES))
def test_a_changed_parent_is_refused_by_name(tmp_path, cohort, export_dir, role) -> None:
    anchor, gain_dir = tmp_path / "anchor", tmp_path / "gain"
    shutil.copytree(alp.ANCHOR_DIR, anchor)
    shutil.copytree(alp.GAIN_DIR, gain_dir)
    out_dir = tmp_path / "out"
    current = {profile: build(cohort, anchor_dir=anchor, gain_dir=gain_dir) for profile, build in
               (("cohort", _cohort_settings), ("calibration", _calibration_settings))}
    for profile, value in current.items():
        lsc.ensure_settings(out_dir / profile, value)
    source, name = alp.PARENT_FILES[role]
    with ((anchor if source == "anchor" else gain_dir) / name).open("a") as handle:
        handle.write("\n")
    calls = (lambda: alp.run_shard(out_dir, current["cohort"], cohort, 0, 1),
             lambda: alp.run_shard(out_dir, current["calibration"], cohort, 0, 1),
             lambda: alp.preflight(out_dir, current["cohort"], cohort),
             lambda: alp.merge(out_dir, "cohort", export_dir),
             lambda: alp.merge(out_dir, "calibration", export_dir),
             lambda: alp.decide(out_dir, export_dir))
    for call in calls:
        with pytest.raises(SystemExit, match=f"Parent {role} .* changed"):
            call()


def test_another_export_is_refused_by_file(tmp_path, export_dir, contract, settings, prepared) -> None:
    out_dir = tmp_path / "out"
    shutil.copytree(prepared, out_dir)
    other = _copy_export(export_dir, tmp_path / "other")
    frame = _block(other, "nmr")
    frame.iloc[0, 0] += 0.5
    _rewrite(other, "nmr.csv", frame, index=True)
    moved = aac.load_cohort(other, contract)
    with pytest.raises(SystemExit, match=r"export\.files\.manifest\.json, export\.files\.nmr\.csv"):
        alp.run_shard(out_dir, _cohort_settings(moved), moved, 0, 1)
    with pytest.raises(SystemExit, match="differ from the settings' for: manifest.json, nmr.csv"):
        alp.run_shard(out_dir, settings, moved, 0, 1)
    for call in (lambda: alp.merge(out_dir, "cohort", other), lambda: alp.merge(out_dir, "calibration", other),
                 lambda: alp.decide(out_dir, other)):
        with pytest.raises(SystemExit, match="differ from the settings' for: manifest.json, nmr.csv"):
            call()


def test_a_contract_violation_is_named_before_any_unit(monkeypatch, tmp_path, export_dir, contract) -> None:
    monkeypatch.setitem(aac.CONTRACTS, "adni", contract)
    broken = _copy_export(export_dir, tmp_path / "broken")
    (broken / "nmr.csv").unlink()
    out_dir = tmp_path / "out"
    for command in ("preflight", "run"):
        for profile in alp.PROFILES:
            with pytest.raises(SystemExit, match="nmr.csv"):
                alp.main([command, "--profile", profile, "--export-dir", str(broken), "--out-dir", str(out_dir)])
    with pytest.raises(SystemExit):  # argparse: --export-dir is required
        alp.main(["run", "--profile", "cohort", "--out-dir", str(out_dir)])
    assert not out_dir.exists()


# ── 1.2 Units ────────────────────────────────────────────────────────────────


def test_a_unit_is_the_anchor_and_gain_drivers_unit(cohort, settings, records) -> None:
    anchor = aac.build_settings(cohort, replicates=2, permutations=9, point_sizes=STAND_IN_POINTS, fixed_ranks=FIXED,
                                cv=SMALL_CV, parent_dir=aac.PARENT_DIR)
    plants = gain.plants_settings(cohort, control_replicates=2, plant_replicates=2, permutations=9,
                                  point_sizes=(STAND_IN_POINTS[1],), cv=SMALL_CV)
    for replicate in (0, 1):
        for point, theirs_settings in (("n40", anchor), ("n88", plants)):
            ours, theirs = (aac.build_unit(cohort, current, point, replicate) for current in (settings,
                                                                                             theirs_settings))
            assert np.array_equal(ours["rows"], theirs["rows"])
            assert ours["metadata"].equals(theirs["metadata"])
    seeds = aac.unit_seeds(anchor, "n40", 0)
    assert {key: records["n40"]["seeds"][key] for key in ("subsample", "labels", "rrpp")} == {
        key: seeds[key] for key in ("subsample", "labels", "rrpp")}


def test_adni_points_and_unit_counts() -> None:
    points = aac.design_points(ADNI_TABLE, alp.POINT_SIZES)
    anchor = json.loads((alp.ANCHOR_DIR / "settings.json").read_text())["settings"]["points"]
    plants = json.loads((alp.GAIN_DIR / "settings.json").read_text())["settings"]["points"]
    committed = {point["name"]: point["table"] for point in (*anchor, *plants)}
    for point in points[:-1]:
        assert point["table"] == committed[point["name"]], point["name"]
    paper = {"profile": "cohort", "points": points, "replicates": alp.REPLICATES}
    assert len(alp.enumerate_units(paper)) == 2000
    assert alp.expected_units(paper) == {"n80": 500, "n200": 500, "n500": 500, "n750": 500}
    calibration = {"profile": "calibration", "points": points, "strengths": list(alp.STRENGTHS),
                   "cohorts": alp.CALIBRATION_COHORTS, "units": alp.CALIBRATION_UNITS}
    units = alp.enumerate_units(calibration)
    assert len(units) == 7200 and len({alp.unit_key(unit) for unit in units}) == 7200
    assert alp.expected_units(calibration) == {name: 1800 for name in ("n80", "n200", "n500", "n750")}


def test_no_committed_driver_uses_the_base_seed(settings) -> None:
    this = Path(alp.__file__).name
    for path in sorted((REPO / "scripts").glob("*.py")):
        if path.name != this:
            assert not re.search(r"\b2041\b", path.read_text(encoding="utf-8")), path.name
    split = alp.split_seeds(settings, alp.Unit("n40", 0))
    assert len(set(split)) == SPLITS
    assert not set(split) & set(aac.unit_seeds(settings, "n40", 0).values())


# ── 1.3 Parts ────────────────────────────────────────────────────────────────


def test_the_parts_partition_the_cohort(cohort, settings) -> None:
    built = aac.build_unit(cohort, settings, "n40", 0)
    remainder = gain.heldout_rows(cohort, built["rows"])
    metadata = cohort.metadata.iloc[remainder].reset_index(drop=True)
    stages = gain.stage_positions(metadata, cohort.contract)
    halves = []
    for seed in alp.split_seeds(settings, alp.Unit("n40", 0)):
        first, second = alp.draw_split(metadata, cohort.contract, seed)
        parts = [set(built["rows"]), set(remainder[first]), set(remainder[second])]
        assert all(not (a & b) for i, a in enumerate(parts) for b in parts[i + 1:])
        assert set().union(*parts) == set(range(len(cohort.metadata)))
        for stage in range(N_STAGES):
            assert abs(int((stages[first] == stage).sum()) - int((stages[second] == stage).sum())) <= 1
        halves.append((remainder[first], remainder[second]))
    assert alp.check_parts(len(cohort.metadata), built["rows"], halves)["passed"]
    assert not np.array_equal(halves[0][0], halves[1][0])


# ── 2.1–2.2 The estimator and the readings ──────────────────────────────────


def test_the_cross_product_is_unbiased_with_the_truth_known() -> None:
    rng = np.random.default_rng(11)
    p, d = 40, 5
    basis = _orthonormal(rng, p, d)
    inside = alp.centred(rng.normal(size=(N_STAGES, d))) @ basis.T
    outside = alp.centred(rng.normal(size=(N_STAGES, p)))
    outside -= outside @ basis @ basis.T
    truth = inside + 0.3 * outside
    share = float(np.sum((truth @ basis) ** 2) / np.sum(truth**2))
    numerator = denominator = 0.0
    naive = []
    for _ in range(4000):
        m1, m2 = (truth + alp.centred(rng.normal(size=(N_STAGES, p))) for _ in range(2))
        products = alp.split_products(m1, m2, {"q": basis})
        numerator += products["spaces"]["q"]["cross"]
        denominator += products["total"]["cross"]
        naive.append(alp.naive_products(m1, {"q": basis}))
    estimate = numerator / denominator
    single = float(np.mean([entry["spaces"]["q"] / entry["total"] for entry in naive]))
    assert abs(estimate - share) < 0.03, (estimate, share)
    assert abs(single - share) > 0.1, (single, share)


def test_equal_halves_recover_the_share_exactly() -> None:
    rng = np.random.default_rng(3)
    basis = _orthonormal(rng, 30, 4)
    m = alp.centred(rng.normal(size=(N_STAGES, 30)))
    products = alp.split_products(m, m, {"q": basis})
    share = float(np.sum((m @ basis) ** 2) / np.sum(m**2))
    assert math.isclose(products["spaces"]["q"]["cross"] / products["total"]["cross"], share, rel_tol=1e-12)
    assert products["total"]["noise"] == 0.0 and products["total"]["correlation"] == pytest.approx(1.0)


def test_the_gram_traces_equal_the_scalars() -> None:
    rng = np.random.default_rng(5)
    bases = {"a": _orthonormal(rng, 25, 3), "b": _orthonormal(rng, 25, 7)}
    products = alp.mean_products([alp.split_products(*(alp.centred(rng.normal(size=(N_STAGES, 25)))
                                                       for _ in range(2)), bases) for _ in range(4)])
    for entry in (products["total"], *products["spaces"].values()):
        assert np.trace(entry["gram"]) == pytest.approx(entry["cross"], abs=1e-12)
    assert alp.check_gram(products)["passed"]


@pytest.mark.parametrize("n_stages", [3, 4])
def test_the_contrast_energies_sum_to_the_total(n_stages) -> None:
    rng = np.random.default_rng(n_stages)
    basis = _orthonormal(rng, 20, 3)
    products = alp.split_products(alp.centred(rng.normal(size=(n_stages, 20))),
                                  alp.centred(rng.normal(size=(n_stages, 20))), {"q": basis})
    projectors = alp.contrast_projectors(n_stages)
    centring = np.eye(n_stages) - 1.0 / n_stages
    assert np.allclose(projectors["net_displacement"] + projectors["interior"], centring)
    assert round(float(np.trace(projectors["interior"]))) == n_stages - 2
    for entry in (products["total"], products["spaces"]["q"]):
        energies = [alp.contrast_energy(entry["gram"], projectors[name]) for name in alp.CONTRASTS]
        assert sum(energies) == pytest.approx(entry["cross"], abs=1e-12)
    record = {"cross_products": {"total": products["total"], "spaces": products["spaces"]}}
    readings = {name: alp.unit_reading(record, "q", name, n_stages) for name in ("total", *alp.CONTRASTS)}
    assert sum(readings[name][0] for name in alp.CONTRASTS) == pytest.approx(readings["total"][0], abs=1e-12)
    assert sum(readings[name][1] for name in alp.CONTRASTS) == pytest.approx(readings["total"][1], abs=1e-12)
    if n_stages == 3:  # the interior is MCI's deviation from the chord
        interior = np.array([1.0, -2.0, 1.0]) / math.sqrt(6.0)
        assert np.allclose(projectors["interior"], np.outer(interior, interior))


def test_the_noise_capture_does_not_depend_on_the_configuration() -> None:
    rng = np.random.default_rng(9)
    basis = _orthonormal(rng, 30, 4)
    e1, e2 = (alp.centred(rng.normal(size=(N_STAGES, 30))) for _ in range(2))
    program = 5.0 * alp.centred(rng.normal(size=(N_STAGES, 30)))
    alone, added = (alp.split_products(e1 + shift, e2 + shift, {"q": basis}) for shift in (0.0, program))
    assert added["total"]["noise"] == pytest.approx(alone["total"]["noise"], rel=1e-12)
    assert added["spaces"]["q"]["noise"] == pytest.approx(alone["spaces"]["q"]["noise"], rel=1e-12)
    assert added["total"]["cross"] != pytest.approx(alone["total"]["cross"])


# ── 2.3 Spaces ───────────────────────────────────────────────────────────────


def test_a_unit_makes_one_cross_validation_and_draws_no_permutation(monkeypatch, cohort, settings) -> None:
    import motco.analysis as analysis_module

    calls = _counting_cv(monkeypatch)

    def refuse(*args, **kwargs):
        raise AssertionError("RRPP was called")

    monkeypatch.setattr(analysis_module, "RRPP", refuse)
    record = alp.evaluate_cohort_unit(alp.Unit("n40", 1), settings, cohort)
    assert record["status"] == "completed", record.get("error_message")
    assert len(calls) == 1
    calibration = alp.evaluate_calibration_unit(alp.Unit("n40", 0, 0.67, 0), _calibration_settings(cohort), cohort)
    assert calibration["status"] == "completed" and len(calls) == 2


def test_the_spaces_bases(cohort, settings) -> None:
    from motco.simulations.diagnostics import LinearLatentMap

    _, production, forced = _fits(cohort, settings)
    model = production.latent_artifacts.model
    rank = int(production.record["latent"]["selected_lv"])
    bases, dimensions, scalings = alp.space_bases(production, forced, settings, N_STAGES)
    assert list(bases) == list(settings["spaces"])
    assert dimensions == {"pls_cv": rank, "predictive": N_STAGES - 1, "pls_rank2": 2, "pls_rank4": 4, "pca": rank}
    for name, basis in bases.items():
        assert basis.shape == (sum(STAND_IN_FEATURES.values()), dimensions[name])
        assert np.allclose(basis.T @ basis, np.eye(basis.shape[1]), atol=1e-10)
    assert np.array_equal(bases["pls_cv"], LinearLatentMap.from_pls(model).basis)
    assert np.array_equal(bases["predictive"], LinearLatentMap.from_pls_predictive(model, aac.PREDICTIVE).basis)
    # The predictive subspace lies inside the fitted span of the same model.
    projected = bases["pls_cv"] @ (bases["pls_cv"].T @ bases["predictive"])
    assert np.allclose(projected, bases["predictive"], atol=1e-10)
    assert alp.check_basis(bases, dimensions, scalings)["passed"]
    assert set(scalings["forced"]) == {"pls_rank2", "pls_rank4"}


def test_the_recorded_dimensions(records) -> None:
    for record in records.values():
        rank = record["production"]["rank"]
        dims = {name: entry["dimension"] for name, entry in record["cross_products"]["spaces"].items()}
        assert dims == {"pls_cv": rank, "predictive": N_STAGES - 1, "pls_rank2": 2, "pls_rank4": 4, "pca": rank}


# ── 2.4 Checks and records ───────────────────────────────────────────────────


def test_every_check_passes_on_stand_in_units(records, cal_records) -> None:
    for record in (*records.values(), *cal_records.values()):
        assert record["status"] == "completed", record.get("error_message") or record.get("inconsistent_checks")
        expected = alp.CHECKS[record["profile"]]
        assert list(record["consistency"]) == list(expected)
        assert all(record["consistency"][name]["passed"] for name in expected)


def test_a_perturbed_part_basis_or_gram_is_flagged() -> None:
    rows = np.arange(10)
    good = [(np.arange(10, 15), np.arange(15, 20))]
    assert alp.check_parts(20, rows, good)["passed"]
    overlap = alp.check_parts(20, rows, [(np.arange(9, 15), np.arange(15, 20))])
    assert not overlap["passed"] and overlap["failing"]["0"]["overlap"] == 1
    uncovered = alp.check_parts(20, rows, [(np.arange(10, 15), np.arange(15, 19))])
    assert not uncovered["passed"] and uncovered["failing"]["0"]["uncovered"] == 1
    rng = np.random.default_rng(1)
    basis = _orthonormal(rng, 12, 3)
    assert alp.check_basis({"a": basis}, {"a": 3})["passed"]
    assert not alp.check_basis({"a": 1.001 * basis}, {"a": 3})["passed"]
    assert alp.check_basis({"a": basis}, {"a": 4})["failing"] == {"a": math.inf}
    products = alp.mean_products([alp.split_products(alp.centred(rng.normal(size=(3, 12))),
                                                     alp.centred(rng.normal(size=(3, 12))), {"a": basis})])
    assert alp.check_gram(products)["passed"]
    perturbed = copy.deepcopy(products)
    perturbed["spaces"]["a"]["gram"][0, 0] += 1e-6
    assert alp.check_gram(perturbed)["failing"] == {"space:a": pytest.approx(1e-6 / max(
        1.0, abs(products["spaces"]["a"]["cross"])))}


def test_an_inconsistent_unit_leaves_the_readings() -> None:
    settings = _paper_settings("cohort")
    good = _record("n80", 0, cross=(1.0, 4.0))
    inconsistent = {**_record("n80", 1, cross=(9.0, 9.5)), "status": "inconsistent", "inconsistent_checks": ["gram"]}
    failed = {"point": "n80", "replicate": 2, "status": "failed", "error_type": "ValueError",
              "error_message": "rank-deficient"}
    rows = alp.capture_rows([good, inconsistent, failed], settings)
    row = next(row for row in rows if row["point"] == "n80" and row["space"] == "pls_cv" and row["reading"] == "total")
    assert row["estimate"] == pytest.approx(0.25) and row["units"] == 1


def test_two_label_permutations_give_equal_readings(cohort, settings) -> None:
    built = aac.build_unit(cohort, settings, "n40", 0)
    seeds = alp.split_seeds(settings, alp.Unit("n40", 0))
    other = dict(built)
    other["metadata"] = built["metadata"].copy()
    stages = other["metadata"]["stage"].astype(str).to_numpy()
    real = other["metadata"]["sex"].astype(str).to_numpy()
    other["metadata"]["sex"] = np.asarray(aac.pc.permute_within_stage(real, stages, cohort.contract.stage_order,
                                                                      np.random.default_rng(99))).astype(str)
    assert not np.array_equal(other["metadata"]["sex"].to_numpy(), built["metadata"]["sex"].to_numpy())
    readings = []
    for current in (built, other):
        parts = alp.evaluate_parts(cohort, settings, current, seeds=seeds, analysis_seed=1, timing={})
        parts.pop("checks")
        readings.append(json.dumps(alp._listed(parts), sort_keys=True))
    assert readings[0] == readings[1]


def test_a_rank_deficient_half_fails_its_unit(monkeypatch, cohort, settings) -> None:
    monkeypatch.setattr(gain, "covariate_columns",
                        lambda metadata, contract: np.column_stack([metadata["age"], metadata["age"]]))
    record = alp.evaluate_cohort_unit(alp.Unit("n40", 0), settings, cohort)
    assert record["status"] == "failed" and "rank-deficient" in record["error_message"]
    assert record["error_type"] == "ValueError"


def test_a_record_holds_no_participant_data(records, cal_records, cohort) -> None:
    ids = cohort.metadata["rid"].tolist()
    features = sum(STAND_IN_FEATURES.values())
    for key, record in (*records.items(), *cal_records.items()):
        text = json.dumps(record)
        assert not any(identifier in text for identifier in ids), key
        sizes = {record["participants"], len(ids), record["halves"]["remainder"], *record["halves"]["sizes"],
                 features}
        for size in sizes:
            assert _walk(record, size) == [], (key, size)
        assert not {"rows", "labels", "configuration", "features", "metadata", "blocks"} & set(record)


# ── 3.1–3.2 The calibration ──────────────────────────────────────────────────


def test_the_generator(cohort) -> None:
    generator = alp.fit_generator(cohort)
    contract = cohort.contract
    clean = alp.synthetic_features(generator, 1.0, 5, noise_scale=0.0)
    recovered = alp.configuration(clean, cohort.metadata, contract)
    assert np.max(np.abs(recovered - generator["centred"])) < 1e-10
    zero = alp.configuration(alp.synthetic_features(generator, 0.0, 5, noise_scale=0.0), cohort.metadata, contract)
    assert np.max(np.abs(zero)) < 1e-10
    first, second = (alp.synthetic_features(generator, 1.0, seed) for seed in (5, 6))
    assert not np.allclose(first, second)
    assert np.array_equal(first, alp.synthetic_features(generator, 1.0, 5))
    synthetic = alp.synthetic_cohort(cohort, first)
    assert synthetic.metadata.equals(cohort.metadata)
    assert [synthetic.blocks[name].shape for name in contract.block_names] == [
        cohort.blocks[name].shape for name in contract.block_names]
    stages = gain.stage_positions(cohort.metadata, contract)
    covariates = gain.covariate_columns(cohort.metadata, contract)
    design = np.hstack([(stages[:, None] == np.arange(N_STAGES)).astype(float), covariates - covariates.mean(0)])
    cross = np.zeros_like(generator["covariance"])
    dof = 0
    for seed in range(100):
        values = alp.synthetic_features(generator, 0.67, 1000 + seed)
        residuals = values - design @ np.linalg.lstsq(design, values, rcond=None)[0]
        cross += residuals.T @ residuals
        dof += values.shape[0] - design.shape[1]
    error = np.linalg.norm(cross / dof - generator["covariance"]) / np.linalg.norm(generator["covariance"])
    assert error < 0.1, error


def test_the_calibration_enumeration_and_seeds(cal_settings) -> None:
    units = alp.enumerate_units(cal_settings)
    assert len(units) == 3 * 2 * 2 * 2 and len({alp.unit_key(unit) for unit in units}) == len(units)
    assert alp.unit_key(units[0]) == "s0/c0/n40/0"
    assert alp.noise_seed(cal_settings, 0.67, 0) != alp.noise_seed(cal_settings, 1.0, 0)
    assert alp.noise_seed(cal_settings, 0.67, 0) != alp.noise_seed(cal_settings, 0.67, 1)
    one, two = alp.Unit("n40", 0, 0.67, 0), alp.Unit("n40", 1, 0.67, 0)
    assert alp.calibration_seed(cal_settings, one, "subsample", 0) != alp.calibration_seed(cal_settings, two,
                                                                                         "subsample", 0)


def test_a_units_truth_is_the_planted_share_in_its_basis(cohort, cal_settings, cal_records) -> None:
    unit = alp.Unit("n40", 0, 1.0, 0)
    record = cal_records[alp.unit_key(unit)]
    generator = alp.generator_for(cohort)
    synthetic = alp.synthetic_cohort(cohort, alp.synthetic_features(
        generator, 1.0, alp.noise_seed(cal_settings, 1.0, 0)))
    built = alp.build_calibration_unit(synthetic, cal_settings, unit)
    dataset = aac.make_dataset(built["blocks"], built["metadata"], cohort.contract)
    production = alp.fit_analysis(dataset, cal_settings, 1, None)
    forced = {rank: alp.fit_analysis(dataset, cal_settings, 1, rank) for rank in cal_settings["fixed_ranks"]}
    bases, _, _ = alp.space_bases(production, forced, cal_settings, N_STAGES)
    scale = gain.joint_scale(gain.sample_scalers(built["blocks"], cohort.contract), cohort.contract)
    planted = generator["centred"] / scale / np.asarray(production.latent_artifacts.model._x_std)
    for space, basis in bases.items():
        share = float(np.sum((planted @ basis) ** 2) / np.sum(planted**2))
        assert alp.ratio(*alp.truth_reading(record, space, "total", N_STAGES)) == pytest.approx(share, rel=1e-10)
    assert record["planted_energy"] == pytest.approx(float(np.sum(planted**2)), rel=1e-10)
    assert cal_records[alp.unit_key(alp.Unit("n88", 1, 0.0, 1))]["truth"] is None
    assert cal_records[alp.unit_key(alp.Unit("n88", 1, 0.0, 1))]["planted_energy"] == 0.0


def test_a_calibration_unit_is_drawn_within_cells(cohort, cal_settings) -> None:
    unit = alp.Unit("n88", 1, 0.67, 1)
    built = alp.build_calibration_unit(cohort, cal_settings, unit)
    assert built["table"] == alp.point_spec(cal_settings, "n88")["table"]
    again = alp.build_calibration_unit(cohort, cal_settings, unit)
    assert np.array_equal(built["rows"], again["rows"])
    other = alp.build_calibration_unit(cohort, cal_settings, alp.Unit("n88", 2, 0.67, 1))
    assert not np.array_equal(built["rows"], other["rows"])


# ── 4.1 Pre-flights ──────────────────────────────────────────────────────────


def test_both_preflights_pass_on_the_stand_in(prepared, settings, cal_settings) -> None:
    record = json.loads((prepared / "cohort" / alp.PREFLIGHT_NAME).read_text())
    assert record["passed"] and record["settings_hash"] == lsc.settings_hash(settings)
    assert [unit["unit"] for unit in record["units"]] == ["cohort/0", "cohort/1"]
    for unit in record["units"]:
        assert unit["fidelity"]["passed"] and unit["forced_rank"] == settings["cohort_rank"]
        assert "cross_validated" not in unit and "forced" not in unit
    reliability = record["reliability"]
    assert reliability["splits"] == SPLITS and -1.0 < reliability["correlation_mean"] < 1.0
    assert len(reliability["cross_gram_eigenvalues"]) == N_STAGES
    calibration = json.loads((prepared / "calibration" / alp.PREFLIGHT_NAME).read_text())
    assert calibration["passed"] and calibration["settings_hash"] == lsc.settings_hash(cal_settings)
    assert [unit["unit"] for unit in calibration["units"]] == ["s1/c0/n40/0", "s1/c0/n88/0"]
    assert all(unit["worst_deviation"] <= 1e-8 for unit in calibration["units"])


@pytest.mark.parametrize("change", ["rank", "auroc"])
def test_a_failed_cohort_preflight_stops_the_run(tmp_path, export_dir, contract, change) -> None:
    wrong = dataclasses.replace(contract, **(
        {"fidelity_rank": contract.fidelity_rank + 1} if change == "rank"
        else {"fidelity_auroc": round(contract.fidelity_auroc + 0.01, 6)}))
    cohort = aac.load_cohort(export_dir, wrong)
    settings = _cohort_settings(cohort)
    result = alp.preflight(tmp_path, settings, cohort)
    assert not result["passed"] and not all(unit["fidelity"]["passed"] for unit in result["units"])
    with pytest.raises(SystemExit, match="pre-flight failed"):
        alp.run_shard(tmp_path, settings, cohort, 0, 1)


def test_a_failed_calibration_preflight_stops_the_run(monkeypatch, tmp_path, cohort, cal_settings) -> None:
    original = alp.truth_products

    def biased(planted, bases):
        truth = original(planted, bases)
        truth["spaces"]["pca"]["cross"] *= 1.01
        return truth

    monkeypatch.setattr(alp, "truth_products", biased)
    result = alp.preflight(tmp_path, cal_settings, cohort)
    assert not result["passed"] and all("pca:total" in unit["failing"] for unit in result["units"])
    with pytest.raises(SystemExit, match="pre-flight failed"):
        alp.run_shard(tmp_path, cal_settings, cohort, 0, 1)


def test_the_run_waits_on_a_preflight_under_its_settings(tmp_path, cohort, settings, cal_settings, prepared) -> None:
    for current, profile in ((settings, "cohort"), (cal_settings, "calibration")):
        with pytest.raises(SystemExit, match=f"run `preflight --profile {profile}` first"):
            alp.run_shard(tmp_path / "fresh", current, cohort, 0, 1)
    out_dir = tmp_path / "out"
    shutil.copytree(prepared, out_dir)
    path = out_dir / "cohort" / alp.PREFLIGHT_NAME
    stored = json.loads(path.read_text())
    stored["settings_hash"] = "0" * 64
    path.write_text(json.dumps(stored))
    with pytest.raises(SystemExit, match="written under other settings"):
        alp.run_shard(out_dir, settings, cohort, 0, 1)


# ── 4.2 Merge (hand-built records) ───────────────────────────────────────────


def _gram(value: float, n_stages: int = N_STAGES) -> list:
    """A centred stage-space Gram matrix with trace ``value`` split equally over the centred directions."""

    return (value / (n_stages - 1) * (np.eye(n_stages) - 1.0 / n_stages)).tolist()


def _record(point: str, replicate: int, *, cross=(1.0, 2.0), strength=None, cohort=None, truth=None,
            rank: int = 2) -> dict:
    inside, total = cross
    spaces = {name: {"dimension": 2, "cross": inside, "gram": _gram(inside), "noise": 0.5, "naive_remainder": 0.4,
                     "naive_sample": 0.8} for name in alp.space_names((2, 8, 19))}
    record = {
        "profile": "cohort" if strength is None else "calibration", "point": point, "replicate": replicate,
        "status": "completed", "production": {"rank": rank, "cv_mean_auroc": 0.6},
        "halves": {"remainder": 100, "sizes": [50, 50]}, "reliability": 0.2,
        "cross_products": {"total": {"cross": total, "gram": _gram(total), "noise": 1.0, "naive_remainder": 1.0,
                                     "naive_sample": 1.0}, "spaces": spaces},
        "runtime_seconds": {"unit": 5.0},
    }
    if strength is not None:
        record.update({"strength": strength, "cohort": cohort})
        if truth is not None:
            record["truth"] = {"total": {"cross": truth[1], "gram": _gram(truth[1])},
                               "spaces": {name: {"cross": truth[0], "gram": _gram(truth[0])} for name in spaces}}
        else:
            record["truth"] = None
    return record


def _paper_settings(profile: str, **changes) -> dict:
    settings = {"profile": profile, "points": aac.design_points(ADNI_TABLE, alp.POINT_SIZES),
                "spaces": alp.space_specs(), "contract": {"stage_order": ["CN", "MCI", "AD"]},
                "replicates": alp.REPLICATES, "strengths": list(alp.STRENGTHS), "cohorts": alp.CALIBRATION_COHORTS,
                "units": alp.CALIBRATION_UNITS}
    settings.update(changes)
    return settings


def _capture(rows, point, reading="total", space="pls_cv", strength=None):
    return next(row for row in rows if row["point"] == point and row["space"] == space and row["reading"] == reading
                and (strength is None or math.isclose(row["strength"], strength)))


def test_the_estimate_is_a_ratio_of_sums() -> None:
    settings = _paper_settings("cohort")
    records = [_record("n80", 0, cross=(1.0, 4.0)), _record("n80", 1, cross=(3.0, 2.0))]
    row = _capture(alp.capture_rows(records, settings), "n80")
    assert row["estimate"] == pytest.approx(4.0 / 6.0)
    assert row["unit_median"] == pytest.approx((0.25 + 1.5) / 2)
    assert row["numerator_sum"] == 4.0 and row["denominator_sum"] == 6.0 and row["denominator_mean"] == 3.0
    for contrast in alp.CONTRASTS:  # half of each energy in each contrast at three stages (``_gram``)
        assert _capture(alp.capture_rows(records, settings), "n80", contrast)["estimate"] == pytest.approx(4.0 / 6.0)


def test_an_estimate_outside_the_unit_interval_is_written_unchanged(tmp_path) -> None:
    settings = _paper_settings("cohort")
    records = [_record("n80", 0, cross=(3.0, 2.0)), _record("n200", 0, cross=(1.0, -4.0))]
    rows = alp.capture_rows(records, settings)
    path = tmp_path / "capture_table.csv"
    lsc._write_csv(path, rows, alp.CAPTURE_COLUMNS)
    written = {(row["point"], row["space"], row["reading"]): row for row in lsc._read_csv(path)}
    assert written[("n80", "pls_cv", "total")]["estimate"] == "1.5"
    assert written[("n200", "pls_cv", "total")]["estimate"] == "-0.25"


def test_the_calibration_table_pools_each_cohort() -> None:
    settings = _paper_settings("calibration", strengths=[0.67])
    records = [
        _record("n80", 0, cross=(1.0, 2.0), strength=0.67, cohort=0, truth=(1.0, 2.0)),
        _record("n80", 1, cross=(1.0, 2.0), strength=0.67, cohort=0, truth=(1.0, 2.0)),
        _record("n80", 0, cross=(3.0, 4.0), strength=0.67, cohort=1, truth=(2.0, 4.0)),
        _record("n80", 1, cross=(3.0, 2.0), strength=0.67, cohort=1, truth=(2.0, 4.0)),
        {**_record("n80", 2, cross=(50.0, 1.0), strength=0.67, cohort=1, truth=(2.0, 4.0)), "status": "inconsistent"},
    ]
    rows = alp.calibration_rows(records, settings)
    row = next(row for row in rows if row["point"] == "n80" and row["space"] == "pls_cv" and row["reading"] == "total")
    assert row["cohorts"] == 2
    assert row["mean_estimate"] == pytest.approx(0.75)  # cohorts pooled to 0.5 and 1.0; all units pooled give 0.8
    assert row["sd_estimate"] == pytest.approx(float(np.std([0.5, 1.0], ddof=1)))
    assert row["mean_truth"] == pytest.approx(0.5) and row["bias"] == pytest.approx(0.25)
    assert row["mean_energy"] == pytest.approx(2.5) and row["sd_energy"] == pytest.approx(float(np.std([2.0, 3.0],
                                                                                                      ddof=1)))
    capture = _capture(alp.capture_rows(records, settings), "n80", strength=0.67)
    assert capture["estimate"] == pytest.approx(0.8) and capture["truth"] == pytest.approx(7.0 / 14.0)


def test_the_point_table_and_accounting() -> None:
    settings = _paper_settings("cohort", replicates=3)
    records = [_record("n80", 0, rank=2), _record("n80", 1, rank=3),
               {**_record("n80", 2), "status": "inconsistent", "inconsistent_checks": ["basis"],
                "consistency": {"basis": {"deviation": 1.0, "failing": {"pca": 1.0}}}},
               {"point": "n200", "replicate": 0, "status": "failed", "error_type": "ValueError",
                "error_message": "rank-deficient", "runtime_seconds": {"unit": 1.0}}]
    rows = alp.point_rows(records, settings)

    def value(quantity, point="n80", **keys):
        return [row["value"] for row in rows if row["point"] == point and row["quantity"] == quantity
                and all(row.get(key) == wanted for key, wanted in keys.items())]

    assert value("cv_rank_units", key="2") == [1] and value("cv_rank_units", key="3") == [1]
    assert value("units_completed") == [2] and value("units_inconsistent") == [1]
    assert value("units_failed", point="n200") == [1]
    assert value("energy_mean", key="total") == [2.0]
    assert value("naive_median", space="pls_cv", key="naive_remainder") == [pytest.approx(0.4)]
    assert value("cross_gram_eigenvalue", key="1") == [pytest.approx(1.0)]
    summary = alp.accounting(records, settings, duplicates=0)
    assert summary["by_point"]["n80"] == {"expected": 3, "complete": 3, "failed": 0, "inconsistent": 1, "usable": 2}
    assert summary["by_point"]["n200"]["failed"] == 1 and summary["failed_unit_errors"][0]["unit"] == "n200/0"
    assert summary["inconsistent_unit_checks"][0]["checks"]["basis"]["failing"] == {"pca": 1.0}


# ── 4.3 Reproductions ────────────────────────────────────────────────────────


def _committed():
    return lsc._read_csv(alp.ANCHOR_DIR / "point_table.csv"), lsc._read_csv(alp.GAIN_DIR / "cell_table.csv")


def _reproducing_rows(anchor_rows, gain_rows) -> list[dict]:
    rows = [{"point": row["point"], "quantity": "cv_rank_units", "space": "", "key": row["key"],
             "value": row["value"]} for row in anchor_rows if row["quantity"] == "cv_rank_units"
            and row["point"] in alp.REPRODUCED_RANK_POINTS]
    for (reading, space), key in alp.NAIVE_PARENT_KEYS.items():
        value = next(row["value"] for row in gain_rows if row["point"] == "n750" and row["cell"] == "anchor"
                     and row["quantity"] == "capture_median" and row["key"] == key)
        rows.append({"point": "n750", "quantity": "naive_median", "space": space, "key": reading, "value": value})
    return rows


def test_the_reproduction_passes_on_the_committed_values() -> None:
    anchor_rows, gain_rows = _committed()
    result = alp.reproduction(_reproducing_rows(anchor_rows, gain_rows), anchor_rows, gain_rows)
    assert result["passed"] and not result["differing"] and len(result["rows"]) == 3 + 4
    naive = {(row["reading"], row["space"]): row["committed"] for row in result["rows"]
             if row["quantity"] == "naive_median"}
    assert naive == {("naive_remainder", "pls_cv"): "0.53603", ("naive_remainder", "predictive"): "0.108803",
                     ("naive_sample", "pls_cv"): "0.925475", ("naive_sample", "predictive"): "0.548427"}


def test_a_differing_count_or_median_is_recorded() -> None:
    anchor_rows, gain_rows = _committed()
    rows = _reproducing_rows(anchor_rows, gain_rows)
    count = next(row for row in rows if row["point"] == "n200" and row["key"] == "2")
    count["value"] = str(int(count["value"]) + 1)
    median = next(row for row in rows if row["quantity"] == "naive_median" and row["key"] == "naive_remainder"
                  and row["space"] == "pls_cv")
    median["value"] = "0.53604"
    result = alp.reproduction(rows, anchor_rows, gain_rows)
    assert not result["passed"] and len(result["differing"]) == 2
    ranks, naive = result["differing"]
    assert ranks["point"] == "n200" and ranks["quantity"] == "cv_rank_units"
    assert ranks["observed"]["2"] == ranks["committed"]["2"] + 1 == 229
    assert naive == {"point": "n750", "quantity": "naive_median", "space": "pls_cv", "reading": "naive_remainder",
                     "parent_key": "heldout:span", "committed": "0.53603", "observed": "0.53604", "equal": False}


def test_the_merge_writes_the_tables_and_the_reproduction(tmp_path, cohort, settings, cal_settings, prepared,
                                                         export_dir) -> None:
    out_dir = tmp_path / "out"
    shutil.copytree(prepared, out_dir)
    assert alp.run_shard(out_dir, settings, cohort, 0, 1) == 4
    assert alp.run_shard(out_dir, cal_settings, cohort, 0, 1) == 24
    before = _snapshot(export_dir)
    summary = alp.merge(out_dir, "cohort", export_dir)
    assert summary["complete_units"] == 4 and summary["failed_units"] == summary["inconsistent_units"] == 0
    assert summary["reproduction_passed"] is False  # the stand-in's points are not ADNI's
    calibration = alp.merge(out_dir, "calibration", export_dir)
    assert calibration["complete_units"] == 24 and calibration["inconsistent_units"] == 0
    rows = lsc._read_csv(out_dir / "calibration" / "calibration_table.csv")
    assert {row["cohorts"] for row in rows if row["reading"] == "total"} == {"2"}
    captures = lsc._read_csv(out_dir / "cohort" / "capture_table.csv")
    assert len(captures) == 2 * len(settings["spaces"]) * len(alp.READINGS)
    assert _snapshot(export_dir) == before
    decision = alp.decide(out_dir, export_dir)
    assert decision["verdict"] == "gate_failed" and "reproduction" in decision["gate"]["causes"]
    assert decision["indicative"]["indicative"]
    assert decision["preflight"]["cohort"]["passed"] and decision["preflight"]["calibration"]["passed"]
    ids = cohort.metadata["rid"].tolist()
    for path in [*out_dir.glob("*/*.csv"), *out_dir.glob("*/*.json"), out_dir / alp.DECISION_NAME]:
        text = path.read_text()
        assert not any(identifier in text for identifier in ids), path


# ── 5.1 Decision (synthetic tables) ──────────────────────────────────────────


POINTS = ("n80", "n200", "n500", "n750")


def _decision_inputs(*, energy=1.0, null_sd=0.1, sd=None, bias=0.0, estimate=0.4, contrasts=None, failed=None,
                     reproduced=True, replicates=500, cohorts=24, units=25):
    """Settings, the cohort's capture table, the calibration table, accountings and the reproduction."""

    sd = sd or {}
    contrasts = contrasts or {}
    points = aac.design_points(ADNI_TABLE, alp.POINT_SIZES)
    settings = {
        "cohort": {"profile": "cohort", "points": points, "replicates": replicates},
        "calibration": {"profile": "calibration", "points": points, "cohorts": cohorts, "units": units,
                        "strengths": list(alp.STRENGTHS)},
    }
    captures, calibration = [], []
    for point in POINTS:
        for reading in ("total", *alp.CONTRASTS):
            point_estimate, point_energy = contrasts.get((point, reading), (estimate, energy / 2.0)) \
                if reading != "total" else (estimate, energy)
            captures.append({"point": point, "space": "pls_cv", "reading": reading, "estimate": str(point_estimate),
                             "denominator_mean": str(point_energy)})
            for strength in alp.STRENGTHS:
                calibration.append({
                    "strength": str(strength), "point": point, "space": "pls_cv", "reading": reading,
                    "sd_estimate": str(sd.get(point, 0.05)),
                    "bias": "" if strength == 0.0 else str(bias if isinstance(bias, float) else bias.get(
                        (strength, point), 0.0)),
                    "sd_energy": str(null_sd),
                })
    failed = failed or {}
    accounting = {
        "cohort": {"by_point": {point: {"expected": replicates, "complete": replicates, "failed": failed.get(point, 0),
                                        "inconsistent": 0} for point in POINTS}, "missing_units": 0},
        "calibration": {"by_point": {point: {"expected": 3 * cohorts * units, "complete": 3 * cohorts * units,
                                             "failed": 0, "inconsistent": 0} for point in POINTS},
                        "missing_units": 0},
    }
    return settings, captures, calibration, accounting, {"passed": reproduced, "differing": [] if reproduced else [
        {"point": "n80", "quantity": "cv_rank_units"}]}


def _decide(**kwargs) -> dict:
    return alp.decision(*_decision_inputs(**kwargs))


def test_the_five_verdicts() -> None:
    assert _decide()["verdict"] == "measured"
    assert _decide(sd={"n750": 0.1})["verdict"] == "partly_measured"
    assert _decide(sd={point: 0.1 for point in POINTS})["verdict"] == "unresolved"
    assert _decide(energy=0.15)["verdict"] == "no_program"
    assert _decide(reproduced=False)["verdict"] == "gate_failed"
    for verdict in alp.VERDICTS:
        assert alp.CONSEQUENCES[verdict]


def test_the_gates_three_causes_separately() -> None:
    units = _decide(failed={"n200": 6})
    assert units["verdict"] == "gate_failed" and units["gate"]["causes"] == ["units"]
    assert not units["gate"]["units"]["cohort/n200"]["passed"] and units["gate"]["units"]["cohort/n80"]["passed"]
    assert _decide(failed={"n200": 5})["gate"]["passed"]
    reproduction = _decide(reproduced=False)
    assert reproduction["gate"]["causes"] == ["reproduction"]
    bias = _decide(bias={(1.0, "n500"): 0.06})
    assert bias["verdict"] == "gate_failed" and bias["gate"]["causes"] == ["bias"]
    assert not bias["gate"]["bias"]["s1/n500"]["passed"] and bias["gate"]["bias"]["s0.67/n500"]["passed"]
    assert _decide(bias={(0.67, "n80"): -0.05})["gate"]["passed"]
    missing = alp.decision(_decision_inputs()[0], [], [], {"cohort": None, "calibration": None}, None)
    assert missing["verdict"] == "gate_failed" and "inputs" in missing["gate"]["causes"]


def test_a_resolution_exactly_at_the_limit_resolves() -> None:
    result = _decide(sd={"n80": 0.075})
    assert result["resolution"]["n80"]["half_width"] == 0.15 and result["resolution"]["n80"]["resolved"]
    assert not _decide(sd={"n80": 0.0751})["resolution"]["n80"]["resolved"]


def test_the_program_test_is_strict() -> None:
    assert _decide(energy=0.2)["verdict"] == "no_program"
    assert _decide(energy=0.2001)["program"]["present"]


def test_each_contrast_grade() -> None:
    result = _decide(contrasts={("n80", "net_displacement"): (0.8, 0.5), ("n80", "interior"): (0.2, 0.5),
                                ("n200", "net_displacement"): (0.5, 0.5), ("n200", "interior"): (0.9, 0.1)})
    grades = result["contrast_grades"]
    assert grades["n80"]["net_displacement"]["grade"] == "captured"
    assert grades["n80"]["interior"]["grade"] == "missed"
    assert grades["n200"]["net_displacement"]["grade"] == "undetermined"
    assert grades["n200"]["interior"] == {**grades["n200"]["interior"], "present": False, "grade": None}
    licenses = result["licenses"]
    assert licenses["follow_up"]["net_displacement"]["licensed"] and licenses["follow_up"]["interior"]["licensed"]
    assert licenses["plants"] == {"magnitude": True, "orientation": True, "shape": True}
    assert [entry["point"] for entry in licenses["stated_captures"]] == list(POINTS)
    only_net = _decide(contrasts={(point, "interior"): (0.4, 0.1) for point in POINTS})["licenses"]
    assert only_net["plants"] == {"magnitude": True, "orientation": False, "shape": False}
    unresolved = _decide(sd={point: 0.1 for point in POINTS})["licenses"]
    assert unresolved["stated_captures"] == [] and not unresolved["plants"]["magnitude"]


@pytest.mark.parametrize("change, reason", [
    ({"replicates": 499}, "499 replicates"),
    ({"cohorts": 23}, "23 cohorts"),
    ({"units": 24}, "24 units"),
])
def test_the_indicative_mark_below_each_committed_count(change, reason) -> None:
    assert not _decide()["indicative"]["indicative"]
    marked = _decide(**change)["indicative"]
    assert marked["indicative"] and any(reason in text for text in marked["reasons"])


# ── The documented commands ──────────────────────────────────────────────────


def test_the_docstrings_commands_run_against_the_stand_in(monkeypatch, tmp_path, export_dir, contract) -> None:
    monkeypatch.setitem(aac.CONTRACTS, "adni", contract)
    monkeypatch.setattr(alp, "POINT_SIZES", STAND_IN_POINTS)
    monkeypatch.setattr(alp, "CV", SMALL_CV)
    monkeypatch.setattr(alp, "FIXED_RANKS", FIXED)
    monkeypatch.setattr(alp, "SPLITS", SPLITS)
    commands = [line.strip() for line in alp.__doc__.splitlines()
                if line.strip().startswith("python scripts/adni_latent_plasmode_check.py")]
    assert [shlex.split(command)[2] for command in commands] == ["preflight", "run", "merge", "preflight", "run",
                                                                 "merge", "decide"]
    out_dir = tmp_path / "smoke"
    before = _snapshot(export_dir)
    for command in commands:
        argv = shlex.split(command.replace("$EXPORT", str(export_dir)).replace("$OUT", str(out_dir)))[2:]
        argv = [value if value != "8" else "2" for value in argv]  # two workers suffice here
        assert alp.main(argv) == 0, command
    decision = json.loads((out_dir / alp.DECISION_NAME).read_text())
    assert decision["complete"]["cohort"]["complete_units"] == 4 * 2
    assert decision["complete"]["calibration"]["complete_units"] == 3 * 2 * 2 * 2
    assert decision["complete"]["cohort"]["inconsistent_units"] == 0
    assert decision["complete"]["calibration"]["inconsistent_units"] == 0
    assert decision["verdict"] == "gate_failed" and decision["indicative"]["indicative"]
    assert _snapshot(export_dir) == before
    merged = [json.loads(line) for line in (out_dir / "calibration" / "merged.jsonl").read_text().splitlines()]
    assert all(_walk(record, sum(STAND_IN_FEATURES.values())) == [] for record in merged)
    assert isinstance(pd.read_csv(out_dir / "calibration" / "calibration_table.csv"), pd.DataFrame)
