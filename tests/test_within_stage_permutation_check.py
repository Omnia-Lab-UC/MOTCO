"""``scripts/within_stage_permutation_check.py``: settings, units, arms, checks, merge, reproduction, decision.

The simulation units run on the pilot's ``seaad`` design at fast CV knobs; the
plasmode profile runs on the plasmode tests' synthetic stand-in export with
stand-in copies of the parents. No real data are read.
"""

from __future__ import annotations

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


wsp = _script("within_stage_permutation_check")
psc, pse, lsc, pc, plsc, grade = wsp.psc, wsp.pse, wsp.lsc, wsp.pc, wsp.plsc, wsp.grade

REPO = Path(__file__).resolve().parents[1]
PERMUTATIONS = 9
FAST_CV = {"cv1_splits": 3, "cv2_splits": 4, "n_repeats": 1, "max_components": 6, "random_state": 1203}


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


# ── Stand-ins ────────────────────────────────────────────────────────────────


@pytest.fixture(scope="module")
def stand_in(tmp_path_factory):
    """The stand-in export, its plasmode and production parents, a stand-in pilot (500/100) and diagnosis."""

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
    shutil.copy(wsp.PILOT_DIR / "plasmode" / "rate_table.csv", pilot_dir / "rate_table.csv")
    diagnosis = pse.build_settings(export, pilot_dir=pilot_dir, parent_dir=parent_dir, production_dir=production_dir,
                                   replicates=2, plant_replicates=1, permutations=PERMUTATIONS)
    diagnosis_dir = tmp_path_factory.mktemp("diagnosis")
    _write_settings(diagnosis_dir, diagnosis)
    shutil.copy(wsp.DIAGNOSIS_DIR / "rate_table.csv", diagnosis_dir / "rate_table.csv")
    return {"directory": directory, "export": export, "parent_dir": parent_dir, "production_dir": production_dir,
            "pilot_dir": pilot_dir, "diagnosis_dir": diagnosis_dir}


def _plasmode(stand_in, **kwargs) -> dict:
    kwargs.setdefault("replicates", 2)
    kwargs.setdefault("fresh_replicates", 2)
    kwargs.setdefault("permutations", PERMUTATIONS)
    for key in ("pilot_dir", "diagnosis_dir", "parent_dir", "production_dir"):
        kwargs.setdefault(key, stand_in[key])
    return wsp.plasmode_settings(stand_in["export"], **kwargs)


@pytest.fixture(scope="module")
def plasmode_settings(stand_in):
    return _plasmode(stand_in)


@pytest.fixture(scope="module")
def anchor_record(stand_in, plasmode_settings):
    return wsp.evaluate_plasmode_unit(lsc.Unit("plasmode", "anchor", 1), plasmode_settings, stand_in["export"])


@pytest.fixture(scope="module")
def fresh_record(stand_in, plasmode_settings):
    return wsp.evaluate_plasmode_unit(lsc.Unit("plasmode", "anchor_fresh", 500), plasmode_settings,
                                      stand_in["export"])


def _fast(settings: dict) -> dict:
    settings = json.loads(json.dumps(settings))
    units = psc._harness_units(settings["pilot"])
    units["pls_cv_params"] = dict(FAST_CV)
    return settings


@pytest.fixture(scope="module")
def small_simulation():
    return _fast(wsp.simulation_settings(replicates=2, permutations=19, designs=["seaad"]))


@pytest.fixture(scope="module")
def harness_record(small_simulation):
    return wsp.evaluate_harness_unit(lsc.Unit("seaad", "orient_e0.25_phi0", 1), small_simulation)


def _counting_cv(monkeypatch) -> list:
    import motco.simulations.evaluation as evaluation_module

    calls: list = []
    original = evaluation_module.plsda_doubleCV

    def counting(*args, **kwargs):
        calls.append(1)
        return original(*args, **kwargs)

    monkeypatch.setattr(evaluation_module, "plsda_doubleCV", counting)
    return calls


# ── 4.1 Settings, parents, resume ────────────────────────────────────────────


def test_settings_bind_the_parents(small_simulation, plasmode_settings) -> None:
    assert {"comparison_settings", "comparison_calibration", "pilot_simulation_settings"} <= set(
        small_simulation["parents"])
    assert {"pilot_plasmode_settings", "pilot_plasmode_rates", "diagnosis_settings", "diagnosis_rates",
            "plasmode_settings", "plasmode_summary"} <= set(plasmode_settings["parents"])
    weak = wsp.weak_size_settings(replicates=2, permutations=PERMUTATIONS)
    assert {"pilot_weak_size_settings", "pilot_weak_size_rates", "span_split_settings"} <= set(weak["parents"])
    for settings in (small_simulation, plasmode_settings, weak):
        assert settings["permutation_strata"] == "stage" and settings["arms"] == list(wsp.ARMS)


def test_the_default_parents_are_the_committed_runs() -> None:
    weak = wsp.weak_size_settings(replicates=2, permutations=PERMUTATIONS)
    assert weak["parents"]["pilot_weak_size_rates"]["path"] == "results/predictive-subspace-2026-11-12/weak_size/" \
                                                                "rate_table.csv"
    assert wsp.DIAGNOSIS_DIR == REPO / "results" / "predictive-size-excess-2026-11-13"


def test_resume_skips_complete_units_and_other_settings_are_refused(tmp_path) -> None:
    settings = wsp.dispersion_settings(replicates=2, permutations=PERMUTATIONS)
    assert wsp.run_shard(tmp_path, settings, 0, 1, limit=3) == 3
    assert wsp.run_shard(tmp_path, settings, 0, 1) == 17
    lines = (tmp_path / "dispersion" / "records" / "shard_000_of_001.jsonl").read_text().splitlines()
    assert len(lines) == 20
    assert wsp.run_shard(tmp_path, settings, 0, 1) == 0
    with pytest.raises(SystemExit, match="replicates"):
        wsp.run_shard(tmp_path, wsp.dispersion_settings(replicates=3, permutations=PERMUTATIONS), 0, 1)


def test_a_changed_parent_is_refused_by_name(tmp_path, stand_in) -> None:
    diagnosis = tmp_path / "diagnosis"
    shutil.copytree(stand_in["diagnosis_dir"], diagnosis)
    settings = _plasmode(stand_in, diagnosis_dir=diagnosis)
    with (diagnosis / "rate_table.csv").open("a") as handle:
        handle.write("\n")
    with pytest.raises(SystemExit, match="diagnosis_rates"):
        wsp.verify(settings)


def test_a_pilot_with_other_units_is_refused(tmp_path, stand_in) -> None:
    pilot = tmp_path / "pilot"
    stored = json.loads((stand_in["pilot_dir"] / "settings.json").read_text())["settings"]
    _write_settings(pilot, {**stored, "permutations": PERMUTATIONS + 1})
    shutil.copy(stand_in["pilot_dir"] / "rate_table.csv", pilot / "rate_table.csv")
    with pytest.raises(SystemExit, match="permutations"):
        _plasmode(stand_in, pilot_dir=pilot)


def test_replicates_outside_the_two_ranges_are_refused(stand_in) -> None:
    with pytest.raises(SystemExit, match="--reps"):
        _plasmode(stand_in, replicates=501)
    with pytest.raises(SystemExit, match="--fresh-reps"):
        _plasmode(stand_in, fresh_replicates=501)


# ── 4.2 Units and seeds ──────────────────────────────────────────────────────


def test_unit_counts(stand_in) -> None:
    assert len(wsp.simulation_units(wsp.simulation_settings())) == 7_800
    assert len(wsp.weak_size_units(wsp.weak_size_settings())) == 1_000
    plasmode = _plasmode(stand_in, replicates=500, fresh_replicates=500)
    units = wsp.plasmode_units(plasmode)
    assert len(units) == 1_000
    assert [u.replicate for u in units if u.cell == "anchor"] == list(range(500))
    assert [u.replicate for u in units if u.cell == "anchor_fresh"] == list(range(500, 1000))
    assert len(wsp.dispersion_units(wsp.dispersion_settings())) == 10_000


def test_a_fresh_simulation_seed_family() -> None:
    ours = {lsc.generator_seed(2038, design, r) for design in psc.SIMULATION_DESIGNS for r in range(100)}
    spent = {
        lsc.generator_seed(base, design, r)
        for base in (2035, 2036, 4242)
        for design in (*lsc.DESIGNS, *psc.SIMULATION_DESIGNS)
        for r in range(500)
    }
    assert len(ours) == len(psc.SIMULATION_DESIGNS) * 100
    assert not ours & spent
    settings = wsp.simulation_settings(replicates=2, permutations=19, designs=["seaad"])
    assert settings["pilot"]["units"]["base_seed"] == 2038


def test_dispersion_seeds_are_the_drivers_own() -> None:
    seeds = {wsp.dispersion_seeds(cell, r)["matrix"] for cell in wsp.dispersion_cells() for r in range(1000)}
    assert len(seeds) == 10_000
    assert 20261116 not in seeds


def test_weak_and_reproduced_plasmode_units_are_the_pilots(stand_in, plasmode_settings, anchor_record) -> None:
    weak = wsp.weak_size_settings()
    assert wsp.weak_size_units(weak) == psc.weak_size_units(psc.weak_size_settings())
    assert anchor_record["seeds"] == pc.unit_seeds(2034, "paper", 1)
    assert plasmode_settings["pilot"]["plasmode"]["stream"] == "paper"


def _paper_stream_replicates(value, found: list) -> None:
    if isinstance(value, dict):
        if value.get("stream") == "paper":
            found.append(value.get("replicates"))
        for item in value.values():
            _paper_stream_replicates(item, found)
    elif isinstance(value, list):
        for item in value:
            _paper_stream_replicates(item, found)


def test_no_committed_run_reads_a_fresh_label_draw() -> None:
    found: list = []
    for path in sorted((REPO / "results").rglob("settings.json")):
        _paper_stream_replicates(json.loads(path.read_text()), found)
    assert found
    assert all(isinstance(count, int) and count <= wsp.FRESH_START for count in found), found


def test_a_changed_export_is_refused(plasmode_settings) -> None:
    class Changed:
        files = {name: "0" * 64 for name in plasmode_settings["pilot"]["export"]["files"]}

    with pytest.raises(SystemExit, match="file hashes differ"):
        wsp.verify(plasmode_settings, Changed())


# ── 4.3 The fitted unit ──────────────────────────────────────────────────────


def test_a_simulation_unit_fits_once(small_simulation, monkeypatch) -> None:
    calls = _counting_cv(monkeypatch)
    record = wsp.evaluate_harness_unit(lsc.Unit("seaad", "anchor", 0), small_simulation)
    assert record["status"] == "completed", record.get("inconsistent_checks")
    assert len(calls) == 1
    assert set(record["arms"]) == set(wsp.ARMS)
    assert {entry["rank"] for entry in record["arms"].values()} == {record["production"]["rank"]}


def test_a_plasmode_unit_makes_no_cross_validation(stand_in, plasmode_settings, monkeypatch) -> None:
    calls = _counting_cv(monkeypatch)
    record = wsp.evaluate_plasmode_unit(lsc.Unit("plasmode", "anchor", 0), plasmode_settings, stand_in["export"])
    assert record["status"] == "completed", record.get("consistency")
    assert calls == []
    for arm, entry in record["arms"].items():
        assert entry["rank"] == plasmode_settings["anchor_rank"], arm


def test_harness_unrestricted_arms_equal_the_pilots(small_simulation, harness_record) -> None:
    theirs = psc.evaluate_harness_unit(lsc.Unit("seaad", "orient_e0.25_phi0", 1), small_simulation["pilot"])
    assert harness_record["status"] == "completed" and theirs["status"] == "completed"
    for arm in wsp.UNRESTRICTED_ARMS:
        assert harness_record["arms"][arm]["pair_statistics"] == theirs["arms"][arm]["pair_statistics"], arm
        assert harness_record["arms"][arm]["p_values"] == theirs["arms"][arm]["p_values"], arm
    assert harness_record["rrpp_seed"] == theirs["rrpp_seed"]


@pytest.mark.parametrize("cell, replicate", [("anchor", 1), ("anchor_fresh", 500)])
def test_plasmode_unrestricted_arms_equal_the_pilots(stand_in, plasmode_settings, anchor_record, fresh_record, cell,
                                                     replicate) -> None:
    ours = anchor_record if cell == "anchor" else fresh_record
    theirs = psc.evaluate_plasmode_unit(lsc.Unit("plasmode", "anchor", replicate), plasmode_settings["pilot"],
                                        stand_in["export"])
    assert ours["status"] == "completed" and theirs["status"] == "completed"
    for arm in wsp.UNRESTRICTED_ARMS:
        assert ours["arms"][arm]["pair_statistics"] == theirs["arms"][arm]["pair_statistics"], arm
        assert ours["arms"][arm]["p_values"] == theirs["arms"][arm]["p_values"], arm
    assert ours["seeds"] == theirs["seeds"]


def test_within_stage_arms_change_only_the_null(anchor_record, harness_record) -> None:
    for record in (anchor_record, harness_record):
        for arm in wsp.UNRESTRICTED_ARMS:
            within = record["arms"][f"{arm}__within_stage"]
            assert within["pair_statistics"] == record["arms"][arm]["pair_statistics"]
            assert within["draws"] == "within_stage" and record["arms"][arm]["draws"] == "unrestricted"
        assert record["consistency"]["same_statistics"]["passed"]
    assert anchor_record["consistency"]["construction_identity"]["passed"]
    assert anchor_record["consistency"]["construction_identity"]["arms"] == list(wsp.WITHIN_STAGE_ARMS)
    assert not harness_record["consistency"]["construction_identity"]["applicable"]


def test_a_record_holds_seeds_and_statistics_only(stand_in, anchor_record) -> None:
    export = stand_in["export"]
    n = len(export.metadata)

    def no_vectors(value) -> bool:
        if isinstance(value, dict):
            return all(no_vectors(item) for item in value.values())
        if isinstance(value, list):
            return len(value) != n and all(no_vectors(item) for item in value)
        return True

    text = json.dumps(anchor_record)
    for donor in export.metadata[export.sample_col]:
        assert donor not in text
    assert "Female" not in text and "Male" not in text
    assert no_vectors(anchor_record)


def test_a_failed_unit_is_recorded_and_the_shard_goes_on(small_simulation, monkeypatch) -> None:
    def broken(*args, **kwargs):
        raise RuntimeError("generator broke")

    monkeypatch.setattr(lsc, "unit_dataset", broken)
    record = wsp.evaluate_harness_unit(lsc.Unit("seaad", "anchor", 0), small_simulation)
    assert record["status"] == "failed" and record["error_message"] == "generator broke"


# ── 4.4 The per-unit checks ──────────────────────────────────────────────────


def test_same_statistics_flags_a_perturbed_arm(anchor_record) -> None:
    arms = json.loads(json.dumps(anchor_record["arms"]))
    assert wsp.same_statistics(arms)["passed"]
    arms["predictive__within_stage"]["pair_statistics"]["delta"] += 1e-12
    check = wsp.same_statistics(arms)
    assert not check["passed"] and check["deviation"] > 0


def test_a_construction_mismatch_is_inconsistent(stand_in, plasmode_settings, monkeypatch) -> None:
    original = pse.infer

    def shifted(*args, **kwargs):
        out = original(*args, **kwargs)
        out["p_values"] = {key: value / 2 for key, value in out["p_values"].items()}
        return out

    monkeypatch.setattr(pse, "infer", shifted)
    record = wsp.evaluate_plasmode_unit(lsc.Unit("plasmode", "anchor", 0), plasmode_settings, stand_in["export"])
    assert record["status"] == "inconsistent"
    assert record["inconsistent_checks"] == ["construction_identity"]
    assert set(record["consistency"]["construction_identity"]["differing"]) == set(wsp.WITHIN_STAGE_ARMS)
    rows = wsp.rate_table([record], plasmode_settings)
    assert all(row["count"] == 0 and row["inconsistent"] == 1 for row in rows)


# ── 4.5 The dispersion profile ───────────────────────────────────────────────


def test_dispersion_settings_and_cells() -> None:
    settings = wsp.dispersion_settings()
    assert settings["cells"] == [f"{p}_d{d}" for p in ("equal", "predictive", "scores", "predictive_doubled",
                                                         "scores_doubled") for d in (2, 5)]
    assert settings["patterns"]["predictive_doubled"] == pytest.approx([1.120**2, 1.278**2, 0.802**2])
    assert settings["patterns"]["scores"] == [0.805, 0.921, 1.137]
    assert settings["replicates"] == 1000 and settings["permutations"] == 999
    assert settings["rejection"] == "p < alpha" and settings["arms"] == ["matrix", "matrix__within_stage"]


def test_dispersion_stage_share_is_the_anchors() -> None:
    settings = wsp.dispersion_settings()
    stacked, stages = [], []
    for replicate in range(400):
        Y, metadata = wsp.dispersion_matrix(settings, "equal_d2", replicate)
        stacked.append(Y)
        stages.append(metadata["stage"].astype(int).to_numpy())
    Y, stage = np.vstack(stacked), np.concatenate(stages)
    means = np.vstack([Y[stage == s].mean(axis=0) for s in range(3)])
    total = np.sum((Y - Y.mean(axis=0)) ** 2)
    between = np.sum((means[stage] - Y.mean(axis=0)) ** 2)
    assert between / total == pytest.approx(0.84, abs=0.01)


@pytest.mark.parametrize("pattern", ["predictive_doubled", "scores"])
def test_dispersion_variances_follow_the_pattern(pattern) -> None:
    settings = wsp.dispersion_settings()
    residuals = {s: [] for s in range(3)}
    for replicate in range(300):
        Y, metadata = wsp.dispersion_matrix(settings, f"{pattern}_d5", replicate)
        stage = metadata["stage"].astype(int).to_numpy()
        for s in range(3):
            rows = Y[stage == s]
            residuals[s].append(rows - rows.mean(axis=0))
    for s, multiplier in enumerate(settings["patterns"][pattern]):
        pooled = np.vstack(residuals[s])
        sizes = {0: 20, 1: 20, 2: 40}[s]
        estimate = np.sum(pooled**2) / (pooled.size * (sizes - 1) / sizes)
        assert estimate == pytest.approx(multiplier, rel=0.05), (pattern, s)


def test_a_dispersion_unit_reads_one_matrix_under_both_nulls(monkeypatch) -> None:
    calls = _counting_cv(monkeypatch)
    settings = wsp.dispersion_settings(replicates=2, permutations=PERMUTATIONS)
    record = wsp.evaluate_dispersion_unit(lsc.Unit("dispersion", "predictive_d5", 1), settings)
    assert record["status"] == "completed" and calls == []
    matrix, within = record["arms"]["matrix"], record["arms"]["matrix__within_stage"]
    assert matrix["pair_statistics"] == within["pair_statistics"]
    assert record["seeds"] == wsp.dispersion_seeds("predictive_d5", 1)
    assert set(record) >= {"pattern", "dimension", "rrpp_seed"} and record["dimension"] == 5
    assert all(not isinstance(value, list) for arm in record["arms"].values() for value in arm.values())


# ── 4.6 Pre-flight ───────────────────────────────────────────────────────────


def test_the_preflight_passes_and_the_run_waits_on_it(tmp_path, stand_in, plasmode_settings) -> None:
    with pytest.raises(SystemExit, match="preflight"):
        wsp.run_shard(tmp_path, plasmode_settings, 0, 1, export=stand_in["export"], limit=1)
    result = wsp.preflight(tmp_path, plasmode_settings, stand_in["export"])
    assert result["passed"], result
    assert [unit["unit"] for unit in result["units"]] == [wsp.unit_key("plasmode", "anchor", 0),
                                                          wsp.unit_key("plasmode", "anchor", 500)]
    assert wsp.run_shard(tmp_path, plasmode_settings, 0, 1, export=stand_in["export"], limit=1) == 1


def test_a_preflight_at_another_rank_fails(tmp_path, stand_in, plasmode_settings) -> None:
    settings = {**plasmode_settings, "anchor_rank": plasmode_settings["anchor_rank"] + 1}
    result = wsp.preflight(tmp_path, settings, stand_in["export"])
    assert not result["passed"]
    assert all(unit["statistics_differing"] or unit["p_values_differing"] for unit in result["units"])
    with pytest.raises(SystemExit, match="pre-flight failed"):
        wsp.run_shard(tmp_path, settings, 0, 1, export=stand_in["export"], limit=1)


# ── 4.7 sbatch ───────────────────────────────────────────────────────────────


@pytest.mark.parametrize("profile", ["simulation", "plasmode"])
def test_sbatch_template_parses_and_dry_runs(tmp_path, profile) -> None:
    script = REPO / "scripts" / "within_stage_permutation_array.sbatch"
    shell = shutil.which("bash") or "bash"
    subprocess.run([shell, "-n", str(script)], check=True)
    env = {"PATH": "/usr/bin:/bin", "WSP_DRY_RUN": "1", "WSP_PROFILE": profile, "WSP_OUT": str(tmp_path / "run"),
           "WSP_EXPORT": str(tmp_path / "export"), "WSP_REPO": str(REPO), "N_SHARDS": "100",
           "SLURM_ARRAY_TASK_ID": "3", "WSP_ARGS": "--limit 1"}
    out = subprocess.run([shell, str(script)], env=env, cwd=tmp_path, check=True, capture_output=True,
                         text=True).stdout
    words = out.split()
    assert words[0] == f"PYTHONPATH={REPO / 'src'}"
    export = ["--export-dir", str(tmp_path / "export")] if profile == "plasmode" else []
    assert words[1:] == ["python", str(REPO / "scripts" / "within_stage_permutation_check.py"), "run", "--profile",
                         profile, *export, "--out-dir", str(tmp_path / "run"), "--shard", "3/100", "--limit", "1"]
    assert not (tmp_path / "run").exists()
    relative = subprocess.run([shell, str(script)], env={**env, "WSP_OUT": "run"}, cwd=tmp_path,
                              capture_output=True, text=True)
    assert relative.returncode == 2 and "absolute" in relative.stderr
    unknown = subprocess.run([shell, str(script)], env={**env, "WSP_PROFILE": "plant"}, cwd=tmp_path,
                             capture_output=True, text=True)
    assert unknown.returncode == 2 and "WSP_PROFILE" in unknown.stderr


# ── 5.1 Merge ────────────────────────────────────────────────────────────────


def test_dispersion_merge_writes_every_table(tmp_path) -> None:
    settings = wsp.dispersion_settings(replicates=2, permutations=PERMUTATIONS)
    wsp.run_shard(tmp_path, settings, 0, 1)
    summary = wsp.merge(tmp_path, "dispersion")
    assert summary["complete_units"] == 20 and summary["failed_units"] == 0 and summary["reproduction_passed"]
    directory = tmp_path / "dispersion"
    rates = lsc._read_csv(directory / "rate_table.csv")
    paired = lsc._read_csv(directory / "paired_table.csv")
    assert len(rates) == 10 * 2 * 3 and len(paired) == 10 * 3
    assert {(row["arm"], row["reference_arm"]) for row in paired} == {("matrix__within_stage", "matrix")}
    for row in paired:
        if int(row["b"]) == 0 and int(row["c"]) == 0:
            assert float(row["p_loss"]) == float(row["p_gain"]) == float(row["p_two_sided"]) == 1.0


def test_plasmode_paired_rows(tmp_path, plasmode_settings, anchor_record, fresh_record) -> None:
    rows = wsp.paired_table([anchor_record, fresh_record], plasmode_settings)
    pairs = {(row["arm"], row["reference_arm"]) for row in rows}
    assert pairs == {
        *((f"{arm}__within_stage", arm) for arm in wsp.UNRESTRICTED_ARMS),
        *((arm, "pls_cv") for arm in wsp.PREDICTIVE_ARMS),
        *((f"{arm}__within_stage", "pls_cv__within_stage") for arm in wsp.PREDICTIVE_ARMS),
    }
    assert len(rows) == 2 * 3 * 10


# ── 5.2 Reproduction ─────────────────────────────────────────────────────────


def _committed(path: Path) -> list[dict]:
    return lsc._read_csv(path)


def test_weak_size_reproduction_reads_every_committed_row() -> None:
    settings = wsp.weak_size_settings(replicates=2, permutations=PERMUTATIONS)
    committed = _committed(wsp.PILOT_DIR / "weak_size" / "rate_table.csv")
    passed = wsp.weak_size_reproduction(committed, settings, committed)
    check = passed["checks"]["pilot_weak_size"]
    assert passed["passed"] and check["rows_compared"] == 24 and check["rows_matching"] == 24
    changed = [dict(row) for row in committed]
    target = next(row for row in changed if row["arm"] == "predictive" and row["statistic"] == "shape")
    target["rejections"] = str(int(target["rejections"]) + 1)
    failed = wsp.weak_size_reproduction(changed, settings, committed)
    assert not failed["passed"] and len(failed["checks"]["pilot_weak_size"]["differing"]) == 1


def test_plasmode_reproduction_reads_both_parents(plasmode_settings) -> None:
    pilot = _committed(wsp.PILOT_DIR / "plasmode" / "rate_table.csv")
    diagnosis = _committed(wsp.DIAGNOSIS_DIR / "rate_table.csv")
    ours = [row for row in pilot if row["arm"] in wsp.UNRESTRICTED_ARMS] + [
        row for row in diagnosis if row["arm"] in wsp.WITHIN_STAGE_ARMS
    ]
    passed = wsp.plasmode_reproduction(ours, plasmode_settings, pilot, diagnosis)
    assert passed["passed"]
    assert passed["checks"]["pilot_anchor"]["rows_compared"] == 12
    assert passed["checks"]["diagnosis_anchor"]["rows_compared"] == 12
    counts = {(row["arm"], row["statistic"]): int(row["rejections"]) for row in ours if row["cell"] == "anchor"}
    assert counts[("pls_cv", "delta")] == 21 and counts[("predictive_orthonormal__within_stage", "delta")] == 29
    changed = [dict(row) for row in ours]
    target = next(row for row in changed if row["cell"] == "anchor" and row["arm"] == "pls_cv__within_stage"
                  and row["statistic"] == "shape")
    target["count"] = str(int(target["count"]) - 1)
    failed = wsp.plasmode_reproduction(changed, plasmode_settings, pilot, diagnosis)
    assert not failed["passed"] and failed["checks"]["pilot_anchor"]["passed"]
    assert len(failed["checks"]["diagnosis_anchor"]["differing"]) == 1


# ── 5.3–5.6 The decision on synthetic tables ─────────────────────────────────

COUNTS = {"simulation": 100, "weak_size": 500, "plasmode": 500, "dispersion": 1000}


def _tables(rejections: dict | None = None, discordant: dict | None = None, *, drop: tuple | None = None,
            rate: float = 0.04) -> dict:
    plan = wsp.criteria_plan()
    readings = sorted({reading for rows in plan.values() for reading in rows})
    tables = {profile: {"rates": [], "paired": [], "cells": []} for profile in wsp.PROFILES}
    for profile, design, cell, statistic in readings:
        n = COUNTS[profile]
        arms = wsp.PROFILE_ARMS[profile]
        for arm in arms:
            key = (profile, design, cell, arm, statistic)
            if key == drop:
                continue
            k = (rejections or {}).get(key, round(rate * n))
            tables[profile]["rates"].append({"design": design, "cell": cell, "arm": arm, "statistic": statistic,
                                             "count": n, "rejections": k, "rate": k / n})
        for arm, reference in wsp.PAIRS:
            if arm not in arms or reference not in arms:
                continue
            b, c = (discordant or {}).get((profile, design, cell, statistic, arm, reference), (0, 0))
            tables[profile]["paired"].append({
                "design": design, "cell": cell, "statistic": statistic, "arm": arm, "reference_arm": reference,
                "n_paired": n, "rate_reference": rate, "rate_arm": rate, "difference": (c - b) / n, "b": b, "c": c,
            })
    return tables


def _decide(tables: dict, *, reproduction: dict | None = None, inconsistent: dict | None = None,
            checked: dict | None = None) -> dict:
    reproduction = reproduction or {profile: {"passed": True, "checks": {}} for profile in wsp.PROFILES}
    summaries = {profile: {"inconsistent_units": (inconsistent or {}).get(profile, 0)} for profile in wsp.PROFILES}
    return wsp.decision(tables, reproduction=reproduction, accounting_summaries=summaries,
                        checked={"passed": True} if checked is None else checked)


NULL, CAND = wsp.NULL_ARM, wsp.CANDIDATE_ARM
STRONG_TARGET = ("simulation", "n200", "orient_e0.25_phi0", "angle")
BAND_TARGET = ("simulation", "n200_delta1", "magnitude_e0.25", "delta")
BAND_OFF = ("simulation", "n200_delta1", "magnitude_e0.25", "angle")


def test_the_plan_has_the_declared_sizes() -> None:
    plan = wsp.criteria_plan()
    sizes = {name: len(rows) for name, rows in plan.items()}
    assert sizes == {"S1_simulation": 42, "S1_weak_size": 6, "S1_plasmode": 3, "S1_dispersion": 30, "S2": 54,
                     "S3_strong": 36, "S3_band": 72, "G": 18, "G_band": 36, "G_weak": 9}
    assert {reading[2] for reading in plan["S1_plasmode"]} == {"anchor_fresh"}


def test_a_null_that_holds_and_a_candidate_without_gain() -> None:
    result = _decide(_tables())
    assert result["verdict"] == "within_stage_validated"
    assert result["null"]["verdict"] == "validated" and result["candidate"]["verdict"] == "no_gain"
    assert result["candidate"]["bears_on_verdict"]
    assert result["consequence"] == wsp.CONSEQUENCES["within_stage_validated"]


def test_a_qualifying_candidate_proposes_a_paper_grade() -> None:
    tables = _tables(discordant={(*STRONG_TARGET, CAND, NULL): (0, 12)})
    result = _decide(tables)
    assert result["candidate"]["verdict"] == "qualifies"
    assert result["verdict"] == "propose_paper_grade"


def test_a_gain_that_is_the_nulls_does_not_count_for_the_candidate() -> None:
    tables = _tables(discordant={(*STRONG_TARGET, NULL, "pls_cv"): (0, 12)})
    result = _decide(tables)
    assert result["candidate"]["verdict"] == "no_gain"
    assert result["verdict"] == "within_stage_validated"
    gains = result["null"]["recorded_gains"]["strong"]
    assert gains["holds"]


@pytest.mark.parametrize(
    "scenario, verdict, failing",
    [
        ({"rejections": {("dispersion", "dispersion", "predictive_d2", "matrix__within_stage", "delta"): 150}},
         "size_excess", "size"),
        ({"rejections": {("plasmode", "plasmode", "anchor_fresh", NULL, "shape"): 80}}, "size_excess", "size"),
        ({"discordant": {(*BAND_TARGET, NULL, "pls_cv"): (12, 0)}}, "power_loss", "power"),
        ({"rejections": {("simulation", "n200", "magnitude_e0.25", NULL, "angle"): 40}}, "specificity_excess",
         "specificity"),
        ({"rejections": {(*BAND_OFF[:3], NULL, BAND_OFF[3]): 12}, "discordant": {(*BAND_OFF, NULL, "pls_cv"): (0, 12)}},
         "specificity_excess", "specificity"),
    ],
)
def test_each_failure_of_the_null(scenario, verdict, failing) -> None:
    tables = _tables(scenario.get("rejections"), scenario.get("discordant"))
    result = _decide(tables)
    assert result["verdict"] == verdict
    assert result["null"]["failing_criteria"][0] == failing
    assert result["consequence"] == wsp.CONSEQUENCES[verdict]


def test_restored_size_is_not_read_as_a_leak() -> None:
    tables = _tables(discordant={(*BAND_OFF, NULL, "pls_cv"): (0, 10)})
    result = _decide(tables)
    band = result["null"]["criteria"]["specificity"]["families"]["band"]
    assert band["holds"] and band["observed"] == 0
    assert result["verdict"] == "within_stage_validated"


def test_an_unavailable_reading_fails_its_criterion() -> None:
    tables = _tables(drop=("weak_size", "seaad_weak", "anchor", NULL, "delta"))
    result = _decide(tables)
    assert result["verdict"] == "size_excess"
    assert result["null"]["criteria"]["size"]["families"]["weak_size"]["n_unavailable"] == 1


@pytest.mark.parametrize(
    "scenario, verdict",
    [
        ({"rejections": {("simulation", "n200", "anchor", CAND, "delta"): 40}}, "size_excess"),
        ({"discordant": {(*BAND_TARGET, CAND, NULL): (12, 0)}}, "power_loss"),
        ({"rejections": {("simulation", "n200", "magnitude_e0.25", CAND, "shape"): 40}}, "specificity_excess"),
    ],
)
def test_each_reading_of_the_candidate(scenario, verdict) -> None:
    discordant = {(*STRONG_TARGET, CAND, NULL): (0, 12), **scenario.get("discordant", {})}
    result = _decide(_tables(scenario.get("rejections"), discordant))
    assert result["candidate"]["verdict"] == verdict
    assert result["verdict"] == "within_stage_validated"


def test_the_candidate_cannot_change_a_failed_null() -> None:
    tables = _tables(
        {("dispersion", "dispersion", "scores_d5", "matrix__within_stage", "shape"): 150},
        {(*STRONG_TARGET, CAND, NULL): (0, 12)},
    )
    result = _decide(tables)
    assert result["verdict"] == "size_excess"
    assert result["candidate"]["verdict"] == "qualifies" and not result["candidate"]["bears_on_verdict"]


@pytest.mark.parametrize(
    "kwargs",
    [
        {"reproduction": {**{p: {"passed": True, "checks": {}} for p in wsp.PROFILES},
                          "weak_size": {"passed": False, "checks": {"pilot_weak_size": {"passed": False}}}}},
        {"inconsistent": {"plasmode": 1}},
        {"checked": {"passed": False}},
    ],
)
def test_the_gate_records_every_criterion(kwargs) -> None:
    result = _decide(_tables(), **kwargs)
    assert result["verdict"] == "reproduction_failed" and not result["gate"]["passed"]
    assert result["null"]["verdict"] == "validated" and result["candidate"]["verdict"] == "no_gain"


def test_a_missing_preflight_fails_the_gate() -> None:
    reproduction = {profile: {"passed": True, "checks": {}} for profile in wsp.PROFILES}
    summaries = {profile: {"inconsistent_units": 0} for profile in wsp.PROFILES}
    result = wsp.decision(_tables(), reproduction=reproduction, accounting_summaries=summaries, checked=None)
    assert result["verdict"] == "reproduction_failed"


def test_descriptive_readings_never_change_the_verdict() -> None:
    plain = _decide(_tables())
    noisy = _decide(_tables({("simulation", "n200", "anchor", "predictive__within_stage", "delta"): 60,
                             ("dispersion", "dispersion", "predictive_d2", "matrix", "delta"): 300}))
    assert plain["verdict"] == noisy["verdict"] == "within_stage_validated"
    recorded = noisy["descriptive"]["recorded_candidates"]["predictive__within_stage"]
    assert recorded["verdict"] == "size_excess"


def test_the_verdicts_and_consequences_are_declared() -> None:
    assert set(wsp.CONSEQUENCES) == set(wsp.VERDICTS)
    assert "500 x 999" in wsp.CONSEQUENCES["propose_paper_grade"]


def _pilot_tables() -> dict:
    tables = {}
    for profile in psc.PROFILES:
        directory = wsp.PILOT_DIR / profile
        tables[profile] = {"rates": lsc._read_csv(directory / "rate_table.csv"),
                           "paired": lsc._read_csv(directory / "paired_table.csv")}
    return tables


@pytest.mark.parametrize("arm", list(psc.CANDIDATE_ARMS))
def test_the_candidate_rule_is_the_pilots_on_its_committed_tables(arm) -> None:
    tables = _pilot_tables()
    rates, paired = psc._index_rates(tables), psc._index_paired(tables)
    plan = psc.criteria_plan()
    ours = wsp.judge_candidate(arm, "pls_cv", plan, rates, paired)
    theirs = psc.judge_arm(arm, plan, rates, paired)
    assert json.loads(json.dumps(ours, default=str)) == json.loads(json.dumps(theirs, default=str))
    assert ours["verdict"] == "size_excess"


def test_judge_band_excess_counts_only_above_the_bound() -> None:
    def row(hits: int, c: int) -> dict:
        return {"n_available": 100, "n_rejected": hits, **plsc.paired_p_values(0, c)}

    within = wsp.judge_band_excess([row(4, 10)] + [row(4, 0)] * 71)
    assert within["holds"] and within["observed"] == 0
    above = wsp.judge_band_excess([row(12, 10)] + [row(4, 0)] * 71)
    assert above["observed"] == 1 and above["holds"] == (1 <= psc.allowed_count(72))
    ceiling = wsp.judge_band_excess([row(12, 12)] + [row(4, 0)] * 71)
    assert not ceiling["holds"] and ceiling["n_below_ceiling"] == 1
    assert not math.isnan(wsp.dispersion_centroid_scale())
