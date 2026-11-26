"""``scripts/candidate_simulation_check.py``: profiles, units, seeds and parents; the arms on one dataset; resume,
merge and reproduction; the decision by design family.

Real units run at fast CV knobs and few permutations; the decision runs on synthetic tables over the real
settings. No real data are read.
"""

from __future__ import annotations

import json
import os
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


csc = _script("candidate_simulation_check")
ubc, lsc, ses, ags, plsc, grade = csc.ubc, csc.lsc, csc.ses, csc.ags, csc.plsc, csc.grade

REPO = Path(__file__).resolve().parents[1]
FAST_CV = {"cv1_splits": 3, "cv2_splits": 4, "n_repeats": 1, "max_components": 6, "random_state": 1203}
PERMS = 19
STATISTICS = csc.STATISTICS


@pytest.fixture(scope="module")
def full():
    return {profile: csc.profile_settings(profile) for profile in csc.PROFILES}


def _small(profile="committed", designs=None, reps=2, control=2, perms=PERMS):
    return csc.profile_settings(profile, replicates=reps, control_replicates=control, permutations=perms,
                                designs=designs)


def _fast(settings):
    for group in settings["groups"]:
        if group["kind"] == "stage_span":
            group["settings"]["units"]["pls_cv_params"] = dict(FAST_CV)
    return settings


# ── 2.1 Settings, units, resume and parents ──────────────────────────────────


def test_unit_counts_per_profile(full) -> None:
    for profile, settings in full.items():
        units = csc.profile_units(settings)
        assert len(units) == csc.EXPECTED_UNITS[profile], profile
        assert len({(u.design, u.cell, u.replicate) for u in units}) == len(units)
        sparse = [u for u in units if csc.arm_applies(settings, u.design, csc.SPARSE_ARM)]
        assert len(sparse) == (9_200 if profile == "committed" else 0), profile
    assert sum(csc.EXPECTED_UNITS.values()) == 58_400


def test_every_profile_binds_its_parents(full) -> None:
    for profile, settings in full.items():
        for name in csc.PROFILE_PARENTS[profile]:
            directory, files, _ = csc.PARENTS[name]
            for file_name in files:
                assert f"results/{directory}/{file_name}" in settings["parent_files"], (profile, name, file_name)
    surgery = full["surgery"]["parent_files"]
    assert any(path.endswith("design_guidance_core_geometry_study.json") for path in surgery)
    adni = full["adni_design"]["parent_files"]
    assert any("adni-real-noise-anchor" in path for path in adni)


def test_parse_arms() -> None:
    assert csc.parse_arms("all") == csc.ARMS
    assert csc.parse_arms("weighted,pls_cv") == ("pls_cv", "weighted")
    with pytest.raises(SystemExit, match="--arms"):
        csc.parse_arms("pls_cv,concat")


def _fake_records(unit, settings, arms):
    base = {**csc.unit_fields(unit, settings), "settings_hash": lsc.settings_hash(settings)}
    return [{**base, "arm": arm, "status": "completed", "rank": 2, "p_values": {s: 0.5 for s in STATISTICS},
             "pair_statistics": {s: 1.0 for s in STATISTICS}, "seconds": {"arm": 0.0, "call_share": 0.0}}
            for arm in arms if csc.arm_applies(settings, unit.design, arm)]


def _lines(path):
    return [json.loads(line) for line in Path(path).read_text().splitlines() if line.strip()]


def test_resume_skips_complete_unit_arm_keys(tmp_path, monkeypatch) -> None:
    settings = _small(designs=["n80"])
    calls = []

    def fake(unit, settings_, arms):
        calls.append((unit.cell, unit.replicate, tuple(arms)))
        return _fake_records(unit, settings_, arms)

    monkeypatch.setattr(csc, "evaluate_unit", fake)
    arms = ("pls_cv", "predictive_orthonormal")
    assert csc.run_shard(tmp_path, settings, 0, 1, arms=arms, limit=3) == 3
    assert csc.run_shard(tmp_path, settings, 0, 1, arms=arms) == len(csc.profile_units(settings)) - 3
    assert csc.run_shard(tmp_path, settings, 0, 1, arms=arms) == 0
    assert len(calls) == len(csc.profile_units(settings))


def test_shards_are_round_robin(tmp_path, monkeypatch) -> None:
    settings = _small(designs=["n80"])
    seen = []
    monkeypatch.setattr(csc, "evaluate_unit", lambda unit, s, arms: seen.append(unit) or _fake_records(unit, s, arms))
    csc.run_shard(tmp_path, settings, 1, 3, arms=("pls_cv",))
    units = csc.profile_units(settings)
    assert seen == [unit for position, unit in enumerate(units) if position % 3 == 1]


def test_an_appended_arm_keeps_the_existing_records(tmp_path, monkeypatch) -> None:
    settings = _small(designs=["n200"])
    monkeypatch.setattr(csc, "evaluate_unit", _fake_records)
    csc.run_shard(tmp_path, settings, 0, 1, arms=("pls_cv", "predictive_orthonormal", "weighted"))
    path = lsc.shard_path(tmp_path / "committed", 0, 1)
    before = _lines(path)
    asked = []

    def only_sparse(unit, settings_, arms):
        asked.append(tuple(arms))
        return _fake_records(unit, settings_, arms)

    monkeypatch.setattr(csc, "evaluate_unit", only_sparse)
    assert csc.run_shard(tmp_path, settings, 0, 1, arms=("sparse_block",)) == len(csc.profile_units(settings))
    assert set(asked) == {("sparse_block",)}
    after = _lines(path)
    assert {csc.record_key(r) for r in before} < {csc.record_key(r) for r in after}
    kept = {csc.record_key(r): r for r in after}
    assert all(kept[csc.record_key(r)] == r for r in before)
    stored = json.loads((tmp_path / "committed" / "settings.json").read_text())
    assert stored["arm_selections"] == [["pls_cv", "predictive_orthonormal", "weighted"], ["sparse_block"]]
    assert stored["hash"] == lsc.settings_hash(settings)
    assert lsc.load_settings(tmp_path / "committed") == json.loads(json.dumps(settings))


def test_a_resume_under_other_settings_is_refused_naming_the_keys(tmp_path, monkeypatch) -> None:
    monkeypatch.setattr(csc, "evaluate_unit", _fake_records)
    csc.run_shard(tmp_path, _small(designs=["n80"]), 0, 1, arms=("pls_cv",), limit=1)
    with pytest.raises(SystemExit, match=r"differing: .*permutations"):
        csc.run_shard(tmp_path, _small(designs=["n80"], perms=29), 0, 1, arms=("pls_cv",), limit=1)


def test_merge_refuses_records_written_under_other_settings(tmp_path, monkeypatch) -> None:
    settings = _small(designs=["n80"])
    monkeypatch.setattr(csc, "evaluate_unit", lambda unit, s, arms: [
        {**record, "settings_hash": "other"} for record in _fake_records(unit, s, arms)])
    csc.run_shard(tmp_path, settings, 0, 1, arms=("pls_cv",), limit=1)
    with pytest.raises(SystemExit, match="other settings"):
        csc.merge(tmp_path, "committed")


def test_a_changed_parent_is_refused_by_name(tmp_path, monkeypatch) -> None:
    monkeypatch.setattr(csc, "evaluate_unit", _fake_records)
    settings = _small(designs=["n80"])
    stored = "results/univariate-baseline-2026-11-23/stage_span/rate_table.csv"
    settings["parent_files"][stored] = "0" * 64
    with pytest.raises(SystemExit, match=r"stage_span/rate_table\.csv has changed"):
        csc.run_shard(tmp_path, settings, 0, 1, arms=("pls_cv",))
    assert not (tmp_path / "committed").exists()


def test_merge_and_decide_refuse_a_changed_parent(tmp_path, monkeypatch) -> None:
    monkeypatch.setattr(csc, "evaluate_unit", _fake_records)
    csc.run_shard(tmp_path, _small(designs=["n80"]), 0, 1, arms=("pls_cv",))
    csc.merge(tmp_path, "committed")
    real = csc.ssc.sha256_file
    monkeypatch.setattr(csc.ssc, "sha256_file", lambda path: "0" * 64 if Path(path).name == "rate_table.csv"
                        else real(path))
    with pytest.raises(SystemExit, match="rate_table.csv has changed"):
        csc.merge(tmp_path, "committed")
    with pytest.raises(SystemExit, match="rate_table.csv has changed"):
        csc.decide(tmp_path)


def test_a_parent_its_builder_no_longer_rebuilds_is_refused(monkeypatch) -> None:
    real = ubc.stage_span_settings

    def drifted(*args, **kwargs):
        settings = real(*args, **kwargs)
        settings["units"]["pls_cv_params"] = {**settings["units"]["pls_cv_params"], "n_repeats": 3}
        return settings

    monkeypatch.setattr(csc.ubc, "stage_span_settings", drifted)
    with pytest.raises(SystemExit, match=r"stage_span/settings\.json differs .*pls_cv_params\.n_repeats"):
        csc.profile_settings("committed", designs=["n80"])


def test_unknown_designs_are_refused() -> None:
    with pytest.raises(SystemExit, match="Unknown conditions design"):
        csc.profile_settings("conditions", designs=["n200"])


# ── 2.2 The committed profile is the parents' units ─────────────────────────


def _committed_units(name):
    settings = csc.committed_settings(name)
    return settings.get("unit_settings", settings)


def test_committed_designs_are_their_parents_units(full) -> None:
    settings = full["committed"]
    sources = {**{design: entry[0] for design, entry in csc.COMMITTED_DESIGNS.items()}, "band_n320": "scaling_band"}
    for design, parent in sources.items():
        group = csc.group_of(settings, design)["settings"]
        theirs = _committed_units(parent)
        assert csc.design_cells(settings, design) == ubc._cell_names(theirs), design
        for cell in csc.design_cells(settings, design):
            for replicate in (0, 137, 199):
                unit = lsc.Unit(design, cell, replicate)
                mine = csc.unit_fields(unit, settings)
                parent_fields = ubc.unit_fields(unit, theirs)
                assert (mine["generator_seed"], mine["rrpp_seed"]) == (
                    parent_fields["generator_seed"], parent_fields["rrpp_seed"]), (design, cell)
                assert ubc.stage_span_generator_params(unit, group) == ubc.stage_span_generator_params(unit, theirs)
                assert ubc.stage_span_params(group, seed=1) == ubc.stage_span_params(theirs, seed=1)


def test_committed_replicates_and_cells(full) -> None:
    counts = csc.expected_counts(full["committed"])
    assert counts[("n80", "anchor")] == counts[("n320_delta1", "translation")] == 500
    assert counts[("n500", "shape_e0.5")] == counts[("band_n320", "magnitude_e2")] == 200
    assert len(csc.design_cells(full["committed"], "band_n320")) == 11
    assert len(csc.design_cells(full["committed"], "n320_delta1")) == 13
    assert full["committed"]["design_index"]["band_n320"]["delta"] == 0.05


# ── 2.3 and 2.4 New conditions and the dispersion designs ────────────────────


def test_condition_designs(full) -> None:
    settings = full["conditions"]
    index = settings["design_index"]
    assert list(index) == list(csc.CONDITION_DESIGNS)
    assert [index[d]["n_stages"] for d in index] == [4, 4, 3, 3, 3]
    assert [index[d]["delta"] for d in index] == [2.0, 0.5, 2.0, 0.5, 2.0]
    assert index["n200_share0.8"]["share"] == 0.8
    for design in index:
        assert len(csc.design_cells(settings, design)) == 13
        group = csc.group_of(settings, design)["settings"]
        assert group["units"]["designs"][design]["layers"] is None
        params = ubc.stage_span_generator_params(lsc.Unit(design, "anchor", 0), group)
        assert params.span_layers is None and params.construction == "stage_span"
        assert params.n_stages == index[design]["n_stages"]
    for design in ("seaad_table", "seaad_table_delta0.5"):
        params = ubc.stage_span_generator_params(lsc.Unit(design, "anchor", 0), csc.group_of(settings, design)[
            "settings"])
        assert params.group_stage_sizes == csc.SEAAD_TABLE
    counts = csc.expected_counts(settings)
    assert counts[("seaad_table", "anchor")] == 500 and counts[("seaad_table", "shape_e0.1")] == 200


def test_condition_designs_measure_three_blocks() -> None:
    from motco.simulations.reference import load_reference
    from motco.simulations.semisynthetic import generate_semisynthetic_trajectory

    settings = _small("conditions", designs=["seaad_table"], reps=1, control=1)
    group = csc.group_of(settings, "seaad_table")["settings"]
    dataset = generate_semisynthetic_trajectory(
        ubc.stage_span_generator_params(lsc.Unit("seaad_table", "anchor", 0), group), reference=load_reference())
    _, params = ubc.unit_inputs(lsc.Unit("seaad_table", "anchor", 0), group, {"rrpp_seed": 1})
    assert "layers" not in params.integration_params
    sizes = dataset.metadata.groupby(["group", "stage"]).size().unstack().to_numpy().tolist()
    assert sizes == [list(row) for row in csc.SEAAD_TABLE]


def test_cells_of_a_design_and_replicate_share_their_null_dataset(full) -> None:
    for profile in ("conditions", "dispersion"):
        settings = full[profile]
        for design in settings["design_index"]:
            seeds = {csc.unit_fields(lsc.Unit(design, cell, 3), settings)["generator_seed"]
                     for cell in csc.design_cells(settings, design)}
            assert seeds == {lsc.generator_seed(csc.BASE_SEED, design, 3)}
            rrpp = {csc.unit_fields(lsc.Unit(design, cell, 3), settings)["rrpp_seed"]
                    for cell in csc.design_cells(settings, design)}
            assert len(rrpp) == len(csc.design_cells(settings, design))


def _committed_design_seeds() -> tuple[set[str], set[int]]:
    """Every design name and every generator seed (replicates 0-499) of every committed settings file.

    This run's own committed aggregates (``results/candidate-simulation-*``) are skipped: they hold the new designs by
    construction, and freshness is asked against the runs before it.
    """

    names: set[str] = set()
    seeds: set[int] = set()

    def walk(node):
        if isinstance(node, dict):
            designs = node.get("designs")
            if "base_seed" in node and isinstance(designs, (dict, list)):
                for name in designs:
                    names.add(str(name))
                    seeds.update(lsc.generator_seed(int(node["base_seed"]), str(name), r) for r in range(500))
            for value in node.values():
                walk(value)
        elif isinstance(node, list):
            for value in node:
                walk(value)

    for path in sorted((REPO / "results").rglob("settings.json")):
        if path.relative_to(REPO / "results").parts[0].startswith("candidate-simulation-"):
            continue
        payload = json.loads(path.read_text(encoding="utf-8"))
        walk(payload)
        text = path.read_text(encoding="utf-8")
        names.update(name for name in [*csc.CONDITION_DESIGNS, *(e["design"] for e in csc.dispersion_entries())]
                     if f'"{name}"' in text)
    return names, seeds


def test_new_names_and_datasets_are_fresh() -> None:
    names, seeds = _committed_design_seeds()
    new = [*csc.CONDITION_DESIGNS, *(entry["design"] for entry in csc.dispersion_entries())]
    assert not set(new) & names
    assert {"n80", "n200", "n320_delta1", "band_n320"} <= names  # the walk sees the committed designs
    mine = {lsc.generator_seed(csc.BASE_SEED, name, r) for name in new for r in range(500)}
    assert not mine & seeds


def test_dispersion_patterns_and_designs(full) -> None:
    assert csc.DISPERSION_PATTERNS["last_narrow"] == (1.058, 1.130, 0.896)
    assert csc.DISPERSION_PATTERNS["last_wide"] == (0.897, 0.960, 1.066)
    assert csc.DISPERSION_PATTERNS["last_narrow_strong"] == (1.120, 1.278, 0.802)
    assert csc.DISPERSION_PATTERNS["last_wide_strong"] == (0.805, 0.921, 1.137)
    settings = full["dispersion"]
    index = settings["design_index"]
    assert len(index) == 10
    assert sum(1 for meta in index.values() if meta["reference_design"] == "n200") == 2
    for design, meta in index.items():
        assert csc.design_cells(settings, design) == ["anchor", "translation"]
        assert csc.expected_counts(settings)[(design, "anchor")] == 500
        mine = csc.group_of(settings, design)["settings"]["units"]["designs"][design]["generator"]
        reference = csc.reference_generator(meta["reference_design"])
        assert {k: v for k, v in mine.items() if k != "stage_dispersion_scale"} == reference
        assert mine["stage_dispersion_scale"] == list(csc.DISPERSION_PATTERNS[meta["pattern"]])
    committed_n200 = csc.committed_settings("univariate_stage_span")["units"]["designs"]["n200"]["generator"]
    assert csc.reference_generator("n200") == committed_n200
    assert csc.reference_generator("seaad_table") == csc.group_of(full["conditions"], "seaad_table")[
        "settings"]["units"]["designs"]["seaad_table"]["generator"]


def test_no_dispersion_unit_plants_a_trajectory_difference(full) -> None:
    settings = full["dispersion"]
    for unit in csc.profile_units(settings)[::97]:
        fields = csc.unit_fields(unit, settings)
        assert fields["trajectory_mode"] in ("none", "translation")


def test_a_dispersion_dataset_is_stretched_by_stage() -> None:
    from motco.simulations.reference import load_reference
    from motco.simulations.semisynthetic import generate_semisynthetic_trajectory

    settings = _small("dispersion", designs=["n200_last_narrow_strong"], control=1)
    group = csc.group_of(settings, "n200_last_narrow_strong")["settings"]
    params = ubc.stage_span_generator_params(lsc.Unit("n200_last_narrow_strong", "anchor", 0), group)
    dataset = generate_semisynthetic_trajectory(params, reference=load_reference())
    assert dataset.truth["stage_dispersion"]["scale"] == list(csc.DISPERSION_PATTERNS["last_narrow_strong"])


# ── 2.5 The surgery and ADNI-design profiles ─────────────────────────────────


def test_surgery_units_are_the_parents(full) -> None:
    settings = full["surgery"]
    theirs = csc.committed_settings("univariate_surgery")
    group = csc.group_of(settings, "n80_three_blocks")["settings"]
    for design in csc.SURGERY_COLUMNS:
        for cell in csc.design_cells(settings, design):
            for replicate in (0, 499):
                unit = lsc.Unit(design, cell, replicate)
                assert ubc.surgery_params(unit, group) == ubc.surgery_params(unit, theirs)
                mine = csc.unit_fields(unit, settings)
                parent = ubc.unit_fields(unit, theirs)
                assert (mine["generator_seed"], mine["rrpp_seed"]) == (parent["generator_seed"], parent["rrpp_seed"])
    assert csc.expected_counts(settings)[("n200_three_blocks", "shape_e1")] == 500


def test_adni_units_are_the_gain_runs(full) -> None:
    settings = full["adni_design"]
    gain = csc.group_of(settings, csc.ADNI_DESIGN)["settings"]
    committed = csc.committed_settings("adni_gain_simulation")
    assert {k: v for k, v in gain.items() if k not in ("arms", "arm_specs")} == {
        k: v for k, v in committed.items() if k not in ("arms", "arm_specs")}
    assert gain["arms"] == ["pls_cv", "predictive_orthonormal"]
    assert gain["arm_specs"] == {arm: committed["arm_specs"][arm] for arm in gain["arms"]}
    mine = {(u.design, u.cell, u.replicate) for u in csc.profile_units(settings)}
    theirs = {(u.point, u.cell, u.replicate) for u in ags.enumerate_units(committed)}
    assert mine == theirs


# ── 2.6 The sbatch array ─────────────────────────────────────────────────────


SBATCH = REPO / "scripts" / "candidate_simulation_array.sbatch"


def test_the_sbatch_array_parses_and_prints_its_command(tmp_path) -> None:
    subprocess.run(["bash", "-n", str(SBATCH)], check=True)
    env = {**os.environ, "CSC_PROFILE": "committed", "CSC_ARMS": "pls_cv,predictive_orthonormal,weighted",
           "CSC_OUT": str(tmp_path), "N_SHARDS": "400", "SLURM_ARRAY_TASK_ID": "7", "CSC_DRY_RUN": "1",
           "SLURM_SUBMIT_DIR": str(REPO)}
    out = subprocess.run(["bash", str(SBATCH)], env=env, capture_output=True, text=True, check=True).stdout
    assert f"PYTHONPATH={REPO}/src" in out
    assert "candidate_simulation_check.py run --profile committed --arms pls_cv\\,predictive_orthonormal\\,weighted" \
        in out or "--arms pls_cv,predictive_orthonormal,weighted" in out
    assert "--shard 7/400" in out and str(tmp_path) in out
    bad = subprocess.run(["bash", str(SBATCH)], env={**env, "CSC_PROFILE": "plants"}, capture_output=True, text=True)
    assert bad.returncode == 2 and "CSC_PROFILE must be" in bad.stderr
    relative = subprocess.run(["bash", str(SBATCH)], env={**env, "CSC_OUT": "run"}, capture_output=True, text=True)
    assert relative.returncode == 2
    assert "never inside --export" in SBATCH.read_text()


# ── 3.1 and 3.2 The unit (real fits, fast knobs) ─────────────────────────────


def _counting_cv():
    import motco.simulations.evaluation as evaluation

    real = evaluation.plsda_doubleCV
    calls = []

    def counting(*args, **kwargs):
        calls.append(1)
        return real(*args, **kwargs)

    return mock.patch.object(evaluation, "plsda_doubleCV", counting), calls


@pytest.fixture(scope="module")
def unit_records():
    settings = _fast(_small(designs=["n200", "n80"], reps=1, control=1))
    unit = lsc.Unit("n200", "orient_e0.5_phi0", 0)
    patch, calls = _counting_cv()
    with patch:
        records = csc.evaluate_unit(unit, settings, csc.ARMS)
    all_arms = len(calls)
    patch, calls = _counting_cv()
    with patch:
        shared = csc.evaluate_unit(unit, settings, ("pls_cv", "predictive_orthonormal"))
    shared_calls = len(calls)
    patch, calls = _counting_cv()
    with patch:
        outside = csc.evaluate_unit(lsc.Unit("n80", "anchor", 0), settings, csc.ARMS)
    outside_calls = len(calls)
    parent = ubc.stage_span_settings(designs=["n200"], replicates=1, control_replicates=1, permutations=PERMS)
    parent["units"]["pls_cv_params"] = dict(FAST_CV)
    return {"settings": settings, "unit": unit, "records": {r["arm"]: r for r in records}, "all_calls": all_arms,
            "shared": {r["arm"]: r for r in shared}, "shared_calls": shared_calls,
            "outside": {r["arm"]: r for r in outside}, "outside_calls": outside_calls,
            "parent": ubc.evaluate_unit(unit, parent)}


def test_every_arm_reads_one_dataset_at_the_units_seed(unit_records) -> None:
    records, unit, settings = unit_records["records"], unit_records["unit"], unit_records["settings"]
    assert set(records) == set(csc.ARMS)
    seeds = {(r["generator_seed"], r["rrpp_seed"]) for r in records.values()}
    assert seeds == {(lsc.generator_seed(csc.BASE_SEED, "n200", 0),
                      lsc.rrpp_seed(csc.BASE_SEED, "n200", unit.cell, 0))}
    for record in records.values():
        assert record["status"] == "completed", record.get("error_message")
        assert record["settings_hash"] == lsc.settings_hash(settings)


def test_the_pls_cv_record_equals_the_parent_drivers(unit_records) -> None:
    mine, parent = unit_records["records"]["pls_cv"], unit_records["parent"]
    assert parent["status"] == "completed"
    production = parent["production"]
    assert mine["rank"] == production["rank"]
    assert mine["cv_mean_auroc"] == production["cv_mean_auroc"]
    assert mine["p_values"] == production["p_values"]
    assert mine["pair_statistics"] == production["pair_statistics"]
    assert (mine["generator_seed"], mine["rrpp_seed"]) == (parent["generator_seed"], parent["rrpp_seed"])


def test_the_predictive_arm_makes_no_cross_validated_fit(unit_records) -> None:
    assert unit_records["shared_calls"] == 1
    assert unit_records["all_calls"] == 3  # production, weighted and sparse_block: one CV each
    predictive, production = unit_records["records"]["predictive_orthonormal"], unit_records["records"]["pls_cv"]
    assert predictive["rank"] == production["rank"]
    assert predictive["cv_mean_auroc"] == production["cv_mean_auroc"]
    assert predictive["consistency"]["rank_equal"] and predictive["consistency"]["auroc_equal"]
    assert predictive["measured_dimension"] == min(production["rank"], 2)
    assert unit_records["shared"]["predictive_orthonormal"]["p_values"] == predictive["p_values"]


def test_own_fit_arms_change_only_their_integration_params(unit_records) -> None:
    sparse, weighted = unit_records["records"]["sparse_block"], unit_records["records"]["weighted"]
    assert sparse["keep_fraction"] in plsc.default_ladder()
    assert sparse["support"] and set(sparse["support"][0]) == {"methylation", "expression", "proteomics"}
    assert weighted["keep_fraction"] is None and weighted["support"] is None


def test_a_design_outside_the_subset_makes_no_sparse_fit(unit_records) -> None:
    assert set(unit_records["outside"]) == {"pls_cv", "predictive_orthonormal", "weighted"}
    assert unit_records["outside_calls"] == 2


def test_a_failed_arm_is_recorded_and_the_other_arms_and_the_next_unit_run(tmp_path, monkeypatch) -> None:
    settings = _fast(_small(designs=["n80"], reps=1, control=1))
    real = csc._own_fit_params

    def broken(params, arm, settings_):
        if arm == "weighted":
            raise RuntimeError("patched weighted failure")
        return real(params, arm, settings_)

    monkeypatch.setattr(csc, "_own_fit_params", broken)
    assert csc.run_shard(tmp_path, settings, 0, 7, arms=("pls_cv", "weighted"), limit=2) == 2
    records = _lines(lsc.shard_path(tmp_path / "committed", 0, 7))
    assert len(records) == 4
    for record in records:
        if record["arm"] == "weighted":
            assert record["status"] == "failed"
            assert (record["error_type"], record["error_message"]) == ("RuntimeError", "patched weighted failure")
        else:
            assert record["status"] == "completed"


@pytest.fixture(scope="module")
def dispersion_records():
    settings = _fast(_small("dispersion", designs=["seaad_table_last_narrow_strong"], control=1))
    return csc.evaluate_unit(lsc.Unit("seaad_table_last_narrow_strong", "anchor", 0), settings, csc.ARMS)


def test_recorded_arms_share_statistics_and_only_p_values_differ(dispersion_records) -> None:
    records = {r["arm"]: r for r in dispersion_records}
    assert set(records) == {"pls_cv", "predictive_orthonormal", "weighted"}
    for arm in ("pls_cv", "predictive_orthonormal"):
        within = records[arm]["within_stage"]
        assert within["same_statistics"] is True
        assert within["pair_statistics"] == records[arm]["pair_statistics"]
        assert set(within["p_values"]) == set(STATISTICS)
    for record in records.values():
        dispersion = record["stage_dispersion"]
        assert {"stage_0", "stage_1", "stage_2", "max_min_ratio"} <= set(dispersion)
    assert "within_stage" not in records["weighted"]


# ── 3.3 Merge (synthetic records) ────────────────────────────────────────────


def _records_for(settings, arms, p_of, status_of=None):
    out = []
    for unit in csc.profile_units(settings):
        for record in _fake_records(unit, settings, arms):
            status = "completed" if status_of is None else status_of(unit, record["arm"])
            record["status"] = status
            if status == "failed":
                record.update(error_type="RuntimeError", error_message="boom")
            record["p_values"] = {s: p_of(unit, record["arm"], s) for s in STATISTICS}
            out.append(record)
    return out


def _write_shard(directory, records, name="shard_000_of_001.jsonl"):
    path = Path(directory) / "records" / name
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("".join(lsc._line(r) for r in records))


def test_merge_counts_failures_handles_duplicates_and_marks_arms_not_run(tmp_path) -> None:
    settings = _small(designs=["n80"], reps=2, control=2)
    directory = tmp_path / "committed"
    csc.ensure_profile_settings(directory, settings, ("pls_cv", "predictive_orthonormal"))

    def p_of(unit, arm, statistic):
        return 0.01 if (arm == "predictive_orthonormal" and unit.replicate == 0) else 0.5

    def status_of(unit, arm):
        return "failed" if (arm == "pls_cv" and unit.cell == "anchor" and unit.replicate == 1) else "completed"

    records = _records_for(settings, ("pls_cv", "predictive_orthonormal"), p_of, status_of)
    _write_shard(directory, records)
    _write_shard(directory, records[:5], "shard_001_of_002.jsonl")  # five duplicates from an old shard count
    summary = csc.merge(tmp_path, "committed")
    assert summary["duplicate_records"] == 5
    accounting = json.loads((directory / "accounting.json").read_text())
    assert accounting["failures_by_cell"] == {"n80/anchor/pls_cv": 1}
    assert accounting["arms"]["pls_cv"]["failed"] == 1 and accounting["arms"]["weighted"]["missing"] == 0
    rates = {(r["cell"], r["arm"], r["statistic"]): r for r in ssc_read(directory / "rate_table.csv")}
    assert len(rates) == 13 * 4 * 3
    anchor = rates[("anchor", "pls_cv", "delta")]
    assert (anchor["count"], anchor["rejections"], anchor["failures"], anchor["status"]) == ("1", "0", "1",
                                                                                              "complete")
    predictive = rates[("magnitude_e0.05", "predictive_orthonormal", "angle")]
    assert (predictive["count"], predictive["rejections"], predictive["status"]) == ("2", "1", "complete")
    assert rates[("anchor", "weighted", "delta")]["status"] == "not_run"
    assert rates[("anchor", "sparse_block", "delta")]["status"] == "not_run"
    assert rates[("anchor", "sparse_block", "delta")]["expected"] == "0"
    paired = {(r["cell"], r["statistic"], r["arm"]): r for r in ssc_read(directory / "paired_table.csv")}
    row = paired[("magnitude_e0.05", "angle", "predictive_orthonormal")]
    assert (row["n_paired"], row["b"], row["c"]) == ("2", "0", "1")
    assert paired[("anchor", "delta", "predictive_orthonormal")]["n_paired"] == "1"
    assert paired[("anchor", "delta", "weighted")]["n_paired"] == "0"
    assert (directory / "cell_table.csv").exists() and (directory / "rank_table.csv").exists()
    reproduction = json.loads((directory / "reproduction.json").read_text())
    assert reproduction["passed"] is False and reproduction["rows_compared"] == 39


def ssc_read(path):
    return csc.ssc._read_csv(path)


def test_a_partially_run_arm_is_marked_partial(tmp_path) -> None:
    settings = _small(designs=["n80"], reps=2, control=2)
    directory = tmp_path / "committed"
    csc.ensure_profile_settings(directory, settings, ("pls_cv",))
    records = [r for r in _records_for(settings, ("pls_cv",), lambda u, a, s: 0.5) if r["replicate"] == 0]
    _write_shard(directory, records)
    csc.merge(tmp_path, "committed")
    statuses = {r["status"] for r in ssc_read(directory / "rate_table.csv") if r["arm"] == "pls_cv"}
    assert statuses == {"partial"}
    accounting = json.loads((directory / "accounting.json").read_text())
    assert accounting["arms"]["pls_cv"]["missing"] == 13


# ── 3.4 Reproduction ─────────────────────────────────────────────────────────


@pytest.mark.parametrize(("profile", "rows"), [("committed", 228), ("surgery", 66), ("adni_design", 60)])
def test_the_committed_tables_resolve_every_compared_row(full, profile, rows) -> None:
    result = csc.reproduction([], full[profile])
    assert result["unresolved_rows"] == []
    assert result["rows_compared"] == rows
    assert result["passed"] is False and len(result["differing"]) == rows


def _rates_from_committed(settings):
    """This driver's rate rows carrying the committed counts (what a faithful rerun would merge)."""

    rows = []
    for entry in settings["reproduction"]:
        for row in csc._committed_rows(entry["parent"]):
            if entry["table"] == "adni":
                if row["point"] == entry["design"] and row["arm"] in entry["arms"]:
                    rows.append({"design": row["point"], "cell": row["cell"], "arm": row["arm"],
                                 "statistic": row["statistic"], "count": row["count"],
                                 "rejections": row["rejections"]})
            elif row["design"] == entry["design"] and row["arm"] in STATISTICS:
                rows.append({"design": row["design"], "cell": row["cell"], "arm": "pls_cv", "statistic": row["arm"],
                             "count": row["count"], "rejections": row["rejections"]})
    return rows


@pytest.mark.parametrize("profile", ["committed", "surgery", "adni_design"])
def test_a_faithful_rerun_passes_and_a_changed_count_fails(full, profile) -> None:
    rates = _rates_from_committed(full[profile])
    assert csc.reproduction(rates, full[profile])["passed"] is True
    changed = [dict(row) for row in rates]
    changed[3]["rejections"] = str(int(changed[3]["rejections"]) + 1)
    result = csc.reproduction(changed, full[profile])
    assert result["passed"] is False
    assert len(result["differing"]) == 1
    differing = result["differing"][0]
    assert differing["reproduced"]["rejections"] == differing["committed"]["rejections"] + 1


# ── 4. The decision (synthetic tables over the real settings) ────────────────


def _table_blocks(full, *, rejected=None, discordant=None, absent=(), skip=(), passed=True):
    """Synthetic merged blocks: ``rejected(...)`` gives a count, ``discordant(...)`` the pair (b, c)."""

    rejected = rejected or (lambda **k: _default_rejected(**k))
    discordant = discordant or (lambda **k: (0, 0))
    blocks = {}
    for profile in csc.PROFILES:
        if profile in absent:
            blocks[profile] = None
            continue
        settings = full[profile]
        expected = csc.expected_counts(settings)
        rates, paired, cells = [], [], []
        for design, cell, fields in csc._cells(settings):
            n = expected[(design, cell)]
            meta = settings["design_index"][design]

            def applies(arm, design=design):
                base = csc._base_arm(arm)
                return csc.arm_applies(settings, design, base) and (design, base) not in skip

            for arm in csc.rate_arms(settings):
                for statistic in STATISTICS:
                    role = csc._role(fields, statistic)
                    k = rejected(design=design, cell=cell, kind=fields["kind"], arm=arm, statistic=statistic,
                                 role=role, n=n, delta=meta["delta"]) if applies(arm) else 0
                    rates.append({**{key: fields[key] for key in csc.CELL_FIELDS}, "arm": arm,
                                  "statistic": statistic, "role": role, "rule": settings["rejection"],
                                  "status": "complete" if applies(arm) else "not_run",
                                  "expected": n if applies(arm) else 0, "records": n if applies(arm) else 0,
                                  "rate": k / n if applies(arm) else None, "se": None,
                                  "count": n if applies(arm) else 0, "rejections": k, "failures": 0,
                                  "inconsistent": 0})
            for statistic in STATISTICS:
                role = csc._role(fields, statistic)
                for arm, reference, comparison in csc.pair_specs(settings):
                    both = applies(arm) and applies(reference)
                    b, c = discordant(design=design, cell=cell, kind=fields["kind"], arm=arm, statistic=statistic,
                                      role=role, n=n, delta=meta["delta"]) if both else (0, 0)
                    paired.append({**{key: fields[key] for key in csc.CELL_FIELDS}, "statistic": statistic,
                                   "role": role, "arm": arm, "reference_arm": reference, "comparison": comparison,
                                   "status": "complete" if both else "not_run", "n_paired": n if both else 0,
                                   "rate_reference": None, "rate_arm": None,
                                   "difference": (c - b) / n if both else None, "b": b, "c": c,
                                   **plsc.paired_p_values(b, c)})
            for arm in csc.ARMS:
                cells.append({**{key: fields[key] for key in csc.CELL_FIELDS}, "arm": arm,
                              "completed": n if applies(arm) else 0, "rank_median": 3, "auroc_median": 0.9,
                              "observed_delta_median": 1.0, "population_delta_median": 1.1,
                              "dispersion_stage_0_median": 1.0})
        block = {"settings": settings, "accounting": {"missing_records": 0}, "rates": rates, "paired": paired,
                 "cells": cells, "ranks": [{"design": "n80", "cell": "shape_e0.5", "rank": 2, "arm": arm,
                                            "statistic": "shape", "count": 10, "rejections": 4}
                                           for arm in ("pls_cv", "predictive_orthonormal")]}
        if profile in csc.REPRODUCED_PROFILES:
            block["reproduction"] = {"passed": passed, "rows_compared": 10, "differing": []}
        blocks[profile] = block
    return blocks


def _default_rejected(*, role, n, **_):
    if role == "control":
        return round(0.05 * n)
    if role == "target":
        return round(0.9 * n)
    return round(0.03 * n)


def _gains(**k):
    return (0, 20) if k["role"] == "target" else (0, 0)


@pytest.fixture(scope="module")
def gaining(full):
    return csc.decision(_table_blocks(full, discordant=_gains))


def test_gains_everywhere_qualify_and_propose_a_default_change(gaining) -> None:
    arms = gaining["arms"]
    for arm in ("predictive_orthonormal", "weighted"):
        assert {f: e["label"] for f, e in arms[arm]["families"].items()} == {f: "qualifies" for f in csc.FAMILIES}
        assert arms[arm]["dispersion"] == "holds"
        assert arms[arm]["verdict"] == "qualifies_everywhere"
    assert arms["sparse_block"]["verdict"] == "scoped_candidate"
    assert arms["sparse_block"]["families"]["core"]["partial"] is True
    assert arms["sparse_block"]["families"]["four_stages"]["label"] == "not_run"
    assert gaining["verdict"] == "propose_default_change"
    assert gaining["named_arm"] == "predictive_orthonormal"  # tied at zero readings above the bound: the order
    assert gaining["rule_version"] == 1 and gaining["gate"]["passed"]


def test_no_gain(full) -> None:
    result = csc.decision(_table_blocks(full))
    for arm in csc.CANDIDATE_ARMS:
        assert {e["label"] for e in result["arms"][arm]["families"].values()} <= {"no_gain", "not_run"}
        assert result["arms"][arm]["verdict"] == "keep_production"
    assert result["verdict"] == "keep_production"


def test_size_excess_in_one_family(full) -> None:
    def rejected(**k):
        if k["arm"] == "weighted" and k["role"] == "control" and k["design"] in csc.FAMILIES["core"]:
            return round(0.12 * k["n"])
        return _default_rejected(**k)

    result = csc.decision(_table_blocks(full, rejected=rejected, discordant=_gains))
    families = result["arms"]["weighted"]["families"]
    assert families["core"]["label"] == "size_excess"
    assert families["weak_stage"]["label"] == "qualifies"
    assert result["arms"]["weighted"]["verdict"] == "scoped_candidate"
    assert result["verdict"] == "propose_default_change" and result["named_arm"] == "predictive_orthonormal"


def test_power_loss(full) -> None:
    def discordant(**k):
        if k["arm"] == "predictive_orthonormal" and k["role"] == "target" and k["design"] == "seaad_table":
            return (30, 0)
        return _gains(**k)

    result = csc.decision(_table_blocks(full, discordant=discordant))
    entry = result["arms"]["predictive_orthonormal"]
    assert entry["families"]["unbalanced"]["label"] == "power_loss"
    assert entry["failing_criteria"]["unbalanced"][0] == "S2"
    assert entry["verdict"] == "scoped_candidate"


def test_specificity_at_the_core_effect_uses_the_bound(full) -> None:
    def rejected(**k):
        if k["arm"] == "weighted" and k["role"] == "off_target" and k["design"] in csc.FAMILIES["core"]:
            return round(0.2 * k["n"])
        return _default_rejected(**k)

    result = csc.decision(_table_blocks(full, rejected=rejected, discordant=_gains))
    core = result["arms"]["weighted"]["criteria"]["core"]["S3"]
    assert core["families"]["band"] is None and not core["families"]["strong"]["holds"]
    assert result["arms"]["weighted"]["families"]["core"]["label"] == "specificity_excess"


def test_specificity_below_the_core_effect_is_paired_excess(full) -> None:
    def rejected(**k):
        # High off-target rates for every arm below delta = 2.0: production's too, so no paired excess.
        if k["role"] == "off_target" and k["delta"] < 2.0:
            return round(0.3 * k["n"])
        return _default_rejected(**k)

    held = csc.decision(_table_blocks(full, rejected=rejected, discordant=_gains))
    weak = held["arms"]["weighted"]["criteria"]["weak_stage"]["S3"]
    assert weak["families"]["strong"] is None and weak["families"]["band"]["holds"]
    assert held["arms"]["weighted"]["families"]["weak_stage"]["label"] == "qualifies"
    below = held["descriptive"]["production_off_target_below_core"]
    assert below and all(row["rate"] == pytest.approx(0.3, abs=0.01) for row in below)

    def discordant(**k):
        if k["arm"] == "weighted" and k["role"] == "off_target" and k["delta"] < 2.0:
            return (0, 20)
        return _gains(**k)

    failed = csc.decision(_table_blocks(full, rejected=rejected, discordant=discordant))
    assert failed["arms"]["weighted"]["families"]["weak_stage"]["label"] == "specificity_excess"


def test_a_family_whose_designs_differ_in_delta(full) -> None:
    def discordant(**k):
        if k["arm"] == "predictive_orthonormal" and k["role"] == "off_target" and k["design"] == \
                "seaad_table_delta0.5":
            return (0, 25)
        return _gains(**k)

    result = csc.decision(_table_blocks(full, discordant=discordant))
    s3 = result["arms"]["predictive_orthonormal"]["criteria"]["unbalanced"]["S3"]
    assert s3["strong_designs"] == ["seaad_table", "n200_share0.8"]
    assert s3["families"]["strong"]["holds"] and not s3["families"]["band"]["holds"]
    assert s3["families"]["band"]["m"] == 18 and s3["families"]["strong"]["m"] == 36
    assert result["arms"]["predictive_orthonormal"]["families"]["unbalanced"]["label"] == "specificity_excess"


def test_a_partial_family_is_judged_on_its_designs(gaining) -> None:
    core = gaining["arms"]["sparse_block"]["families"]["core"]
    assert core["designs"] == {"n80": "not_run", "n200": "complete", "n500": "not_run"}
    assert core["label"] == "qualifies" and core["partial"] is True
    assert gaining["arms"]["sparse_block"]["criteria"]["core"]["S1"]["m"] == 6


def test_an_arm_that_qualifies_wherever_run_with_a_family_not_run_is_scoped(full) -> None:
    result = csc.decision(_table_blocks(full, discordant=_gains, absent=("adni_design",)))
    for arm in ("predictive_orthonormal", "weighted"):
        assert result["arms"][arm]["families"]["adni_design"]["label"] == "not_run"
        assert result["arms"][arm]["verdict"] == "scoped_candidate"
    assert result["verdict"] == "scoped_candidates"
    assert result["scoped_arms"] == list(csc.CANDIDATE_ARMS)
    assert result["gate"]["profiles"]["adni_design"] == {"status": "not_run"}


def test_the_tie_rule_names_the_arm_with_fewest_readings_above_the_bound(full) -> None:
    n = 200
    threshold = grade.exceedance_threshold(n)

    def rejected(**k):
        if (k["arm"] == "predictive_orthonormal" and k["design"] == "n200" and k["cell"] == "shape_e0.5"
                and k["statistic"] == "delta"):
            return threshold + 1
        return _default_rejected(**k)

    result = csc.decision(_table_blocks(full, rejected=rejected, discordant=_gains))
    assert result["arms"]["predictive_orthonormal"]["verdict"] == "qualifies_everywhere"
    assert result["arms"]["predictive_orthonormal"]["off_target_above_bound"] == 1
    assert result["named_arm"] == "weighted"


def test_a_failed_gate_still_records_the_criteria(full) -> None:
    result = csc.decision(_table_blocks(full, discordant=_gains, passed=False))
    assert result["verdict"] == "reproduction_failed"
    assert result["named_arm"] is None
    assert result["arms"]["predictive_orthonormal"]["verdict"] == "qualifies_everywhere"
    assert result["arms"]["weighted"]["criteria"]["core"]["S1"]["holds"]
    assert result["consequence"] == csc.CONSEQUENCES["reproduction_failed"]


DISPERSION_DESIGNS = [entry["design"] for entry in csc.dispersion_entries()]


def _dispersion_rejected(production_rate, candidate_rate):
    def rejected(**k):
        if k["design"] in DISPERSION_DESIGNS and k["role"] == "control":
            return round((production_rate if k["arm"] == "pls_cv" else candidate_rate) * k["n"])
        return _default_rejected(**k)

    return rejected


def test_dispersion_production_holds_and_a_candidate_fails(full) -> None:
    result = csc.decision(_table_blocks(full, rejected=_dispersion_rejected(0.05, 0.12), discordant=_gains))
    assert result["dispersion"]["production_holds"] is True
    assert result["dispersion"]["production_finding"] is None
    assert result["arms"]["weighted"]["dispersion"] == "fails"
    assert result["arms"]["weighted"]["verdict"] == "scoped_candidate"
    assert result["dispersion"]["arms"]["weighted"]["judged"]["m"] == 60


def test_dispersion_production_fails_and_a_candidate_not_in_excess_holds(full) -> None:
    result = csc.decision(_table_blocks(full, rejected=_dispersion_rejected(0.12, 0.12), discordant=_gains))
    assert result["dispersion"]["production_holds"] is False
    assert result["dispersion"]["production_finding"]["finding"] == "production_dispersion_excess"
    assert result["dispersion"]["arms"]["pls_cv"]["label"] == "production_dispersion_excess"
    assert result["arms"]["weighted"]["dispersion"] == "holds"
    assert result["dispersion"]["arms"]["weighted"]["rule"] == "count_and_ceiling_on_paired_excess"


def test_dispersion_production_fails_and_a_candidate_in_excess_fails(full) -> None:
    def discordant(**k):
        if k["arm"] == "weighted" and k["design"] in DISPERSION_DESIGNS:
            return (0, 25)
        return _gains(**k)

    result = csc.decision(_table_blocks(full, rejected=_dispersion_rejected(0.12, 0.12), discordant=discordant))
    assert result["dispersion"]["production_holds"] is False
    assert result["arms"]["weighted"]["dispersion"] == "fails"
    assert result["arms"]["predictive_orthonormal"]["dispersion"] == "holds"


def _incidental(b, c):
    def discordant(**k):
        if k["arm"] == "predictive_orthonormal" and k["role"] == "target" and k["kind"] in ("surgery",
                                                                                           "out_of_span"):
            return (b, c)
        return _gains(**k)

    return discordant


def test_incidental_labels(full) -> None:
    loses = csc.decision(_table_blocks(full, discordant=_incidental(80, 0)))["arms"]["predictive_orthonormal"]
    labels = loses["incidental_labels"]
    assert labels["m"] == 18 + 20
    assert labels["counts"]["candidate_loses"] == 38
    assert "loses at" in loses["consequence"] and "n80_three_blocks/" in loses["consequence"]
    gains = csc.decision(_table_blocks(full, discordant=_incidental(0, 80)))["arms"]["predictive_orthonormal"]
    assert gains["incidental_labels"]["counts"]["candidate_gains"] == 38
    # Significant at 0.025 / m but below the margin of 0.10: unresolved.
    below = csc.decision(_table_blocks(full, discordant=_incidental(30, 0)))["arms"]["predictive_orthonormal"]
    surgery = [e for e in below["incidental_labels"]["labels"] if e["kind"] == "surgery"]
    assert surgery and all(e["label"] == "unresolved" and e["p_loss"] < 1e-6 for e in surgery)


def test_incidental_rates_change_no_family_label_or_verdict(full, gaining) -> None:
    changed = csc.decision(_table_blocks(full, discordant=_incidental(80, 0)))
    assert changed["verdict"] == gaining["verdict"] and changed["named_arm"] == gaining["named_arm"]
    for arm in csc.CANDIDATE_ARMS:
        assert changed["arms"][arm]["families"] == gaining["arms"][arm]["families"]
        assert changed["arms"][arm]["verdict"] == gaining["arms"][arm]["verdict"]


def test_descriptive_tables_have_their_rows(full, gaining) -> None:
    tables = gaining["descriptive"]
    n_cells = sum(len(list(csc._cells(full[p]))) for p in csc.PROFILES)
    assert len(tables["fit"]) == len(tables["statistics"]) == n_cells * len(csc.ARMS)
    below_designs = {design for designs in csc.FAMILIES.values() for design in designs
                     if [full[p] for p in csc.PROFILES if design in full[p]["design_index"]][0]["design_index"][
                         design]["delta"] < 2.0}
    assert {row["design"] for row in tables["production_off_target_below_core"]} == below_designs
    assert len(tables["production_off_target_below_core"]) == 5 * 18 + 16
    assert len(tables["dispersion"]) == 10 * 2 * len(csc.ARMS)
    assert tables["rank_strata"] == []


def test_rank_strata_name_the_cells_an_arm_loses(full) -> None:
    def discordant(**k):
        if k["arm"] == "predictive_orthonormal" and k["role"] == "target" and k["design"] == "n80" \
                and k["cell"] == "shape_e0.5":
            return (12, 0)
        return _gains(**k)

    result = csc.decision(_table_blocks(full, discordant=discordant))
    strata = result["descriptive"]["rank_strata"]
    assert [(r["arm"], r["design"], r["cell"], r["statistic"]) for r in strata] == [
        ("predictive_orthonormal", "n80", "shape_e0.5", "shape")]
    assert strata[0]["strata"] == [{"rank": 2, "pls_cv": {"count": 10, "rejections": 4},
                                    "predictive_orthonormal": {"count": 10, "rejections": 4}}]


def test_decide_writes_the_decision_from_merged_profiles(tmp_path, monkeypatch) -> None:
    monkeypatch.setattr(csc, "evaluate_unit", _fake_records)
    csc.run_shard(tmp_path, _small(designs=["n80"]), 0, 1, arms=("pls_cv", "predictive_orthonormal", "weighted"))
    csc.merge(tmp_path, "committed")
    result = csc.decide(tmp_path)
    assert (tmp_path / csc.DECISION_NAME).exists()
    assert result["verdict"] == "reproduction_failed"  # two replicates do not reproduce the committed counts
    assert result["profiles"]["conditions"] == "not_run"
    assert result["arms"]["weighted"]["families"]["core"]["partial"] is True
