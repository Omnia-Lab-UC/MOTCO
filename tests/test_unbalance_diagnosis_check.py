"""``scripts/unbalance_diagnosis_check.py``: designs, seeds, the parent binding, the unit and its corrected arms,
merge, reproduction, the decision and the re-read of the group-dispersion families.

Units run at fast CV knobs and few permutations; the decision runs on synthetic tables over the real settings and
the committed group-dispersion decision. No real data are read.
"""

from __future__ import annotations

import json
import os
import shutil
import subprocess
import sys
from pathlib import Path
from types import SimpleNamespace
from unittest import mock

import numpy as np
import pandas as pd
import pytest


def _script(name: str):
    sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
    try:
        return __import__(name)
    finally:
        sys.path.pop(0)


ud = _script("unbalance_diagnosis_check")
ses = _script("stage_effect_scaling_check")
gd = ud.gd
ubc = ud.ubc
lsc = ud.lsc
ssc = ud.ssc

FAST_CV = {"cv1_splits": 3, "cv2_splits": 4, "n_repeats": 1, "max_components": 6, "random_state": 1203}
PERMS = 39
NEW = tuple(ud.NEW_DESIGNS)


@pytest.fixture(scope="module")
def settings():
    return ud.build_settings()


@pytest.fixture(scope="module")
def committed():
    _, stored = ud.parent_binding()
    return stored


@pytest.fixture(scope="module")
def parent_decision():
    return json.loads((ud.PARENT_DIR / gd.DECISION_NAME).read_text(encoding="utf-8"))


def _params(settings, design, cell="anchor", replicate=0):
    return ubc.stage_span_generator_params(lsc.Unit(design, cell, replicate), settings["unit_settings"])


def _cell_sizes(params) -> tuple[list[int], list[int]]:
    from motco.simulations import semisynthetic

    return semisynthetic._group_stage_sizes(semisynthetic._stage_sizes(params), params)


# ── 1.1 Settings, units and the parent ───────────────────────────────────────


def test_unit_counts(settings) -> None:
    units = ud.run_units(settings)
    assert len(units) == 7_500
    assert len({ubc.unit_key(unit) for unit in units}) == 7_500
    assert settings["replicates"] == {ud.RERUN_DESIGN: 500, **{design: 1000 for design in NEW}}
    assert sum(1 for unit in units if unit.design == ud.RERUN_DESIGN) == 500
    assert {unit.cell for unit in units} == {"anchor"}
    assert settings["unit_settings"]["permutations"] == 999
    assert settings["unit_settings"]["recorded_spaces"] == ["pca_matched"]
    assert settings["unit_settings"]["designs"] == [ud.RERUN_DESIGN, *NEW]
    assert settings["min_available"] == 900


def test_the_parent_files_are_bound(settings) -> None:
    expected = {f"results/group-dispersion-2026-11-24/{name}": ssc.sha256_file(ud.PARENT_DIR / name)
                for name in ud.PARENT_FILE_NAMES}
    assert settings["parent_files"] == expected
    assert len(expected) == 5
    assert settings["parent"]["dir"] == "results/group-dispersion-2026-11-24"
    stored = json.loads((ud.PARENT_DIR / "columns" / "settings.json").read_text())
    assert settings["parent"]["columns_settings_hash"] == stored["hash"]


def _fake_record(unit, units, *, rejections=None, status="completed", rank=2, within="consistent",
                 corrected="consistent", corrected_rejects=None):
    fields = ubc.unit_fields(unit, units)
    base = {**fields, "settings_hash": lsc.settings_hash(units), "status": status,
            "runtime_seconds": {"unit": 1.0, "fit": 0.5, "within_group": 0.2, "corrected": 0.3}}
    if status != "completed":
        return {**base, "error_type": "RuntimeError", "error_message": "patched",
                "within_group": {"status": "skipped", "reason": "the parent's unit failed"},
                "corrected": {"status": "skipped", "reason": "the parent's unit failed"}}
    rejections = rejections or {}
    arms = {arm: {"p_value": 0.01 if rejections.get(arm) else 0.5, "rejects": bool(rejections.get(arm, False))}
            for arm in ubc.ARMS}
    within_block = {"status": "completed", "strata": "group", "rank": rank, "same_statistics": within == "consistent",
                    "rejects": {s: False for s in gd.STATISTICS}}
    if corrected == "failed":
        corrected_block = {"status": "failed", "error_type": "RuntimeError", "error_message": "patched"}
    else:
        corrected_rejects = corrected_rejects or {}
        corrected_block = {"status": "completed", "rank": rank, "same_statistics": corrected != "inconsistent",
                           "pair_statistics": {"delta": 1.0, "angle": 2.0, "shape": 0.1}}
        for correction in ud.CORRECTIONS:
            if corrected == f"{correction}_failed":
                corrected_block[correction] = {"status": "failed", "error_type": "ValueError", "error_message": "x"}
            else:
                corrected_block[correction] = {
                    "status": "completed", "scale": {"min": 1.0, "max": 1.0}, "p_values": {},
                    "rejects": {s: bool(corrected_rejects.get(f"{correction}_{s}", False)) for s in gd.STATISTICS}}
    return {
        **base,
        "blocks": ["methylation", "expression", "proteomics"],
        "n_features": 658,
        "n_population_null": 100,
        "production": {"rank": rank, "cv_mean_auroc": 0.9,
                       "span_capture": {"group_difference": 0.5, "pooled_stage": 0.9, "random_share": 0.01}},
        "scan": {"hits": 0, "hits_population_null": 0},
        "fdr": {"discoveries": 0, "discoveries_population_null": 0},
        "spaces": {"pca_matched": {"dimension": rank}},
        "arms": arms,
        "within_group": within_block,
        "corrected": corrected_block,
    }


def _small(**overrides):
    overrides.setdefault("designs", [ud.RERUN_DESIGN, "n42_r0.5"])
    overrides.setdefault("replicates", 2)
    overrides.setdefault("permutations", PERMS)
    return ud.build_settings(**overrides)


@pytest.fixture()
def fake_units(monkeypatch):
    monkeypatch.setattr(ud, "evaluate_unit", lambda unit, units: _fake_record(unit, units))


def test_resume_skips_complete_units(tmp_path, monkeypatch) -> None:
    settings = _small()
    calls = []

    def fake(unit, units):
        calls.append(ubc.unit_key(unit))
        return _fake_record(unit, units)

    monkeypatch.setattr(ud, "evaluate_unit", fake)
    assert ud.run_shard(tmp_path, settings, 0, 1, limit=3) == 3
    assert ud.run_shard(tmp_path, settings, 0, 1) == 4 - 3
    assert len(calls) == 4 and len(set(calls)) == 4
    assert ud.run_shard(tmp_path, settings, 0, 1) == 0
    summary = ud.merge(tmp_path)
    assert summary["complete_units"] == 4 and summary["missing_units"] == 0


def test_an_interrupted_line_is_rerun_not_duplicated(tmp_path, fake_units) -> None:
    settings = _small()
    ud.run_shard(tmp_path, settings, 0, 1, limit=2)
    path = lsc.shard_path(tmp_path, 0, 1)
    with path.open("a", encoding="utf-8") as handle:
        handle.write('{"design": "n42_r0.5", "cell": "anc')
    assert ud.run_shard(tmp_path, settings, 0, 1) == 2
    assert len(path.read_text().splitlines()) == 4


def test_shards_are_round_robin(tmp_path, fake_units) -> None:
    settings = _small(replicates=3)
    assert [ud.run_shard(tmp_path, settings, index, 4) for index in range(4)] == [2, 2, 1, 1]
    assert ud.merge(tmp_path)["complete_units"] == 6


def test_a_resume_under_other_settings_is_refused_naming_the_keys(tmp_path, fake_units) -> None:
    ud.run_shard(tmp_path, _small(), 0, 1, limit=1)
    with pytest.raises(SystemExit, match=r"differing: unit_settings.permutations"):
        ud.run_shard(tmp_path, _small(permutations=19), 0, 1)
    with pytest.raises(SystemExit, match=r"replicates.n42_r0.5"):
        ud.run_shard(tmp_path, _small(replicates=3), 0, 1)


def test_merge_refuses_records_written_under_other_settings(tmp_path, monkeypatch) -> None:
    monkeypatch.setattr(ud, "evaluate_unit",
                        lambda unit, units: {**_fake_record(unit, units), "settings_hash": "other"})
    ud.run_shard(tmp_path, _small(), 0, 1, limit=2)
    with pytest.raises(SystemExit, match=r"written under other settings"):
        ud.merge(tmp_path)


def _parent_copy(tmp_path: Path) -> Path:
    target = tmp_path / "parent"
    (target / "columns").mkdir(parents=True)
    for name in ud.PARENT_FILE_NAMES:
        shutil.copy(ud.PARENT_DIR / name, target / name)
    return target


@pytest.mark.parametrize("name", ud.PARENT_FILE_NAMES)
def test_a_changed_parent_file_is_refused_by_name(tmp_path, fake_units, name) -> None:
    parent = _parent_copy(tmp_path)
    settings = _small(parent_dir=parent)
    assert str(parent / name) in settings["parent_files"]
    out = tmp_path / "run"
    ud.run_shard(out, settings, 0, 1)
    ud.merge(out)
    ud.decide(out)
    with (parent / name).open("a") as handle:
        handle.write("\n")
    for call in (lambda: ud.run_shard(out, settings, 0, 1), lambda: ud.merge(out), lambda: ud.decide(out)):
        with pytest.raises(SystemExit, match=rf"parent/{name} has changed since the settings were written"):
            call()


def test_a_parent_settings_file_failing_its_hash_is_refused(tmp_path) -> None:
    parent = _parent_copy(tmp_path)
    path = parent / "columns" / "settings.json"
    stored = json.loads(path.read_text())
    stored["hash"] = "0" * 64
    path.write_text(json.dumps(stored))
    with pytest.raises(SystemExit, match=r"settings.json fails its hash"):
        _small(parent_dir=parent)


def _edit_parent_settings(parent: Path, edit) -> None:
    path = parent / "columns" / "settings.json"
    stored = json.loads(path.read_text())
    edit(stored["settings"])
    stored["hash"] = lsc.settings_hash(stored["settings"])
    path.write_text(json.dumps(stored))


def test_parent_units_other_than_the_group_dispersion_drivers_are_refused(tmp_path) -> None:
    parent = _parent_copy(tmp_path)
    _edit_parent_settings(parent, lambda s: s["unit_settings"]["rules"].update({"delta": "p <= alpha"}))
    with pytest.raises(SystemExit, match=r"differs from the units group_dispersion_check writes in rules.delta"):
        _small(parent_dir=parent)


def test_a_parent_rerun_design_other_than_the_committed_one_is_refused(tmp_path) -> None:
    parent = _parent_copy(tmp_path)

    def edit(settings):
        settings["unit_settings"]["units"]["designs"]["n80_r0.8"]["generator"]["group_ratio"] = 0.75

    _edit_parent_settings(parent, edit)
    with pytest.raises(SystemExit, match=r"'n80_r0.8' is not the equal-dispersion three-block design"):
        _small(parent_dir=parent)


def test_unknown_designs_are_refused() -> None:
    with pytest.raises(SystemExit, match=r"Unknown design\(s\) \['n80_r0.2_s2'\]"):
        ud.build_settings(designs=["n80_r0.2_s2"])


# ── 1.2 The design builder ───────────────────────────────────────────────────


def test_the_rerun_design_is_the_parents(settings, committed) -> None:
    parent_units = committed["unit_settings"]
    units = settings["unit_settings"]
    assert units["units"]["designs"][ud.RERUN_DESIGN] == parent_units["units"]["designs"][ud.RERUN_DESIGN]
    assert units["units"]["cells"] == parent_units["units"]["cells"]
    for replicate in (0, 1, 137, 250, 499):
        unit = lsc.Unit(ud.RERUN_DESIGN, "anchor", replicate)
        assert ubc.unit_fields(unit, units) == ubc.unit_fields(unit, parent_units)
        assert ubc.stage_span_generator_params(unit, units) == ubc.stage_span_generator_params(unit, parent_units)
    assert settings["design_index"][ud.RERUN_DESIGN]["family"] is None
    assert settings["rerun_design"] == ud.RERUN_DESIGN


def test_the_unit_settings_differ_from_the_parents_only_in_designs_and_grade(settings, committed) -> None:
    mine = ubc._canonical(settings["unit_settings"])
    theirs = ubc._canonical(committed["unit_settings"])
    differing = lsc.settings_difference(theirs, mine)
    assert differing
    assert all(key.startswith(("units.designs.", "designs", "replicates.")) for key in differing), differing


def test_each_new_design_differs_from_the_rerun_only_in_size_and_balance(settings) -> None:
    designs = settings["unit_settings"]["units"]["designs"]
    reference = designs[ud.RERUN_DESIGN]["generator"]
    for design in NEW:
        generator = designs[design]["generator"]
        assert designs[design]["layers"] is None
        assert set(generator) == set(reference)
        assert {key for key in generator if generator[key] != reference[key]} <= {"n_samples", "group_ratio"}
        n_samples, share, *_ = ud.NEW_DESIGNS[design]
        assert (generator["n_samples"], generator["group_ratio"]) == (n_samples, share)
        assert "group_dispersion_scale" not in generator
    assert designs["n80_r0.8_b"]["generator"] == reference


def test_the_realized_cells_are_the_declared_ones(settings) -> None:
    expected = {"n80_r0.8": (22, 5), "n80_r0.8_b": (22, 5), "n200_r0.925": (62, 5), "n42_r0.5": (7, 7),
                "n80_r0.74": (20, 7), "n200_r0.895": (60, 7), "n80_r0.65": (18, 9), "n80_r0.5_b": (14, 13)}
    for design, (larger, smaller) in expected.items():
        meta = settings["design_index"][design]
        assert (meta["larger_per_stage"], meta["smaller_per_stage"]) == (larger, smaller)
        sizes_a, sizes_b = _cell_sizes(_params(settings, design))
        assert len(sizes_a) == len(sizes_b) == 3
        assert sizes_b == [smaller] * 3, design
        assert max(sizes_a) == larger and min(sizes_a) >= larger - 1, design
        assert sum(sizes_a) + sum(sizes_b) == meta["n_samples"]


def test_a_new_design_measures_the_three_blocks(settings) -> None:
    from motco.simulations.semisynthetic import generate_semisynthetic_trajectory

    dataset = generate_semisynthetic_trajectory(_params(settings, "n42_r0.5"))
    assert sum(getattr(dataset, layer).shape[1] for layer in ("methylation", "expression", "proteomics")) == 658
    assert len(dataset.metadata) == 42
    assert "group_dispersion" not in dataset.truth or dataset.truth["group_dispersion"] is None


def test_new_names_and_datasets_are_fresh(settings) -> None:
    group_dispersion = {entry["design"] for entry in (*gd.core_entries(), *gd.column_entries())}
    committed_names = (group_dispersion | set(ubc.STAGE_SPAN_DESIGNS) | set(ses.SATURATED_DESIGNS)
                       | set(ses.BAND_DESIGNS))
    assert ud.RERUN_DESIGN in group_dispersion
    assert not set(NEW) & committed_names
    new_seeds = {lsc.generator_seed(2043, design, replicate) for design in NEW for replicate in range(1000)}
    assert len(new_seeds) == 7_000
    # The fresh replicates share a committed design's generator, so they must share none of its seeds. (A 31-bit
    # seed may coincide across designs whose generators differ; their datasets are not shared.)
    designs = settings["unit_settings"]["units"]["designs"]
    pairs = {"n80_r0.8_b": (ud.RERUN_DESIGN, designs[ud.RERUN_DESIGN]["generator"]),
             "n80_r0.5_b": ("n80", ubc.stage_span_generator(*ubc.STAGE_SPAN_DESIGNS["n80"]))}
    for fresh, (committed, generator) in pairs.items():
        assert designs[fresh]["generator"] == generator
        theirs = {lsc.generator_seed(2043, committed, replicate) for replicate in range(500)}
        mine = {lsc.generator_seed(2043, fresh, replicate) for replicate in range(1000)}
        assert not theirs & mine, fresh


def test_families_take_every_new_design_once(settings) -> None:
    members = settings["families"]
    assert {family: len(designs) for family, designs in members.items()} == {
        "balanced_small": 1, "unbalanced_small": 4, "reference": 2}
    flat = [design for designs in members.values() for design in designs]
    assert sorted(flat) == sorted(NEW)
    assert ud.RERUN_DESIGN not in flat


# ── 1.3 and 1.4 The unit and the corrected arms (real units, fast knobs) ─────


def _fast(settings):
    settings["unit_settings"]["units"]["pls_cv_params"] = dict(FAST_CV)
    return settings


UNIT = lsc.Unit(ud.RERUN_DESIGN, "anchor", 1)


@pytest.fixture(scope="module")
def unit_records():
    mine = _fast(_small())
    parent = _fast(gd.profile_settings("columns", designs=[ud.RERUN_DESIGN], replicates=2, permutations=PERMS))
    return ud.evaluate_unit(UNIT, mine["unit_settings"]), gd.evaluate_unit(UNIT, parent["unit_settings"]), mine


def test_a_unit_record_equals_the_parent_drivers(unit_records) -> None:
    record, parent, mine = unit_records
    assert record["status"] == "completed", record.get("error_message")
    assert record["settings_hash"] == lsc.settings_hash(mine["unit_settings"])

    def strip(r):
        out = {key: value for key, value in r.items() if key not in ("settings_hash", "runtime_seconds", "corrected")}
        return out

    assert strip(record) == strip(parent)
    assert record["within_group"]["status"] == "completed"
    assert record["generator_seed"] == lsc.generator_seed(2043, UNIT.design, UNIT.replicate)
    assert {"group_dispersion_unit", "parent_unit", "within_group", "corrected", "unit"} <= set(
        record["runtime_seconds"])


def test_the_refit_reproduces_the_production_statistics(unit_records) -> None:
    record, _, _ = unit_records
    corrected = record["corrected"]
    assert corrected["status"] == "completed", corrected.get("error_message")
    assert corrected["same_statistics"] is True
    assert corrected["pair_statistics"] == record["production"]["pair_statistics"]
    assert corrected["rank"] == record["production"]["rank"]
    for correction in ud.CORRECTIONS:
        block = corrected[correction]
        assert block["status"] == "completed", block.get("error_message")
        assert set(block["p_values"]) == set(ud.STATISTICS)
        assert set(block["rejects"]) == set(ud.STATISTICS)
    # Without covariates the reduced model is intercept, group and two stage columns: rank 4.
    assert corrected["df"]["scale"]["min"] == pytest.approx((80 / 76) ** 0.5)
    assert corrected["df"]["scale"]["min"] == corrected["df"]["scale"]["max"]
    assert corrected["leverage"]["scale"]["min"] >= 1.0


def _forced_inputs(unit, units, rank):
    from motco.simulations.evaluation import integrate_semisynthetic_dataset

    fields = ubc.unit_fields(unit, units)
    dataset, params = ubc.unit_inputs(unit, units, fields)
    return dataset, params, integrate_semisynthetic_dataset(dataset, gd.forced_rank_params(params, rank))


def test_no_correction_reproduces_the_production_p_values(unit_records) -> None:
    record, _, mine = unit_records
    dataset, params, latent = _forced_inputs(UNIT, mine["unit_settings"], record["production"]["rank"])
    plain = ud.corrected_reading(dataset, latent.matrix, params, residual_correction=None)
    assert plain["pair_statistics"] == record["production"]["pair_statistics"]
    assert plain["p_values"] == record["production"]["p_values"]
    assert "scale" not in plain
    for correction in ud.CORRECTIONS:
        reading = ud.corrected_reading(dataset, latent.matrix, params, residual_correction=correction)
        assert reading["p_values"] == record["corrected"][correction]["p_values"]


def test_a_unit_scale_of_one_reproduces_the_uncorrected_nulls(unit_records, monkeypatch) -> None:
    from motco.stats import permutation

    record, _, mine = unit_records
    dataset, params, latent = _forced_inputs(UNIT, mine["unit_settings"], record["production"]["rank"])
    plain = ud.corrected_reading(dataset, latent.matrix, params, residual_correction=None, return_nulls=True)
    monkeypatch.setattr(permutation, "residual_correction_scale", lambda model, correction: np.ones(len(model)))
    for correction in ud.CORRECTIONS:
        unit_scale = ud.corrected_reading(dataset, latent.matrix, params, residual_correction=correction,
                                          return_nulls=True)
        assert unit_scale["nulls"] == plain["nulls"]
        assert unit_scale["p_values"] == plain["p_values"]


def test_each_correction_uses_the_uncorrected_permutation_indices(monkeypatch) -> None:
    from motco.simulations.evaluation import SimulationEvaluationParams, build_simulation_trajectory_design
    from motco.stats import permutation
    from tests.test_residual_correction import _reduced_fit

    sizes = ((6, 5, 9), (2, 3, 2))
    groups = np.concatenate([[g] * sizes[g][s] for g in range(2) for s in range(3)])
    stages = np.concatenate([[s] * sizes[g][s] for g in range(2) for s in range(3)])
    metadata = pd.DataFrame({"group": np.where(groups == 0, "A", "B"), "stage": stages.astype(str)})
    rng = np.random.default_rng(4)
    Y = 2.0 * rng.normal(size=(3, 3))[stages] + rng.normal(size=(groups.size, 3))
    params = SimulationEvaluationParams(permutations=PERMS, seed=5, orientation_estimator="net_displacement")
    dataset = SimpleNamespace(metadata=metadata)
    design = build_simulation_trajectory_design(metadata, group_col="group", stage_col="stage")
    reduced = np.asarray(design.model_reduced, dtype=float)
    y_hat, y_res = _reduced_fit(Y, reduced)

    def draws_under(correction):
        captured = []
        real = permutation.estimate_difference

        def capture(y_random, *args, **kwargs):
            captured.append(np.asarray(y_random, dtype=float).copy())
            return real(y_random, *args, **kwargs)

        with monkeypatch.context() as patch:
            patch.setattr(permutation, "estimate_difference", capture)
            ud.corrected_reading(dataset, Y, params, residual_correction=correction)
        return captured

    plain = draws_under(None)
    assert len(plain) == PERMS
    indices = []
    for y_random in plain:
        distances = (((y_random - y_hat)[:, None, :] - y_res[None, :, :]) ** 2).sum(axis=2)
        index = distances.argmin(axis=1)
        assert sorted(index.tolist()) == list(range(groups.size))
        np.testing.assert_allclose(y_random, y_hat + y_res[index], atol=1e-12)
        indices.append(index)
    assert any(np.any(groups[index] != groups) for index in indices)  # unrestricted draws
    for correction in ud.CORRECTIONS:
        scale = permutation.residual_correction_scale(reduced, correction)
        assert not np.allclose(scale, 1.0)
        corrected = draws_under(correction)
        assert len(corrected) == PERMS
        for y_random, index in zip(corrected, indices, strict=True):
            np.testing.assert_allclose(y_random, y_hat + (y_res * scale[:, None])[index], atol=1e-12)


def test_an_inconsistent_refit_is_marked_for_these_arms_only(unit_records, monkeypatch) -> None:
    record, _, mine = unit_records
    real = ud.corrected_reading

    def shifted(*args, **kwargs):
        reading = real(*args, **kwargs)
        statistics = dict(reading["pair_statistics"])
        statistics["delta"] = statistics["delta"] + 1e-9
        return {**reading, "pair_statistics": statistics}

    monkeypatch.setattr(ud, "corrected_reading", shifted)
    corrected = ud.corrected_arms(UNIT, mine["unit_settings"], record)
    assert corrected["status"] == "completed" and corrected["same_statistics"] is False
    marked = {**record, "corrected": corrected}
    rates, table = ud.corrected_rows([marked], mine["unit_settings"])
    row = next(r for r in table if r["design"] == UNIT.design)
    assert (row["refit_completed"], row["inconsistent"]) == (1, 1)
    assert all(r["count"] == 0 and r["failures"] == 1 for r in rates if r["design"] == UNIT.design)
    index = ubc._index(ud.rate_rows([marked], mine["unit_settings"]), "design", "cell", "arm")
    assert index[(UNIT.design, "anchor", "delta")]["count"] == 1
    assert index[(UNIT.design, "anchor", "within_group_delta")]["count"] == 1
    assert index[(UNIT.design, "anchor", "df_delta")]["count"] == 0


def test_a_failed_unit_is_recorded_and_the_shard_continues(tmp_path) -> None:
    settings = _fast(_small())

    def broken(*args, **kwargs):
        raise RuntimeError("patched scan failure")

    with mock.patch("motco.stats.interaction.interaction_feature_scan", broken):
        assert ud.run_shard(tmp_path, settings, 0, 4, limit=2) == 1
        assert ud.run_shard(tmp_path, settings, 1, 4, limit=2) == 1
    records = [json.loads(line) for path in sorted((tmp_path / "records").glob("*.jsonl")) for line in path.open()]
    assert len(records) == 2
    for record in records:
        assert record["status"] == "failed"
        assert (record["error_type"], record["error_message"]) == ("RuntimeError", "patched scan failure")
        assert record["within_group"]["status"] == "skipped"
        assert record["corrected"]["status"] == "skipped"


def test_a_failed_correction_is_recorded_for_that_correction_only(unit_records, monkeypatch) -> None:
    record, _, mine = unit_records
    real = ud.corrected_reading

    def leverage_fails(*args, residual_correction, **kwargs):
        if residual_correction == "leverage":
            raise ValueError("patched leverage failure")
        return real(*args, residual_correction=residual_correction, **kwargs)

    monkeypatch.setattr(ud, "corrected_reading", leverage_fails)
    corrected = ud.corrected_arms(UNIT, mine["unit_settings"], record)
    assert corrected["same_statistics"] is True
    assert corrected["df"]["status"] == "completed"
    assert corrected["leverage"] == {"status": "failed", "error_type": "ValueError",
                                     "error_message": "patched leverage failure"}


# ── 1.5 The sbatch array ─────────────────────────────────────────────────────


SBATCH = Path(__file__).resolve().parents[1] / "scripts" / "unbalance_diagnosis_array.sbatch"


def test_the_sbatch_parses_and_prints_its_command(tmp_path) -> None:
    subprocess.run(["bash", "-n", str(SBATCH)], check=True)
    env = {**os.environ, "SLURM_ARRAY_TASK_ID": "3", "N_SHARDS": "150", "UD_OUT": str(tmp_path / "run"),
           "UD_DRY_RUN": "1", "UD_ARGS": "--limit 1", "SLURM_SUBMIT_DIR": str(SBATCH.parents[1])}
    result = subprocess.run(["bash", str(SBATCH)], env=env, check=True, capture_output=True, text=True)
    line = result.stdout.strip()
    assert line.startswith(f"PYTHONPATH={SBATCH.parents[1] / 'src'} python ")
    assert line.endswith(f"run --out-dir {tmp_path / 'run'} --shard 3/150 --limit 1")
    assert not (tmp_path / "run").exists()
    relative = subprocess.run(["bash", str(SBATCH)], env={**env, "UD_OUT": "relative"}, capture_output=True,
                              text=True)
    assert relative.returncode == 2 and "must be an absolute path" in relative.stderr


# ── 2.1 Merge ────────────────────────────────────────────────────────────────


def _write_shard(out_dir: Path, settings, records, name="shard_000_of_001.jsonl") -> None:
    lsc.ensure_settings(Path(out_dir), settings)
    (Path(out_dir) / "records").mkdir(parents=True, exist_ok=True)
    with (Path(out_dir) / "records" / name).open("a", encoding="utf-8") as handle:
        for record in records:
            handle.write(lsc._line(record))


def test_merge_writes_one_row_per_design_and_arm_and_counts_failures(tmp_path) -> None:
    settings = _small(replicates=4)
    units = settings["unit_settings"]
    records = []
    special = {
        0: dict(status="failed"),
        1: dict(rejections={"delta": True}, corrected="inconsistent"),
        2: dict(corrected="failed"),
        3: dict(corrected="leverage_failed", corrected_rejects={"df_angle": True}),
    }
    for unit in ud.run_units(settings):
        if unit.design == "n42_r0.5":
            records.append(_fake_record(unit, units, **special[unit.replicate]))
        else:
            records.append(_fake_record(unit, units, corrected_rejects={"leverage_shape": True}))
    _write_shard(tmp_path, settings, records)
    # A duplicate written to a second shard file is counted once.
    _write_shard(tmp_path, settings, records[:1], name="shard_001_of_002.jsonl")
    summary = ud.merge(tmp_path)
    assert summary["complete_units"] == 8 and summary["failed_units"] == 1
    assert summary["corrected_inconsistent"] == 1
    rates = ubc._read_csv(tmp_path / "rate_table.csv")
    assert len(rates) == 2 * len(ud.RATE_ARMS)
    keys = [(r["design"], r["cell"], r["arm"]) for r in rates]
    assert len(set(keys)) == len(keys)
    assert [r["arm"] for r in rates[:len(ud.RATE_ARMS)]] == list(ud.RATE_ARMS)
    index = ubc._index(rates, "design", "cell", "arm")
    delta = index[("n42_r0.5", "anchor", "delta")]
    assert (int(delta["count"]), int(delta["rejections"]), int(delta["failures"])) == (3, 1, 1)
    df_angle = index[("n42_r0.5", "anchor", "df_angle")]
    assert (int(df_angle["count"]), int(df_angle["rejections"]), int(df_angle["failures"])) == (1, 1, 3)
    leverage = index[("n42_r0.5", "anchor", "leverage_angle")]
    assert (int(leverage["count"]), int(leverage["failures"])) == (0, 4)
    assert index[("n42_r0.5", "anchor", "df_angle")]["rule"] == "p < alpha"
    rerun = index[(ud.RERUN_DESIGN, "anchor", "leverage_shape")]
    assert (int(rerun["count"]), int(rerun["rejections"]), int(rerun["failures"])) == (4, 4, 0)
    table = ubc._index(ubc._read_csv(tmp_path / "corrected_table.csv"), "design", "cell")
    row = table[("n42_r0.5", "anchor")]
    assert {key: int(row[key]) for key in ud.CORRECTED_COLUMNS[len(ubc.CELL_FIELDS):]} == {
        "units": 4, "parent_completed": 3, "refit_completed": 2, "refit_failed": 1, "inconsistent": 1,
        "df_failed": 0, "leverage_failed": 1}
    accounting = json.loads((tmp_path / "accounting.json").read_text())
    assert accounting["duplicate_records"] == 1
    assert accounting["failures_by_cell"] == {"n42_r0.5/anchor": 1}
    assert accounting["corrected"]["inconsistent_by_design"] == {"n42_r0.5/anchor": 1}
    assert accounting["corrected"]["refit_failed_by_design"] == {"n42_r0.5/anchor": 1}
    assert accounting["corrected"]["leverage_failed_by_design"] == {"n42_r0.5/anchor": 1}
    assert accounting["expected_by_design"] == {ud.RERUN_DESIGN: 4, "n42_r0.5": 4}
    for name in ("cell_table.csv", "within_group_table.csv", "reproduction.json"):
        assert (tmp_path / name).exists()


# ── 2.2 Reproduction ─────────────────────────────────────────────────────────


@pytest.fixture(scope="module")
def committed_tables():
    return (lsc._read_csv(ud.PARENT_DIR / "columns" / "rate_table.csv"),
            lsc._read_csv(ud.PARENT_DIR / "columns" / "within_group_table.csv"))


def _as_mine(committed_tables):
    rates, within = committed_tables
    mine_rates = [dict(row) for row in rates if row["design"] == ud.RERUN_DESIGN]
    mine_rates += [{**mine_rates[0], "arm": arm, "count": "7", "rejections": "1"} for arm in ud.CORRECTED_ARMS]
    return mine_rates, [dict(row) for row in within if row["design"] == ud.RERUN_DESIGN]


def test_the_committed_tables_resolve_every_compared_row(settings, committed_tables) -> None:
    rates, within = _as_mine(committed_tables)
    result = ud.rerun_reproduction(rates, within, *committed_tables, settings)
    assert result["passed"] is True
    assert result["rows_compared"] == len(gd.RATE_ARMS) + 1 == 22
    assert result["unresolved_rows"] == [] and result["differing"] == []


def test_the_reproduction_fails_and_lists_every_differing_row(settings, committed_tables) -> None:
    rates, within = _as_mine(committed_tables)
    for row in rates:
        if row["arm"] in ("delta", "within_group_shape"):
            row["rejections"] = str(int(row["rejections"]) + 1)
    within[0]["inconsistent"] = "2"
    result = ud.rerun_reproduction(rates, within, *committed_tables, settings)
    assert result["passed"] is False
    assert result["rows_compared"] == 22
    assert {(d["table"], d["arm"]) for d in result["differing"]} == {
        ("rate_table", "delta"), ("rate_table", "within_group_shape"), ("within_group_table", "within_group")}
    for entry in result["differing"]:
        if entry["table"] == "rate_table":
            assert entry["reproduced"]["rejections"] == entry["committed"]["rejections"] + 1
        else:
            assert (entry["committed"]["inconsistent"], entry["reproduced"]["inconsistent"]) == (0, 2)
    missing = ud.rerun_reproduction([row for row in rates if row["arm"] != "fdr"], within, *committed_tables,
                                    settings)
    assert missing["passed"] is False
    assert missing["unresolved_rows"] == [
        {"table": "rate_table", "cell": "anchor", "arm": "fdr", "missing_from": "this run"}]


def test_a_run_without_the_rerun_design_does_not_pass(committed_tables) -> None:
    settings = ud.build_settings(designs=["n42_r0.5"], replicates=2)
    result = ud.rerun_reproduction([], [], *committed_tables, settings)
    assert result["passed"] is False and result["rows_compared"] == 0


# ── 2.3 to 2.5 Decision and re-read (synthetic tables over the real settings) ─


def _rates(settings, rejections=50):
    units = settings["unit_settings"]
    return [{"design": design, "cell": "anchor", "arm": arm, "count": settings["replicates"][design],
             "rejections": rejections * settings["replicates"][design] // 1000, "rate": rejections / 1000}
            for design in units["designs"] for arm in ud.RATE_ARMS]


def _cells(settings):
    return [{"design": design, "cell": "anchor", "completed": settings["replicates"][design], "rank_q10": 2,
             "rank_median": 3, "rank_q90": 4, "auroc_median": 0.99} for design in settings["unit_settings"]["designs"]]


@pytest.fixture()
def block(settings):
    return {"settings": settings,
            "accounting": {"missing_units": 0, "failed_units": 0, "within_group": {}, "corrected": {}},
            "rates": _rates(settings), "cells": _cells(settings),
            "reproduction": {"passed": True, "rows_compared": 22}}


def _set(block, design, arm, rejections, count=None):
    for row in block["rates"]:
        if (row["design"], row["arm"]) == (design, arm):
            row["rejections"] = rejections
            if count is not None:
                row["count"] = count
            return
    raise KeyError((design, arm))


VERDICT_CASES = {
    "not_replicated": [],
    "small_cells": [("n42_r0.5", "delta", 70), ("n80_r0.65", "shape", 80)],
    "small_unbalanced_cells": [("n200_r0.925", "angle", 80)],
    "moderate_excess": [("n80_r0.5_b", "shape", 80)],
}


def test_families_have_their_sizes_and_the_rerun_enters_none(block, parent_decision) -> None:
    result = ud.decision(block, parent_decision)
    assert {family: entry["m"] for family, entry in result["families"].items()} == ud.FAMILY_M_EXPECTED
    assert {family: (entry["allowed_exceedances"], entry["ceiling_at_median"])
            for family, entry in result["families"].items()} == {
        "balanced_small": (1, 0.065), "unbalanced_small": (1, 0.069), "reference": (1, 0.067)}
    observations = [row for family in ud.FAMILIES for row in result["criteria"]["size"][family]["observations"]]
    assert not [row for row in observations if row["design"] == ud.RERUN_DESIGN]
    keys = [(row["design"], row["arm"]) for row in observations]
    assert len(keys) == len(set(keys)) == 21
    assert {row["design"] for row in observations} == set(NEW)


@pytest.mark.parametrize("verdict", sorted(VERDICT_CASES))
def test_each_verdict_and_every_family_recorded(block, parent_decision, verdict) -> None:
    for design, arm, rejections in VERDICT_CASES[verdict]:
        _set(block, design, arm, rejections)
    result = ud.decision(block, parent_decision)
    assert result["verdict"] == verdict
    assert result["consequence"] == ud.CONSEQUENCES[verdict]
    assert set(result["families"]) == set(ud.FAMILIES)
    assert result["rule_version"] == 1
    failing = {family for family in ud.FAMILIES if not result["families"][family]["holds"]}
    assert failing == set(result["failing"])
    for design, arm, _ in VERDICT_CASES[verdict]:
        family = ud.NEW_DESIGNS[design][4]
        assert (design, arm) in {(row["design"], row["arm"]) for row in result["failing"][family]}


def test_the_gate_comes_first(block, parent_decision) -> None:
    block["reproduction"] = {"passed": False, "rows_compared": 22,
                             "differing": [{"table": "rate_table", "cell": "anchor", "arm": "delta"}]}
    _set(block, "n42_r0.5", "delta", 70)
    result = ud.decision(block, parent_decision)
    assert result["verdict"] == "reproduction_failed"
    assert result["families"]["balanced_small"]["holds"] is False
    assert set(result["families"]) == set(ud.FAMILIES)
    assert result["failing"]["reproduction"]["differing"][0]["arm"] == "delta"


def test_an_exceedance_within_the_allowed_count_holds(block, parent_decision) -> None:
    _set(block, "n80_r0.74", "delta", 64)  # above the bound, within the ceiling
    result = ud.decision(block, parent_decision)
    assert result["families"]["unbalanced_small"]["observed_exceedances"] == 1
    assert result["verdict"] == "not_replicated"
    _set(block, "n200_r0.895", "angle", 64)
    assert ud.decision(block, parent_decision)["verdict"] == "small_unbalanced_cells"


@pytest.mark.parametrize("count", [899, 0])
def test_a_design_short_of_units_fails_its_family(block, parent_decision, count) -> None:
    _set(block, "n80_r0.65", "angle", 40, count=count)
    result = ud.decision(block, parent_decision)
    assert result["families"]["reference"]["holds"] is False
    assert result["families"]["reference"]["n_unavailable"] == 1
    assert result["verdict"] == "moderate_excess"
    [row] = result["failing"]["reference"]
    assert (row["design"], row["arm"], row["status"]) == ("n80_r0.65", "angle", "unavailable")
    _set(block, "n80_r0.65", "angle", 40, count=900)
    assert ud.decision(block, parent_decision)["verdict"] == "not_replicated"


def test_the_pooled_delta_count(block, parent_decision) -> None:
    _set(block, ud.RERUN_DESIGN, "delta", 44)
    _set(block, "n80_r0.8_b", "delta", 61)
    pooled = ud.decision(block, parent_decision)["pooled_delta"]
    assert (pooled["n_rejected"], pooled["n_available"]) == (105, 1500)
    assert pooled["rate"] == pytest.approx(0.07)
    assert pooled["by_design"] == {ud.RERUN_DESIGN: {"n_available": 500, "n_rejected": 44},
                                   "n80_r0.8_b": {"n_available": 1000, "n_rejected": 61}}


def test_the_committed_re_read_with_nothing_changed_reproduces_the_committed_decision(parent_decision) -> None:
    judged = ud.reread_families(parent_decision)
    assert judged["verdict"] == parent_decision["verdict"] == "unbalance_excess"
    keys = ("holds", "m", "median_n_available", "allowed_exceedances", "observed_exceedances", "ceiling_at_median",
            "n_above_ceiling", "n_unavailable")
    for family in gd.FAMILIES:
        assert {key: judged["size"][family][key] for key in keys} == {
            key: parent_decision["families"][family][key] for key in keys}
        assert judged["size"][family]["observations"] == parent_decision["criteria"]["size"][family]["observations"]


@pytest.mark.parametrize("verdict", ["small_cells", "small_unbalanced_cells"])
def test_small_cell_verdicts_set_the_five_designs_aside(block, parent_decision, verdict) -> None:
    for design, arm, rejections in VERDICT_CASES[verdict]:
        _set(block, design, arm, rejections)
    result = ud.decision(block, parent_decision)
    assert result["verdict"] == verdict
    reread = result["reread"]
    assert reread["performed"] is True and reread["branch"] == "set_aside"
    assert reread["set_aside"] == list(ud.SET_ASIDE_DESIGNS) and reread["replaced"] is None
    assert {family: entry["m"] for family, entry in reread["families"].items()} == {
        "unbalanced_equal": 21, "balanced_unequal": 42, "wider_smaller": 54, "wider_larger": 54}
    assert all(entry["m"] == entry["m_expected"] for entry in reread["families"].values())
    designs = {row["design"] for family in gd.FAMILIES for row in reread["criteria"]["size"][family]["observations"]}
    assert not designs & set(ud.SET_ASIDE_DESIGNS)
    assert reread["families"]["unbalanced_equal"]["holds"] is True
    assert reread["verdict"] == "dispersion_excess"
    assert reread["consequence"] == gd.CONSEQUENCES["dispersion_excess"]
    assert reread["scope"] == ud.REREAD_SCOPES["set_aside"]
    assert reread["committed"]["verdict"] == "unbalance_excess"


def test_not_replicated_replaces_the_committed_reading(block, parent_decision) -> None:
    _set(block, "n80_r0.8_b", "delta", 52)
    result = ud.decision(block, parent_decision)
    assert result["verdict"] == "not_replicated"
    reread = result["reread"]
    assert reread["performed"] is True and reread["branch"] == "replaced"
    assert reread["set_aside"] == []
    assert reread["replaced"]["design"] == ud.RERUN_DESIGN and reread["replaced"]["by"] == "n80_r0.8_b"
    assert {family: entry["m"] for family, entry in reread["families"].items()} == gd.FAMILY_M_EXPECTED
    rows = reread["criteria"]["size"]["unbalanced_equal"]["observations"]
    assert ud.RERUN_DESIGN not in {row["design"] for row in rows}
    fresh = [row for row in rows if row["design"] == "n80_r0.8_b"]
    assert [(row["arm"], row["n_rejected"], row["n_available"]) for row in fresh] == [
        ("delta", 52, 1000), ("angle", 50, 1000), ("shape", 50, 1000)]
    committed_rows = parent_decision["criteria"]["size"]["unbalanced_equal"]["observations"]
    position = next(i for i, row in enumerate(committed_rows) if row["design"] == ud.RERUN_DESIGN)
    assert [row["design"] for row in rows[position:position + 3]] == ["n80_r0.8_b"] * 3
    assert reread["families"]["unbalanced_equal"]["holds"] is True
    assert reread["verdict"] == "dispersion_excess"
    for family in ("balanced_unequal", "wider_smaller", "wider_larger"):
        assert reread["families"][family]["observed_exceedances"] == parent_decision["families"][family][
            "observed_exceedances"]
    # A replacement above its ceiling, or one short of units, keeps unbalanced_equal failing.
    for rejected, available in ((90, 1000), (40, 899)):
        failing = ud.size_observation({"count": available, "rejections": rejected}, profile=ud.PROFILE,
                                      design="n80_r0.8_b", cell="anchor", arm="delta")
        judged = ud.reread_families(parent_decision, replacement=[failing, dict(fresh[1]), dict(fresh[2])])
        assert judged["holds"]["unbalanced_equal"] is False, (rejected, available)
        assert judged["verdict"] == "unbalance_excess"


@pytest.mark.parametrize("verdict", ["moderate_excess", "reproduction_failed"])
def test_no_re_read_under_the_other_verdicts(block, parent_decision, verdict) -> None:
    if verdict == "reproduction_failed":
        block["reproduction"] = {"passed": False, "rows_compared": 22, "differing": [{"arm": "delta"}]}
    else:
        for design, arm, rejections in VERDICT_CASES[verdict]:
            _set(block, design, arm, rejections)
    result = ud.decision(block, parent_decision)
    assert result["verdict"] == verdict
    reread = result["reread"]
    assert reread["performed"] is False and reread["branch"] is None
    assert reread["verdict"] == "unbalance_excess"
    assert reread["set_aside"] == [] and reread["replaced"] is None


def test_recorded_arms_are_labelled_and_leave_the_verdict_and_re_read(block, parent_decision) -> None:
    _set(block, "n42_r0.5", "delta", 70)
    before = ud.decision(block, parent_decision)
    for row in block["rates"]:
        if row["arm"] in ud.RECORDED_ARMS:
            row["rejections"] = row["count"] // 2
    after = ud.decision(block, parent_decision)
    assert after["verdict"] == before["verdict"] == "small_cells"
    assert after["families"] == before["families"]
    assert after["reread"] == before["reread"]
    assert set(after["recorded_arms"]) == set(ud.RECORDED_ARMS)
    assert {"df_delta", "leverage_shape", "within_group_angle", "production", "pca_matched"} <= set(
        after["recorded_arms"])
    for arm in ud.RECORDED_ARMS:
        assert set(after["recorded_arms"][arm]["families"]) == set(ud.FAMILIES)
        assert before["recorded_arms"][arm]["families"]["reference"]["label"] == "holds"
        assert after["recorded_arms"][arm]["families"]["reference"]["label"] == "exceeds"
    assert after["recorded_arms"]["df_angle"]["rule"] == "p < alpha"
    assert after["recorded_arms"]["production"]["rule"] == "min p <= alpha / 3"
    assert after["recorded_arms"]["production"]["families"]["unbalanced_small"]["m"] == 4


def test_counts_fit_and_rates_are_reported_per_design(block, parent_decision) -> None:
    _set(block, "n42_r0.5", "shape", 61)
    result = ud.decision(block, parent_decision)
    assert [entry["design"] for entry in result["counts"]] == [ud.RERUN_DESIGN, *NEW]
    entry = next(e for e in result["counts"] if e["design"] == "n42_r0.5")
    assert entry["counts"]["shape"]["n_rejected"] == 61 and entry["counts"]["delta"]["n_available"] == 1000
    assert (entry["larger_per_stage"], entry["smaller_per_stage"], entry["family"]) == (7, 7, "balanced_small")
    fit = next(e for e in result["fit"] if e["design"] == "n200_r0.925")
    assert (fit["rank_median"], fit["auroc_median"], fit["share"]) == (3.0, 0.99, 0.925)
    rates = next(e for e in result["rates"] if e["design"] == ud.RERUN_DESIGN)
    assert set(rates["rates"]) == set(ud.RATE_ARMS)


# ── End to end on fake units: the parent files are only read ────────────────


def test_decide_runs_from_the_aggregates_alone_and_leaves_the_parent_unchanged(tmp_path, fake_units) -> None:
    before = {name: ssc.sha256_file(ud.PARENT_DIR / name) for name in ud.PARENT_FILE_NAMES}
    settings = _small()
    ud.run_shard(tmp_path, settings, 0, 1)
    summary = ud.merge(tmp_path)
    assert summary["reproduction_passed"] is False
    first = ud.decide(tmp_path)
    assert first["verdict"] == "reproduction_failed"
    assert first["reread"]["performed"] is False
    written = (tmp_path / ud.DECISION_NAME).read_bytes()
    for path in (tmp_path / "merged.jsonl", *(tmp_path / "records").glob("*.jsonl")):
        path.unlink()
    ud.decide(tmp_path)
    assert (tmp_path / ud.DECISION_NAME).read_bytes() == written
    assert {name: ssc.sha256_file(ud.PARENT_DIR / name) for name in ud.PARENT_FILE_NAMES} == before


def test_decide_requires_the_merge(tmp_path, fake_units) -> None:
    ud.run_shard(tmp_path, _small(), 0, 1, limit=1)
    with pytest.raises(SystemExit, match=r"run `merge` first"):
        ud.decide(tmp_path)
