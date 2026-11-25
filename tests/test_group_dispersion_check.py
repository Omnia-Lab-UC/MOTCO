"""``scripts/group_dispersion_check.py``: designs, seeds, parents, the unit and its within-group arm, merge,
reproduction and decision.

Units run at fast CV knobs and few permutations; the decision runs on synthetic tables over the real settings.
No real data are read.
"""

from __future__ import annotations

import json
import shutil
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


gd = _script("group_dispersion_check")
ses = _script("stage_effect_scaling_check")
ubc = gd.ubc
lsc = gd.lsc
ssc = gd.ssc

FAST_CV = {"cv1_splits": 3, "cv2_splits": 4, "n_repeats": 1, "max_components": 6, "random_state": 1203}
PERMS = 39
COMPARED = ("n_stages", "n_samples", "p_dmp", "baseline_continuity", "magnitude_kind", "construction",
            "delta_methyl", "delta_expr", "delta_protein")


@pytest.fixture(scope="module")
def core():
    return gd.profile_settings("core")


@pytest.fixture(scope="module")
def columns():
    return gd.profile_settings("columns")


@pytest.fixture(scope="module")
def committed():
    _, settings = gd.parent_binding()
    return settings


def _without_balance(generator):
    return {key: value for key, value in generator.items() if key not in ("group_ratio", "group_dispersion_scale")}


def _params(settings, design, cell="anchor", replicate=0):
    return ubc.stage_span_generator_params(lsc.Unit(design, cell, replicate), settings["unit_settings"])


# ── 2.1 Settings, units and parents ──────────────────────────────────────────


def test_unit_counts_per_profile(core, columns) -> None:
    assert len(core["unit_settings"]["designs"]) == 18
    assert len(gd.profile_units(core)) == 18_000
    assert len(columns["unit_settings"]["designs"]) == 28
    assert len(gd.profile_units(columns)) == 14_000
    for settings in (core, columns):
        assert settings["unit_settings"]["permutations"] == 999
        assert set(settings["unit_settings"]["replicates"].values()) == {500}
        assert settings["unit_settings"]["recorded_spaces"] == ["pca_matched"]


def test_both_parents_are_bound_in_both_profiles(core, columns) -> None:
    expected = {
        f"results/{run}/{name}": ssc.sha256_file(gd.REPO / "results" / run / name)
        for run in ("univariate-baseline-2026-11-23/stage_span", "stage-effect-scaling-2026-11-24/band")
        for name in gd.PARENT_FILE_NAMES
    }
    assert core["parent_files"] == expected
    assert columns["parent_files"] == expected
    assert core["parents"]["univariate_baseline"]["dir"] == "results/univariate-baseline-2026-11-23/stage_span"
    assert core["parents"]["stage_effect_band"]["dir"] == "results/stage-effect-scaling-2026-11-24/band"


def _fake_record(unit, units, *, rejections=None, status="completed", rank=2, within="consistent", within_rejects=None):
    fields = ubc.unit_fields(unit, units)
    base = {**fields, "settings_hash": lsc.settings_hash(units), "status": status,
            "runtime_seconds": {"unit": 1.0, "fit": 0.5, "within_group": 0.2}}
    if status != "completed":
        return {**base, "error_type": "RuntimeError", "error_message": "patched",
                "within_group": {"status": "skipped", "reason": "the parent's unit failed"}}
    rejections = rejections or {}
    arms = {arm: {"p_value": 0.01 if rejections.get(arm) else 0.5, "rejects": bool(rejections.get(arm, False))}
            for arm in ubc.ARMS}
    if within == "failed":
        within_block = {"status": "failed", "error_type": "RuntimeError", "error_message": "patched"}
    else:
        within_rejects = within_rejects or {}
        within_block = {"status": "completed", "strata": "group", "rank": rank,
                        "same_statistics": within == "consistent",
                        "rejects": {s: bool(within_rejects.get(s, False)) for s in gd.STATISTICS}}
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
    }


def _small(profile="core", **overrides):
    overrides.setdefault("designs", [gd.RERUN_DESIGN, "n200_r0.8_s2"] if profile == "core"
                         else ["n80_r0.8", "n80_r0.2_s2"])
    return gd.profile_settings(profile, replicates=2, permutations=PERMS, **overrides)


@pytest.fixture()
def fake_units(monkeypatch):
    monkeypatch.setattr(gd, "evaluate_unit", lambda unit, units: _fake_record(unit, units))


def test_resume_skips_complete_units(tmp_path, monkeypatch) -> None:
    settings = _small()
    calls = []

    def fake(unit, units):
        calls.append(ubc.unit_key(unit))
        return _fake_record(unit, units)

    monkeypatch.setattr(gd, "evaluate_unit", fake)
    assert gd.run_shard(tmp_path, settings, 0, 1, limit=3) == 3
    assert gd.run_shard(tmp_path, settings, 0, 1) == 8 - 3
    assert len(calls) == 8 and len(set(calls)) == 8
    assert gd.run_shard(tmp_path, settings, 0, 1) == 0
    summary = gd.merge(tmp_path, "core")
    assert summary["complete_units"] == 8 and summary["missing_units"] == 0


def test_shards_are_round_robin(tmp_path, fake_units) -> None:
    settings = _small()
    assert [gd.run_shard(tmp_path, settings, index, 3) for index in range(3)] == [3, 3, 2]
    assert gd.merge(tmp_path, "core")["complete_units"] == 8


def test_a_resume_under_other_settings_is_refused_naming_the_keys(tmp_path, fake_units) -> None:
    gd.run_shard(tmp_path, _small(), 0, 1, limit=1)
    with pytest.raises(SystemExit, match=r"differing: unit_settings.permutations"):
        gd.run_shard(tmp_path, gd.profile_settings("core", designs=[gd.RERUN_DESIGN, "n200_r0.8_s2"], replicates=2,
                                                   permutations=19), 0, 1)
    with pytest.raises(SystemExit, match=r"unit_settings.replicates.control"):
        gd.run_shard(tmp_path, gd.profile_settings("core", designs=[gd.RERUN_DESIGN, "n200_r0.8_s2"], replicates=3,
                                                   permutations=PERMS), 0, 1)


def test_merge_refuses_records_written_under_other_settings(tmp_path, monkeypatch) -> None:
    monkeypatch.setattr(gd, "evaluate_unit",
                        lambda unit, units: {**_fake_record(unit, units), "settings_hash": "other"})
    gd.run_shard(tmp_path, _small(), 0, 1, limit=2)
    with pytest.raises(SystemExit, match=r"written under other settings"):
        gd.merge(tmp_path, "core")


def _parent_copies(tmp_path: Path) -> dict[str, Path]:
    out = {}
    for name, (directory, _) in gd.PARENTS.items():
        target = tmp_path / "parents" / name
        target.mkdir(parents=True)
        for file_name in gd.PARENT_FILE_NAMES:
            shutil.copy(directory / file_name, target / file_name)
        out[name] = target
    return out


@pytest.mark.parametrize("parent", sorted(gd.PARENTS))
@pytest.mark.parametrize("name", ["rate_table.csv", "cell_table.csv"])
def test_a_changed_parent_table_is_refused_by_name(tmp_path, monkeypatch, parent, name) -> None:
    monkeypatch.setattr(gd, "evaluate_unit", lambda unit, units: _fake_record(unit, units))
    dirs = _parent_copies(tmp_path)
    settings = _small(parent_dirs=dirs)
    assert {str(dirs[parent] / n) for n in gd.PARENT_FILE_NAMES} <= set(settings["parent_files"])
    out = tmp_path / "run"
    gd.run_shard(out, settings, 0, 1)
    gd.merge(out, "core")
    with (dirs[parent] / name).open("a") as handle:
        handle.write("\n")
    with pytest.raises(SystemExit, match=rf"{parent}/{name} has changed since the settings were written"):
        gd.run_shard(out, settings, 0, 1)
    with pytest.raises(SystemExit, match=rf"{name} has changed"):
        gd.merge(out, "core")


def test_decide_refuses_a_changed_parent(tmp_path, monkeypatch) -> None:
    monkeypatch.setattr(gd, "evaluate_unit", lambda unit, units: _fake_record(unit, units))
    dirs = _parent_copies(tmp_path)
    out = tmp_path / "run"
    for profile in gd.PROFILES:
        gd.run_shard(out, _small(profile, parent_dirs=dirs), 0, 1)
        gd.merge(out, profile)
    assert gd.decide(out)["verdict"] == "reproduction_failed"
    with (dirs["stage_effect_band"] / "cell_table.csv").open("a") as handle:
        handle.write("\n")
    with pytest.raises(SystemExit, match=r"stage_effect_band/cell_table.csv has changed"):
        gd.decide(out)


def test_a_changed_parent_settings_file_is_refused(tmp_path, fake_units) -> None:
    dirs = _parent_copies(tmp_path)
    settings = _small(parent_dirs=dirs)
    path = dirs["univariate_baseline"] / "settings.json"
    stored = json.loads(path.read_text())
    stored["settings"]["alpha"] = 0.1
    stored["hash"] = lsc.settings_hash(stored["settings"])
    path.write_text(json.dumps(stored))
    with pytest.raises(SystemExit, match=r"univariate_baseline/settings.json has changed since the settings"):
        gd.run_shard(tmp_path / "run", settings, 0, 1, limit=0)
    stored["hash"] = "0" * 64
    path.write_text(json.dumps(stored))
    with pytest.raises(SystemExit, match=r"settings.json fails its hash"):
        _small(parent_dirs=dirs)


def test_a_band_settings_file_with_another_stage_effect_is_refused(tmp_path) -> None:
    dirs = _parent_copies(tmp_path)
    path = dirs["stage_effect_band"] / "settings.json"
    stored = json.loads(path.read_text())
    generator = stored["settings"]["unit_settings"]["units"]["designs"]["band_n320"]["generator"]
    for key in ubc.STAGE_EFFECT_KEYS:
        generator[key] = 0.06
    stored["hash"] = lsc.settings_hash(stored["settings"])
    path.write_text(json.dumps(stored))
    with pytest.raises(SystemExit, match=r"band_n320 has stage effect \[0.06\], not the committed 0.05"):
        gd.profile_settings("columns", parent_dirs=dirs)


def test_unknown_designs_are_refused() -> None:
    with pytest.raises(SystemExit, match=r"Unknown core design\(s\) \['n80_r0.8'\]"):
        gd.profile_settings("core", designs=["n80_r0.8"])


# ── 2.2 The core builder ─────────────────────────────────────────────────────


def test_the_rerun_design_is_the_parents(core, committed) -> None:
    parent = ubc.stage_span_settings(designs=[gd.RERUN_DESIGN])
    units = core["unit_settings"]
    assert units["units"]["designs"][gd.RERUN_DESIGN] == committed["univariate_baseline"]["units"]["designs"]["n200"]
    assert [list(cell) for cell in units["units"]["cells"]] == [list(cell) for cell in gd.CORE_CELLS]
    for cell in gd.CORE_CELLS:
        assert list(cell) in committed["univariate_baseline"]["units"]["cells"]
    for unit in (lsc.Unit("n200", "anchor", 0), lsc.Unit("n200", "translation", 499), lsc.Unit("n200", "anchor", 250)):
        assert ubc.unit_fields(unit, units) == ubc.unit_fields(unit, parent)
        assert ubc.stage_span_generator_params(unit, units) == ubc.stage_span_generator_params(unit, parent)
    assert core["design_index"][gd.RERUN_DESIGN]["family"] is None


def test_core_has_the_seventeen_new_designs(core) -> None:
    index = core["design_index"]
    new = {design: meta for design, meta in index.items() if design != gd.RERUN_DESIGN}
    assert len(new) == 17
    expected = {("n200", share, None) for share in (0.65, 0.8)} | {
        ("n200", share, scale) for share in (0.5, 0.65, 0.8, 0.35, 0.2) for scale in (1.25, 1.5, 2.0)}
    assert {(meta["column"], meta["share"], meta["scale"]) for meta in new.values()} == expected
    assert "n200_r0.65" in new and "n200_r0.5_s1.25" in new and "n200_r0.2_s2" in new


def test_each_core_design_differs_from_n200_only_in_balance_and_dispersion(core, committed) -> None:
    reference = _without_balance(committed["univariate_baseline"]["units"]["designs"]["n200"]["generator"])
    for design, meta in core["design_index"].items():
        generator = core["unit_settings"]["units"]["designs"][design]["generator"]
        assert {key: value for key, value in generator.items()
                if key not in ("group_ratio", "group_dispersion_scale")} == reference
        assert generator["group_ratio"] == meta["share"]
        assert generator.get("group_dispersion_scale") == meta["scale"]


def test_a_core_design_measures_the_three_blocks(core) -> None:
    from motco.simulations.semisynthetic import generate_semisynthetic_trajectory

    dataset = generate_semisynthetic_trajectory(_params(core, "n200_r0.65_s1.25"))
    assert sum(getattr(dataset, layer).shape[1] for layer in ("methylation", "expression", "proteomics")) == 658
    assert dataset.truth["group_dispersion"] == {"scale": 1.25, "covariance_multiplier": 1.5625, "group": "B"}


def _cell_sizes(params) -> tuple[list[int], list[int]]:
    from motco.simulations import semisynthetic

    return semisynthetic._group_stage_sizes(semisynthetic._stage_sizes(params), params)


def test_the_share_sets_which_group_is_wider(core, columns) -> None:
    for settings in (core, columns):
        for design, meta in settings["design_index"].items():
            if meta["scale"] is None:
                continue
            sizes_a, sizes_b = _cell_sizes(_params(settings, design))
            if meta["share"] > 0.5:
                assert sum(sizes_b) < sum(sizes_a) and meta["family"] == "wider_smaller", design
            elif meta["share"] < 0.5:
                assert sum(sizes_b) > sum(sizes_a) and meta["family"] == "wider_larger", design
            else:
                # The parent's balanced designs round each stage's half (n200: 101 against 99).
                assert abs(sum(sizes_b) - sum(sizes_a)) <= len(sizes_a), design
                assert meta["family"] == "balanced_unequal", design


def test_new_names_and_datasets_are_fresh(core, columns) -> None:
    committed_names = set(ubc.STAGE_SPAN_DESIGNS) | set(ses.SATURATED_DESIGNS) | set(ses.BAND_DESIGNS)
    new = [design for settings in (core, columns) for design in settings["design_index"]
           if design != gd.RERUN_DESIGN]
    assert not set(new) & committed_names
    assert len(new) == len(set(new)) == 45
    committed_seeds = {lsc.generator_seed(2043, design, replicate)
                       for design in committed_names for replicate in range(500)}
    new_seeds = {lsc.generator_seed(2043, design, replicate) for design in new for replicate in range(500)}
    assert not committed_seeds & new_seeds


def test_the_two_cells_of_a_design_and_replicate_share_their_null_dataset(core) -> None:
    units = core["unit_settings"]
    for design in ("n200_r0.8_s2", "n200_r0.35_s1.5"):
        for replicate in (0, 7, 499):
            anchor = ubc.unit_fields(lsc.Unit(design, "anchor", replicate), units)
            translation = ubc.unit_fields(lsc.Unit(design, "translation", replicate), units)
            assert anchor["generator_seed"] == translation["generator_seed"]
            assert anchor["rrpp_seed"] != translation["rrpp_seed"]


# ── 2.3 The columns builder ──────────────────────────────────────────────────


def test_columns_designs_and_their_committed_generators(columns, committed) -> None:
    index = columns["design_index"]
    assert len(index) == 28
    for column in gd.COLUMN_SOURCES:
        mine = {design: meta for design, meta in index.items() if meta["column"] == column}
        assert len(mine) == 7
        assert {(meta["share"], meta["scale"]) for meta in mine.values()} == {(0.8, None)} | {
            (share, scale) for share in (0.5, 0.8, 0.2) for scale in (1.5, 2.0)}
        source = committed[gd.COLUMN_SOURCES[column]]
        source_units = source["unit_settings"] if "unit_settings" in source else source
        reference = _without_balance(source_units["units"]["designs"][column]["generator"])
        assert source_units["units"]["designs"][column]["generator"]["group_ratio"] == 0.5
        for design, meta in mine.items():
            generator = columns["unit_settings"]["units"]["designs"][design]["generator"]
            assert _without_balance(generator) == reference
            assert generator["group_ratio"] == meta["share"]
            assert generator.get("group_dispersion_scale") == meta["scale"]
    band = columns["unit_settings"]["units"]["designs"]["band_n320_r0.8"]["generator"]
    assert {band[key] for key in ubc.STAGE_EFFECT_KEYS} == {0.05}
    assert band["n_samples"] == 320


def test_every_cell_of_every_design_has_at_least_two_samples(core, columns) -> None:
    smallest = min(min(*_cell_sizes(_params(settings, design))[0], *_cell_sizes(_params(settings, design))[1])
                   for settings in (core, columns) for design in settings["design_index"])
    assert smallest >= 2
    sizes_a, sizes_b = _cell_sizes(_params(columns, "n80_r0.8_s2"))
    assert min(sizes_b) == 5 and sum(sizes_a) + sum(sizes_b) == 80


def test_columns_cells_are_the_anchor(columns) -> None:
    assert [list(cell) for cell in columns["unit_settings"]["units"]["cells"]] == [list(gd.COLUMN_CELLS[0])]
    assert gd.COLUMN_CELLS[0][0] == "anchor"
    assert columns["rerun_design"] is None


# ── 2.4 and 2.5 The unit and the within-group arm (real units, fast knobs) ───


def _fast(settings):
    settings["unit_settings"]["units"]["pls_cv_params"] = dict(FAST_CV)
    return settings


UNIT = lsc.Unit(gd.RERUN_DESIGN, "anchor", 1)


@pytest.fixture(scope="module")
def unit_records():
    import motco.simulations.evaluation as evaluation

    mine = _fast(_small())
    parent = ubc.stage_span_settings(designs=[gd.RERUN_DESIGN], replicates=2, control_replicates=2,
                                     permutations=PERMS)
    parent["units"]["pls_cv_params"] = dict(FAST_CV)
    real = evaluation.plsda_doubleCV
    calls = []

    def counting(*args, **kwargs):
        calls.append(1)
        return real(*args, **kwargs)

    with mock.patch.object(evaluation, "plsda_doubleCV", counting):
        record = gd.evaluate_unit(UNIT, mine["unit_settings"])
    return record, ubc.evaluate_unit(UNIT, parent), len(calls), mine


def test_a_unit_record_equals_the_parent_drivers(unit_records) -> None:
    record, parent, _, mine = unit_records
    assert record["status"] == "completed", record.get("error_message")
    assert record["settings_hash"] == lsc.settings_hash(mine["unit_settings"])

    def strip(r):
        return {key: value for key, value in r.items()
                if key not in ("settings_hash", "runtime_seconds", "within_group")}

    assert strip(record) == strip(parent)
    assert record["generator_seed"] == lsc.generator_seed(2043, UNIT.design, UNIT.replicate)
    assert set(record["spaces"]) == {"pca_matched"}
    assert {"parent_unit", "within_group", "unit"} <= set(record["runtime_seconds"])


def test_one_cross_validated_fit_per_unit(unit_records) -> None:
    assert unit_records[2] == 1


def test_a_failed_unit_is_recorded_and_the_shard_continues(tmp_path) -> None:
    settings = _fast(_small())

    def broken(*args, **kwargs):
        raise RuntimeError("patched scan failure")

    with mock.patch("motco.stats.interaction.interaction_feature_scan", broken):
        assert gd.run_shard(tmp_path, settings, 0, 4, limit=2) == 2
    records = [json.loads(line) for line in next((tmp_path / "core" / "records").glob("*.jsonl")).open()]
    assert len(records) == 2
    for record in records:
        assert record["status"] == "failed"
        assert (record["error_type"], record["error_message"]) == ("RuntimeError", "patched scan failure")
        assert record["within_group"]["status"] == "skipped"


@pytest.mark.parametrize(("profile", "unit"), [
    ("core", lsc.Unit("n200_r0.8_s2", "translation", 0)),
    ("columns", lsc.Unit("band_n320_r0.2_s1.5", "anchor", 1)),
])
def test_the_refit_reproduces_the_production_statistics(profile, unit) -> None:
    settings = _fast(gd.profile_settings(profile, replicates=2, permutations=PERMS, designs=[unit.design]))
    record = gd.evaluate_unit(unit, settings["unit_settings"])
    assert record["status"] == "completed", record.get("error_message")
    within = record["within_group"]
    assert within["status"] == "completed", within.get("error_message")
    assert within["same_statistics"] is True
    assert within["pair_statistics"] == record["production"]["pair_statistics"]
    assert within["rank"] == record["production"]["rank"]
    assert set(within["p_values"]) == set(gd.STATISTICS)


def _forced_inputs(unit, units, rank):
    from motco.simulations.evaluation import integrate_semisynthetic_dataset

    fields = ubc.unit_fields(unit, units)
    dataset, params = ubc.unit_inputs(unit, units, fields)
    return dataset, params, integrate_semisynthetic_dataset(dataset, gd.forced_rank_params(params, rank))


def test_one_stratum_reproduces_the_production_p_values(unit_records) -> None:
    record, _, _, mine = unit_records
    dataset, params, latent = _forced_inputs(UNIT, mine["unit_settings"], record["production"]["rank"])
    one = gd.strata_reading(dataset, latent.matrix, params, strata=np.zeros(len(dataset.metadata), dtype=int))
    assert one["pair_statistics"] == record["production"]["pair_statistics"]
    assert one["p_values"] == record["production"]["p_values"]
    grouped = gd.strata_reading(dataset, latent.matrix, params,
                                strata=gd.group_strata(dataset.metadata, params.group_col, params.stage_col))
    assert grouped["pair_statistics"] == record["production"]["pair_statistics"]
    assert grouped["p_values"] == record["within_group"]["p_values"]


def test_every_draw_keeps_each_residual_inside_its_group(monkeypatch) -> None:
    from motco.simulations.evaluation import SimulationEvaluationParams, build_simulation_trajectory_design
    from motco.stats import permutation
    from tests.test_residual_correction import _reduced_fit

    sizes = ((6, 5, 9), (4, 7, 5))
    groups = np.concatenate([[g] * sizes[g][s] for g in range(2) for s in range(3)])
    stages = np.concatenate([[s] * sizes[g][s] for g in range(2) for s in range(3)])
    metadata = pd.DataFrame({"group": np.where(groups == 0, "A", "B"), "stage": stages.astype(str)})
    rng = np.random.default_rng(4)
    Y = 2.0 * rng.normal(size=(3, 3))[stages] + rng.normal(size=(groups.size, 3)) * np.where(groups == 1, 2.0, 1.0)[
        :, None]
    params = SimulationEvaluationParams(permutations=PERMS, seed=5, orientation_estimator="net_displacement")
    dataset = SimpleNamespace(metadata=metadata)
    design = build_simulation_trajectory_design(metadata, group_col="group", stage_col="stage")
    y_hat, y_res = _reduced_fit(Y, np.asarray(design.model_reduced, dtype=float))

    def draws_under(strata):
        captured = []
        real = permutation.estimate_difference

        def capture(y_random, *args, **kwargs):
            captured.append(np.asarray(y_random, dtype=float).copy())
            return real(y_random, *args, **kwargs)

        with monkeypatch.context() as patch:
            patch.setattr(permutation, "estimate_difference", capture)
            gd.strata_reading(dataset, Y, params, strata=strata)
        out = []
        for y_random in captured:
            distances = (((y_random - y_hat)[:, None, :] - y_res[None, :, :]) ** 2).sum(axis=2)
            out.append(distances.argmin(axis=1))
        return out

    strata = gd.group_strata(metadata, "group", "stage")
    np.testing.assert_array_equal(strata, groups)
    within = draws_under(strata)
    assert len(within) == PERMS
    for draw in within:
        assert sorted(draw.tolist()) == list(range(groups.size))
        np.testing.assert_array_equal(groups[draw], groups)
    unrestricted = draws_under(np.zeros(groups.size, dtype=int))
    assert any(np.any(groups[draw] != groups) for draw in unrestricted)


def test_an_inconsistent_refit_is_marked_for_this_arm_only(unit_records, monkeypatch) -> None:
    record, _, _, mine = unit_records
    real = gd.strata_reading

    def shifted(*args, **kwargs):
        reading = real(*args, **kwargs)
        statistics = dict(reading["pair_statistics"])
        statistics["delta"] = statistics["delta"] + 1e-9
        return {**reading, "pair_statistics": statistics}

    monkeypatch.setattr(gd, "strata_reading", shifted)
    within = gd.within_group_reading(UNIT, mine["unit_settings"], record)
    assert within["status"] == "completed" and within["same_statistics"] is False
    marked = {**record, "within_group": within}
    rates, table = gd.within_group_rows([marked], mine["unit_settings"])
    row = next(r for r in table if (r["design"], r["cell"]) == (UNIT.design, UNIT.cell))
    assert (row["within_completed"], row["inconsistent"]) == (1, 1)
    assert all(r["count"] == 0 for r in rates if (r["design"], r["cell"]) == (UNIT.design, UNIT.cell))
    parent_rates = ubc.rate_rows([marked], mine["unit_settings"])
    assert next(r for r in parent_rates if (r["design"], r["cell"], r["arm"]) == (UNIT.design, UNIT.cell, "delta"))[
        "count"] == 1


# ── 3.1 Merge ────────────────────────────────────────────────────────────────


def _write_shard(out_dir: Path, settings, records, name="shard_000_of_001.jsonl") -> None:
    directory = gd.profile_dir(out_dir, settings["profile"])
    lsc.ensure_settings(directory, settings)
    (directory / "records").mkdir(parents=True, exist_ok=True)
    with (directory / "records" / name).open("a", encoding="utf-8") as handle:
        for record in records:
            handle.write(lsc._line(record))


def test_merge_writes_one_row_per_design_cell_and_arm_and_counts_failures(tmp_path) -> None:
    settings = _small()
    units = settings["unit_settings"]
    records = []
    for unit in gd.profile_units(settings):
        if (unit.design, unit.cell, unit.replicate) == ("n200_r0.8_s2", "anchor", 0):
            records.append(_fake_record(unit, units, status="failed"))
        elif (unit.design, unit.cell, unit.replicate) == ("n200_r0.8_s2", "anchor", 1):
            records.append(_fake_record(unit, units, rejections={"delta": True}, within="inconsistent"))
        elif (unit.design, unit.cell, unit.replicate) == ("n200_r0.8_s2", "translation", 0):
            records.append(_fake_record(unit, units, within="failed"))
        else:
            records.append(_fake_record(unit, units, within_rejects={"angle": True}))
    _write_shard(tmp_path, settings, records)
    # A duplicate written to a second shard file is counted once.
    _write_shard(tmp_path, settings, records[:1], name="shard_001_of_002.jsonl")
    summary = gd.merge(tmp_path, "core")
    assert summary["complete_units"] == 8 and summary["failed_units"] == 1
    directory = tmp_path / "core"
    rates = ubc._read_csv(directory / "rate_table.csv")
    assert len(rates) == 2 * 2 * len(gd.RATE_ARMS)
    keys = [(r["design"], r["cell"], r["arm"]) for r in rates]
    assert len(set(keys)) == len(keys)
    index = ubc._index(rates, "design", "cell", "arm")
    anchor = index[("n200_r0.8_s2", "anchor", "delta")]
    assert (int(anchor["count"]), int(anchor["rejections"]), int(anchor["failures"])) == (1, 1, 1)
    within = index[("n200_r0.8_s2", "anchor", "within_group_angle")]
    assert (int(within["count"]), int(within["failures"])) == (0, 2)
    translation = index[("n200_r0.8_s2", "translation", "within_group_angle")]
    assert (int(translation["count"]), int(translation["rejections"]), int(translation["failures"])) == (1, 1, 1)
    assert index[("n200", "anchor", "within_group_angle")]["rule"] == "p < alpha"
    table = ubc._index(ubc._read_csv(directory / "within_group_table.csv"), "design", "cell")
    row = table[("n200_r0.8_s2", "anchor")]
    assert {key: int(row[key]) for key in ("units", "parent_completed", "within_completed", "within_failed",
                                           "inconsistent")} == {
        "units": 2, "parent_completed": 1, "within_completed": 1, "within_failed": 0, "inconsistent": 1}
    assert int(table[("n200_r0.8_s2", "translation")]["within_failed"]) == 1
    accounting = json.loads((directory / "accounting.json").read_text())
    assert accounting["duplicate_records"] == 1
    assert accounting["within_group"]["inconsistent_by_cell"] == {"n200_r0.8_s2/anchor": 1}
    assert accounting["within_group"]["failed_by_cell"] == {"n200_r0.8_s2/translation": 1}
    assert accounting["failures_by_cell"] == {"n200_r0.8_s2/anchor": 1}
    assert (directory / "cell_table.csv").exists() and (directory / "reproduction.json").exists()


def test_the_columns_merge_writes_no_reproduction(tmp_path, fake_units) -> None:
    gd.run_shard(tmp_path, _small("columns"), 0, 1)
    assert gd.merge(tmp_path, "columns")["reproduction_passed"] is None
    assert not (tmp_path / "columns" / "reproduction.json").exists()


# ── 3.2 Reproduction ─────────────────────────────────────────────────────────


@pytest.fixture(scope="module")
def committed_rates():
    return lsc._read_csv(gd.UNIVARIATE_DIR / "rate_table.csv")


def _as_mine(committed_rates):
    return [dict(row) for row in committed_rates if row["design"] == "n200" and row["cell"] in gd.CONTROL_CELLS]


def test_the_committed_table_resolves_every_compared_row(core, committed_rates) -> None:
    result = gd.rerun_reproduction(_as_mine(committed_rates), committed_rates, core)
    assert result["passed"] is True
    assert result["rows_compared"] == 2 * len(ubc.ARMS)
    assert result["unresolved_rows"] == [] and result["differing"] == []


def test_the_reproduction_fails_and_lists_every_differing_row(core, committed_rates) -> None:
    mine = _as_mine(committed_rates)
    for row in mine:
        if (row["cell"], row["arm"]) in (("anchor", "delta"), ("translation", "scan")):
            row["rejections"] = str(int(row["rejections"]) + 1)
    mine.extend({**row, "arm": arm} for row in mine[:1] for arm in gd.WITHIN_ARMS)
    result = gd.rerun_reproduction(mine, committed_rates, core)
    assert result["passed"] is False
    assert result["rows_compared"] == 2 * len(ubc.ARMS)
    assert {(d["cell"], d["arm"]) for d in result["differing"]} == {("anchor", "delta"), ("translation", "scan")}
    for entry in result["differing"]:
        assert entry["reproduced"]["rejections"] == entry["committed"]["rejections"] + 1
    missing = gd.rerun_reproduction([row for row in mine if row["arm"] != "fdr"], committed_rates, core)
    assert missing["passed"] is False
    assert {(r["cell"], r["arm"], r["missing_from"]) for r in missing["unresolved_rows"]} == {
        ("anchor", "fdr", "this run"), ("translation", "fdr", "this run")}


# ── 3.3 and 3.4 Decision (synthetic tables over the real settings) ───────────


def _rates(settings, rejections=25, count=500):
    units = settings["unit_settings"]
    return [{"design": design, "cell": cell, "arm": arm, "count": count, "rejections": rejections}
            for design in units["designs"] for cell in ubc._cell_names(units) for arm in gd.RATE_ARMS]


def _cells(settings):
    units = settings["unit_settings"]
    return [{"design": design, "cell": cell, "completed": 500, "rank_q10": 2, "rank_median": 3, "rank_q90": 4,
             "auroc_median": 0.99} for design in units["designs"] for cell in ubc._cell_names(units)]


def _block(settings, **extra):
    return {"settings": settings, "accounting": {"missing_units": 0, "failed_units": 0, "within_group": {}},
            "rates": _rates(settings), "cells": _cells(settings), **extra}


@pytest.fixture()
def blocks(core, columns):
    return {"core": _block(core, reproduction={"passed": True, "rows_compared": 2 * len(ubc.ARMS)}),
            "columns": _block(columns)}


def _set(block, design, cell, arm, rejections, count=None):
    for row in block["rates"]:
        if (row["design"], row["cell"], row["arm"]) == (design, cell, arm):
            row["rejections"] = rejections
            if count is not None:
                row["count"] = count
            return
    raise KeyError((design, cell, arm))


def _decide(blocks):
    return gd.decision(blocks["core"], blocks["columns"])


def test_families_have_their_sizes_and_the_rerun_enters_none(blocks) -> None:
    result = _decide(blocks)
    assert {family: entry["m"] for family, entry in result["families"].items()} == gd.FAMILY_M_EXPECTED
    assert {family: (entry["allowed_exceedances"], entry["ceiling_at_median"])
            for family, entry in result["families"].items()} == {
        "unbalanced_equal": (2, 0.08), "balanced_unequal": (3, 0.082), "wider_smaller": (4, 0.084),
        "wider_larger": (4, 0.084)}
    observations = [row for family in gd.FAMILIES for row in result["criteria"]["size"][family]["observations"]]
    assert not [row for row in observations if row["design"] == gd.RERUN_DESIGN]
    keys = [(row["design"], row["cell"], row["arm"]) for row in observations]
    assert len(keys) == len(set(keys)) == 186
    assert {row["arm"] for row in observations} == set(gd.STATISTICS)


def test_no_excess_reads_size_holds(blocks) -> None:
    result = _decide(blocks)
    assert result["verdict"] == "size_holds"
    assert result["consequence"] == gd.CONSEQUENCES["size_holds"]
    assert all(entry["holds"] for entry in result["families"].values())
    assert result["conservative"] == []
    assert all(entry["none_failed"] for entry in result["map"])


def test_the_gate_comes_first_and_every_family_is_recorded(blocks) -> None:
    blocks["core"]["reproduction"] = {"passed": False, "rows_compared": 46,
                                      "differing": [{"cell": "anchor", "arm": "delta"}]}
    _set(blocks["core"], "n200_r0.8", "anchor", "delta", 60)
    result = _decide(blocks)
    assert result["verdict"] == "reproduction_failed"
    assert set(result["families"]) == set(gd.FAMILIES)
    assert result["families"]["unbalanced_equal"]["holds"] is False
    assert result["failing"]["reproduction"]["differing"] == [{"cell": "anchor", "arm": "delta"}]


def test_unbalance_alone_reads_unbalance_excess(blocks) -> None:
    _set(blocks["columns"], "n500_r0.8", "anchor", "angle", 60)
    _set(blocks["core"], "n200_r0.5_s2", "anchor", "delta", 60)
    result = _decide(blocks)
    assert result["verdict"] == "unbalance_excess"
    assert result["families"]["unbalanced_equal"]["holds"] is False
    assert result["families"]["balanced_unequal"]["holds"] is False
    assert set(result["failing"]) == {"unbalanced_equal", "balanced_unequal"}


def test_balanced_designs_losing_size_read_dispersion_excess(blocks) -> None:
    _set(blocks["core"], "n200_r0.5_s2", "translation", "delta", 60)
    _set(blocks["core"], "n200_r0.8_s2", "anchor", "delta", 80)
    result = _decide(blocks)
    assert result["verdict"] == "dispersion_excess"
    assert result["families"]["wider_smaller"]["holds"] is False


@pytest.mark.parametrize("design", ["n200_r0.8_s1.5", "band_n320_r0.2_s2"])
def test_only_unbalanced_designs_losing_size_read_unbalanced_dispersion_excess(blocks, design) -> None:
    block = blocks["core"] if design.startswith("n200_") else blocks["columns"]
    _set(block, design, "anchor", "shape", 60)
    result = _decide(blocks)
    assert result["verdict"] == "unbalanced_dispersion_excess"
    family = "wider_smaller" if "_r0.8_" in design else "wider_larger"
    assert result["families"][family]["holds"] is False
    assert [row["design"] for row in result["failing"][family]] == [design]


def test_exceedances_within_the_allowed_count_hold(blocks) -> None:
    for design in ("n200_r0.5_s1.25", "n200_r0.5_s1.5", "n200_r0.5_s2"):
        _set(blocks["core"], design, "anchor", "delta", 36)
    assert _decide(blocks)["families"]["balanced_unequal"]["holds"] is True
    _set(blocks["columns"], "n80_r0.5_s2", "anchor", "delta", 36)
    result = _decide(blocks)
    assert result["families"]["balanced_unequal"]["observed_exceedances"] == 4
    assert result["verdict"] == "dispersion_excess"


def test_an_unavailable_reading_fails_its_family(blocks) -> None:
    _set(blocks["columns"], "n80_r0.2_s2", "anchor", "delta", 0, count=0)
    result = _decide(blocks)
    assert result["families"]["wider_larger"]["holds"] is False
    assert result["families"]["wider_larger"]["n_unavailable"] == 1
    assert result["verdict"] == "unbalanced_dispersion_excess"
    assert [row["status"] for row in result["failing"]["wider_larger"]] == ["unavailable"]


def _map_entry(result, statistic, column, share):
    return next(e for e in result["map"] if (e["statistic"], e["column"], e["share"]) == (statistic, column, share))


def test_the_map_gives_the_smallest_failing_scale(blocks) -> None:
    _set(blocks["core"], "n200_r0.5_s1.5", "anchor", "delta", 60)
    _set(blocks["core"], "n200_r0.5_s2", "translation", "delta", 60)
    _set(blocks["core"], "n200_r0.5_s1.25", "anchor", "delta", 36)  # exceeds the bound, within the ceiling
    result = _decide(blocks)
    entry = _map_entry(result, "delta", "n200", 0.5)
    assert entry["first_failing_scale"] == 1.5 and entry["none_failed"] is False
    assert entry["by_cell"] == {"anchor": 1.5, "translation": 2.0}
    assert entry["scales_read"] == [1.25, 1.5, 2.0] and entry["family"] == "balanced_unequal"
    quiet = _map_entry(result, "angle", "n200", 0.5)
    assert quiet["first_failing_scale"] is None and quiet["none_failed"] is True
    assert quiet["largest_scale_read"] == 2.0


def test_the_map_covers_every_statistic_column_and_share(blocks) -> None:
    result = _decide(blocks)
    cells = {(e["statistic"], e["column"], e["share"]) for e in result["map"]}
    expected = {(s, "n200", share) for s in gd.STATISTICS for share in gd.CORE_SHARES} | {
        (s, column, share) for s in gd.STATISTICS for column in gd.COLUMN_SOURCES for share in gd.COLUMN_SHARES}
    assert cells == expected
    entry = _map_entry(result, "shape", "n80", 0.8)
    assert entry["scales_read"] == [1.5, 2.0]
    assert entry["equal_dispersion"] == {"anchor": "holds"}
    assert _map_entry(result, "shape", "n80", 0.5)["equal_dispersion"] is None


def test_the_map_order_survives_stored_settings(blocks) -> None:
    stored = {name: {**block, "settings": json.loads(json.dumps(block["settings"], sort_keys=True))}
              for name, block in blocks.items()}
    result = _decide(stored)
    order = list(dict.fromkeys(entry["column"] for entry in result["map"] if entry["statistic"] == "delta"))
    assert order == ["n200", "n80", "n500", "n200_delta0.5", "band_n320"]
    shares = [entry["share"] for entry in result["map"] if (entry["statistic"], entry["column"]) == ("delta", "n200")]
    assert shares == [0.5, 0.65, 0.8, 0.35, 0.2]
    assert result["map"] == _decide(blocks)["map"]


def test_conservative_readings_are_flagged_without_failing(blocks) -> None:
    _set(blocks["core"], "n200_r0.2_s2", "anchor", "delta", 15)
    _set(blocks["core"], "n200_r0.2_s2", "translation", "delta", 16)
    result = _decide(blocks)
    assert result["verdict"] == "size_holds"
    flagged = [(row["design"], row["cell"], row["arm"]) for row in result["conservative"]]
    assert flagged == [("n200_r0.2_s2", "anchor", "delta")]
    assert gd.conservative_limit(500) == pytest.approx(0.05 - 2 * (0.05 * 0.95 / 500) ** 0.5)


def test_recorded_arms_are_labelled_and_leave_the_verdict_and_map(blocks) -> None:
    _set(blocks["core"], "n200_r0.5_s2", "anchor", "delta", 60)
    before = _decide(blocks)
    for block in blocks.values():
        for row in block["rates"]:
            if row["arm"] in (*gd.WITHIN_ARMS, "production", "scan", "fdr", "total", "span", "complement",
                              "pca_matched"):
                row["rejections"] = 200
    after = _decide(blocks)
    assert after["verdict"] == before["verdict"] == "dispersion_excess"
    assert after["map"] == before["map"]
    assert after["families"] == before["families"]
    assert set(after["recorded_arms"]) == set(gd.RECORDED_ARMS)
    for arm in gd.RECORDED_ARMS:
        assert before["recorded_arms"][arm]["families"]["wider_smaller"]["label"] == "holds"
        assert after["recorded_arms"][arm]["families"]["wider_smaller"]["label"] == "exceeds"
        assert set(after["recorded_arms"][arm]["families"]) == set(gd.FAMILIES)
    assert after["recorded_arms"]["production"]["families"]["unbalanced_equal"]["m"] == 8
    assert after["recorded_arms"]["within_group_delta"]["rule"] == "p < alpha"


def test_fit_and_rates_are_reported_per_design_and_cell(blocks) -> None:
    result = _decide(blocks)
    assert len(result["fit"]) == 18 * 2 + 28
    entry = next(e for e in result["fit"] if (e["design"], e["cell"]) == ("n200_r0.8_s2", "translation"))
    assert (entry["rank_median"], entry["auroc_median"], entry["share"], entry["scale"]) == (3.0, 0.99, 0.8, 2.0)
    rates = next(e for e in result["rates"] if (e["design"], e["cell"]) == ("band_n320_r0.8", "anchor"))
    assert set(rates["rates"]) == set(gd.RATE_ARMS)
    assert rates["rates"]["within_group_shape"] == {"rate": 0.05, "count": 500, "rejections": 25}
