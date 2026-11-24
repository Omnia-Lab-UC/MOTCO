"""``scripts/univariate_baseline_check.py``: profiles, units, arms, merge, reproduction and the decision.

Units run at fast CV knobs and few permutations on the smallest designs; the
decision runs on synthetic tables. No real data are read.
"""

from __future__ import annotations

import dataclasses
import json
import os
import shutil
import subprocess
import sys
from pathlib import Path
from unittest import mock

import numpy as np
import pytest


def _script(name: str):
    sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
    try:
        return __import__(name)
    finally:
        sys.path.pop(0)


ubc = _script("univariate_baseline_check")
isb = _script("interaction_scan_blocks_check")
psc = _script("predictive_subspace_check")
isc = ubc.isc
lsc = ubc.lsc
ssc = ubc.ssc
slc = ubc.slc

REPO = Path(__file__).resolve().parents[1]
FAST_CV = {"cv1_splits": 3, "cv2_splits": 4, "n_repeats": 1, "max_components": 6, "random_state": 1203}
PERMS = 39


# ── Profiles and unit counts ─────────────────────────────────────────────────


@pytest.fixture(scope="module")
def full_settings():
    return {
        "parent": ubc.parent_settings(),
        "stage_span": ubc.stage_span_settings(),
        "surgery": ubc.surgery_settings(),
    }


def test_unit_counts_per_profile(full_settings) -> None:
    counts = {name: len(ubc.profile_units(settings)) for name, settings in full_settings.items()}
    assert counts == {"parent": 4_400, "stage_span": 12_800, "surgery": 16_500}
    signal = {name: len(ubc.signal_cells(settings)) for name, settings in full_settings.items()}
    assert signal["stage_span"] + signal["surgery"] == ubc.SIGNAL_CELLS_EXPECTED == 71
    for name, settings in full_settings.items():
        units = ubc.profile_units(settings)
        assert len({ubc.unit_key(unit) for unit in units}) == len(units)
        assert settings["permutations"] == 999


def test_parent_units_are_the_scan_validations(full_settings) -> None:
    settings = full_settings["parent"]
    committed = json.loads((ubc.SCAN_DIR / "settings.json").read_text())["settings"]
    comparison = settings["scan_validation"]["span_split"]["comparison"]
    stored = committed["span_split"]["comparison"]
    assert settings["designs"] == ["n200", "seaad"]
    for design in settings["designs"]:
        assert comparison["designs"][design] == stored["designs"][design]
    assert comparison["cells"] == stored["cells"]
    assert comparison["base_seed"] == stored["base_seed"] == 2035
    assert settings["scan_validation"]["replicates"] == committed["replicates"] == {"control": 500, "other": 100}
    units = ubc.profile_units(settings)
    parent_units = [unit for unit in isc.simulation_units(committed) if unit.design in settings["designs"]]
    assert units == parent_units
    for unit in (units[0], units[1234], units[-1]):
        fields = ubc.unit_fields(unit, settings)
        assert fields["generator_seed"] == lsc.generator_seed(stored["base_seed"], unit.design, unit.replicate)
        assert fields["rrpp_seed"] == lsc.rrpp_seed(stored["base_seed"], unit.design, unit.cell, unit.replicate)
        cell = lsc.cell_spec(stored, unit.cell)
        assert (fields["trajectory_mode"], fields["effect_size"], fields["span_alignment"], fields["span_support"]) == (
            cell["mode"], float(cell["effect"]), cell["phi"], cell["support"])
    assert set(settings["parent_files"]) == {
        "results/interaction-scan-2026-11-09/simulation/settings.json",
        "results/interaction-scan-2026-11-09/simulation/comparator_table.csv",
    }


def _scan_copy(tmp_path: Path) -> Path:
    target = tmp_path / "scan" / "simulation"
    target.mkdir(parents=True)
    for name in ubc.PARENT_FILE_NAMES:
        shutil.copy(ubc.SCAN_DIR / name, target / name)
    return target


def test_a_changed_parent_settings_file_is_refused_by_name(tmp_path) -> None:
    scan = _scan_copy(tmp_path)
    stored = json.loads((scan / "settings.json").read_text())
    stored["settings"]["alpha"] = 0.1
    stored["hash"] = lsc.settings_hash(stored["settings"])
    (scan / "settings.json").write_text(json.dumps(stored))
    with pytest.raises(SystemExit, match=r"scan/simulation/settings.json differs .*differing: alpha\."):
        ubc.parent_settings(scan_dir=scan)
    stored["hash"] = "0" * 64
    (scan / "settings.json").write_text(json.dumps(stored))
    with pytest.raises(SystemExit, match=r"settings.json fails its hash"):
        ubc.parent_settings(scan_dir=scan)


def test_a_changed_comparator_table_is_refused_by_name(tmp_path) -> None:
    scan = _scan_copy(tmp_path)
    settings = ubc.parent_settings(scan_dir=scan, designs=["seaad"], replicates=1, control_replicates=1,
                                   permutations=9)
    ubc.check_parent_files(settings)
    with (scan / "comparator_table.csv").open("a") as handle:
        handle.write("\n")
    with pytest.raises(SystemExit, match=r"comparator_table.csv has changed since the settings were written"):
        ubc.check_parent_files(settings)
    with pytest.raises(SystemExit, match=r"comparator_table.csv has changed"):
        ubc.run_shard(tmp_path / "run", settings, 0, 1, limit=0)
    out = tmp_path / "run2"
    (scan / "comparator_table.csv").write_bytes((ubc.SCAN_DIR / "comparator_table.csv").read_bytes())
    ubc.run_shard(out, settings, 0, 1, limit=0)
    (scan / "comparator_table.csv").write_text("changed\n")
    with pytest.raises(SystemExit, match=r"comparator_table.csv has changed"):
        ubc.merge(out, "parent")


def test_stage_span_designs_differ_only_in_size_and_stage_effect(full_settings) -> None:
    settings = full_settings["stage_span"]
    designs = settings["units"]["designs"]
    assert list(designs) == ["n80", "n200", "n500", "n200_delta0.5"]
    allowed = {"n_samples", *ubc.STAGE_EFFECT_KEYS}
    for first in designs.values():
        assert first["layers"] is None and "span_layers" not in first["generator"]
        assert first["generator"]["construction"] == "stage_span" and first["generator"]["n_stages"] == 3
        for second in designs.values():
            a, b = first["generator"], second["generator"]
            assert {key for key in set(a) | set(b) if a.get(key) != b.get(key)} <= allowed
    assert {name: (spec["generator"]["n_samples"], spec["generator"]["delta_expr"]) for name, spec in designs.items()} \
        == {"n80": (80, 2.0), "n200": (200, 2.0), "n500": (500, 2.0), "n200_delta0.5": (200, 0.5)}
    assert all(spec["generator"]["delta_methyl"] == spec["generator"]["delta_protein"]
               == spec["generator"]["delta_expr"] for spec in designs.values())


def test_stage_span_cells_and_scales(full_settings) -> None:
    settings = full_settings["stage_span"]
    cells = {cell[0]: cell[1:] for cell in settings["units"]["cells"]}
    assert len(cells) == 13
    assert cells["anchor"] == ["none", 0.0, None, None, "planted"]
    assert cells["translation"] == ["translation", 1.0, None, None, "planted"]
    for effect in (0.05, 0.1, 0.25):
        assert cells[f"magnitude_e{effect:g}"] == ["magnitude", effect, None, None, "planted"]
    for effect in (0.1, 0.25, 0.5):
        assert cells[f"orient_e{effect:g}_phi0"] == ["orientation", effect, 0.0, None, "planted"]
        assert cells[f"shape_e{effect:g}"] == ["shape", effect, None, None, "planted"]
    assert cells["orient_e0.5_phi90"] == ["orientation", 0.5, 90.0, None, None]
    assert cells["orient_e0.5_phi90_k20"] == ["orientation", 0.5, 90.0, 20, None]
    units = ubc.profile_units(settings)
    for design in settings["designs"]:
        for cell in cells:
            replicates = sorted(u.replicate for u in units if (u.design, u.cell) == (design, cell))
            assert replicates == list(range(500 if cell in ("anchor", "translation") else 200))
    unit = lsc.Unit("n80", "orient_e0.5_phi90", 3)
    params = ubc.stage_span_generator_params(unit, settings)
    assert (params.span_scale, params.span_alignment, params.seed) == (None, 90.0, lsc.generator_seed(2043, "n80", 3))
    params = ubc.stage_span_generator_params(lsc.Unit("n80", "shape_e0.25", 3), settings)
    assert (params.span_scale, params.group_effect_size, params.trajectory_mode) == ("planted", 0.25, "shape")


def _seed_set(base_seed: int, designs, replicates: int) -> set[int]:
    return {lsc.generator_seed(base_seed, design, r) for design in designs for r in range(replicates)}


def test_stage_span_is_a_fresh_seed_family(full_settings) -> None:
    settings = full_settings["stage_span"]
    ours = _seed_set(2043, settings["designs"], 500)
    assert len(ours) == len(settings["designs"]) * 500
    # The comparison, the split-test validation and the scan validation share base seed 2035 (replicates 0-499);
    # the predictive-subspace pilot is base seed 2036.
    comparison = _seed_set(lsc.BASE_SEED, lsc.DESIGNS, 500)
    predictive = _seed_set(psc.BASE_SEED, psc.SIMULATION_DESIGNS, 500)
    assert lsc.BASE_SEED == 2035 and psc.BASE_SEED == 2036
    assert not ours & comparison
    assert not ours & predictive
    assert not ours & {lsc.calibration_seed(lsc.BASE_SEED, r) for r in range(500)}


def test_cells_of_one_design_and_replicate_share_their_null_dataset(full_settings) -> None:
    from motco.simulations.reference import load_reference
    from motco.simulations.semisynthetic import generate_semisynthetic_trajectory

    settings = full_settings["stage_span"]
    reference = load_reference()

    def dataset(cell: str):
        return generate_semisynthetic_trajectory(
            ubc.stage_span_generator_params(lsc.Unit("n80", cell, 7), settings), reference=reference
        )

    anchor = dataset("anchor")
    for cell in ("magnitude_e0.1", "orient_e0.5_phi90"):
        other = dataset(cell)
        assert anchor.metadata.equals(other.metadata)
        group_a = (anchor.metadata["group"] == "A").to_numpy()
        for layer in ("methylation", "expression", "proteomics"):
            np.testing.assert_array_equal(getattr(anchor, layer).to_numpy()[group_a],
                                          getattr(other, layer).to_numpy()[group_a])
            assert not np.array_equal(getattr(anchor, layer).to_numpy()[~group_a],
                                      getattr(other, layer).to_numpy()[~group_a])
    assert anchor.truth["span_layers"] == ["methylation", "expression", "proteomics"]


def _study_cell(config: Path, column: str, mode: str, effect: float):
    from motco.simulations.study import enumerate_study, load_study_config

    for cell in enumerate_study(load_study_config(config)).cells:
        if cell.phase not in ("power_primary", "power_design") or slc.column_of(cell) != column:
            continue
        cell_effect = float(cell.generator_params.group_effect_size)
        cell_mode = "none" if cell_effect == 0 else cell.generator_params.trajectory_mode
        if cell_mode == mode and cell_effect == float(effect):
            return cell
    raise AssertionError(f"no study cell for {column} {mode} {effect}")


def test_surgery_units_are_the_study_engines(full_settings) -> None:
    from motco.simulations.grid import derive_replicate_seed

    settings = full_settings["surgery"]
    assert settings["designs"] == list(ubc.SURGERY_COLUMNS)
    assert settings["n_samples"] == {"n80_three_blocks": 80, "n200_three_blocks": 200, "n500_three_blocks": 500}
    cells = {cell[0]: cell[1:] for cell in settings["cells"]}
    assert len(cells) == 11 and cells["anchor"] == ["none", 0.0, "geometry"]
    assert cells["translation"] == ["translation", 1.0, "geometry"]
    assert {name for name, (_, _, source) in cells.items() if source == "magnitude"} == {
        "magnitude_e0.02", "magnitude_e0.05", "magnitude_e0.1"}
    sampled = (("n80_three_blocks", "anchor", 0), ("n200_three_blocks", "magnitude_e0.05", 17),
               ("n500_three_blocks", "orientation_e1", 499), ("n80_three_blocks", "shape_e0.25", 3),
               ("n200_three_blocks", "translation", 250))
    for column, cell, replicate in sampled:
        mode, effect, source = cells[cell]
        study = _study_cell(ubc.SURGERY_SOURCES[source]["config"], column, mode, effect)
        generator, evaluation = ubc.surgery_params(lsc.Unit(column, cell, replicate), settings)
        replicate_seed = derive_replicate_seed(study, replicate)
        assert generator == dataclasses.replace(study.generator_params, seed=replicate_seed)
        expected_seed = study.evaluation_params.seed if study.evaluation_params.seed is not None else replicate_seed
        assert evaluation == dataclasses.replace(study.evaluation_params, seed=expected_seed, span_capture=True,
                                                 include_null_distributions=True)
        assert evaluation.integration_params.get("layers") is None
        fields = ubc.unit_fields(lsc.Unit(column, cell, replicate), settings)
        assert (fields["generator_seed"], fields["rrpp_seed"]) == (replicate_seed, expected_seed)
    anchor = _study_cell(ubc.SURGERY_SOURCES["geometry"]["config"], "n80_three_blocks", "none", 0.0)
    generator, _ = ubc.surgery_params(lsc.Unit("n80_three_blocks", "anchor", 4), settings)
    assert generator == dataclasses.replace(anchor.generator_params, seed=derive_replicate_seed(anchor, 4))


def test_a_changed_surgery_parent_is_refused_by_name(tmp_path) -> None:
    sources = {}
    for name, paths in ubc.SURGERY_SOURCES.items():
        sources[name] = {}
        for kind, path in paths.items():
            target = tmp_path / name / path.name
            target.parent.mkdir(parents=True, exist_ok=True)
            shutil.copy(path, target)
            sources[name][kind] = target
    settings = ubc.surgery_settings(sources=sources, columns=["n80_three_blocks"], replicates=1, control_replicates=1)
    ubc.check_parent_files(settings)
    with sources["geometry"]["operating"].open("a") as handle:
        handle.write("\n")
    with pytest.raises(SystemExit, match=r"geometry/design_point_operating.csv has changed since the settings"):
        ubc.check_parent_files(settings)


# ── Resume and settings binding (fast fake units) ────────────────────────────


def _fake_record(unit, settings, *, rejections=None, status="completed", **extra):
    fields = ubc.unit_fields(unit, settings)
    base = {**fields, "settings_hash": lsc.settings_hash(settings), "status": status,
            "runtime_seconds": {"unit": 1.0, "fit": 0.5}}
    if status != "completed":
        return {**base, "error_type": "RuntimeError", "error_message": "patched"}
    rejections = rejections or {}
    arms = {arm: {"p_value": 0.01 if rejections.get(arm) else 0.5, "rejects": bool(rejections.get(arm, False))}
            for arm in ubc.ARMS}
    return {
        **base,
        "blocks": ["methylation", "expression", "proteomics"],
        "n_features": 658,
        "n_population_null": extra.get("n_population_null", 100),
        "production": {"rank": extra.get("rank", 2), "cv_mean_auroc": 0.9,
                       "span_capture": {"group_difference": 0.5, "pooled_stage": 0.9, "random_share": 0.01}},
        "scan": {"hits": extra.get("hits", 0), "hits_population_null": extra.get("hits_null", 0)},
        "fdr": {"discoveries": extra.get("discoveries", 0), "discoveries_population_null": extra.get("disc_null", 0)},
        "spaces": {"concat": {"dimension": 658}, "pca_matched": {"dimension": extra.get("rank", 2)}},
        "arms": arms,
    }


def _small_stage_span(**overrides):
    overrides.setdefault("recorded_spaces", ubc.RECORDED_SPACES)
    return ubc.stage_span_settings(designs=["n80"], replicates=2, control_replicates=2, permutations=PERMS,
                                   **overrides)


def test_resume_skips_complete_units(tmp_path, monkeypatch) -> None:
    settings = _small_stage_span()
    calls = []

    def fake(unit, settings):
        calls.append(ubc.unit_key(unit))
        return _fake_record(unit, settings)

    monkeypatch.setattr(ubc, "evaluate_unit", fake)
    assert ubc.run_shard(tmp_path, settings, 0, 1, limit=5) == 5
    assert ubc.run_shard(tmp_path, settings, 0, 1) == 26 - 5
    assert len(calls) == 26 and len(set(calls)) == 26
    assert ubc.run_shard(tmp_path, settings, 0, 1) == 0
    summary = ubc.merge(tmp_path, "stage_span")
    assert summary["complete_units"] == 26 and summary["missing_units"] == 0
    merged = (tmp_path / "stage_span" / "merged.jsonl").read_text().splitlines()
    assert len(merged) == 26


def test_shards_are_round_robin(tmp_path, monkeypatch) -> None:
    settings = _small_stage_span()
    monkeypatch.setattr(ubc, "evaluate_unit", lambda unit, settings: _fake_record(unit, settings))
    counts = [ubc.run_shard(tmp_path, settings, index, 3) for index in range(3)]
    assert counts == [9, 9, 8]
    assert ubc.merge(tmp_path, "stage_span")["complete_units"] == 26


def test_a_resume_under_other_settings_is_refused_naming_the_keys(tmp_path, monkeypatch) -> None:
    monkeypatch.setattr(ubc, "evaluate_unit", lambda unit, settings: _fake_record(unit, settings))
    ubc.run_shard(tmp_path, _small_stage_span(), 0, 1, limit=1)
    with pytest.raises(SystemExit, match=r"differing: permutations"):
        ubc.run_shard(tmp_path, ubc.stage_span_settings(designs=["n80"], replicates=2, control_replicates=2,
                                                        permutations=19,
                                                        recorded_spaces=ubc.RECORDED_SPACES), 0, 1)
    with pytest.raises(SystemExit, match=r"differing: replicates.other"):
        ubc.run_shard(tmp_path, ubc.stage_span_settings(designs=["n80"], replicates=3, control_replicates=2,
                                                        permutations=PERMS,
                                                        recorded_spaces=ubc.RECORDED_SPACES), 0, 1)


def test_merge_refuses_records_written_under_other_settings(tmp_path, monkeypatch) -> None:
    settings = _small_stage_span()
    monkeypatch.setattr(ubc, "evaluate_unit", lambda unit, settings: {**_fake_record(unit, settings),
                                                                      "settings_hash": "other"})
    ubc.run_shard(tmp_path, settings, 0, 1, limit=2)
    with pytest.raises(SystemExit, match=r"written under other settings"):
        ubc.merge(tmp_path, "stage_span")


# ── The unit and its arms (real units, fast knobs) ───────────────────────────


def _fast_stage_span():
    settings = _small_stage_span()
    settings["units"]["pls_cv_params"] = dict(FAST_CV)
    return settings


def _fast_parent():
    settings = ubc.parent_settings(designs=["seaad"], replicates=2, control_replicates=2, permutations=PERMS,
                                   recorded_spaces=ubc.RECORDED_SPACES)
    settings["scan_validation"]["span_split"]["comparison"]["pls_cv_params"] = dict(FAST_CV)
    return settings


SPAN_UNIT = lsc.Unit("n80", "orient_e0.5_phi90_k20", 1)


@pytest.fixture(scope="module")
def stage_span_unit():
    import motco.simulations.evaluation as evaluation

    settings = _fast_stage_span()
    real = evaluation.plsda_doubleCV
    calls = []

    def counting(*args, **kwargs):
        calls.append(1)
        return real(*args, **kwargs)

    with mock.patch.object(evaluation, "plsda_doubleCV", counting):
        record = ubc.evaluate_unit(SPAN_UNIT, settings)
    return settings, record, len(calls)


@pytest.fixture(scope="module")
def parent_unit():
    settings = _fast_parent()
    return settings, ubc.evaluate_unit(lsc.Unit("seaad", "magnitude", 0), settings)


def test_one_cross_validated_fit_per_unit(stage_span_unit) -> None:
    _, record, calls = stage_span_unit
    assert record["status"] == "completed", record.get("error_message")
    assert calls == 1


def test_the_production_record_equals_a_plain_evaluation(stage_span_unit) -> None:
    from motco.simulations.evaluation import evaluate_semisynthetic_trajectory

    settings, record, _ = stage_span_unit
    fields = ubc.unit_fields(SPAN_UNIT, settings)
    dataset, params = ubc.unit_inputs(SPAN_UNIT, settings, fields)
    plain = evaluate_semisynthetic_trajectory(dataset, dataclasses.replace(params, include_null_distributions=False))
    production = record["production"]
    assert production["rank"] == int(plain.latent_matrix_metadata["selected_lv"])
    assert production["cv_mean_auroc"] == pytest.approx(plain.latent_matrix_metadata["cv_mean_auroc"])
    assert production["pair_statistics"] == pytest.approx(lsc._finite(plain.pair_statistics))
    assert production["p_values"] == lsc._finite(plain.p_values)
    assert record["rrpp_seed"] == lsc.rrpp_seed(2043, *SPAN_UNIT.key) == params.seed


def test_the_baseline_does_not_depend_on_the_fit(stage_span_unit) -> None:
    from motco.stats.interaction import interaction_feature_scan

    settings, record, _ = stage_span_unit
    fields = ubc.unit_fields(SPAN_UNIT, settings)
    dataset, params = ubc.unit_inputs(SPAN_UNIT, settings, fields)
    features, full, reduced = isb.simulation_inputs(dataset, params)
    scan = interaction_feature_scan(features.to_numpy(dtype=float), full, reduced, permutations=PERMS,
                                    seed=record["rrpp_seed"])
    assert record["scan"]["max_f"] == pytest.approx(scan.max_f, rel=1e-9)
    assert record["scan"]["p_value_global"] == scan.p_value_global
    assert record["scan"]["n_features"] == scan.n_features == 658
    fdr = ubc.fdr_reading(scan.f, scan.degenerate, scan.df_interaction, scan.df_residual)
    assert record["fdr"]["discoveries"] == fdr["discoveries"]
    assert record["fdr"]["min_adjusted"] == pytest.approx(fdr["min_adjusted"], rel=1e-6)


def test_the_scan_arm_is_the_packages_scan_on_the_split_inputs(parent_unit) -> None:
    from motco.simulations.evaluation import integrate_semisynthetic_dataset
    from motco.stats.interaction import interaction_feature_scan

    settings, record = parent_unit
    assert record["status"] == "completed", record.get("error_message")
    unit = lsc.Unit("seaad", "magnitude", 0)
    dataset, params = ubc.unit_inputs(unit, settings, ubc.unit_fields(unit, settings))
    latent = integrate_semisynthetic_dataset(dataset, params)
    features, _, full, reduced = ssc.split_inputs(dataset, latent, params)
    scan = interaction_feature_scan(features, full, reduced, permutations=PERMS, seed=record["rrpp_seed"])
    assert record["scan"]["max_f"] == scan.max_f
    assert record["scan"]["p_value_global"] == scan.p_value_global
    assert record["arms"]["scan"]["rejects"] == (scan.p_value_global <= 0.05)


def _arrays(value, path=""):
    if isinstance(value, dict):
        for key, item in value.items():
            yield from _arrays(item, f"{path}.{key}")
    elif isinstance(value, list):
        yield path, len(value)
        for item in value:
            yield from _arrays(item, path)


@pytest.mark.parametrize("which", ["stage_span", "parent"])
def test_the_record_carries_every_field_and_no_feature_array(which, stage_span_unit, parent_unit) -> None:
    record = stage_span_unit[1] if which == "stage_span" else parent_unit[1]
    assert record["status"] == "completed", record.get("error_message")
    for key in ("profile", "design", "cell", "trajectory_mode", "effect_size", "replicate", "generator_seed",
                "rrpp_seed", "settings_hash", "blocks", "n_population_null", "production", "split", "scan", "fdr",
                "spaces", "arms", "runtime_seconds"):
        assert key in record, key
    production = record["production"]
    assert {"rank", "cv_mean_auroc", "p_values", "p_min", "p_minp", "span_capture"} <= set(production)
    assert set(production["span_capture"]) == {"group_difference", "pooled_stage", "random_share"}
    assert {"hits", "hits_population_null", "p_value_global", "max_f"} <= set(record["scan"])
    assert {"discoveries", "discoveries_population_null", "n_degenerate", "min_adjusted"} <= set(record["fdr"])
    assert set(record["arms"]) == set(ubc.ARMS)
    for arm, reading in record["arms"].items():
        assert reading["rejects"] == ubc.rejects(reading["p_value"], ubc.RULES[arm]), arm
    assert record["arms"]["production"]["p_value"] == min(
        value for value in production["p_values"].values() if value is not None)
    assert all(length < 20 for _, length in _arrays(record)), list(_arrays(record))
    assert set(record["runtime_seconds"]) == {"fit", "evaluation", "split", "scan", "concat", "pca_matched", "unit"}


def test_the_recorded_spaces_measure_the_units_blocks(stage_span_unit, parent_unit) -> None:
    for (_, record, *_), blocks, features in ((stage_span_unit, ["methylation", "expression", "proteomics"], 658),
                                              (parent_unit, ["methylation", "expression"], 498)):
        assert record["blocks"] == blocks
        assert record["n_features"] == features
        assert record["spaces"]["concat"]["dimension"] == features
        assert record["spaces"]["pca_matched"]["dimension"] == record["production"]["rank"]
        assert record["split"]["total"]["dimension"] == features


def test_recorded_spaces_are_bound_and_a_left_out_space_is_unavailable() -> None:
    both = _small_stage_span()
    lean = _small_stage_span(recorded_spaces=["pca_matched"])
    assert both["recorded_spaces"] == ["concat", "pca_matched"] and lean["recorded_spaces"] == ["pca_matched"]
    assert lsc.settings_difference(both, lean) == ["recorded_spaces"]
    readings = ubc.arm_readings(
        {"p_values": {"delta": 0.01, "angle": 0.5, "shape": 0.5}, "p_min": 0.01, "p_minp": 0.02},
        {part: {"p_value_f": 0.5} for part in ("total", "span", "complement")},
        {"p_value_global": 0.5}, {"min_adjusted": 0.5},
        {"pca_matched": {"p_values": {"delta": 0.001}, "p_min": 0.001}},
    )
    assert all(readings[arm] == {"p_value": None, "rejects": None}
               for arm in ("concat", "concat_delta", "concat_angle", "concat_shape"))
    assert readings["pca_matched"] == {"p_value": 0.001, "rejects": True}
    assert readings["production"] == {"p_value": 0.01, "rejects": True}
    with pytest.raises(SystemExit, match=r"Unknown recorded space"):
        _small_stage_span(recorded_spaces=["pls_residual_pca"])
    args = ubc.build_parser().parse_args(["run", "--profile", "stage_span", "--out-dir", "x", "--designs", "n80",
                                          "--recorded-spaces", "none", "--reps", "2", "--control-reps", "2"])
    settings = ubc.settings_from_args(args)
    assert settings["recorded_spaces"] == [] and settings["designs"] == ["n80"]


def test_concat_is_left_out_by_default() -> None:
    """The approval gate dropped ``concat``: every profile and the CLI record ``pca_matched`` alone."""
    assert ubc.stage_span_settings(designs=["n80"], replicates=2, control_replicates=2)["recorded_spaces"] == [
        "pca_matched"]
    assert ubc.parent_settings(designs=["seaad"])["recorded_spaces"] == ["pca_matched"]
    assert ubc.surgery_settings(columns=["n80_three_blocks"])["recorded_spaces"] == ["pca_matched"]
    args = ubc.build_parser().parse_args(["run", "--profile", "surgery", "--out-dir", "x"])
    assert ubc.settings_from_args(args)["recorded_spaces"] == ["pca_matched"]


def test_recorded_space_parameters_change_only_the_method() -> None:
    from motco.simulations.evaluation import SimulationEvaluationParams

    params = SimulationEvaluationParams(integration_method="pls", integration_params={**FAST_CV, "layers": ["a"]},
                                        permutations=9, seed=4, span_capture=True, include_null_distributions=True)
    concat = ubc.recorded_space_params(params, "concat", rank=3)
    pca = ubc.recorded_space_params(params, "pca_matched", rank=3)
    assert dataclasses.replace(concat, integration_method="pls", integration_params=params.integration_params) == params
    assert concat.integration_params == {"layers": ["a"]}
    assert pca.integration_method == "pca" and pca.integration_params == {"n_components": 3, "layers": ["a"]}
    three = dataclasses.replace(params, integration_params=dict(FAST_CV))
    assert ubc.recorded_space_params(three, "concat", rank=3).integration_params == {}


def test_a_failed_arm_is_recorded_and_the_shard_continues(tmp_path) -> None:
    settings = _fast_stage_span()

    def broken(*args, **kwargs):
        raise RuntimeError("patched scan failure")

    with mock.patch("motco.stats.interaction.interaction_feature_scan", broken):
        assert ubc.run_shard(tmp_path, settings, 0, 13, limit=2) == 2
    records = [json.loads(line) for line in next((tmp_path / "stage_span" / "records").glob("*.jsonl")).open()]
    assert len(records) == 2
    for record in records:
        assert record["status"] == "failed"
        assert (record["error_type"], record["error_message"]) == ("RuntimeError", "patched scan failure")
        assert record["generator_seed"] == lsc.generator_seed(2043, record["design"], record["replicate"])
        assert {"fit", "evaluation", "split", "unit"} <= set(record["runtime_seconds"])


# ── Production readings ──────────────────────────────────────────────────────


def test_the_omnibus_rejects_at_alpha_over_three() -> None:
    assert ubc.rejects(ubc.omnibus_p({"delta": 0.012, "angle": 0.4, "shape": 0.7}), ubc.RULES["production"])
    p_values = {"delta": 0.030, "angle": 0.400, "shape": 0.700}
    assert not ubc.rejects(ubc.omnibus_p(p_values), ubc.RULES["production"])
    assert ubc.rejects(p_values["delta"], ubc.RULES["delta"])
    assert ubc.rejects(0.05 / 3, ubc.RULES["production"]) and not ubc.rejects(0.017, ubc.RULES["production"])
    assert ubc.omnibus_p({"delta": None, "angle": float("nan"), "shape": 0.2}) == 0.2
    assert ubc.omnibus_p({"delta": None}) is None
    assert ubc.rejects(None, ubc.RULES["production"]) is None


def test_each_arm_keeps_its_validated_rule() -> None:
    assert ubc.RULES["scan"] == "p <= alpha" and ubc.rejects(0.05, ubc.RULES["scan"])
    for arm in ("total", "span", "complement", "delta", "angle", "shape", "concat_delta", "pca_matched_shape"):
        assert ubc.RULES[arm] == "p < alpha" and not ubc.rejects(0.05, ubc.RULES[arm]), arm
    assert ubc.rejects(0.05, ubc.RULES["fdr"]) and not ubc.rejects(0.0501, ubc.RULES["fdr"])
    assert ubc.rejects(0.05, ubc.RULES["production_minp"])
    for arm in ("production", "concat", "pca_matched"):
        assert ubc.RULES[arm] == "min p <= alpha / 3", arm


def _minp_direct(observed, nulls):
    names = [s for s in ("delta", "angle", "shape") if s in nulls and np.isfinite(observed.get(s, np.nan))]
    columns = {s: [observed[s], *nulls[s]] for s in names}
    size = len(next(iter(columns.values())))

    def p_of(s, j):
        values = [v if np.isfinite(v) else -np.inf for v in columns[s]]
        return sum(1 for v in values if v >= values[j]) / size

    smallest = [min(p_of(s, j) for s in names) for j in range(size)]
    return sum(1 for value in smallest if value <= smallest[0]) / size


def test_minp_equals_a_direct_rank_computation() -> None:
    rng = np.random.default_rng(5)
    nulls: dict[str, list[float]] = {}
    for trial in range(5):
        nulls = {s: list(rng.gamma(2.0, size=49)) for s in ("delta", "angle", "shape")}
        nulls["angle"][3] = float("nan")
        observed = {s: float(rng.gamma(2.0) * (1.0 + trial / 2)) for s in ("delta", "angle", "shape")}
        assert ubc.minp_p_value(observed, nulls) == pytest.approx(_minp_direct(observed, nulls))
    observed = {"delta": 1.0, "angle": float("nan"), "shape": 2.0}
    assert ubc.minp_p_value(observed, nulls) == pytest.approx(_minp_direct(observed, nulls))
    assert ubc.minp_p_value({"delta": float("nan")}, nulls) is None


def test_minp_is_never_below_one_over_draws_plus_one() -> None:
    nulls = {s: [0.1] * 99 for s in ("delta", "angle", "shape")}
    assert ubc.minp_p_value({"delta": 50.0, "angle": 50.0, "shape": 50.0}, nulls) == pytest.approx(1 / 100)
    assert ubc.minp_p_value({"delta": 0.0, "angle": 0.0, "shape": 0.0}, nulls) == 1.0


def test_minp_on_one_statistic_is_its_empirical_p_value() -> None:
    from motco.simulations.evaluation import _empirical_p_value

    rng = np.random.default_rng(2)
    null = list(rng.normal(size=199))
    for observed in (-1.0, 0.5, 2.5):
        assert ubc.minp_p_value({"delta": observed}, {"delta": null}) == pytest.approx(
            _empirical_p_value(null, observed))


def test_minp_record_matches_the_units_draws(stage_span_unit) -> None:
    _, record, _ = stage_span_unit
    production = record["production"]
    assert 1 / (PERMS + 1) <= production["p_minp"] <= 1.0
    assert production["p_minp"] >= production["p_min"]


# ── The false-discovery arm ──────────────────────────────────────────────────


def test_the_fdr_arm_on_a_hand_computed_case() -> None:
    from scipy.stats import f as f_distribution

    f_values = np.array([12.0, 0.4, 7.5, float("nan"), 3.1, 5.2, 1.0, 9.8])
    degenerate = np.array([False, False, False, True, False, False, False, False])
    df1, df2 = 2, 50
    tested, p, adjusted = ubc.bh_adjusted(f_values, degenerate, df1, df2)
    assert tested.tolist() == [True, True, True, False, True, True, True, True]
    expected_p = f_distribution.sf(f_values[tested], df1, df2)
    np.testing.assert_allclose(p[tested], expected_p)
    m = len(expected_p)
    order = np.argsort(expected_p)
    ranked = expected_p[order] * m / np.arange(1, m + 1)
    hand = np.minimum.accumulate(ranked[::-1])[::-1].clip(max=1.0)
    expected_adjusted = np.empty(m)
    expected_adjusted[order] = hand
    np.testing.assert_allclose(adjusted[tested], expected_adjusted)
    assert np.isnan(adjusted[3])
    reading = ubc.fdr_reading(f_values, degenerate, df1, df2)
    assert reading["n_tested"] == 7 and reading["n_degenerate"] == 1
    assert reading["min_adjusted"] == pytest.approx(expected_adjusted.min())
    assert reading["discoveries"] == int((expected_adjusted <= 0.05).sum())
    assert ubc.rejects(reading["min_adjusted"], ubc.RULES["fdr"]) == bool(expected_adjusted.min() <= 0.05)
    mask = np.array([True, True, False, False, False, True, True, True])
    with_mask = ubc.fdr_reading(f_values, degenerate, df1, df2, null_mask=mask)
    assert with_mask["discoveries_population_null"] == int(((adjusted <= 0.05) & mask & tested).sum())


def test_fdr_rejects_where_the_maximum_does_not() -> None:
    """A hundred near-copies of one moderately significant feature among 900 null ones.

    Benjamini–Hochberg discovers every copy (each p-value is below 100 · 0.05 / 1000), while the maximum's
    permutation p-value pays for the search over the 900 null features.
    """

    import pandas as pd

    from motco.stats.design import get_model_matrix
    from motco.stats.interaction import interaction_feature_scan

    frame = pd.DataFrame(
        [(group, stage) for group in ("A", "B") for stage in range(3) for _ in range(10)], columns=["group", "stage"]
    ).astype(str)
    full = np.asarray(get_model_matrix(frame, "group", "stage", full=True), dtype=float)
    reduced = np.asarray(get_model_matrix(frame, "group", "stage", full=False), dtype=float)
    interaction = ((frame["group"] == "B") & (frame["stage"] == "2")).to_numpy(dtype=float)
    rng = np.random.default_rng(0)
    common, noise, jitter = rng.normal(size=60), rng.normal(size=(60, 900)), rng.normal(size=(60, 100))
    copies = (2.4 * interaction + common)[:, None] + 0.05 * jitter
    scan = interaction_feature_scan(np.hstack([copies, noise]), full, reduced, permutations=199, seed=3)
    reading = ubc.fdr_reading(scan.f, scan.degenerate, scan.df_interaction, scan.df_residual)
    assert scan.p_value_global > 0.2 and not ubc.rejects(scan.p_value_global, ubc.RULES["scan"])
    assert reading["min_adjusted"] < 0.03 and ubc.rejects(reading["min_adjusted"], ubc.RULES["fdr"])
    assert reading["discoveries"] >= 100


# ── Merge on synthetic records ───────────────────────────────────────────────


def _write_shard(out_dir: Path, settings, records, name="shard_000_of_001.jsonl") -> None:
    directory = ubc.profile_dir(out_dir, settings["profile"])
    lsc.ensure_settings(directory, settings)
    (directory / "records").mkdir(parents=True, exist_ok=True)
    with (directory / "records" / name).open("a", encoding="utf-8") as handle:
        for record in records:
            handle.write(lsc._line(record))


def test_merge_writes_the_rate_paired_and_cell_tables(tmp_path) -> None:
    settings = _small_stage_span()
    units = ubc.profile_units(settings)
    records = []
    for unit in units:
        if (unit.cell, unit.replicate) == ("magnitude_e0.1", 1):
            records.append(_fake_record(unit, settings, status="failed"))
        elif unit.cell == "magnitude_e0.1":
            records.append(_fake_record(unit, settings, rejections={"production": True, "fdr": True}))
        elif unit.cell == "shape_e0.5":
            records.append(_fake_record(unit, settings, rejections={"production": unit.replicate == 0, "scan": True},
                                        hits=4, hits_null=1, discoveries=8, disc_null=2, rank=3))
        else:
            records.append(_fake_record(unit, settings))
    _write_shard(tmp_path, settings, records)
    _write_shard(tmp_path, settings, records[:3], name="shard_dup.jsonl")
    summary = ubc.merge(tmp_path, "stage_span")
    directory = tmp_path / "stage_span"
    assert summary["failed_units"] == 1 and summary["reproduction_passed"] is None
    assert not (directory / "reproduction.json").exists()
    accounting = json.loads((directory / "accounting.json").read_text())
    assert accounting["failures_by_cell"] == {"n80/magnitude_e0.1": 1}
    assert accounting["duplicate_records"] == 3 and accounting["complete_units"] == 26
    rates = lsc._read_csv(directory / "rate_table.csv")
    assert list(rates[0]) == list(ubc.RATE_COLUMNS)
    assert len(rates) == 13 * len(ubc.ARMS)
    by = {(r["cell"], r["arm"]): r for r in rates}
    assert (by["magnitude_e0.1", "production"]["count"], by["magnitude_e0.1", "production"]["rejections"],
            by["magnitude_e0.1", "production"]["failures"]) == ("1", "1", "1")
    assert by["magnitude_e0.1", "fdr"]["rule"] == "min BH-adjusted p <= alpha"
    paired = lsc._read_csv(directory / "paired_table.csv")
    assert list(paired[0]) == list(ubc.PAIRED_COLUMNS)
    assert len(paired) == 13 * len(ubc.BASELINE_ARMS)
    pairs = {(r["cell"], r["arm"]): r for r in paired}
    none = pairs["anchor", "scan"]
    assert (none["b"], none["c"], none["p_gain"], none["p_loss"]) == ("0", "0", "1", "1")
    gain = pairs["magnitude_e0.1", "scan"]
    assert (gain["n_paired"], gain["b"], gain["c"]) == ("1", "0", "1")
    mixed = pairs["shape_e0.5", "scan"]
    assert (mixed["n_paired"], mixed["b"], mixed["c"], mixed["rate_production"], mixed["rate_arm"]) == (
        "2", "1", "0", "0.5", "1")
    cells = {r["cell"]: r for r in lsc._read_csv(directory / "cell_table.csv")}
    assert list(next(iter(cells.values()))) == list(ubc.CELL_COLUMNS)
    shape = cells["shape_e0.5"]
    assert (shape["rank_median"], shape["scan_hits_median"], shape["fdr_discoveries_median"]) == ("3", "4", "8")
    assert (shape["scan_fdp_mean"], shape["fdr_fdp_mean"], shape["pca_matched_rank_equal"]) == ("0.25", "0.25", "1")
    assert (cells["magnitude_e0.1"]["completed"], cells["magnitude_e0.1"]["failures"]) == ("1", "1")
    assert cells["anchor"]["fdr_fdp_mean"] == "0"


def test_paired_row_with_no_discordant_unit() -> None:
    row = ubc.paired_row([(True, True), (False, False), (None, True), (True, None)])
    assert (row["n_paired"], row["b"], row["c"], row["p_gain"], row["p_loss"]) == (2, 0, 0, 1.0, 1.0)
    row = ubc.paired_row([(True, False)] * 6 + [(False, True)])
    assert (row["b"], row["c"]) == (1, 6)
    assert row["difference"] == pytest.approx(5 / 7)
    assert row["p_gain"] < row["p_loss"]


# ── Reproduction ─────────────────────────────────────────────────────────────


def _parent_rates_from_committed(settings):
    committed = lsc._read_csv(ubc.SCAN_DIR / "comparator_table.csv")
    rows = []
    for row in committed:
        if row["design"] in settings["designs"] and row["test"] in ubc.PARENT_TESTS:
            rows.append({"design": row["design"], "cell": row["cell"], "arm": ubc.PARENT_TESTS[row["test"]],
                         "count": row["count"], "rejections": row["rejections"]})
    return committed, rows


def test_parent_reproduction_both_outcomes(full_settings) -> None:
    settings = full_settings["parent"]
    committed, rates = _parent_rates_from_committed(settings)
    result = ubc.parent_reproduction(rates, settings, committed)
    assert result["passed"] and result["rows_compared"] == 2 * 14 * 7 and result["differing"] == []
    changed = [dict(row) for row in rates]
    target = next(row for row in changed if (row["design"], row["cell"], row["arm"]) == ("seaad", "anchor", "total"))
    target["rejections"] = str(int(target["rejections"]) + 1)
    result = ubc.parent_reproduction(changed, settings, committed)
    assert not result["passed"]
    assert result["differing"] == [{
        "cell": ["seaad", "anchor"], "test": "split_total",
        "committed": {"count": int(target["count"]), "rejections": int(target["rejections"]) - 1},
        "reproduced": {"count": int(target["count"]), "rejections": int(target["rejections"])},
    }]
    missing = [row for row in rates if row["arm"] != "scan"]
    result = ubc.parent_reproduction(missing, settings, committed)
    assert not result["passed"] and len(result["differing"]) == 28
    assert all(entry["reproduced"] is None for entry in result["differing"])


def _surgery_rates_from_committed(settings):
    operating = {source: lsc._read_csv(paths["operating"]) for source, paths in ubc.SURGERY_SOURCES.items()}
    rows = []
    for column in settings["designs"]:
        n_samples = str(settings["n_samples"][column])
        for cell, mode, effect, source in settings["cells"]:
            for row in operating[source]:
                if (row["evaluation.integration_params.layers"], row["generator.n_samples"], row["trajectory_mode"]) \
                        == ("", n_samples, mode) and float(row["effect_size"]) == float(effect):
                    rows.append({"design": column, "cell": cell, "arm": row["statistic"],
                                 "count": row["n_available"], "rejections": row["n_rejected"]})
    return operating, rows


def test_surgery_reproduction_both_outcomes(full_settings) -> None:
    settings = full_settings["surgery"]
    operating, rates = _surgery_rates_from_committed(settings)
    assert len(rates) == 3 * 11 * 3
    result = ubc.surgery_reproduction(rates, settings, operating)
    assert result["passed"] and result["rows_compared"] == 99
    assert result["differing"] == [] and result["unresolved_rows"] == []
    changed = [dict(row) for row in rates]
    changed[5]["count"] = str(int(changed[5]["count"]) - 1)
    result = ubc.surgery_reproduction(changed, settings, operating)
    assert not result["passed"] and len(result["differing"]) == 1
    entry = result["differing"][0]
    assert entry["cell"] == [changed[5]["design"], changed[5]["cell"]] and entry["test"] == changed[5]["arm"]
    assert entry["committed"]["count"] == entry["reproduced"]["count"] + 1
    anchor = [row for row in rates if row["cell"] == "anchor"]
    geometry_anchor = [row for row in operating["geometry"] if row["trajectory_mode"] == "none"
                       and row["evaluation.integration_params.layers"] == "" and row["generator.n_samples"] == "80"]
    assert {(r["arm"], r["count"], r["rejections"]) for r in anchor if r["design"] == "n80_three_blocks"} == {
        (r["statistic"], r["n_available"], r["n_rejected"]) for r in geometry_anchor}


def test_the_surgery_reproduction_reads_the_merged_rate_rows(tmp_path, monkeypatch) -> None:
    settings = ubc.surgery_settings(columns=["n80_three_blocks"], replicates=1, control_replicates=1)
    monkeypatch.setattr(ubc, "evaluate_unit", lambda unit, settings: _fake_record(unit, settings))
    ubc.run_shard(tmp_path, settings, 0, 1)
    summary = ubc.merge(tmp_path, "surgery")
    assert summary["reproduction_passed"] is False
    reproduction = json.loads((tmp_path / "surgery" / "reproduction.json").read_text())
    assert reproduction["rows_compared"] == 33 and len(reproduction["differing"]) == 33


# ── The decision on synthetic tables ─────────────────────────────────────────


def _rate(design, cell, arm, count, rejections):
    return {"design": design, "cell": cell, "arm": arm, "count": str(count), "rejections": str(rejections)}


def _block(settings, *, reproduction=True):
    rates, paired = [], []
    for design in settings["designs"]:
        for cell in ubc._cell_names(settings):
            control = cell in ubc.CONTROL_CELLS
            count = 500 if control else 200
            for arm in ubc.ARMS:
                rates.append(_rate(design, cell, arm, count, 25 if control else 100))
            for arm in ubc.BASELINE_ARMS:
                paired.append({"design": design, "cell": cell, "arm": arm, "n_paired": str(count),
                               "rate_production": "0.5", "rate_arm": "0.5", "b": "0", "c": "0"})
    block = {
        "settings": settings,
        "accounting": {"missing_units": 0, "failed_units": 0},
        "rates": rates,
        "paired": paired,
        "cells": [],
    }
    if reproduction:
        block["reproduction"] = {"passed": True, "rows_compared": 1, "differing": []}
    return block


@pytest.fixture()
def blocks(full_settings):
    return {
        "parent": _block(full_settings["parent"]),
        "stage_span": _block(full_settings["stage_span"], reproduction=False),
        "surgery": _block(full_settings["surgery"]),
    }


def _set_rate(block, design, cell, arm, rejections, count=None):
    for row in block["rates"]:
        if (row["design"], row["cell"], row["arm"]) == (design, cell, arm):
            row["rejections"] = str(rejections)
            if count is not None:
                row["count"] = str(count)
            return
    raise AssertionError((design, cell, arm))


def _set_pair(block, design, cell, arm, *, n=200, b=0, c=0):
    for row in block["paired"]:
        if (row["design"], row["cell"], row["arm"]) == (design, cell, arm):
            base = 100
            row.update({"n_paired": str(n), "b": str(b), "c": str(c),
                        "rate_production": str((base + c) / n) if n else "",
                        "rate_arm": str((base + b) / n) if n else ""})
            return
    raise AssertionError((design, cell, arm))


def _decide(blocks):
    return ubc.decision(blocks["parent"], blocks["stage_span"], blocks["surgery"])


def _label(result, profile, design, cell):
    return next(entry for entry in result["criteria"]["labels"]
                if (entry["profile"], entry["design"], entry["cell"]) == (profile, design, cell))


def test_every_cell_unresolved_reads_unresolved(blocks) -> None:
    result = _decide(blocks)
    assert result["verdict"] == "unresolved" and result["m_signal_cells"] == 71
    assert result["admissible_baseline_arms"] == ["scan", "fdr", "total"]
    assert len(result["criteria"]["labels"]) == 71
    assert result["levels"] == {"production_leads": pytest.approx(0.025 / 71),
                                "baseline_leads": pytest.approx(0.025 / 213)}
    assert {arm: reading["admissible"] for arm, reading in result["criteria"]["size"].items()} == {
        "production": True, "scan": True, "fdr": True, "total": True}
    families = result["criteria"]["size"]["production"]["families"]
    assert (families["stage_span"]["m"], families["surgery"]["m"]) == (8, 6)
    assert result["consequence"] == ubc.CONSEQUENCES["unresolved"]


def test_production_leads(blocks) -> None:
    for arm in ubc.BASELINE_ARMS:
        _set_pair(blocks["stage_span"], "n200", "magnitude_e0.05", arm, c=60)
    result = _decide(blocks)
    assert result["verdict"] == "production_leads"
    assert result["criteria"]["label_cells"]["production_leads"] == [["stage_span", "n200", "magnitude_e0.05"]]


def test_production_must_beat_every_admissible_arm(blocks) -> None:
    for arm in ("scan", "fdr"):
        _set_pair(blocks["stage_span"], "n200", "magnitude_e0.05", arm, c=60)
    entry = _label(_decide(blocks), "stage_span", "n200", "magnitude_e0.05")
    assert entry["label"] == "unresolved"
    assert entry["against"]["scan"]["gain"] and not entry["against"]["total"]["gain"]


def test_one_arm_is_enough_for_the_baseline(blocks) -> None:
    _set_pair(blocks["surgery"], "n80_three_blocks", "orientation_e0.5", "scan", b=60)
    result = _decide(blocks)
    assert result["verdict"] == "baseline_leads"
    entry = _label(result, "surgery", "n80_three_blocks", "orientation_e0.5")
    assert entry["label"] == "baseline_leads" and entry["leading_arms"] == ["scan"]


def test_complementary_names_the_cells_of_each_label(blocks) -> None:
    for arm in ubc.BASELINE_ARMS:
        _set_pair(blocks["stage_span"], "n500", "magnitude_e0.25", arm, c=60)
    _set_pair(blocks["stage_span"], "n80", "orient_e0.5_phi90", "fdr", b=80)
    result = _decide(blocks)
    assert result["verdict"] == "complementary"
    assert result["criteria"]["label_cells"]["production_leads"] == [["stage_span", "n500", "magnitude_e0.25"]]
    assert result["criteria"]["label_cells"]["baseline_leads"] == [["stage_span", "n80", "orient_e0.5_phi90"]]
    by_map = {(row["family"], row["design"]): row for row in result["descriptive"]["map"]}
    assert by_map["stage_span_magnitude", "n500"]["production_leads"] == 1
    assert by_map["stage_span_out_of_span", "n80"]["baseline_leads"] == 1
    assert sum(row[label] for row in by_map.values() for label in ubc.LABELS) == 71


def test_a_significant_difference_below_the_margin_is_unresolved(blocks) -> None:
    for arm in ubc.BASELINE_ARMS:
        _set_pair(blocks["stage_span"], "n200", "shape_e0.5", arm, n=500, c=45)
    entry = _label(_decide(blocks), "stage_span", "n200", "shape_e0.5")
    assert entry["label"] == "unresolved"
    assert all(reading["p_gain"] < 0.025 / 71 for reading in entry["against"].values())
    assert all(reading["difference"] == pytest.approx(0.09) for reading in entry["against"].values())
    for arm in ubc.BASELINE_ARMS:
        _set_pair(blocks["stage_span"], "n200", "shape_e0.5", arm, n=500, c=50)
    assert _label(_decide(blocks), "stage_span", "n200", "shape_e0.5")["label"] == "production_leads"


def test_an_unavailable_reading_is_unresolved_and_named(blocks) -> None:
    for arm in ubc.BASELINE_ARMS:
        _set_pair(blocks["stage_span"], "n200", "shape_e0.5", arm, c=60)
    _set_pair(blocks["stage_span"], "n200", "shape_e0.5", "total", n=0)
    entry = _label(_decide(blocks), "stage_span", "n200", "shape_e0.5")
    assert entry["label"] == "unresolved" and entry["unavailable"] == ["total"] and entry["reason"] == "unavailable"


def test_an_inadmissible_baseline_arm_leaves_the_labels(blocks) -> None:
    _set_rate(blocks["surgery"], "n200_three_blocks", "anchor", "fdr", 120)
    _set_pair(blocks["stage_span"], "n200", "orient_e0.5_phi90", "fdr", b=100)
    for arm in ("scan", "total"):
        _set_pair(blocks["stage_span"], "n200", "magnitude_e0.1", arm, c=60)
    result = _decide(blocks)
    assert result["admissible_baseline_arms"] == ["scan", "total"]
    assert result["levels"]["baseline_leads"] == pytest.approx(0.025 / (71 * 2))
    assert "fdr" in result["inadmissible_baseline_arms"]
    failing = result["inadmissible_baseline_arms"]["fdr"]["surgery"]
    assert [(row["design"], row["cell"], row["n_rejected"]) for row in failing] == [
        ("n200_three_blocks", "anchor", 120)]
    assert _label(result, "stage_span", "n200", "orient_e0.5_phi90")["label"] == "unresolved"
    assert "fdr" not in _label(result, "stage_span", "n200", "orient_e0.5_phi90")["against"]
    assert _label(result, "stage_span", "n200", "magnitude_e0.1")["label"] == "production_leads"
    assert result["verdict"] == "production_leads"


@pytest.mark.parametrize("failing_arms, a", [((), 3), (("total",), 2), (("total", "scan"), 1)])
def test_the_levels_follow_the_admissible_arms(blocks, failing_arms, a) -> None:
    for arm in failing_arms:
        _set_rate(blocks["stage_span"], "n80", "translation", arm, 150)
    result = _decide(blocks)
    assert len(result["admissible_baseline_arms"]) == a
    assert result["levels"]["baseline_leads"] == pytest.approx(0.025 / (71 * a))
    assert result["levels"]["production_leads"] == pytest.approx(0.025 / 71)


def test_production_size_excess_names_the_family_and_rates(blocks) -> None:
    _set_rate(blocks["stage_span"], "n500", "anchor", "production", 80)
    result = _decide(blocks)
    assert result["verdict"] == "production_size_excess"
    failing = result["failing"]["production_size"]
    assert list(failing) == ["stage_span"]
    assert [(row["design"], row["cell"], row["n_rejected"]) for row in failing["stage_span"]] == [
        ("n500", "anchor", 80)]


def test_no_admissible_baseline(blocks) -> None:
    for arm in ubc.BASELINE_ARMS:
        _set_rate(blocks["surgery"], "n80_three_blocks", "translation", arm, 200)
    result = _decide(blocks)
    assert result["verdict"] == "no_admissible_baseline" and result["admissible_baseline_arms"] == []
    assert result["levels"]["baseline_leads"] is None
    assert all(entry["label"] == "unresolved" for entry in result["criteria"]["labels"])


def test_the_gate_still_records_size_labels_and_grade(blocks) -> None:
    for arm in ubc.BASELINE_ARMS:
        _set_pair(blocks["stage_span"], "n200", "magnitude_e0.05", arm, c=60)
    blocks["surgery"]["reproduction"] = {"passed": False, "rows_compared": 99, "differing": [{"cell": ["x"]}]}
    result = _decide(blocks)
    assert result["verdict"] == "reproduction_failed"
    assert result["failing"]["reproduction"] == [{"profile": "surgery", "differing": [{"cell": ["x"]}]}]
    assert result["criteria"]["size"]["production"]["admissible"]
    assert result["criteria"]["label_cells"]["production_leads"] == [["stage_span", "n200", "magnitude_e0.05"]]
    assert result["characterization"]["grade"] in ubc.GRADES
    blocks["surgery"]["reproduction"]["passed"] = True
    blocks["parent"]["reproduction"] = {"passed": False, "rows_compared": 0, "differing": []}
    assert _decide(blocks)["verdict"] == "reproduction_failed"


# ── The characterization grade ───────────────────────────────────────────────


def test_the_off_target_statistics_per_mode() -> None:
    assert ubc.OFF_TARGET == {"magnitude": ("angle", "shape"), "orientation": ("delta", "shape"),
                              "shape": ("delta", "angle")}


def _detect(block, design, cell, mode, *, off_target=10):
    _set_rate(block, design, cell, "production", 180)
    for statistic in ubc.OFF_TARGET[mode]:
        _set_rate(block, design, cell, statistic, off_target)


def test_no_detected_cell(blocks) -> None:
    result = _decide(blocks)
    assert result["characterization"]["grade"] == "no_detected_cell"
    assert result["characterization"]["considered_cells"] == ubc.GRADED_CELLS_EXPECTED == 36


def test_names_kind_and_detects_only(blocks) -> None:
    _detect(blocks["stage_span"], "n200", "orient_e0.25_phi0", "orientation")
    _detect(blocks["stage_span"], "n500", "magnitude_e0.25", "magnitude")
    result = _decide(blocks)
    grade = result["characterization"]
    assert grade["grade"] == "names_kind"
    assert [(cell["design"], cell["cell"]) for cell in grade["detected_cells"]] == [
        ("n200", "orient_e0.25_phi0"), ("n500", "magnitude_e0.25")]
    assert [(row["cell"], row["statistic"]) for row in grade["family"]["observations"]] == [
        ("orient_e0.25_phi0", "delta"), ("orient_e0.25_phi0", "shape"),
        ("magnitude_e0.25", "angle"), ("magnitude_e0.25", "shape")]
    verdict = result["verdict"]
    _set_rate(blocks["stage_span"], "n500", "magnitude_e0.25", "angle", 60)
    result = _decide(blocks)
    assert result["characterization"]["grade"] == "detects_only"
    assert [(row["cell"], row["statistic"]) for row in result["characterization"]["failing"]] == [
        ("magnitude_e0.25", "angle")]
    assert result["verdict"] == verdict


def test_the_floor_is_production_detecting_in_eighty_percent(blocks) -> None:
    _detect(blocks["stage_span"], "n80", "shape_e0.5", "shape", off_target=60)
    _set_rate(blocks["stage_span"], "n80", "shape_e0.5", "production", 159)
    assert _decide(blocks)["characterization"]["grade"] == "no_detected_cell"
    _set_rate(blocks["stage_span"], "n80", "shape_e0.5", "production", 160)
    assert _decide(blocks)["characterization"]["grade"] == "detects_only"


def test_surgery_and_out_of_span_cells_are_not_graded(blocks) -> None:
    _detect(blocks["stage_span"], "n200", "shape_e0.25", "shape")
    _detect(blocks["stage_span"], "n200", "orient_e0.5_phi90", "orientation", off_target=120)
    _detect(blocks["surgery"], "n500_three_blocks", "orientation_e1", "orientation", off_target=150)
    result = _decide(blocks)
    grade = result["characterization"]
    assert grade["grade"] == "names_kind"
    assert [(cell["design"], cell["cell"]) for cell in grade["detected_cells"]] == [("n200", "shape_e0.25")]
    reported = {(row["profile"], row["design"], row["cell"]): row for row in grade["off_target_reported"]}
    assert len(reported) == 8 + 27  # the eight out-of-span cells and the 27 surgery magnitude/geometry cells
    assert reported["stage_span", "n200", "orient_e0.5_phi90"]["delta"]["rejections"] == 120
    assert reported["surgery", "n500_three_blocks", "orientation_e1"]["shape"]["rejections"] == 150


# ── Descriptive tables ───────────────────────────────────────────────────────


def test_descriptive_tables_have_one_row_per_cell(blocks) -> None:
    result = _decide(blocks)
    descriptive = result["descriptive"]
    assert len(descriptive["cell_rates"]) == 4 * 13 + 3 * 11
    for row in descriptive["cell_rates"]:
        assert set(row["rates"]) == set(ubc.ARMS)
        for arm in ("production", "scan", "fdr", "total", "concat", "pca_matched"):
            assert row["rates"][arm] is not None
    assert len(descriptive["discoveries"]) == len(descriptive["fit"]) == 4 * 13 + 3 * 11
    assert len(descriptive["parent"]) == 2 * 14
    assert set(descriptive["parent"][0]) >= {"production", "fdr", "scan", "total", "span", "complement", "delta"}
    families = {(row["family"], row["design"]) for row in descriptive["map"]}
    assert len(families) == 4 * 4 + 3 * 3


def test_recorded_arms_never_change_the_verdict_or_grade(blocks) -> None:
    for arm in ubc.BASELINE_ARMS:
        _set_pair(blocks["stage_span"], "n200", "magnitude_e0.05", arm, c=60)
    _detect(blocks["stage_span"], "n200", "orient_e0.25_phi0", "orientation")
    before = _decide(blocks)
    for block in (blocks["stage_span"], blocks["surgery"]):
        for row in block["rates"]:
            if row["arm"] in ("span", "complement", "production_minp") or row["arm"].startswith(("concat", "pca")):
                row["rejections"] = str(int(row["count"]))
    after = _decide(blocks)
    assert (after["verdict"], after["characterization"]["grade"]) == (before["verdict"],
                                                                      before["characterization"]["grade"])
    assert after["criteria"]["labels"] == before["criteria"]["labels"]


# ── The sbatch array ─────────────────────────────────────────────────────────


SBATCH = REPO / "scripts" / "univariate_baseline_array.sbatch"


def test_the_sbatch_array_parses() -> None:
    subprocess.run(["bash", "-n", str(SBATCH)], check=True)


@pytest.mark.parametrize("profile", ubc.PROFILES)
def test_the_sbatch_dry_run_prints_the_run_command(tmp_path, profile) -> None:
    env = {**os.environ, "UBC_PROFILE": profile, "UBC_OUT": str(tmp_path / "run"), "N_SHARDS": "40",
           "SLURM_ARRAY_TASK_ID": "7", "UBC_DRY_RUN": "1", "UBC_ARGS": "--limit 1"}
    result = subprocess.run(["bash", str(SBATCH)], cwd=tmp_path, env=env, capture_output=True, text=True, check=True)
    command = result.stdout.strip()
    assert command.startswith("python scripts/univariate_baseline_check.py run")
    assert f"--profile {profile}" in command and "--shard 7/40" in command and "--limit 1" in command
    assert f"--out-dir {tmp_path / 'run'}" in command
    bad = {**env, "UBC_PROFILE": "plasmode"}
    failed = subprocess.run(["bash", str(SBATCH)], cwd=tmp_path, env=bad, capture_output=True, text=True)
    assert failed.returncode == 2 and "UBC_PROFILE must be" in failed.stderr
