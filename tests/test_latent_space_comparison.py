"""Tests for ``scripts/latent_space_comparison.py``: arms on matched data, resume, tables, calibration, reading."""

from __future__ import annotations

import importlib.util
import json
import math
import sys
from pathlib import Path

import numpy as np
import pytest

_SCRIPT = Path(__file__).resolve().parents[1] / "scripts" / "latent_space_comparison.py"
_spec = importlib.util.spec_from_file_location("latent_space_comparison", _SCRIPT)
assert _spec is not None and _spec.loader is not None
lsc = importlib.util.module_from_spec(_spec)
sys.modules[_spec.name] = lsc  # dataclasses resolve their module through sys.modules
_spec.loader.exec_module(lsc)

FAST_CV = {"cv1_splits": 3, "cv2_splits": 4, "n_repeats": 1, "max_components": 6, "random_state": 1203}


def _small_settings(**overrides):
    settings = lsc.build_settings(weak_delta=0.15, replicates=2, permutations=19, designs=["seaad"])
    settings["pls_cv_params"] = dict(FAST_CV)
    settings.update(overrides)
    return settings


# ── 4.1 Arms on matched data ─────────────────────────────────────────────────


@pytest.fixture(scope="module")
def unit_records():
    settings = _small_settings()
    return settings, lsc.evaluate_unit(lsc.Unit("seaad", "orient_e0.5_phi90", 1), settings)


def test_every_arm_shares_the_unit_seeds(unit_records) -> None:
    settings, records = unit_records
    assert [r["arm"] for r in records] == list(lsc.ARMS)
    assert all(r["status"] == "completed" for r in records), [r.get("error_message") for r in records]
    assert {(r["generator_seed"], r["rrpp_seed"]) for r in records} == {
        (lsc.generator_seed(2035, "seaad", 1), lsc.rrpp_seed(2035, "seaad", "orient_e0.5_phi90", 1))
    }
    for record in records:
        assert record["span_capture"] is not None and record["population_latent"] is not None
        assert record["design"] == "seaad" and record["cell"] == "orient_e0.5_phi90" and record["replicate"] == 1
        assert record["out_of_span"]["support_size"] == record["out_of_span"]["available_features"]


def test_production_arm_equals_a_direct_evaluation(unit_records) -> None:
    from motco.simulations.evaluation import SimulationEvaluationParams, evaluate_semisynthetic_trajectory

    settings, records = unit_records
    record = next(r for r in records if r["arm"] == "pls_cv")
    dataset = lsc.unit_dataset(lsc.Unit("seaad", "orient_e0.5_phi90", 1), settings)
    direct = evaluate_semisynthetic_trajectory(
        dataset,
        SimulationEvaluationParams(
            integration_method="pls",
            integration_params={**FAST_CV, "layers": ["methylation", "expression"]},
            permutations=19,
            seed=record["rrpp_seed"],
            span_capture=True,
        ),
    )
    assert record["retained_dimensionality"] == direct.latent_matrix_metadata["selected_lv"]
    assert record["pair_statistics"] == pytest.approx(direct.pair_statistics, abs=0, rel=0)
    assert record["p_values"] == direct.p_values
    assert record["span_capture"] == direct.realized_geometry["span_capture"]


def test_matched_pca_uses_the_production_rank(unit_records) -> None:
    _, records = unit_records
    by_arm = {r["arm"]: r for r in records}
    assert by_arm["pca_matched"]["retained_dimensionality"] == by_arm["pls_cv"]["retained_dimensionality"]
    assert by_arm["pca_fixed"]["retained_dimensionality"] == 4
    assert by_arm["pls_residual_pca"]["retained_dimensionality"] == by_arm["pls_cv"]["retained_dimensionality"] + 2
    assert by_arm["concat"]["retained_dimensionality"] == 367 + 131
    assert by_arm["concat"]["span_capture"]["group_difference"] == 1.0


def test_cells_at_a_replicate_share_the_null_dataset() -> None:
    settings = _small_settings()
    anchor = lsc.unit_dataset(lsc.Unit("seaad", "anchor", 0), settings)
    rotated = lsc.unit_dataset(lsc.Unit("seaad", "orient_e0.25_phi60", 0), settings)
    other = lsc.unit_dataset(lsc.Unit("seaad", "anchor", 1), settings)
    group_a = anchor.metadata["group"] == "A"
    assert anchor.metadata.equals(rotated.metadata)
    for layer in ("methylation", "expression", "proteomics"):
        assert getattr(anchor, layer).loc[group_a].equals(getattr(rotated, layer).loc[group_a])
    assert anchor.proteomics.equals(rotated.proteomics)
    assert not anchor.methylation.equals(other.methylation)


def test_one_pls_model_is_fitted_per_unit(monkeypatch) -> None:
    import motco.simulations.evaluation as evaluation_module

    calls = {"fit": 0, "cv": 0}
    fit, cv = evaluation_module.fit_plsda_model, evaluation_module.plsda_doubleCV

    def counting_fit(*args, **kwargs):
        calls["fit"] += 1
        return fit(*args, **kwargs)

    def counting_cv(*args, **kwargs):
        calls["cv"] += 1
        return cv(*args, **kwargs)

    monkeypatch.setattr(evaluation_module, "fit_plsda_model", counting_fit)
    monkeypatch.setattr(evaluation_module, "plsda_doubleCV", counting_cv)
    records = lsc.evaluate_unit(lsc.Unit("seaad", "anchor", 0), _small_settings())
    assert all(r["status"] == "completed" for r in records)
    assert calls == {"fit": 1, "cv": 1}


def test_an_arm_failure_is_recorded_without_aborting_the_unit(monkeypatch) -> None:
    import motco.simulations.evaluation as evaluation_module

    def broken(*args, **kwargs):
        raise evaluation_module.SimulationEvaluationError("no residual space today")

    monkeypatch.setattr(evaluation_module, "pls_residual_pca_from_pls", broken)
    records = lsc.evaluate_unit(lsc.Unit("seaad", "anchor", 0), _small_settings())
    by_arm = {r["arm"]: r for r in records}
    assert by_arm["pls_residual_pca"]["status"] == "failed"
    assert by_arm["pls_residual_pca"]["error_message"] == "no residual space today"
    assert all(by_arm[arm]["status"] == "completed" for arm in lsc.ARMS if arm != "pls_residual_pca")


def test_the_design_grid() -> None:
    settings = lsc.build_settings(weak_delta=0.15)
    assert list(settings["designs"]) == ["n200", "seaad", "seaad_weak"]
    assert len(settings["cells"]) == 14
    assert len(lsc.enumerate_units(settings)) == 3 * 14 * 100
    assert settings["permutations"] == 199 and settings["base_seed"] == 2035
    weak = settings["designs"]["seaad_weak"]["generator"]
    assert weak["delta_methyl"] == weak["delta_expr"] == weak["delta_protein"] == 0.15
    assert "delta_methyl" not in settings["designs"]["seaad"]["generator"]
    assert lsc.cell_spec(settings, "orient_e0.5_phi90_k20") == {
        "cell": "orient_e0.5_phi90_k20", "mode": "orientation", "effect": 0.5, "phi": 90.0, "support": 20,
    }
    assert lsc.cell_spec(settings, "orient_e0.25_phi0")["phi"] == 0.0


# ── 4.2 Resume, settings binding and tables ──────────────────────────────────


def _fake_unit(unit, settings):
    """Deterministic stand-in for a unit, derived from its seeds only."""

    perm = lsc.rrpp_seed(settings["base_seed"], unit.design, unit.cell, unit.replicate)
    cell = lsc.cell_spec(settings, unit.cell)
    records = []
    for position, arm in enumerate(settings["arms"]):
        rng = np.random.default_rng([perm, position])
        records.append(
            {
                "design": unit.design,
                "cell": unit.cell,
                "trajectory_mode": cell["mode"],
                "effect_size": cell["effect"],
                "span_alignment": cell["phi"],
                "span_support": cell["support"],
                "replicate": unit.replicate,
                "rrpp_seed": perm,
                "settings_hash": lsc.settings_hash(settings),
                "arm": arm,
                "status": "completed",
                "retained_dimensionality": 2 + position,
                "p_values": {statistic: float(rng.uniform()) for statistic in lsc.STATISTICS},
                "span_capture": {"group_difference": float(rng.uniform()), "pooled_stage": 0.9, "random_share": 0.01},
                "population_latent": {"delta": 0.0, "angle": 40.0, "shape": 0.0},
                "intended_geometry": {"delta": 0.0, "angle": 45.0, "shape": 0.0},
                "runtime_seconds": {"pls_fit": 0.0, "arm": 1.0},
                "unit_runtime_seconds": 5.0,
            }
        )
    return records


def test_resume_skips_completed_units(tmp_path, monkeypatch) -> None:
    monkeypatch.setattr(lsc, "evaluate_unit", _fake_unit)
    settings = _small_settings()
    lsc.ensure_settings(tmp_path, settings)
    assert lsc.run_shard(tmp_path, settings, 0, 2, limit=3) == 3
    path = lsc.shard_path(tmp_path, 0, 2)
    # A torn final line, as an interrupted write leaves it.
    with path.open("a", encoding="utf-8") as handle:
        handle.write('{"design": "seaad", "cell": "anc')
    total = len(lsc.shard_units(lsc.enumerate_units(settings), 0, 2))
    assert lsc.run_shard(tmp_path, settings, 0, 2) == total - 3
    assert lsc.run_shard(tmp_path, settings, 0, 2) == 0
    lsc.run_shard(tmp_path, settings, 1, 2)
    outcome = lsc.merge(tmp_path, settings)
    assert outcome["complete_units"] == 28 and outcome["missing_units"] == 0
    assert outcome["duplicate_records"] == 0 and outcome["failed_records"] == 0
    merged = [json.loads(line) for line in (tmp_path / "merged.jsonl").read_text().splitlines()]
    keys = [(r["design"], r["cell"], r["replicate"], r["arm"]) for r in merged]
    assert len(keys) == len(set(keys)) == 28 * len(lsc.ARMS)
    accounting = json.loads((tmp_path / "accounting.json").read_text())
    assert accounting["recorded_core_hours"] == pytest.approx(28 * 5.0 / 3600)


def test_settings_mismatch_is_refused(tmp_path, monkeypatch) -> None:
    settings = _small_settings()
    lsc.ensure_settings(tmp_path, settings)
    with pytest.raises(SystemExit, match="differing: permutations, replicates"):
        lsc.ensure_settings(tmp_path, {**settings, "replicates": 3, "permutations": 99})
    with pytest.raises(SystemExit, match=r"differing: .*weak_delta\.value"):
        lsc.main(["run", "--out-dir", str(tmp_path), "--weak-delta", "0.2", "--reps", "2", "--perms", "19",
                  "--designs", "seaad"])
    # A record written under other settings is refused at the merge.
    monkeypatch.setattr(lsc, "evaluate_unit", _fake_unit)
    lsc.run_shard(tmp_path, settings, 0, 1, limit=1)
    stray = _fake_unit(lsc.Unit("seaad", "magnitude", 0), {**settings, "permutations": 99})
    with lsc.shard_path(tmp_path, 0, 1).open("a", encoding="utf-8") as handle:
        handle.write("".join(json.dumps(record) + "\n" for record in stray))
    with pytest.raises(SystemExit, match="refusing to merge"):
        lsc.merge(tmp_path, settings)


def test_table_columns_and_failure_counts() -> None:
    settings = _small_settings()
    records = [r for unit in lsc.enumerate_units(settings) for r in _fake_unit(unit, settings)]
    failed = dict(records[0], status="failed", error_type="ValueError", error_message="x")
    for key in ("p_values", "span_capture", "population_latent", "retained_dimensionality"):
        failed.pop(key)
    records[0] = failed
    rates = lsc.rate_table(records, alpha=0.05)
    assert tuple(rates[0]) == lsc.RATE_COLUMNS
    assert {(r["arm"], r["cell"]) for r in rates} == {(arm, cell[0]) for arm in lsc.ARMS for cell in lsc.CELLS}
    assert len(rates) == len(lsc.ARMS) * 14 * 3
    first = next(r for r in rates if (r["design"], r["cell"], r["arm"]) == (failed["design"], failed["cell"],
                                                                            failed["arm"]))
    assert first["count"] == 1 and first["failures"] == 1
    for row in rates:
        if row["rate"] is not None:
            assert row["se"] == pytest.approx(math.sqrt(row["rate"] * (1 - row["rate"]) / row["count"]))
    captures = lsc.capture_table(records)
    assert tuple(captures[0]) == lsc.CAPTURE_COLUMNS
    assert len(captures) == len(lsc.ARMS) * 14
    row = next(r for r in captures if r["cell"] == "orient_e0.5_phi90" and r["arm"] == "pls_cv")
    assert row["target_statistic"] == "angle"
    assert row["population_latent_target_median"] == 40.0 and row["intended_target"] == 45.0
    assert row["population_latent_angle_median"] == 40.0
    assert row["population_latent_delta_median"] == 0.0 and row["population_latent_shape_median"] == 0.0
    assert row["n"] == 2 and row["retained_rank_median"] == 2.0
    anchor = next(r for r in captures if r["cell"] == "anchor" and r["arm"] == "concat")
    assert anchor["target_statistic"] is None and anchor["population_latent_target_median"] is None
    summary = lsc.accounting(records, settings, duplicates=0, missing=[])
    assert summary["failures_by_arm_and_cell"] == {failed["arm"]: {f"{failed['design']}/{failed['cell']}": 1}}


# ── 4.3 Calibration ──────────────────────────────────────────────────────────


LADDER = (0.1, 0.15, 0.2, 0.25, 0.3, 0.4, 0.5)


def _steps(aurocs):
    return [{"delta": delta, "auroc_median": auroc} for delta, auroc in zip(LADDER, aurocs)]


def test_the_amended_ladder() -> None:
    assert lsc.CALIBRATION_LADDER == (0.1, 0.11, 0.12, 0.13, 0.14, 0.15, 0.2, 0.25, 0.3, 0.4, 0.5)
    assert list(lsc.CALIBRATION_LADDER) == sorted(lsc.CALIBRATION_LADDER)


def test_first_step_in_the_band_is_adopted() -> None:
    adopted = lsc.adopt_calibration(_steps([0.52, 0.58, 0.61, 0.64, 0.70, 0.80, 0.90]))
    assert adopted == {"delta": 0.2, "rule": "first_in_band", "auroc_median": 0.61}
    # The ladder order decides, not the input order.
    shuffled = lsc.adopt_calibration(list(reversed(_steps([0.52, 0.58, 0.61, 0.64, 0.70, 0.80, 0.90]))))
    assert shuffled["delta"] == 0.2


def test_the_closest_step_is_adopted_when_the_band_is_missed() -> None:
    adopted = lsc.adopt_calibration(_steps([0.52, 0.55, 0.59, 0.69, 0.75, 0.80, 0.90]))
    assert adopted == {"delta": 0.2, "rule": "band_missed", "auroc_median": 0.59}
    tie = lsc.adopt_calibration([{"delta": 0.3, "auroc_median": 0.686}, {"delta": 0.2, "auroc_median": 0.586}])
    assert tie["delta"] == 0.2 and tie["rule"] == "band_missed"
    with pytest.raises(SystemExit, match="No calibration step"):
        lsc.adopt_calibration([{"delta": 0.1, "auroc_median": None}])


def test_run_refuses_to_start_without_a_strength(tmp_path) -> None:
    with pytest.raises(SystemExit, match="--calibration <calibration.json>.*--weak-delta"):
        lsc.main(["run", "--out-dir", str(tmp_path / "run")])
    assert not (tmp_path / "run").exists()
    with pytest.raises(SystemExit, match="does not exist; run `calibrate` first"):
        lsc.main(["run", "--out-dir", str(tmp_path / "run"), "--calibration", str(tmp_path / "missing.json")])
    with pytest.raises(SystemExit, match="not both"):
        lsc.build_settings(weak_delta=0.1, calibration=tmp_path / "x.json")


def test_run_reads_the_adopted_delta(tmp_path, monkeypatch) -> None:
    settings = lsc.calibration_settings(replicates=2)

    curve = {0.1: 0.52, 0.11: 0.53, 0.12: 0.54, 0.13: 0.545, 0.14: 0.55, 0.15: 0.55, 0.2: 0.58, 0.25: 0.625,
             0.3: 0.66, 0.4: 0.72, 0.5: 0.8}

    def fake_unit(delta, replicate, settings):
        return {"delta": delta, "replicate": replicate, "status": "completed", "auroc": curve[delta],
                "rank": 2, "seconds": 1.0}

    monkeypatch.setattr(lsc, "calibration_unit", fake_unit)
    payload = lsc.calibrate(tmp_path, settings)
    assert [step["delta"] for step in payload["steps"]] == list(lsc.CALIBRATION_LADDER)
    assert all(step["replicates"] == 2 for step in payload["steps"])
    assert payload["adopted"] == {"delta": 0.25, "rule": "first_in_band", "auroc_median": 0.625}
    run = lsc.build_settings(calibration=tmp_path / "calibration.json")
    assert run["weak_delta"]["value"] == 0.25 and run["weak_delta"]["source"] == "calibration"
    assert run["designs"]["seaad_weak"]["generator"]["delta_expr"] == 0.25
    assert (tmp_path / "calibration_units.csv").read_text().splitlines()[0].startswith("delta,replicate")


def test_calibration_unit_fits_production_pls_on_the_seaad_anchor() -> None:
    settings = {**lsc.calibration_settings(replicates=1), "pls_cv_params": dict(FAST_CV)}
    entry = lsc.calibration_unit(0.15, 0, settings)
    assert entry["status"] == "completed", entry
    assert 0.0 <= entry["auroc"] <= 1.0 and entry["rank"] >= 2
    assert entry["generator_seed"] == lsc.calibration_seed(2035, 0)
    assert entry["generator_seed"] != lsc.generator_seed(2035, "seaad", 0)


# ── 5.1 The advisory reading ─────────────────────────────────────────────────


def _rates(overrides=None, count=100):
    """A synthetic rate table: controls at 0.05, every other rate at 0.5, then the overrides."""

    rows = []
    for design in lsc.DESIGNS:
        for cell, *_ in lsc.CELLS:
            for arm in lsc.ARMS:
                for statistic in lsc.STATISTICS:
                    rate = 0.05 if cell in ("anchor", "translation") else 0.5
                    rate = (overrides or {}).get((design, cell, arm, statistic), rate)
                    rows.append(
                        {
                            "design": design, "cell": cell, "arm": arm, "statistic": statistic, "rate": rate,
                            "se": math.sqrt(rate * (1 - rate) / count), "count": count,
                        }
                    )
    return rows


def _gain(arms, rate=0.9):
    return {(design, "orient_e0.5_phi90", arm, "angle"): rate for design in ("n200", "seaad") for arm in arms}


def test_reading_qualifies_and_proposes_a_paper_grade() -> None:
    result = lsc.reading(_rates(_gain(["pls_residual_pca"])), [], alpha=0.05)
    assert result["verdict"] == "propose_paper_grade"
    assert result["qualifying_arms"] == ["pls_residual_pca"]
    assert result["arms"]["pls_residual_pca"]["verdict"] == "qualifies"
    assert result["arms"]["pls_residual_pca"]["failing"] == {}
    assert result["arms"]["pca_matched"]["verdict"] == "no_gain"
    assert result["arms"]["pca_matched"]["failing"]["gain"] == [
        {"design": "n200", "cell": "orient_e0.5_phi90", "statistic": "angle"},
        {"design": "seaad", "cell": "orient_e0.5_phi90", "statistic": "angle"},
    ]
    assert "concat" not in result["arms"] and "concat" in result["ineligible_arms"]
    assert "pls_cv" not in result["arms"]


def test_reading_documents_a_tradeoff() -> None:
    overrides = {
        **_gain(["pca_fixed", "pca_matched"]),
        ("seaad", "shape_e0.5", "pca_fixed", "shape"): 0.2,
        ("seaad_weak", "translation", "pca_matched", "angle"): 0.12,
        # A gain at one strong design only is no gain.
        ("seaad", "orient_e0.5_phi90", "pls_residual_pca", "angle"): 0.9,
    }
    result = lsc.reading(_rates(overrides), [], alpha=0.05)
    assert result["verdict"] == "document_tradeoff"
    assert result["qualifying_arms"] == []
    fixed = result["arms"]["pca_fixed"]
    assert fixed["verdict"] == "trade_off"
    assert fixed["failing"] == {"protected": [{"design": "seaad", "cell": "shape_e0.5", "statistic": "shape"}]}
    matched = result["arms"]["pca_matched"]
    assert matched["verdict"] == "size_excess"
    assert matched["failing"]["size"] == [{"design": "seaad_weak", "cell": "translation", "statistic": "angle"}]
    residual = result["arms"]["pls_residual_pca"]
    assert residual["verdict"] == "no_gain"
    assert residual["failing"]["gain"] == [{"design": "n200", "cell": "orient_e0.5_phi90", "statistic": "angle"}]


def test_reading_thresholds_and_missing_rates() -> None:
    bound = lsc.size_bound(0.05, 100)
    assert bound == pytest.approx(0.0936, abs=5e-5)
    at_bound = {("n200", "anchor", "pca_fixed", "delta"): 0.09}
    result = lsc.reading(_rates({**at_bound, **_gain(["pca_fixed"])}), [], alpha=0.05)
    assert result["arms"]["pca_fixed"]["verdict"] == "qualifies"
    # A gain within 2 pooled SE does not count; a missing rate fails its check.
    rows = [r for r in _rates(_gain(["pca_fixed"], rate=0.6)) if not (r["cell"] == "magnitude" and
                                                                     r["arm"] == "pca_fixed")]
    result = lsc.reading(rows, [], alpha=0.05)
    gain_checks = result["arms"]["pca_fixed"]["criteria"]["gain"]["checks"]
    assert all(not check["pass"] for check in gain_checks)
    margin = 2 * math.sqrt(0.6 * 0.4 / 100 + 0.5 * 0.5 / 100)
    assert all(check["threshold"] == pytest.approx(margin) for check in gain_checks)
    assert all(check["difference"] == pytest.approx(0.1) for check in gain_checks)
    assert margin > 0.1
    protected = result["arms"]["pca_fixed"]["criteria"]["protected"]["checks"]
    assert [check.get("missing") for check in protected if check["cell"] == "magnitude"] == [True, True]


def test_reading_descriptive_tables() -> None:
    captures = [
        {"design": "n200", "cell": "orient_e0.5_phi90_k20", "arm": "pls_cv", "n": 100,
         "group_difference_mean": 0.05, "group_difference_sd": 0.02, "group_difference_median": 0.05},
        {"design": "n200", "cell": "orient_e0.5_phi90", "arm": "pls_cv", "n": 100,
         "group_difference_mean": 0.01, "group_difference_sd": 0.01, "group_difference_median": 0.01},
    ]
    overrides = {("n200", "orient_e0.5_phi90_k20", "pls_cv", "angle"): 0.3,
                 ("n200", "orient_e0.5_phi90", "pls_cv", "angle"): 0.1}
    tables = lsc.reading(_rates(overrides), captures, alpha=0.05)["descriptive"]
    contrast = next(r for r in tables["sparse_minus_dense"] if r["design"] == "n200" and r["arm"] == "pls_cv")
    assert contrast["capture_difference"] == pytest.approx(0.04)
    assert contrast["capture_pooled_se"] == pytest.approx(math.sqrt(0.02**2 / 100 + 0.01**2 / 100))
    assert contrast["angle_difference"] == pytest.approx(0.2)
    assert contrast["angle_pooled_se"] == pytest.approx(math.sqrt(0.3 * 0.7 / 100 + 0.1 * 0.9 / 100))
    curve = [r for r in tables["visibility_against_phi"] if r["design"] == "n200" and r["arm"] == "pls_cv"]
    assert len(curve) == 8
    assert next(r for r in curve if r["phi"] == 90 and r["effect"] == 0.5)["capture"] == 0.01
    controls = [r for r in tables["pls_cv_controls"] if r["design"] == "seaad"]
    assert len(controls) == 6 and all(r["holds"] for r in controls)
    assert {r["arm"] for r in tables["weak_design"]} == set(lsc.ARMS)
