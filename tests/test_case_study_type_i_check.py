"""``scripts/case_study_type_i_check.py``: reproducibility, arms, resume, merge, calibration and the decision rule."""

from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import pytest


def _script():
    sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
    try:
        return __import__("case_study_type_i_check")
    finally:
        sys.path.pop(0)


check = _script()
SMALL_CV = {"cv1_splits": 2, "cv2_splits": 3, "n_repeats": 2, "max_components": 5, "random_state": 1203}
SIX = ["z", "x1", "x2", "x3", "x4", "x5"]
STATISTICS = ("delta", "angle", "shape")


def _settings(out_dir: Path) -> dict:
    args = check.build_parser().parse_args(["run", "--reps", "2", "--perms", "9", "--out-dir", str(out_dir)])
    settings = check.build_settings(args)
    settings["cv"] = dict(SMALL_CV)
    return settings


def _records(out_dir: Path, index: int = 0, count: int = 1) -> list[dict]:
    lines = check.shard_path(out_dir, index, count).read_text(encoding="utf-8").splitlines()
    records = sorted((json.loads(line) for line in lines), key=lambda record: record["key"])
    for record in records:
        record.pop("runtime_seconds")
    return records


@pytest.fixture(scope="module")
def runs(tmp_path_factory):
    outputs = []
    for name in ("first", "second"):
        out_dir = tmp_path_factory.mktemp(name)
        settings = _settings(out_dir)
        check.ensure_settings(out_dir, settings)
        assert check.run_shard(out_dir, settings, 0, 1) == 4  # 2 cells x 2 replicates
        outputs.append((out_dir, settings))
    return outputs


# ── Construction ─────────────────────────────────────────────────────────────


def test_units_are_reproducible(runs) -> None:
    (first, _), (second, _) = runs
    assert _records(first) == _records(second)


def test_controls_share_the_baseline_draw_and_the_covariates(runs) -> None:
    out_dir, settings = runs[0]
    records = _records(out_dir)
    assert {r["mode"] for r in records} == {"anchor", "translation"}
    by_replicate: dict[int, dict[str, dict]] = {}
    for record in records:
        by_replicate.setdefault(record["replicate"], {})[record["mode"]] = record
    assert set(by_replicate) == {0, 1}
    for replicate, per_mode in by_replicate.items():
        anchor, translation = per_mode["anchor"], per_mode["translation"]
        assert anchor["seeds"] == translation["seeds"]
        assert set(anchor["seeds"]) == {"generator", "covariates", "rrpp", "noise_covariates"}
        built_anchor = check.build_unit("anchor", 0.0, replicate, settings)
        built_translation = check.build_unit("translation", 1.0, replicate, settings)
        assert built_anchor["seeds"] == anchor["seeds"]
        for name in SIX:
            np.testing.assert_array_equal(built_anchor["metadata"][name], built_translation["metadata"][name])
        # The surgery is the only difference: group B's features move, the build is deterministic.
        assert not built_anchor["blocks"]["methylation"].equals(built_translation["blocks"]["methylation"])
        rebuilt = check.build_unit("translation", 1.0, replicate, settings)
        for name in settings["layers"]:
            np.testing.assert_array_equal(
                built_translation["blocks"][name].to_numpy(), rebuilt["blocks"][name].to_numpy()
            )


def test_settings_record_six_covariates_and_no_categorical_injection(runs) -> None:
    out_dir, settings = runs[0]
    assert settings["check"] == "case_study_type_i"
    assert check.declared_covariates(settings) == SIX
    assert settings["covariates"]["declared"] == SIX
    assert settings["covariates"]["noise_numeric"]["count"] == 5
    assert settings["injection"]["s_c"] == 0.0
    assert settings["injection"]["s_z"] == 3.0
    assert settings["min_components"] == 2
    assert settings["generator"]["p_dmp"] == 0.1
    assert {settings["generator"][f] for f in ("delta_methyl", "delta_expr", "delta_protein")} == {0.1}
    assert settings["generator"]["group_stage_sizes"] == [[11, 10, 28], [9, 10, 12]]
    assert settings["controls"] == {"anchor": 0.0, "translation": 1.0}
    assert settings["calibration"]["auroc_band"] == [0.6, 0.75]
    assert settings["calibration"]["rho_threshold"] == 0.8
    stored = json.loads((out_dir / "settings.json").read_text(encoding="utf-8"))
    assert stored["hash"] == check.settings_hash(settings)
    assert stored["settings"]["covariates"]["declared"] == SIX


def test_headroom_report_covers_the_translation_control(runs) -> None:
    _, settings = runs[0]
    report = check.headroom_report(settings)
    assert set(report) == {'["translation", 1.0]'}
    assert all(entry["fits"] for entry in report.values())


# ── Arms ─────────────────────────────────────────────────────────────────────


def test_three_arms_on_matched_data(runs) -> None:
    out_dir, settings = runs[0]
    for record in _records(out_dir):
        arms = record["arms"]
        assert set(arms) == {"adjusted", "unadjusted", "latent_stage"}
        for arm, entry in arms.items():
            assert entry["error"] is None
            assert entry["component_selection"] == "cv"
            assert set(entry["p_values"]) == set(STATISTICS)
            assert set(entry["statistics"]) == set(STATISTICS)
            assert set(entry["null_summary"]) == set(STATISTICS)
            assert isinstance(entry["selected_lv"], int) and entry["selected_lv"] >= 2
            assert 0.0 <= entry["cv_mean_auroc"] <= 1.0
            assert -1.0 <= entry["first_component_spearman_z"] <= 1.0
            assert entry["pooled_spectrum"] is not None
            assert set(entry) == set(check.ARM_FIELDS)
        adjusted, unadjusted, latent_stage = arms["adjusted"], arms["unadjusted"], arms["latent_stage"]
        assert adjusted["covariate_design_columns"] == SIX
        assert unadjusted["covariate_design_columns"] is None
        assert latent_stage["covariate_design_columns"] == SIX
        # The model-only arms share the covariate-blind basis.
        assert adjusted["selected_lv"] == unadjusted["selected_lv"]
        assert adjusted["cv_mean_auroc"] == unadjusted["cv_mean_auroc"]
        assert adjusted["first_component_spearman_z"] == unadjusted["first_component_spearman_z"]
        assert adjusted["latent_adjustment"] is None and unadjusted["latent_adjustment"] is None
        # The latent-adjusted arm fits its own basis on the stage conditioning plus the six covariates.
        block = latent_stage["latent_adjustment"]
        assert block["conditioning"] == "stage"
        assert block["covariates"] == SIX
        assert len(block["design_columns"]) == 2 + 6  # two stage indicators (drop-first) + six covariates
        n = sum(sum(row) for row in settings["generator"]["group_stage_sizes"])
        assert block["residual_df"] == n - 1 - 2 - 6
        assert block["cv_min_training_df"] is not None


def test_failed_arm_is_unavailable_and_the_others_proceed(runs, monkeypatch) -> None:
    import motco.analysis as analysis

    _, settings = runs[0]
    real = analysis.analyze_trajectories

    def raising(dataset, params):
        if params.covariate_latent_adjustment == "stage":
            raise analysis.TrajectoryAnalysisError("synthetic latent_stage failure")
        return real(dataset, params)

    monkeypatch.setattr(analysis, "analyze_trajectories", raising)
    built = check.build_unit("anchor", 0.0, 0, settings)
    arms = check.analyze_arms(built, settings)
    failed = arms["latent_stage"]
    assert failed["error"] == "synthetic latent_stage failure"
    assert failed["p_values"] == {} and failed["statistics"] == {} and failed["null_summary"] == {}
    assert failed["selected_lv"] is None and failed["latent_adjustment"] is None
    assert set(failed) == set(check.ARM_FIELDS)
    for arm in ("adjusted", "unadjusted"):
        assert arms[arm]["error"] is None
        assert set(arms[arm]["p_values"]) == set(STATISTICS)
    # A failed arm is unavailable in the rates: R drops for that arm alone.
    record = {"key": "k", "mode": "anchor", "effect": 0.0, "replicate": 0, "arms": arms}
    rows = check.p_value_rows({"k": record})
    rates = {(r["mode"], r["arm"], r["statistic"]): r for r in check.rejection_rates(rows, alpha=0.05, k=2.0)}
    assert rates[("anchor", "latent_stage", "delta")]["replicates"] == 0
    assert rates[("anchor", "latent_stage", "delta")]["verdict"] == "unavailable"
    assert rates[("anchor", "adjusted", "delta")]["replicates"] == 1
    assert rates[("anchor", "unadjusted", "delta")]["replicates"] == 1


# ── Resume and merge ─────────────────────────────────────────────────────────


def test_torn_line_is_re_evaluated(runs) -> None:
    out_dir, settings = runs[1]
    before = _records(out_dir)
    path = check.shard_path(out_dir, 0, 1)
    lines = path.read_text(encoding="utf-8").splitlines()
    assert len(lines) == 4
    path.write_text("\n".join(lines[:3]) + "\n" + lines[3][: len(lines[3]) // 2] + "\n", encoding="utf-8")
    assert check.run_shard(out_dir, settings, 0, 1) == 1
    assert _records(out_dir) == before


def test_resume_skips_recorded_units_under_another_shard_count(runs) -> None:
    out_dir, settings = runs[0]
    before = _records(out_dir)
    assert check.run_shard(out_dir, settings, 0, 1) == 0
    assert check.run_shard(out_dir, settings, 1, 2) == 0
    assert check.shard_path(out_dir, 1, 2).read_text(encoding="utf-8") == ""
    assert _records(out_dir) == before


def test_changed_settings_are_refused(runs) -> None:
    out_dir, settings = runs[0]
    with pytest.raises(SystemExit, match="different settings"):
        check.ensure_settings(out_dir, {**settings, "permutations": 19})


def test_merge_writes_the_tables_and_the_decision(runs) -> None:
    out_dir, settings = runs[0]
    summary = check.merge(out_dir, settings)
    assert summary["units"] == 4 and summary["missing"] == []
    assert summary["decision"] in check.DECISIONS
    assert summary["failures"] == {"adjusted": 0, "unadjusted": 0, "latent_stage": 0}
    assert len((out_dir / "merged.jsonl").read_text(encoding="utf-8").splitlines()) == 4
    pvalues = (out_dir / "case_study_type_i_pvalues.csv").read_text(encoding="utf-8").splitlines()
    assert pvalues[0] == "mode,effect,replicate,arm,statistic,p_value"
    assert len(pvalues) == 1 + 4 * 3 * 3
    replicates = (out_dir / "case_study_type_i_replicates.csv").read_text(encoding="utf-8").splitlines()
    assert replicates[0].split(",") == list(check.REPLICATE_COLUMNS)
    assert len(replicates) == 1 + 4 * 3
    rates = (out_dir / "case_study_type_i_rates.csv").read_text(encoding="utf-8").splitlines()
    assert rates[0].split(",") == list(check.RATE_COLUMNS)
    assert len(rates) == 1 + 2 * 3 * 3
    decision = json.loads((out_dir / "case_study_type_i_decision.json").read_text(encoding="utf-8"))
    assert decision["decision"] == summary["decision"]
    assert decision["judged_arms"] == ["adjusted", "latent_stage"]
    assert set(decision["controls"]) == {"anchor", "translation"}
    assert decision["controls"]["translation"]["effect"] == 1.0
    assert set(decision["verdicts"]) == {"adjusted", "latent_stage"}
    assert set(decision["verdicts"]["adjusted"]["anchor"]) == set(STATISTICS)
    assert decision["controls"]["anchor"]["unadjusted"]["delta"]["verdict"] == "descriptive"
    assert set(decision["calibration"]["criteria"]) == {"adjusted_auroc", "unadjusted_rho", "unadjusted_rejects"}
    assert decision["calibration"]["criteria"]["adjusted_auroc"]["replicates"] == 2
    assert decision["calibration"]["criteria"]["unadjusted_rho"]["replicates"] == 2
    assert decision["construction"]["strengths"] == {
        "delta": 0.1, "delta_methyl": 0.1, "delta_expr": 0.1, "delta_protein": 0.1, "p_dmp": 0.1, "s_z": 3.0
    }
    assert decision["construction"]["declared_covariates"] == SIX
    assert decision["settings_hash"] == check.settings_hash(settings)
    assert set(decision["latent"]) == set(check.ARMS)
    assert set(decision["latent"]["adjusted"]) == {'["anchor", 0.0]', '["translation", 1.0]'}
    assert decision["units"] == 4


# ── Decision rule ────────────────────────────────────────────────────────────


def _uniform(n: int, rejections: int) -> list[float]:
    return [0.01] * rejections + [0.5] * (n - rejections)


def _rows(rates: dict[tuple[str, float], dict[str, dict[str, int]]], n: int = 500) -> list[dict]:
    """``{(mode, effect): {arm: {statistic: rejections}}}`` → p-value rows."""

    rows = []
    for (mode, effect), per_arm in rates.items():
        for arm, per_statistic in per_arm.items():
            for statistic, rejections in per_statistic.items():
                for replicate, p_value in enumerate(_uniform(n, rejections)):
                    rows.append({"mode": mode, "effect": effect, "replicate": replicate, "arm": arm,
                                 "statistic": statistic, "p_value": p_value})
    return rows


HOLD = {"delta": 25, "angle": 30, "shape": 34}  # 0.068 <= 0.0695 at R = 500
EXCESS_SHAPE = {"delta": 25, "angle": 30, "shape": 36}  # 0.072 > 0.0695
EXCESS_ANGLE = {"delta": 25, "angle": 40, "shape": 30}  # 0.080 > 0.0695
INFLATED = {"delta": 200, "angle": 200, "shape": 200}
ANCHOR, TRANSLATION = ("anchor", 0.0), ("translation", 1.0)


def _both(adjusted: dict, latent_stage: dict, unadjusted: dict = INFLATED) -> dict:
    per_arm = {"adjusted": adjusted, "unadjusted": unadjusted, "latent_stage": latent_stage}
    return {ANCHOR: per_arm, TRANSLATION: per_arm}


def test_adjusted_holding_everywhere_validates_the_primary() -> None:
    decision = check.case_study_type_i_decision(_rows(_both(HOLD, HOLD)), alpha=0.05, k=2.0)
    assert decision["controls"]["anchor"]["adjusted"]["shape"]["bound"] == pytest.approx(0.0695, abs=1e-4)
    assert decision["verdicts"]["adjusted"] == {
        "anchor": dict.fromkeys(STATISTICS, "holds"), "translation": dict.fromkeys(STATISTICS, "holds")
    }
    assert decision["admissible"] == ["adjusted", "latent_stage"]
    assert decision["exceeding"] == [] and decision["unavailable"] == []
    assert decision["decision"] == "primary_validated"


def test_adjusted_shape_excess_with_latent_stage_holding_adopts_latent_stage() -> None:
    rows = _rows({
        ANCHOR: {"adjusted": EXCESS_SHAPE, "unadjusted": INFLATED, "latent_stage": HOLD},
        TRANSLATION: {"adjusted": HOLD, "unadjusted": INFLATED, "latent_stage": HOLD},
    })
    decision = check.case_study_type_i_decision(rows, alpha=0.05, k=2.0)
    assert decision["verdicts"]["adjusted"]["anchor"]["shape"] == "exceeds"
    assert decision["verdicts"]["adjusted"]["translation"] == dict.fromkeys(STATISTICS, "holds")
    assert decision["exceeding"] == [{"arm": "adjusted", "mode": "anchor", "effect": 0.0, "statistic": "shape",
                                      "rate": pytest.approx(0.072), "bound": pytest.approx(0.0695, abs=1e-4)}]
    assert decision["admissible"] == ["latent_stage"]
    assert decision["decision"] == "adopt_latent_stage"


def test_both_judged_arms_exceeding_is_diagnosed_naming_each_triple() -> None:
    rows = _rows({
        ANCHOR: {"adjusted": EXCESS_SHAPE, "unadjusted": INFLATED, "latent_stage": HOLD},
        TRANSLATION: {"adjusted": HOLD, "unadjusted": INFLATED, "latent_stage": EXCESS_ANGLE},
    })
    decision = check.case_study_type_i_decision(rows, alpha=0.05, k=2.0)
    assert decision["admissible"] == []
    assert [(e["arm"], e["mode"], e["statistic"]) for e in decision["exceeding"]] == [
        ("adjusted", "anchor", "shape"), ("latent_stage", "translation", "angle")
    ]
    assert decision["decision"] == "diagnose_configuration_excess"


def test_unavailable_p_values_lower_r_and_raise_the_bound() -> None:
    rows = _rows(_both(HOLD, HOLD))
    full = check.case_study_type_i_decision(rows, alpha=0.05, k=2.0)
    for replicate in range(10):
        rows.append({"mode": "anchor", "effect": 0.0, "replicate": replicate, "arm": "adjusted",
                     "statistic": "delta", "p_value": None})
    rows = [r for r in rows if not (r["mode"] == "anchor" and r["arm"] == "adjusted" and r["statistic"] == "delta"
                                    and r["replicate"] >= 490 and r["p_value"] is not None)]
    reduced = check.case_study_type_i_decision(rows, alpha=0.05, k=2.0)
    cell = reduced["controls"]["anchor"]["adjusted"]["delta"]
    assert cell["replicates"] == 490
    assert cell["bound"] > full["controls"]["anchor"]["adjusted"]["delta"]["bound"]
    assert cell["bound"] == pytest.approx(0.05 + 2 * (0.05 * 0.95 / 490) ** 0.5)
    assert reduced["controls"]["anchor"]["latent_stage"]["delta"]["replicates"] == 500
    absent = check.case_study_type_i_decision(_rows({ANCHOR: _both(HOLD, HOLD)[ANCHOR]}), alpha=0.05, k=2.0)
    assert set(absent["verdicts"]["adjusted"]["translation"].values()) == {"unavailable"}
    assert absent["controls"]["translation"]["effect"] is None
    assert len(absent["unavailable"]) == 6
    assert absent["decision"] == "diagnose_configuration_excess"


def test_unadjusted_excess_alone_validates_the_primary() -> None:
    inflated = check.case_study_type_i_decision(_rows(_both(HOLD, HOLD, INFLATED)), alpha=0.05, k=2.0)
    held = check.case_study_type_i_decision(_rows(_both(HOLD, HOLD, HOLD)), alpha=0.05, k=2.0)
    assert inflated["controls"]["anchor"]["unadjusted"]["delta"]["rate"] == pytest.approx(0.4)
    assert inflated["controls"]["anchor"]["unadjusted"]["delta"]["verdict"] == "descriptive"
    assert "unadjusted" not in inflated["verdicts"]
    assert inflated["decision"] == held["decision"] == "primary_validated"
    rates = check.rejection_rates(_rows(_both(HOLD, HOLD, INFLATED)), alpha=0.05, k=2.0)
    assert {r["verdict"] for r in rates if r["arm"] == "unadjusted"} == {"descriptive"}
    assert {r["verdict"] for r in rates if r["arm"] != "unadjusted"} == {"holds"}
    assert len(rates) == 2 * 3 * 3


def test_combine_decisions_takes_the_conjunction() -> None:
    validated = check.case_study_type_i_decision(_rows(_both(HOLD, HOLD)), alpha=0.05, k=2.0)
    fallback = check.case_study_type_i_decision(_rows(_both(EXCESS_SHAPE, HOLD)), alpha=0.05, k=2.0)
    diagnosed = check.case_study_type_i_decision(_rows(_both(EXCESS_SHAPE, EXCESS_ANGLE)), alpha=0.05, k=2.0)
    assert (validated["decision"], fallback["decision"], diagnosed["decision"]) == check.DECISIONS
    assert check.combine_decisions([validated, validated])["decision"] == "primary_validated"
    combined = check.combine_decisions([validated, fallback])
    assert combined["decision"] == "adopt_latent_stage"
    assert combined["admissible"] == ["latent_stage"]
    assert [run["decision"] for run in combined["runs"]] == ["primary_validated", "adopt_latent_stage"]
    assert check.combine_decisions([validated, diagnosed])["decision"] == "diagnose_configuration_excess"
    assert check.combine_decisions([fallback, diagnosed])["decision"] == "diagnose_configuration_excess"
    with pytest.raises(ValueError):
        check.combine_decisions([])


# ── Calibration ──────────────────────────────────────────────────────────────


def _calibration_settings(**overrides) -> dict:
    args = check.build_parser().parse_args(["calibrate", "--reps", "4", "--perms", "99", "--out-dir", "unused"])
    settings = check.build_settings(args)
    settings["calibration"].update(overrides)
    return settings


def _anchor_records(aurocs: list[float], rhos: list[float], unadjusted_p: list[dict]) -> list[dict]:
    records = []
    for replicate, (auroc, rho, p_values) in enumerate(zip(aurocs, rhos, unadjusted_p)):
        records.append({
            "key": check.unit_key("anchor", 0.0, replicate),
            "mode": "anchor",
            "effect": 0.0,
            "replicate": replicate,
            "arms": {
                "adjusted": {"cv_mean_auroc": auroc, "first_component_spearman_z": rho, "selected_lv": 2,
                             "p_values": dict.fromkeys(STATISTICS, 0.5), "error": None},
                "unadjusted": {"cv_mean_auroc": auroc, "first_component_spearman_z": rho, "selected_lv": 2,
                               "p_values": p_values, "error": None},
                "latent_stage": {"cv_mean_auroc": 0.5, "first_component_spearman_z": 0.1, "selected_lv": 3,
                                 "p_values": dict.fromkeys(STATISTICS, 0.5), "error": None},
            },
        })
    return records


REJECT_ALL = [dict.fromkeys(STATISTICS, 0.01)] * 4
REJECT_NONE = [dict.fromkeys(STATISTICS, 0.5)] * 4


def test_calibration_criteria_are_independent_and_met_is_their_conjunction() -> None:
    settings = _calibration_settings()
    met = check.calibration_criteria(_anchor_records([0.62, 0.70, 0.66, 0.68], [0.9, -0.85, 0.82, 0.95], REJECT_ALL),
                                     settings)
    assert met["criteria"]["adjusted_auroc"]["mean"] == pytest.approx(0.665)
    assert met["criteria"]["adjusted_auroc"]["quantiles"]["min"] == pytest.approx(0.62)
    assert met["criteria"]["unadjusted_rho"]["median_abs"] == pytest.approx(0.875)
    assert met["criteria"]["unadjusted_rejects"]["rates"] == dict.fromkeys(STATISTICS, 1.0)
    assert met["criteria"]["unadjusted_rejects"]["replicates"] == dict.fromkeys(STATISTICS, 4)
    assert all(entry["met"] for entry in met["criteria"].values()) and met["met"] is True
    # Each criterion fails on its own and takes the conjunction down with it.
    high_auroc = check.calibration_criteria(
        _anchor_records([0.9, 0.9, 0.9, 0.9], [0.9, 0.9, 0.9, 0.9], REJECT_ALL), settings
    )
    assert high_auroc["criteria"]["adjusted_auroc"]["met"] is False
    assert high_auroc["criteria"]["unadjusted_rho"]["met"] and high_auroc["criteria"]["unadjusted_rejects"]["met"]
    assert high_auroc["met"] is False
    low_rho = check.calibration_criteria(
        _anchor_records([0.65, 0.65, 0.65, 0.65], [0.5, 0.6, 0.7, 0.75], REJECT_ALL), settings
    )
    assert low_rho["criteria"]["unadjusted_rho"]["met"] is False and low_rho["met"] is False
    vacuous = check.calibration_criteria(
        _anchor_records([0.65, 0.65, 0.65, 0.65], [0.9, 0.9, 0.9, 0.9], REJECT_NONE), settings
    )
    assert vacuous["criteria"]["unadjusted_rejects"]["rates"] == dict.fromkeys(STATISTICS, 0.0)
    assert vacuous["criteria"]["unadjusted_rejects"]["met"] is False and vacuous["met"] is False
    one_statistic = check.calibration_criteria(
        _anchor_records([0.65] * 4, [0.9] * 4, [{"delta": 0.5, "angle": 0.01, "shape": 0.5}] * 4), settings
    )
    assert one_statistic["criteria"]["unadjusted_rejects"]["met"] is True
    # The band and the threshold are read from the settings.
    narrow = check.calibration_criteria(
        _anchor_records([0.62, 0.70, 0.66, 0.68], [0.9] * 4, REJECT_ALL),
        _calibration_settings(auroc_band=[0.6, 0.65], rho_threshold=0.95),
    )
    assert narrow["criteria"]["adjusted_auroc"]["met"] is False
    assert narrow["criteria"]["unadjusted_rho"]["met"] is False
    # Translation records never enter the calibration.
    with_translation = _anchor_records([0.65] * 4, [0.9] * 4, REJECT_ALL)
    with_translation.append({**with_translation[0], "mode": "translation", "effect": 1.0,
                             "arms": {arm: {**entry, "cv_mean_auroc": 0.99} for arm, entry in
                                      with_translation[0]["arms"].items()}})
    assert check.calibration_criteria(with_translation, settings)["criteria"]["adjusted_auroc"]["replicates"] == 4


def test_calibration_entries_append(tmp_path) -> None:
    settings = _calibration_settings()
    records = _anchor_records([0.62, 0.70, 0.66, 0.68], [0.9, -0.85, 0.82, 0.95], REJECT_ALL)
    entry = check.calibration_entry(records, settings)
    assert entry["strengths"] == {
        "delta": 0.1, "delta_methyl": 0.1, "delta_expr": 0.1, "delta_protein": 0.1, "p_dmp": 0.1, "s_z": 3.0
    }
    assert entry["replicates"] == 4 and entry["permutations"] == 99 and entry["met"] is True
    assert entry["declared_covariates"] == SIX
    assert entry["selected_rank"] == {"adjusted": {"2": 4}, "unadjusted": {"2": 4}, "latent_stage": {"3": 4}}
    assert [row["replicate"] for row in entry["per_replicate"]] == [0, 1, 2, 3]
    first = check.record_calibration(tmp_path, entry)
    second = check.record_calibration(tmp_path, {**entry, "met": False})
    assert len(first) == 1 and len(second) == 2
    history = json.loads((tmp_path / "calibration.json").read_text(encoding="utf-8"))
    assert [item["met"] for item in history] == [True, False]
