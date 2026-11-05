"""``scripts/stage_span_paper_grade_decision.py`` on synthetic paper-grade reports."""

from __future__ import annotations

import json
import sys
from pathlib import Path

import pandas as pd
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
try:
    import stage_span_paper_grade_decision as decision
finally:
    sys.path.pop(0)

N, SIZES, STAGES, LAYERS = (
    "generator.n_samples",
    "generator.group_stage_sizes",
    "generator.n_stages",
    "evaluation.integration_params.layers",
)
TWO = "('methylation', 'expression')"
SEAAD = "((11, 10, 28), (9, 10, 12))"
AXES = {
    "n200": {N: 200, SIZES: None, STAGES: 3, LAYERS: TWO},
    "n80": {N: 80, SIZES: None, STAGES: 3, LAYERS: TWO},
    "seaad": {N: 600, SIZES: SEAAD, STAGES: 3, LAYERS: TWO},
    "n200_three_blocks": {N: 200, SIZES: None, STAGES: 3, LAYERS: None},
    "n200_four_stages": {N: 200, SIZES: None, STAGES: 4, LAYERS: TWO},
}
PROFILE_CELLS = {
    "span_geometry": (
        decision.COLUMNS,
        [("orientation", e) for e in (0.1, 0.25, 0.5)]
        + [("shape", e) for e in (0.1, 0.25, 0.5)]
        + [("translation", e) for e in (0.1, 0.25, 0.5)],
    ),
    "span_magnitude": (decision.COLUMNS, [("magnitude", e) for e in (0.02, 0.05, 0.1)]),
    "surgery_geometry": (
        decision.SURGERY_COLUMNS,
        [(mode, e) for mode in ("orientation", "shape", "translation") for e in (0.25, 0.5, 1.0)],
    ),
}
R = 500


def _default_rejected(profile, column, mode, effect, statistic) -> int:
    if mode in decision.TARGETS and statistic == decision.TARGETS[mode]:
        return int(R * min(0.99, 0.3 + 2 * effect))
    return 25


def _write_report(root: Path, profile: str, *, rejected=None, capture=None, available=None, drop=()) -> None:
    rejected = rejected or (lambda *key: _default_rejected(profile, *key))
    capture = capture or (lambda column, mode, effect: 0.95)
    available = available or (lambda column, mode, effect: R)
    report = root / profile / "report"
    report.mkdir(parents=True)
    columns, cells = PROFILE_CELLS[profile]
    operating, shares = [], []
    for column in columns:
        phase = "power_primary" if column == ("seaad" if profile == "surgery_geometry" else "n200") else "power_design"
        for mode, effect in [("none", 0.0), *cells]:
            n_available = available(column, mode, effect)
            for statistic in decision.STATISTICS:
                if (column, mode, effect, statistic) in drop:
                    continue
                count = min(rejected(column, mode, effect, statistic), n_available)
                rate = count / n_available
                operating.append(
                    {
                        **AXES[column],
                        "phase": phase,
                        "trajectory_mode": mode,
                        "effect_size": effect,
                        "statistic": statistic,
                        "n_replicates": R,
                        "n_available": n_available,
                        "n_rejected": count,
                        "rejection_rate": rate,
                        "monte_carlo_se": (rate * (1 - rate) / n_available) ** 0.5,
                        "median_selected_lv": 3,
                    }
                )
            constant = mode in ("none", "translation")
            shares.append(
                {
                    **AXES[column],
                    "phase": phase,
                    "cell_id": f"{profile}-{column}-{mode}-{effect}",
                    "trajectory_mode": mode,
                    "effect_size": effect,
                    "construction": "stage_span" if profile.startswith("span") and mode != "none" else "surgery",
                    "n_replicates": R,
                    "n_group_difference_available": 0 if constant else R,
                    "q10_group_difference": None if constant else capture(column, mode, effect) - 0.02,
                    "median_group_difference": None if constant else capture(column, mode, effect),
                    "q90_group_difference": None if constant else capture(column, mode, effect) + 0.01,
                    "median_pooled_stage": 0.97,
                    "median_rank": 3,
                    "median_population_latent_target": 0.4,
                    "median_intended_target": 0.4,
                }
            )
    shares.append({**shares[0], N: None, SIZES: None, "phase": "type_i_baseline", "cell_id": "type-i"})
    pd.DataFrame(operating).to_csv(report / "design_point_operating.csv", index=False)
    pd.DataFrame(shares).to_csv(report / "span_capture.csv", index=False)


def _results(tmp_path: Path, **overrides) -> Path:
    root = tmp_path / "results"
    for profile in decision.PROFILES:
        _write_report(root, profile, **overrides.get(profile, {}))
    return root


def _run(tmp_path: Path, **overrides) -> dict:
    root = _results(tmp_path, **overrides)
    return decision.run(root, slice_report=tmp_path / "absent", pilot_capture=tmp_path / "absent.csv")


# --------------------------------------------------------------------------- #
# 5.1 The bound rule
# --------------------------------------------------------------------------- #


def test_bound_rule_matches_the_design() -> None:
    p500, p100 = decision.exceedance_probability(500), decision.exceedance_probability(100)
    assert p500 == pytest.approx(0.0303, abs=5e-4)
    assert p100 == pytest.approx(0.0282, abs=5e-4)
    assert decision.allowed_exceedances(60, p500) == 4
    assert decision.allowed_exceedances(90, p500) == 6
    assert decision.allowed_exceedances(60, p100) == 4
    assert decision.bonferroni_ceiling(500, 60) == pytest.approx(0.084)
    assert decision.bonferroni_ceiling(500, 90) == pytest.approx(0.084)
    assert decision.bonferroni_ceiling(100, 60) == pytest.approx(0.130)
    # The bound at R = 500 is 0.0695: 34 rejections hold it, 35 exceed it.
    assert decision.exceedance_threshold(500) == 34


def test_judge_rates_counts_exceedances_against_the_allowed_number() -> None:
    def rows(counts):
        return [{"n_rejected": count, "n_available": 500, "statistic": "delta"} for count in counts]

    four = decision.judge_rates(rows([36] * 4 + [25] * 56))
    assert four["holds"] and four["observed_exceedances"] == 4 and four["allowed_exceedances"] == 4
    five = decision.judge_rates(rows([36] * 5 + [25] * 55))
    assert not five["holds"] and five["observed_exceedances"] == 5 and five["n_above_ceiling"] == 0
    ceiling = decision.judge_rates(rows([43] + [25] * 59))
    assert not ceiling["holds"] and ceiling["n_above_ceiling"] == 1 and ceiling["observed_exceedances"] == 1
    at_ceiling = decision.judge_rates(rows([42] + [25] * 59))
    assert at_ceiling["holds"]
    unavailable = decision.judge_rates([*rows([25] * 59), {"statistic": "delta", "status": "unavailable"}])
    assert not unavailable["holds"] and unavailable["n_unavailable"] == 1 and unavailable["m"] == 60


def test_each_row_uses_its_own_available_replicates() -> None:
    rows = [{"n_rejected": 35, "n_available": 491}] + [{"n_rejected": 25, "n_available": 500}] * 59
    judged = decision.judge_rates(rows)
    first = judged["observations"][0]
    assert first["bound"] == pytest.approx(0.05 + 2 * (0.05 * 0.95 / 491) ** 0.5)
    assert first["status"] == "exceeds"
    assert judged["median_n_available"] == 500


@pytest.mark.parametrize(
    ("row", "expected"),
    [
        ({N: 200, SIZES: None, STAGES: 3, LAYERS: TWO}, "n200"),
        ({N: "80.0", SIZES: float("nan"), STAGES: "3", LAYERS: TWO}, "n80"),
        ({N: 600, SIZES: SEAAD, STAGES: 3, LAYERS: TWO}, "seaad"),
        ({N: 200, SIZES: None, STAGES: 3, LAYERS: float("nan")}, "n200_three_blocks"),
        ({N: 200, SIZES: None, STAGES: 4, LAYERS: TWO}, "n200_four_stages"),
        ({N: 200, SIZES: ""}, "n200"),
        ({N: 80, STAGES: 3, LAYERS: ""}, "n80_three_blocks"),
    ],
)
def test_column_names_are_read_from_any_axis_set(row, expected) -> None:
    assert decision.column_name(row) == expected


# --------------------------------------------------------------------------- #
# 5.2 Criteria and verdicts
# --------------------------------------------------------------------------- #


def test_clean_reports_establish_sensitivity(tmp_path: Path) -> None:
    result = _run(tmp_path)
    assert result["verdict"] == "sensitivity_established"
    assert result["failing_criteria"] == []
    assert result["failures"] == []
    assert result["criteria"]["P1"]["m"] == 60
    assert result["criteria"]["P2"]["m"] == 45
    assert result["criteria"]["P3"]["m"] == 90
    assert result["criteria"]["P1"]["allowed_exceedances"] == 4
    assert result["criteria"]["P3"]["allowed_exceedances"] == 6
    assert result["criteria"]["P1"]["ceiling_at_median"] == pytest.approx(0.084)
    written = json.loads((tmp_path / "results" / decision.OUT_NAME).read_text())
    assert written["verdict"] == "sensitivity_established"
    assert written["rule_version"] == decision.DECISION_RULE_VERSION


def _geometry(rejected) -> dict:
    return {"span_geometry": {"rejected": rejected}}


def test_an_anchor_above_the_ceiling_is_a_size_excess(tmp_path: Path) -> None:
    def rejected(column, mode, effect, statistic):
        if (column, mode, statistic) == ("n80", "none", "delta"):
            return 60
        return _default_rejected("span_geometry", column, mode, effect, statistic)

    result = _run(tmp_path, **_geometry(rejected))
    assert result["verdict"] == "size_excess"
    named = [f for f in result["failures"] if f["criterion"] == "P1"]
    assert [(f["column"], f["trajectory_mode"], f["statistic"], f["status"]) for f in named] == [
        ("n80", "none", "delta", "above_ceiling")
    ]


@pytest.mark.parametrize(("n_exceeding", "verdict"), [(4, "sensitivity_established"), (5, "size_excess")])
def test_translation_exceedances_are_counted(tmp_path, n_exceeding, verdict) -> None:
    exceeding = [(column, 0.1, "angle") for column in decision.COLUMNS][:n_exceeding]

    def rejected(column, mode, effect, statistic):
        if mode == "translation" and (column, effect, statistic) in exceeding:
            return 36
        return _default_rejected("span_geometry", column, mode, effect, statistic)

    result = _run(tmp_path, **_geometry(rejected))
    assert result["verdict"] == verdict
    assert result["criteria"]["P1"]["observed_exceedances"] == n_exceeding
    assert sum(f["status"] == "exceeds" for f in result["failures"]) == n_exceeding


def test_a_capture_shortfall_is_named(tmp_path: Path) -> None:
    def capture(column, mode, effect):
        return 0.85 if (column, mode, effect) == ("seaad", "shape", 0.5) else 0.95

    result = _run(tmp_path, span_geometry={"capture": capture})
    assert result["verdict"] == "capture_shortfall"
    assert [(f["criterion"], f["column"], f["trajectory_mode"], f["effect_size"]) for f in result["failures"]] == [
        ("P2", "seaad", "shape", 0.5)
    ]


def test_an_off_target_rate_above_the_ceiling_is_a_specificity_excess(tmp_path: Path) -> None:
    def rejected(column, mode, effect, statistic):
        if (column, mode, effect, statistic) == ("n200_four_stages", "shape", 0.5, "delta"):
            return 50
        return _default_rejected("span_geometry", column, mode, effect, statistic)

    result = _run(tmp_path, **_geometry(rejected))
    assert result["verdict"] == "specificity_excess"
    assert result["failing_criteria"] == ["P3"]
    [failure] = result["failures"]
    assert (failure["column"], failure["trajectory_mode"], failure["statistic"]) == (
        "n200_four_stages",
        "shape",
        "delta",
    )


def test_magnitude_off_target_rates_enter_p3(tmp_path: Path) -> None:
    def rejected(column, mode, effect, statistic):
        if (column, mode, effect, statistic) == ("n80", "magnitude", 0.1, "shape"):
            return 60
        return _default_rejected("span_magnitude", column, mode, effect, statistic)

    result = _run(tmp_path, span_magnitude={"rejected": rejected})
    assert result["verdict"] == "specificity_excess"
    assert result["failures"][0]["profile"] == "span_magnitude"


def test_the_first_failing_criterion_sets_the_verdict_and_all_failures_are_listed(tmp_path: Path) -> None:
    def rejected(column, mode, effect, statistic):
        if (column, mode, statistic) == ("n200", "none", "angle"):
            return 60
        if (column, mode, effect, statistic) == ("n80", "orientation", 0.25, "shape"):
            return 60
        return _default_rejected("span_geometry", column, mode, effect, statistic)

    result = _run(
        tmp_path,
        span_geometry={"rejected": rejected, "capture": lambda column, mode, effect: 0.5},
    )
    assert result["verdict"] == "size_excess"
    assert result["failing_criteria"] == ["P1", "P2", "P3"]
    assert {f["criterion"] for f in result["failures"]} == {"P1", "P2", "P3"}


def test_a_missing_row_is_unavailable_and_fails_its_criterion(tmp_path: Path) -> None:
    result = _run(tmp_path, span_geometry={"drop": {("seaad", "translation", 0.1, "angle")}})
    assert result["verdict"] == "size_excess"
    [failure] = result["failures"]
    assert failure["status"] == "unavailable"
    assert (failure["column"], failure["trajectory_mode"], failure["statistic"]) == ("seaad", "translation", "angle")
    assert result["criteria"]["P1"]["m"] == 60


def test_a_missing_report_is_an_error(tmp_path: Path) -> None:
    root = tmp_path / "results"
    _write_report(root, "span_geometry")
    with pytest.raises(decision.DecisionError, match="Missing"):
        decision.run(root)


# --------------------------------------------------------------------------- #
# 5.3 Descriptive outputs
# --------------------------------------------------------------------------- #


def test_target_power_informative_effect_and_incomplete_cells(tmp_path: Path) -> None:
    def available(column, mode, effect):
        return 491 if (column, mode, effect) == ("n200_four_stages", "orientation", 0.5) else R

    result = _run(tmp_path, span_geometry={"available": available})
    descriptive = result["descriptive"]
    power = {(r["column"], r["trajectory_mode"], r["effect_size"]): r for r in descriptive["target_power_and_capture"]}
    assert len(power) == 45
    shape = power[("n80", "shape", 0.25)]
    assert shape["target_statistic"] == "shape" and shape["rate"] == pytest.approx(0.8)
    assert shape["q10_group_difference"] == pytest.approx(0.93) and shape["median_rank"] == 3
    informative = {(r["column"], r["trajectory_mode"]): r["informative_from"] for r in descriptive["informative_from"]}
    assert informative[("n80", "shape")] == 0.25
    assert informative[("n80", "magnitude")] is None
    assert descriptive["incomplete_cells"] == [
        {
            "profile": "span_geometry",
            "column": "n200_four_stages",
            "trajectory_mode": "orientation",
            "effect_size": 0.5,
            "n_replicates": 500,
            "n_available": 491,
            "n_missing": 9,
        }
    ]
    assert result["verdict"] == "sensitivity_established"


def test_missing_contrast_sources_are_unavailable_with_their_labels(tmp_path: Path) -> None:
    result = _run(tmp_path)
    contrast = {
        (r["column"], r["trajectory_mode"], r["effect_size"]): r for r in result["descriptive"]["visibility_contrast"]
    }
    here = contrast[("seaad", "orientation", 0.5)]
    assert here["power"]["status"] == "available" and here["power"]["source"] == decision.RUN_SOURCE["source"]
    assert here["capture"]["status"] == "available"
    elsewhere = contrast[("n200", "shape", 1.0)]
    assert elsewhere["power"] == {
        **decision.SLICE_SOURCE,
        "path": str(tmp_path / "absent" / "design_point_operating.csv"),
        "status": "unavailable",
    }
    assert elsewhere["capture"]["status"] == "unavailable"
    assert elsewhere["capture"]["source"] == decision.PILOT_SOURCE["source"]
    assert result["inputs"][str(tmp_path / "absent.csv")] == "unavailable"


def test_present_contrast_sources_are_read_and_three_block_capture_is_unavailable(tmp_path: Path) -> None:
    root = _results(tmp_path)
    slice_report = tmp_path / "slice"
    slice_report.mkdir()
    rows = []
    for column, axes in (("n200", (200, TWO)), ("n80", (80, TWO)), ("n200_three_blocks", (200, None))):
        for mode in ("orientation", "shape"):
            for effect in (0.25, 0.5, 1.0):
                for statistic in decision.STATISTICS:
                    rows.append(
                        {
                            LAYERS: axes[1],
                            N: axes[0],
                            STAGES: 3,
                            "phase": "power_design",
                            "trajectory_mode": mode,
                            "effect_size": effect,
                            "statistic": statistic,
                            "n_replicates": 500,
                            "n_available": 500,
                            "n_rejected": 100,
                            "rejection_rate": 0.2,
                            "monte_carlo_se": 0.018,
                        }
                    )
    pd.DataFrame(rows).to_csv(slice_report / "design_point_operating.csv", index=False)
    pilot = tmp_path / "pilot_capture.csv"
    pd.DataFrame(
        [
            {
                SIZES: None,
                N: n,
                "phase": "power_design",
                "cell_id": f"pilot-{n}-{mode}-{effect}",
                "trajectory_mode": mode,
                "effect_size": effect,
                "construction": "surgery",
                "n_replicates": 100,
                "n_group_difference_available": 100,
                "q10_group_difference": 0.001,
                "median_group_difference": 0.01,
                "q90_group_difference": 0.02,
                "median_pooled_stage": 0.95,
                "median_rank": 3,
                "median_population_latent_target": None,
                "median_intended_target": None,
            }
            for n in (200, 80)
            for mode in ("orientation", "shape")
            for effect in (0.25, 0.5, 1.0)
        ]
    ).to_csv(pilot, index=False)
    result = decision.run(root, slice_report=slice_report, pilot_capture=pilot)
    contrast = {
        (r["column"], r["trajectory_mode"], r["effect_size"]): r for r in result["descriptive"]["visibility_contrast"]
    }
    assert contrast[("n200", "orientation", 0.5)]["power"]["rate"] == pytest.approx(0.2)
    assert contrast[("n200", "orientation", 0.5)]["power"]["grade"] == decision.SLICE_SOURCE["grade"]
    assert contrast[("n80", "shape", 1.0)]["capture"]["median_group_difference"] == pytest.approx(0.01)
    assert contrast[("n80", "shape", 1.0)]["capture"]["grade"] == decision.PILOT_SOURCE["grade"]
    assert contrast[("n200_three_blocks", "shape", 1.0)]["power"]["status"] == "available"
    assert contrast[("n200_three_blocks", "shape", 1.0)]["capture"]["status"] == "unavailable"
