"""Design-point reporting: operating table, continuity rekeying, decision rule, files.

Records are built by hand with known p-values, eigengaps, null widths, and
selected ranks so every expected number is constructed, not inferred.
"""

from __future__ import annotations

import json
from pathlib import Path

import matplotlib
import numpy as np
import pandas as pd
import pytest

matplotlib.use("Agg")

from motco.simulations.grid import SimulationReplicateResult
from motco.simulations.study import (
    DesignPointDecisionRule,
    build_report_frames,
    evaluate_design_point_decision,
    render_design_point_power,
    resolve_operating_by_design_point,
    resolve_orientation_by_continuity,
    resolve_stage_signal_by_design_point,
    write_design_point_decision,
    write_report_csvs,
)
from motco.simulations.study.enumerate import DESIGN_AXIS_MARKER, DESIGN_PHASE, DESIGN_POINT_KEY

RHO = "generator.baseline_continuity"
N = "generator.n_samples"


def _spectrum(gap: float) -> dict:
    return {
        "version": 1,
        "pooled": {"n_points": 4, "n_dimensions": 3, "total_variance": 4.0, "spectrum": [], "relative_eigengap": gap},
        "groups": {},
    }


def record(
    *,
    point: dict | None,
    phase: str,
    mode: str,
    effect: float,
    index: int,
    p_angle: float,
    p_delta: float = 0.5,
    gap: float = 0.05,
    null_q95: float = 90.0,
    selected_lv: int | None = 3,
    continuity: float | None = None,
    cv_auroc: float | None = None,
    integration: dict | None = None,
) -> SimulationReplicateResult:
    meta: dict = {"trajectory_mode": mode, "effect_size": effect}
    if point is not None:
        meta[DESIGN_POINT_KEY] = dict(point)
        meta["varied_axis"] = DESIGN_AXIS_MARKER if phase == DESIGN_PHASE else None
    cell = f"{phase}-{mode}-{effect}-" + "-".join(f"{k.split('.')[-1]}{v}" for k, v in (point or {}).items())
    rho = continuity if continuity is not None else (point or {}).get(RHO, 0.0)
    if integration is None:
        integration = {"integration_method": "pls", "component_selection": "cv"}
        if selected_lv is not None:
            integration["selected_lv"] = selected_lv
        if cv_auroc is not None:
            integration["cv_mean_auroc"] = cv_auroc
    return SimulationReplicateResult(
        cell_id=cell,
        phase=phase,
        replicate_index=index,
        replicate_seed=index,
        generator_seed=index,
        evaluation_seed=index,
        parameter_signature="sig",
        status="completed",
        p_values={"angle": p_angle, "delta": p_delta, "shape": 0.5},
        cell_metadata=meta,
        truth_metadata={"baseline_continuity": rho},
        integration_metadata=integration,
        null_summary={"angle": {"q95": null_q95}},
        config_spectrum=_spectrum(gap),
    )


def _column(
    point: dict,
    phase: str,
    *,
    angle_power: float,
    n: int = 20,
    gap: float = 0.05,
    q95: float = 90.0,
    auroc: float | None = None,
):
    """One design column: an anchor (never rejects) + orientation at effects 0.5 and 1.0.

    With ``auroc`` every record carries a CV stage AUROC of ``auroc + 0.01·index``.
    """

    rows = []
    n_reject = int(round(angle_power * n))
    for index in range(n):
        signal = None if auroc is None else auroc + 0.01 * index
        common = {"point": point, "phase": phase, "index": index, "gap": gap, "cv_auroc": signal}
        rows.append(record(mode="none", effect=0.0, p_angle=0.5, **common))
        rows.append(record(mode="orientation", effect=0.5, p_angle=0.5, **common))
        rows.append(
            record(
                point=point,
                phase=phase,
                mode="orientation",
                effect=1.0,
                index=index,
                p_angle=0.01 if index < n_reject else 0.5,
                gap=gap,
                null_q95=q95,
                cv_auroc=signal,
            )
        )
    return rows


def _grid_records() -> list[SimulationReplicateResult]:
    base = {RHO: 0.0, N: 300}
    return (
        _column(base, "power_primary", angle_power=0.60, gap=0.04, q95=120.0)
        + _column({RHO: 0.0, N: 600}, DESIGN_PHASE, angle_power=0.75, gap=0.05, q95=100.0)
        + _column({RHO: 0.5, N: 300}, DESIGN_PHASE, angle_power=0.85, gap=0.10, q95=60.0)
        + _column({RHO: 0.5, N: 600}, DESIGN_PHASE, angle_power=1.00, gap=0.12, q95=40.0)
    )


# ── Operating table ───────────────────────────────────────────────────────────


def test_design_point_table_resolves_every_coordinate_and_the_anchors() -> None:
    frame = resolve_operating_by_design_point(_grid_records(), alpha=0.05)
    assert list(frame.columns[:2]) == [RHO, N]
    assert len(frame) == 4 * 3 * 3  # 4 points × (anchor, e=0.5, e=1.0) × 3 statistics
    angle_top = frame[(frame["statistic"] == "angle") & (frame["effect_size"] == 1.0)].set_index([RHO, N])
    assert angle_top.loc[(0.0, 300), "rejection_rate"] == pytest.approx(0.60)
    assert angle_top.loc[(0.5, 600), "rejection_rate"] == pytest.approx(1.00)
    assert bool(angle_top.loc[(0.0, 300), "is_baseline"]) is True
    assert bool(angle_top.loc[(0.5, 300), "is_baseline"]) is False
    assert angle_top.loc[(0.5, 300), "median_eigengap"] == pytest.approx(0.10)
    assert angle_top.loc[(0.5, 300), "median_angle_null_q95"] == pytest.approx(60.0)
    assert angle_top.loc[(0.5, 300), "median_selected_lv"] == 3
    anchors = frame[(frame["trajectory_mode"] == "none") & (frame["effect_size"] == 0.0)]
    assert len(anchors) == 4 * 3
    assert (anchors["rejection_rate"] == 0.0).all()


def test_design_point_table_is_empty_without_design_cells() -> None:
    baseline_only = _column({RHO: 0.0, N: 300}, "power_primary", angle_power=0.6)
    assert resolve_operating_by_design_point(baseline_only).empty
    # Primary cells without any design_point metadata (pre-grid studies) too.
    plain = [
        record(point=None, phase="power_primary", mode="orientation", effect=1.0, index=i, p_angle=0.01)
        for i in range(5)
    ]
    assert resolve_operating_by_design_point(plain).empty


# ── Stage-signal table ────────────────────────────────────────────────────────


def _auroc_records() -> list[SimulationReplicateResult]:
    return (
        _column({RHO: 0.0, N: 300}, "power_primary", angle_power=0.60, auroc=0.60)
        + _column({RHO: 0.0, N: 600}, DESIGN_PHASE, angle_power=0.75, auroc=0.70)
        + _column({RHO: 0.5, N: 300}, DESIGN_PHASE, angle_power=0.85, auroc=0.80)
        + _column({RHO: 0.5, N: 600}, DESIGN_PHASE, angle_power=1.00, auroc=0.90)
    )


def test_stage_signal_table_has_a_row_per_coordinate_mode_and_effect_with_the_anchors() -> None:
    records = _auroc_records()
    frame = resolve_stage_signal_by_design_point(records)
    assert list(frame.columns[:2]) == [RHO, N]
    assert len(frame) == 4 * 3  # 4 points × (anchor, e=0.5, e=1.0); AUROC does not vary by statistic
    # Rows line up with the operating table's (point, mode, effect) keys, in the same order.
    keys = [RHO, N, "trajectory_mode", "effect_size"]
    operating = resolve_operating_by_design_point(records)[keys].drop_duplicates()
    assert frame[keys].to_records(index=False).tolist() == operating.to_records(index=False).tolist()

    anchors = frame[(frame["trajectory_mode"] == "none") & (frame["effect_size"] == 0.0)]
    assert len(anchors) == 4 and (anchors["n_auroc_available"] == 20).all()
    baseline = frame[frame["is_baseline"]]
    assert len(baseline) == 3 and set(baseline[N]) == {300} and set(baseline[RHO]) == {0.0}
    assert set(frame.loc[~frame["is_baseline"], "phase"]) == {DESIGN_PHASE}

    row = frame[(frame[RHO] == 0.5) & (frame[N] == 300) & (frame["effect_size"] == 1.0)].iloc[0]
    values = [0.80 + 0.01 * index for index in range(20)]
    assert row["n_replicates"] == 20 and row["n_cells"] == 1 and row["n_auroc_available"] == 20
    assert row["mean_cv_auroc"] == pytest.approx(float(np.mean(values)))
    assert row["q33_cv_auroc"] == pytest.approx(float(np.quantile(values, 1 / 3)))
    assert row["median_cv_auroc"] == pytest.approx(float(np.quantile(values, 0.5)))
    assert row["q67_cv_auroc"] == pytest.approx(float(np.quantile(values, 2 / 3)))


def test_stage_signal_reports_columns_without_a_cv_auroc_as_unavailable() -> None:
    """A forced rank records a null AUROC; concat and pca record none at all."""

    forced = {"integration_method": "pls", "component_selection": "forced", "selected_lv": 4, "cv_mean_auroc": None}
    no_cv = {
        (0.0, 600): forced,
        (0.5, 300): {"integration_method": "concat"},
        (0.5, 600): {"integration_method": "pca"},
    }
    records = _column({RHO: 0.0, N: 300}, "power_primary", angle_power=0.60, auroc=0.60) + [
        record(point={RHO: rho, N: n}, phase=DESIGN_PHASE, mode="none", effect=0.0, index=i, p_angle=0.5, integration=m)
        for (rho, n), m in no_cv.items()
        for i in range(5)
    ]
    frame = resolve_stage_signal_by_design_point(records)
    assert len(frame) == 3 + 3
    unavailable = frame[~frame["is_baseline"]]
    assert (unavailable["n_replicates"] == 5).all()
    assert (unavailable["n_auroc_available"] == 0).all()
    for column in ("mean_cv_auroc", "q33_cv_auroc", "median_cv_auroc", "q67_cv_auroc"):
        assert unavailable[column].isna().all()
    assert (frame.loc[frame["is_baseline"], "n_auroc_available"] == 20).all()


def test_stage_signal_table_is_empty_without_design_cells() -> None:
    baseline_only = _column({RHO: 0.0, N: 300}, "power_primary", angle_power=0.6, auroc=0.6)
    assert resolve_stage_signal_by_design_point(baseline_only).empty
    plain = [
        record(point=None, phase="power_primary", mode="orientation", effect=1.0, index=i, p_angle=0.01, cv_auroc=0.7)
        for i in range(5)
    ]
    assert resolve_stage_signal_by_design_point(plain).empty


# ── Continuity table rekeyed on other design coordinates ─────────────────────


def test_continuity_table_does_not_pool_across_sample_size() -> None:
    frame = resolve_orientation_by_continuity(_grid_records(), alpha=0.05)
    assert N in frame.columns
    assert list(frame.columns[:2]) == ["baseline_continuity", N]
    angle = frame[(frame["statistic"] == "angle") & (frame["effect_size"] == 1.0)]
    assert len(angle) == 4
    by = angle.set_index(["baseline_continuity", N])["rejection_rate"]
    assert by.loc[(0.0, 300)] == pytest.approx(0.60)
    assert by.loc[(0.0, 600)] == pytest.approx(0.75)
    assert by.loc[(0.5, 300)] == pytest.approx(0.85)
    assert (angle["n_cells"] == 1).all()


def test_continuity_table_without_a_design_grid_keeps_its_columns() -> None:
    records = [
        record(point=None, phase="power_primary", mode="orientation", effect=1.0, index=i, p_angle=0.01, continuity=rho)
        for rho in (0.0, 0.9)
        for i in range(5)
    ]
    frame = resolve_orientation_by_continuity(records)
    assert list(frame.columns[:3]) == ["baseline_continuity", "trajectory_mode", "effect_size"]
    assert N not in frame.columns and RHO not in frame.columns


# ── Decision rule ─────────────────────────────────────────────────────────────


def _rule(**overrides) -> DesignPointDecisionRule:
    defaults = {"trajectory_mode": "orientation", "statistic": "angle", "min_power_at_top": 0.8, "prefer": (N, RHO)}
    defaults.update(overrides)
    return DesignPointDecisionRule(**defaults)


def test_decision_picks_the_first_confirmed_column_in_preference_order() -> None:
    decision = evaluate_design_point_decision(_grid_records(), _rule(), alpha=0.05)
    # Preference is n ascending then rho ascending. At n=300: rho=0 fails (0.60),
    # rho=0.5 is 0.85 with SE 0.08 → lower bound 0.77 < 0.80 → marginal. At
    # n=600: rho=0 is 0.75 → fails; rho=0.5 is 1.00 → meets.
    assert decision.verdict == "chosen"
    assert decision.chosen == {RHO: 0.5, N: 600}
    statuses = {tuple(sorted(c.design_point.items())): c.status for c in decision.columns}
    assert statuses[tuple(sorted({RHO: 0.0, N: 300}.items()))] == "fails"
    assert statuses[tuple(sorted({RHO: 0.5, N: 300}.items()))] == "marginal"
    assert statuses[tuple(sorted({RHO: 0.0, N: 600}.items()))] == "fails"
    assert statuses[tuple(sorted({RHO: 0.5, N: 600}.items()))] == "meets"
    order = [(c.design_point[N], c.design_point[RHO]) for c in decision.columns]
    assert order == [(300, 0.0), (300, 0.5), (600, 0.0), (600, 0.5)]
    for column in decision.columns:
        assert column.anchor_rates == {"delta": 0.0, "angle": 0.0, "shape": 0.0}
        assert column.anchor_available == 20
        assert column.top_effect == 1.0


def test_decision_thresholds_come_from_the_rule() -> None:
    lenient = evaluate_design_point_decision(_grid_records(), _rule(confirmation_se_threshold=0.0))
    assert lenient.chosen == {RHO: 0.5, N: 300}  # 0.85 ≥ 0.80 without confirmation
    # 0.60 at n=20 has SE 0.11, so a 0.45 floor is confirmed at the baseline column.
    lower_floor = evaluate_design_point_decision(_grid_records(), _rule(min_power_at_top=0.45))
    assert lower_floor.chosen == {RHO: 0.0, N: 300}
    by_rho_first = evaluate_design_point_decision(_grid_records(), _rule(prefer=(RHO, N)))
    assert by_rho_first.chosen == {RHO: 0.5, N: 600}
    assert [c.design_point[RHO] for c in by_rho_first.columns] == [0.0, 0.0, 0.5, 0.5]


def test_decision_revises_the_claim_when_no_column_is_confirmed() -> None:
    # Drop the saturated column so the best remaining rate (0.85) is only marginal.
    records = [r for r in _grid_records() if r.cell_metadata.get(DESIGN_POINT_KEY) != {RHO: 0.5, N: 600}]
    decision = evaluate_design_point_decision(records, _rule())
    assert decision.verdict == "revise_claim"
    assert decision.chosen is None
    assert "revise the method or the claim" in decision.rationale
    assert {c.status for c in decision.columns} <= {"marginal", "fails"}


def test_decision_marks_columns_without_the_target_mode_unavailable() -> None:
    records = _grid_records()
    records = [
        r
        for r in records
        if not (
            r.cell_metadata.get(DESIGN_POINT_KEY) == {RHO: 0.5, N: 600}
            and r.cell_metadata["trajectory_mode"] == "orientation"
        )
    ]
    decision = evaluate_design_point_decision(records, _rule())
    column = next(c for c in decision.columns if c.design_point == {RHO: 0.5, N: 600})
    assert column.status == "unavailable"
    assert column.anchor_available == 20
    assert decision.verdict == "revise_claim"


def test_decision_without_design_records_has_nothing_to_evaluate() -> None:
    decision = evaluate_design_point_decision(_column({RHO: 0.0, N: 300}, "power_primary", angle_power=1.0), _rule())
    assert decision.verdict == "no_design_grid"
    assert decision.columns == ()


def test_decision_files_are_written(tmp_path: Path) -> None:
    decision = evaluate_design_point_decision(_grid_records(), _rule())
    paths = write_design_point_decision(decision, tmp_path)
    payload = json.loads(paths["design_point_decision"].read_text(encoding="utf-8"))
    assert payload["verdict"] == "chosen"
    assert payload["chosen"] == {RHO: 0.5, N: 600}
    assert payload["rule"]["min_power_at_top"] == 0.8
    assert len(payload["columns"]) == 4
    frame = pd.read_csv(paths["design_point_decision_csv"])
    assert set(frame["status"]) == {"fails", "marginal", "meets"}
    assert frame["chosen"].sum() == 1
    assert "anchor_angle_rate" in frame.columns


# ── Report wiring ─────────────────────────────────────────────────────────────


def test_report_writes_design_outputs_only_for_design_grid_studies(tmp_path: Path) -> None:
    grid = build_report_frames([], [], _grid_records())
    assert not grid.design_point_operating.empty
    assert not grid.design_point_stage_signal.empty
    paths = write_report_csvs(grid, tmp_path / "grid")
    assert paths["design_point_operating"].exists()
    # The stage-signal table sits beside the operating table (no AUROC recorded here: zero available).
    assert paths["design_point_stage_signal"] == tmp_path / "grid" / "design_point_stage_signal.csv"
    signal = pd.read_csv(paths["design_point_stage_signal"])
    assert len(signal) == 4 * 3 and (signal["n_auroc_available"] == 0).all()
    figure = render_design_point_power(grid.design_point_operating, tmp_path / "grid" / "design_point_power.png")
    assert figure.exists() and figure.stat().st_size > 0

    plain = build_report_frames([], [], _column({RHO: 0.0, N: 300}, "power_primary", angle_power=0.6))
    assert plain.design_point_operating.empty
    assert plain.design_point_stage_signal.empty
    plain_paths = write_report_csvs(plain, tmp_path / "plain")
    assert "design_point_operating" not in plain_paths
    assert "design_point_stage_signal" not in plain_paths
    assert not (tmp_path / "plain" / "design_point_operating.csv").exists()
    assert not (tmp_path / "plain" / "design_point_stage_signal.csv").exists()
    assert {
        "specificity_matrix",
        "power_curves",
        "type_i_table",
        "config_spectrum",
        "eigengap_stratified_power",
        "realized_surgery",
    } <= set(plain_paths)
    placeholder = render_design_point_power(plain.design_point_operating, tmp_path / "plain" / "empty.png")
    assert placeholder.exists()


def test_design_cells_are_invisible_to_baseline_readers() -> None:
    from motco.simulations.study.summary import summarize_combined_rule, summarize_study

    records = _grid_records()
    per_stat = summarize_study(records, alpha=0.05)
    combined = summarize_combined_rule(records, alpha=0.05)
    frames = build_report_frames(per_stat, combined, records)
    assert set(frames.power_curves["phase"]) == {"power_primary"}
    assert set(frames.specificity_matrix["phase"]) == {"power_primary"}
    assert frames.type_i_table.empty


def test_list_valued_design_axes_survive_a_jsonl_round_trip() -> None:
    """Size tables and layer lists come back from JSONL as nested lists; they must still key a point."""

    sizes = "generator.group_stage_sizes"
    layers = "evaluation.integration_params.layers"
    baseline = {sizes: [[11, 10, 28], [9, 10, 12]], layers: ["methylation", "expression"]}
    # Same size table, layers differing only by list vs None: the three-block column.
    three_blocks = {sizes: [[11, 10, 28], [9, 10, 12]], layers: None}
    smaller = {sizes: [[10, 9, 25], [9, 9, 12]], layers: ["methylation", "expression"]}
    records = [
        record(point=baseline, phase="power_primary", mode="orientation", effect=1.0, index=i, p_angle=0.01)
        for i in range(4)
    ] + [
        record(point=point, phase=DESIGN_PHASE, mode="orientation", effect=1.0, index=i, p_angle=0.9)
        for point in (three_blocks, smaller)
        for i in range(4)
    ]
    frame = resolve_operating_by_design_point(records)
    angle = frame[frame["statistic"] == "angle"]
    assert len(angle) == 3
    primary = angle[angle["phase"] == "power_primary"].iloc[0]
    assert primary["rejection_rate"] == 1.0
    assert primary[sizes] == ((11, 10, 28), (9, 10, 12))
    assert primary[layers] == ("methylation", "expression")
    design = angle[angle["phase"] == DESIGN_PHASE]
    assert design["rejection_rate"].tolist() == [0.0, 0.0]
    assert design[layers].tolist().count(None) == 1
