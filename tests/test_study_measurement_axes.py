"""Orientation estimator and integration method as design axes, end to end.

Tiny real evaluations (no mock): the estimator axis must change only ``angle``
and leave the baseline column's records equal to a configuration without the
axis; the integration-method axis must evaluate ``concat`` and ``pca`` columns;
and a forced-rank PLS column declared as an explicit point must pass the
production component-selection guard.
"""

from __future__ import annotations

from dataclasses import asdict

import pytest

from motco.simulations import SemiSyntheticTrajectoryParams, SimulationEvaluationParams
from motco.simulations.grid import run_simulation_replicate
from motco.simulations.study import DesignGrid, MatchedSeedPolicy, StudyConfig, enumerate_study
from motco.simulations.study.enumerate import DESIGN_PHASE, DESIGN_POINT_KEY
from motco.simulations.study.report import StudyReportError, assert_production_component_selection

METHOD = "evaluation.integration_method"
FORCED = "evaluation.integration_params.forced_components"
NCOMP = "evaluation.integration_params.n_components"
EST = "evaluation.orientation_estimator"

#: Fields that legitimately differ between the same dataset evaluated under two
#: cells: identity, the stamped design point (hashed into the signature), timing.
_IDENTITY_FIELDS = {"cell_id", "cell_metadata", "parameter_signature", "runtime_metadata"}


def _config(evaluation: SimulationEvaluationParams, design_grid: DesignGrid | None = None) -> StudyConfig:
    return StudyConfig(
        generator=SemiSyntheticTrajectoryParams(
            seed=2, trajectory_mode="magnitude", n_samples=60, n_stages=3, p_dmp=0.1
        ),
        evaluation=evaluation,
        trajectory_modes=("translation",),
        effect_sizes=(0.0, 1.0),
        n_replicates=2,
        base_seed=11,
        matched_seeds=MatchedSeedPolicy(enabled=True, primary_family="axes"),
        design_grid=design_grid or DesignGrid(),
    )


def _anchor(config: StudyConfig, *, point: dict | None = None):
    for cell in enumerate_study(config).cells:
        if cell.metadata.get("trajectory_mode") != "none" or cell.metadata.get("effect_size") != 0.0:
            continue
        if cell.phase == "power_primary" and point is None:
            return cell
        if cell.phase == DESIGN_PHASE and point is not None and cell.metadata[DESIGN_POINT_KEY] == point:
            return cell
    raise AssertionError(f"no anchor cell for point {point!r}")


def _comparable(record) -> dict:
    payload = asdict(record)
    runtime = {k: v for k, v in payload["runtime_metadata"].items() if k != "runtime_seconds"}
    return {**{k: v for k, v in payload.items() if k not in _IDENTITY_FIELDS}, "runtime": runtime}


def test_estimator_axis_changes_only_angle_and_leaves_the_baseline_column_unchanged() -> None:
    evaluation = SimulationEvaluationParams(
        integration_method="pca", integration_params={"n_components": 3}, permutations=19, seed=3
    )
    with_axis = _config(evaluation, DesignGrid(axes={EST: (None, "net_displacement")}))
    without_axis = _config(evaluation)

    baseline_cell = _anchor(with_axis)
    plain_cell = _anchor(without_axis)
    nd_cell = _anchor(with_axis, point={EST: "net_displacement"})
    assert nd_cell.evaluation_params.orientation_estimator == "net_displacement"
    assert baseline_cell.evaluation_params == plain_cell.evaluation_params

    for index in range(2):
        baseline = run_simulation_replicate(baseline_cell, index)
        plain = run_simulation_replicate(plain_cell, index)
        nd = run_simulation_replicate(nd_cell, index)

        assert "orientation_estimator" not in baseline.runtime_metadata
        assert _comparable(baseline) == _comparable(plain)

        assert nd.runtime_metadata["orientation_estimator"] == "net_displacement"
        assert nd.generator_seed == baseline.generator_seed
        # Same data, same draws: delta and shape are untouched, angle is re-measured.
        for statistic in ("delta", "shape"):
            assert nd.pair_statistics[statistic] == baseline.pair_statistics[statistic]
            assert nd.p_values[statistic] == baseline.p_values[statistic]
        assert nd.pair_statistics["angle"] != pytest.approx(baseline.pair_statistics["angle"])


PLS_PARAMS = {"cv1_splits": 3, "cv2_splits": 3, "n_repeats": 2, "max_components": 6, "random_state": 1}
POINTS = (
    {METHOD: "pls", FORCED: None, NCOMP: None, EST: None},
    {METHOD: "pls", FORCED: 3, NCOMP: None, EST: "net_displacement"},
    {METHOD: "pca", FORCED: None, NCOMP: 3, EST: "net_displacement"},
    {METHOD: "concat", FORCED: None, NCOMP: None, EST: "net_displacement"},
)


def test_integration_method_axis_evaluates_concat_pca_and_a_forced_pls_point() -> None:
    evaluation = SimulationEvaluationParams(
        integration_method="pls", integration_params=PLS_PARAMS, permutations=9, seed=3
    )
    config = _config(evaluation, DesignGrid(points=POINTS))
    records = {}
    for point in POINTS[1:]:
        cell = _anchor(config, point=point)
        assert cell.evaluation_params.integration_method == point[METHOD]
        record = run_simulation_replicate(cell, 0)
        assert record.status == "completed"
        assert set(record.p_values) >= {"delta", "angle", "shape"}
        records[point[METHOD] if point[FORCED] is None else "pls_forced"] = record

    forced = records["pls_forced"]
    assert forced.integration_metadata["component_selection"] == "forced"
    assert forced.integration_metadata["selected_lv"] == 3
    for method in ("concat", "pca"):
        assert records[method].integration_metadata.get("component_selection") is None
        assert records[method].integration_metadata.get("selected_lv") is None

    # The guard admits the declared forced column, and the non-PLS columns.
    assert_production_component_selection(list(records.values()))
    # ... but not a forced record whose design point declares a different rank.
    mismatched = forced.__class__(
        **{**asdict(forced), "cell_metadata": {**forced.cell_metadata, DESIGN_POINT_KEY: {**POINTS[1], FORCED: 4}}}
    )
    with pytest.raises(StudyReportError, match="undeclared forced latent rank"):
        assert_production_component_selection([mismatched])


# ── Reporting over non-PLS columns (hand-built records) ───────────────────────


def _record(point: dict, phase: str, mode: str, effect: float, index: int, *, method: str, forced: int | None):
    from motco.simulations.grid import SimulationReplicateResult

    integration: dict = {"integration_method": method}
    if method == "pls":
        integration["component_selection"] = "forced" if forced is not None else "cv"
        integration["selected_lv"] = forced if forced is not None else 2
    return SimulationReplicateResult(
        cell_id=f"{phase}-{mode}-{effect}-{point[METHOD]}-{point[FORCED]}-{point[EST]}",
        phase=phase,
        replicate_index=index,
        replicate_seed=index,
        generator_seed=index,
        evaluation_seed=3,
        parameter_signature="sig",
        status="completed",
        p_values={"delta": 0.5, "angle": 0.01 if index == 0 else 0.5, "shape": 0.5},
        cell_metadata={
            "trajectory_mode": mode,
            "effect_size": effect,
            DESIGN_POINT_KEY: dict(point),
            "varied_axis": "design_grid" if phase == DESIGN_PHASE else None,
        },
        integration_metadata=integration,
        null_summary={"angle": {"q95": 50.0}},
        config_spectrum={"version": 1, "pooled": {"relative_eigengap": 0.1}, "groups": {}},
    )


def _five_column_records() -> list:
    points = [
        ({METHOD: "pls", FORCED: None, NCOMP: None, EST: None}, "power_primary", "pls", None),
        ({METHOD: "pls", FORCED: None, NCOMP: None, EST: "net_displacement"}, DESIGN_PHASE, "pls", None),
        ({METHOD: "pls", FORCED: 4, NCOMP: None, EST: "net_displacement"}, DESIGN_PHASE, "pls", 4),
        ({METHOD: "pca", FORCED: None, NCOMP: 4, EST: "net_displacement"}, DESIGN_PHASE, "pca", None),
        ({METHOD: "concat", FORCED: None, NCOMP: None, EST: "net_displacement"}, DESIGN_PHASE, "concat", None),
    ]
    return [
        _record(point, phase, mode, effect, index, method=method, forced=forced)
        for point, phase, method, forced in points
        for mode, effect in (("none", 0.0), ("translation", 1.0))
        for index in range(4)
    ]


def test_operating_table_reports_non_pls_columns_with_unavailable_selection(tmp_path) -> None:
    import pandas as pd

    from motco.simulations.study import build_report_frames, write_report_csvs
    from motco.simulations.study.summary import summarize_combined_rule, summarize_study

    records = _five_column_records()
    frames = build_report_frames(
        summarize_study(records, alpha=0.05), summarize_combined_rule(records, alpha=0.05), records
    )
    path = write_report_csvs(frames, tmp_path)["design_point_operating"]
    table = pd.read_csv(path)
    anchors = table[(table["trajectory_mode"] == "none") & (table["statistic"] == "angle")]
    assert len(anchors) == 5
    for method in ("concat", "pca"):
        row = anchors[anchors[METHOD] == method].iloc[0]
        assert row["rejection_rate"] == pytest.approx(0.25)
        assert row["monte_carlo_se"] > 0
        assert row["median_eigengap"] == pytest.approx(0.1)
        assert row["median_angle_null_q95"] == pytest.approx(50.0)
        assert pd.isna(row["component_selection"]) and pd.isna(row["median_selected_lv"])
    assert set(anchors.loc[anchors[METHOD] == "pls", "component_selection"]) == {"cv", "forced"}


def test_rank_ladder_is_unchanged_by_columns_without_component_selection(tmp_path) -> None:
    # ``concat``/``pca`` rows have a null rank and no selected rank, so they have
    # no position on the ladder: the figure must equal the PLS-only figure.
    from motco.simulations.study import render_rank_ladder, resolve_operating_by_design_point

    frame = resolve_operating_by_design_point(_five_column_records())
    pls_only = frame[frame[METHOD] == "pls"]
    mixed = render_rank_ladder(frame, tmp_path / "mixed.png")
    reference = render_rank_ladder(pls_only, tmp_path / "pls.png")
    assert mixed is not None and reference is not None
    assert mixed.read_bytes() == reference.read_bytes()
