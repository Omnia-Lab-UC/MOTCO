"""Predictive-subspace coordinates in the real-data entry point (``TrajectoryAnalysisParams.latent_coordinates``)."""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from motco.analysis import (
    TrajectoryAnalysisError,
    TrajectoryAnalysisParams,
    analyze_trajectories,
    write_analysis_outputs,
)
from motco.analysis import _json_safe as json_safe
from motco.simulations.evaluation import (
    AttributionDiagnosticSettings,
    SimulationEvaluationParams,
    evaluate_semisynthetic_trajectory,
    integrate_semisynthetic_dataset,
)
from tests.test_analysis import (
    COVARIATES,
    SMALL_CV,
    _as_trajectory_dataset,
    make_dataset,
    make_real_like,
    synthetic,  # noqa: F401 - the module-scoped fixture
    with_covariates,
)

PREDICTIVE = ["predictive", "predictive_orthonormal", "predictive_fitted"]


def _no_fit(*args, **kwargs):
    raise AssertionError("the setting must be rejected before any fitting")


# ── 2.1 The setting and its rejections ───────────────────────────────────────


@pytest.mark.parametrize("value", PREDICTIVE)
def test_predictive_values_are_accepted_and_recorded(value: str) -> None:
    params = TrajectoryAnalysisParams(latent_coordinates=value)
    assert params.latent_coordinates == value
    assert params.to_dict()["latent_coordinates"] == value


def test_absent_setting_is_omitted_from_to_dict() -> None:
    assert "latent_coordinates" not in TrajectoryAnalysisParams().to_dict()


@pytest.mark.parametrize("value", ["orthonormal", "scores", "whitened"])
def test_other_values_are_rejected(value: str) -> None:
    allowed = r"\['predictive', 'predictive_orthonormal', 'predictive_fitted'\] or None"
    with pytest.raises(TrajectoryAnalysisError, match=rf"Unknown latent_coordinates '{value}'; allowed: {allowed}"):
        TrajectoryAnalysisParams(latent_coordinates=value)


@pytest.mark.parametrize("value", PREDICTIVE)
def test_attribution_is_rejected(value: str) -> None:
    with pytest.raises(
        TrajectoryAnalysisError, match=rf"latent_coordinates='{value}' cannot be combined with attribution"
    ):
        TrajectoryAnalysisParams(
            latent_coordinates=value, attribution=AttributionDiagnosticSettings(enabled=True, bootstrap_replicates=0)
        )


@pytest.mark.parametrize("value", PREDICTIVE)
def test_split_test_is_rejected(value: str) -> None:
    with pytest.raises(
        TrajectoryAnalysisError,
        match=rf"latent_coordinates='{value}' cannot be combined with the split interaction test \(span_split_test\)",
    ):
        TrajectoryAnalysisParams(latent_coordinates=value, span_split_test=True)


@pytest.mark.parametrize("value", PREDICTIVE)
def test_two_stages_are_rejected_before_fitting(value: str, monkeypatch) -> None:
    import motco.analysis as analysis_module

    monkeypatch.setattr(analysis_module, "fit_pls_latent_space", _no_fit)
    monkeypatch.setattr(analysis_module, "_trajectory_inputs", _no_fit)
    metadata, blocks = make_real_like(stages=["Low", "High"])
    dataset = make_dataset(metadata, blocks, stage_order=["Low", "High"])
    with pytest.raises(
        TrajectoryAnalysisError, match=rf"latent_coordinates='{value}' needs at least three stages.*2 stage\(s\)"
    ):
        analyze_trajectories(dataset, TrajectoryAnalysisParams(latent_coordinates=value))


# ── 2.2 The record and the outputs ───────────────────────────────────────────


@pytest.fixture(scope="module")
def forced_pair():
    dataset = make_dataset()
    plain = analyze_trajectories(dataset, TrajectoryAnalysisParams(forced_components=5, permutations=19, seed=3))
    predictive = analyze_trajectories(
        dataset,
        TrajectoryAnalysisParams(forced_components=5, permutations=19, seed=3, latent_coordinates="predictive"),
    )
    return plain, predictive


def test_record_block_and_its_absence(forced_pair) -> None:
    plain, predictive = forced_pair
    assert predictive.record["measurement"]["latent_coordinates"] == {"value": "predictive", "measured_dimension": 2}
    assert predictive.record["measurement"]["selected_lv"] == 5
    assert predictive.record["latent"]["integration_params"]["latent_coordinates"] == "predictive"
    assert predictive.record["params"]["latent_coordinates"] == "predictive"
    assert "latent_coordinates" not in plain.record["measurement"]
    assert "latent_coordinates" not in plain.record["latent"]["integration_params"]
    assert "latent_coordinates" not in plain.record["params"]


def test_outputs_at_a_forced_rank_of_five(forced_pair, tmp_path) -> None:
    plain, predictive = forced_pair
    assert predictive.latent_scores.shape[1] == 2
    assert list(predictive.latent_scores.columns) == ["pred_0", "pred_1"]
    assert predictive.stage_means.shape[1] == 2 + 2
    assert list(predictive.pls_loadings.columns) == ["block", "feature", *(f"pls_{i}" for i in range(5))]
    assert predictive.pls_loadings.equals(plain.pls_loadings)
    assert predictive.record["latent"]["selected_lv"] == plain.record["latent"]["selected_lv"] == 5
    paths = write_analysis_outputs(predictive, tmp_path)
    assert pd.read_csv(paths["latent_scores"], index_col=0).shape[1] == 2
    assert pd.read_csv(paths["pls_loadings"]).shape[1] == 2 + 5
    assert list(pd.read_csv(paths["stage_means"]).columns) == ["group", "stage", "pred_0", "pred_1"]


@pytest.mark.parametrize("value", PREDICTIVE)
def test_rank_and_auroc_equal_the_analysis_without_the_setting(value: str) -> None:
    dataset = make_dataset()
    plain = analyze_trajectories(dataset, TrajectoryAnalysisParams(**SMALL_CV))
    measured = analyze_trajectories(dataset, TrajectoryAnalysisParams(**SMALL_CV, latent_coordinates=value))
    assert measured.record["latent"]["selected_lv"] == plain.record["latent"]["selected_lv"]
    assert measured.record["latent"]["cv_mean_auroc"] == plain.record["latent"]["cv_mean_auroc"]
    assert measured.pls_loadings.equals(plain.pls_loadings)
    expected_columns = 3 if value == "predictive_fitted" else 2
    assert measured.latent_scores.shape[1] == expected_columns


@pytest.mark.parametrize("value", PREDICTIVE)
def test_matches_harness_predictive(synthetic, value: str) -> None:  # noqa: F811 - the fixture
    harness_params = SimulationEvaluationParams(
        integration_method="pls", integration_params=SMALL_CV, permutations=49, seed=11, latent_coordinates=value
    )
    harness = evaluate_semisynthetic_trajectory(synthetic, harness_params)
    analysis = analyze_trajectories(
        _as_trajectory_dataset(synthetic),
        TrajectoryAnalysisParams(permutations=49, seed=11, **SMALL_CV, latent_coordinates=value),
    )
    record = analysis.record
    assert record["statistics"] == json_safe(harness.pair_statistics)
    assert record["p_values"] == json_safe(harness.p_values)
    assert record["null_summary"] == json_safe(harness.null_summary)
    assert record["config_spectrum"] == json_safe(harness.config_spectrum)
    assert record["latent"]["selected_lv"] == harness.latent_matrix_metadata["selected_lv"]
    assert record["measurement"]["latent_coordinates"]["measured_dimension"] == (
        harness.latent_matrix_metadata["measured_dimension"]
    )
    latent = integrate_semisynthetic_dataset(synthetic, harness_params)
    assert np.array_equal(latent.matrix.to_numpy(), analysis.latent_scores.to_numpy())


@pytest.mark.parametrize(
    "options",
    [
        {"residual_correction": "df"},
        {"forced_components": 4},
        {"covariate_latent_adjustment": "stage"},
    ],
    ids=["residual_correction", "forced_rank", "latent_adjustment"],
)
def test_options_run_under_the_setting(options) -> None:
    metadata, blocks = with_covariates()
    dataset = make_dataset(
        metadata, blocks, covariate_cols=list(COVARIATES), categorical_covariates=["education_level"]
    )
    result = analyze_trajectories(
        dataset,
        TrajectoryAnalysisParams(**SMALL_CV, permutations=9, seed=2, latent_coordinates="predictive", **options),
    )
    assert result.record["measurement"]["latent_coordinates"]["measured_dimension"] == 2
    assert result.latent_scores.shape[1] == 2
    assert set(result.record["p_values"]) == {"delta", "angle", "shape"}
    assert "covariates" in result.record
    if "residual_correction" in options:
        assert "residual_correction" in result.record["measurement"]
    if "forced_components" in options:
        assert result.record["latent"]["selected_lv"] == 4
