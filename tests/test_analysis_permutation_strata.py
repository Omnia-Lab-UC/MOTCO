"""Within-stage permutation draws in the real-data entry point (``TrajectoryAnalysisParams.permutation_strata``)."""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pytest

import motco.analysis as analysis_module
from motco.analysis import (
    _STAGE,
    SparsitySettings,
    TrajectoryAnalysisError,
    TrajectoryAnalysisParams,
    _trajectory_inputs,
    analyze_trajectories,
)
from motco.analysis import _json_safe as json_safe
from motco.simulations.evaluation import (
    AttributionDiagnosticSettings,
    SimulationEvaluationParams,
    _empirical_p_value,
    _extract_null_distributions,
    _summarize_null_distributions,
    evaluate_semisynthetic_trajectory,
    stage_strata,
)
from motco.stats.permutation import RRPP
from tests.test_analysis import (
    SMALL_CV,
    _as_trajectory_dataset,
    _assert_equivalent,
    make_dataset,
    make_real_like,
    synthetic,  # noqa: F401 - the module-scoped fixture
    with_covariates,
)

PERMUTATIONS = 19


def _size_excess_check():
    sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
    try:
        return __import__("predictive_size_excess_check")
    finally:
        sys.path.pop(0)


def _no_fit(*args, **kwargs):
    raise AssertionError("the setting must be rejected before any fitting")


def _params(**kwargs) -> TrajectoryAnalysisParams:
    kwargs.setdefault("permutations", PERMUTATIONS)
    kwargs.setdefault("seed", 3)
    kwargs.setdefault("forced_components", 3)
    return TrajectoryAnalysisParams(**kwargs)


def _covariate_dataset():
    metadata, blocks = with_covariates()
    return make_dataset(metadata, blocks, covariate_cols=["age", "pmi"])


def _direct_null(dataset, params: TrajectoryAnalysisParams, latent_scores, *, strata: bool):
    """The analysis's null recomputed by calling RRPP on its latent scores and covariate-adjusted models."""

    inputs = _trajectory_inputs(dataset, params)
    deltas, angles, shapes, _ = RRPP(
        latent_scores.to_numpy(dtype=float), inputs.model_full, inputs.model_reduced, inputs.ls_means,
        inputs.design.contrast, permutations=params.permutations, progress=False, seed=params.seed,
        return_eigengaps=True, orientation_estimator=params.orientation_estimator,
        residual_correction=params.residual_correction,
        strata=stage_strata(inputs.working, _STAGE, inputs.design.stage_levels) if strata else None,
    )
    return _extract_null_distributions(deltas, angles, shapes, shape_available=len(inputs.stage_levels) >= 3)


def _assert_null_is(result, null) -> None:
    record = result.record
    expected = {
        statistic: _empirical_p_value(null[statistic], observed)
        for statistic, observed in record["statistics"].items()
        if statistic in null and observed is not None and np.isfinite(observed)
    }
    assert record["p_values"] == json_safe(expected)
    assert record["null_summary"] == json_safe(_summarize_null_distributions(null))


# ── The setting and its rejections ───────────────────────────────────────────


def test_stage_is_accepted_and_recorded_in_the_params() -> None:
    params = TrajectoryAnalysisParams(permutation_strata="stage")
    assert params.to_dict()["permutation_strata"] == "stage"


def test_absent_setting_is_omitted_from_to_dict() -> None:
    assert "permutation_strata" not in TrajectoryAnalysisParams().to_dict()


@pytest.mark.parametrize("value", ["group", "cell", ""])
def test_unknown_values_are_rejected(value: str) -> None:
    with pytest.raises(
        TrajectoryAnalysisError, match=rf"Unknown permutation_strata '{value}'; allowed: \['stage'\] or None"
    ):
        TrajectoryAnalysisParams(permutation_strata=value)


def test_split_test_is_rejected() -> None:
    with pytest.raises(
        TrajectoryAnalysisError,
        match=r"permutation_strata='stage' cannot be combined with the split interaction test \(span_split_test\)",
    ):
        TrajectoryAnalysisParams(permutation_strata="stage", span_split_test=True)


def test_interaction_scan_is_rejected() -> None:
    with pytest.raises(
        TrajectoryAnalysisError,
        match=r"permutation_strata='stage' cannot be combined with the interaction scan \(interaction_scan\)",
    ):
        TrajectoryAnalysisParams(permutation_strata="stage", interaction_scan=True)


def test_rejections_come_before_any_fitting(monkeypatch) -> None:
    monkeypatch.setattr(analysis_module, "fit_pls_latent_space", _no_fit)
    monkeypatch.setattr(analysis_module, "_trajectory_inputs", _no_fit)
    with pytest.raises(TrajectoryAnalysisError, match="cannot be combined with the interaction scan"):
        analyze_trajectories(make_dataset(), TrajectoryAnalysisParams(permutation_strata="stage",
                                                                      interaction_scan=True))


# ── Only the null moves ──────────────────────────────────────────────────────


@pytest.fixture(scope="module")
def pair():
    dataset = make_dataset()
    plain = analyze_trajectories(dataset, _params())
    within = analyze_trajectories(dataset, _params(permutation_strata="stage"))
    return dataset, plain, within


def test_only_the_null_moves(pair) -> None:
    dataset, plain, within = pair
    assert within.record["statistics"] == plain.record["statistics"]
    assert within.record["latent"] == plain.record["latent"]
    assert within.latent_scores.equals(plain.latent_scores)
    assert within.pls_loadings.equals(plain.pls_loadings)
    assert within.stage_means.equals(plain.stage_means)
    assert within.record["null_summary"] != plain.record["null_summary"]
    _assert_null_is(within, _direct_null(dataset, _params(permutation_strata="stage"), within.latent_scores,
                                         strata=True))
    _assert_null_is(plain, _direct_null(dataset, _params(), plain.latent_scores, strata=False))


def test_record_block_and_its_absence(pair) -> None:
    _, plain, within = pair
    assert within.record["measurement"]["permutation_strata"] == {"value": "stage", "sizes": [16, 16, 16]}
    assert within.record["params"]["permutation_strata"] == "stage"
    assert "permutation_strata" not in plain.record["measurement"]
    assert "permutation_strata" not in plain.record["params"]


def test_cross_validated_fit_is_unchanged() -> None:
    dataset = make_dataset()
    plain = analyze_trajectories(dataset, TrajectoryAnalysisParams(**SMALL_CV, permutations=9, seed=1))
    within = analyze_trajectories(
        dataset, TrajectoryAnalysisParams(**SMALL_CV, permutations=9, seed=1, permutation_strata="stage")
    )
    assert within.record["latent"]["selected_lv"] == plain.record["latent"]["selected_lv"]
    assert within.record["latent"]["cv_mean_auroc"] == plain.record["latent"]["cv_mean_auroc"]
    assert within.latent_scores.equals(plain.latent_scores)


def test_record_sizes_of_unequal_stages() -> None:
    metadata, blocks = make_real_like(n_per_cell=20)
    replicate = metadata["donor"].str.rsplit("_", n=1).str[-1].astype(int)
    keep = (metadata["ADNC"] == "High") | (replicate < 10)
    metadata = metadata[keep].reset_index(drop=True)
    blocks = {name: frame.loc[metadata["donor"]] for name, frame in blocks.items()}
    result = analyze_trajectories(make_dataset(metadata, blocks), _params(permutation_strata="stage"))
    assert result.record["measurement"]["permutation_strata"] == {"value": "stage", "sizes": [20, 20, 40]}


def test_no_permutations_records_no_block() -> None:
    result = analyze_trajectories(make_dataset(), _params(permutations=0, permutation_strata="stage"))
    assert result.record["params"]["permutation_strata"] == "stage"
    assert "permutation_strata" not in result.record["measurement"]
    assert result.record["p_values"] == {}


# ── Declared stage order ─────────────────────────────────────────────────────


def test_four_stages_in_a_non_sorted_declared_order() -> None:
    stages = ["Mid", "Early", "Late", "Advanced"]
    metadata, blocks = make_real_like(stages=stages, n_per_cell=5)
    first = metadata["ADNC"] == "Mid"
    metadata = metadata[~first | (metadata["donor"].str.endswith(("_0", "_1", "_2")))].reset_index(drop=True)
    blocks = {name: frame.loc[metadata["donor"]] for name, frame in blocks.items()}
    dataset = make_dataset(metadata, blocks, stage_order=stages)
    params = _params(permutation_strata="stage")
    result = analyze_trajectories(dataset, params)
    assert result.record["design"]["stage_levels"] == stages
    assert result.record["measurement"]["permutation_strata"]["sizes"] == [6, 10, 10, 10]
    inputs = _trajectory_inputs(dataset, params)
    positions = stage_strata(inputs.working, _STAGE, inputs.design.stage_levels)
    declared = metadata["ADNC"].map({stage: index for index, stage in enumerate(stages)}).to_numpy()
    np.testing.assert_array_equal(positions, declared)
    _assert_null_is(result, _direct_null(dataset, params, result.latent_scores, strata=True))


def test_eleven_stages_follow_the_declared_order() -> None:
    stages = [f"S{index}" for index in range(11)]
    metadata, blocks = make_real_like(stages=stages, n_per_cell=3)
    dataset = make_dataset(metadata, blocks, stage_order=stages)
    params = _params(permutation_strata="stage", permutations=9)
    inputs = _trajectory_inputs(dataset, params)
    assert inputs.design.stage_levels[:3] == ["0", "1", "2"] and inputs.design.stage_levels[-1] == "10"
    positions = stage_strata(inputs.working, _STAGE, inputs.design.stage_levels)
    np.testing.assert_array_equal(positions, metadata["ADNC"].str[1:].astype(int).to_numpy())
    result = analyze_trajectories(dataset, params)
    assert result.record["measurement"]["permutation_strata"]["sizes"] == [6] * 11
    _assert_null_is(result, _direct_null(dataset, params, result.latent_scores, strata=True))


# ── Equivalences ─────────────────────────────────────────────────────────────


def test_matches_harness_under_within_stage_draws(synthetic) -> None:  # noqa: F811 - the fixture
    harness = evaluate_semisynthetic_trajectory(
        synthetic,
        SimulationEvaluationParams(
            integration_method="pls", integration_params=dict(SMALL_CV), permutations=49, seed=11,
            permutation_strata="stage",
        ),
    )
    analysis = analyze_trajectories(
        _as_trajectory_dataset(synthetic),
        TrajectoryAnalysisParams(permutations=49, seed=11, permutation_strata="stage", **SMALL_CV),
    )
    _assert_equivalent(harness, analysis)
    assert analysis.record["measurement"]["permutation_strata"]["sizes"] == [20, 20, 40]


def test_case_study_configuration_equals_the_diagnosis_construction() -> None:
    check = _size_excess_check()
    dataset = _covariate_dataset()
    params = _params(permutation_strata="stage", residual_correction="df", permutations=29, seed=5)
    result = analyze_trajectories(dataset, params)
    inputs = _trajectory_inputs(dataset, params)
    indices = check.within_stage_indices(inputs.working[_STAGE], inputs.design.stage_levels, 29, 5)
    inferred = check.infer(result.latent_scores, inputs, permutations=29, seed=5, indices=indices,
                           orientation_estimator="net_displacement", residual_correction="df")
    assert result.record["p_values"] == json_safe(inferred["p_values"])


# ── Combinations ─────────────────────────────────────────────────────────────


@pytest.mark.parametrize(
    "extra",
    [
        {"attribution": AttributionDiagnosticSettings(enabled=True, bootstrap_replicates=0)},
        {"block_weighting": "equal_variance"},
        {"sparsity": SparsitySettings(scope="block", keep_fractions=(0.5,))},
        {"latent_coordinates": "predictive_orthonormal"},
        {"residual_correction": "leverage"},
    ],
    ids=["attribution", "block_weighting", "sparsity", "predictive_orthonormal", "leverage"],
)
def test_combines_with_other_settings(extra) -> None:
    dataset = make_dataset()
    params = _params(permutation_strata="stage", permutations=9, **extra)
    result = analyze_trajectories(dataset, params)
    assert result.record["measurement"]["permutation_strata"]["value"] == "stage"
    _assert_null_is(result, _direct_null(dataset, params, result.latent_scores, strata=True))


def test_combines_with_latent_adjustment() -> None:
    dataset = _covariate_dataset()
    params = _params(permutation_strata="stage", permutations=9, covariate_latent_adjustment="stage")
    result = analyze_trajectories(dataset, params)
    assert result.record["measurement"]["permutation_strata"]["value"] == "stage"
    _assert_null_is(result, _direct_null(dataset, params, result.latent_scores, strata=True))
