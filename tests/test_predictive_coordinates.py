"""Predictive-subspace coordinates of the PLS measurement (``latent_coordinates``)."""

from __future__ import annotations

import json
from dataclasses import replace

import numpy as np
import pandas as pd
import pytest

import motco.simulations.evaluation as evaluation_module
from motco.simulations.diagnostics import (
    POPULATION_LATENT_CHECKPOINT,
    PREDICTIVE_COORDINATES,
    LinearLatentMap,
    calculate_realized_geometry,
    geometry_from_means,
    predictive_directions,
)
from motco.simulations.evaluation import (
    LATENT_COORDINATES,
    AttributionDiagnosticSettings,
    SimulationEvaluationError,
    SimulationEvaluationParams,
    build_simulation_trajectory_design,
    evaluate_semisynthetic_trajectory,
    integrate_semisynthetic_dataset,
    pls_coordinates_from_pls,
    pls_orthonormal_from_pls,
)
from motco.simulations.grid import _to_jsonable, make_simulation_cell, parameter_signature
from motco.simulations.preprocessing import concatenate_blocks, fit_omics_preprocessor
from motco.simulations.reference import load_reference
from motco.simulations.semisynthetic import SemiSyntheticTrajectoryParams, generate_semisynthetic_trajectory
from motco.stats.trajectory import estimate_difference
from tests.orthonormal_coordinates_fixture_capture import DATASET, FAST_CV, evaluation_summary
from tests.test_latent_adjustment import _fit_core, _with_nuisance
from tests.test_orthonormal_coordinates import FIXTURE, _assert_matches

PREDICTIVE = list(PREDICTIVE_COORDINATES)


@pytest.fixture(scope="module")
def dataset():
    return generate_semisynthetic_trajectory(DATASET, reference=load_reference())


@pytest.fixture(scope="module")
def two_stage_dataset():
    params = SemiSyntheticTrajectoryParams(
        seed=3, n_samples=60, n_stages=2, p_dmp=0.1, trajectory_mode="magnitude", group_effect_size=0.2
    )
    return generate_semisynthetic_trajectory(params, reference=load_reference())


def _params(**overrides) -> SimulationEvaluationParams:
    base = {"integration_method": "pls", "integration_params": FAST_CV, "permutations": 0}
    return SimulationEvaluationParams(**{**base, **overrides})


def _forced(rank: int, **overrides) -> SimulationEvaluationParams:
    return _params(integration_params={**FAST_CV, "forced_components": rank}, **overrides)


@pytest.fixture(scope="module", params=[2, 5], ids=["rank2", "rank5"])
def forced(request, dataset):
    """The score integration at a forced rank and the three predictive integrations built from it."""

    rank = request.param
    scores = integrate_semisynthetic_dataset(dataset, _forced(rank))
    built = {
        value: pls_coordinates_from_pls(dataset, scores, _forced(rank, latent_coordinates=value))
        for value in PREDICTIVE
    }
    return rank, scores, built


@pytest.fixture(scope="module")
def cross_validated(dataset):
    scores = integrate_semisynthetic_dataset(dataset, _params())
    built = {value: integrate_semisynthetic_dataset(dataset, _params(latent_coordinates=value)) for value in PREDICTIVE}
    return scores, built


def _inputs(latent) -> np.ndarray:
    model = latent.artifacts.model
    return (latent.artifacts.features.to_numpy(dtype=float) - model._x_mean) / model._x_std


def _strip(metadata: dict) -> dict:
    """The metadata without the entries a coordinate system adds or describes."""

    described = ("measured_dimension", "shape", "n_features")
    stripped = {key: value for key, value in metadata.items() if key not in described}
    stripped["integration_params"] = {
        key: value for key, value in metadata["integration_params"].items() if key != "latent_coordinates"
    }
    return stripped


def _no_preprocessing(*args, **kwargs):
    raise AssertionError("the setting must be rejected before any preprocessing")


# --------------------------------------------------------------------------- #
# 1.1 The setting
# --------------------------------------------------------------------------- #


@pytest.mark.parametrize("value", PREDICTIVE)
def test_predictive_values_are_accepted(value: str) -> None:
    params = SimulationEvaluationParams(latent_coordinates=value)
    assert params.latent_coordinates == value
    assert params.resolved_latent_coordinates == value


def test_allowed_values_and_scores_normalization() -> None:
    assert LATENT_COORDINATES == ("scores", "orthonormal", "predictive", "predictive_orthonormal", "predictive_fitted")
    assert SimulationEvaluationParams(latent_coordinates="scores") == SimulationEvaluationParams()


def test_unknown_value_names_the_allowed_values() -> None:
    allowed = r"\['scores', 'orthonormal', 'predictive', 'predictive_orthonormal', 'predictive_fitted'\]"
    with pytest.raises(SimulationEvaluationError, match=rf"Unknown latent_coordinates 'predicted'; allowed: {allowed}"):
        SimulationEvaluationParams(latent_coordinates="predicted")


@pytest.mark.parametrize("value", PREDICTIVE)
def test_cell_identity_moves_only_when_set(value: str) -> None:
    plain = make_simulation_cell(phase="power_primary", generator_params=DATASET)
    predictive = make_simulation_cell(
        phase="power_primary",
        generator_params=DATASET,
        evaluation_params=SimulationEvaluationParams(latent_coordinates=value),
    )
    assert predictive.cell_id != plain.cell_id
    assert parameter_signature(predictive) != parameter_signature(plain)
    assert _to_jsonable(predictive.evaluation_params)["latent_coordinates"] == value
    assert "latent_coordinates" not in _to_jsonable(SimulationEvaluationParams())


# --------------------------------------------------------------------------- #
# 1.2 The predictive step in the core
# --------------------------------------------------------------------------- #


def test_column_counts_and_ranks(forced) -> None:
    rank, scores, built = forced
    assert scores.matrix.shape[1] == rank
    for value, columns in (("predictive", 2), ("predictive_orthonormal", 2), ("predictive_fitted", 3)):
        matrix = built[value].matrix.to_numpy()
        assert matrix.shape == (scores.matrix.shape[0], columns), value
        assert np.linalg.matrix_rank(matrix, tol=1e-8 * np.abs(matrix).max()) == 2, value
        assert built[value].metadata["measured_dimension"] == 2
        assert built[value].metadata["shape"] == matrix.shape
        assert built[value].metadata["n_features"] == columns
        assert built[value].matrix.index.equals(scores.matrix.index)
    assert list(built["predictive"].matrix.columns) == ["pred_0", "pred_1"]
    assert list(built["predictive_fitted"].matrix.columns) == ["fitted_0", "fitted_1", "fitted_2"]


def test_the_fitted_stage_labels_are_kept(forced) -> None:
    _, scores, built = forced
    model = scores.artifacts.model
    T = scores.matrix.to_numpy()
    Q = np.asarray(model.y_loadings_)
    V = predictive_directions(model)
    np.testing.assert_allclose(T @ V @ V.T @ Q.T, T @ Q.T, atol=1e-8)
    # The predictive coordinates carry the fitted labels: T·V·(Q·V)ᵀ = T·Qᵀ.
    np.testing.assert_allclose(built["predictive"].matrix.to_numpy() @ (Q @ V).T, T @ Q.T, atol=1e-8)
    predictions = model.predict(scores.artifacts.features.to_numpy(dtype=float))
    np.testing.assert_allclose(
        built["predictive_fitted"].matrix.to_numpy(), (predictions - model._y_mean) / model._y_std, atol=1e-8
    )


def test_the_three_matrices_are_linear_images_of_one_another(forced) -> None:
    _, _, built = forced
    reference = built["predictive"].matrix.to_numpy()
    for value in ("predictive_orthonormal", "predictive_fitted"):
        other = built[value].matrix.to_numpy()
        mapping, *_ = np.linalg.lstsq(reference, other, rcond=None)
        np.testing.assert_allclose(reference @ mapping, other, atol=1e-8)
        assert np.linalg.matrix_rank(mapping, tol=1e-8) == 2, value
        back, *_ = np.linalg.lstsq(other, reference, rcond=None)
        np.testing.assert_allclose(other @ back, reference, atol=1e-8)


def test_predictive_orthonormal_is_the_inputs_on_the_subspace_basis(forced) -> None:
    _, scores, built = forced
    model = scores.artifacts.model
    basis, _ = np.linalg.qr(np.asarray(model.x_rotations_) @ predictive_directions(model))
    latent = built["predictive_orthonormal"]
    np.testing.assert_allclose(latent.matrix.to_numpy(), _inputs(latent) @ basis, atol=1e-12)
    np.testing.assert_allclose(latent.linear_map.basis.T @ latent.linear_map.basis, np.eye(2), atol=1e-12)


def test_rank_and_cv_are_those_of_the_score_measurement(cross_validated) -> None:
    scores, built = cross_validated
    for value, latent in built.items():
        assert latent.metadata["integration_params"]["latent_coordinates"] == value
        assert _strip(latent.metadata) == _strip(scores.metadata)
        assert latent.metadata["selected_lv"] == scores.metadata["selected_lv"]
        assert latent.metadata["cv_mean_auroc"] == scores.metadata["cv_mean_auroc"]
        assert latent.artifacts.model.x_rotations_.tobytes() == scores.artifacts.model.x_rotations_.tobytes()
    assert "latent_coordinates" not in scores.metadata["integration_params"]
    assert "measured_dimension" not in scores.metadata


def test_the_score_map_is_unchanged_and_the_maps_measure_the_matrices(cross_validated) -> None:
    scores, built = cross_validated
    rows = scores.artifacts.features.to_numpy(dtype=float)
    assert scores.matrix.to_numpy().tobytes() == np.asarray(scores.artifacts.model.x_scores_).tobytes()
    for value, latent in built.items():
        np.testing.assert_allclose(latent.linear_map.transform(rows), latent.matrix.to_numpy(), atol=1e-10)
        for attribute in ("center", "scale"):
            assert getattr(latent.linear_map, attribute).tobytes() == getattr(scores.linear_map, attribute).tobytes()
        assert latent.linear_map.basis.tobytes() == built["predictive"].linear_map.basis.tobytes(), value


def test_latent_adjusted_path_measures_the_raw_features() -> None:
    inputs = _with_nuisance()
    scores = _fit_core(inputs, {"forced_components": 4}, nuisance=inputs["nuisance"])
    model = scores.artifacts.model
    assert model.scale is False
    V = predictive_directions(model)
    X = inputs["X"].to_numpy()
    for value in PREDICTIVE:
        latent = _fit_core(inputs, {"forced_components": 4}, nuisance=inputs["nuisance"], coordinates=value)
        assert latent.metadata["integration_params"]["latent_adjustment"] == "stage"
        if value == "predictive":
            np.testing.assert_allclose(latent.matrix.to_numpy(), model.transform(X) @ V, atol=1e-10)
        elif value == "predictive_orthonormal":
            basis, _ = np.linalg.qr(np.asarray(model.x_rotations_) @ V)
            np.testing.assert_allclose(latent.matrix.to_numpy(), (X - model._x_mean) @ basis, atol=1e-10)
        else:
            np.testing.assert_allclose(latent.matrix.to_numpy(), model.transform(X) @ model.y_loadings_.T, atol=1e-10)


def test_a_degenerate_fit_names_the_rank(dataset, monkeypatch) -> None:
    monkeypatch.setattr(
        evaluation_module, "predictive_directions", lambda model: np.zeros((model.x_rotations_.shape[1], 1))
    )
    with pytest.raises(SimulationEvaluationError, match=r"have rank 1, below min\(rank, n_stages − 1\) = 2"):
        integrate_semisynthetic_dataset(dataset, _forced(3, latent_coordinates="predictive"))


def test_the_core_rejects_two_stages() -> None:
    inputs = _with_nuisance()
    keep = inputs["y"] != inputs["y"].iloc[-1]
    two = {
        "X": inputs["X"][keep.to_numpy()].reset_index(drop=True),
        "y": inputs["y"][keep].reset_index(drop=True),
        "sample_ids": [s for s, k in zip(inputs["sample_ids"], keep) if k],
    }
    with pytest.raises(SimulationEvaluationError, match="'predictive_fitted' needs at least three stages"):
        _fit_core(two, coordinates="predictive_fitted")


# --------------------------------------------------------------------------- #
# 1.3 Rejections before any fit
# --------------------------------------------------------------------------- #


@pytest.mark.parametrize("value", PREDICTIVE)
@pytest.mark.parametrize(
    ("method", "integration_params"),
    [
        ("pca", {"n_components": 3}),
        ("concat", {}),
        ("snf", {}),
        ("pls_within_group", {}),
        ("pls_residual_pca", {}),
    ],
)
def test_other_integration_methods_are_rejected(dataset, monkeypatch, value, method, integration_params) -> None:
    monkeypatch.setattr(evaluation_module, "fit_omics_preprocessor", _no_preprocessing)
    params = _params(
        integration_method=method,
        integration_params={"layers": FAST_CV["layers"], **integration_params},
        latent_coordinates=value,
    )
    match = rf"latent_coordinates='{value}' is specific to integration_method='pls'; got it with '{method}'"
    with pytest.raises(SimulationEvaluationError, match=match):
        evaluate_semisynthetic_trajectory(dataset, params)
    with pytest.raises(SimulationEvaluationError, match=match):
        integrate_semisynthetic_dataset(dataset, params)


@pytest.mark.parametrize("value", PREDICTIVE)
def test_attribution_is_rejected(dataset, monkeypatch, value) -> None:
    monkeypatch.setattr(evaluation_module, "fit_omics_preprocessor", _no_preprocessing)
    params = _params(
        latent_coordinates=value, attribution=AttributionDiagnosticSettings(enabled=True, bootstrap_replicates=0)
    )
    with pytest.raises(
        SimulationEvaluationError, match=rf"latent_coordinates='{value}' cannot be combined with attribution"
    ):
        evaluate_semisynthetic_trajectory(dataset, params)


@pytest.mark.parametrize("value", PREDICTIVE)
def test_two_stages_are_rejected(two_stage_dataset, monkeypatch, value) -> None:
    monkeypatch.setattr(evaluation_module, "fit_omics_preprocessor", _no_preprocessing)
    params = _params(integration_params={}, latent_coordinates=value)
    match = rf"latent_coordinates='{value}' needs at least three stages.*has 2 stage\(s\)"
    with pytest.raises(SimulationEvaluationError, match=match):
        evaluate_semisynthetic_trajectory(two_stage_dataset, params)
    with pytest.raises(SimulationEvaluationError, match=match):
        integrate_semisynthetic_dataset(two_stage_dataset, params)


def test_orthonormal_coordinates_still_run_on_two_stages(two_stage_dataset) -> None:
    params = _params(integration_params={"forced_components": 2}, latent_coordinates="orthonormal")
    assert integrate_semisynthetic_dataset(two_stage_dataset, params).matrix.shape[1] == 2


# --------------------------------------------------------------------------- #
# 1.4 The measurement in predictive coordinates
# --------------------------------------------------------------------------- #


def test_predictive_equals_the_scores_at_rank_k_minus_1(dataset) -> None:
    scores_params = _forced(2, permutations=49, seed=11)
    scores = evaluate_semisynthetic_trajectory(dataset, scores_params)
    predictive = evaluate_semisynthetic_trajectory(dataset, replace(scores_params, latent_coordinates="predictive"))
    for statistic, value in scores.pair_statistics.items():
        assert predictive.pair_statistics[statistic] == pytest.approx(value, abs=1e-10)
    assert predictive.p_values == scores.p_values


def test_statistics_do_not_depend_on_the_orthonormal_basis(dataset, cross_validated) -> None:
    _, built = cross_validated
    latent = built["predictive_orthonormal"]
    dimension = latent.matrix.shape[1]
    rotation, _ = np.linalg.qr(np.random.default_rng(7).normal(size=(dimension, dimension)))
    rotated = replace(
        latent,
        matrix=pd.DataFrame(
            latent.matrix.to_numpy() @ rotation, index=latent.matrix.index, columns=latent.matrix.columns
        ),
    )
    params = _params(latent_coordinates="predictive_orthonormal", permutations=49, seed=11)
    plain = evaluate_semisynthetic_trajectory(dataset, params, latent=latent)
    turned = evaluate_semisynthetic_trajectory(dataset, params, latent=rotated)
    for statistic, value in plain.pair_statistics.items():
        assert turned.pair_statistics[statistic] == pytest.approx(value, abs=1e-10)
    assert turned.p_values == plain.p_values


def test_span_capture_and_the_population_checkpoint(dataset, forced) -> None:
    rank, scores, built = forced
    on_scores = evaluate_semisynthetic_trajectory(dataset, _forced(rank, span_capture=True), latent=scores)
    preprocessor = fit_omics_preprocessor(dataset, tuple(FAST_CV["layers"]))
    population = concatenate_blocks(preprocessor.transform_population(dataset.population_trajectories))
    groups, stages = ["A", "B"], sorted(dataset.metadata["stage"].astype(str).unique())
    means = np.vstack([population.loc[(group, stage)].to_numpy(dtype=float) for group in groups for stage in stages])
    model = scores.artifacts.model
    V = predictive_directions(model)
    through = {
        "predictive": np.asarray(model.transform(means)) @ V,
        "predictive_orthonormal": built["predictive_orthonormal"].linear_map.inputs(means)
        @ built["predictive_orthonormal"].linear_map.basis,
        "predictive_fitted": np.asarray(model.transform(means)) @ model.y_loadings_.T,
    }
    shares = {}
    for value, latent in built.items():
        result = evaluate_semisynthetic_trajectory(
            dataset, _forced(rank, span_capture=True, latent_coordinates=value), latent=latent
        )
        shares[value] = result.realized_geometry["span_capture"]
        expected = geometry_from_means(through[value], groups, stages, orientation_estimator="net_displacement")
        _assert_matches(
            result.realized_geometry["checkpoints"][POPULATION_LATENT_CHECKPOINT]["joint"],
            json.loads(json.dumps(expected.to_dict())),
        )
    assert shares["predictive"] == shares["predictive_orthonormal"] == shares["predictive_fitted"]
    assert shares["predictive"]["rank"] == 2
    for key in ("group_difference", "pooled_stage"):
        assert shares["predictive"][key] <= on_scores.realized_geometry["span_capture"][key] + 1e-12


def test_latent_checkpoint_and_spectrum_come_from_the_predictive_matrix(dataset, cross_validated) -> None:
    scores, built = cross_validated
    latent = built["predictive"]
    result = evaluate_semisynthetic_trajectory(dataset, _params(latent_coordinates="predictive"), latent=latent)
    preprocessor = fit_omics_preprocessor(dataset, tuple(FAST_CV["layers"]))
    expected = calculate_realized_geometry(
        dataset,
        preprocessor=preprocessor,
        observed_blocks=preprocessor.transform_dataset(dataset),
        latent_matrix=latent.matrix,
        orientation_estimator="net_displacement",
    ).to_dict()
    assert result.realized_geometry["checkpoints"]["pls_latent"] == expected["checkpoints"]["pls_latent"]
    design = build_simulation_trajectory_design(dataset.metadata)
    *_, spectra = estimate_difference(
        latent.matrix,
        design.model_full,
        design.ls_means,
        design.contrast,
        return_spectra=True,
        orientation_estimator="net_displacement",
    )
    assert result.config_spectrum["pooled"] == spectra["pooled"]
    assert result.runtime_metadata["latent_coordinates"] == "predictive"
    assert result.latent_matrix_metadata["integration_params"]["latent_coordinates"] == "predictive"


@pytest.mark.parametrize(
    "integration_params",
    [
        {**FAST_CV, "forced_components": 4},
        {**FAST_CV, "block_weighting": "equal_variance"},
        {**FAST_CV, "sparsity": {"scope": "joint", "keep_fractions": [0.5]}, "forced_components": 3},
    ],
    ids=["forced_rank", "block_weighting", "sparsity"],
)
def test_options_run_under_the_setting(dataset, integration_params) -> None:
    params = _params(integration_params=integration_params, latent_coordinates="predictive", permutations=9, seed=2)
    result = evaluate_semisynthetic_trajectory(dataset, params)
    assert result.latent_matrix_metadata["measured_dimension"] == 2
    assert set(result.p_values) == {"delta", "angle", "shape"}
    if "forced_components" in integration_params:
        assert result.latent_matrix_metadata["selected_lv"] == integration_params["forced_components"]


@pytest.mark.parametrize("name", ["pls", "pls_forced"])
def test_evaluations_without_the_setting_match_pre_change(name: str) -> None:
    _assert_matches(json.loads(json.dumps(evaluation_summary(name))), FIXTURE["evaluations"][name])


# --------------------------------------------------------------------------- #
# 1.5 The prebuilt fit
# --------------------------------------------------------------------------- #


@pytest.mark.parametrize("value", PREDICTIVE)
def test_prebuilt_fit_equals_the_fit_from_scratch(dataset, cross_validated, monkeypatch, value) -> None:
    scores, built = cross_validated

    def no_fit(*args, **kwargs):
        raise AssertionError("the prebuilt path must not fit a PLS model")

    monkeypatch.setattr(evaluation_module, "fit_plsda_model", no_fit)
    monkeypatch.setattr(evaluation_module, "plsda_doubleCV", no_fit)
    reused = pls_coordinates_from_pls(dataset, scores, _params(latent_coordinates=value))
    assert reused.matrix.to_numpy().tobytes() == built[value].matrix.to_numpy().tobytes()
    assert reused.matrix.index.equals(built[value].matrix.index)
    assert list(reused.matrix.columns) == list(built[value].matrix.columns)
    assert reused.metadata == built[value].metadata
    assert json.dumps(_to_jsonable(reused.metadata)) == json.dumps(_to_jsonable(built[value].metadata))
    assert reused.artifacts is scores.artifacts


def test_three_values_from_one_fit(dataset, monkeypatch) -> None:
    fits = []
    original = evaluation_module.fit_plsda_model

    def counting(*args, **kwargs):
        fits.append(1)
        return original(*args, **kwargs)

    monkeypatch.setattr(evaluation_module, "fit_plsda_model", counting)
    scores = integrate_semisynthetic_dataset(dataset, _params())
    for value in PREDICTIVE:
        latent = pls_coordinates_from_pls(dataset, scores, _params(latent_coordinates=value))
        evaluate_semisynthetic_trajectory(dataset, _params(latent_coordinates=value, span_capture=True), latent=latent)
    evaluate_semisynthetic_trajectory(dataset, _params(span_capture=True), latent=scores)
    assert len(fits) == 1


def test_from_scratch_fits_one_model(dataset, monkeypatch) -> None:
    fits = []
    original = evaluation_module.fit_plsda_model

    def counting(*args, **kwargs):
        fits.append(1)
        return original(*args, **kwargs)

    monkeypatch.setattr(evaluation_module, "fit_plsda_model", counting)
    evaluate_semisynthetic_trajectory(dataset, _params(latent_coordinates="predictive_fitted", span_capture=True))
    assert len(fits) == 1


def test_orthonormal_is_the_generic_helpers_orthonormal_case(dataset, cross_validated) -> None:
    scores, _ = cross_validated
    params = _params(latent_coordinates="orthonormal")
    direct = pls_orthonormal_from_pls(dataset, scores, params)
    generic = pls_coordinates_from_pls(dataset, scores, params)
    assert direct.matrix.to_numpy().tobytes() == generic.matrix.to_numpy().tobytes()
    assert direct.metadata == generic.metadata


def test_prebuilt_fit_refusals(dataset, cross_validated) -> None:
    scores, built = cross_validated
    params = _params(latent_coordinates="predictive")
    with pytest.raises(SimulationEvaluationError, match="already in 'predictive' coordinates"):
        pls_coordinates_from_pls(dataset, built["predictive"], params)
    with pytest.raises(SimulationEvaluationError, match="without artifacts"):
        pls_coordinates_from_pls(dataset, replace(scores, artifacts=None), params)
    concat = integrate_semisynthetic_dataset(
        dataset, _params(integration_method="concat", integration_params={"layers": FAST_CV["layers"]})
    )
    with pytest.raises(SimulationEvaluationError, match="'pls' integration with fitted artifacts; got 'concat'"):
        pls_coordinates_from_pls(dataset, concat, params)
    with pytest.raises(SimulationEvaluationError, match="a latent_coordinates value other than 'scores'"):
        pls_coordinates_from_pls(dataset, scores, _params())
    with pytest.raises(SimulationEvaluationError, match="measures layers"):
        pls_coordinates_from_pls(dataset, scores, _params(latent_coordinates="predictive", integration_params={}))


def test_evaluation_refuses_mismatched_prebuilt_coordinates(dataset, cross_validated) -> None:
    scores, built = cross_validated
    with pytest.raises(
        SimulationEvaluationError, match="in 'scores' coordinates, not latent_coordinates='predictive'"
    ):
        evaluate_semisynthetic_trajectory(dataset, _params(latent_coordinates="predictive"), latent=scores)
    with pytest.raises(
        SimulationEvaluationError,
        match="in 'predictive_fitted' coordinates, not latent_coordinates='predictive_orthonormal'",
    ):
        evaluate_semisynthetic_trajectory(
            dataset, _params(latent_coordinates="predictive_orthonormal"), latent=built["predictive_fitted"]
        )


def test_predictive_map_constructor(cross_validated) -> None:
    scores, _ = cross_validated
    model = scores.artifacts.model
    rows = scores.artifacts.features.to_numpy(dtype=float)[:5]
    V = predictive_directions(model)
    maps = {value: LinearLatentMap.from_pls_predictive(model, value) for value in PREDICTIVE}
    np.testing.assert_allclose(maps["predictive"].transform(rows), model.transform(rows) @ V, atol=1e-12)
    orthonormal = maps["predictive_orthonormal"]
    np.testing.assert_allclose(orthonormal.transform(rows), orthonormal.inputs(rows) @ orthonormal.basis, atol=1e-12)
    np.testing.assert_allclose(
        maps["predictive_fitted"].transform(rows), model.transform(rows) @ model.y_loadings_.T, atol=1e-12
    )
    assert all(built.method == "pls" for built in maps.values())
    with pytest.raises(ValueError, match="Unknown predictive coordinates 'orthonormal'"):
        LinearLatentMap.from_pls_predictive(model, "orthonormal")
