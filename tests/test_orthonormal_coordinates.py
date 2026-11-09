"""Orthonormal span coordinates of the PLS measurement (``latent_coordinates``)."""

from __future__ import annotations

import json
import math
from dataclasses import replace
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

import motco.simulations.evaluation as evaluation_module
from motco.simulations.diagnostics import (
    POPULATION_LATENT_CHECKPOINT,
    LinearLatentMap,
    calculate_realized_geometry,
    geometry_from_means,
)
from motco.simulations.evaluation import (
    AttributionDiagnosticSettings,
    LatentIntegrationResult,
    SimulationEvaluationError,
    SimulationEvaluationParams,
    build_simulation_trajectory_design,
    evaluate_semisynthetic_trajectory,
    integrate_semisynthetic_dataset,
    pls_orthonormal_from_pls,
    pls_residual_pca_from_pls,
)
from motco.simulations.grid import _to_jsonable, make_simulation_cell, parameter_signature
from motco.simulations.preprocessing import concatenate_blocks, fit_omics_preprocessor
from motco.simulations.reference import load_reference
from motco.simulations.semisynthetic import generate_semisynthetic_trajectory
from motco.stats.trajectory import estimate_difference
from tests.orthonormal_coordinates_fixture_capture import DATASET, EVALUATIONS, FAST_CV, evaluation_summary
from tests.test_latent_adjustment import _fit_core, _with_nuisance
from tests.test_phase6_byte_identity import _assert_close

FIXTURE = json.loads(
    (Path(__file__).parent / "data" / "orthonormal_coordinates_prechange_fixture.json").read_text()
)

#: Degrees. An angle that is exactly zero in the population (the shape plant keeps
#: the net displacement) is ill-conditioned: arccos near 1 turns rounding of order
#: 1e-16 into about 1e-6 degrees, which differs between BLAS thread counts.
ANGLE_ABS_TOL = 1e-5


def _assert_matches(actual, expected, path: str = "") -> None:
    """``_assert_close``, with angles compared at ``ANGLE_ABS_TOL``."""

    if isinstance(expected, dict):
        assert isinstance(actual, dict) and sorted(actual) == sorted(expected), path
        for key in expected:
            _assert_matches(actual[key], expected[key], f"{path}.{key}")
    elif isinstance(expected, list):
        assert isinstance(actual, list) and len(actual) == len(expected), path
        for index, (a, e) in enumerate(zip(actual, expected, strict=True)):
            _assert_matches(a, e, f"{path}[{index}]")
    elif isinstance(expected, float) and path.endswith(".angle"):
        assert math.isclose(actual, expected, rel_tol=1e-9, abs_tol=ANGLE_ABS_TOL), (path, actual, expected)
    else:
        _assert_close(actual, expected, path)


@pytest.fixture(scope="module")
def dataset():
    return generate_semisynthetic_trajectory(DATASET, reference=load_reference())


def _params(**overrides) -> SimulationEvaluationParams:
    base = {"integration_method": "pls", "integration_params": FAST_CV, "permutations": 0}
    return SimulationEvaluationParams(**{**base, **overrides})


@pytest.fixture(scope="module")
def integrations(dataset):
    scores = integrate_semisynthetic_dataset(dataset, _params())
    orthonormal = integrate_semisynthetic_dataset(dataset, _params(latent_coordinates="orthonormal"))
    return scores, orthonormal


def _inputs(latent: LatentIntegrationResult) -> np.ndarray:
    model = latent.artifacts.model
    return (latent.artifacts.features.to_numpy(dtype=float) - model._x_mean) / model._x_std


# --------------------------------------------------------------------------- #
# 1.1 The setting and its identity
# --------------------------------------------------------------------------- #


def test_scores_normalize_to_absent() -> None:
    assert SimulationEvaluationParams(latent_coordinates="scores").latent_coordinates is None
    assert SimulationEvaluationParams(latent_coordinates="scores") == SimulationEvaluationParams()
    assert SimulationEvaluationParams().resolved_latent_coordinates == "scores"
    assert SimulationEvaluationParams(latent_coordinates="orthonormal").resolved_latent_coordinates == "orthonormal"


def test_unknown_value_is_rejected() -> None:
    with pytest.raises(
        SimulationEvaluationError, match=r"Unknown latent_coordinates 'whitened'; allowed: \['scores', 'orthonormal'\]"
    ):
        SimulationEvaluationParams(latent_coordinates="whitened")


def test_absent_setting_keeps_the_cell_identity() -> None:
    plain = make_simulation_cell(phase="power_primary", generator_params=DATASET)
    scores = make_simulation_cell(
        phase="power_primary",
        generator_params=DATASET,
        evaluation_params=SimulationEvaluationParams(latent_coordinates="scores"),
    )
    orthonormal = make_simulation_cell(
        phase="power_primary",
        generator_params=DATASET,
        evaluation_params=SimulationEvaluationParams(latent_coordinates="orthonormal"),
    )
    assert scores.cell_id == plain.cell_id
    assert parameter_signature(scores) == parameter_signature(plain)
    assert "latent_coordinates" not in _to_jsonable(SimulationEvaluationParams())
    assert orthonormal.cell_id != plain.cell_id
    assert _to_jsonable(orthonormal.evaluation_params)["latent_coordinates"] == "orthonormal"


# --------------------------------------------------------------------------- #
# 1.2 The coordinates in the core
# --------------------------------------------------------------------------- #


def test_orthonormal_matrix_is_the_inputs_on_the_weights(integrations) -> None:
    scores, orthonormal = integrations
    model = orthonormal.artifacts.model
    inputs = _inputs(orthonormal)
    by_weights = inputs @ np.asarray(model.x_weights_)
    signs = np.sign(np.sum(orthonormal.matrix.to_numpy() * by_weights, axis=0))
    np.testing.assert_allclose(orthonormal.matrix.to_numpy(), by_weights * signs, atol=1e-8)
    # The scores are the orthonormal matrix times one square matrix (same span).
    basis = orthonormal.linear_map.basis
    np.testing.assert_allclose(basis.T @ basis, np.eye(basis.shape[1]), atol=1e-10)
    np.testing.assert_allclose(
        scores.matrix.to_numpy(), orthonormal.matrix.to_numpy() @ (basis.T @ model.x_rotations_), atol=1e-8
    )
    assert list(orthonormal.matrix.columns) == list(scores.matrix.columns)
    assert orthonormal.matrix.index.equals(scores.matrix.index)


def test_fit_and_rank_are_unchanged(integrations) -> None:
    scores, orthonormal = integrations
    assert orthonormal.artifacts.model.x_rotations_.tobytes() == scores.artifacts.model.x_rotations_.tobytes()
    assert orthonormal.metadata["integration_params"]["latent_coordinates"] == "orthonormal"
    stripped = dict(orthonormal.metadata)
    stripped["integration_params"] = {
        key: value for key, value in orthonormal.metadata["integration_params"].items() if key != "latent_coordinates"
    }
    assert stripped == scores.metadata
    assert "latent_coordinates" not in scores.metadata["integration_params"]


def test_scores_map_is_unchanged(integrations) -> None:
    scores, orthonormal = integrations
    rows = scores.artifacts.features.to_numpy(dtype=float)
    assert scores.matrix.to_numpy().tobytes() == np.asarray(scores.artifacts.model.x_scores_).tobytes()
    np.testing.assert_allclose(scores.linear_map.transform(rows), scores.matrix.to_numpy(), atol=1e-10)
    np.testing.assert_allclose(orthonormal.linear_map.transform(rows), orthonormal.matrix.to_numpy(), atol=1e-12)
    for attribute in ("center", "scale", "basis"):
        assert getattr(orthonormal.linear_map, attribute).tobytes() == getattr(scores.linear_map, attribute).tobytes()


def test_latent_adjusted_path_projects_the_raw_features() -> None:
    inputs = _with_nuisance()
    scores = _fit_core(inputs, nuisance=inputs["nuisance"])
    orthonormal = _fit_core(inputs, nuisance=inputs["nuisance"], coordinates="orthonormal")
    model = orthonormal.artifacts.model
    assert model.scale is False
    basis, _ = np.linalg.qr(np.asarray(model.x_rotations_))
    X = inputs["X"].to_numpy()
    np.testing.assert_allclose(orthonormal.matrix.to_numpy(), (X - model._x_mean) @ basis, atol=1e-12)
    np.testing.assert_allclose(
        scores.matrix.to_numpy(), orthonormal.matrix.to_numpy() @ (basis.T @ model.x_rotations_), atol=1e-8
    )
    assert orthonormal.metadata["integration_params"]["latent_adjustment"] == "stage"
    assert orthonormal.metadata["selected_lv"] == scores.metadata["selected_lv"]


def test_core_rejects_an_unknown_coordinate_system() -> None:
    inputs = _with_nuisance()
    with pytest.raises(SimulationEvaluationError, match="unknown coordinates 'whitened'"):
        _fit_core(inputs, coordinates="whitened")


# --------------------------------------------------------------------------- #
# 1.3 Rejections before any fit
# --------------------------------------------------------------------------- #


def _no_preprocessing(*args, **kwargs):
    raise AssertionError("the setting must be rejected before any preprocessing")


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
def test_other_integration_methods_are_rejected(dataset, monkeypatch, method, integration_params) -> None:
    monkeypatch.setattr(evaluation_module, "fit_omics_preprocessor", _no_preprocessing)
    params = _params(
        integration_method=method,
        integration_params={"layers": FAST_CV["layers"], **integration_params},
        latent_coordinates="orthonormal",
    )
    match = rf"latent_coordinates='orthonormal' is specific to integration_method='pls'; got it with '{method}'"
    with pytest.raises(SimulationEvaluationError, match=match):
        evaluate_semisynthetic_trajectory(dataset, params)
    with pytest.raises(SimulationEvaluationError, match=match):
        integrate_semisynthetic_dataset(dataset, params)


def test_attribution_is_rejected(dataset, monkeypatch) -> None:
    monkeypatch.setattr(evaluation_module, "fit_omics_preprocessor", _no_preprocessing)
    params = _params(
        latent_coordinates="orthonormal",
        attribution=AttributionDiagnosticSettings(enabled=True, bootstrap_replicates=0),
    )
    with pytest.raises(SimulationEvaluationError, match="attribution is defined only in score coordinates"):
        evaluate_semisynthetic_trajectory(dataset, params)


def test_an_unknown_value_set_after_construction_is_rejected(dataset, monkeypatch) -> None:
    monkeypatch.setattr(evaluation_module, "fit_omics_preprocessor", _no_preprocessing)
    params = _params()
    object.__setattr__(params, "latent_coordinates", "whitened")
    with pytest.raises(SimulationEvaluationError, match="Unknown latent_coordinates 'whitened'"):
        evaluate_semisynthetic_trajectory(dataset, params)


@pytest.mark.parametrize("name", sorted(EVALUATIONS))
def test_evaluations_without_the_setting_match_pre_change(name: str) -> None:
    _assert_matches(json.loads(json.dumps(evaluation_summary(name))), FIXTURE["evaluations"][name])


def test_explicit_scores_match_pre_change() -> None:
    explicit = {**EVALUATIONS, "pls": replace(EVALUATIONS["pls"], latent_coordinates="scores")}
    assert explicit["pls"] == EVALUATIONS["pls"]
    _assert_matches(json.loads(json.dumps(evaluation_summary("pls"))), FIXTURE["evaluations"]["pls"])


# --------------------------------------------------------------------------- #
# 1.4 The measurement in orthonormal coordinates
# --------------------------------------------------------------------------- #


def test_statistics_do_not_depend_on_the_orthonormal_basis(dataset, integrations) -> None:
    _, orthonormal = integrations
    rank = orthonormal.matrix.shape[1]
    rotation, _ = np.linalg.qr(np.random.default_rng(4).normal(size=(rank, rank)))
    rotated = replace(
        orthonormal,
        matrix=pd.DataFrame(
            orthonormal.matrix.to_numpy() @ rotation, index=orthonormal.matrix.index, columns=orthonormal.matrix.columns
        ),
    )
    params = _params(latent_coordinates="orthonormal", permutations=49, seed=11)
    plain = evaluate_semisynthetic_trajectory(dataset, params, latent=orthonormal)
    turned = evaluate_semisynthetic_trajectory(dataset, params, latent=rotated)
    for statistic, value in plain.pair_statistics.items():
        assert turned.pair_statistics[statistic] == pytest.approx(value, abs=1e-10)
    assert turned.p_values == plain.p_values


def test_population_checkpoint_and_span_capture(dataset, integrations) -> None:
    scores, orthonormal = integrations
    on_scores = evaluate_semisynthetic_trajectory(dataset, _params(span_capture=True), latent=scores)
    on_span = evaluate_semisynthetic_trajectory(
        dataset, _params(span_capture=True, latent_coordinates="orthonormal"), latent=orthonormal
    )
    assert on_span.realized_geometry["span_capture"] == on_scores.realized_geometry["span_capture"]

    preprocessor = fit_omics_preprocessor(dataset, tuple(FAST_CV["layers"]))
    population = concatenate_blocks(preprocessor.transform_population(dataset.population_trajectories))
    groups, stages = ["A", "B"], sorted(dataset.metadata["stage"].astype(str).unique())
    means = np.vstack([population.loc[(group, stage)].to_numpy(dtype=float) for group in groups for stage in stages])
    basis = orthonormal.linear_map.basis
    expected = geometry_from_means(
        orthonormal.linear_map.inputs(means) @ basis, groups, stages, orientation_estimator="net_displacement"
    )
    measured = on_span.realized_geometry["checkpoints"][POPULATION_LATENT_CHECKPOINT]["joint"]
    _assert_close(measured, json.loads(json.dumps(expected.to_dict())))
    assert measured != on_scores.realized_geometry["checkpoints"][POPULATION_LATENT_CHECKPOINT]["joint"]


def test_latent_checkpoint_and_spectrum_come_from_the_orthonormal_matrix(dataset, integrations) -> None:
    scores, orthonormal = integrations
    params = _params(latent_coordinates="orthonormal")
    result = evaluate_semisynthetic_trajectory(dataset, params, latent=orthonormal)
    preprocessor = fit_omics_preprocessor(dataset, tuple(FAST_CV["layers"]))
    expected = calculate_realized_geometry(
        dataset,
        preprocessor=preprocessor,
        observed_blocks=preprocessor.transform_dataset(dataset),
        latent_matrix=orthonormal.matrix,
        orientation_estimator="net_displacement",
    ).to_dict()
    assert result.realized_geometry["checkpoints"]["pls_latent"] == expected["checkpoints"]["pls_latent"]
    design = build_simulation_trajectory_design(dataset.metadata)
    *_, spectra = estimate_difference(
        orthonormal.matrix,
        design.model_full,
        design.ls_means,
        design.contrast,
        return_spectra=True,
        orientation_estimator="net_displacement",
    )
    assert result.config_spectrum["pooled"] == spectra["pooled"]
    on_scores = evaluate_semisynthetic_trajectory(dataset, _params(), latent=scores)
    assert on_scores.config_spectrum["pooled"] != result.config_spectrum["pooled"]
    assert result.runtime_metadata["latent_coordinates"] == "orthonormal"
    assert "latent_coordinates" not in on_scores.runtime_metadata
    assert result.latent_matrix_metadata["integration_params"]["latent_coordinates"] == "orthonormal"


def test_forced_rank_works_under_the_setting(dataset) -> None:
    params = _params(integration_params={**FAST_CV, "forced_components": 3}, latent_coordinates="orthonormal")
    latent = integrate_semisynthetic_dataset(dataset, params)
    assert latent.metadata["component_selection"] == "forced"
    assert latent.matrix.shape[1] == 3
    basis = latent.linear_map.basis
    np.testing.assert_allclose(latent.matrix.to_numpy(), _inputs(latent) @ basis, atol=1e-12)
    result = evaluate_semisynthetic_trajectory(dataset, replace(params, permutations=9, seed=2), latent=latent)
    assert set(result.p_values) == {"delta", "angle", "shape"}


# --------------------------------------------------------------------------- #
# 1.5 The prebuilt fit
# --------------------------------------------------------------------------- #


def test_prebuilt_fit_equals_the_fit_from_scratch(dataset, integrations, monkeypatch) -> None:
    scores, orthonormal = integrations

    def no_fit(*args, **kwargs):
        raise AssertionError("the prebuilt path must not fit a PLS model")

    monkeypatch.setattr(evaluation_module, "fit_plsda_model", no_fit)
    monkeypatch.setattr(evaluation_module, "plsda_doubleCV", no_fit)
    reused = pls_orthonormal_from_pls(dataset, scores, _params(latent_coordinates="orthonormal"))
    assert reused.matrix.to_numpy().tobytes() == orthonormal.matrix.to_numpy().tobytes()
    assert reused.matrix.index.equals(orthonormal.matrix.index)
    assert reused.metadata == orthonormal.metadata
    assert json.dumps(_to_jsonable(reused.metadata)) == json.dumps(_to_jsonable(orthonormal.metadata))
    assert reused.artifacts is scores.artifacts


def test_from_scratch_fits_one_model(dataset, monkeypatch) -> None:
    fits = []
    original = evaluation_module.fit_plsda_model

    def counting(*args, **kwargs):
        fits.append(1)
        return original(*args, **kwargs)

    monkeypatch.setattr(evaluation_module, "fit_plsda_model", counting)
    evaluate_semisynthetic_trajectory(dataset, _params(latent_coordinates="orthonormal", span_capture=True))
    assert len(fits) == 1


def test_prebuilt_fit_refusals(dataset, integrations) -> None:
    scores, orthonormal = integrations
    params = _params(latent_coordinates="orthonormal")
    with pytest.raises(SimulationEvaluationError, match="already in 'orthonormal' coordinates"):
        pls_orthonormal_from_pls(dataset, orthonormal, params)
    with pytest.raises(SimulationEvaluationError, match="without artifacts"):
        pls_orthonormal_from_pls(dataset, replace(scores, artifacts=None), params)
    concat = integrate_semisynthetic_dataset(
        dataset, _params(integration_method="concat", integration_params={"layers": FAST_CV["layers"]})
    )
    with pytest.raises(SimulationEvaluationError, match="'pls' integration with fitted artifacts; got 'concat'"):
        pls_orthonormal_from_pls(dataset, concat, params)
    with pytest.raises(SimulationEvaluationError, match="latent_coordinates='orthonormal'"):
        pls_orthonormal_from_pls(dataset, scores, _params())
    with pytest.raises(SimulationEvaluationError, match="measures layers"):
        pls_orthonormal_from_pls(dataset, scores, _params(latent_coordinates="orthonormal", integration_params={}))
    with pytest.raises(SimulationEvaluationError, match="built on the PLS scores"):
        pls_residual_pca_from_pls(
            dataset, orthonormal, _params(integration_method="pls_residual_pca", integration_params=FAST_CV)
        )


def test_evaluation_refuses_mismatched_prebuilt_coordinates(dataset, integrations) -> None:
    scores, orthonormal = integrations
    with pytest.raises(
        SimulationEvaluationError, match="in 'scores' coordinates, not latent_coordinates='orthonormal'"
    ):
        evaluate_semisynthetic_trajectory(dataset, _params(latent_coordinates="orthonormal"), latent=scores)
    with pytest.raises(
        SimulationEvaluationError, match="in 'orthonormal' coordinates, not latent_coordinates='scores'"
    ):
        evaluate_semisynthetic_trajectory(dataset, _params(), latent=orthonormal)


def test_orthonormal_map_constructor(integrations) -> None:
    _, orthonormal = integrations
    model = orthonormal.artifacts.model
    built = LinearLatentMap.from_pls_orthonormal(model)
    rows = orthonormal.artifacts.features.to_numpy(dtype=float)[:5]
    np.testing.assert_allclose(built.transform(rows), built.inputs(rows) @ built.basis, atol=1e-14)
    assert built.method == "pls"
