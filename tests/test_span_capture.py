"""Span-capture diagnostics and the stage-span guard in the evaluation harness."""

from __future__ import annotations

from dataclasses import replace

import numpy as np
import pytest

import motco.simulations.evaluation as evaluation_module
from motco.simulations.diagnostics import (
    POPULATION_LATENT_CHECKPOINT,
    LinearLatentMap,
    SpanCapture,
    geometry_from_means,
    measure_span_capture,
    span_share,
)
from motco.simulations.evaluation import (
    SimulationEvaluationError,
    SimulationEvaluationParams,
    evaluate_semisynthetic_trajectory,
)
from motco.simulations.grid import make_simulation_cell, parameter_signature
from motco.simulations.preprocessing import concatenate_blocks, fit_omics_preprocessor
from motco.simulations.reference import load_reference
from motco.simulations.semisynthetic import (
    SemiSyntheticTrajectoryParams,
    generate_semisynthetic_trajectory,
)
from motco.stats.pls import fit_plsda_model

TWO_BLOCKS = ["methylation", "expression"]
FAST_CV = {"layers": TWO_BLOCKS, "n_repeats": 1, "max_components": 6}


@pytest.fixture(scope="module")
def reference():
    return load_reference()


def _params(**overrides) -> SimulationEvaluationParams:
    base = {
        "integration_method": "pls",
        "integration_params": FAST_CV,
        "permutations": 0,
        "span_capture": True,
    }
    return SimulationEvaluationParams(**{**base, **overrides})


def _dataset(reference, *, construction=None, mode="none", effect=0.0, n=160, seed=8):
    return generate_semisynthetic_trajectory(
        SemiSyntheticTrajectoryParams(
            seed=seed,
            n_samples=n,
            p_dmp=0.1,
            trajectory_mode=mode,
            group_effect_size=effect,
            construction=construction,
            span_layers=tuple(TWO_BLOCKS) if construction else None,
        ),
        reference=reference,
    )


# --------------------------------------------------------------------------- #
# 3.1 The setting, its identity, and statistics that do not depend on it
# --------------------------------------------------------------------------- #


def test_false_span_capture_normalizes_to_absent() -> None:
    assert SimulationEvaluationParams(span_capture=False).span_capture is None
    assert SimulationEvaluationParams(span_capture=False) == SimulationEvaluationParams()
    with pytest.raises(SimulationEvaluationError, match="span_capture must be true or false"):
        SimulationEvaluationParams(span_capture="yes")  # type: ignore[arg-type]


def test_absent_span_capture_keeps_the_cell_identity() -> None:
    generator = SemiSyntheticTrajectoryParams(seed=1)
    plain = make_simulation_cell(phase="power_primary", generator_params=generator)
    off = make_simulation_cell(
        phase="power_primary",
        generator_params=generator,
        evaluation_params=SimulationEvaluationParams(span_capture=False),
    )
    on = make_simulation_cell(
        phase="power_primary",
        generator_params=generator,
        evaluation_params=SimulationEvaluationParams(span_capture=True),
    )
    assert off.cell_id == plain.cell_id
    assert parameter_signature(off) == parameter_signature(plain)
    assert on.cell_id != plain.cell_id


def test_span_capture_reuses_the_single_fit_and_changes_no_statistic(reference, monkeypatch) -> None:
    dataset = _dataset(reference, mode="orientation", effect=0.5)
    fits = []
    original = evaluation_module.fit_plsda_model

    def counting(*args, **kwargs):
        fits.append(1)
        return original(*args, **kwargs)

    monkeypatch.setattr(evaluation_module, "fit_plsda_model", counting)
    with_capture = evaluate_semisynthetic_trajectory(dataset, _params(permutations=19, seed=5))
    assert len(fits) == 1
    without = evaluate_semisynthetic_trajectory(dataset, _params(permutations=19, seed=5, span_capture=None))

    assert with_capture.pair_statistics == without.pair_statistics
    assert with_capture.p_values == without.p_values
    assert with_capture.latent_matrix_metadata == without.latent_matrix_metadata
    assert with_capture.null_summary == without.null_summary
    assert with_capture.runtime_metadata["span_capture"] is True
    assert "span_capture" not in without.runtime_metadata
    assert "span_capture" not in without.realized_geometry
    assert POPULATION_LATENT_CHECKPOINT not in without.realized_geometry["checkpoints"]
    extended = dict(with_capture.realized_geometry)
    extended.pop("span_capture")
    extended["checkpoints"] = {
        key: value for key, value in extended["checkpoints"].items() if key != POPULATION_LATENT_CHECKPOINT
    }
    assert extended == without.realized_geometry


def test_latent_matrix_is_identical_with_and_without_span_capture(reference) -> None:
    dataset = _dataset(reference, mode="magnitude", effect=0.1)
    results = []
    for capture in (True, None):
        params = _params(span_capture=capture)
        latent = evaluation_module.integrate_semisynthetic_dataset(dataset, params)
        results.append(latent.matrix)
    assert results[0].equals(results[1])


# --------------------------------------------------------------------------- #
# 3.2 The population_latent checkpoint
# --------------------------------------------------------------------------- #


@pytest.mark.parametrize(
    ("method", "integration_params"),
    [("snf", {"layers": TWO_BLOCKS}), ("pls_within_group", {"layers": TWO_BLOCKS})],
)
def test_span_capture_is_inert_for_other_integrations(reference, method, integration_params) -> None:
    dataset = _dataset(reference, mode="magnitude", effect=0.1)
    result = evaluate_semisynthetic_trajectory(
        dataset, _params(integration_method=method, integration_params=integration_params)
    )
    assert POPULATION_LATENT_CHECKPOINT not in result.realized_geometry["checkpoints"]
    assert "span_capture" not in result.realized_geometry


def test_checkpoint_is_absent_without_population_trajectories(reference) -> None:
    dataset = replace(_dataset(reference, mode="magnitude", effect=0.1), population_trajectories=None)
    result = evaluate_semisynthetic_trajectory(dataset, _params())
    assert POPULATION_LATENT_CHECKPOINT not in result.realized_geometry["checkpoints"]
    assert "span_capture" not in result.realized_geometry


def test_population_latent_measures_the_truth_through_the_fitted_basis(reference) -> None:
    dataset = _dataset(reference, mode="orientation", effect=0.5)
    result = evaluate_semisynthetic_trajectory(dataset, _params(orientation_estimator="net_displacement"))
    checkpoint = result.realized_geometry["checkpoints"][POPULATION_LATENT_CHECKPOINT]["joint"]
    assert set(checkpoint) == {"path_lengths", "delta", "angle", "shape", "availability"}
    assert set(checkpoint["path_lengths"]) == {"A", "B"}
    assert checkpoint["availability"] == {"delta": True, "angle": True, "shape": True}


@pytest.fixture(scope="module")
def large_n(reference):
    """Stage-span magnitude, orientation and shape evaluated at n = 800."""

    effects = {"magnitude": 0.1, "orientation": 0.5, "shape": 0.5}
    out = {}
    for mode, effect in effects.items():
        dataset = _dataset(reference, construction="stage_span", mode=mode, effect=effect, n=800)
        params = _params()
        out[mode] = (dataset, evaluate_semisynthetic_trajectory(dataset, params))
    return out


def test_stage_span_orientation_is_recovered_at_large_n(large_n) -> None:
    dataset, result = large_n["orientation"]
    checkpoint = result.realized_geometry["checkpoints"][POPULATION_LATENT_CHECKPOINT]["joint"]
    assert dataset.truth["stage_span"]["intended_geometry"]["angle"] == pytest.approx(45.0)
    assert checkpoint["angle"] == pytest.approx(45.0, abs=3.0)


def test_stage_span_magnitude_ratio_is_recovered_at_large_n(large_n) -> None:
    _, result = large_n["magnitude"]
    lengths = result.realized_geometry["checkpoints"][POPULATION_LATENT_CHECKPOINT]["joint"]["path_lengths"]
    assert lengths["B"] / lengths["A"] == pytest.approx(1.1, abs=1e-6)


# --------------------------------------------------------------------------- #
# 3.3 The span_capture block
# --------------------------------------------------------------------------- #


@pytest.mark.parametrize(
    ("construction", "mode", "effect", "constant"),
    [
        (None, "none", 0.0, True),
        (None, "translation", 1.0, True),
        ("stage_span", "translation", 1.0, True),
        (None, "orientation", 1.0, False),
        (None, "magnitude", 0.1, False),
        ("stage_span", "shape", 0.5, False),
    ],
)
def test_span_capture_shares(reference, construction, mode, effect, constant) -> None:
    dataset = _dataset(reference, construction=construction, mode=mode, effect=effect)
    capture = evaluate_semisynthetic_trajectory(dataset, _params()).realized_geometry["span_capture"]
    assert set(capture) == {"group_difference", "pooled_stage", "rank", "n_features", "random_share"}
    assert capture["n_features"] == 367 + 131
    assert capture["random_share"] == pytest.approx(capture["rank"] / capture["n_features"])
    assert 0.0 <= capture["pooled_stage"] <= 1.0
    if constant:
        assert capture["group_difference"] is None
    else:
        assert 0.0 <= capture["group_difference"] <= 1.0


@pytest.mark.parametrize("mode", ["magnitude", "orientation", "shape"])
def test_stage_span_plant_is_captured_at_large_n(large_n, mode) -> None:
    _, result = large_n[mode]
    capture = result.realized_geometry["span_capture"]
    assert capture["group_difference"] > 0.9
    assert capture["pooled_stage"] > 0.9


def test_share_reproduces_the_plasmode_latent_capture_definition() -> None:
    rng = np.random.default_rng(0)
    X = rng.normal(size=(60, 40))
    y = np.repeat(["0", "1", "2"], 20)
    X[y == "1", :5] += 1.0
    X[y == "2", 5:10] += 1.0
    model = fit_plsda_model(X, y, n_components=3)
    direction = rng.normal(size=40)
    # results/plasmode-2026-11-02/diagnostic/latent_capture.py:
    basis, _ = np.linalg.qr(np.asarray(model.x_rotations_))
    plasmode_share = float(np.linalg.norm(basis.T @ direction) / np.linalg.norm(direction))
    assert span_share(direction, basis) == pytest.approx(plasmode_share**2)
    # Through measure_span_capture: B − A varies along ±direction across two stages.
    means = np.vstack([np.zeros(40), np.zeros(40), -direction, direction]) * np.asarray(model._x_std)
    means = means + np.asarray(model._x_mean)
    capture = measure_span_capture(model, means)
    assert capture.group_difference == pytest.approx(plasmode_share**2)
    assert capture.rank == 3
    assert capture.random_share == pytest.approx(3 / 40)


# --------------------------------------------------------------------------- #
# 3.4 The stage-span guard
# --------------------------------------------------------------------------- #


def test_stage_span_dataset_on_other_blocks_is_refused_before_any_fit(reference, monkeypatch) -> None:
    dataset = _dataset(reference, construction="stage_span", mode="orientation", effect=0.5)

    def no_fit(*args, **kwargs):
        raise AssertionError("no fit may be attempted")

    for name in ("fit_omics_preprocessor", "plsda_doubleCV", "fit_plsda_model"):
        monkeypatch.setattr(evaluation_module, name, no_fit)
    params = _params(integration_params={key: value for key, value in FAST_CV.items() if key != "layers"})
    with pytest.raises(SimulationEvaluationError) as excinfo:
        evaluate_semisynthetic_trajectory(dataset, params)
    message = str(excinfo.value)
    assert "['methylation', 'expression']" in message
    assert "['methylation', 'expression', 'proteomics']" in message
    with pytest.raises(SimulationEvaluationError, match="span layers"):
        evaluation_module.integrate_semisynthetic_dataset(dataset, params)


def test_stage_span_dataset_on_its_span_layers_is_accepted(reference) -> None:
    dataset = _dataset(reference, construction="stage_span", mode="orientation", effect=0.5)
    result = evaluate_semisynthetic_trajectory(dataset, _params(span_capture=None))
    assert result.truth_metadata["construction"] == "stage_span"


# --------------------------------------------------------------------------- #
# Span capture on every linear space (latent-space comparison, 0.16.0)
# --------------------------------------------------------------------------- #

LINEAR_METHODS = {
    "pls": FAST_CV,
    "pca": {"layers": TWO_BLOCKS, "n_components": 3},
    "pls_residual_pca": FAST_CV,
    "concat": {"layers": TWO_BLOCKS},
}


def _population_rows(dataset, preprocessor) -> np.ndarray:
    population = concatenate_blocks(preprocessor.transform_population(dataset.population_trajectories))
    return np.vstack(
        [np.asarray(population.loc[(group, stage)], dtype=float) for group in ("A", "B") for stage in ("0", "1", "2")]
    )


def _frozen_0150_span_capture(model, means) -> SpanCapture:
    """``measure_span_capture`` as released in 0.15.0, kept to pin the PLS records."""

    rotations = np.asarray(model.x_rotations_, dtype=float)
    basis, _ = np.linalg.qr(rotations)
    inputs = (np.asarray(means, dtype=float) - np.asarray(model._x_mean)) / np.asarray(model._x_std)
    n_stages = inputs.shape[0] // 2
    difference = inputs[n_stages:] - inputs[:n_stages]
    varying = difference - difference.mean(axis=0)
    pooled = (inputs[:n_stages] + inputs[n_stages:]) / 2.0
    pooled = pooled - pooled.mean(axis=0)
    scale = max(float(np.linalg.norm(difference)), float(np.linalg.norm(pooled)), np.finfo(float).tiny)
    constant = float(np.linalg.norm(varying)) <= 1e-9 * scale
    return SpanCapture(
        group_difference=None if constant else span_share(varying, basis),
        pooled_stage=span_share(pooled, basis),
        rank=int(rotations.shape[1]),
        n_features=int(rotations.shape[0]),
        random_share=rotations.shape[1] / rotations.shape[0],
    )


def test_pls_map_reproduces_the_model(reference) -> None:
    dataset = _dataset(reference, mode="orientation", effect=0.5)
    latent = evaluation_module.integrate_semisynthetic_dataset(dataset, _params())
    model, linear_map = latent.artifacts.model, latent.linear_map
    assert linear_map.method == "pls"
    rows = latent.artifacts.features.to_numpy()
    assert linear_map.transform(rows).tobytes() == model.transform(rows).tobytes()
    np.testing.assert_allclose(linear_map.transform(rows), latent.matrix.to_numpy(), atol=1e-10)
    np.testing.assert_allclose(linear_map.basis.T @ linear_map.basis, np.eye(model.x_rotations_.shape[1]), atol=1e-12)


def test_pca_map_reproduces_its_latent_matrix(reference) -> None:
    dataset = _dataset(reference, mode="orientation", effect=0.5)
    params = _params(integration_method="pca", integration_params=LINEAR_METHODS["pca"])
    latent = evaluation_module.integrate_semisynthetic_dataset(dataset, params)
    linear_map = latent.linear_map
    preprocessor = fit_omics_preprocessor(dataset, TWO_BLOCKS)
    rows = concatenate_blocks(preprocessor.transform_dataset(dataset)).to_numpy()
    np.testing.assert_allclose(linear_map.transform(rows), latent.matrix.to_numpy(), atol=1e-9)
    np.testing.assert_allclose(linear_map.basis.T @ linear_map.basis, np.eye(3), atol=1e-12)
    assert np.all(linear_map.scale == 1.0)


@pytest.mark.parametrize(("construction", "mode", "effect"), [
    (None, "orientation", 1.0), (None, "none", 0.0), ("stage_span", "shape", 0.5), ("stage_span", "orientation", 0.5),
])
def test_pls_span_capture_records_are_unchanged(reference, construction, mode, effect) -> None:
    """The PLS checkpoint and shares equal the 0.15.0 computation exactly, bit for bit."""

    dataset = _dataset(reference, construction=construction, mode=mode, effect=effect)
    params = _params()
    result = evaluate_semisynthetic_trajectory(dataset, params)
    latent = evaluation_module.integrate_semisynthetic_dataset(dataset, params)
    model = latent.artifacts.model
    rows = _population_rows(dataset, fit_omics_preprocessor(dataset, TWO_BLOCKS))
    expected_capture = _frozen_0150_span_capture(model, rows).to_dict()
    expected_checkpoint = geometry_from_means(
        np.asarray(model.transform(rows), dtype=float), ["A", "B"], ["0", "1", "2"],
        orientation_estimator="net_displacement",
    ).to_dict()
    assert result.realized_geometry["span_capture"] == expected_capture
    assert result.realized_geometry["checkpoints"][POPULATION_LATENT_CHECKPOINT]["joint"] == expected_checkpoint
    # The model-shorthand entry point reads the same map.
    assert measure_span_capture(model, rows) == measure_span_capture(LinearLatentMap.from_pls(model), rows)


def test_pca_shares_match_the_projection_definition(reference) -> None:
    dataset = _dataset(reference, construction="stage_span", mode="orientation", effect=0.5)
    params = _params(integration_method="pca", integration_params=LINEAR_METHODS["pca"])
    result = evaluate_semisynthetic_trajectory(dataset, params)
    latent = evaluation_module.integrate_semisynthetic_dataset(dataset, params)
    preprocessor = fit_omics_preprocessor(dataset, TWO_BLOCKS)
    rows = _population_rows(dataset, preprocessor)
    # The known basis: the top three principal axes of the standardized matrix.
    values = concatenate_blocks(preprocessor.transform_dataset(dataset)).to_numpy()
    _, _, right = np.linalg.svd(values - values.mean(axis=0), full_matrices=False)
    axes = right[:3].T
    difference = rows[3:] - rows[:3]
    varying = difference - difference.mean(axis=0)
    pooled = (rows[:3] + rows[3:]) / 2.0
    pooled = pooled - pooled.mean(axis=0)
    capture = result.realized_geometry["span_capture"]
    expected = np.sum((varying @ axes) ** 2) / np.sum(varying**2)
    assert capture["group_difference"] == pytest.approx(expected, rel=1e-9)
    assert capture["pooled_stage"] == pytest.approx(np.sum((pooled @ axes) ** 2) / np.sum(pooled**2), rel=1e-9)
    assert capture["rank"] == 3 and capture["n_features"] == 498
    checkpoint = result.realized_geometry["checkpoints"][POPULATION_LATENT_CHECKPOINT]["joint"]
    expected_geometry = geometry_from_means(
        latent.linear_map.transform(rows), ["A", "B"], ["0", "1", "2"], orientation_estimator="net_displacement"
    ).to_dict()
    assert checkpoint == expected_geometry


@pytest.mark.parametrize(("mode", "effect"), [("orientation", 1.0), ("translation", 1.0)])
def test_concat_keeps_everything(reference, mode, effect) -> None:
    dataset = _dataset(reference, mode=mode, effect=effect)
    params = _params(integration_method="concat", integration_params=LINEAR_METHODS["concat"])
    geometry = evaluate_semisynthetic_trajectory(dataset, params).realized_geometry
    capture = geometry["span_capture"]
    assert capture["pooled_stage"] == 1.0
    assert capture["random_share"] == 1.0
    assert capture["rank"] == capture["n_features"] == 498
    if mode == "translation":
        assert capture["group_difference"] is None
    else:
        assert capture["group_difference"] == 1.0
    # Through the identity map, the checkpoint is the standardized population joint scope.
    assert (
        geometry["checkpoints"][POPULATION_LATENT_CHECKPOINT]["joint"]
        == geometry["checkpoints"]["population_standardized"]["joint"]
    )


def test_unstandardized_concat_records_no_map(reference) -> None:
    dataset = _dataset(reference, mode="orientation", effect=1.0)
    params = _params(integration_method="concat", integration_params={"layers": TWO_BLOCKS, "standardize": False})
    geometry = evaluate_semisynthetic_trajectory(dataset, params).realized_geometry
    assert "span_capture" not in geometry
    assert POPULATION_LATENT_CHECKPOINT not in geometry["checkpoints"]


@pytest.mark.parametrize("method", sorted([*LINEAR_METHODS, "pls_within_group"]))
def test_statistics_do_not_depend_on_span_capture(reference, method) -> None:
    dataset = _dataset(reference, construction="stage_span", mode="orientation", effect=0.5)
    integration_params = LINEAR_METHODS.get(method, {"layers": TWO_BLOCKS})
    results = [
        evaluate_semisynthetic_trajectory(
            dataset,
            _params(integration_method=method, integration_params=integration_params, permutations=19, seed=3,
                    span_capture=capture),
        )
        for capture in (True, None)
    ]
    on, off = results
    assert on.pair_statistics == off.pair_statistics
    assert on.p_values == off.p_values
    assert on.null_summary == off.null_summary
    assert on.latent_matrix_metadata == off.latent_matrix_metadata
    assert on.config_spectrum == off.config_spectrum
    stripped = dict(on.realized_geometry)
    stripped.pop("span_capture", None)
    stripped["checkpoints"] = {
        key: value for key, value in stripped["checkpoints"].items() if key != POPULATION_LATENT_CHECKPOINT
    }
    assert stripped == off.realized_geometry
