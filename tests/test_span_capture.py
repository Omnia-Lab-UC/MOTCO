"""Span-capture diagnostics and the stage-span guard in the evaluation harness."""

from __future__ import annotations

from dataclasses import replace

import numpy as np
import pytest

import motco.simulations.evaluation as evaluation_module
from motco.simulations.diagnostics import POPULATION_LATENT_CHECKPOINT, measure_span_capture, span_share
from motco.simulations.evaluation import (
    SimulationEvaluationError,
    SimulationEvaluationParams,
    evaluate_semisynthetic_trajectory,
)
from motco.simulations.grid import make_simulation_cell, parameter_signature
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
    [("concat", {"layers": TWO_BLOCKS}), ("pca", {"layers": TWO_BLOCKS, "n_components": 3})],
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
