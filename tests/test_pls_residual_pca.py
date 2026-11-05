"""The diagnostic ``pls_residual_pca`` integration: production PLS scores plus stage-residual components."""

from __future__ import annotations

from dataclasses import replace

import numpy as np
import pytest

import motco.simulations.evaluation as evaluation_module
from motco.simulations.diagnostics import POPULATION_LATENT_CHECKPOINT
from motco.simulations.evaluation import (
    AttributionDiagnosticSettings,
    SimulationEvaluationError,
    SimulationEvaluationParams,
    evaluate_semisynthetic_trajectory,
    integrate_semisynthetic_dataset,
    pls_residual_pca_from_pls,
)
from motco.simulations.reference import load_reference
from motco.simulations.semisynthetic import SemiSyntheticTrajectoryParams, generate_semisynthetic_trajectory

TWO_BLOCKS = ["methylation", "expression"]
FAST_CV = {"layers": TWO_BLOCKS, "n_repeats": 1, "max_components": 6}


@pytest.fixture(scope="module")
def reference():
    return load_reference()


@pytest.fixture(scope="module")
def dataset(reference):
    return generate_semisynthetic_trajectory(
        SemiSyntheticTrajectoryParams(
            seed=8,
            n_samples=160,
            p_dmp=0.1,
            trajectory_mode="orientation",
            group_effect_size=0.5,
            construction="stage_span",
            span_layers=tuple(TWO_BLOCKS),
            span_alignment=90.0,
        ),
        reference=reference,
    )


def _params(method: str = "pls_residual_pca", **overrides) -> SimulationEvaluationParams:
    base = {"integration_method": method, "integration_params": FAST_CV, "permutations": 0}
    return SimulationEvaluationParams(**{**base, **overrides})


@pytest.fixture(scope="module")
def integrations(dataset):
    return integrate_semisynthetic_dataset(dataset, _params("pls")), integrate_semisynthetic_dataset(dataset, _params())


# --------------------------------------------------------------------------- #
# 3.1 The space
# --------------------------------------------------------------------------- #


def test_leading_columns_are_the_production_scores(integrations) -> None:
    pls, residual = integrations
    rank = pls.matrix.shape[1]
    assert residual.matrix.iloc[:, :rank].to_numpy().tobytes() == pls.matrix.to_numpy().tobytes()
    assert list(residual.matrix.columns[:rank]) == list(pls.matrix.columns)
    assert list(residual.matrix.columns[rank:]) == ["resid_0", "resid_1"]
    assert residual.matrix.index.equals(pls.matrix.index)


def test_residual_directions_are_orthogonal_to_the_pls_span(integrations) -> None:
    pls, residual = integrations
    rank = pls.matrix.shape[1]
    pls_basis, _ = np.linalg.qr(np.asarray(pls.artifacts.model.x_rotations_))
    axes = residual.linear_map.basis[:, rank:]
    np.testing.assert_allclose(pls_basis.T @ axes, 0.0, atol=1e-10)
    np.testing.assert_allclose(axes.T @ axes, np.eye(axes.shape[1]), atol=1e-10)
    # The map reproduces the residual scores from the standardized joint rows.
    rows = pls.artifacts.features.to_numpy()
    np.testing.assert_allclose(residual.linear_map.transform(rows), residual.matrix.to_numpy(), atol=1e-9)


def test_basis_does_not_use_the_group_labels(dataset, integrations) -> None:
    _, residual = integrations
    metadata = dataset.metadata.copy()
    rng = np.random.default_rng(0)
    for stage in metadata["stage"].unique():
        rows = metadata.index[metadata["stage"] == stage]
        metadata.loc[rows, "group"] = rng.permutation(metadata.loc[rows, "group"].to_numpy())
    assert not metadata["group"].equals(dataset.metadata["group"])
    relabelled = integrate_semisynthetic_dataset(replace(dataset, metadata=metadata), _params())
    assert relabelled.matrix.to_numpy().tobytes() == residual.matrix.to_numpy().tobytes()


def test_prebuilt_pls_fit_equals_the_fit_from_scratch(dataset, integrations, monkeypatch) -> None:
    pls, residual = integrations

    def no_fit(*args, **kwargs):
        raise AssertionError("the prebuilt path must not fit a PLS model")

    monkeypatch.setattr(evaluation_module, "fit_plsda_model", no_fit)
    monkeypatch.setattr(evaluation_module, "plsda_doubleCV", no_fit)
    reused = pls_residual_pca_from_pls(dataset, pls, _params())
    assert reused.matrix.to_numpy().tobytes() == residual.matrix.to_numpy().tobytes()
    assert reused.metadata == residual.metadata
    result = evaluate_semisynthetic_trajectory(dataset, _params(span_capture=True), latent=reused)
    assert result.latent_matrix_metadata == residual.metadata


@pytest.mark.parametrize("value", [0, -1, 10_000, 1.5, "two"])
def test_infeasible_residual_rank_is_rejected(dataset, integrations, value) -> None:
    pls, _ = integrations
    params = _params(integration_params={**FAST_CV, "residual_components": value})
    match = r"residual_components=-?\d+ is outside the feasible range \[1, \d+\]" if isinstance(value, int) else (
        "residual_components must be an integer"
    )
    with pytest.raises(SimulationEvaluationError, match=match):
        pls_residual_pca_from_pls(dataset, pls, params)


def test_residual_rank_is_configurable(dataset, integrations) -> None:
    pls, _ = integrations
    reused = pls_residual_pca_from_pls(dataset, pls, _params(integration_params={**FAST_CV, "residual_components": 4}))
    assert reused.matrix.shape[1] == pls.matrix.shape[1] + 4
    assert reused.metadata["residual_components"] == 4


def test_prebuilt_fit_must_match(dataset, integrations) -> None:
    pls, residual = integrations
    with pytest.raises(SimulationEvaluationError, match="built on a 'pls' integration"):
        pls_residual_pca_from_pls(dataset, residual, _params())
    with pytest.raises(SimulationEvaluationError, match="measures layers"):
        pls_residual_pca_from_pls(dataset, pls, _params(integration_params={"n_repeats": 1}))
    with pytest.raises(SimulationEvaluationError, match="integration_method='pls_residual_pca'"):
        pls_residual_pca_from_pls(dataset, pls, _params("pls"))


# --------------------------------------------------------------------------- #
# 3.2 Diagnostic label, metadata, checkpoint and rejections
# --------------------------------------------------------------------------- #


def test_metadata_records_the_ranks(integrations) -> None:
    pls, residual = integrations
    metadata = residual.metadata
    rank = pls.metadata["selected_lv"]
    assert metadata["integration_method"] == "pls_residual_pca"
    assert metadata["integration_role"] == "diagnostic"
    assert metadata["pls_rank"] == rank
    assert metadata["residual_components"] == 2
    assert metadata["retained_dimensionality"] == rank + 2
    assert metadata["integration_params"] == {**pls.metadata["integration_params"], "residual_components": 2}
    assert metadata["cv_mean_auroc"] == pls.metadata["cv_mean_auroc"]


def test_checkpoint_and_span_capture(dataset, integrations) -> None:
    pls, residual = integrations
    result = evaluate_semisynthetic_trajectory(dataset, _params(span_capture=True), latent=residual)
    checkpoints = result.realized_geometry["checkpoints"]
    assert "pls_residual_pca_latent" in checkpoints and "pls_latent" not in checkpoints
    assert POPULATION_LATENT_CHECKPOINT in checkpoints
    capture = result.realized_geometry["span_capture"]
    assert capture["rank"] == residual.metadata["retained_dimensionality"]
    production = evaluate_semisynthetic_trajectory(dataset, _params("pls", span_capture=True), latent=pls)
    # The residual space contains the PLS span, so it captures at least as much.
    assert capture["group_difference"] >= production.realized_geometry["span_capture"]["group_difference"] - 1e-12
    assert capture["pooled_stage"] >= production.realized_geometry["span_capture"]["pooled_stage"] - 1e-12


def test_attribution_and_forced_components_are_rejected(dataset) -> None:
    with pytest.raises(SimulationEvaluationError, match="require integration_method='pls'"):
        evaluate_semisynthetic_trajectory(
            dataset, _params(attribution=AttributionDiagnosticSettings(enabled=True, bootstrap_replicates=0))
        )
    with pytest.raises(SimulationEvaluationError, match="forced_components"):
        evaluate_semisynthetic_trajectory(dataset, _params(integration_params={**FAST_CV, "forced_components": 3}))


def test_default_residual_rank_follows_the_stage_count(reference) -> None:
    four = generate_semisynthetic_trajectory(
        SemiSyntheticTrajectoryParams(seed=3, n_samples=160, n_stages=4, p_dmp=0.1), reference=reference
    )
    latent = integrate_semisynthetic_dataset(four, _params(integration_params={**FAST_CV, "layers": None}))
    assert latent.metadata["residual_components"] == 3
