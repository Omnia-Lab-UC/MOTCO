"""Diagnostic ``pca`` / ``pls_within_group`` integration and the orientation-estimator parameter.

See ``openspec/changes/archive/2026-10-12-orientation-measurement-probe/`` (the
``simulation-evaluation-harness`` delta). Default-path byte identity is pinned
separately by ``tests/test_orientation_probe_byte_identity.py``.
"""

from __future__ import annotations

from dataclasses import replace

import numpy as np
import pytest

from motco.simulations import (
    SemiSyntheticTrajectoryParams,
    SimulationEvaluationError,
    SimulationEvaluationParams,
    evaluate_semisynthetic_trajectory,
    generate_semisynthetic_trajectory,
    integrate_semisynthetic_dataset,
    load_reference,
)
from motco.simulations.evaluation import (
    AttributionDiagnosticSettings,
    build_simulation_trajectory_design,
    within_group_pls_basis,
)
from motco.simulations.grid import _to_jsonable, make_simulation_cell, parameter_signature
from motco.simulations.preprocessing import concatenate_blocks, fit_omics_preprocessor
from motco.stats.trajectory import estimate_difference

LAYERS = {"layers": ["methylation", "expression"]}


@pytest.fixture(scope="module")
def dataset():
    params = SemiSyntheticTrajectoryParams(
        seed=8, trajectory_mode="orientation", n_samples=72, n_stages=3, group_effect_size=0.5, p_dmp=0.1
    )
    return generate_semisynthetic_trajectory(params, reference=load_reference())


def _concat_matrix(dataset):
    preprocessor = fit_omics_preprocessor(dataset, ("methylation", "expression", "proteomics"))
    return concatenate_blocks(preprocessor.transform_dataset(dataset))


# ── Orientation estimator parameter ──────────────────────────────────────────


def test_default_estimator_is_invisible() -> None:
    generator = SemiSyntheticTrajectoryParams(seed=1)
    default = SimulationEvaluationParams(integration_method="pls")
    explicit = SimulationEvaluationParams(integration_method="pls", orientation_estimator="pc1")
    assert explicit.orientation_estimator is None
    assert explicit == default
    assert "orientation_estimator" not in _to_jsonable(default)
    left = make_simulation_cell(phase="p", generator_params=generator, evaluation_params=default)
    right = make_simulation_cell(phase="p", generator_params=generator, evaluation_params=explicit)
    assert left.cell_id == right.cell_id
    assert parameter_signature(left) == parameter_signature(right)


def test_net_displacement_estimator_is_recorded_and_used(dataset) -> None:
    params = SimulationEvaluationParams(
        integration_method="concat", permutations=9, seed=3, orientation_estimator="net_displacement"
    )
    result = evaluate_semisynthetic_trajectory(dataset, params)
    assert result.runtime_metadata["orientation_estimator"] == "net_displacement"
    design = build_simulation_trajectory_design(dataset.metadata)
    _, angles, _ = estimate_difference(
        _concat_matrix(dataset),
        design.model_full,
        design.ls_means,
        design.contrast,
        orientation_estimator="net_displacement",
    )
    assert result.pair_statistics["angle"] == pytest.approx(float(angles[0, 1]))
    default = evaluate_semisynthetic_trajectory(dataset, replace(params, orientation_estimator=None))
    assert "orientation_estimator" not in default.runtime_metadata
    assert result.pair_statistics["delta"] == default.pair_statistics["delta"]
    assert result.pair_statistics["angle"] != pytest.approx(default.pair_statistics["angle"])
    observed = result.realized_geometry["checkpoints"]["observed_standardized"]["joint"]["angle"]
    assert observed != pytest.approx(
        default.realized_geometry["checkpoints"]["observed_standardized"]["joint"]["angle"]
    )
    # The cell identity carries a non-default estimator.
    generator = SemiSyntheticTrajectoryParams(seed=1)
    assert (
        make_simulation_cell(phase="p", generator_params=generator, evaluation_params=params).cell_id
        != make_simulation_cell(
            phase="p", generator_params=generator, evaluation_params=replace(params, orientation_estimator=None)
        ).cell_id
    )


def test_unknown_estimator_is_rejected() -> None:
    with pytest.raises(SimulationEvaluationError, match="pc1.*net_displacement"):
        SimulationEvaluationParams(orientation_estimator="folded")


# ── pca ──────────────────────────────────────────────────────────────────────


def test_pca_is_label_free_principal_component_scores(dataset) -> None:
    params = SimulationEvaluationParams(integration_method="pca", integration_params={"n_components": 4})
    latent = integrate_semisynthetic_dataset(dataset, params)
    X = _concat_matrix(dataset).to_numpy()
    centered = X - X.mean(axis=0)
    U, S, _ = np.linalg.svd(centered, full_matrices=False)
    np.testing.assert_allclose(latent.matrix.to_numpy(), U[:, :4] * S[:4])
    rng = np.random.default_rng(0)
    metadata = dataset.metadata.copy()
    metadata["group"] = rng.permutation(metadata["group"].to_numpy())
    metadata["stage"] = rng.permutation(metadata["stage"].to_numpy())
    permuted = integrate_semisynthetic_dataset(replace(dataset, metadata=metadata), params)
    np.testing.assert_array_equal(permuted.matrix.to_numpy(), latent.matrix.to_numpy())
    assert latent.metadata["integration_role"] == "diagnostic"
    assert latent.metadata["integration_method"] == "pca"
    assert latent.metadata["retained_dimensionality"] == 4


@pytest.mark.parametrize(
    "integration_params", [{}, {"n_components": 0}, {"n_components": 10_000}, {"n_components": 2.5}]
)
def test_pca_without_a_feasible_count_is_rejected(dataset, integration_params) -> None:
    params = SimulationEvaluationParams(integration_method="pca", integration_params=integration_params)
    with pytest.raises(SimulationEvaluationError, match="n_components"):
        integrate_semisynthetic_dataset(dataset, params)


# ── pls_within_group ─────────────────────────────────────────────────────────


def test_within_group_basis_spans_orthogonal_group_progressions() -> None:
    rng = np.random.default_rng(1)
    n_per_cell, n_features = 15, 30
    groups = np.repeat(["A", "B"], 3 * n_per_cell)
    stages = np.tile(np.repeat(["0", "1", "2"], n_per_cell), 2)
    # Each group varies only within its own progression plane: group A along
    # features 0 then 1, group B along 5 then 6.
    X = np.zeros((groups.size, n_features))
    X[groups == "A", :2] = rng.normal(scale=0.3, size=((groups == "A").sum(), 2))
    X[groups == "B", 5:7] = rng.normal(scale=0.3, size=((groups == "B").sum(), 2))
    stage_index = stages.astype(int)
    X[groups == "A", 0] += 3.0 * (stage_index[groups == "A"] >= 1)
    X[groups == "A", 1] += 3.0 * (stage_index[groups == "A"] >= 2)
    X[groups == "B", 5] += 3.0 * (stage_index[groups == "B"] >= 1)
    X[groups == "B", 6] += 3.0 * (stage_index[groups == "B"] >= 2)
    basis = within_group_pls_basis(X, groups, stages, components_per_group=2)
    assert basis.shape == (n_features, 4)
    np.testing.assert_allclose(basis.T @ basis, np.eye(4), atol=1e-10)
    projector = basis @ basis.T
    for group in ("A", "B"):
        means = np.vstack([X[(groups == group) & (stages == s)].mean(axis=0) for s in ("0", "1", "2")])
        displacements = np.diff(means, axis=0)
        residual = displacements - displacements @ projector
        assert np.linalg.norm(residual) < 1e-8 * max(1.0, np.linalg.norm(displacements))


def test_within_group_pls_requires_two_samples_per_stage(dataset) -> None:
    keep = ~(
        (dataset.metadata["group"].astype(str) == "B") & (dataset.metadata["stage"].astype(str) == "2")
    ).to_numpy()
    first_dropped = int(np.flatnonzero(~keep)[0])
    keep[first_dropped] = True  # leave exactly one sample in (B, 2)
    subset = replace(
        dataset,
        metadata=dataset.metadata.loc[keep].reset_index(drop=True),
        methylation=dataset.methylation.loc[keep],
        expression=dataset.expression.loc[keep],
        proteomics=dataset.proteomics.loc[keep],
    )
    params = SimulationEvaluationParams(integration_method="pls_within_group")
    with pytest.raises(SimulationEvaluationError, match=r"group 'B'.*\['2'\]"):
        integrate_semisynthetic_dataset(subset, params)


def test_within_group_pls_metadata(dataset) -> None:
    params = SimulationEvaluationParams(integration_method="pls_within_group", integration_params=LAYERS)
    latent = integrate_semisynthetic_dataset(dataset, params)
    metadata = latent.metadata
    assert metadata["integration_role"] == "diagnostic"
    assert metadata["components_per_group"] == 2
    assert 1 <= metadata["basis_rank"] <= 4
    assert metadata["retained_dimensionality"] == metadata["basis_rank"] == latent.matrix.shape[1]
    assert metadata["integration_params"]["layers"] == ["methylation", "expression"]
    # The projection is pooled-centered.
    np.testing.assert_allclose(latent.matrix.to_numpy().mean(axis=0), 0.0, atol=1e-10)


# ── Guards ───────────────────────────────────────────────────────────────────


@pytest.mark.parametrize("method", ["pca", "pls_within_group"])
def test_attribution_rejects_diagnostic_methods(dataset, method) -> None:
    params = SimulationEvaluationParams(
        integration_method=method,
        integration_params={"n_components": 4},
        attribution=AttributionDiagnosticSettings(enabled=True),
    )
    with pytest.raises(SimulationEvaluationError, match="attribution"):
        evaluate_semisynthetic_trajectory(dataset, params)


@pytest.mark.parametrize("method", ["pca", "pls_within_group"])
def test_forced_components_rejects_diagnostic_methods(dataset, method) -> None:
    params = SimulationEvaluationParams(
        integration_method=method, integration_params={"n_components": 4, "forced_components": 4}
    )
    with pytest.raises(SimulationEvaluationError, match="forced_components"):
        evaluate_semisynthetic_trajectory(dataset, params)


# ── Latent checkpoint ────────────────────────────────────────────────────────


@pytest.mark.parametrize(
    ("method", "integration_params", "checkpoint"),
    [
        ("pca", {"n_components": 4}, "pca_latent"),
        ("pls_within_group", {}, "pls_within_group_latent"),
        ("pls", {"n_repeats": 1}, "pls_latent"),
    ],
)
def test_latent_checkpoint_is_named_for_its_method(dataset, method, integration_params, checkpoint) -> None:
    params = SimulationEvaluationParams(integration_method=method, integration_params=integration_params)
    result = evaluate_semisynthetic_trajectory(dataset, params)
    checkpoints = result.realized_geometry["checkpoints"]
    latent_names = {name for name in checkpoints if name.endswith("_latent")}
    assert latent_names == {checkpoint}
    assert checkpoints[checkpoint]["joint"]["delta"] == pytest.approx(result.pair_statistics["delta"])


@pytest.mark.parametrize("method", ["concat", "snf"])
def test_concat_and_snf_carry_no_latent_checkpoint(dataset, method) -> None:
    result = evaluate_semisynthetic_trajectory(dataset, SimulationEvaluationParams(integration_method=method))
    assert not [name for name in result.realized_geometry["checkpoints"] if name.endswith("_latent")]
