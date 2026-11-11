"""Sparse stage-PLS in the production core and the harness (``integration_params["sparsity"]``).

The core's metadata, the forced-rank rule, composition with the latent
adjustment, the orthonormal coordinates and block weighting; the harness's
rejections before any preprocessing, the PLS-derived spaces, span capture and
attribution on a sparse fit; and invisibility while the key is absent.
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

import motco.simulations.evaluation as evaluation_module
from motco.simulations.evaluation import (
    AttributionDiagnosticSettings,
    SimulationEvaluationError,
    SimulationEvaluationParams,
    evaluate_semisynthetic_trajectory,
    fit_pls_latent_space,
    integrate_semisynthetic_dataset,
    pls_orthonormal_from_pls,
    pls_residual_pca_from_pls,
)
from motco.simulations.grid import make_simulation_cell, parameter_signature
from motco.simulations.preprocessing import concatenate_blocks, fit_omics_preprocessor
from motco.simulations.reference import load_reference
from motco.simulations.semisynthetic import SemiSyntheticTrajectoryParams, generate_semisynthetic_trajectory
from motco.stats.pls import NuisanceDesign
from motco.stats.sparse_pls import DEFAULT_KEEP_FRACTIONS, SparsePLSDA, SparsityDesign
from tests.orthonormal_coordinates_fixture_capture import evaluation_summary
from tests.test_orthonormal_coordinates import _assert_matches

TWO_BLOCKS = ("methylation", "expression")
COUNTS = {"methylation": 367, "expression": 131}
FAST_CV = {"layers": list(TWO_BLOCKS), "n_repeats": 1, "max_components": 6}
BLOCK = {"scope": "block"}
PRECHANGE_FIXTURE = json.loads(
    (Path(__file__).parent / "data" / "orthonormal_coordinates_prechange_fixture.json").read_text()
)


@pytest.fixture(scope="module")
def dataset():
    return generate_semisynthetic_trajectory(
        SemiSyntheticTrajectoryParams(seed=8, n_samples=120, p_dmp=0.1, trajectory_mode="orientation",
                                      group_effect_size=0.5),
        reference=load_reference(),
    )


def _params(sparsity=None, **overrides) -> SimulationEvaluationParams:
    integration_params = dict(overrides.pop("integration_params", FAST_CV))
    if sparsity is not None:
        integration_params["sparsity"] = sparsity
    base = {"integration_method": "pls", "integration_params": integration_params, "permutations": 0}
    return SimulationEvaluationParams(**{**base, **overrides})


def _no_preprocessing(*args, **kwargs):
    raise AssertionError("the setting must be rejected before any preprocessing")


@pytest.fixture(scope="module")
def sparse_pls(dataset):
    return integrate_semisynthetic_dataset(dataset, _params(BLOCK))


def _features(dataset, weighting=None) -> tuple[pd.DataFrame, pd.Series]:
    preprocessor = fit_omics_preprocessor(dataset, TWO_BLOCKS, weighting=weighting)
    X = concatenate_blocks(preprocessor.transform_dataset(dataset)).reset_index(drop=True)
    return X, dataset.metadata["stage"].astype(str).reset_index(drop=True)


def _core(X, y, integration_params, **kwargs):
    return fit_pls_latent_space(
        X, y, integration_params, stage_col="stage", sample_ids=[str(i) for i in range(len(y))],
        layer_feature_counts=COUNTS, layers_metadata={}, original_scales=np.ones(X.shape[1]),
        methylation_units="mvalue", **kwargs,
    )


# --------------------------------------------------------------------------- #
# 3.1 The production core
# --------------------------------------------------------------------------- #


def test_sparse_metadata_records_setting_selection_and_support(sparse_pls) -> None:
    metadata = sparse_pls.metadata
    params = metadata["integration_params"]
    assert params["sparsity"] == {"scope": "block", "keep_fractions": list(DEFAULT_KEEP_FRACTIONS)}
    fraction = params["selected_keep_fraction"]
    assert fraction in DEFAULT_KEEP_FRACTIONS
    block = metadata["sparsity"]
    assert block["scope"] == "block" and block["keep_fractions"] == list(DEFAULT_KEEP_FRACTIONS)
    assert block["selected_keep_fraction"] == fraction
    expected = dict(zip(TWO_BLOCKS, SparsityDesign("block", (fraction,), (367, 131)).counts(fraction)))
    assert block["support"] == [expected] * metadata["selected_lv"]
    model = sparse_pls.artifacts.model
    assert isinstance(model, SparsePLSDA)
    assert model.n_components == metadata["selected_lv"]
    # The latent matrix is the sparse model's scores.
    np.testing.assert_array_equal(sparse_pls.matrix.to_numpy(), model.x_scores_)


def test_cv_selects_over_the_requested_ladder(dataset) -> None:
    latent = integrate_semisynthetic_dataset(dataset, _params({"scope": "joint", "keep_fractions": [0.1, 0.05]}))
    block = latent.metadata["sparsity"]
    assert block["keep_fractions"] == [0.1, 0.05]
    fraction = block["selected_keep_fraction"]
    assert fraction in (0.1, 0.05)
    total = round(fraction * 498)
    assert all(sum(row.values()) == total for row in block["support"])


def test_forced_rank_fits_at_a_one_element_ladder(dataset) -> None:
    latent = integrate_semisynthetic_dataset(
        dataset, _params({"scope": "block", "keep_fractions": [0.2]},
                         integration_params={**FAST_CV, "forced_components": 3})
    )
    assert latent.metadata["component_selection"] == "forced"
    assert latent.metadata["selected_lv"] == 3
    assert latent.metadata["sparsity"]["selected_keep_fraction"] == 0.2
    assert latent.metadata["sparsity"]["support"] == [{"methylation": 73, "expression": 26}] * 3


@pytest.mark.parametrize("entry", ["evaluate", "integrate", "core"])
def test_forced_rank_with_a_ladder_is_rejected_naming_both(dataset, monkeypatch, entry) -> None:
    integration_params = {**FAST_CV, "forced_components": 3, "sparsity": BLOCK}
    match = r"forced_components'\]=3 .*ladder \[1\.0, 0\.5, 0\.2, 0\.1, 0\.05\]"
    if entry == "core":
        X, y = _features(dataset)
        monkeypatch.setattr(evaluation_module, "fit_plsda_model", _no_preprocessing)
        with pytest.raises(SimulationEvaluationError, match=match):
            _core(X, y, integration_params)
        return
    monkeypatch.setattr(evaluation_module, "fit_omics_preprocessor", _no_preprocessing)
    call = evaluate_semisynthetic_trajectory if entry == "evaluate" else integrate_semisynthetic_dataset
    with pytest.raises(SimulationEvaluationError, match=match):
        call(dataset, _params(integration_params=integration_params))


def test_membership_must_cover_the_columns(dataset) -> None:
    X, y = _features(dataset)
    with pytest.raises(SimulationEvaluationError, match=r"cover 400 features, but X has 498"):
        fit_pls_latent_space(
            X, y, {"sparsity": BLOCK, "n_repeats": 1}, stage_col="stage", sample_ids=[str(i) for i in range(len(y))],
            layer_feature_counts={"methylation": 367, "expression": 33}, layers_metadata={},
            original_scales=np.ones(X.shape[1]), methylation_units="mvalue",
        )


def _stage_nuisance(dataset, seed: int = 5) -> NuisanceDesign:
    stages = dataset.metadata["stage"].astype(str).to_numpy()
    levels = sorted(set(stages))
    conditioning = np.column_stack([(stages == level).astype(float) for level in levels[1:]])
    z = np.random.default_rng(seed).normal(size=stages.size)
    return NuisanceDesign(conditioning, z[:, None], tuple(f"stage[{v}]" for v in levels[1:]), ("z",), "stage")


def test_composes_with_the_latent_adjustment(dataset) -> None:
    X, y = _features(dataset)
    nuisance = _stage_nuisance(dataset)
    latent = _core(X, y, {"sparsity": BLOCK, "n_repeats": 1, "max_components": 6}, nuisance=nuisance)
    model = latent.artifacts.model
    assert isinstance(model, SparsePLSDA) and not model.scale
    np.testing.assert_array_equal(np.asarray(model._x_std), np.ones(X.shape[1]))
    # The basis is fitted on the residuals; the raw features are projected onto it.
    np.testing.assert_allclose(latent.matrix.to_numpy(), model.transform(X.to_numpy(dtype=float)), atol=1e-12)
    assert not np.allclose(latent.matrix.to_numpy(), model.x_scores_)
    assert latent.metadata["integration_params"]["latent_adjustment"] == "stage"
    assert latent.metadata["sparsity"]["scope"] == "block"


def test_composes_with_orthonormal_coordinates(dataset, sparse_pls) -> None:
    orthonormal = integrate_semisynthetic_dataset(dataset, _params(BLOCK, latent_coordinates="orthonormal"))
    model = orthonormal.artifacts.model
    basis = orthonormal.linear_map.basis
    rotations = np.asarray(model.x_rotations_)
    np.testing.assert_allclose(basis @ (basis.T @ rotations), rotations, atol=1e-10)
    inputs = (orthonormal.artifacts.features.to_numpy(dtype=float) - model._x_mean) / model._x_std
    np.testing.assert_allclose(orthonormal.matrix.to_numpy(), inputs @ basis, atol=1e-10)
    # Same fit as the score-coordinate space, and the derived builder agrees with it.
    np.testing.assert_array_equal(model.x_weights_, sparse_pls.artifacts.model.x_weights_)
    derived = pls_orthonormal_from_pls(dataset, sparse_pls, _params(BLOCK, latent_coordinates="orthonormal"))
    np.testing.assert_array_equal(derived.matrix.to_numpy(), orthonormal.matrix.to_numpy())


def test_composes_with_block_weighting(dataset) -> None:
    latent = integrate_semisynthetic_dataset(
        dataset, _params(BLOCK, integration_params={**FAST_CV, "block_weighting": "equal_variance"})
    )
    model = latent.artifacts.model
    assert isinstance(model, SparsePLSDA) and not model.scale
    assert latent.metadata["integration_params"]["block_weighting"] == "equal_variance"
    assert latent.metadata["integration_params"]["sparsity"] == {
        "scope": "block", "keep_fractions": list(DEFAULT_KEEP_FRACTIONS)
    }


def test_full_budget_ladder_is_the_production_space(dataset) -> None:
    production = integrate_semisynthetic_dataset(dataset, _params())
    full = integrate_semisynthetic_dataset(dataset, _params({"scope": "joint", "keep_fractions": [1.0]}))
    np.testing.assert_array_equal(full.matrix.to_numpy(), production.matrix.to_numpy())
    assert full.metadata["selected_lv"] == production.metadata["selected_lv"]
    assert full.metadata["cv_mean_auroc"] == production.metadata["cv_mean_auroc"]


# --------------------------------------------------------------------------- #
# 3.2 The harness
# --------------------------------------------------------------------------- #


@pytest.mark.parametrize(
    ("method", "integration_params"),
    [
        ("pca", {"n_components": 3}),
        ("concat", {}),
        ("snf", {"spectral_components": 4}),
        ("pls_within_group", {}),
    ],
)
@pytest.mark.parametrize("entry", ["evaluate", "integrate"])
def test_other_methods_reject_the_setting_before_preprocessing(dataset, monkeypatch, method, integration_params,
                                                                entry) -> None:
    monkeypatch.setattr(evaluation_module, "fit_omics_preprocessor", _no_preprocessing)
    params = _params(BLOCK, integration_method=method,
                     integration_params={"layers": list(TWO_BLOCKS), **integration_params})
    call = evaluate_semisynthetic_trajectory if entry == "evaluate" else integrate_semisynthetic_dataset
    with pytest.raises(SimulationEvaluationError, match=rf"'pls' or 'pls_residual_pca'; got it with '{method}'"):
        call(dataset, params)


@pytest.mark.parametrize(
    ("value", "match"),
    [
        ("block", r"must be a mapping"),
        ({"scope": "block", "ladder": [0.5]}, r"unknown key\(s\) \['ladder'\]"),
        ({"keep_fractions": [0.5]}, r"unknown scope None"),
        ({"scope": "blocks"}, r"unknown scope 'blocks'"),
        ({"scope": "joint", "keep_fractions": [1.0, 0.0]}, r"keep fraction 0\.0 is outside \(0, 1\]"),
        ({"scope": "joint", "keep_fractions": [0.5, 0.5]}, r"repeat a fraction"),
        ({"scope": "joint", "keep_fractions": []}, r"at least one fraction"),
    ],
)
@pytest.mark.parametrize("entry", ["evaluate", "integrate"])
def test_malformed_setting_is_rejected_before_preprocessing(dataset, monkeypatch, value, match, entry) -> None:
    monkeypatch.setattr(evaluation_module, "fit_omics_preprocessor", _no_preprocessing)
    call = evaluate_semisynthetic_trajectory if entry == "evaluate" else integrate_semisynthetic_dataset
    with pytest.raises(SimulationEvaluationError, match=match):
        call(dataset, _params(value))


def test_residual_pca_accepts_the_setting(dataset, sparse_pls) -> None:
    params = _params(BLOCK, integration_method="pls_residual_pca")
    built = integrate_semisynthetic_dataset(dataset, params)
    derived = pls_residual_pca_from_pls(dataset, sparse_pls, params)
    np.testing.assert_array_equal(built.matrix.to_numpy(), derived.matrix.to_numpy())
    assert derived.metadata["sparsity"] == sparse_pls.metadata["sparsity"]
    assert derived.metadata["integration_params"]["sparsity"]["scope"] == "block"
    pls_columns = sparse_pls.matrix.shape[1]
    np.testing.assert_array_equal(derived.matrix.to_numpy()[:, :pls_columns], sparse_pls.matrix.to_numpy())


def test_prebuilt_and_derived_spaces_refuse_another_sparsity(dataset, sparse_pls) -> None:
    plain = integrate_semisynthetic_dataset(dataset, _params())
    with pytest.raises(SimulationEvaluationError, match=r"sparsity=\{'scope': 'block'.*not the requested None"):
        evaluate_semisynthetic_trajectory(dataset, _params(), latent=sparse_pls)
    with pytest.raises(SimulationEvaluationError, match=r"sparsity=None, not the requested \{'scope': 'joint'"):
        evaluate_semisynthetic_trajectory(dataset, _params({"scope": "joint"}), latent=plain)
    with pytest.raises(SimulationEvaluationError, match=r"fitted with sparsity=None, but this evaluation"):
        pls_residual_pca_from_pls(dataset, plain, _params(BLOCK, integration_method="pls_residual_pca"))
    with pytest.raises(SimulationEvaluationError, match=r"fitted with sparsity=\{'scope': 'block'"):
        pls_orthonormal_from_pls(dataset, sparse_pls, _params(latent_coordinates="orthonormal"))
    reused = evaluate_semisynthetic_trajectory(dataset, _params(BLOCK), latent=sparse_pls)
    built = evaluate_semisynthetic_trajectory(dataset, _params(BLOCK))
    assert reused.pair_statistics == built.pair_statistics


def test_span_capture_and_attribution_run_on_a_sparse_fit(dataset) -> None:
    params = _params(
        BLOCK, span_capture=True, permutations=9, seed=3,
        attribution=AttributionDiagnosticSettings(enabled=True, bootstrap_replicates=3, top_k=10),
    )
    result = evaluate_semisynthetic_trajectory(dataset, params)
    capture = result.realized_geometry["span_capture"]
    assert 0.0 <= capture["pooled_stage"] <= 1.0
    assert capture["rank"] == result.latent_matrix_metadata["selected_lv"]
    attribution = result.attribution_diagnostics
    assert attribution["status"] == "computed"
    assert attribution["model"]["n_features"] == 367 + 131
    assert set(result.p_values) >= {"delta", "angle"}
    assert result.latent_matrix_metadata["sparsity"]["scope"] == "block"


def test_a_set_key_is_part_of_the_cell_identity() -> None:
    generator = SemiSyntheticTrajectoryParams(seed=1, n_samples=60)
    plain = make_simulation_cell(phase="power_primary", generator_params=generator,
                                 evaluation_params=_params(integration_params={"layers": list(TWO_BLOCKS)}))
    sparse = make_simulation_cell(phase="power_primary", generator_params=generator,
                                  evaluation_params=_params(BLOCK, integration_params={"layers": list(TWO_BLOCKS)}))
    assert parameter_signature(plain) != parameter_signature(sparse)
    assert plain.cell_id != sparse.cell_id
    assert "sparsity" not in json.dumps(parameter_signature(plain), default=str)


@pytest.mark.parametrize("name", ["pls", "pls_forced", "pca", "concat"])
def test_absent_sparsity_is_the_prechange_evaluation(name) -> None:
    summary = json.loads(json.dumps(evaluation_summary(name)))
    _assert_matches(summary, PRECHANGE_FIXTURE["evaluations"][name])
    assert "sparsity" not in json.dumps(summary)
