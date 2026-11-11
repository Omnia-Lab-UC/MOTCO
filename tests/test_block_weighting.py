"""The label-blind block weighting of the shared preprocessing and the harness (``block_weighting``)."""

from __future__ import annotations

import json
import math
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

import motco.simulations.evaluation as evaluation_module
import motco.simulations.preprocessing as preprocessing_module
from motco.simulations.diagnostics import POPULATION_LATENT_CHECKPOINT, geometry_from_means
from motco.simulations.evaluation import (
    SimulationEvaluationError,
    SimulationEvaluationParams,
    evaluate_semisynthetic_trajectory,
    fit_pls_latent_space,
    integrate_semisynthetic_dataset,
    pls_orthonormal_from_pls,
    pls_residual_pca_from_pls,
)
from motco.simulations.grid import make_simulation_cell, parameter_signature
from motco.simulations.preprocessing import (
    BLOCK_WEIGHTINGS,
    concatenate_blocks,
    fit_block_scalers,
    fit_omics_preprocessor,
    transform_blocks,
)
from motco.simulations.reference import load_reference
from motco.simulations.semisynthetic import SemiSyntheticTrajectoryParams, generate_semisynthetic_trajectory
from motco.stats.pls import fit_plsda_model, plsda_doubleCV
from tests.orthonormal_coordinates_fixture_capture import evaluation_summary
from tests.test_orthonormal_coordinates import _assert_matches
from tests.test_span_capture import _frozen_0150_span_capture, _population_rows

TWO_BLOCKS = ("methylation", "expression")
FAST_CV = {"layers": list(TWO_BLOCKS), "n_repeats": 1, "max_components": 6}


PRECHANGE_FIXTURE = json.loads(
    (Path(__file__).parent / "data" / "orthonormal_coordinates_prechange_fixture.json").read_text()
)


def _two_blocks() -> dict[str, pd.DataFrame]:
    """Beta-valued methylation (37 features + one constant) and correlated expression (13 features)."""

    rng = np.random.default_rng(20261110)
    n = 40
    index = [f"s{i}" for i in range(n)]
    beta = pd.DataFrame(rng.uniform(0.05, 0.95, size=(n, 37)), index=index, columns=[f"cg{j}" for j in range(37)])
    beta["cg_const"] = 0.5
    expression = pd.DataFrame(
        rng.normal(size=(n, 13)) @ rng.normal(size=(13, 13)), index=index, columns=[f"g{j}" for j in range(13)]
    )
    return {"methylation": beta, "expression": expression}


def _sized_blocks(sizes=(367, 131), n=60, seed=4) -> dict[str, pd.DataFrame]:
    """Blocks of the simulations' feature counts, the second with a dominant direction."""

    rng = np.random.default_rng(seed)
    index = [f"s{i}" for i in range(n)]
    first = rng.normal(size=(n, sizes[0]))
    second = rng.normal(size=(n, sizes[1])) + 3.0 * rng.normal(size=(n, 1))
    return {
        "a": pd.DataFrame(first, index=index, columns=[f"a{j}" for j in range(sizes[0])]),
        "b": pd.DataFrame(second, index=index, columns=[f"b{j}" for j in range(sizes[1])]),
    }


TRANSFORMS = {"methylation": "beta_to_mvalue"}


@pytest.fixture(scope="module")
def reference():
    return load_reference()


@pytest.fixture(scope="module")
def dataset(reference):
    return generate_semisynthetic_trajectory(
        SemiSyntheticTrajectoryParams(seed=8, n_samples=120, p_dmp=0.1, trajectory_mode="orientation",
                                      group_effect_size=0.5),
        reference=reference,
    )


def _params(**overrides) -> SimulationEvaluationParams:
    integration_params = overrides.pop("integration_params", FAST_CV)
    weighting = overrides.pop("weighting", None)
    if weighting is not None:
        integration_params = {**integration_params, "block_weighting": weighting}
    base = {"integration_method": "pls", "integration_params": integration_params, "permutations": 0}
    return SimulationEvaluationParams(**{**base, **overrides})


def _no_preprocessing(*args, **kwargs):
    raise AssertionError("the setting must be rejected before any preprocessing")


# --------------------------------------------------------------------------- #
# 1.1 The weights
# --------------------------------------------------------------------------- #


def test_kinds() -> None:
    assert BLOCK_WEIGHTINGS == ("equal_variance", "mfa")


def test_equal_variance_gives_unit_total_variance() -> None:
    blocks = _sized_blocks()
    scalers = fit_block_scalers(blocks, {}, weighting="equal_variance")
    transformed = transform_blocks(blocks, scalers, {})
    for name, size in (("a", 367), ("b", 131)):
        assert scalers[name].weight == 1.0 / math.sqrt(size)
        assert transformed[name].to_numpy().var(axis=0).sum() == pytest.approx(1.0, rel=1e-12)


def test_equal_variance_counts_constant_features() -> None:
    blocks = _two_blocks()
    scalers = fit_block_scalers(blocks, TRANSFORMS, weighting="equal_variance")
    transformed = transform_blocks(blocks, scalers, TRANSFORMS)
    assert scalers["methylation"].weight == 1.0 / math.sqrt(38)
    # The constant feature counts against the block's 38 features but carries no variance.
    assert transformed["methylation"].to_numpy().var(axis=0).sum() == pytest.approx(37 / 38, rel=1e-12)


def test_mfa_gives_unit_variance_to_the_leading_component() -> None:
    blocks = _sized_blocks()
    scalers = fit_block_scalers(blocks, {}, weighting="mfa")
    transformed = transform_blocks(blocks, scalers, {})
    for name in ("a", "b"):
        values = transformed[name].to_numpy()
        leading = np.linalg.svd(values - values.mean(axis=0), compute_uv=False)[0]
        assert leading**2 / values.shape[0] == pytest.approx(1.0, rel=1e-12)
    # λ₁ ≤ p, so an mfa weight is never below equal variance's; relative to the isotropic block, the block
    # with a dominant direction loses weight.
    equal = fit_block_scalers(blocks, {}, weighting="equal_variance")
    assert all(scalers[name].weight >= equal[name].weight for name in ("a", "b"))
    assert scalers["b"].weight / scalers["a"].weight < equal["b"].weight / equal["a"].weight


def test_mfa_refuses_a_block_without_variance() -> None:
    blocks = {"flat": pd.DataFrame(np.ones((10, 4)), columns=list("wxyz"))}
    with pytest.raises(ValueError, match="Block 'flat' has no variance"):
        fit_block_scalers(blocks, {}, weighting="mfa")


@pytest.mark.parametrize("kind", ["pareto", "MFA", "", 1])
def test_unknown_kind_is_rejected_before_any_block_is_transformed(monkeypatch, kind) -> None:
    def no_transform(*args, **kwargs):
        raise AssertionError("no block may be transformed")

    monkeypatch.setattr(preprocessing_module, "apply_unit_transform", no_transform)
    allowed = r"allowed: \['equal_variance', 'mfa'\] or None"
    with pytest.raises(ValueError, match=rf"Unknown block weighting {kind!r}; {allowed}"):
        fit_block_scalers(_two_blocks(), TRANSFORMS, weighting=kind)


def _frozen_0200_standardize(blocks, transforms) -> dict[str, np.ndarray]:
    """``fit_block_scalers`` + ``transform_blocks`` as released in 0.20.0, kept to pin the unweighted path.

    A digest of the output is not portable (``log`` and BLAS differ in the last bit between machines),
    so the pre-change computation runs beside the current one on the same machine.
    """

    standardized = {}
    for name, matrix in blocks.items():
        values = preprocessing_module.apply_unit_transform(
            name, matrix.astype(float).to_numpy(dtype=float), transforms.get(name)
        )
        scale = values.std(axis=0)
        scale[scale < 1e-10] = 1.0
        standardized[name] = (values - values.mean(axis=0)) / scale
    return standardized


def test_unweighted_transform_is_byte_identical_to_the_prechange_output() -> None:
    blocks = _two_blocks()
    scalers = fit_block_scalers(blocks, TRANSFORMS)
    assert all(scaler.weight == 1.0 for scaler in scalers.values())
    out = transform_blocks(blocks, scalers, TRANSFORMS)
    expected = _frozen_0200_standardize(blocks, TRANSFORMS)
    for name in ("methylation", "expression"):
        assert out[name].to_numpy(dtype=float).tobytes() == expected[name].tobytes()


def test_weights_are_label_blind_row_order_invariant() -> None:
    blocks = _sized_blocks()
    order = np.random.default_rng(0).permutation(60)
    shuffled = {name: frame.iloc[order] for name, frame in blocks.items()}
    for kind in BLOCK_WEIGHTINGS:
        original = fit_block_scalers(blocks, {}, weighting=kind)
        permuted = fit_block_scalers(shuffled, {}, weighting=kind)
        for name in blocks:
            assert permuted[name].weight == pytest.approx(original[name].weight, rel=1e-12)


# --------------------------------------------------------------------------- #
# 1.2 The fitted preprocessor
# --------------------------------------------------------------------------- #


@pytest.mark.parametrize("kind", BLOCK_WEIGHTINGS)
def test_population_means_are_weighted_like_the_observed_rows(dataset, kind) -> None:
    plain = fit_omics_preprocessor(dataset, TWO_BLOCKS)
    weighted = fit_omics_preprocessor(dataset, TWO_BLOCKS, weighting=kind)
    assert plain.block_weighting is None and weighted.block_weighting == kind
    population = dataset.population_trajectories
    plain_population = plain.transform_population(population)
    weighted_population = weighted.transform_population(population)
    plain_observed = plain.transform_dataset(dataset)
    weighted_observed = weighted.transform_dataset(dataset)
    for layer in TWO_BLOCKS:
        weight = weighted.scalers[layer].weight
        assert weight != 1.0
        np.testing.assert_array_equal(weighted.scalers[layer].scale, plain.scalers[layer].scale)
        pairs = ((weighted_population, plain_population), (weighted_observed, plain_observed))
        for weighted_frame, plain_frame in pairs:
            np.testing.assert_array_equal(weighted_frame[layer].to_numpy(), plain_frame[layer].to_numpy() * weight)


def test_require_layers_is_unchanged_and_weighting_is_required(dataset) -> None:
    weighted = fit_omics_preprocessor(dataset, TWO_BLOCKS, weighting="mfa")
    assert weighted.layers == TWO_BLOCKS
    weighted.require_layers(TWO_BLOCKS)
    with pytest.raises(ValueError, match=r"fitted on layers \['methylation', 'expression'\]"):
        weighted.require_layers(("methylation", "expression", "proteomics"))
    weighted.require_block_weighting("mfa")
    with pytest.raises(ValueError, match="block_weighting='mfa', not the requested None"):
        weighted.require_block_weighting(None)
    with pytest.raises(ValueError, match="Unknown block weighting 'pareto'"):
        fit_omics_preprocessor(dataset, TWO_BLOCKS, weighting="pareto")


# --------------------------------------------------------------------------- #
# 2.1 The harness setting
# --------------------------------------------------------------------------- #


@pytest.mark.parametrize("entry", ["evaluate", "integrate"])
def test_unknown_kind_is_rejected_before_any_preprocessing(dataset, monkeypatch, entry) -> None:
    monkeypatch.setattr(evaluation_module, "fit_omics_preprocessor", _no_preprocessing)
    params = _params(weighting="pareto")
    call = evaluate_semisynthetic_trajectory if entry == "evaluate" else integrate_semisynthetic_dataset
    with pytest.raises(
        SimulationEvaluationError,
        match=r"block_weighting'\]: Unknown block weighting 'pareto'; allowed: \['equal_variance', 'mfa'\]",
    ):
        call(dataset, params)


@pytest.mark.parametrize(
    ("method", "integration_params", "match"),
    [
        ("snf", {"spectral_components": 4}, "cannot be used with integration_method='snf'"),
        ("concat", {"standardize": False}, r"unstandardized concat baseline \(standardize=False\)"),
    ],
)
@pytest.mark.parametrize("entry", ["evaluate", "integrate"])
def test_methods_the_weights_cannot_reach_are_refused(dataset, monkeypatch, method, integration_params, match,
                                                      entry) -> None:
    monkeypatch.setattr(evaluation_module, "fit_omics_preprocessor", _no_preprocessing)
    params = _params(
        integration_method=method,
        integration_params={"layers": list(TWO_BLOCKS), **integration_params},
        weighting="equal_variance",
    )
    call = evaluate_semisynthetic_trajectory if entry == "evaluate" else integrate_semisynthetic_dataset
    with pytest.raises(SimulationEvaluationError, match=match):
        call(dataset, params)


@pytest.mark.parametrize(("fitted", "requested"), [(None, "equal_variance"), ("mfa", "equal_variance"), ("mfa", None)])
def test_a_prebuilt_preprocessor_with_another_weighting_is_refused(dataset, fitted, requested) -> None:
    preprocessor = fit_omics_preprocessor(dataset, TWO_BLOCKS, weighting=fitted)
    with pytest.raises(
        SimulationEvaluationError, match=rf"block_weighting={fitted!r}, not the requested {requested!r}"
    ):
        integrate_semisynthetic_dataset(dataset, _params(weighting=requested), preprocessor=preprocessor)


@pytest.fixture(scope="module")
def weighted_pls(dataset):
    return integrate_semisynthetic_dataset(dataset, _params(weighting="equal_variance"))


@pytest.fixture(scope="module")
def plain_pls(dataset):
    return integrate_semisynthetic_dataset(dataset, _params())


def test_a_prebuilt_latent_with_another_weighting_is_refused(dataset, weighted_pls, plain_pls) -> None:
    with pytest.raises(SimulationEvaluationError, match="block_weighting='equal_variance', not the requested None"):
        evaluate_semisynthetic_trajectory(dataset, _params(), latent=weighted_pls)
    with pytest.raises(SimulationEvaluationError, match="block_weighting=None, not the requested 'mfa'"):
        evaluate_semisynthetic_trajectory(dataset, _params(weighting="mfa"), latent=plain_pls)
    # The matching prebuilt latent is accepted and measured as if built here.
    reused = evaluate_semisynthetic_trajectory(dataset, _params(weighting="equal_variance"), latent=weighted_pls)
    built = evaluate_semisynthetic_trajectory(dataset, _params(weighting="equal_variance"))
    assert reused.pair_statistics == built.pair_statistics


def test_derived_spaces_refuse_a_pls_fit_with_another_weighting(dataset, weighted_pls, plain_pls) -> None:
    with pytest.raises(SimulationEvaluationError, match="fitted with block_weighting=None, but this evaluation"):
        pls_orthonormal_from_pls(dataset, plain_pls, _params(weighting="mfa", latent_coordinates="orthonormal"))
    with pytest.raises(SimulationEvaluationError, match="fitted with block_weighting='equal_variance'"):
        pls_residual_pca_from_pls(dataset, weighted_pls, _params(integration_method="pls_residual_pca"))
    derived = pls_residual_pca_from_pls(
        dataset, weighted_pls, _params(integration_method="pls_residual_pca", weighting="equal_variance")
    )
    assert derived.metadata["integration_params"]["block_weighting"] == "equal_variance"
    orthonormal = pls_orthonormal_from_pls(
        dataset, weighted_pls, _params(weighting="equal_variance", latent_coordinates="orthonormal")
    )
    assert orthonormal.metadata["integration_params"]["block_weighting"] == "equal_variance"


def _per_block_total_variance(matrix: pd.DataFrame) -> dict[str, float]:
    variances = matrix.to_numpy().var(axis=0)
    blocks = np.array([column.split("__", 1)[0] for column in matrix.columns])
    return {layer: float(variances[blocks == layer].sum()) for layer in TWO_BLOCKS}


def test_concat_measures_the_weighted_blocks(dataset) -> None:
    params = _params(integration_method="concat", integration_params={"layers": list(TWO_BLOCKS)},
                     weighting="equal_variance")
    latent = integrate_semisynthetic_dataset(dataset, params)
    assert latent.metadata["integration_params"] == {
        "standardize": True, "layers": list(TWO_BLOCKS), "block_weighting": "equal_variance"
    }
    counts = latent.metadata["layer_feature_counts"]
    constant = {
        layer: int((getattr(dataset, layer).to_numpy().std(axis=0) < 1e-10).sum()) for layer in TWO_BLOCKS
    }
    for layer, total in _per_block_total_variance(latent.matrix).items():
        assert total == pytest.approx((counts[layer] - constant[layer]) / counts[layer], rel=1e-10)


def test_pca_measures_the_weighted_blocks(dataset) -> None:
    params = _params(integration_method="pca", integration_params={"layers": list(TWO_BLOCKS), "n_components": 3},
                     weighting="mfa")
    latent = integrate_semisynthetic_dataset(dataset, params)
    assert latent.metadata["integration_params"]["block_weighting"] == "mfa"
    preprocessor = fit_omics_preprocessor(dataset, TWO_BLOCKS, weighting="mfa")
    weighted = concatenate_blocks(preprocessor.transform_dataset(dataset))
    values = weighted.to_numpy()
    singular = np.linalg.svd(values - values.mean(axis=0), compute_uv=False)[:3]
    np.testing.assert_allclose(np.linalg.norm(latent.matrix.to_numpy(), axis=0), singular, rtol=1e-10)


def test_pls_measures_the_weighted_blocks(dataset, weighted_pls) -> None:
    assert weighted_pls.metadata["integration_params"]["block_weighting"] == "equal_variance"
    expected = concatenate_blocks(
        fit_omics_preprocessor(dataset, TWO_BLOCKS, weighting="equal_variance").transform_dataset(dataset)
    ).reset_index(drop=True)
    pd.testing.assert_frame_equal(weighted_pls.artifacts.features, expected)
    assert _per_block_total_variance(weighted_pls.artifacts.features) == pytest.approx(
        {layer: 1.0 for layer in TWO_BLOCKS}, rel=1e-10
    )


def test_a_set_key_is_part_of_the_cell_identity() -> None:
    generator = SemiSyntheticTrajectoryParams(seed=1, n_samples=60)
    plain = make_simulation_cell(phase="power_primary", generator_params=generator,
                                 evaluation_params=_params(integration_params={"layers": list(TWO_BLOCKS)}))
    weighted = make_simulation_cell(phase="power_primary", generator_params=generator,
                                    evaluation_params=_params(integration_params={"layers": list(TWO_BLOCKS)},
                                                              weighting="equal_variance"))
    assert parameter_signature(plain) != parameter_signature(weighted)
    assert plain.cell_id != weighted.cell_id


# --------------------------------------------------------------------------- #
# 2.2 The PLS metric
# --------------------------------------------------------------------------- #


def test_weighted_fit_has_unit_input_scaling_and_keeps_the_weights(dataset, weighted_pls) -> None:
    model = weighted_pls.artifacts.model
    np.testing.assert_array_equal(np.asarray(model._x_std), np.ones(model.x_weights_.shape[0]))
    stages = dataset.metadata["stage"].astype(str).reset_index(drop=True)
    rescaled = fit_plsda_model(weighted_pls.artifacts.features, stages, n_components=model.n_components, scale=True)
    first, other = model.x_weights_[:, 0], rescaled.x_weights_[:, 0]
    assert min(np.linalg.norm(first - other), np.linalg.norm(first + other)) > 0.05
    # scale=True would have undone the weights: its weights equal the unweighted production fit's.
    plain = fit_plsda_model(
        concatenate_blocks(fit_omics_preprocessor(dataset, TWO_BLOCKS).transform_dataset(dataset)).reset_index(
            drop=True
        ),
        stages,
        n_components=model.n_components,
        scale=True,
    )
    np.testing.assert_allclose(np.abs(rescaled.x_weights_[:, 0]), np.abs(plain.x_weights_[:, 0]), atol=1e-8)


def _core(X, y, integration_params):
    return fit_pls_latent_space(
        X, y, integration_params, stage_col="stage", sample_ids=[str(i) for i in range(len(y))],
        layer_feature_counts={}, layers_metadata={}, original_scales=np.ones(X.shape[1]), methylation_units="mvalue",
    )


def test_unit_weights_fit_production_at_a_fixed_rank(dataset) -> None:
    """At a fixed rank the unweighted metric (``scale=False``) is production's up to a uniform factor."""

    X = concatenate_blocks(fit_omics_preprocessor(dataset, TWO_BLOCKS).transform_dataset(dataset)).reset_index(
        drop=True
    )
    y = dataset.metadata["stage"].astype(str).reset_index(drop=True)
    production = _core(X, y, {"forced_components": 3})
    # The core with the key set on unweighted features: every weight is 1, only the metric switch differs.
    unit = _core(X, y, {"forced_components": 3, "block_weighting": "equal_variance"})
    n = X.shape[0]
    np.testing.assert_allclose(
        np.abs(unit.artifacts.model.x_weights_), np.abs(production.artifacts.model.x_weights_), atol=1e-10
    )
    ratio = np.abs(production.matrix.to_numpy()) / np.abs(unit.matrix.to_numpy())
    np.testing.assert_allclose(ratio, math.sqrt((n - 1) / n), rtol=1e-8)


def test_unit_weights_select_productions_rank_on_a_pinned_dataset(reference) -> None:
    """The CV's per-fold rescaling is the only other difference; on this dataset it changes nothing."""

    pinned = generate_semisynthetic_trajectory(
        SemiSyntheticTrajectoryParams(seed=31, n_samples=90, p_dmp=0.1, trajectory_mode="shape",
                                      group_effect_size=0.5),
        reference=reference,
    )
    X = concatenate_blocks(fit_omics_preprocessor(pinned, TWO_BLOCKS).transform_dataset(pinned)).reset_index(drop=True)
    y = pinned.metadata["stage"].astype(str).reset_index(drop=True)
    knobs = {"cv1_splits": 3, "cv2_splits": 4, "n_repeats": 5, "max_components": 20, "random_state": 1203,
             "min_components": 2, "progress": False}
    production = plsda_doubleCV(X, y, **knobs)
    unit = plsda_doubleCV(X, y, scale=False, **knobs)
    pd.testing.assert_frame_equal(unit["table"], production["table"])
    production_core = _core(X, y, {"max_components": 20})
    unit_core = _core(X, y, {"max_components": 20, "block_weighting": "equal_variance"})
    assert unit_core.metadata["selected_lv"] == production_core.metadata["selected_lv"]
    assert unit_core.metadata["cv_mean_auroc"] == production_core.metadata["cv_mean_auroc"]


def test_weighted_metadata_and_original_scales(dataset, weighted_pls, plain_pls) -> None:
    params = weighted_pls.metadata["integration_params"]
    assert params["block_weighting"] == "equal_variance"
    assert "block_weighting" not in plain_pls.metadata["integration_params"]
    preprocessor = fit_omics_preprocessor(dataset, TWO_BLOCKS, weighting="equal_variance")
    expected = np.concatenate(
        [preprocessor.scalers[layer].scale / preprocessor.scalers[layer].weight for layer in TWO_BLOCKS]
    )
    np.testing.assert_array_equal(weighted_pls.artifacts.original_scales, expected)
    plain_expected = np.concatenate([preprocessor.scalers[layer].scale for layer in TWO_BLOCKS])
    np.testing.assert_array_equal(plain_pls.artifacts.original_scales, plain_expected)
    # A forced rank records the weighting too.
    forced = integrate_semisynthetic_dataset(
        dataset, _params(integration_params={**FAST_CV, "forced_components": 2}, weighting="mfa")
    )
    assert forced.metadata["integration_params"]["block_weighting"] == "mfa"
    assert forced.metadata["component_selection"] == "forced"


# --------------------------------------------------------------------------- #
# 2.3 Span capture and identity
# --------------------------------------------------------------------------- #


@pytest.mark.parametrize("kind", BLOCK_WEIGHTINGS)
def test_span_capture_reads_the_weighted_inputs(dataset, kind) -> None:
    params = _params(weighting=kind, span_capture=True)
    result = evaluate_semisynthetic_trajectory(dataset, params)
    latent = integrate_semisynthetic_dataset(dataset, params)
    model = latent.artifacts.model
    rows = _population_rows(dataset, fit_omics_preprocessor(dataset, TWO_BLOCKS, weighting=kind))
    assert result.realized_geometry["span_capture"] == _frozen_0150_span_capture(model, rows).to_dict()
    expected_checkpoint = geometry_from_means(
        np.asarray(model.transform(rows), dtype=float), ["A", "B"], ["0", "1", "2"],
        orientation_estimator="net_displacement",
    ).to_dict()
    assert result.realized_geometry["checkpoints"][POPULATION_LATENT_CHECKPOINT]["joint"] == expected_checkpoint
    plain = evaluate_semisynthetic_trajectory(dataset, _params(span_capture=True))
    assert (
        result.realized_geometry["checkpoints"][POPULATION_LATENT_CHECKPOINT]
        != plain.realized_geometry["checkpoints"][POPULATION_LATENT_CHECKPOINT]
    )


@pytest.mark.parametrize("name", ["pls", "pls_forced", "pca", "concat"])
def test_absent_weighting_is_the_prechange_evaluation(name) -> None:
    summary = json.loads(json.dumps(evaluation_summary(name)))
    _assert_matches(summary, PRECHANGE_FIXTURE["evaluations"][name])
    assert "block_weighting" not in json.dumps(summary)
