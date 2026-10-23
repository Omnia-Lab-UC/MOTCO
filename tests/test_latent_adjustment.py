"""Covariate-residualized latent space: the shared PLS core and the conditioning design."""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
import pytest

from motco.analysis import _covariate_adjusted_design
from motco.simulations.evaluation import (
    SimulationEvaluationError,
    build_simulation_trajectory_design,
    fit_pls_latent_space,
)
from motco.stats.design import build_conditioning_matrix, build_covariate_matrix
from motco.stats.permutation import RRPP
from motco.stats.pls import (
    NuisanceDesign,
    fit_nuisance_coefficients,
    fit_plsda_model,
    nuisance_partial_r2,
    plsda_doubleCV,
    residualize,
)
from motco.stats.trajectory import estimate_difference
from tests.test_phase6_byte_identity import _assert_close

#: ``fit_pls_latent_space`` output on a fixed input, captured before the
#: ``nuisance`` argument existed (``uv run python -m tests.test_latent_adjustment``).
CORE_IDENTITY = Path(__file__).parent / "data" / "pls_latent_space_identity.json"
STAGE_LABELS = ("0", "1", "2")


def _core_input(seed: int = 20261020, n_per_cell: int = 8, n_features: int = 14) -> dict[str, Any]:
    rng = np.random.default_rng(seed)
    groups = np.repeat([0, 1], 3 * n_per_cell)
    stages = np.tile(np.repeat([0, 1, 2], n_per_cell), 2)
    n = groups.size
    loading = rng.normal(size=n_features)
    X = rng.normal(size=(n, n_features)) + np.outer(stages * (1.0 + 0.3 * groups), loading)
    X = (X - X.mean(axis=0)) / X.std(axis=0, ddof=0)
    return {
        "X": pd.DataFrame(X, columns=[f"rna__g{j}" for j in range(n_features)]),
        "groups": groups,
        "stages": stages,
        "y": pd.Series([STAGE_LABELS[s] for s in stages]),
        "sample_ids": [f"s{i}" for i in range(n)],
        "rng": rng,
    }


def _fit_core(
    inputs: dict[str, Any],
    integration_params: dict[str, Any] | None = None,
    *,
    min_components: int = 1,
    **kwargs: Any,
):
    """The core on the module's fixed input.

    ``min_components`` defaults to 1 here (the rule the committed identity
    fixture and the nuisance tests were captured under); the core's own default
    is 2, exercised by the rank-floor tests below.
    """

    params = {"cv1_splits": 2, "cv2_splits": 3, "n_repeats": 2, "max_components": 5, "random_state": 3}
    params.update(integration_params or {})
    X = inputs["X"]
    return fit_pls_latent_space(
        X,
        inputs["y"],
        params,
        stage_col="stage",
        sample_ids=inputs["sample_ids"],
        layer_feature_counts={"rna": int(X.shape[1])},
        layers_metadata={},
        original_scales=np.ones(X.shape[1]),
        methylation_units="not_applicable",
        min_components=min_components,
        **kwargs,
    )


def _core_identity_payload() -> dict[str, Any]:
    latent = _fit_core(_core_input())
    return json.loads(
        json.dumps(
            {
                "matrix": latent.matrix.to_numpy().tolist(),
                "metadata": latent.metadata,
                "x_rotations": np.asarray(latent.artifacts.model.x_rotations_).tolist(),
            },
            default=list,
        )
    )


def test_core_without_nuisance_reproduces_stored_output() -> None:
    expected = json.loads(CORE_IDENTITY.read_text(encoding="utf-8"))
    _assert_close(_core_identity_payload(), expected, "fit_pls_latent_space")
    assert expected["metadata"]["selected_lv"] == 1  # the fixture is a rank-1 selection under floor 1
    assert "min_components" not in expected["metadata"]["integration_params"]


# ── Rank floor in the core ──────────────────────────────────────────────────


def _core_default_floor(inputs: dict[str, Any], integration_params: dict[str, Any] | None = None, **kwargs: Any):
    """The core at its own default floor (no ``min_components`` passed)."""

    params = {"cv1_splits": 2, "cv2_splits": 3, "n_repeats": 2, "max_components": 5, "random_state": 3}
    params.update(integration_params or {})
    X = inputs["X"]
    return fit_pls_latent_space(
        X,
        inputs["y"],
        params,
        stage_col="stage",
        sample_ids=inputs["sample_ids"],
        layer_feature_counts={"rna": int(X.shape[1])},
        layers_metadata={},
        original_scales=np.ones(X.shape[1]),
        methylation_units="not_applicable",
        **kwargs,
    )


def test_default_floor_lifts_the_identity_input_to_rank_two() -> None:
    latent = _core_default_floor(_core_input())
    assert latent.metadata["selected_lv"] >= 2
    assert latent.matrix.shape[1] >= 2
    assert latent.metadata["integration_params"]["min_components"] == 2
    assert latent.metadata["component_selection"] == "cv"
    explicit = _fit_core(_core_input(), min_components=2)
    assert explicit.metadata == latent.metadata
    np.testing.assert_array_equal(explicit.matrix.to_numpy(), latent.matrix.to_numpy())


def test_named_floor_is_honoured() -> None:
    latent = _fit_core(_core_input(), min_components=3)
    assert latent.metadata["selected_lv"] >= 3
    assert latent.metadata["integration_params"]["min_components"] == 3


def test_maximum_at_the_floor_is_raised_to_one_candidate() -> None:
    latent = _core_default_floor(_core_input(), {"max_components": 2})
    assert latent.metadata["integration_params"]["max_components"] == 3
    assert latent.metadata["selected_lv"] == 2


def test_sample_too_small_for_the_floor_is_rejected() -> None:
    inputs = _core_input(n_features=2)
    with pytest.raises(SimulationEvaluationError, match=r"min_components=2 needs at least 3 features.*2 feature"):
        _core_default_floor(inputs)
    # The unfloored rule still fits a two-feature input (candidates {1}).
    assert _fit_core(inputs).metadata["selected_lv"] == 1


def test_forced_components_ignore_the_floor() -> None:
    latent = _fit_core(_core_input(), {"forced_components": 2}, min_components=3)
    assert latent.metadata["component_selection"] == "forced"
    assert latent.metadata["selected_lv"] == 2
    assert "min_components" not in latent.metadata["integration_params"]


@pytest.mark.parametrize("floor", [0, 2.5, -1, "2"])
def test_invalid_floor_is_rejected(floor) -> None:
    with pytest.raises(SimulationEvaluationError, match="min_components must be an integer of at least 1"):
        _fit_core(_core_input(), min_components=floor)



# ── Residualized basis in the core ──────────────────────────────────────────


def _with_nuisance(kind: str = "stage", seed: int = 20261020) -> dict[str, Any]:
    """The core input plus a stage-correlated covariate ``z`` that loads on most features."""

    inputs = _core_input(seed=seed)
    rng = inputs["rng"]
    n, p = inputs["X"].shape
    z = inputs["stages"] + 0.5 * inputs["groups"] + rng.normal(scale=0.7, size=n)
    loading = np.where(np.arange(p) < p - 2, 1.5, 0.0) * rng.choice([-1.0, 1.0], size=p)
    X = inputs["X"].to_numpy() + np.outer(z - z.mean(), loading)
    X = (X - X.mean(axis=0)) / X.std(axis=0, ddof=0)
    inputs["X"] = pd.DataFrame(X, columns=inputs["X"].columns)
    frame = pd.DataFrame({"group": inputs["groups"].astype(str), "stage": inputs["stages"].astype(str)})
    conditioning = build_conditioning_matrix(frame, "group", "stage", kind)
    inputs["z"] = z
    inputs["frame"] = frame
    inputs["nuisance"] = NuisanceDesign(conditioning.matrix, z[:, None], conditioning.column_names, ("z",), kind)
    return inputs


@pytest.fixture(scope="module")
def adjusted_core():
    inputs = _with_nuisance()
    return inputs, _fit_core(inputs, nuisance=inputs["nuisance"])


def test_nuisance_scores_project_the_raw_features(adjusted_core) -> None:
    inputs, latent = adjusted_core
    model = latent.artifacts.model
    X = inputs["X"].to_numpy()
    assert model.scale is False
    np.testing.assert_allclose(latent.matrix.to_numpy(), model.transform(X), atol=1e-12)
    difference = latent.matrix.to_numpy() - np.asarray(model.x_scores_)
    assert np.linalg.norm(difference) > 1e-3
    span = np.column_stack([np.ones(X.shape[0]), inputs["nuisance"].covariates])
    leftover = difference - span @ np.linalg.lstsq(span, difference, rcond=None)[0]
    np.testing.assert_allclose(leftover, 0.0, atol=1e-10)


def test_nuisance_basis_is_fitted_on_residualized_features_at_the_pooled_scale(adjusted_core) -> None:
    inputs, latent = adjusted_core
    X, nuisance = inputs["X"], inputs["nuisance"]
    coefficients, _ = fit_nuisance_coefficients(X, nuisance)
    X_res = residualize(X.to_numpy(), nuisance.covariates, coefficients)
    expected = fit_plsda_model(X_res, inputs["y"], latent.metadata["selected_lv"], scale=False)
    np.testing.assert_allclose(latent.artifacts.model.x_rotations_, expected.x_rotations_, atol=1e-10)
    np.testing.assert_allclose(latent.artifacts.model.x_scores_, expected.x_scores_, atol=1e-10)
    # A feature that was mostly covariate keeps a residual variance below one; it is not rescaled.
    r2 = nuisance_partial_r2(X, nuisance)
    heavy = int(np.argmax(r2))
    assert r2[heavy] > 0.5
    assert X_res[:, heavy].var() < 0.5


def test_nuisance_metadata_and_summary_match_a_direct_computation(adjusted_core) -> None:
    inputs, latent = adjusted_core
    X, nuisance = inputs["X"], inputs["nuisance"]
    assert latent.metadata["integration_params"]["latent_adjustment"] == "stage"
    assert latent.metadata["component_selection"] == "cv"
    summary = latent.artifacts.latent_adjustment
    n = X.shape[0]
    assert summary.conditioning == "stage"
    assert summary.design_columns == ("stage[1]", "stage[2]", "z")
    assert summary.n_samples == n
    assert summary.residual_df == n - 1 - 2 - 1
    assert summary.feature_scale == "pooled_standardized"
    cv = plsda_doubleCV(
        X, inputs["y"], cv1_splits=2, cv2_splits=3, n_repeats=2, max_components=5, random_state=3,
        progress=False, scale=False, nuisance=nuisance,
    )
    assert summary.cv_min_training_df == cv["nuisance_min_training_df"]
    assert latent.metadata["selected_lv"] == int(pd.Series(cv["table"]["LV"]).mode().min())
    r2 = nuisance_partial_r2(X, nuisance)
    expected = {f"q{round(100 * q):02d}": float(np.quantile(r2, q)) for q in (0.05, 0.25, 0.5, 0.75, 0.95)}
    expected.update({"mean": float(r2.mean()), "fraction_above_0.5": float(np.mean(r2 > 0.5))})
    assert summary.feature_r2 == pytest.approx(expected)
    assert list(summary.feature_r2) == ["q05", "q25", "q50", "q75", "q95", "mean", "fraction_above_0.5"]


def test_without_nuisance_there_is_no_adjustment_entry() -> None:
    latent = _fit_core(_with_nuisance())
    assert "latent_adjustment" not in latent.metadata["integration_params"]
    assert latent.artifacts.latent_adjustment is None
    assert latent.artifacts.model.scale is True


def test_forced_components_are_honoured_with_nuisance() -> None:
    inputs = _with_nuisance(kind="group_stage")
    latent = _fit_core(inputs, {"forced_components": 3}, nuisance=inputs["nuisance"])
    assert latent.metadata["component_selection"] == "forced"
    assert latent.metadata["selected_lv"] == 3
    assert latent.metadata["integration_params"]["latent_adjustment"] == "group_stage"
    summary = latent.artifacts.latent_adjustment
    assert summary.cv_min_training_df is None
    assert summary.residual_df == inputs["X"].shape[0] - 1 - 5 - 1
    np.testing.assert_allclose(latent.matrix.to_numpy(), latent.artifacts.model.transform(inputs["X"].to_numpy()))


def test_group_stage_cross_validation_runs_on_group_sorted_rows() -> None:
    inputs = _with_nuisance(kind="group_stage")
    assert np.all(np.diff(inputs["groups"]) >= 0)  # rows are sorted by group
    latent = _fit_core(inputs, nuisance=inputs["nuisance"])
    assert latent.metadata["component_selection"] == "cv"
    assert latent.metadata["integration_params"]["latent_adjustment"] == "group_stage"
    assert latent.artifacts.latent_adjustment.cv_min_training_df > 0

def test_rank_deficient_full_sample_design_is_rejected() -> None:
    inputs = _with_nuisance()
    stage_index = inputs["stages"].astype(float)[:, None]
    nuisance = NuisanceDesign(
        inputs["nuisance"].conditioning, stage_index, inputs["nuisance"].conditioning_names, ("stage_index",), "stage"
    )
    with pytest.raises(SimulationEvaluationError, match="latent adjustment failed.*stage_index\\] is not of full rank"):
        _fit_core(inputs, nuisance=nuisance)


# ── The equivalence the raw-feature projection relies on ────────────────────


@pytest.mark.parametrize("estimator", ["pc1", "net_displacement"])
def test_raw_and_residualized_projections_give_identical_inference(adjusted_core, estimator: str) -> None:
    inputs, latent = adjusted_core
    model = latent.artifacts.model
    Y_raw = latent.matrix.to_numpy()
    Y_res = np.asarray(model.x_scores_)
    design = build_simulation_trajectory_design(inputs["frame"], group_col="group", stage_col="stage")
    covariates = build_covariate_matrix(pd.DataFrame({"z": inputs["z"]}), ["z"])
    full, reduced, ls_means = _covariate_adjusted_design(
        design.model_full, design.model_reduced, design.ls_means, covariates
    )
    observed = [
        estimate_difference(Y, full, ls_means, design.contrast, orientation_estimator=estimator)
        for Y in (Y_raw, Y_res)
    ]
    for a, b in zip(*observed, strict=True):
        np.testing.assert_allclose(a, b, atol=1e-9)
    nulls = [
        RRPP(Y, full, reduced, ls_means, design.contrast, permutations=29, seed=11, orientation_estimator=estimator,
             progress=False)
        for Y in (Y_raw, Y_res)
    ]
    for a, b in zip(*nulls, strict=True):
        np.testing.assert_allclose(np.asarray(a, dtype=float), np.asarray(b, dtype=float), atol=1e-9)
    assert not np.allclose(Y_raw, Y_res)

if __name__ == "__main__":
    CORE_IDENTITY.write_text(json.dumps(_core_identity_payload(), indent=2) + "\n", encoding="utf-8")
    print(f"wrote {CORE_IDENTITY}")
