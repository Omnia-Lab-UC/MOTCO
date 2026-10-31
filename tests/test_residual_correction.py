"""Residual-corrected RRPP: the multipliers, the primitive's null and the real-data analysis."""

from __future__ import annotations

import math
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from motco.analysis import (
    TrajectoryAnalysisError,
    TrajectoryAnalysisParams,
    _covariate_adjusted_design,
    analyze_trajectories,
)
from motco.simulations.evaluation import (
    _empirical_p_value,
    _extract_null_distributions,
    build_simulation_trajectory_design,
)
from motco.stats import permutation
from motco.stats.design import build_covariate_matrix
from motco.stats.permutation import RESIDUAL_CORRECTIONS, RRPP, _RRPPWorker, residual_correction_scale
from motco.stats.trajectory import estimate_betas, estimate_difference
from tests.test_analysis import SMALL_CV, make_dataset, with_covariates

ESTIMATOR = "net_displacement"


def _diagnosis():
    sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
    try:
        return __import__("case_study_diagnosis")
    finally:
        sys.path.pop(0)


def _hat_diagonal(X: np.ndarray) -> np.ndarray:
    return np.diag(X @ np.linalg.inv(X.T @ X) @ X.T)


def _unbalanced_inputs(seed: int = 0):
    """Unequal group × stage cells, one numeric covariate, a 3-column latent with group structure."""

    rng = np.random.default_rng(seed)
    sizes = [[6, 5, 9], [4, 7, 5]]
    groups = np.concatenate([[g] * sizes[g][s] for g in range(2) for s in range(3)])
    stages = np.concatenate([[s] * sizes[g][s] for g in range(2) for s in range(3)])
    design = build_simulation_trajectory_design(
        pd.DataFrame({"group": groups.astype(str), "stage": stages.astype(str)}), group_col="group", stage_col="stage"
    )
    z = 0.7 * groups * stages + rng.normal(size=groups.size)
    covariates = build_covariate_matrix(pd.DataFrame({"z": z}), ["z"])
    full, reduced, ls_means = _covariate_adjusted_design(
        design.model_full, design.model_reduced, design.ls_means, covariates
    )
    Y = rng.normal(size=(3, 3))[stages] + np.outer(z, rng.normal(size=3)) + rng.normal(size=(groups.size, 3))
    return Y, np.asarray(full, dtype=float), np.asarray(reduced, dtype=float), np.asarray(ls_means), design.contrast


def _reduced_fit(Y: np.ndarray, reduced: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """The reduced fitted values and residuals, computed exactly as RRPP computes them."""

    frame = pd.DataFrame(Y)
    y_hat = np.matmul(reduced, estimate_betas(reduced, frame))
    return np.asarray(y_hat, dtype=float), np.asarray(frame - y_hat, dtype=float)


def _manual_null(Y, full, reduced, ls_means, contrast, scale, *, permutations, seed, center=False):
    y_hat, y_res = _reduced_fit(Y, reduced)
    residuals = y_res * scale[:, None]
    if center:
        residuals = residuals - residuals.mean(axis=0)
    rng = np.random.default_rng(seed)
    draws = []
    for _ in range(permutations):
        perm = rng.permutation(Y.shape[0])
        draws.append(
            estimate_difference(y_hat + residuals[perm], full, ls_means, contrast, orientation_estimator=ESTIMATOR)
        )
    return tuple(list(column) for column in zip(*draws))


def _assert_nulls_equal(a, b, *, atol: float = 0.0) -> None:
    assert len(a) == len(b)
    for left, right in zip(a, b):
        np.testing.assert_allclose(np.asarray(left), np.asarray(right), rtol=0.0, atol=atol)


def _rrpp(Y, full, reduced, ls_means, contrast, **kwargs):
    kwargs.setdefault("permutations", 25)
    kwargs.setdefault("seed", 7)
    return RRPP(Y, full, reduced, ls_means, contrast, n_jobs=kwargs.pop("n_jobs", 1), progress=False,
                orientation_estimator=ESTIMATOR, **kwargs)


# ── Multipliers ──────────────────────────────────────────────────────────────


def test_corrections_are_named() -> None:
    assert RESIDUAL_CORRECTIONS == ("df", "leverage")


def test_df_multiplier_at_the_case_study_model() -> None:
    rng = np.random.default_rng(1)
    X = np.column_stack([np.ones(80), rng.normal(size=(80, 9))])
    scale = residual_correction_scale(X, "df")
    assert scale.shape == (80,)
    np.testing.assert_allclose(scale, math.sqrt(80 / 70), rtol=0, atol=1e-15)


def test_df_multiplier_uses_the_rank() -> None:
    rng = np.random.default_rng(2)
    X = np.column_stack([np.ones(30), rng.normal(size=(30, 4))])
    duplicated = np.column_stack([X, X[:, 2]])
    np.testing.assert_allclose(residual_correction_scale(duplicated, "df"), math.sqrt(30 / 25))


def test_leverage_multipliers_match_the_hat_matrix() -> None:
    _, _, reduced, _, _ = _unbalanced_inputs()
    h = _hat_diagonal(reduced)
    scale = residual_correction_scale(reduced, "leverage")
    np.testing.assert_allclose(scale, 1.0 / np.sqrt(1.0 - h), rtol=1e-12)
    assert np.sum(1.0 - scale**-2.0) == pytest.approx(np.linalg.matrix_rank(reduced), abs=1e-10)
    assert np.ptp(scale) > 0.01  # unequal cells: the leverages differ


def test_balanced_one_factor_leverage_equals_kappa() -> None:
    factor = np.repeat([0, 1, 2], 7)
    X = np.column_stack([np.ones(21), factor == 1, factor == 2]).astype(float)
    np.testing.assert_allclose(residual_correction_scale(X, "leverage"), residual_correction_scale(X, "df"))
    np.testing.assert_allclose(residual_correction_scale(X, "df"), math.sqrt(21 / 18))


def test_helper_is_importable_from_the_permutation_module() -> None:
    from motco.stats.permutation import residual_correction_scale as helper

    assert helper(np.ones((5, 1)), "df").shape == (5,)


def test_unknown_correction_is_rejected() -> None:
    with pytest.raises(ValueError, match=r"'hc3'.*\['df', 'leverage'\]"):
        residual_correction_scale(np.ones((5, 1)), "hc3")


def test_saturated_model_is_rejected() -> None:
    with pytest.raises(ValueError, match="n = 5 rows and rank 5"):
        residual_correction_scale(np.eye(5), "df")


def test_leverage_one_names_the_row() -> None:
    X = np.column_stack([np.ones(6), np.eye(6)[:, 3]])
    with pytest.raises(ValueError, match="row 3 of the reduced model has leverage 1"):
        residual_correction_scale(X, "leverage")


# ── RRPP ─────────────────────────────────────────────────────────────────────


def test_absent_correction_equals_explicit_none() -> None:
    Y, full, reduced, ls_means, contrast = _unbalanced_inputs()
    default = RRPP(Y, full, reduced, ls_means, contrast, permutations=25, n_jobs=1, progress=False, seed=7,
                   orientation_estimator=ESTIMATOR)
    explicit = _rrpp(Y, full, reduced, ls_means, contrast, residual_correction=None)
    _assert_nulls_equal(default, explicit)


def test_df_equals_the_diagnosis_construction() -> None:
    Y, full, reduced, ls_means, contrast = _unbalanced_inputs()
    kappa = residual_correction_scale(reduced, "df")[0]
    corrected = _rrpp(Y, full, reduced, ls_means, contrast, residual_correction="df")
    rescaled = _diagnosis().rescaled_residuals(Y, reduced, kappa)
    y_hat, y_res = _reduced_fit(Y, reduced)
    np.testing.assert_allclose(rescaled, y_hat + kappa * y_res, atol=1e-12)
    _assert_nulls_equal(corrected, _rrpp(rescaled, full, reduced, ls_means, contrast), atol=1e-10)
    uncorrected = _rrpp(Y, full, reduced, ls_means, contrast)
    assert not np.allclose(np.asarray(corrected[0]), np.asarray(uncorrected[0]))


def test_leverage_permutes_the_rescaled_residuals_with_the_same_indices() -> None:
    Y, full, reduced, ls_means, contrast = _unbalanced_inputs()
    scale = residual_correction_scale(reduced, "leverage")
    corrected = _rrpp(Y, full, reduced, ls_means, contrast, residual_correction="leverage")
    manual = _manual_null(Y, full, reduced, ls_means, contrast, scale, permutations=25, seed=7)
    _assert_nulls_equal(corrected, manual, atol=1e-12)
    # Unit multipliers through the same loop reproduce the uncorrected null: the indices are RRPP's.
    plain = _manual_null(Y, full, reduced, ls_means, contrast, np.ones(Y.shape[0]), permutations=25, seed=7)
    _assert_nulls_equal(_rrpp(Y, full, reduced, ls_means, contrast), plain, atol=1e-12)


@pytest.mark.parametrize("correction", [None, "df", "leverage"])
def test_statistics_do_not_depend_on_the_residual_mean(correction) -> None:
    Y, full, reduced, ls_means, contrast = _unbalanced_inputs()
    scale = np.ones(Y.shape[0]) if correction is None else residual_correction_scale(reduced, correction)
    raw = _manual_null(Y, full, reduced, ls_means, contrast, scale, permutations=25, seed=7)
    centred = _manual_null(Y, full, reduced, ls_means, contrast, scale, permutations=25, seed=7, center=True)
    if correction == "leverage":
        y_res = _reduced_fit(Y, reduced)[1] * scale[:, None]
        assert np.abs(y_res.mean(axis=0)).max() > 1e-6  # the leverage-rescaled residuals are not centred
    _assert_nulls_equal(raw, centred, atol=1e-10)


@pytest.mark.parametrize("correction", ["df", "leverage"])
def test_parallel_path_draws_from_the_corrected_residuals(correction) -> None:
    Y, full, reduced, ls_means, contrast = _unbalanced_inputs()
    first = _rrpp(Y, full, reduced, ls_means, contrast, permutations=9, n_jobs=2, residual_correction=correction)
    second = _rrpp(Y, full, reduced, ls_means, contrast, permutations=9, n_jobs=2, residual_correction=correction)
    _assert_nulls_equal(first, second)
    # The workers' seeds and counts as RRPP spawns them, run in process on the corrected residuals.
    seeds = [int(s.generate_state(1)[0]) for s in np.random.SeedSequence(7).spawn(2)]
    y_hat, y_res = _reduced_fit(Y, reduced)
    worker = _RRPPWorker(y_hat, y_res * residual_correction_scale(reduced, correction)[:, None], full, ls_means,
                         contrast, orientation_estimator=ESTIMATOR)
    parts = [worker(count, seed) for count, seed in zip([5, 4], seeds)]
    expected = tuple([draw for part in parts for draw in part[index]] for index in range(3))
    _assert_nulls_equal(first, expected, atol=1e-12)


@pytest.mark.parametrize(
    ("correction", "inputs", "message"),
    [
        ("hc3", "unbalanced", "'hc3'"),
        ("df", "saturated", "n = 6 rows and rank 6"),
        ("leverage", "singleton", "row 0 of the reduced model"),
    ],
)
def test_undefined_correction_raises_before_any_permutation(monkeypatch, correction, inputs, message) -> None:
    def forbidden(*args, **kwargs):
        raise AssertionError("no permutation may run")

    if inputs == "unbalanced":
        Y, full, reduced, ls_means, contrast = _unbalanced_inputs()
    else:
        Y = np.random.default_rng(3).normal(size=(6, 2))
        full = np.eye(6)
        ls_means = np.eye(6)
        contrast = [[0, 1, 2], [3, 4, 5]]
        reduced = np.eye(6) if inputs == "saturated" else np.column_stack([np.ones(6), np.eye(6)[:, 0]])
    monkeypatch.setattr(permutation, "estimate_betas", forbidden)
    monkeypatch.setattr(permutation, "estimate_difference", forbidden)
    with pytest.raises(ValueError, match=message):
        _rrpp(Y, full, reduced, ls_means, contrast, residual_correction=correction)


# ── Analysis ─────────────────────────────────────────────────────────────────


def _covariate_dataset():
    metadata, blocks = with_covariates()
    return make_dataset(metadata, blocks, covariate_cols=["age", "pmi"])


def _params(**kwargs) -> TrajectoryAnalysisParams:
    kwargs.setdefault("permutations", 49)
    kwargs.setdefault("seed", 3)
    return TrajectoryAnalysisParams(**SMALL_CV, **kwargs)


@pytest.fixture(scope="module")
def covariate_runs():
    dataset = _covariate_dataset()
    return dataset, {
        correction: analyze_trajectories(dataset, _params(residual_correction=correction))
        for correction in (None, "df", "leverage")
    }


def test_correction_moves_only_the_null(covariate_runs) -> None:
    _, runs = covariate_runs
    base = runs[None]
    for correction in ("df", "leverage"):
        corrected = runs[correction]
        pd.testing.assert_frame_equal(base.latent_scores, corrected.latent_scores, check_exact=True)
        pd.testing.assert_frame_equal(base.stage_means, corrected.stage_means, check_exact=True)
        assert base.record["latent"] == corrected.record["latent"]
        assert base.record["statistics"] == corrected.record["statistics"]
        assert base.record["covariates"] == corrected.record["covariates"]
        assert base.record["null_summary"] != corrected.record["null_summary"]


def test_analysis_uses_the_primitives_leverage_null(covariate_runs) -> None:
    dataset, runs = covariate_runs
    result = runs["leverage"]
    working = pd.DataFrame(
        {
            "group": dataset.metadata[dataset.group_col].map({g: str(i) for i, g in enumerate(dataset.groups)}),
            "stage": dataset.metadata[dataset.stage_col].map({s: str(i) for i, s in enumerate(dataset.stages)}),
        }
    )
    design = build_simulation_trajectory_design(working, group_col="group", stage_col="stage")
    full, reduced, ls_means = _covariate_adjusted_design(
        design.model_full, design.model_reduced, design.ls_means, dataset.covariate_design
    )
    nulls = RRPP(result.latent_scores, full, reduced, ls_means, design.contrast, permutations=49, n_jobs=1,
                 progress=False, seed=3, orientation_estimator=ESTIMATOR, residual_correction="leverage")
    null_values = _extract_null_distributions(*nulls, shape_available=True)
    expected = {
        statistic: _empirical_p_value(null_values[statistic], observed)
        for statistic, observed in result.record["statistics"].items()
    }
    assert result.record["p_values"] == expected


def test_df_record_states_kind_n_rank_and_kappa(covariate_runs) -> None:
    _, runs = covariate_runs
    record = runs["df"].record
    block = record["measurement"]["residual_correction"]
    assert block == {"kind": "df", "n": 48, "reduced_rank": 6, "kappa": pytest.approx(math.sqrt(48 / 42))}
    assert record["params"]["residual_correction"] == "df"
    keys = list(record["measurement"])
    assert keys.index("residual_correction") == keys.index("n_jobs") + 1


def test_leverage_record_states_the_range(covariate_runs) -> None:
    _, runs = covariate_runs
    block = runs["leverage"].record["measurement"]["residual_correction"]
    assert set(block) == {"kind", "n", "reduced_rank", "leverage_min", "leverage_max", "scale_min", "scale_max"}
    assert block["kind"] == "leverage" and block["n"] == 48 and block["reduced_rank"] == 6
    assert 0 < block["leverage_min"] < block["leverage_max"] < 1
    assert block["scale_min"] == pytest.approx(1 / math.sqrt(1 - block["leverage_min"]))
    assert block["scale_max"] == pytest.approx(1 / math.sqrt(1 - block["leverage_max"]))


def test_no_correction_carries_no_key(covariate_runs) -> None:
    _, runs = covariate_runs
    record = runs[None].record
    assert "residual_correction" not in record["params"]
    assert "residual_correction" not in record["measurement"]


def test_unknown_correction_is_rejected_before_fitting(monkeypatch) -> None:
    from motco import analysis

    def forbidden(*args, **kwargs):
        raise AssertionError("the latent space must not be fitted")

    monkeypatch.setattr(analysis, "fit_pls_latent_space", forbidden)
    with pytest.raises(TrajectoryAnalysisError, match=r"'uniform'.*\['df', 'leverage'\]"):
        analyze_trajectories(_covariate_dataset(), _params(residual_correction="uniform"))


def test_correction_without_covariates_uses_the_group_stage_model() -> None:
    record = analyze_trajectories(make_dataset(), _params(permutations=19, residual_correction="df")).record
    block = record["measurement"]["residual_correction"]
    assert block["reduced_rank"] == 4 and block["n"] == 48
    assert block["kappa"] == pytest.approx(math.sqrt(48 / 44))


def test_correction_without_permutations_is_not_reported_as_applied() -> None:
    record = analyze_trajectories(
        _covariate_dataset(), _params(permutations=0, forced_components=2, residual_correction="df")
    ).record
    assert record["params"]["residual_correction"] == "df"
    assert "residual_correction" not in record["measurement"]
    assert record["p_values"] == {}
