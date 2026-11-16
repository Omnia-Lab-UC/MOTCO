"""Stratified RRPP: the draws stay in their strata, identity while absent, the diagnosis's construction, parallel."""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest
from tqdm import tqdm

from motco.simulations.evaluation import build_simulation_trajectory_design
from motco.stats import permutation
from motco.stats.permutation import RRPP, residual_correction_scale, strata_groups
from motco.stats.trajectory import estimate_betas, estimate_difference, pooled_relative_eigengap
from tests.test_residual_correction import _reduced_fit, _unbalanced_inputs

ESTIMATOR = "net_displacement"
PERMUTATIONS = 40


def _size_excess_check():
    sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
    try:
        return __import__("predictive_size_excess_check")
    finally:
        sys.path.pop(0)


def _stage_design(sizes=((6, 5, 9), (4, 7, 5)), seed: int = 1):
    """Unequal group × stage cells and a 3-column response with stage structure; stage positions as strata."""

    rng = np.random.default_rng(seed)
    groups = np.concatenate([[g] * sizes[g][s] for g in range(2) for s in range(3)])
    stages = np.concatenate([[s] * sizes[g][s] for g in range(2) for s in range(3)])
    design = build_simulation_trajectory_design(
        pd.DataFrame({"group": groups.astype(str), "stage": stages.astype(str)}), group_col="group", stage_col="stage"
    )
    Y = 2.0 * rng.normal(size=(3, 3))[stages] + rng.normal(size=(groups.size, 3))
    return Y, design, groups, stages


def _rrpp(Y, design, **kwargs):
    return RRPP(
        Y, design.model_full, design.model_reduced, design.ls_means, design.contrast,
        permutations=kwargs.pop("permutations", PERMUTATIONS), progress=False, seed=kwargs.pop("seed", 11),
        orientation_estimator=kwargs.pop("orientation_estimator", ESTIMATOR), **kwargs,
    )


def _bytes(null) -> list[bytes]:
    return [np.asarray(part, dtype=float).tobytes() for part in null]


def _captured_draws(monkeypatch, Y, design, strata, *, n_jobs=None):
    """Every permuted response RRPP evaluates, mapped back to the residual index each row received."""

    captured: list[np.ndarray] = []

    def capture(y_random, *args, **kwargs):
        captured.append(np.asarray(y_random, dtype=float).copy())
        return np.zeros(1), np.zeros(1), np.zeros(1)

    monkeypatch.setattr(permutation, "estimate_difference", capture)
    _rrpp(Y, design, strata=strata, n_jobs=n_jobs)
    y_hat, y_res = _reduced_fit(Y, np.asarray(design.model_reduced, dtype=float))
    draws = []
    for y_random in captured:
        permuted = y_random - y_hat
        distances = ((permuted[:, None, :] - y_res[None, :, :]) ** 2).sum(axis=2)
        draws.append(distances.argmin(axis=1))
    return draws


# ── Validation and the strata ────────────────────────────────────────────────


def test_strata_groups_are_in_ascending_label_order() -> None:
    groups = strata_groups([2, 0, 1, 0, 2, 2], 6)
    assert [rows.tolist() for rows in groups] == [[1, 3], [2], [0, 4, 5]]


def test_wrong_length_is_rejected_before_any_draw(monkeypatch) -> None:
    Y, design, _, stages = _stage_design()
    monkeypatch.setattr(permutation, "estimate_betas", lambda *a, **k: pytest.fail("fitted before validating"))
    with pytest.raises(ValueError, match=rf"strata has {Y.shape[0] - 1} labels but Y has {Y.shape[0]} rows"):
        _rrpp(Y, design, strata=stages[:-1])


@pytest.mark.parametrize("missing", [None, np.nan])
def test_missing_label_is_rejected_naming_the_row(monkeypatch, missing) -> None:
    Y, design, _, stages = _stage_design()
    labels = stages.astype(object)
    labels[4] = missing
    monkeypatch.setattr(permutation, "estimate_betas", lambda *a, **k: pytest.fail("fitted before validating"))
    with pytest.raises(ValueError, match="row 4 is missing"):
        _rrpp(Y, design, strata=labels)


def test_two_dimensional_strata_are_rejected() -> None:
    Y, design, _, stages = _stage_design()
    with pytest.raises(ValueError, match="one-dimensional"):
        _rrpp(Y, design, strata=np.column_stack([stages, stages]))


# ── The draws ────────────────────────────────────────────────────────────────


def test_every_draw_stays_in_its_stratum_and_is_a_permutation(monkeypatch) -> None:
    Y, design, _, stages = _stage_design()
    draws = _captured_draws(monkeypatch, Y, design, stages)
    assert len(draws) == PERMUTATIONS
    moved = 0
    for idx in draws:
        assert sorted(idx.tolist()) == list(range(Y.shape[0]))
        np.testing.assert_array_equal(stages[idx], stages)
        moved += int(np.any(idx != np.arange(Y.shape[0])))
    assert moved == PERMUTATIONS


def test_a_one_row_stratum_keeps_its_residual(monkeypatch) -> None:
    Y, design, _, stages = _stage_design()
    strata = stages.copy()
    strata[7] = 9
    for idx in _captured_draws(monkeypatch, Y, design, strata):
        assert idx[7] == 7
        np.testing.assert_array_equal(strata[idx], strata)


def test_stratum_means_stay_at_the_reduced_fit(monkeypatch) -> None:
    Y, design, _, stages = _stage_design()
    captured: list[np.ndarray] = []

    def capture(y_random, *args, **kwargs):
        captured.append(np.asarray(y_random, dtype=float).copy())
        return np.zeros(1), np.zeros(1), np.zeros(1)

    monkeypatch.setattr(permutation, "estimate_difference", capture)
    _rrpp(Y, design, strata=stages)
    y_hat, _ = _reduced_fit(Y, np.asarray(design.model_reduced, dtype=float))
    for y_random in captured:
        for stage in range(3):
            rows = stages == stage
            np.testing.assert_allclose(y_random[rows].mean(axis=0), y_hat[rows].mean(axis=0), atol=1e-12)


def _manual_draws(strata: np.ndarray, visit: list[int], seed: int, count: int) -> list[np.ndarray]:
    rng = np.random.default_rng(seed)
    draws = []
    for _ in range(count):
        draw = np.empty(strata.size, dtype=np.intp)
        for label in visit:
            rows = np.flatnonzero(strata == label)
            draw[rows] = rows[rng.permutation(rows.size)]
        draws.append(draw)
    return draws


def test_strata_are_visited_in_label_order_not_row_order(monkeypatch) -> None:
    """Rows reordered so stage 2 comes first: the draws still visit stage 0, then 1, then 2."""

    Y, _, groups, stages = _stage_design()
    order = np.argsort(-stages, kind="stable")
    Y, groups, stages = Y[order], groups[order], stages[order]
    design = build_simulation_trajectory_design(
        pd.DataFrame({"group": groups.astype(str), "stage": stages.astype(str)}), group_col="group", stage_col="stage"
    )
    assert stages[0] == 2
    draws = _captured_draws(monkeypatch, Y, design, stages)
    ascending = _manual_draws(stages, [0, 1, 2], 11, PERMUTATIONS)
    first_seen = _manual_draws(stages, [2, 1, 0], 11, PERMUTATIONS)
    for ours, expected in zip(draws, ascending):
        np.testing.assert_array_equal(ours, expected)
    assert any(not np.array_equal(ours, other) for ours, other in zip(draws, first_seen))


def test_one_stratum_equals_the_unrestricted_draws() -> None:
    Y, design, _, _ = _stage_design()
    assert _bytes(_rrpp(Y, design, strata=np.zeros(Y.shape[0], dtype=int))) == _bytes(_rrpp(Y, design))


def test_strata_change_the_null() -> None:
    Y, design, _, stages = _stage_design()
    assert _bytes(_rrpp(Y, design, strata=stages)) != _bytes(_rrpp(Y, design))


# ── Identity while absent: a frozen copy of the 0.24.0 serial loop ───────────


def _frozen_0240_rrpp_serial(Y, model_full, model_reduced, LS_means, contrast, permutations, seed,
                             return_eigengaps=False, orientation_estimator="pc1", residual_correction=None):
    correction_scale = (
        None if residual_correction is None
        else residual_correction_scale(np.asarray(model_reduced, dtype=float), residual_correction)
    )
    Y = pd.DataFrame(Y)
    betas_red = estimate_betas(model_reduced, Y)
    y_hat = np.matmul(model_reduced, betas_red)
    y_hat.index = Y.index
    y_res = Y - y_hat
    y_hat_np = np.asarray(y_hat)
    y_res_np = np.asarray(y_res)
    if correction_scale is not None:
        y_res_np = y_res_np * correction_scale[:, None]
    deltas, angles, shapes, eigengaps = [], [], [], []
    n = y_res_np.shape[0]
    rng = np.random.default_rng(seed)
    for _ in tqdm(range(permutations), desc="RRPP", unit="perm", disable=True):
        idx = rng.permutation(n)
        y_random = y_hat_np + y_res_np[idx, :]
        if return_eigengaps:
            d, a, s, spectra = estimate_difference(
                y_random, model_full, LS_means, contrast, return_spectra=True,
                orientation_estimator=orientation_estimator,
            )
            eigengaps.append(pooled_relative_eigengap(spectra))
        else:
            d, a, s = estimate_difference(
                y_random, model_full, LS_means, contrast, orientation_estimator=orientation_estimator
            )
        deltas.append(d)
        angles.append(a)
        shapes.append(s)
    if return_eigengaps:
        return deltas, angles, shapes, eigengaps
    return deltas, angles, shapes


@pytest.mark.parametrize("estimator", ["pc1", "net_displacement"])
@pytest.mark.parametrize("eigengaps", [False, True])
@pytest.mark.parametrize("correction", [None, "df"])
def test_absent_strata_are_byte_identical_to_0240(estimator, eigengaps, correction) -> None:
    Y, design, _, _ = _stage_design(seed=5)
    frozen = _frozen_0240_rrpp_serial(
        Y, design.model_full, design.model_reduced, design.ls_means, design.contrast, PERMUTATIONS, 23,
        return_eigengaps=eigengaps, orientation_estimator=estimator, residual_correction=correction,
    )
    kwargs = dict(seed=23, orientation_estimator=estimator, return_eigengaps=eigengaps, residual_correction=correction)
    implicit = _rrpp(Y, design, **dict(kwargs))
    explicit = _rrpp(Y, design, strata=None, **dict(kwargs))
    assert _bytes(implicit) == _bytes(frozen)
    assert _bytes(explicit) == _bytes(frozen)


# ── The diagnosis's construction and the residual correction ─────────────────


def _unbalanced_stages() -> np.ndarray:
    """The stage of each row of ``_unbalanced_inputs`` (cells group-major, stage-minor, sizes 6/5/9 and 4/7/5)."""

    return np.concatenate([[s] * size for sizes in ((6, 5, 9), (4, 7, 5)) for s, size in enumerate(sizes)])


@pytest.mark.parametrize("correction", [None, "df"])
def test_strata_equal_the_diagnosis_within_stage_loop(correction) -> None:
    check = _size_excess_check()
    Y, full, reduced, ls_means, contrast = _unbalanced_inputs(seed=2)
    stages = _unbalanced_stages()
    ours = RRPP(
        Y, full, reduced, ls_means, contrast, permutations=PERMUTATIONS, progress=False, seed=31,
        orientation_estimator=ESTIMATOR, residual_correction=correction, strata=stages,
    )
    indices = check.within_stage_indices(stages.astype(str), ["0", "1", "2"], PERMUTATIONS, 31)
    theirs = check.restricted_null(
        Y, full, reduced, ls_means, contrast, indices, orientation_estimator=ESTIMATOR, residual_correction=correction
    )
    assert _bytes(ours) == _bytes(theirs)


def test_df_under_strata_equals_kappa_scaled_residuals() -> None:
    Y, full, reduced, ls_means, contrast = _unbalanced_inputs(seed=4)
    stages = _unbalanced_stages()
    kappa = residual_correction_scale(reduced, "df")
    y_hat, y_res = _reduced_fit(Y, reduced)
    corrected = RRPP(
        Y, full, reduced, ls_means, contrast, permutations=PERMUTATIONS, progress=False, seed=41,
        orientation_estimator=ESTIMATOR, residual_correction="df", strata=stages,
    )
    manual = RRPP(
        y_hat + kappa[:, None] * y_res, full, reduced, ls_means, contrast, permutations=PERMUTATIONS,
        progress=False, seed=41, orientation_estimator=ESTIMATOR, strata=stages,
    )
    for ours, theirs in zip(corrected, manual):
        np.testing.assert_allclose(np.asarray(ours, dtype=float), np.asarray(theirs, dtype=float), rtol=1e-9,
                                   atol=1e-9)


def test_eigengap_return_does_not_change_the_stratified_draws() -> None:
    Y, design, _, stages = _stage_design()
    plain = _rrpp(Y, design, strata=stages)
    with_gaps = _rrpp(Y, design, strata=stages, return_eigengaps=True)
    assert len(with_gaps) == 4
    assert _bytes(with_gaps[:3]) == _bytes(plain)


# ── The parallel path ────────────────────────────────────────────────────────


def test_parallel_stratified_runs_are_reproducible() -> None:
    Y, design, _, stages = _stage_design()
    first = _rrpp(Y, design, strata=stages, n_jobs=2, seed=7)
    second = _rrpp(Y, design, strata=stages, n_jobs=2, seed=7)
    assert _bytes(first) == _bytes(second)
    assert _bytes(first) != _bytes(_rrpp(Y, design, n_jobs=2, seed=7))


@pytest.mark.parametrize("n_jobs", [None, 2])
def test_residuals_constant_within_strata_never_move(n_jobs) -> None:
    """Rows identical within each group × stage cell: a within-cell draw is the unpermuted response."""

    _, design, groups, stages = _stage_design()
    cells = groups * 3 + stages
    rng = np.random.default_rng(13)
    Y = rng.normal(size=(6, 4))[cells]
    intercept_only = np.ones((Y.shape[0], 1))
    null = RRPP(
        Y, design.model_full, intercept_only, design.ls_means, design.contrast, permutations=12, progress=False,
        seed=3, n_jobs=n_jobs, orientation_estimator=ESTIMATOR, strata=cells,
    )
    y_hat, y_res = _reduced_fit(Y, intercept_only)
    observed = estimate_difference(
        y_hat + y_res, design.model_full, design.ls_means, design.contrast, orientation_estimator=ESTIMATOR
    )
    for statistic, draws in zip(observed, null):
        assert len(draws) == 12
        for draw in draws:
            assert np.asarray(draw, dtype=float).tobytes() == np.asarray(statistic, dtype=float).tobytes()
    unrestricted = RRPP(
        Y, design.model_full, intercept_only, design.ls_means, design.contrast, permutations=12, progress=False,
        seed=3, n_jobs=n_jobs, orientation_estimator=ESTIMATOR,
    )
    assert any(
        np.asarray(draw).tobytes() != np.asarray(observed[0]).tobytes() for draw in unrestricted[0]
    )
