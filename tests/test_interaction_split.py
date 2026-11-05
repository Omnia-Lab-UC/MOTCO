"""The split interaction test: definitions, decomposition, permutation path and validation."""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from motco.simulations.evaluation import build_simulation_trajectory_design
from motco.stats import interaction, permutation
from motco.stats.interaction import (
    INTERACTION_SPLIT_PARTS,
    InteractionSplitResult,
    interaction_split_test,
)
from motco.stats.permutation import RRPP

REL = 1e-10


def _design(sizes: list[list[int]]):
    groups = np.concatenate([[g] * sizes[g][s] for g in range(2) for s in range(len(sizes[0]))])
    stages = np.concatenate([[s] * sizes[g][s] for g in range(2) for s in range(len(sizes[0]))])
    design = build_simulation_trajectory_design(
        pd.DataFrame({"group": groups.astype(str), "stage": stages.astype(str)}), group_col="group", stage_col="stage"
    )
    return design, groups, stages


def _inputs(seed: int = 0, p: int = 12, r: int = 3, sizes: list[list[int]] | None = None, covariate: bool = True):
    """Unequal cells, an optional numeric covariate in both models, a random orthonormal basis."""

    rng = np.random.default_rng(seed)
    design, groups, stages = _design(sizes or [[6, 5, 9], [4, 7, 5]])
    full = np.asarray(design.model_full, dtype=float)
    reduced = np.asarray(design.model_reduced, dtype=float)
    if covariate:
        z = rng.normal(size=(groups.size, 1))
        full = np.hstack([full, z])
        reduced = np.hstack([reduced, z])
    Y = rng.normal(size=(3, p))[stages] + 0.5 * np.outer(groups * stages, rng.normal(size=p))
    Y = Y + rng.normal(size=Y.shape)
    basis, _ = np.linalg.qr(rng.normal(size=(p, r)))
    return Y, full, reduced, basis


def _hat(X: np.ndarray) -> np.ndarray:
    return X @ np.linalg.pinv(X)


def _direct(Y: np.ndarray, full: np.ndarray, reduced: np.ndarray, basis: np.ndarray) -> dict[str, tuple]:
    n = Y.shape[0]
    H_full, H_reduced = _hat(full), _hat(reduced)
    projector = basis @ basis.T
    df_int = np.linalg.matrix_rank(full) - np.linalg.matrix_rank(reduced)
    df_res = n - np.linalg.matrix_rank(full)
    out = {}
    for name, part in (("total", Y), ("span", Y @ projector), ("complement", Y - Y @ projector)):
        ss = float(np.sum(((H_full - H_reduced) @ part) ** 2))
        rss = float(np.sum(((np.eye(n) - H_full) @ part) ** 2))
        out[name] = (ss, rss, (ss / df_int) / (rss / df_res), int(df_int), int(df_res))
    return out


def _brute_force_null(Y, full, reduced, basis, permutations, seed, scale=None):
    """Refit Ŷ_R + s·E[π] per draw, with π from ``np.random.default_rng(seed)``."""

    n = Y.shape[0]
    H_reduced = _hat(reduced)
    fitted = H_reduced @ Y
    residuals = Y - fitted
    if scale is not None:
        residuals = residuals * scale[:, None]
    rng = np.random.default_rng(seed)
    nulls = {name: [] for name in INTERACTION_SPLIT_PARTS}
    for _ in range(permutations):
        idx = rng.permutation(n)
        for name, value in _direct(fitted + residuals[idx], full, reduced, basis).items():
            nulls[name].append(value)
    return {name: np.array(values) for name, values in nulls.items()}


# ── Observed quantities ──────────────────────────────────────────────────────


def test_definitions_match_direct_hat_matrices():
    Y, full, reduced, basis = _inputs()
    result = interaction_split_test(Y, full, reduced, basis)
    direct = _direct(Y, full, reduced, basis)
    for name in INTERACTION_SPLIT_PARTS:
        part = result.parts[name]
        ss, rss, f_value, df_int, df_res = direct[name]
        assert part.ss_interaction == pytest.approx(ss, rel=REL)
        assert part.rss == pytest.approx(rss, rel=REL)
        assert part.f == pytest.approx(f_value, rel=REL)
        assert (part.df_interaction, part.df_residual) == (df_int, df_res)
    assert (result.total.df_interaction, result.total.df_residual) == (2, 36 - 7)


def test_dimensions_are_p_r_and_p_minus_r():
    Y, full, reduced, basis = _inputs(p=12, r=3)
    result = interaction_split_test(Y, full, reduced, basis)
    assert [result.parts[name].dimension for name in INTERACTION_SPLIT_PARTS] == [12, 3, 9]


def test_span_and_complement_add_to_total():
    for seed in range(4):
        Y, full, reduced, basis = _inputs(seed=seed)
        result = interaction_split_test(Y, full, reduced, basis)
        for attribute in ("ss_interaction", "rss"):
            total = getattr(result.total, attribute)
            parts = getattr(result.span, attribute) + getattr(result.complement, attribute)
            assert parts == pytest.approx(total, rel=REL)
        assert 0.0 <= result.span_share <= 1.0
        assert result.span_share == pytest.approx(result.span.ss_interaction / result.total.ss_interaction, rel=REL)


def test_result_depends_on_the_basis_only_through_its_span():
    Y, full, reduced, basis = _inputs()
    rotation, _ = np.linalg.qr(np.random.default_rng(9).normal(size=(basis.shape[1], basis.shape[1])))
    first = interaction_split_test(Y, full, reduced, basis, permutations=30, seed=4)
    second = interaction_split_test(Y, full, reduced, basis @ rotation, permutations=30, seed=4)
    for name in INTERACTION_SPLIT_PARTS:
        a, b = first.parts[name], second.parts[name]
        for attribute in ("ss_interaction", "rss", "f"):
            assert getattr(b, attribute) == pytest.approx(getattr(a, attribute), rel=REL)
        assert (b.p_value_f, b.p_value_ss) == (a.p_value_f, a.p_value_ss)
    assert second.span_share == pytest.approx(first.span_share, rel=REL)


def test_noise_free_interaction_inside_the_span():
    rng = np.random.default_rng(2)
    design, groups, stages = _design([[5, 5, 5], [5, 5, 5]])
    full = np.asarray(design.model_full, dtype=float)
    reduced = np.asarray(design.model_reduced, dtype=float)
    basis, _ = np.linalg.qr(rng.normal(size=(10, 2)))
    coefficients = rng.normal(size=(full.shape[1], 2))
    Y = (full @ coefficients) @ basis.T
    result = interaction_split_test(Y, full, reduced, basis)
    assert result.span_share == pytest.approx(1.0, abs=1e-12)
    assert result.complement.ss_interaction == pytest.approx(0.0, abs=1e-20 * result.total.ss_interaction + 1e-20)
    assert result.total.ss_interaction > 0


def test_share_is_absent_without_interaction_energy():
    Y, full, reduced, basis = _inputs()
    result = interaction_split_test(np.zeros_like(Y), full, reduced, basis, permutations=5, seed=0)
    assert result.total.ss_interaction == 0.0
    assert result.span_share is None
    # 0 / 0: the pseudo-F is undefined and has no p-value; the SS p-value is still reported.
    assert np.isnan(result.total.f) and result.total.p_value_f is None
    assert result.total.p_value_ss == 1.0


# ── Permutation path ─────────────────────────────────────────────────────────


def test_null_matches_brute_force_refits():
    Y, full, reduced, basis = _inputs(seed=3)
    result = interaction_split_test(Y, full, reduced, basis, permutations=25, seed=11, return_null=True)
    brute = _brute_force_null(Y, full, reduced, basis, 25, 11)
    for name in INTERACTION_SPLIT_PARTS:
        part = result.parts[name]
        np.testing.assert_allclose(part.null_ss, brute[name][:, 0], rtol=REL)
        np.testing.assert_allclose(part.null_f, brute[name][:, 2], rtol=REL)
        expected_f = (1 + np.sum(part.null_f >= part.f)) / 26
        expected_ss = (1 + np.sum(part.null_ss >= part.ss_interaction)) / 26
        assert (part.p_value_f, part.p_value_ss) == (expected_f, expected_ss)


def test_permuted_outcomes_are_rrpps_serial_draws(monkeypatch):
    Y, full, reduced, basis = _inputs(seed=5, p=4, r=1)
    ls_means = np.zeros((6, full.shape[1]))
    ls_means[:, 0] = 1.0
    captured: list[np.ndarray] = []

    def capture(y_random, *args, **kwargs):
        captured.append(np.array(y_random, dtype=float))
        return np.zeros((2, 2)), np.zeros((2, 2)), np.zeros((2, 2))

    monkeypatch.setattr(permutation, "estimate_difference", capture)
    RRPP(Y, full, reduced, ls_means, [[0, 1, 2], [3, 4, 5]], permutations=15, progress=False, seed=21)
    result = interaction_split_test(Y, full, reduced, basis, permutations=15, seed=21, return_null=True)
    assert len(captured) == 15
    for draw, outcome in enumerate(captured):
        direct = _direct(outcome, full, reduced, basis)
        for name in INTERACTION_SPLIT_PARTS:
            assert result.parts[name].null_ss[draw] == pytest.approx(direct[name][0], rel=REL)
            assert result.parts[name].null_f[draw] == pytest.approx(direct[name][2], rel=REL)


def test_parts_share_each_draw():
    Y, full, reduced, basis = _inputs(seed=6)
    result = interaction_split_test(Y, full, reduced, basis, permutations=40, seed=2, return_null=True)
    np.testing.assert_allclose(result.span.null_ss + result.complement.null_ss, result.total.null_ss, rtol=REL)


def test_df_correction_leaves_the_f_null_and_scales_the_ss_null():
    Y, full, reduced, basis = _inputs(seed=7)
    plain = interaction_split_test(Y, full, reduced, basis, permutations=60, seed=8, return_null=True)
    corrected = interaction_split_test(
        Y, full, reduced, basis, permutations=60, seed=8, residual_correction="df", return_null=True
    )
    n = Y.shape[0]
    kappa_squared = n / (n - np.linalg.matrix_rank(reduced))
    for name in INTERACTION_SPLIT_PARTS:
        a, b = plain.parts[name], corrected.parts[name]
        assert b.p_value_f == a.p_value_f
        np.testing.assert_allclose(b.null_f, a.null_f, rtol=REL)
        np.testing.assert_allclose(b.null_ss, kappa_squared * a.null_ss, rtol=REL)
        assert (b.ss_interaction, b.rss, b.f) == (a.ss_interaction, a.rss, a.f)
    assert corrected.residual_correction == "df"


def test_leverage_correction_scales_rows_before_permutation():
    Y, full, reduced, basis = _inputs(seed=8)
    plain = interaction_split_test(Y, full, reduced, basis)
    corrected = interaction_split_test(
        Y, full, reduced, basis, permutations=20, seed=5, residual_correction="leverage", return_null=True
    )
    leverage = np.diag(_hat(reduced))
    brute = _brute_force_null(Y, full, reduced, basis, 20, 5, scale=1.0 / np.sqrt(1.0 - leverage))
    for name in INTERACTION_SPLIT_PARTS:
        part = corrected.parts[name]
        np.testing.assert_allclose(part.null_ss, brute[name][:, 0], rtol=REL)
        np.testing.assert_allclose(part.null_f, brute[name][:, 2], rtol=REL)
        observed = plain.parts[name]
        assert (part.ss_interaction, part.rss, part.f) == (observed.ss_interaction, observed.rss, observed.f)


def test_unknown_correction_is_rejected_before_any_draw(monkeypatch):
    Y, full, reduced, basis = _inputs()

    def no_draws(*args, **kwargs):  # pragma: no cover - must not be reached
        raise AssertionError("a permutation was drawn")

    monkeypatch.setattr(np.random, "default_rng", no_draws)
    with pytest.raises(ValueError, match=r"Unknown residual_correction 'uniform'.*'df', 'leverage'"):
        interaction_split_test(Y, full, reduced, basis, permutations=10, residual_correction="uniform")


def test_zero_permutations_give_no_p_values():
    Y, full, reduced, basis = _inputs()
    result = interaction_split_test(Y, full, reduced, basis, permutations=0, return_null=True)
    for name in INTERACTION_SPLIT_PARTS:
        part = result.parts[name]
        assert part.p_value_f is None and part.p_value_ss is None
        assert part.null_f is None and part.null_ss is None
        assert "p_value_f" not in part.to_dict() and "p_value_ss" not in part.to_dict()
    assert result.permutations == 0


def test_result_dictionary_holds_parts_and_settings():
    Y, full, reduced, basis = _inputs()
    result = interaction_split_test(Y, full, reduced, basis, permutations=9, seed=3, return_null=True)
    assert isinstance(result, InteractionSplitResult)
    payload = result.to_dict()
    assert list(payload["parts"]) == list(INTERACTION_SPLIT_PARTS)
    assert set(payload["parts"]["total"]) == {
        "dimension",
        "ss_interaction",
        "rss",
        "df_interaction",
        "df_residual",
        "f",
        "p_value_f",
        "p_value_ss",
    }
    assert (payload["permutations"], payload["seed"], payload["residual_correction"]) == (9, 3, None)


def test_high_dimensional_case_never_forms_a_permuted_outcome(monkeypatch):
    rng = np.random.default_rng(0)
    design, groups, stages = _design([[14, 13, 13], [13, 14, 13]])
    full = np.asarray(design.model_full, dtype=float)
    reduced = np.asarray(design.model_reduced, dtype=float)
    n, p = full.shape[0], 10_000
    Y = rng.normal(size=(n, p))
    basis, _ = np.linalg.qr(rng.normal(size=(p, 5)))
    shapes: set[tuple[tuple[int, ...], tuple[int, ...]]] = set()
    sums = interaction._sums

    def watch(gram, weights, df_interaction, trace):
        shapes.add((gram.shape, weights.shape))
        return sums(gram, weights, df_interaction, trace)

    # Every per-draw quantity goes through ``_sums``: it must only ever see n × n Gram matrices.
    monkeypatch.setattr(interaction, "_sums", watch)
    result = interaction_split_test(Y, full, reduced, basis, permutations=999, seed=1, return_null=True)
    assert shapes == {((n, n), (n, full.shape[1]))}
    assert result.total.null_ss.shape == (999,)
    assert [result.parts[name].dimension for name in INTERACTION_SPLIT_PARTS] == [p, 5, p - 5]


# ── Validation ───────────────────────────────────────────────────────────────


@pytest.mark.parametrize(
    "mutate, message",
    [
        (lambda Y, F, R, Q: (Y[:-1], F, R, Q), "Y has 35 rows but model_full has 36 rows"),
        (lambda Y, F, R, Q: (Y, F, R[:-1], Q), "Y has 36 rows but model_reduced has 35 rows"),
        (lambda Y, F, R, Q: (Y, F, R, Q[:-1]), "Y has 12 columns but basis has 11 rows"),
        (lambda Y, F, R, Q: (np.where(np.arange(Y.size).reshape(Y.shape) == 3, np.nan, Y), F, R, Q), "Y contains NaN"),
        (lambda Y, F, R, Q: (Y, F, R, Q * 1.01), "not orthonormal: the largest deviation .* is 0.0201"),
        (lambda Y, F, R, Q: (Y, F, R, np.eye(12)), "The complement is empty"),
        (lambda Y, F, R, Q: (Y, F, R, np.zeros((12, 0))), "The basis has no columns"),
        (
            lambda Y, F, R, Q: (Y, F, np.hstack([R, np.random.default_rng(0).normal(size=(36, 1))]), Q),
            "The models are not nested: column 4 of model_reduced",
        ),
        (lambda Y, F, R, Q: (Y, R, R, Q), "df_int must be positive"),
    ],
)
def test_invalid_inputs_are_rejected(mutate, message):
    Y, full, reduced, basis = _inputs(covariate=False)
    with pytest.raises(ValueError, match=message):
        interaction_split_test(*mutate(Y, full, reduced, basis))


def test_non_positive_residual_df_is_rejected():
    design, groups, stages = _design([[1, 1, 1], [1, 1, 1]])
    full = np.asarray(design.model_full, dtype=float)
    reduced = np.asarray(design.model_reduced, dtype=float)
    basis, _ = np.linalg.qr(np.random.default_rng(0).normal(size=(4, 2)))
    with pytest.raises(ValueError, match="df_res must be positive"):
        interaction_split_test(np.random.default_rng(1).normal(size=(6, 4)), full, reduced, basis)


def test_negative_permutations_are_rejected():
    Y, full, reduced, basis = _inputs()
    with pytest.raises(ValueError, match="permutations must be greater than or equal to 0"):
        interaction_split_test(Y, full, reduced, basis, permutations=-1)


def test_flat_and_direct_imports_are_the_same_object():
    from motco.stats import InteractionSplitResult as flat_result
    from motco.stats import interaction_split_test as flat
    from motco.stats.interaction import InteractionSplitResult as direct_result
    from motco.stats.interaction import interaction_split_test as direct

    assert flat is direct
    assert flat_result is direct_result
