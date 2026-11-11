"""The per-feature interaction scan: definitions, degenerate features, permutation path, step-down and validation."""

from __future__ import annotations

from time import perf_counter

import numpy as np
import pandas as pd
import pytest

from motco.simulations.evaluation import build_simulation_trajectory_design
from motco.stats import interaction, permutation
from motco.stats.interaction import (
    INTERACTION_SCAN_ALPHA,
    InteractionScanResult,
    interaction_feature_scan,
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


def _inputs(seed: int = 0, p: int = 12, sizes: list[list[int]] | None = None, covariate: bool = True):
    """Unequal cells, an optional numeric covariate in both models, an interaction on a few features."""

    rng = np.random.default_rng(seed)
    design, groups, stages = _design(sizes or [[6, 5, 9], [4, 7, 5]])
    full = np.asarray(design.model_full, dtype=float)
    reduced = np.asarray(design.model_reduced, dtype=float)
    if covariate:
        z = rng.normal(size=(groups.size, 1))
        full = np.hstack([full, z])
        reduced = np.hstack([reduced, z])
    loading = np.zeros(p)
    loading[: max(1, p // 4)] = rng.normal(scale=1.5, size=max(1, p // 4))
    Y = rng.normal(size=(3, p))[stages] + np.outer(groups * stages, loading)
    Y = Y + rng.normal(size=Y.shape) * rng.uniform(0.5, 3.0, size=p)
    return Y, full, reduced


def _hat(X: np.ndarray) -> np.ndarray:
    return X @ np.linalg.pinv(X)


def _direct(Y: np.ndarray, full: np.ndarray, reduced: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    n = Y.shape[0]
    H_full, H_reduced = _hat(full), _hat(reduced)
    df_int = np.linalg.matrix_rank(full) - np.linalg.matrix_rank(reduced)
    df_res = n - np.linalg.matrix_rank(full)
    ss = np.sum(((H_full - H_reduced) @ Y) ** 2, axis=0)
    rss = np.sum(((np.eye(n) - H_full) @ Y) ** 2, axis=0)
    return ss, rss, (ss / df_int) / (rss / df_res)


def _brute_force_null(Y, full, reduced, permutations, seed, scale=None) -> np.ndarray:
    """Refit Ŷ_R + s·E[π] per draw, with π from ``np.random.default_rng(seed)``; draws × columns F."""

    n = Y.shape[0]
    fitted = _hat(reduced) @ Y
    residuals = Y - fitted
    if scale is not None:
        residuals = residuals * scale[:, None]
    rng = np.random.default_rng(seed)
    return np.array([_direct(fitted + residuals[rng.permutation(n)], full, reduced)[2] for _ in range(permutations)])


def _step_down(f: np.ndarray, null_f: np.ndarray, valid: np.ndarray) -> np.ndarray:
    """A direct, loop-based Westfall–Young step-down max-T on the returned null."""

    columns = sorted(valid.tolist(), key=lambda column: (-f[column], column))
    B = null_f.shape[0]
    adjusted = np.full(f.size, np.nan)
    running = 0.0
    for rank, column in enumerate(columns):
        tail = columns[rank:]
        count = sum(1 for b in range(B) if max(null_f[b, c] for c in tail) >= f[column])
        running = max(running, (1 + count) / (1 + B))
        adjusted[column] = running
    return adjusted


# ── Observed quantities ──────────────────────────────────────────────────────


def test_definitions_match_direct_hat_matrices():
    Y, full, reduced = _inputs()
    result = interaction_feature_scan(Y, full, reduced)
    ss, rss, f_values = _direct(Y, full, reduced)
    np.testing.assert_allclose(result.ss_interaction, ss, rtol=REL)
    np.testing.assert_allclose(result.rss, rss, rtol=REL)
    np.testing.assert_allclose(result.f, f_values, rtol=REL)
    assert (result.df_interaction, result.df_residual) == (2, 36 - 7)
    assert result.max_column == int(np.argmax(f_values))
    assert result.max_f == pytest.approx(float(f_values.max()), rel=REL)
    assert result.n_features == 12 and result.n_degenerate == 0
    assert not result.degenerate.any()


def test_each_feature_ignores_its_location_and_scale():
    Y, full, reduced = _inputs(seed=1)
    base = interaction_feature_scan(Y, full, reduced)
    moved = Y.copy()
    moved[:, 3] = -7.5 * Y[:, 3] + 40.0
    result = interaction_feature_scan(moved, full, reduced)
    np.testing.assert_allclose(result.f, base.f, rtol=REL)


def test_feature_sums_add_to_the_split_tests_total():
    for seed in range(3):
        Y, full, reduced = _inputs(seed=seed)
        basis, _ = np.linalg.qr(np.random.default_rng(seed).normal(size=(Y.shape[1], 2)))
        scan = interaction_feature_scan(Y, full, reduced)
        split = interaction_split_test(Y, full, reduced, basis)
        assert float(scan.ss_interaction.sum()) == pytest.approx(split.total.ss_interaction, rel=REL)
        assert float(scan.rss.sum()) == pytest.approx(split.total.rss, rel=REL)


def test_constant_feature_is_degenerate_and_ignored():
    Y, full, reduced = _inputs(seed=2)
    with_constant = np.hstack([Y[:, :5], np.full((Y.shape[0], 1), 0.1), Y[:, 5:]])
    base = interaction_feature_scan(Y, full, reduced, permutations=40, seed=3)
    result = interaction_feature_scan(with_constant, full, reduced, permutations=40, seed=3)
    assert result.degenerate[5] and result.n_degenerate == 1
    assert list(result.degenerate_columns) == [5]
    assert np.isnan(result.f[5]) and np.isnan(result.p_value_adjusted[5])
    keep = np.r_[0:5, 6:13]
    np.testing.assert_allclose(result.f[keep], base.f, rtol=REL)
    np.testing.assert_array_equal(result.p_value_adjusted[keep], base.p_value_adjusted)
    # A product over a different column set can move the last bit, so the F are compared to tolerance.
    assert result.max_f == pytest.approx(base.max_f, rel=REL)
    assert result.max_column == int(keep[base.max_column])
    assert result.p_value_global == base.p_value_global
    np.testing.assert_allclose(result.null_max_f, base.null_max_f, rtol=REL)
    assert 5 not in set(result.null_max_column.tolist())


def test_feature_fitted_exactly_by_the_full_model_is_degenerate():
    Y, full, reduced = _inputs(seed=3)
    exact = full @ np.random.default_rng(0).normal(size=full.shape[1])
    result = interaction_feature_scan(np.column_stack([Y, exact]), full, reduced)
    assert result.degenerate[-1] and not result.degenerate[:-1].any()


def test_every_feature_degenerate_is_rejected_before_any_draw(monkeypatch):
    _, full, reduced = _inputs()

    def no_draws(*args, **kwargs):  # pragma: no cover - must not be reached
        raise AssertionError("a permutation was drawn")

    monkeypatch.setattr(np.random, "default_rng", no_draws)
    constant = np.ones((full.shape[0], 3))
    with pytest.raises(ValueError, match="No feature has a defined F: all 3 outcome column"):
        interaction_feature_scan(constant, full, reduced, permutations=10, seed=0)


# ── Permutation path ─────────────────────────────────────────────────────────


def test_null_f_match_brute_force_refits():
    Y, full, reduced = _inputs(seed=4)
    result = interaction_feature_scan(Y, full, reduced, permutations=25, seed=11, return_null=True)
    brute = _brute_force_null(Y, full, reduced, 25, 11)
    assert result.null_f.shape == (25, Y.shape[1])
    np.testing.assert_allclose(result.null_f, brute, rtol=REL)


def test_permuted_outcomes_are_rrpps_serial_draws(monkeypatch):
    Y, full, reduced = _inputs(seed=5, p=4)
    ls_means = np.zeros((6, full.shape[1]))
    ls_means[:, 0] = 1.0
    captured: list[np.ndarray] = []

    def capture(y_random, *args, **kwargs):
        captured.append(np.array(y_random, dtype=float))
        return np.zeros((2, 2)), np.zeros((2, 2)), np.zeros((2, 2))

    monkeypatch.setattr(permutation, "estimate_difference", capture)
    RRPP(Y, full, reduced, ls_means, [[0, 1, 2], [3, 4, 5]], permutations=15, progress=False, seed=21)
    result = interaction_feature_scan(Y, full, reduced, permutations=15, seed=21, return_null=True)
    assert len(captured) == 15
    for draw, outcome in enumerate(captured):
        np.testing.assert_allclose(result.null_f[draw], _direct(outcome, full, reduced)[2], rtol=REL)


def test_global_p_value_is_the_permutation_p_value_of_the_maximum():
    Y, full, reduced = _inputs(seed=6)
    B = 60
    result = interaction_feature_scan(Y, full, reduced, permutations=B, seed=2, return_null=True)
    maxima = np.nanmax(result.null_f, axis=1)
    np.testing.assert_array_equal(result.null_max_f, maxima)
    assert result.p_value_global == (1 + np.sum(maxima >= result.max_f)) / (1 + B)


def test_adjusted_p_values_match_a_direct_step_down():
    for seed in range(3):
        Y, full, reduced = _inputs(seed=seed + 10, p=9)
        result = interaction_feature_scan(Y, full, reduced, permutations=40, seed=seed, return_null=True)
        expected = _step_down(result.f, result.null_f, np.flatnonzero(~result.degenerate))
        np.testing.assert_array_equal(result.p_value_adjusted, expected)


def test_step_down_breaks_tied_f_by_column_index(monkeypatch):
    pseudo_f = interaction._pseudo_f

    def rounded(*args):
        # A duplicated column's F can differ in the last bit; rounding makes the tie exact.
        return np.round(pseudo_f(*args), 8)

    monkeypatch.setattr(interaction, "_pseudo_f", rounded)
    Y, full, reduced = _inputs(seed=7, p=6)
    tied = np.column_stack([Y, Y[:, 2]])
    result = interaction_feature_scan(tied, full, reduced, permutations=40, seed=1, return_null=True)
    assert result.f[6] == result.f[2]
    expected = _step_down(result.f, result.null_f, np.arange(7))
    np.testing.assert_array_equal(result.p_value_adjusted, expected)
    assert result.p_value_adjusted[2] <= result.p_value_adjusted[6]


def test_adjusted_p_values_are_monotone_and_the_top_one_is_global():
    Y, full, reduced = _inputs(seed=8, p=20)
    result = interaction_feature_scan(Y, full, reduced, permutations=99, seed=4)
    order = np.argsort(-result.f, kind="stable")
    assert np.all(np.diff(result.p_value_adjusted[order]) >= 0)
    assert result.p_value_adjusted[result.max_column] == result.p_value_global
    assert np.all(result.p_value_adjusted >= result.p_value_global)


def test_null_maxima_name_columns_attaining_them():
    Y, full, reduced = _inputs(seed=9)
    result = interaction_feature_scan(Y, full, reduced, permutations=50, seed=6, return_null=True)
    assert result.null_max_f.shape == (50,) and result.null_max_column.shape == (50,)
    for draw in range(50):
        column = int(result.null_max_column[draw])
        assert not result.degenerate[column]
        assert result.null_f[draw, column] == result.null_max_f[draw]


def test_null_f_are_kept_only_on_request():
    Y, full, reduced = _inputs()
    result = interaction_feature_scan(Y, full, reduced, permutations=10, seed=0)
    assert result.null_f is None
    assert result.null_max_f.shape == (10,) and result.null_max_column.shape == (10,)


def test_hits_are_ordered_by_f_at_the_fixed_level():
    rng = np.random.default_rng(12)
    design, groups, stages = _design([[10, 10, 10], [10, 10, 10]])
    full = np.asarray(design.model_full, dtype=float)
    reduced = np.asarray(design.model_reduced, dtype=float)
    Y = rng.normal(size=(60, 30))
    Y[:, [4, 17]] += np.outer(groups * stages, [3.0, 2.0])
    result = interaction_feature_scan(Y, full, reduced, permutations=199, seed=3)
    hits = result.hits()
    assert INTERACTION_SCAN_ALPHA == 0.05
    assert set(hits.tolist()) >= {4, 17}
    assert np.all(np.diff(result.f[hits]) <= 0)
    assert np.all(result.p_value_adjusted[hits] <= 0.05)
    payload = result.to_dict()
    assert [hit["column"] for hit in payload["hits"]] == hits.tolist()
    assert payload["hit_alpha"] == 0.05 and payload["p_value_global"] == result.p_value_global


def test_df_correction_leaves_every_f_and_p_value_unchanged():
    Y, full, reduced = _inputs(seed=13)
    plain = interaction_feature_scan(Y, full, reduced, permutations=60, seed=8, return_null=True)
    corrected = interaction_feature_scan(
        Y, full, reduced, permutations=60, seed=8, residual_correction="df", return_null=True
    )
    np.testing.assert_array_equal(corrected.f, plain.f)
    np.testing.assert_array_equal(corrected.null_f, plain.null_f)
    np.testing.assert_array_equal(corrected.p_value_adjusted, plain.p_value_adjusted)
    assert corrected.p_value_global == plain.p_value_global
    assert corrected.residual_correction == "df"


def test_leverage_correction_scales_rows_before_permutation():
    Y, full, reduced = _inputs(seed=14)
    plain = interaction_feature_scan(Y, full, reduced)
    corrected = interaction_feature_scan(
        Y, full, reduced, permutations=20, seed=5, residual_correction="leverage", return_null=True
    )
    leverage = np.diag(_hat(reduced))
    assert np.ptp(leverage) > 0.01  # unequal cells and a covariate: the multipliers vary
    brute = _brute_force_null(Y, full, reduced, 20, 5, scale=1.0 / np.sqrt(1.0 - leverage))
    np.testing.assert_allclose(corrected.null_f, brute, rtol=REL)
    np.testing.assert_array_equal(corrected.f, plain.f)


def test_unknown_correction_is_rejected_before_any_draw(monkeypatch):
    Y, full, reduced = _inputs()

    def no_draws(*args, **kwargs):  # pragma: no cover - must not be reached
        raise AssertionError("a permutation was drawn")

    monkeypatch.setattr(np.random, "default_rng", no_draws)
    with pytest.raises(ValueError, match=r"Unknown residual_correction 'uniform'.*'df', 'leverage'"):
        interaction_feature_scan(Y, full, reduced, permutations=10, residual_correction="uniform")


def test_zero_permutations_give_no_p_values():
    Y, full, reduced = _inputs()
    result = interaction_feature_scan(Y, full, reduced, permutations=0, return_null=True)
    assert isinstance(result, InteractionScanResult)
    assert result.p_value_global is None and result.p_value_adjusted is None
    assert result.null_max_f is None and result.null_max_column is None and result.null_f is None
    assert result.permutations == 0
    payload = result.to_dict()
    assert {"p_value_global", "hits", "hit_alpha"}.isdisjoint(payload)
    assert payload["max_column"] == result.max_column and payload["n_features"] == 12
    with pytest.raises(ValueError, match="without permutations"):
        result.hits()


def test_high_dimensional_case_never_forms_a_permuted_outcome(monkeypatch, record_property):
    rng = np.random.default_rng(0)
    design, groups, stages = _design([[14, 13, 13], [13, 14, 13]])
    covariates = rng.normal(size=(80, 6))
    full = np.hstack([np.asarray(design.model_full, dtype=float), covariates])
    reduced = np.hstack([np.asarray(design.model_reduced, dtype=float), covariates])
    n, p = full.shape[0], 10_000
    assert full.shape[1] == 12
    Y = rng.normal(size=(n, p))
    calls: list[tuple[int, tuple[int, ...], tuple[int, ...]]] = []
    sums = interaction._feature_sums

    def watch(residuals, weights, df_interaction, norms):
        calls.append((id(residuals), residuals.shape, weights.shape))
        return sums(residuals, weights, df_interaction, norms)

    # Every per-draw quantity goes through ``_feature_sums``: the n × p residuals are one fixed matrix
    # and only the n × rank(F) weights are permuted.
    monkeypatch.setattr(interaction, "_feature_sums", watch)
    start = perf_counter()
    result = interaction_feature_scan(Y, full, reduced, permutations=999, seed=1, residual_correction="leverage")
    elapsed = perf_counter() - start
    record_property("scan_seconds_10000x80x999", elapsed)
    print(f"interaction_feature_scan: p = {p}, n = {n}, 12 model columns, 999 draws in {elapsed:.2f} s")
    draws = calls[1:]
    assert len(draws) == 999
    assert {(shape, weights) for _, shape, weights in draws} == {((n, p), (n, 12))}
    assert len({identifier for identifier, _, _ in draws}) == 1
    assert result.null_max_f.shape == (999,)


# ── By-block combination ─────────────────────────────────────────────────────


def _labels(p: int, sizes: tuple[int, ...] = (4, 8), names: tuple[str, ...] = ("rna", "atac")) -> list[str]:
    """Contiguous block labels of the given sizes, which must sum to ``p``."""

    assert sum(sizes) == p
    return [name for name, size in zip(names, sizes) for _ in range(size)]


def _block_step_down(f: np.ndarray, null_f: np.ndarray, degenerate: np.ndarray, labels: list[str]) -> np.ndarray:
    """min(1, K × the direct step-down applied to each block alone), K the blocks with a non-degenerate column."""

    by_block = {}
    for name in dict.fromkeys(labels):
        valid = np.array([c for c, label in enumerate(labels) if label == name and not degenerate[c]], dtype=int)
        if valid.size:
            by_block[name] = valid
    expected = np.full(f.size, np.nan)
    for valid in by_block.values():
        within = _step_down(f, null_f, valid)
        expected[valid] = np.minimum(1.0, len(by_block) * within[valid])
    return expected


def test_by_block_adjusted_p_values_match_a_brute_force_per_block_step_down():
    for seed in range(3):
        Y, full, reduced = _inputs(seed=seed + 20, p=12)
        labels = _labels(12)
        result = interaction_feature_scan(Y, full, reduced, permutations=40, seed=seed, return_null=True, blocks=labels)
        expected = _block_step_down(result.f, result.null_f, result.degenerate, labels)
        np.testing.assert_array_equal(result.p_value_adjusted, expected)
        assert result.combination == "block_bonferroni" and result.n_blocks_tested == 2


def test_by_block_handles_interleaved_labels_and_names_blocks_by_first_appearance():
    Y, full, reduced = _inputs(seed=24, p=10)
    labels = ["b", "a", "b", "c", "a", "a", "b", "c", "c", "a"]
    result = interaction_feature_scan(Y, full, reduced, permutations=40, seed=3, return_null=True, blocks=labels)
    assert result.block_names == ("b", "a", "c")
    assert [block.columns.tolist() for block in result.blocks] == [[0, 2, 6], [1, 4, 5, 9], [3, 7, 8]]
    expected = _block_step_down(result.f, result.null_f, result.degenerate, labels)
    np.testing.assert_array_equal(result.p_value_adjusted, expected)
    assert result.n_blocks_tested == 3


def test_by_block_global_p_value_is_the_smallest_adjusted_p_value():
    Y, full, reduced = _inputs(seed=25, p=16)
    labels = _labels(16, (6, 10))
    result = interaction_feature_scan(Y, full, reduced, permutations=99, seed=4, blocks=labels)
    block_p = [block.p_value for block in result.blocks]
    assert result.p_value_global == float(np.nanmin(result.p_value_adjusted))
    assert result.p_value_global == min(1.0, 2 * min(block_p))
    for block in result.blocks:
        # The block's own p-value is its top feature's within-block value: K times it, capped, is that feature's.
        assert result.p_value_adjusted[block.max_column] == min(1.0, 2 * block.p_value)
        assert block.p_value == (1 + np.sum(block.null_max_f >= block.max_f)) / (1 + 99)


def test_by_block_adjusted_p_values_are_monotone_in_f_within_a_block():
    Y, full, reduced = _inputs(seed=26, p=20)
    labels = _labels(20, (5, 15))
    result = interaction_feature_scan(Y, full, reduced, permutations=99, seed=5, blocks=labels)
    for block in result.blocks:
        columns = block.columns[np.argsort(-result.f[block.columns], kind="stable")]
        assert np.all(np.diff(result.p_value_adjusted[columns]) >= 0)


def test_one_block_label_reproduces_the_joint_scan_exactly():
    for seed in range(3):
        Y, full, reduced = _inputs(seed=seed + 30, p=14)
        joint = interaction_feature_scan(Y, full, reduced, permutations=60, seed=seed)
        single = interaction_feature_scan(Y, full, reduced, permutations=60, seed=seed, blocks=["all"] * 14)
        assert single.p_value_global == joint.p_value_global
        np.testing.assert_array_equal(single.p_value_adjusted, joint.p_value_adjusted)
        np.testing.assert_array_equal(single.hits(), joint.hits())
        np.testing.assert_array_equal(single.null_max_f, joint.null_max_f)
        assert single.n_blocks_tested == 1 and single.blocks[0].p_value == joint.p_value_global
        np.testing.assert_array_equal(single.null_block_max_f[:, 0], joint.null_max_f)


def test_by_block_keeps_every_joint_quantity():
    Y, full, reduced = _inputs(seed=33, p=12)
    joint = interaction_feature_scan(Y, full, reduced, permutations=50, seed=9, return_null=True)
    by_block = interaction_feature_scan(Y, full, reduced, permutations=50, seed=9, return_null=True, blocks=_labels(12))
    for name in ("ss_interaction", "rss", "f", "degenerate", "null_max_f", "null_max_column", "null_f"):
        np.testing.assert_array_equal(getattr(by_block, name), getattr(joint, name))
    assert (by_block.max_f, by_block.max_column) == (joint.max_f, joint.max_column)


def test_joint_scan_has_no_block_fields():
    Y, full, reduced = _inputs(seed=34)
    result = interaction_feature_scan(Y, full, reduced, permutations=20, seed=1)
    assert result.combination == "joint"
    assert result.blocks is None and result.block_names is None
    assert result.n_blocks_tested is None and result.null_block_max_f is None
    assert {"combination", "n_blocks_tested", "blocks"}.isdisjoint(result.to_dict())


def test_per_block_null_maxima_are_attained_by_a_column_of_the_block():
    Y, full, reduced = _inputs(seed=35, p=12)
    result = interaction_feature_scan(Y, full, reduced, permutations=50, seed=6, return_null=True, blocks=_labels(12))
    assert result.null_block_max_f.shape == (50, 2)
    for slot, block in enumerate(result.blocks):
        assert block.null_max_f.shape == (50,)
        np.testing.assert_array_equal(block.null_max_f, result.null_block_max_f[:, slot])
        np.testing.assert_array_equal(block.null_max_f, np.nanmax(result.null_f[:, block.columns], axis=1))
        assert block.max_f == float(np.nanmax(result.f[block.columns]))
        assert block.max_column in set(block.columns.tolist())
    np.testing.assert_array_equal(result.null_max_f, result.null_block_max_f.max(axis=1))


def test_a_wholly_degenerate_block_is_not_tested():
    Y, full, reduced = _inputs(seed=36, p=10)
    constant = np.full((Y.shape[0], 3), 2.5)
    outcome = np.hstack([Y[:, :4], constant, Y[:, 4:]])
    labels = ["x"] * 4 + ["dead"] * 3 + ["y"] * 6
    result = interaction_feature_scan(outcome, full, reduced, permutations=40, seed=2, return_null=True, blocks=labels)
    dead = result.blocks[1]
    assert result.block_names == ("x", "dead", "y")
    assert dead.name == "dead" and not dead.tested and dead.n_degenerate == 3 and dead.n_features == 3
    assert dead.max_f is None and dead.max_column is None and dead.p_value is None and dead.null_max_f is None
    assert result.n_blocks_tested == 2 and result.null_block_max_f.shape == (40, 2)
    assert np.all(np.isnan(result.p_value_adjusted[4:7]))
    expected = _block_step_down(result.f, result.null_f, result.degenerate, labels)
    np.testing.assert_array_equal(result.p_value_adjusted, expected)
    assert result.to_dict()["blocks"][1] == {
        "name": "dead",
        "n_features": 3,
        "n_degenerate": 3,
        "max_f": None,
        "max_column": None,
        "p_value": None,
    }


def test_by_block_without_permutations_reports_the_block_maxima_only():
    Y, full, reduced = _inputs(seed=37, p=12)
    result = interaction_feature_scan(Y, full, reduced, blocks=_labels(12))
    assert result.combination == "block_bonferroni" and result.n_blocks_tested == 2
    assert result.p_value_global is None and result.p_value_adjusted is None and result.null_block_max_f is None
    for block in result.blocks:
        assert block.p_value is None and block.null_max_f is None
        assert block.max_f == float(np.nanmax(result.f[block.columns]))
    payload = result.to_dict()
    assert payload["combination"] == "block_bonferroni" and payload["n_blocks_tested"] == 2
    assert all("p_value" not in entry for entry in payload["blocks"])


def test_by_block_serialization_adds_the_block_fields():
    Y, full, reduced = _inputs(seed=38, p=12)
    result = interaction_feature_scan(Y, full, reduced, permutations=99, seed=7, blocks=_labels(12))
    payload = result.to_dict()
    joint = interaction_feature_scan(Y, full, reduced, permutations=99, seed=7).to_dict()
    assert set(payload) - set(joint) == {"combination", "n_blocks_tested", "blocks"}
    assert [entry["name"] for entry in payload["blocks"]] == ["rna", "atac"]
    for entry, block in zip(payload["blocks"], result.blocks):
        assert entry == {
            "name": block.name,
            "n_features": block.n_features,
            "n_degenerate": 0,
            "max_f": block.max_f,
            "max_column": block.max_column,
            "p_value": block.p_value,
        }
    assert payload["p_value_global"] == result.p_value_global
    assert [hit["column"] for hit in payload["hits"]] == result.hits().tolist()


def test_by_block_hits_read_the_by_block_p_values():
    rng = np.random.default_rng(40)
    design, groups, stages = _design([[10, 10, 10], [10, 10, 10]])
    full = np.asarray(design.model_full, dtype=float)
    reduced = np.asarray(design.model_reduced, dtype=float)
    Y = rng.normal(size=(60, 40))
    Y[:, [3, 30]] += np.outer(groups * stages, [3.0, 3.0])
    labels = _labels(40, (5, 35), ("small", "large"))
    result = interaction_feature_scan(Y, full, reduced, permutations=199, seed=3, blocks=labels)
    hits = result.hits()
    assert set(hits.tolist()) >= {3, 30}
    called = np.flatnonzero(result.p_value_adjusted <= 0.05)
    np.testing.assert_array_equal(hits, called[np.lexsort((called, -result.f[called]))])


@pytest.mark.parametrize("labels, message", [(["a"] * 9, "blocks has 9 label.*Y has 10 column"), ("rna", "single str")])
def test_malformed_block_labels_are_rejected_before_any_draw(monkeypatch, labels, message):
    Y, full, reduced = _inputs(p=10)

    def no_draws(*args, **kwargs):  # pragma: no cover - must not be reached
        raise AssertionError("a permutation was drawn")

    monkeypatch.setattr(np.random, "default_rng", no_draws)
    with pytest.raises(ValueError, match=message):
        interaction_feature_scan(Y, full, reduced, permutations=10, seed=0, blocks=labels)


def test_block_labels_accept_arrays_and_indexes():
    Y, full, reduced = _inputs(seed=41, p=12)
    labels = _labels(12)
    listed = interaction_feature_scan(Y, full, reduced, permutations=30, seed=2, blocks=labels)
    for other in (np.asarray(labels), pd.Index(labels), tuple(labels)):
        result = interaction_feature_scan(Y, full, reduced, permutations=30, seed=2, blocks=other)
        assert result.block_names == ("rna", "atac")
        np.testing.assert_array_equal(result.p_value_adjusted, listed.p_value_adjusted)


# ── Validation ───────────────────────────────────────────────────────────────


@pytest.mark.parametrize(
    "mutate, message",
    [
        (lambda Y, F, R: (Y[:-1], F, R), "Y has 35 rows but model_full has 36 rows"),
        (lambda Y, F, R: (Y, F, R[:-1]), "Y has 36 rows but model_reduced has 35 rows"),
        (lambda Y, F, R: (Y[:, 0], F, R), "Y must be a 2-D matrix"),
        (lambda Y, F, R: (np.where(np.arange(Y.size).reshape(Y.shape) == 3, np.nan, Y), F, R), "Y contains NaN"),
        (lambda Y, F, R: (Y, np.where(np.arange(F.size).reshape(F.shape) == 3, np.inf, F), R), "model_full contains"),
        (
            lambda Y, F, R: (Y, F, np.hstack([R, np.random.default_rng(0).normal(size=(36, 1))])),
            "The models are not nested: column 4 of model_reduced",
        ),
        (lambda Y, F, R: (Y, R, R), "df_int must be positive"),
    ],
)
def test_invalid_inputs_are_rejected(mutate, message):
    Y, full, reduced = _inputs(covariate=False)
    with pytest.raises(ValueError, match=message):
        interaction_feature_scan(*mutate(Y, full, reduced))


def test_non_positive_residual_df_is_rejected():
    design, groups, stages = _design([[1, 1, 1], [1, 1, 1]])
    full = np.asarray(design.model_full, dtype=float)
    reduced = np.asarray(design.model_reduced, dtype=float)
    with pytest.raises(ValueError, match="df_res must be positive"):
        interaction_feature_scan(np.random.default_rng(1).normal(size=(6, 4)), full, reduced)


def test_negative_permutations_are_rejected():
    Y, full, reduced = _inputs()
    with pytest.raises(ValueError, match="permutations must be greater than or equal to 0"):
        interaction_feature_scan(Y, full, reduced, permutations=-1)


def test_dataframe_outcome_is_accepted():
    Y, full, reduced = _inputs()
    frame = pd.DataFrame(Y, columns=[f"f{i}" for i in range(Y.shape[1])])
    np.testing.assert_array_equal(
        interaction_feature_scan(frame, full, reduced).f, interaction_feature_scan(Y, full, reduced).f
    )


def test_flat_and_direct_imports_are_the_same_object():
    from motco.stats import InteractionScanResult as flat_result
    from motco.stats import interaction_feature_scan as flat
    from motco.stats.interaction import InteractionScanResult as direct_result
    from motco.stats.interaction import interaction_feature_scan as direct

    assert flat is direct
    assert flat_result is direct_result
