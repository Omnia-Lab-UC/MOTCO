"""Joint selection of rank and keep fraction in ``plsda_doubleCV`` (regularized-stage-pls).

The selected pair on a scripted inner-AUROC pattern, both tie rules, the
``keep_fraction`` column, a ``[1.0]`` ladder reproducing the run without
sparsity, composition with the nuisance residualization and the strata, and
rejection before any fit. The default output's byte identity is pinned by
``tests/test_pls.py`` against ``tests/data/plsda_doublecv_identity.json``.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest
from sklearn.cross_decomposition import PLSRegression

from motco.stats import pls
from motco.stats.pls import NuisanceDesign, plsda_doubleCV
from motco.stats.sparse_pls import SparsePLSDA, SparsityDesign

CV = dict(cv1_splits=2, cv2_splits=3, n_repeats=2, max_components=5, random_state=11, progress=False)


def _three_stage(n_per_class: int = 12, p: int = 30, seed: int = 20261110) -> tuple[pd.DataFrame, pd.Series]:
    rng = np.random.default_rng(seed)
    labels = np.repeat(["Low", "Mid", "High"], n_per_class)
    X = rng.standard_normal((labels.size, p))
    X[:, :6] += 1.1 * (labels == "Mid")[:, None]
    X[:, 6:12] += 0.8 * (labels == "High")[:, None]
    return pd.DataFrame(X, columns=[f"f{i}" for i in range(p)]), pd.Series(labels)


def _scripted(monkeypatch, scores: dict[tuple[int, float], float], calls: list | None = None) -> None:
    """Replace the AUROC of every candidate by a scripted value (0.5 unless listed)."""

    def scripted(n_components, X_train, Y_train, X_test, Y_test, return_full=False, *, scale=True,
                 sparsity=None, keep_fraction=None):
        if calls is not None:
            calls.append((int(n_components), keep_fraction, len(X_train)))
        return scores.get((int(n_components), keep_fraction), 0.5)

    monkeypatch.setattr(pls, "_plsda_auroc", scripted)


def test_selected_pair_is_evaluated_and_recorded(monkeypatch):
    X, y = _three_stage()
    calls: list = []
    _scripted(monkeypatch, {(3, 0.2): 0.9, (2, 1.0): 0.8}, calls)
    result = plsda_doubleCV(X, y, sparsity=SparsityDesign("joint", (1.0, 0.5, 0.2)), **CV)
    table = result["table"]
    assert list(table.columns) == ["rep", "LV", "keep_fraction", "AUROC", "AUROC_std"]
    assert table["LV"].tolist() == [3, 3]
    assert table["keep_fraction"].tolist() == [0.2, 0.2]
    assert table["keep_fraction"].dtype == float
    # The outer evaluations (on the larger rest sets) are at the selected pair only.
    outer_size = max(size for _, _, size in calls)
    assert {(n, f) for n, f, size in calls if size == outer_size} == {(3, 0.2)}
    # The inner loop scored every pair of the grid.
    assert {(n, f) for n, f, _ in calls} == {(n, f) for n in range(1, 5) for f in (1.0, 0.5, 0.2)}
    for model in result["models"]:
        assert isinstance(model, SparsePLSDA)
        assert model.n_components == 3
        assert model.support_sizes().tolist() == [[6]] * 3  # round(0.2 · 30)


def test_tie_goes_to_the_fuller_budget(monkeypatch):
    X, y = _three_stage()
    _scripted(monkeypatch, {(2, 1.0): 0.8, (2, 0.5): 0.8, (3, 0.2): 0.7})
    result = plsda_doubleCV(X, y, sparsity=SparsityDesign("joint", (0.5, 0.2, 1.0)), **CV)
    assert result["table"]["LV"].tolist() == [2, 2]
    assert result["table"]["keep_fraction"].tolist() == [1.0, 1.0]


def test_tie_goes_to_the_smaller_rank(monkeypatch):
    X, y = _three_stage()
    _scripted(monkeypatch, {(3, 0.5): 0.8, (2, 0.2): 0.8})
    result = plsda_doubleCV(X, y, sparsity=SparsityDesign("joint", (1.0, 0.5, 0.2)), **CV)
    assert result["table"]["LV"].tolist() == [2, 2]
    assert result["table"]["keep_fraction"].tolist() == [0.2, 0.2]


def test_modal_fraction_breaks_ties_toward_the_fuller_budget():
    assert pls._modal_fraction_fullest([0.2, 0.5, 0.2, 0.5]) == 0.5
    assert pls._modal_fraction_fullest([0.2, 0.2, 0.5]) == 0.2
    assert pls._modal_fraction_fullest([0.05]) == 0.05
    with pytest.raises(ValueError):
        pls._modal_fraction_fullest([])


@pytest.mark.parametrize("scope_groups", [("joint", None), ("block", (20, 10))])
def test_ladder_of_one_full_budget_reproduces_production(scope_groups):
    X, y = _three_stage()
    scope, groups = scope_groups
    default = plsda_doubleCV(X, y, **CV)
    sparse = plsda_doubleCV(X, y, sparsity=SparsityDesign(scope, (1.0,), groups), **CV)
    # The full budget is PLSRegression.fit, so the agreement is exact (the contract is 1e-6 up to signs).
    pd.testing.assert_frame_equal(sparse["table"].drop(columns="keep_fraction"), default["table"])
    assert sparse["table"]["keep_fraction"].tolist() == [1.0, 1.0]
    for ours, theirs in zip(sparse["models"], default["models"], strict=True):
        np.testing.assert_array_equal(ours.x_scores_, theirs.x_scores_)


def test_column_is_absent_without_sparsity():
    X, y = _three_stage()
    result = plsda_doubleCV(X, y, **CV)
    assert list(result["table"].columns) == ["rep", "LV", "AUROC", "AUROC_std"]
    assert all(type(model) is PLSRegression for model in result["models"])


def test_sparse_cv_runs_end_to_end_in_both_scopes():
    X, y = _three_stage()
    for design in (SparsityDesign("joint"), SparsityDesign("block", feature_groups=(18, 12))):
        result = plsda_doubleCV(X, y, sparsity=design, **CV)
        table = result["table"]
        assert set(table["keep_fraction"]) <= set(design.keep_fractions)
        for (_, row), model in zip(table.iterrows(), result["models"], strict=True):
            assert model.n_components == row["LV"]
            expected = design.resolve(30).counts(row["keep_fraction"])
            assert model.keep_.tolist() == [list(expected)] * int(row["LV"])
            assert 0.0 <= row["AUROC"] <= 1.0


def test_parallel_inner_loop_equals_serial():
    X, y = _three_stage()
    design = SparsityDesign("block", (1.0, 0.2), (18, 12))
    serial = plsda_doubleCV(X, y, sparsity=design, **CV)
    parallel = plsda_doubleCV(X, y, sparsity=design, n_jobs=2, **CV)
    pd.testing.assert_frame_equal(serial["table"], parallel["table"])


def _stage_nuisance(X: pd.DataFrame, y: pd.Series, seed: int = 3) -> NuisanceDesign:
    rng = np.random.default_rng(seed)
    stages = pd.Categorical(y, categories=["Low", "Mid", "High"]).codes
    conditioning = np.column_stack([(stages == 1).astype(float), (stages == 2).astype(float)])
    z = stages + rng.normal(size=stages.size)
    return NuisanceDesign(conditioning, z[:, None], ("stage[1]", "stage[2]"), ("z",), "stage")


def test_nuisance_and_strata_compose_with_the_grid():
    X, y = _three_stage()
    nuisance = _stage_nuisance(X, y)
    options = dict(scale=False, nuisance=nuisance, strata=nuisance.cells(), **CV)
    plain = plsda_doubleCV(X, y, **options)
    full = plsda_doubleCV(X, y, sparsity=SparsityDesign("joint", (1.0,)), **options)
    pd.testing.assert_frame_equal(full["table"].drop(columns="keep_fraction"), plain["table"])
    assert full["nuisance_min_training_df"] == plain["nuisance_min_training_df"]
    for ours, theirs in zip(full["models"], plain["models"], strict=True):
        np.testing.assert_array_equal(ours.x_scores_, theirs.x_scores_)
    sparse = plsda_doubleCV(X, y, sparsity=SparsityDesign("block", (0.5, 0.1), (18, 12)), **options)
    assert "keep_fraction" in sparse["table"].columns
    assert all(not model.scale for model in sparse["models"])


def test_membership_mismatch_is_rejected_before_any_fit(monkeypatch):
    X, y = _three_stage()

    def forbidden(*args, **kwargs):
        raise AssertionError("no fold may be fitted when the sparsity design is rejected")

    monkeypatch.setattr(pls, "_plsda_auroc", forbidden)
    with pytest.raises(ValueError, match=r"cover 29 features, but X has 30"):
        plsda_doubleCV(X, y, sparsity=SparsityDesign("block", (0.5,), (20, 9)), **CV)
    with pytest.raises(ValueError, match=r"sparsity must be a SparsityDesign"):
        plsda_doubleCV(X, y, sparsity={"scope": "joint"}, **CV)  # type: ignore[arg-type]
