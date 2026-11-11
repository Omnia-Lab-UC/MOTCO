"""The sparse stage-PLS-DA estimator (regularized-stage-pls).

Budget semantics in both scopes, the soft threshold, determinism, rejection of
malformed budgets before any fit, equality with ``PLSRegression`` at the full
budget, agreement with the committed mixOmics fixture, and the module
contract.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest
from sklearn.base import clone
from sklearn.cross_decomposition import PLSRegression

from motco.simulations.diagnostics import LinearLatentMap
from motco.stats.pls import calculate_vips
from motco.stats.sparse_pls import (
    DEFAULT_KEEP_FRACTIONS,
    SPARSITY_SCOPES,
    SparsePLSDA,
    SparsityDesign,
    keep_count,
    soft_threshold,
)

FIXTURE = Path(__file__).parent / "data" / "splsda_reference_fixture.npz"


def _three_class(n_features: int, *, n_per_class: int = 20, seed: int = 20261110) -> tuple[np.ndarray, np.ndarray]:
    rng = np.random.default_rng(seed)
    labels = np.repeat([0, 1, 2], n_per_class)
    X = rng.standard_normal((labels.size, n_features))
    X[:, :8] += 0.9 * (labels == 1)[:, None]
    X[:, 8:16] += 0.7 * (labels == 2)[:, None]
    return X, np.eye(3)[labels]


def _column_signs(a: np.ndarray, b: np.ndarray) -> np.ndarray:
    return np.where((a * b).sum(axis=0) < 0.0, -1.0, 1.0)


# ── Module contract ──────────────────────────────────────────────────────────


def test_flat_namespace_returns_the_submodule_objects():
    import motco.stats as stats
    from motco.stats import sparse_pls

    for name in ("SparsePLSDA", "SparsityDesign", "SPARSITY_SCOPES", "DEFAULT_KEEP_FRACTIONS"):
        assert getattr(stats, name) is getattr(sparse_pls, name)
        assert name in stats.__all__
    assert SPARSITY_SCOPES == ("joint", "block")
    assert DEFAULT_KEEP_FRACTIONS == (1.0, 0.5, 0.2, 0.1, 0.05)
    assert isinstance(SparsePLSDA(2, keep=3), PLSRegression)


def test_clone_keeps_the_budget():
    model = SparsePLSDA(3, keep=(4, 2), groups=(10, 5), scale=False)
    copy = clone(model)
    assert copy.get_params() == model.get_params()
    assert copy.get_params()["keep"] == (4, 2) and copy.get_params()["groups"] == (10, 5)


# ── Budget semantics ─────────────────────────────────────────────────────────


def test_joint_budget_keeps_fifty_of_five_hundred():
    X, Y = _three_class(500)
    model = SparsityDesign("joint", (0.1,)).resolve(500).estimator(3, 0.1).fit(X, Y)
    assert (model.support_.sum(axis=0) == 50).all()
    assert model.keep_.tolist() == [[50]] * 3
    assert model.support_sizes().tolist() == [[50]] * 3


def test_block_budget_keeps_each_block_share():
    X, Y = _three_class(367 + 131)
    design = SparsityDesign("block", (0.1,), (367, 131))
    assert design.counts(0.1) == (37, 13)
    model = design.estimator(3, 0.1).fit(X, Y)
    assert model.support_sizes().tolist() == [[37, 13]] * 3
    assert model.keep_.tolist() == [[37, 13]] * 3


def test_smallest_block_keeps_one_feature():
    X, Y = _three_class(4 + 120)
    design = SparsityDesign("block", (0.05,), (4, 120))
    assert design.counts(0.05) == (1, 6)
    model = design.estimator(2, 0.05).fit(X, Y)
    assert model.support_sizes().tolist() == [[1, 6]] * 2
    assert keep_count(0.001, 10) == 1


def test_counts_round_to_nearest_with_halves_up():
    assert keep_count(0.1, 367) == 37
    assert keep_count(0.1, 131) == 13
    assert keep_count(0.5, 5) == 3
    assert keep_count(1.0, 9) == 9


def test_soft_threshold_shrinks_survivors_by_the_largest_discarded_magnitude():
    weights = np.array([0.5, -2.0, 0.1, 3.0, -0.7, 0.2])
    out = soft_threshold(weights.copy(), (6,), (2,))
    # Kept 3.0 and −2.0; the largest discarded magnitude is 0.7.
    np.testing.assert_allclose(out, [0.0, -1.3, 0.0, 2.3, 0.0, 0.0])
    grouped = soft_threshold(weights.copy(), (3, 3), (1, 1))
    np.testing.assert_allclose(grouped, [0.0, -1.5, 0.0, 2.3, 0.0, 0.0])
    unchanged = soft_threshold(weights.copy(), (6,), (6,))
    np.testing.assert_array_equal(unchanged, weights)


def test_single_outcome_first_component_is_the_soft_thresholded_covariance():
    # PLS1 converges in one step: the first weight is soft(Xcᵀ yc), normalized.
    X, Y = _three_class(40)
    y = Y[:, 1]
    model = SparsePLSDA(1, keep=5, scale=False).fit(X, y)
    expected = soft_threshold((X - X.mean(axis=0)).T @ (y - y.mean()), (40,), (5,))
    expected /= np.linalg.norm(expected)
    np.testing.assert_allclose(model.x_weights_[:, 0] * _column_signs(model.x_weights_, expected[:, None])[0],
                               expected, atol=1e-12)


def test_fit_is_deterministic():
    X, Y = _three_class(200)
    design = SparsityDesign("block", (0.2,), (120, 80))
    first = design.estimator(3, 0.2).fit(X, Y)
    second = design.estimator(3, 0.2).fit(X, Y)
    for attribute in ("x_weights_", "x_loadings_", "y_loadings_", "x_rotations_", "x_scores_", "coef_"):
        np.testing.assert_array_equal(getattr(first, attribute), getattr(second, attribute))


# ── Rejections before fitting ────────────────────────────────────────────────


@pytest.mark.parametrize("fraction", [0, 0.0, 1.5, -0.1, float("nan"), True, "0.5"])
def test_fraction_out_of_range_is_rejected(fraction):
    with pytest.raises(ValueError, match=r"keep fraction"):
        SparsityDesign("joint", (1.0, fraction))


@pytest.mark.parametrize(
    ("kwargs", "message"),
    [
        ({"scope": "blocks"}, r"Unknown sparsity scope 'blocks'"),
        ({"scope": "joint", "keep_fractions": ()}, r"at least one fraction"),
        ({"scope": "joint", "keep_fractions": (0.5, 0.5)}, r"repeat a fraction"),
        ({"scope": "joint", "keep_fractions": "0.5"}, r"sequence of fractions"),
        ({"scope": "block"}, r"needs feature_groups"),
        ({"scope": "block", "feature_groups": (10, 0)}, r"group size 0"),
        ({"scope": "block", "feature_groups": (10, 2.5)}, r"group size 2.5"),
    ],
)
def test_malformed_design_is_rejected(kwargs, message):
    with pytest.raises(ValueError, match=message):
        SparsityDesign(**kwargs)


def test_design_ladder_is_held_fullest_first():
    design = SparsityDesign("joint", (0.1, 1.0, 0.5))
    assert design.keep_fractions == (1.0, 0.5, 0.1)
    assert design.to_dict() == {"scope": "joint", "keep_fractions": [1.0, 0.5, 0.1]}


def test_membership_must_cover_every_feature():
    with pytest.raises(ValueError, match=r"cover 15 features, but X has 16"):
        SparsityDesign("block", (0.5,), (10, 5)).resolve(16)


@pytest.mark.parametrize(
    ("keep", "groups", "message"),
    [
        (0, None, r"keep count 0"),
        (31, None, r"keep count 31"),
        ((3, 2), None, r"2 count\(s\) for 1 feature group"),
        (3, (20, 10), r"1 count\(s\) for 2 feature group"),
        ((3, 11), (20, 10), r"keep count 11"),
        ((3, 2), (20, 9), r"cover 29 features, but X has 30"),
        ((3, 2), (20, -10), r"group size -10"),
        (2.5, None, r"keep 2.5 must be a count"),
    ],
)
def test_malformed_budget_is_rejected_before_fitting(monkeypatch, keep, groups, message):
    X, Y = _three_class(30)

    def forbidden(*args, **kwargs):
        raise AssertionError("no component may be fitted when the budget is rejected")

    monkeypatch.setattr("motco.stats.sparse_pls._center_scale_xy", forbidden)
    monkeypatch.setattr(PLSRegression, "fit", forbidden)
    model = SparsePLSDA(2, keep=keep, groups=groups)
    with pytest.raises(ValueError, match=message):
        model.fit(X, Y)
    assert not hasattr(model, "x_weights_")


# ── Full budget is ordinary PLS ──────────────────────────────────────────────


@pytest.mark.parametrize("scale", [True, False])
@pytest.mark.parametrize("n_components", [2, 3, 4, 5])
def test_full_budget_equals_plsregression(n_components, scale):
    X, Y = _three_class(150)
    reference = PLSRegression(n_components=n_components, scale=scale, max_iter=1000).fit(X, Y)
    for model in (
        SparsePLSDA(n_components, keep=150, scale=scale),
        SparsityDesign("block", (1.0,), (100, 50)).estimator(n_components, 1.0, scale=scale),
    ):
        model.fit(X, Y)
        # The contract is 1e-6 up to column signs; the full budget delegates to
        # PLSRegression.fit, so the agreement is exact.
        for attribute in ("x_scores_", "x_rotations_", "x_weights_", "x_loadings_", "y_loadings_"):
            np.testing.assert_array_equal(getattr(model, attribute), getattr(reference, attribute))
        np.testing.assert_array_equal(model.predict(X), reference.predict(X))
        np.testing.assert_array_equal(model.transform(X), reference.transform(X))
        np.testing.assert_allclose(calculate_vips(model), calculate_vips(reference), rtol=1e-12)
        assert model.support_.all()


def test_sparse_fit_exposes_the_pls_surface():
    X, Y = _three_class(90)
    model = SparsityDesign("block", (0.2,), (60, 30)).estimator(3, 0.2).fit(X, Y)
    scores = model.transform(X)
    np.testing.assert_allclose(scores, model.x_scores_, atol=1e-10)
    np.testing.assert_allclose(scores, ((X - model._x_mean) / model._x_std) @ model.x_rotations_, atol=1e-12)
    assert model.inverse_transform(scores).shape == X.shape
    assert model.predict(X).shape == Y.shape
    # Rows of the rotations are zero exactly where no component loads on the feature.
    unused = ~model.support_.any(axis=1)
    assert unused.any()
    assert np.all(model.x_rotations_[unused] == 0.0)
    latent_map = LinearLatentMap.from_pls(model)
    projection = latent_map.basis @ latent_map.basis.T
    np.testing.assert_allclose(projection @ model.x_rotations_, model.x_rotations_, atol=1e-10)
    vips = calculate_vips(model)
    assert vips.shape == (90,) and np.all(vips[unused] == 0.0)


# ── Reference implementation ─────────────────────────────────────────────────


def test_reference_fixture_is_reproduced():
    fixture = np.load(FIXTURE)
    weights, scores = fixture["weights"], fixture["scores"]
    # mixOmics ran to a squared weight change below 1e-20; fit to the same fixed point.
    model = SparsePLSDA(int(fixture["ncomp"]), keep=int(fixture["keep"]), scale=False, tol=1e-24, max_iter=10000)
    model.fit(fixture["X"], fixture["Y"])
    np.testing.assert_array_equal(model.support_, weights != 0.0)
    signs = _column_signs(model.x_weights_, weights)
    np.testing.assert_allclose(model.x_weights_ * signs, weights, atol=1e-6)
    np.testing.assert_allclose(model.x_scores_ * signs, scores, atol=1e-6)
    assert str(fixture["mixomics_version"]) and len(str(fixture["script_sha256"])) == 64


def test_reference_fixture_records_the_script_that_made_it():
    import hashlib

    fixture = np.load(FIXTURE)
    script = Path(__file__).parent / "data" / "reference" / "splsda_fixture.R"
    assert str(fixture["script_sha256"]) == hashlib.sha256(script.read_bytes()).hexdigest()


def test_default_tolerance_keeps_the_reference_support():
    fixture = np.load(FIXTURE)
    model = SparsePLSDA(int(fixture["ncomp"]), keep=int(fixture["keep"]), scale=False).fit(
        fixture["X"], fixture["Y"]
    )
    np.testing.assert_array_equal(model.support_, fixture["weights"] != 0.0)
