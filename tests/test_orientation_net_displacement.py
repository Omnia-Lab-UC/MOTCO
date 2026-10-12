"""Contract for the opt-in net-displacement orientation estimator.

See ``openspec/changes/orientation-measurement-probe/`` (the
``trajectory-orientation-invariance`` delta). The PC1 estimator stays the
default; these tests pin the alternative and its plumbing through RRPP and the
realized-geometry checkpoints.
"""

from __future__ import annotations

import numpy as np
import pytest

from motco.simulations.diagnostics import geometry_from_means
from motco.stats.permutation import RRPP
from motco.stats.trajectory import _estimate_net_displacement, estimate_difference

ANGLE_TOL_DEG = 1e-6


def _estimate(left: np.ndarray, right: np.ndarray, estimator: str) -> tuple[np.ndarray, ...]:
    vectors = np.vstack([left, right])
    n_levels = left.shape[0]
    contrast = [list(range(n_levels)), list(range(n_levels, 2 * n_levels))]
    identity = np.eye(vectors.shape[0])
    return estimate_difference(vectors, identity, identity, contrast, orientation_estimator=estimator)


def _angle(left: np.ndarray, right: np.ndarray, estimator: str) -> float:
    return float(_estimate(left, right, estimator)[1][0, 1])


def _degrees(u: np.ndarray, v: np.ndarray) -> float:
    cosine = float(u @ v / (np.linalg.norm(u) * np.linalg.norm(v)))
    return float(np.degrees(np.arccos(np.clip(cosine, -1.0, 1.0))))


EQUILATERAL = np.array([[0.0, 0.0], [1.0, 0.0], [0.5, np.sqrt(3.0) / 2.0]]) + np.array([3.0, -2.0])


def test_net_displacement_is_the_first_to_last_direction():
    rng = np.random.default_rng(0)
    a = rng.normal(size=(4, 6))
    b = rng.normal(size=(4, 6))
    expected = _degrees(a[-1] - a[0], b[-1] - b[0])
    assert _angle(a, b, "net_displacement") == pytest.approx(expected, abs=ANGLE_TOL_DEG)


def test_two_stage_estimators_agree():
    rng = np.random.default_rng(1)
    for _ in range(20):
        a = rng.normal(size=(2, 5))
        b = rng.normal(size=(2, 5))
        assert _angle(a, b, "net_displacement") == pytest.approx(_angle(a, b, "pc1"), abs=1e-6)


def test_equilateral_configuration_has_zero_angle():
    assert _angle(EQUILATERAL, EQUILATERAL + 1.5, "net_displacement") == pytest.approx(0.0, abs=ANGLE_TOL_DEG)


def test_reversed_progression_is_180_degrees():
    assert _angle(EQUILATERAL, EQUILATERAL[::-1], "net_displacement") == pytest.approx(180.0, abs=ANGLE_TOL_DEG)


def test_only_angle_depends_on_the_estimator():
    rng = np.random.default_rng(2)
    a = rng.normal(size=(3, 7))
    b = rng.normal(size=(3, 7))
    pc1 = _estimate(a, b, "pc1")
    net = _estimate(a, b, "net_displacement")
    np.testing.assert_array_equal(pc1[0], net[0])
    np.testing.assert_array_equal(pc1[2], net[2])
    assert not np.allclose(pc1[1], net[1])


def test_default_is_pc1():
    rng = np.random.default_rng(3)
    vectors = rng.normal(size=(6, 4))
    contrast = [[0, 1, 2], [3, 4, 5]]
    identity = np.eye(6)
    default = estimate_difference(vectors, identity, identity, contrast)
    explicit = estimate_difference(vectors, identity, identity, contrast, orientation_estimator="pc1")
    for left, right in zip(default, explicit, strict=True):
        np.testing.assert_array_equal(left, right)


def test_unknown_estimator_is_rejected():
    with pytest.raises(ValueError, match="orientation_estimator"):
        _estimate(EQUILATERAL, EQUILATERAL, "pc2")


def test_zero_net_displacement_is_unavailable():
    closed = np.array([[0.0, 0.0], [1.0, 1.0], [0.0, 0.0]])
    assert np.all(np.isnan(_estimate_net_displacement(closed, [0, 1, 2])))
    assert np.isnan(_angle(closed, EQUILATERAL, "net_displacement"))
    # The checkpoint reports it as unavailable, never as a non-finite value.
    diagnostic = geometry_from_means(
        np.vstack([closed, EQUILATERAL]), ["A", "B"], ["0", "1", "2"], orientation_estimator="net_displacement"
    )
    assert diagnostic.angle is None
    assert diagnostic.availability["angle"] is False


def test_checkpoint_angle_follows_the_estimator():
    rng = np.random.default_rng(4)
    means = rng.normal(size=(6, 5))
    pc1 = geometry_from_means(means, ["A", "B"], ["0", "1", "2"])
    net = geometry_from_means(means, ["A", "B"], ["0", "1", "2"], orientation_estimator="net_displacement")
    assert net.angle == pytest.approx(_degrees(means[2] - means[0], means[5] - means[3]), abs=ANGLE_TOL_DEG)
    assert net.angle != pytest.approx(pc1.angle)
    assert (net.delta, net.shape) == (pc1.delta, pc1.shape)


def _rrpp_problem():
    import pandas as pd

    from motco.stats.design import build_ls_means, get_model_matrix

    rng = np.random.default_rng(5)
    frame = pd.DataFrame(
        {"group": np.repeat(["A", "B"], 24), "stage": np.tile(np.repeat(["0", "1", "2"], 8), 2)}
    )
    Y = rng.normal(size=(48, 4)) + np.where(frame["stage"].astype(int).to_numpy()[:, None] > 0, 0.5, 0.0)
    full = get_model_matrix(frame, group_col="group", level_col="stage", full=True)
    reduced = get_model_matrix(frame, group_col="group", level_col="stage", full=False)
    ls_means = build_ls_means(["A", "B"], ["0", "1", "2"], full=True)
    return Y, full, reduced, ls_means, [[0, 1, 2], [3, 4, 5]]


@pytest.mark.parametrize("n_jobs", [1, 2])
def test_rrpp_uses_the_selected_estimator(n_jobs: int):
    Y, full, reduced, ls_means, contrast = _rrpp_problem()
    common = dict(permutations=12, n_jobs=n_jobs, progress=False, seed=9)
    _, pc1_angles, pc1_shapes = RRPP(Y, full, reduced, ls_means, contrast, **common)
    net_deltas, net_angles, net_shapes = RRPP(
        Y, full, reduced, ls_means, contrast, orientation_estimator="net_displacement", **common
    )
    # Recompute every permuted angle directly with the net-displacement estimator.
    from motco.stats.trajectory import estimate_betas

    y_hat = np.asarray(reduced) @ np.asarray(estimate_betas(reduced, Y))
    residuals = Y - y_hat
    if n_jobs == 1:
        rng = np.random.default_rng(9)
        for observed in net_angles:
            idx = rng.permutation(Y.shape[0])
            _, expected, _ = estimate_difference(
                y_hat + residuals[idx], full, ls_means, contrast, orientation_estimator="net_displacement"
            )
            np.testing.assert_allclose(observed, expected)
    for left, right in zip(pc1_shapes, net_shapes, strict=True):
        np.testing.assert_array_equal(left, right)
    assert not all(np.allclose(a, b) for a, b in zip(pc1_angles, net_angles, strict=True))
    assert len(net_deltas) == 12
