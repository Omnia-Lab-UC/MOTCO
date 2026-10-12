from __future__ import annotations

import logging
from typing import Any, Sequence, Union

import numpy as np
import pandas as pd

from motco.stats.design import _sort_levels, build_ls_means, get_model_matrix

logger = logging.getLogger(__name__)

#: Version of the configuration-spectrum contract produced by
#: :func:`configuration_spectra`. Bumped whenever the recorded fields or their
#: definitions change, so a persisted record can never be read under the wrong
#: contract.
CONFIG_SPECTRUM_VERSION = 1

#: Version of the trajectory ``shape`` statistic contract implemented by
#: :func:`_estimate_shape`. Bumped whenever the statistic's definition changes,
#: so records produced under different contracts can never be mixed (the study
#: parameter signature carries this key; see ``simulations/grid.py``).
#: 1 = proper-rotation alignment (reflections retained where the configuration
#: spanned the ambient space, aligned away otherwise);
#: 2 = full orthogonal alignment at every ambient dimension.
SHAPE_STATISTIC_VERSION = 2

#: Orientation estimators :func:`estimate_difference` can measure ``angle`` with.
#: ``"pc1"`` (the default) is the signed leading principal axis of the stage
#: configuration; ``"net_displacement"`` is the unit first-to-last stage vector.
ORIENTATION_ESTIMATORS: tuple[str, ...] = ("pc1", "net_displacement")

#: Net displacements shorter than this have no defined direction.
_NET_DISPLACEMENT_TOLERANCE = 1e-15


def estimate_betas(
    X: Union[pd.DataFrame, np.ndarray], Y: Union[pd.DataFrame, np.ndarray]
) -> Union[pd.DataFrame, np.ndarray]:
    """
    Estimate the beta coefficients between an outcome matrix
    and a model matrix

    Parameters
    ----------
    X: Union[pd.DataFrame, np.ndarray]
        Model matrix with intercept.
    Y: Union[pd.DataFrame, np.ndarray]
        Outcome matrix.

    Returns
    -------
    betas: Union[pd.DataFrame, np.ndarray]
        Beta coefficients
    """
    # Convert inputs to arrays for linear algebra while preserving Y's metadata
    X_arr = np.asarray(X, dtype=float)
    if isinstance(Y, pd.DataFrame):
        Y_arr = Y.to_numpy(dtype=float)
        y_cols = Y.columns
    else:
        Y_arr = np.asarray(Y, dtype=float)
        y_cols = None

    # Solve normal equations using factorization with robust fallbacks
    XtX = X_arr.T @ X_arr
    XtY = X_arr.T @ Y_arr
    try:
        # Cholesky is fastest and most stable for SPD XtX
        L = np.linalg.cholesky(XtX)
        tmp = np.linalg.solve(L, XtY)
        betas_arr = np.linalg.solve(L.T, tmp)
    except np.linalg.LinAlgError:
        logger.warning("Cholesky decomposition failed; falling back to direct solve. Check for near-singular XtX.")
        try:
            # Fall back to a direct solve of the normal equations
            betas_arr = np.linalg.solve(XtX, XtY)
        except np.linalg.LinAlgError:
            logger.warning("Direct solve failed; falling back to lstsq. Model matrix may be rank-deficient.")
            # Final fallback: least-squares without forming normal equations
            # This handles rank deficiency and ill-conditioning better.
            betas_arr, *_ = np.linalg.lstsq(X_arr, Y_arr, rcond=None)

    if y_cols is not None:
        # Return a DataFrame so downstream matmul with numpy yields a DataFrame
        # and index/column handling stays consistent with previous behavior.
        return pd.DataFrame(betas_arr, columns=y_cols)
    return betas_arr


def pair_difference(
    dat: pd.DataFrame,
    group_col: str,
    level_col: str,
    groups: tuple[str, str] | None = None,
    levels: tuple[str, str] | None = None,
    feature_cols: Sequence[str] | None = None,
) -> tuple[float, float]:
    """
    Estimate difference in direction (angle, degrees) and magnitude (delta)
    between two groups across two levels.

    The change vector for a group is defined as `level1 - level2` over the
    selected feature columns.

    Parameters
    ----------
    dat: pd.DataFrame
        DataFrame containing features plus `group_col` and `level_col`.
    group_col: str
        Column with groups (between-subject factor).
    level_col: str
        Column with levels/states (within-subject factor).
    groups: tuple[str, str] | None
        Pair of group labels to compare. If None, infer and require exactly two.
    levels: tuple[str, str] | None
        Pair of level labels to use for the change vector. If None, infer and
        require exactly two.
    feature_cols: Sequence[str] | None
        Feature columns to use. If None, all numeric columns except `group_col`
        and `level_col` are used.

    Returns
    -------
    tuple[float, float]
        (angle_degrees, delta_magnitude_difference)

    Notes
    -----
    See [1]_ for more information on two-state comparisons.

    References
    ----------
    .. [1] Collyer, Michael L., and Dean C. Adams. "Analysis of two‐state
           multivariate phenotypic change in ecological studies." Ecology 88.3
           (2007): 683-692. https://doi.org/10.1890/06-0727
    """
    if feature_cols is None:
        feature_cols = [
            c
            for c in dat.select_dtypes(include=[np.number]).columns.tolist()
            if c not in {group_col, level_col}
        ]
    if not feature_cols:
        raise ValueError("No feature columns provided or detected.")

    g_vals = _sort_levels(pd.unique(dat[group_col].astype(str)).tolist())
    l_vals = _sort_levels(pd.unique(dat[level_col].astype(str)).tolist())
    if groups is None:
        if len(g_vals) != 2:
            raise ValueError(
                f"Expected exactly 2 groups, found {len(g_vals)}: {g_vals}"
            )
        groups = (g_vals[0], g_vals[1])
    if levels is None:
        if len(l_vals) != 2:
            raise ValueError(
                f"Expected exactly 2 levels, found {len(l_vals)}: {l_vals}"
            )
        levels = (l_vals[0], l_vals[1])

    # Compute per (group, level) means
    means = (
        dat.assign(
            __g=dat[group_col].astype(str), __l=dat[level_col].astype(str)
        )
        .groupby(["__g", "__l"])[feature_cols]
        .mean()
    )
    try:
        y_g1 = np.asarray(means.loc[(groups[0], levels[0])], dtype=float) - np.asarray(
            means.loc[(groups[0], levels[1])], dtype=float
        )
        y_g2 = np.asarray(means.loc[(groups[1], levels[0])], dtype=float) - np.asarray(
            means.loc[(groups[1], levels[1])], dtype=float
        )
    except KeyError as e:
        raise ValueError(
            "Missing combinations for the requested groups/levels in the data."
        ) from e

    d1 = float(np.linalg.norm(y_g1))
    d2 = float(np.linalg.norm(y_g2))
    delta = abs(d1 - d2)
    if d1 == 0 or d2 == 0:
        raise ValueError(
            "Zero-magnitude change vector for at least one group; angle is undefined."
        )
    angle = float(np.degrees(np.arccos(np.inner(y_g1 / d1, y_g2 / d2))))
    return angle, delta


def configuration_spectrum(configuration: Union[pd.DataFrame, np.ndarray]) -> dict[str, Any]:
    """Describe the eigenspectrum of one centered stage configuration.

    The configuration is a ``k x d`` matrix of stage points in the measurement
    space. It is column-centered, and the squared singular values of the centered
    matrix — the eigenvalues of its scatter matrix — are normalized by their sum.

    Returns
    -------
    dict[str, Any]
        ``n_points`` / ``n_dimensions`` of the configuration, its
        ``total_variance`` (the untransformed eigenvalue sum), the normalized
        ``spectrum``, and the **relative eigengap** ``(l1 - l2) / sum(l)``.

    Notes
    -----
    The relative eigengap measures how strongly a single axis dominates the
    configuration: it is the observable that predicts how well an orientation
    (PC1) can be resolved from noise, since PC1's estimator variance scales
    like noise over the eigengap. See ``docs/reports/geometry-audit-2026-09-01.md``.

    Two degeneracies are recorded rather than hidden:

    - A configuration of ``k = 2`` stages has rank at most 1 after centering, so
      its relative eigengap is identically ``1.0``. It is recorded as computed
      and is uninformative by construction.
    - A configuration with zero total variance (every stage at one point, or
      fewer than two stages) has no defined spectrum. It records an empty
      ``spectrum`` and ``relative_eigengap = None`` — never a non-finite float,
      which JSON cannot represent conformingly.
    """

    X_raw = np.asarray(configuration, dtype=float)
    if X_raw.ndim != 2:
        raise ValueError(f"configuration must be a 2-D matrix; got shape {X_raw.shape}.")
    n_points, n_dimensions = X_raw.shape
    entry: dict[str, Any] = {
        "n_points": int(n_points),
        "n_dimensions": int(n_dimensions),
    }
    if n_points == 0 or n_dimensions == 0:
        entry.update({"total_variance": 0.0, "spectrum": [], "relative_eigengap": None})
        return entry

    X = X_raw - X_raw.mean(axis=0, keepdims=True)
    eigenvalues = np.linalg.svd(X, compute_uv=False) ** 2
    total = float(eigenvalues.sum())
    if not np.isfinite(total) or total <= 0.0:
        entry.update({"total_variance": 0.0, "spectrum": [], "relative_eigengap": None})
        return entry

    normalized = eigenvalues / total
    second = float(normalized[1]) if normalized.size > 1 else 0.0
    entry.update(
        {
            "total_variance": total,
            "spectrum": [float(value) for value in normalized],
            "relative_eigengap": float(normalized[0] - second),
        }
    )
    return entry


def configuration_spectra(
    obs_vect: Union[pd.DataFrame, np.ndarray],
    contrast: list[list[int]],
) -> dict[str, Any]:
    """Spectra of the pooled and per-group stage-mean configurations.

    ``obs_vect`` is the LS-mean matrix the trajectory statistics are measured
    from (``LS_means @ betas``); ``contrast`` enumerates, per group, the rows
    belonging to that group's trajectory in stage order.

    The **per-group** configurations are those rows; the **pooled**
    configuration averages the groups stage by stage. Pooling requires every
    group to carry the same number of stages — the design the study builds — and
    is reported as ``None`` when they do not, rather than silently pooling
    mismatched stages.

    Because the LS-mean rows exist per permutation, this description is
    computable inside RRPP without a second beta solve. Under a balanced design
    the pooled LS-mean configuration coincides with the pooled sample stage
    means.
    """

    V = np.asarray(obs_vect, dtype=float)
    groups = [V[np.asarray(levels, dtype=int), :] for levels in contrast]
    sizes = {group.shape[0] for group in groups}
    pooled = (
        configuration_spectrum(np.mean(np.stack(groups, axis=0), axis=0))
        if groups and len(sizes) == 1
        else None
    )
    return {
        "version": CONFIG_SPECTRUM_VERSION,
        "pooled": pooled,
        "groups": [configuration_spectrum(group) for group in groups],
    }


def pooled_relative_eigengap(spectra: dict[str, Any] | None) -> float | None:
    """Read the pooled relative eigengap out of a spectrum block, or ``None``."""

    pooled = (spectra or {}).get("pooled")
    if not pooled:
        return None
    value = pooled.get("relative_eigengap")
    if value is None:
        return None
    value = float(value)
    return value if np.isfinite(value) else None


def principal_orientation(
    configuration: Union[pd.DataFrame, np.ndarray],
    levels: Sequence[int] | None = None,
) -> np.ndarray:
    """Signed leading principal axis of one stage configuration.

    This is the public entry point to the orientation functional the tested
    ``angle`` statistic is built from: PC1 of the column-centered configuration,
    signed so it points along the net displacement from the first stage row to
    the last. It is exposed so that any consumer explaining an ``angle`` result
    — attribution in particular — decomposes *the same* quantity the test
    measured, rather than a second implementation of the convention.

    Parameters
    ----------
    configuration:
        ``k x d`` matrix of stage points, in stage order, in the measurement
        space. Rows are the trajectory's stage means.
    levels:
        Optional row indices selecting (and ordering) the stages to use. When
        omitted every row is used in the order given.

    Returns
    -------
    np.ndarray
        Unit vector of length ``d`` along the trajectory's direction of
        progression.

    Notes
    -----
    See :func:`_estimate_orientation` for the sign convention, its deliberate
    deviation from the reference supplement, and the one accepted degeneracy —
    a closed trajectory whose net displacement vanishes, for which no direction
    is defined and the returned sign is arbitrary. Callers that must not report
    an arbitrary sign SHOULD test the net displacement against a tolerance
    before calling.
    """

    values = np.asarray(configuration, dtype=float)
    if values.ndim != 2:
        raise ValueError(f"configuration must be a 2-D matrix; got shape {values.shape}.")
    indices = list(range(values.shape[0])) if levels is None else [int(level) for level in levels]
    if len(indices) < 2:
        raise ValueError("configuration must contain at least two stages.")
    return _estimate_orientation(values, indices)


def estimate_difference(
    Y: Union[pd.DataFrame, np.ndarray],
    model_matrix: Union[pd.DataFrame, np.ndarray],
    LS_means: Union[pd.DataFrame, np.ndarray],
    contrast: list[list[int]],
    *,
    return_spectra: bool = False,
    orientation_estimator: str = "pc1",
) -> tuple:
    """
    Estimate parameters angle, delta, and shape given an outcome
    matrix, model matrix, and contrast to compare. This is a comparison
    of more than two states.

    Parameters
    ----------
    Y: Union[pd.DataFrame, np.ndarray]
        Outcome matrix.
    model_matrix: Union[pd.DataFrame, np.ndarray]
        Model matrix with intercept.
    LS_means: Union[pd.DataFrame, np.ndarray]
        Least-squares means to estimate.
    contrast: list[list[int]]
        Indices indicating the groups to compare based on LS means.
        Each list must contain the cohorts that belong to the same group.
    return_spectra: bool
        When ``True``, additionally return the pooled and per-group stage-mean
        configuration spectra (:func:`configuration_spectra`) computed from the
        LS-mean vectors this call already fits. Off by default: the returned
        tuple is then exactly the historical three-element one.
    orientation_estimator: str
        How each group's orientation is estimated for ``angles``: ``"pc1"``
        (default; :func:`_estimate_orientation`) or ``"net_displacement"``
        (:func:`_estimate_net_displacement`). Only ``angles`` depends on it. A
        net displacement of zero length has no direction, and its pairwise
        angles are ``NaN``.

    Returns
    -------
    deltas: np.ndarray
        Symmetric matrix (n_groups x n_groups) with differences in magnitude.
    angles: np.ndarray
        Symmetric matrix (n_groups x n_groups) with differences in direction (degrees).
    shapes: np.ndarray
        Symmetric matrix (n_groups x n_groups) with Procrustes shape distances
        after removing translation, uniform scale, and any orthogonal
        transformation. Reflections are aligned away at every ambient
        dimension; a mirror difference is an orientation difference and
        surfaces in ``angles``.
    spectra: dict[str, Any]
        Only when ``return_spectra=True``. A recorded covariate — it enters no
        statistic and no decision rule.

    Notes
    -----
    See [1]_ for more information on trajectory analysis.

    References
    ----------
    .. [1] Adams, Dean C., and Michael L. Collyer.
           "A general framework for the analysis of phenotypic trajectories in
           evolutionary studies."
           Evolution: International Journal of Organic Evolution 63.5 (2009):
           1143-1154.
           https://doi.org/10.1111/j.1558-5646.2009.00649.x
    """
    # --- Input validation ---
    _Y = np.asarray(Y, dtype=float)
    _X = np.asarray(model_matrix, dtype=float)
    _LS = np.asarray(LS_means, dtype=float)
    if _Y.shape[0] != _X.shape[0]:
        raise ValueError(
            f"Y has {_Y.shape[0]} rows but model_matrix has {_X.shape[0]} rows — "
            "number of rows must match."
        )
    if _LS.shape[1] != _X.shape[1]:
        raise ValueError(
            f"LS_means has {_LS.shape[1]} columns but model_matrix has {_X.shape[1]} columns — "
            "number of columns must match."
        )
    _n_ls = _LS.shape[0]
    for _gi, _group in enumerate(contrast):
        for _idx in _group:
            if not (0 <= _idx < _n_ls):
                raise ValueError(
                    f"contrast[{_gi}] contains index {_idx}, but LS_means only has {_n_ls} rows "
                    f"(valid indices: 0–{_n_ls - 1})."
                )
    if not np.all(np.isfinite(_Y)):
        raise ValueError("Y contains NaN or Inf values.")
    if not np.all(np.isfinite(_X)):
        raise ValueError("model_matrix contains NaN or Inf values.")
    if orientation_estimator not in ORIENTATION_ESTIMATORS:
        raise ValueError(
            f"Unknown orientation_estimator {orientation_estimator!r}; "
            f"expected one of {list(ORIENTATION_ESTIMATORS)}."
        )
    # --- End validation ---
    estimate_orientation = (
        _estimate_net_displacement if orientation_estimator == "net_displacement" else _estimate_orientation
    )
    n_groups = len(contrast)
    betas = estimate_betas(model_matrix, Y)
    # Compute LS-mean vectors; keep as DataFrame to minimize behavioral drift
    obs_vect = pd.DataFrame(
        np.matmul(np.asarray(LS_means, dtype=float), np.asarray(betas, dtype=float))
    )
    ys = []
    des = []
    angles = np.zeros((n_groups, n_groups))
    deltas = np.zeros((n_groups, n_groups))
    for i in range(n_groups):
        y = estimate_orientation(obs_vect, contrast[i])
        d = _estimate_size(obs_vect, contrast[i])
        des.append(d)
        ys.append(y)
    shapes = _estimate_shape(obs_vect, contrast)
    for i in range(n_groups):
        comp = i + 1
        while comp < n_groups:
            delta = np.abs(des[i] - des[comp])
            # When using SVD, no need to divide by size
            dot = np.clip(np.inner(ys[i], ys[comp]), -1.0, 1.0)
            angle = np.arccos(dot) * 180 / np.pi
            deltas[i, comp] = delta
            deltas[comp, i] = delta
            angles[i, comp] = angle
            angles[comp, i] = angle
            comp += 1
    if return_spectra:
        return deltas, angles, shapes, configuration_spectra(obs_vect, contrast)
    return deltas, angles, shapes


def get_observed_vectors(
    X: pd.DataFrame,
    Y: Union[pd.DataFrame, np.ndarray],
    group_col: str,
    level_col: str,
    full: bool = True,
) -> pd.DataFrame:
    """
    Get LS-mean vectors for each group × level cell.

    Parameters
    ----------
    X: pd.DataFrame
        DataFrame containing factors `group_col` and `level_col` for each row
        corresponding to `Y`.
    Y: Union[pd.DataFrame, np.ndarray]
        Outcome matrix (n_samples × n_features).
    group_col: str
        Group column name in `X`.
    level_col: str
        Level/state column name in `X`.
    full: bool
        Whether to include interactions in the model.

    Returns
    -------
    pd.DataFrame
        LS means arranged with a MultiIndex (group, level). Columns follow `Y`.
    """
    model_full = get_model_matrix(X[[group_col, level_col]], group_col, level_col, full)
    betas = estimate_betas(model_full, Y)

    g_levels = _sort_levels(pd.unique(X[group_col].astype(str)).tolist())
    l_levels = _sort_levels(pd.unique(X[level_col].astype(str)).tolist())
    ls_matrix = build_ls_means(g_levels, l_levels, full)
    means = np.matmul(ls_matrix, np.asarray(betas, dtype=float))

    # Build a clear index and columns
    idx = pd.MultiIndex.from_product([g_levels, l_levels], names=[group_col, level_col])
    if isinstance(Y, pd.DataFrame):
        cols = Y.columns
    else:
        cols = [f"f{i}" for i in range(means.shape[1])]  # type: ignore[assignment]
    return pd.DataFrame(means, index=idx, columns=cols)


def _estimate_size(obs_vect: pd.DataFrame | np.ndarray, levels: list[int]) -> float:
    """
    Estimate the size of a trajectory of two or more levels.

    Parameters
    ----------
    obs_vect: pd.DataFrame
        Matrix of observed mean vectors.
    levels: list[int]
        List of indices indicating the levels to consider.

    Returns
    -------
    size: float
        Size of the trajectory.
    """
    # Use a fully vectorized NumPy implementation to avoid Python-loop overhead
    X = np.asarray(obs_vect, dtype=float)[levels, :]
    if X.shape[0] < 2:
        return 0.0
    diffs = X[:-1] - X[1:]
    size = np.linalg.norm(diffs, axis=1).sum()
    return float(size)


def _estimate_orientation(
    obs_vect: pd.DataFrame | np.ndarray,
    levels: list[int],
) -> np.ndarray:
    """
    Estimate the orientation of a trajectory of two or more levels.

    Parameters
    ----------
    obs_vect: pd.DataFrame
        Matrix of observed mean vectors.
    levels: list[int]
        List of indices indicating the levels to consider.

    Returns
    -------
    orientation: np.ndarray
        Unit vector along the trajectory's direction of progression.

    Notes
    -----
    The estimator is PC1 of the centered stage configuration, as in the
    reference implementation (``tests/data/reference/evo_649_sm_suppmat.r:57``).

    The **sign convention deviates deliberately** from that reference. Its
    ``#check startingpoint location`` line (``:64``) resolves PC1's inherent
    sign ambiguity by anchoring on the raw first-stage row, ``M[1, ] . PC1``.
    That anchor is unusable here for two reasons:

    - Taken literally on raw coordinates, the sign depends on where the
      trajectory sits relative to the coordinate origin, so a pure translation
      can reverse it — which would break MOTCO's translation null control.
    - Taken on the centered configuration, the anchor is
      ``PC1 . (stage_first - centroid)``, which vanishes whenever a trajectory
      departs and returns laterally along its own principal axis. The sign is
      then decided by noise, and two trajectories carrying no orientation
      difference get reported as near-antiparallel.

    Anchoring on **net displacement** instead states what the reference line is
    groping toward — orient PC1 along the direction of progression. It is
    intrinsic to the trajectory, hence invariant to translation and to uniform
    scale; it stays well away from zero for any real progression; and for two
    stages it reduces exactly to the transition direction.

    Fixture effect: ``results_example2.csv`` (5 levels) is reproduced exactly by
    all three anchors. On ``results_example1.csv`` (2 levels) only the raw
    anchor matches R — pairs ``t1/t3`` and ``t2/t3`` sit on opposite sides of
    the PCA origin, so R reports 105.30/103.51 where both the previous centered
    anchor and this one report the supplements 74.70/76.49. This convention
    therefore leaves the committed fixture outputs unchanged; it does not
    introduce the example1 deviation, which predates it. See
    ``tests/test_trajectory_orientation.py`` for the invariance contract.

    The one accepted degeneracy is a closed trajectory whose last stage returns
    to its first: net displacement vanishes and no direction is defined.
    """
    X_raw = np.asarray(obs_vect, dtype=float)[levels, :]
    # Center rows
    X = X_raw - X_raw.mean(axis=0, keepdims=True)
    # The leading eigenvector of X.T @ X is the right singular vector of X
    # for the largest singular value; svd on the (k x p) data matrix avoids
    # ever forming the (p x p) covariance matrix.
    _, _, Vt = np.linalg.svd(X, full_matrices=False)
    orientation = Vt[0, :]
    # Sign PC1 along the trajectory's direction of progression (see Notes).
    c1 = float(orientation @ (X_raw[-1, :] - X_raw[0, :]))
    if c1 < 0:
        orientation = -orientation
    return orientation


def _estimate_net_displacement(
    obs_vect: pd.DataFrame | np.ndarray,
    levels: list[int],
) -> np.ndarray:
    """
    Estimate the orientation of a trajectory as its unit net displacement.

    Parameters
    ----------
    obs_vect: pd.DataFrame
        Matrix of observed mean vectors.
    levels: list[int]
        List of indices indicating the levels to consider, in stage order.

    Returns
    -------
    orientation: np.ndarray
        Unit vector from the first stage's mean to the last stage's mean, or a
        vector of ``NaN`` when that displacement has zero length.

    Notes
    -----
    For two stages this equals the signed PC1 of :func:`_estimate_orientation`.
    Unlike PC1 it stays well defined when the stage configuration is
    near-isotropic, where PC1 can swap with PC2 or flip sign; it ignores the
    interior stages.
    """
    X = np.asarray(obs_vect, dtype=float)[levels, :]
    displacement = X[-1, :] - X[0, :]
    length = float(np.linalg.norm(displacement))
    if not np.isfinite(length) or length <= _NET_DISPLACEMENT_TOLERANCE:
        return np.full(displacement.shape, np.nan)
    return displacement / length


def _estimate_shape(
    vectors: Union[pd.DataFrame, np.ndarray], contrast: list[list[int]]
) -> np.ndarray:
    """
    Estimate pairwise trajectory shape distances after Procrustes alignment.

    Each group trajectory is centered and scaled to unit centroid size, then
    each pair is aligned over the full orthogonal group: reflections are
    aligned away, at every ambient dimension, so a mirror pair has zero shape
    distance and the statistic never depends on whether the stage
    configuration spans the ambient space. A trajectory of ``k`` stages has
    rank at most ``k - 1``, so at every pre-integration checkpoint the
    configuration is rank-deficient and a proper-rotation constraint would be
    vacuous anyway (it could be satisfied for free in the null space); this
    policy makes the statistic identical across checkpoints by construction.
    A genuine mirror difference is an orientation difference and surfaces in
    ``angle`` instead. See ``SHAPE_STATISTIC_VERSION``.

    Parameters
    ----------
    vectors: Union[pd.DataFrame, np.ndarray]
        A n by k point matrix with the vectors to align,
        where n is the number of points, and k the number of dimensions.
    contrast: list[list[int]]
        Indices indicating the groups to compare based on LS means.
        Each list within the list must contain the cohorts that belong
        to the same group.

    Returns
    -------
    shape_distance: np.ndarray
        Matrix with shape distances.
    """
    V = np.asarray(vectors, dtype=float)
    n_groups = len(contrast)
    shapes = [_center_scale_unit(V[np.asarray(levels, dtype=int), :]) for levels in contrast]
    distances = np.zeros((n_groups, n_groups), dtype=float)

    for i in range(n_groups):
        for j in range(i + 1, n_groups):
            distance = _proper_procrustes_distance(shapes[i], shapes[j])
            distances[i, j] = distance
            distances[j, i] = distance

    return distances


def _center_scale_unit(A: np.ndarray) -> np.ndarray:
    """Center a trajectory and scale it to unit centroid size."""
    Z = A - A.mean(axis=0, keepdims=True)
    cs = float(np.linalg.norm(Z))
    if not np.isfinite(cs) or cs <= 1e-15:
        cs = 1.0
    return Z / cs


def _proper_procrustes_distance(reference: np.ndarray, target: np.ndarray) -> float:
    """
    Return the residual norm after aligning ``target`` to ``reference``.

    The alignment optimizes over the full orthogonal group ``O(k)``, not the
    rotation subgroup ``SO(k)``: reflections are aligned away rather than
    retained. This is the single reflection policy of the ``shape`` statistic
    and it holds at every ambient dimension, so the returned distance does not
    depend on whether the configuration spans the ambient space. See
    :func:`_estimate_shape` and ``SHAPE_STATISTIC_VERSION``.
    """
    H = target.T @ reference
    U, _, Vt = np.linalg.svd(H, full_matrices=False)
    R = U @ Vt
    return float(np.linalg.norm(reference - target @ R))
