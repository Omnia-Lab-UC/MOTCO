"""Sparse stage-PLS-DA: PLS-DA whose components load on a budget of features.

This module provides:
- `SparsePLSDA`: a :class:`sklearn.cross_decomposition.PLSRegression` whose X
  weight vector is soft-thresholded to a budget of features per component (the
  sPLS-DA of Lê Cao et al. 2008, 2011, in regression mode), with one budget over
  all features or one per feature group (block). It inherits every fitted
  attribute and method of ``PLSRegression`` (``transform``,
  ``inverse_transform``, ``predict``, ``x_rotations_`` …), so every consumer of
  a fitted PLS-DA model accepts it unchanged. At the full budget it *is*
  ``PLSRegression.fit``.
- `SparsityDesign`: the budget a cross-validation selects over — a scope
  (``joint`` or ``block``), a ladder of keep fractions and the block sizes.
- `SPARSITY_SCOPES`, `DEFAULT_KEEP_FRACTIONS`.

The sparse fit reuses scikit-learn's PLS machinery (input validation,
``_center_scale_xy``, the sign flip, the rotations and coefficients) and
changes one step: inside the NIPALS power iteration the X weight vector is
soft-thresholded before it is normalized. The rest of the iteration — its
initialization from the first non-constant Y column, its updates and its
convergence rule — is scikit-learn's, so a full budget, where the threshold is
zero, gives the ordinary fit bit for bit.
"""

from __future__ import annotations

import math
import warnings
from dataclasses import dataclass
from numbers import Integral, Real
from typing import Any, Sequence

import numpy as np
from scipy.linalg import pinv
from sklearn.base import _fit_context
from sklearn.cross_decomposition import PLSRegression
from sklearn.cross_decomposition._pls import _center_scale_xy, _svd_flip_1d
from sklearn.exceptions import ConvergenceWarning
from sklearn.utils import check_array, check_consistent_length
from sklearn.utils.validation import validate_data

#: Where a sparsity budget is counted: over all features, or within each block.
SPARSITY_SCOPES: tuple[str, ...] = ("joint", "block")

#: The default ladder of keep fractions: production (1.0) down to one feature in twenty.
DEFAULT_KEEP_FRACTIONS: tuple[float, ...] = (1.0, 0.5, 0.2, 0.1, 0.05)


def keep_count(fraction: float, n_features: int) -> int:
    """The number of features a keep fraction retains: ``max(1, round(f·p))``, halves rounded up."""

    return max(1, int(math.floor(float(fraction) * int(n_features) + 0.5)))


def _validate_fraction(fraction: Any, *, where: str) -> float:
    if isinstance(fraction, bool) or not isinstance(fraction, Real):
        raise ValueError(f"{where}: keep fraction {fraction!r} is not a number.")
    value = float(fraction)
    if not (math.isfinite(value) and 0.0 < value <= 1.0):
        raise ValueError(f"{where}: keep fraction {fraction!r} is outside (0, 1].")
    return value


def _validate_groups(groups: Any, *, where: str) -> tuple[int, ...]:
    try:
        values = tuple(groups)
    except TypeError:
        raise ValueError(f"{where}: feature groups {groups!r} must be a sequence of block sizes.") from None
    if not values:
        raise ValueError(f"{where}: feature groups must name at least one block.")
    for size in values:
        if isinstance(size, bool) or not isinstance(size, Integral) or size < 1:
            raise ValueError(f"{where}: feature group size {size!r} must be a positive integer.")
    return tuple(int(size) for size in values)


@dataclass(frozen=True)
class SparsityDesign:
    """A sparsity budget to select over: scope, keep-fraction ladder and block sizes.

    ``scope`` is ``"joint"`` (one budget over all features) or ``"block"``
    (one budget per block, so every block keeps its share whatever its size).
    ``keep_fractions`` are fractions of the scope's feature count, distinct and
    in (0, 1]; they are held fullest first. ``feature_groups`` are the block
    sizes in column order; the block scope requires them, the joint scope uses
    them only to report support per block. A joint design without groups is
    completed by :meth:`resolve` once the feature count is known.
    """

    scope: str
    keep_fractions: tuple[float, ...] = DEFAULT_KEEP_FRACTIONS
    feature_groups: tuple[int, ...] | None = None

    def __post_init__(self) -> None:
        if self.scope not in SPARSITY_SCOPES:
            raise ValueError(f"Unknown sparsity scope {self.scope!r}; allowed: {list(SPARSITY_SCOPES)}.")
        if isinstance(self.keep_fractions, (str, bytes)):
            raise ValueError(f"keep_fractions {self.keep_fractions!r} must be a sequence of fractions.")
        try:
            raw = tuple(self.keep_fractions)
        except TypeError:
            raise ValueError(f"keep_fractions {self.keep_fractions!r} must be a sequence of fractions.") from None
        if not raw:
            raise ValueError("keep_fractions must hold at least one fraction.")
        fractions = tuple(_validate_fraction(value, where="SparsityDesign") for value in raw)
        if len(set(fractions)) != len(fractions):
            raise ValueError(f"keep_fractions {list(raw)} repeat a fraction; they must be distinct.")
        object.__setattr__(self, "keep_fractions", tuple(sorted(fractions, reverse=True)))
        if self.feature_groups is None:
            if self.scope == "block":
                raise ValueError("The block sparsity scope needs feature_groups (the block sizes in column order).")
        else:
            object.__setattr__(self, "feature_groups", _validate_groups(self.feature_groups, where="SparsityDesign"))

    @property
    def n_features(self) -> int | None:
        """The total feature count the groups cover (``None`` while unresolved)."""

        return None if self.feature_groups is None else sum(self.feature_groups)

    def resolve(self, n_features: int) -> SparsityDesign:
        """This design checked against ``n_features`` columns (a joint design without groups gets one group)."""

        if self.feature_groups is None:
            return SparsityDesign(self.scope, self.keep_fractions, (int(n_features),))
        if sum(self.feature_groups) != int(n_features):
            raise ValueError(
                f"The sparsity feature groups {list(self.feature_groups)} cover {sum(self.feature_groups)} "
                f"features, but X has {n_features}."
            )
        return self

    @property
    def budget_groups(self) -> tuple[int, ...]:
        """The groups a budget is counted in: the whole feature set (joint) or every block (block)."""

        if self.feature_groups is None:
            raise ValueError("The sparsity design is unresolved; call resolve(n_features) first.")
        return (sum(self.feature_groups),) if self.scope == "joint" else self.feature_groups

    def counts(self, fraction: float) -> tuple[int, ...]:
        """The features kept per budget group at ``fraction``: ``max(1, round(f·p_g))``."""

        value = _validate_fraction(fraction, where="SparsityDesign.counts")
        return tuple(keep_count(value, size) for size in self.budget_groups)

    def estimator(self, n_components: int, fraction: float, *, scale: bool = True) -> SparsePLSDA:
        """An unfitted :class:`SparsePLSDA` at ``n_components`` and ``fraction``."""

        return SparsePLSDA(
            n_components, keep=self.counts(fraction), groups=self.budget_groups, scale=scale, max_iter=1000
        )

    def to_dict(self) -> dict[str, Any]:
        """The setting as recorded: scope and ladder (fullest first)."""

        return {"scope": self.scope, "keep_fractions": list(self.keep_fractions)}


def _resolve_budget(keep: Any, groups: Any, n_features: int) -> tuple[tuple[int, ...], tuple[int, ...]]:
    """Validate ``keep`` against ``groups`` and return ``(group sizes, kept counts)``."""

    sizes = (int(n_features),) if groups is None else _validate_groups(groups, where="SparsePLSDA")
    if sum(sizes) != n_features:
        raise ValueError(f"SparsePLSDA: groups {list(sizes)} cover {sum(sizes)} features, but X has {n_features}.")
    if keep is None:
        return sizes, sizes
    if isinstance(keep, Integral) and not isinstance(keep, bool):
        counts: tuple[Any, ...] = (keep,)
    else:
        try:
            counts = tuple(keep)
        except TypeError:
            raise ValueError(f"SparsePLSDA: keep {keep!r} must be a count or one count per group.") from None
    if len(counts) != len(sizes):
        raise ValueError(
            f"SparsePLSDA: keep {list(counts)} gives {len(counts)} count(s) for {len(sizes)} feature group(s)."
        )
    for count, size in zip(counts, sizes):
        if isinstance(count, bool) or not isinstance(count, Integral) or not 1 <= int(count) <= size:
            raise ValueError(f"SparsePLSDA: keep count {count!r} must be an integer in [1, {size}] for its group.")
    return sizes, tuple(int(count) for count in counts)


def soft_threshold(weights: np.ndarray, groups: Sequence[int], counts: Sequence[int]) -> np.ndarray:
    """Keep the ``counts[g]`` largest absolute weights of every group, shrunk by the largest discarded one.

    Every other weight becomes zero. Ties are broken by column order (the
    earlier column is kept). A group whose count equals its size is unchanged.
    Operates in place and returns ``weights``.
    """

    start = 0
    for size, count in zip(groups, counts):
        stop = start + size
        if count < size:
            segment = weights[start:stop]
            magnitude = np.abs(segment)
            order = np.argsort(-magnitude, kind="stable")
            threshold = magnitude[order[count]]
            kept = np.zeros(size, dtype=bool)
            kept[order[:count]] = True
            segment[:] = np.where(kept, np.sign(segment) * (magnitude - threshold), 0.0)
        start = stop
    return weights


def _sparse_power_method(
    X: np.ndarray,
    y: np.ndarray,
    groups: Sequence[int],
    counts: Sequence[int],
    *,
    max_iter: int,
    tol: float,
) -> tuple[np.ndarray, np.ndarray, int]:
    """scikit-learn's NIPALS power method (mode A, unnormalized Y weights) with soft-thresholded X weights."""

    eps = np.finfo(X.dtype).eps
    try:
        y_score = next(col for col in y.T if np.any(np.abs(col) > eps))
    except StopIteration as e:
        raise StopIteration("y residual is constant") from e

    x_weights_old: Any = 100  # init to big value for first convergence check
    for i in range(max_iter):
        x_weights = np.dot(X.T, y_score) / np.dot(y_score, y_score)
        soft_threshold(x_weights, groups, counts)
        x_weights /= np.sqrt(np.dot(x_weights, x_weights)) + eps
        x_score = np.dot(X, x_weights)
        y_weights = np.dot(y.T, x_score) / np.dot(x_score.T, x_score)
        y_score = np.dot(y, y_weights) / (np.dot(y_weights, y_weights) + eps)
        x_weights_diff = x_weights - x_weights_old
        if np.dot(x_weights_diff, x_weights_diff) < tol or y.shape[1] == 1:
            break
        x_weights_old = x_weights

    n_iter = i + 1
    if n_iter == max_iter:
        warnings.warn("Maximum number of iterations reached", ConvergenceWarning)
    return x_weights, y_weights, n_iter


class SparsePLSDA(PLSRegression):
    """PLS regression (PLS-DA on one-hot labels) with a feature budget per component.

    Per component the power iteration's X weight vector is soft-thresholded in
    every group: the ``keep`` largest absolute entries are kept and shrunk
    toward zero by the largest discarded absolute entry, every other entry is
    zero, and the vector is then normalized. Deflation is regression mode, so
    the fitted attributes and methods are those of ``PLSRegression`` with the
    same meanings (``x_rotations_ = W·(PᵀW)⁺``, ``transform``,
    ``inverse_transform``, ``predict``). When every count equals its group's
    size the fit is ``PLSRegression.fit``.

    Parameters
    ----------
    n_components : int, default=2
        Number of components.
    keep : int | sequence of int | None, default=None
        Features each component's X weight vector may load on: one count
        (one group) or one count per group of ``groups``, each in
        ``[1, group size]``. ``None`` keeps every feature, i.e. ordinary PLS.
    groups : sequence of int | None, default=None
        Sizes of consecutive column groups (blocks) the budget is counted in.
        ``None`` is one group of all features.
    scale : bool, default=True
        As :class:`~sklearn.cross_decomposition.PLSRegression`.
    max_iter : int, default=1000
        As ``PLSRegression`` (MOTCO's setting rather than sklearn's 500).
    tol : float, default=1e-06
        As ``PLSRegression``: the iteration stops when the squared change of
        the X weight vector is below ``tol``.
    copy : bool, default=True
        As ``PLSRegression``.

    Attributes
    ----------
    keep_ : ndarray of shape (n_components, n_groups)
        The resolved count per component and group.
    support_ : ndarray of shape (n_features, n_components)
        Boolean; the features each component's weight vector loads on
        (its non-zero entries).
    """

    _parameter_constraints: dict = {
        **PLSRegression._parameter_constraints,
        "keep": "no_validation",
        "groups": "no_validation",
    }

    def __init__(
        self,
        n_components: int = 2,
        keep: int | Sequence[int] | None = None,
        *,
        groups: Sequence[int] | None = None,
        scale: bool = True,
        max_iter: int = 1000,
        tol: float = 1e-06,
        copy: bool = True,
    ) -> None:
        super().__init__(n_components=n_components, scale=scale, max_iter=max_iter, tol=tol, copy=copy)
        self.keep = keep
        self.groups = groups

    @_fit_context(prefer_skip_nested_validation=True)
    def fit(self, X: Any, y: Any) -> SparsePLSDA:
        """Fit the sparse model; the budget is validated before any component is fitted."""

        n_features = int(np.shape(X)[1]) if np.ndim(X) == 2 else -1
        if n_features < 1:
            raise ValueError("SparsePLSDA: X must be a 2-D matrix with at least one column.")
        sizes, counts = _resolve_budget(self.keep, self.groups, n_features)
        if counts == sizes:
            # The full budget: the threshold is zero, so this is ordinary PLS.
            super().fit(X, y)
        else:
            self._fit_sparse(X, y, sizes, counts)
        n_components = self.x_weights_.shape[1]
        self.keep_ = np.tile(np.asarray(counts, dtype=int), (n_components, 1))
        self.support_ = self.x_weights_ != 0.0
        return self

    def _fit_sparse(self, X: Any, y: Any, groups: tuple[int, ...], counts: tuple[int, ...]) -> None:
        """``sklearn.cross_decomposition._PLS.fit`` (regression mode, NIPALS) with the thresholded power method."""

        check_consistent_length(X, y)
        X = validate_data(self, X, dtype=np.float64, force_writeable=True, copy=self.copy, ensure_min_samples=2)
        y = check_array(y, input_name="y", dtype=np.float64, force_writeable=True, copy=self.copy, ensure_2d=False)
        if y.ndim == 1:
            self._predict_1d = True
            y = y.reshape(-1, 1)
        else:
            self._predict_1d = False

        n, p = X.shape
        q = y.shape[1]
        n_components = self.n_components
        rank_upper_bound = min(n, p)
        if n_components > rank_upper_bound:
            raise ValueError(
                f"`n_components` upper bound is {rank_upper_bound}. "
                f"Got {n_components} instead. Reduce `n_components`."
            )
        self._norm_y_weights = False
        Xk, yk, self._x_mean, self._y_mean, self._x_std, self._y_std = _center_scale_xy(X, y, self.scale)

        self.x_weights_ = np.zeros((p, n_components))
        self.y_weights_ = np.zeros((q, n_components))
        self._x_scores = np.zeros((n, n_components))
        self._y_scores = np.zeros((n, n_components))
        self.x_loadings_ = np.zeros((p, n_components))
        self.y_loadings_ = np.zeros((q, n_components))
        self.n_iter_ = []

        y_eps = np.finfo(yk.dtype).eps
        for k in range(n_components):
            yk_mask = np.all(np.abs(yk) < 10 * y_eps, axis=0)
            yk[:, yk_mask] = 0.0
            try:
                x_weights, y_weights, n_iter_ = _sparse_power_method(
                    Xk, yk, groups, counts, max_iter=self.max_iter, tol=self.tol
                )
            except StopIteration as e:
                if str(e) != "y residual is constant":
                    raise
                warnings.warn(f"y residual is constant at iteration {k}")
                break
            self.n_iter_.append(n_iter_)
            _svd_flip_1d(x_weights, y_weights)

            x_scores = np.dot(Xk, x_weights)
            y_ss = np.dot(y_weights, y_weights)
            y_scores = np.dot(yk, y_weights) / y_ss
            x_loadings = np.dot(x_scores, Xk) / np.dot(x_scores, x_scores)
            Xk -= np.outer(x_scores, x_loadings)
            y_loadings = np.dot(x_scores, yk) / np.dot(x_scores, x_scores)
            yk -= np.outer(x_scores, y_loadings)

            self.x_weights_[:, k] = x_weights
            self.y_weights_[:, k] = y_weights
            self._x_scores[:, k] = x_scores
            self._y_scores[:, k] = y_scores
            self.x_loadings_[:, k] = x_loadings
            self.y_loadings_[:, k] = y_loadings

        self.x_rotations_ = np.dot(
            self.x_weights_, pinv(np.dot(self.x_loadings_.T, self.x_weights_), check_finite=False)
        )
        self.y_rotations_ = np.dot(
            self.y_weights_, pinv(np.dot(self.y_loadings_.T, self.y_weights_), check_finite=False)
        )
        self.coef_ = np.dot(self.x_rotations_, self.y_loadings_.T)
        self.coef_ = (self.coef_ * self._y_std).T / self._x_std
        self.intercept_ = self._y_mean
        self._n_features_out = self.x_rotations_.shape[1]
        self.x_scores_ = self._x_scores
        self.y_scores_ = self._y_scores

    def support_sizes(self, groups: Sequence[int] | None = None) -> np.ndarray:
        """Non-zero weights per component (rows) and group (columns); ``groups`` default to the fit's."""

        requested = groups if groups is not None else self.groups
        sizes = (self.support_.shape[0],) if requested is None else tuple(int(size) for size in requested)
        if sum(sizes) != self.support_.shape[0]:
            raise ValueError(f"groups {list(sizes)} do not cover the model's {self.support_.shape[0]} features.")
        bounds = np.cumsum((0,) + tuple(sizes))
        return np.array(
            [[int(self.support_[bounds[g] : bounds[g + 1], k].sum()) for g in range(len(sizes))]
             for k in range(self.support_.shape[1])],
            dtype=int,
        )
