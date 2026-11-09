"""Group × stage interaction tests: split into a span and its complement, and feature by feature.

:func:`interaction_split_test` measures the interaction of a multivariate
outcome ``Y`` — the variation the full model (group × stage) explains beyond
the nested reduced model (group + stage) — in three parts:

- ``total``: every coordinate of ``Y``;
- ``span``: the projection ``Y·Q·Qᵀ`` onto the columns of an orthonormal basis ``Q``;
- ``complement``: the projection ``Y·(I − Q·Qᵀ)`` onto its orthogonal complement.

Per part it reports the interaction sum of squares ``‖(H_F − H_R)·Y_part‖²``,
the full model's residual sum of squares ``‖(I − H_F)·Y_part‖²`` and the trace
pseudo-F ``(SS_int / df_int) / (RSS / df_res)``. Orthogonal projectors split
both sums exactly, so ``span`` + ``complement`` = ``total``.

With permutations, the reduced model's residuals are permuted by
Freedman–Lane with :func:`motco.stats.permutation.RRPP`'s serial index sequence
at the same seed, optionally rescaled by one of its residual corrections. The
pseudo-F is the test statistic: its numerator and denominator come from the
same permuted residuals, so a residual scale bias cancels, whereas an energy
summed over many coordinates turns it into a liberal null. The
interaction-SS p-value is returned as a diagnostic.

The permutation path never forms a permuted ``n × p`` outcome. A permuted
outcome's interaction and residual sums depend only on its permuted reduced
residuals, so each part's ``n × n`` residual Gram matrix is formed once and
every draw costs ``O(n²)`` per part, whatever the number of features.

:func:`interaction_feature_scan` tests the same interaction one outcome column
at a time. Each feature gets the one-column F ``(SS_int / df_int) / (RSS /
df_res)``, and the features are combined by their maximum: the global p-value
is the permutation p-value of the largest F, and each feature gets a
Westfall–Young step-down max-T adjusted p-value, which controls the
familywise error rate over the scanned features. The draws, models and
residual corrections are the split test's. A draw costs one ``rank(F) × n``
by ``n × p`` product, and no permuted ``n × p`` outcome is formed.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any, Optional, Union

import numpy as np
import pandas as pd

from motco.stats.permutation import residual_correction_scale

#: The parts the test reports, in order.
INTERACTION_SPLIT_PARTS: tuple[str, ...] = ("total", "span", "complement")

# Largest admissible deviation of the basis's Gram matrix from the identity.
_ORTHONORMAL_TOLERANCE = 1e-8

# Relative residual above which a reduced-model column is outside the full model's column space.
_NESTED_TOLERANCE = 1e-8

#: Adjusted p-value at or below which the feature scan calls a feature a hit.
INTERACTION_SCAN_ALPHA = 0.05

# Full-model RSS, relative to the centred sum of squares, at or below which a feature is degenerate.
_DEGENERATE_TOLERANCE = 1e-12


@dataclass(frozen=True)
class InteractionSplitPart:
    """Interaction and residual energy of one part.

    ``p_value_f`` and ``p_value_ss`` are ``(1 + #{null ≥ observed}) / (1 +
    permutations)`` for the pseudo-F and the interaction sum of squares, and
    ``None`` without permutations (or for an undefined pseudo-F). ``null_f``
    and ``null_ss`` hold the per-draw null values when they were requested.
    """

    dimension: int
    ss_interaction: float
    rss: float
    df_interaction: int
    df_residual: int
    f: float
    p_value_f: float | None = None
    p_value_ss: float | None = None
    null_f: np.ndarray | None = field(default=None, repr=False, compare=False)
    null_ss: np.ndarray | None = field(default=None, repr=False, compare=False)

    def to_dict(self) -> dict[str, Any]:
        """Observed quantities and p-values; null distributions are never included."""

        payload: dict[str, Any] = {
            "dimension": self.dimension,
            "ss_interaction": self.ss_interaction,
            "rss": self.rss,
            "df_interaction": self.df_interaction,
            "df_residual": self.df_residual,
            "f": self.f,
        }
        if self.p_value_f is not None or self.p_value_ss is not None:
            payload["p_value_f"] = self.p_value_f
            payload["p_value_ss"] = self.p_value_ss
        return payload


@dataclass(frozen=True)
class InteractionSplitResult:
    """The split interaction test: one :class:`InteractionSplitPart` per part.

    ``span_share`` is the observed ``span`` interaction SS over the ``total``
    one, ``None`` when the total interaction SS is zero. ``permutations``,
    ``seed`` and ``residual_correction`` are the settings the null was drawn
    with.
    """

    parts: dict[str, InteractionSplitPart]
    span_share: float | None
    permutations: int
    seed: int | None
    residual_correction: str | None

    @property
    def total(self) -> InteractionSplitPart:
        return self.parts["total"]

    @property
    def span(self) -> InteractionSplitPart:
        return self.parts["span"]

    @property
    def complement(self) -> InteractionSplitPart:
        return self.parts["complement"]

    def to_dict(self) -> dict[str, Any]:
        return {
            "parts": {name: self.parts[name].to_dict() for name in INTERACTION_SPLIT_PARTS},
            "span_share": self.span_share,
            "permutations": self.permutations,
            "seed": self.seed,
            "residual_correction": self.residual_correction,
        }


@dataclass(frozen=True, eq=False)
class InteractionScanResult:
    """The per-feature interaction scan: one F per outcome column and their maximum.

    ``ss_interaction``, ``rss`` and ``f`` hold one value per outcome column;
    ``f`` is ``NaN`` where ``degenerate`` is true (zero full-model residual
    variance, a constant column included). ``max_f`` and ``max_column`` are the
    largest F over the non-degenerate columns and its column (ties go to the
    lower column index).

    With permutations, ``p_value_global`` is ``(1 + #{draws whose largest null F
    ≥ max_f}) / (1 + permutations)``, ``p_value_adjusted`` holds each column's
    Westfall–Young step-down max-T adjusted p-value (``NaN`` where degenerate),
    and ``null_max_f`` / ``null_max_column`` hold each draw's largest null F and
    the column attaining it. ``null_f`` (draws × columns, ``NaN`` where
    degenerate) is kept only when requested. Without permutations all of these
    are ``None``.
    """

    ss_interaction: np.ndarray
    rss: np.ndarray
    f: np.ndarray
    degenerate: np.ndarray
    df_interaction: int
    df_residual: int
    max_f: float
    max_column: int
    permutations: int
    seed: int | None
    residual_correction: str | None
    p_value_global: float | None = None
    p_value_adjusted: np.ndarray | None = None
    null_max_f: np.ndarray | None = field(default=None, repr=False)
    null_max_column: np.ndarray | None = field(default=None, repr=False)
    null_f: np.ndarray | None = field(default=None, repr=False)

    @property
    def n_features(self) -> int:
        return int(self.f.size)

    @property
    def n_degenerate(self) -> int:
        return int(np.count_nonzero(self.degenerate))

    @property
    def degenerate_columns(self) -> np.ndarray:
        return np.flatnonzero(self.degenerate)

    def hits(self, alpha: float = INTERACTION_SCAN_ALPHA) -> np.ndarray:
        """Columns with adjusted p-value ≤ ``alpha``, by F descending (ties by column index).

        Raises
        ------
        ValueError
            When the scan ran without permutations.
        """

        if self.p_value_adjusted is None:
            raise ValueError("The scan ran without permutations; it has no adjusted p-values and no hits.")
        columns = np.flatnonzero(~self.degenerate & (self.p_value_adjusted <= alpha))
        return columns[np.lexsort((columns, -self.f[columns]))]

    def to_dict(self) -> dict[str, Any]:
        """Counts, df, the maximum and, with permutations, the global p-value and the hits by column index."""

        payload: dict[str, Any] = {
            "n_features": self.n_features,
            "n_degenerate": self.n_degenerate,
            "degenerate_columns": [int(column) for column in self.degenerate_columns],
            "df_interaction": self.df_interaction,
            "df_residual": self.df_residual,
            "max_f": self.max_f,
            "max_column": self.max_column,
            "permutations": self.permutations,
            "seed": self.seed,
            "residual_correction": self.residual_correction,
        }
        if self.p_value_adjusted is not None:
            assert self.p_value_global is not None
            payload["p_value_global"] = self.p_value_global
            payload["hit_alpha"] = INTERACTION_SCAN_ALPHA
            adjusted = self.p_value_adjusted
            payload["hits"] = [
                {"column": int(column), "f": float(self.f[column]), "p_value_adjusted": float(adjusted[column])}
                for column in self.hits()
            ]
        return payload


def interaction_split_test(
    Y: Union[pd.DataFrame, np.ndarray],
    model_full: Union[pd.DataFrame, np.ndarray],
    model_reduced: Union[pd.DataFrame, np.ndarray],
    basis: np.ndarray,
    *,
    permutations: int = 0,
    seed: Optional[int] = None,
    residual_correction: Optional[str] = None,
    return_null: bool = False,
) -> InteractionSplitResult:
    """Test the group × stage interaction inside a span and in its orthogonal complement.

    Parameters
    ----------
    Y: Union[pd.DataFrame, np.ndarray]
        Outcome matrix, ``n × p``, in the coordinates the basis is orthonormal in.
    model_full: Union[pd.DataFrame, np.ndarray]
        Model matrix of the full model (group × stage, plus any covariates), including intercept.
    model_reduced: Union[pd.DataFrame, np.ndarray]
        Model matrix of the reduced model (group + stage, plus the same covariates), nested in the full one.
    basis: np.ndarray
        ``p × r`` matrix with orthonormal columns (tolerance 1e-8), ``1 ≤ r < p``.
    permutations: int
        Number of Freedman–Lane draws; 0 reports the observed quantities only.
    seed: Optional[int]
        Seed of ``np.random.default_rng``; each draw is one ``rng.permutation(n)``,
        RRPP's serial sequence, so the test pairs with a trajectory analysis at the same seed.
    residual_correction: Optional[str]
        ``"df"``, ``"leverage"`` or ``None``, applied to the reduced residuals before any draw
        (see :func:`motco.stats.permutation.residual_correction_scale`). The observed quantities
        do not depend on it; under ``"df"`` the pseudo-F null is unchanged and the SS null scales by κ².
    return_null: bool
        Keep the per-draw null pseudo-F and SS of every part.

    Returns
    -------
    InteractionSplitResult

    Raises
    ------
    ValueError
        Before any draw, for inputs whose rows or columns disagree, non-finite values, a basis
        whose columns are not orthonormal, an empty basis or complement, models that are not
        nested, a non-positive df_int or df_res, a negative permutation count, or an unknown or
        undefined residual correction.
    """

    outcome = np.asarray(Y, dtype=float)
    full = np.asarray(model_full, dtype=float)
    reduced = np.asarray(model_reduced, dtype=float)
    Q = np.asarray(basis, dtype=float)
    _validate_shapes(outcome, full, reduced, Q)
    _check_permutations(permutations)
    n, p = outcome.shape
    r = Q.shape[1]
    spaces = _model_spaces(full, reduced)
    df_interaction = spaces.df_interaction
    df_residual = spaces.df_residual
    U_reduced = spaces.U_reduced
    scale = None if residual_correction is None else residual_correction_scale(reduced, residual_correction)

    residuals = outcome - U_reduced @ (U_reduced.T @ outcome)
    in_span = residuals @ Q
    # Formed directly: G_total − G_span would cancel badly when the span holds most of the energy.
    out_of_span = residuals - in_span @ Q.T
    grams = {
        "total": residuals @ residuals.T,
        "span": in_span @ in_span.T,
        "complement": out_of_span @ out_of_span.T,
    }
    dimensions = {"total": p, "span": r, "complement": p - r}
    weights = spaces.weights

    observed: dict[str, tuple[float, float]] = {
        name: _sums(gram, weights, df_interaction, float(np.trace(gram))) for name, gram in grams.items()
    }
    null_ss: dict[str, np.ndarray] = {}
    null_rss: dict[str, np.ndarray] = {}
    if permutations > 0:
        scaled = grams if scale is None else {name: gram * np.outer(scale, scale) for name, gram in grams.items()}
        traces = {name: float(np.trace(gram)) for name, gram in scaled.items()}
        for name in INTERACTION_SPLIT_PARTS:
            null_ss[name] = np.empty(permutations)
            null_rss[name] = np.empty(permutations)
        rng = np.random.default_rng(seed)
        for draw in range(permutations):
            idx = rng.permutation(n)
            # Row i of the permuted residuals is row idx[i], and Σ_ij w_i·G[idx_i, idx_j]·w_j = w'ᵀ·G·w'
            # with w' = w[argsort(idx)], so the permuted Gram matrix is never formed.
            permuted_weights = weights[np.argsort(idx)]
            for name, gram in scaled.items():
                ss, rss = _sums(gram, permuted_weights, df_interaction, traces[name])
                null_ss[name][draw] = ss
                null_rss[name][draw] = rss

    parts: dict[str, InteractionSplitPart] = {}
    for name in INTERACTION_SPLIT_PARTS:
        ss, rss = observed[name]
        f_value = _pseudo_f(np.array([ss]), np.array([rss]), df_interaction, df_residual)[0]
        p_value_f: float | None = None
        p_value_ss: float | None = None
        part_null_f: np.ndarray | None = None
        if permutations > 0:
            part_null_f = _pseudo_f(null_ss[name], null_rss[name], df_interaction, df_residual)
            p_value_ss = _p_value(null_ss[name], ss)
            p_value_f = _p_value(part_null_f, f_value) if np.isfinite(f_value) else None
        parts[name] = InteractionSplitPart(
            dimension=dimensions[name],
            ss_interaction=ss,
            rss=rss,
            df_interaction=df_interaction,
            df_residual=df_residual,
            f=float(f_value),
            p_value_f=p_value_f,
            p_value_ss=p_value_ss,
            null_f=part_null_f if return_null else None,
            null_ss=null_ss[name] if (return_null and permutations > 0) else None,
        )
    total_ss = parts["total"].ss_interaction
    span_share = None if total_ss <= 0.0 else float(np.clip(parts["span"].ss_interaction / total_ss, 0.0, 1.0))
    return InteractionSplitResult(
        parts=parts,
        span_share=span_share,
        permutations=int(permutations),
        seed=seed,
        residual_correction=residual_correction,
    )


def interaction_feature_scan(
    Y: Union[pd.DataFrame, np.ndarray],
    model_full: Union[pd.DataFrame, np.ndarray],
    model_reduced: Union[pd.DataFrame, np.ndarray],
    *,
    permutations: int = 0,
    seed: Optional[int] = None,
    residual_correction: Optional[str] = None,
    return_null: bool = False,
) -> InteractionScanResult:
    """Test the group × stage interaction feature by feature, combined by the maximum F.

    For each outcome column y, the interaction sum of squares is
    ``‖(H_F − H_R)·y‖²``, the residual sum of squares ``‖(I − H_F)·y‖²`` and
    ``F = (SS_int / df_int) / (RSS / df_res)``, with ``df_int = rank(F) −
    rank(R)`` and ``df_res = n − rank(F)``. Each F ignores its column's location
    and scale, and the per-feature sums add to
    :func:`interaction_split_test`'s ``total`` sums.

    A column whose RSS is at most 1e-12 of its centred sum of squares, or that
    is constant, is degenerate: it gets no F and no p-value and never enters a
    maximum.

    With permutations, the reduced residuals are permuted by Freedman–Lane with
    RRPP's serial index sequence (one ``rng.permutation(n)`` per draw from
    ``np.random.default_rng(seed)``), and every feature's F is evaluated at
    every draw from the same permuted residuals. The global p-value is the
    permutation p-value of the largest F. Each feature's adjusted p-value is
    Westfall–Young step-down max-T: with features ordered by observed F
    (descending, ties by column index), rank j's raw value is ``(1 + #{draws
    whose largest null F over ranks ≥ j is ≥ the j-th observed F}) / (1 + B)``,
    and the adjusted p-value is the running maximum down the order. The top
    feature's adjusted p-value is the global p-value.

    Parameters
    ----------
    Y: Union[pd.DataFrame, np.ndarray]
        Outcome matrix, ``n × p``.
    model_full: Union[pd.DataFrame, np.ndarray]
        Model matrix of the full model (group × stage, plus any covariates), including intercept.
    model_reduced: Union[pd.DataFrame, np.ndarray]
        Model matrix of the reduced model (group + stage, plus the same covariates), nested in the full one.
    permutations: int
        Number of Freedman–Lane draws; 0 reports the observed F only.
    seed: Optional[int]
        Seed of ``np.random.default_rng``, so the scan pairs with a trajectory analysis at the same seed.
    residual_correction: Optional[str]
        ``"df"``, ``"leverage"`` or ``None``, applied to the reduced residuals before any draw
        (see :func:`motco.stats.permutation.residual_correction_scale`). The observed F do not depend on
        it, and a constant multiplier (``"df"``, or ``"leverage"`` at equal leverages) cancels in every F,
        so it leaves every null F and p-value unchanged.
    return_null: bool
        Keep the draws × columns null F.

    Returns
    -------
    InteractionScanResult

    Raises
    ------
    ValueError
        Before any draw, for inputs whose rows disagree, non-finite values, models that are not nested,
        a non-positive df_int or df_res, a negative permutation count, an unknown or undefined residual
        correction, or an outcome whose every column is degenerate.
    """

    outcome = np.asarray(Y, dtype=float)
    full = np.asarray(model_full, dtype=float)
    reduced = np.asarray(model_reduced, dtype=float)
    named = (("Y", outcome), ("model_full", full), ("model_reduced", reduced))
    _check_two_dimensional(named)
    _check_rows(outcome, full, reduced)
    _check_finite(named)
    _check_permutations(permutations)
    n, p = outcome.shape
    spaces = _model_spaces(full, reduced)
    df_interaction = spaces.df_interaction
    df_residual = spaces.df_residual
    scale = None if residual_correction is None else residual_correction_scale(reduced, residual_correction)

    residuals = outcome - spaces.U_reduced @ (spaces.U_reduced.T @ outcome)
    ss, rss = _feature_sums(residuals, spaces.weights, df_interaction, np.sum(residuals**2, axis=0))
    centred_ss = np.sum((outcome - outcome.mean(axis=0)) ** 2, axis=0)
    constant = np.all(outcome == outcome[:1], axis=0)
    degenerate = constant | (rss <= _DEGENERATE_TOLERANCE * centred_ss)
    valid = np.flatnonzero(~degenerate)
    if valid.size == 0:
        raise ValueError(
            f"No feature has a defined F: all {p} outcome column(s) are degenerate (zero full-model residual "
            "variance)."
        )
    f_values = np.full(p, np.nan)
    f_values[valid] = _pseudo_f(ss[valid], rss[valid], df_interaction, df_residual)
    # Step-down order: F descending, ties by column index.
    order = valid[np.lexsort((valid, -f_values[valid]))]
    observed_sorted = f_values[order]
    max_column = int(order[0])
    max_f = float(observed_sorted[0])
    if permutations == 0:
        return InteractionScanResult(
            ss_interaction=ss,
            rss=rss,
            f=f_values,
            degenerate=degenerate,
            df_interaction=df_interaction,
            df_residual=df_residual,
            max_f=max_f,
            max_column=max_column,
            permutations=int(permutations),
            seed=seed,
            residual_correction=residual_correction,
        )

    permuted = residuals[:, order]
    if scale is not None and not np.all(scale == scale[0]):
        # A constant multiplier cancels in every F, so only a varying one is applied.
        permuted = permuted * scale[:, None]
    norms = np.sum(permuted**2, axis=0)
    exceed = np.zeros(order.size, dtype=np.int64)
    null_max_f = np.empty(permutations)
    null_max_column = np.empty(permutations, dtype=np.int64)
    null_f = np.full((permutations, p), np.nan) if return_null else None
    rng = np.random.default_rng(seed)
    for draw in range(permutations):
        idx = rng.permutation(n)
        # Row i of the permuted residuals is row idx[i], so Wᵀ·E[π] = W[argsort(idx)]ᵀ·E.
        draw_ss, draw_rss = _feature_sums(permuted, spaces.weights[np.argsort(idx)], df_interaction, norms)
        draw_f = _pseudo_f(draw_ss, draw_rss, df_interaction, df_residual)
        # 0 / 0 (a permuted column inside col(reduced)) carries no interaction.
        draw_f = np.where(np.isnan(draw_f), 0.0, draw_f)
        suffix_max = np.maximum.accumulate(draw_f[::-1])[::-1]
        exceed += suffix_max >= observed_sorted
        position = int(np.argmax(draw_f))
        null_max_f[draw] = draw_f[position]
        null_max_column[draw] = order[position]
        if null_f is not None:
            null_f[draw, order] = draw_f
    adjusted = np.full(p, np.nan)
    adjusted[order] = np.maximum.accumulate((1.0 + exceed) / (1.0 + permutations))
    return InteractionScanResult(
        ss_interaction=ss,
        rss=rss,
        f=f_values,
        degenerate=degenerate,
        df_interaction=df_interaction,
        df_residual=df_residual,
        max_f=max_f,
        max_column=max_column,
        p_value_global=_p_value(null_max_f, max_f),
        p_value_adjusted=adjusted,
        null_max_f=null_max_f,
        null_max_column=null_max_column,
        null_f=null_f,
        permutations=int(permutations),
        seed=seed,
        residual_correction=residual_correction,
    )


def _validate_shapes(outcome: np.ndarray, full: np.ndarray, reduced: np.ndarray, Q: np.ndarray) -> None:
    named = (("Y", outcome), ("model_full", full), ("model_reduced", reduced), ("basis", Q))
    _check_two_dimensional(named)
    _check_rows(outcome, full, reduced)
    p = outcome.shape[1]
    if Q.shape[0] != p:
        raise ValueError(f"Y has {p} columns but basis has {Q.shape[0]} rows — they must match.")
    _check_finite(named)
    r = Q.shape[1]
    if r == 0:
        raise ValueError("The basis has no columns; the span is empty.")
    if r >= p:
        raise ValueError(f"The complement is empty: the basis has {r} columns and spans all {p} coordinates of Y.")
    deviation = np.abs(Q.T @ Q - np.eye(r))
    largest = float(deviation.max())
    if largest > _ORTHONORMAL_TOLERANCE:
        raise ValueError(
            f"The basis columns are not orthonormal: the largest deviation of basisᵀ·basis from the identity is "
            f"{largest:.3g} (tolerance {_ORTHONORMAL_TOLERANCE:g})."
        )


def _check_two_dimensional(named: tuple[tuple[str, np.ndarray], ...]) -> None:
    for name, matrix in named:
        if matrix.ndim != 2:
            raise ValueError(f"{name} must be a 2-D matrix, got {matrix.ndim} dimension(s).")


def _check_rows(outcome: np.ndarray, full: np.ndarray, reduced: np.ndarray) -> None:
    n = outcome.shape[0]
    for name, matrix in (("model_full", full), ("model_reduced", reduced)):
        if matrix.shape[0] != n:
            raise ValueError(f"Y has {n} rows but {name} has {matrix.shape[0]} rows — number of rows must match.")


def _check_finite(named: tuple[tuple[str, np.ndarray], ...]) -> None:
    for name, matrix in named:
        if not np.all(np.isfinite(matrix)):
            raise ValueError(f"{name} contains NaN or Inf values.")


def _check_permutations(permutations: int) -> None:
    if permutations < 0:
        raise ValueError(f"permutations must be greater than or equal to 0, got {permutations}.")


@dataclass(frozen=True)
class _ModelSpaces:
    """Column spaces of a nested pair of models.

    ``weights`` is an orthonormal basis of col(full) with the interaction
    directions (col(full) ⊖ col(reduced)) first and the reduced model's after.
    """

    U_reduced: np.ndarray
    weights: np.ndarray
    df_interaction: int
    df_residual: int


def _model_spaces(full: np.ndarray, reduced: np.ndarray) -> _ModelSpaces:
    """Orthonormal bases and df of a full and a reduced model, which must be nested with positive df."""

    n = full.shape[0]
    U_full = _column_basis(full)
    U_reduced = _column_basis(reduced)
    rank_full = U_full.shape[1]
    rank_reduced = U_reduced.shape[1]
    nested_residual = reduced - U_full @ (U_full.T @ reduced)
    column_norms = np.linalg.norm(reduced, axis=0)
    relative = np.linalg.norm(nested_residual, axis=0) / np.where(column_norms > 0, column_norms, 1.0)
    if relative.size and float(relative.max()) > _NESTED_TOLERANCE:
        column = int(np.argmax(relative))
        raise ValueError(
            f"The models are not nested: column {column} of model_reduced lies outside the column space of "
            f"model_full (relative residual {float(relative[column]):.3g})."
        )
    df_interaction = rank_full - rank_reduced
    df_residual = n - rank_full
    if df_interaction <= 0:
        raise ValueError(
            f"df_int must be positive, but model_full has rank {rank_full} and model_reduced rank {rank_reduced}."
        )
    if df_residual <= 0:
        raise ValueError(f"df_res must be positive, but n = {n} and model_full has rank {rank_full}.")
    U_interaction = _interaction_basis(U_full, U_reduced, df_interaction)
    return _ModelSpaces(
        U_reduced=U_reduced,
        weights=np.hstack([U_interaction, U_reduced]),
        df_interaction=df_interaction,
        df_residual=df_residual,
    )


def _column_basis(matrix: np.ndarray) -> np.ndarray:
    """Orthonormal basis of the column space, at numpy's ``matrix_rank`` tolerance."""

    U, s, _ = np.linalg.svd(matrix, full_matrices=False)
    tolerance = s.max() * max(matrix.shape) * np.finfo(float).eps if s.size else 0.0
    return U[:, : int(np.sum(s > tolerance))]


def _interaction_basis(U_full: np.ndarray, U_reduced: np.ndarray, df_interaction: int) -> np.ndarray:
    """Orthonormal basis of col(full) ⊖ col(reduced), so that H_F − H_R = U·Uᵀ."""

    beyond = U_full - U_reduced @ (U_reduced.T @ U_full)
    U, _, _ = np.linalg.svd(beyond, full_matrices=False)
    return U[:, :df_interaction]


def _sums(gram: np.ndarray, weights: np.ndarray, df_interaction: int, trace: float) -> tuple[float, float]:
    """Interaction SS ⟨H_int, G⟩ and residual SS tr(G) − ⟨H_F, G⟩ from the weight columns' quadratic forms."""

    quadratic = np.sum((gram @ weights) * weights, axis=0)
    ss = max(float(np.sum(quadratic[:df_interaction])), 0.0)
    rss = max(trace - float(np.sum(quadratic)), 0.0)
    return ss, rss


def _feature_sums(
    residuals: np.ndarray, weights: np.ndarray, df_interaction: int, norms: np.ndarray
) -> tuple[np.ndarray, np.ndarray]:
    """Per-column interaction SS ‖U_intᵀ·e‖² and residual SS ‖e‖² − ‖Wᵀ·e‖² (clipped at 0)."""

    projected = weights.T @ residuals
    squared = projected**2
    ss = np.sum(squared[:df_interaction], axis=0)
    rss = np.maximum(norms - np.sum(squared, axis=0), 0.0)
    return ss, rss


def _pseudo_f(ss: np.ndarray, rss: np.ndarray, df_interaction: int, df_residual: int) -> np.ndarray:
    with np.errstate(divide="ignore", invalid="ignore"):
        return (ss / df_interaction) / (rss / df_residual)


def _p_value(null: np.ndarray, observed: float) -> float:
    return float((1.0 + np.sum(null >= observed)) / (1.0 + null.size))


__all__ = [
    "INTERACTION_SCAN_ALPHA",
    "INTERACTION_SPLIT_PARTS",
    "InteractionScanResult",
    "InteractionSplitPart",
    "InteractionSplitResult",
    "interaction_feature_scan",
    "interaction_split_test",
]
