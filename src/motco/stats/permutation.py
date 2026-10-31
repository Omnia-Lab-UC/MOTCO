from __future__ import annotations

import logging
import multiprocessing
import os
from typing import Optional, Union

import numpy as np
import pandas as pd
from tqdm import tqdm

from motco.stats.trajectory import (
    estimate_betas,
    estimate_difference,
    pooled_relative_eigengap,
)

logger = logging.getLogger(__name__)

#: Residual corrections :func:`RRPP` accepts; ``None`` is the uncorrected Freedman–Lane null.
RESIDUAL_CORRECTIONS: tuple[str, ...] = ("df", "leverage")

# Smallest admissible ``1 − h_ii`` under the leverage correction.
_LEVERAGE_TOLERANCE = 1e-10


def residual_correction_scale(
    model_reduced: Union[pd.DataFrame, np.ndarray],
    correction: str,
) -> np.ndarray:
    """Per-sample multipliers a residual correction applies to the reduced model's residuals.

    The reduced model's mean squared residual is RSS/n, not the unbiased
    RSS/(n − r), so permuting the raw residuals draws a null that is too
    narrow for scale-dependent statistics. The corrections restore that scale:

    - ``"df"``: every residual is multiplied by κ = √(n / (n − r));
    - ``"leverage"``: residual i is divided by √(1 − h_ii), with h the reduced
      model's hat matrix.

    r is the numerical rank of ``model_reduced`` (numpy's ``matrix_rank``
    tolerance on its thin SVD), so a rank-deficient model uses its rank, not
    its column count. A balanced one-factor model gives every sample the same
    leverage r/n, where both corrections coincide.

    Parameters
    ----------
    model_reduced: Union[pd.DataFrame, np.ndarray]
        Model matrix of the reduced model, including intercept.
    correction: str
        One of :data:`RESIDUAL_CORRECTIONS`.

    Returns
    -------
    np.ndarray
        One multiplier per row of ``model_reduced``.

    Raises
    ------
    ValueError
        For an unknown correction, for ``"df"`` when n − r ≤ 0, and for
        ``"leverage"`` when a row has leverage 1 to numerical tolerance.
    """
    if correction not in RESIDUAL_CORRECTIONS:
        raise ValueError(
            f"Unknown residual_correction {correction!r}; allowed: {list(RESIDUAL_CORRECTIONS)} or None."
        )
    X = np.asarray(model_reduced, dtype=float)
    n, p = X.shape
    U, s, _ = np.linalg.svd(X, full_matrices=False)
    tolerance = s.max() * max(n, p) * np.finfo(float).eps if s.size else 0.0
    rank = int(np.sum(s > tolerance))
    if correction == "df":
        if n - rank <= 0:
            raise ValueError(
                f"The df residual correction needs n - rank > 0, but the reduced model has n = {n} rows "
                f"and rank {rank}."
            )
        return np.full(n, np.sqrt(n / (n - rank)))
    complement = 1.0 - np.sum(U[:, :rank] ** 2, axis=1)
    saturated = np.flatnonzero(complement <= _LEVERAGE_TOLERANCE)
    if saturated.size:
        raise ValueError(
            f"The leverage residual correction is undefined: row {int(saturated[0])} of the reduced model has "
            f"leverage 1 (1 - h_ii = {complement[saturated[0]]:.3g})."
        )
    return 1.0 / np.sqrt(complement)


class _RRPPWorker:
    """Picklable callable that runs a chunk of RRPP permutations."""

    def __init__(
        self,
        y_hat: np.ndarray,
        y_res: np.ndarray,
        model_full: Union[pd.DataFrame, np.ndarray],
        ls_means: Union[pd.DataFrame, np.ndarray],
        contrast: list[list[int]],
        return_eigengaps: bool = False,
        orientation_estimator: str = "pc1",
    ) -> None:
        self.y_hat = y_hat
        self.y_res = y_res
        self.model_full = model_full
        self.ls_means = ls_means
        self.contrast = contrast
        self.return_eigengaps = return_eigengaps
        self.orientation_estimator = orientation_estimator

    def __call__(self, n_iters: int, seed: int):
        rng = np.random.default_rng(seed)
        n = self.y_res.shape[0]
        out_d, out_a, out_s = [], [], []
        out_g: list[float | None] = []
        for _ in range(n_iters):
            # The permutation index is drawn identically whether or not the
            # spectrum is recorded: the spectrum path consumes no randomness.
            idx = rng.permutation(n)
            y_random = self.y_hat + self.y_res[idx, :]
            if self.return_eigengaps:
                d, a, s, spectra = estimate_difference(
                    y_random,
                    self.model_full,
                    self.ls_means,
                    self.contrast,
                    return_spectra=True,
                    orientation_estimator=self.orientation_estimator,
                )
                out_g.append(pooled_relative_eigengap(spectra))
            else:
                d, a, s = estimate_difference(
                    y_random,
                    self.model_full,
                    self.ls_means,
                    self.contrast,
                    orientation_estimator=self.orientation_estimator,
                )
            out_d.append(d)
            out_a.append(a)
            out_s.append(s)
        if self.return_eigengaps:
            return out_d, out_a, out_s, out_g
        return out_d, out_a, out_s


def RRPP(
    Y: Union[pd.DataFrame, np.ndarray],
    model_full: Union[pd.DataFrame, np.ndarray],
    model_reduced: Union[pd.DataFrame, np.ndarray],
    LS_means: Union[pd.DataFrame, np.ndarray],
    contrast: list[list[int]],
    permutations: int = 999,
    n_jobs: Optional[int] = None,
    progress: bool = True,
    seed: Optional[int] = None,
    return_eigengaps: bool = False,
    orientation_estimator: str = "pc1",
    residual_correction: Optional[str] = None,
) -> tuple:
    """
    Residual Randomization in a Permutation Procedure to evaluate
    linear models.

    Parameters
    ----------
    Y: Union[pd.DataFrame, np.ndarray]
        Outcome matrix.
    model_full: Union[pd.DataFrame, np.ndarray]
        Model matrix for full model, including intercept.
    model_reduced: Union[pd.DataFrame, np.ndarray]
        Model matrix for reduced model, including intercept.
    LS_means: Union[pd.DataFrame, np.ndarray]
        Least-squares means to estimate.
    contrast: list[list[int]]
        Indices indicating the groups to compare based on LS means.
        Each list must contain the cohorts that belong to the same group.
    permutations: int
        Number of permutations.
    n_jobs: Optional[int]
        If provided and > 1, run permutations in parallel using multiple
        worker processes. Use -1 to use all available CPUs. When None or 1,
        runs single-threaded (backward-compatible default).
    seed: Optional[int]
        Optional seed for reproducible residual randomization.
    return_eigengaps: bool
        When ``True``, additionally return the pooled stage-mean configuration
        relative eigengap of each permutation draw, so a replicate can locate its
        observed eigengap against its own permutation distribution. The spectrum
        consumes no randomness, so the permutation draws and the three returned
        null distributions are identical with the flag on or off. Off by default:
        the returned tuple is then exactly the historical three-element one.
    orientation_estimator: str
        Orientation estimator for every permuted ``angle`` (``"pc1"``, the
        default, or ``"net_displacement"``); see
        :func:`motco.stats.trajectory.estimate_difference`.
    residual_correction: Optional[str]
        Rescale the reduced model's residuals once before they are permuted:
        ``"df"`` multiplies each by κ = √(n / (n − r)) and ``"leverage"``
        divides each by √(1 − h_ii) (see :func:`residual_correction_scale`).
        The permutation indices, the reduced fitted values and the full-model
        estimation are unchanged, so only the null distributions move. ``None``
        (the default) permutes the raw residuals, exactly as before the
        setting existed. An unknown or undefined correction raises before any
        permutation is drawn.

    Returns
    -------
    dist_delta: list[float]
        Distribution of deltas.
    dist_angle: list[float]
        Distribution of angles.
    """
    # --- Input validation ---
    _Y = np.asarray(Y, dtype=float)
    _Xf = np.asarray(model_full, dtype=float)
    _Xr = np.asarray(model_reduced, dtype=float)
    _LS = np.asarray(LS_means, dtype=float)
    if _Y.shape[0] != _Xf.shape[0]:
        raise ValueError(
            f"Y has {_Y.shape[0]} rows but model_full has {_Xf.shape[0]} rows — "
            "number of rows must match."
        )
    if _Y.shape[0] != _Xr.shape[0]:
        raise ValueError(
            f"Y has {_Y.shape[0]} rows but model_reduced has {_Xr.shape[0]} rows — "
            "number of rows must match."
        )
    if _LS.shape[1] != _Xf.shape[1]:
        raise ValueError(
            f"LS_means has {_LS.shape[1]} columns but model_full has {_Xf.shape[1]} columns — "
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
    if not np.all(np.isfinite(_Xf)):
        raise ValueError("model_full contains NaN or Inf values.")
    if not np.all(np.isfinite(_Xr)):
        raise ValueError("model_reduced contains NaN or Inf values.")
    correction_scale = None if residual_correction is None else residual_correction_scale(_Xr, residual_correction)
    # --- End validation ---
    # Set Y to be pandas df
    Y = pd.DataFrame(Y)
    # Set-up permutation procedure
    betas_red = estimate_betas(model_reduced, Y)
    # Predicted values from reduced model
    y_hat = np.matmul(model_reduced, betas_red)
    y_hat.index = Y.index
    # Resdiuals of reduced mode (these are the permuted units)
    y_res = Y - y_hat
    # Prepare NumPy views to avoid pandas inside the permutation loop
    y_hat_np = np.asarray(y_hat)
    y_res_np = np.asarray(y_res)
    if correction_scale is not None:
        # Rescaled once, before either path: the draws below are untouched.
        y_res_np = y_res_np * correction_scale[:, None]
    # Prepare containers
    deltas: list[np.ndarray] = []
    angles: list[np.ndarray] = []
    shapes: list[np.ndarray] = []
    eigengaps: list[float | None] = []

    # Serial path
    if n_jobs in (None, 1):
        n = y_res_np.shape[0]
        rng = np.random.default_rng(seed)
        for _ in tqdm(range(permutations), desc="RRPP", unit="perm", disable=not progress):
            # Drawn before any spectrum work, and identically either way.
            idx = rng.permutation(n)
            y_random = y_hat_np + y_res_np[idx, :]
            if return_eigengaps:
                d, a, s, spectra = estimate_difference(
                    y_random,
                    model_full,
                    LS_means,
                    contrast,
                    return_spectra=True,
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

    # Parallel path
    n_workers = (os.cpu_count() or 1) if n_jobs == -1 else max(1, n_jobs or 1)
    n_workers = min(n_workers, max(1, permutations))
    logger.info("Running %d permutations across %d workers", permutations, n_workers)
    base = permutations // n_workers
    rem = permutations % n_workers
    counts = [base + (1 if i < rem else 0) for i in range(n_workers)]

    ss = np.random.SeedSequence(seed)
    seeds = [int(s.generate_state(1)[0]) for s in ss.spawn(n_workers)]

    worker = _RRPPWorker(
        y_hat_np,
        y_res_np,
        model_full,
        LS_means,
        contrast,
        return_eigengaps=return_eigengaps,
        orientation_estimator=orientation_estimator,
    )
    with multiprocessing.Pool(processes=n_workers) as pool:
        parts = pool.starmap(worker, zip(counts, seeds))

    for part in parts:
        deltas.extend(part[0])
        angles.extend(part[1])
        shapes.extend(part[2])
        if return_eigengaps:
            eigengaps.extend(part[3])

    if return_eigengaps:
        return deltas, angles, shapes, eigengaps
    return deltas, angles, shapes
