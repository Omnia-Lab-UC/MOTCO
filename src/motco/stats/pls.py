"""Partial Least Squares utilities (PLS-DA, VIP computation).

This module provides:
- `plsda_doubleCV`: double cross-validation over number of components
  for PLS-DA classification using AUROC as criterion, optionally residualizing
  the predictors on nuisance covariates inside every training set.
- `NuisanceDesign`, `fit_nuisance_coefficients`, `residualize`: the per-fold
  covariate residualizer used by that option.
- `calculate_vips`: compute Variable Importance in Projection (VIP) scores
  from a fitted `sklearn.cross_decomposition.PLSRegression` model.
"""

import logging
import multiprocessing
from dataclasses import dataclass
from functools import partial
from itertools import repeat
from typing import NotRequired, Sequence, TypedDict, Union, overload

import numpy as np
import pandas as pd
from sklearn.cross_decomposition import PLSRegression
from sklearn.metrics import roc_auc_score
from sklearn.model_selection import RepeatedStratifiedKFold, StratifiedKFold
from sklearn.preprocessing import OneHotEncoder
from tqdm import tqdm

logger = logging.getLogger(__name__)


class PLSDAResult(TypedDict):
    table: pd.DataFrame
    models: list[PLSRegression]
    nuisance_min_training_df: NotRequired[int]


@dataclass(frozen=True)
class NuisanceDesign:
    """Row-aligned nuisance design for covariate residualization.

    ``conditioning`` (``n × k``) holds the design the covariate coefficients
    are estimated conditional on (for example drop-first stage indicators),
    without an intercept; ``covariates`` (``n × q``) holds the covariate
    columns whose fitted part is removed. The residualization regresses each
    predictor on ``[1, conditioning, covariates]`` and subtracts only the
    covariate part, so signal collinear with the covariates but explained by
    the conditioning design is kept (Frisch–Waugh–Lovell). ``kind`` optionally
    names the conditioning design for callers that record it.
    """

    conditioning: np.ndarray
    covariates: np.ndarray
    conditioning_names: tuple[str, ...]
    covariate_names: tuple[str, ...]
    kind: str | None = None

    def __post_init__(self) -> None:
        conditioning = np.asarray(self.conditioning, dtype=float)
        covariates = np.asarray(self.covariates, dtype=float)
        if conditioning.ndim != 2 or covariates.ndim != 2:
            raise ValueError("NuisanceDesign conditioning and covariates must be 2-D matrices.")
        if conditioning.shape[0] != covariates.shape[0]:
            raise ValueError(
                f"NuisanceDesign conditioning has {conditioning.shape[0]} rows but covariates have "
                f"{covariates.shape[0]}."
            )
        if covariates.shape[1] == 0:
            raise ValueError("NuisanceDesign needs at least one covariate column.")
        if conditioning.shape[1] != len(self.conditioning_names) or covariates.shape[1] != len(
            self.covariate_names
        ):
            raise ValueError("NuisanceDesign column names do not match its matrices.")
        object.__setattr__(self, "conditioning", conditioning)
        object.__setattr__(self, "covariates", covariates)
        object.__setattr__(self, "conditioning_names", tuple(self.conditioning_names))
        object.__setattr__(self, "covariate_names", tuple(self.covariate_names))

    @property
    def n_samples(self) -> int:
        return int(self.covariates.shape[0])

    @property
    def column_names(self) -> tuple[str, ...]:
        """Design column names beyond the intercept: conditioning, then covariates."""

        return self.conditioning_names + self.covariate_names

    def matrix(self) -> np.ndarray:
        """The full residualization design ``[1, conditioning, covariates]``."""

        return np.column_stack([np.ones(self.n_samples), self.conditioning, self.covariates])

    def cells(self) -> np.ndarray:
        """One label per row naming its conditioning cell (its distinct conditioning row).

        Under stage indicators these are the stages; under group × stage
        indicators, the group × stage cells. Used to stratify cross-validation
        folds so that every training set contains every conditioning cell.
        """

        return np.array(["|".join(f"{value:g}" for value in row) for row in self.conditioning], dtype=object)

    def rows(self, index: Union[np.ndarray, Sequence[int]]) -> "NuisanceDesign":
        """The design restricted to ``index`` (positional rows)."""

        index = np.asarray(index)
        return NuisanceDesign(
            self.conditioning[index], self.covariates[index], self.conditioning_names, self.covariate_names, self.kind
        )


def fit_nuisance_coefficients(
    X: Union[np.ndarray, pd.DataFrame],
    design: NuisanceDesign,
) -> tuple[np.ndarray, int]:
    """Estimate the covariate coefficients of ``X`` conditional on the design.

    Regresses every column of ``X`` on ``[1, conditioning, covariates]`` by
    least squares and returns the covariate block of the coefficients
    (``q × p``) and the residual degrees of freedom ``n − 1 − k − q``.

    Raises
    ------
    ValueError
        If ``X`` and the design differ in rows, or the design does not have
        full column rank; the message names the design columns.
    """

    X_arr = np.asarray(X, dtype=float)
    if X_arr.shape[0] != design.n_samples:
        raise ValueError(f"X has {X_arr.shape[0]} rows but the nuisance design has {design.n_samples}.")
    M = design.matrix()
    rank = int(np.linalg.matrix_rank(M))
    if rank < M.shape[1]:
        raise ValueError(
            f"The residualization design [intercept, {', '.join(design.column_names)}] is not of full rank "
            f"(rank {rank} of {M.shape[1]} columns on {M.shape[0]} samples); a covariate is collinear with "
            "the conditioning design or with another covariate, or a factor level is absent."
        )
    coefficients = np.linalg.lstsq(M, X_arr, rcond=None)[0]
    k = design.conditioning.shape[1]
    return coefficients[1 + k :, :], int(M.shape[0] - M.shape[1])


@overload
def residualize(X: pd.DataFrame, covariates: np.ndarray, coefficients: np.ndarray) -> pd.DataFrame: ...


@overload
def residualize(X: np.ndarray, covariates: np.ndarray, coefficients: np.ndarray) -> np.ndarray: ...


def residualize(
    X: Union[np.ndarray, pd.DataFrame],
    covariates: np.ndarray,
    coefficients: np.ndarray,
) -> Union[np.ndarray, pd.DataFrame]:
    """Subtract the fitted covariate part: ``X − covariates @ coefficients``.

    Only the covariate columns are read, so held-out rows can be residualized
    with training coefficients without their labels. A DataFrame keeps its
    index and columns.
    """

    shift = np.asarray(covariates, dtype=float) @ np.asarray(coefficients, dtype=float)
    if isinstance(X, pd.DataFrame):
        return X - pd.DataFrame(shift, index=X.index, columns=X.columns)
    return np.asarray(X, dtype=float) - shift


def nuisance_partial_r2(X: Union[np.ndarray, pd.DataFrame], design: NuisanceDesign) -> np.ndarray:
    """Per-feature partial R² of the covariates given the intercept and conditioning.

    ``1 − RSS([1, D, C]) / RSS([1, D])`` for every column of ``X``: the share
    of the variance left after the conditioning design that the covariates
    explain, i.e. what the residualization removes. A column with no variance
    left after the conditioning design gets 0.
    """

    X_arr = np.asarray(X, dtype=float)
    full = design.matrix()
    reduced = full[:, : 1 + design.conditioning.shape[1]]
    rss_full = _residual_sum_of_squares(full, X_arr)
    rss_reduced = _residual_sum_of_squares(reduced, X_arr)
    with np.errstate(divide="ignore", invalid="ignore"):
        r2 = np.where(rss_reduced > 0.0, 1.0 - rss_full / rss_reduced, 0.0)
    return np.clip(r2, 0.0, 1.0)


def _residual_sum_of_squares(M: np.ndarray, X: np.ndarray) -> np.ndarray:
    coefficients = np.linalg.lstsq(M, X, rcond=None)[0]
    residuals = X - M @ coefficients
    return np.einsum("ij,ij->j", residuals, residuals)


def _fit_and_residualize(
    X_train: pd.DataFrame,
    X_held_out: pd.DataFrame,
    design: NuisanceDesign,
    train: np.ndarray,
    held_out: np.ndarray,
    where: str,
) -> tuple[pd.DataFrame, pd.DataFrame, int]:
    """Residualize a training set and its held-out rows with the training coefficients."""

    train_design = design.rows(train)
    try:
        coefficients, residual_df = fit_nuisance_coefficients(X_train, train_design)
    except ValueError as exc:
        raise ValueError(f"Nuisance residualization failed in {where}: {exc}") from exc
    return (
        residualize(X_train, train_design.covariates, coefficients),
        residualize(X_held_out, design.covariates[held_out], coefficients),
        residual_df,
    )


def plsda_doubleCV(
    X: pd.DataFrame,
    y: Union[pd.DataFrame, pd.Series],
    cv1_splits: int = 7,
    cv2_splits: int = 8,
    n_repeats: int = 30,
    max_components: int = 50,
    random_state: int = 1203,
    n_jobs: int = 1,
    progress: bool = True,
    *,
    scale: bool = True,
    nuisance: NuisanceDesign | None = None,
    strata: Union[np.ndarray, pd.Series, Sequence[str], None] = None,
    min_components: int = 1,
) -> PLSDAResult:
    """
    Run canonical double-nested cross-validation for PLS-DA.

    For each outer fold, the inner CV averages AUROC across its V folds and
    selects the n_LV with the highest mean. The outer fold's selected n_LV is
    then evaluated on its held-out test set. Per repeat, the K outer-fold test
    AUROCs are aggregated by mean (and sample std), and the per-fold n_LV
    choices are aggregated by mode (parsimony tie-break). One final model is
    refit per repeat on the full input using the modal n_LV.

    Parameters
    ----------
    X: pd.DataFrame
        The predictor variables.
    y: Union[pd.DataFrame, pd.Series]
        The outcome variable.
    cv1_splits: int
        Number of folds in the CV1 (inner) loop. Default: 7.
    cv2_splits: int
        Number of folds in the CV2 (outer) loop. Default: 8.
    n_repeats: int
        Number of repeats of the K-fold outer CV. Default: 30.
    max_components: int
        Maximum number of LV to test (candidates are
        ``min_components..max_components-1``). Default: 50.
    random_state: int
        For reproducibility. Default: 1203.
    n_jobs: int
        Number of parallel workers for the inner CV loop. Use -1 for all
        available CPUs. Default: 1 (serial).
    progress: bool
        Whether to display a tqdm progress bar over outer folds. Default: True.
    scale: bool
        Keyword-only. Whether every PLS fit (inner, outer and the per-repeat
        refits) rescales the predictors to unit variance. ``False`` fits in the
        metric ``X`` is already in. Default: True.
    nuisance: NuisanceDesign | None
        Keyword-only. When given (row-aligned to ``X``), every training set —
        each outer rest set and each inner training set — estimates the
        covariate coefficients of ``X`` on ``[1, conditioning, covariates]``
        from its own rows (:func:`fit_nuisance_coefficients`), residualizes its
        own rows, and residualizes the held-out rows with the same coefficients
        and the held-out covariate columns alone, so held-out labels never
        enter the residualization. The per-repeat refits use coefficients
        estimated on all rows. Every training design must have full column
        rank; otherwise a ``ValueError`` names the repeat and fold. Default:
        None (no residualization).
    strata: array-like | None
        Keyword-only. Row-aligned labels the outer and inner folds are
        stratified on instead of ``y`` (for example the conditioning cells of
        ``nuisance``, so that every training set contains every cell). Folds
        depend only on the partition the labels define, so labels that
        partition the rows like ``y`` reproduce the default folds. Default:
        None (stratify on ``y``).
    min_components: int
        Keyword-only. The smallest number of latent variables the inner CV
        considers, so the candidates are ``min_components..max_components-1``
        and every outer-fold selection, every modal ``LV`` and every refit
        model has at least this many components. Must satisfy
        ``1 <= min_components < max_components``; otherwise a ``ValueError``
        names both values before any fold is fitted. Default: 1 (the
        unfloored rule; the trajectory latent space applies its own floor,
        see :func:`motco.simulations.evaluation.fit_pls_latent_space`).

    Returns
    -------
    dict
        With two keys:

        - ``"table"``: pd.DataFrame with one row per repeat and four columns —
          ``"rep"`` (int), ``"LV"`` (int, modal across K outer folds with
          parsimony tie-break), ``"AUROC"`` (float, mean across K outer folds),
          ``"AUROC_std"`` (float, sample std across K outer folds; NaN if K < 2).
        - ``"models"``: list of length ``n_repeats``. Each entry is a
          PLSRegression refit on the full input (X, one-hot(y)) with the
          corresponding row's modal n_LV.

        With ``nuisance`` given, a third key ``"nuisance_min_training_df"``
        holds the smallest residual degrees of freedom of any training-set
        residualization fit.
    """
    _X_arr = np.asarray(X, dtype=float)
    _y_arr = np.asarray(y)
    _n_x = _X_arr.shape[0]
    _n_y = _y_arr.shape[0]
    if _n_x != _n_y:
        raise ValueError(
            f"X has {_n_x} rows but y has {_n_y} rows — number of rows must match."
        )
    _n_classes = len(np.unique(_y_arr))
    if _n_classes < 2:
        raise ValueError(
            f"y has {_n_classes} unique class(es); at least 2 are required."
        )
    if max_components > _X_arr.shape[1]:
        raise ValueError(
            f"max_components={max_components} exceeds the number of features in X "
            f"({_X_arr.shape[1]}); reduce max_components."
        )
    try:
        floor_ok = int(min_components) == min_components and 1 <= min_components < max_components
    except (TypeError, ValueError):
        floor_ok = False
    if not floor_ok:
        raise ValueError(
            f"min_components={min_components!r} must be an integer with "
            f"1 <= min_components < max_components={max_components}; the candidates are "
            "min_components..max_components-1."
        )
    min_components = int(min_components)
    if not np.all(np.isfinite(_X_arr)):
        raise ValueError("X contains NaN or Inf values.")
    if nuisance is not None and nuisance.n_samples != _n_x:
        raise ValueError(f"X has {_n_x} rows but the nuisance design has {nuisance.n_samples}.")
    folds_on: Union[pd.DataFrame, pd.Series] = y
    if strata is not None:
        folds_on = pd.Series(np.asarray(strata))
        if folds_on.shape[0] != _n_x:
            raise ValueError(f"X has {_n_x} rows but strata has {folds_on.shape[0]}.")
    encoder = OneHotEncoder(sparse_output=False)
    yd_arr = encoder.fit_transform(np.array(y).reshape(-1, 1))
    yd = pd.DataFrame(yd_arr)
    cv2 = RepeatedStratifiedKFold(
        n_splits=cv2_splits, n_repeats=n_repeats, random_state=random_state
    )
    cv1 = StratifiedKFold(n_splits=cv1_splits)

    n_lv_candidates = list(range(min_components, max_components))
    n_candidates = len(n_lv_candidates)

    rep_means: list[float] = []
    rep_stds: list[float] = []
    rep_lvs: list[int] = []

    outer_aurocs: list[float] = []
    outer_n_lv: list[int] = []
    fold_in_repeat = 0
    training_dfs: list[int] = []
    auroc = partial(_plsda_auroc, scale=scale)

    cv2_iter = tqdm(
        cv2.split(X, folds_on),
        total=cv2_splits * n_repeats,
        desc="PLS-DA CV2",
        unit="fold",
        disable=not progress,
    )
    for rest, test in cv2_iter:
        X_rest = X.iloc[rest, :]
        X_test = X.iloc[test, :]
        folds_on_rest = folds_on.iloc[rest]
        yd_rest = yd.iloc[rest, :]
        yd_test = yd.iloc[test, :]
        outer_where = f"repeat {len(rep_lvs) + 1}, outer fold {fold_in_repeat + 1}"
        X_rest_fit, X_test_fit = X_rest, X_test
        if nuisance is not None:
            X_rest_fit, X_test_fit, rest_df = _fit_and_residualize(
                X_rest, X_test, nuisance, rest, test, outer_where
            )
            training_dfs.append(rest_df)

        # Inner CV: average AUROC across V folds for each candidate n_LV.
        inner_aurocs = np.zeros((cv1_splits, n_candidates))
        for v_idx, (train, validation) in enumerate(cv1.split(X_rest, folds_on_rest)):
            X_train = X_rest.iloc[train, :]
            yd_train = yd_rest.iloc[train, :]
            X_val = X_rest.iloc[validation, :]
            yd_val = yd_rest.iloc[validation, :]
            if nuisance is not None:
                # Coefficients from the inner training rows of the raw predictors.
                X_train, X_val, train_df = _fit_and_residualize(
                    X_train,
                    X_val,
                    nuisance,
                    rest[train],
                    rest[validation],
                    f"{outer_where}, inner fold {v_idx + 1}",
                )
                training_dfs.append(train_df)
            args = zip(
                n_lv_candidates,
                repeat(X_train),
                repeat(yd_train),
                repeat(X_val),
                repeat(yd_val),
            )
            fold_aurocs: list[float]
            if n_jobs == 1:
                fold_aurocs = [auroc(*a) for a in args]  # type: ignore[misc]
            else:
                n_workers = (
                    (multiprocessing.cpu_count() or 1) if n_jobs == -1 else max(1, n_jobs)
                )
                with multiprocessing.Pool(processes=n_workers) as pool:
                    fold_aurocs = pool.starmap(auroc, args)  # type: ignore[arg-type]
            inner_aurocs[v_idx, :] = fold_aurocs

        mean_per_n_lv = inner_aurocs.mean(axis=0)
        # np.argmax returns the first occurrence on ties → parsimony. Select
        # by value so the floor is a property of the candidate set.
        n_lv_star = n_lv_candidates[int(np.argmax(mean_per_n_lv))]

        outer_score = float(
            auroc(n_lv_star, X_rest_fit, yd_rest, X_test_fit, yd_test)  # type: ignore[arg-type]
        )
        outer_aurocs.append(outer_score)
        outer_n_lv.append(n_lv_star)

        fold_in_repeat += 1
        if fold_in_repeat == cv2_splits:
            mean_auroc = float(np.mean(outer_aurocs))
            std_auroc = (
                float(np.std(outer_aurocs, ddof=1))
                if len(outer_aurocs) > 1
                else float("nan")
            )
            mode_lv = _modal_int_with_parsimony(outer_n_lv)
            rep_means.append(mean_auroc)
            rep_stds.append(std_auroc)
            rep_lvs.append(mode_lv)
            outer_aurocs = []
            outer_n_lv = []
            fold_in_repeat = 0

    model_table = pd.DataFrame(
        {
            "rep": list(range(1, n_repeats + 1)),
            "LV": rep_lvs,
            "AUROC": rep_means,
            "AUROC_std": rep_stds,
        }
    )
    model_table["rep"] = model_table["rep"].astype(int)
    model_table["LV"] = model_table["LV"].astype(int)
    model_table["AUROC"] = model_table["AUROC"].astype(float)
    model_table["AUROC_std"] = model_table["AUROC_std"].astype(float)

    X_full = _X_arr
    if nuisance is not None:
        try:
            full_coefficients, _ = fit_nuisance_coefficients(_X_arr, nuisance)
        except ValueError as exc:
            raise ValueError(f"Nuisance residualization failed on the full input: {exc}") from exc
        X_full = residualize(_X_arr, nuisance.covariates, full_coefficients)

    best_models: list[PLSRegression] = []
    for lv in rep_lvs:
        model = PLSRegression(
            n_components=lv, scale=scale, max_iter=1000
        ).fit(X_full, yd_arr)
        best_models.append(model)

    if nuisance is not None:
        return {"models": best_models, "table": model_table, "nuisance_min_training_df": min(training_dfs)}
    return {"models": best_models, "table": model_table}


def _modal_int_with_parsimony(values: list[int]) -> int:
    """Return the most frequent integer; ties broken by smaller value."""
    if not values:
        raise ValueError("values must be non-empty")
    counts: dict[int, int] = {}
    for v in values:
        counts[v] = counts.get(v, 0) + 1
    max_count = max(counts.values())
    return min(v for v, c in counts.items() if c == max_count)


def _plsda_auroc(
    n_components: int,
    X_train: pd.DataFrame,
    Y_train: Union[pd.Series, pd.DataFrame],
    X_test: pd.DataFrame,
    Y_test: Union[pd.Series, pd.DataFrame],
    return_full: bool = False,
    *,
    scale: bool = True,
) -> Union[float, dict[str, Union[PLSRegression, float]]]:
    """
    Estimate a partial least squares regression and return the AUROC value
    and the model, or just the AUROC value.

    Parameters
    ----------
    n_components: int
        Number of components to use.
    X_train: pd.DataFrame
        The predictors to use for training.
    Y_train: Union[pd.Series, pd.DataFrame]
        The outcome to use for training.
    X_test: pd.DataFrame,
        The predictors to use for testing.
    Y_test: Union[pd.Series, pd.DataFrame]
        The outcomes to use for testing.
    return_full: bool
        Whether to return the model and the auroc score.
        If False, returns only the auroc score. Default: False.
    scale: bool
        Keyword-only. Whether the PLS fit rescales the predictors to unit
        variance. Default: True.

    Returns
    -------
    auroc: Union[float,
                 dict[str,
                      Union[PLSRegression,
                            float]]]
        Return the auroc score, and optionally the regression model.
    """
    pls = PLSRegression(n_components=n_components, scale=scale, max_iter=1000).fit(
        X=X_train, y=Y_train
    )
    y_pred = pls.predict(X_test)
    score = roc_auc_score(Y_test, y_pred)
    if return_full:
        auroc = {"model": pls, "score": score}
    else:
        auroc = score
    return auroc


def fit_plsda_model(
    X: Union[np.ndarray, "pd.DataFrame"],
    y: Union["pd.Series", np.ndarray],
    n_components: int,
    *,
    scale: bool = True,
) -> PLSRegression:
    """Fit a PLS-DA model on the full dataset and return the fitted estimator.

    No cross-validation: the model is fit once on all rows with a fixed number
    of components. Use this when you need the projector itself (e.g. to call
    ``.transform()`` on new vectors), rather than just the training scores.

    Parameters
    ----------
    X:
        Feature matrix of shape (n_samples, n_features).
    y:
        Class labels of shape (n_samples,). String or numeric; one-hot encoded
        internally to match the behaviour of ``plsda_doubleCV``.
    n_components:
        Number of latent variables.
    scale:
        Whether ``PLSRegression`` rescales each feature to unit variance before
        fitting (the default). Pass ``False`` to fit in the metric ``X`` is
        already in.

    Returns
    -------
    PLSRegression
        Fitted estimator whose ``.transform()`` projects feature vectors onto
        the ``n_components`` latent axes.
    """
    X_arr = np.asarray(X, dtype=float)
    encoder = OneHotEncoder(sparse_output=False)
    y_enc = encoder.fit_transform(np.asarray(y).reshape(-1, 1))
    return PLSRegression(n_components=n_components, scale=scale, max_iter=1000).fit(X_arr, y_enc)


def fit_plsda_transform(
    X: Union[np.ndarray, "pd.DataFrame"],
    y: Union["pd.Series", np.ndarray],
    n_components: int,
) -> np.ndarray:
    """Fit a PLS model on the full dataset and return the X score matrix.

    Parameters
    ----------
    X:
        Feature matrix of shape (n_samples, n_features).
    y:
        Class labels of shape (n_samples,). String or numeric; one-hot encoded
        internally to match the behaviour of ``plsda_doubleCV``.
    n_components:
        Number of latent variables.

    Returns
    -------
    scores : np.ndarray of shape (n_samples, n_components)
    """
    return np.asarray(fit_plsda_model(X, y, n_components).x_scores_)


def calculate_vips(
    model,
    components: Union[None, list[int]] = None,
) -> np.ndarray:
    """
    Estimates Variable Importance in Projection (VIP)
    in Partial Least Squares (PLS)

    Parameters
    ----------
    model: PLSRegression
        model generated from the PLSRegression function
    components: Union[None, list[int]]
        if not None, a list of integers indicating the components to compute
        the VIPs from. If None, all components are taken into account.
        Default None.

    Returns
    -------
    vips: np.array
        variable importance in projection for each variable
    """
    if components is not None:
        t = model.x_scores_[:, components]
        w = model.x_weights_[:, components]
        q = model.y_loadings_[:, components]
        p, h = w.shape
    else:
        w = model.x_weights_
        t = model.x_scores_
        q = model.y_loadings_
        p, h = w.shape

    vips = np.zeros((p,))
    s = np.diag(np.matmul(np.matmul(np.matmul(t.T, t), q.T), q)).reshape(h, -1)
    total_s = np.sum(s)
    for i in range(p):
        weight = np.array([(w[i, j] / np.linalg.norm(w[:, j])) ** 2 for j in range(h)])
        vips[i] = np.sqrt(p * np.matmul(s.T, weight)[0] / total_s)
    return vips
