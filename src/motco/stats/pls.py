"""Partial Least Squares utilities (PLS-DA, VIP computation).

This module provides:
- `plsda_doubleCV`: double cross-validation over number of components
  for PLS-DA classification using AUROC as criterion.
- `calculate_vips`: compute Variable Importance in Projection (VIP) scores
  from a fitted `sklearn.cross_decomposition.PLSRegression` model.
"""

import logging
import multiprocessing
from itertools import repeat
from typing import TypedDict, Union

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
        Maximum number of LV to test (candidates are 1..max_components-1). Default: 50.
    random_state: int
        For reproducibility. Default: 1203.
    n_jobs: int
        Number of parallel workers for the inner CV loop. Use -1 for all
        available CPUs. Default: 1 (serial).
    progress: bool
        Whether to display a tqdm progress bar over outer folds. Default: True.

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
    if not np.all(np.isfinite(_X_arr)):
        raise ValueError("X contains NaN or Inf values.")
    encoder = OneHotEncoder(sparse_output=False)
    yd_arr = encoder.fit_transform(np.array(y).reshape(-1, 1))
    yd = pd.DataFrame(yd_arr)
    cv2 = RepeatedStratifiedKFold(
        n_splits=cv2_splits, n_repeats=n_repeats, random_state=random_state
    )
    cv1 = StratifiedKFold(n_splits=cv1_splits)

    n_lv_candidates = list(range(1, max_components))
    n_candidates = len(n_lv_candidates)

    rep_means: list[float] = []
    rep_stds: list[float] = []
    rep_lvs: list[int] = []

    outer_aurocs: list[float] = []
    outer_n_lv: list[int] = []
    fold_in_repeat = 0

    cv2_iter = tqdm(
        cv2.split(X, y),
        total=cv2_splits * n_repeats,
        desc="PLS-DA CV2",
        unit="fold",
        disable=not progress,
    )
    for rest, test in cv2_iter:
        X_rest = X.iloc[rest, :]
        X_test = X.iloc[test, :]
        y_rest = y.iloc[rest]
        yd_rest = yd.iloc[rest, :]
        yd_test = yd.iloc[test, :]

        # Inner CV: average AUROC across V folds for each candidate n_LV.
        inner_aurocs = np.zeros((cv1_splits, n_candidates))
        for v_idx, (train, validation) in enumerate(cv1.split(X_rest, y_rest)):
            X_train = X_rest.iloc[train, :]
            yd_train = yd_rest.iloc[train, :]
            X_val = X_rest.iloc[validation, :]
            yd_val = yd_rest.iloc[validation, :]
            args = zip(
                n_lv_candidates,
                repeat(X_train),
                repeat(yd_train),
                repeat(X_val),
                repeat(yd_val),
            )
            fold_aurocs: list[float]
            if n_jobs == 1:
                fold_aurocs = [_plsda_auroc(*a) for a in args]  # type: ignore[misc]
            else:
                n_workers = (
                    (multiprocessing.cpu_count() or 1) if n_jobs == -1 else max(1, n_jobs)
                )
                with multiprocessing.Pool(processes=n_workers) as pool:
                    fold_aurocs = pool.starmap(_plsda_auroc, args)  # type: ignore[arg-type]
            inner_aurocs[v_idx, :] = fold_aurocs

        mean_per_n_lv = inner_aurocs.mean(axis=0)
        # np.argmax returns the first occurrence on ties → parsimony.
        n_lv_star = int(np.argmax(mean_per_n_lv)) + 1

        outer_score = float(
            _plsda_auroc(n_lv_star, X_rest, yd_rest, X_test, yd_test)  # type: ignore[arg-type]
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

    best_models: list[PLSRegression] = []
    for lv in rep_lvs:
        model = PLSRegression(
            n_components=lv, scale=True, max_iter=1000
        ).fit(_X_arr, yd_arr)
        best_models.append(model)

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

    Returns
    -------
    auroc: Union[float,
                 dict[str,
                      Union[PLSRegression,
                            float]]]
        Return the auroc score, and optionally the regression model.
    """
    pls = PLSRegression(n_components=n_components, scale=True, max_iter=1000).fit(
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
