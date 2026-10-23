# tests/test_pls.py
from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
import pytest
from sklearn.model_selection import RepeatedStratifiedKFold

from motco.stats import pls
from motco.stats.pls import (
    NuisanceDesign,
    calculate_vips,
    fit_nuisance_coefficients,
    fit_plsda_model,
    fit_plsda_transform,
    nuisance_partial_r2,
    plsda_doubleCV,
    residualize,
)
from tests.test_phase6_byte_identity import _assert_close

#: ``plsda_doubleCV`` output on a fixed input and ``random_state``, captured
#: before the keyword-only ``scale`` / ``nuisance`` arguments existed
#: (``uv run python -m tests.test_pls``).
DEFAULT_IDENTITY = Path(__file__).parent / "data" / "plsda_doublecv_identity.json"


def _synthetic_data(n: int = 30, p: int = 10, seed: int = 0) -> tuple[pd.DataFrame, pd.Series]:
    rng = np.random.default_rng(seed)
    X = pd.DataFrame(rng.standard_normal((n, p)), columns=[f"f{i}" for i in range(p)])
    y = pd.Series(["A"] * (n // 2) + ["B"] * (n - n // 2))
    return X, y


def test_plsda_returns_required_keys():
    X, y = _synthetic_data()
    result = plsda_doubleCV(X, y, cv1_splits=3, cv2_splits=3, n_repeats=2, max_components=3)
    assert set(result.keys()) >= {"table", "models"}


def test_plsda_table_has_correct_shape():
    n_repeats = 3
    X, y = _synthetic_data()
    result = plsda_doubleCV(X, y, cv1_splits=3, cv2_splits=3, n_repeats=n_repeats, max_components=3)
    assert result["table"].shape == (n_repeats, 4)
    assert list(result["table"].columns) == ["rep", "LV", "AUROC", "AUROC_std"]


def test_plsda_auroc_values_in_unit_interval():
    X, y = _synthetic_data()
    result = plsda_doubleCV(X, y, cv1_splits=3, cv2_splits=3, n_repeats=2, max_components=3)
    aurocs = result["table"].iloc[:, 2].values  # AUROC is third column
    assert np.all(aurocs >= 0.0) and np.all(aurocs <= 1.0)


def test_plsda_auroc_std_non_negative():
    X, y = _synthetic_data()
    result = plsda_doubleCV(X, y, cv1_splits=3, cv2_splits=3, n_repeats=2, max_components=3)
    stds = result["table"]["AUROC_std"].values
    finite = stds[~np.isnan(stds)]
    assert np.all(finite >= 0.0)


def test_plsda_doubleCV_mean_aggregation_is_honest():
    """Synthetic data with moderate SNR: mean AUROC should be < 0.99 with K outer folds.

    Group A ~ N(0, I_20); Group B ~ N(mu, I_20) with mu = 0.5 in 20 features.
    Bayes-optimal AUROC ≈ 0.94; small-n CV variance produces visible per-fold spread.
    Under the old max-of-K aggregation, at least one fold per repeat would saturate
    at AUROC = 1.0, pinning the reported value. Mean-of-K reporting produces an
    honest summary.
    """
    rng = np.random.default_rng(42)
    n_per = 30
    X_A = rng.standard_normal((n_per, 20))
    X_B = rng.standard_normal((n_per, 20)) + 0.5
    X = pd.DataFrame(np.vstack([X_A, X_B]), columns=[f"f{i}" for i in range(20)])
    y = pd.Series(["A"] * n_per + ["B"] * n_per)
    result = plsda_doubleCV(X, y, cv1_splits=3, cv2_splits=4, n_repeats=3, max_components=5)
    aurocs = result["table"]["AUROC"].values
    auroc_stds = result["table"]["AUROC_std"].values
    assert np.mean(aurocs) < 0.99, (
        f"Mean AUROC {np.mean(aurocs):.4f} is suspiciously close to 1.0; "
        "max-of-K aggregation may have leaked through."
    )
    assert np.mean(auroc_stds) > 0.0, (
        "AUROC_std is zero across repeats; expected fold-to-fold variance "
        "with moderate-SNR data."
    )


def test_plsda_doubleCV_models_refit_on_full_data():
    """Models stored per repeat are full-data refits, not CV-fold models."""
    n_repeats = 2
    X, y = _synthetic_data(n=30, p=10)
    result = plsda_doubleCV(
        X, y, cv1_splits=3, cv2_splits=3, n_repeats=n_repeats, max_components=3
    )
    for i in range(n_repeats):
        model = result["models"][i]
        assert model.x_scores_.shape[0] == X.shape[0], (
            f"models[{i}].x_scores_.shape[0]={model.x_scores_.shape[0]}, "
            f"expected full-data refit on {X.shape[0]} samples"
        )
        # Per-repeat n_components should match the modal LV recorded in the table.
        assert model.n_components == int(result["table"].iloc[i, 1])


def test_plsda_model_count_matches_repeats():
    n_repeats = 2
    X, y = _synthetic_data()
    result = plsda_doubleCV(X, y, cv1_splits=3, cv2_splits=3, n_repeats=n_repeats, max_components=3)
    assert len(result["models"]) == n_repeats


def test_plsda_lv_values_are_positive_integers():
    X, y = _synthetic_data()
    result = plsda_doubleCV(X, y, cv1_splits=3, cv2_splits=3, n_repeats=2, max_components=3)
    lvs = result["table"].iloc[:, 1].values  # LV is second column
    assert np.all(lvs >= 1)
    assert np.all(lvs == lvs.astype(int))


def test_calculate_vips_shape():
    X, y = _synthetic_data()
    result = plsda_doubleCV(X, y, cv1_splits=3, cv2_splits=3, n_repeats=2, max_components=3)
    model = result["models"][0]
    vips = calculate_vips(model)
    assert vips.shape == (X.shape[1],)  # one VIP per feature


def test_calculate_vips_non_negative():
    X, y = _synthetic_data()
    result = plsda_doubleCV(X, y, cv1_splits=3, cv2_splits=3, n_repeats=2, max_components=3)
    model = result["models"][0]
    vips = calculate_vips(model)
    assert np.all(vips >= 0)


def test_fit_plsda_transform_shape():
    X, y = _synthetic_data(n=30, p=10)
    scores = fit_plsda_transform(X, y, n_components=3)
    assert scores.shape == (30, 3)


def test_fit_plsda_transform_single_component():
    X, y = _synthetic_data(n=30, p=10)
    scores = fit_plsda_transform(X, y, n_components=1)
    assert scores.shape == (30, 1)


def test_fit_plsda_transform_returns_ndarray():
    X, y = _synthetic_data(n=30, p=10)
    scores = fit_plsda_transform(X, y, n_components=2)
    assert isinstance(scores, np.ndarray)


def test_fit_plsda_transform_numeric_labels():
    rng = np.random.default_rng(0)
    X = pd.DataFrame(rng.standard_normal((30, 8)), columns=[f"f{i}" for i in range(8)])
    y = pd.Series([0, 1, 2] * 10)
    scores = fit_plsda_transform(X, y, n_components=2)
    assert scores.shape == (30, 2)


def test_fit_plsda_model_returns_fitted_projector():
    from sklearn.cross_decomposition import PLSRegression

    X, y = _synthetic_data(n=30, p=10)
    model = fit_plsda_model(X, y, n_components=2)
    assert isinstance(model, PLSRegression)
    # The projector transforms new feature vectors into the 2-D latent space.
    projected = model.transform(X.to_numpy())
    assert projected.shape == (30, 2)
    # fit_plsda_transform is a thin wrapper returning the training scores.
    np.testing.assert_allclose(projected, fit_plsda_transform(X, y, n_components=2))


def test_stats_top_level_imports():
    from motco.stats import (  # noqa: F401
        RRPP,
        SNF,
        build_ls_means,
        calculate_vips,
        center_matrix,
        estimate_betas,
        estimate_difference,
        fit_plsda_model,
        fit_plsda_transform,
        get_affinity_matrix,
        get_model_matrix,
        get_observed_vectors,
        get_spectral,
        pair_difference,
        plsda_doubleCV,
    )


def _default_identity_payload() -> dict[str, Any]:
    rng = np.random.default_rng(20261020)
    labels = np.repeat(["Low", "Mid", "High"], 12)
    shift = np.repeat([0.0, 0.8, 1.6], 12)[:, None]
    X = pd.DataFrame(rng.standard_normal((36, 9)) + shift * np.linspace(0, 1, 9), columns=[f"f{i}" for i in range(9)])
    result = plsda_doubleCV(
        X, pd.Series(labels), cv1_splits=2, cv2_splits=3, n_repeats=2, max_components=4, random_state=11,
        progress=False,
    )
    return json.loads(
        json.dumps(
            {
                "table": result["table"].to_dict(orient="list"),
                "x_scores": [np.asarray(model.x_scores_).tolist() for model in result["models"]],
                "x_rotations": [np.asarray(model.x_rotations_).tolist() for model in result["models"]],
            }
        )
    )


def test_plsda_doubleCV_defaults_reproduce_stored_output():
    expected = json.loads(DEFAULT_IDENTITY.read_text(encoding="utf-8"))
    _assert_close(_default_identity_payload(), expected, "plsda_doubleCV")



# ── Nuisance residualization ─────────────────────────────────────────────────


def _stage_indicators(stages: np.ndarray) -> np.ndarray:
    return np.column_stack([(stages == level).astype(float) for level in (1, 2)])


def _nuisance_problem(seed: int = 0, n_per_stage: int = 12, p: int = 8):
    """Features carrying a stage effect plus a covariate effect, the covariate correlated with stage."""

    rng = np.random.default_rng(seed)
    stages = np.repeat([0, 1, 2], n_per_stage)
    z = stages + rng.normal(scale=0.8, size=stages.size)
    stage_effect = np.array([0.0, 1.0, 2.5])[stages]
    v = rng.normal(size=p)
    X = np.outer(stage_effect, np.linspace(0.5, 1.5, p)) + np.outer(z, v) + rng.normal(scale=0.5, size=(stages.size, p))
    design = NuisanceDesign(_stage_indicators(stages), z[:, None], ("stage[1]", "stage[2]"), ("z",))
    labels = pd.Series(np.array(["Low", "Mid", "High"])[stages])
    return pd.DataFrame(X, columns=[f"f{i}" for i in range(p)]), labels, design, stages, z


def test_stage_collinear_covariate_keeps_the_stage_effect():
    rng = np.random.default_rng(3)
    stages = np.repeat([0, 1, 2], 15)
    z = stages + rng.normal(scale=0.6, size=stages.size)
    stage_effect = np.array([0.0, 1.0, 2.5])[stages]
    feature = stage_effect + 2.0 * z
    design = NuisanceDesign(_stage_indicators(stages), z[:, None], ("stage[1]", "stage[2]"), ("z",))
    coefficients, residual_df = fit_nuisance_coefficients(feature[:, None], design)
    assert coefficients.shape == (1, 1)
    assert coefficients[0, 0] == pytest.approx(2.0, abs=1e-10)  # within-stage effect
    marginal = np.polyfit(z, feature, 1)[0]
    assert abs(marginal - 2.0) > 0.3  # the marginal slope absorbs the stage effect
    assert residual_df == stages.size - 1 - 2 - 1
    kept = residualize(feature[:, None], design.covariates, coefficients)[:, 0]
    np.testing.assert_allclose(kept, stage_effect, atol=1e-10)


def test_covariate_effect_is_removed_within_the_design():
    X, _, design, stages, z = _nuisance_problem()
    coefficients, _ = fit_nuisance_coefficients(X, design)
    X_res = np.asarray(residualize(X.to_numpy(), design.covariates, coefficients))
    base = np.column_stack([np.ones(stages.size), design.conditioning])
    within = np.eye(stages.size) - base @ np.linalg.pinv(base)
    z_within = within @ z
    X_within = within @ X_res
    correlation = z_within @ X_within / (np.linalg.norm(z_within) * np.linalg.norm(X_within, axis=0))
    np.testing.assert_allclose(correlation, 0.0, atol=1e-10)
    assert np.all(np.abs(np.corrcoef(z, X.to_numpy(), rowvar=False)[0, 1:]) > 0.1)


def test_residualize_keeps_dataframe_labels():
    X, _, design, _, _ = _nuisance_problem()
    coefficients, _ = fit_nuisance_coefficients(X, design)
    out = residualize(X, design.covariates, coefficients)
    assert isinstance(out, pd.DataFrame)
    pd.testing.assert_index_equal(out.columns, X.columns)
    pd.testing.assert_index_equal(out.index, X.index)


def test_rank_deficient_nuisance_design_names_the_columns():
    stages = np.repeat([0, 1, 2], 6)
    design = NuisanceDesign(
        _stage_indicators(stages), stages[:, None].astype(float), ("stage[1]", "stage[2]"), ("stage_index",)
    )
    X = np.random.default_rng(0).normal(size=(stages.size, 3))
    with pytest.raises(ValueError, match=r"\[intercept, stage\[1\], stage\[2\], stage_index\] is not of full rank"):
        fit_nuisance_coefficients(X, design)


def test_partial_r2_is_the_covariate_share_after_the_conditioning():
    X, _, design, stages, z = _nuisance_problem()
    r2 = nuisance_partial_r2(X, design)
    base = np.column_stack([np.ones(stages.size), design.conditioning])
    full = np.column_stack([base, z])
    for j in range(X.shape[1]):
        x = X.iloc[:, j].to_numpy()
        rss_base = np.sum((x - base @ np.linalg.lstsq(base, x, rcond=None)[0]) ** 2)
        rss_full = np.sum((x - full @ np.linalg.lstsq(full, x, rcond=None)[0]) ** 2)
        assert r2[j] == pytest.approx(1.0 - rss_full / rss_base, abs=1e-12)
    assert np.all((r2 >= 0.0) & (r2 <= 1.0))


def test_plsda_doubleCV_default_scale_is_true():
    X, y = _synthetic_data()
    result = plsda_doubleCV(X, y, cv1_splits=3, cv2_splits=3, n_repeats=2, max_components=3, progress=False)
    assert all(model.scale is True for model in result["models"])
    assert "nuisance_min_training_df" not in result


def test_plsda_doubleCV_scale_false_and_nuisance_refits():
    X, y, design, _, _ = _nuisance_problem()
    result = plsda_doubleCV(
        X, y, cv1_splits=2, cv2_splits=3, n_repeats=2, max_components=4, progress=False, scale=False, nuisance=design
    )
    assert all(model.scale is False for model in result["models"])
    coefficients, _ = fit_nuisance_coefficients(X, design)
    X_res = np.asarray(residualize(X.to_numpy(), design.covariates, coefficients))
    for model, lv in zip(result["models"], result["table"]["LV"]):
        assert model.n_components == lv
        expected = fit_plsda_model(X_res, y, n_components=int(lv), scale=False)
        np.testing.assert_allclose(model.x_scores_, expected.x_scores_, atol=1e-10)
    # Smallest training set: an inner training set of an outer rest set (36 → 24 → 12 rows).
    assert result["nuisance_min_training_df"] == 12 - 1 - 2 - 1


def _recording_residualizer(calls: list):
    original = pls.residualize

    def recording(X, covariates, coefficients):
        out = original(X, covariates, coefficients)
        calls.append((np.asarray(X).copy(), np.asarray(coefficients).copy(), np.asarray(out).copy()))
        return out

    return recording


def test_nuisance_cv_never_reads_held_out_labels(monkeypatch):
    """Permuting one outer fold's held-out conditioning rows leaves that fold's fits and matrices unchanged."""

    X, y, design, _, _ = _nuisance_problem()
    cv1_splits, cv2_splits, random_state = 2, 3, 5
    kwargs = dict(cv1_splits=cv1_splits, cv2_splits=cv2_splits, n_repeats=1, max_components=4,
                  random_state=random_state, progress=False, scale=False)

    def run(nuisance: NuisanceDesign) -> list:
        calls: list = []
        monkeypatch.setattr(pls, "residualize", _recording_residualizer(calls))
        plsda_doubleCV(X, y, nuisance=nuisance, **kwargs)
        monkeypatch.undo()
        return calls

    baseline = run(design)
    per_fold = 2 + 2 * cv1_splits  # rest + test, then train + validation per inner fold
    assert len(baseline) == cv2_splits * per_fold + 1  # + the full-input refit
    folds = list(RepeatedStratifiedKFold(n_splits=cv2_splits, n_repeats=1, random_state=random_state).split(X, y))
    for f, (_, test) in enumerate(folds):
        conditioning = design.conditioning.copy()
        conditioning[test] = conditioning[np.roll(test, 1)]
        assert not np.array_equal(conditioning, design.conditioning)
        permuted = NuisanceDesign(conditioning, design.covariates, design.conditioning_names, design.covariate_names)
        calls = run(permuted)
        own = slice(f * per_fold, (f + 1) * per_fold)
        for (x_a, coef_a, out_a), (x_b, coef_b, out_b) in zip(baseline[own], calls[own], strict=True):
            assert np.array_equal(x_a, x_b)
            assert np.array_equal(coef_a, coef_b)
            assert np.array_equal(out_a, out_b)
        # The permuted rows do train the other folds, so the comparison has teeth.
        assert not all(np.array_equal(a[1], b[1]) for a, b in zip(baseline, calls))


def test_fold_coefficients_differ_from_full_sample_coefficients(monkeypatch):
    X, y, design, _, _ = _nuisance_problem()
    calls: list = []
    monkeypatch.setattr(pls, "residualize", _recording_residualizer(calls))
    plsda_doubleCV(X, y, cv1_splits=2, cv2_splits=3, n_repeats=1, max_components=4, progress=False, nuisance=design)
    full, _ = fit_nuisance_coefficients(X, design)
    np.testing.assert_allclose(calls[-1][1], full)
    assert all(not np.allclose(coef, full) for _, coef, _ in calls[:-1])


def test_fold_missing_a_factor_level_names_the_fold():
    X, y, design, _, z = _nuisance_problem()
    rare = np.zeros(z.size)
    rare[5] = 1.0  # one sample carries this level: some training set lacks it
    covariates = np.column_stack([z, rare])
    nuisance = NuisanceDesign(design.conditioning, covariates, design.conditioning_names, ("z", "batch[rare]"))
    fit_nuisance_coefficients(X, nuisance)  # the full sample is of full rank
    with pytest.raises(ValueError, match=r"repeat 1, outer fold \d+(, inner fold \d+)?: .*not of full rank"):
        plsda_doubleCV(X, y, cv1_splits=2, cv2_splits=3, n_repeats=1, max_components=4, progress=False,
                       nuisance=nuisance)


def test_nuisance_rows_must_match_X():
    X, y, design, _, _ = _nuisance_problem()
    with pytest.raises(ValueError, match="nuisance design has 35"):
        plsda_doubleCV(X, y, cv1_splits=2, cv2_splits=3, n_repeats=1, max_components=4, progress=False,
                       nuisance=design.rows(np.arange(35)))


def test_strata_with_the_same_partition_reproduce_the_default_folds():
    X, y, _, _, _ = _nuisance_problem()
    kwargs = dict(cv1_splits=2, cv2_splits=3, n_repeats=2, max_components=4, random_state=7, progress=False)
    default = plsda_doubleCV(X, y, **kwargs)
    relabelled = plsda_doubleCV(X, y, strata=y.map({"Low": "c", "Mid": "a", "High": "b"}).to_numpy(), **kwargs)
    pd.testing.assert_frame_equal(default["table"], relabelled["table"])


def _group_sorted_nuisance():
    """Rows sorted by group within the stage-major generator order: one group × stage cell is contiguous."""

    rng = np.random.default_rng(4)
    groups = np.repeat([0, 1], 18)
    stages = np.tile(np.repeat([0, 1, 2], 6), 2)
    order = np.lexsort((groups, stages))
    groups, stages = groups[order], stages[order]
    g = groups.astype(float)
    s1, s2 = (stages == 1).astype(float), (stages == 2).astype(float)
    conditioning = np.column_stack([g, s1, s2, g * s1, g * s2])
    z = stages + rng.normal(size=stages.size)
    X = pd.DataFrame(rng.normal(size=(stages.size, 6)) + np.outer(stages + 0.5 * z, np.ones(6)))
    names = ("group[1]", "stage[1]", "stage[2]", "group[1]:stage[1]", "group[1]:stage[2]")
    nuisance = NuisanceDesign(conditioning, z[:, None], names, ("z",), "group_stage")
    labels = pd.Series(np.array(["Low", "Mid", "High"])[stages])
    return X, labels, nuisance, groups, stages


def test_conditioning_cells_partition_the_rows():
    _, _, nuisance, groups, stages = _group_sorted_nuisance()
    cells = nuisance.cells()
    assert len(set(cells)) == 6
    for label in set(cells):
        members = cells == label
        assert len(set(zip(groups[members], stages[members]))) == 1


def test_strata_on_conditioning_cells_keep_every_cell_in_training():
    X, y, nuisance, _, _ = _group_sorted_nuisance()
    kwargs = dict(cv1_splits=3, cv2_splits=3, n_repeats=1, max_components=4, progress=False, scale=False)
    with pytest.raises(ValueError, match=r"repeat 1, outer fold \d+, inner fold \d+: .*not of full rank"):
        plsda_doubleCV(X, y, nuisance=nuisance, **kwargs)  # stage-only folds empty a group × stage cell
    result = plsda_doubleCV(X, y, nuisance=nuisance, strata=nuisance.cells(), **kwargs)
    assert result["table"].shape[0] == 1


def test_strata_rows_must_match_X():
    X, y = _synthetic_data()
    with pytest.raises(ValueError, match="strata has 29"):
        plsda_doubleCV(X, y, cv1_splits=3, cv2_splits=3, n_repeats=1, max_components=3, progress=False,
                       strata=y.to_numpy()[:29])

# ── Candidate floor (min_components) ─────────────────────────────────────────


def _floor_data(n_per_class: int = 12, p: int = 9, seed: int = 20261020) -> tuple[pd.DataFrame, pd.Series]:
    rng = np.random.default_rng(seed)
    labels = np.repeat(["Low", "Mid", "High"], n_per_class)
    shift = np.repeat([0.0, 0.8, 1.6], n_per_class)[:, None]
    X = pd.DataFrame(rng.standard_normal((3 * n_per_class, p)) + shift * np.linspace(0, 1, p),
                     columns=[f"f{i}" for i in range(p)])
    return X, pd.Series(labels)


def test_min_components_restricts_the_candidates_fitted(monkeypatch):
    X, y = _floor_data()
    requested: list[int] = []
    original = pls._plsda_auroc

    def recording(n_components, *args, **kwargs):
        requested.append(int(n_components))
        return original(n_components, *args, **kwargs)

    monkeypatch.setattr(pls, "_plsda_auroc", recording)
    result = plsda_doubleCV(X, y, cv1_splits=2, cv2_splits=3, n_repeats=2, max_components=5, random_state=11,
                            progress=False, min_components=2)
    assert set(requested) == {2, 3, 4}
    assert (result["table"]["LV"] >= 2).all()
    assert all(model.n_components >= 2 for model in result["models"])
    assert all(np.asarray(model.x_scores_).shape[1] >= 2 for model in result["models"])


@pytest.mark.parametrize("floor", [5, 0, 2.5, -1])
def test_infeasible_min_components_is_rejected_before_any_fit(monkeypatch, floor):
    X, y = _floor_data()

    def forbidden(*args, **kwargs):
        raise AssertionError("no fold may be fitted when the floor is rejected")

    monkeypatch.setattr(pls, "_plsda_auroc", forbidden)
    with pytest.raises(ValueError, match=rf"min_components={floor!r}.*max_components=5"):
        plsda_doubleCV(X, y, cv1_splits=2, cv2_splits=3, n_repeats=1, max_components=5, progress=False,
                       min_components=floor)


def test_min_components_of_one_is_the_default_rule():
    X, y = _floor_data()
    kwargs = dict(cv1_splits=2, cv2_splits=3, n_repeats=2, max_components=4, random_state=11, progress=False)
    default = plsda_doubleCV(X, y, **kwargs)
    explicit = plsda_doubleCV(X, y, min_components=1, **kwargs)
    pd.testing.assert_frame_equal(default["table"], explicit["table"])
    for a, b in zip(default["models"], explicit["models"]):
        np.testing.assert_array_equal(a.x_scores_, b.x_scores_)


if __name__ == "__main__":
    DEFAULT_IDENTITY.write_text(json.dumps(_default_identity_payload(), indent=2) + "\n", encoding="utf-8")
    print(f"wrote {DEFAULT_IDENTITY}")
