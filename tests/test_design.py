from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from motco.stats.design import (
    build_conditioning_matrix,
    build_covariate_matrix,
    build_ls_means,
    center_matrix,
    get_model_matrix,
)


def _feature_columns(df: pd.DataFrame, group_col: str, level_col: str) -> list[str]:
    return [
        c for c in df.select_dtypes(include=[np.number]).columns if c not in {group_col, level_col}
    ]


# ── Smoke tests ───────────────────────────────────────────────────────────────

def test_center_matrix_group_means_zero(example_df: pd.DataFrame, group_col: str, level_col: str):
    feat_cols = _feature_columns(example_df, group_col, level_col)
    centered = center_matrix(example_df, group_col=group_col, level_col=level_col, feature_cols=feat_cols)

    grp_means = centered.groupby(group_col)[feat_cols].mean()
    assert np.allclose(grp_means.values, 0.0, atol=1e-10)


def test_model_matrix_and_ls_means_shapes(example_df: pd.DataFrame, group_col: str, level_col: str):
    X = example_df[[group_col, level_col]].copy()
    g_levels = sorted(pd.unique(X[group_col].astype(str)).tolist())
    l_levels = sorted(pd.unique(X[level_col].astype(str)).tolist())

    M = get_model_matrix(X, group_col=group_col, level_col=level_col, full=True)
    Gm1 = max(len(g_levels) - 1, 0)
    Lm1 = max(len(l_levels) - 1, 0)
    expected_cols = 1 + Gm1 + Lm1 + (Gm1 * Lm1)
    assert M.shape[0] == len(X)
    assert M.shape[1] == expected_cols

    LS = build_ls_means(g_levels, l_levels, full=True)
    assert LS.shape[0] == len(g_levels) * len(l_levels)
    assert LS.shape[1] == expected_cols


# ── Validation tests ──────────────────────────────────────────────────────────

def test_get_model_matrix_missing_group_col():
    X = pd.DataFrame({"group": ["A", "B", "A", "B"], "level": ["t0", "t1", "t0", "t1"]})
    with pytest.raises(ValueError, match="'missing'"):
        get_model_matrix(X, group_col="missing", level_col="level")


def test_get_model_matrix_missing_level_col():
    X = pd.DataFrame({"group": ["A", "B", "A", "B"], "level": ["t0", "t1", "t0", "t1"]})
    with pytest.raises(ValueError, match="'missing'"):
        get_model_matrix(X, group_col="group", level_col="missing")


def test_get_model_matrix_single_group():
    X = pd.DataFrame({"group": ["A", "A", "A"], "level": ["t0", "t1", "t0"]})
    with pytest.raises(ValueError, match="unique"):
        get_model_matrix(X, group_col="group", level_col="level")


def test_get_model_matrix_single_level():
    X = pd.DataFrame({"group": ["A", "B", "A"], "level": ["t0", "t0", "t0"]})
    with pytest.raises(ValueError, match="unique"):
        get_model_matrix(X, group_col="group", level_col="level")


def test_center_matrix_missing_feature_col():
    dat = pd.DataFrame({
        "group": ["A", "A", "B", "B"],
        "level": ["t0", "t1", "t0", "t1"],
        "f1": [1.0, 2.0, 3.0, 4.0],
    })
    with pytest.raises(ValueError, match="'ghost_col'"):
        center_matrix(dat, group_col="group", level_col="level", feature_cols=["ghost_col"])


# ── Level ordering (audit D8: "10" must not sort before "2") ─────────────────

def test_sort_levels_numeric_order_for_integer_labels():
    from motco.stats.design import _sort_levels

    labels = [str(i) for i in range(12)]
    assert _sort_levels(list(reversed(labels))) == labels


def test_sort_levels_lexicographic_for_non_numeric_labels():
    from motco.stats.design import _sort_levels

    assert _sort_levels(["B", "A"]) == ["A", "B"]
    # Mixed alphanumeric labels stay lexicographic — no general natural sort.
    assert _sort_levels(["stage10", "stage2"]) == ["stage10", "stage2"]


def test_sort_levels_lexicographic_on_numeric_collision():
    from motco.stats.design import _sort_levels

    # "01" and "1" collide numerically; a numeric key would be ill-defined.
    assert _sort_levels(["1", "01"]) == ["01", "1"]


def test_model_matrix_orders_ten_plus_integer_levels_numerically():
    n_stages = 12
    frame = pd.DataFrame(
        {
            "group": ["A", "B"] * n_stages,
            "stage": [str(s) for s in range(n_stages) for _ in (0, 1)],
        }
    )
    model = get_model_matrix(frame, group_col="group", level_col="stage", full=True)
    # Intercept + (G-1) + (L-1) + (G-1)(L-1) columns.
    assert model.shape == (2 * n_stages, 1 + 1 + (n_stages - 1) + (n_stages - 1))
    # Level dummy columns are 2..2+(L-1): drop-first over the ordered levels,
    # so the row for stage "2" must activate the dummy at ordinal position 2
    # (numeric order), not position 4 (lexicographic "0","1","10","11","2").
    row_stage2 = frame.index[frame["stage"] == "2"][0]
    level_dummies = model[row_stage2, 2 : 2 + n_stages - 1]
    assert level_dummies[1] == 1.0  # dummy for level ordinal 2 (first was dropped)
    assert level_dummies.sum() == 1.0


def test_simulation_design_orders_ten_plus_stages_numerically():
    from motco.simulations.evaluation import build_simulation_trajectory_design

    n_stages = 12
    metadata = pd.DataFrame(
        {
            "group": ["A", "B"] * n_stages,
            "stage": [str(s) for s in range(n_stages) for _ in (0, 1)],
        }
    )
    design = build_simulation_trajectory_design(metadata)
    assert design.stage_levels == [str(s) for s in range(n_stages)]
    assert design.group_levels == ["A", "B"]
    # Contrast indices are group-major, level-minor over the numeric order.
    assert design.contrast == [
        list(range(n_stages)),
        list(range(n_stages, 2 * n_stages)),
    ]
    assert design.ls_means.shape[0] == 2 * n_stages


def test_observed_vectors_index_orders_integer_levels_numerically():
    from motco.stats.trajectory import get_observed_vectors

    n_stages = 11
    rng = np.random.default_rng(1203)
    frame = pd.DataFrame(
        {
            "group": ["A", "B"] * n_stages,
            "stage": [str(s) for s in range(n_stages) for _ in (0, 1)],
        }
    )
    Y = pd.DataFrame(rng.normal(size=(2 * n_stages, 3)), columns=["f1", "f2", "f3"])
    means = get_observed_vectors(frame, Y, "group", "stage", full=True)
    stage_order = means.index.get_level_values("stage")[:n_stages].tolist()
    assert stage_order == [str(s) for s in range(n_stages)]


# ── Covariate design columns ─────────────────────────────────────────────────


def _covariate_frame() -> pd.DataFrame:
    return pd.DataFrame(
        {
            "age": [70.0, 80.0, 90.0, 60.0, 75.0, 85.0],
            "education_level": ["HighSchool", "College", "Graduate", "College", "HighSchool", "College"],
            "batch": [3, 1, 2, 2, 1, 3],
        }
    )


def test_covariate_matrix_numeric_and_categorical() -> None:
    design = build_covariate_matrix(_covariate_frame(), ["age", "education_level"])
    assert design.column_names == ("age", "education_level[Graduate]", "education_level[HighSchool]")
    assert design.kinds == {"age": "numeric", "education_level": "categorical"}
    assert design.levels == {"education_level": ["College", "Graduate", "HighSchool"]}
    assert design.reference_levels == {"education_level": "College"}
    np.testing.assert_array_equal(design.matrix[:, 0], _covariate_frame()["age"].to_numpy())
    np.testing.assert_array_equal(design.matrix[:, 1], [0, 0, 1, 0, 0, 0])
    np.testing.assert_array_equal(design.matrix[:, 2], [1, 0, 0, 0, 1, 0])
    np.testing.assert_allclose(design.reference_row, [76.666666666666667, 1 / 6, 2 / 6])


def test_covariate_matrix_follows_declared_order() -> None:
    design = build_covariate_matrix(_covariate_frame(), ["education_level", "age"])
    assert design.column_names == ("education_level[Graduate]", "education_level[HighSchool]", "age")
    assert design.columns == ("education_level", "age")


def test_numeric_coded_factor_declared_categorical() -> None:
    design = build_covariate_matrix(_covariate_frame(), ["batch"], categorical=["batch"])
    assert design.kinds == {"batch": "categorical"}
    assert design.levels == {"batch": ["1", "2", "3"]}
    assert design.column_names == ("batch[2]", "batch[3]")
    assert design.matrix.shape == (6, 2)


def test_boolean_covariate_is_categorical() -> None:
    frame = pd.DataFrame({"apoe4": [True, False, True, False]})
    design = build_covariate_matrix(frame, ["apoe4"])
    assert design.kinds == {"apoe4": "categorical"}
    assert design.column_names == ("apoe4[True]",)


@pytest.mark.parametrize(
    ("frame", "categorical", "message"),
    [
        (pd.DataFrame({"pmi": [5.0, 5.0, 5.0]}), [], "'pmi' is constant"),
        (pd.DataFrame({"pmi": ["a", "a", "a"]}), [], "'pmi' has 1 observed level"),
        (pd.DataFrame({"pmi": [1, 1, 1]}), ["pmi"], "'pmi' has 1 observed level"),
        (pd.DataFrame({"pmi": [1.0, np.inf, 3.0]}), [], "'pmi' contains non-finite"),
        (pd.DataFrame({"pmi": [1.0, np.nan, 3.0]}), [], "'pmi' contains missing"),
    ],
    ids=["constant", "one-level", "one-level-coded", "infinite", "missing"],
)
def test_covariate_matrix_rejections(frame: pd.DataFrame, categorical: list[str], message: str) -> None:
    with pytest.raises(ValueError, match=message):
        build_covariate_matrix(frame, ["pmi"], categorical=categorical)


# ── Conditioning design for covariate residualization ───────────────────────


def _conditioning_frame() -> pd.DataFrame:
    groups = np.repeat(["0", "1"], 6)
    stages = np.tile(["0", "1", "2"], 4)
    return pd.DataFrame({"group": groups, "stage": stages})


def test_stage_conditioning_is_drop_first_stage_indicators() -> None:
    frame = _conditioning_frame()
    design = build_conditioning_matrix(frame, "group", "stage", "stage")
    assert design.kind == "stage"
    assert design.column_names == ("stage[1]", "stage[2]")
    expected = np.column_stack([(frame["stage"] == "1").astype(float), (frame["stage"] == "2").astype(float)])
    np.testing.assert_array_equal(design.matrix, expected)
    reduced = get_model_matrix(frame, group_col="group", level_col="stage", full=False)
    np.testing.assert_array_equal(design.matrix, reduced[:, 2:])


def test_group_stage_conditioning_matches_the_full_model_without_intercept() -> None:
    frame = _conditioning_frame()
    design = build_conditioning_matrix(frame, "group", "stage", "group_stage")
    assert design.column_names == (
        "group[1]",
        "stage[1]",
        "stage[2]",
        "group[1]:stage[1]",
        "group[1]:stage[2]",
    )
    full = get_model_matrix(frame, group_col="group", level_col="stage", full=True)
    np.testing.assert_array_equal(design.matrix, full[:, 1:])


def test_conditioning_orders_levels_numerically() -> None:
    frame = pd.DataFrame({"group": ["a", "b"] * 6, "stage": [str(level) for level in (10, 2, 1) * 4]})
    design = build_conditioning_matrix(frame, "group", "stage", "stage")
    assert design.column_names == ("stage[2]", "stage[10]")


def test_conditioning_rejects_unknown_kind_and_single_level() -> None:
    frame = _conditioning_frame()
    with pytest.raises(ValueError, match="Unknown conditioning 'group'"):
        build_conditioning_matrix(frame, "group", "stage", "group")
    with pytest.raises(ValueError, match="'stage' has 1 observed level"):
        build_conditioning_matrix(frame.assign(stage="0"), "group", "stage", "stage")
