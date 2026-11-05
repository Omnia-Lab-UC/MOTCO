"""The split interaction test as an opt-in of the real-data analysis."""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from motco import analysis
from motco.analysis import TrajectoryAnalysisError, TrajectoryAnalysisParams, analyze_trajectories
from motco.simulations.diagnostics import LinearLatentMap
from motco.simulations.evaluation import build_simulation_trajectory_design
from motco.stats.design import build_covariate_matrix
from motco.stats.interaction import interaction_split_test
from tests.test_analysis import SMALL_CV, make_dataset, make_real_like

SIX = ("age", "pmi", "rin", "log10_n_rna", "log10_n_atac", "batch_score")


def six_covariate_dataset(seed: int = 4):
    metadata, blocks = make_real_like()
    rng = np.random.default_rng(seed)
    for name in SIX:
        metadata[name] = rng.normal(size=metadata.shape[0])
    return make_dataset(metadata, blocks, covariate_cols=list(SIX))


def _params(**kwargs) -> TrajectoryAnalysisParams:
    kwargs.setdefault("permutations", 99)
    kwargs.setdefault("seed", 3)
    return TrajectoryAnalysisParams(**SMALL_CV, **kwargs)


@pytest.fixture(scope="module")
def pair():
    dataset = six_covariate_dataset()
    without = analyze_trajectories(dataset, _params(residual_correction="df"))
    with_split = analyze_trajectories(dataset, _params(residual_correction="df", span_split_test=True))
    return dataset, without, with_split


def test_default_is_off_and_omitted_from_resolved_parameters() -> None:
    params = TrajectoryAnalysisParams()
    assert params.span_split_test is False
    assert "span_split_test" not in params.to_dict()
    assert TrajectoryAnalysisParams(span_split_test=True).to_dict()["span_split_test"] is True


def test_rejected_with_latent_adjustment_before_fitting(monkeypatch) -> None:
    def no_fit(*args, **kwargs):  # pragma: no cover - must not be reached
        raise AssertionError("the latent space was fitted")

    monkeypatch.setattr(analysis, "fit_pls_latent_space", no_fit)
    params = _params(covariate_latent_adjustment="stage", span_split_test=True)
    with pytest.raises(TrajectoryAnalysisError, match=r"span_split_test.*covariate_latent_adjustment='stage'"):
        analyze_trajectories(six_covariate_dataset(), params)


def test_other_outputs_are_unchanged(pair) -> None:
    _, without, with_split = pair
    pd.testing.assert_frame_equal(with_split.latent_scores, without.latent_scores)
    pd.testing.assert_frame_equal(with_split.stage_means, without.stage_means)
    pd.testing.assert_frame_equal(with_split.pls_loadings, without.pls_loadings)
    for key in ("statistics", "p_values", "null_summary", "config_spectrum", "latent", "design", "covariates"):
        assert with_split.record[key] == without.record[key]
    measurement = {k: v for k, v in with_split.record["measurement"].items() if k != "runtime_seconds"}
    assert measurement == {k: v for k, v in without.record["measurement"].items() if k != "runtime_seconds"}
    assert "span_split" not in without.record
    assert "span_split_test" not in without.record["params"]
    assert with_split.record["params"]["span_split_test"] is True


def test_block_follows_measurement(pair) -> None:
    _, _, with_split = pair
    keys = list(with_split.record)
    assert keys[keys.index("measurement") + 1] == "span_split"
    block = with_split.record["span_split"]
    assert set(block) == {"parts", "span_share", "permutations", "seed", "residual_correction", "draws_match_rrpp"}
    assert (block["permutations"], block["seed"], block["residual_correction"]) == (99, 3, "df")
    assert block["draws_match_rrpp"] is True


def test_block_equals_the_primitive_on_the_same_inputs(pair) -> None:
    dataset, _, with_split = pair
    artifacts = with_split.latent_artifacts
    linear_map = LinearLatentMap.from_pls(artifacts.model)
    features = (artifacts.features.to_numpy(dtype=float) - artifacts.model._x_mean) / artifacts.model._x_std
    basis, _ = np.linalg.qr(np.asarray(artifacts.model.x_rotations_, dtype=float))
    np.testing.assert_array_equal(basis, linear_map.basis)
    groups = {label: str(index) for index, label in enumerate(dataset.groups)}
    stages = {label: str(index) for index, label in enumerate(dataset.stages)}
    working = pd.DataFrame(
        {
            "group": dataset.metadata["Sex"].astype(str).map(groups).to_numpy(),
            "stage": dataset.metadata["ADNC"].astype(str).map(stages).to_numpy(),
        }
    )
    design = build_simulation_trajectory_design(working, group_col="group", stage_col="stage")
    covariates = build_covariate_matrix(dataset.metadata[list(SIX)], list(SIX))
    full = np.hstack([np.asarray(design.model_full, dtype=float), covariates.matrix])
    reduced = np.hstack([np.asarray(design.model_reduced, dtype=float), covariates.matrix])
    expected = interaction_split_test(
        features, full, reduced, basis, permutations=99, seed=3, residual_correction="df"
    ).to_dict()
    block = with_split.record["span_split"]
    assert block["parts"] == analysis._json_safe(expected["parts"])
    assert block["span_share"] == expected["span_share"]
    assert block["parts"]["total"]["dimension"] == 25
    assert block["parts"]["span"]["dimension"] == with_split.record["latent"]["selected_lv"]


def test_zero_permutations_record_observed_values_only() -> None:
    result = analyze_trajectories(
        six_covariate_dataset(), _params(permutations=0, residual_correction="leverage", span_split_test=True)
    )
    block = result.record["span_split"]
    for values in block["parts"].values():
        assert "p_value_f" not in values and "p_value_ss" not in values
        assert values["f"] is not None and values["ss_interaction"] > 0
    assert block["permutations"] == 0
    assert block["residual_correction"] is None
    assert 0.0 <= block["span_share"] <= 1.0


def test_two_jobs_do_not_match_rrpps_draws() -> None:
    dataset = make_dataset()
    serial = analyze_trajectories(dataset, _params(permutations=19, forced_components=2, span_split_test=True))
    parallel = analyze_trajectories(
        dataset, _params(permutations=19, forced_components=2, n_jobs=2, span_split_test=True)
    )
    assert serial.record["span_split"]["draws_match_rrpp"] is True
    assert parallel.record["span_split"]["draws_match_rrpp"] is False
    # The split test itself always draws the serial sequence at the seed.
    assert parallel.record["span_split"]["parts"] == serial.record["span_split"]["parts"]


def test_high_dimensional_block_dimensions() -> None:
    rng = np.random.default_rng(0)
    metadata, _ = make_real_like()
    n = metadata.shape[0]
    stage_index = metadata["ADNC"].map({"Low": 0, "Intermediate": 1, "High": 2}).to_numpy()
    ids = metadata["donor"].tolist()
    blocks = {
        name: pd.DataFrame(
            np.outer(stage_index, rng.normal(size=width)) + rng.normal(size=(n, width)),
            index=ids,
            columns=[f"{name}{j}" for j in range(width)],
        )
        for name, width in (("rna", 6_000), ("atac", 4_000))
    }
    result = analyze_trajectories(
        make_dataset(metadata, blocks), _params(permutations=999, forced_components=5, span_split_test=True)
    )
    block = result.record["span_split"]
    assert [block["parts"][name]["dimension"] for name in ("total", "span", "complement")] == [10_000, 5, 9_995]
    for values in block["parts"].values():
        assert {"p_value_f", "p_value_ss"} <= set(values)
    assert result.record["params"]["span_split_test"] is True
