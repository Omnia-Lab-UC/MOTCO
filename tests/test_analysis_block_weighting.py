"""Block weighting as an opt-in of the real-data analysis (``TrajectoryAnalysisParams.block_weighting``)."""

from __future__ import annotations

import math

import numpy as np
import pandas as pd
import pytest

import motco.analysis as analysis_module
from motco.analysis import TrajectoryAnalysisError, TrajectoryAnalysisParams, analyze_trajectories
from motco.simulations.evaluation import AttributionDiagnosticSettings
from tests.test_analysis import SMALL_CV, STAGES, make_dataset, make_real_like

ATTRIBUTION = AttributionDiagnosticSettings(enabled=True, bootstrap_replicates=0, top_k=5)


def _params(**kwargs) -> TrajectoryAnalysisParams:
    kwargs.setdefault("permutations", 19)
    kwargs.setdefault("seed", 3)
    return TrajectoryAnalysisParams(**SMALL_CV, **kwargs)


@pytest.fixture(scope="module")
def pair():
    dataset = make_dataset()
    plain = analyze_trajectories(dataset, _params(attribution=ATTRIBUTION))
    weighted = analyze_trajectories(dataset, _params(attribution=ATTRIBUTION, block_weighting="equal_variance"))
    return dataset, plain, weighted


# --------------------------------------------------------------------------- #
# 3.1 The parameter
# --------------------------------------------------------------------------- #


def test_default_is_absent_and_omitted_from_resolved_parameters() -> None:
    assert TrajectoryAnalysisParams().block_weighting is None
    assert "block_weighting" not in TrajectoryAnalysisParams().to_dict()
    for kind in ("equal_variance", "mfa"):
        assert TrajectoryAnalysisParams(block_weighting=kind).to_dict()["block_weighting"] == kind


def test_unknown_kind_is_rejected_before_fitting(monkeypatch) -> None:
    def no_fit(*args, **kwargs):
        raise AssertionError("no fit may be attempted")

    monkeypatch.setattr(analysis_module, "fit_pls_latent_space", no_fit)
    with pytest.raises(
        TrajectoryAnalysisError, match=r"Unknown block_weighting 'pareto'; allowed: \['equal_variance', 'mfa'\]"
    ):
        analyze_trajectories(make_dataset(), _params(block_weighting="pareto"))


def test_the_space_is_fitted_on_weighted_blocks(pair) -> None:
    _, _, weighted = pair
    artifacts = weighted.latent_artifacts
    features = artifacts.features
    variances = features.to_numpy().var(axis=0)
    blocks = np.array([column.split("__", 1)[0] for column in features.columns])
    for block in ("rna", "atac"):
        assert variances[blocks == block].sum() == pytest.approx(1.0, rel=1e-10)
    np.testing.assert_array_equal(np.asarray(artifacts.model._x_std), np.ones(features.shape[1]))


def test_original_scales_divide_by_the_weights(pair) -> None:
    _, plain, weighted = pair
    weights = {block["name"]: block["weight"] for block in weighted.record["block_weighting"]["blocks"]}
    per_feature = np.array([weights[column.split("__", 1)[0]] for column in weighted.latent_artifacts.features.columns])
    np.testing.assert_allclose(
        weighted.latent_artifacts.original_scales, plain.latent_artifacts.original_scales / per_feature, rtol=1e-15
    )


def _raw_transitions(group: str) -> np.ndarray:
    """Original-unit stage-to-stage changes of one group's cell means, straight from the blocks."""

    metadata, blocks = make_real_like()
    joint = pd.concat([blocks["rna"], blocks["atac"]], axis=1).to_numpy()
    means = [joint[((metadata["Sex"] == group) & (metadata["ADNC"] == stage)).to_numpy()].mean(axis=0)
             for stage in STAGES]
    return np.vstack([later - earlier for earlier, later in zip(means, means[1:])])


@pytest.mark.parametrize("which", ["plain", "weighted"])
def test_attribution_maps_back_to_original_units(pair, which) -> None:
    """Observed transitions times the exposed original scales are the raw cell-mean changes, weighted or not."""

    _, plain, weighted = pair
    result = {"plain": plain, "weighted": weighted}[which]
    scales = result.latent_artifacts.original_scales
    for group_index, group in enumerate(("F", "M")):
        transitions = result.attribution.transitions
        mapped = np.vstack([transition.observed.group_transitions[group_index] * scales for transition in transitions])
        np.testing.assert_allclose(mapped, _raw_transitions(group), rtol=1e-9, atol=1e-12)
    effects = result.attribution.feature_effects.dropna(subset=["effect_original"])
    index = {name: position for position, name in enumerate(result.latent_artifacts.features.columns)}
    expected = effects["effect_standardized"].to_numpy() * scales[[index[name] for name in effects["feature"]]]
    np.testing.assert_allclose(effects["effect_original"].to_numpy(), expected, rtol=1e-12)


def test_angle_shape_and_rank_are_measured_in_the_weighted_metric(pair) -> None:
    _, plain, weighted = pair
    assert weighted.record["latent"]["integration_params"]["block_weighting"] == "equal_variance"
    assert set(weighted.record["statistics"]) == set(plain.record["statistics"])
    assert weighted.record["statistics"] != plain.record["statistics"]


# --------------------------------------------------------------------------- #
# 3.2 The record
# --------------------------------------------------------------------------- #


@pytest.mark.parametrize("kind", ["equal_variance", "mfa"])
def test_record_block(kind) -> None:
    dataset = make_dataset()
    result = analyze_trajectories(dataset, _params(permutations=0, block_weighting=kind))
    record = result.record
    keys = list(record)
    assert keys.index("block_weighting") == keys.index("blocks") + 1
    block = record["block_weighting"]
    assert block["kind"] == kind
    assert [entry["name"] for entry in block["blocks"]] == ["rna", "atac"]
    assert [entry["n_features"] for entry in block["blocks"]] == [15, 10]
    if kind == "equal_variance":
        assert [entry["weight"] for entry in block["blocks"]] == [1 / math.sqrt(15), 1 / math.sqrt(10)]
    else:
        features = result.latent_artifacts.features
        for entry in block["blocks"]:
            columns = [column for column in features.columns if column.startswith(f"{entry['name']}__")]
            values = features[columns].to_numpy()
            leading = np.linalg.svd(values - values.mean(axis=0), compute_uv=False)[0]
            assert leading**2 / values.shape[0] == pytest.approx(1.0, rel=1e-10)
    assert record["latent"]["integration_params"]["block_weighting"] == kind
    assert record["params"]["block_weighting"] == kind


def test_absent_weighting_writes_nothing(pair) -> None:
    _, plain, _ = pair
    assert "block_weighting" not in plain.record
    assert "block_weighting" not in plain.record["latent"]["integration_params"]
    assert "block_weighting" not in plain.record["params"]


def test_the_interaction_scan_reads_the_unweighted_features() -> None:
    dataset = make_dataset()
    plain = analyze_trajectories(dataset, _params(permutations=19, interaction_scan=True, forced_components=2))
    weighted = analyze_trajectories(
        dataset, _params(permutations=19, interaction_scan=True, forced_components=2, block_weighting="mfa")
    )
    pd.testing.assert_frame_equal(weighted.interaction_scan, plain.interaction_scan)
    assert weighted.record["interaction_scan"] == plain.record["interaction_scan"]


def test_latent_adjustment_runs_on_the_weighted_blocks() -> None:
    from tests.test_analysis_span_split import six_covariate_dataset

    result = analyze_trajectories(
        six_covariate_dataset(),
        _params(permutations=0, covariate_latent_adjustment="stage", block_weighting="equal_variance"),
    )
    assert result.record["latent"]["integration_params"]["latent_adjustment"] == "stage"
    assert result.record["latent"]["integration_params"]["block_weighting"] == "equal_variance"
    assert result.record["latent_adjustment"]["feature_scale"] == "pooled_standardized_block_weighted"
    np.testing.assert_array_equal(np.asarray(result.latent_artifacts.model._x_std), np.ones(25))


def test_the_split_test_runs_in_the_weighted_metric() -> None:
    dataset = make_dataset()
    result = analyze_trajectories(dataset, _params(span_split_test=True, block_weighting="mfa"))
    split = result.record["span_split"]["parts"]
    assert split["total"]["dimension"] == 25
    assert split["span"]["dimension"] == result.record["latent"]["selected_lv"]
    assert split["total"]["ss_interaction"] == pytest.approx(
        split["span"]["ss_interaction"] + split["complement"]["ss_interaction"], rel=1e-9
    )
