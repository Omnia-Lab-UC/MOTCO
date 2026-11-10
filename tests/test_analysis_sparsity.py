"""Sparse stage-PLS as an opt-in of the real-data analysis (``TrajectoryAnalysisParams.sparsity``)."""

from __future__ import annotations

import json

import numpy as np
import pandas as pd
import pytest

import motco.analysis as analysis_module
from motco.analysis import (
    SparsitySettings,
    TrajectoryAnalysisError,
    TrajectoryAnalysisParams,
    analyze_trajectories,
    write_analysis_outputs,
)
from motco.simulations.evaluation import AttributionDiagnosticSettings
from motco.stats.sparse_pls import DEFAULT_KEEP_FRACTIONS, SparsePLSDA
from tests.test_analysis import SMALL_CV, make_dataset

ATTRIBUTION = AttributionDiagnosticSettings(enabled=True, bootstrap_replicates=0, top_k=5)
BLOCK = SparsitySettings("block", (0.2,))


def _params(**kwargs) -> TrajectoryAnalysisParams:
    kwargs.setdefault("permutations", 19)
    kwargs.setdefault("seed", 3)
    return TrajectoryAnalysisParams(**SMALL_CV, **kwargs)


@pytest.fixture(scope="module")
def pair():
    dataset = make_dataset()
    plain = analyze_trajectories(dataset, _params(attribution=ATTRIBUTION))
    sparse = analyze_trajectories(dataset, _params(attribution=ATTRIBUTION, sparsity=BLOCK))
    return dataset, plain, sparse


# --------------------------------------------------------------------------- #
# 4.1 The parameter and the record
# --------------------------------------------------------------------------- #


def test_default_is_absent_and_omitted_from_resolved_parameters() -> None:
    assert TrajectoryAnalysisParams().sparsity is None
    assert "sparsity" not in TrajectoryAnalysisParams().to_dict()
    resolved = TrajectoryAnalysisParams(sparsity=SparsitySettings("joint")).to_dict()
    assert resolved["sparsity"] == {"scope": "joint", "keep_fractions": DEFAULT_KEEP_FRACTIONS}
    assert SparsitySettings("block", (0.1, 1.0)).keep_fractions == (1.0, 0.1)


@pytest.mark.parametrize(
    ("make", "match"),
    [
        (lambda: SparsitySettings("blocks"), r"Unknown sparsity scope 'blocks'"),
        (lambda: SparsitySettings("joint", (1.0, 1.5)), r"keep fraction 1\.5 is outside \(0, 1\]"),
        (lambda: SparsitySettings("joint", (0.5, 0.5)), r"repeat a fraction"),
        (lambda: SparsitySettings("joint", ()), r"at least one fraction"),
        (lambda: TrajectoryAnalysisParams(sparsity="block"), r"must be a SparsitySettings"),
        (lambda: TrajectoryAnalysisParams(forced_components=2, sparsity=SparsitySettings("block")),
         r"forced_components=2 .*ladder \[1\.0, 0\.5, 0\.2, 0\.1, 0\.05\]"),
    ],
)
def test_malformed_setting_is_rejected_before_fitting(monkeypatch, make, match) -> None:
    def no_fit(*args, **kwargs):
        raise AssertionError("no fit may be attempted")

    monkeypatch.setattr(analysis_module, "fit_pls_latent_space", no_fit)
    with pytest.raises(TrajectoryAnalysisError, match=match):
        make()


def test_every_component_loads_on_the_selected_fraction_of_each_block(pair) -> None:
    _, _, sparse = pair
    model = sparse.latent_artifacts.model
    assert isinstance(model, SparsePLSDA)
    lv = sparse.record["latent"]["selected_lv"]
    # rna has 15 features, atac 10: round(0.2 · 15) = 3, round(0.2 · 10) = 2.
    assert model.support_sizes((15, 10)).tolist() == [[3, 2]] * lv


def test_record_contents(pair) -> None:
    _, _, sparse = pair
    record = sparse.record
    keys = list(record)
    assert keys.index("sparsity") == keys.index("latent") + 1
    block = record["sparsity"]
    assert block["scope"] == "block"
    assert block["keep_fractions"] == [0.2]
    assert block["selected_keep_fraction"] == 0.2
    assert block["support"] == [{"rna": 3, "atac": 2}] * record["latent"]["selected_lv"]
    params = record["latent"]["integration_params"]
    assert params["sparsity"] == {"scope": "block", "keep_fractions": [0.2]}
    assert params["selected_keep_fraction"] == 0.2
    assert record["params"]["sparsity"] == {"scope": "block", "keep_fractions": [0.2]}
    json.dumps(record)  # JSON-safe


def test_absent_sparsity_writes_nothing(pair) -> None:
    _, plain, _ = pair
    assert "sparsity" not in plain.record
    assert "sparsity" not in plain.record["latent"]
    assert "sparsity" not in plain.record["latent"]["integration_params"]
    assert "sparsity" not in plain.record["params"]


def test_written_loadings_show_dropped_features_at_zero(pair, tmp_path) -> None:
    _, _, sparse = pair
    write_analysis_outputs(sparse, tmp_path)
    loadings = pd.read_csv(tmp_path / "pls_loadings.csv")
    components = [column for column in loadings.columns if column not in ("block", "feature")]
    support = sparse.latent_artifacts.model.support_
    dropped = ~support.any(axis=1)
    assert dropped.sum() >= 25 - 5 * len(components)
    assert (loadings.loc[dropped, components].to_numpy() == 0.0).all()
    assert (loadings.loc[~dropped, components].abs().to_numpy().sum(axis=1) > 0.0).all()


def test_attribution_runs_on_the_sparse_space(pair) -> None:
    _, _, sparse = pair
    assert sparse.attribution is not None
    assert sparse.record["attribution"]["status"] == "computed"


def test_the_split_test_and_the_scan_run_on_the_sparse_space() -> None:
    dataset = make_dataset()
    result = analyze_trajectories(dataset, _params(span_split_test=True, interaction_scan=True, sparsity=BLOCK))
    split = result.record["span_split"]["parts"]
    assert split["span"]["dimension"] == result.record["latent"]["selected_lv"]
    assert split["total"]["ss_interaction"] == pytest.approx(
        split["span"]["ss_interaction"] + split["complement"]["ss_interaction"], rel=1e-9
    )
    plain = analyze_trajectories(dataset, _params(interaction_scan=True))
    # The scan uses neither the latent basis nor the fit, so it is unchanged.
    pd.testing.assert_frame_equal(result.interaction_scan, plain.interaction_scan)


def test_forced_rank_with_one_fraction_and_latent_adjustment() -> None:
    from tests.test_analysis_span_split import six_covariate_dataset

    forced = analyze_trajectories(
        make_dataset(), _params(permutations=0, forced_components=2, sparsity=SparsitySettings("joint", (0.4,)))
    )
    # A joint budget of round(0.4 · 25) = 10 features per component, split across the blocks by the fit.
    support = forced.record["sparsity"]["support"]
    assert len(support) == 2 and all(sum(row.values()) == 10 for row in support)
    assert forced.record["latent"]["component_selection"] == "forced"
    adjusted = analyze_trajectories(
        six_covariate_dataset(), _params(permutations=0, covariate_latent_adjustment="stage", sparsity=BLOCK)
    )
    assert adjusted.record["latent"]["integration_params"]["latent_adjustment"] == "stage"
    assert adjusted.record["sparsity"]["selected_keep_fraction"] == 0.2
    np.testing.assert_array_equal(np.asarray(adjusted.latent_artifacts.model._x_std), np.ones(25))


def test_full_budget_ladder_is_the_ordinary_analysis(pair) -> None:
    dataset, plain, _ = pair
    full = analyze_trajectories(dataset, _params(attribution=ATTRIBUTION, sparsity=SparsitySettings("joint", (1.0,))))
    pd.testing.assert_frame_equal(full.latent_scores, plain.latent_scores)
    assert full.record["statistics"] == plain.record["statistics"]
    assert full.record["p_values"] == plain.record["p_values"]
