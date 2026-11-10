"""The per-feature interaction scan by block in the real-data analysis, and the scan's input helper."""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from motco import analysis
from motco.analysis import TrajectoryAnalysisError, TrajectoryAnalysisParams, analyze_trajectories
from motco.stats.interaction import interaction_feature_scan
from tests.test_analysis import SMALL_CV, make_dataset
from tests.test_analysis_span_split import six_covariate_dataset

COMMON = {"residual_correction": "df", "span_split_test": True}


def _params(**kwargs) -> TrajectoryAnalysisParams:
    kwargs.setdefault("permutations", 99)
    kwargs.setdefault("seed", 3)
    return TrajectoryAnalysisParams(**SMALL_CV, **kwargs)


def _block_labels(dataset) -> list[str]:
    return [block.name for block in dataset.blocks for _ in range(block.values.shape[1])]


@pytest.fixture(scope="module")
def pair():
    dataset = six_covariate_dataset()
    joint = analyze_trajectories(dataset, _params(**COMMON, interaction_scan=True))
    by_block = analyze_trajectories(dataset, _params(**COMMON, interaction_scan=True, interaction_scan_by_block=True))
    return dataset, joint, by_block


# ── The input helper ─────────────────────────────────────────────────────────


@pytest.mark.parametrize("weighting", [None, "equal_variance", "mfa"])
def test_helper_inputs_give_the_analysis_scan_bit_for_bit(weighting) -> None:
    dataset = six_covariate_dataset()
    params = _params(**COMMON, interaction_scan=True, block_weighting=weighting)
    result = analyze_trajectories(dataset, params)
    inputs = analysis._interaction_scan_inputs(dataset, params)
    record, table = analysis._interaction_scan_outputs(inputs, params)
    assert analysis._json_safe(record) == result.record["interaction_scan"]
    pd.testing.assert_frame_equal(table, result.interaction_scan, check_exact=True)
    scan = interaction_feature_scan(
        inputs.features.to_numpy(dtype=float),
        inputs.model_full,
        inputs.model_reduced,
        permutations=99,
        seed=3,
        residual_correction="df",
    )
    np.testing.assert_array_equal(scan.f, result.interaction_scan["f"].to_numpy())
    np.testing.assert_array_equal(scan.p_value_adjusted, result.interaction_scan["p_value_adjusted"].to_numpy())
    assert list(inputs.features.columns) == list(result.latent_artifacts.features.columns)
    assert inputs.stage_levels == list(dataset.stages)


def test_helper_features_are_unweighted() -> None:
    dataset = six_covariate_dataset()
    plain = analysis._interaction_scan_inputs(dataset, _params(interaction_scan=True))
    weighted = analysis._interaction_scan_inputs(dataset, _params(interaction_scan=True, block_weighting="mfa"))
    pd.testing.assert_frame_equal(weighted.features, plain.features, check_exact=True)
    fitted = analyze_trajectories(dataset, _params(permutations=0, interaction_scan=True))
    pd.testing.assert_frame_equal(plain.features, fitted.latent_artifacts.features, check_exact=True)


def test_helper_fits_no_latent_space(monkeypatch) -> None:
    def forbidden(*args, **kwargs):  # pragma: no cover - must not be reached
        raise AssertionError("the latent space must not be fitted")

    monkeypatch.setattr(analysis, "fit_pls_latent_space", forbidden)
    inputs = analysis._interaction_scan_inputs(six_covariate_dataset(), _params(interaction_scan=True))
    assert inputs.model_full.shape[1] == inputs.model_reduced.shape[1] + 2
    assert inputs.ls_means.shape == (6, inputs.model_full.shape[1])
    assert len(inputs.contrast) == 2


# ── The setting ──────────────────────────────────────────────────────────────


def test_default_is_off_and_omitted_from_resolved_parameters() -> None:
    assert TrajectoryAnalysisParams().interaction_scan_by_block is False
    assert "interaction_scan_by_block" not in TrajectoryAnalysisParams(interaction_scan=True).to_dict()
    params = TrajectoryAnalysisParams(interaction_scan=True, interaction_scan_by_block=True)
    assert params.to_dict()["interaction_scan_by_block"] is True


def test_by_block_without_the_scan_is_rejected_before_fitting(monkeypatch) -> None:
    def forbidden(*args, **kwargs):  # pragma: no cover - must not be reached
        raise AssertionError("the latent space must not be fitted")

    monkeypatch.setattr(analysis, "fit_pls_latent_space", forbidden)
    with pytest.raises(TrajectoryAnalysisError, match="interaction_scan_by_block requires interaction_scan"):
        TrajectoryAnalysisParams(interaction_scan_by_block=True)


def test_by_block_equals_the_primitive_with_the_declared_blocks(pair) -> None:
    dataset, _, by_block = pair
    inputs = analysis._interaction_scan_inputs(dataset, _params(interaction_scan=True))
    expected = interaction_feature_scan(
        inputs.features.to_numpy(dtype=float),
        inputs.model_full,
        inputs.model_reduced,
        permutations=99,
        seed=3,
        residual_correction="df",
        blocks=_block_labels(dataset),
    )
    block = by_block.record["interaction_scan"]
    assert block["combination"] == "block_bonferroni" and block["n_blocks_tested"] == 2
    assert block["p_value_global"] == expected.p_value_global
    np.testing.assert_array_equal(by_block.interaction_scan["p_value_adjusted"].to_numpy(), expected.p_value_adjusted)
    columns = list(inputs.features.columns)
    assert [f"{hit['block']}__{hit['feature']}" for hit in block["hits"]] == [columns[c] for c in expected.hits()]
    for entry in expected.blocks:
        recorded = block["blocks"][entry.name]
        assert recorded["max_f"] == entry.max_f
        assert f"{entry.name}__{recorded['max_feature']}" == columns[entry.max_column]
        assert recorded["p_value"] == entry.p_value
    assert list(block["blocks"]) == ["rna", "atac"]


def test_other_outputs_and_the_f_are_unchanged(pair) -> None:
    _, joint, by_block = pair
    pd.testing.assert_frame_equal(by_block.latent_scores, joint.latent_scores, check_exact=True)
    pd.testing.assert_frame_equal(by_block.stage_means, joint.stage_means, check_exact=True)
    pd.testing.assert_frame_equal(by_block.pls_loadings, joint.pls_loadings, check_exact=True)
    for key in ("statistics", "p_values", "null_summary", "config_spectrum", "latent", "design", "covariates"):
        assert by_block.record[key] == joint.record[key]
    assert by_block.record["span_split"] == joint.record["span_split"]
    assert by_block.record["attribution"] == joint.record["attribution"]
    strip = lambda block: {k: v for k, v in block.items() if k != "runtime_seconds"}  # noqa: E731
    assert strip(by_block.record["measurement"]) == strip(joint.record["measurement"])
    table, reference = by_block.interaction_scan, joint.interaction_scan
    assert list(table.columns) == list(reference.columns)
    unchanged = [column for column in reference.columns if column != "p_value_adjusted"]
    pd.testing.assert_frame_equal(table[unchanged], reference[unchanged], check_exact=True)
    scan, joint_scan = by_block.record["interaction_scan"], joint.record["interaction_scan"]
    for key in ("n_features", "n_degenerate", "df_interaction", "df_residual", "max_f", "max_feature", "max_block"):
        assert scan[key] == joint_scan[key]
    for name, values in joint_scan["blocks"].items():
        assert scan["blocks"][name]["null_max_share"] == values["null_max_share"]
    assert by_block.record["params"]["interaction_scan_by_block"] is True


def test_record_block_order_and_contents(pair) -> None:
    _, _, by_block = pair
    block = by_block.record["interaction_scan"]
    assert list(block) == [
        "n_features",
        "n_degenerate",
        "df_interaction",
        "df_residual",
        "max_f",
        "max_feature",
        "max_block",
        "combination",
        "n_blocks_tested",
        "p_value_global",
        "hit_alpha",
        "hits",
        "blocks",
        "permutations",
        "seed",
        "residual_correction",
        "draws_match_rrpp",
    ]
    entry_keys = ["n_features", "n_degenerate", "max_f", "max_feature", "hits", "null_max_share", "p_value"]
    for values in block["blocks"].values():
        assert list(values) == entry_keys
    smallest = min(values["p_value"] for values in block["blocks"].values())
    assert block["p_value_global"] == min(1.0, 2 * smallest)


def test_the_joint_scan_carries_no_by_block_key(pair) -> None:
    _, joint, _ = pair
    block = joint.record["interaction_scan"]
    assert {"combination", "n_blocks_tested"}.isdisjoint(block)
    for values in block["blocks"].values():
        assert set(values) == {"n_features", "n_degenerate", "hits", "null_max_share"}
    assert "interaction_scan_by_block" not in joint.record["params"]


def test_zero_permutations_record_the_block_maxima_only() -> None:
    result = analyze_trajectories(
        six_covariate_dataset(), _params(permutations=0, interaction_scan=True, interaction_scan_by_block=True)
    )
    block = result.record["interaction_scan"]
    assert block["combination"] == "block_bonferroni" and block["n_blocks_tested"] == 2
    assert {"p_value_global", "hits", "hit_alpha"}.isdisjoint(block)
    for values in block["blocks"].values():
        assert set(values) == {"n_features", "n_degenerate", "max_f", "max_feature"}
        assert values["max_f"] > 0


def test_a_wholly_degenerate_block_is_recorded_untested() -> None:
    from tests.test_analysis import make_real_like

    metadata, blocks = make_real_like()
    blocks["atac"] = blocks["atac"] * 0.0 + 1.5
    result = analyze_trajectories(
        make_dataset(metadata, blocks),
        _params(permutations=49, forced_components=2, interaction_scan=True, interaction_scan_by_block=True),
    )
    block = result.record["interaction_scan"]
    assert block["n_blocks_tested"] == 1
    atac = block["blocks"]["atac"]
    assert atac["n_degenerate"] == atac["n_features"]
    assert atac["max_f"] is None and atac["max_feature"] is None and atac["p_value"] is None
    assert atac["hits"] == 0
    # One tested block: its p-values are the joint scan's.
    joint = analyze_trajectories(
        make_dataset(metadata, blocks), _params(permutations=49, forced_components=2, interaction_scan=True)
    )
    assert block["p_value_global"] == joint.record["interaction_scan"]["p_value_global"]
    np.testing.assert_array_equal(
        result.interaction_scan["p_value_adjusted"].to_numpy(), joint.interaction_scan["p_value_adjusted"].to_numpy()
    )
