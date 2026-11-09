"""The per-feature interaction scan as an opt-in of the real-data analysis."""

from __future__ import annotations

import json

import numpy as np
import pandas as pd
import pytest

from motco import analysis
from motco.analysis import TrajectoryAnalysisParams, analyze_trajectories, write_analysis_outputs
from motco.simulations.evaluation import AttributionDiagnosticSettings, build_simulation_trajectory_design
from motco.stats.design import build_covariate_matrix
from motco.stats.interaction import interaction_feature_scan
from tests.test_analysis import SMALL_CV, make_dataset, make_real_like
from tests.test_analysis_span_split import SIX, six_covariate_dataset

STAGES = ["Low", "Intermediate", "High"]


def _params(**kwargs) -> TrajectoryAnalysisParams:
    kwargs.setdefault("permutations", 99)
    kwargs.setdefault("seed", 3)
    return TrajectoryAnalysisParams(**SMALL_CV, **kwargs)


def _models(dataset):
    """The analysis's covariate-adjusted full and reduced models, rebuilt independently."""

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
    return full, reduced, covariates


@pytest.fixture(scope="module")
def pair():
    dataset = six_covariate_dataset()
    common = {"residual_correction": "df", "span_split_test": True}
    without = analyze_trajectories(dataset, _params(**common))
    with_scan = analyze_trajectories(dataset, _params(**common, interaction_scan=True))
    return dataset, without, with_scan


def test_default_is_off_and_omitted_from_resolved_parameters() -> None:
    params = TrajectoryAnalysisParams()
    assert params.interaction_scan is False
    assert "interaction_scan" not in params.to_dict()
    assert TrajectoryAnalysisParams(interaction_scan=True).to_dict()["interaction_scan"] is True


def test_other_outputs_are_unchanged(pair) -> None:
    _, without, with_scan = pair
    pd.testing.assert_frame_equal(with_scan.latent_scores, without.latent_scores)
    pd.testing.assert_frame_equal(with_scan.stage_means, without.stage_means)
    pd.testing.assert_frame_equal(with_scan.pls_loadings, without.pls_loadings)
    for key in ("statistics", "p_values", "null_summary", "config_spectrum", "latent", "design", "covariates"):
        assert with_scan.record[key] == without.record[key]
    assert with_scan.record["span_split"] == without.record["span_split"]
    assert with_scan.record["attribution"] == without.record["attribution"]
    measurement = {k: v for k, v in with_scan.record["measurement"].items() if k != "runtime_seconds"}
    assert measurement == {k: v for k, v in without.record["measurement"].items() if k != "runtime_seconds"}
    assert "interaction_scan" not in without.record
    assert without.interaction_scan is None
    assert "interaction_scan" not in without.record["params"]
    assert with_scan.record["params"]["interaction_scan"] is True


def test_attribution_outputs_are_unchanged() -> None:
    dataset = six_covariate_dataset()
    settings = AttributionDiagnosticSettings(enabled=True, bootstrap_replicates=4, bootstrap_seed=1, top_k=5)
    common = {"permutations": 19, "attribution": settings}
    without = analyze_trajectories(dataset, _params(**common))
    with_scan = analyze_trajectories(dataset, _params(**common, interaction_scan=True))
    strip = lambda record: {k: v for k, v in record.items() if k != "runtime"}  # noqa: E731
    assert strip(with_scan.record["attribution"]) == strip(without.record["attribution"])
    assert with_scan.record["attribution"]["status"] == "computed"
    assert with_scan.record["p_values"] == without.record["p_values"]
    assert with_scan.interaction_scan is not None


def test_block_follows_measurement_and_span_split(pair) -> None:
    _, _, with_scan = pair
    keys = list(with_scan.record)
    assert keys[keys.index("measurement") + 1 : keys.index("measurement") + 3] == ["span_split", "interaction_scan"]
    block = with_scan.record["interaction_scan"]
    assert list(block) == [
        "n_features",
        "n_degenerate",
        "df_interaction",
        "df_residual",
        "max_f",
        "max_feature",
        "max_block",
        "p_value_global",
        "hit_alpha",
        "hits",
        "blocks",
        "permutations",
        "seed",
        "residual_correction",
        "draws_match_rrpp",
    ]
    assert (block["n_features"], block["n_degenerate"]) == (25, 0)
    assert (block["df_interaction"], block["df_residual"]) == (2, 48 - 12)
    assert (block["permutations"], block["seed"], block["residual_correction"]) == (99, 3, "df")
    assert block["draws_match_rrpp"] is True
    assert block["hit_alpha"] == 0.05
    shares = [values["null_max_share"] for values in block["blocks"].values()]
    assert sum(shares) == pytest.approx(1.0)
    assert {name: values["n_features"] for name, values in block["blocks"].items()} == {"rna": 15, "atac": 10}
    assert sum(values["hits"] for values in block["blocks"].values()) == len(block["hits"])


def test_block_and_table_equal_the_primitive_on_the_same_inputs(pair) -> None:
    dataset, _, with_scan = pair
    full, reduced, _ = _models(dataset)
    X = with_scan.latent_artifacts.features.to_numpy(dtype=float)
    expected = interaction_feature_scan(X, full, reduced, permutations=99, seed=3, residual_correction="df")
    table = with_scan.interaction_scan
    np.testing.assert_allclose(table["f"].to_numpy(), expected.f, rtol=1e-10)
    np.testing.assert_allclose(table["ss_interaction"].to_numpy(), expected.ss_interaction, rtol=1e-10)
    np.testing.assert_array_equal(table["p_value_adjusted"].to_numpy(), expected.p_value_adjusted)
    np.testing.assert_array_equal(table["degenerate"].to_numpy(), expected.degenerate)
    block = with_scan.record["interaction_scan"]
    assert block["p_value_global"] == expected.p_value_global
    assert block["max_f"] == pytest.approx(expected.max_f, rel=1e-10)
    columns = list(with_scan.latent_artifacts.features.columns)
    assert f"{block['max_block']}__{block['max_feature']}" == columns[expected.max_column]
    assert [f"{hit['block']}__{hit['feature']}" for hit in block["hits"]] == [columns[c] for c in expected.hits()]
    assert len(block["hits"]) > 0
    hit_f = [hit["f"] for hit in block["hits"]]
    assert hit_f == sorted(hit_f, reverse=True)


def test_table_rows_follow_the_joint_feature_order(pair) -> None:
    _, _, with_scan = pair
    table = with_scan.interaction_scan
    assert list(table.columns) == [
        "block",
        "feature",
        "ss_interaction",
        "rss",
        "f",
        "p_value_adjusted",
        "degenerate",
        *(f"difference[{stage}]" for stage in STAGES),
    ]
    joint = [f"{block}__{feature}" for block, feature in zip(table["block"], table["feature"])]
    assert joint == list(with_scan.latent_artifacts.features.columns)
    assert list(table["feature"][:3]) == ["g0", "g1", "g2"] and table["block"].iloc[-1] == "atac"


def test_differences_are_the_full_models_ls_mean_differences(pair) -> None:
    dataset, _, with_scan = pair
    full, _, covariates = _models(dataset)
    X = with_scan.latent_artifacts.features.to_numpy(dtype=float)
    coefficients, *_ = np.linalg.lstsq(full, X, rcond=None)
    reference = covariates.matrix.mean(axis=0)
    q = covariates.matrix.shape[1]
    sex = dataset.metadata["Sex"].astype(str).to_numpy()
    stage = dataset.metadata["ADNC"].astype(str).to_numpy()
    first, second = dataset.groups
    for name in STAGES:
        rows = []
        for group in (first, second):
            row = full[np.flatnonzero((sex == group) & (stage == name))[0]].copy()
            row[-q:] = reference  # the cell's design columns at the pooled covariate means
            rows.append(row)
        expected = (rows[1] - rows[0]) @ coefficients
        np.testing.assert_allclose(with_scan.interaction_scan[f"difference[{name}]"].to_numpy(), expected, atol=1e-10)


def test_latent_adjustment_does_not_change_the_scan() -> None:
    dataset = six_covariate_dataset()
    plain = analyze_trajectories(dataset, _params(interaction_scan=True))
    adjusted = analyze_trajectories(dataset, _params(interaction_scan=True, covariate_latent_adjustment="stage"))
    assert "latent_adjustment" in adjusted.record
    assert adjusted.record["interaction_scan"] == plain.record["interaction_scan"]
    pd.testing.assert_frame_equal(adjusted.interaction_scan, plain.interaction_scan)


def test_zero_permutations_record_observed_values_only() -> None:
    result = analyze_trajectories(
        six_covariate_dataset(), _params(permutations=0, residual_correction="leverage", interaction_scan=True)
    )
    block = result.record["interaction_scan"]
    assert {"p_value_global", "hits", "hit_alpha", "draws_match_rrpp"}.isdisjoint(block)
    assert block["n_features"] == 25 and block["max_f"] > 0 and block["max_feature"]
    assert (block["permutations"], block["residual_correction"]) == (0, None)
    for values in block["blocks"].values():
        assert set(values) == {"n_features", "n_degenerate"}
    assert "p_value_adjusted" not in result.interaction_scan.columns


def test_two_jobs_do_not_match_rrpps_draws() -> None:
    dataset = make_dataset()
    serial = analyze_trajectories(dataset, _params(permutations=19, forced_components=2, interaction_scan=True))
    parallel = analyze_trajectories(
        dataset, _params(permutations=19, forced_components=2, n_jobs=2, interaction_scan=True)
    )
    assert serial.record["interaction_scan"]["draws_match_rrpp"] is True
    assert parallel.record["interaction_scan"]["draws_match_rrpp"] is False
    # The scan itself always draws the serial sequence at the seed.
    strip = lambda block: {k: v for k, v in block.items() if k != "draws_match_rrpp"}  # noqa: E731
    assert strip(parallel.record["interaction_scan"]) == strip(serial.record["interaction_scan"])


def test_constant_feature_is_flagged_and_ignored() -> None:
    metadata, blocks = make_real_like()
    blocks["rna"]["g3"] = 2.5
    result = analyze_trajectories(make_dataset(metadata, blocks), _params(permutations=49, interaction_scan=True))
    table = result.interaction_scan
    row = table.loc[(table["block"] == "rna") & (table["feature"] == "g3")].iloc[0]
    assert bool(row["degenerate"]) and np.isnan(row["f"]) and np.isnan(row["p_value_adjusted"])
    block = result.record["interaction_scan"]
    assert block["n_degenerate"] == 1 and block["blocks"]["rna"]["n_degenerate"] == 1
    assert "g3" not in {hit["feature"] for hit in block["hits"] if hit["block"] == "rna"}
    assert block["max_feature"] != "g3" or block["max_block"] != "rna"


def test_scan_errors_are_analysis_errors(monkeypatch) -> None:
    def failing(*args, **kwargs):
        raise ValueError("No feature has a defined F")

    monkeypatch.setattr(analysis, "interaction_feature_scan", failing)
    with pytest.raises(analysis.TrajectoryAnalysisError, match="Interaction scan: No feature has a defined F"):
        analyze_trajectories(make_dataset(), _params(permutations=0, forced_components=2, interaction_scan=True))


def test_write_outputs_add_the_scan_table_only_when_it_ran(pair, tmp_path) -> None:
    _, without, with_scan = pair
    paths = write_analysis_outputs(with_scan, tmp_path / "scan")
    written = pd.read_csv(paths["interaction_scan"])
    assert paths["interaction_scan"].name == "interaction_scan.csv"
    assert list(written.columns) == list(with_scan.interaction_scan.columns)
    assert written.shape[0] == 25
    assert list(written["block"] + "__" + written["feature"]) == list(with_scan.latent_artifacts.features.columns)
    np.testing.assert_allclose(written["f"].to_numpy(), with_scan.interaction_scan["f"].to_numpy(), rtol=1e-12)
    record = json.loads((tmp_path / "scan" / "result.json").read_text())
    assert record["interaction_scan"]["hits"] == with_scan.record["interaction_scan"]["hits"]

    plain = write_analysis_outputs(without, tmp_path / "plain")
    assert "interaction_scan" not in plain
    assert not (tmp_path / "plain" / "interaction_scan.csv").exists()
    assert sorted(path.name for path in (tmp_path / "plain").iterdir()) == [
        "latent_scores.csv",
        "pls_loadings.csv",
        "result.json",
        "stage_means.csv",
    ]


def test_high_dimensional_blocks() -> None:
    rng = np.random.default_rng(0)
    metadata, _ = make_real_like()
    n = metadata.shape[0]
    stage_index = metadata["ADNC"].map({"Low": 0, "Intermediate": 1, "High": 2}).to_numpy()
    ids = metadata["donor"].tolist()
    blocks = {
        name: pd.DataFrame(
            np.outer(stage_index, rng.normal(size=5_000)) + rng.normal(size=(n, 5_000)),
            index=ids,
            columns=[f"{name}{j}" for j in range(5_000)],
        )
        for name in ("rna", "atac")
    }
    result = analyze_trajectories(
        make_dataset(metadata, blocks), _params(permutations=999, forced_components=5, interaction_scan=True)
    )
    block = result.record["interaction_scan"]
    assert block["n_features"] == 10_000 and block["n_degenerate"] == 0
    assert block["max_block"] in {"rna", "atac"}
    assert result.interaction_scan.shape[0] == 10_000
    assert sum(values["null_max_share"] for values in block["blocks"].values()) == pytest.approx(1.0)
    assert result.record["params"]["interaction_scan"] is True
