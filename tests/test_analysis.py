"""Real-data entry point: input contract, measurement record, harness equivalence, artifact."""

from __future__ import annotations

import json
import os

import numpy as np
import pandas as pd
import pytest

from motco.analysis import (
    OmicsBlock,
    TrajectoryAnalysisError,
    TrajectoryAnalysisParams,
    TrajectoryDataset,
    analyze_trajectories,
    write_analysis_outputs,
)
from motco.analysis import _json_safe as json_safe
from motco.simulations.diagnostics import geometry_from_means
from motco.simulations.evaluation import (
    AttributionDiagnosticSettings,
    SimulationEvaluationParams,
    evaluate_semisynthetic_trajectory,
)
from motco.simulations.reference import load_reference
from motco.simulations.semisynthetic import SemiSyntheticTrajectoryParams, generate_semisynthetic_trajectory

STAGES = ["Low", "Intermediate", "High"]
SMALL_CV = {"cv2_splits": 3, "cv1_splits": 2, "n_repeats": 2, "max_components": 5}


def make_real_like(
    *,
    stages: list[str] = STAGES,
    n_per_cell: int = 8,
    seed: int = 0,
    groups: tuple[str, str] = ("F", "M"),
) -> tuple[pd.DataFrame, dict[str, pd.DataFrame]]:
    """Two blocks (``rna``, ``atac``) with a stage signal and a group-specific bend."""

    rng = np.random.default_rng(seed)
    rows = []
    rna, atac = [], []
    for stage_index, stage in enumerate(stages):
        for group_index, group in enumerate(groups):
            for replicate in range(n_per_cell):
                rows.append((f"{group}_{stage}_{replicate}", group, stage))
                shift = float(stage_index) * (1.0 + 0.5 * group_index)
                rna.append(shift + rng.normal(scale=0.6, size=15) + np.linspace(0, group_index * stage_index, 15))
                atac.append(0.5 * shift + rng.normal(scale=0.6, size=10))
    ids = [row[0] for row in rows]
    metadata = pd.DataFrame(
        {"donor": ids, "Sex": [row[1] for row in rows], "ADNC": [row[2] for row in rows], "age": 70}
    )
    blocks = {
        "rna": pd.DataFrame(np.asarray(rna), index=ids, columns=[f"g{j}" for j in range(15)]),
        "atac": pd.DataFrame(np.asarray(atac), index=ids, columns=[f"p{j}" for j in range(10)]),
    }
    return metadata, blocks


def make_dataset(metadata=None, blocks=None, **kwargs) -> TrajectoryDataset:
    if metadata is None or blocks is None:
        metadata, blocks = make_real_like()
    transforms = kwargs.pop("transforms", {})
    return TrajectoryDataset(
        blocks=[OmicsBlock(name, frame, transforms.get(name, "identity")) for name, frame in blocks.items()],
        metadata=metadata,
        sample_col="donor",
        group_col="Sex",
        stage_col="ADNC",
        stage_order=kwargs.pop("stage_order", STAGES),
        **kwargs,
    )


# ── Input contract ───────────────────────────────────────────────────────────


def test_misaligned_block_is_rejected() -> None:
    metadata, blocks = make_real_like()
    blocks["atac"] = blocks["atac"].iloc[::-1]
    with pytest.raises(TrajectoryAnalysisError, match="'atac' rows are not aligned"):
        make_dataset(metadata, blocks)


def test_row_count_mismatch_is_rejected() -> None:
    metadata, blocks = make_real_like()
    blocks["rna"] = blocks["rna"].iloc[:-1]
    with pytest.raises(TrajectoryAnalysisError, match="'rna' has 47 rows"):
        make_dataset(metadata, blocks)


@pytest.mark.parametrize("bad", [np.nan, np.inf])
def test_non_finite_value_is_rejected(bad: float) -> None:
    metadata, blocks = make_real_like()
    blocks["rna"].iloc[3, 2] = bad
    with pytest.raises(TrajectoryAnalysisError, match="'rna' contains missing or infinite"):
        make_dataset(metadata, blocks)


def test_duplicate_block_names_are_rejected() -> None:
    metadata, blocks = make_real_like()
    with pytest.raises(TrajectoryAnalysisError, match="repeated: \\['rna'\\]"):
        TrajectoryDataset(
            blocks=[OmicsBlock("rna", blocks["rna"]), OmicsBlock("rna", blocks["atac"])],
            metadata=metadata,
            sample_col="donor",
            group_col="Sex",
            stage_col="ADNC",
            stage_order=STAGES,
        )


@pytest.mark.parametrize(
    "order",
    [["Low", "High"], ["Low", "Intermediate", "High", "Severe"], ["Low", "Low", "Intermediate", "High"]],
    ids=["omits", "extra", "repeats"],
)
def test_incomplete_stage_order_is_rejected(order: list[str]) -> None:
    with pytest.raises(TrajectoryAnalysisError, match="declared .* observed \\['High', 'Intermediate', 'Low'\\]"):
        make_dataset(stage_order=order)


def test_three_groups_are_rejected() -> None:
    metadata, blocks = make_real_like()
    metadata.loc[0, "Sex"] = "X"
    with pytest.raises(TrajectoryAnalysisError, match="Exactly two groups.*\\['F', 'M', 'X'\\]"):
        make_dataset(metadata, blocks)


def test_missing_cell_is_rejected() -> None:
    metadata, blocks = make_real_like()
    keep = ~((metadata["Sex"] == "M") & (metadata["ADNC"] == "High")).to_numpy()
    metadata = metadata[keep].reset_index(drop=True)
    blocks = {name: frame[keep] for name, frame in blocks.items()}
    with pytest.raises(TrajectoryAnalysisError, match="\\('M', 'High'\\)"):
        make_dataset(metadata, blocks)


def test_unknown_transform_is_rejected() -> None:
    with pytest.raises(TrajectoryAnalysisError, match="'rna' declares unknown transform 'log'"):
        make_dataset(transforms={"rna": "log"})


def test_beta_out_of_range_is_rejected() -> None:
    with pytest.raises(TrajectoryAnalysisError, match="'rna' declares beta_to_mvalue.*\\(0, 1\\)"):
        make_dataset(transforms={"rna": "beta_to_mvalue"})


def test_bad_group_order_and_estimator_are_rejected() -> None:
    with pytest.raises(TrajectoryAnalysisError, match="group_order"):
        make_dataset(group_order=["F", "X"])
    with pytest.raises(TrajectoryAnalysisError, match="orientation_estimator"):
        TrajectoryAnalysisParams(orientation_estimator="procrustes")


# ── Measurement record ───────────────────────────────────────────────────────


@pytest.fixture(scope="module")
def result():
    return analyze_trajectories(
        make_dataset(), TrajectoryAnalysisParams(permutations=19, seed=3, **SMALL_CV)
    )


def test_record_layout(result) -> None:
    record = result.record
    assert set(record) == {
        "schema_version", "statistics", "p_values", "null_summary", "config_spectrum", "latent",
        "design", "blocks", "measurement", "params", "software", "attribution",
    }
    assert set(record["p_values"]) == {"delta", "angle", "shape"}
    assert set(record["null_summary"]["angle"]) >= {"count", "mean", "sd", "q50", "q90", "q95", "q99"}
    assert "permutation_pooled_eigengap" in record["config_spectrum"]
    assert set(record["config_spectrum"]["groups"]) == {"F", "M"}
    assert record["latent"]["component_selection"] == "cv"
    assert record["latent"]["cv_mean_auroc"] is not None
    assert record["latent"]["layer_feature_counts"] == {"rna": 15, "atac": 10}
    assert record["measurement"]["orientation_estimator"] == "net_displacement"
    assert record["measurement"]["selected_lv"] == record["latent"]["selected_lv"]
    assert record["design"]["cell_sizes"] == {g: {s: 8 for s in STAGES} for g in ("F", "M")}
    assert record["blocks"][0] == {"name": "rna", "n_features": 15, "transform": "identity", "units": "as_supplied"}
    assert set(record["software"]) >= {"motco", "numpy", "scipy", "sklearn", "pandas"}
    assert record["attribution"]["status"] == "not_requested"
    json.dumps(record, allow_nan=False)


def test_stage_order_is_declared_not_sorted(result) -> None:
    assert result.record["design"]["stage_levels"] == STAGES
    assert result.stage_means["stage"].tolist() == STAGES * 2
    assert result.stage_means["group"].tolist() == ["F"] * 3 + ["M"] * 3


def test_angle_is_net_displacement_on_declared_order(result) -> None:
    means = result.stage_means.drop(columns=["group", "stage"]).to_numpy()
    expected = geometry_from_means(means, ["F", "M"], STAGES, orientation_estimator="net_displacement")
    assert result.record["statistics"]["angle"] == pytest.approx(expected.angle, abs=1e-9)
    # By hand: High minus Low per group.
    d_f, d_m = means[2] - means[0], means[5] - means[3]
    cos = np.dot(d_f, d_m) / (np.linalg.norm(d_f) * np.linalg.norm(d_m))
    assert result.record["statistics"]["angle"] == pytest.approx(np.degrees(np.arccos(cos)), abs=1e-6)


def test_joint_features_are_block_prefixed_in_declared_order(result) -> None:
    loadings = result.pls_loadings
    assert loadings["block"].tolist() == ["rna"] * 15 + ["atac"] * 10
    assert loadings["feature"].tolist()[:2] == ["g0", "g1"]
    assert loadings.shape == (25, 2 + result.record["latent"]["selected_lv"])


def test_forced_rank_skips_cv() -> None:
    forced = analyze_trajectories(make_dataset(), TrajectoryAnalysisParams(forced_components=2))
    assert forced.record["latent"]["component_selection"] == "forced"
    assert forced.record["latent"]["cv_repeats_completed"] == 0
    assert forced.latent_scores.shape == (48, 2)
    assert forced.record["p_values"] == {}


def test_same_seed_gives_identical_p_values(result) -> None:
    again = analyze_trajectories(make_dataset(), TrajectoryAnalysisParams(permutations=19, seed=3, **SMALL_CV))
    assert again.record["p_values"] == result.record["p_values"]
    assert again.record["null_summary"] == result.record["null_summary"]


def test_shape_unavailable_with_two_stages() -> None:
    metadata, blocks = make_real_like(stages=["Low", "High"])
    out = analyze_trajectories(
        make_dataset(metadata, blocks, stage_order=["Low", "High"]),
        TrajectoryAnalysisParams(permutations=9, seed=1, **SMALL_CV),
    )
    assert out.record["statistics"]["shape"] is None
    assert set(out.record["p_values"]) == {"delta", "angle"}
    assert out.record["measurement"]["shape_available"] is False


def test_pc1_is_selectable_and_moves_only_angle(result) -> None:
    pc1 = analyze_trajectories(
        make_dataset(),
        TrajectoryAnalysisParams(permutations=19, seed=3, orientation_estimator="pc1", **SMALL_CV),
    )
    assert pc1.record["measurement"]["orientation_estimator"] == "pc1"
    assert pc1.record["statistics"]["delta"] == result.record["statistics"]["delta"]
    assert pc1.record["statistics"]["shape"] == result.record["statistics"]["shape"]


# ── Attribution ──────────────────────────────────────────────────────────────


def test_attribution_reuses_fitted_estimator(monkeypatch) -> None:
    import motco.simulations.attribution_diagnostics as diag

    seen = {}
    original = diag.analyze_orientation_attribution

    def spy(X, metadata, model, **kwargs):
        seen["model"] = model
        seen["groups"], seen["stages"] = kwargs["groups"], kwargs["stages"]
        return original(X, metadata, model, **kwargs)

    monkeypatch.setattr(diag, "analyze_orientation_attribution", spy)
    import motco.analysis as analysis

    fitted = {}
    core = analysis.fit_pls_latent_space

    def capture(*args, **kwargs):
        out = core(*args, **kwargs)
        fitted["model"] = out.artifacts.model
        return out

    monkeypatch.setattr(analysis, "fit_pls_latent_space", capture)
    out = analyze_trajectories(
        make_dataset(),
        TrajectoryAnalysisParams(
            attribution=AttributionDiagnosticSettings(enabled=True, bootstrap_replicates=5, bootstrap_seed=2, top_k=5),
            **SMALL_CV,
        ),
    )
    attribution = out.record["attribution"]
    assert attribution["status"] == "computed"
    assert attribution["model"]["model_components"] == out.record["latent"]["selected_lv"]
    assert attribution["model"]["stages"] == STAGES
    assert attribution["truth"]["available"] is False
    assert "real data carry no generator truth" in attribution["truth"]["definition"]
    assert all(row["available"] is False for row in attribution["truth_recovery"])
    assert attribution["units"]["original_basis"] == {"rna": "as_supplied", "atac": "as_supplied"}
    assert seen["model"] is fitted["model"]
    assert list(seen["stages"]) == STAGES
    assert out.attribution is not None


# ── Artifact ─────────────────────────────────────────────────────────────────


def test_write_outputs_round_trip(result, tmp_path) -> None:
    out_dir = tmp_path / "analysis"
    paths = write_analysis_outputs(result, out_dir)
    assert {p.name for p in out_dir.iterdir()} == {
        "result.json", "latent_scores.csv", "pls_loadings.csv", "stage_means.csv"
    }
    loaded = json.loads(paths["result"].read_text())
    assert set(loaded["p_values"]) == {"delta", "angle", "shape"}
    assert loaded["latent"]["selected_lv"] == result.record["latent"]["selected_lv"]
    rank = loaded["latent"]["selected_lv"]
    scores = pd.read_csv(paths["latent_scores"], index_col=0)
    assert scores.shape == (48, rank) and scores.index.name == "donor"
    loadings = pd.read_csv(paths["pls_loadings"])
    assert list(loadings.columns) == ["block", "feature"] + [f"pls_{i}" for i in range(rank)]
    means = pd.read_csv(paths["stage_means"])
    assert means.shape == (6, 2 + rank) and means["stage"].tolist() == STAGES * 2

    mtimes = {p: os.stat(p).st_mtime_ns for p in paths.values()}
    with pytest.raises(FileExistsError):
        write_analysis_outputs(result, out_dir)
    assert {p: os.stat(p).st_mtime_ns for p in paths.values()} == mtimes
    write_analysis_outputs(result, out_dir, overwrite=True)


def test_write_outputs_with_attribution(tmp_path) -> None:
    out = analyze_trajectories(
        make_dataset(),
        TrajectoryAnalysisParams(
            attribution=AttributionDiagnosticSettings(enabled=True, bootstrap_replicates=3, top_k=5), **SMALL_CV
        ),
    )
    paths = write_analysis_outputs(out, tmp_path)
    attribution_files = [key for key in paths if key.startswith("attribution/")]
    assert attribution_files and all(paths[key].exists() for key in attribution_files)


# ── Equivalence with the simulation harness ──────────────────────────────────


@pytest.fixture(scope="module")
def synthetic():
    params = SemiSyntheticTrajectoryParams(
        seed=4,
        trajectory_mode="orientation",
        group_effect_size=0.3,
        p_dmp=0.1,
        group_stage_sizes=((11, 10, 28), (9, 10, 12)),
    )
    return generate_semisynthetic_trajectory(params, reference=load_reference())


def _as_trajectory_dataset(dataset, layers=("methylation", "expression", "proteomics")) -> TrajectoryDataset:
    return TrajectoryDataset(
        blocks=[
            OmicsBlock(layer, getattr(dataset, layer), "beta_to_mvalue" if layer == "methylation" else "identity")
            for layer in layers
        ],
        metadata=dataset.metadata,
        sample_col="sample_id",
        group_col="group",
        stage_col="stage",
        stage_order=["0", "1", "2"],
    )


def _assert_equivalent(harness_result, analysis_result) -> None:
    record = analysis_result.record
    assert record["statistics"] == json_safe(harness_result.pair_statistics)
    assert record["p_values"] == json_safe(harness_result.p_values)
    assert record["null_summary"] == json_safe(harness_result.null_summary)
    assert record["config_spectrum"] == json_safe(harness_result.config_spectrum)
    assert record["latent"]["selected_lv"] == harness_result.latent_matrix_metadata["selected_lv"]
    assert record["latent"]["feature_order_signature"] == harness_result.latent_matrix_metadata[
        "feature_order_signature"
    ]


def _run_pair(dataset, *, estimator=None, layers=None):
    integration_params = dict(SMALL_CV)
    if layers is not None:
        integration_params["layers"] = list(layers)
    harness_kwargs = {} if estimator is None else {"orientation_estimator": estimator}
    harness = evaluate_semisynthetic_trajectory(
        dataset,
        SimulationEvaluationParams(
            integration_method="pls",
            integration_params=integration_params,
            permutations=49,
            seed=11,
            **harness_kwargs,
        ),
    )
    analysis_kwargs = {} if estimator is None else {"orientation_estimator": estimator}
    analysis = analyze_trajectories(
        _as_trajectory_dataset(dataset, layers or ("methylation", "expression", "proteomics")),
        TrajectoryAnalysisParams(permutations=49, seed=11, **SMALL_CV, **analysis_kwargs),
    )
    return harness, analysis


def test_matches_harness_default_estimator(synthetic) -> None:
    from motco.simulations.evaluation import integrate_semisynthetic_dataset

    harness, analysis = _run_pair(synthetic)
    _assert_equivalent(harness, analysis)
    latent = integrate_semisynthetic_dataset(
        synthetic, SimulationEvaluationParams(integration_method="pls", integration_params=SMALL_CV)
    )
    assert np.array_equal(latent.matrix.to_numpy(), analysis.latent_scores.to_numpy())
    assert analysis.record["measurement"]["orientation_estimator"] == "net_displacement"


def test_matches_harness_pc1(synthetic) -> None:
    harness, analysis = _run_pair(synthetic, estimator="pc1")
    _assert_equivalent(harness, analysis)


def test_matches_harness_layer_selection(synthetic) -> None:
    from motco.simulations.evaluation import integrate_semisynthetic_dataset

    layers = ("methylation", "expression")
    harness, analysis = _run_pair(synthetic, layers=layers)
    _assert_equivalent(harness, analysis)
    latent = integrate_semisynthetic_dataset(
        synthetic,
        SimulationEvaluationParams(integration_method="pls", integration_params={**SMALL_CV, "layers": list(layers)}),
    )
    assert np.array_equal(latent.matrix.to_numpy(), analysis.latent_scores.to_numpy())
