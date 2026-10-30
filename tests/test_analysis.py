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


# ── Rank floor ───────────────────────────────────────────────────────────────


@pytest.mark.parametrize("floor", [0, 2.5, -1, True])
def test_invalid_rank_floor_is_rejected_before_fitting(monkeypatch, floor) -> None:
    from motco import analysis

    def forbidden(*args, **kwargs):
        raise AssertionError("the latent space must not be fitted when the floor is invalid")

    monkeypatch.setattr(analysis, "fit_pls_latent_space", forbidden)
    with pytest.raises(TrajectoryAnalysisError, match="min_components must be an integer of at least 1"):
        analyze_trajectories(make_dataset(), TrajectoryAnalysisParams(min_components=floor, **SMALL_CV))


def test_rank_floor_spelling() -> None:
    assert TrajectoryAnalysisParams().min_components == 2
    assert TrajectoryAnalysisParams().resolved_min_components == 2
    one = TrajectoryAnalysisParams(min_components=1)
    assert one.min_components is None and one.resolved_min_components == 1
    assert "min_components" not in one.to_dict()
    assert TrajectoryAnalysisParams(min_components=3).to_dict()["min_components"] == 3


def test_default_floor_is_recorded_and_selects_rank_two(result) -> None:
    record = result.record
    assert record["params"]["min_components"] == 2
    assert record["latent"]["integration_params"]["min_components"] == 2
    assert record["latent"]["selected_lv"] >= 2
    assert result.latent_scores.shape[1] >= 2


def test_floor_of_one_omits_the_key_from_params_and_latent() -> None:
    out = analyze_trajectories(make_dataset(), TrajectoryAnalysisParams(min_components=1, **SMALL_CV))
    assert "min_components" not in out.record["params"]
    assert "min_components" not in out.record["latent"]["integration_params"]
    assert out.record["latent"]["selected_lv"] >= 1


def test_named_floor_is_honoured_and_forced_rank_ignores_it() -> None:
    three = analyze_trajectories(make_dataset(), TrajectoryAnalysisParams(min_components=3, **SMALL_CV))
    assert three.record["latent"]["selected_lv"] >= 3
    assert three.record["latent"]["integration_params"]["min_components"] == 3
    forced = analyze_trajectories(make_dataset(), TrajectoryAnalysisParams(forced_components=2, min_components=3))
    assert forced.record["latent"]["selected_lv"] == 2
    assert forced.record["latent"]["component_selection"] == "forced"
    assert forced.record["params"]["min_components"] == 3


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


def _run_pair(dataset, *, estimator=None, layers=None, floor=None):
    """The harness and the analysis on one dataset; ``floor=None`` leaves both at their defaults."""

    integration_params = dict(SMALL_CV)
    if layers is not None:
        integration_params["layers"] = list(layers)
    harness_kwargs = {} if estimator is None else {"orientation_estimator": estimator}
    floor_kwargs = {} if floor is None else {"min_components": floor}
    harness = evaluate_semisynthetic_trajectory(
        dataset,
        SimulationEvaluationParams(
            integration_method="pls",
            integration_params=integration_params,
            permutations=49,
            seed=11,
            **harness_kwargs,
            **floor_kwargs,
        ),
    )
    analysis_kwargs = {} if estimator is None else {"orientation_estimator": estimator}
    analysis = analyze_trajectories(
        _as_trajectory_dataset(dataset, layers or ("methylation", "expression", "proteomics")),
        TrajectoryAnalysisParams(permutations=49, seed=11, **SMALL_CV, **analysis_kwargs, **floor_kwargs),
    )
    return harness, analysis


FLOORS = pytest.mark.parametrize("floor", [None, 1], ids=["default_floor", "floor_1"])


@FLOORS
def test_matches_harness_default_estimator(synthetic, floor) -> None:
    from motco.simulations.evaluation import integrate_semisynthetic_dataset

    harness, analysis = _run_pair(synthetic, floor=floor)
    _assert_equivalent(harness, analysis)
    floor_kwargs = {} if floor is None else {"min_components": floor}
    latent = integrate_semisynthetic_dataset(
        synthetic,
        SimulationEvaluationParams(integration_method="pls", integration_params=SMALL_CV, **floor_kwargs),
    )
    assert np.array_equal(latent.matrix.to_numpy(), analysis.latent_scores.to_numpy())
    assert analysis.record["measurement"]["orientation_estimator"] == "net_displacement"
    if floor is None:
        assert analysis.record["latent"]["selected_lv"] >= 2
        assert harness.runtime_metadata["min_components"] == 2
    else:
        assert "min_components" not in harness.runtime_metadata


@FLOORS
def test_matches_harness_pc1(synthetic, floor) -> None:
    harness, analysis = _run_pair(synthetic, estimator="pc1", floor=floor)
    _assert_equivalent(harness, analysis)


@FLOORS
def test_matches_harness_layer_selection(synthetic, floor) -> None:
    from motco.simulations.evaluation import integrate_semisynthetic_dataset

    layers = ("methylation", "expression")
    harness, analysis = _run_pair(synthetic, layers=layers, floor=floor)
    _assert_equivalent(harness, analysis)
    floor_kwargs = {} if floor is None else {"min_components": floor}
    latent = integrate_semisynthetic_dataset(
        synthetic,
        SimulationEvaluationParams(
            integration_method="pls", integration_params={**SMALL_CV, "layers": list(layers)}, **floor_kwargs
        ),
    )
    assert np.array_equal(latent.matrix.to_numpy(), analysis.latent_scores.to_numpy())


# ── Covariates ───────────────────────────────────────────────────────────────


def with_covariates(seed: int = 5, **kwargs) -> tuple[pd.DataFrame, dict[str, pd.DataFrame]]:
    metadata, blocks = make_real_like(**kwargs)
    rng = np.random.default_rng(seed)
    n = metadata.shape[0]
    metadata["age"] = 75.0 + 5.0 * rng.normal(size=n)
    metadata["pmi"] = rng.uniform(2.0, 12.0, size=n)
    metadata["education_level"] = rng.choice(["College", "Graduate", "HighSchool"], size=n)
    return metadata, blocks


COVARIATES = ("age", "education_level")


def test_covariate_naming_the_group_column_is_rejected() -> None:
    metadata, blocks = with_covariates()
    with pytest.raises(TrajectoryAnalysisError, match="must differ from the sample, group and stage columns"):
        make_dataset(metadata, blocks, covariate_cols=["age", "Sex"])


@pytest.mark.parametrize(
    ("kwargs", "message"),
    [
        ({"covariate_cols": ["rin"]}, "missing covariate column.*'rin'"),
        ({"covariate_cols": ["age", "age"]}, "declared once; repeated: \\['age'\\]"),
        ({"covariate_cols": ["age"], "categorical_covariates": ["pmi"]}, "not declared: \\['pmi'\\]"),
        ({"covariate_cols": ["age"], "missing_covariates": "impute"}, "Unknown missing_covariates policy 'impute'"),
    ],
    ids=["absent", "repeated", "categorical-undeclared", "policy"],
)
def test_covariate_declaration_rejections(kwargs: dict, message: str) -> None:
    metadata, blocks = with_covariates()
    with pytest.raises(TrajectoryAnalysisError, match=message):
        make_dataset(metadata, blocks, **kwargs)


def test_constant_covariate_is_rejected() -> None:
    metadata, blocks = make_real_like()  # age is 70 for everyone
    with pytest.raises(TrajectoryAnalysisError, match="'age' is constant"):
        make_dataset(metadata, blocks, covariate_cols=["age"])


def test_missing_covariate_is_an_error_by_default() -> None:
    metadata, blocks = with_covariates()
    metadata.loc[4, "pmi"] = np.nan
    with pytest.raises(TrajectoryAnalysisError, match=f"'pmi' has missing values.*'{metadata.loc[4, 'donor']}'"):
        make_dataset(metadata, blocks, covariate_cols=["age", "pmi"])


def test_complete_case_drops_the_sample_everywhere() -> None:
    metadata, blocks = with_covariates(n_per_cell=14)
    metadata = metadata.iloc[4:].reset_index(drop=True)
    blocks = {name: frame.iloc[4:] for name, frame in blocks.items()}
    assert metadata.shape[0] == 80
    metadata["rin"] = np.linspace(5.0, 9.0, 80)
    dropped = [metadata.loc[10, "donor"], metadata.loc[50, "donor"]]
    metadata.loc[[10, 50], "rin"] = np.nan
    dataset = make_dataset(metadata, blocks, covariate_cols=["rin"], missing_covariates="complete_case")
    assert dataset.dropped_sample_ids == tuple(dropped)
    assert len(dataset.sample_ids) == 78
    assert all(block.values.shape[0] == 78 for block in dataset.blocks)
    out = analyze_trajectories(dataset, TrajectoryAnalysisParams(forced_components=2))
    assert out.latent_scores.shape[0] == 78
    assert not set(dropped) & set(out.latent_scores.index)
    cells = out.record["design"]["cell_sizes"]
    assert sum(size for stages in cells.values() for size in stages.values()) == 78
    assert out.record["design"]["n_samples"] == 78
    assert out.record["covariates"]["n_dropped"] == 2
    assert out.record["covariates"]["dropped_sample_ids"] == dropped
    assert out.record["covariates"]["missing_policy"] == "complete_case"


def test_complete_case_that_empties_a_cell_is_rejected() -> None:
    metadata, blocks = with_covariates()
    metadata.loc[(metadata["Sex"] == "M") & (metadata["ADNC"] == "High"), "pmi"] = np.nan
    with pytest.raises(TrajectoryAnalysisError, match="emptied group/stage combination.*\\('M', 'High'\\)"):
        make_dataset(metadata, blocks, covariate_cols=["pmi"], missing_covariates="complete_case")


def test_attribution_with_covariates_is_rejected_before_fitting(monkeypatch) -> None:
    import motco.analysis as analysis

    def forbidden(*args, **kwargs):
        raise AssertionError("the latent space must not be fitted")

    monkeypatch.setattr(analysis, "fit_pls_latent_space", forbidden)
    metadata, blocks = with_covariates()
    with pytest.raises(TrajectoryAnalysisError, match="Attribution is not available with covariates"):
        analyze_trajectories(
            make_dataset(metadata, blocks, covariate_cols=COVARIATES),
            TrajectoryAnalysisParams(attribution=AttributionDiagnosticSettings(enabled=True)),
        )


def test_collinear_covariate_is_rejected() -> None:
    metadata, blocks = with_covariates()
    metadata["is_m"] = (metadata["Sex"] == "M").astype(float)
    with pytest.raises(TrajectoryAnalysisError, match="covariates \\['age', 'is_m'\\] is not of full rank"):
        analyze_trajectories(
            make_dataset(metadata, blocks, covariate_cols=["age", "is_m"]),
            TrajectoryAnalysisParams(forced_components=2),
        )


@pytest.fixture(scope="module")
def adjusted():
    metadata, blocks = with_covariates()
    return analyze_trajectories(
        make_dataset(metadata, blocks, covariate_cols=COVARIATES),
        TrajectoryAnalysisParams(permutations=19, seed=3, **SMALL_CV),
    )


def test_latent_space_is_unchanged_by_covariates(adjusted) -> None:
    metadata, blocks = with_covariates()
    plain = analyze_trajectories(make_dataset(metadata, blocks), TrajectoryAnalysisParams(**SMALL_CV))
    assert np.array_equal(plain.latent_scores.to_numpy(), adjusted.latent_scores.to_numpy())
    pd.testing.assert_frame_equal(plain.pls_loadings, adjusted.pls_loadings)
    assert plain.record["latent"]["selected_lv"] == adjusted.record["latent"]["selected_lv"]
    assert "covariates" not in plain.record


def test_covariate_record_layout(adjusted, tmp_path) -> None:
    record = adjusted.record
    keys = list(record)
    assert keys.index("covariates") == keys.index("design") + 1
    block = record["covariates"]
    assert block["columns"] == ["age", "education_level"]
    assert block["kinds"] == {"age": "numeric", "education_level": "categorical"}
    assert block["levels"] == {"education_level": ["College", "Graduate", "HighSchool"]}
    assert block["reference_levels"] == {"education_level": "College"}
    assert block["design_columns"] == ["age", "education_level[Graduate]", "education_level[HighSchool]"]
    assert set(block["reference_row"]) == set(block["design_columns"])
    assert block["reference_row"]["age"] == pytest.approx(with_covariates()[0]["age"].mean())
    assert block["lsmean_reference"] == "pooled_mean"
    assert block["missing_policy"] == "error"
    assert block["n_dropped"] == 0 and block["dropped_sample_ids"] == []
    components = list(adjusted.latent_scores.columns)
    assert list(block["coefficients"]) == block["design_columns"]
    assert all(list(row) == components for row in block["coefficients"].values())
    assert record["schema_version"] == 1
    loaded = json.loads(write_analysis_outputs(adjusted, tmp_path)["result"].read_text())
    assert loaded["covariates"]["design_columns"] == block["design_columns"]


def test_adjusted_p_values_are_reproducible(adjusted) -> None:
    metadata, blocks = with_covariates()
    again = analyze_trajectories(
        make_dataset(metadata, blocks, covariate_cols=COVARIATES),
        TrajectoryAnalysisParams(permutations=19, seed=3, **SMALL_CV),
    )
    assert again.record["p_values"] == adjusted.record["p_values"]
    assert again.record["null_summary"] == adjusted.record["null_summary"]
    assert set(adjusted.record["p_values"]) == {"delta", "angle", "shape"}


def _cell_design(groups: np.ndarray, stages: np.ndarray):
    from motco.simulations.evaluation import build_simulation_trajectory_design

    return build_simulation_trajectory_design(
        pd.DataFrame({"group": groups.astype(str), "stage": stages.astype(str)}), group_col="group", stage_col="stage"
    )


def _confounded_latent(seed: int = 0):
    """Pure cell means plus ``z·v`` with ``z`` depending on group × stage."""

    rng = np.random.default_rng(seed)
    groups = np.repeat([0, 1], 30)
    stages = np.tile(np.repeat([0, 1, 2], 10), 2)
    cell_means = rng.normal(size=(2, 3, 4))
    Y0 = cell_means[groups, stages]
    z = 0.8 * groups * stages + 0.4 * stages + rng.normal(size=groups.size)
    v = rng.normal(size=4)
    return groups, stages, Y0, z, v


def test_statistics_do_not_depend_on_the_reference_row() -> None:
    from dataclasses import replace

    from motco.analysis import _covariate_adjusted_design
    from motco.stats.design import build_covariate_matrix
    from motco.stats.trajectory import estimate_betas, estimate_difference

    groups, stages, Y0, z, v = _confounded_latent()
    Y = Y0 + np.outer(z, v) + np.random.default_rng(1).normal(scale=0.3, size=Y0.shape)
    design = _cell_design(groups, stages)
    covariates = build_covariate_matrix(pd.DataFrame({"z": z}), ["z"])
    shifted = replace(covariates, reference_row=covariates.reference_row + 3.7)
    outputs = []
    for cov in (covariates, shifted):
        full, _, ls = _covariate_adjusted_design(design.model_full, design.model_reduced, design.ls_means, cov)
        stats = estimate_difference(Y, full, ls, design.contrast, orientation_estimator="net_displacement")
        means = ls @ np.asarray(estimate_betas(full, Y))
        outputs.append((stats, means))
    (stats_a, means_a), (stats_b, means_b) = outputs
    for a, b in zip(stats_a, stats_b):
        np.testing.assert_allclose(a, b, atol=1e-9)
    shift = means_b - means_a
    np.testing.assert_allclose(shift, np.tile(shift[0], (shift.shape[0], 1)), atol=1e-12)
    assert np.linalg.norm(shift[0]) > 0


def test_rrpp_p_values_do_not_depend_on_the_reference_row() -> None:
    from dataclasses import replace

    from motco.analysis import _covariate_adjusted_design
    from motco.simulations.evaluation import _empirical_p_value
    from motco.stats.design import build_covariate_matrix
    from motco.stats.permutation import RRPP
    from motco.stats.trajectory import estimate_difference

    groups, stages, _, z, v = _confounded_latent()
    rng = np.random.default_rng(2)
    # No group difference beyond the covariate, so the p-values sit inside (1/(B + 1), 1].
    Y = rng.normal(size=(3, 4))[stages] + np.outer(z, v) + rng.normal(scale=0.5, size=(groups.size, 4))
    design = _cell_design(groups, stages)
    covariates = build_covariate_matrix(pd.DataFrame({"z": z}), ["z"])
    shifted = replace(covariates, reference_row=covariates.reference_row + 3.7)
    kwargs = {"orientation_estimator": "net_displacement"}
    outputs = []
    for cov in (covariates, shifted):
        full, reduced, ls = _covariate_adjusted_design(design.model_full, design.model_reduced, design.ls_means, cov)
        observed = estimate_difference(Y, full, ls, design.contrast, **kwargs)
        nulls = RRPP(Y, full, reduced, ls, design.contrast, permutations=99, n_jobs=1, progress=False, seed=11,
                     **kwargs)
        p_values = [
            _empirical_p_value([draw[0, 1] for draw in null], statistic[0, 1])
            for statistic, null in zip(observed, nulls)
        ]
        outputs.append((observed, nulls, p_values))
    (observed_a, nulls_a, p_a), (observed_b, nulls_b, p_b) = outputs
    for a, b in zip(observed_a, observed_b):
        np.testing.assert_allclose(a, b, atol=1e-9)
    for null_a, null_b in zip(nulls_a, nulls_b):
        np.testing.assert_allclose(np.asarray(null_a), np.asarray(null_b), atol=1e-9)
    assert p_a == p_b
    assert any(0.01 < p < 1.0 for p in p_a)


def test_adjustment_removes_a_known_covariate_effect() -> None:
    from motco.analysis import _covariate_adjusted_design
    from motco.stats.design import build_covariate_matrix
    from motco.stats.trajectory import estimate_difference

    groups, stages, Y0, z, v = _confounded_latent()
    Y = Y0 + np.outer(z, v)
    design = _cell_design(groups, stages)
    covariates = build_covariate_matrix(pd.DataFrame({"z": z}), ["z"])
    full, _, ls = _covariate_adjusted_design(design.model_full, design.model_reduced, design.ls_means, covariates)
    kwargs = {"orientation_estimator": "net_displacement"}
    truth = estimate_difference(Y0, design.model_full, design.ls_means, design.contrast, **kwargs)
    adjusted_stats = estimate_difference(Y, full, ls, design.contrast, **kwargs)
    unadjusted = estimate_difference(Y, design.model_full, design.ls_means, design.contrast, **kwargs)
    for expected, got in zip(truth, adjusted_stats):
        np.testing.assert_allclose(got, expected, atol=1e-8)
    assert not np.allclose(unadjusted[0], truth[0], atol=1e-3)
    assert not np.allclose(unadjusted[1], truth[1], atol=1e-3)


# ── Latent adjustment ────────────────────────────────────────────────────────


def test_unknown_latent_adjustment_is_rejected() -> None:
    with pytest.raises(TrajectoryAnalysisError, match="Unknown covariate_latent_adjustment 'group'"):
        TrajectoryAnalysisParams(covariate_latent_adjustment="group")


def test_latent_adjustment_is_omitted_from_params_while_none() -> None:
    assert "covariate_latent_adjustment" not in TrajectoryAnalysisParams().to_dict()
    params = TrajectoryAnalysisParams(covariate_latent_adjustment="group_stage").to_dict()
    assert params["covariate_latent_adjustment"] == "group_stage"


def test_latent_adjustment_without_covariates_is_rejected_before_fitting(monkeypatch) -> None:
    import motco.analysis as analysis

    def forbidden(*args, **kwargs):
        raise AssertionError("the latent space must not be fitted")

    monkeypatch.setattr(analysis, "fit_pls_latent_space", forbidden)
    with pytest.raises(TrajectoryAnalysisError, match="requires declared covariates"):
        analyze_trajectories(make_dataset(), TrajectoryAnalysisParams(covariate_latent_adjustment="stage"))


@pytest.fixture(scope="module")
def latent_adjusted():
    metadata, blocks = with_covariates()
    return analyze_trajectories(
        make_dataset(metadata, blocks, covariate_cols=COVARIATES),
        TrajectoryAnalysisParams(permutations=19, seed=3, covariate_latent_adjustment="stage", **SMALL_CV),
    )


def test_latent_adjustment_changes_the_basis_only(monkeypatch) -> None:
    import motco.analysis as analysis

    designs: list[tuple] = []
    original = analysis.RRPP

    def recording(Y, model_full, model_reduced, LS_means, contrast, **kwargs):
        designs.append((np.array(model_full), np.array(model_reduced), np.array(LS_means), contrast))
        return original(Y, model_full, model_reduced, LS_means, contrast, **kwargs)

    monkeypatch.setattr(analysis, "RRPP", recording)
    metadata, blocks = with_covariates()
    dataset = make_dataset(metadata, blocks, covariate_cols=COVARIATES)
    plain = analyze_trajectories(dataset, TrajectoryAnalysisParams(permutations=9, seed=3, **SMALL_CV))
    adjusted = analyze_trajectories(
        dataset,
        TrajectoryAnalysisParams(permutations=9, seed=3, covariate_latent_adjustment="stage", **SMALL_CV),
    )
    (full_a, reduced_a, ls_a, contrast_a), (full_b, reduced_b, ls_b, contrast_b) = designs
    assert np.array_equal(full_a, full_b)
    assert np.array_equal(reduced_a, reduced_b)
    assert np.array_equal(ls_a, ls_b)
    assert contrast_a == contrast_b
    assert full_a.shape[1] == 6 + 3  # group × stage plus age and two education indicators
    assert not np.allclose(plain.latent_scores.to_numpy()[:, :1], adjusted.latent_scores.to_numpy()[:, :1])
    assert "latent_adjustment" not in plain.record["latent"]["integration_params"]
    assert adjusted.record["latent"]["integration_params"]["latent_adjustment"] == "stage"


def test_latent_adjustment_record_layout(latent_adjusted) -> None:
    record = latent_adjusted.record
    keys = list(record)
    assert keys.index("latent_adjustment") == keys.index("covariates") + 1
    block = record["latent_adjustment"]
    assert list(block) == [
        "conditioning",
        "covariates",
        "design_columns",
        "n_samples",
        "residual_df",
        "cv_min_training_df",
        "feature_scale",
        "feature_r2",
    ]
    assert block["conditioning"] == "stage"
    assert block["covariates"] == ["age", "education_level"]
    assert block["design_columns"] == [
        "stage[1]",
        "stage[2]",
        "age",
        "education_level[Graduate]",
        "education_level[HighSchool]",
    ]
    n, k, q = 48, 2, 3
    assert block["n_samples"] == n
    assert block["residual_df"] == n - 1 - k - q
    assert isinstance(block["cv_min_training_df"], int) and 0 < block["cv_min_training_df"] < block["residual_df"]
    assert block["feature_scale"] == "pooled_standardized"
    assert list(block["feature_r2"]) == ["q05", "q25", "q50", "q75", "q95", "mean", "fraction_above_0.5"]
    assert all(0.0 <= value <= 1.0 for value in block["feature_r2"].values())
    assert record["params"]["covariate_latent_adjustment"] == "stage"
    assert record["schema_version"] == 1
    # The covariate block keeps its 0.8.0 meaning: covariate effects on the projected raw features.
    coefficients = np.array([list(row.values()) for row in record["covariates"]["coefficients"].values()])
    assert np.any(np.abs(coefficients) > 1e-6)


def test_latent_adjustment_block_absent_without_the_setting(adjusted) -> None:
    assert "latent_adjustment" not in adjusted.record
    assert "latent_adjustment" not in adjusted.record["latent"]["integration_params"]
    assert "covariate_latent_adjustment" not in adjusted.record["params"]


def test_latent_adjustment_spec_design_on_80_samples() -> None:
    metadata, blocks = with_covariates(n_per_cell=14)
    metadata = metadata.iloc[4:].reset_index(drop=True)
    blocks = {name: frame.iloc[4:] for name, frame in blocks.items()}
    metadata["rin"] = np.random.default_rng(9).uniform(5.0, 9.0, size=80)
    out = analyze_trajectories(
        make_dataset(metadata, blocks, covariate_cols=["age", "pmi", "rin"]),
        TrajectoryAnalysisParams(covariate_latent_adjustment="stage", **SMALL_CV),
    )
    block = out.record["latent_adjustment"]
    assert block["design_columns"] == ["stage[1]", "stage[2]", "age", "pmi", "rin"]
    assert block["residual_df"] == 74


def test_latent_adjusted_p_values_are_reproducible(latent_adjusted) -> None:
    metadata, blocks = with_covariates()
    again = analyze_trajectories(
        make_dataset(metadata, blocks, covariate_cols=COVARIATES),
        TrajectoryAnalysisParams(permutations=19, seed=3, covariate_latent_adjustment="stage", **SMALL_CV),
    )
    assert again.record["p_values"] == latent_adjusted.record["p_values"]
    assert again.record["null_summary"] == latent_adjusted.record["null_summary"]
    assert np.array_equal(again.latent_scores.to_numpy(), latent_adjusted.latent_scores.to_numpy())


def test_latent_adjustment_collinear_with_stage_names_the_covariate(monkeypatch) -> None:
    import motco.analysis as analysis

    def forbidden(*args, **kwargs):
        raise AssertionError("the latent space must not be fitted")

    monkeypatch.setattr(analysis, "fit_pls_latent_space", forbidden)
    metadata, blocks = with_covariates()
    metadata["stage_index"] = metadata["ADNC"].map({stage: index for index, stage in enumerate(STAGES)}).astype(float)
    with pytest.raises(
        TrajectoryAnalysisError, match="residualization design with covariates \\['age', 'stage_index'\\].*not of full"
    ):
        analyze_trajectories(
            make_dataset(metadata, blocks, covariate_cols=["age", "stage_index"]),
            TrajectoryAnalysisParams(covariate_latent_adjustment="stage", **SMALL_CV),
        )


def test_latent_adjustment_fold_missing_a_level_names_the_fold() -> None:
    metadata, blocks = with_covariates()
    metadata["batch"] = "common"
    metadata.loc[5, "batch"] = "rare"
    with pytest.raises(TrajectoryAnalysisError, match="repeat 1, outer fold \\d+(, inner fold \\d+)?: .*not of full"):
        analyze_trajectories(
            make_dataset(metadata, blocks, covariate_cols=["age", "batch"]),
            TrajectoryAnalysisParams(covariate_latent_adjustment="group_stage", **SMALL_CV),
        )
