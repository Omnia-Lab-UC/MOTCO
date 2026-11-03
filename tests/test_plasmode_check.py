"""``scripts/plasmode_check.py``: contract, null-making step, plant, units, metrics, merge, calibration, decision.

Every test runs on a synthetic stand-in export written to ``tmp_path`` (the
astrocyte ``cohort_30`` cell table, 30 features per block); no real data.
"""

from __future__ import annotations

import dataclasses
import hashlib
import json
import math
import shutil
import sys
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pandas as pd
import pytest


def _script(name: str):
    sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
    try:
        return __import__(name)
    finally:
        sys.path.pop(0)


check = _script("plasmode_check")
TABLE = [[11, 10, 28], [9, 10, 12]]
N_FEATURES = 30
SMALL_CV = {"cv1_splits": 2, "cv2_splits": 3, "n_repeats": 2, "max_components": 5, "random_state": 1203}


# ── The stand-in export ──────────────────────────────────────────────────────


def write_export(directory: Path, seed: int = 7) -> Path:
    """A ``sead_sexdiff``-shaped export: manifest, rna/atac blocks, metadata and covariates."""

    rng = np.random.default_rng(seed)
    cells = [(group, stage) for g, group in enumerate(check.GROUPS) for s, stage in enumerate(check.STAGE_ORDER)
             for _ in range(TABLE[g][s])]
    order = rng.permutation(len(cells))
    donors = [f"H{index:03d}" for index in range(len(cells))]
    metadata = pd.DataFrame(
        {
            "donor_id": donors,
            "sex": [cells[i][0] for i in order],
            "stage": [cells[i][1] for i in order],
            "stage_label": [f"label {cells[i][1]}" for i in order],
        }
    )
    stage = metadata["stage"].astype(int).to_numpy()
    directory.mkdir(parents=True, exist_ok=True)
    blocks = {}
    for name in check.BLOCKS:
        loadings = rng.normal(size=N_FEATURES) * (rng.random(N_FEATURES) < 0.5)
        values = 5 + rng.normal(size=(len(donors), N_FEATURES)) + 1.2 * stage[:, None] * loadings
        blocks[name] = pd.DataFrame(values, index=pd.Index(donors, name="donor_id"),
                                    columns=[f"{name}_g{j}" for j in range(N_FEATURES)])
        blocks[name].to_csv(directory / f"{name}.csv", float_format="%.10g")
    metadata.to_csv(directory / "metadata.csv", index=False)
    covariates = pd.DataFrame(
        {
            "donor_id": donors,
            "age_at_death": rng.normal(85, 6, len(donors)).round(),
            "age_capped": False,
            "pmi": rng.uniform(3, 10, len(donors)),
            "rin": rng.uniform(6, 9, len(donors)),
            "years_education": rng.integers(10, 20, len(donors)).astype(float),
            "education_level": "Bachelors",
            "n_rna": rng.integers(100, 3000, len(donors)),
            "n_atac": rng.integers(30, 900, len(donors)),
        }
    ).sample(frac=1.0, random_state=3)
    covariates.to_csv(directory / "covariates.csv", index=False)
    manifest = {
        "cohort": "stand_in",
        "n_donors": len(donors),
        "blocks": {
            name: {"file": f"{name}.csv", "units": "log2-CPM", "motco_transform": "identity", "n_features": N_FEATURES}
            for name in check.BLOCKS
        },
        "metadata": {
            "file": "metadata.csv",
            "sample_col": "donor_id",
            "group_col": "sex",
            "group_labels": list(check.GROUPS),
            "stage_col": "stage",
            "stage_order": list(check.STAGE_ORDER),
        },
        "covariates": {"file": "covariates.csv"},
    }
    (directory / "manifest.json").write_text(json.dumps(manifest, indent=1), encoding="utf-8")
    return directory


@pytest.fixture(scope="module")
def export_dir(tmp_path_factory) -> Path:
    return write_export(tmp_path_factory.mktemp("export"))


@pytest.fixture(scope="module")
def export(export_dir):
    return check.load_export(export_dir)


@pytest.fixture(scope="module")
def stand_in_fidelity(export) -> dict:
    """The stand-in's own label-free latent space at the small CV (the real gate is rank 5, AUROC 0.635556)."""

    settings = {**check.base_settings(export), "cv": dict(SMALL_CV), "attribution": None}
    built = check.build_unit(0.0, 0, {**settings, "stream": "fidelity"}, export)
    _, result = check.analyse(
        built["blocks"], built["metadata"], export, check.analysis_params(settings, None, permutations=0)
    )
    latent = result.record["latent"]
    return {"selected_lv": int(latent["selected_lv"]), "cv_mean_auroc": round(float(latent["cv_mean_auroc"]), 6)}


def small_settings(export, fidelity, *, reps: int = 2, perms: int = 9, bootstrap: int = 5, **strength) -> dict:
    strength = strength or {"e_star": 1.0}
    settings = check.build_settings(
        export, fidelity, replicates=reps, permutations=perms, bootstrap=bootstrap, **strength
    )
    settings["cv"] = dict(SMALL_CV)
    return settings


def _copy_export(source: Path, target: Path) -> Path:
    shutil.copytree(source, target)
    return target


# ── Settings and CLI ─────────────────────────────────────────────────────────


def test_help_lists_the_five_subcommands() -> None:
    text = check.build_parser().format_help()
    for command in ("fidelity", "calibrate", "run", "merge", "decide"):
        assert command in text


def test_settings_hold_every_constant_and_the_export_hashes(export) -> None:
    settings = check.build_settings(export, {"selected_lv": 5, "cv_mean_auroc": 0.635556}, e_star=1.0)
    assert settings["base_seed"] == 2034 and settings["ladder"] == [0.5, 1.0, 1.5, 2.0, 3.0]
    assert settings["attribution"]["bootstrap_replicates"] == 200 and settings["attribution"]["top_k"] == 20
    assert settings["plant"]["per_block"] == 10 and settings["plant"]["group"] == "Male"
    assert settings["covariates"]["declared"] == [
        "age_at_death", "pmi", "rin", "years_education", "log10_n_rna", "log10_n_atac"
    ]
    assert settings["cv"]["random_state"] == 1203 and settings["residual_correction"] == "df"
    assert settings["min_components"] == 2 and settings["orientation_estimator"] == "net_displacement"
    assert settings["fidelity_gate"] == {"selected_lv": 5, "cv_mean_auroc": 0.635556, "decimals": 6}
    assert settings["rule"]["thresholds"] == {"planted_precision": 0.5, "sign_stability": 0.9, "top_k_frequency": 0.5}
    assert settings["rule"]["detection_floor"] == 0.8 and settings["rule"]["se_multiplier"] == 2.0
    assert settings["replicates"] == 500 and settings["permutations"] == 999
    assert settings["cells"] == [0.0, 1.0] and settings["strength"]["source"] == "explicit"
    assert settings["export"]["group_stage_sizes"] == TABLE
    files = {"manifest.json", "metadata.csv", "covariates.csv", "rna.csv", "atac.csv"}
    assert set(settings["export"]["files"]) == files
    with pytest.raises(SystemExit, match="exactly one"):
        check.build_settings(export, {}, e_star=None)


# ── Input contract ───────────────────────────────────────────────────────────


def test_derived_covariates_are_log10_of_the_counts(export, export_dir) -> None:
    raw = pd.read_csv(export_dir / "covariates.csv").set_index("donor_id").loc[export.metadata["donor_id"]]
    np.testing.assert_array_equal(export.metadata["log10_n_rna"], np.log10(raw["n_rna"].to_numpy(dtype=float)))
    np.testing.assert_array_equal(export.metadata["log10_n_atac"], np.log10(raw["n_atac"].to_numpy(dtype=float)))
    assert list(export.metadata.columns[3:]) == list(check.COVARIATES)
    assert export.metadata["donor_id"].tolist() == sorted(export.metadata["donor_id"])
    for name in check.BLOCKS:
        assert export.blocks[name].index.tolist() == export.metadata["donor_id"].tolist()


def test_recorded_hashes_are_the_files_sha256(export, export_dir) -> None:
    for name, digest in export.files.items():
        assert digest == hashlib.sha256((export_dir / name).read_bytes()).hexdigest()


def _break(directory: Path, violation: str) -> None:
    if violation == "donor_mismatch":
        frame = pd.read_csv(directory / "rna.csv", index_col=0)
        frame.iloc[1:].to_csv(directory / "rna.csv")
    elif violation in ("unknown_stage", "group_count"):
        metadata = pd.read_csv(directory / "metadata.csv", dtype=str)
        column, value = ("stage", "3") if violation == "unknown_stage" else ("sex", "Other")
        metadata.loc[0, column] = value
        metadata.to_csv(directory / "metadata.csv", index=False)
    elif violation == "missing_covariate":
        covariates = pd.read_csv(directory / "covariates.csv")
        covariates.loc[4, "pmi"] = np.nan
        covariates.to_csv(directory / "covariates.csv", index=False)
    elif violation == "missing_column":
        covariates = pd.read_csv(directory / "covariates.csv").drop(columns=["rin"])
        covariates.to_csv(directory / "covariates.csv", index=False)


@pytest.mark.parametrize(
    ("violation", "error", "match"),
    [
        ("donor_mismatch", "DonorMismatchError", "different donors"),
        ("unknown_stage", "UnknownStageError", r"stage label\(s\) \['3'\]"),
        ("group_count", "GroupCountError", "Other"),
        ("missing_covariate", "MissingCovariateError", "'pmi' is missing"),
        ("missing_column", "MissingCovariateError", "rin"),
    ],
)
def test_each_contract_violation_raises_its_named_error(export_dir, tmp_path, violation, error, match) -> None:
    directory = _copy_export(export_dir, tmp_path / "export")
    _break(directory, violation)
    with pytest.raises(getattr(check, error), match=match):
        check.load_export(directory)
    with pytest.raises(SystemExit, match="export contract violation"):
        check.main(["fidelity", "--export-dir", str(directory), "--out-dir", str(tmp_path / "out")])
    assert not (tmp_path / "out").exists()


def test_manifest_must_declare_the_primary_blocks_and_design(export_dir, tmp_path) -> None:
    directory = _copy_export(export_dir, tmp_path / "export")
    manifest = json.loads((directory / "manifest.json").read_text())
    manifest["metadata"]["stage_order"] = ["2", "1", "0"]
    (directory / "manifest.json").write_text(json.dumps(manifest))
    with pytest.raises(check.ManifestError, match="stage order"):
        check.load_export(directory)


# ── Fidelity gate ────────────────────────────────────────────────────────────


def _fake_analysis(rank: int, auroc: float):
    record = {"latent": {"selected_lv": rank, "cv_mean_auroc": auroc, "component_selection": "cv"}, "software": {}}
    return lambda export, n_jobs=1: SimpleNamespace(record=record)


def test_fidelity_match_writes_the_gate_result(export, tmp_path) -> None:
    payload = check.run_fidelity(export, tmp_path, analysis=_fake_analysis(5, 0.63555612))
    stored = json.loads((tmp_path / "fidelity.json").read_text())
    assert stored["passed"] is True and stored["selected_lv"] == 5 and stored["cv_mean_auroc_rounded"] == 0.635556
    assert stored["export_files"] == export.files and payload["gate"] == check.FIDELITY_GATE
    assert check.check_fidelity(tmp_path / "fidelity.json", export) == {"selected_lv": 5, "cv_mean_auroc": 0.635556}


@pytest.mark.parametrize(("rank", "auroc"), [(4, 0.635556), (5, 0.6355549), (3, 0.61)])
def test_fidelity_mismatch_names_both_values_and_writes_nothing(export, export_dir, tmp_path, monkeypatch, rank, auroc):
    with pytest.raises(check.FidelityError, match=rf"rank {rank} \(gate 5\), AUROC {auroc:.6f} \(gate 0.635556\)"):
        check.run_fidelity(export, tmp_path, analysis=_fake_analysis(rank, auroc))
    assert not (tmp_path / "fidelity.json").exists()
    monkeypatch.setattr(check, "fidelity_analysis", _fake_analysis(rank, auroc))
    with pytest.raises(SystemExit, match="fidelity gate failed"):
        check.main(["fidelity", "--export-dir", str(export_dir), "--out-dir", str(tmp_path)])
    with pytest.raises(SystemExit, match="run fidelity first"):
        check.main(["calibrate", "--export-dir", str(export_dir), "--fidelity", str(tmp_path / "fidelity.json"),
                    "--out-dir", str(tmp_path)])
    assert sorted(path.name for path in tmp_path.iterdir()) == []


def test_fidelity_from_another_export_is_refused(export, tmp_path) -> None:
    check.run_fidelity(export, tmp_path, analysis=_fake_analysis(5, 0.635556))
    other = check.load_export(write_export(tmp_path / "other", seed=8))
    with pytest.raises(SystemExit, match="another export"):
        check.check_fidelity(tmp_path / "fidelity.json", other)


# ── Null-making step ─────────────────────────────────────────────────────────


def test_permutation_keeps_the_cell_table(export) -> None:
    settings = small_settings(export, {"selected_lv": 2, "cv_mean_auroc": 0.5})
    for replicate in range(5):
        built = check.build_unit(1.0, replicate, settings, export)
        labels = built["metadata"]["sex"].to_numpy()
        assert export.group_stage_sizes(labels) == TABLE
        assert not np.array_equal(labels, export.real_labels)


def test_covariates_and_blocks_stay_on_their_donors(export) -> None:
    settings = small_settings(export, {"selected_lv": 2, "cv_mean_auroc": 0.5})
    built = check.build_unit(0.0, 3, settings, export)
    pd.testing.assert_frame_equal(
        built["metadata"].drop(columns=["sex"]), export.metadata.drop(columns=["sex"]), check_exact=True
    )
    for name in check.BLOCKS:
        assert built["blocks"][name] is export.blocks[name]


def test_cells_of_a_replicate_are_matched(export) -> None:
    settings = small_settings(export, {"selected_lv": 2, "cv_mean_auroc": 0.5})
    built = [check.build_unit(effect, 4, settings, export) for effect in (0.0, 1.0, 3.0)]
    for other in built[1:]:
        np.testing.assert_array_equal(other["metadata"]["sex"], built[0]["metadata"]["sex"])
        assert other["seeds"] == built[0]["seeds"] and other["plant"] == built[0]["plant"]
    assert set(built[0]["seeds"]) == {"labels", "plant", "rrpp", "bootstrap"}


def test_paper_and_calibration_streams_differ(export) -> None:
    paper = check.unit_seeds(2034, "paper", 0)
    calibration = check.unit_seeds(2034, "calibration", 0)
    assert not set(paper.values()) & set(calibration.values())
    assert not set(paper.values()) & set(check.unit_seeds(2034, "paper", 1).values())
    settings = small_settings(export, {"selected_lv": 2, "cv_mean_auroc": 0.5})
    a = check.build_unit(1.0, 0, settings, export)
    b = check.build_unit(1.0, 0, {**settings, "stream": "calibration"}, export)
    assert a["plant"] != b["plant"] or not np.array_equal(a["metadata"]["sex"], b["metadata"]["sex"])


def test_a_draw_equal_to_the_real_labels_is_redrawn_and_counted(export, monkeypatch) -> None:
    original = check.permute_within_stage
    calls = []

    def first_is_real(labels, stages, stage_order, rng):
        calls.append(1)
        drawn = original(labels, stages, stage_order, rng)
        return np.asarray(labels, dtype=object).copy() if len(calls) == 1 else drawn

    monkeypatch.setattr(check, "permute_within_stage", first_is_real)
    settings = small_settings(export, {"selected_lv": 2, "cv_mean_auroc": 0.5})
    built = check.build_unit(1.0, 2, settings, export)
    assert built["label_redraws"] == 1 and len(calls) == 2
    rng = np.random.default_rng(check.labels_seed(2034, "paper", 2, 1))
    expected = original(export.real_labels, export.stages, export.stage_order, rng)
    np.testing.assert_array_equal(built["metadata"]["sex"].to_numpy(), expected.astype(str))


# ── Plant ────────────────────────────────────────────────────────────────────


def test_pooled_within_stage_sd_matches_hand_computation(export) -> None:
    frame = export.blocks["rna"]
    stages = export.stages
    column = frame.iloc[:, 3].to_numpy()
    squares = sum(((column[stages == s] - column[stages == s].mean()) ** 2).sum() for s in check.STAGE_ORDER)
    assert export.sigma()["rna"][3] == pytest.approx(math.sqrt(squares / (len(column) - 3)), rel=1e-12)
    assert check.stage_weights(check.STAGE_ORDER) == {"0": 0.0, "1": 0.5, "2": 1.0}


def test_shift_formula_on_a_middle_stage_planted_donor(export) -> None:
    settings = small_settings(export, {"selected_lv": 2, "cv_mean_auroc": 0.5})
    built = check.build_unit(2.0, 1, settings, export)
    labels = built["metadata"]["sex"].to_numpy()
    for name in check.BLOCKS:
        indices = built["plant"][name]["indices"]
        signs = np.asarray(built["plant"][name]["signs"], dtype=float)
        planted, original = built["blocks"][name].to_numpy(), export.blocks[name].to_numpy()
        middle = np.flatnonzero((labels == "Male") & (export.stages == "1"))
        assert middle.size == 10
        expected = original.copy()
        expected[np.ix_(middle, indices)] = original[np.ix_(middle, indices)] + export.sigma()[name][indices] * signs
        last = np.flatnonzero((labels == "Male") & (export.stages == "2"))
        expected[np.ix_(last, indices)] = original[np.ix_(last, indices)] + 2.0 * export.sigma()[name][indices] * signs
        np.testing.assert_array_equal(planted, expected)
        untouched = np.ones(planted.shape, dtype=bool)
        untouched[np.ix_(np.concatenate([middle, last]), indices)] = False
        np.testing.assert_array_equal(planted[untouched], original[untouched])


def test_zero_strength_is_the_export_byte_for_byte(export) -> None:
    settings = small_settings(export, {"selected_lv": 2, "cv_mean_auroc": 0.5})
    built = check.build_unit(0.0, 1, settings, export)
    for name in check.BLOCKS:
        pd.testing.assert_frame_equal(built["blocks"][name], export.blocks[name], check_exact=True)


def test_plant_is_identical_across_strengths_and_drawn_per_block(export) -> None:
    settings = small_settings(export, {"selected_lv": 2, "cv_mean_auroc": 0.5})
    plants = [check.build_unit(effect, 6, settings, export)["plant"] for effect in (0.0, 1.5, 3.0)]
    assert plants[0] == plants[1] == plants[2]
    for name in check.BLOCKS:
        assert len(set(plants[0][name]["features"])) == 10 and set(plants[0][name]["signs"]) <= {-1, 1}
        assert all(feature.startswith(f"{name}_") for feature in plants[0][name]["features"])


def test_plant_cosine_matches_hand_computation(export) -> None:
    plant = check.draw_plant({name: list(export.blocks[name].columns) for name in check.BLOCKS}, 11, 10)
    values = np.concatenate([export.blocks[name].to_numpy() for name in check.BLOCKS], axis=1)
    z = (values - values.mean(axis=0)) / values.std(axis=0)
    d = z[export.stages == "2"].mean(axis=0) - z[export.stages == "0"].mean(axis=0)
    w = np.zeros(values.shape[1])
    for offset, name in ((0, "rna"), (N_FEATURES, "atac")):
        idx = np.asarray(plant[name]["indices"])
        scale = values.std(axis=0)[offset + idx]
        w[offset + idx] = export.sigma()[name][idx] * np.asarray(plant[name]["signs"]) / scale
    assert check.plant_cosine(export, plant) == pytest.approx(w @ d / (np.linalg.norm(w) * np.linalg.norm(d)))


# ── Metrics ──────────────────────────────────────────────────────────────────


def _bootstrap(sign, frequency, valid=10):
    return SimpleNamespace(
        valid_replicates=valid,
        sign_stability=np.asarray(sign, dtype=float),
        top_k_selection_frequency=np.asarray(frequency, dtype=float),
    )


NAMES = [f"f{j}" for j in range(6)]
CONTRAST = np.array([0.9, -0.8, 0.1, 0.7, 0.0, -0.05])
SIGN = [1.0, 0.8, 0.5, 0.9, 0.7, 0.6]
FREQUENCY = [1.0, 0.5, 0.0, 0.7, 0.2, 0.1]


def _planted(*indices: int, signs: tuple[int, ...]) -> tuple[np.ndarray, np.ndarray]:
    mask = np.zeros(6, dtype=bool)
    sign_vector = np.zeros(6, dtype=int)
    mask[list(indices)] = True
    sign_vector[list(indices)] = signs
    return mask, sign_vector


def test_contrast_metrics_match_hand_counts() -> None:
    planted, signs = _planted(0, 1, 4, signs=(1, 1, -1))
    metrics = check.contrast_metrics(
        CONTRAST, planted=planted, signs=signs, bootstrap=_bootstrap(SIGN, FREQUENCY), feature_names=NAMES, top_k=3
    )
    # top-3 = f0, f1, f3; K = {f0, f1, f4}; hits f0 (+, planted +) and f1 (−, planted +).
    assert metrics["planted_precision"] == pytest.approx(2 / 3)
    assert metrics["sign_agreement"] == pytest.approx(0.5)
    assert metrics["sign_stability"] == pytest.approx((1.0 + 0.8 + 0.7) / 3)
    assert metrics["top_k_frequency"] == pytest.approx((1.0 + 0.5 + 0.2) / 3)
    assert metrics["valid_replicates"] == 10 and metrics["planted_size"] == 3 and metrics["available"] is True


def test_planted_set_outside_the_top_k_keeps_stability_defined() -> None:
    planted, signs = _planted(2, 5, signs=(1, -1))
    metrics = check.contrast_metrics(
        CONTRAST, planted=planted, signs=signs, bootstrap=_bootstrap(SIGN, FREQUENCY), feature_names=NAMES, top_k=3
    )
    assert metrics["planted_precision"] == 0.0 and metrics["sign_agreement"] is None
    assert metrics["sign_stability"] == pytest.approx((0.5 + 0.6) / 2)
    assert metrics["top_k_frequency"] == pytest.approx((0.0 + 0.1) / 2)


def test_precision_divides_by_the_smaller_of_k_and_the_planted_set() -> None:
    planted, signs = _planted(0, 1, 3, signs=(1, -1, 1))
    metrics = check.contrast_metrics(
        CONTRAST, planted=planted, signs=signs, bootstrap=_bootstrap(SIGN, FREQUENCY), feature_names=NAMES, top_k=2
    )
    assert metrics["planted_precision"] == 1.0 and metrics["sign_agreement"] == 1.0


def test_unavailable_contrast_keeps_the_bootstrap_readings() -> None:
    planted, signs = _planted(0, 1, signs=(1, 1))
    metrics = check.contrast_metrics(
        None, planted=planted, signs=signs, bootstrap=_bootstrap(SIGN, FREQUENCY, valid=0), feature_names=NAMES, top_k=3
    )
    assert metrics["available"] is False and metrics["planted_precision"] is None
    assert metrics["sign_stability"] == pytest.approx(0.9) and metrics["valid_replicates"] == 0


def test_planted_vectors_name_the_joint_features() -> None:
    plant = {
        "rna": {"features": ["a"], "signs": [-1], "indices": [0]},
        "atac": {"features": ["b"], "signs": [1], "indices": [1]},
    }
    mask, signs = check.planted_vectors(["rna__a", "rna__b", "atac__a", "atac__b"], plant)
    assert mask.tolist() == [True, False, False, True] and signs.tolist() == [-1, 0, 0, 1]
    with pytest.raises(ValueError, match="not all among"):
        check.planted_vectors(["rna__a"], plant)


# ── Units and consistency checks ─────────────────────────────────────────────


@pytest.fixture(scope="module")
def analysed(export, stand_in_fidelity):
    settings = small_settings(export, stand_in_fidelity)
    built = check.build_unit(1.0, 0, settings, export)
    params = check.analysis_params(settings, built["seeds"], permutations=9)
    return check.analyse(built["blocks"], built["metadata"], export, params)


def test_analysis_runs_the_case_study_configuration(analysed) -> None:
    dataset, result = analysed
    assert list(dataset.covariate_cols) == list(check.COVARIATES)
    params = result.record["params"]
    assert params["min_components"] == 2 and params["residual_correction"] == "df"
    assert params["orientation_estimator"] == "net_displacement"
    assert result.record["measurement"]["residual_correction"]["kind"] == "df"
    assert result.attribution.config.mean_source == "covariate_adjusted"
    check.check_consistency(result)


@pytest.mark.parametrize(
    ("breakage", "match"),
    [
        ("latent_scores", "latent_scores"),
        ("stage_means", "projection"),
        ("rank", "anchor_latent"),
    ],
)
def test_a_failed_check_raises_and_writes_no_record(
    analysed, export, stand_in_fidelity, tmp_path, monkeypatch, breakage, match
) -> None:
    dataset, result = analysed
    if breakage == "latent_scores":
        broken = dataclasses.replace(result, latent_scores=result.latent_scores + 1.0)
    elif breakage == "stage_means":
        means = result.stage_means.copy()
        means[means.columns[2:]] = means[means.columns[2:]] * (1 + 1e-6)
        broken = dataclasses.replace(result, stage_means=means)
    else:
        record = json.loads(json.dumps(result.record))
        record["latent"]["selected_lv"] = stand_in_fidelity["selected_lv"] + 1
        broken = dataclasses.replace(result, record=record)
    monkeypatch.setattr(check, "analyse", lambda *args, **kwargs: (dataset, broken))
    settings = small_settings(export, stand_in_fidelity)
    with pytest.raises(check.ConsistencyError, match=match):
        check.run_shard(tmp_path, settings, export, 0, 1, limit=1)
    assert check.shard_path(tmp_path, 0, 1).read_text() == ""


def _strings(value) -> list:
    if isinstance(value, dict):
        return [*value.keys(), *(item for child in value.values() for item in _strings(child))]
    if isinstance(value, list):
        return [item for child in value for item in _strings(child)]
    return [value]


def test_end_to_end_unit_record_holds_no_donor_level_value(export, stand_in_fidelity) -> None:
    settings = small_settings(export, stand_in_fidelity)
    record = json.loads(json.dumps(check.evaluate_unit(1.0, 0, settings, export), allow_nan=False))
    leaves = _strings(record)
    donors = set(export.metadata["donor_id"])
    assert not donors & {leaf for leaf in leaves if isinstance(leaf, str)}

    def lists(value):
        if isinstance(value, dict):
            return [item for child in value.values() for item in lists(child)]
        if isinstance(value, list):
            return [value, *(item for child in value for item in lists(child))]
        return []

    assert all(len(items) < len(donors) for items in lists(record))
    assert not {"group_labels", "sex", "metadata", *check.GROUPS} & {leaf for leaf in leaves if isinstance(leaf, str)}
    assert set(record["arms"]) == {"adjusted", "unadjusted"}
    assert set(record["seeds"]) == {"labels", "plant", "rrpp", "bootstrap"} and record["label_redraws"] == 0
    for name in check.BLOCKS:
        assert len(record["planted"][name]["features"]) == 10 and len(record["planted"][name]["signs"]) == 10
    assert -1.0 <= record["plant_cosine"] <= 1.0
    analysis = record["analysis"]
    assert set(analysis["p_values"]) == {"delta", "angle", "shape"} and analysis["selected_lv"] >= 2
    assert isinstance(analysis["degenerate"], bool) and "cv_mean_auroc" in analysis
    principal = record["arms"]["adjusted"]["blocks"]["principal"]["observed"]
    assert principal["valid_replicates"] == 5 and principal["planted_size"] == 20
    assert set(record["arms"]["adjusted"]["blocks"]) == {"0->1", "1->2", "principal"}


# ── Run, resume and merge ────────────────────────────────────────────────────


def _records(out_dir: Path) -> list[dict]:
    lines = [
        line
        for path in sorted((out_dir / "records").glob("shard_*.jsonl"))
        for line in path.read_text(encoding="utf-8").splitlines()
    ]
    return sorted((json.loads(line) for line in lines), key=lambda record: record["key"])


@pytest.fixture(scope="module")
def sharded_run(tmp_path_factory, export, stand_in_fidelity):
    out_dir = tmp_path_factory.mktemp("plasmode-run")
    settings = small_settings(export, stand_in_fidelity, descriptive=0.5, e_star=1.0)
    check.ensure_settings(out_dir, settings)
    evaluated = [check.run_shard(out_dir, settings, export, index, 2) for index in (0, 1)]
    return out_dir, settings, evaluated


def test_two_shards_cover_every_unit_once(sharded_run) -> None:
    out_dir, settings, evaluated = sharded_run
    assert evaluated == [3, 3] and settings["cells"] == [0.0, 1.0, 0.5]
    keys = [record["key"] for record in _records(out_dir)]
    assert sorted(keys) == sorted(check.unit_key(*unit) for unit in check.units(settings))
    assert len(keys) == len(set(keys))


def test_resume_on_a_complete_shard_appends_nothing(sharded_run, export) -> None:
    out_dir, settings, _ = sharded_run
    before = check.shard_path(out_dir, 0, 2).read_bytes()
    assert check.run_shard(out_dir, settings, export, 0, 2) == 0
    assert check.shard_path(out_dir, 0, 2).read_bytes() == before


def test_resume_under_a_changed_threshold_is_refused(sharded_run, export, stand_in_fidelity, monkeypatch) -> None:
    out_dir, _, _ = sharded_run
    before = {path.name: path.read_bytes() for path in (out_dir / "records").iterdir()}
    monkeypatch.setattr(check, "THRESHOLDS", {**check.THRESHOLDS, "sign_stability": 0.85})
    with pytest.raises(SystemExit, match="different settings"):
        check.ensure_settings(out_dir, small_settings(export, stand_in_fidelity, descriptive=0.5, e_star=1.0))
    assert {path.name: path.read_bytes() for path in (out_dir / "records").iterdir()} == before


def test_resume_on_a_changed_export_is_refused(export_dir, tmp_path) -> None:
    directory = _copy_export(export_dir, tmp_path / "export")
    export = check.load_export(directory)
    check.run_fidelity(export, tmp_path, analysis=_fake_analysis(5, 0.635556))
    argv = ["run", "--export-dir", str(directory), "--fidelity", str(tmp_path / "fidelity.json"),
            "--e-star", "1.0", "--reps", "2", "--out-dir", str(tmp_path / "run"), "--limit", "0"]
    assert check.main(argv) == 0
    shard = check.shard_path(tmp_path / "run", 0, 1)
    assert shard.read_text() == ""
    text = (directory / "rna.csv").read_text().splitlines()
    text[1] = text[1].rsplit(",", 1)[0] + ",9.5"
    (directory / "rna.csv").write_text("\n".join(text) + "\n")
    with pytest.raises(SystemExit, match="another export"):
        check.main(argv)
    stale = check.load_export(directory)
    assert stale.files["rna.csv"] != export.files["rna.csv"]
    settings = check.load_settings(tmp_path / "run")
    with pytest.raises(SystemExit, match="differ from the settings"):
        check.run_shard(tmp_path / "run", settings, stale, 0, 1)
    assert shard.read_text() == ""


def test_merge_writes_the_summary_and_is_byte_identical(sharded_run) -> None:
    out_dir, settings, _ = sharded_run
    first = check.merge(out_dir, settings)
    names = ("plasmode_summary.csv", "merge.json", "plasmode_units.csv", "plasmode_metrics.csv", "merged.jsonl")
    snapshot = {name: (out_dir / name).read_bytes() for name in names}
    second = check.merge(out_dir, settings)
    assert {name: (out_dir / name).read_bytes() for name in names} == snapshot and first == second
    assert first["units"] == first["expected"] == 6 and first["missing"] == []
    assert first["units_by_cell"] == {"anchor": 2, "e=1": 2, "e=0.5": 2}
    summary = pd.read_csv(out_dir / "plasmode_summary.csv", keep_default_na=False)
    reference = summary[summary["role"] == "false_recovery_reference"]
    assert len(reference) == 1
    row = reference.iloc[0]
    assert (row["cell"], row["arm"], row["block"], row["component"], row["subset"]) == (
        "anchor", "adjusted", "principal", "observed", "all"
    )
    metrics = pd.read_csv(out_dir / "plasmode_metrics.csv")
    anchor = metrics[(metrics["effect"] == 0.0) & (metrics["arm"] == "adjusted") & (metrics["block"] == "principal")
                     & (metrics["component"] == "observed")]
    assert float(row["mean"]) == pytest.approx(anchor["planted_precision"].mean())
    assert first["false_recovery_reference"]["mean"] == pytest.approx(anchor["planted_precision"].mean())
    assert set(summary["subset"]) == {"all", "angle_rejects"}
    assert len(summary[summary["role"] == "criterion_i"]) == 3
    assert len(summary[summary["role"] == "criterion_iii"]) == 3


def test_smoke_decision_is_indicative_and_reports_the_descriptive_readings(sharded_run) -> None:
    out_dir, settings, _ = sharded_run
    check.merge(out_dir, settings)
    decision = check.decide(out_dir, settings)
    assert decision["indicative"] is True and decision["verdict"] in check.VERDICTS
    assert any("permutations" in reason for reason in decision["indicative_reasons"])
    assert any("explicitly" in reason for reason in decision["indicative_reasons"])
    assert decision["replicates"] == {"anchor": 2, "e=1": 2, "e=0.5": 2}
    assert decision["settings_hash"] == check.settings_hash(settings)
    descriptive = decision["descriptive"]
    assert set(descriptive["transitions"]) == {"0->1", "1->2"}
    assert descriptive["descriptive_cell"]["cell"] == "e=0.5"
    assert descriptive["false_recovery_reference"]["n"] == 2
    assert json.loads((out_dir / "plasmode_decision.json").read_text())["verdict"] == decision["verdict"]


# ── Calibration ──────────────────────────────────────────────────────────────


def _calibrate(tmp_path, export, rates: dict[float, float]):
    tried = []

    def run_once(settings, export, workers=1):
        effect = settings["cells"][0]
        tried.append(effect)
        assert settings["stream"] == "calibration" and settings["attribution"] is None
        assert settings["replicates"] == 40 and settings["permutations"] == 99
        return {"rate": rates[effect], "rates": {}, "units": 40}

    settings = check.calibration_settings(export, {"selected_lv": 5, "cv_mean_auroc": 0.635556})
    return check.calibrate(tmp_path, settings, export, run_once=run_once), tried


def test_first_qualifying_step_is_adopted(tmp_path, export) -> None:
    entry, tried = _calibrate(tmp_path, export, {0.5: 0.1, 1.0: 0.4, 1.5: 0.80, 2.0: 0.925, 3.0: 1.0})
    assert tried == [0.5, 1.0, 1.5, 2.0]
    assert (entry["outcome"], entry["e_star"], entry["descriptive"]) == ("adopt", 2.0, 1.5)
    history = json.loads((tmp_path / "calibration.json").read_text())
    assert [item["effect"] for item in history] == [0.5, 1.0, 1.5, 2.0]
    assert [item["outcome"] for item in history] == ["climb", "climb", "climb", "adopt"]
    assert [item["rate"] for item in history] == [0.1, 0.4, 0.80, 0.925]


def test_cap_reached_proceeds_at_the_cap(tmp_path, export) -> None:
    entry, tried = _calibrate(tmp_path, export, dict.fromkeys(check.LADDER, 0.5))
    assert tried == list(check.LADDER)
    assert (entry["outcome"], entry["e_star"], entry["descriptive"]) == ("cap_reached", 3.0, 2.0)
    assert json.loads((tmp_path / "calibration.json").read_text())[-1]["outcome"] == "cap_reached"


def test_lowest_step_qualifying_halves_for_the_descriptive_cell(tmp_path, export) -> None:
    entry, tried = _calibrate(tmp_path, export, dict.fromkeys(check.LADDER, 0.95))
    assert tried == [0.5] and (entry["e_star"], entry["descriptive"]) == (0.5, 0.25)


def test_calibration_resumes_after_its_last_iteration(tmp_path, export) -> None:
    rates = {0.5: 0.1, 1.0: 0.95}

    def interrupted(settings, export, workers=1):
        if settings["cells"][0] == 1.0:
            raise KeyboardInterrupt
        return {"rate": rates[settings["cells"][0]]}

    settings = check.calibration_settings(export, {"selected_lv": 5, "cv_mean_auroc": 0.635556})
    with pytest.raises(KeyboardInterrupt):
        check.calibrate(tmp_path, settings, export, run_once=interrupted)
    entry, tried = _calibrate(tmp_path, export, rates)
    assert tried == [1.0] and entry["e_star"] == 1.0
    again, tried = _calibrate(tmp_path, export, rates)
    assert tried == [] and again["iteration"] == 1
    assert len(json.loads((tmp_path / "calibration.json").read_text())) == 2


def test_settings_name_the_adopted_iteration(tmp_path, export) -> None:
    _calibrate(tmp_path, export, {0.5: 0.1, 1.0: 0.4, 1.5: 0.80, 2.0: 0.925})
    path = tmp_path / "calibration.json"
    settings = check.build_settings(export, {"selected_lv": 5, "cv_mean_auroc": 0.635556}, calibration=path)
    assert settings["adopted_calibration"]["iteration"] == 3 and settings["adopted_calibration"]["e_star"] == 2.0
    assert settings["cells"] == [0.0, 2.0, 1.5] and settings["strength"]["source"] == "calibration"
    other = check.load_export(write_export(tmp_path / "other", seed=9))
    with pytest.raises(SystemExit, match="another export"):
        check.build_settings(other, {"selected_lv": 5, "cv_mean_auroc": 0.635556}, calibration=path)


# ── Decision ─────────────────────────────────────────────────────────────────


def _decision_settings(export, *, e_star: float = 2.0, descriptive: float = 1.5, **scale) -> dict:
    settings = check.build_settings(
        export, {"selected_lv": 5, "cv_mean_auroc": 0.635556}, e_star=e_star, descriptive=descriptive, **scale
    )
    settings["strength"]["source"] = "calibration"
    return settings


def _summary(settings, overrides: dict | None = None, *, n: int = 500) -> list[dict]:
    """A summary where (i), (ii) and (iii) hold; ``overrides`` maps a row key to (mean, se)."""

    e_star, descriptive = settings["strength"]["e_star"], settings["strength"]["descriptive"]
    values: dict[tuple, tuple[float, float | None]] = {}
    for effect in (0.0, e_star, descriptive):
        values[(effect, "analysis", "", "", "all", "selected_lv")] = (5.0, 0.01)
    for statistic, rate in zip(check.STATISTICS, (0.046, 0.058, 0.040)):
        values[(0.0, "analysis", "", "", "all", f"reject_{statistic}")] = (rate, 0.01)
    values[(e_star, "analysis", "", "", "all", "reject_angle")] = (0.92, 0.012)
    for quantity, (mean, se) in {"planted_precision": (0.71, 0.01), "sign_stability": (0.95, 0.005),
                                 "top_k_frequency": (0.62, 0.01)}.items():
        values[(e_star, "adjusted", "principal", "observed", "angle_rejects", quantity)] = (mean, se)
    values.update(overrides or {})
    return [
        {"effect": key[0], "cell": check.cell_label(key[0]), "arm": key[1], "block": key[2], "component": key[3],
         "subset": key[4], "quantity": key[5], "n": n, "mean": mean, "se": se, "role": ""}
        for key, (mean, se) in values.items()
    ]


def _criterion(decision: dict, criterion: str, quantity: str) -> dict:
    return next(entry for entry in decision["criteria"][criterion] if entry["quantity"] == quantity)


RECOVERY = (2.0, "adjusted", "principal", "observed", "angle_rejects")


def test_size_holds_at_the_anchor(export) -> None:
    settings = _decision_settings(export)
    decision = check.plasmode_decision(_summary(settings), settings)
    for statistic in check.STATISTICS:
        entry = _criterion(decision, "size", statistic)
        assert entry["verdict"] == "holds" and entry["bound"] == pytest.approx(0.0695, abs=5e-5) and entry["n"] == 500


def test_size_fails_on_shape(export) -> None:
    settings = _decision_settings(export)
    decision = check.plasmode_decision(
        _summary(settings, {(0.0, "analysis", "", "", "all", "reject_shape"): (0.072, 0.012)}), settings
    )
    assert decision["verdict"] == "attribution_leg_failed"
    assert [(item["criterion"], item["quantity"]) for item in decision["failing"]] == [("size", "shape")]


def test_detection_fails_at_the_cap(export) -> None:
    settings = _decision_settings(export, e_star=3.0, descriptive=2.0)
    overrides = {(3.0, "analysis", "", "", "all", "reject_angle"): (0.62, 0.02)}
    decision = check.plasmode_decision(_summary(settings, overrides), settings)
    assert decision["failing"] == [
        {"criterion": "detection", "cell": "e=3", "quantity": "angle", "value": 0.62, "bound": 0.8, "verdict": "fails"}
    ]


def test_recovery_holds_and_the_leg_is_lifted(export) -> None:
    settings = _decision_settings(export)
    decision = check.plasmode_decision(_summary(settings), settings)
    assert decision["verdict"] == "attribution_leg_lifted" and decision["failing"] == []
    assert decision["indicative"] is False and decision["indicative_reasons"] == []
    precision = _criterion(decision, "recovery", "planted_precision")
    assert precision["value"] == pytest.approx(0.69) and precision["bound"] == 0.5 and precision["se"] == 0.01
    assert _criterion(decision, "recovery", "sign_stability")["value"] == pytest.approx(0.94)
    assert decision["settings_hash"] == check.settings_hash(settings)


def test_stability_fails_and_names_top_k_frequency(export) -> None:
    settings = _decision_settings(export)
    decision = check.plasmode_decision(_summary(settings, {RECOVERY + ("top_k_frequency",): (0.51, 0.01)}), settings)
    assert decision["verdict"] == "attribution_leg_failed"
    (failing,) = decision["failing"]
    assert (failing["criterion"], failing["quantity"]) == ("recovery", "top_k_frequency")
    assert failing["value"] == pytest.approx(0.49) and failing["bound"] == 0.5


def test_failed_on_two_criteria_names_both(export) -> None:
    settings = _decision_settings(export)
    decision = check.plasmode_decision(
        _summary(settings, {(2.0, "analysis", "", "", "all", "reject_angle"): (0.62, 0.02),
                            RECOVERY + ("sign_stability",): (0.89, 0.005)}),
        settings,
    )
    assert decision["verdict"] == "attribution_leg_failed"
    assert {(item["criterion"], item["quantity"]) for item in decision["failing"]} == {
        ("detection", "angle"), ("recovery", "sign_stability")
    }


def test_a_missing_reading_fails_its_criterion(export) -> None:
    settings = _decision_settings(export)
    summary = [row for row in _summary(settings) if row["quantity"] != "sign_stability"]
    decision = check.plasmode_decision(summary, settings)
    assert decision["failing"] == [
        {"criterion": "recovery", "cell": "e=2", "quantity": "sign_stability", "value": None, "bound": 0.9,
         "verdict": "unavailable"}
    ]


def test_smoke_decision_is_indicative(export) -> None:
    settings = _decision_settings(export, replicates=20, permutations=99, bootstrap=50)
    decision = check.plasmode_decision(_summary(settings, n=20), settings)
    assert decision["indicative"] is True and len(decision["indicative_reasons"]) == 3
    settings["strength"]["source"] = "explicit"
    assert len(check.plasmode_decision(_summary(settings, n=20), settings)["indicative_reasons"]) == 4


def test_descriptive_readings_do_not_move_the_verdict(export) -> None:
    settings = _decision_settings(export)
    poor = (0.01, 0.01)
    overrides = {}
    for quantity in check.THRESHOLDS:
        overrides[(2.0, "unadjusted", "principal", "observed", "angle_rejects", quantity)] = poor
        overrides[(2.0, "adjusted", "principal", "pls_captured", "angle_rejects", quantity)] = poor
        overrides[(2.0, "adjusted", "principal", "residual", "angle_rejects", quantity)] = poor
        overrides[(2.0, "adjusted", "0->1", "observed", "angle_rejects", quantity)] = poor
        overrides[(2.0, "adjusted", "principal", "observed", "all", quantity)] = poor
        overrides[(1.5, "adjusted", "principal", "observed", "angle_rejects", quantity)] = poor
    overrides[(2.0, "analysis", "", "", "all", "reject_delta")] = (1.0, 0.0)
    overrides[(2.0, "analysis", "", "", "all", "reject_shape")] = (0.0, 0.0)
    overrides[(1.5, "analysis", "", "", "all", "reject_angle")] = (0.1, 0.01)
    overrides[(0.0, "adjusted", "principal", "observed", "all", "planted_precision")] = (0.9, 0.01)
    decision = check.plasmode_decision(_summary(settings, overrides), settings)
    assert decision["verdict"] == "attribution_leg_lifted" and decision["failing"] == []
    descriptive = decision["descriptive"]
    assert all(entry["verdict"] == "fails" for entry in descriptive["unadjusted_arm"])
    assert all(entry["verdict"] == "fails" for entry in descriptive["components"]["residual"])
    assert all(entry["verdict"] == "fails" for entry in descriptive["transitions"]["0->1"])
    assert all(entry["verdict"] == "fails" for entry in descriptive["unconditional"])
    assert descriptive["e_star_rates"]["delta"]["rate"] == 1.0 and descriptive["e_star_rates"]["shape"]["rate"] == 0.0
    assert descriptive["descriptive_cell"]["rates"]["angle"]["rate"] == 0.1
    assert descriptive["false_recovery_reference"]["mean"] == 0.9
