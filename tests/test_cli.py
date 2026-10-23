# tests/test_cli.py
from __future__ import annotations

import json
from pathlib import Path
from unittest.mock import patch

import numpy as np
import pandas as pd
import pytest

from motco.cli import main
from motco.simulations.reference import ReferenceCacheMissingError
from motco.simulations.semisynthetic import SemiSyntheticTrajectoryDataset
from motco.stats.design import build_ls_means, get_model_matrix

# ── Fixtures ──────────────────────────────────────────────────────────────────

@pytest.fixture()
def plsr_csv(tmp_path: Path) -> Path:
    rng = np.random.default_rng(0)
    n, p = 30, 5
    df = pd.DataFrame(rng.standard_normal((n, p)), columns=[f"f{i}" for i in range(p)])
    df["label"] = ["A"] * 15 + ["B"] * 15
    path = tmp_path / "plsr_data.csv"
    df.to_csv(path, index=False)
    return path


@pytest.fixture()
def snf_csvs(tmp_path: Path) -> list[Path]:
    rng = np.random.default_rng(0)
    paths = []
    for i in range(2):
        df = pd.DataFrame(rng.standard_normal((15, 4)))
        p = tmp_path / f"omics{i}.csv"
        df.to_csv(p, index=False)
        paths.append(p)
    return paths


@pytest.fixture()
def de_files(tmp_path: Path) -> dict[str, Path]:
    rng = np.random.default_rng(0)
    n = 20
    Y = pd.DataFrame(rng.standard_normal((n, 3)))
    factors = pd.DataFrame({
        "group": ["A"] * 10 + ["B"] * 10,
        "level": ["t0", "t1"] * 10,
    })
    M_full = get_model_matrix(factors, group_col="group", level_col="level", full=True)
    M_red = get_model_matrix(factors, group_col="group", level_col="level", full=False)
    LS = build_ls_means(["A", "B"], ["t0", "t1"], full=True)
    contrast = [[0, 1], [2, 3]]

    paths = {}
    paths["Y"] = tmp_path / "Y.csv"
    paths["model"] = tmp_path / "model.csv"
    paths["model_full"] = tmp_path / "model_full.csv"
    paths["model_red"] = tmp_path / "model_red.csv"
    paths["ls"] = tmp_path / "ls.csv"
    paths["contrast"] = tmp_path / "contrast.json"

    Y.to_csv(paths["Y"], index=False)
    pd.DataFrame(M_full).to_csv(paths["model"], index=False)
    pd.DataFrame(M_full).to_csv(paths["model_full"], index=False)
    pd.DataFrame(M_red).to_csv(paths["model_red"], index=False)
    pd.DataFrame(LS).to_csv(paths["ls"], index=False)
    paths["contrast"].write_text(json.dumps(contrast))
    return paths


# ── plsr subcommand ───────────────────────────────────────────────────────────

def test_plsr_saves_csv(tmp_path: Path, plsr_csv: Path) -> None:
    out = tmp_path / "table.csv"
    main([
        "plsr", "--data", str(plsr_csv), "--label-col", "label",
        "--cv1-splits", "3", "--cv2-splits", "3", "--n-repeats", "2",
        "--max-components", "3", "--out-table", str(out),
    ])
    assert out.exists()
    df = pd.read_csv(out)
    assert df.shape[0] == 2  # n_repeats rows


def test_plsr_prints_to_stdout(capsys: pytest.CaptureFixture, plsr_csv: Path) -> None:
    main([
        "plsr", "--data", str(plsr_csv), "--label-col", "label",
        "--cv1-splits", "3", "--cv2-splits", "3", "--n-repeats", "2",
        "--max-components", "3",
    ])
    out = capsys.readouterr().out
    assert "AUROC" in out


def test_plsr_bad_label_col_exits(plsr_csv: Path) -> None:
    with pytest.raises(SystemExit):
        main([
            "plsr", "--data", str(plsr_csv), "--label-col", "nonexistent",
            "--cv1-splits", "3", "--cv2-splits", "3", "--n-repeats", "2",
            "--max-components", "3",
        ])


# ── snf subcommand ────────────────────────────────────────────────────────────

def test_snf_saves_fused_matrix(tmp_path: Path, snf_csvs: list[Path]) -> None:
    out = tmp_path / "fused.csv"
    main([
        "snf",
        "--input", str(snf_csvs[0]), "--input", str(snf_csvs[1]),
        "--K", "5", "--k", "5", "--t", "3", "--out-fused", str(out),
    ])
    assert out.exists()
    df = pd.read_csv(out)
    assert df.shape == (15, 15)


def test_snf_saves_embedding(tmp_path: Path, snf_csvs: list[Path]) -> None:
    out_fused = tmp_path / "fused.csv"
    out_emb = tmp_path / "emb.csv"
    main([
        "snf",
        "--input", str(snf_csvs[0]), "--input", str(snf_csvs[1]),
        "--K", "5", "--k", "5", "--t", "3",
        "--out-fused", str(out_fused), "--out-embedding", str(out_emb),
    ])
    assert out_emb.exists()
    emb = pd.read_csv(out_emb)
    assert emb.shape == (15, 10)


def test_snf_requires_two_inputs() -> None:
    with pytest.raises(SystemExit):
        main(["snf", "--input", "only_one.csv"])


# ── de subcommand — estimate path ─────────────────────────────────────────────

def test_de_estimate_saves_json(tmp_path: Path, de_files: dict[str, Path]) -> None:
    out = tmp_path / "result.json"
    main([
        "de",
        "--Y", str(de_files["Y"]),
        "--model-matrix", str(de_files["model"]),
        "--ls-means", str(de_files["ls"]),
        "--contrast", str(de_files["contrast"]),
        "--out-json", str(out),
    ])
    assert out.exists()
    result = json.loads(out.read_text())
    assert {"deltas", "angles", "shapes"} <= result.keys()
    # 2 groups → 2×2 matrices
    assert len(result["deltas"]) == 2
    assert len(result["deltas"][0]) == 2


# ── de subcommand — RRPP path ─────────────────────────────────────────────────

def test_de_rrpp_saves_json(tmp_path: Path, de_files: dict[str, Path]) -> None:
    out = tmp_path / "rrpp.json"
    main([
        "de",
        "--Y", str(de_files["Y"]),
        "--model-full", str(de_files["model_full"]),
        "--model-reduced", str(de_files["model_red"]),
        "--ls-means", str(de_files["ls"]),
        "--contrast", str(de_files["contrast"]),
        "--rrpp-permutations", "2",
        "--out-json", str(out),
    ])
    assert out.exists()
    result = json.loads(out.read_text())
    assert {"deltas", "angles", "shapes"} <= result.keys()
    assert len(result["deltas"]) == 2  # 2 permutations


# ── --version flag ────────────────────────────────────────────────────────────

def test_version_flag_prints_version(capsys: pytest.CaptureFixture) -> None:
    with pytest.raises(SystemExit) as exc:
        main(["--version"])
    assert exc.value.code == 0
    # argparse prints version to stdout
    assert "motco" in capsys.readouterr().out


# ── --out-observed flag ───────────────────────────────────────────────────────

def test_de_out_observed_saves_csv(tmp_path: Path, de_files: dict[str, Path]) -> None:
    out_json = tmp_path / "result.json"
    out_obs = tmp_path / "observed.csv"
    main([
        "de",
        "--Y", str(de_files["Y"]),
        "--model-matrix", str(de_files["model"]),
        "--ls-means", str(de_files["ls"]),
        "--contrast", str(de_files["contrast"]),
        "--out-json", str(out_json),
        "--out-observed", str(out_obs),
    ])
    assert out_obs.exists()
    obs = pd.read_csv(out_obs)
    # LS has 4 rows (2 groups × 2 levels), Y has 3 features
    assert obs.shape == (4, 3)


@pytest.fixture()
def multi_omics_csvs(tmp_path: Path) -> tuple[list[Path], Path]:
    """Two small omics CSVs and a metadata CSV with a label column."""
    rng = np.random.default_rng(42)
    n = 30
    paths = []
    for i, name in enumerate(["methyl", "expr"]):
        df = pd.DataFrame(rng.standard_normal((n, 6)), columns=[f"{name}_{j}" for j in range(6)])
        p = tmp_path / f"{name}.csv"
        df.to_csv(p, index=False)
        paths.append(p)
    meta = pd.DataFrame({"sample_id": range(n), "stage": [0, 1, 2] * 10})
    meta_path = tmp_path / "metadata.csv"
    meta.to_csv(meta_path, index=False)
    return paths, meta_path


def test_plsr_multi_input_runs(tmp_path: Path, multi_omics_csvs: tuple) -> None:
    input_paths, meta_path = multi_omics_csvs
    out = tmp_path / "table.csv"
    main([
        "plsr",
        "--input", str(input_paths[0]), "--input", str(input_paths[1]),
        "--metadata", str(meta_path), "--label-col", "stage",
        "--cv1-splits", "3", "--cv2-splits", "3", "--n-repeats", "2",
        "--max-components", "3", "--out-table", str(out),
    ])
    assert out.exists()
    df = pd.read_csv(out)
    assert df.shape[0] == 2


def test_plsr_out_scores_correct_shape(tmp_path: Path, plsr_csv: Path) -> None:
    out_table = tmp_path / "table.csv"
    out_scores = tmp_path / "scores.csv"
    main([
        "plsr", "--data", str(plsr_csv), "--label-col", "label",
        "--cv1-splits", "3", "--cv2-splits", "3", "--n-repeats", "2",
        "--max-components", "3",
        "--out-table", str(out_table),
        "--out-scores", str(out_scores),
    ])
    assert out_scores.exists()
    scores = pd.read_csv(out_scores)
    assert scores.shape[0] == 30  # n_samples from plsr_csv fixture
    assert scores.shape[1] >= 1
    assert all(c.startswith("lv_") for c in scores.columns)


@pytest.fixture()
def plsr_wide_csv(tmp_path: Path) -> Path:
    rng = np.random.default_rng(3)
    n, p = 36, 8
    labels = np.repeat(["A", "B", "C"], n // 3)
    shift = np.repeat([0.0, 0.8, 1.6], n // 3)[:, None]
    df = pd.DataFrame(rng.standard_normal((n, p)) + shift * np.linspace(0, 1, p), columns=[f"f{i}" for i in range(p)])
    df["label"] = labels
    path = tmp_path / "plsr_wide.csv"
    df.to_csv(path, index=False)
    return path


def _plsr_floor_args(csv: Path, out_table: Path, *extra: str) -> list[str]:
    return [
        "plsr", "--data", str(csv), "--label-col", "label",
        "--cv1-splits", "2", "--cv2-splits", "3", "--n-repeats", "2",
        "--out-table", str(out_table), *extra,
    ]


def test_plsr_min_components_floors_the_table_and_the_scores(tmp_path: Path, plsr_wide_csv: Path) -> None:
    out_table = tmp_path / "table.csv"
    out_scores = tmp_path / "scores.csv"
    main(_plsr_floor_args(plsr_wide_csv, out_table, "--min-components", "2", "--max-components", "6",
                          "--out-scores", str(out_scores)))
    table = pd.read_csv(out_table)
    assert (table["LV"] >= 2).all()
    assert pd.read_csv(out_scores).shape[1] >= 2


def test_plsr_infeasible_min_components_exits(tmp_path: Path, plsr_wide_csv: Path,
                                              capsys: pytest.CaptureFixture) -> None:
    with pytest.raises(SystemExit) as exc_info:
        main(_plsr_floor_args(plsr_wide_csv, tmp_path / "t.csv", "--min-components", "6", "--max-components", "6"))
    assert exc_info.value.code != 0
    message = str(exc_info.value)
    assert "min_components=6" in message and "max_components=6" in message


def test_plsr_default_min_components_is_unchanged(tmp_path: Path, plsr_wide_csv: Path) -> None:
    without = tmp_path / "without.csv"
    explicit = tmp_path / "explicit.csv"
    main(_plsr_floor_args(plsr_wide_csv, without, "--max-components", "5"))
    main(_plsr_floor_args(plsr_wide_csv, explicit, "--max-components", "5", "--min-components", "1"))
    pd.testing.assert_frame_equal(pd.read_csv(without), pd.read_csv(explicit))


def test_plsr_n_components_override(tmp_path: Path, plsr_csv: Path) -> None:
    out_scores = tmp_path / "scores.csv"
    main([
        "plsr", "--data", str(plsr_csv), "--label-col", "label",
        "--cv1-splits", "3", "--cv2-splits", "3", "--n-repeats", "2",
        "--max-components", "3",
        "--out-scores", str(out_scores),
        "--n-components", "2",
    ])
    assert out_scores.exists()
    scores = pd.read_csv(out_scores)
    assert scores.shape == (30, 2)


def test_plsr_input_and_data_exclusive(tmp_path: Path, plsr_csv: Path, multi_omics_csvs: tuple) -> None:
    input_paths, meta_path = multi_omics_csvs
    with pytest.raises(SystemExit):
        main([
            "plsr",
            "--input", str(input_paths[0]),
            "--data", str(plsr_csv),
            "--label-col", "stage",
            "--metadata", str(meta_path),
            "--cv1-splits", "3", "--cv2-splits", "3", "--n-repeats", "2",
            "--max-components", "3",
        ])


def test_plsr_input_without_metadata_exits(tmp_path: Path, multi_omics_csvs: tuple) -> None:
    input_paths, _ = multi_omics_csvs
    with pytest.raises(SystemExit):
        main([
            "plsr",
            "--input", str(input_paths[0]),
            "--label-col", "stage",
            "--cv1-splits", "3", "--cv2-splits", "3", "--n-repeats", "2",
            "--max-components", "3",
        ])


def test_plsr_missing_label_col_in_metadata_exits(tmp_path: Path, multi_omics_csvs: tuple) -> None:
    input_paths, meta_path = multi_omics_csvs
    with pytest.raises(SystemExit):
        main([
            "plsr",
            "--input", str(input_paths[0]),
            "--metadata", str(meta_path),
            "--label-col", "nonexistent_col",
            "--cv1-splits", "3", "--cv2-splits", "3", "--n-repeats", "2",
            "--max-components", "3",
        ])


def test_plsr_out_vips_saves_csv(tmp_path: Path, plsr_csv: Path) -> None:
    out_table = tmp_path / "table.csv"
    out_vips = tmp_path / "vips.csv"
    main([
        "plsr", "--data", str(plsr_csv), "--label-col", "label",
        "--cv1-splits", "3", "--cv2-splits", "3", "--n-repeats", "2",
        "--max-components", "3",
        "--out-table", str(out_table),
        "--out-vips", str(out_vips),
    ])
    assert out_vips.exists()
    vips_df = pd.read_csv(out_vips)
    # 5 features (from plsr_csv fixture), 2 repeats
    assert vips_df.shape == (5, 2)
    assert list(vips_df.columns) == ["rep_1", "rep_2"]


def test_snf_custom_spectral_components(tmp_path: Path, snf_csvs: list[Path]) -> None:
    out_fused = tmp_path / "fused.csv"
    out_emb = tmp_path / "emb.csv"
    main([
        "snf",
        "--input", str(snf_csvs[0]), "--input", str(snf_csvs[1]),
        "--K", "5", "--k", "5", "--t", "3",
        "--out-fused", str(out_fused),
        "--out-embedding", str(out_emb),
        "--spectral-components", "5",
    ])
    assert out_emb.exists()
    emb = pd.read_csv(out_emb)
    assert emb.shape == (15, 5)


# ── simulate subcommand helpers ───────────────────────────────────────────────

_SIMULATE_EXPECTED_FILES = {
    "methylation.csv", "expression.csv", "proteomics.csv",
    "metadata.csv", "model_full.csv", "model_reduced.csv",
    "ls_means.csv", "contrast.json", "truth.json",
}


def _make_fake_simulate_dataset(n: int = 30, seed: int = 0) -> SemiSyntheticTrajectoryDataset:
    """Minimal SemiSyntheticTrajectoryDataset for testing cmd_simulate."""
    rng = np.random.default_rng(seed)
    sample_ids = [f"S{i:03d}" for i in range(n)]
    # Interleave groups so all (group, stage) cells are populated
    groups = (["A", "B"] * n)[:n]
    stages = ([0] * (n // 3) + [1] * (n // 3) + [2] * (n // 3))[:n]
    metadata = pd.DataFrame({
        "sample_id": sample_ids,
        "group": groups,
        "stage": stages,
        "cluster": stages,
    })
    methyl = pd.DataFrame(
        rng.standard_normal((n, 5)), index=sample_ids,
        columns=[f"cg{i}" for i in range(5)],
    )
    expr = pd.DataFrame(
        rng.standard_normal((n, 4)), index=sample_ids,
        columns=[f"gene{i}" for i in range(4)],
    )
    prot = pd.DataFrame(
        rng.standard_normal((n, 3)), index=sample_ids,
        columns=[f"prot{i}" for i in range(3)],
    )
    truth = {
        "trajectory_mode": "orientation",
        "group_effect_size": 1.0,
        "group_labels": ["A", "B"],
        "group_ratio": 0.5,
        "n_stages": 3,
        "p_dmp": 0.2,
        "seed": seed,
        "stage_assumption": "clusters-as-stages",
        "deltas": {"A": [2.0, 2.0, 2.0], "B": [2.0, 2.0, 2.0]},
        "indicator_counts": {"A": {"methylation": [1, 2, 1]}, "B": {"methylation": [1, 1, 2]}},
        "location_offset": None,
    }
    return SemiSyntheticTrajectoryDataset(
        methylation=methyl, expression=expr, proteomics=prot,
        metadata=metadata, truth=truth,
    )


@pytest.fixture()
def mock_simulate_env():
    """Patch the reference loader and dataset generator for simulate tests."""
    fake_dataset = _make_fake_simulate_dataset(n=30, seed=0)
    with patch("motco.simulations.reference.load_reference", return_value=object()), \
         patch("motco.simulations.semisynthetic.generate_semisynthetic_trajectory",
               return_value=fake_dataset):
        yield fake_dataset


# ── simulate subcommand tests ─────────────────────────────────────────────────

def test_simulate_produces_all_files(tmp_path: Path, mock_simulate_env: SemiSyntheticTrajectoryDataset) -> None:
    main(["simulate", "--seed", "0", "--out-dir", str(tmp_path)])
    produced = {f.name for f in tmp_path.iterdir()}
    assert _SIMULATE_EXPECTED_FILES <= produced


def test_simulate_reproducible(tmp_path: Path, mock_simulate_env: SemiSyntheticTrajectoryDataset) -> None:
    out1 = tmp_path / "run1"
    out2 = tmp_path / "run2"
    main(["simulate", "--seed", "42", "--out-dir", str(out1)])
    main(["simulate", "--seed", "42", "--out-dir", str(out2)])
    for fname in _SIMULATE_EXPECTED_FILES:
        assert (out1 / fname).read_bytes() == (out2 / fname).read_bytes()


def test_simulate_missing_reference_cache_exits(tmp_path: Path) -> None:
    def _raise(*args, **kwargs):
        raise ReferenceCacheMissingError("reference cache not found at intersim_reference.npz")

    with patch("motco.simulations.reference.load_reference", side_effect=_raise):
        with pytest.raises(SystemExit) as exc:
            main(["simulate", "--seed", "0", "--out-dir", str(tmp_path)])
    assert "intersim_reference.npz" in str(exc.value)


def test_simulate_wires_trajectory_params(tmp_path: Path) -> None:
    fake_dataset = _make_fake_simulate_dataset(n=30, seed=0)
    with patch("motco.simulations.reference.load_reference", return_value=object()), \
         patch("motco.simulations.semisynthetic.generate_semisynthetic_trajectory",
               return_value=fake_dataset) as generate_mock:
        main([
            "simulate",
            "--seed", "0",
            "--out-dir", str(tmp_path),
            "--trajectory-mode", "magnitude",
            "--n-samples", "300",
            "--effect-size", "0.5",
            "--magnitude-kind", "extremes",
        ])
    (traj_params,) = generate_mock.call_args.args
    assert traj_params.trajectory_mode == "magnitude"
    assert traj_params.n_samples == 300
    assert traj_params.group_effect_size == 0.5
    assert traj_params.magnitude_kind == "extremes"


def test_simulate_magnitude_kind_defaults_to_joint(tmp_path: Path) -> None:
    fake_dataset = _make_fake_simulate_dataset(n=30, seed=0)
    with patch("motco.simulations.reference.load_reference", return_value=object()), \
         patch("motco.simulations.semisynthetic.generate_semisynthetic_trajectory",
               return_value=fake_dataset) as generate_mock:
        main(["simulate", "--seed", "0", "--out-dir", str(tmp_path)])
    (traj_params,) = generate_mock.call_args.args
    assert traj_params.magnitude_kind == "joint"


@pytest.mark.parametrize(("flag", "expected"), [([], "joint"), (["--magnitude-kind", "all"], "all")])
def test_simulate_truth_records_magnitude_kind(tmp_path: Path, flag: list[str], expected: str) -> None:
    """Unmocked: the written truth follows the generator default unless the flag overrides it."""

    main([
        "simulate", "--seed", "0", "--out-dir", str(tmp_path),
        "--trajectory-mode", "magnitude", "--n-samples", "60", "--effect-size", "0.5",
        *flag,
    ])
    truth = json.loads((tmp_path / "truth.json").read_text())
    assert truth["magnitude_kind"] == expected


def test_simulate_cluster_mean_shift_fans_out_to_all_omics(tmp_path: Path) -> None:
    fake_dataset = _make_fake_simulate_dataset(n=30, seed=0)
    with patch("motco.simulations.reference.load_reference", return_value=object()), \
         patch("motco.simulations.semisynthetic.generate_semisynthetic_trajectory",
               return_value=fake_dataset) as generate_mock:
        main([
            "simulate",
            "--seed", "0",
            "--out-dir", str(tmp_path),
            "--cluster-mean-shift", "0.7",
        ])
    (traj_params,) = generate_mock.call_args.args
    assert traj_params.delta_methyl == 0.7
    assert traj_params.delta_expr == 0.7
    assert traj_params.delta_protein == 0.7


def test_simulate_per_omic_delta_overrides_cluster_mean_shift(tmp_path: Path) -> None:
    fake_dataset = _make_fake_simulate_dataset(n=30, seed=0)
    with patch("motco.simulations.reference.load_reference", return_value=object()), \
         patch("motco.simulations.semisynthetic.generate_semisynthetic_trajectory",
               return_value=fake_dataset) as generate_mock:
        main([
            "simulate",
            "--seed", "0",
            "--out-dir", str(tmp_path),
            "--cluster-mean-shift", "0.7",
            "--delta-expr", "1.2",
        ])
    (traj_params,) = generate_mock.call_args.args
    assert traj_params.delta_methyl == 0.7
    assert traj_params.delta_expr == 1.2
    assert traj_params.delta_protein == 0.7


def test_simulate_without_delta_flags_uses_default_shift(tmp_path: Path) -> None:
    fake_dataset = _make_fake_simulate_dataset(n=30, seed=0)
    with patch("motco.simulations.reference.load_reference", return_value=object()), \
         patch("motco.simulations.semisynthetic.generate_semisynthetic_trajectory",
               return_value=fake_dataset) as generate_mock:
        main(["simulate", "--seed", "0", "--out-dir", str(tmp_path)])
    (traj_params,) = generate_mock.call_args.args
    assert traj_params.delta_methyl == 2.0
    assert traj_params.delta_expr == 2.0
    assert traj_params.delta_protein == 2.0


def test_simulate_negative_delta_exits_before_load(tmp_path: Path) -> None:
    with patch("motco.simulations.reference.load_reference") as load_mock:
        with pytest.raises(SystemExit) as exc:
            main([
                "simulate",
                "--seed", "0",
                "--out-dir", str(tmp_path),
                "--delta-methyl", "-0.1",
            ])
    assert "-0.1" in str(exc.value)
    load_mock.assert_not_called()


def test_simulate_design_files_compatible_with_de(
    tmp_path: Path, mock_simulate_env: SemiSyntheticTrajectoryDataset
) -> None:
    sim_dir = tmp_path / "sim"
    main(["simulate", "--seed", "0", "--out-dir", str(sim_dir)])

    # Build a small Y matrix with correct row count
    n_samples = mock_simulate_env.metadata.shape[0]
    rng = np.random.default_rng(0)
    Y = pd.DataFrame(rng.standard_normal((n_samples, 5)))
    Y_path = tmp_path / "Y.csv"
    Y.to_csv(Y_path, index=False)

    out_json = tmp_path / "de_result.json"
    main([
        "de",
        "--Y", str(Y_path),
        "--model-matrix", str(sim_dir / "model_full.csv"),
        "--ls-means", str(sim_dir / "ls_means.csv"),
        "--contrast", str(sim_dir / "contrast.json"),
        "--out-json", str(out_json),
    ])
    assert out_json.exists()
    result = json.loads(out_json.read_text())
    assert {"deltas", "angles", "shapes"} <= result.keys()
