"""`motco analyze` and the `motco de --orientation-estimator` flag."""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from motco.cli import main
from motco.stats.design import build_ls_means, get_model_matrix
from motco.stats.trajectory import estimate_difference
from tests.test_analysis import STAGES, make_real_like


@pytest.fixture()
def analyze_files(tmp_path: Path) -> dict[str, Path]:
    metadata, blocks = make_real_like()
    paths = {"meta": tmp_path / "meta.csv"}
    metadata.to_csv(paths["meta"], index=False)
    for name, frame in blocks.items():
        paths[name] = tmp_path / f"{name}.csv"
        frame.rename_axis("donor").to_csv(paths[name])
    return paths


def _analyze_args(files: dict[str, Path], out_dir: Path, *extra: str) -> list[str]:
    return [
        "analyze",
        "--block", f"rna={files['rna']}",
        "--block", f"atac={files['atac']}",
        "--metadata", str(files["meta"]),
        "--sample-col", "donor",
        "--group-col", "Sex",
        "--stage-col", "ADNC",
        "--out-dir", str(out_dir),
        "--cv2-splits", "3", "--cv1-splits", "2", "--n-repeats", "2", "--max-components", "5",
        *extra,
    ]


def test_analyze_two_blocks(tmp_path: Path, analyze_files, capsys: pytest.CaptureFixture) -> None:
    out_dir = tmp_path / "out"
    main(_analyze_args(analyze_files, out_dir, "--stage-order", ",".join(STAGES), "--permutations", "19",
                       "--seed", "1"))
    record = json.loads((out_dir / "result.json").read_text())
    assert set(record["p_values"]) == {"delta", "angle", "shape"}
    assert record["design"]["stage_levels"] == STAGES
    printed = capsys.readouterr().out
    assert f"{record['latent']['selected_lv']} component(s) (cv)" in printed
    assert "net_displacement" in printed
    assert str(out_dir / "result.json") in printed
    for name in ("latent_scores.csv", "pls_loadings.csv", "stage_means.csv"):
        assert (out_dir / name).exists()


def test_analyze_refuses_existing_result(tmp_path: Path, analyze_files) -> None:
    out_dir = tmp_path / "out"
    args = _analyze_args(analyze_files, out_dir, "--stage-order", ",".join(STAGES), "--forced-components", "2")
    main(args)
    with pytest.raises(SystemExit, match="already exists"):
        main(args)
    main([*args, "--overwrite"])


def test_analyze_undeclared_transform_block(tmp_path: Path, analyze_files) -> None:
    args = _analyze_args(analyze_files, tmp_path / "out", "--stage-order", ",".join(STAGES),
                         "--transform", "methyl=beta_to_mvalue")
    with pytest.raises(SystemExit, match="undeclared block\\(s\\) \\['methyl'\\]"):
        main(args)


def test_analyze_missing_stage_order(tmp_path: Path, analyze_files) -> None:
    with pytest.raises(SystemExit, match="--stage-order is required"):
        main(_analyze_args(analyze_files, tmp_path / "out"))


def test_analyze_validation_error_is_clean(tmp_path: Path, analyze_files) -> None:
    with pytest.raises(SystemExit, match="Error: stage_order must list every observed stage"):
        main(_analyze_args(analyze_files, tmp_path / "out", "--stage-order", "Low,High"))


def test_analyze_malformed_block(tmp_path: Path, analyze_files) -> None:
    args = _analyze_args(analyze_files, tmp_path / "out", "--stage-order", ",".join(STAGES))
    args[2] = str(analyze_files["rna"])
    with pytest.raises(SystemExit, match="--block expects NAME=VALUE"):
        main(args)


# ── de --orientation-estimator ───────────────────────────────────────────────


@pytest.fixture()
def de_three_stage(tmp_path: Path) -> dict[str, Path]:
    rng = np.random.default_rng(5)
    factors = pd.DataFrame({"group": ["A"] * 15 + ["B"] * 15, "level": ["t0", "t1", "t2"] * 10})
    Y = rng.standard_normal((30, 4)) + np.outer(factors["level"].map({"t0": 0, "t1": 1, "t2": 3}), [1, 0.5, 0, 0])
    model = get_model_matrix(factors, group_col="group", level_col="level", full=True)
    reduced = get_model_matrix(factors, group_col="group", level_col="level", full=False)
    paths = {name: tmp_path / f"{name}.csv" for name in ("Y", "model", "reduced", "ls")}
    pd.DataFrame(Y).to_csv(paths["Y"], index=False)
    pd.DataFrame(model).to_csv(paths["model"], index=False)
    pd.DataFrame(reduced).to_csv(paths["reduced"], index=False)
    pd.DataFrame(build_ls_means(["A", "B"], ["t0", "t1", "t2"], full=True)).to_csv(paths["ls"], index=False)
    paths["contrast"] = tmp_path / "contrast.json"
    paths["contrast"].write_text(json.dumps([[0, 1, 2], [3, 4, 5]]))
    return paths


def _de(files: dict[str, Path], out: Path, *extra: str) -> dict:
    main([
        "de", "--Y", str(files["Y"]), "--model-matrix", str(files["model"]), "--ls-means", str(files["ls"]),
        "--contrast", str(files["contrast"]), "--out-json", str(out), *extra,
    ])
    return json.loads(out.read_text())


def test_de_default_output_is_unchanged(tmp_path: Path, de_three_stage) -> None:
    out = _de(de_three_stage, tmp_path / "default.json")
    assert set(out) == {"deltas", "angles", "shapes"}
    Y = pd.read_csv(de_three_stage["Y"]).values
    X = pd.read_csv(de_three_stage["model"]).values
    LS = pd.read_csv(de_three_stage["ls"]).values
    deltas, angles, shapes = estimate_difference(Y=Y, model_matrix=X, LS_means=LS, contrast=[[0, 1, 2], [3, 4, 5]])
    assert out == {
        "deltas": np.asarray(deltas).tolist(),
        "angles": np.asarray(angles).tolist(),
        "shapes": np.asarray(shapes).tolist(),
    }


def test_de_net_displacement_changes_only_angles(tmp_path: Path, de_three_stage) -> None:
    default = _de(de_three_stage, tmp_path / "default.json")
    net = _de(de_three_stage, tmp_path / "net.json", "--orientation-estimator", "net_displacement")
    assert net["orientation_estimator"] == "net_displacement"
    assert net["deltas"] == default["deltas"]
    assert net["shapes"] == default["shapes"]
    assert net["angles"] != default["angles"]
    explicit_pc1 = _de(de_three_stage, tmp_path / "pc1.json", "--orientation-estimator", "pc1")
    assert explicit_pc1["orientation_estimator"] == "pc1"
    assert explicit_pc1["angles"] == default["angles"]


def test_de_rrpp_uses_net_displacement(tmp_path: Path, de_three_stage) -> None:
    out = tmp_path / "rrpp.json"
    main([
        "de", "--Y", str(de_three_stage["Y"]), "--model-full", str(de_three_stage["model"]),
        "--model-reduced", str(de_three_stage["reduced"]), "--ls-means", str(de_three_stage["ls"]),
        "--contrast", str(de_three_stage["contrast"]), "--rrpp-permutations", "3",
        "--orientation-estimator", "net_displacement", "--out-json", str(out),
    ])
    result = json.loads(out.read_text())
    assert result["orientation_estimator"] == "net_displacement"
    assert len(result["angles"]) == 3


def test_de_unknown_estimator_exits(tmp_path: Path, de_three_stage, capsys: pytest.CaptureFixture) -> None:
    with pytest.raises(SystemExit) as exc:
        _de(de_three_stage, tmp_path / "bad.json", "--orientation-estimator", "procrustes")
    assert exc.value.code != 0
    err = capsys.readouterr().err
    assert "pc1" in err and "net_displacement" in err


@pytest.fixture()
def covariate_files(tmp_path: Path) -> dict[str, Path]:
    metadata, blocks = make_real_like()
    rng = np.random.default_rng(9)
    n = metadata.shape[0]
    metadata["age_at_death"] = 80.0 + 6.0 * rng.normal(size=n)
    metadata["pmi"] = rng.uniform(2.0, 12.0, size=n)
    metadata.loc[[3, 17], "pmi"] = np.nan
    metadata["education_level"] = rng.choice(["College", "Graduate", "HighSchool"], size=n)
    paths = {"meta": tmp_path / "meta.csv"}
    metadata.to_csv(paths["meta"], index=False)
    for name, frame in blocks.items():
        paths[name] = tmp_path / f"{name}.csv"
        frame.rename_axis("donor").to_csv(paths[name])
    return paths


def test_analyze_with_covariates(tmp_path: Path, covariate_files, capsys: pytest.CaptureFixture) -> None:
    out_dir = tmp_path / "out"
    main(_analyze_args(covariate_files, out_dir, "--stage-order", ",".join(STAGES), "--forced-components", "2",
                       "--covariate", "age_at_death", "--covariate", "pmi",
                       "--categorical-covariate", "education_level", "--missing-covariates", "complete_case"))
    block = json.loads((out_dir / "result.json").read_text())["covariates"]
    assert block["columns"] == ["age_at_death", "pmi", "education_level"]
    assert block["kinds"] == {"age_at_death": "numeric", "pmi": "numeric", "education_level": "categorical"}
    assert block["n_dropped"] == 2 and block["missing_policy"] == "complete_case"
    printed = capsys.readouterr().out
    assert "Covariates: age_at_death (numeric), pmi (numeric), education_level (categorical)" in printed
    assert "dropped samples: 2" in printed


def test_analyze_missing_covariate_column(tmp_path: Path, analyze_files) -> None:
    args = _analyze_args(analyze_files, tmp_path / "out", "--stage-order", ",".join(STAGES), "--covariate", "rin")
    with pytest.raises(SystemExit, match="'rin'"):
        main(args)


def test_analyze_missing_covariate_value_errors_by_default(tmp_path: Path, covariate_files) -> None:
    args = _analyze_args(covariate_files, tmp_path / "out", "--stage-order", ",".join(STAGES), "--covariate", "pmi")
    with pytest.raises(SystemExit, match="'pmi' has missing values"):
        main(args)
