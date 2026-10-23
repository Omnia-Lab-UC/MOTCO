"""End-to-end smoke test for the trajectory power study pipeline.

Uses a tiny config and a mock evaluator to validate runner → merge → report
without invoking the real InterSIM + RRPP pipeline.
"""

from __future__ import annotations

import dataclasses
import hashlib
import json
import shutil
import sys
from pathlib import Path
from unittest.mock import patch

import matplotlib
import numpy as np
import pandas as pd

matplotlib.use("Agg")

from motco.simulations import (
    SemiSyntheticTrajectoryParams,
    SimulationEvaluationParams,
    SimulationEvaluationResult,
)
from motco.simulations.study import (
    AcceptanceTargets,
    PowerMonotonicityTarget,
    SpecificityTarget,
    StudyConfig,
    TypeIControlTarget,
    dump_study_config,
    enumerate_study,
)
from motco.simulations.study.merge import discover_shard_paths, merge_shards
from motco.simulations.study.report import (
    build_report_frames,
    write_report_csvs,
)
from motco.simulations.study.sharding import run_shard
from motco.simulations.study.summary import (
    summarize_combined_rule,
    summarize_study,
)
from motco.simulations.study.targets import evaluate_targets, write_target_report


def _smoke_config() -> StudyConfig:
    return StudyConfig(
        generator=SemiSyntheticTrajectoryParams(seed=2, trajectory_mode="magnitude", group_effect_size=0.1),
        # PC1 and the unfloored rank rule pinned: the byte-identity checks below
        # compare against pre-flip outputs, whose cell ids hash neither default.
        evaluation=SimulationEvaluationParams(
            integration_method="concat", permutations=0, seed=3, orientation_estimator="pc1", min_components=1
        ),
        trajectory_modes=("magnitude",),
        effect_sizes=(0.1, 0.5),
        n_replicates=2,
        base_seed=42,
        alpha=0.05,
        acceptance=AcceptanceTargets(
            type_i=(TypeIControlTarget(alpha=0.05, se_tolerance=2.0),),
            power=(PowerMonotonicityTarget(trajectory_mode="magnitude", statistic="delta", min_power_at_top=0.3),),
        ),
    )


def _mock_evaluator(generator_params, evaluation_params) -> SimulationEvaluationResult:
    """A mock evaluator: nulls reject ~alpha, power cells reject more as effect grows."""

    rng = np.random.default_rng(generator_params.seed)
    matrix = np.zeros((2, 2), dtype=float)
    mode = generator_params.trajectory_mode
    effect = generator_params.group_effect_size
    if mode == "none" or effect == 0.0:
        delta_p = float(rng.uniform(0.0, 1.0))
        angle_p = float(rng.uniform(0.0, 1.0))
    else:
        # higher effect → smaller p-values on average
        bias = max(0.0, 0.6 - 1.2 * effect)
        delta_p = float(np.clip(rng.uniform(0.0, 1.0) - 0.5 + bias, 0.0, 1.0))
        angle_p = float(np.clip(rng.uniform(0.0, 1.0) - 0.3 + bias, 0.0, 1.0))
    return SimulationEvaluationResult(
        observed_deltas=matrix,
        observed_angles=matrix,
        observed_shapes=matrix,
        pair_statistics={"delta": 1.0, "angle": 2.0, "shape": float("nan")},
        p_values={"delta": delta_p, "angle": angle_p},
        latent_matrix_metadata={"integration_method": "concat"},
        truth_metadata={"trajectory_mode": mode, "effect_size": effect},
        runtime_metadata={"runtime_seconds": 0.01},
        evaluation_params=evaluation_params,
        group_levels=["A", "B"],
        stage_levels=["0", "1"],
        contrast=[[0, 1], [2, 3]],
    )


def test_end_to_end_smoke(tmp_path: Path) -> None:
    config = _smoke_config()
    config_path = tmp_path / "smoke.json"
    dump_study_config(config, config_path)

    grid = enumerate_study(config)
    n_shards = 3
    for shard_index in range(n_shards):
        run_shard(
            grid,
            shard_index=shard_index,
            n_shards=n_shards,
            out_dir=tmp_path,
            evaluator=_mock_evaluator,
        )

    shards = discover_shard_paths(tmp_path)
    assert len(shards) == n_shards
    merged_path = tmp_path / "merged.jsonl"
    records = merge_shards(shards, out_path=merged_path)
    assert merged_path.exists()
    total_units = sum(c.n_replicates for c in grid.cells)
    assert len(records) == total_units

    per_stat = summarize_study(records, alpha=config.alpha)
    combined = summarize_combined_rule(records, alpha=config.alpha)
    frames = build_report_frames(per_stat, combined, records, alpha=config.alpha)
    report_dir = tmp_path / "report"
    csv_paths = write_report_csvs(frames, report_dir)
    for path in csv_paths.values():
        assert path.exists()
        assert pd.read_csv(path) is not None

    evaluations = evaluate_targets(config.acceptance, per_stat, records)
    paths = write_target_report(evaluations, report_dir)
    payload = json.loads(paths["json"].read_text(encoding="utf-8"))
    assert any(ev["target_kind"] == "type_i_control" for ev in payload)


def test_run_study_shard_script_invocation(tmp_path: Path, monkeypatch) -> None:
    """The runner script wires through to run_shard and emits a JSONL."""

    config = _smoke_config()
    config_path = tmp_path / "smoke.json"
    dump_study_config(config, config_path)

    sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
    try:
        import run_study_shard  # type: ignore[import-not-found]
    finally:
        sys.path.pop(0)

    with patch("run_study_shard.run_shard") as mock_run:
        mock_run.return_value = []
        rc = run_study_shard.main(
            [
                "--config",
                str(config_path),
                "--out-dir",
                str(tmp_path / "shards"),
                "--shard-index",
                "0",
                "--n-shards",
                "1",
                "--n-jobs",
                "1",
            ]
        )
    assert rc == 0
    mock_run.assert_called_once()


def test_motco_study_script_merge_and_report(tmp_path: Path) -> None:
    """The interactive CLI runs merge then report end-to-end."""

    config = _smoke_config()
    config_path = tmp_path / "smoke.json"
    dump_study_config(config, config_path)
    grid = enumerate_study(config)
    for shard_index in range(2):
        run_shard(
            grid,
            shard_index=shard_index,
            n_shards=2,
            out_dir=tmp_path,
            evaluator=_mock_evaluator,
        )

    sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
    try:
        import motco_study  # type: ignore[import-not-found]
    finally:
        sys.path.pop(0)

    assert motco_study.main(["merge", "--out-dir", str(tmp_path)]) == 0
    assert (tmp_path / "merged.jsonl").exists()

    assert (
        motco_study.main(
            ["report", "--config", str(config_path), "--out-dir", str(tmp_path)]
        )
        == 0
    )
    report_dir = tmp_path / "report"
    for name in (
        "specificity_matrix.csv",
        "power_curves.csv",
        "type_i_table.csv",
        "specificity_matrix.png",
        "power_curves.png",
        "type_i.png",
        "acceptance_report.csv",
        "acceptance_report.json",
    ):
        assert (report_dir / name).exists(), f"missing {name}"


def test_end_to_end_design_grid(tmp_path: Path) -> None:
    """Shards, merge, resume, and report over a tiny crossed design grid."""

    from motco.simulations.study import DesignGrid, DesignPointDecisionRule, MatchedSeedPolicy
    from motco.simulations.study.enumerate import DESIGN_PHASE
    from motco.simulations.study.targets import evaluate_design_point_decision

    config = StudyConfig(
        generator=SemiSyntheticTrajectoryParams(seed=2, trajectory_mode="magnitude", n_samples=60, p_dmp=0.1),
        evaluation=SimulationEvaluationParams(integration_method="concat", permutations=0, seed=3),
        trajectory_modes=("magnitude", "orientation"),
        effect_sizes=(0.0, 0.5, 1.0),
        n_replicates=2,
        base_seed=7,
        alpha=0.05,
        matched_seeds=MatchedSeedPolicy(enabled=True, primary_family="fam"),
        design_grid=DesignGrid(
            axes={"generator.baseline_continuity": (0.0, 0.5), "generator.n_samples": (60, 120)}
        ),
        acceptance=AcceptanceTargets(
            type_i=(TypeIControlTarget(alpha=0.05),),
            design_point=DesignPointDecisionRule(
                trajectory_mode="orientation",
                statistic="angle",
                min_power_at_top=0.5,
                prefer=("generator.n_samples", "generator.baseline_continuity"),
            ),
        ),
    )
    grid = enumerate_study(config)
    assert sum(cell.phase == DESIGN_PHASE for cell in grid.cells) == 3 * (1 + 2 * 2)

    for shard_index in range(2):
        run_shard(grid, shard_index=shard_index, n_shards=2, out_dir=tmp_path, evaluator=_mock_evaluator)
    records = merge_shards(discover_shard_paths(tmp_path), out_path=tmp_path / "merged.jsonl")
    assert len(records) == sum(c.n_replicates for c in grid.cells)
    assert {r.phase for r in records} >= {"power_primary", DESIGN_PHASE, "type_i_baseline"}

    # Resume: nothing left to do, no duplicates.
    run_shard(grid, shard_index=0, n_shards=2, out_dir=tmp_path, evaluator=_mock_evaluator)
    records = merge_shards(discover_shard_paths(tmp_path), out_path=tmp_path / "merged.jsonl")
    assert len(records) == sum(c.n_replicates for c in grid.cells)

    per_stat = summarize_study(records, alpha=config.alpha)
    combined = summarize_combined_rule(records, alpha=config.alpha)
    frames = build_report_frames(per_stat, combined, records, alpha=config.alpha)
    paths = write_report_csvs(frames, tmp_path / "report")
    design = pd.read_csv(paths["design_point_operating"])
    assert set(design["generator.n_samples"]) == {60, 120}
    assert set(design["generator.baseline_continuity"]) == {0.0, 0.5}
    # Baseline readers never see design cells (the translation control is a
    # baseline Type I cell and has always been part of the curves frame).
    assert DESIGN_PHASE not in set(frames.power_curves["phase"])
    assert DESIGN_PHASE not in set(frames.specificity_matrix["phase"])

    decision = evaluate_design_point_decision(records, config.acceptance.design_point, alpha=config.alpha)
    assert decision.verdict in {"chosen", "revise_claim"}
    assert len(decision.columns) == 4


def _mock_pls_evaluator(generator_params, evaluation_params) -> SimulationEvaluationResult:
    """Mock evaluator that honours ``forced_components`` the way ``_pls_integration`` does."""

    result = _mock_evaluator(generator_params, evaluation_params)
    forced = dict(evaluation_params.integration_params).get("forced_components")
    metadata = {
        "integration_method": "pls",
        "component_selection": "forced" if forced is not None else "cv",
        "selected_lv": int(forced) if forced is not None else 3,
    }
    return SimulationEvaluationResult(
        **{**result.__dict__, "latent_matrix_metadata": metadata, "evaluation_params": evaluation_params}
    )


def test_end_to_end_rank_ladder(tmp_path: Path) -> None:
    """Shards, merge, resume, and report over a tiny retained-rank ladder."""

    from motco.simulations.study import (
        DesignGrid,
        MatchedSeedPolicy,
        RankDecisionRule,
        StatisticPair,
        render_rank_ladder,
    )
    from motco.simulations.study.enumerate import DESIGN_PHASE
    from motco.simulations.study.spectrum import RANK_AXIS
    from motco.simulations.study.targets import evaluate_rank_decision, write_rank_decision

    config = StudyConfig(
        generator=SemiSyntheticTrajectoryParams(seed=2, trajectory_mode="magnitude", n_samples=60, p_dmp=0.1),
        evaluation=SimulationEvaluationParams(
            integration_method="pls",
            integration_params={"cv1_splits": 3, "cv2_splits": 4, "n_repeats": 5, "max_components": 20},
            permutations=9,
            seed=3,
        ),
        trajectory_modes=("magnitude", "orientation", "shape", "translation"),
        effect_sizes=(0.0, 0.5, 1.0),
        n_replicates=1,
        base_seed=7,
        alpha=0.05,
        matched_seeds=MatchedSeedPolicy(enabled=True, primary_family="ladder"),
        design_grid=DesignGrid(axes={RANK_AXIS: (None, 3)}),
        acceptance=AcceptanceTargets(
            type_i=(TypeIControlTarget(alpha=0.05),),
            rank_decision=RankDecisionRule(
                axis=RANK_AXIS,
                target=StatisticPair("orientation", "angle"),
                protected=(StatisticPair("magnitude", "delta"), StatisticPair("shape", "shape")),
            ),
        ),
    )
    grid = enumerate_study(config)
    forced_cells = [cell for cell in grid.cells if cell.phase == DESIGN_PHASE]
    assert len(forced_cells) == 1 + 4 * 2
    assert all(cell.evaluation_params.integration_params["forced_components"] == 3 for cell in forced_cells)

    for shard_index in range(2):
        run_shard(grid, shard_index=shard_index, n_shards=2, out_dir=tmp_path, evaluator=_mock_pls_evaluator)
    records = merge_shards(discover_shard_paths(tmp_path), out_path=tmp_path / "merged.jsonl")
    assert len(records) == sum(c.n_replicates for c in grid.cells)
    assert {r.phase for r in records} >= {"power_primary", DESIGN_PHASE, "type_i_baseline"}

    forced = [r for r in records if r.phase == DESIGN_PHASE]
    assert forced and all(r.integration_metadata["component_selection"] == "forced" for r in forced)
    assert all(r.integration_metadata["selected_lv"] == 3 for r in forced)
    cv = [r for r in records if r.phase == "power_primary"]
    assert cv and all(r.integration_metadata["component_selection"] == "cv" for r in cv)
    # Same dataset, different measurement: paired generator seeds across the two columns.
    by_key: dict[tuple, set[int]] = {}
    for r in cv + forced:
        key = (r.cell_metadata["trajectory_mode"], r.cell_metadata["effect_size"], r.replicate_index)
        by_key.setdefault(key, set()).add(r.generator_seed)
    assert all(len(seeds) == 1 for seeds in by_key.values())

    # Resume: nothing left to do, no duplicates.
    run_shard(grid, shard_index=0, n_shards=2, out_dir=tmp_path, evaluator=_mock_pls_evaluator)
    records = merge_shards(discover_shard_paths(tmp_path), out_path=tmp_path / "merged.jsonl")
    assert len(records) == sum(c.n_replicates for c in grid.cells)

    per_stat = summarize_study(records, alpha=config.alpha)
    combined = summarize_combined_rule(records, alpha=config.alpha)
    frames = build_report_frames(per_stat, combined, records, alpha=config.alpha)
    paths = write_report_csvs(frames, tmp_path / "report")
    design = pd.read_csv(paths["design_point_operating"])
    assert set(design["component_selection"]) == {"cv", "forced"}
    assert set(design.loc[design["component_selection"] == "forced", RANK_AXIS]) == {3}
    assert DESIGN_PHASE not in set(frames.power_curves["phase"])
    assert render_rank_ladder(frames.design_point_operating, tmp_path / "report" / "rank_ladder.png") is not None

    decision = evaluate_rank_decision(records, config.acceptance.rank_decision, alpha=config.alpha)
    assert decision.verdict in {"keep_cv", "adopt_fixed_rank"}
    assert [c.rank for c in decision.columns] == [None, 3]
    out = write_rank_decision(decision, tmp_path / "report")
    assert out["rank_decision"].exists() and out["rank_decision_csv"].exists()


# --- acceptance target bound form: default-form byte identity --------------------

# sha256 of the smoke study's acceptance_report.* captured on the code *before*
# the ``bound`` key existed (Type I + specificity targets, no ``bound`` key).
_PRE_BOUND_ACCEPTANCE_SHA256 = {
    "acceptance_report.csv": "a86151a9d876a20056ba361f8bfad7369f50ff20be94b3083421e2f1b12a2bd2",
    "acceptance_report.json": "8c61de31aba5a8a642eb250428dcf890aeef71005d79b2299147e23edba0028d",
}


def test_default_bound_form_regenerates_acceptance_report_byte_identically(tmp_path: Path) -> None:
    """A keyless config and an explicit ``two_sided`` config report the pre-change bytes."""

    config = _smoke_config()
    config = dataclasses.replace(
        config,
        acceptance=dataclasses.replace(
            config.acceptance,
            specificity=(SpecificityTarget(trajectory_mode="magnitude", statistic="angle", alpha=0.05),),
        ),
    )
    run_dir = tmp_path / "run"
    dump_study_config(config, run_dir / "smoke.json")
    grid = enumerate_study(config)
    for shard_index in range(3):
        run_shard(grid, shard_index=shard_index, n_shards=3, out_dir=run_dir, evaluator=_mock_evaluator)

    sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
    try:
        import motco_study  # type: ignore[import-not-found]
    finally:
        sys.path.pop(0)
    assert motco_study.main(["merge", "--out-dir", str(run_dir)]) == 0

    explicit = json.loads((run_dir / "smoke.json").read_text(encoding="utf-8"))
    targets = explicit["acceptance"]["type_i"] + explicit["acceptance"]["specificity"]
    assert targets and all(t["bound"] == "two_sided" for t in targets)
    keyless = json.loads(json.dumps(explicit))
    for entry in keyless["acceptance"]["type_i"] + keyless["acceptance"]["specificity"]:
        del entry["bound"]

    reports: dict[str, dict[str, bytes]] = {}
    for label, payload in (("keyless", keyless), ("explicit", explicit)):
        out_dir = tmp_path / label
        out_dir.mkdir()
        shutil.copy(run_dir / "merged.jsonl", out_dir / "merged.jsonl")
        config_path = out_dir / "config.json"
        config_path.write_text(json.dumps(payload), encoding="utf-8")
        assert motco_study.main(["report", "--config", str(config_path), "--out-dir", str(out_dir)]) == 0
        reports[label] = {
            name: (out_dir / "report" / name).read_bytes() for name in _PRE_BOUND_ACCEPTANCE_SHA256
        }

    assert reports["keyless"] == reports["explicit"]
    for name, digest in _PRE_BOUND_ACCEPTANCE_SHA256.items():
        assert hashlib.sha256(reports["keyless"][name]).hexdigest() == digest, name


# --- report contract: CLI wiring and runner refusal -----------------------------


def _import_script(name: str):
    sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
    try:
        return __import__(name)
    finally:
        sys.path.pop(0)


def _contract_config(n_jobs_override: str = "forbid") -> StudyConfig:
    from motco.simulations.study import ReportContract

    return StudyConfig(
        generator=SemiSyntheticTrajectoryParams(seed=2, trajectory_mode="magnitude", group_effect_size=0.1),
        evaluation=SimulationEvaluationParams(integration_method="concat", permutations=0, seed=3, n_jobs=1),
        trajectory_modes=("magnitude",),
        effect_sizes=(0.1, 0.5),
        n_replicates=1,
        base_seed=42,
        acceptance=AcceptanceTargets(type_i=(TypeIControlTarget(alpha=0.05),)),
        report_contract=ReportContract(driver_component="observed", n_jobs_override=n_jobs_override),
    )


def test_report_writes_contract_files_only_when_declared(tmp_path: Path) -> None:
    motco_study = _import_script("motco_study")

    plain = _smoke_config()
    plain_dir = tmp_path / "plain"
    dump_study_config(plain, plain_dir / "config.json")
    grid = enumerate_study(plain)
    run_shard(grid, shard_index=0, n_shards=1, out_dir=plain_dir, evaluator=_mock_evaluator)
    assert motco_study.main(["merge", "--out-dir", str(plain_dir)]) == 0
    assert motco_study.main(["report", "--config", str(plain_dir / "config.json"), "--out-dir", str(plain_dir)]) == 0
    plain_files = {p.name for p in (plain_dir / "report").iterdir()}
    assert "driver_report.csv" not in plain_files and "report_contract.json" not in plain_files

    contracted = _contract_config()
    c_dir = tmp_path / "contract"
    dump_study_config(contracted, c_dir / "config.json")
    grid = enumerate_study(contracted)
    run_shard(grid, shard_index=0, n_shards=1, out_dir=c_dir, evaluator=_mock_evaluator)
    assert motco_study.main(["merge", "--out-dir", str(c_dir)]) == 0
    assert motco_study.main(["report", "--config", str(c_dir / "config.json"), "--out-dir", str(c_dir)]) == 0
    files = {p.name for p in (c_dir / "report").iterdir()}
    assert {"driver_report.csv", "report_contract.json"} <= files
    # Same study apart from the contract: the shared files are the same set.
    assert files - {"driver_report.csv", "report_contract.json"} == plain_files
    echo = json.loads((c_dir / "report" / "report_contract.json").read_text(encoding="utf-8"))
    assert echo["driver_component"] == "observed" and echo["n_jobs_override"] == "forbid"
    assert echo["config_n_jobs"] == 1
    driver = pd.read_csv(c_dir / "report" / "driver_report.csv")
    assert list(driver.columns)[:3] == ["trajectory_mode", "effect_size", "transition_id"]


def test_runner_refuses_a_differing_n_jobs_under_forbid(tmp_path: Path, capsys) -> None:
    run_study_shard = _import_script("run_study_shard")
    config_path = tmp_path / "forbid.json"
    dump_study_config(_contract_config("forbid"), config_path)
    out_dir = tmp_path / "shards"
    common = ["--config", str(config_path), "--out-dir", str(out_dir), "--shard-index", "0", "--n-shards", "1"]

    with patch("run_study_shard.enumerate_study") as mock_enumerate, patch("run_study_shard.run_shard") as mock_run:
        rc = run_study_shard.main([*common, "--n-jobs", "4"])
    assert rc == 2
    mock_enumerate.assert_not_called()
    mock_run.assert_not_called()
    assert not out_dir.exists()
    err = capsys.readouterr().err
    assert "--n-jobs 4" in err and "evaluation.n_jobs 1" in err
    assert "n_jobs_override: forbid" in err and "parameter signature" in err

    # Equal to the config value: accepted, no warning.
    with patch("run_study_shard.run_shard") as mock_run:
        mock_run.return_value = []
        assert run_study_shard.main([*common, "--n-jobs", "1"]) == 0
    mock_run.assert_called_once()
    assert capsys.readouterr().err == ""

    # Absent: runs the shard for real.
    with patch("run_study_shard.run_shard") as mock_run:
        mock_run.return_value = []
        assert run_study_shard.main(common) == 0
    mock_run.assert_called_once()


def test_runner_warns_and_proceeds_under_warn_or_no_contract(tmp_path: Path, capsys) -> None:
    run_study_shard = _import_script("run_study_shard")
    for name, config in (("warn.json", _contract_config("warn")), ("plain.json", _smoke_config())):
        config_path = tmp_path / name
        dump_study_config(config, config_path)
        common = ["--config", str(config_path), "--out-dir", str(tmp_path / "out"), "--shard-index", "0"]
        with patch("run_study_shard.run_shard") as mock_run:
            mock_run.return_value = []
            assert run_study_shard.main([*common, "--n-shards", "1", "--n-jobs", "4"]) == 0
        mock_run.assert_called_once()
        assert mock_run.call_args.args[0].cells[0].evaluation_params.n_jobs == 4
        assert "WARNING: --n-jobs 4 overrides" in capsys.readouterr().err
