"""``scripts/case_study_diagnosis.py``: parent binding, oracle, factorial, correction, latent arms, verdicts."""

from __future__ import annotations

import itertools
import json
import math
import shutil
import sys
from dataclasses import replace
from pathlib import Path

import numpy as np
import pytest


def _script():
    sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
    try:
        return __import__("case_study_diagnosis")
    finally:
        sys.path.pop(0)


diag = _script()
cs = diag.cs
SMALL_CV = {"cv1_splits": 2, "cv2_splits": 3, "n_repeats": 2, "max_components": 5, "random_state": 1203}
SIX = ["z", "x1", "x2", "x3", "x4", "x5"]
STATISTICS = ("delta", "angle", "shape")
P_RED = {"oracle_noise0": 4, "oracle_noise2": 6, "oracle_noise5": 9, "z_noise0": 5, "z_noise2": 7, "z_noise5": 10}


@pytest.fixture(scope="module")
def small_cv():
    """The parent check's CV knobs at test size, for its own runs and for the diagnosis's rebuild check."""

    with pytest.MonkeyPatch.context() as patch:
        patch.setattr(cs.d10, "PLS_CV_PARAMS", dict(SMALL_CV))
        yield


@pytest.fixture(scope="module")
def parent(tmp_path_factory, small_cv) -> Path:
    out_dir = tmp_path_factory.mktemp("parent")
    assert cs.main(["run", "--reps", "2", "--perms", "9", "--out-dir", str(out_dir)]) == 0
    return out_dir


@pytest.fixture(scope="module")
def diagnosis(tmp_path_factory, parent) -> Path:
    out_dir = tmp_path_factory.mktemp("diagnosis")
    assert diag.main(["run", "--parent", str(parent), "--out-dir", str(out_dir)]) == 0
    return out_dir


@pytest.fixture(scope="module")
def smoke(tmp_path_factory, parent) -> Path:
    out_dir = tmp_path_factory.mktemp("smoke")
    assert diag.main(["run", "--parent", str(parent), "--reps", "1", "--perms", "5", "--out-dir", str(out_dir)]) == 0
    return out_dir


@pytest.fixture(scope="module")
def unit(diagnosis) -> dict:
    """The anchor unit of replicate 0: its build, datasets, covariate-blind analysis and oracle scores."""

    from motco import analysis

    settings = cs.load_settings(diagnosis)
    parent_settings = settings["parent"]["settings"]
    built = diag.build_unit("anchor", 0.0, 0, settings)
    metadata = built["metadata"]
    common = {
        "blocks": [analysis.OmicsBlock(name, built["blocks"][name], "identity") for name in parent_settings["layers"]],
        "metadata": metadata,
        "sample_col": "sample_id",
        "group_col": "group",
        "stage_col": "stage",
        "stage_order": [str(stage) for stage in sorted(set(metadata["stage"].tolist()))],
    }
    six = analysis.TrajectoryDataset(**common, covariate_cols=SIX)
    seed = int(built["seeds"]["rrpp"])
    blind = analysis.analyze_trajectories(six, diag.analysis_params(settings, seed))
    Y0, checks = diag.oracle_latent(built, blind, parent_settings["layers"])
    return {
        "settings": settings,
        "built": built,
        "common": common,
        "six": six,
        "seed": seed,
        "blind": blind,
        "Y0": Y0,
        "checks": checks,
        "design": diag.trajectory_design(six),
        "layers": list(parent_settings["layers"]),
    }


def _records(out_dir: Path, index: int = 0, count: int = 1) -> list[dict]:
    lines = cs.shard_path(out_dir, index, count).read_text(encoding="utf-8").splitlines()
    records = sorted((json.loads(line) for line in lines), key=lambda record: record["key"])
    for record in records:
        record.pop("runtime_seconds")
    return records


# ── Parent binding ───────────────────────────────────────────────────────────


def _copy_parent(parent: Path, target: Path, edit) -> Path:
    shutil.copytree(parent, target)
    path = target / "settings.json"
    stored = json.loads(path.read_text(encoding="utf-8"))
    edit(stored)
    path.write_text(json.dumps(stored, indent=1, sort_keys=True) + "\n", encoding="utf-8")
    return target


def test_tampered_parent_is_refused_and_writes_nothing(parent, tmp_path) -> None:
    def tamper(stored):
        stored["settings"]["alpha"] = 0.1  # the recorded hash no longer matches

    tampered = _copy_parent(parent, tmp_path / "tampered", tamper)
    with pytest.raises(SystemExit, match="fails its recorded hash"):
        diag.load_parent(tampered)
    out_dir = tmp_path / "out"
    with pytest.raises(SystemExit, match="tampered"):
        diag.main(["run", "--parent", str(tampered), "--out-dir", str(out_dir)])
    assert not out_dir.exists()


def test_parent_the_current_code_does_not_rebuild_is_refused(parent, tmp_path) -> None:
    def drift(stored):
        # A construction value no recorded argument sets: the hash is re-sealed, the rebuild differs.
        stored["settings"]["generator"]["p_dmp"] = 0.2
        stored["hash"] = cs.settings_hash(stored["settings"])

    drifted = _copy_parent(parent, tmp_path / "drifted", drift)
    with pytest.raises(SystemExit, match="rebuilds its settings to hash .*differing: generator"):
        diag.load_parent(drifted)


def test_resume_under_another_parent_is_refused(parent, diagnosis, tmp_path) -> None:
    other = tmp_path / "other"
    assert cs.main(["run", "--reps", "2", "--perms", "19", "--out-dir", str(other), "--limit", "0"]) == 0
    shard = cs.shard_path(diagnosis, 0, 1)
    before = shard.read_text(encoding="utf-8")
    with pytest.raises(SystemExit, match="written against parent"):
        diag.main(["run", "--parent", str(other), "--out-dir", str(diagnosis)])
    assert shard.read_text(encoding="utf-8") == before
    # The same parent at other diagnosis settings is refused too.
    with pytest.raises(SystemExit, match="different settings"):
        diag.main(["run", "--parent", str(parent), "--perms", "5", "--out-dir", str(diagnosis)])


def test_reps_and_perms_are_bounded_by_the_parent(parent) -> None:
    bound = diag.load_parent(parent)
    with pytest.raises(SystemExit, match="prefix"):
        diag.build_settings(parent, bound, reps=3)
    with pytest.raises(SystemExit, match="permutations|--perms"):
        diag.build_settings(parent, bound, perms=19)
    prefix = diag.build_settings(parent, bound, reps=1, perms=5)
    assert prefix["smoke"] is True and prefix["replicates"] == 1
    assert diag.build_settings(parent, bound)["smoke"] is False
    assert diag.units(prefix) == [("anchor", 0.0, 0), ("translation", 1.0, 0)]


def test_units_regenerate_the_parent(parent, diagnosis) -> None:
    theirs = {
        record["key"]: record
        for record in map(json.loads, (parent / "merged.jsonl").read_text(encoding="utf-8").splitlines())
    }
    ours = _records(diagnosis)
    assert sorted(record["key"] for record in ours) == sorted(theirs)
    for record in ours:
        assert record["seeds"] == theirs[record["key"]]["seeds"]
        assert record["covariate_cells"] == theirs[record["key"]]["covariate_cells"]
        assert record["group_labels"] == theirs[record["key"]]["group_labels"]
    settings = cs.load_settings(diagnosis)
    built = diag.build_unit("translation", 1.0, 1, settings)
    reference = cs.build_unit("translation", 1.0, 1, settings["parent"]["settings"])
    assert built["seeds"] == reference["seeds"]
    assert built["metadata"].equals(reference["metadata"])
    for name in settings["parent"]["settings"]["layers"]:
        np.testing.assert_array_equal(built["blocks"][name].to_numpy(), reference["blocks"][name].to_numpy())


def test_judged_arms_reproduce_the_parent(diagnosis) -> None:
    report = json.loads((diagnosis / "diagnosis_reproduction.json").read_text(encoding="utf-8"))
    assert report["status"] == "holds"
    assert report["scale"] == 10
    assert report["arms"] == {
        "z_noise5": {"parent_arm": "adjusted", "compared": 12, "mismatches": 0},
        "latent_stage_cv": {"parent_arm": "latent_stage", "compared": 12, "mismatches": 0},
    }
    assert report["mismatching_units"] == []


# ── Oracle ───────────────────────────────────────────────────────────────────


def test_reconstruction_reproduces_the_latent_scores(diagnosis, unit) -> None:
    for record in _records(diagnosis):
        checks = record["checks"]
        assert checks["reconstruction_max_abs_error"] <= diag.RECONSTRUCTION_TOLERANCE
        assert checks["latent_shift_rank_one_residual"] <= diag.RANK_ONE_TOLERANCE
        assert max(checks["injection_rank_one_residual"].values()) <= diag.RANK_ONE_TOLERANCE
        assert set(checks["injection_rank_one_residual"]) == set(unit["layers"])
    assert unit["checks"]["reconstruction_max_abs_error"] <= diag.RECONSTRUCTION_TOLERANCE


def test_oracle_removes_exactly_the_injected_effect(unit) -> None:
    built, blind = unit["built"], unit["blind"]
    parent = unit["settings"]["parent"]["settings"]
    z = built["metadata"]["z"].to_numpy(dtype=float)
    shift = blind.latent_scores.to_numpy(dtype=float) - unit["Y0"].to_numpy(dtype=float)
    # Direct construction: D10's injection c = s_z · w · sd (pre-injection SD) through the fitted linear map.
    project = diag.CovariateBlindProjection(built["blocks"], blind.pls_loadings, unit["layers"])
    clean = {name: built["blocks_without_injection"][name].to_numpy(dtype=float) for name in unit["layers"]}
    n_features = sum(values.shape[1] for values in clean.values())
    w, _ = cs.d10.feature_directions(n_features, parent)
    injection = np.concatenate(
        [
            parent["injection"]["s_z"] * w[offset : offset + values.shape[1]] * values.std(axis=0, ddof=1)
            for offset, values in zip(np.cumsum([0] + [v.shape[1] for v in clean.values()]), clean.values())
        ]
    )
    scales = np.concatenate([project.scalers[name].scale for name in unit["layers"]])
    direct = (injection / scales / project.x_std) @ project.rotations
    np.testing.assert_allclose(shift, np.outer(z, direct), atol=1e-9)
    assert np.linalg.norm(direct) > 0


def test_wrong_scaling_makes_the_unit_raise(unit, monkeypatch) -> None:
    monkeypatch.setattr(diag, "PLS_SCALE_DDOF", 0)
    with pytest.raises(RuntimeError, match="reconstructed covariate-blind projection"):
        diag.oracle_latent(unit["built"], unit["blind"], unit["layers"])


# ── Factorial and correction ─────────────────────────────────────────────────


def test_identity_assertion_catches_a_perturbed_design(unit, monkeypatch) -> None:
    record = unit["blind"].record
    arm = {"statistics": dict(record["statistics"]), "p_values": dict(record["p_values"])}
    diag.check_identity(arm, record)
    with pytest.raises(RuntimeError, match="differs from analyze_trajectories"):
        diag.check_identity({**arm, "p_values": {**arm["p_values"], "delta": 1.0}}, record)
    # A design that drops one noise column from the case-study model no longer matches the analysis.
    monkeypatch.setitem(diag.FACTORIAL, "z_noise5", (False, 4))
    with pytest.raises(RuntimeError, match="differs from analyze_trajectories"):
        diag.evaluate_unit("anchor", 0.0, 0, unit["settings"])


def test_noise_covariates_are_nested(diagnosis, unit) -> None:
    from motco.analysis import TrajectoryDataset

    expected = {
        "oracle_noise0": [], "oracle_noise2": ["x1", "x2"], "oracle_noise5": SIX[1:],
        "z_noise0": ["z"], "z_noise2": ["z", "x1", "x2"], "z_noise5": SIX,
    }
    for record in _records(diagnosis):
        for arm, columns in expected.items():
            assert record["arms"][arm]["covariate_design_columns"] == columns
            assert record["arms"][f"{arm}_corrected"]["covariate_design_columns"] == columns
        assert record["arms"]["latent_stage_cv"]["covariate_design_columns"] == SIX
    designs = {
        arm: TrajectoryDataset(**unit["common"], covariate_cols=diag.arm_covariates(arm)).covariate_design
        for arm in ("z_noise2", "z_noise5", "oracle_noise2", "oracle_noise5")
    }
    np.testing.assert_array_equal(designs["z_noise2"].matrix, designs["z_noise5"].matrix[:, :3])
    np.testing.assert_array_equal(designs["oracle_noise2"].matrix, designs["oracle_noise5"].matrix[:, :2])
    np.testing.assert_array_equal(designs["z_noise5"].matrix[:, 1:], designs["oracle_noise5"].matrix)
    # The noise columns are the parent's own draw.
    parent_unit = cs.build_unit("anchor", 0.0, 0, unit["settings"]["parent"]["settings"])
    np.testing.assert_array_equal(designs["z_noise5"].matrix[:, 1:], parent_unit["metadata"][SIX[1:]].to_numpy())


def test_arms_share_the_basis_and_the_permutation_draws(unit, monkeypatch) -> None:
    calls = []
    real = diag.RRPP

    def spy(Y, full, reduced, ls_means, contrast, **kwargs):
        calls.append({"Y": np.asarray(Y, dtype=float).copy(), "seed": kwargs["seed"],
                      "permutations": kwargs["permutations"], "n_jobs": kwargs["n_jobs"], "p_red": reduced.shape[1]})
        return real(Y, full, reduced, ls_means, contrast, **kwargs)

    monkeypatch.setattr(diag, "RRPP", spy)
    blind = unit["blind"]
    arms = diag.factorial_arms(
        blind.latent_scores,
        unit["Y0"],
        unit["common"],
        unit["design"],
        permutations=9,
        seed=unit["seed"],
        orientation_estimator="net_displacement",
        latent_fields=diag._latent_fields(blind.record),
    )
    assert len(calls) == 12  # six arms, each uncorrected then corrected
    assert {call["seed"] for call in calls} == {unit["seed"]}
    assert {call["permutations"] for call in calls} == {9}
    assert {call["n_jobs"] for call in calls} == {1}
    assert {call["Y"].shape for call in calls} == {blind.latent_scores.shape}
    for (arm, (oracle, _)), call in zip(diag.FACTORIAL.items(), calls[::2]):
        expected = unit["Y0"] if oracle else blind.latent_scores
        np.testing.assert_array_equal(call["Y"], expected.to_numpy(dtype=float))
        assert call["p_red"] == P_RED[arm]
    assert {arms[arm]["selected_lv"] for arm in arms} == {blind.record["latent"]["selected_lv"]}
    # One seed and one n: every arm draws the same permutation indices.
    n = blind.latent_scores.shape[0]
    first = np.random.default_rng(unit["seed"]).permutation(n)
    np.testing.assert_array_equal(np.random.default_rng(calls[-1]["seed"]).permutation(n), first)


def test_p_red_and_kappa_per_arm(diagnosis) -> None:
    for record in _records(diagnosis):
        arms = record["arms"]
        n = 80
        assert {arm: arms[arm]["p_red"] for arm in P_RED} == P_RED
        for arm, p_red in P_RED.items():
            assert arms[arm]["kappa"] is None
            assert arms[f"{arm}_corrected"]["p_red"] == p_red
            assert arms[f"{arm}_corrected"]["kappa"] == pytest.approx(math.sqrt(n / (n - p_red)), rel=1e-12)
        assert arms["z_noise5_corrected"]["kappa"] > arms["oracle_noise0_corrected"]["kappa"]
        assert arms["latent_stage_cv"]["p_red"] == 10 and arms["latent_stage_cv"]["kappa"] is None


def test_corrected_arms_keep_the_observed_statistics(diagnosis) -> None:
    for record in _records(diagnosis):
        arms = record["arms"]
        for arm in diag.FACTORIAL_ARMS:
            assert arms[f"{arm}_corrected"]["statistics"] == arms[arm]["statistics"]
            assert set(arms[f"{arm}_corrected"]["p_values"]) == set(STATISTICS)
        # The script's case-study cell is the analysis's (asserted per unit while running).
        assert arms["z_noise5"]["selected_lv"] == arms["oracle_noise0"]["selected_lv"]


def test_kappa_one_reproduces_the_uncorrected_p_values(unit) -> None:
    from motco.analysis import TrajectoryDataset

    for arm, (oracle, _) in diag.FACTORIAL.items():
        columns = diag.arm_covariates(arm)
        covariates = TrajectoryDataset(**unit["common"], covariate_cols=columns).covariate_design if columns else None
        latent = unit["Y0"] if oracle else unit["blind"].latent_scores
        kwargs = {"permutations": 49, "seed": unit["seed"], "orientation_estimator": "net_displacement"}
        uncorrected = diag.infer(latent, unit["design"], covariates, **kwargs)
        at_one = diag.infer(latent, unit["design"], covariates, kappa=1.0, **kwargs)
        assert at_one["statistics"] == uncorrected["statistics"]
        assert at_one["p_values"] == uncorrected["p_values"]
        assert at_one["null_q95"] == pytest.approx(uncorrected["null_q95"], rel=1e-9)


def test_rescaled_residuals_keep_the_reduced_fit() -> None:
    from motco.stats.trajectory import estimate_betas

    rng = np.random.default_rng(0)
    reduced = np.column_stack([np.ones(30), rng.normal(size=(30, 3))])
    Y = rng.normal(size=(30, 2))
    rescaled = diag.rescaled_residuals(Y, reduced, 1.5)
    fitted = reduced @ np.asarray(estimate_betas(reduced, Y))
    np.testing.assert_allclose(reduced @ np.asarray(estimate_betas(reduced, rescaled)), fitted, atol=1e-12)
    np.testing.assert_allclose(rescaled - fitted, 1.5 * (Y - fitted), atol=1e-12)


# ── Latent arms ──────────────────────────────────────────────────────────────


def test_latent_arms_record_their_rank_rules(diagnosis) -> None:
    for record in _records(diagnosis):
        arms = record["arms"]
        assert arms["latent_stage_rank2"]["component_selection"] == "forced"
        assert arms["latent_stage_rank2"]["selected_lv"] == 2
        assert arms["latent_stage_rank2"]["cv_mean_auroc"] is None
        assert arms["latent_stage_cv"]["component_selection"] == "cv"
        assert arms["latent_stage_cv"]["selected_lv"] >= 2
        assert arms["latent_stage_floor1"]["component_selection"] == "cv"
        for arm in diag.LATENT_ARMS:
            assert arms[arm]["error"] is None
            assert set(arms[arm]["p_values"]) == set(STATISTICS)


def test_floor1_arm_is_unfloored_and_records_rank_one(unit, monkeypatch) -> None:
    from motco import analysis

    real = analysis.analyze_trajectories
    seen = []

    def narrowed(dataset, params):
        seen.append(params)
        # With max_components 2 the unfloored CV has the single candidate rank 1.
        return real(dataset, replace(params, max_components=2))

    monkeypatch.setattr(analysis, "analyze_trajectories", narrowed)
    entry = diag.latent_arm(unit["six"], unit["settings"], unit["seed"], "latent_stage_floor1", 10)
    assert seen[0].resolved_min_components == 1 and seen[0].forced_components is None
    assert seen[0].covariate_latent_adjustment == "stage"
    assert entry["selected_lv"] == 1 and entry["component_selection"] == "cv"
    assert set(entry["p_values"]) == set(STATISTICS)  # angle at rank 1 is recorded
    record = {"key": "k", "mode": "anchor", "effect": 0.0, "replicate": 0,
              "arms": {**{arm: _synthetic_arm(0.5) for arm in diag.ARMS}, "latent_stage_floor1": entry}}
    table = diag.rank_rows({"k": record}, "run", alpha=0.05)
    rank1 = [row for row in table if row["arm"] == "latent_stage_floor1" and row["rank_bucket"] == "1"]
    assert {row["statistic"]: row["replicates"] for row in rank1} == dict.fromkeys(STATISTICS, 1)
    rates = diag.rate_rows(diag.p_value_rows({"k": record}, "run"), alpha=0.05, k=2.0)
    angle = [r for r in rates if (r["mode"], r["arm"], r["statistic"]) == ("anchor", "latent_stage_floor1", "angle")]
    assert angle[0]["replicates"] == 1


def test_failed_latent_arm_leaves_the_other_arms(unit, monkeypatch) -> None:
    from motco import analysis

    real = analysis.analyze_trajectories

    def raising(dataset, params):
        if params.forced_components == 2:
            raise analysis.TrajectoryAnalysisError("synthetic rank-2 failure")
        return real(dataset, params)

    monkeypatch.setattr(analysis, "analyze_trajectories", raising)
    record = diag.evaluate_unit("anchor", 0.0, 0, unit["settings"])
    failed = record["arms"]["latent_stage_rank2"]
    assert failed["error"] == "synthetic rank-2 failure"
    assert failed["p_values"] == {} and failed["statistics"] == {} and failed["selected_lv"] is None
    assert set(failed) == set(diag.ARM_FIELDS)
    for arm in diag.ARMS:
        if arm != "latent_stage_rank2":
            assert record["arms"][arm]["error"] is None
            assert set(record["arms"][arm]["p_values"]) == set(STATISTICS)
    rates = diag.rate_rows(diag.p_value_rows({record["key"]: record}, "run"), alpha=0.05, k=2.0)
    by_arm = {(r["mode"], r["arm"], r["statistic"]): r for r in rates}
    assert by_arm[("anchor", "latent_stage_rank2", "delta")]["verdict"] == "unavailable"
    assert by_arm[("anchor", "latent_stage_cv", "delta")]["replicates"] == 1


# ── Records, resume, merge ───────────────────────────────────────────────────


def test_torn_line_is_re_evaluated(smoke) -> None:
    settings = cs.load_settings(smoke)
    before = _records(smoke)
    path = cs.shard_path(smoke, 0, 1)
    lines = path.read_text(encoding="utf-8").splitlines()
    assert len(lines) == 2
    path.write_text(lines[0] + "\n" + lines[1][: len(lines[1]) // 2] + "\n", encoding="utf-8")
    assert diag.run_shard(smoke, settings, 0, 1) == 1
    assert _records(smoke) == before


def test_resume_skips_recorded_units_under_another_shard_count(diagnosis) -> None:
    settings = cs.load_settings(diagnosis)
    before = _records(diagnosis)
    assert diag.run_shard(diagnosis, settings, 0, 1) == 0
    assert diag.run_shard(diagnosis, settings, 1, 2) == 0
    assert cs.shard_path(diagnosis, 1, 2).read_text(encoding="utf-8") == ""
    assert _records(diagnosis) == before


def test_merge_writes_the_tables(diagnosis) -> None:
    summary = diag.merge(diagnosis, cs.load_settings(diagnosis))
    assert summary["units"] == 4 and summary["missing"] == []
    assert summary["failures"] == dict.fromkeys(diag.ARMS, 0)
    assert summary["reproduction"]["status"] == "holds"
    assert summary["max_reconstruction_error"] <= diag.RECONSTRUCTION_TOLERANCE
    pvalues = (diagnosis / "diagnosis_pvalues.csv").read_text(encoding="utf-8").splitlines()
    assert pvalues[0].split(",") == list(diag.PVALUE_COLUMNS)
    assert len(pvalues) == 1 + 4 * len(diag.ARMS) * 3
    rates = (diagnosis / "diagnosis_rates.csv").read_text(encoding="utf-8").splitlines()
    assert rates[0].split(",") == list(diag.RATE_COLUMNS)
    assert len(rates) == 1 + 2 * len(diag.ARMS) * 3
    ranks = (diagnosis / "diagnosis_rank_table.csv").read_text(encoding="utf-8").splitlines()
    assert ranks[0].split(",") == list(diag.RANK_COLUMNS)
    assert len(ranks) == 1 + len(diag.RANK_TABLE_ARMS) * len(diag.RANK_BUCKETS) * 3


def test_smoke_run_is_not_compared(smoke) -> None:
    report = json.loads((smoke / "diagnosis_reproduction.json").read_text(encoding="utf-8"))
    assert report["status"] == "not_compared"
    assert report["permutations"] == 5 and report["parent_permutations"] == 9
    assert all(entry["compared"] == 0 for entry in report["arms"].values())


def test_decide_on_merged_runs(diagnosis, smoke, tmp_path) -> None:
    decision = diag.decide([diagnosis], tmp_path / "full")
    assert decision["decision"] == "diagnosed" and decision["indicative"] is False
    assert decision["verdicts"]["adjusted_delta"]["verdict"] in diag.ADJUSTED_VERDICTS
    assert decision["verdicts"]["latent_stage"]["verdict"] in diag.LATENT_VERDICTS
    settings = cs.load_settings(diagnosis)
    assert decision["run_settings"][0]["settings_hash"] == json.loads(
        (diagnosis / "settings.json").read_text(encoding="utf-8"))["hash"] == cs.settings_hash(settings)
    assert decision["run_settings"][0]["parent_hash"] == settings["parent"]["hash"]
    written = json.loads((tmp_path / "full" / "case_study_diagnosis_decision.json").read_text(encoding="utf-8"))
    assert written["decision"] == "diagnosed"
    contrasts = (tmp_path / "full" / "diagnosis_contrasts.csv").read_text(encoding="utf-8").splitlines()
    assert contrasts[0].split(",") == list(diag.CONTRAST_COLUMNS)
    assert len(contrasts) == 1 + 3 * 4 * (1 + 2)  # statistics × contrasts × (clustered + two readings)
    indicative = diag.decide([smoke], tmp_path / "smoke")
    assert indicative["indicative"] is True
    assert indicative["verdicts"]["adjusted_delta"]["indicative"] is True
    with pytest.raises(SystemExit, match="distinct parents"):
        diag.decide([diagnosis, smoke], tmp_path / "both")


# ── Tables on synthetic records ──────────────────────────────────────────────


def _synthetic_arm(p_value: float, lv: int = 2, error: str | None = None) -> dict:
    if error is not None:
        return {**dict.fromkeys(diag.ARM_FIELDS), "statistics": {}, "p_values": {}, "null_q95": {}, "error": error}
    return {**dict.fromkeys(diag.ARM_FIELDS), "p_values": dict.fromkeys(STATISTICS, p_value),
            "statistics": dict.fromkeys(STATISTICS, 1.0), "selected_lv": lv, "error": None}


def _synthetic_rows(rejecting: dict, *, replicates: int = 100, runs=("delta0.1", "delta0.05")) -> list[dict]:
    """``{arm: {statistic: set of rejecting replicates}}`` → p-value rows; other arms reject replicates 0–4."""

    default = set(range(5))
    rows = []
    for run, (mode, effect), replicate, arm, statistic in itertools.product(
        runs, (("anchor", 0.0), ("translation", 1.0)), range(replicates), diag.ARMS, STATISTICS
    ):
        rejected = replicate in rejecting.get(arm, {}).get(statistic, default)
        rows.append({"run": run, "mode": mode, "effect": effect, "replicate": replicate, "arm": arm,
                     "statistic": statistic, "p_value": 0.01 if rejected else 0.5})
    return rows


def test_bound_and_verdict_boundary_at_500_replicates() -> None:
    rows = []
    for replicate in range(500):
        for arm, rejections in (("z_noise5", 34), ("z_noise0", 35)):
            rows.append({"run": "r", "mode": "anchor", "effect": 0.0, "replicate": replicate, "arm": arm,
                         "statistic": "delta", "p_value": 0.01 if replicate < rejections else 0.5})
    rates = {r["arm"]: r for r in diag.rate_rows(rows, alpha=0.05, k=2.0)
             if r["mode"] == "anchor" and r["statistic"] == "delta" and r["arm"] in ("z_noise5", "z_noise0")}
    assert rates["z_noise5"]["bound"] == pytest.approx(0.0695, abs=1e-4)
    assert rates["z_noise5"]["replicates"] == 500
    assert (rates["z_noise5"]["rejections"], rates["z_noise5"]["verdict"]) == (34, "holds")
    assert (rates["z_noise0"]["rejections"], rates["z_noise0"]["verdict"]) == (35, "exceeds")


def test_every_anchor_unit_falls_in_one_rank_bucket() -> None:
    records = {}
    for replicate, lv in enumerate([1, 2, 2, 3, 4, 7, 11, None]):
        for mode, effect in (("anchor", 0.0), ("translation", 1.0)):
            arms = {arm: _synthetic_arm(0.01 if replicate % 2 else 0.5, lv=lv or 2) for arm in diag.ARMS}
            if lv is None:
                arms["latent_stage_cv"] = _synthetic_arm(0.5, error="failed")
            key = cs.unit_key(mode, effect, replicate)
            records[key] = {"key": key, "mode": mode, "effect": effect, "replicate": replicate, "arms": arms}
    table = diag.rank_rows(records, "run", alpha=0.05)
    for arm in diag.RANK_TABLE_ARMS:
        delta = {row["rank_bucket"]: row for row in table if row["arm"] == arm and row["statistic"] == "delta"}
        assert set(delta) == set(diag.RANK_BUCKETS)
        available = 7 if arm == "latent_stage_cv" else 8  # the failed unit is unavailable; translation never enters
        assert sum(row["units"] for row in delta.values()) == available
    cv = {row["rank_bucket"]: row for row in table if row["arm"] == "latent_stage_cv" and row["statistic"] == "angle"}
    assert {bucket: row["units"] for bucket, row in cv.items()} == {"1": 1, "2": 2, "3": 1, ">=4": 3}
    assert cv[">=4"]["rejections"] == 1 and cv[">=4"]["rate"] == pytest.approx(1 / 3)  # replicates 4, 5, 6
    assert diag.rank_bucket(None) is None and diag.rank_bucket("4") == ">=4"


def test_altered_parent_p_value_is_a_mismatch(parent, diagnosis) -> None:
    settings = cs.load_settings(diagnosis)
    records = diag.read_records([cs.shard_path(diagnosis, 0, 1)])
    parent_p_values = diag.read_parent_p_values(parent)
    assert diag.reproduction_report(records, settings, parent_p_values)["status"] == "holds"
    unit = ("translation", 1.0, 1, "latent_stage", "shape")
    altered = {**parent_p_values, unit: 1.0 if parent_p_values[unit] != 1.0 else 0.1}
    report = diag.reproduction_report(records, settings, altered)
    assert report["status"] == "failed"
    assert report["arms"]["latent_stage_cv"]["mismatches"] == 1 and report["arms"]["z_noise5"]["mismatches"] == 0
    assert report["mismatching_units"] == [cs.unit_key("translation", 1.0, 1)]
    smoke = diag.reproduction_report(records, {**settings, "smoke": True}, None)
    assert smoke["status"] == "not_compared" and smoke["mismatching_units"] == []


# ── Contrasts and verdicts ───────────────────────────────────────────────────


def test_sign_flip_matches_brute_force() -> None:
    for differences in ([2, -1, 0, 3, 1, -2, 4], [1, 1, 1], [-1, -2, 0], [0, 0], [4, -4, 3, 3, -1, 2, 2, 0, 1]):
        observed = sum(differences)
        signs = list(itertools.product((-1, 1), repeat=len(differences)))
        brute = sum(sum(s * d for s, d in zip(sign, differences)) >= observed for sign in signs) / len(signs)
        assert diag.sign_flip_p_value(differences) == pytest.approx(brute, abs=1e-12)
    assert diag.sign_flip_p_value([]) == 1.0


def test_one_replicate_contributes_one_cluster() -> None:
    rows = _synthetic_rows({"z_noise5": {"delta": {0}}, "z_noise0": {"delta": set()}}, replicates=3)
    rejected = diag._rejections(rows, "delta", alpha=0.05)
    result = diag.paired_contrast(rejected, "z_noise5", "z_noise0")
    assert (result["readings"], result["clusters"]) == (12, 3)
    assert (result["b_only"], result["a_only"], result["net"]) == (4, 0, 4)
    # Four discordant readings of one replicate flip together: P(±4 >= 4) = 1/2, not 1/16.
    assert result["p_value"] == pytest.approx(0.5)
    one = diag.paired_contrast(rejected, "z_noise5", "z_noise0", readings={("delta0.1", "anchor", 0.0)})
    assert (one["readings"], one["b_only"], one["p_value"]) == (3, 1, pytest.approx(0.5))


def _reproductions(status: str = "holds", mismatching: tuple = ()) -> list[dict]:
    return [
        {"parent": {"label": "delta0.1"}, "status": status, "arms": {}, "mismatching_units": list(mismatching)},
        {"parent": {"label": "delta0.05"}, "status": "holds" if status == "failed" else status, "arms": {},
         "mismatching_units": []},
    ]


def _decide(rejecting: dict, reproductions: list[dict] | None = None) -> dict:
    return diag.diagnosis_decision(
        _synthetic_rows(rejecting), [], reproductions or _reproductions(), alpha=0.05, k=2.0
    )


HIGH, BASE = set(range(15)), set(range(5))


def test_shrinkage_verdict() -> None:
    decision = _decide({"z_noise5": {"delta": HIGH}})
    verdict = decision["verdicts"]["adjusted_delta"]
    assert decision["decision"] == "diagnosed" and decision["indicative"] is False
    assert verdict["verdict"] == "residual_df_shrinkage"
    assert verdict["conditions"]["column_contrast_significant"] == {"C1": False, "C2": True}
    assert verdict["conditions"]["corrected_case_study_model_holds"] is True
    assert verdict["conditions"]["p_values"]["C2"] == pytest.approx(2.0 ** -10)
    assert "corrected permutation" in verdict["consequence"]
    clustered = {row["contrast"]: row for row in decision["contrasts"] if row["statistic"] == "delta"}
    assert (clustered["C2"]["b_only"], clustered["C2"]["a_only"], clustered["C2"]["readings"]) == (40, 0, 400)


def test_collinear_verdict() -> None:
    # Estimating z costs size, and the column contrast is flat.
    flat = _decide({"z_noise5": {"delta": HIGH}, "z_noise0": {"delta": HIGH}})
    assert flat["verdicts"]["adjusted_delta"]["verdict"] == "collinear_covariate"
    # Shrinkage fails because the corrected case-study model still exceeds; the collinear contrast decides.
    exceeding = _decide({
        "z_noise5": {"delta": set(range(20))},
        "z_noise0": {"delta": set(range(10))},
        "z_noise5_corrected": {"delta": HIGH},
    })
    verdict = exceeding["verdicts"]["adjusted_delta"]
    assert verdict["conditions"]["column_contrast_significant"]["C2"] is True
    assert verdict["conditions"]["corrected_case_study_model_holds"] is False
    assert {(row["run"], row["mode"], row["statistic"]) for row in verdict["conditions"]["corrected_not_holding"]} == {
        (run, mode, "delta") for run in ("delta0.1", "delta0.05") for mode in ("anchor", "translation")
    }
    assert verdict["verdict"] == "collinear_covariate"


def test_unresolved_verdict_names_the_replication() -> None:
    decision = _decide({})
    verdict = decision["verdicts"]["adjusted_delta"]
    assert verdict["verdict"] == "unresolved"
    assert "R = 1,000" in verdict["consequence"] and "R = 1,000" in verdict["next_step"]
    # A significant corrected column contrast also blocks the shrinkage verdict.
    corrected_trend = _decide({"z_noise5": {"delta": HIGH}, "z_noise5_corrected": {"delta": set(range(9))},
                               "z_noise0_corrected": {"delta": set()}})
    assert corrected_trend["verdicts"]["adjusted_delta"]["conditions"]["corrected_column_contrast_significant"]
    assert corrected_trend["verdicts"]["adjusted_delta"]["verdict"] == "unresolved"


def test_latent_verdicts() -> None:
    selected = _decide({"latent_stage_cv": {"angle": HIGH, "shape": HIGH}})
    latent = selected["verdicts"]["latent_stage"]
    assert latent["verdict"] == "rank_selection" and latent["rank2_not_holding"] == []
    assert {(row["mode"], row["statistic"]) for row in latent["cv_exceeding"]} == {
        (mode, statistic) for mode in ("anchor", "translation") for statistic in ("angle", "shape")
    }
    assert "fixed rank of 2" in latent["consequence"]
    basis = _decide({"latent_stage_rank2": {"angle": HIGH}})
    latent = basis["verdicts"]["latent_stage"]
    assert latent["verdict"] == "residualized_basis"
    assert {(row["run"], row["mode"], row["statistic"]) for row in latent["rank2_not_holding"]} == {
        (run, mode, "angle") for run in ("delta0.1", "delta0.05") for mode in ("anchor", "translation")
    }
    # The floor-1 arm never enters the verdict.
    floor1 = _decide({"latent_stage_floor1": {"angle": HIGH, "shape": HIGH}})
    assert floor1["verdicts"]["latent_stage"]["verdict"] == "rank_selection"


def test_reproduction_failure_blocks_every_verdict() -> None:
    decision = _decide({"z_noise5": {"delta": HIGH}}, _reproductions("failed", ('["anchor", 0.0, 3]',)))
    assert decision["decision"] == "reproduction_failed"
    assert decision["verdicts"] is None
    assert decision["reproduction"]["status"] == "failed"
    assert decision["reproduction"]["mismatching_units"] == ['["anchor", 0.0, 3]']


def test_smoke_verdicts_are_indicative() -> None:
    decision = _decide({"z_noise5": {"delta": HIGH}}, _reproductions("not_compared"))
    assert decision["decision"] == "diagnosed" and decision["indicative"] is True
    assert decision["reproduction"]["status"] == "not_compared"
    assert decision["verdicts"]["adjusted_delta"]["indicative"] is True
    assert decision["verdicts"]["latent_stage"]["indicative"] is True


def test_floor1_summary_reports_the_rank1_fraction() -> None:
    table = [
        {"run": "delta0.1", "arm": "latent_stage_floor1", "rank_bucket": bucket, "statistic": "delta", "units": units}
        for bucket, units in (("1", "3"), ("2", "5"), ("3", "2"), (">=4", "0"))
    ]
    summary = diag.floor1_summary(table)
    assert summary["delta0.1"]["units"] == 10
    assert summary["delta0.1"]["rank1_fraction"] == pytest.approx(0.3)
