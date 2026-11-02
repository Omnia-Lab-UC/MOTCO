"""``scripts/adni_design_check.py``: the ADNI-design construction, the units, merge, calibration and the decision."""

from __future__ import annotations

import csv
import json
import math
import sys
from pathlib import Path
from typing import Any

import numpy as np
import pytest


def _script(name: str):
    sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
    try:
        return __import__(name)
    finally:
        sys.path.pop(0)


check = _script("adni_design_check")
cs = _script("case_study_type_i_check")
d10 = _script("covariate_type_i_check")
SMALL_CV = {"cv1_splits": 2, "cv2_splits": 3, "n_repeats": 2, "max_components": 5, "random_state": 1203}
SMALL = "test36"
#: A small registered table (n = 36) so the run path stays fast; the shares are the table's own.
check.TABLES.setdefault(
    SMALL,
    {
        "group_stage_sizes": [[6, 6, 6], [6, 6, 6]],
        "phase_shares": [[0.6, 0.6, 0.4], [0.55, 0.5, 0.5]],
        "source": {"table": "tests/test_adni_design_check.py"},
        "sha256": {},
        "reference": None,
    },
)
ADNI = [[237, 310, 125], [223, 448, 158]]
SIX = ["z", "category[c1]", "category[c2]", "phase", "x1", "x2"]
#: Base seeds of the committed checks and studies (design D7).
COMMITTED_BASE_SEEDS = (100, 400, 500, 600, 800, 900, 1028, 2026, 2031)


def _args(*argv: str, command: str = "run", out_dir: str = "unused"):
    return check.build_parser().parse_args([command, *argv, "--out-dir", out_dir])


def _settings(*argv: str, reps: int = 2, perms: int = 9, table: str = SMALL) -> dict:
    settings = check.build_settings(_args("--table", table, "--reps", str(reps), "--perms", str(perms), *argv))
    settings["cv"] = dict(SMALL_CV)
    return settings


def _few_cells(settings: dict) -> dict:
    settings["effects"] = {"magnitude": [0.5], "orientation": [0.5], "shape": [0.5]}
    return settings


# ── Settings ─────────────────────────────────────────────────────────────────


def test_adni_settings_carry_the_design() -> None:
    settings = check.build_settings(_args())
    assert settings["check"] == "adni_design" and settings["base_seed"] == 2032
    assert settings["generator"]["group_stage_sizes"] == ADNI
    assert sum(map(sum, ADNI)) == 1501
    assert settings["table"]["name"] == "adni" and len(settings["table"]["sha256"]["result.json"]) == 64
    shares = settings["covariates"]["phase"]["shares"]
    np.testing.assert_allclose(shares, [[0.603, 0.623, 0.424], [0.561, 0.507, 0.481]], atol=5e-4)
    assert settings["covariates"]["declared"] == ["z", "category", "phase", "x1", "x2"]
    assert settings["covariates"]["categorical_declared"] == ["category"]
    assert settings["covariates"]["noise_numeric"]["count"] == 2
    assert {key: settings["injection"][key] for key in ("s_z", "s_c", "s_p")} == {"s_z": 0.1, "s_c": 0.0, "s_p": 1.75}
    assert settings["injection"]["phase"]["block"] == "expression"
    assert {settings["generator"][f] for f in ("delta_methyl", "delta_expr", "delta_protein")} == {0.02}
    assert settings["generator"]["p_dmp"] == 0.1 and settings["layers"] == ["methylation", "expression"]
    assert settings["min_components"] == 2 and settings["orientation_estimator"] == "net_displacement"
    assert settings["replicates"] == 500 and settings["permutations"] == 999
    assert settings["alpha"] == 0.05 and settings["k"] == 2.0
    assert settings["rule"]["power_floor"] == 0.8 and settings["rule"]["k"] == 2.0
    assert settings["calibration"]["auroc_band"] == [0.55, 0.65] and settings["calibration"]["rho_threshold"] == 0.4
    assert set(settings["arms"]) == {"adjusted", "unadjusted"}


def test_adni_run_enumerates_fourteen_cells_replicate_major() -> None:
    settings = check.build_settings(_args())
    cells = check.cells(settings)
    assert len(cells) == 14
    assert cells[:2] == [("anchor", 0.0), ("translation", 1.0)]
    assert [e for m, e in cells if m == "magnitude"] == [0.1, 0.25, 0.5, 1.0]
    assert [e for m, e in cells if m == "orientation"] == [0.25, 0.5, 0.75, 1.0]
    assert [e for m, e in cells if m == "shape"] == [0.25, 0.5, 0.75, 1.0]
    units = check.units(settings)
    assert len(units) == 14 * 500
    assert [replicate for _, _, replicate in units[:14]] == [0] * 14
    assert [(mode, effect) for mode, effect, _ in units[14:28]] == cells


@pytest.mark.parametrize(
    "argv",
    [
        ("--table", SMALL),
        ("--delta", "0.05"),
        ("--s-z", "1.0"),
        ("--s-p", "2.0"),
        ("--base-seed", "2033"),
        ("--perms", "99"),
    ],
)
def test_settings_hash_changes(argv) -> None:
    base = check.build_settings(_args())
    assert check.settings_hash(check.build_settings(_args(*argv))) != check.settings_hash(base)


def test_settings_hash_changes_with_the_cells() -> None:
    base = check.build_settings(_args())
    changed = check.build_settings(_args())
    changed["effects"]["magnitude"] = [0.1, 0.25, 0.5]
    assert check.settings_hash(changed) != check.settings_hash(base)


@pytest.mark.parametrize("argv", [("--s-p", "2.0"), ("--delta", "0.05"), ("--perms", "19")])
def test_resume_under_other_settings_is_refused(tmp_path, argv) -> None:
    check.ensure_settings(tmp_path, check.build_settings(_args("--table", SMALL)))
    with pytest.raises(SystemExit, match="different settings"):
        check.ensure_settings(tmp_path, check.build_settings(_args("--table", SMALL, *argv)))
    assert not (tmp_path / "records").exists()


# ── Construction ─────────────────────────────────────────────────────────────


def test_phase_offset_lands_on_the_expression_block_only() -> None:
    settings = _settings()
    null = check.build_unit("anchor", 0.0, 0, {**settings, "injection": {**settings["injection"], "s_p": 0.0}})
    built = check.build_unit("anchor", 0.0, 0, settings)
    np.testing.assert_array_equal(null["blocks"]["methylation"].to_numpy(), built["blocks"]["methylation"].to_numpy())
    np.testing.assert_array_equal(null["metadata"]["phase"], built["metadata"]["phase"])
    phase = built["metadata"]["phase"].to_numpy()
    assert set(np.unique(phase)) == {0.0, 1.0}
    difference = built["blocks"]["expression"].to_numpy() - null["blocks"]["expression"].to_numpy()
    np.testing.assert_allclose(difference, np.outer(phase, built["phase_offset"]), atol=1e-12)
    assert not np.allclose(built["phase_offset"], 0.0) and np.allclose(null["phase_offset"], 0.0)
    realized = built["phase_realized"]
    assert set(realized) == {"methylation", "expression"}
    assert realized["expression"]["n_features"] == 131 and realized["methylation"]["n_features"] == 367
    assert realized["expression"]["median_abs_d"] > null["phase_realized"]["expression"]["median_abs_d"]
    assert 0.0 <= realized["expression"]["abs_corr_pc1_phase"] <= 1.0


def test_build_draws_phase_by_the_table_cells() -> None:
    settings = _settings()
    built = check.build_unit("anchor", 0.0, 1, settings)
    metadata = built["metadata"]
    group, stage = check.table_cells(settings["generator"]["group_stage_sizes"])
    labels = ["A", "B"]
    np.testing.assert_array_equal([labels.index(g) for g in metadata["group"]], group)
    np.testing.assert_array_equal(metadata["stage"].to_numpy(), stage)
    rng = np.random.default_rng(check.phase_seed(2032, 1))
    expected = check.draw_phase(group, stage, settings["covariates"]["phase"]["shares"], rng)
    np.testing.assert_array_equal(metadata["phase"].to_numpy(), expected)
    assert built["seeds"]["phase"] == check.phase_seed(2032, 1)


def test_phase_draw_matches_the_adni_shares() -> None:
    settings = check.build_settings(_args())
    shares = np.asarray(settings["covariates"]["phase"]["shares"])
    group, stage = check.table_cells(ADNI)
    totals = np.zeros_like(shares)
    replicates = 2000
    for replicate in range(replicates):
        phase = check.draw_phase(group, stage, shares, np.random.default_rng(check.phase_seed(2032, replicate)))
        for g in range(2):
            for s in range(3):
                totals[g, s] += phase[(group == g) & (stage == s)].mean()
    np.testing.assert_allclose(totals / replicates, shares, atol=0.03)


def test_cells_of_a_replicate_share_their_seeds_and_covariates() -> None:
    settings = _settings()
    for replicate in (0, 1):
        anchor = check.build_unit("anchor", 0.0, replicate, settings)
        other = check.build_unit("orientation", 0.5, replicate, settings)
        assert anchor["seeds"] == other["seeds"]
        assert set(anchor["seeds"]) == {"generator", "covariates", "rrpp", "noise_covariates", "phase"}
        for name in ("z", "category", "phase", "x1", "x2"):
            np.testing.assert_array_equal(anchor["metadata"][name], other["metadata"][name])
        np.testing.assert_array_equal(anchor["phase_offset"] != 0, other["phase_offset"] != 0)


def test_seeds_are_fresh() -> None:
    settings = _settings()
    for replicate in (0, 1, 2):
        seeds = check.build_unit("anchor", 0.0, replicate, settings)["seeds"]
        assert {key: seeds[key] for key in ("generator", "covariates", "rrpp")} == d10.replicate_seeds(2032, replicate)
        committed: set[int] = set()
        for base in COMMITTED_BASE_SEEDS:
            committed |= set(d10.replicate_seeds(base, replicate).values())
            committed.add(cs._seed(base, "noise_covariates", replicate))
        assert not set(seeds.values()) & committed
    direction = check.phase_direction(131, settings)
    assert not np.allclose(direction, check.phase_direction(131, {**settings, "base_seed": 2031}))


# ── Units ────────────────────────────────────────────────────────────────────


@pytest.fixture(scope="module")
def unit_record():
    return check.evaluate_unit("anchor", 0.0, 0, _settings())


def test_unit_record_holds_every_field(unit_record) -> None:
    record = unit_record
    assert record["key"] == check.unit_key("anchor", 0.0, 0)
    assert record["settings_hash"] == check.settings_hash(_settings())
    assert set(record["seeds"]) == {"generator", "covariates", "rrpp", "noise_covariates", "phase"}
    assert set(record["phase_realized"]) == {"methylation", "expression"}
    assert len(record["phase_cell_shares"]) == 2 and len(record["phase_cell_shares"][0]) == 3
    assert record["runtime_seconds"] > 0
    assert set(record["arms"]) == {"adjusted"}
    arm = record["arms"]["adjusted"]
    assert set(arm) == set(check.ARM_FIELDS)
    assert arm["error"] is None
    assert set(arm["statistics"]) == set(arm["p_values"]) == {"delta", "angle", "shape"}
    assert set(arm["null_summary"]) == {"delta", "angle", "shape"}
    assert arm["covariate_design_columns"] == SIX
    assert arm["selected_lv"] >= 2 and arm["component_selection"] == "cv"
    assert 0.0 <= arm["cv_mean_auroc"] <= 1.0
    assert abs(arm["first_component_spearman_phase"]) <= 1 and abs(arm["first_component_spearman_z"]) <= 1
    assert 0.0 <= arm["pooled_relative_eigengap"] <= 1.0
    assert len(arm["latent_digest"]) == 64


def test_units_are_reproducible(unit_record) -> None:
    again = check.evaluate_unit("anchor", 0.0, 0, _settings())
    drop = {"runtime_seconds"}
    assert {k: v for k, v in again.items() if k not in drop} == {k: v for k, v in unit_record.items() if k not in drop}


def test_analysis_error_is_recorded_without_p_values(monkeypatch) -> None:
    def refuse(dataset, params):
        raise check.TrajectoryAnalysisError("the covariate design is rank deficient")

    monkeypatch.setattr(check, "analyze_trajectories", refuse)
    record = check.evaluate_unit("anchor", 0.0, 0, _settings())
    arm = record["arms"]["adjusted"]
    assert arm["error"] == "the covariate design is rank deficient"
    assert arm["p_values"] == {} and arm["statistics"] == {} and arm["selected_lv"] is None
    rows = check.p_value_rows({record["key"]: record})
    assert [row["p_value"] for row in rows] == [None, None, None]
    assert check.rejection_tallies(rows, alpha=0.05) == {}


def test_calibration_arms_share_the_latent_space() -> None:
    record = check.evaluate_unit("anchor", 0.0, 0, _settings(), arms=check.CALIBRATION_ARMS)
    adjusted, unadjusted = record["arms"]["adjusted"], record["arms"]["unadjusted"]
    assert adjusted["latent_digest"] == unadjusted["latent_digest"]
    assert adjusted["selected_lv"] == unadjusted["selected_lv"]
    assert adjusted["covariate_design_columns"] == SIX and unadjusted["covariate_design_columns"] is None


def test_calibration_arms_that_differ_raise(monkeypatch) -> None:
    real = check.analyze_trajectories

    def perturbed(dataset, params):
        result = real(dataset, params)
        if not dataset.covariate_cols:
            result.latent_scores.iloc[0, 0] += 1e-12
        return result

    monkeypatch.setattr(check, "analyze_trajectories", perturbed)
    with pytest.raises(RuntimeError, match="latent spaces differ"):
        check.evaluate_unit("anchor", 0.0, 0, _settings(), arms=check.CALIBRATION_ARMS)


def _records(out_dir: Path) -> list[dict]:
    lines = [
        line
        for path in sorted((out_dir / "records").glob("shard_*.jsonl"))
        for line in path.read_text(encoding="utf-8").splitlines()
    ]
    return sorted((json.loads(line) for line in lines), key=lambda record: record["key"])


@pytest.fixture(scope="module")
def sharded_run(tmp_path_factory):
    out_dir = tmp_path_factory.mktemp("adni-design")
    settings = _few_cells(_settings())
    check.ensure_settings(out_dir, settings)
    evaluated = [check.run_shard(out_dir, settings, index, 2) for index in (0, 1)]
    return out_dir, settings, evaluated


def test_two_shards_cover_every_unit_once(sharded_run) -> None:
    out_dir, settings, evaluated = sharded_run
    assert evaluated == [5, 5]
    keys = [record["key"] for record in _records(out_dir)]
    assert sorted(keys) == sorted(check.unit_key(*unit) for unit in check.units(settings))
    assert len(keys) == len(set(keys)) == 10


def test_resume_skips_recorded_units(sharded_run) -> None:
    out_dir, settings, _ = sharded_run
    before = _records(out_dir)
    assert check.run_shard(out_dir, settings, 0, 1) == 0
    assert check.run_shard(out_dir, settings, 1, 2) == 0
    assert _records(out_dir) == before


def test_torn_line_is_rerun(tmp_path, sharded_run) -> None:
    out_dir, settings, _ = sharded_run
    lines = check.shard_path(out_dir, 0, 2).read_text(encoding="utf-8").splitlines()
    copy = tmp_path / "run"
    check.ensure_settings(copy, settings)
    torn = check.shard_path(copy, 0, 2)
    torn.parent.mkdir(parents=True)
    torn.write_text("\n".join(lines[:-1]) + "\n" + lines[-1][: len(lines[-1]) // 2], encoding="utf-8")
    assert check.run_shard(copy, settings, 0, 2) == 1
    assert len(_records(copy)) == 5


def test_over_headroom_cell_refuses_to_start(tmp_path) -> None:
    settings = _settings()
    settings["effects"]["orientation"] = [0.25, 3.0]
    with pytest.raises(SystemExit, match=r'\["orientation", 3.0\]'):
        check.run_shard(tmp_path, settings, 0, 1)
    assert check.shard_path(tmp_path, 0, 1).read_text() == ""


def test_merge_writes_the_tables_and_is_idempotent(sharded_run) -> None:
    out_dir, settings, _ = sharded_run
    names = ("adni_design_pvalues.csv", "adni_design_replicates.csv", "adni_design_rates.csv",
             "adni_design_rank_strata.csv")
    first = check.merge(out_dir, settings)
    contents = {name: (out_dir / name).read_text(encoding="utf-8") for name in names}
    second = check.merge(out_dir, settings)
    assert first == second
    for name, text in contents.items():
        assert (out_dir / name).read_text(encoding="utf-8") == text
    assert first["units"] == first["expected"] == 10 and first["missing"] == [] and first["failures"] == 0
    assert len(contents["adni_design_pvalues.csv"].splitlines()) == 1 + 10 * 3
    assert len(contents["adni_design_replicates.csv"].splitlines()) == 1 + 10
    assert len(contents["adni_design_rates.csv"].splitlines()) == 1 + 5 * 3


def test_merge_refuses_records_from_other_settings(tmp_path, sharded_run) -> None:
    out_dir, settings, _ = sharded_run
    copy = tmp_path / "run"
    (copy / "records").mkdir(parents=True)
    path = check.shard_path(out_dir, 0, 2)
    (copy / "records" / path.name).write_text(path.read_text(encoding="utf-8"), encoding="utf-8")
    other = {**settings, "permutations": 19}
    with pytest.raises(SystemExit, match="other settings"):
        check.merge(copy, other)


# ── Merge on a fixed record set ──────────────────────────────────────────────


def _fake_record(settings: dict, mode: str, effect: float, replicate: int, p: dict, rank: int | None) -> dict:
    arm: dict[str, Any] = {field: None for field in check.ARM_FIELDS}
    arm.update(
        p_values=p,
        statistics={s: 1.0 for s in p},
        null_summary={"angle": {"q95": 30.0, "sd": 7.0}},
        selected_lv=rank,
        component_selection="cv",
        cv_mean_auroc=0.6,
        first_component_spearman_phase=0.5,
        first_component_spearman_z=0.1,
        pooled_relative_eigengap=0.7,
    )
    if rank is None:
        arm.update(p_values={}, statistics={}, null_summary={}, component_selection=None, error="failed")
    return {
        "key": check.unit_key(mode, effect, replicate),
        "mode": mode,
        "effect": float(effect),
        "replicate": replicate,
        "settings_hash": check.settings_hash(settings),
        "phase_realized": {
            "expression": {"median_abs_d": 1.0, "abs_corr_pc1_phase": 0.9},
            "methylation": {"median_abs_d": 0.03, "abs_corr_pc1_phase": 0.02},
        },
        "arms": {"adjusted": arm},
        "runtime_seconds": 1.0,
    }


def test_merge_rates_match_hand_counts(tmp_path) -> None:
    settings = _few_cells(_settings(reps=6))
    # anchor: delta rejects at replicates 0 and 1; translation: nothing; magnitude: delta at 0-3; replicate 5 fails.
    plan = {
        ("anchor", 0.0): lambda r: {"delta": 0.01 if r < 2 else 0.5, "angle": 0.04 if r == 3 else 0.6, "shape": 0.9},
        ("translation", 1.0): lambda r: {"delta": 0.3, "angle": 0.3, "shape": 0.3},
        ("magnitude", 0.5): lambda r: {"delta": 0.001 if r < 4 else 0.2, "angle": 0.05, "shape": 0.051},
        ("orientation", 0.5): lambda r: {"delta": 0.5, "angle": 0.01, "shape": 0.5},
        ("shape", 0.5): lambda r: {"delta": 0.5, "angle": 0.5, "shape": 0.02 if r % 2 else 0.5},
    }
    ranks = [2, 3, 5, 8, 12, None]
    records = [
        _fake_record(settings, mode, effect, replicate, plan[(mode, effect)](replicate), ranks[replicate])
        for replicate in range(6)
        for mode, effect in check.cells(settings)
    ]
    check.ensure_settings(tmp_path, settings)
    path = check.shard_path(tmp_path, 0, 1)
    path.parent.mkdir(parents=True)
    path.write_text("".join(json.dumps(record) + "\n" for record in records), encoding="utf-8")
    summary = check.merge(tmp_path, settings)
    assert summary["units"] == 30 and summary["failures"] == 5
    with (tmp_path / "adni_design_rates.csv").open(encoding="utf-8") as handle:
        rates = {(row["mode"], float(row["effect"]), row["statistic"]): row for row in csv.DictReader(handle)}
    expected = {
        ("anchor", 0.0, "delta"): (5, 2),
        ("anchor", 0.0, "angle"): (5, 1),
        ("anchor", 0.0, "shape"): (5, 0),
        ("translation", 1.0, "delta"): (5, 0),
        ("magnitude", 0.5, "delta"): (5, 4),
        ("magnitude", 0.5, "angle"): (5, 5),
        ("magnitude", 0.5, "shape"): (5, 0),
        ("orientation", 0.5, "angle"): (5, 5),
        ("shape", 0.5, "shape"): (5, 2),
    }
    for key, (replicates, rejections) in expected.items():
        row = rates[key]
        assert (int(row["replicates"]), int(row["rejections"])) == (replicates, rejections), key
        assert float(row["rate"]) == pytest.approx(rejections / replicates)
        assert int(row["units"]) == 6 and int(row["failures"]) == 1
    assert rates[("magnitude", 0.5, "delta")]["role"] == "target"
    assert rates[("magnitude", 0.5, "angle")]["role"] == "off_target"
    assert rates[("anchor", 0.0, "delta")]["role"] == "size"
    assert json.loads(rates[("anchor", 0.0, "delta")]["selected_rank"]) == {"2": 1, "3": 1, "5": 1, "8": 1, "12": 1}
    with (tmp_path / "adni_design_rank_strata.csv").open(encoding="utf-8") as handle:
        strata = [row for row in csv.DictReader(handle) if row["mode"] == "magnitude" and row["statistic"] == "delta"]
    by_stratum = {row["stratum"]: (int(row["replicates"]), int(row["rejections"])) for row in strata}
    assert by_stratum == {"2-3": (2, 2), "4-7": (1, 1), "8+": (2, 1)}
    assert {int(row["pooled_replicates"]) for row in strata} == {5}


# ── Calibration ──────────────────────────────────────────────────────────────


def _calibration_record(replicate: int, auroc: float, rho: float, unadjusted_p: float) -> dict:
    return {
        "mode": "anchor",
        "replicate": replicate,
        "runtime_seconds": 1.0,
        "phase_realized": {"expression": {"median_abs_d": 1.0, "abs_corr_pc1_phase": 0.9, "abs_d_gt_1": 60}},
        "arms": {
            "adjusted": {
                "cv_mean_auroc": auroc,
                "first_component_spearman_phase": rho,
                "first_component_spearman_z": 0.2,
                "pooled_relative_eigengap": 0.7,
                "selected_lv": 4,
                "null_summary": {"angle": {"q95": 30.0}},
                "p_values": {"delta": 0.5, "angle": 0.5, "shape": 0.5},
                "error": None,
            },
            "unadjusted": {
                "cv_mean_auroc": auroc,
                "first_component_spearman_phase": rho,
                "selected_lv": 4,
                "p_values": {"delta": unadjusted_p, "angle": 0.5, "shape": 0.5},
                "error": None,
            },
        },
    }


@pytest.mark.parametrize(
    ("auroc", "rho", "unadjusted_p", "met"),
    [
        (0.60, -0.5, 0.01, {"adjusted_auroc": True, "adjusted_rho_phase": True, "unadjusted_rejects": True}),
        (0.70, -0.5, 0.01, {"adjusted_auroc": False, "adjusted_rho_phase": True, "unadjusted_rejects": True}),
        (0.50, 0.5, 0.01, {"adjusted_auroc": False, "adjusted_rho_phase": True, "unadjusted_rejects": True}),
        (0.60, 0.3, 0.01, {"adjusted_auroc": True, "adjusted_rho_phase": False, "unadjusted_rejects": True}),
        (0.60, 0.5, 0.5, {"adjusted_auroc": True, "adjusted_rho_phase": True, "unadjusted_rejects": False}),
    ],
)
def test_calibration_criteria_are_evaluated_on_their_own(auroc, rho, unadjusted_p, met) -> None:
    settings = check.build_settings(_args(command="calibrate"))
    records = [_calibration_record(r, auroc, rho, unadjusted_p) for r in range(10)]
    result = check.calibration_criteria(records, settings)
    assert {name: entry["met"] for name, entry in result["criteria"].items()} == met
    assert result["met"] is all(met.values())
    assert result["criteria"]["adjusted_rho_phase"]["median_abs"] == pytest.approx(abs(rho))


def test_calibration_entry_carries_the_descriptives() -> None:
    settings = check.build_settings(_args(command="calibrate"))
    entry = check.calibration_entry([_calibration_record(r, 0.6, 0.5, 0.01) for r in range(4)], settings)
    assert entry["strengths"]["s_p"] == 1.75 and entry["strengths"]["delta"] == 0.02 and entry["table"] == "adni"
    descriptives = entry["descriptives"]
    assert descriptives["selected_rank"] == {"adjusted": {"4": 4}, "unadjusted": {"4": 4}}
    assert descriptives["first_component_spearman_z"]["median_abs"] == pytest.approx(0.2)
    assert descriptives["phase_offset"]["expression"]["median_abs_d"] == 1.0
    assert descriptives["pooled_relative_eigengap"]["median"] == 0.7
    assert descriptives["angle_null_q95"]["median"] == 30.0
    assert entry["adni_reference"]["selected_lv"] == 8 and len(entry["per_replicate"]) == 4


def test_two_calibrations_leave_two_entries(tmp_path) -> None:
    settings = _settings(reps=1)
    check.calibrate(tmp_path, settings)
    check.calibrate(tmp_path, {**settings, "injection": {**settings["injection"], "s_p": 2.0}})
    history = json.loads((tmp_path / "calibration.json").read_text(encoding="utf-8"))
    assert len(history) == 2
    assert [entry["strengths"]["s_p"] for entry in history] == [1.75, 2.0]
    for entry in history:
        assert set(entry["criteria"]) == {"adjusted_auroc", "adjusted_rho_phase", "unadjusted_rejects"}
        assert isinstance(entry["met"], bool) and entry["replicates"] == 1


# ── Decision (pure) ──────────────────────────────────────────────────────────


def _decision_rows(settings: dict, rates: dict, *, replicates: int = 500, ranks=(3, 5, 9)):
    """p-value rows with ``round(rate·R)`` rejections per (mode, effect, statistic); default 0.05 off target."""

    rows, replicate_rows = [], []
    for mode, effect in check.cells(settings):
        for replicate in range(replicates):
            replicate_rows.append(
                {"mode": mode, "effect": effect, "replicate": replicate, "selected_lv": ranks[replicate % len(ranks)]}
            )
        for statistic in check.STATISTICS:
            rejections = round(rates.get((mode, effect, statistic), 0.05) * replicates)
            rows.extend(
                {
                    "mode": mode,
                    "effect": effect,
                    "replicate": replicate,
                    "statistic": statistic,
                    "p_value": 0.01 if replicate < rejections else 0.5,
                }
                for replicate in range(replicates)
            )
    return rows, replicate_rows


def _curve(mode: str, statistic: str, values) -> dict:
    return {(mode, effect, statistic): rate for effect, rate in zip(check.EFFECTS[mode], values)}


def test_bound_at_the_default_replicates() -> None:
    assert check.anchor_bound(0.05, 500, 2.0) == pytest.approx(0.05 + 2 * math.sqrt(0.0475 / 500))
    assert check.anchor_bound(0.05, 500, 2.0) == pytest.approx(0.0695, abs=5e-5)


def test_single_exceeding_reading_diagnoses_and_issues_no_statement() -> None:
    settings = check.build_settings(_args())
    rates = {("translation", 1.0, "delta"): 0.072, **_curve("magnitude", "delta", (0.5, 0.9, 1.0, 1.0))}
    decision = check.adni_design_decision(*_decision_rows(settings, rates), settings)
    assert decision["verdict"] == "diagnose_design_excess"
    assert decision["size"]["holds"] is False
    assert decision["size"]["exceeding"] == [
        {
            "mode": "translation",
            "effect": 1.0,
            "statistic": "delta",
            "replicates": 500,
            "rate": 0.072,
            "bound": pytest.approx(0.0695, abs=5e-5),
        }
    ]
    assert decision["statements"] == {}
    assert "translation `delta` 0.072" in decision["size"]["sentence"]
    assert decision["consequence"] == check.CONSEQUENCES["diagnose_design_excess"]
    verdicts = {(row["mode"], row["statistic"]): row["verdict"] for row in decision["size"]["readings"]}
    assert list(verdicts.values()).count("exceeds") == 1 and len(verdicts) == 6


def test_informative_effect_found_and_absent() -> None:
    settings = check.build_settings(_args())
    rates = {
        **_curve("orientation", "angle", (0.42, 0.71, 0.86, 0.95)),
        **_curve("shape", "shape", (0.1, 0.2, 0.31, 0.3)),
    }
    decision = check.adni_design_decision(*_decision_rows(settings, rates), settings)
    angle, shape = decision["power"]["angle"], decision["power"]["shape"]
    assert angle["informative_from"] == 0.75 and angle["floor_held_from_informative"] is True
    assert [entry["rate"] for entry in angle["effects"]] == [0.42, 0.71, 0.86, 0.95]
    assert angle["effects"][2]["monte_carlo_se"] == pytest.approx(math.sqrt(0.86 * 0.14 / 500))
    assert angle["effects"][2]["selected_rank"] == {"3": 167, "5": 167, "9": 166}
    assert shape["informative_from"] is None and shape["max_rate"] == 0.31 and shape["effect_at_max"] == 0.75
    assert decision["verdict"] == "power_statement"


def test_power_statement_wording() -> None:
    settings = check.build_settings(_args())
    rates = _curve("magnitude", "delta", (0.3, 0.6, 0.85, 1.0))
    decision = check.adni_design_decision(*_decision_rows(settings, rates), settings)
    assert decision["verdict"] == "power_statement"
    assert decision["power"]["delta"]["informative_from"] == 0.5
    delta = decision["statements"]["delta"]
    assert "rejects a size difference of e = 0.5 (the magnitude construction) with power 0.850" in delta
    assert "a non-rejection argues against size differences of e ≥ 0.5 at this design" in delta
    assert "CV-selected rank 3–9 (median 5)" in delta
    for statistic in ("angle", "shape"):
        text = decision["statements"][statistic]
        assert "reaches at most 0.050" in text
        assert "a non-rejection is uninformative up to the construction maximum (e = 1" in text
    assert decision["consequence"] == check.CONSEQUENCES["power_statement"]


def test_non_rejections_uninformative() -> None:
    settings = check.build_settings(_args())
    rates = {**_curve("magnitude", "delta", (0.1, 0.2, 0.4, 0.79)), **_curve("orientation", "angle", (0.3,) * 4)}
    decision = check.adni_design_decision(*_decision_rows(settings, rates), settings)
    assert decision["verdict"] == "non_rejections_uninformative"
    assert decision["size"]["holds"] is True and decision["indicative"] is False
    assert all(curve["informative_from"] is None for curve in decision["power"].values())
    assert "reaches at most 0.790" in decision["statements"]["delta"]
    assert "orientation construction" in decision["statements"]["angle"]
    assert decision["size"]["sentence"].startswith("Size holds at the ADNI design")
    assert decision["settings_hash"] == check.settings_hash(settings)


def test_low_permutation_run_is_indicative() -> None:
    settings = check.build_settings(_args("--perms", "99"))
    decision = check.adni_design_decision(*_decision_rows(settings, {}), settings)
    assert decision["indicative"] is True
    assert decision["indicative_reasons"]["below_default_permutations"] == 99


def test_unavailable_size_reading_is_indicative() -> None:
    settings = check.build_settings(_args())
    rows, replicate_rows = _decision_rows(settings, {})
    rows = [row for row in rows if not (row["mode"] == "anchor" and row["statistic"] == "shape")]
    decision = check.adni_design_decision(rows, replicate_rows, settings)
    assert decision["indicative"] is True and decision["size"]["holds"] is False
    assert decision["indicative_reasons"]["unavailable_size_readings"] == ["anchor/shape"]
    assert decision["verdict"] == "non_rejections_uninformative"


def test_rank_strata_exclude_nothing_from_the_pooled_rate() -> None:
    settings = check.build_settings(_args())
    rates = _curve("orientation", "angle", (0.42, 0.71, 0.86, 0.95))
    rows, replicate_rows = _decision_rows(settings, rates, ranks=(2, 3, 4, 7, 8, 20))
    decision = check.adni_design_decision(rows, replicate_rows, settings)
    cells = {(row["mode"], row["effect"], row["statistic"]): row for row in decision["descriptive"]["cells"]}
    strata: dict[tuple, list[dict]] = {}
    for row in decision["descriptive"]["rank_strata"]:
        strata.setdefault((row["mode"], row["effect"], row["statistic"]), []).append(row)
    assert set(strata) == set(cells) and len(cells) == 14 * 3
    for key, rows_ in strata.items():
        assert [row["stratum"] for row in rows_] == ["2-3", "4-7", "8+"]
        assert sum(row["replicates"] for row in rows_) == cells[key]["replicates"] == rows_[0]["pooled_replicates"]
        assert sum(row["rejections"] for row in rows_) == cells[key]["rejections"]
        assert rows_[0]["pooled_rate"] == pytest.approx(cells[key]["rate"])


def test_decide_reads_the_merged_tables(tmp_path, sharded_run) -> None:
    out_dir, settings, _ = sharded_run
    check.merge(out_dir, settings)
    decision = check.decide(out_dir, tmp_path / "decision")
    assert decision["indicative"] is True and decision["verdict"] in check.VERDICTS
    assert decision["settings_hash"] == check.settings_hash(settings)
    assert decision["units"] == {"recorded": 10, "expected": 10, "failures": 0}
    written = json.loads((tmp_path / "decision" / "adni_design_decision.json").read_text(encoding="utf-8"))
    assert written["verdict"] == decision["verdict"]
    with (tmp_path / "decision" / "size_readings.csv").open(encoding="utf-8") as handle:
        assert len(list(csv.DictReader(handle))) == 6


def test_decide_refuses_an_unmerged_run(tmp_path) -> None:
    check.ensure_settings(tmp_path, _settings())
    with pytest.raises(SystemExit, match="run merge first"):
        check.decide(tmp_path)
