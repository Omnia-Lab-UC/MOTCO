"""``scripts/corrected_permutation_check.py``: construction, the three-arm units, merge and the decision rule."""

from __future__ import annotations

import json
import math
import sys
from pathlib import Path

import numpy as np
import pytest


def _script(name: str):
    sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
    try:
        return __import__(name)
    finally:
        sys.path.pop(0)


check = _script("corrected_permutation_check")
cs = _script("case_study_type_i_check")
d10 = _script("covariate_type_i_check")
SMALL_CV = {"cv1_splits": 2, "cv2_splits": 3, "n_repeats": 2, "max_components": 5, "random_state": 1203}
SIX = ["z", "x1", "x2", "x3", "x4", "x5"]
ASTROCYTE = [[11, 10, 28], [9, 10, 12]]
MICROGLIA = [[11, 10, 24], [9, 11, 11]]


def _args(*argv: str, out_dir: str = "unused"):
    return check.build_parser().parse_args(["run", *argv, "--out-dir", out_dir])


def _settings(*argv: str, reps: int = 2, perms: int = 9) -> dict:
    settings = check.build_settings(_args(*argv, "--reps", str(reps), "--perms", str(perms)))
    settings["cv"] = dict(SMALL_CV)
    return settings


TYPE_I_RUNS = [
    ("astrocyte", "0.1", ASTROCYTE),
    ("astrocyte", "0.05", ASTROCYTE),
    ("microglia", "0.1", MICROGLIA),
    ("microglia", "0.05", MICROGLIA),
]


# ── Settings and construction ────────────────────────────────────────────────


@pytest.mark.parametrize(("table", "delta", "sizes"), TYPE_I_RUNS)
def test_type_i_runs_have_the_controls_and_the_table(table, delta, sizes) -> None:
    settings = check.build_settings(_args("--profile", "type_i", "--table", table, "--delta", delta))
    assert check.run_label(settings) == f"type-i-{table}-delta{delta}"
    assert check.cells(settings) == [("anchor", 0.0), ("translation", 1.0)]
    assert settings["generator"]["group_stage_sizes"] == sizes
    assert {settings["generator"][f] for f in ("delta_methyl", "delta_expr", "delta_protein")} == {float(delta)}
    assert settings["generator"]["p_dmp"] == 0.1
    assert settings["table"]["name"] == table and len(settings["table"]["sha256"]) == 64
    assert settings["base_seed"] == 2031 and settings["injection"]["s_z"] == 10.0
    assert settings["covariates"]["declared"] == SIX and settings["min_components"] == 2
    assert settings["orientation_estimator"] == "net_displacement"
    assert settings["replicates"] == 500 and settings["permutations"] == 999
    assert settings["arms"] == {
        "none": {"residual_correction": None},
        "df": {"residual_correction": "df"},
        "leverage": {"residual_correction": "leverage"},
    }
    assert "k" not in settings
    assert settings["rule"]["reading_k"] == 3.0 and settings["rule"]["pooled_k"] == 2.0


def test_power_run_has_fourteen_cells_on_the_astrocyte_table() -> None:
    settings = check.build_settings(_args("--profile", "power"))
    assert check.run_label(settings) == "power"
    cells = check.cells(settings)
    assert len(cells) == 14
    assert cells[:2] == [("anchor", 0.0), ("translation", 1.0)]
    assert [e for m, e in cells if m == "magnitude"] == [0.02, 0.05, 0.1, 0.25]
    assert [e for m, e in cells if m == "orientation"] == [0.25, 0.5, 0.75, 1.0]
    assert [e for m, e in cells if m == "shape"] == [0.25, 0.5, 0.75, 1.0]
    assert settings["generator"]["group_stage_sizes"] == ASTROCYTE
    assert {settings["generator"][f] for f in ("delta_methyl", "delta_expr", "delta_protein")} == {2.0}
    assert len(check.units(settings)) == 14 * 500
    assert set(check.expected_runs()) == {
        "type-i-astrocyte-delta0.1",
        "type-i-astrocyte-delta0.05",
        "type-i-microglia-delta0.1",
        "type-i-microglia-delta0.05",
        "power",
    }


@pytest.mark.parametrize("argv", [("--table", "microglia"), ("--delta", "0.1")])
def test_power_profile_is_fixed(argv) -> None:
    with pytest.raises(SystemExit, match="power profile runs the astrocyte table at delta 2.0"):
        check.build_settings(_args("--profile", "power", *argv))


def test_seeds_are_fresh_and_matched() -> None:
    settings = _settings("--profile", "power")
    for replicate in (0, 1):
        first = cs.build_unit("anchor", 0.0, replicate, settings)
        other = cs.build_unit("orientation", 0.5, replicate, settings)
        assert first["seeds"] == other["seeds"]
        assert set(first["seeds"]) == {"generator", "covariates", "rrpp", "noise_covariates"}
        for name in SIX:
            np.testing.assert_array_equal(first["metadata"][name], other["metadata"][name])
        assert {key: first["seeds"][key] for key in ("generator", "covariates", "rrpp")} == d10.replicate_seeds(
            2031, replicate
        )
        parent = d10.replicate_seeds(2026, replicate)
        parent["noise_covariates"] = cs._seed(2026, "noise_covariates", replicate)
        assert not set(first["seeds"].values()) & set(parent.values())
    w_fresh, _ = d10.feature_directions(100, settings)
    w_parent, _ = d10.feature_directions(100, {**settings, "base_seed": 2026})
    assert not np.allclose(w_fresh, w_parent)


@pytest.mark.parametrize(
    "argv",
    [
        ("--table", "microglia"),
        ("--delta", "0.05"),
        ("--base-seed", "2032"),
        ("--perms", "19"),
    ],
)
def test_resume_under_other_settings_is_refused(tmp_path, argv) -> None:
    base = ["--profile", "type_i", "--table", "astrocyte", "--delta", "0.1", "--reps", "2", "--perms", "9",
            "--base-seed", "2031"]
    check.ensure_settings(tmp_path, check.build_settings(_args(*base)))
    changed = list(base)
    for flag, value in zip(argv[::2], argv[1::2]):
        changed[changed.index(flag) + 1] = value
    with pytest.raises(SystemExit, match="different settings"):
        check.ensure_settings(tmp_path, check.build_settings(_args(*changed)))
    assert not (tmp_path / "records").exists()


def test_over_headroom_cell_refuses_to_start(tmp_path) -> None:
    settings = _settings("--profile", "power")
    settings["effects"]["orientation"] = [0.25, 3.0]
    with pytest.raises(SystemExit, match=r'\["orientation", 3.0\]'):
        check.run_shard(tmp_path, settings, 0, 1)
    assert check.shard_path(tmp_path, 0, 1).read_text() == ""


# ── Units ────────────────────────────────────────────────────────────────────


def _records(out_dir: Path) -> list[dict]:
    lines = [
        line
        for path in sorted((out_dir / "records").glob("shard_*.jsonl"))
        for line in path.read_text(encoding="utf-8").splitlines()
    ]
    return sorted((json.loads(line) for line in lines), key=lambda record: record["key"])


@pytest.fixture(scope="module")
def sharded_run(tmp_path_factory):
    out_dir = tmp_path_factory.mktemp("type-i")
    settings = _settings("--profile", "type_i", "--table", "astrocyte", "--delta", "0.1")
    check.ensure_settings(out_dir, settings)
    evaluated = [check.run_shard(out_dir, settings, index, 2) for index in (0, 1)]
    return out_dir, settings, evaluated


def test_two_shards_cover_every_unit_once(sharded_run) -> None:
    out_dir, settings, evaluated = sharded_run
    assert evaluated == [2, 2]
    keys = [record["key"] for record in _records(out_dir)]
    assert sorted(keys) == sorted(check.unit_key(*unit) for unit in check.units(settings))
    assert len(keys) == len(set(keys)) == 4


def test_arms_share_the_measurement(sharded_run) -> None:
    out_dir, settings, _ = sharded_run
    for record in _records(out_dir):
        assert record["run"] == "type-i-astrocyte-delta0.1"
        assert record["settings_hash"] == check.settings_hash(settings)
        arms = record["arms"]
        assert set(arms) == {"none", "df", "leverage"}
        none, df, leverage = arms["none"], arms["df"], arms["leverage"]
        for other in (df, leverage):
            assert other["error"] is None
            differing = {field for field in check.ARM_FIELDS if other[field] != none[field]}
            assert differing <= {"p_values", "null_summary", "residual_correction"}
            assert {"null_summary", "residual_correction"} <= differing
        assert none["latent_digest"] == df["latent_digest"] == leverage["latent_digest"]
        assert none["residual_correction"] is None
        assert df["residual_correction"]["kind"] == "df"
        assert df["residual_correction"]["n"] == 80 and df["residual_correction"]["reduced_rank"] == 10
        assert df["residual_correction"]["kappa"] == pytest.approx(math.sqrt(80 / 70))
        assert leverage["residual_correction"]["kind"] == "leverage"
        assert none["selected_lv"] >= 2 and none["cv_mean_auroc"] is not None
        assert abs(none["first_component_spearman_z"]) <= 1


def test_resume_skips_recorded_units(sharded_run) -> None:
    out_dir, settings, _ = sharded_run
    before = _records(out_dir)
    assert check.run_shard(out_dir, settings, 0, 1) == 0
    assert check.run_shard(out_dir, settings, 1, 2) == 0
    assert _records(out_dir) == before


def test_torn_line_is_rerun(tmp_path, sharded_run) -> None:
    out_dir, settings, _ = sharded_run
    path = check.shard_path(out_dir, 0, 2)
    lines = path.read_text(encoding="utf-8").splitlines()
    copy = tmp_path / "run"
    check.ensure_settings(copy, settings)
    torn = check.shard_path(copy, 0, 2)
    torn.parent.mkdir(parents=True)
    torn.write_text(lines[0] + "\n" + lines[1][: len(lines[1]) // 2], encoding="utf-8")
    assert check.run_shard(copy, settings, 0, 2) == 1
    assert len(_records(copy)) == 2


def test_arm_disagreement_raises_and_writes_nothing(tmp_path, monkeypatch) -> None:
    settings = _settings("--profile", "type_i", "--table", "astrocyte", "--delta", "0.1", reps=1)
    real = check.analyze_trajectories

    def perturbed(dataset, params):
        result = real(dataset, params)
        if params.residual_correction == "leverage":
            result.latent_scores.iloc[0, 0] += 1e-12
        return result

    monkeypatch.setattr(check, "analyze_trajectories", perturbed)
    check.ensure_settings(tmp_path, settings)
    with pytest.raises(RuntimeError, match=r"\['leverage'\] differ from arm 'none'"):
        check.run_shard(tmp_path, settings, 0, 1)
    assert check.shard_path(tmp_path, 0, 1).read_text(encoding="utf-8") == ""


def test_merge_writes_the_tables_and_is_idempotent(sharded_run) -> None:
    out_dir, settings, _ = sharded_run
    first = check.merge(out_dir, settings)
    contents = {
        name: (out_dir / name).read_text(encoding="utf-8")
        for name in ("corrected_permutation_pvalues.csv", "corrected_permutation_rates.csv")
    }
    second = check.merge(out_dir, settings)
    assert first == second
    for name, text in contents.items():
        assert (out_dir / name).read_text(encoding="utf-8") == text
    assert not (out_dir / "corrected_permutation_contrasts.csv").exists()  # Type I runs have no power cells
    assert first["units"] == first["expected"] == 4 and first["missing"] == []
    assert len(contents["corrected_permutation_pvalues.csv"].splitlines()) == 1 + 4 * 3 * 3


def test_decide_over_one_run_is_indicative(tmp_path, sharded_run) -> None:
    out_dir, settings, _ = sharded_run
    check.merge(out_dir, settings)
    decision = check.decide([out_dir], tmp_path / "decision")
    assert decision["indicative"] is True
    assert decision["indicative_reasons"]["below_default_permutations"] == ["type-i-astrocyte-delta0.1"]
    assert "power" in decision["indicative_reasons"]["missing_runs"]
    assert decision["runs"]["type-i-astrocyte-delta0.1"]["settings_hash"] == check.settings_hash(settings)
    assert (tmp_path / "decision" / "size_readings.csv").exists()
    stored = json.loads((tmp_path / "decision" / "corrected_permutation_decision.json").read_text())
    assert stored["verdict"] in check.VERDICTS
    assert set(stored["calibration"]["type-i-astrocyte-delta0.1"]) >= {"auroc", "rho", "selected_rank"}


# ── Merge tables on synthetic records ────────────────────────────────────────


def _record(mode: str, effect: float, replicate: int, p_values: dict[str, dict[str, float]], settings) -> dict:
    arms = {}
    for arm in check.ARMS:
        arms[arm] = {field: None for field in check.ARM_FIELDS}
        arms[arm].update(
            {
                "p_values": p_values[arm],
                "statistics": {},
                "null_summary": {},
                "selected_lv": 2,
                "cv_mean_auroc": 0.7,
                "first_component_spearman_z": 0.9,
            }
        )
    return {
        "key": check.unit_key(mode, effect, replicate),
        "run": check.run_label(settings),
        "mode": mode,
        "effect": effect,
        "replicate": replicate,
        "settings_hash": check.settings_hash(settings),
        "arms": arms,
        "runtime_seconds": 1.0,
    }


def _p(reject: bool) -> float:
    return 0.01 if reject else 0.5


def test_rates_and_contrasts_match_hand_counts(tmp_path) -> None:
    settings = _settings("--profile", "power", reps=3)
    rejects = {  # (replicate, effect) → arms rejecting orientation/angle
        (0, 0.25): {"leverage"},
        (0, 0.5): {"leverage", "df"},
        (0, 1.0): {"leverage"},
        (1, 0.75): {"df", "none"},
        (2, 1.0): {"df", "leverage", "none"},
    }
    records = []
    for replicate in range(3):
        for mode, effect in check.cells(settings):
            arms = {}
            for arm in check.ARMS:
                reject = mode == "orientation" and arm in rejects.get((replicate, effect), set())
                arms[arm] = {"delta": 0.5, "angle": _p(reject), "shape": 0.5}
            records.append(_record(mode, effect, replicate, arms, settings))
    check.ensure_settings(tmp_path, settings)
    path = check.shard_path(tmp_path, 0, 1)
    path.parent.mkdir(parents=True)
    path.write_text("".join(json.dumps(r) + "\n" for r in records), encoding="utf-8")
    first = check.merge(tmp_path, settings)
    assert first["units"] == 42
    rates = check.diagnosis._read_csv(tmp_path / "corrected_permutation_rates.csv")
    rate = {(r["mode"], float(r["effect"]), r["arm"], r["statistic"]): r for r in rates}
    assert rate[("orientation", 1.0, "leverage", "angle")]["rejections"] == "2"
    assert float(rate[("orientation", 1.0, "leverage", "angle")]["rate"]) == pytest.approx(2 / 3)
    assert rate[("orientation", 0.75, "none", "angle")]["rejections"] == "1"
    assert rate[("anchor", 0.0, "df", "angle")]["rejections"] == "0"
    assert rate[("anchor", 0.0, "df", "angle")]["selected_rank"] == '{"2":3}'
    contrasts = check.diagnosis._read_csv(tmp_path / "corrected_permutation_contrasts.csv")
    row = {(c["scope"], c["mode"], c["effect"], c["statistic"], c["arm_b"], c["arm_a"]): c for c in contrasts}
    pooled = row[("mode", "orientation", None, "angle", "leverage", "df")]
    # leverage-only: (0, 0.25), (0, 1.0); df-only: (1, 0.75). Replicate sums d = (+2, -1, 0).
    assert (pooled["b_only"], pooled["a_only"], pooled["net"], pooled["clusters"]) == ("2", "1", "1", "3")
    assert pooled["role"] == "tie_break"
    assert float(pooled["p_b_greater"]) == pytest.approx(check.diagnosis.sign_flip_p_value([2, -1, 0]))
    assert float(pooled["p_b_greater"]) == pytest.approx(0.5)
    cell = row[("cell", "orientation", "1.0", "angle", "leverage", "none")]
    assert (cell["b_only"], cell["a_only"], cell["readings"]) == ("1", "0", "3")
    assert float(cell["b_rate"]) == pytest.approx(2 / 3) and float(cell["a_rate"]) == pytest.approx(1 / 3)
    assert row[("mode", "orientation", None, "angle", "df", "none")]["role"] == "descriptive"
    assert check.merge(tmp_path, settings) == first


def test_clustered_p_value_pools_a_replicates_cells() -> None:
    rows = []
    for replicate in range(2):
        for effect in (0.25, 0.5):
            for arm in check.ARMS:
                reject = replicate == 0 and arm == "leverage"
                rows.append(
                    {"run": "power", "mode": "orientation", "effect": effect, "replicate": replicate, "arm": arm,
                     "statistic": "angle", "p_value": _p(reject)}
                )
    contrasts = check.contrast_rows(rows, {"orientation": [0.25, 0.5]}, alpha=0.05)
    pooled = next(
        c for c in contrasts
        if c["scope"] == "mode" and c["statistic"] == "angle" and (c["arm_b"], c["arm_a"]) == ("leverage", "df")
    )
    # Two leverage-only readings in one replicate form one cluster (d = 2): p = 1/2, not the 1/4 of two clusters.
    assert pooled["b_only"] == 2 and pooled["clusters"] == 2
    assert pooled["p_b_greater"] == pytest.approx(0.5)
    assert check.diagnosis.sign_flip_p_value([1, 1]) == pytest.approx(0.25)


# ── Decision ─────────────────────────────────────────────────────────────────


R = 500


def _decision_settings(perms: int = 999) -> dict[str, dict]:
    out = {}
    for table in ("astrocyte", "microglia"):
        for delta in ("0.1", "0.05"):
            settings = check.build_settings(
                _args("--profile", "type_i", "--table", table, "--delta", delta, "--perms", str(perms))
            )
            out[check.run_label(settings)] = settings
    power = check.build_settings(_args("--profile", "power", "--perms", str(perms)))
    out["power"] = power
    return out


def _rows(run: str, mode: str, effect: float, arm: str, statistic: str, rejected: set[int]) -> list[dict]:
    return [
        {"run": run, "mode": mode, "effect": effect, "replicate": r, "arm": arm, "statistic": statistic,
         "p_value": _p(r in rejected)}
        for r in range(R)
    ]


def _size_rows(settings: dict[str, dict], overrides: dict | None = None) -> list[dict]:
    """Every control reading at 25/500 unless ``overrides[(run or table, arm, statistic)]`` names a count."""

    overrides = overrides or {}
    rows = []
    for run, run_settings in settings.items():
        table = run_settings["table"]["name"]
        for mode in ("anchor", "translation"):
            effect = float(run_settings["controls"][mode])
            for arm in check.ARMS:
                for statistic in check.STATISTICS:
                    count = overrides.get((run, arm, statistic), overrides.get((table, arm, statistic), 25))
                    if run == "power" and (table, arm, statistic) in overrides:
                        count = 25
                    rows.extend(_rows(run, mode, effect, arm, statistic, set(range(count))))
    return rows


def _power_rows(gains: dict[tuple[str, str], int]) -> list[dict]:
    """Target cells where df and leverage agree, except ``gains[(mode, statistic)]`` replicates at the top effect.

    A positive count makes leverage reject where df does not; a negative one the reverse.
    """

    rows = []
    for mode, statistic in check.TARGETS:
        effects = check.POWER_EFFECTS[mode]
        for effect in effects:
            shared = set(range(100))
            extra = gains.get((mode, statistic), 0) if effect == effects[-1] else 0
            extra_set = set(range(100, 100 + abs(extra)))
            for arm in check.ARMS:
                rejected = set(shared)
                if (arm == "leverage" and extra > 0) or (arm == "df" and extra < 0):
                    rejected |= extra_set
                rows.extend(_rows("power", mode, effect, arm, statistic, rejected))
    return rows


def _decide(overrides=None, gains=None, perms: int = 999) -> dict:
    settings = _decision_settings(perms)
    return check.corrected_permutation_decision(_size_rows(settings, overrides), _power_rows(gains or {}), settings)


def test_bounds_at_the_default_replicates() -> None:
    assert check.size_bound(0.05, 500, 3.0) == pytest.approx(0.0792, abs=5e-5)
    assert check.size_bound(0.05, 2000, 2.0) == pytest.approx(0.0597, abs=5e-5)
    assert math.floor(check.size_bound(0.05, 500, 3.0) * 500) == 39
    assert math.floor(check.size_bound(0.05, 2000, 2.0) * 2000) == 119


def test_only_df_holds() -> None:
    decision = _decide({("type-i-microglia-delta0.1", "leverage", "angle"): 40})
    assert decision["verdict"] == "adopt_df" and decision["tie_break"] is None
    assert decision["holds"] == {"df": True, "leverage": False}
    assert decision["indicative"] is False
    assert any(f["tier"] == "reading" and f["statistic"] == "angle" for f in decision["failing"]["leverage"])
    assert "--residual-correction df" in decision["consequence"]


def test_only_leverage_holds() -> None:
    decision = _decide({("type-i-astrocyte-delta0.05", "df", "shape"): 40})
    assert decision["verdict"] == "adopt_leverage" and decision["tie_break"] is None


def test_leverage_wins_the_tie_break_on_a_gain_without_loss() -> None:
    decision = _decide(gains={("orientation", "angle"): 20})
    assert decision["verdict"] == "adopt_leverage"
    pairs = {(p["mode"], p["statistic"]): p for p in decision["tie_break"]["pairs"]}
    assert pairs[("orientation", "angle")]["gain"] and pairs[("orientation", "angle")]["p_gain"] == pytest.approx(
        2.0**-20
    )
    assert not any(p["loss"] for p in pairs.values())


def test_no_significant_difference_keeps_df() -> None:
    decision = _decide()
    assert decision["verdict"] == "adopt_df"
    assert decision["tie_break"]["status"] == "computed"
    assert all(p["p_gain"] == 1.0 and p["p_loss"] == 1.0 for p in decision["tie_break"]["pairs"])


def test_gain_and_loss_keep_df() -> None:
    decision = _decide(gains={("orientation", "angle"): 20, ("shape", "shape"): -20})
    assert decision["verdict"] == "adopt_df"
    pairs = {(p["mode"], p["statistic"]): p for p in decision["tie_break"]["pairs"]}
    assert pairs[("orientation", "angle")]["gain"] and pairs[("shape", "shape")]["loss"]


def test_neither_holds_lists_the_pooled_microglia_delta() -> None:
    # 31/500 = 0.062 in each microglia reading: under the 0.0792 reading bound, over the 0.0597 pooled bound.
    overrides = {("microglia", arm, "delta"): 31 for arm in check.ARMS}
    decision = _decide(overrides)
    assert decision["verdict"] == "no_correction_holds"
    for arm in ("df", "leverage"):
        failing = decision["failing"][arm]
        assert [f["tier"] for f in failing] == ["pooled"]
        assert failing[0]["table"] == "microglia" and failing[0]["statistic"] == "delta"
        assert failing[0]["rate"] == pytest.approx(0.062) and failing[0]["replicates"] == 2000
    assert "keeps reading covariate-adjusted delta against its realized size" in decision["consequence"]
    readings = [r for r in decision["readings"] if r["table"] == "microglia" and r["statistic"] == "delta"]
    assert all(r["verdict"] == "holds" for r in readings if r["arm"] != "none")


def test_uncorrected_excess_is_descriptive() -> None:
    decision = _decide({("astrocyte", "none", "delta"): 33})
    assert decision["verdict"] == "adopt_df"
    assert decision["uncorrected"]["uncorrected_excess_reproduced"] is True
    pooled = {(p["table"], p["statistic"]): p for p in decision["uncorrected"]["pooled"]}
    assert pooled[("astrocyte", "delta")]["rate"] == pytest.approx(0.066)
    assert _decide()["uncorrected"]["uncorrected_excess_reproduced"] is False


def test_smoke_decision_is_indicative() -> None:
    decision = _decide(perms=99)
    assert decision["indicative"] is True
    assert decision["indicative_reasons"]["below_default_permutations"] == sorted(check.expected_runs())
    assert decision["indicative_reasons"]["missing_runs"] == []


def test_decision_records_every_settings_hash_and_the_power_cost() -> None:
    settings = _decision_settings()
    decision = check.corrected_permutation_decision(_size_rows(settings), _power_rows({}), settings)
    assert {run: entry["settings_hash"] for run, entry in decision["runs"].items()} == {
        run: check.settings_hash(s) for run, s in settings.items()
    }
    assert len(decision["readings"]) == 5 * 2 * 3 * 3
    assert len(decision["pooled"]) == 2 * 3 * 3
    cost = {(c["correction"], c["mode"], c["statistic"]): c for c in decision["power_cost"]}
    assert cost[("df", "orientation", "angle")]["uncorrected_rate"] == pytest.approx(0.2)
    assert cost[("leverage", "shape", "shape")]["p_loss"] == 1.0
