"""``scripts/orthonormal_coordinates_check.py``: units, the two-arm unit, resume, merge, reproduction, decision.

The unit runs on a small stage-span design at fast CV knobs; the profiles'
units are checked against the study engine and the span-split check without
evaluating them. No real data are read.
"""

from __future__ import annotations

import dataclasses
import json
import shutil
import sys
from pathlib import Path

import numpy as np
import pytest
from scipy.stats import binom


def _script(name: str):
    sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
    try:
        return __import__(name)
    finally:
        sys.path.pop(0)


occ = _script("orthonormal_coordinates_check")
lsc, ssc, leak, grade = occ.lsc, occ.ssc, occ.leak, occ.grade

REPO = Path(__file__).resolve().parents[1]
FAST_CV = {"layers": ["methylation", "expression"], "n_repeats": 1, "max_components": 6}


@pytest.fixture(scope="module")
def paper_settings():
    return occ.paper_grade_settings()


@pytest.fixture(scope="module")
def weak_settings():
    return occ.weak_settings()


# ── 2.1 Units, resume and settings ───────────────────────────────────────────


def test_unit_counts(paper_settings, weak_settings) -> None:
    paper = occ.enumerate_units(paper_settings)
    assert len(paper) == 65 * 500 == 32_500
    cells = {(column, mode, effect) for column, mode, effect, _ in paper}
    assert len(cells) == 65
    for cell in cells:
        assert sorted(r for c, m, e, r in paper if (c, m, e) == cell) == list(range(500))
    weak = occ.enumerate_units(weak_settings)
    assert len(weak) == 2 * 500 == 1_000
    assert {(c, m, e) for c, m, e, _ in weak} == {("seaad_weak", "none", 0.0), ("seaad_weak", "translation", 1.0)}
    assert paper_settings["permutations"] == weak_settings["permutations"] == 999
    assert paper_settings["span_scale"] == "planted" and paper_settings["span_capture"] is True


def _fake_arm(p: float, *, share: float = 0.01) -> dict:
    return {
        "rank": 2,
        "cv_mean_auroc": 0.9,
        "pair_statistics": {"delta": 1.0, "angle": 10.0, "shape": 0.1},
        "p_values": {statistic: p for statistic in occ.STATISTICS},
        "pooled_eigengap": 0.2,
        "angle_null_q95": 30.0,
        "span_capture": {"group_difference": 0.9, "pooled_stage": 0.95},
        "group_levels": ["A", "B"],
        "population_latent": {"path_lengths": {"A": 10.0, "B": 10.0 * (1 + share)}},
    }


def _fake_record(key, settings, *, p=0.5, p_orthonormal=None, status="completed", consistent=True):
    column, mode, effect, replicate = key
    record = {
        "profile": settings["profile"],
        "column": column,
        "cell": occ._cell_label(settings, mode, effect),
        "trajectory_mode": mode,
        "effect_size": float(effect),
        "replicate": int(replicate),
        "settings_hash": lsc.settings_hash(settings),
        "generator_seed": 1,
        "rrpp_seed": 2,
        "status": status,
        "runtime_seconds": {"unit": 2.0, "generate": 0.5, "fit": 1.0, "scores": 0.2, "orthonormal": 0.2},
    }
    if status != "completed":
        return {**record, "error_type": "RuntimeError", "error_message": "boom"}
    deviation = 0.0 if consistent else 0.5
    record["consistency"] = {
        "tolerance": 1e-8,
        "deviations": {"orthonormal_basis": deviation, "same_span": 0.0},
        "consistent": consistent,
        "failing": [] if consistent else [{"check": "orthonormal_basis", "deviation": deviation}],
    }
    record["arms"] = {"scores": _fake_arm(p), "orthonormal": _fake_arm(p if p_orthonormal is None else p_orthonormal)}
    return record


def test_resume_skips_completed_units(tmp_path, monkeypatch) -> None:
    settings = occ.paper_grade_settings(replicates=2, permutations=39)
    calls: list = []

    def fake(key, settings):
        calls.append(key)
        return _fake_record(key, settings)

    monkeypatch.setattr(occ, "evaluate_unit", fake)
    assert occ.run_shard(tmp_path, settings, 0, 2, limit=3) == 3
    with lsc.shard_path(tmp_path / "paper_grade", 0, 2).open("a") as handle:
        handle.write('{"torn": ')  # an interrupted write
    first = list(calls)
    total = len(occ.enumerate_units(settings))
    assert occ.run_shard(tmp_path, settings, 0, 2) == total // 2 - 3
    assert not set(first) & set(calls[3:])
    occ.run_shard(tmp_path, settings, 1, 2)
    summary = occ.merge(tmp_path, "paper_grade")
    assert summary["complete_units"] == total == 130 and summary["missing_units"] == 0
    assert len((tmp_path / "paper_grade" / "merged.jsonl").read_text().splitlines()) == total


def test_settings_mismatch_is_refused_naming_the_keys(tmp_path, monkeypatch) -> None:
    monkeypatch.setattr(occ, "evaluate_unit", _fake_record)
    occ.run_shard(tmp_path, occ.paper_grade_settings(replicates=2, permutations=39), 0, 1, limit=1)
    with pytest.raises(SystemExit, match="differing: .*permutations"):
        occ.run_shard(tmp_path, occ.paper_grade_settings(replicates=2, permutations=49), 0, 1, limit=1)


def _profile_copies(tmp_path: Path) -> dict[str, Path]:
    copies = {}
    for name, path in occ.PAPER_GRADE_PROFILES.items():
        copies[name] = tmp_path / "profiles" / path.name
        copies[name].parent.mkdir(parents=True, exist_ok=True)
        shutil.copy(path, copies[name])
    return copies


def test_a_changed_profile_is_refused(tmp_path, monkeypatch) -> None:
    calls: list = []
    monkeypatch.setattr(occ, "evaluate_unit", lambda key, settings: calls.append(key) or _fake_record(key, settings))
    profiles = _profile_copies(tmp_path)
    out = tmp_path / "run"
    occ.run_shard(out, occ.paper_grade_settings(replicates=2, profiles=profiles), 0, 1, limit=1)
    assert len(calls) == 1
    changed = profiles["span_geometry"]
    changed.write_text(changed.read_text() + "\n")
    with pytest.raises(SystemExit, match=r"Parent profile:span_geometry .* changed"):
        occ.run_shard(out, occ.paper_grade_settings(replicates=2, profiles=profiles), 0, 1, limit=1)
    with pytest.raises(SystemExit, match=r"Parent profile:span_geometry .* changed"):
        occ.merge(out, "paper_grade")
    assert len(calls) == 1


# ── 2.2 The profiles' units are their parents' ───────────────────────────────


SAMPLED_PAPER_KEYS = [
    ("seaad", "orientation", 0.5, 7),
    ("n200_four_stages", "magnitude", 0.05, 3),
    ("n80", "none", 0.0, 499),
    ("n200_three_blocks", "translation", 0.1, 0),
    ("n200", "shape", 0.25, 250),
]


@pytest.mark.parametrize("key", SAMPLED_PAPER_KEYS)
def test_paper_grade_units_are_the_paper_grades(paper_settings, key) -> None:
    from motco.simulations.grid import derive_replicate_seed
    from motco.simulations.study import enumerate_study, load_study_config

    column, mode, effect, replicate = key
    profile = "span_magnitude" if mode == "magnitude" else "span_geometry"
    study = enumerate_study(load_study_config(occ.PAPER_GRADE_PROFILES[profile]))
    matches = [
        cell for cell in study.cells
        if cell.phase in ("power_primary", "power_design")
        and leak.column_of(cell) == column
        and float(cell.generator_params.group_effect_size) == effect
        and (effect == 0.0 or cell.generator_params.trajectory_mode == mode)
    ]
    assert matches, key
    cell = matches[0]
    seed = derive_replicate_seed(cell, replicate)
    expected_generator = dataclasses.replace(cell.generator_params, seed=seed)
    expected_rrpp = cell.evaluation_params.seed if cell.evaluation_params.seed is not None else seed

    generator, evaluation = occ.unit_params(key, paper_settings)
    assert generator == expected_generator
    assert evaluation.seed == expected_rrpp
    assert evaluation.permutations == 999 and evaluation.n_jobs == 1 and evaluation.span_capture is True
    assert evaluation.latent_coordinates is None
    theirs_generator, theirs_evaluation = leak.unit_params(
        leak.column_cells(occ.PAPER_GRADE_PROFILES[profile]), key, span_scale="planted", permutations=999
    )
    assert generator == theirs_generator
    assert dataclasses.replace(evaluation, n_jobs=theirs_evaluation.n_jobs) == theirs_evaluation


@pytest.mark.parametrize("cell", ["anchor", "translation"])
def test_weak_units_are_the_span_split_units(weak_settings, cell) -> None:
    spec = lsc.cell_spec(weak_settings["units"]["comparison"], cell)
    key = ("seaad_weak", spec["mode"], float(spec["effect"]), 17)
    generator, evaluation = occ.unit_params(key, weak_settings)
    assert generator.seed == lsc.generator_seed(2035, "seaad_weak", 17)
    assert evaluation.seed == lsc.rrpp_seed(2035, "seaad_weak", cell, 17)
    committed = ssc.load_parent(occ.SPAN_SPLIT_DIR)["settings"]
    expected = ssc.production_params(committed, seed=evaluation.seed, n_stages=3)
    assert evaluation == dataclasses.replace(expected, include_null_distributions=False)
    mine = occ.generate_dataset(generator)
    theirs = lsc.unit_dataset(lsc.Unit("seaad_weak", cell, 17), committed["comparison"])
    for name in ("methylation", "expression", "proteomics", "metadata"):
        assert getattr(mine, name).equals(getattr(theirs, name))


def test_weak_delta_is_the_calibration(weak_settings) -> None:
    calibration = json.loads((ssc.COMPARISON_DIR / "calibration" / "calibration.json").read_text())
    assert weak_settings["units"]["comparison"]["weak_delta"]["value"] == calibration["adopted"]["delta"] == 0.11
    generator, _ = occ.unit_params(("seaad_weak", "none", 0.0, 0), weak_settings)
    assert generator.delta_methyl == generator.delta_expr == generator.delta_protein == 0.11


# ── 2.3 The unit ─────────────────────────────────────────────────────────────


def _small_generator():
    from motco.simulations.semisynthetic import SemiSyntheticTrajectoryParams

    return SemiSyntheticTrajectoryParams(
        seed=31, n_samples=90, p_dmp=0.1, trajectory_mode="orientation", group_effect_size=0.5,
        construction="stage_span", span_layers=("methylation", "expression"), span_scale="planted",
    )


def _small_params():
    from motco.simulations.evaluation import SimulationEvaluationParams

    return SimulationEvaluationParams(
        integration_method="pls", integration_params=FAST_CV, permutations=19, n_jobs=1, seed=5, span_capture=True,
    )


KEY = ("n200", "orientation", 0.5, 3)


@pytest.fixture(scope="module")
def unit_record(paper_settings):
    with pytest.MonkeyPatch.context() as patch:
        patch.setattr(occ, "unit_params", lambda key, settings: (_small_generator(), _small_params()))
        return occ.evaluate_unit(KEY, paper_settings)


def test_record_carries_every_field(unit_record, paper_settings) -> None:
    record = unit_record
    assert record["status"] == "completed", record.get("error_message")
    assert (record["profile"], record["column"], record["cell"]) == ("paper_grade", "n200", "orientation_e0.5")
    assert (record["trajectory_mode"], record["effect_size"], record["replicate"]) == ("orientation", 0.5, 3)
    assert record["generator_seed"] == 31 and record["rrpp_seed"] == 5
    assert record["settings_hash"] == lsc.settings_hash(paper_settings)
    for arm in occ.ARMS:
        entry = record["arms"][arm]
        assert set(entry["pair_statistics"]) == set(entry["p_values"]) == set(occ.STATISTICS)
        assert entry["pooled_eigengap"] is not None and entry["angle_null_q95"] is not None
        assert set(entry["population_latent"]["path_lengths"]) == {"A", "B"}
    assert record["consistency"]["consistent"] is True
    for name in occ.CONSISTENCY_CHECKS:
        assert record["consistency"]["deviations"][name] <= 1e-8
    assert {"generate", "fit", "scores", "orthonormal", "unit"} <= set(record["runtime_seconds"])


def test_both_arms_share_the_fit(unit_record) -> None:
    scores, orthonormal = unit_record["arms"]["scores"], unit_record["arms"]["orthonormal"]
    assert scores["rank"] == orthonormal["rank"]
    assert scores["cv_mean_auroc"] == orthonormal["cv_mean_auroc"]
    assert scores["span_capture"] == orthonormal["span_capture"]
    assert scores["population_latent"] != orthonormal["population_latent"]


def test_arms_equal_plain_evaluations(unit_record) -> None:
    from motco.simulations.evaluation import evaluate_semisynthetic_trajectory

    dataset = occ.generate_dataset(_small_generator())
    params = _small_params()
    plain = occ.arm_record(evaluate_semisynthetic_trajectory(dataset, params))
    assert unit_record["arms"]["scores"] == json.loads(json.dumps(plain))
    scratch = occ.arm_record(
        evaluate_semisynthetic_trajectory(dataset, dataclasses.replace(params, latent_coordinates="orthonormal"))
    )
    assert unit_record["arms"]["orthonormal"] == json.loads(json.dumps(scratch))


def test_an_inconsistent_unit_is_flagged(monkeypatch) -> None:
    import motco.simulations.evaluation as evaluation_module

    original = evaluation_module.pls_orthonormal_from_pls

    def patched(dataset, pls, params):
        latent = original(dataset, pls, params)
        skewed = dataclasses.replace(latent.linear_map, basis=latent.linear_map.basis * 1.01)
        return dataclasses.replace(latent, linear_map=skewed)

    monkeypatch.setattr(evaluation_module, "pls_orthonormal_from_pls", patched)
    measured = occ.measure_unit(occ.generate_dataset(_small_generator()), dataclasses.replace(_small_params(),
                                                                                             permutations=0))
    consistency = measured["consistency"]
    assert consistency["consistent"] is False
    failing = {row["check"]: row["deviation"] for row in consistency["failing"]}
    assert failing["orthonormal_basis"] == pytest.approx(1.01**2 - 1, rel=1e-6)


def test_a_failing_unit_is_recorded_without_aborting_the_shard(tmp_path, monkeypatch) -> None:
    def broken(key, settings):
        raise RuntimeError(f"no data for {key[0]}")

    monkeypatch.setattr(occ, "unit_params", broken)
    settings = occ.paper_grade_settings(replicates=2)
    assert occ.run_shard(tmp_path, settings, 0, 1, limit=2) == 2
    records = [json.loads(line) for line in lsc.shard_path(tmp_path / "paper_grade", 0, 1).read_text().splitlines()]
    assert [r["status"] for r in records] == ["failed", "failed"]
    assert records[0]["error_type"] == "RuntimeError" and "no data" in records[0]["error_message"]


# ── 3.1 Merge ────────────────────────────────────────────────────────────────


def _rows(rows, **match):
    return [row for row in rows if all(row[k] == v for k, v in match.items())]


def test_tables_count_failures_inconsistencies_and_pairs() -> None:
    settings = occ.paper_grade_settings(replicates=4)
    keys = [key for key in occ.enumerate_units(settings) if key[:3] == ("n80", "orientation", 0.5)]
    records = [
        _fake_record(keys[0], settings, p=0.01, p_orthonormal=0.01),  # both reject
        _fake_record(keys[1], settings, p=0.01, p_orthonormal=0.5),  # only scores rejects (b)
        _fake_record(keys[2], settings, p=0.5, p_orthonormal=0.5, consistent=False),
        _fake_record(keys[3], settings, status="failed"),
    ]
    tables = occ.tables(records, settings)
    scores = _rows(tables["rates"], column="n80", statistic="angle", arm="scores")[0]
    assert (scores["n_units"], scores["n_failed"], scores["n_available"], scores["n_rejected"]) == (4, 1, 3, 2)
    orthonormal = _rows(tables["rates"], column="n80", statistic="angle", arm="orthonormal")[0]
    assert (orthonormal["n_available"], orthonormal["n_rejected"]) == (3, 1)
    assert orthonormal["se"] == pytest.approx(np.sqrt((1 / 3) * (2 / 3) / 3))
    paired = _rows(tables["paired"], column="n80", statistic="angle")[0]
    assert (paired["n_paired"], paired["b"], paired["c"]) == (3, 1, 0)
    assert paired["difference"] == pytest.approx(-1 / 3)
    assert paired["p_two_sided"] == 1.0 and paired["p_loss"] == 0.5
    cell = _rows(tables["cells"], column="n80", arm="orthonormal")[0]
    assert cell["n"] == 3 and cell["median_population_share"] == pytest.approx(0.01)
    summary = occ.accounting(records, settings, expected=4, duplicates=0)
    assert summary["failed_units"] == 1 and summary["failures_by_cell"] == {"n80/orientation_e0.5": 1}
    assert summary["inconsistent_units"] == 1
    assert summary["inconsistent"][0]["failing"][0]["check"] == "orthonormal_basis"
    assert summary["max_consistency_deviation"]["orthonormal_basis"] == 0.5


def test_merge_counts_duplicates(tmp_path, monkeypatch) -> None:
    monkeypatch.setattr(occ, "evaluate_unit", _fake_record)
    settings = occ.paper_grade_settings(replicates=2)
    occ.run_shard(tmp_path, settings, 0, 1, limit=4)
    directory = tmp_path / "paper_grade"
    shutil.copy(lsc.shard_path(directory, 0, 1), lsc.shard_path(directory, 0, 2))
    summary = occ.merge(tmp_path, "paper_grade")
    accounting = json.loads((directory / "accounting.json").read_text())
    assert summary["complete_units"] == 4 and accounting["duplicate_records"] == 4
    assert summary["reproduction_passed"] is False
    for name in ("rate_table.csv", "paired_table.csv", "cell_table.csv", "reproduction.json", "merged.jsonl"):
        assert (directory / name).exists()


def test_merge_refuses_records_written_under_other_settings(tmp_path, monkeypatch) -> None:
    monkeypatch.setattr(occ, "evaluate_unit", _fake_record)
    settings = occ.paper_grade_settings(replicates=2)
    occ.run_shard(tmp_path, settings, 0, 1, limit=1)
    path = lsc.shard_path(tmp_path / "paper_grade", 0, 1)
    record = json.loads(path.read_text())
    record["settings_hash"] = "other"
    path.write_text(json.dumps(record) + "\n")
    with pytest.raises(SystemExit, match="written under other settings"):
        occ.merge(tmp_path, "paper_grade")


# ── 3.2 Reproduction ─────────────────────────────────────────────────────────


def _committed_rate_rows(settings) -> list[dict]:
    rows = []
    for (column, mode, effect), statistics in occ.committed_operating(settings).items():
        for statistic, entry in statistics.items():
            rows.append({"column": column, "trajectory_mode": mode, "effect_size": effect, "statistic": statistic,
                         "arm": "scores", "n_available": entry["n_available"], "n_rejected": entry["n_rejected"]})
    return rows


def test_committed_operating_resolves_every_row(paper_settings) -> None:
    committed = occ.committed_operating(paper_settings)
    assert len(committed) == 65
    assert sum(len(statistics) for statistics in committed.values()) == 65 * 3
    assert {column for column, _, _ in committed} == set(leak.column_cells(occ.PAPER_GRADE_PROFILES["span_geometry"]))
    # The anchors of the two siblings are the same counts.
    magnitude = grade.load_rates(occ._resolve(paper_settings["parents"]["operating:span_magnitude"]["path"]))
    for column in occ.PAPER_GRADE_COLUMNS:
        assert magnitude[(column, "none", 0.0)] == committed[(column, "none", 0.0)]
    # The four-stage orientation cells carry the failed units.
    assert committed[("n200_four_stages", "orientation", 0.5)]["angle"]["n_available"] < 500


def test_paper_grade_reproduction_passes_and_fails(paper_settings) -> None:
    rows = _committed_rate_rows(paper_settings)
    result = occ.paper_grade_reproduction(rows, paper_settings)
    assert result["passed"] and result["rows_compared"] == 195 and result["differing"] == []
    rows[7] = {**rows[7], "n_rejected": rows[7]["n_rejected"] + 1}
    rows.pop(0)
    failed = occ.paper_grade_reproduction(rows, paper_settings)
    assert not failed["passed"] and len(failed["differing"]) == 2
    changed = [row for row in failed["differing"] if row["reproduced"] is not None][0]
    assert changed["reproduced"]["n_rejected"] == changed["committed"]["n_rejected"] + 1
    assert [row for row in failed["differing"] if row["reproduced"] is None]


def test_weak_reproduction_passes_and_fails(weak_settings) -> None:
    committed = lsc._read_csv(occ.SPAN_SPLIT_DIR / "production_table.csv")
    rows = [
        {"column": row["design"], "cell": row["cell"], "statistic": row["statistic"], "arm": "scores",
         "n_available": int(row["count"]), "n_rejected": int(row["rejections"])}
        for row in committed
        if row["design"] == "seaad_weak" and row["cell"] in ("anchor", "translation")
    ]
    result = occ.weak_reproduction(rows, weak_settings, committed)
    assert result["passed"] and result["rows_compared"] == 6
    rows[0] = {**rows[0], "n_rejected": rows[0]["n_rejected"] + 1}
    failed = occ.weak_reproduction(rows, weak_settings, committed)
    assert not failed["passed"] and len(failed["differing"]) == 1


# ── 4.1 Decision ─────────────────────────────────────────────────────────────


def _synthetic_tables(*, overrides=None, paired_overrides=None, weak_overrides=None):
    """Rate, paired and cell rows at the null for controls and off-targets and high target power."""

    overrides, paired_overrides, weak_overrides = overrides or {}, paired_overrides or {}, weak_overrides or {}
    rates, paired, cells = [], [], []
    keys = [(column, "none", 0.0) for column in occ.PAPER_GRADE_COLUMNS]
    keys += [(column, "translation", effect) for column in occ.PAPER_GRADE_COLUMNS for effect in occ.GEOMETRY_EFFECTS]
    keys += occ._target_cells()
    for column, mode, effect in keys:
        for statistic in occ.STATISTICS:
            target = occ.TARGETS.get(mode) == statistic
            for arm in occ.ARMS:
                n_rejected = overrides.get((column, mode, effect, statistic, arm), 450 if target else 25)
                rates.append({"column": column, "trajectory_mode": mode, "effect_size": effect,
                              "statistic": statistic, "arm": arm, "n_available": 500, "n_rejected": n_rejected})
            b, c = paired_overrides.get((column, mode, effect, statistic), (0, 0))
            two_sided, loss = occ.paired_p_values(b, c)
            paired.append({"column": column, "trajectory_mode": mode, "effect_size": effect, "statistic": statistic,
                           "n_paired": 500, "rate_scores": 0.9, "rate_orthonormal": 0.9 + (c - b) / 500,
                           "difference": (c - b) / 500, "b": b, "c": c, "p_two_sided": two_sided, "p_loss": loss})
        for arm in occ.ARMS:
            cells.append({"column": column, "trajectory_mode": mode, "effect_size": effect, "arm": arm,
                          "median_population_share": -0.005 if arm == "scores" else 0.0005,
                          "median_pooled_eigengap": 0.2, "median_angle_null_q95": 30.0})
    weak = []
    for column, mode, effect in occ.WEAK_CONTROLS:
        for statistic in occ.STATISTICS:
            for arm in occ.ARMS:
                weak.append({"column": column, "trajectory_mode": mode, "effect_size": effect, "statistic": statistic,
                             "arm": arm, "n_available": 500,
                             "n_rejected": weak_overrides.get((mode, statistic, arm), 25)})
    return {"rates": rates, "paired": paired, "cells": cells}, {"rates": weak, "paired": [], "cells": []}


PASSING = {"paper_grade": {"passed": True}, "weak": {"passed": True}}
CONSISTENT = {"paper_grade": {"inconsistent_units": 0}, "weak": {"inconsistent_units": 0}}


def _decide(paper, weak, reproduction=PASSING, accounting=CONSISTENT):
    return occ.decision(paper, weak, reproduction=reproduction, accounting_by_profile=accounting)


def test_adopt_orthonormal() -> None:
    result = _decide(*_synthetic_tables())
    assert result["verdict"] == "adopt_orthonormal" and result["failing_criteria"] == []
    assert result["criteria"]["O1"]["families"]["paper_grade"]["m"] == 60
    assert result["criteria"]["O1"]["families"]["weak"]["m"] == 6
    assert result["criteria"]["O2"]["m"] == 90
    assert result["criteria"]["O3"]["m"] == 45 and result["criteria"]["O3"]["allowed_losses"] == 3
    assert result["consequence"] == occ.CONSEQUENCES["adopt_orthonormal"]
    watched = result["descriptive"]["watched_cells"]
    assert [(row["column"], row["statistic"]) for row in watched] == [("n200", "shape"), ("n200_four_stages", "shape")]
    assert all("p_two_sided" in row for row in watched)
    residual = result["descriptive"]["population_residual_shape_e0.5"]
    assert len(residual) == 5 and residual[0]["share_scores"] == -0.005 and residual[0]["share_orthonormal"] == 0.0005


def test_size_excess_names_the_reading() -> None:
    paper, weak = _synthetic_tables(overrides={("seaad", "none", 0.0, "angle", "orthonormal"): 100})
    result = _decide(paper, weak)
    assert result["verdict"] == "size_excess" and result["failing_criteria"] == ["O1"]
    named = [row for row in result["failures"] if row["criterion"] == "O1"]
    assert named[0]["column"] == "seaad" and named[0]["statistic"] == "angle" and named[0]["family"] == "paper_grade"


def test_the_weak_family_alone_fails_size() -> None:
    paper, weak = _synthetic_tables(weak_overrides={("translation", "delta", "orthonormal"): 100})
    result = _decide(paper, weak)
    assert result["verdict"] == "size_excess"
    assert not result["criteria"]["O1"]["families"]["weak"]["holds"]
    assert result["criteria"]["O1"]["families"]["paper_grade"]["holds"]


def test_specificity_excess() -> None:
    paper, weak = _synthetic_tables(overrides={("n200", "orientation", 0.5, "shape", "orthonormal"): 100})
    result = _decide(paper, weak)
    assert result["verdict"] == "specificity_excess" and result["failing_criteria"] == ["O2"]


def test_power_loss_names_the_cell() -> None:
    paper, weak = _synthetic_tables(paired_overrides={("n80", "shape", 0.25, "shape"): (12, 0)})
    result = _decide(paper, weak)
    assert result["verdict"] == "power_loss" and result["failing_criteria"] == ["O3"]
    named = [row for row in result["failures"] if row["criterion"] == "O3"]
    assert len(named) == 1 and named[0]["column"] == "n80" and (named[0]["b"], named[0]["c"]) == (12, 0)
    assert named[0]["status"] == "below_ceiling"
    assert {"rate_scores", "rate_orthonormal"} <= set(named[0])


def test_losses_within_the_allowed_count_hold() -> None:
    # p_loss = 0.5**6 ≈ 0.016: below 0.025, above 0.05/45.
    cells = {key: (6, 0) for key in [("n80", "shape", 0.25, "shape"), ("n200", "shape", 0.25, "shape"),
                                     ("seaad", "shape", 0.25, "shape")]}
    assert _decide(*_synthetic_tables(paired_overrides=cells))["verdict"] == "adopt_orthonormal"
    cells[("n80", "shape", 0.5, "shape")] = (6, 0)
    assert _decide(*_synthetic_tables(paired_overrides=cells))["verdict"] == "power_loss"


def test_size_excess_takes_precedence_over_a_failing_power() -> None:
    paper, weak = _synthetic_tables(
        overrides={("n80", "none", 0.0, "delta", "orthonormal"): 100},
        paired_overrides={("n80", "shape", 0.25, "shape"): (12, 0)},
    )
    result = _decide(paper, weak)
    assert result["verdict"] == "size_excess" and result["failing_criteria"] == ["O1", "O3"]


def test_the_gate_still_records_every_criterion() -> None:
    paper, weak = _synthetic_tables(overrides={("seaad", "none", 0.0, "angle", "orthonormal"): 100})
    failed = {**PASSING, "weak": {"passed": False}}
    result = _decide(paper, weak, reproduction=failed)
    assert result["verdict"] == "reproduction_failed"
    assert set(result["criteria"]) == {"O1", "O2", "O3"} and result["failing_criteria"] == ["O1"]
    inconsistent = {**CONSISTENT, "paper_grade": {"inconsistent_units": 2}}
    gated = _decide(*_synthetic_tables(), accounting=inconsistent)
    assert gated["verdict"] == "reproduction_failed"
    assert gated["gate"]["failures"] == [{"profile": "paper_grade", "reason": "2 inconsistent unit(s)"}]


def test_no_discordant_units_give_p_one() -> None:
    assert occ.paired_p_values(0, 0) == (1.0, 1.0)
    two_sided, loss = occ.paired_p_values(3, 1)
    assert loss == pytest.approx(binom.sf(2, 4, 0.5)) and two_sided == pytest.approx(min(1.0, 2 * binom.cdf(1, 4, 0.5)))
    assert occ.paired_p_values(0, 5)[1] == 1.0


def test_allowed_count_at_45() -> None:
    assert occ.loss_allowed(45) == int(binom.ppf(0.95, 45, 0.025)) == 3


def test_an_unavailable_reading_fails_its_criterion() -> None:
    paper, weak = _synthetic_tables()
    paper["paired"] = [row for row in paper["paired"]
                       if (row["column"], row["trajectory_mode"], row["effect_size"]) != ("n80", "magnitude", 0.05)]
    result = _decide(paper, weak)
    assert result["verdict"] == "power_loss"
    unavailable = [row for row in result["failures"] if row.get("status") == "unavailable"]
    assert unavailable == [{"criterion": "O3", "column": "n80", "trajectory_mode": "magnitude", "effect_size": 0.05,
                            "statistic": "delta", "status": "unavailable"}]
    paper, weak = _synthetic_tables()
    weak["rates"] = [row for row in weak["rates"] if row["trajectory_mode"] != "translation"]
    assert _decide(paper, weak)["verdict"] == "size_excess"


def test_decide_reads_both_merged_profiles(tmp_path, monkeypatch) -> None:
    monkeypatch.setattr(occ, "evaluate_unit", _fake_record)
    for settings in (occ.paper_grade_settings(replicates=1), occ.weak_settings(replicates=1)):
        occ.run_shard(tmp_path, settings, 0, 1)
        occ.merge(tmp_path, settings["profile"])
    result = occ.decide(tmp_path)
    assert result["verdict"] == "reproduction_failed"
    assert result["complete"]["paper_grade"]["complete_units"] == 65
    assert result["complete"]["weak"]["complete_units"] == 2
    written = json.loads((tmp_path / occ.DECISION_NAME).read_text())
    assert written["rule_version"] == 1 and written["verdict"] == "reproduction_failed"
    assert set(written["inputs"]) >= {"paper_grade/rate_table.csv", "weak/reproduction.json"}
