"""``scripts/adni_latent_plasmode_diagnosis.py``: arithmetic, bindings, the reproduction and the reading.

No export is read. Records are hand-built in the latent-plasmode driver's own format from its own product
functions (Gaussian half configurations around a known stage configuration, with cohort-level and unit-level
noise), and the parent's committed files are written from them by the parent's own table builders.
"""

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


ald = _script("adni_latent_plasmode_diagnosis")
alp, lsc = ald.alp, ald.lsc

STAGES = ["CN", "MCI", "AD"]
FEATURES = 6
SPACE_DIMENSIONS = {"pls_cv": 2, "predictive": 1}
POINTS = [{"name": "n80", "participants": 80, "subsampled": True, "table": [[1, 1, 1], [1, 1, 1]]},
          {"name": "n200", "participants": 200, "subsampled": True, "table": [[1, 1, 1], [1, 1, 1]]}]
COHORTS, UNITS, REPLICATES = 3, 4, 6


# ── Hand-built records ───────────────────────────────────────────────────────


def _settings(profile: str) -> dict:
    base = {"profile": profile, "spaces": sorted(SPACE_DIMENSIONS), "contract": {"stage_order": STAGES},
            "points": POINTS}
    if profile == alp.CALIBRATION:
        return {**base, "strengths": [0.0, 0.67, 1.0], "cohorts": COHORTS, "units": UNITS}
    return {**base, "replicates": REPLICATES}


def _centred(rng: np.random.Generator) -> np.ndarray:
    z = rng.standard_normal((len(STAGES), FEATURES))
    return z - z.mean(axis=0)


def _record(unit, settings: dict, configuration: np.ndarray, cohort_noise: tuple[np.ndarray, np.ndarray],
            rng: np.random.Generator) -> dict:
    bases = {}
    for name, dimension in SPACE_DIMENSIONS.items():
        q, _ = np.linalg.qr(rng.standard_normal((FEATURES, dimension)))
        bases[name] = q
    m1 = configuration + cohort_noise[0] + 0.3 * rng.standard_normal(configuration.shape)
    m2 = configuration + cohort_noise[1] + 0.3 * rng.standard_normal(configuration.shape)
    products = alp.mean_products([alp.split_products(m1, m2, bases)])
    naive = alp.naive_products(m1, bases)
    record = {
        "status": "completed", "point": unit.point, "replicate": int(unit.replicate),
        "settings_hash": lsc.settings_hash(settings),
        "cross_products": {
            "total": {**{key: products["total"][key] for key in ("cross", "gram", "noise")},
                      "naive_remainder": naive["total"], "naive_sample": naive["total"]},
            "spaces": {name: {"dimension": SPACE_DIMENSIONS[name], **products["spaces"][name],
                              "naive_remainder": naive["spaces"][name], "naive_sample": naive["spaces"][name]}
                       for name in bases},
        },
    }
    if unit.strength is not None:
        record.update({"strength": float(unit.strength), "cohort": int(unit.cohort),
                       "truth": alp.truth_products(configuration, bases)})
    return json.loads(json.dumps(record, default=lambda value: np.asarray(value).tolist()))


def _records(profile: str, settings: dict, seed: int = 7) -> list[dict]:
    rng = np.random.default_rng(seed)
    program = _centred(rng)
    noise = {cohort: (0.5 * rng.standard_normal(program.shape), 0.5 * rng.standard_normal(program.shape))
             for cohort in range(COHORTS + 1)}
    out = []
    for unit in alp.enumerate_units(settings):
        strength = 0.8 if unit.strength is None else float(unit.strength)
        cohort = COHORTS if unit.cohort is None else int(unit.cohort)
        out.append(_record(unit, settings, strength * program, noise[cohort], rng))
    return out


def _write_lines(path: Path, records: list[dict]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("".join(json.dumps(record) + "\n" for record in records), encoding="utf-8")


@pytest.fixture()
def world(tmp_path: Path) -> dict:
    """A parent results directory, the run's records directory and the diagnosis's output directory."""

    parent, records_dir, out = tmp_path / "parent", tmp_path / "records", tmp_path / "out"
    profiles = {profile: _settings(profile) for profile in alp.PROFILES}
    records = {profile: _records(profile, profiles[profile]) for profile in alp.PROFILES}
    for profile in alp.PROFILES:
        lsc.ensure_settings(parent / profile, profiles[profile])
        _write_lines(records_dir / profile / "merged.jsonl", records[profile])
    lsc._write_csv(parent / "calibration" / "calibration_table.csv",
                   alp.calibration_rows(records[alp.CALIBRATION], profiles[alp.CALIBRATION]), alp.CALIBRATION_COLUMNS)
    lsc._write_csv(parent / "cohort" / "capture_table.csv",
                   alp.capture_rows(records[alp.COHORT], profiles[alp.COHORT]), alp.CAPTURE_COLUMNS)
    (parent / alp.DECISION_NAME).write_text(json.dumps({"verdict": "gate_failed"}), encoding="utf-8")
    return {"parent": parent, "records": records_dir, "out": out, "profiles": profiles, "data": records}


# ── Arithmetic ───────────────────────────────────────────────────────────────


def test_variance_components_match_the_hand_computation() -> None:
    rng = np.random.default_rng(1)
    groups = np.repeat(np.arange(5), 8)
    values = rng.standard_normal((5, 2))[groups] + 0.4 * rng.standard_normal((40, 2))
    parts = ald.variance_components(values, groups)
    within = sum(np.cov(values[groups == g].T, ddof=1) * 7 for g in range(5)) / (40 - 5)
    means = np.stack([values[groups == g].mean(axis=0) for g in range(5)])
    np.testing.assert_allclose(parts["W"], within)
    np.testing.assert_allclose(parts["C"], np.cov(means.T, ddof=1))
    np.testing.assert_allclose(parts["B"], np.cov(means.T, ddof=1) - within / 8)
    assert parts["m"] == 8


def test_delta_ratio_formula() -> None:
    cov = np.array([[0.04, 0.01], [0.01, 0.09]])
    sd, bias = ald.delta_ratio([1.0, 2.0], cov)
    r = 0.5
    assert sd == pytest.approx(math.sqrt((0.04 - 2 * r * 0.01 + r * r * 0.09) / 4))
    assert bias == pytest.approx((r * 0.09 - 0.01) / 4)
    assert ald.delta_ratio([1.0, 2.0], np.zeros((2, 2)))[0] is None


def test_fieller_kinds() -> None:
    small = np.diag([0.01, 0.01])
    bounded = ald.fieller(1.0, 2.0, small)
    assert bounded["kind"] == "bounded" and bounded["lower"] < 0.5 < bounded["upper"]
    assert bounded["t"] == pytest.approx(20.0)
    assert ald.fieller(3.0, 0.1, np.diag([0.01, 0.25]))["kind"] == "exclusive"
    assert ald.fieller(0.1, 0.1, np.diag([0.25, 0.25]))["kind"] == "whole_line"


def test_the_verdict_order() -> None:
    def found(**passed: bool) -> dict:
        return {name: {"passed": passed.get(name, True)} for name in ("F1", "F2", "F3", "F4")}

    assert ald.verdict_of(False, found()) == "reproduction_failed"
    assert ald.verdict_of(True, found(F1=False)) == "estimator_biased"
    assert ald.verdict_of(True, found(F2=False)) == "estimator_biased"
    assert ald.verdict_of(True, found(F4=False, F3=False)) == "resolvable"
    assert ald.verdict_of(True, found(F3=False)) == "units_limited"
    assert ald.verdict_of(True, found()) == "information_limited"
    assert set(ald.CONSEQUENCES) == set(ald.VERDICTS)


# ── Run and decide ───────────────────────────────────────────────────────────


def test_run_reproduces_and_decides(world: dict) -> None:
    summary = ald.run(world["records"], world["out"], world["parent"])
    assert summary["reproduction_passed"]
    cells = 3 * len(POINTS) * len(SPACE_DIMENSIONS) * len(ald.READINGS)
    assert summary["spread_rows"] == cells
    assert summary["component_rows"] == 2 * len(POINTS) * len(SPACE_DIMENSIONS) * len(ald.READINGS)  # no truth at 0
    assert summary["adni_rows"] == len(POINTS) * len(SPACE_DIMENSIONS) * len(ald.READINGS) * len(ald.STAND_INS)
    payload = ald.decide(world["out"])
    assert payload["verdict"] in ald.VERDICTS and payload["verdict"] != "reproduction_failed"
    assert payload["parent_verdict"] == "gate_failed"
    assert set(payload["findings"]) == {"F1", "F2", "F3", "F4"}
    assert set(payload["findings"]["F4"]["cells"]) == {"n80", "n200"}
    for cell in payload["findings"]["F4"]["cells"].values():  # every stand-in found, strength 0 included
        assert set(cell["stand_ins"]) == {"s0/raw", "s0/scaled", "s0.67/raw", "s0.67/scaled"}
        assert all(entry is not None and entry["fieller"] in ("bounded", "exclusive", "whole_line")
                   for entry in cell["stand_ins"].values())


def test_component_rows_are_differences_of_means(world: dict) -> None:
    ald.run(world["records"], world["out"], world["parent"])
    row = next(r for r in lsc._read_csv(world["out"] / ald.COMPONENT_NAME)
               if r["space"] == "pls_cv" and r["reading"] == "total" and r["point"] == "n80" and r["strength"] == "1")
    members = [r for r in world["data"][alp.CALIBRATION] if r["point"] == "n80" and r["strength"] == 1.0]
    observed = np.array([alp.unit_reading(r, "pls_cv", "total", 3) for r in members])
    truth = np.array([alp.truth_reading(r, "pls_cv", "total", 3) for r in members])
    assert float(row["bias_N"]) == pytest.approx(float(np.mean(observed[:, 0] - truth[:, 0])), rel=1e-5)
    assert float(row["grand_truth"]) == pytest.approx(truth[:, 0].sum() / truth[:, 1].sum(), rel=1e-5)


def test_cohort_bias_equals_the_parents_calibration_table(world: dict) -> None:
    ald.run(world["records"], world["out"], world["parent"])
    ours = {(r["strength"], r["point"], r["space"], r["reading"]): r["cohort_bias"]
            for r in lsc._read_csv(world["out"] / ald.SPREAD_NAME)}
    theirs = lsc._read_csv(world["parent"] / "calibration" / "calibration_table.csv")
    compared = 0
    for row in theirs:
        key = (row["strength"], row["point"], row["space"], row["reading"])
        if key in ours and row["bias"] != "":
            assert ours[key] == row["bias"]
            compared += 1
    assert compared > 0


def test_a_changed_parent_file_is_refused_by_name(world: dict) -> None:
    ald.run(world["records"], world["out"], world["parent"])
    table = world["parent"] / "calibration" / "calibration_table.csv"
    table.write_text(table.read_text(encoding="utf-8") + "\n", encoding="utf-8")
    with pytest.raises(SystemExit, match="calibration_table.csv"):
        ald.decide(world["out"])
    with pytest.raises(SystemExit, match="parents.calibration_table.sha256"):
        ald.run(world["records"], world["out"], world["parent"])


def test_records_under_other_settings_are_refused(world: dict) -> None:
    path = world["records"] / "cohort" / "merged.jsonl"
    records = world["data"][alp.COHORT]
    records[0] = {**records[0], "settings_hash": "0" * 64}
    _write_lines(path, records)
    with pytest.raises(SystemExit, match="other settings"):
        ald.run(world["records"], world["out"], world["parent"])


def test_a_missing_unit_is_refused(world: dict) -> None:
    _write_lines(world["records"] / "calibration" / "merged.jsonl", world["data"][alp.CALIBRATION][1:])
    with pytest.raises(SystemExit, match="expected units"):
        ald.run(world["records"], world["out"], world["parent"])


def test_records_that_do_not_reproduce_fail_the_reading(world: dict) -> None:
    records = world["data"][alp.COHORT]
    first = records[0]
    first["cross_products"]["total"]["cross"] = float(first["cross_products"]["total"]["cross"]) + 1.0
    _write_lines(world["records"] / "cohort" / "merged.jsonl", records)
    assert not ald.run(world["records"], world["out"], world["parent"])["reproduction_passed"]
    assert ald.decide(world["out"])["verdict"] == "reproduction_failed"


def test_outputs_hold_aggregates_only(world: dict) -> None:
    ald.run(world["records"], world["out"], world["parent"])
    ald.decide(world["out"])
    for name, columns in ((ald.COMPONENT_NAME, ald.COMPONENT_COLUMNS), (ald.SPREAD_NAME, ald.SPREAD_COLUMNS),
                          (ald.ADNI_NAME, ald.ADNI_COLUMNS)):
        rows = lsc._read_csv(world["out"] / name)
        assert list(rows[0]) == list(columns)  # aggregate columns only: no per-unit or per-cohort field

    def lists(value) -> list:
        if isinstance(value, list):
            return [value] + [item for entry in value for item in lists(entry)]
        if isinstance(value, dict):
            return [item for entry in value.values() for item in lists(entry)]
        return []

    for name in (ald.DECISION_NAME, ald.REPRODUCTION_NAME):
        payload = json.loads((world["out"] / name).read_text(encoding="utf-8"))
        assert all(len(entry) <= 4 for entry in lists(payload))
