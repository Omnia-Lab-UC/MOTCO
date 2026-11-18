#!/usr/bin/env python3
"""Diagnosis of the ADNI latent plasmode's ``gate_failed`` verdict, from the run's own records (no new units).

The driver behind ``openspec/changes/adni-latent-plasmode-diagnosis/`` (spec
``stage-program-capture-diagnosis``). The latent-plasmode run
(``adni_latent_plasmode_check.py``, 2026-11-18) read ``gate_failed``: the bias of
``pls_cv``'s total capture exceeded ±0.05 in its calibration. Its rule names the
consequence: a diagnosis of the estimator comes before any reading. This driver
re-reads the run's merged records and asks four questions of the designated
reading (``pls_cv``, total capture), recording the same quantities for every
space and stage contrast:

- **F1, the parts.** Are a unit's numerator ``N`` and denominator ``D`` unbiased
  for their truths ``TN``, ``TD`` in the calibration? The mean of ``N − TN`` and
  ``D − TD`` over a cell's units, with a standard error from the synthetic
  cohorts' means.
- **F2, the aggregate.** Is the ratio of sums over all of a cell's units (24
  cohorts) on its truth ratio? A cluster bootstrap over cohorts gives its
  standard error and interval.
- **F3, where the spread lives.** The covariance of ``(N, D)`` within cohorts
  ``W`` (pooled over cohorts) and between them ``B`` (the covariance of the cohort
  means minus ``W / m``); the between share of the variance of a cohort's mean
  ``D``. When ``B`` dominates, more units or splits of one cohort do not narrow it.
- **F4, ADNI's own reading.** For each point, ADNI's mean ``(N, D)`` over its 500
  units and a covariance ``V = f·B + W_ADNI / 500``, with ``B`` from a calibration
  strength (a stand-in: ADNI's between-cohort covariance cannot be observed from
  one cohort) and ``f`` 1 (raw) or ADNI's within-cohort variance of ``D`` over the
  calibration's (scaled). Fieller's 95 % interval for the ratio is bounded only
  when ``t = D̄ / √V_DD`` exceeds 1.96.

Beside them, per cell: the per-cohort ratio's mean error (the gate's bias), its
median error, SD and robust SD (IQR / 1.349), the cohorts whose mean ``D`` is
below a quarter of the cell's, and the delta-method bias and SD at the
calibration's units per cohort, at 500 units and as units grow without bound.

**Bindings.** The parent's five committed files (both profiles' settings, the
calibration table, the cohort's capture table and the decision) are bound by
sha256 in this driver's ``settings.json`` and re-checked by ``run`` and
``decide``; a changed file is refused by name. Every record must carry its
profile's committed settings hash, and the unit counts must equal the
settings'. The merged records stay where the run wrote them (on ing) and their
sha256 is recorded.

**Reproduction.** Before any reading, the parent driver's own table builders
regenerate ``calibration/calibration_table.csv`` and
``cohort/capture_table.csv`` from the records; each must equal the committed
file byte for byte, or the verdict is ``reproduction_failed``.

**The reading** (rule version 1; written after an exploratory pass on the same
records, recorded in the change's ``notes.md``, so it is a reading, not a
predeclared test): ``reproduction_failed``; else ``estimator_biased`` (F1 or F2
fails: |z| ≥ 3 at a gate cell); else ``resolvable`` (F4 fails: ADNI's interval is
bounded at some point under the primary stand-in, strength 0.67 raw); else
``units_limited`` (F3 fails: the between share of ``var(D̄)`` at strength 0.67 is
below 0.5 at some point); else ``information_limited``. The parent's verdict is
not changed.

Outputs under ``<out>``: ``settings.json``, ``component_table.csv`` (F1, F2),
``spread_table.csv`` (F3 and the per-cohort spread), ``adni_table.csv`` (F4),
``reproduction.json`` and, from ``decide``, ``adni_latent_plasmode_diagnosis.json``.
Every output holds aggregates only: no participant identifier, no per-unit and no
per-cohort value.

Example (on ing, where the run's records are):

    RECORDS=/home1/tgonza/MOTCO-adni-latent-run/out
    OUT=$HOME/motco-runs/adni-latent-plasmode-diagnosis
    python scripts/adni_latent_plasmode_diagnosis.py run --records-dir $RECORDS --out-dir $OUT
    python scripts/adni_latent_plasmode_diagnosis.py decide --out-dir $OUT
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import sys
import tempfile
from pathlib import Path
from typing import Any, Mapping, Sequence

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import adni_latent_plasmode_check as alp  # noqa: E402  (readings, table builders, profiles)
import latent_space_comparison as lsc  # noqa: E402  (settings hashing and files, JSON and CSV helpers)
import production_latent_space_check as plsc  # noqa: E402  (sha256, path display and resolution)

sys.path.pop(0)

DRIVER_NAME = "adni-latent-plasmode-diagnosis"
DRIVER_VERSION = 1
RULE_VERSION = 1
BASE_SEED = 2042

REPO = Path(__file__).resolve().parents[1]
PARENT_DIR = REPO / "results" / "adni-latent-plasmode-2026-11-18"
PARENT_FILES: dict[str, str] = {
    "calibration_settings": "calibration/settings.json",
    "calibration_table": "calibration/calibration_table.csv",
    "cohort_settings": "cohort/settings.json",
    "cohort_captures": "cohort/capture_table.csv",
    "decision": alp.DECISION_NAME,
}
REPRODUCED: dict[str, tuple[str, str]] = {
    # name → (profile, committed file inside the parent directory)
    "calibration_table": (alp.CALIBRATION, "calibration/calibration_table.csv"),
    "cohort_captures": (alp.COHORT, "cohort/capture_table.csv"),
}

DESIGNATED = alp.DESIGNATED
READINGS: tuple[str, ...] = alp.TRUE_READINGS
GATE_STRENGTHS: tuple[float, ...] = alp.GATE_STRENGTHS
NULL_STRENGTH = alp.NULL_STRENGTH
INTERVAL_STRENGTH = alp.INTERVAL_STRENGTH
STAND_INS: tuple[tuple[float, str], ...] = (
    (NULL_STRENGTH, "raw"), (NULL_STRENGTH, "scaled"), (INTERVAL_STRENGTH, "raw"), (INTERVAL_STRENGTH, "scaled"))
PRIMARY_STAND_IN: tuple[float, str] = (INTERVAL_STRENGTH, "raw")
ADNI_UNITS_FOR_V = 500
BOOTSTRAP = 2000
Z95 = 1.959963984540054
BIAS_Z_LIMIT = 3.0
BETWEEN_SHARE_LIMIT = 0.5
LOW_DENOMINATOR_SHARE = 0.25
VERDICTS: tuple[str, ...] = (
    "reproduction_failed", "estimator_biased", "resolvable", "units_limited", "information_limited")

SETTINGS_NAME = "settings.json"
COMPONENT_NAME = "component_table.csv"
SPREAD_NAME = "spread_table.csv"
ADNI_NAME = "adni_table.csv"
REPRODUCTION_NAME = "reproduction.json"
DECISION_NAME = "adni_latent_plasmode_diagnosis.json"

COMPONENT_COLUMNS: tuple[str, ...] = (
    "strength", "point", "space", "reading", "cohorts", "units", "mean_N", "mean_TN", "bias_N", "bias_N_se", "z_N",
    "mean_D", "mean_TD", "bias_D", "bias_D_se", "z_D", "grand_ratio", "grand_truth", "grand_bias", "grand_bias_se",
    "z_grand", "grand_bias_lower", "grand_bias_upper")
SPREAD_COLUMNS: tuple[str, ...] = (
    "strength", "point", "space", "reading", "cohorts", "units_per_cohort", "mean_N", "mean_D", "W_NN", "W_ND",
    "W_DD", "B_NN", "B_ND", "B_DD", "between_share_N", "between_share_D", "cohort_D_sd", "cohort_D_cv",
    "low_denominator_cohorts", "estimate_sd", "estimate_robust_sd", "cohort_bias", "cohort_bias_se",
    "cohort_error_median", "delta_bias_m", "delta_bias_m500", "delta_sd_m", "delta_sd_m500", "delta_sd_inf")
ADNI_COLUMNS: tuple[str, ...] = (
    "point", "space", "reading", "units", "mean_N", "mean_D", "estimate", "W_NN", "W_ND", "W_DD", "stand_in_strength",
    "stand_in", "within_ratio_D", "V_NN", "V_ND", "V_DD", "t", "fieller", "lower", "upper")


# ── Settings and bindings ────────────────────────────────────────────────────


def parents(parent_dir: Path | None = None) -> dict[str, dict[str, str]]:
    directory = PARENT_DIR if parent_dir is None else Path(parent_dir)
    return {name: {"path": plsc._display_path(directory / relative),
                   "sha256": plsc.sha256_file(directory / relative)}
            for name, relative in PARENT_FILES.items()}


def build_settings(parent_dir: Path | None = None) -> dict[str, Any]:
    return {
        "driver": DRIVER_NAME, "driver_version": DRIVER_VERSION, "rule_version": RULE_VERSION,
        "base_seed": BASE_SEED, "bootstrap": BOOTSTRAP, "designated_space": DESIGNATED, "readings": list(READINGS),
        "gate_strengths": list(GATE_STRENGTHS), "stand_ins": [list(entry) for entry in STAND_INS],
        "primary_stand_in": list(PRIMARY_STAND_IN), "adni_units_for_v": ADNI_UNITS_FOR_V, "z": Z95,
        "bias_z_limit": BIAS_Z_LIMIT, "between_share_limit": BETWEEN_SHARE_LIMIT,
        "low_denominator_share": LOW_DENOMINATOR_SHARE, "parents": parents(parent_dir),
    }


def verify_parents(settings: Mapping[str, Any]) -> None:
    """Refuse a parent file that is missing or whose hash differs from the settings', naming it."""

    for name, entry in settings["parents"].items():
        path = plsc._resolve(entry["path"])
        if not path.exists():
            raise SystemExit(f"Parent file {entry['path']} ({name}) is missing; refusing.")
        if plsc.sha256_file(path) != entry["sha256"]:
            raise SystemExit(f"Parent file {entry['path']} ({name}) differs from the hash bound in the settings; "
                             "refusing.")


def parent_settings(settings: Mapping[str, Any], profile: str) -> dict[str, Any]:
    path = plsc._resolve(settings["parents"][f"{profile}_settings"]["path"])
    stored = json.loads(path.read_text(encoding="utf-8"))
    if stored.get("hash") != lsc.settings_hash(stored["settings"]):
        raise SystemExit(f"{path} is corrupt: its hash does not match its settings.")
    return dict(stored["settings"])


def _seed(*parts: Any) -> int:
    encoded = json.dumps([DRIVER_NAME, *parts], separators=(",", ":")).encode("utf-8")
    return int(hashlib.sha256(encoded).hexdigest()[:8], 16) & 0x7FFFFFFF


def load_records(records_dir: Path, profile: str, profile_settings: Mapping[str, Any]
                 ) -> tuple[list[dict[str, Any]], dict[str, Any]]:
    """A profile's merged records, refused unless every record carries the committed settings hash and every
    expected unit is present once."""

    path = Path(records_dir) / profile / "merged.jsonl"
    if not path.exists():
        raise SystemExit(f"{path} does not exist; point --records-dir at the run's output directory.")
    records = alp.read_records([path])
    expected_hash = lsc.settings_hash(profile_settings)
    stale = sorted({alp.record_key(record) for record in records if record.get("settings_hash") != expected_hash})
    if stale:
        raise SystemExit(f"{len(stale)} record(s) of {path} were written under other settings than the committed "
                         f"{profile} settings, e.g. {stale[:3]}; refusing.")
    keys = [alp.record_key(record) for record in records]
    expected = {alp.unit_key(unit) for unit in alp.enumerate_units(profile_settings)}
    if len(keys) != len(set(keys)) or set(keys) != expected:
        raise SystemExit(f"{path} holds {len(set(keys))} distinct of {len(expected)} expected units "
                         f"({len(keys) - len(set(keys))} duplicates); refusing.")
    return records, {"path": str(path), "sha256": plsc.sha256_file(path), "records": len(records)}


# ── Reproduction ─────────────────────────────────────────────────────────────


def reproduce(records: Mapping[str, Sequence[Mapping[str, Any]]], profiles: Mapping[str, Mapping[str, Any]],
              settings: Mapping[str, Any]) -> dict[str, Any]:
    """The parent's table builders on the records, against the committed files byte for byte."""

    builders = {"calibration_table": (alp.calibration_rows, alp.CALIBRATION_COLUMNS),
                "cohort_captures": (alp.capture_rows, alp.CAPTURE_COLUMNS)}
    rows = []
    with tempfile.TemporaryDirectory() as scratch:
        for name, (profile, _) in REPRODUCED.items():
            build, columns = builders[name]
            target = Path(scratch) / f"{name}.csv"
            lsc._write_csv(target, build(records[profile], profiles[profile]), columns)
            observed = plsc.sha256_file(target)
            committed = settings["parents"][name]["sha256"]
            rows.append({"file": settings["parents"][name]["path"], "committed": committed, "observed": observed,
                         "equal": observed == committed})
    return {"rule": "the parent driver's table builders, run on the merged records, write the committed "
                    "calibration table and the cohort's capture table byte for byte",
            "rows": rows, "passed": all(row["equal"] for row in rows)}


# ── Arithmetic ───────────────────────────────────────────────────────────────


def variance_components(values: np.ndarray, groups: np.ndarray) -> dict[str, Any]:
    """``values``: units × k; ``groups``: one label per unit. ``W`` the pooled within-group covariance (ddof 1 per
    group), ``C`` the covariance of the group means, ``m`` the mean group size and ``B = C − W / m``."""

    values = np.asarray(values, dtype=float)
    groups = np.asarray(groups)
    labels = np.unique(groups)
    means = np.stack([values[groups == label].mean(axis=0) for label in labels])
    centred = np.concatenate([values[groups == label] - values[groups == label].mean(axis=0) for label in labels])
    within = centred.T @ centred / float(len(values) - len(labels))
    between_means = np.atleast_2d(np.cov(means.T, ddof=1)) if len(labels) > 1 else np.full_like(within, np.nan)
    m = float(np.mean([np.sum(groups == label) for label in labels]))
    return {"W": within, "C": between_means, "B": between_means - within / m, "m": m, "means": means,
            "mu": values.mean(axis=0), "labels": labels}


def delta_ratio(mu: Sequence[float], cov: np.ndarray) -> tuple[float | None, float | None]:
    """First-order SD and bias of ``N̄ / D̄`` with mean ``mu`` and covariance ``cov`` of ``(N̄, D̄)``."""

    mean_n, mean_d = float(mu[0]), float(mu[1])
    if mean_d == 0:
        return None, None
    r = mean_n / mean_d
    variance = (cov[0, 0] - 2.0 * r * cov[0, 1] + r * r * cov[1, 1]) / mean_d**2
    bias = (r * cov[1, 1] - cov[0, 1]) / mean_d**2
    return (math.sqrt(variance) if variance > 0 else None), float(bias)


def fieller(mean_n: float, mean_d: float, cov: np.ndarray, z: float = Z95) -> dict[str, Any]:
    """Fieller's interval for ``N/D``: ``bounded`` (lower, upper), ``exclusive`` (the line outside (lower, upper))
    or ``whole_line``; ``t`` is ``D̄ / √V_DD``."""

    a = mean_d**2 - z**2 * cov[1, 1]
    b = -2.0 * (mean_n * mean_d - z**2 * cov[0, 1])
    c = mean_n**2 - z**2 * cov[0, 0]
    disc = b * b - 4.0 * a * c
    t = mean_d / math.sqrt(cov[1, 1]) if cov[1, 1] > 0 else None
    if disc < 0:
        return {"kind": "whole_line", "lower": None, "upper": None, "t": t}
    roots = sorted(((-b - math.sqrt(disc)) / (2 * a), (-b + math.sqrt(disc)) / (2 * a))) if a != 0 else [None, None]
    return {"kind": "bounded" if a > 0 else "exclusive", "lower": roots[0], "upper": roots[1], "t": t}


def _pairs(records: Sequence[Mapping[str, Any]], space: str, reading: str, n_stages: int, *, truth: bool = False
           ) -> np.ndarray:
    reader = alp.truth_reading if truth else alp.unit_reading
    out = [reader(record, space, reading, n_stages) for record in records]
    if any(pair is None for pair in out):
        raise ValueError(f"A unit has no {'truth' if truth else 'reading'} for {space}/{reading}.")
    return np.asarray(out, dtype=float)


def _robust_sd(values: np.ndarray) -> float:
    q1, q3 = np.quantile(values, [0.25, 0.75])
    return float((q3 - q1) / 1.349)


def calibration_cell(records: Sequence[Mapping[str, Any]], space: str, reading: str, n_stages: int,
                     rng: np.random.Generator | None) -> tuple[dict[str, Any], dict[str, Any] | None]:
    """One strength × point cell: its spread row and, when the cell has a truth, its component row."""

    groups = np.asarray([int(record["cohort"]) for record in records])
    observed = _pairs(records, space, reading, n_stages)
    parts = variance_components(observed, groups)
    labels, mu, w, b, c = parts["labels"], parts["mu"], parts["W"], parts["B"], parts["C"]
    m = parts["m"]
    numerators, denominators = parts["means"][:, 0], parts["means"][:, 1]
    estimates = numerators / denominators
    spread: dict[str, Any] = {
        "cohorts": int(len(labels)), "units_per_cohort": m, "mean_N": float(mu[0]), "mean_D": float(mu[1]),
        "W_NN": float(w[0, 0]), "W_ND": float(w[0, 1]), "W_DD": float(w[1, 1]),
        "B_NN": float(b[0, 0]), "B_ND": float(b[0, 1]), "B_DD": float(b[1, 1]),
        "between_share_N": float(b[0, 0] / c[0, 0]) if c[0, 0] > 0 else None,
        "between_share_D": float(b[1, 1] / c[1, 1]) if c[1, 1] > 0 else None,
        "cohort_D_sd": float(math.sqrt(c[1, 1])),
        "cohort_D_cv": float(math.sqrt(c[1, 1]) / abs(mu[1])) if mu[1] != 0 else None,
        "low_denominator_cohorts": int(np.sum(denominators < LOW_DENOMINATOR_SHARE * mu[1])),
        "estimate_sd": float(np.std(estimates, ddof=1)), "estimate_robust_sd": _robust_sd(estimates),
    }
    spread["delta_sd_m"], spread["delta_bias_m"] = delta_ratio(mu, w / m + b)
    spread["delta_sd_m500"], spread["delta_bias_m500"] = delta_ratio(mu, w / 500.0 + b)
    spread["delta_sd_inf"], _ = delta_ratio(mu, b)
    if not all(record.get("truth") for record in records):
        return spread, None
    truths = _pairs(records, space, reading, n_stages, truth=True)
    truth_means = np.stack([truths[groups == label].mean(axis=0) for label in labels])
    errors = estimates - truth_means[:, 0] / truth_means[:, 1]
    spread.update({"cohort_bias": float(np.mean(errors)),
                   "cohort_bias_se": float(np.std(errors, ddof=1) / math.sqrt(len(labels))),
                   "cohort_error_median": float(np.median(errors))})
    component: dict[str, Any] = {"cohorts": int(len(labels)), "units": int(len(records))}
    for index, name in enumerate(("N", "D")):
        difference = observed[:, index] - truths[:, index]
        cohort_means = np.asarray([difference[groups == label].mean() for label in labels])
        se = float(np.std(cohort_means, ddof=1) / math.sqrt(len(labels)))
        component.update({f"mean_{name}": float(observed[:, index].mean()),
                          f"mean_T{name}": float(truths[:, index].mean()), f"bias_{name}": float(difference.mean()),
                          f"bias_{name}_se": se, f"z_{name}": float(difference.mean() / se) if se > 0 else None})
    grand = float(observed[:, 0].sum() / observed[:, 1].sum())
    grand_truth = float(truths[:, 0].sum() / truths[:, 1].sum())
    draws = []
    if rng is not None:
        members = {label: np.flatnonzero(groups == label) for label in labels}
        for _ in range(BOOTSTRAP):
            pick = np.concatenate([members[label] for label in rng.choice(labels, size=len(labels), replace=True)])
            draws.append(observed[pick, 0].sum() / observed[pick, 1].sum()
                         - truths[pick, 0].sum() / truths[pick, 1].sum())
    se = float(np.std(draws, ddof=1)) if len(draws) > 1 else None
    component.update({
        "grand_ratio": grand, "grand_truth": grand_truth, "grand_bias": grand - grand_truth, "grand_bias_se": se,
        "z_grand": (grand - grand_truth) / se if se else None,
        "grand_bias_lower": float(np.quantile(draws, 0.025)) if draws else None,
        "grand_bias_upper": float(np.quantile(draws, 0.975)) if draws else None,
    })
    return spread, component


def calibration_tables(records: Sequence[Mapping[str, Any]], profile_settings: Mapping[str, Any]
                       ) -> tuple[list[dict[str, Any]], list[dict[str, Any]]]:
    n_stages = len(profile_settings["contract"]["stage_order"])
    spread_rows, component_rows = [], []
    for strength, point, members in alp._cells(alp.usable(records), profile_settings):
        for space in profile_settings["spaces"]:
            for reading in READINGS:
                rng = np.random.default_rng(_seed(BASE_SEED, "bootstrap", strength, point, space, reading))
                spread, component = calibration_cell(members, space, reading, n_stages, rng)
                keys = {"strength": strength, "point": point, "space": space, "reading": reading}
                spread_rows.append({**keys, **spread})
                if component is not None:
                    component_rows.append({**keys, **component})
    return component_rows, spread_rows


def _row(rows: Sequence[Mapping[str, Any]], **keys: Any) -> Mapping[str, Any] | None:
    for row in rows:
        if all((math.isclose(float(row[key]), float(value)) if key == "strength" else str(row[key]) == str(value))
               for key, value in keys.items()):
            return row
    return None


def adni_tables(records: Sequence[Mapping[str, Any]], profile_settings: Mapping[str, Any],
                spread_rows: Sequence[Mapping[str, Any]]) -> list[dict[str, Any]]:
    """Per point, space, reading and stand-in: ADNI's means, within covariance, ``V`` and Fieller's interval."""

    n_stages = len(profile_settings["contract"]["stage_order"])
    rows = []
    for _, point, members in alp._cells(alp.usable(records), profile_settings):
        for space in profile_settings["spaces"]:
            for reading in READINGS:
                pairs = _pairs(members, space, reading, n_stages)
                mean_n, mean_d = (float(value) for value in pairs.mean(axis=0))
                within = np.cov(pairs.T, ddof=1)
                base = {"point": point, "space": space, "reading": reading, "units": int(len(pairs)),
                        "mean_N": mean_n, "mean_D": mean_d, "estimate": float(pairs[:, 0].sum() / pairs[:, 1].sum()),
                        "W_NN": float(within[0, 0]), "W_ND": float(within[0, 1]), "W_DD": float(within[1, 1])}
                for strength, kind in STAND_INS:
                    spread = _row(spread_rows, strength=strength, point=point, space=space, reading=reading)
                    if spread is None:
                        raise ValueError(f"No calibration row at strength {strength}, {point}, {space}, {reading}.")
                    between = np.array([[spread["B_NN"], spread["B_ND"]], [spread["B_ND"], spread["B_DD"]]])
                    ratio = float(within[1, 1] / spread["W_DD"])
                    factor = 1.0 if kind == "raw" else ratio
                    cov = factor * between + within / float(ADNI_UNITS_FOR_V)
                    interval = fieller(mean_n, mean_d, cov)
                    rows.append({**base, "stand_in_strength": strength, "stand_in": kind, "within_ratio_D": ratio,
                                 "V_NN": float(cov[0, 0]), "V_ND": float(cov[0, 1]), "V_DD": float(cov[1, 1]),
                                 "t": interval["t"], "fieller": interval["kind"], "lower": interval["lower"],
                                 "upper": interval["upper"]})
    return rows


# ── The reading (rule version 1) ─────────────────────────────────────────────


RULE: dict[str, Any] = {
    "version": RULE_VERSION,
    "written": "after an exploratory pass on the same records (notes.md, 2026-11-18); a reading, not a "
               "predeclared test, and it does not change the parent's verdict",
    "designated": f"{DESIGNATED}, total capture",
    "F1": f"|z| of the mean N - TN and of D - TD below {BIAS_Z_LIMIT:g} at every gate strength and point",
    "F2": f"|z| of the grand ratio's bias (cluster bootstrap over cohorts) below {BIAS_Z_LIMIT:g} at every gate "
          "strength and point",
    "F3": f"the between share of var(cohort mean D) at strength {INTERVAL_STRENGTH:g} at least "
          f"{BETWEEN_SHARE_LIMIT:g} at every point",
    "F4": f"ADNI's 95 % Fieller interval unbounded at every point under the primary stand-in (strength "
          f"{PRIMARY_STAND_IN[0]:g}, {PRIMARY_STAND_IN[1]}); the other stand-ins are recorded",
    "verdict": "reproduction_failed; else estimator_biased (F1 or F2 fails); else resolvable (F4 fails); else "
               "units_limited (F3 fails); else information_limited",
}
CONSEQUENCES: dict[str, str] = {
    "reproduction_failed": "The records do not regenerate the committed tables; nothing is read until they do.",
    "estimator_biased": "The estimator's parts or its aggregate are biased; the estimator needs repair before any "
                        "rule version 2.",
    "resolvable": "A rule version 2 with an interval for a ratio (Fieller's) could read a capture at the points "
                  "where ADNI's interval is bounded; it is its own change, written before it reads.",
    "units_limited": "The spread is mostly within cohorts; a run with more units per point could narrow it.",
    "information_limited": "The estimator is sound, but its spread is a property of the one cohort and ADNI's "
                           "program energy does not clear it: no 95 % interval for the capture is bounded at any "
                           "point. No rule version 2 is proposed on ADNI; the capture waits for a larger cohort "
                           "or several cohorts, and the parent's gate_failed stands.",
}


def findings(component_rows: Sequence[Mapping[str, Any]], spread_rows: Sequence[Mapping[str, Any]],
             adni_rows: Sequence[Mapping[str, Any]], points: Sequence[str]) -> dict[str, Any]:
    def designated(rows: Sequence[Mapping[str, Any]], **keys: Any) -> Mapping[str, Any] | None:
        return _row(rows, space=DESIGNATED, reading="total", **keys)

    f1, f2, f3, f4 = {}, {}, {}, {}
    for strength in GATE_STRENGTHS:
        for point in points:
            row = designated(component_rows, strength=strength, point=point)
            z = [None if row is None else _num(row.get(f"z_{name}")) for name in ("N", "D")]
            f1[f"s{alp.strength_label(strength)}/{point}"] = {
                "z_N": z[0], "z_D": z[1], "passed": all(v is not None and abs(v) < BIAS_Z_LIMIT for v in z)}
            zg = None if row is None else _num(row.get("z_grand"))
            f2[f"s{alp.strength_label(strength)}/{point}"] = {
                "grand_bias": None if row is None else _num(row.get("grand_bias")), "z_grand": zg,
                "passed": zg is not None and abs(zg) < BIAS_Z_LIMIT}
    for point in points:
        row = designated(spread_rows, strength=INTERVAL_STRENGTH, point=point)
        share = None if row is None else _num(row.get("between_share_D"))
        f3[point] = {"between_share_D": share, "between_share_N": None if row is None else
                     _num(row.get("between_share_N")),
                     "passed": share is not None and share >= BETWEEN_SHARE_LIMIT}
        by_stand_in = {}
        for strength, kind in STAND_INS:
            entry = _row(adni_rows, point=point, space=DESIGNATED, reading="total", stand_in_strength=strength,
                         stand_in=kind)
            by_stand_in[f"s{alp.strength_label(strength)}/{kind}"] = None if entry is None else {
                "t": _num(entry.get("t")), "fieller": str(entry.get("fieller"))}
        primary = by_stand_in[f"s{alp.strength_label(PRIMARY_STAND_IN[0])}/{PRIMARY_STAND_IN[1]}"]
        f4[point] = {"stand_ins": by_stand_in,
                     "passed": primary is not None and primary["fieller"] != "bounded",
                     "stand_ins_agree": all(entry is not None and entry["fieller"] != "bounded"
                                            for entry in by_stand_in.values())}
    out = {}
    for name, cells in (("F1", f1), ("F2", f2), ("F3", f3), ("F4", f4)):
        out[name] = {"rule": RULE[name], "cells": cells, "passed": bool(cells) and all(c["passed"]
                                                                                        for c in cells.values())}
    return out


def _num(value: Any) -> float | None:
    if value is None or value == "":
        return None
    return float(value)


def verdict_of(reproduced: bool, found: Mapping[str, Mapping[str, Any]]) -> str:
    if not reproduced:
        return "reproduction_failed"
    if not (found["F1"]["passed"] and found["F2"]["passed"]):
        return "estimator_biased"
    if not found["F4"]["passed"]:
        return "resolvable"
    if not found["F3"]["passed"]:
        return "units_limited"
    return "information_limited"


# ── Commands ─────────────────────────────────────────────────────────────────


def run(records_dir: Path, out_dir: Path, parent_dir: Path | None = None) -> dict[str, Any]:
    settings = build_settings(parent_dir)
    lsc.ensure_settings(out_dir, settings)
    verify_parents(settings)
    profiles = {profile: parent_settings(settings, profile) for profile in alp.PROFILES}
    records, sources = {}, {}
    for profile in alp.PROFILES:
        records[profile], sources[profile] = load_records(records_dir, profile, profiles[profile])
    reproduced = reproduce(records, profiles, settings)
    reproduced["records"] = {profile: {"sha256": entry["sha256"], "records": entry["records"]}
                             for profile, entry in sources.items()}
    lsc._write_json(out_dir / REPRODUCTION_NAME, reproduced)
    component_rows, spread_rows = calibration_tables(records[alp.CALIBRATION], profiles[alp.CALIBRATION])
    adni_rows = adni_tables(records[alp.COHORT], profiles[alp.COHORT], spread_rows)
    lsc._write_csv(out_dir / COMPONENT_NAME, component_rows, COMPONENT_COLUMNS)
    lsc._write_csv(out_dir / SPREAD_NAME, spread_rows, SPREAD_COLUMNS)
    lsc._write_csv(out_dir / ADNI_NAME, adni_rows, ADNI_COLUMNS)
    return {"reproduction_passed": reproduced["passed"], "component_rows": len(component_rows),
            "spread_rows": len(spread_rows), "adni_rows": len(adni_rows)}


def decide(out_dir: Path) -> dict[str, Any]:
    settings = lsc.load_settings(out_dir)
    verify_parents(settings)
    reproduced = json.loads((out_dir / REPRODUCTION_NAME).read_text(encoding="utf-8"))
    tables = {name: lsc._read_csv(out_dir / name) for name in (COMPONENT_NAME, SPREAD_NAME, ADNI_NAME)}
    points = alp.point_names(parent_settings(settings, alp.COHORT))
    found = findings(tables[COMPONENT_NAME], tables[SPREAD_NAME], tables[ADNI_NAME], points)
    verdict = verdict_of(bool(reproduced.get("passed")), found)
    parent = json.loads(plsc._resolve(settings["parents"]["decision"]["path"]).read_text(encoding="utf-8"))
    payload = {
        "verdict": verdict, "consequence": CONSEQUENCES[verdict], "rule": RULE, "rule_version": RULE_VERSION,
        "parent_verdict": parent.get("verdict"), "reproduction": {"passed": reproduced.get("passed")},
        "findings": found, "settings_hash": lsc.settings_hash(settings),
        "inputs": {name: plsc.sha256_file(out_dir / name) for name in
                   (REPRODUCTION_NAME, COMPONENT_NAME, SPREAD_NAME, ADNI_NAME)},
    }
    lsc._write_json(out_dir / DECISION_NAME, payload)
    return payload


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    commands = parser.add_subparsers(dest="command", required=True)
    run_parser = commands.add_parser("run", help="re-read the run's records and write the tables")
    run_parser.add_argument("--records-dir", type=Path, required=True,
                            help="the latent-plasmode run's output directory (holds <profile>/merged.jsonl)")
    run_parser.add_argument("--out-dir", type=Path, required=True)
    decide_parser = commands.add_parser("decide", help="apply the reading rule to the tables")
    decide_parser.add_argument("--out-dir", type=Path, required=True)
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    if args.command == "run":
        print(json.dumps(run(args.records_dir, args.out_dir), indent=1))
    else:
        payload = decide(args.out_dir)
        print(json.dumps({key: payload[key] for key in ("verdict", "consequence", "parent_verdict")}, indent=1))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
