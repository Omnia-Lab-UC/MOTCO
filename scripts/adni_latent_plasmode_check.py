#!/usr/bin/env python3
"""The latent plasmode on ADNI: how much of the cohort's real stage program the production span holds (no R).

The driver behind ``openspec/changes/adni-latent-plasmode-capture/`` (spec
``cohort-stage-program-capture``, roadmap open item 2). It estimates the share
of ADNI's covariate-adjusted, group-blind stage program that lies inside a
fitted measurement space, from two independent held-out halves, calibrates the
estimator on synthetic cohorts with the truth known, and reads both under a
decision fixed before any unit runs. No statistic of group difference is read:
no RRPP draw, no plant and no analysis under the real sex labels. No package
code changes; every fit is one ``analyze_trajectories`` call with
``permutations=0``.

**The estimand** (design D1). For a space with orthonormal basis ``Q`` fitted
on the unit's sample, in the fitted model's input coordinates, and ``M`` the
cohort's covariate-adjusted stage configuration in those coordinates, centred
over stages (stages × features), the capture is ``tr(M·Q·Qᵀ·Mᵀ) / tr(M·Mᵀ)``.
The configuration of a set of participants is the stage coefficients of
``features ~ stage indicators + centred covariate columns`` in the sample's
pooled-standardized coordinates (``adni_gain_specificity_check``'s
``adjusted_configuration``, imported); the group label is in no fit.

**Three disjoint parts** (design D2). The sample ``A`` is the anchor driver's
subsample (``adni_anchor_check``, imported unchanged: base seed 2039, within
sex × stage cells, labels permuted within stage), so a unit at ``n80``,
``n200`` or ``n500`` is the anchor run's unit and one at ``n750`` the gain
run's. The participants outside ``A`` are split within stage into two halves
(``split_within_stage``, imported) 20 times; each half's configuration ``m1``,
``m2`` is fitted on that half alone. With the halves independent of each other
and of ``Q``, ``tr(m1·Q·Qᵀ·m2ᵀ)`` and ``tr(m1·m2ᵀ)`` are unbiased for the
program's energy inside the space and in total.

**The estimate is a ratio of sums** (design D4): per point, space and reading,
the sum over units of the numerator over the sum of the denominator, never
clipped; per-unit ratios are quantiles for description only.

**Readings** (design D5), per space: ``total``; the stage contrasts
``net_displacement`` (the unit contrast last stage minus first) and
``interior`` (the rest of the centred stage space), whose energies sum to the
total, from the K × K cross-Gram matrices ``(m1·m2ᵀ + m2·m1ᵀ)/2`` and their
in-space counterparts; ``noise``, the share of the halves' difference
``m1 − m2`` (no program, only estimation noise); and the naive shares of the
whole remainder's configuration (``naive_remainder``) and of the sample's own
(``naive_sample``), item 22's two readings. Beside them: the program's energy,
the eigenvalues of the pooled cross-Gram matrix and the halves' correlation.

**Spaces** (design D6): ``pls_cv`` (the designated space: production's fitted
span at the case study's configuration, cross-validated once per unit),
``predictive`` (the predictive subspace of the same model), ``pls_rank2``,
``pls_rank8`` and ``pls_rank19`` (forced fits) and ``pca`` (the leading right
singular vectors of the sample's model inputs at the unit's CV rank).

**Profiles** (``--profile``; each writes under ``<out>/<profile>/``):

- ``cohort`` — ADNI's blocks: 500 units at each of ``n80``, ``n200``,
  ``n500`` and ``n750`` (2,000).
- ``calibration`` (design D7) — synthetic cohorts that keep every
  participant's stage, sex and covariates and draw features as the full
  cohort's grand mean + ``s`` × its centred stage configuration + its
  covariate part + Gaussian noise with its residual covariance, at strengths
  ``s`` 0, 0.67 and 1; 24 cohorts per strength and 25 units per cohort and
  point (7,200), each evaluated as a cohort unit is, with its true capture
  recorded. It gives the estimator's bias, the cohort-level standard deviation
  of the pooled estimate and, at ``s = 0``, the null level of each energy.

**Own seeds** (sha256 of ``[driver name, 2041, …]``, truncated as
``plasmode_check._seed`` does): ``reference_split`` per unit and split index;
in the calibration also ``cohort_noise`` per strength and cohort, and
``subsample`` and ``labels`` per unit (with the redraw attempt).

**Checks on every unit** (design D9; a failing unit is ``inconsistent`` and
enters no reading): ``cell_table`` (the anchor driver's, imported), ``parts``
(the sample and each split's halves are pairwise disjoint and cover the
cohort), ``basis`` (each basis orthonormal to 1e-8 at its expected dimension,
every fit in the same input coordinates), ``gram`` (each Gram matrix's trace
equals its scalar cross-product to 1e-8) and, in the calibration, ``truth``
(the planted configuration is recovered exactly from the noise-free features).

**Parents** (sha256 in both profiles' settings; every subcommand refuses a
changed file by name): the ADNI anchor run's ``point_table.csv`` (the CV-rank
counts at 80, 200 and 500) and the gain run's plants ``cell_table.csv`` (the
naive captures' medians at 750). The export's five files are bound the same
way, and ``merge`` and ``decide`` re-hash them.

Subcommands (``--export-dir`` required by every subcommand, no default):

- ``preflight`` — ``settings.json`` and ``preflight.json``. Cohort: the anchor
  driver's pre-flight units (cohort replicates 0 and 1 select rank 8 at AUROC
  0.599514 and equal the forced call), and the full cohort configuration's
  split-half reliability recorded. Calibration: at ``s = 1`` with the noise
  scaled to zero, every space's estimated capture equals its true capture to
  1e-8. ``run`` refuses to start without a pass.
- ``run`` — ``records/shard_<i>_of_<N>.jsonl``; round-robin shards, resume of
  complete units, refusal of other settings naming the keys.
- ``merge`` — ``merged.jsonl`` (keep off the repository), ``capture_table.csv``,
  ``point_table.csv``, ``accounting.json``; cohort: ``reproduction.json``;
  calibration: ``calibration_table.csv`` (one pooled estimate per synthetic
  cohort behind each standard deviation).
- ``decide`` — ``adni_latent_plasmode_decision.json`` (rule version 1) in
  ``<out>``: the verdict ``gate_failed`` | ``no_program`` | ``unresolved`` |
  ``partly_measured`` | ``measured`` (design D8), the contrast grades and what
  the verdict licenses.

Records hold seeds, ranks, scalars and stage × stage matrices only: no
participant identifier, label vector, subsample or split membership, stage
configuration or synthetic feature. The export is read in place, never
copied; there is no sbatch array. With ``--workers`` above 1, pin BLAS in the
shell.

Examples (the smoke: four units per point; two cohorts of two units per
strength and point):

    export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
    EXPORT=../adni_sexdiff/results/motco_input
    OUT=$HOME/motco-runs/adni-latent-plasmode-smoke
    python scripts/adni_latent_plasmode_check.py preflight --profile cohort --export-dir $EXPORT --out-dir $OUT --reps 4
    python scripts/adni_latent_plasmode_check.py run --profile cohort --export-dir $EXPORT --out-dir $OUT --reps 4 --workers 8
    python scripts/adni_latent_plasmode_check.py merge --profile cohort --export-dir $EXPORT --out-dir $OUT
    python scripts/adni_latent_plasmode_check.py preflight --profile calibration --export-dir $EXPORT --out-dir $OUT --cohorts 2 --units 2
    python scripts/adni_latent_plasmode_check.py run --profile calibration --export-dir $EXPORT --out-dir $OUT --cohorts 2 --units 2 --workers 8
    python scripts/adni_latent_plasmode_check.py merge --profile calibration --export-dir $EXPORT --out-dir $OUT
    python scripts/adni_latent_plasmode_check.py decide --export-dir $EXPORT --out-dir $OUT
"""  # noqa: E501

from __future__ import annotations

import argparse
import dataclasses
import hashlib
import json
import math
import sys
import warnings
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path
from time import perf_counter
from typing import Any, Iterable, Mapping, NamedTuple, Sequence

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
import adni_anchor_check as aac  # noqa: E402  (contract, units, seeds, analysis parameters, the cell-table check)
import adni_gain_specificity_check as gain  # noqa: E402  (scalers, the adjusted configuration, the split)
import latent_space_comparison as lsc  # noqa: E402  (settings hashing and files, JSON and CSV helpers)
import plasmode_check as pc  # noqa: E402  (the within-stage label permutation)
import predictive_subspace_check as psc  # noqa: E402  (record reading, check and failure records)
import production_latent_space_check as plsc  # noqa: E402  (parents, export binding, quantiles)

sys.path.pop(0)

DRIVER_NAME = "adni-latent-plasmode"
DRIVER_VERSION = 1
RULE_VERSION = 1
#: This driver's own seed family; no committed driver uses it.
BASE_SEED = 2041
PROFILES: tuple[str, ...] = ("cohort", "calibration")
COHORT, CALIBRATION = PROFILES

#: Design points (participants), drawn as the anchor driver draws them; the cohort point is the pre-flight's only.
POINT_SIZES: tuple[int, ...] = (80, 200, 500, 750)
REPLICATES = 500
SPLITS = 20
FIXED_RANKS: tuple[int, ...] = (2, 8, 19)
PREFLIGHT_REPLICATES: tuple[int, ...] = aac.PREFLIGHT_REPLICATES
CV: dict[str, Any] = dict(aac.CV)

STRENGTHS: tuple[float, ...] = (0.0, 0.67, 1.0)
NULL_STRENGTH = 0.0
INTERVAL_STRENGTH = 0.67
GATE_STRENGTHS: tuple[float, ...] = (0.67, 1.0)
CALIBRATION_COHORTS = 24
CALIBRATION_UNITS = 25
PREFLIGHT_STRENGTH = 1.0

TOLERANCE = 1e-8
SPACE_CV, SPACE_PREDICTIVE, SPACE_PCA = "pls_cv", "predictive", "pca"
DESIGNATED = SPACE_CV
CONTRASTS: tuple[str, ...] = ("net_displacement", "interior")
READINGS: tuple[str, ...] = ("total", *CONTRASTS, "noise", "naive_remainder", "naive_sample")
#: Readings with a true capture in the calibration.
TRUE_READINGS: tuple[str, ...] = ("total", *CONTRASTS)
CHECKS: dict[str, tuple[str, ...]] = {
    COHORT: ("cell_table", "parts", "basis", "gram"),
    CALIBRATION: ("cell_table", "parts", "basis", "gram", "truth"),
}

# ── Decision constants (design D8, rule version 1) ──────────────────────────

FAILURE_SHARE = 0.01
BIAS_LIMIT = 0.05
HALF_WIDTH_LIMIT = 0.15
INTERVAL_MULTIPLIER = 2.0
PRESENCE_MULTIPLIER = 2.0
CAPTURE_THRESHOLD = 0.5
PROGRAM_POINT = "n80"

REPO = Path(__file__).resolve().parents[1]
#: The ADNI anchor run (item 21) and the gain run's plants profile (item 22), bound by sha256.
ANCHOR_DIR = REPO / "results" / "adni-real-noise-anchor-2026-11-17"
GAIN_DIR = REPO / "results" / "adni-gain-specificity-2026-11-17" / "plants"
PARENT_FILES: dict[str, tuple[str, str]] = {
    "anchor_points": ("anchor", "point_table.csv"),
    "gain_cells": ("gain", "cell_table.csv"),
}
REPRODUCED_RANK_POINTS: tuple[str, ...] = ("n80", "n200", "n500")
REPRODUCED_NAIVE_POINT = "n750"
#: (reading, space) → the gain run's ``capture_median`` key.
NAIVE_PARENT_KEYS: dict[tuple[str, str], str] = {
    ("naive_remainder", SPACE_CV): "heldout:span",
    ("naive_remainder", SPACE_PREDICTIVE): "heldout:predictive",
    ("naive_sample", SPACE_CV): "sample:span",
    ("naive_sample", SPACE_PREDICTIVE): "sample:predictive",
}

PREFLIGHT_NAME = "preflight.json"
REPRODUCTION_NAME = "reproduction.json"
DECISION_NAME = "adni_latent_plasmode_decision.json"


# ── Spaces, settings and seeds ───────────────────────────────────────────────


def rank_space(rank: int) -> str:
    return f"pls_rank{int(rank)}"


def space_names(fixed_ranks: Sequence[int] = FIXED_RANKS) -> list[str]:
    return [SPACE_CV, SPACE_PREDICTIVE, *(rank_space(rank) for rank in fixed_ranks), SPACE_PCA]


def space_specs(fixed_ranks: Sequence[int] = FIXED_RANKS) -> dict[str, dict[str, Any]]:
    """Each space's fit, basis and dimension (design D6), in reading order."""

    specs: dict[str, dict[str, Any]] = {
        SPACE_CV: {"fit": "cross_validated", "rank": None, "dimension": "cv_rank",
                   "basis": "qr(x_rotations_) of the cross-validated fit (LinearLatentMap.from_pls)"},
        SPACE_PREDICTIVE: {"fit": "cross_validated", "rank": None, "dimension": "min(cv_rank, n_stages - 1)",
                           "basis": "qr(R·V) of the same model (LinearLatentMap.from_pls_predictive)"},
    }
    for rank in fixed_ranks:
        specs[rank_space(rank)] = {"fit": "forced", "rank": int(rank), "dimension": int(rank),
                                   "basis": "qr(x_rotations_) of the fit forced at the rank"}
    specs[SPACE_PCA] = {"fit": "none", "rank": None, "dimension": "cv_rank",
                        "basis": "the leading right singular vectors of the sample's model inputs"}
    return specs


def expected_dimension(spec: Mapping[str, Any], cv_rank: int, n_stages: int) -> int:
    if spec["dimension"] == "cv_rank":
        return int(cv_rank)
    if spec["dimension"] == "min(cv_rank, n_stages - 1)":
        return min(int(cv_rank), int(n_stages) - 1)
    return int(spec["dimension"])


def _seed(*parts: Any) -> int:
    encoded = json.dumps([DRIVER_NAME, *parts], separators=(",", ":")).encode("utf-8")
    return int(hashlib.sha256(encoded).hexdigest()[:8], 16) & 0x7FFFFFFF


def own_seed(settings: Mapping[str, Any], *parts: Any) -> int:
    return _seed(int(settings["own_base_seed"]), *parts)


def strength_label(strength: float) -> str:
    return f"{float(strength):g}"


def parents(anchor_dir: Path | None = None, gain_dir: Path | None = None) -> dict[str, dict[str, str]]:
    directories = {"anchor": Path(ANCHOR_DIR if anchor_dir is None else anchor_dir),
                   "gain": Path(GAIN_DIR if gain_dir is None else gain_dir)}
    return {role: plsc._parent(directories[source] / name) for role, (source, name) in PARENT_FILES.items()}


def _common_settings(
    cohort: aac.Cohort, profile: str, *, splits: int | None, point_sizes: Sequence[int] | None,
    fixed_ranks: Sequence[int] | None, cv: Mapping[str, Any] | None, anchor_dir: Path | None, gain_dir: Path | None,
) -> dict[str, Any]:
    contract = cohort.contract
    if contract.transform != "identity":
        raise SystemExit(f"The calibration draws raw features; the contract's transform is {contract.transform!r}, "
                         "not 'identity'.")
    splits = SPLITS if splits is None else int(splits)
    if splits < 1:
        raise SystemExit("--splits must be at least 1.")
    fixed_ranks = FIXED_RANKS if fixed_ranks is None else tuple(int(rank) for rank in fixed_ranks)
    return {
        "driver": DRIVER_NAME,
        "driver_version": DRIVER_VERSION,
        "rule_version": RULE_VERSION,
        "profile": profile,
        "contract": contract.to_dict(),
        "export": {"files": dict(cohort.files), "participants": int(len(cohort.metadata)),
                   "features": {name: int(frame.shape[1]) for name, frame in cohort.blocks.items()}},
        "parents": parents(anchor_dir, gain_dir),
        "own_base_seed": BASE_SEED,
        "max_attempts": aac.MAX_ATTEMPTS,
        # The cohort point is listed for the cohort pre-flight's units only; no unit of either profile is drawn there.
        "points": aac.design_points(cohort.table(), POINT_SIZES if point_sizes is None else tuple(point_sizes)),
        "splits": splits,
        "fixed_ranks": list(fixed_ranks),
        "spaces": space_specs(fixed_ranks),
        "designated_space": DESIGNATED,
        "configuration": "the stage coefficients of features ~ stage indicators + covariate columns (coded as "
                         "analyze_trajectories codes them, centred on the set's means), fitted on each set alone in "
                         "the sample's pooled-standardized coordinates, centred over stages and mapped to the "
                         "cross-validated model's input coordinates; the group label is in no fit",
        "split": "the participants outside the sample, split within stage into two halves "
                 "(adni_gain_specificity_check.split_within_stage; the first half takes each stage's floor)",
        "contrasts": {"net_displacement": "the unit contrast (last stage - first stage) / sqrt(2)",
                      "interior": "the centred stage space orthogonal to the net displacement"},
        "estimate": "per point, space and reading, the sum over units of the numerator over the sum of the "
                    "denominator, never clipped; each unit's numerator and denominator are means over its splits",
        "cohort_rank": int(contract.fidelity_rank),
        "permutations": 0,
        "analysis": {
            "cv": dict(CV if cv is None else cv),
            "min_components": aac.MIN_COMPONENTS,
            "orientation_estimator": aac.ORIENTATION_ESTIMATOR,
            "missing_covariates": aac.MISSING_COVARIATES,
            "permutations": 0,
            "residual_correction": None,
            "n_jobs": 1,
        },
        "fit": "one cross-validated fit per unit (pls_cv; the predictive subspace is the same model's); every other "
               "fit is forced at its rank; each is one analyze_trajectories call with permutations=0",
        "consistency_tolerance": TOLERANCE,
    }


def cohort_settings(
    cohort: aac.Cohort,
    *,
    replicates: int | None = None,
    splits: int | None = None,
    point_sizes: Sequence[int] | None = None,
    fixed_ranks: Sequence[int] | None = None,
    cv: Mapping[str, Any] | None = None,
    anchor_dir: Path | None = None,
    gain_dir: Path | None = None,
) -> dict[str, Any]:
    """The cohort profile's settings; unset arguments take the module's values at call time."""

    replicates = REPLICATES if replicates is None else int(replicates)
    if replicates < 1:
        raise SystemExit("--reps must be at least 1.")
    settings = _common_settings(cohort, COHORT, splits=splits, point_sizes=point_sizes, fixed_ranks=fixed_ranks,
                                cv=cv, anchor_dir=anchor_dir, gain_dir=gain_dir)
    settings.update({
        "base_seed": aac.BASE_SEED,
        "seed_derivation": f"the anchor driver's: sha256 of [{aac.DRIVER_NAME!r}, {aac.BASE_SEED}, point, replicate, "
                           "purpose(, attempt)]; purposes subsample, labels (the fit's seed is the unit's rrpp seed; "
                           "no permutation is drawn)",
        "own_seed_derivation": f"sha256 of [{DRIVER_NAME!r}, {BASE_SEED}, point, replicate, 'reference_split', split]",
        "replicates": replicates,
        "preflight_replicates": list(PREFLIGHT_REPLICATES),
        "preflight_reliability": f"the full cohort's split-half reliability: splits of sha256 of [{DRIVER_NAME!r}, "
                                 f"{BASE_SEED}, 'cohort', 0, 'reference_split', split], recorded and not judged",
    })
    return settings


def calibration_settings(
    cohort: aac.Cohort,
    *,
    cohorts: int | None = None,
    units: int | None = None,
    splits: int | None = None,
    strengths: Sequence[float] | None = None,
    point_sizes: Sequence[int] | None = None,
    fixed_ranks: Sequence[int] | None = None,
    cv: Mapping[str, Any] | None = None,
    anchor_dir: Path | None = None,
    gain_dir: Path | None = None,
) -> dict[str, Any]:
    """The calibration profile's settings; unset arguments take the module's values at call time."""

    cohorts = CALIBRATION_COHORTS if cohorts is None else int(cohorts)
    units = CALIBRATION_UNITS if units is None else int(units)
    if cohorts < 1 or units < 1:
        raise SystemExit("--cohorts and --units must be at least 1.")
    settings = _common_settings(cohort, CALIBRATION, splits=splits, point_sizes=point_sizes, fixed_ranks=fixed_ranks,
                                cv=cv, anchor_dir=anchor_dir, gain_dir=gain_dir)
    settings.update({
        "own_seed_derivation": f"sha256 of [{DRIVER_NAME!r}, {BASE_SEED}, 'calibration', strength, cohort, …]: "
                               "'cohort_noise' per strength and cohort; per unit [point, unit, purpose(, attempt or "
                               "split)] with purposes subsample, labels, reference_split (the fit's seed is the "
                               "subsample's first attempt; no permutation is drawn)",
        "strengths": [float(value) for value in (STRENGTHS if strengths is None else strengths)],
        "cohorts": cohorts,
        "units": units,
        "interval_strength": INTERVAL_STRENGTH,
        "null_strength": NULL_STRENGTH,
        "generator": {
            "fit": "least squares of the full cohort's raw joint features on stage indicators and covariate columns "
                   "(coded as analyze_trajectories codes them, centred on the cohort's means)",
            "features": "grand mean (the stage coefficients' mean) + strength × the centred stage coefficients at "
                        "the participant's stage + the participant's centred covariates × their coefficients + "
                        "N(0, S); every participant keeps their stage, sex and covariates",
            "covariance": "S = the residuals' cross-product / (participants - design rank), drawn through its "
                          "Cholesky factor",
            "truth": "strength × the centred stage coefficients (raw units), in each unit's sample scaler and model "
                     "input coordinates; no truth at strength 0",
        },
        "preflight": {"strength": PREFLIGHT_STRENGTH, "noise_scale": 0.0, "cohort": 0, "unit": 0,
                      "rule": "at strength 1 with the noise scaled to zero, every space's estimated capture equals "
                              "its true capture (total and both contrasts) to the tolerance on one unit per point",
                      "tolerance": TOLERANCE},
    })
    return settings


def point_names(settings: Mapping[str, Any]) -> list[str]:
    """The points units are drawn at (the cohort point is the pre-flight's only)."""

    return [str(point["name"]) for point in settings["points"] if point["subsampled"]]


def point_spec(settings: Mapping[str, Any], name: str) -> dict[str, Any]:
    return aac.point_spec(settings, name)


class Unit(NamedTuple):
    point: str
    replicate: int
    strength: float | None = None
    cohort: int | None = None


def unit_key(unit: Unit) -> str:
    if unit.strength is None:
        return f"{unit.point}/{int(unit.replicate)}"
    return f"s{strength_label(unit.strength)}/c{int(unit.cohort or 0)}/{unit.point}/{int(unit.replicate)}"


def record_unit(record: Mapping[str, Any]) -> Unit:
    strength = record.get("strength")
    return Unit(str(record["point"]), int(record["replicate"]), None if strength is None else float(strength),
                None if record.get("cohort") is None else int(record["cohort"]))


def record_key(record: Mapping[str, Any]) -> str:
    return unit_key(record_unit(record))


def enumerate_units(settings: Mapping[str, Any]) -> list[Unit]:
    """Replicate-major (cohort) or unit-major (calibration), so a partial run covers every cell."""

    points = point_names(settings)
    if settings["profile"] == COHORT:
        return [Unit(point, replicate) for replicate in range(int(settings["replicates"])) for point in points]
    return [Unit(point, unit, float(strength), cohort)
            for unit in range(int(settings["units"]))
            for cohort in range(int(settings["cohorts"]))
            for strength in settings["strengths"]
            for point in points]


def expected_units(settings: Mapping[str, Any]) -> dict[str, int]:
    counts: dict[str, int] = {}
    for unit in enumerate_units(settings):
        counts[unit.point] = counts.get(unit.point, 0) + 1
    return counts


def calibration_seed(settings: Mapping[str, Any], unit: Unit, purpose: str, *extra: Any) -> int:
    return own_seed(settings, "calibration", strength_label(float(unit.strength or 0.0)), int(unit.cohort or 0),
                    unit.point, int(unit.replicate), purpose, *extra)


def noise_seed(settings: Mapping[str, Any], strength: float, cohort: int) -> int:
    return own_seed(settings, "calibration", strength_label(strength), int(cohort), "cohort_noise")


def split_seeds(settings: Mapping[str, Any], unit: Unit) -> list[int]:
    count = int(settings["splits"])
    if unit.strength is None:
        return [own_seed(settings, unit.point, int(unit.replicate), "reference_split", index) for index in range(count)]
    return [calibration_seed(settings, unit, "reference_split", index) for index in range(count)]


# ── Configurations and parts (designs D1, D2) ───────────────────────────────


def centred(z: np.ndarray) -> np.ndarray:
    z = np.asarray(z, dtype=float)
    return z - z.mean(axis=0)


def configuration(features_std: np.ndarray, metadata: pd.DataFrame, contract: aac.CohortContract) -> np.ndarray:
    """A set's centred covariate-adjusted stage configuration (stages × features); raises when rank-deficient."""

    stages = gain.stage_positions(metadata, contract)
    return centred(gain.adjusted_configuration(features_std, stages, gain.covariate_columns(metadata, contract),
                                               len(contract.stage_order)))


def draw_split(remainder_metadata: pd.DataFrame, contract: aac.CohortContract, seed: int) -> tuple[np.ndarray, ...]:
    """Positions (into the remainder) of the two halves of one within-stage split."""

    stages = gain.stage_positions(remainder_metadata, contract)
    return gain.split_within_stage(stages, int(seed), len(contract.stage_order))


def _named_check(deviations: Mapping[str, float], tolerance: float) -> dict[str, Any]:
    """The worst deviation against the tolerance, naming every failing entry."""

    worst = max(deviations.values(), default=None)
    check = psc._check(worst, tolerance, applicable=bool(deviations))
    failing = {name: value for name, value in deviations.items()
               if not (math.isfinite(float(value)) and float(value) <= tolerance)}
    return {**check, "failing": failing}


def _renamed(check: Mapping[str, Any]) -> dict[str, Any]:
    """The anchor driver's check with its failing entries under ``failing``."""

    out = {key: value for key, value in check.items() if key != "failing_arms"}
    out["failing"] = check.get("failing_arms", {})
    return out


def check_parts(n_total: int, sample_rows: np.ndarray, halves: Sequence[tuple[np.ndarray, np.ndarray]]) -> dict:
    """Per split: the sample and the two halves are pairwise disjoint and cover the cohort's rows."""

    sample = np.asarray(sample_rows, dtype=int)
    failing: dict[str, Any] = {}
    for index, (first, second) in enumerate(halves):
        parts = [sample, np.asarray(first, dtype=int), np.asarray(second, dtype=int)]
        sizes = sum(part.size for part in parts)
        union = np.unique(np.concatenate(parts))
        overlap = int(sizes - union.size)
        uncovered = int(n_total - np.intersect1d(union, np.arange(int(n_total))).size)
        outside = int(np.setdiff1d(union, np.arange(int(n_total))).size)
        if overlap or uncovered or outside:
            failing[str(index)] = {"overlap": overlap, "uncovered": uncovered, "outside": outside}
    deviation = float(len(failing)) if halves else math.inf
    return {**psc._check(deviation, 0.0, applicable=True), "failing": failing}


# ── Spaces (design D6) ───────────────────────────────────────────────────────


def _filtered(call: Any, *args: Any) -> Any:
    """A fit with sklearn's 'residual is constant' warnings silenced (a noise-free cohort exhausts its rank)."""

    with warnings.catch_warnings():
        warnings.filterwarnings("ignore", message=".*residual is constant.*")
        return call(*args)


def fit_analysis(dataset: Any, settings: Mapping[str, Any], seed: int, rank: int | None) -> Any:
    """One ``analyze_trajectories`` call at the case study's configuration with no permutations."""

    params = aac.analysis_params(settings, int(seed), forced=rank, permutations=0, correction=None)
    return _filtered(aac._analyse, dataset, params)


def pca_basis(features: np.ndarray, model: Any, dimension: int) -> np.ndarray:
    """The leading right singular vectors of the sample's rows in the model's input coordinates."""

    inputs = (np.asarray(features, dtype=float) - np.asarray(model._x_mean)) / np.asarray(model._x_std)
    _, _, right = np.linalg.svd(inputs, full_matrices=False)
    return np.ascontiguousarray(right[: int(dimension)].T)


def space_bases(
    production: Any, forced: Mapping[int, Any], settings: Mapping[str, Any], n_stages: int
) -> tuple[dict[str, np.ndarray], dict[str, int], dict[str, Any]]:
    """Each space's orthonormal basis in the cross-validated model's input coordinates, its expected dimension,
    and each forced model's (centre, scale) for the coordinate check."""

    from motco.simulations.diagnostics import LinearLatentMap

    model = production.latent_artifacts.model
    cv_rank = int(production.record["latent"]["selected_lv"])
    features = production.latent_artifacts.features.to_numpy(dtype=float)
    specs = settings["spaces"]
    bases: dict[str, np.ndarray] = {}
    scalings: dict[str, Any] = {}
    for name, spec in specs.items():
        if name == SPACE_CV:
            bases[name] = LinearLatentMap.from_pls(model).basis
        elif name == SPACE_PREDICTIVE:
            bases[name] = LinearLatentMap.from_pls_predictive(model, aac.PREDICTIVE).basis
        elif name == SPACE_PCA:
            bases[name] = pca_basis(features, model, cv_rank)
        else:
            fitted = forced[int(spec["rank"])].latent_artifacts.model
            bases[name] = LinearLatentMap.from_pls(fitted).basis
            scalings[name] = fitted
    dimensions = {name: expected_dimension(spec, cv_rank, n_stages) for name, spec in specs.items()}
    return bases, dimensions, {"reference": model, "forced": scalings}


def check_basis(bases: Mapping[str, np.ndarray], dimensions: Mapping[str, int],
                scalings: Mapping[str, Any] | None = None) -> dict[str, Any]:
    """Each basis orthonormal at its expected dimension; each forced fit in the reference model's coordinates."""

    deviations: dict[str, float] = {}
    for name, basis in bases.items():
        basis = np.asarray(basis, dtype=float)
        if basis.ndim != 2 or basis.shape[1] != int(dimensions[name]):
            deviations[name] = math.inf
            continue
        gram = basis.T @ basis
        deviations[name] = float(np.max(np.abs(gram - np.eye(gram.shape[0])))) if gram.size else 0.0
    if scalings:
        reference = scalings["reference"]
        centre, scale = np.asarray(reference._x_mean, dtype=float), np.asarray(reference._x_std, dtype=float)
        for name, model in (scalings.get("forced") or {}).items():
            shift = float(np.max(np.abs((np.asarray(model._x_mean, dtype=float) - centre) / scale)))
            ratio = float(np.max(np.abs(np.asarray(model._x_std, dtype=float) / scale - 1.0)))
            deviations[name] = max(deviations.get(name, 0.0), shift, ratio)
    return _named_check(deviations, TOLERANCE)


# ── Cross-products (designs D4, D5) ──────────────────────────────────────────


def _sym(left: np.ndarray, right: np.ndarray) -> np.ndarray:
    return (left @ right.T + right @ left.T) / 2.0


def split_products(m1: np.ndarray, m2: np.ndarray, bases: Mapping[str, np.ndarray]) -> dict[str, Any]:
    """One split's cross-products: total and inside each space (scalars and stage × stage Gram matrices)."""

    m1, m2 = np.asarray(m1, dtype=float), np.asarray(m2, dtype=float)
    difference = m1 - m2
    norms = float(np.sum(m1**2)), float(np.sum(m2**2))
    total = {"cross": float(np.sum(m1 * m2)), "gram": _sym(m1, m2), "noise": float(np.sum(difference**2)),
             "correlation": float(np.sum(m1 * m2) / math.sqrt(norms[0] * norms[1])) if min(norms) > 0 else None}
    spaces = {}
    for name, basis in bases.items():
        a1, a2, d = m1 @ basis, m2 @ basis, difference @ basis
        spaces[name] = {"cross": float(np.sum(a1 * a2)), "gram": _sym(a1, a2), "noise": float(np.sum(d**2))}
    return {"total": total, "spaces": spaces}


def mean_products(products: Sequence[Mapping[str, Any]]) -> dict[str, Any]:
    """The products averaged over splits (the unit's numerators and denominators)."""

    def mean(values: Iterable[Any]) -> Any:
        values = list(values)
        if isinstance(values[0], np.ndarray):
            return np.mean(np.stack(values), axis=0)
        finite = [float(value) for value in values if value is not None]
        return float(np.mean(finite)) if finite else None

    total = {key: mean(entry["total"][key] for entry in products) for key in ("cross", "gram", "noise")}
    total["correlation"] = mean(entry["total"]["correlation"] for entry in products)
    spaces = {
        name: {key: mean(entry["spaces"][name][key] for entry in products) for key in ("cross", "gram", "noise")}
        for name in products[0]["spaces"]
    }
    return {"total": total, "spaces": spaces}


def naive_products(vectors: np.ndarray, bases: Mapping[str, np.ndarray]) -> dict[str, Any]:
    """One configuration's energy in total and inside each space (``span_share``'s arithmetic)."""

    rows = np.atleast_2d(np.asarray(vectors, dtype=float))
    return {"total": float(np.sum(rows**2)),
            "spaces": {name: float(np.sum((rows @ basis) ** 2)) for name, basis in bases.items()}}


def check_gram(products: Mapping[str, Any]) -> dict[str, Any]:
    """Each Gram matrix's trace against its scalar cross-product (relative to max(1, |scalar|))."""

    entries = {"total": products["total"], **{f"space:{name}": entry for name, entry in products["spaces"].items()}}
    deviations = {}
    for name, entry in entries.items():
        gram, cross = np.asarray(entry["gram"], dtype=float), float(entry["cross"])
        deviations[name] = float(abs(np.trace(gram) - cross) / max(1.0, abs(cross)))
    return _named_check(deviations, TOLERANCE)


def truth_products(planted: np.ndarray, bases: Mapping[str, np.ndarray]) -> dict[str, Any] | None:
    """The planted configuration's energy and Gram matrices, in total and inside each space; ``None`` when zero."""

    planted = np.asarray(planted, dtype=float)
    if not float(np.sum(planted**2)) > 0:
        return None
    spaces = {}
    for name, basis in bases.items():
        inside = planted @ basis
        spaces[name] = {"cross": float(np.sum(inside**2)), "gram": inside @ inside.T}
    return {"total": {"cross": float(np.sum(planted**2)), "gram": planted @ planted.T}, "spaces": spaces}


# ── Stage contrasts (design D5) ──────────────────────────────────────────────


def contrast_projectors(n_stages: int) -> dict[str, np.ndarray]:
    """The net displacement's and the interior's projectors in stage space; they sum to the centring projector."""

    k = int(n_stages)
    vector = np.zeros(k)
    vector[0], vector[-1] = -1.0, 1.0
    vector /= math.sqrt(2.0)
    net = np.outer(vector, vector)
    centring = np.eye(k) - np.full((k, k), 1.0 / k)
    return {"net_displacement": net, "interior": centring - net}


def contrast_energy(gram: Any, projector: np.ndarray) -> float:
    return float(np.trace(projector @ np.asarray(gram, dtype=float)))


# ── One unit (designs D2–D6, D9) ─────────────────────────────────────────────


def evaluate_parts(
    cohort: aac.Cohort, settings: Mapping[str, Any], built: Mapping[str, Any], *, seeds: Sequence[int],
    analysis_seed: int, timing: dict[str, float], planted_raw: np.ndarray | None = None,
    noise_free: np.ndarray | None = None,
) -> dict[str, Any]:
    """The unit's readings from its built sample: configurations, one CV fit, forced fits, cross-products, checks.

    Raises on a rank-deficient configuration or a failed fit (the caller records the unit as failed).
    """

    contract = cohort.contract
    n_stages = len(contract.stage_order)
    names = contract.block_names
    mark = perf_counter()
    scalers = gain.sample_scalers(built["blocks"], contract)
    scale = gain.joint_scale(scalers, contract)
    sample_std = gain.standardized(built["blocks"], scalers, contract)
    sample_configuration = configuration(sample_std, built["metadata"], contract)
    remainder = gain.heldout_rows(cohort, built["rows"])
    if remainder.size == 0:
        raise ValueError("The remainder is empty: the point holds every participant.")
    remainder_metadata = cohort.metadata.iloc[remainder].reset_index(drop=True)
    remainder_std = gain.standardized({name: cohort.blocks[name].iloc[remainder] for name in names}, scalers,
                                      contract)
    remainder_configuration = configuration(remainder_std, remainder_metadata, contract)
    halves: list[tuple[np.ndarray, np.ndarray]] = []
    pairs: list[tuple[np.ndarray, np.ndarray]] = []
    for seed in seeds:
        first, second = draw_split(remainder_metadata, contract, seed)
        halves.append((remainder[first], remainder[second]))
        pairs.append(tuple(configuration(remainder_std[part], remainder_metadata.iloc[part].reset_index(drop=True),
                                         contract) for part in (first, second)))  # type: ignore[misc]
    timing["configurations"] = perf_counter() - mark

    dataset = aac.make_dataset(built["blocks"], built["metadata"], contract)
    mark = perf_counter()
    production = fit_analysis(dataset, settings, analysis_seed, None)
    timing["cv_fit"] = perf_counter() - mark
    mark = perf_counter()
    forced = {int(rank): fit_analysis(dataset, settings, analysis_seed, int(rank)) for rank in settings["fixed_ranks"]}
    timing["forced_fits"] = perf_counter() - mark

    mark = perf_counter()
    model = production.latent_artifacts.model
    bases, dimensions, scalings = space_bases(production, forced, settings, n_stages)

    def inputs(z: np.ndarray) -> np.ndarray:
        # The gain driver's arithmetic: standardized → raw → the model's input coordinates.
        return gain.model_inputs(np.asarray(z, dtype=float) * scale, scale, model)

    products = mean_products([split_products(inputs(m1), inputs(m2), bases) for m1, m2 in pairs])
    naive = {"naive_remainder": naive_products(inputs(remainder_configuration), bases),
             "naive_sample": naive_products(inputs(sample_configuration), bases)}
    latent = production.record["latent"]
    realized = {"labels": built["table"], SPACE_CV: aac.recorded_table(production, contract)}
    checks: dict[str, dict[str, Any]] = {
        "cell_table": _renamed(aac.check_cell_table(point_spec(settings, built["point"])["table"], realized)),
        "parts": check_parts(len(cohort.metadata), built["rows"], halves),
        "basis": check_basis(bases, dimensions, scalings),
        "gram": check_gram(products),
    }
    out: dict[str, Any] = {}
    if planted_raw is not None:
        truth_inputs = gain.model_inputs(np.asarray(planted_raw, dtype=float), scale, model)
        out["truth"] = truth_products(truth_inputs, bases)
        out["planted_energy"] = float(np.sum(truth_inputs**2))
        checks["truth"] = check_truth(cohort, built, remainder, planted_raw, noise_free, scale)
    timing["readings"] = perf_counter() - mark

    total = products["total"]
    out.update({
        "production": {"rank": int(latent["selected_lv"]), "component_selection": latent.get("component_selection"),
                       "cv_mean_auroc": latent.get("cv_mean_auroc")},
        "halves": {"remainder": int(remainder.size), "sizes": [int(halves[0][0].size), int(halves[0][1].size)]},
        "reliability": total["correlation"],
        "cross_products": {
            "total": {"cross": total["cross"], "gram": total["gram"], "noise": total["noise"],
                      **{reading: entry["total"] for reading, entry in naive.items()}},
            "spaces": {
                name: {"dimension": int(np.asarray(bases[name]).shape[1]), **products["spaces"][name],
                       **{reading: entry["spaces"][name] for reading, entry in naive.items()}}
                for name in bases
            },
        },
        "checks": checks,
    })
    return out


def check_truth(cohort: aac.Cohort, built: Mapping[str, Any], remainder: np.ndarray, planted_raw: np.ndarray,
                noise_free: np.ndarray | None, scale: np.ndarray) -> dict[str, Any]:
    """The planted configuration recovered from the noise-free features of the sample and of the remainder."""

    if noise_free is None:
        return {**psc._check(None, TOLERANCE), "failing": {"noise_free": "absent"}}
    contract = cohort.contract
    target = np.asarray(planted_raw, dtype=float) / scale
    deviations = {}
    for name, rows in (("sample", np.asarray(built["rows"], dtype=int)), ("remainder", remainder)):
        metadata = cohort.metadata.iloc[rows].reset_index(drop=True)
        recovered = configuration(np.asarray(noise_free, dtype=float)[rows] / scale, metadata, contract)
        deviations[name] = float(np.max(np.abs(recovered - target)))
    return _named_check(deviations, TOLERANCE)


def _finish(base: dict[str, Any], parts: dict[str, Any], names: Sequence[str], start: float,
            timing: dict[str, float]) -> dict[str, Any]:
    checks = parts.pop("checks")
    failing = [name for name in names if not checks[name]["passed"]]
    timing["unit"] = perf_counter() - start
    return lsc._jsonable(_listed({
        **base,
        **parts,
        "status": "inconsistent" if failing else "completed",
        "inconsistent_checks": failing,
        "consistency": checks,
        "runtime_seconds": timing,
    }))


def _listed(value: Any) -> Any:
    if isinstance(value, Mapping):
        return {key: _listed(item) for key, item in value.items()}
    if isinstance(value, np.ndarray):
        return value.tolist()
    if isinstance(value, list | tuple):
        return [_listed(item) for item in value]
    return value


def _failed_unit(base: Mapping[str, Any], exc: Exception, start: float, timing: dict[str, float]) -> dict[str, Any]:
    return lsc._jsonable(psc._failed_unit(base, exc, start, timing))


def evaluate_cohort_unit(unit: Unit, settings: Mapping[str, Any], cohort: aac.Cohort) -> dict[str, Any]:
    """The anchor driver's unit at the point and replicate, with its splits; one record of seeds and statistics."""

    start = perf_counter()
    point = point_spec(settings, unit.point)
    seeds = aac.unit_seeds(settings, unit.point, unit.replicate)
    splits = split_seeds(settings, unit)
    base: dict[str, Any] = {
        "profile": COHORT, "point": unit.point, "replicate": int(unit.replicate),
        "participants": int(point["participants"]),
        "seeds": {**{key: seeds[key] for key in ("subsample", "labels", "rrpp")}, "reference_split": splits},
        "settings_hash": lsc.settings_hash(settings),
    }
    timing: dict[str, float] = {}
    try:
        mark = perf_counter()
        built = aac.build_unit(cohort, settings, unit.point, unit.replicate)
        base.update({"subsample_redraws": built["subsample_redraws"], "label_redraws": built["label_redraws"],
                     "cell_table": built["table"]})
        timing["build"] = perf_counter() - mark
        parts = evaluate_parts(cohort, settings, built, seeds=splits, analysis_seed=int(seeds["rrpp"]), timing=timing)
    except Exception as exc:  # recorded, never aborts the shard
        return _failed_unit(base, exc, start, timing)
    return _finish(base, parts, CHECKS[COHORT], start, timing)


# ── The calibration (design D7) ─────────────────────────────────────────────


def raw_features(cohort: aac.Cohort) -> np.ndarray:
    return np.hstack([cohort.blocks[name].to_numpy(dtype=float) for name in cohort.contract.block_names])


def fit_generator(cohort: aac.Cohort) -> dict[str, Any]:
    """The full cohort's covariate-adjusted stage fit and residual covariance; in memory only."""

    contract = cohort.contract
    n_stages = len(contract.stage_order)
    features = raw_features(cohort)
    stages = gain.stage_positions(cohort.metadata, contract)
    covariates = np.asarray(gain.covariate_columns(cohort.metadata, contract), dtype=float)
    covariates = covariates.reshape(features.shape[0], -1)
    centre = covariates.mean(axis=0)
    indicators = (stages[:, None] == np.arange(n_stages)[None, :]).astype(float)
    design = np.hstack([indicators, covariates - centre])
    rank = int(np.linalg.matrix_rank(design))
    if rank < design.shape[1]:
        raise ValueError(f"The generator's design is rank-deficient: rank {rank} of {design.shape[1]} columns.")
    coefficients, *_ = np.linalg.lstsq(design, features, rcond=None)
    residuals = features - design @ coefficients
    covariance = residuals.T @ residuals / float(features.shape[0] - rank)
    stage_coefficients = coefficients[:n_stages]
    grand_mean = stage_coefficients.mean(axis=0)
    return {
        "grand_mean": grand_mean,
        "centred": stage_coefficients - grand_mean,
        "covariate_part": (covariates - centre) @ coefficients[n_stages:],
        "stages": stages,
        "covariance": covariance,
        "factor": np.linalg.cholesky(covariance),
        "dof": int(features.shape[0] - rank),
    }


_GENERATORS: dict[tuple[tuple[str, str], ...], dict[str, Any]] = {}


def generator_for(cohort: aac.Cohort) -> dict[str, Any]:
    key = tuple(sorted(cohort.files.items()))
    if key not in _GENERATORS:
        _GENERATORS[key] = fit_generator(cohort)
    return _GENERATORS[key]


def synthetic_features(generator: Mapping[str, Any], strength: float, seed: int, *, noise_scale: float = 1.0
                       ) -> np.ndarray:
    """One synthetic cohort's raw joint features (participants × features), in the cohort's row order."""

    stages = np.asarray(generator["stages"], dtype=int)
    signal = generator["grand_mean"] + float(strength) * generator["centred"][stages] + generator["covariate_part"]
    factor = np.asarray(generator["factor"], dtype=float)
    noise = np.random.default_rng(int(seed)).standard_normal(signal.shape) @ factor.T
    return signal + float(noise_scale) * noise


def synthetic_cohort(cohort: aac.Cohort, features: np.ndarray) -> aac.Cohort:
    """The cohort with its blocks replaced by ``features`` (same participants, stages, sexes and covariates)."""

    blocks: dict[str, pd.DataFrame] = {}
    start = 0
    for name in cohort.contract.block_names:
        frame = cohort.blocks[name]
        width = frame.shape[1]
        blocks[name] = pd.DataFrame(np.asarray(features)[:, start:start + width], index=frame.index,
                                    columns=frame.columns)
        start += width
    return dataclasses.replace(cohort, blocks=blocks, _cache={})


def build_calibration_unit(synthetic: aac.Cohort, settings: Mapping[str, Any], unit: Unit) -> dict[str, Any]:
    """The anchor driver's unit construction under this driver's seeds: a subsample within cells, labels permuted."""

    contract = synthetic.contract
    table = point_spec(settings, unit.point)["table"]
    labels, stages = synthetic.real_labels, synthetic.stages
    rows = None
    subsample_redraws = 0
    for attempt in range(int(settings["max_attempts"])):
        rng = np.random.default_rng(calibration_seed(settings, unit, "subsample", attempt))
        chosen = []
        for g, group in enumerate(contract.groups):
            for s, stage in enumerate(contract.stage_order):
                cell = np.flatnonzero((labels == group) & (stages == stage))
                chosen.append(rng.choice(cell, size=int(table[g][s]), replace=False))
        candidate = np.sort(np.concatenate(chosen))
        if aac._levels_present(synthetic.metadata.iloc[candidate], contract):
            rows, subsample_redraws = candidate, attempt
            break
    if rows is None:
        raise RuntimeError(f"{settings['max_attempts']} subsamples all missed a declared covariate level.")
    metadata = synthetic.metadata.iloc[rows].reset_index(drop=True)
    real = metadata[contract.group_col].astype(str).to_numpy()
    unit_stages = metadata[contract.stage_col].astype(str).to_numpy()
    drawn = None
    label_redraws = 0
    for attempt in range(int(settings["max_attempts"])):
        rng = np.random.default_rng(calibration_seed(settings, unit, "labels", attempt))
        candidate_labels = pc.permute_within_stage(real, unit_stages, contract.stage_order, rng)
        if not np.array_equal(candidate_labels, np.asarray(real, dtype=object)):
            drawn, label_redraws = candidate_labels, attempt
            break
    if drawn is None:
        raise RuntimeError(f"{settings['max_attempts']} within-stage permutations all reproduced the labels.")
    metadata[contract.group_col] = np.asarray(drawn).astype(str)
    return {
        "point": unit.point,
        "replicate": int(unit.replicate),
        "subsample_redraws": int(subsample_redraws),
        "label_redraws": int(label_redraws),
        "rows": rows,
        "metadata": metadata,
        "blocks": {name: frame.iloc[rows] for name, frame in synthetic.blocks.items()},
        "table": aac.cell_table(metadata[contract.group_col], unit_stages, contract),
    }


def evaluate_calibration_unit(unit: Unit, settings: Mapping[str, Any], cohort: aac.Cohort, *,
                              noise_scale: float = 1.0) -> dict[str, Any]:
    """A synthetic cohort's unit through the cohort profile's pipeline, with its true capture."""

    start = perf_counter()
    point = point_spec(settings, unit.point)
    strength = float(unit.strength or 0.0)
    noise = noise_seed(settings, strength, int(unit.cohort or 0))
    splits = split_seeds(settings, unit)
    analysis_seed = calibration_seed(settings, unit, "subsample", 0)
    base: dict[str, Any] = {
        "profile": CALIBRATION, "point": unit.point, "replicate": int(unit.replicate), "strength": strength,
        "cohort": int(unit.cohort or 0), "participants": int(point["participants"]),
        "seeds": {"cohort_noise": noise, "subsample": analysis_seed, "reference_split": splits},
        "settings_hash": lsc.settings_hash(settings),
    }
    if float(noise_scale) != 1.0:
        base["noise_scale"] = float(noise_scale)
    timing: dict[str, float] = {}
    try:
        mark = perf_counter()
        generator = generator_for(cohort)
        synthetic = synthetic_cohort(cohort, synthetic_features(generator, strength, noise, noise_scale=noise_scale))
        noise_free = synthetic_features(generator, strength, noise, noise_scale=0.0)
        built = build_calibration_unit(synthetic, settings, unit)
        base.update({"subsample_redraws": built["subsample_redraws"], "label_redraws": built["label_redraws"],
                     "cell_table": built["table"]})
        timing["build"] = perf_counter() - mark
        parts = evaluate_parts(synthetic, settings, built, seeds=splits, analysis_seed=analysis_seed, timing=timing,
                               planted_raw=strength * generator["centred"], noise_free=noise_free)
    except Exception as exc:  # recorded, never aborts the shard
        return _failed_unit(base, exc, start, timing)
    return _finish(base, parts, CHECKS[CALIBRATION], start, timing)


def evaluate_unit(unit: Unit, settings: Mapping[str, Any], cohort: aac.Cohort) -> dict[str, Any]:
    if settings["profile"] == COHORT:
        return evaluate_cohort_unit(unit, settings, cohort)
    return evaluate_calibration_unit(unit, settings, cohort)


# ── Unit readings (used by the pre-flight, merge and tests) ─────────────────


def unit_reading(record: Mapping[str, Any], space: str, reading: str, n_stages: int) -> tuple[float, float] | None:
    """A unit's (numerator, denominator) for a space and reading; ``None`` when absent."""

    products = record.get("cross_products") or {}
    total, inside = products.get("total"), (products.get("spaces") or {}).get(space)
    if total is None or inside is None:
        return None
    if reading == "total":
        return float(inside["cross"]), float(total["cross"])
    if reading in CONTRASTS:
        projector = contrast_projectors(n_stages)[reading]
        return contrast_energy(inside["gram"], projector), contrast_energy(total["gram"], projector)
    if reading in ("noise", "naive_remainder", "naive_sample"):
        return float(inside[reading]), float(total[reading])
    raise ValueError(f"Unknown reading {reading!r}.")


def truth_reading(record: Mapping[str, Any], space: str, reading: str, n_stages: int) -> tuple[float, float] | None:
    truth = record.get("truth")
    if not truth or reading not in TRUE_READINGS:
        return None
    total, inside = truth["total"], truth["spaces"][space]
    if reading == "total":
        return float(inside["cross"]), float(total["cross"])
    projector = contrast_projectors(n_stages)[reading]
    return contrast_energy(inside["gram"], projector), contrast_energy(total["gram"], projector)


def ratio(numerator: float, denominator: float) -> float | None:
    return float(numerator) / float(denominator) if denominator != 0 and math.isfinite(denominator) else None


# ── Pre-flights (design D9) ──────────────────────────────────────────────────


def profile_dir(out_dir: Path, profile: str) -> Path:
    if profile not in PROFILES:
        raise SystemExit(f"--profile must be one of {list(PROFILES)}, got {profile!r}.")
    return Path(out_dir) / profile


def verify(settings: Mapping[str, Any], cohort: aac.Cohort | None = None) -> None:
    """Refuse a changed parent file, or an export whose files differ from the settings', by name."""

    plsc.verify_parents(settings)
    if cohort is not None:
        plsc.verify_export(cohort, settings)


def _prepare(directory: Path, settings: Mapping[str, Any], cohort: aac.Cohort) -> None:
    if (directory / "settings.json").exists():
        verify(lsc.load_settings(directory))
    verify(settings, cohort)
    lsc.ensure_settings(directory, settings)


def load_preflight(directory: Path, settings: Mapping[str, Any]) -> dict[str, Any] | None:
    path = Path(directory) / PREFLIGHT_NAME
    if not path.exists():
        return None
    record = json.loads(path.read_text(encoding="utf-8"))
    return record if record.get("settings_hash") == lsc.settings_hash(settings) else None


def cohort_reliability(cohort: aac.Cohort, settings: Mapping[str, Any]) -> dict[str, Any]:
    """The full cohort configuration's split-half agreement (its own scaler), over the settings' splits."""

    contract = cohort.contract
    features = gain.standardized(cohort.blocks, gain.sample_scalers(cohort.blocks, contract), contract)
    correlations, grams = [], []
    for index in range(int(settings["splits"])):
        seed = own_seed(settings, aac.COHORT_POINT, 0, "reference_split", index)
        halves = draw_split(cohort.metadata, contract, seed)
        first, second = (configuration(features[part], cohort.metadata.iloc[part].reset_index(drop=True), contract)
                         for part in halves)
        correlations.append(float(np.sum(first * second) / (np.linalg.norm(first) * np.linalg.norm(second))))
        grams.append(_sym(first, second))
    correlation = float(np.mean(correlations))
    full = 2.0 * correlation / (1.0 + correlation) if correlation > -1.0 else None
    eigenvalues = np.sort(np.linalg.eigvalsh(np.mean(np.stack(grams), axis=0)))[::-1]
    return {
        "splits": int(settings["splits"]),
        "halves": [int(part.size) for part in draw_split(cohort.metadata, contract, own_seed(
            settings, aac.COHORT_POINT, 0, "reference_split", 0))],
        "correlation_mean": correlation,
        "correlation_range": [float(min(correlations)), float(max(correlations))],
        "spearman_brown_full_cohort": full,
        "implied_strength": math.sqrt(full) if full is not None and full > 0 else None,
        "settings_interval_strength": INTERVAL_STRENGTH,
        "cross_gram_eigenvalues": [float(value) for value in eigenvalues],
        "note": "recorded, not judged: the calibration's strength 0.67 extrapolates item 22's split-half 0.166 at "
                "375 participants per half (Spearman-Brown 0.443 at the full cohort, square root 0.666)",
    }


def cohort_preflight(out_dir: Path, settings: Mapping[str, Any], cohort: aac.Cohort) -> dict[str, Any]:
    """The anchor driver's pre-flight units under this driver's settings, and the cohort's reliability."""

    directory = profile_dir(out_dir, COHORT)
    _prepare(directory, settings, cohort)
    units = []
    for replicate in settings["preflight_replicates"]:
        try:
            entry = _filtered(aac.preflight_unit, int(replicate), settings, cohort)
            # The fit's equality is kept by name; the statistics (of permuted labels) are not recorded.
            units.append({key: value for key, value in entry.items() if key not in ("cross_validated", "forced")})
        except Exception as exc:  # a pre-flight that raises fails
            units.append({"unit": aac.unit_key(aac.COHORT_POINT, int(replicate)), "passed": False,
                          "error": f"{type(exc).__name__}: {exc}"})
    try:
        reliability: dict[str, Any] = cohort_reliability(cohort, settings)
    except Exception as exc:  # recorded, never decisive
        reliability = {"error": f"{type(exc).__name__}: {exc}"}
    result = {
        "profile": COHORT,
        "settings_hash": lsc.settings_hash(settings),
        "rule": "the cross-validated pls_cv on the cohort anchor's units selects the contract's rank and AUROC "
                "(rounded), and its statistics equal the call forced at that rank exactly (no permutation drawn)",
        "expected": dict(settings["contract"]["fidelity"]),
        "units": units,
        "reliability": reliability,
        "passed": bool(units) and all(unit["passed"] for unit in units),
    }
    lsc._write_json(directory / PREFLIGHT_NAME, result)
    return result


def calibration_preflight(out_dir: Path, settings: Mapping[str, Any], cohort: aac.Cohort) -> dict[str, Any]:
    """At strength 1 with the noise scaled to zero, each point's unit estimates its true capture exactly."""

    directory = profile_dir(out_dir, CALIBRATION)
    _prepare(directory, settings, cohort)
    spec = settings["preflight"]
    n_stages = len(cohort.contract.stage_order)
    units = []
    for point in point_names(settings):
        unit = Unit(point, int(spec["unit"]), float(spec["strength"]), int(spec["cohort"]))
        record = evaluate_calibration_unit(unit, settings, cohort, noise_scale=float(spec["noise_scale"]))
        entry: dict[str, Any] = {"unit": unit_key(unit), "status": record["status"],
                                 "inconsistent_checks": record.get("inconsistent_checks", []),
                                 "seconds": (record.get("runtime_seconds") or {}).get("unit")}
        if record["status"] == "failed":
            entry.update({"error": f"{record.get('error_type')}: {record.get('error_message')}", "passed": False})
            units.append(entry)
            continue
        deviations: dict[str, float] = {}
        for space in settings["spaces"]:
            for reading in TRUE_READINGS:
                estimated, truth = unit_reading(record, space, reading, n_stages), truth_reading(
                    record, space, reading, n_stages)
                if estimated is None or truth is None:
                    deviations[f"{space}:{reading}"] = math.inf
                    continue
                left, right = ratio(*estimated), ratio(*truth)
                if left is None and right is None and abs(estimated[1]) <= TOLERANCE and abs(truth[1]) <= TOLERANCE:
                    deviations[f"{space}:{reading}"] = 0.0  # an empty contrast (two stages): nothing to estimate
                    continue
                deviations[f"{space}:{reading}"] = (abs(left - right) if left is not None and right is not None
                                                    else math.inf)
        worst = max(deviations.values(), default=math.inf)
        entry.update({"rank": (record.get("production") or {}).get("rank"), "worst_deviation": worst,
                      "failing": {key: value for key, value in deviations.items() if not value <= TOLERANCE},
                      "passed": record["status"] == "completed" and worst <= TOLERANCE})
        units.append(entry)
    result = {
        "profile": CALIBRATION,
        "settings_hash": lsc.settings_hash(settings),
        "rule": spec["rule"],
        "tolerance": TOLERANCE,
        "units": units,
        "passed": bool(units) and all(unit["passed"] for unit in units),
    }
    lsc._write_json(directory / PREFLIGHT_NAME, result)
    return result


def preflight(out_dir: Path, settings: Mapping[str, Any], cohort: aac.Cohort) -> dict[str, Any]:
    if settings["profile"] == COHORT:
        return cohort_preflight(out_dir, settings, cohort)
    return calibration_preflight(out_dir, settings, cohort)


# ── Shards, resume and the run ───────────────────────────────────────────────


def read_records(paths: Iterable[Path]) -> list[dict[str, Any]]:
    """Every parsable unit line; a line torn by an interrupted write is dropped, so its unit is re-run."""

    return psc.read_records(paths)


_WORKER_COHORT: aac.Cohort | None = None


def _init_worker(directory: str, contract: aac.CohortContract, files: Mapping[str, str]) -> None:
    global _WORKER_COHORT
    cohort = aac.load_cohort(Path(directory), contract)
    if cohort.files != dict(files):
        raise SystemExit(f"{directory} changed after the run started (file hashes differ).")
    _WORKER_COHORT = cohort


def _evaluate_pending(unit: Unit, settings: Mapping[str, Any]) -> dict[str, Any]:
    assert _WORKER_COHORT is not None
    return evaluate_unit(unit, settings, _WORKER_COHORT)


def run_shard(
    out_dir: Path,
    settings: Mapping[str, Any],
    cohort: aac.Cohort,
    index: int,
    count: int,
    *,
    limit: int | None = None,
    workers: int = 1,
) -> int:
    """Evaluate this shard's missing units of the settings' profile; returns how many were evaluated."""

    profile = str(settings["profile"])
    directory = profile_dir(out_dir, profile)
    _prepare(directory, settings, cohort)
    checked = load_preflight(directory, settings)
    if checked is None:
        raise SystemExit(f"{directory / PREFLIGHT_NAME} is missing or was written under other settings; run "
                         f"`preflight --profile {profile}` first.")
    if not checked.get("passed"):
        raise SystemExit(f"The pre-flight failed ({directory / PREFLIGHT_NAME}); the run stops.")
    path = lsc.shard_path(directory, index, count)
    path.parent.mkdir(parents=True, exist_ok=True)
    done = {record_key(record): record for record in (read_records([path]) if path.exists() else [])}
    # Rewritten with complete units only, so a torn line is re-run, not duplicated.
    with path.open("w", encoding="utf-8") as handle:
        handle.write("".join(lsc._line(done[key]) for key in sorted(done)))
    pending = [unit for position, unit in enumerate(enumerate_units(settings))
               if position % count == index and unit_key(unit) not in done]
    if limit is not None:
        pending = pending[:limit]

    def write(position: int, unit: Unit, record: Mapping[str, Any]) -> None:
        with path.open("a", encoding="utf-8") as handle:
            handle.write(lsc._line(record))
            handle.flush()
        print(f"[{profile} {index}/{count}] {position}/{len(pending)} {unit_key(unit)} {record['status']}",
              file=sys.stderr, flush=True)

    if workers <= 1:
        for position, unit in enumerate(pending, start=1):
            write(position, unit, evaluate_unit(unit, settings, cohort))
        return len(pending)
    with ProcessPoolExecutor(max_workers=workers, initializer=_init_worker,
                             initargs=(str(cohort.directory), cohort.contract, dict(cohort.files))) as pool:
        futures = {pool.submit(_evaluate_pending, unit, settings): unit for unit in pending}
        for position, future in enumerate(as_completed(futures), start=1):
            write(position, futures[future], future.result())
    return len(pending)


# ── Merge: capture, point and calibration tables, accounting, reproduction ──


CAPTURE_COLUMNS: tuple[str, ...] = (
    "profile", "strength", "point", "participants", "space", "reading", "estimate", "numerator_sum",
    "denominator_sum", "denominator_mean", "units", "unit_q05", "unit_q25", "unit_median", "unit_q75", "unit_q95",
    "truth",
)
POINT_COLUMNS: tuple[str, ...] = ("profile", "strength", "point", "quantity", "space", "key", "value")
CALIBRATION_COLUMNS: tuple[str, ...] = (
    "strength", "point", "space", "reading", "cohorts", "mean_estimate", "mean_truth", "bias", "sd_estimate",
    "mean_energy", "sd_energy",
)
QUANTILES: tuple[tuple[str, float], ...] = (("unit_q05", 0.05), ("unit_q25", 0.25), ("unit_median", 0.5),
                                            ("unit_q75", 0.75), ("unit_q95", 0.95))


def _median(values: Iterable[Any]) -> float | None:
    return plsc._quantile(values, 0.5)


def usable(records: Iterable[Mapping[str, Any]]) -> list[Mapping[str, Any]]:
    """Completed units only: failed and inconsistent units enter no reading."""

    return [record for record in records if record.get("status") == "completed"]


def _cells(records: Iterable[Mapping[str, Any]], settings: Mapping[str, Any]
           ) -> list[tuple[float | None, str, list[Mapping[str, Any]]]]:
    """(strength, point, records) in settings order; the cohort profile has one strength, ``None``."""

    records = list(records)
    strengths: list[float | None] = ([None] if settings["profile"] == COHORT
                                     else [float(value) for value in settings["strengths"]])
    out = []
    for strength in strengths:
        for point in point_names(settings):
            members = [record for record in records if str(record["point"]) == point and (
                strength is None or math.isclose(float(record.get("strength", math.nan)), strength))]
            out.append((strength, point, members))
    return out


def pooled(records: Sequence[Mapping[str, Any]], space: str, reading: str, n_stages: int, *, truth: bool = False
           ) -> dict[str, Any]:
    """The ratio of sums over usable units, with sums, the mean denominator and per-unit ratios."""

    reader = truth_reading if truth else unit_reading
    pairs = [pair for pair in (reader(record, space, reading, n_stages) for record in usable(records)) if pair]
    numerator = float(sum(pair[0] for pair in pairs))
    denominator = float(sum(pair[1] for pair in pairs))
    return {
        "units": len(pairs),
        "numerator_sum": numerator if pairs else None,
        "denominator_sum": denominator if pairs else None,
        "denominator_mean": denominator / len(pairs) if pairs else None,
        "estimate": ratio(numerator, denominator) if pairs else None,
        "unit_ratios": [ratio(*pair) for pair in pairs],
    }


def capture_rows(records: Sequence[Mapping[str, Any]], settings: Mapping[str, Any]) -> list[dict[str, Any]]:
    n_stages = len(settings["contract"]["stage_order"])
    rows = []
    for strength, point, members in _cells(records, settings):
        participants = point_spec(settings, point)["participants"]
        for space in settings["spaces"]:
            for reading in READINGS:
                entry = pooled(members, space, reading, n_stages)
                truth = (pooled(members, space, reading, n_stages, truth=True)["estimate"]
                         if settings["profile"] == CALIBRATION else None)
                rows.append({
                    "profile": settings["profile"], "strength": strength, "point": point,
                    "participants": participants, "space": space, "reading": reading,
                    **{key: entry[key] for key in ("estimate", "numerator_sum", "denominator_sum", "denominator_mean",
                                                   "units")},
                    **{name: plsc._quantile(entry["unit_ratios"], q) for name, q in QUANTILES},
                    "truth": truth,
                })
    return rows


def point_rows(records: Sequence[Mapping[str, Any]], settings: Mapping[str, Any]) -> list[dict[str, Any]]:
    """The point table in long form: one quantity per row; blank keys are not applicable."""

    n_stages = len(settings["contract"]["stage_order"])
    rows: list[dict[str, Any]] = []
    profile = settings["profile"]

    for strength, point, members in _cells(records, settings):
        def add(quantity: str, value: Any, **keys: Any) -> None:
            rows.append({"profile": profile, "strength": strength, "point": point, "quantity": quantity, **keys,
                         "value": value})

        units = usable(members)
        add("participants", point_spec(settings, point)["participants"])
        add("units_completed", len(units))
        add("units_failed", sum(1 for record in members if record.get("status") == "failed"))
        add("units_inconsistent", sum(1 for record in members if record.get("status") == "inconsistent"))
        ranks = [int(record["production"]["rank"]) for record in units]
        for rank in sorted(set(ranks)):
            add("cv_rank_units", ranks.count(rank), key=str(rank))
        add("cv_rank_median", _median(ranks))
        add("cv_auroc_median", _median(record["production"].get("cv_mean_auroc") for record in units))
        if units:
            for position, size in zip(("first", "second"), units[0]["halves"]["sizes"]):
                add("half_size", size, key=position)
        add("reliability_median", _median(record.get("reliability") for record in units))
        add("reliability_mean", float(np.mean([record["reliability"] for record in units
                                               if record.get("reliability") is not None]))
            if any(record.get("reliability") is not None for record in units) else None)
        for reading in READINGS:
            add("energy_mean", pooled(units, SPACE_CV, reading, n_stages)["denominator_mean"], key=reading)
        if units:
            gram = np.mean([np.asarray(record["cross_products"]["total"]["gram"], dtype=float) for record in units],
                           axis=0)
            for index, value in enumerate(np.sort(np.linalg.eigvalsh(gram))[::-1], start=1):
                add("cross_gram_eigenvalue", float(value), key=str(index))
        for space in settings["spaces"]:
            add("dimension_median", _median(record["cross_products"]["spaces"][space]["dimension"]
                                            for record in units), space=space)
            for reading in ("naive_remainder", "naive_sample"):
                shares = [ratio(*pair) for pair in (unit_reading(record, space, reading, n_stages)
                                                    for record in units) if pair]
                add("naive_median", _median(shares), space=space, key=reading)
        add("subsample_redraws", sum(int(record.get("subsample_redraws") or 0) for record in members))
        add("label_redraws", sum(int(record.get("label_redraws") or 0) for record in members))
        add("unit_seconds_median", _median((record.get("runtime_seconds") or {}).get("unit") for record in members))
    return rows


def calibration_rows(records: Sequence[Mapping[str, Any]], settings: Mapping[str, Any]) -> list[dict[str, Any]]:
    """Per strength, point, space and reading: one pooled estimate (and truth, and energy) per synthetic cohort."""

    n_stages = len(settings["contract"]["stage_order"])
    rows = []
    for strength, point, members in _cells(records, settings):
        by_cohort: dict[int, list[Mapping[str, Any]]] = {}
        for record in members:
            by_cohort.setdefault(int(record["cohort"]), []).append(record)
        for space in settings["spaces"]:
            for reading in READINGS:
                estimates, truths, energies, biases = [], [], [], []
                for cohort in sorted(by_cohort):
                    estimate = pooled(by_cohort[cohort], space, reading, n_stages)
                    truth = pooled(by_cohort[cohort], space, reading, n_stages, truth=True)["estimate"]
                    if estimate["estimate"] is not None:
                        estimates.append(estimate["estimate"])
                        if truth is not None:
                            biases.append(estimate["estimate"] - truth)
                    if truth is not None:
                        truths.append(truth)
                    if estimate["denominator_mean"] is not None:
                        energies.append(estimate["denominator_mean"])
                rows.append({
                    "strength": strength, "point": point, "space": space, "reading": reading,
                    "cohorts": len(estimates),
                    "mean_estimate": float(np.mean(estimates)) if estimates else None,
                    "mean_truth": float(np.mean(truths)) if truths else None,
                    "bias": float(np.mean(biases)) if biases else None,
                    "sd_estimate": float(np.std(estimates, ddof=1)) if len(estimates) > 1 else None,
                    "mean_energy": float(np.mean(energies)) if energies else None,
                    "sd_energy": float(np.std(energies, ddof=1)) if len(energies) > 1 else None,
                })
    return rows


def accounting(records: Sequence[Mapping[str, Any]], settings: Mapping[str, Any], *, duplicates: int) -> dict:
    """Units, failures and inconsistent units by point (and strength), with errors, checks and core-hours."""

    expected = expected_units(settings)
    by_point: dict[str, dict[str, int]] = {point: {"expected": count, "complete": 0, "failed": 0, "inconsistent": 0}
                                           for point, count in expected.items()}
    by_cell: dict[str, dict[str, int]] = {}
    failed, inconsistent = [], []
    seconds: list[float] = []
    for record in records:
        point = str(record["point"])
        cell = point if record.get("strength") is None else f"s{strength_label(record['strength'])}/{point}"
        counts = by_point[point]
        cell_counts = by_cell.setdefault(cell, {"complete": 0, "failed": 0, "inconsistent": 0})
        counts["complete"] += 1
        cell_counts["complete"] += 1
        status = record.get("status")
        if status in ("failed", "inconsistent"):
            counts[status] += 1
            cell_counts[status] += 1
        if status == "failed":
            failed.append({"unit": record_key(record), "error_type": record.get("error_type"),
                           "error_message": record.get("error_message")})
        elif status == "inconsistent":
            inconsistent.append({"unit": record_key(record), "checks": {
                name: {"deviation": record["consistency"][name].get("deviation"),
                       "failing": record["consistency"][name].get("failing")}
                for name in record["inconsistent_checks"]}})
        unit_seconds = (record.get("runtime_seconds") or {}).get("unit")
        if unit_seconds is not None:
            seconds.append(float(unit_seconds))
    for counts in by_point.values():
        counts["usable"] = counts["complete"] - counts["failed"] - counts["inconsistent"]
    total_expected = sum(expected.values())
    complete = sum(counts["complete"] for counts in by_point.values())
    return {
        "settings_hash": lsc.settings_hash(settings),
        "profile": settings["profile"],
        "expected_units": total_expected,
        "complete_units": complete,
        "missing_units": total_expected - complete,
        "failed_units": len(failed),
        "inconsistent_units": len(inconsistent),
        "by_point": by_point,
        "by_cell": by_cell,
        "failed_unit_errors": failed,
        "inconsistent_unit_checks": inconsistent,
        "duplicate_records": int(duplicates),
        "recorded_core_hours": sum(seconds) / 3600.0,
        "unit_seconds_median": _median(seconds),
    }


def _value(text: Any) -> str:
    return lsc._csv_value(text) if isinstance(text, float) else ("" if text is None else str(text))


def reproduction(points_rows: Sequence[Mapping[str, Any]], anchor_rows: Sequence[Mapping[str, Any]],
                 gain_rows: Sequence[Mapping[str, Any]]) -> dict[str, Any]:
    """The CV-rank counts at 80, 200 and 500 against the anchor run, the naive medians at 750 against the gain run.

    Values are compared as the tables write them (six significant digits).
    """

    compared: list[dict[str, Any]] = []
    for point in REPRODUCED_RANK_POINTS:
        ours = {str(row["key"]): int(float(row["value"])) for row in points_rows
                if str(row["point"]) == point and str(row["quantity"]) == "cv_rank_units"}
        theirs = {str(row["key"]): int(float(row["value"])) for row in anchor_rows
                  if str(row["point"]) == point and str(row["quantity"]) == "cv_rank_units"}
        compared.append({"point": point, "quantity": "cv_rank_units", "committed": theirs, "observed": ours,
                         "equal": bool(theirs) and theirs == ours})
    for (reading, space), key in NAIVE_PARENT_KEYS.items():
        observed = next((row["value"] for row in points_rows if str(row["point"]) == REPRODUCED_NAIVE_POINT
                         and str(row["quantity"]) == "naive_median" and str(row["space"]) == space
                         and str(row["key"]) == reading), None)
        committed = next((row["value"] for row in gain_rows if str(row["point"]) == REPRODUCED_NAIVE_POINT
                          and str(row["cell"]) == "anchor" and str(row["quantity"]) == "capture_median"
                          and str(row["key"]) == key), None)
        observed_text, committed_text = _value(observed), _value(committed)
        compared.append({"point": REPRODUCED_NAIVE_POINT, "quantity": "naive_median", "space": space,
                         "reading": reading, "parent_key": key, "committed": committed_text,
                         "observed": observed_text,
                         "equal": committed_text != "" and committed_text == observed_text})
    differing = [entry for entry in compared if not entry["equal"]]
    return {"rule": "the cohort profile's CV-rank counts at n80, n200 and n500 equal the anchor run's committed "
                    "point table, and its medians of the naive captures at n750 (remainder and sample, fitted span "
                    "and predictive subspace) equal the gain run's committed cell table, as written",
            "rows": compared, "differing": differing, "passed": not differing}


def merge(out_dir: Path, profile: str, export_dir: Path) -> dict[str, Any]:
    """Merge a profile's shards and write its tables, the accounting and (cohort) the reproduction."""

    directory = profile_dir(out_dir, profile)
    settings = lsc.load_settings(directory)
    verify(settings)
    aac.verify_export_files(export_dir, settings)
    raw = read_records(sorted((directory / "records").glob("shard_*.jsonl")))
    stale = sorted({record_key(r) for r in raw if r.get("settings_hash") != lsc.settings_hash(settings)})
    if stale:
        raise SystemExit(f"{len(stale)} unit(s) were written under other settings, e.g. {stale[:3]}; refusing.")
    seen: dict[str, int] = {}
    by_key: dict[str, dict[str, Any]] = {}
    for record in raw:
        key = record_key(record)
        seen[key] = seen.get(key, 0) + 1
        by_key[key] = record
    records = [by_key[unit_key(unit)] for unit in enumerate_units(settings) if unit_key(unit) in by_key]
    with (directory / "merged.jsonl").open("w", encoding="utf-8") as handle:
        handle.write("".join(lsc._line(record) for record in records))
    lsc._write_csv(directory / "capture_table.csv", capture_rows(records, settings), CAPTURE_COLUMNS)
    lsc._write_csv(directory / "point_table.csv", point_rows(records, settings), POINT_COLUMNS)
    if profile == CALIBRATION:
        lsc._write_csv(directory / "calibration_table.csv", calibration_rows(records, settings), CALIBRATION_COLUMNS)
    summary = accounting(records, settings, duplicates=sum(n - 1 for n in seen.values()))
    lsc._write_json(directory / "accounting.json", summary)
    out = {key: summary[key] for key in ("expected_units", "complete_units", "missing_units", "failed_units",
                                         "inconsistent_units", "duplicate_records", "recorded_core_hours")}
    if profile == COHORT:
        reproduced = reproduction(lsc._read_csv(directory / "point_table.csv"),
                                  lsc._read_csv(plsc._resolve(settings["parents"]["anchor_points"]["path"])),
                                  lsc._read_csv(plsc._resolve(settings["parents"]["gain_cells"]["path"])))
        reproduced["settings_hash"] = lsc.settings_hash(settings)
        lsc._write_json(directory / REPRODUCTION_NAME, reproduced)
        out["reproduction_passed"] = reproduced["passed"]
    return out


# ── Decision (pure functions; design D8) ─────────────────────────────────────


VERDICTS: tuple[str, ...] = ("gate_failed", "no_program", "unresolved", "partly_measured", "measured")
CONSEQUENCES: dict[str, str] = {
    "gate_failed": "Nothing is read. A diagnosis of the estimator, of the failing units or of the reproduction comes "
                   "before any reading.",
    "no_program": "ADNI's blocks carry no stage program the halves agree on, and stage-span plants on ADNI are edits "
                  "of noise. The roadmap records it and item 22's held-out readings are annotated accordingly.",
    "unresolved": "The driver's estimator stands, and the number waits for ROSMAP or a larger cohort. No capture is "
                  "stated.",
    "partly_measured": "The resolved captures are stated on the measurement-visibility page as the real-noise capture "
                       "at ADNI's design points, beside the simulated 0.10; the unresolved points are not. A "
                       "follow-up with plants along the cross-validated program is proposed for the contrasts graded "
                       "present; a rotation or a bend needs both.",
    "measured": "The captures are stated on the measurement-visibility page as the real-noise capture at ADNI's design "
                "points, beside the simulated 0.10. A follow-up with plants along the cross-validated program is "
                "proposed for the contrasts graded present; a rotation or a bend needs both.",
}
RULE: dict[str, Any] = {
    "version": RULE_VERSION,
    "gate": f"failed plus inconsistent units at most {FAILURE_SHARE:g} of each point's expected units in each profile "
            "(the calibration's points pool their strengths); both reproductions pass; at strengths "
            f"{', '.join(f'{s:g}' for s in GATE_STRENGTHS)} the absolute bias of pls_cv's total capture is at most "
            f"{BIAS_LIMIT:g} at every point; both profiles' tables present",
    "program": f"the cohort's total energy at {PROGRAM_POINT} exceeds {PRESENCE_MULTIPLIER:g} × its standard deviation "
               f"across synthetic cohorts at strength {NULL_STRENGTH:g}",
    "interval": f"estimate ± {INTERVAL_MULTIPLIER:g} × the standard deviation of the pooled estimate across synthetic "
                f"cohorts at strength {INTERVAL_STRENGTH:g}, at the same point, space and reading",
    "resolved": f"a point whose pls_cv total-capture half-width is at most {HALF_WIDTH_LIMIT:g}",
    "verdict": "gate_failed; else no_program; else unresolved (no point resolved), partly_measured (one to three), "
               "measured (all four)",
    "contrast": f"present when its energy exceeds {PRESENCE_MULTIPLIER:g} × its standard deviation at strength "
                f"{NULL_STRENGTH:g}; then captured, missed or undetermined as pls_cv's interval lies above "
                f"{CAPTURE_THRESHOLD:g}, below it or across it",
    "licenses": "under measured or partly_measured, a contrast licenses a follow-up with plants when graded present "
                "at a resolved point; a stretch needs one present contrast, a rotation or a bend needs both",
    "indicative": f"fewer than {REPLICATES} units at a cohort point, or fewer than {CALIBRATION_COHORTS} cohorts or "
                  f"{CALIBRATION_UNITS} units per cohort and point in the calibration",
}


def _number(value: Any) -> float | None:
    return plsc._number(value)


def _matches(row: Mapping[str, Any], **keys: Any) -> bool:
    for key, wanted in keys.items():
        value = row.get(key)
        if key == "strength":
            if value in (None, "") or not math.isclose(float(value), float(wanted)):
                return False
        elif str(value) != str(wanted):
            return False
    return True


def _row(rows: Sequence[Mapping[str, Any]], **keys: Any) -> Mapping[str, Any] | None:
    return next((row for row in rows if _matches(row, **keys)), None)


def gate(
    *, accounting_by_profile: Mapping[str, Mapping[str, Any] | None], reproduced: Mapping[str, Any] | None,
    calibration: Sequence[Mapping[str, Any]], points: Sequence[str],
) -> dict[str, Any]:
    """The three causes (units, reproduction, bias) and the inputs, each with what failed."""

    causes: list[str] = []
    inputs = {profile: accounting_by_profile.get(profile) is not None for profile in PROFILES}
    if not all(inputs.values()) or reproduced is None:
        causes.append("inputs")
    units: dict[str, Any] = {}
    for profile in PROFILES:
        summary = accounting_by_profile.get(profile)
        if summary is None:
            continue
        for point, counts in summary["by_point"].items():
            bad = int(counts.get("failed", 0)) + int(counts.get("inconsistent", 0))
            allowed = FAILURE_SHARE * int(counts.get("expected", 0))
            units[f"{profile}/{point}"] = {"failed_or_inconsistent": bad, "allowed": allowed, "passed": bad <= allowed}
    if any(not entry["passed"] for entry in units.values()):
        causes.append("units")
    reproduction_passed = bool(reproduced and reproduced.get("passed"))
    if reproduced is not None and not reproduction_passed:
        causes.append("reproduction")
    bias: dict[str, Any] = {}
    for strength in GATE_STRENGTHS:
        for point in points:
            row = _row(calibration, strength=strength, point=point, space=DESIGNATED, reading="total")
            value = None if row is None else _number(row.get("bias"))
            bias[f"s{strength_label(strength)}/{point}"] = {
                "bias": value, "passed": value is not None and abs(value) <= BIAS_LIMIT}
    if any(not entry["passed"] for entry in bias.values()):
        causes.append("bias")
    return {"passed": not causes, "causes": causes, "inputs": inputs, "units": units,
            "reproduction": None if reproduced is None else {"passed": reproduction_passed,
                                                             "differing": reproduced.get("differing", [])},
            "bias": bias}


def _null_sd(calibration: Sequence[Mapping[str, Any]], point: str, reading: str) -> float | None:
    row = _row(calibration, strength=NULL_STRENGTH, point=point, space=DESIGNATED, reading=reading)
    return None if row is None else _number(row.get("sd_energy"))


def _energy(captures: Sequence[Mapping[str, Any]], point: str, reading: str) -> float | None:
    row = _row(captures, point=point, space=DESIGNATED, reading=reading)
    return None if row is None else _number(row.get("denominator_mean"))


def program_present(captures: Sequence[Mapping[str, Any]], calibration: Sequence[Mapping[str, Any]]) -> dict:
    energy, sd = _energy(captures, PROGRAM_POINT, "total"), _null_sd(calibration, PROGRAM_POINT, "total")
    threshold = None if sd is None else PRESENCE_MULTIPLIER * sd
    return {"point": PROGRAM_POINT, "energy": energy, "null_sd": sd, "threshold": threshold,
            "present": energy is not None and threshold is not None and energy > threshold}


def intervals(captures: Sequence[Mapping[str, Any]], calibration: Sequence[Mapping[str, Any]]) -> list[dict]:
    rows = []
    for row in captures:
        estimate = _number(row.get("estimate"))
        reference = _row(calibration, strength=INTERVAL_STRENGTH, point=row["point"], space=row["space"],
                         reading=row["reading"])
        sd = None if reference is None else _number(reference.get("sd_estimate"))
        half = None if sd is None else INTERVAL_MULTIPLIER * sd
        rows.append({"point": str(row["point"]), "space": str(row["space"]), "reading": str(row["reading"]),
                     "estimate": estimate, "sd": sd, "half_width": half,
                     "lower": None if estimate is None or half is None else estimate - half,
                     "upper": None if estimate is None or half is None else estimate + half})
    return rows


def _interval(table: Sequence[Mapping[str, Any]], point: str, reading: str, space: str = DESIGNATED
              ) -> Mapping[str, Any]:
    return next((row for row in table if row["point"] == point and row["space"] == space
                 and row["reading"] == reading), {"estimate": None, "half_width": None, "lower": None, "upper": None})


def resolution(table: Sequence[Mapping[str, Any]], points: Sequence[str]) -> dict[str, Any]:
    out = {}
    for point in points:
        half = _interval(table, point, "total")["half_width"]
        out[point] = {"half_width": half, "resolved": half is not None and half <= HALF_WIDTH_LIMIT}
    return out


def contrast_grades(captures: Sequence[Mapping[str, Any]], calibration: Sequence[Mapping[str, Any]],
                    table: Sequence[Mapping[str, Any]], points: Sequence[str]) -> dict[str, Any]:
    out: dict[str, Any] = {}
    for point in points:
        out[point] = {}
        for contrast in CONTRASTS:
            energy, sd = _energy(captures, point, contrast), _null_sd(calibration, point, contrast)
            present = energy is not None and sd is not None and energy > PRESENCE_MULTIPLIER * sd
            interval = _interval(table, point, contrast)
            grade: str | None = None
            if present:
                lower, upper = interval["lower"], interval["upper"]
                if lower is None or upper is None:
                    grade = "undetermined"
                elif lower > CAPTURE_THRESHOLD:
                    grade = "captured"
                elif upper < CAPTURE_THRESHOLD:
                    grade = "missed"
                else:
                    grade = "undetermined"
            out[point][contrast] = {"energy": energy, "null_sd": sd, "present": present,
                                    "estimate": interval["estimate"], "lower": interval["lower"],
                                    "upper": interval["upper"], "grade": grade}
    return out


def verdict_of(gated: Mapping[str, Any], program: Mapping[str, Any], resolved: Mapping[str, Any]) -> str:
    if not gated["passed"]:
        return "gate_failed"
    if not program["present"]:
        return "no_program"
    count = sum(1 for entry in resolved.values() if entry["resolved"])
    if count == 0:
        return "unresolved"
    return "measured" if count == len(resolved) else "partly_measured"


def licenses(verdict: str, table: Sequence[Mapping[str, Any]], resolved: Mapping[str, Any],
             grades: Mapping[str, Any]) -> dict[str, Any]:
    stating = verdict in ("measured", "partly_measured")
    resolved_points = [point for point, entry in resolved.items() if entry["resolved"]]
    stated = [{key: _interval(table, point, "total")[key] for key in ("estimate", "lower", "upper")} | {
        "point": point, "space": DESIGNATED} for point in resolved_points] if stating else []
    contrasts = {}
    for contrast in CONTRASTS:
        present_at = [point for point in resolved_points if grades[point][contrast]["present"]]
        contrasts[contrast] = {"present_at_resolved_points": present_at, "licensed": stating and bool(present_at)}
    any_present = any(entry["licensed"] for entry in contrasts.values())
    both = all(entry["licensed"] for entry in contrasts.values())
    return {
        "stated_captures": stated,
        "follow_up": contrasts,
        "plants": {"magnitude": any_present, "orientation": both, "shape": both},
        "statement": CONSEQUENCES[verdict],
    }


def indicative(settings_by_profile: Mapping[str, Mapping[str, Any] | None],
               accounting_by_profile: Mapping[str, Mapping[str, Any] | None]) -> dict[str, Any]:
    reasons: list[str] = []
    cohort, calibration = settings_by_profile.get(COHORT), settings_by_profile.get(CALIBRATION)
    if cohort is None or calibration is None:
        reasons.append("a profile is missing")
    if cohort is not None:
        if int(cohort["replicates"]) < REPLICATES:
            reasons.append(f"cohort: {cohort['replicates']} replicates per point, below {REPLICATES}")
        summary = accounting_by_profile.get(COHORT)
        for point, counts in ((summary or {}).get("by_point") or {}).items():
            if int(counts.get("complete", 0)) < REPLICATES:
                reasons.append(f"cohort: {counts.get('complete', 0)} complete units at {point}, below {REPLICATES}")
    if calibration is not None:
        if int(calibration["cohorts"]) < CALIBRATION_COHORTS:
            reasons.append(f"calibration: {calibration['cohorts']} cohorts, below {CALIBRATION_COHORTS}")
        if int(calibration["units"]) < CALIBRATION_UNITS:
            reasons.append(f"calibration: {calibration['units']} units per cohort and point, below "
                           f"{CALIBRATION_UNITS}")
        summary = accounting_by_profile.get(CALIBRATION)
        if summary is not None and int(summary.get("missing_units", 0)) > 0:
            reasons.append(f"calibration: {summary['missing_units']} units missing")
    return {"indicative": bool(reasons), "reasons": reasons}


def decision(
    settings_by_profile: Mapping[str, Mapping[str, Any] | None],
    captures: Sequence[Mapping[str, Any]],
    calibration: Sequence[Mapping[str, Any]],
    accounting_by_profile: Mapping[str, Mapping[str, Any] | None],
    reproduced: Mapping[str, Any] | None,
) -> dict[str, Any]:
    """The decision from the cohort's capture table, the calibration table, both accountings and the reproduction."""

    reference = settings_by_profile.get(COHORT) or settings_by_profile.get(CALIBRATION)
    points = point_names(reference) if reference is not None else [f"n{size}" for size in POINT_SIZES]
    gated = gate(accounting_by_profile=accounting_by_profile, reproduced=reproduced, calibration=calibration,
                 points=points)
    program = program_present(captures, calibration)
    table = intervals(captures, calibration)
    resolved = resolution(table, points)
    grades = contrast_grades(captures, calibration, table, points)
    verdict = verdict_of(gated, program, resolved)
    return {
        "rule": RULE,
        "rule_version": RULE_VERSION,
        "verdict": verdict,
        "consequence": CONSEQUENCES[verdict],
        "gate": gated,
        "program": program,
        "resolution": resolved,
        "designated": [row for row in table if row["space"] == DESIGNATED],
        "intervals": table,
        "contrast_grades": grades,
        "licenses": licenses(verdict, table, resolved, grades),
        "indicative": indicative(settings_by_profile, accounting_by_profile),
    }


def _profile_inputs(directory: Path, tables: Sequence[str]) -> tuple[Any, dict[str, Any], Any, Any, dict]:
    """A merged profile's settings, tables, accounting and pre-flight (``None`` where absent or stale)."""

    hashes: dict[str, str] = {}
    if not (directory / "settings.json").exists():
        return None, {}, None, None, hashes
    settings = lsc.load_settings(directory)
    verify(settings)
    digest = lsc.settings_hash(settings)
    summary = None
    path = directory / "accounting.json"
    if path.exists():
        stored = json.loads(path.read_text(encoding="utf-8"))
        if stored.get("settings_hash") == digest:
            summary = stored
            hashes[f"{directory.name}/accounting.json"] = plsc.sha256_file(path)
    loaded: dict[str, Any] = {}
    if summary is not None:
        for name in tables:
            table_path = directory / name
            if table_path.exists():
                loaded[name] = lsc._read_csv(table_path)
                hashes[f"{directory.name}/{name}"] = plsc.sha256_file(table_path)
    checked = load_preflight(directory, settings)
    return settings, loaded, summary, checked, hashes


def decide(out_dir: Path, export_dir: Path) -> dict[str, Any]:
    """Read both merged profiles (a missing one fails the gate) and write the decision."""

    out_dir = Path(out_dir)
    settings_by_profile: dict[str, Any] = {}
    summaries: dict[str, Any] = {}
    preflights: dict[str, Any] = {}
    tables: dict[str, Any] = {}
    inputs: dict[str, str] = {}
    wanted = {COHORT: ("capture_table.csv", "point_table.csv"), CALIBRATION: ("calibration_table.csv",)}
    for profile in PROFILES:
        settings, loaded, summary, checked, hashes = _profile_inputs(profile_dir(out_dir, profile), wanted[profile])
        if settings is not None:
            aac.verify_export_files(export_dir, settings)
        settings_by_profile[profile], summaries[profile], preflights[profile] = settings, summary, checked
        tables[profile] = loaded
        inputs.update(hashes)
    reproduced = None
    path = profile_dir(out_dir, COHORT) / REPRODUCTION_NAME
    if path.exists() and settings_by_profile.get(COHORT) is not None and summaries.get(COHORT) is not None:
        stored = json.loads(path.read_text(encoding="utf-8"))
        if stored.get("settings_hash") == lsc.settings_hash(settings_by_profile[COHORT]):
            reproduced = stored
            inputs[f"{COHORT}/{REPRODUCTION_NAME}"] = plsc.sha256_file(path)
    result = decision(settings_by_profile, tables[COHORT].get("capture_table.csv", []),
                      tables[CALIBRATION].get("calibration_table.csv", []), summaries, reproduced)
    result["complete"] = {profile: None if summaries.get(profile) is None else {
        key: summaries[profile][key] for key in ("expected_units", "complete_units", "missing_units", "failed_units",
                                                 "inconsistent_units")} for profile in PROFILES}
    result["preflight"] = {profile: None if preflights.get(profile) is None else {
        "passed": preflights[profile].get("passed"), "units": preflights[profile].get("units"),
        **({"reliability": preflights[profile]["reliability"]} if "reliability" in preflights[profile] else {})}
        for profile in PROFILES}
    result["inputs"] = inputs
    lsc._write_json(out_dir / DECISION_NAME, result)
    return result


# ── CLI ──────────────────────────────────────────────────────────────────────


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = parser.add_subparsers(dest="command", required=True)

    def common(command: argparse.ArgumentParser, *, profile: bool = True) -> None:
        if profile:
            command.add_argument("--profile", choices=PROFILES, required=True)
        command.add_argument("--export-dir", type=Path, required=True,
                             help="The ADNI export, read in place (adni_sexdiff/results/motco_input)")
        command.add_argument("--out-dir", type=Path, required=True)

    def counts(command: argparse.ArgumentParser) -> None:
        command.add_argument("--reps", type=int, default=None,
                             help=f"Cohort: units per point (default {REPLICATES})")
        command.add_argument("--cohorts", type=int, default=None,
                             help=f"Calibration: synthetic cohorts per strength (default {CALIBRATION_COHORTS})")
        command.add_argument("--units", type=int, default=None,
                             help=f"Calibration: units per cohort and point (default {CALIBRATION_UNITS})")
        command.add_argument("--splits", type=int, default=None, help=f"Splits per unit (default {SPLITS})")

    pre = sub.add_parser("preflight", help="Write a profile's settings and run its pre-flight")
    common(pre)
    counts(pre)
    run = sub.add_parser("run", help="Evaluate (a shard of) a profile's units")
    common(run)
    counts(run)
    run.add_argument("--shard", type=lsc.parse_shard, default=(0, 1), help="i/N, round-robin over the units")
    run.add_argument("--limit", type=int, default=None, help="Evaluate at most this many pending units")
    run.add_argument("--workers", type=int, default=1, help="Units evaluated in parallel (not part of the settings)")
    merge_parser = sub.add_parser("merge", help="Merge a profile's shards and write its tables")
    common(merge_parser)
    decide_parser = sub.add_parser("decide", help=f"Write {DECISION_NAME}")
    common(decide_parser, profile=False)
    return parser


def _load(directory: Path) -> aac.Cohort:
    try:
        return aac.load_cohort(directory, aac.CONTRACTS["adni"])
    except pc.ExportContractError as exc:
        raise aac.contract_exit(exc) from exc


def settings_from_args(args: argparse.Namespace, cohort: aac.Cohort) -> dict[str, Any]:
    if args.profile == COHORT:
        if args.cohorts is not None or args.units is not None:
            raise SystemExit("--cohorts and --units are the calibration's; the cohort profile takes --reps.")
        return cohort_settings(cohort, replicates=args.reps, splits=args.splits)
    if args.reps is not None:
        raise SystemExit("--reps is the cohort profile's; the calibration takes --cohorts and --units.")
    return calibration_settings(cohort, cohorts=args.cohorts, units=args.units, splits=args.splits)


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    if args.command in ("preflight", "run"):
        cohort = _load(args.export_dir)
        settings = settings_from_args(args, cohort)
        if args.command == "preflight":
            result = preflight(args.out_dir, settings, cohort)
            print(json.dumps({"profile": args.profile, "passed": result["passed"], "units": [
                {key: unit.get(key) for key in ("unit", "passed", "fidelity", "forced_rank", "statistics_differing",
                                                "rank", "worst_deviation", "failing", "error") if key in unit}
                for unit in result["units"]
            ], **({"reliability": result["reliability"]} if "reliability" in result else {})}, indent=1))
            return 0 if result["passed"] else 1
        index, count = args.shard
        evaluated = run_shard(args.out_dir, settings, cohort, index, count, limit=args.limit, workers=args.workers)
        print(f"evaluated {evaluated} unit(s)", file=sys.stderr)
        return 0
    if args.command == "merge":
        print(json.dumps(merge(args.out_dir, args.profile, args.export_dir), indent=1))
        return 0
    result = decide(args.out_dir, args.export_dir)
    print(json.dumps({
        "verdict": result["verdict"],
        "gate": {key: result["gate"][key] for key in ("passed", "causes")},
        "program": result["program"],
        "resolution": result["resolution"],
        "indicative": result["indicative"],
        "complete": result["complete"],
    }, indent=1))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
