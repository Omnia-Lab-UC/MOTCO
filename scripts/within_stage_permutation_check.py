#!/usr/bin/env python3
"""Within-stage permutation draws against RRPP's unrestricted null: a paired pilot (no R).

The driver behind ``openspec/changes/within-stage-permutation/`` (spec
``within-stage-permutation-validation``, roadmap open item 1). RRPP permutes the
reduced model's residuals across all samples; the opt-in permutation strata
(``RRPP(..., strata=)``, ``permutation_strata="stage"`` in the harness and in
``motco.analysis``, 0.25.0) exchange them only within each stage. This driver
reads the restricted null against the unrestricted one on matched units under a
decision fixed before any unit runs, and re-reads ``predictive_orthonormal``
against production under the same draws.

**Arms** (one production fit per fitted unit; every arm at the unit's RRPP
seed): ``pls_cv`` (the scores) and the three predictive-subspace coordinate
systems of the same fit (``predictive``, ``predictive_orthonormal``,
``predictive_fitted``), each under the unrestricted draws and, with the suffix
``__within_stage``, under the package's within-stage draws. In simulation each
arm is one harness evaluation on the prebuilt fit; on the plasmode each arm is
one ``analyze_trajectories`` call at the case-study configuration (six
covariates, the ``df`` correction), forced at the plasmode's fidelity rank. The
``dispersion`` profile has two arms, ``matrix`` and ``matrix__within_stage``,
read by the RRPP primitive on one Gaussian matrix.

**Checks on every fitted unit.** The pilot's five (``rank``, ``fitted_labels``,
``same_fit``, ``one_subspace``, ``identity_at_low_rank``) on the unrestricted
arms; ``same_statistics`` (each within-stage arm's observed statistics equal its
counterpart's exactly); and, on the plasmode, ``construction_identity`` (each
within-stage arm's p-values equal the predictive size-excess diagnosis's
within-stage loop on that arm's latent matrix, models, correction and seed,
exactly). A failing unit is ``inconsistent`` and leaves every rate.

**Profiles** (output under ``--out-dir``/``<profile>/``):

- ``simulation`` — the predictive-subspace pilot's seven stage-span designs,
  eleven cells and surgery cell (``predictive_subspace_check.simulation_settings``,
  imported unchanged) at a fresh base seed 2038, 100 × 199, p < α; 7,800 units.
- ``weak_size`` — the pilot's 1,000 weak-design control units (500 × 999,
  p < α). The four unrestricted arms must reproduce the pilot's committed
  ``weak_size/rate_table.csv`` exactly.
- ``plasmode`` — the SEA-AD astrocyte export through the plasmode's contract,
  the anchor at the ``paper`` stream's replicates 0–499 (cell ``anchor``) and
  500–999 (cell ``anchor_fresh``, label draws no committed run has read), 999
  permutations, p ≤ α, every arm at the committed fidelity rank. On ``anchor``
  the unrestricted arms must reproduce the pilot's committed plasmode counts and
  the within-stage arms the diagnosis's, exactly. ``preflight`` checks the
  cross-validated fit on replicates 0 and 500 against the forced call first.
- ``dispersion`` — Gaussian latent matrices on the SEA-AD cell table (F 11 / 10
  / 28, M 9 / 10 / 12) with no group difference and a stage-dependent noise
  variance: five patterns (``equal``; the diagnosis's median relative dispersion
  by stage in ``predictive_orthonormal`` and in the scores; each with its
  multipliers squared) at 2 and 5 dimensions, 1,000 × 999, p < α; 10,000 units
  with no fit.

**Parents** (sha256 in the settings; ``run``, ``merge`` and ``decide`` refuse a
changed file by name): the pilot's three ``settings.json`` and its
``weak_size`` and ``plasmode`` ``rate_table.csv``; the diagnosis's
``settings.json`` and ``rate_table.csv``; and, through the pilot's settings
builders, the comparison's, the span-split run's, the plasmode's and the
production selection's files.

Subcommands:

- ``preflight`` (plasmode) — writes ``plasmode/settings.json`` and
  ``plasmode/preflight.json``; ``run --profile plasmode`` waits on its pass.
- ``run --profile …`` — evaluates (a shard of) the units into
  ``records/shard_<i>_of_<N>.jsonl``; round-robin shards, resume of complete
  units, refusal of a resume under other settings naming the keys.
- ``merge --profile …`` — ``merged.jsonl`` (keep off the repository), the rate,
  paired and cell tables, ``accounting.json`` and ``reproduction.json``.
- ``decide`` — ``within_stage_permutation_decision.json`` (rule version 1, design
  D12): the gate, then the null on the production scores (size in four families,
  protected power, specificity), then ``predictive_orthonormal__within_stage``
  against ``pls_cv__within_stage`` by the pilot's S1–S3 and G; verdict
  ``reproduction_failed`` | ``size_excess`` | ``power_loss`` |
  ``specificity_excess`` | ``within_stage_validated`` | ``propose_paper_grade``.

With ``--workers`` above 1, pin BLAS to one thread per process
(``OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1``).

Examples:

    export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
    EXPORT=../sead_sexdiff/results/motco_input/cohort_30
    python scripts/within_stage_permutation_check.py preflight --export-dir $EXPORT --out-dir /tmp/wsp
    python scripts/within_stage_permutation_check.py run --profile plasmode --export-dir $EXPORT --reps 2 \\
        --fresh-reps 2 --out-dir /tmp/wsp --workers 7
    python scripts/within_stage_permutation_check.py run --profile dispersion --reps 20 --out-dir /tmp/wsp
    python scripts/within_stage_permutation_check.py run --profile simulation --reps 2 --out-dir /tmp/wsp
    python scripts/within_stage_permutation_check.py merge --profile plasmode --out-dir /tmp/wsp
    python scripts/within_stage_permutation_check.py decide --out-dir /tmp/wsp

    # one SLURM array task per shard (scripts/within_stage_permutation_array.sbatch), then merge and decide
"""

from __future__ import annotations

import argparse
import dataclasses
import hashlib
import json
import math
import sys
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path
from time import perf_counter
from typing import Any, Iterable, Mapping, Sequence

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
import predictive_size_excess_check as pse  # noqa: E402  (the within-stage loop, the inference path, pre-flight)
import predictive_subspace_check as psc  # noqa: E402  (settings, units, arm records, checks, decision helpers)

sys.path.pop(0)

lsc, pc, plsc, grade = psc.lsc, psc.pc, psc.plsc, psc.grade

DRIVER_NAME = "within-stage-permutation"
DRIVER_VERSION = 1
RULE_VERSION = 1
PROFILES: tuple[str, ...] = ("simulation", "weak_size", "plasmode", "dispersion")
FITTED_PROFILES: tuple[str, ...] = ("simulation", "weak_size", "plasmode")

WITHIN_STAGE = "__within_stage"
STRATA = "stage"
REFERENCE_ARM = psc.REFERENCE_ARM
PREDICTIVE_ARMS: tuple[str, ...] = psc.CANDIDATE_ARMS
UNRESTRICTED_ARMS: tuple[str, ...] = (REFERENCE_ARM, *PREDICTIVE_ARMS)
WITHIN_STAGE_ARMS: tuple[str, ...] = tuple(f"{arm}{WITHIN_STAGE}" for arm in UNRESTRICTED_ARMS)
ARMS: tuple[str, ...] = (*UNRESTRICTED_ARMS, *WITHIN_STAGE_ARMS)
#: The null on the production scores (layer 1) and the candidate under the same draws (layer 2).
NULL_ARM = f"{REFERENCE_ARM}{WITHIN_STAGE}"
CANDIDATE_ARM = f"predictive_orthonormal{WITHIN_STAGE}"
RECORDED_CANDIDATE_ARMS: tuple[str, ...] = (f"predictive{WITHIN_STAGE}", f"predictive_fitted{WITHIN_STAGE}")
DISPERSION_ARM = "matrix"
DISPERSION_WITHIN_ARM = f"{DISPERSION_ARM}{WITHIN_STAGE}"
DISPERSION_ARMS: tuple[str, ...] = (DISPERSION_ARM, DISPERSION_WITHIN_ARM)
PROFILE_ARMS: dict[str, tuple[str, ...]] = {
    "simulation": ARMS, "weak_size": ARMS, "plasmode": ARMS, "dispersion": DISPERSION_ARMS,
}
#: Paired comparisons written by ``merge``: (arm, reference arm).
PAIRS: tuple[tuple[str, str], ...] = (
    *((f"{arm}{WITHIN_STAGE}", arm) for arm in UNRESTRICTED_ARMS),
    *((arm, REFERENCE_ARM) for arm in PREDICTIVE_ARMS),
    *((f"{arm}{WITHIN_STAGE}", NULL_ARM) for arm in PREDICTIVE_ARMS),
    (DISPERSION_WITHIN_ARM, DISPERSION_ARM),
)

STATISTICS: tuple[str, ...] = psc.STATISTICS
TARGETS: dict[str, str] = psc.TARGETS
ALPHA = psc.ALPHA
REJECTION_RULES: dict[str, str] = {**psc.REJECTION_RULES, "dispersion": "p < alpha"}
CONSISTENCY_TOLERANCE = psc.CONSISTENCY_TOLERANCE
CONSISTENCY_CHECKS: tuple[str, ...] = (*psc.CONSISTENCY_CHECKS, "same_statistics", "construction_identity")
ORIENTATION_ESTIMATOR = "net_displacement"

REPO = Path(__file__).resolve().parents[1]
PILOT_DIR = REPO / "results" / "predictive-subspace-2026-11-12"
DIAGNOSIS_DIR = REPO / "results" / "predictive-size-excess-2026-11-13"

# Simulation profile (design D9).
BASE_SEED = 2038
SIMULATION_REPLICATES = psc.SIMULATION_REPLICATES
SIMULATION_PERMUTATIONS = psc.SIMULATION_PERMUTATIONS
STRONG_DESIGNS: tuple[str, ...] = psc.STRONG_DESIGNS
BAND_DESIGNS: tuple[str, ...] = psc.BAND_DESIGNS
WEAK_DESIGN = psc.WEAK_DESIGN
SURGERY_DESIGN = psc.SURGERY_DESIGN
STAGE_SPAN_DESIGNS: tuple[str, ...] = psc.STAGE_SPAN_DESIGNS
CONTROL_CELLS: tuple[str, ...] = psc.CONTROL_CELLS
TARGET_CELLS: tuple[str, ...] = psc.TARGET_CELLS

# Weak-size profile.
WEAK_REPLICATES = psc.WEAK_REPLICATES
WEAK_PERMUTATIONS = psc.WEAK_PERMUTATIONS

# Plasmode profile (design D10).
PLASMODE_DESIGN = psc.PLASMODE_DESIGN
ANCHOR_CELL = psc.ANCHOR_CELL
FRESH_CELL = "anchor_fresh"
PLASMODE_CELLS: tuple[str, ...] = (ANCHOR_CELL, FRESH_CELL)
#: The committed runs read the ``paper`` stream's replicates 0–499; the fresh draws start here.
FRESH_START = psc.PLASMODE_ANCHOR_REPLICATES
PLASMODE_REPLICATES = psc.PLASMODE_ANCHOR_REPLICATES
PLASMODE_FRESH_REPLICATES = 500
PLASMODE_PERMUTATIONS = psc.PLASMODE_PERMUTATIONS
RESIDUAL_CORRECTION = "df"
#: The pilot's plasmode settings this run must share (replicates aside).
_PILOT_SHARED_KEYS: tuple[str, ...] = tuple(key for key in pse._PILOT_SHARED_KEYS if key != "arms")
PREFLIGHT_NAME = "preflight.json"

# Dispersion profile (design D11).
DISPERSION_DESIGN = "dispersion"
DISPERSION_TABLE: list[list[int]] = [[11, 10, 28], [9, 10, 12]]
#: Variance multipliers by stage: the diagnosis's median relative dispersion by stage (its report, section 9).
_MEASURED_PATTERNS: dict[str, tuple[float, float, float]] = {
    "predictive": (1.120, 1.278, 0.802),
    "scores": (0.805, 0.921, 1.137),
}
DISPERSION_PATTERNS: dict[str, tuple[float, ...]] = {
    "equal": (1.0, 1.0, 1.0),
    "predictive": _MEASURED_PATTERNS["predictive"],
    "scores": _MEASURED_PATTERNS["scores"],
    "predictive_doubled": tuple(value**2 for value in _MEASURED_PATTERNS["predictive"]),
    "scores_doubled": tuple(value**2 for value in _MEASURED_PATTERNS["scores"]),
}
DISPERSION_DIMENSIONS: tuple[int, ...] = (2, 5)
#: Population stage share of variance at d = 2 under equal dispersion (the anchor's in-sample stage R²).
DISPERSION_STAGE_SHARE = 0.84
#: Unit centroids before scaling: stage 0, 1 and 2 in the first two dimensions.
DISPERSION_CENTROID_SHAPE: tuple[tuple[float, float], ...] = ((0.0, 0.0), (1.0, 0.0), (0.5, 0.9))
DISPERSION_REPLICATES = 1000
DISPERSION_PERMUTATIONS = 999


def dispersion_cells() -> list[str]:
    return [f"{pattern}_d{dimension}" for pattern in DISPERSION_PATTERNS for dimension in DISPERSION_DIMENSIONS]


# Decision (design D12).
LOSS_ALPHA = psc.LOSS_ALPHA
FAMILY_ALPHA = psc.FAMILY_ALPHA
NULL_CRITERIA: tuple[str, ...] = ("size", "power", "specificity")
_NULL_VERDICT_OF = {"size": "size_excess", "power": "power_loss", "specificity": "specificity_excess"}
VERDICTS: tuple[str, ...] = (
    "reproduction_failed", "size_excess", "power_loss", "specificity_excess", "within_stage_validated",
    "propose_paper_grade",
)
CONSEQUENCES: dict[str, str] = {
    "reproduction_failed": (
        "No reading is interpreted. The cause is fixed and the affected units are rerun under the same rule."
    ),
    "size_excess": (
        "The within-stage null does not hold size in the named family. The option stays an opt-in diagnostic, "
        "documented as such; no case study reports it as a test, and predictive_orthonormal is not reopened."
    ),
    "power_loss": (
        "The null holds size and costs power on the production scores. It stays an opt-in null for data whose stage "
        "dispersion is unequal, documented with the cells it loses; predictive_orthonormal is not reopened."
    ),
    "specificity_excess": (
        "The null holds size and reads more off-target rejections than the bound allows. It stays opt-in, "
        "documented with those readings; predictive_orthonormal is not reopened."
    ),
    "within_stage_validated": (
        "The within-stage null is a validated opt-in null for the production measurement. A separate change may "
        "apply it to the SEA-AD case study and decide the default. predictive_orthonormal stays an opt-in diagnostic "
        "for the named reason."
    ),
    "propose_paper_grade": (
        "As within_stage_validated, and a separate change runs predictive_orthonormal under within-stage draws at "
        "500 x 999 with a rank rule for a fixed measured dimension. Only its verdict may change a default."
    ),
}
DECISION_NAME = "within_stage_permutation_decision.json"


# ── Shared helpers ───────────────────────────────────────────────────────────


def _canonical(value: Any) -> Any:
    return json.loads(json.dumps(value, sort_keys=True))


def _seed(*parts: Any) -> int:
    encoded = json.dumps([DRIVER_NAME, *parts], separators=(",", ":")).encode("utf-8")
    return int(hashlib.sha256(encoded).hexdigest()[:8], 16) & 0x7FFFFFFF


def profile_dir(out_dir: Path, profile: str) -> Path:
    return Path(out_dir) / profile


def unit_key(design: str, cell: str, replicate: int) -> str:
    return psc.unit_key(design, cell, replicate)


def record_key(record: Mapping[str, Any]) -> str:
    return unit_key(record["design"], record["cell"], record["replicate"])


def _arm_spec(arm: str) -> dict[str, str]:
    base = arm.removesuffix(WITHIN_STAGE)
    coordinates = "scores" if base == REFERENCE_ARM else base
    return {"coordinates": coordinates, "draws": "within_stage" if arm.endswith(WITHIN_STAGE) else "unrestricted"}


def _counterpart(arm: str) -> str:
    return arm.removesuffix(WITHIN_STAGE)


# ── Settings and units ───────────────────────────────────────────────────────


def _common(profile: str) -> dict[str, Any]:
    return {
        "driver": DRIVER_NAME,
        "driver_version": DRIVER_VERSION,
        "profile": profile,
        "arms": list(PROFILE_ARMS[profile]),
        "arm_specs": {arm: _arm_spec(arm) for arm in PROFILE_ARMS[profile]},
        "permutation_strata": STRATA,
        "alpha": ALPHA,
        "rejection": REJECTION_RULES[profile],
        "consistency_tolerance": CONSISTENCY_TOLERANCE,
    }


def simulation_settings(
    *,
    replicates: int = SIMULATION_REPLICATES,
    permutations: int = SIMULATION_PERMUTATIONS,
    designs: Sequence[str] = psc.SIMULATION_DESIGNS,
    comparison_dir: Path = psc.COMPARISON_DIR,
    pilot_dir: Path = PILOT_DIR / "simulation",
) -> dict[str, Any]:
    """The pilot's designs, cells and seeds at the fresh base seed 2038, refused if the pilot's would differ."""

    pilot = psc.simulation_settings(parent_dir=comparison_dir, replicates=replicates, permutations=permutations,
                                    designs=designs, base_seed=BASE_SEED)
    committed = plsc.load_parent(Path(pilot_dir))["settings"]
    reference = psc.simulation_settings(parent_dir=comparison_dir, replicates=int(committed["replicates"]),
                                        permutations=int(committed["permutations"]))
    if _canonical(reference) != _canonical(committed):
        differing = sorted(key for key in set(reference) | set(committed)
                           if _canonical(reference.get(key)) != _canonical(committed.get(key)))
        raise SystemExit(f"The pilot's simulation settings ({pilot_dir}/settings.json) no longer rebuild from its "
                         f"code and parents; differing: {', '.join(differing)}.")
    return {
        **_common("simulation"),
        "pilot": pilot,
        "parents": {**pilot["parents"], "pilot_simulation_settings": psc._parent(Path(pilot_dir) / "settings.json")},
        "base_seed": BASE_SEED,
        "replicates": int(replicates),
        "permutations": int(permutations),
    }


def weak_size_settings(
    *,
    replicates: int = WEAK_REPLICATES,
    permutations: int = WEAK_PERMUTATIONS,
    pilot_dir: Path = PILOT_DIR / "weak_size",
    span_split_dir: Path = psc.SPAN_SPLIT_DIR,
    comparison_dir: Path = psc.COMPARISON_DIR,
) -> dict[str, Any]:
    """The pilot's weak-design control units, refused if they would not be the pilot's committed units."""

    pilot = psc.weak_size_settings(replicates=replicates, permutations=permutations, parent_dir=span_split_dir,
                                   comparison_dir=comparison_dir)
    committed = plsc.load_parent(Path(pilot_dir))["settings"]
    reference = psc.weak_size_settings(replicates=int(committed["replicates"]),
                                       permutations=int(committed["permutations"]), parent_dir=span_split_dir,
                                       comparison_dir=comparison_dir)
    if _canonical(reference) != _canonical(committed):
        raise SystemExit(f"The pilot's weak-size settings ({pilot_dir}/settings.json) no longer rebuild from its code "
                         "and parents; refusing.")
    return {
        **_common("weak_size"),
        "pilot": pilot,
        "parents": {
            **pilot["parents"],
            "pilot_weak_size_settings": psc._parent(Path(pilot_dir) / "settings.json"),
            "pilot_weak_size_rates": psc._parent(Path(pilot_dir) / "rate_table.csv"),
        },
        "replicates": int(replicates),
        "permutations": int(permutations),
    }


def plasmode_settings(
    export: Any,
    *,
    replicates: int = PLASMODE_REPLICATES,
    fresh_replicates: int = PLASMODE_FRESH_REPLICATES,
    permutations: int = PLASMODE_PERMUTATIONS,
    pilot_dir: Path = PILOT_DIR / "plasmode",
    diagnosis_dir: Path = DIAGNOSIS_DIR,
    parent_dir: Path = psc.PLASMODE_DIR,
    production_dir: Path = psc.PRODUCTION_PLASMODE_DIR,
) -> dict[str, Any]:
    """The anchor on reproduced (0–499) and fresh (500–999) label draws, refused on another export or pilot."""

    if not 0 <= int(replicates) <= FRESH_START or not 0 <= int(fresh_replicates) <= PLASMODE_FRESH_REPLICATES:
        raise SystemExit(f"--reps must be in 0..{FRESH_START} and --fresh-reps in 0..{PLASMODE_FRESH_REPLICATES}.")
    pilot_dir, diagnosis_dir = Path(pilot_dir), Path(diagnosis_dir)
    pilot = psc.plasmode_settings(export, parent_dir=parent_dir, production_dir=production_dir,
                                  replicates=FRESH_START, plant_replicates=0, permutations=permutations)
    committed = plsc.load_parent(pilot_dir)["settings"]
    differing = [key for key in _PILOT_SHARED_KEYS if _canonical(pilot.get(key)) != _canonical(committed.get(key))]
    if differing:
        raise SystemExit(f"The plasmode units differ from the pilot's ({pilot_dir}/settings.json) in: "
                         f"{', '.join(differing)}; they would not be the units the excess was read on.")
    diagnosis = plsc.load_parent(diagnosis_dir)["settings"]
    expected = {"residual_correction": RESIDUAL_CORRECTION, "orientation_estimator": ORIENTATION_ESTIMATOR,
                "anchor_rank": int(pilot["fidelity"]["selected_lv"])}
    mismatched = [key for key, value in expected.items() if diagnosis.get(key) != value]
    if _canonical((diagnosis.get("pilot") or {}).get("plasmode")) != _canonical(pilot["plasmode"]):
        mismatched.append("pilot.plasmode")
    if mismatched:
        raise SystemExit(f"The diagnosis's settings ({diagnosis_dir}/settings.json) differ in: {', '.join(mismatched)}"
                         "; its within-stage counts would not be these units'.")
    return {
        **_common("plasmode"),
        "pilot": pilot,
        "parents": {
            **pilot["parents"],
            "pilot_plasmode_settings": psc._parent(pilot_dir / "settings.json"),
            "pilot_plasmode_rates": psc._parent(pilot_dir / "rate_table.csv"),
            "diagnosis_settings": psc._parent(diagnosis_dir / "settings.json"),
            "diagnosis_rates": psc._parent(diagnosis_dir / "rate_table.csv"),
        },
        "design": PLASMODE_DESIGN,
        "cells": {ANCHOR_CELL: [0, int(replicates)], FRESH_CELL: [FRESH_START, FRESH_START + int(fresh_replicates)]},
        "anchor_rank": int(pilot["fidelity"]["selected_lv"]),
        "fit": "every arm at forced_components = the plasmode's fidelity rank (no cross-validation)",
        "permutations": int(permutations),
        "residual_correction": RESIDUAL_CORRECTION,
        "orientation_estimator": ORIENTATION_ESTIMATOR,
        "preflight_replicates": [0] + ([FRESH_START] if int(fresh_replicates) > 0 else []),
    }


def dispersion_centroid_scale(
    table: Sequence[Sequence[int]] = DISPERSION_TABLE, share: float = DISPERSION_STAGE_SHARE, dimension: int = 2,
) -> float:
    """The centroid scale s giving the population stage share ``share`` at ``dimension`` with unit noise variance.

    The stage weights are the table's stage sizes; the between-stage variance of
    the unit centroids grows as s², and the within-stage variance is
    ``dimension``.
    """

    sizes = np.asarray(table, dtype=float).sum(axis=0)
    weights = sizes / sizes.sum()
    centroids = np.asarray(DISPERSION_CENTROID_SHAPE, dtype=float)
    centre = weights @ centroids
    between = float(weights @ np.sum((centroids - centre) ** 2, axis=1))
    return math.sqrt(share * dimension / ((1.0 - share) * between))


def dispersion_settings(
    *, replicates: int = DISPERSION_REPLICATES, permutations: int = DISPERSION_PERMUTATIONS,
) -> dict[str, Any]:
    scale = dispersion_centroid_scale()
    return {
        **_common("dispersion"),
        "design": DISPERSION_DESIGN,
        "cell_table": DISPERSION_TABLE,
        "patterns": {name: list(values) for name, values in DISPERSION_PATTERNS.items()},
        "dimensions": list(DISPERSION_DIMENSIONS),
        "cells": dispersion_cells(),
        "centroids": {
            "shape": [list(row) for row in DISPERSION_CENTROID_SHAPE],
            "scale": scale,
            "stage_share_at_d2": DISPERSION_STAGE_SHARE,
            "rule": "stage s has centroid scale * shape[s] in the first two dimensions and zero elsewhere",
        },
        "noise": "N(0, v_stage I_d), v_stage the pattern's multiplier",
        "base_seed": BASE_SEED,
        "seed_derivation": f"sha256 of [{DRIVER_NAME!r}, {BASE_SEED}, 'dispersion', pattern, dimension, replicate, "
                           "purpose] (purpose 'matrix' or 'rrpp')",
        "orientation_estimator": ORIENTATION_ESTIMATOR,
        "residual_correction": None,
        "replicates": int(replicates),
        "permutations": int(permutations),
    }


def simulation_units(settings: Mapping[str, Any]) -> list[lsc.Unit]:
    return psc.simulation_units(settings["pilot"])


def weak_size_units(settings: Mapping[str, Any]) -> list[lsc.Unit]:
    return psc.weak_size_units(settings["pilot"])


def plasmode_units(settings: Mapping[str, Any]) -> list[lsc.Unit]:
    return [
        lsc.Unit(PLASMODE_DESIGN, cell, replicate)
        for cell in PLASMODE_CELLS
        for replicate in range(*settings["cells"][cell])
    ]


def dispersion_units(settings: Mapping[str, Any]) -> list[lsc.Unit]:
    return [
        lsc.Unit(DISPERSION_DESIGN, cell, replicate)
        for replicate in range(int(settings["replicates"]))
        for cell in settings["cells"]
    ]


def enumerate_units(settings: Mapping[str, Any]) -> list[lsc.Unit]:
    profile = settings["profile"]
    if profile == "simulation":
        return simulation_units(settings)
    if profile == "weak_size":
        return weak_size_units(settings)
    if profile == "plasmode":
        return plasmode_units(settings)
    return dispersion_units(settings)


def verify(settings: Mapping[str, Any], export: Any = None) -> None:
    """Refuse a changed parent file, or an export whose files differ from the settings', by name."""

    plsc.verify_parents(settings)
    if export is not None and settings["profile"] == "plasmode":
        plsc.verify_export(export, settings["pilot"])


# ── Per-unit checks ──────────────────────────────────────────────────────────


def _statistics_deviation(left: Mapping[str, Any], right: Mapping[str, Any]) -> float:
    worst = 0.0
    for statistic in STATISTICS:
        a, b = left.get(statistic), right.get(statistic)
        if a is None and b is None:
            continue
        if a is None or b is None:
            return math.inf
        worst = max(worst, abs(float(a) - float(b)))
    return worst


def same_statistics(arms: Mapping[str, Mapping[str, Any]]) -> dict[str, Any]:
    """Each within-stage arm's observed statistics against its counterpart's, exactly."""

    pairs = [
        (arm, _counterpart(arm)) for arm in WITHIN_STAGE_ARMS
        if (arms.get(arm) or {}).get("status") == "completed"
        and (arms.get(_counterpart(arm)) or {}).get("status") == "completed"
    ]
    deviation = max(
        (_statistics_deviation(arms[arm]["pair_statistics"], arms[source]["pair_statistics"]) for arm, source in pairs),
        default=None,
    )
    return psc._check(deviation, 0.0, applicable=bool(pairs))


def _finish(base: dict[str, Any], arms: dict[str, Any], checks: dict[str, Any], production: dict[str, Any],
            start: float, timing: dict[str, float]) -> dict[str, Any]:
    failing = [name for name in CONSISTENCY_CHECKS if not checks[name]["passed"]]
    timing["unit"] = perf_counter() - start
    return {
        **base,
        "status": "inconsistent" if failing else "completed",
        "inconsistent_checks": failing,
        "consistency": checks,
        "production": production,
        "arms": {arm: arms[arm] for arm in ARMS if arm in arms},
        "runtime_seconds": timing,
    }


def _with_spec(arm: str, record: dict[str, Any]) -> dict[str, Any]:
    return {**record, **_arm_spec(arm)}


# ── One harness unit (simulation, weak_size) ─────────────────────────────────


def evaluate_harness_unit(unit: lsc.Unit, settings: Mapping[str, Any]) -> dict[str, Any]:
    """One dataset and one production fit; the four coordinate systems under both draws at the unit's seed."""

    from motco.simulations.evaluation import (
        evaluate_semisynthetic_trajectory,
        integrate_semisynthetic_dataset,
        pls_coordinates_from_pls,
    )

    start = perf_counter()
    pilot = settings["pilot"]
    units = psc._harness_units(pilot)
    cell = lsc.cell_spec(units, unit.cell)
    base_seed = int(units["base_seed"])
    seed = lsc.rrpp_seed(base_seed, unit.design, unit.cell, unit.replicate)
    n_stages = int(units["designs"][unit.design]["generator"]["n_stages"])
    base = {
        "profile": settings["profile"],
        "design": unit.design,
        "cell": unit.cell,
        "trajectory_mode": cell["mode"],
        "effect_size": float(cell["effect"]),
        "replicate": unit.replicate,
        "generator_seed": lsc.generator_seed(base_seed, unit.design, unit.replicate),
        "rrpp_seed": seed,
        "settings_hash": lsc.settings_hash(settings),
    }
    timing: dict[str, float] = {}
    try:
        mark = perf_counter()
        dataset = lsc.unit_dataset(unit, units, reference=psc._reference())
        timing["generate"] = perf_counter() - mark
        params = lsc.arm_params(lsc.REFERENCE_ARM, units, seed=seed, n_stages=n_stages)
        mark = perf_counter()
        scores = integrate_semisynthetic_dataset(dataset, params)
        timing["fit"] = perf_counter() - mark
        model = scores.artifacts.model
        production = {
            "rank": int(scores.metadata["selected_lv"]),
            "cv_mean_auroc": scores.metadata.get("cv_mean_auroc"),
            "loadings_rank": psc.loadings_rank(model),
            "stage_r2": psc.stage_r2(
                model, scores.artifacts.features.to_numpy(dtype=float), dataset.metadata[params.stage_col].tolist()
            ),
        }
    except Exception as exc:  # recorded, never aborts the shard
        return psc._failed_unit(base, exc, start, timing)

    arms: dict[str, Any] = {}
    candidates: dict[str, np.ndarray] = {}
    for source in UNRESTRICTED_ARMS:
        source_params = params if source == REFERENCE_ARM else dataclasses.replace(params, latent_coordinates=source)
        try:
            latent = scores if source == REFERENCE_ARM else pls_coordinates_from_pls(dataset, scores, source_params)
        except Exception as exc:  # both draws of this coordinate system fail on the unit
            for arm in (source, f"{source}{WITHIN_STAGE}"):
                arms[arm] = _with_spec(arm, psc._failed_arm(exc, 0.0))
            continue
        if source != REFERENCE_ARM:
            candidates[source] = latent.matrix.to_numpy(dtype=float)
        for strata in (None, STRATA):
            arm = source if strata is None else f"{source}{WITHIN_STAGE}"
            mark = perf_counter()
            try:
                result = evaluate_semisynthetic_trajectory(
                    dataset, dataclasses.replace(source_params, permutation_strata=strata), latent=latent
                )
                arms[arm] = _with_spec(arm, psc.harness_arm_record(result, perf_counter() - mark))
            except Exception as exc:  # recorded on the arm, the unit goes on
                arms[arm] = _with_spec(arm, psc._failed_arm(exc, perf_counter() - mark))
    checks = psc.consistency_checks(
        model=model,
        scores=scores.matrix.to_numpy(dtype=float),
        n_stages=n_stages,
        candidates=candidates,
        arms={arm: arms[arm] for arm in UNRESTRICTED_ARMS if arm in arms},
        with_capture=bool(params.span_capture),
    )
    checks["same_statistics"] = same_statistics(arms)
    checks["construction_identity"] = psc._check(None, 0.0, applicable=False)
    return _finish(base, arms, checks, production, start, timing)


# ── One plasmode unit ────────────────────────────────────────────────────────


def construction_identity(
    results: Mapping[str, Any], arms: Mapping[str, Mapping[str, Any]], dataset: Any, params: Any, *,
    permutations: int, seed: int,
) -> dict[str, Any]:
    """Each within-stage arm's p-values against the diagnosis's within-stage loop on its latent matrix, exactly."""

    from motco import analysis

    inputs = analysis._trajectory_inputs(dataset, params)
    order = [str(index) for index in range(len(inputs.stage_levels))]
    strata = inputs.working[analysis._STAGE].astype(str).to_numpy()
    indices = pse.within_stage_indices(strata, order, int(permutations), int(seed))
    differing: dict[str, list[str]] = {}
    compared = []
    for source in UNRESTRICTED_ARMS:
        arm = f"{source}{WITHIN_STAGE}"
        if source not in results or (arms.get(arm) or {}).get("status") != "completed":
            continue
        inferred = pse.infer(results[source].latent_scores, inputs, permutations=int(permutations), seed=int(seed),
                             indices=indices, orientation_estimator=ORIENTATION_ESTIMATOR,
                             residual_correction=RESIDUAL_CORRECTION)
        theirs, ours = lsc._finite(inferred["p_values"]), arms[arm]["p_values"]
        compared.append(arm)
        names = sorted(name for name in set(theirs) | set(ours) if theirs.get(name) != ours.get(name))
        if names:
            differing[arm] = names
    check = psc._check(math.inf if differing else 0.0, 0.0, applicable=bool(compared))
    return {**check, "arms": compared, "differing": differing}


def evaluate_plasmode_unit(unit: lsc.Unit, settings: Mapping[str, Any], export: Any) -> dict[str, Any]:
    """Build the unit's label draw once; analyse each coordinate system under both draws at the fidelity rank."""

    start = perf_counter()
    pilot = settings["pilot"]
    parent = pilot["plasmode"]
    permutations = int(settings["permutations"])
    seeds = pc.unit_seeds(int(parent["base_seed"]), str(parent["stream"]), unit.replicate)
    rank = int(settings["anchor_rank"])
    base: dict[str, Any] = {
        "profile": "plasmode",
        "design": unit.design,
        "cell": unit.cell,
        "trajectory_mode": "none",
        "effect_size": 0.0,
        "replicate": unit.replicate,
        "stream": parent["stream"],
        "seeds": seeds,
        "rrpp_seed": seeds["rrpp"],
        "settings_hash": lsc.settings_hash(settings),
    }
    timing: dict[str, float] = {}
    try:
        mark = perf_counter()
        built = plsc.build_plasmode_unit(lsc.Unit(PLASMODE_DESIGN, ANCHOR_CELL, unit.replicate), pilot, export)
        base["label_redraws"] = built["label_redraws"]
        timing["generate"] = perf_counter() - mark
        production_params = psc.plasmode_production_params(pilot, built["seeds"])
        n_stages = len(export.stage_order)
    except Exception as exc:  # recorded, never aborts the shard
        return lsc._jsonable(psc._failed_unit(base, exc, start, timing))

    arms: dict[str, Any] = {}
    results: dict[str, Any] = {}
    dataset = None
    for source in UNRESTRICTED_ARMS:
        source_params = (
            dataclasses.replace(production_params, forced_components=rank)
            if source == REFERENCE_ARM
            else psc.plasmode_arm_params(source, production_params, n_stages=n_stages, production_rank=rank)
        )
        for strata in (None, STRATA):
            arm = source if strata is None else f"{source}{WITHIN_STAGE}"
            mark = perf_counter()
            try:
                built_dataset, result = pc.analyse(
                    built["blocks"], built["metadata"], export,
                    dataclasses.replace(source_params, permutation_strata=strata),
                )
                record = psc.analysis_arm_record(result, perf_counter() - mark)
                record.pop("span_capture", None)
                record.pop("population_latent", None)
                arms[arm] = _with_spec(arm, {**record, "fit": "forced_at_fidelity_rank"})
                if strata is None:
                    results[source] = result
                    dataset = built_dataset
            except Exception as exc:  # recorded on the arm, the unit goes on
                arms[arm] = _with_spec(arm, psc._failed_arm(exc, perf_counter() - mark))
    if REFERENCE_ARM not in results:
        error = arms.get(REFERENCE_ARM, {})
        return lsc._jsonable(psc._failed_unit(
            base, RuntimeError(f"the production analysis failed: {error.get('error_message')}"), start, timing
        ))
    production_result = results[REFERENCE_ARM]
    artifacts = production_result.latent_artifacts
    model = artifacts.model
    production = {
        "rank": int(production_result.record["latent"]["selected_lv"]),
        "cv_mean_auroc": production_result.record["latent"].get("cv_mean_auroc"),
        "loadings_rank": psc.loadings_rank(model),
        "stage_r2": psc.stage_r2(
            model, artifacts.features.to_numpy(dtype=float), built["metadata"][export.stage_col].tolist()
        ),
    }
    checks = psc.consistency_checks(
        model=model,
        scores=production_result.latent_scores.to_numpy(dtype=float),
        n_stages=n_stages,
        candidates={arm: results[arm].latent_scores.to_numpy(dtype=float) for arm in PREDICTIVE_ARMS if arm in results},
        arms={arm: arms[arm] for arm in UNRESTRICTED_ARMS if arm in arms},
        with_capture=False,
        models={arm: results[arm].latent_artifacts.model for arm in PREDICTIVE_ARMS if arm in results},
    )
    checks["same_statistics"] = same_statistics(arms)
    mark = perf_counter()
    try:
        checks["construction_identity"] = construction_identity(
            results, arms, dataset, production_params, permutations=permutations, seed=int(seeds["rrpp"])
        )
    except Exception as exc:  # a construction that raises fails its own identity
        checks["construction_identity"] = {**psc._check(math.inf, 0.0), "error": f"{type(exc).__name__}: {exc}"}
    timing["construction_identity"] = perf_counter() - mark
    return lsc._jsonable(_finish(base, arms, checks, production, start, timing))


# ── One dispersion unit ──────────────────────────────────────────────────────


def _dispersion_cell(cell: str) -> tuple[str, int]:
    pattern, _, dimension = cell.rpartition("_d")
    if pattern not in DISPERSION_PATTERNS or int(dimension) not in DISPERSION_DIMENSIONS:
        raise KeyError(f"Unknown dispersion cell {cell!r}.")
    return pattern, int(dimension)


def dispersion_seeds(cell: str, replicate: int) -> dict[str, int]:
    pattern, dimension = _dispersion_cell(cell)
    return {
        purpose: _seed(BASE_SEED, "dispersion", pattern, dimension, int(replicate), purpose)
        for purpose in ("matrix", "rrpp")
    }


def dispersion_metadata(table: Sequence[Sequence[int]] = DISPERSION_TABLE) -> pd.DataFrame:
    """One row per sample, group-major and stage-minor, with string labels."""

    rows = [(str(group), str(stage)) for group, sizes in enumerate(table) for stage, size in enumerate(sizes)
            for _ in range(int(size))]
    return pd.DataFrame(rows, columns=["group", "stage"])


def dispersion_matrix(settings: Mapping[str, Any], cell: str, replicate: int) -> tuple[np.ndarray, pd.DataFrame]:
    """The unit's Gaussian latent matrix: stage centroids plus noise of the pattern's variance by stage."""

    pattern, dimension = _dispersion_cell(cell)
    metadata = dispersion_metadata(settings["cell_table"])
    stages = metadata["stage"].astype(int).to_numpy()
    centroids = np.zeros((len(settings["cell_table"][0]), dimension))
    shape = np.asarray(settings["centroids"]["shape"], dtype=float)
    centroids[:, :2] = float(settings["centroids"]["scale"]) * shape
    variances = np.asarray(settings["patterns"][pattern], dtype=float)
    rng = np.random.default_rng(dispersion_seeds(cell, replicate)["matrix"])
    noise = rng.standard_normal((stages.size, dimension)) * np.sqrt(variances[stages])[:, None]
    return centroids[stages] + noise, metadata


def evaluate_dispersion_unit(unit: lsc.Unit, settings: Mapping[str, Any]) -> dict[str, Any]:
    """Both nulls by the RRPP primitive on one matrix at one seed; no fit."""

    from motco.simulations.evaluation import (
        _empirical_p_value,
        _extract_null_distributions,
        _extract_pair_statistics,
        _summarize_null_distributions,
        build_simulation_trajectory_design,
        stage_strata,
    )
    from motco.stats.permutation import RRPP
    from motco.stats.trajectory import estimate_difference

    start = perf_counter()
    pattern, dimension = _dispersion_cell(unit.cell)
    seeds = dispersion_seeds(unit.cell, unit.replicate)
    base = {
        "profile": "dispersion",
        "design": unit.design,
        "cell": unit.cell,
        "pattern": pattern,
        "dimension": dimension,
        "trajectory_mode": "none",
        "effect_size": 0.0,
        "replicate": unit.replicate,
        "seeds": seeds,
        "rrpp_seed": seeds["rrpp"],
        "settings_hash": lsc.settings_hash(settings),
    }
    timing: dict[str, float] = {}
    try:
        Y, metadata = dispersion_matrix(settings, unit.cell, unit.replicate)
        design = build_simulation_trajectory_design(metadata, group_col="group", stage_col="stage")
        estimator = str(settings["orientation_estimator"])
        observed = _extract_pair_statistics(
            *estimate_difference(Y, design.model_full, design.ls_means, design.contrast,
                                 orientation_estimator=estimator),
            shape_available=True,
        )
        arms: dict[str, Any] = {}
        for arm, strata in ((DISPERSION_ARM, None),
                            (DISPERSION_WITHIN_ARM, stage_strata(metadata, "stage", design.stage_levels))):
            mark = perf_counter()
            null = _extract_null_distributions(
                *RRPP(Y, design.model_full, design.model_reduced, design.ls_means, design.contrast,
                      permutations=int(settings["permutations"]), n_jobs=1, progress=False, seed=seeds["rrpp"],
                      orientation_estimator=estimator, strata=strata),
                shape_available=True,
            )
            p_values = {statistic: _empirical_p_value(null[statistic], value)
                        for statistic, value in observed.items() if statistic in null and np.isfinite(value)}
            arms[arm] = {
                "status": "completed",
                "draws": "within_stage" if strata is not None else "unrestricted",
                "pair_statistics": lsc._finite(observed),
                "p_values": lsc._finite(p_values),
                "angle_null_q95": (_summarize_null_distributions(null).get("angle") or {}).get("q95"),
                "seconds": perf_counter() - mark,
            }
    except Exception as exc:  # recorded, never aborts the shard
        return lsc._jsonable(psc._failed_unit(base, exc, start, timing))
    timing["unit"] = perf_counter() - start
    return lsc._jsonable({**base, "status": "completed", "inconsistent_checks": [], "arms": arms,
                          "runtime_seconds": timing})


# ── Pre-flight (plasmode) ────────────────────────────────────────────────────


def preflight(out_dir: Path, settings: Mapping[str, Any], export: Any) -> dict[str, Any]:
    """Write the plasmode settings and ``preflight.json``: the CV against the forced fit on replicates 0 and 500."""

    directory = profile_dir(out_dir, "plasmode")
    if (directory / "settings.json").exists():
        verify(lsc.load_settings(directory))
    verify(settings, export)
    lsc.ensure_settings(directory, settings)
    view = {"pilot": settings["pilot"], "anchor_rank": settings["anchor_rank"]}
    units = []
    for replicate in settings["preflight_replicates"]:
        try:
            units.append(pse.preflight_unit(int(replicate), view, export))
        except Exception as exc:  # a pre-flight that raises fails
            units.append({"unit": unit_key(PLASMODE_DESIGN, ANCHOR_CELL, int(replicate)), "passed": False,
                          "error": f"{type(exc).__name__}: {exc}"})
    result = {
        "settings_hash": lsc.settings_hash(settings),
        "rule": "the cross-validated pls_cv selects the fidelity rank and AUROC, and its statistics and p-values "
                "equal the forced call's exactly",
        "expected": {"rank": int(settings["anchor_rank"]), **settings["pilot"]["fidelity"]},
        "units": units,
        "passed": bool(units) and all(unit["passed"] for unit in units),
    }
    lsc._write_json(directory / PREFLIGHT_NAME, result)
    return result


def load_preflight(out_dir: Path, settings: Mapping[str, Any]) -> dict[str, Any] | None:
    path = profile_dir(out_dir, "plasmode") / PREFLIGHT_NAME
    if not path.exists():
        return None
    record = json.loads(path.read_text(encoding="utf-8"))
    return record if record.get("settings_hash") == lsc.settings_hash(settings) else None


# ── Shards, resume and the run ───────────────────────────────────────────────


def read_records(paths: Iterable[Path]) -> list[dict[str, Any]]:
    return psc.read_records(paths)


_WORKER_EXPORT: Any = None


def _init_worker(directory: str, files: Mapping[str, str]) -> None:
    global _WORKER_EXPORT
    export = pc.load_export(Path(directory))
    if export.files != dict(files):
        raise SystemExit(f"{directory} changed after the run started (file hashes differ).")
    _WORKER_EXPORT = export


def evaluate_unit(unit: lsc.Unit, settings: Mapping[str, Any], export: Any = None) -> dict[str, Any]:
    profile = settings["profile"]
    if profile in ("simulation", "weak_size"):
        return evaluate_harness_unit(unit, settings)
    if profile == "plasmode":
        return evaluate_plasmode_unit(unit, settings, export if export is not None else _WORKER_EXPORT)
    return evaluate_dispersion_unit(unit, settings)


def run_shard(
    out_dir: Path,
    settings: Mapping[str, Any],
    index: int,
    count: int,
    *,
    export: Any = None,
    limit: int | None = None,
    workers: int = 1,
) -> int:
    """Evaluate this shard's missing units into its file; returns how many were evaluated."""

    profile = settings["profile"]
    directory = profile_dir(out_dir, profile)
    if (directory / "settings.json").exists():
        verify(lsc.load_settings(directory))
    if profile == "plasmode" and export is None:
        raise SystemExit("The plasmode profile needs the export.")
    verify(settings, export)
    lsc.ensure_settings(directory, settings)
    if profile == "plasmode":
        checked = load_preflight(out_dir, settings)
        if checked is None:
            raise SystemExit(f"{directory / PREFLIGHT_NAME} is missing or was written under other settings; run "
                             "`preflight` first.")
        if not checked.get("passed"):
            raise SystemExit(f"The pre-flight failed ({directory / PREFLIGHT_NAME}); the run stops.")
    path = lsc.shard_path(directory, index, count)
    path.parent.mkdir(parents=True, exist_ok=True)
    done = {record_key(record): record for record in (read_records([path]) if path.exists() else [])}
    # Rewritten with complete units only, so a torn line is re-run, not duplicated.
    with path.open("w", encoding="utf-8") as handle:
        handle.write("".join(lsc._line(done[key]) for key in sorted(done)))
    pending = [
        unit
        for position, unit in enumerate(enumerate_units(settings))
        if position % count == index and unit_key(*unit.key) not in done
    ]
    if limit is not None:
        pending = pending[:limit]

    def write(position: int, unit: lsc.Unit, record: Mapping[str, Any]) -> None:
        with path.open("a", encoding="utf-8") as handle:
            handle.write(lsc._line(record))
            handle.flush()
        print(f"[{profile} {index}/{count}] {position}/{len(pending)} {unit.key} {record['status']}",
              file=sys.stderr, flush=True)

    if workers <= 1:
        for position, unit in enumerate(pending, start=1):
            write(position, unit, evaluate_unit(unit, settings, export))
        return len(pending)
    initializer: Any = None
    initargs: tuple[Any, ...] = ()
    if profile == "plasmode":
        initializer, initargs = _init_worker, (str(export.directory), dict(export.files))
    with ProcessPoolExecutor(max_workers=workers, initializer=initializer, initargs=initargs) as pool:
        futures = {pool.submit(evaluate_unit, unit, settings): unit for unit in pending}
        for position, future in enumerate(as_completed(futures), start=1):
            write(position, futures[future], future.result())
    return len(pending)


# ── Merge: rate, paired and cell tables, accounting, reproduction ────────────


RATE_COLUMNS = psc.RATE_COLUMNS
PAIRED_COLUMNS = psc.PAIRED_COLUMNS
CELL_COLUMNS = psc.CELL_COLUMNS


def rate_table(records: Sequence[Mapping[str, Any]], settings: Mapping[str, Any]) -> list[dict[str, Any]]:
    return psc.rate_table(records, settings)


def paired_table(records: Sequence[Mapping[str, Any]], settings: Mapping[str, Any]) -> list[dict[str, Any]]:
    """Every pair of :data:`PAIRS` present in the profile, per cell and statistic, over units available in both."""

    alpha, rule = float(settings["alpha"]), str(settings["rejection"])
    arms = set(settings["arms"])
    rows = []
    for _, group in sorted(psc._by_cell(records).items()):
        fields = psc._cell_fields(group[0])
        ordered = sorted(group, key=lambda record: int(record["replicate"]))
        for statistic in STATISTICS:
            for arm, reference in PAIRS:
                if arm not in arms or reference not in arms:
                    continue
                decisions = [
                    (psc._rejects(psc._p(record, reference, statistic), alpha, rule),
                     psc._rejects(psc._p(record, arm, statistic), alpha, rule))
                    for record in ordered
                ]
                both = [(ref, cand) for ref, cand in decisions if ref is not None and cand is not None]
                n = len(both)
                b = sum(1 for ref, cand in both if ref and not cand)
                c = sum(1 for ref, cand in both if cand and not ref)
                rows.append({
                    **fields, "statistic": statistic, "arm": arm, "reference_arm": reference, "n_paired": n,
                    "rate_reference": sum(ref for ref, _ in both) / n if n else None,
                    "rate_arm": sum(cand for _, cand in both) / n if n else None,
                    "difference": (c - b) / n if n else None,
                    "b": b, "c": c, **plsc.paired_p_values(b, c),
                })
    return rows


def cell_table(records: Sequence[Mapping[str, Any]], settings: Mapping[str, Any]) -> list[dict[str, Any]]:
    return psc.cell_table(records, settings)


def accounting(
    records: Sequence[Mapping[str, Any]], settings: Mapping[str, Any], *, expected: int, duplicates: int
) -> dict[str, Any]:
    summary = psc.accounting(records, settings, expected=expected, duplicates=duplicates)
    summary["settings_hash"] = lsc.settings_hash(settings)
    return summary


def _counts(rows: Sequence[Mapping[str, Any]], *, design: str, cells: Sequence[str],
            arms: Sequence[str]) -> dict[tuple[str, ...], dict[str, int]]:
    return {
        (str(row["cell"]), str(row["arm"]), str(row["statistic"])): {
            "count": int(float(row["count"])), "rejections": int(float(row["rejections"]))
        }
        for row in rows
        if str(row["design"]) == design and str(row["cell"]) in cells and str(row["arm"]) in arms
    }


def _compare(label: str, keys: Sequence[tuple[str, ...]], committed: Mapping[tuple[str, ...], Any],
             ours: Mapping[tuple[str, ...], Any], *, source: Mapping[str, Any], rule: str) -> dict[str, Any]:
    return psc._compare(label, keys, committed, ours, source=source, rule=rule)


def weak_size_reproduction(rates: Sequence[Mapping[str, Any]], settings: Mapping[str, Any],
                           committed_rows: Sequence[Mapping[str, Any]]) -> dict[str, Any]:
    """The four unrestricted arms' counts against the pilot's committed weak-size rate table, exact."""

    cells = list(CONTROL_CELLS)
    keys = [(cell, arm, statistic) for cell in cells for arm in UNRESTRICTED_ARMS for statistic in STATISTICS]
    check = _compare(
        "pilot_weak_size", keys,
        _counts(committed_rows, design=WEAK_DESIGN, cells=cells, arms=UNRESTRICTED_ARMS),
        _counts(rates, design=WEAK_DESIGN, cells=cells, arms=UNRESTRICTED_ARMS),
        source=settings["parents"]["pilot_weak_size_rates"],
        rule="each unrestricted arm's count and rejections (p < alpha) equal the pilot's weak_size/rate_table.csv",
    )
    return {"profile": "weak_size", "passed": check["passed"], "checks": {"pilot_weak_size": check}}


def plasmode_reproduction(rates: Sequence[Mapping[str, Any]], settings: Mapping[str, Any],
                          pilot_rows: Sequence[Mapping[str, Any]],
                          diagnosis_rows: Sequence[Mapping[str, Any]]) -> dict[str, Any]:
    """Cell ``anchor``: the unrestricted arms against the pilot, the within-stage arms against the diagnosis."""

    cells = [ANCHOR_CELL]
    pilot = _compare(
        "pilot_anchor",
        [(ANCHOR_CELL, arm, statistic) for arm in UNRESTRICTED_ARMS for statistic in STATISTICS],
        _counts(pilot_rows, design=PLASMODE_DESIGN, cells=cells, arms=UNRESTRICTED_ARMS),
        _counts(rates, design=PLASMODE_DESIGN, cells=cells, arms=UNRESTRICTED_ARMS),
        source=settings["parents"]["pilot_plasmode_rates"],
        rule="each unrestricted arm's anchor count and rejections (p <= alpha) equal the pilot's plasmode rate table",
    )
    diagnosis = _compare(
        "diagnosis_anchor",
        [(ANCHOR_CELL, arm, statistic) for arm in WITHIN_STAGE_ARMS for statistic in STATISTICS],
        _counts(diagnosis_rows, design=PLASMODE_DESIGN, cells=cells, arms=WITHIN_STAGE_ARMS),
        _counts(rates, design=PLASMODE_DESIGN, cells=cells, arms=WITHIN_STAGE_ARMS),
        source=settings["parents"]["diagnosis_rates"],
        rule="each within-stage arm's anchor count and rejections (p <= alpha) equal the diagnosis's rate table",
    )
    checks = {"pilot_anchor": pilot, "diagnosis_anchor": diagnosis}
    return {"profile": "plasmode", "passed": all(check["passed"] for check in checks.values()), "checks": checks}


def merge(out_dir: Path, profile: str) -> dict[str, Any]:
    """Merge one profile's shards and write its tables, accounting and reproduction."""

    directory = profile_dir(out_dir, profile)
    settings = lsc.load_settings(directory)
    verify(settings)
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
    expected = [unit_key(*unit.key) for unit in enumerate_units(settings)]
    records = [by_key[key] for key in expected if key in by_key]
    with (directory / "merged.jsonl").open("w", encoding="utf-8") as handle:
        handle.write("".join(lsc._line(record) for record in records))
    rates = rate_table(records, settings)
    lsc._write_csv(directory / "rate_table.csv", rates, RATE_COLUMNS)
    lsc._write_csv(directory / "paired_table.csv", paired_table(records, settings), PAIRED_COLUMNS)
    lsc._write_csv(directory / "cell_table.csv", cell_table(records, settings), CELL_COLUMNS)
    summary = accounting(records, settings, expected=len(expected), duplicates=sum(n - 1 for n in seen.values()))
    rates_read = lsc._read_csv(directory / "rate_table.csv")
    if profile == "weak_size":
        reproduction = weak_size_reproduction(
            rates_read, settings, lsc._read_csv(plsc._resolve(settings["parents"]["pilot_weak_size_rates"]["path"]))
        )
    elif profile == "plasmode":
        reproduction = plasmode_reproduction(
            rates_read,
            settings,
            lsc._read_csv(plsc._resolve(settings["parents"]["pilot_plasmode_rates"]["path"])),
            lsc._read_csv(plsc._resolve(settings["parents"]["diagnosis_rates"]["path"])),
        )
    else:
        reproduction = {"profile": profile, "passed": True, "checks": {},
                        "note": f"a fresh seed family (base seed {BASE_SEED}): there is no committed run to reproduce"}
    lsc._write_json(directory / "accounting.json", summary)
    lsc._write_json(directory / "reproduction.json", reproduction)
    return {
        "profile": profile,
        **{key: summary[key] for key in ("complete_units", "missing_units", "failed_units", "inconsistent_units",
                                         "duplicate_records", "recorded_core_hours")},
        "reproduction_passed": reproduction["passed"],
        "reproduction_differing": {name: len(check.get("differing", []))
                                   for name, check in reproduction["checks"].items()},
    }


# ── Decision (pure functions, design D12) ────────────────────────────────────

#: A reading: (profile, design, cell, statistic), as in the pilot.
Reading = tuple[str, str, str, str]


def criteria_plan() -> dict[str, list[Reading]]:
    """The readings of every criterion, fixed before any unit.

    Size: the anchor and translation of the seven stage-span designs (m = 42),
    ``weak_size`` (m = 6), the plasmode's ``anchor_fresh`` (m = 3) and the ten
    dispersion cells (m = 30). Protected power: the target cells of the six
    designs at δ ≥ 0.25 (m = 54). Specificity: the target cells' off-target
    statistics at δ = 2.0 (m = 36) and in the band (m = 72). Gain: the target
    cells at δ = 2.0 (m = 18); recorded only: the band (36) and the weak design
    (9).
    """

    pilot = psc.criteria_plan()
    return {
        "S1_simulation": pilot["S1_simulation"],
        "S1_weak_size": pilot["S1_weak_size"],
        "S1_plasmode": [("plasmode", PLASMODE_DESIGN, FRESH_CELL, statistic) for statistic in STATISTICS],
        "S1_dispersion": [
            ("dispersion", DISPERSION_DESIGN, cell, statistic)
            for cell in dispersion_cells()
            for statistic in STATISTICS
        ],
        "S2": [reading for reading in pilot["S2"] if reading[0] == "simulation"],
        "S3_strong": pilot["S3_strong"],
        "S3_band": pilot["S3_band"],
        "G": pilot["G"],
        "G_band": pilot["G_band"],
        "G_weak": pilot["G_weak"],
    }


def _index_rates(tables: Mapping[str, Mapping[str, Sequence[Mapping[str, Any]]]]) -> dict[tuple[str, ...], Any]:
    return psc._index_rates(tables)


def _index_paired(tables: Mapping[str, Mapping[str, Sequence[Mapping[str, Any]]]]) -> dict[tuple[str, ...], Any]:
    return psc._index_paired(tables)


def judge_band_excess(rows: Sequence[Mapping[str, Any]]) -> dict[str, Any]:
    """The band's specificity for a null: a reading counts only above the size bound with a paired excess.

    A reading counts when its rejected count exceeds the stage-span bound
    α + 2·SE at its own count and its one-sided exact paired excess p-value
    against the unrestricted null is below 0.025. The family holds when no more
    readings count than the 95th percentile of Binomial(m, 0.025), no reading
    above the bound has an excess p-value below 0.05/m, and none is unavailable.
    """

    m = len(rows)
    allowed = psc.allowed_count(m) if m else None
    ceiling = FAMILY_ALPHA / m if m else None
    out, counted, over_ceiling, unavailable = [], 0, 0, 0
    for row in rows:
        entry = dict(row)
        if row.get("status") == "unavailable" or row.get("n_available") is None:
            unavailable += 1
            entry["status"] = "unavailable"
            out.append(entry)
            continue
        n, hits = int(row["n_available"]), int(row["n_rejected"])
        above = hits > grade.exceedance_threshold(n)
        p = float(row["p_gain"])
        counts = above and p < LOSS_ALPHA
        beyond = above and ceiling is not None and p < ceiling
        counted += counts
        over_ceiling += beyond
        entry.update({"rate": hits / n, "bound": grade.bound(n), "above_bound": above,
                      "status": "below_ceiling" if beyond else ("excess" if counts else "holds")})
        out.append(entry)
    holds = unavailable == 0 and allowed is not None and counted <= allowed and over_ceiling == 0
    return {
        "rule": "count_and_ceiling_on_paired_excess_above_bound",
        "holds": holds,
        "m": m,
        "alpha": LOSS_ALPHA,
        "allowed": allowed,
        "observed": counted,
        "ceiling": ceiling,
        "n_below_ceiling": over_ceiling,
        "n_unavailable": unavailable,
        "observations": out,
    }


def judge_null(plan: Mapping[str, Sequence[Reading]], rates: Mapping[tuple[str, ...], Any],
               paired: Mapping[tuple[str, ...], Any]) -> dict[str, Any]:
    """Layer 1: size, protected power and specificity of ``pls_cv__within_stage`` against ``pls_cv``."""

    size = {
        family: grade.judge_rates([psc._observation(rates, reading, NULL_ARM) for reading in plan[f"S1_{family}"]])
        for family in ("simulation", "weak_size", "plasmode")
    }
    size["dispersion"] = grade.judge_rates(
        [psc._observation(rates, reading, DISPERSION_WITHIN_ARM) for reading in plan["S1_dispersion"]]
    )
    power = psc.judge_paired(
        [psc._paired(paired, reading, NULL_ARM, REFERENCE_ARM) for reading in plan["S2"]],
        p_field="p_loss", rule="count_and_ceiling_on_paired_loss",
    )
    strong = grade.judge_rates([psc._observation(rates, reading, NULL_ARM) for reading in plan["S3_strong"]])
    band_rows = []
    for reading in plan["S3_band"]:
        observation = psc._observation(rates, reading, NULL_ARM)
        pair = psc._paired(paired, reading, NULL_ARM, REFERENCE_ARM)
        unavailable = observation.get("status") == "unavailable" or pair.get("status") == "unavailable"
        band_rows.append({**pair, **({"status": "unavailable"} if unavailable else {
            "n_available": observation["n_available"], "n_rejected": observation["n_rejected"]})})
    band = judge_band_excess(band_rows)
    criteria = {
        "size": {"holds": all(family["holds"] for family in size.values()), "families": size},
        "power": power,
        "specificity": {"holds": strong["holds"] and band["holds"], "families": {"strong": strong, "band": band}},
    }
    gains = {
        family: psc.judge_gain([psc._paired(paired, reading, NULL_ARM, REFERENCE_ARM) for reading in plan[key]])
        for family, key in (("strong", "G"), ("band", "G_band"), ("weak", "G_weak"))
    }
    failing = [name for name in NULL_CRITERIA if not criteria[name]["holds"]]
    return {
        "arm": NULL_ARM,
        "reference_arm": REFERENCE_ARM,
        "dispersion_arm": DISPERSION_WITHIN_ARM,
        "verdict": _NULL_VERDICT_OF[failing[0]] if failing else "validated",
        "failing_criteria": failing,
        "criteria": criteria,
        "recorded_gains": gains,
    }


def judge_candidate(
    arm: str,
    reference: str,
    plan: Mapping[str, Sequence[Reading]],
    rates: Mapping[tuple[str, ...], Any],
    paired: Mapping[tuple[str, ...], Any],
) -> dict[str, Any]:
    """The pilot's S1, S2, S3 and G for ``arm`` against ``reference``; every criterion is evaluated.

    With ``reference = "pls_cv"`` and the pilot's plan this is
    ``predictive_subspace_check.judge_arm`` (a test pins it on the pilot's
    committed tables). The S1 families are those of the plan's ``S1_<family>``
    keys other than ``dispersion``; the recorded gains those of its ``G_<family>``
    keys.
    """

    families = [key.removeprefix("S1_") for key in plan if key.startswith("S1_") and key != "S1_dispersion"]
    s1 = {
        family: grade.judge_rates([psc._observation(rates, reading, arm) for reading in plan[f"S1_{family}"]])
        for family in families
    }
    strong = grade.judge_rates([psc._observation(rates, reading, arm) for reading in plan["S3_strong"]])
    band_rows = []
    for reading in plan["S3_band"]:
        row = psc._paired(paired, reading, arm, reference)
        production = psc._observation(rates, reading, reference)
        row["production_rate"] = (
            None if production.get("status") == "unavailable"
            else production["n_rejected"] / production["n_available"]
        )
        band_rows.append(row)
    band = psc.judge_paired(band_rows, p_field="p_gain", rule="count_and_ceiling_on_paired_excess")
    off_target = [psc._observation(rates, reading, arm) for reading in (*plan["S3_strong"], *plan["S3_band"])]
    above = psc._above_bound(off_target)
    criteria = {
        "S1": {"holds": all(family["holds"] for family in s1.values()), "families": s1},
        "S2": psc.judge_paired(
            [psc._paired(paired, reading, arm, reference) for reading in plan["S2"]],
            p_field="p_loss",
            rule="count_and_ceiling_on_paired_loss",
        ),
        "S3": {
            "holds": strong["holds"] and band["holds"],
            "families": {"strong": strong, "band": band},
            "off_target_above_bound": len(above),
            "off_target_above_bound_readings": above,
        },
        "G": psc.judge_gain([psc._paired(paired, reading, arm, reference) for reading in plan["G"]]),
    }
    recorded_gains = {
        key.removeprefix("G_"): psc.judge_gain([psc._paired(paired, reading, arm, reference) for reading in plan[key]])
        for key in plan if key.startswith("G_")
    }
    failing_readings = {
        "S1": [
            {"family": family, **row}
            for family, judged in s1.items()
            for row in judged["observations"]
            if row.get("status") != "holds"
        ],
        "S2": [row for row in criteria["S2"]["observations"] if row.get("status") != "holds"],
        "S3": [
            {"family": family, **row}
            for family, judged in (("strong", strong), ("band", band))
            for row in judged["observations"]
            if row.get("status") != "holds"
        ],
        "G": [] if criteria["G"]["holds"] else criteria["G"]["observations"],
    }
    failing = [name for name in ("S1", "S2", "S3", "G") if not criteria[name]["holds"]]
    return {
        "verdict": psc._VERDICT_OF[failing[0]] if failing else "qualifies",
        "failing_criteria": failing,
        "failing_readings": {name: rows for name, rows in failing_readings.items() if name in failing},
        "criteria": criteria,
        "recorded_gains": recorded_gains,
    }


def gate(
    reproduction: Mapping[str, Mapping[str, Any]],
    accounting_summaries: Mapping[str, Mapping[str, Any]],
    checked: Mapping[str, Any] | None,
) -> dict[str, Any]:
    """``reproduction_failed`` on a failed reproduction or pre-flight, or any inconsistent unit."""

    failures = []
    for profile in PROFILES:
        entry = reproduction.get(profile)
        if entry is None:
            failures.append({"profile": profile, "reason": "the profile was not merged"})
            continue
        for name, check in (entry.get("checks") or {}).items():
            if not check.get("passed"):
                failures.append({"profile": profile, "reason": f"{name} failed"})
        if not entry.get("passed", False) and not any(f["profile"] == profile for f in failures):
            failures.append({"profile": profile, "reason": "reproduction failed"})
        inconsistent = int((accounting_summaries.get(profile) or {}).get("inconsistent_units", 0) or 0)
        if inconsistent:
            failures.append({"profile": profile, "reason": f"{inconsistent} inconsistent unit(s)"})
    if checked is None:
        failures.append({"profile": "plasmode", "reason": "the pre-flight record is missing"})
    elif not checked.get("passed"):
        failures.append({"profile": "plasmode", "reason": "the pre-flight did not pass"})
    return {"passed": not failures, "failures": failures}


def decision(
    tables: Mapping[str, Mapping[str, Sequence[Mapping[str, Any]]]],
    *,
    reproduction: Mapping[str, Mapping[str, Any]],
    accounting_summaries: Mapping[str, Mapping[str, Any]],
    checked: Mapping[str, Any] | None,
    diagnosis_rows: Sequence[Mapping[str, Any]] = (),
) -> dict[str, Any]:
    """The predeclared verdict (rule version 1); the descriptive readings beside it never change it."""

    rates, paired = _index_rates(tables), _index_paired(tables)
    plan = criteria_plan()
    gated = gate(reproduction, accounting_summaries, checked)
    null = judge_null(plan, rates, paired)
    candidate = judge_candidate(CANDIDATE_ARM, NULL_ARM, plan, rates, paired)
    if not gated["passed"]:
        verdict = "reproduction_failed"
    elif null["failing_criteria"]:
        verdict = null["verdict"]
    elif candidate["verdict"] == "qualifies":
        verdict = "propose_paper_grade"
    else:
        verdict = "within_stage_validated"
    return {
        "rule_version": RULE_VERSION,
        "verdict": verdict,
        "consequence": CONSEQUENCES[verdict],
        "consequences": dict(CONSEQUENCES),
        "gate": gated,
        "null": null,
        "candidate": {
            "arm": CANDIDATE_ARM,
            "reference_arm": NULL_ARM,
            "bears_on_verdict": gated["passed"] and not null["failing_criteria"],
            **candidate,
        },
        "alpha": ALPHA,
        "rejection": dict(REJECTION_RULES),
        "plan": {name: [list(reading) for reading in readings] for name, readings in plan.items()},
        "descriptive": descriptive_tables(tables, plan, diagnosis_rows=diagnosis_rows),
        "notes": [
            "Layer 1 judges the within-stage null on the production scores; layer 2 judges predictive_orthonormal "
            "under within-stage draws against production under the same draws, so a gain is the coordinates' and "
            "not the null's. Layer 2 bears on the verdict only when layer 1 holds.",
            "A band off-target reading counts against the null only above the size bound and with a paired excess "
            "against the unrestricted null, so restored size is not read as a leak.",
            "Every table under 'descriptive' is a reading beside the verdict and never changes it.",
        ],
    }


def _rate(rates: Mapping[tuple[str, ...], Any], profile: str, design: str, cell: str, arm: str,
          statistic: str) -> float | None:
    row = rates.get((profile, design, cell, arm, statistic))
    return None if row is None else plsc._number(row.get("rate"))


def _count(rates: Mapping[tuple[str, ...], Any], profile: str, design: str, cell: str, arm: str,
           statistic: str) -> tuple[int, int] | None:
    row = rates.get((profile, design, cell, arm, statistic))
    if row is None:
        return None
    return int(float(row["count"])), int(float(row["rejections"]))


def descriptive_tables(
    tables: Mapping[str, Mapping[str, Sequence[Mapping[str, Any]]]],
    plan: Mapping[str, Sequence[Reading]],
    *,
    diagnosis_rows: Sequence[Mapping[str, Any]] = (),
) -> dict[str, Any]:
    """The other candidates, the dispersion cells, the anchor draws, the null's gains, null widths, plant rows."""

    rates, paired = _index_rates(tables), _index_paired(tables)
    recorded = {arm: judge_candidate(arm, NULL_ARM, plan, rates, paired) for arm in RECORDED_CANDIDATE_ARMS}
    unrestricted_candidate = judge_candidate("predictive_orthonormal", REFERENCE_ARM, plan, rates, paired)
    dispersion = [
        {
            "cell": cell, "pattern": _dispersion_cell(cell)[0], "dimension": _dispersion_cell(cell)[1],
            "statistic": statistic,
            **{f"rate_{arm}": _rate(rates, "dispersion", DISPERSION_DESIGN, cell, arm, statistic)
               for arm in DISPERSION_ARMS},
            "paired": psc._paired(paired, ("dispersion", DISPERSION_DESIGN, cell, statistic), DISPERSION_WITHIN_ARM,
                                  DISPERSION_ARM),
        }
        for cell in dispersion_cells() for statistic in STATISTICS
    ]
    anchor = []
    for arm in ARMS:
        for statistic in STATISTICS:
            counts = {cell: _count(rates, "plasmode", PLASMODE_DESIGN, cell, arm, statistic) for cell in PLASMODE_CELLS}
            present = [value for value in counts.values() if value is not None]
            pooled = (sum(n for n, _ in present), sum(k for _, k in present)) if present else None
            anchor.append({
                "arm": arm, "statistic": statistic,
                **{cell: None if value is None else {"count": value[0], "rejections": value[1]}
                   for cell, value in counts.items()},
                "pooled": None if pooled is None else {"count": pooled[0], "rejections": pooled[1],
                                                       "rate": pooled[1] / pooled[0] if pooled[0] else None},
                "fresh_size_rule": grade.judge_rates([psc._observation(rates, (
                    "plasmode", PLASMODE_DESIGN, FRESH_CELL, statistic), arm)]),
            })
    targets = [*plan["G"], *plan["G_band"], *plan["G_weak"]]
    null_gains = [psc._paired(paired, reading, NULL_ARM, REFERENCE_ARM) for reading in targets]
    surgery = [
        {"arm": arm, "statistic": statistic,
         "rate": _rate(rates, "simulation", SURGERY_DESIGN, psc.SURGERY_CELL[0], arm, statistic)}
        for arm in ARMS for statistic in STATISTICS
    ]
    null_widths = [
        {"profile": profile, "design": str(row["design"]), "cell": str(row["cell"]), "arm": str(row["arm"]),
         "angle_null_q95_median": plsc._number(row.get("angle_null_q95_median"))}
        for profile, table in tables.items() for row in table.get("cells", [])
    ]
    plant_rows = [
        {key: row.get(key) for key in ("cell", "arm", "statistic", "count", "rejections", "rate")}
        for row in diagnosis_rows
        if str(row.get("cell", "")).startswith("span_e")
        and str(row.get("arm")) in (*UNRESTRICTED_ARMS, *WITHIN_STAGE_ARMS)
    ]
    return {
        "recorded_candidates": recorded,
        "predictive_orthonormal_unrestricted_against_pls_cv": unrestricted_candidate,
        "dispersion": dispersion,
        "plasmode_anchor_draws": anchor,
        "null_gains": null_gains,
        "surgery_cell": surgery,
        "null_widths": null_widths,
        "diagnosis_plant_rows": plant_rows,
    }


_TABLE_FILES: dict[str, str] = {"rates": "rate_table.csv", "paired": "paired_table.csv", "cells": "cell_table.csv"}


def _read_json(path: Path) -> dict[str, Any]:
    if not path.exists():
        raise SystemExit(f"{path} does not exist; run `merge` first.")
    return json.loads(path.read_text(encoding="utf-8"))


def decide(out_dir: Path) -> dict[str, Any]:
    """Read every profile's merged tables, reproduction and accounting and the pre-flight; write the decision."""

    out_dir = Path(out_dir)
    tables: dict[str, dict[str, list[dict[str, Any]]]] = {}
    reproduction: dict[str, Any] = {}
    summaries: dict[str, Any] = {}
    inputs: dict[str, str] = {}
    settings_by_profile: dict[str, Any] = {}
    for profile in PROFILES:
        directory = profile_dir(out_dir, profile)
        settings = lsc.load_settings(directory)
        verify(settings)
        settings_by_profile[profile] = settings
        summary = _read_json(directory / "accounting.json")
        if summary.get("settings_hash") != lsc.settings_hash(settings):
            raise SystemExit(f"{directory}/accounting.json was merged under other settings; run `merge` again.")
        summaries[profile] = summary
        reproduction[profile] = _read_json(directory / "reproduction.json")
        tables[profile] = {}
        for name, filename in _TABLE_FILES.items():
            path = directory / filename
            if not path.exists():
                raise SystemExit(f"{path} does not exist; run `merge` first.")
            tables[profile][name] = lsc._read_csv(path)
        for filename in ("settings.json", *_TABLE_FILES.values(), "accounting.json", "reproduction.json"):
            inputs[f"{profile}/{filename}"] = plsc.sha256_file(directory / filename)
    plasmode = settings_by_profile["plasmode"]
    checked = load_preflight(out_dir, plasmode)
    if (profile_dir(out_dir, "plasmode") / PREFLIGHT_NAME).exists():
        inputs[f"plasmode/{PREFLIGHT_NAME}"] = plsc.sha256_file(profile_dir(out_dir, "plasmode") / PREFLIGHT_NAME)
    diagnosis_rows = lsc._read_csv(plsc._resolve(plasmode["parents"]["diagnosis_rates"]["path"]))
    result = decision(tables, reproduction=reproduction, accounting_summaries=summaries, checked=checked,
                      diagnosis_rows=diagnosis_rows)
    completeness = ("expected_units", "complete_units", "missing_units", "failed_units", "inconsistent_units")
    result["complete"] = {
        profile: {key: summaries[profile][key] for key in completeness} for profile in PROFILES
    }
    result["reproduction"] = {
        profile: {name: {key: check.get(key) for key in ("passed", "rows_compared", "rows_matching")}
                  for name, check in (reproduction[profile].get("checks") or {}).items()}
        for profile in PROFILES
    }
    result["preflight"] = None if checked is None else {
        "passed": checked.get("passed"),
        "units": [{key: unit.get(key) for key in ("unit", "passed", "forced_rank", "statistics_differing",
                                                  "p_values_differing")} for unit in checked.get("units", [])],
    }
    result["inputs"] = inputs
    lsc._write_json(out_dir / DECISION_NAME, result)
    return result


# ── CLI ──────────────────────────────────────────────────────────────────────


_DEFAULT_REPS = {
    "simulation": SIMULATION_REPLICATES, "weak_size": WEAK_REPLICATES, "plasmode": PLASMODE_REPLICATES,
    "dispersion": DISPERSION_REPLICATES,
}
_DEFAULT_PERMS = {
    "simulation": SIMULATION_PERMUTATIONS, "weak_size": WEAK_PERMUTATIONS, "plasmode": PLASMODE_PERMUTATIONS,
    "dispersion": DISPERSION_PERMUTATIONS,
}


def _unit_arguments(parser: argparse.ArgumentParser, *, plasmode_only: bool = False) -> None:
    if not plasmode_only:
        parser.add_argument("--profile", choices=PROFILES, required=True)
    parser.add_argument("--out-dir", type=Path, required=True)
    parser.add_argument("--export-dir", type=Path, default=None, help="The SEA-AD cohort_30 export (plasmode)")
    parser.add_argument("--reps", type=int, default=None,
                        help="Replicates (default: simulation 100, weak_size 500, plasmode 500 reproduced, "
                             "dispersion 1000)")
    parser.add_argument("--fresh-reps", type=int, default=PLASMODE_FRESH_REPLICATES,
                        help=f"Plasmode replicates from {FRESH_START} on (default {PLASMODE_FRESH_REPLICATES})")
    parser.add_argument("--perms", type=int, default=None, help="Permutations (default: simulation 199, else 999)")
    parser.add_argument("--designs", type=str, default=None, help="Simulation designs, comma-separated (default all)")


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = parser.add_subparsers(dest="command", required=True)
    pre = sub.add_parser("preflight", help="Write the plasmode settings and check the forced anchor fit")
    _unit_arguments(pre, plasmode_only=True)
    run = sub.add_parser("run", help="Evaluate (a shard of) a profile's units")
    _unit_arguments(run)
    run.add_argument("--shard", type=lsc.parse_shard, default=(0, 1), help="i/N, round-robin over the units")
    run.add_argument("--limit", type=int, default=None, help="Evaluate at most this many pending units")
    run.add_argument("--workers", type=int, default=1, help="Units evaluated in parallel (not part of the settings)")
    merge_parser = sub.add_parser("merge", help="Merge one profile's shards")
    merge_parser.add_argument("--profile", choices=PROFILES, required=True)
    merge_parser.add_argument("--out-dir", type=Path, required=True)
    decide_parser = sub.add_parser("decide", help=f"Write {DECISION_NAME}")
    decide_parser.add_argument("--out-dir", type=Path, required=True)
    return parser


def settings_from_args(args: argparse.Namespace, profile: str) -> tuple[dict[str, Any], Any]:
    reps = _DEFAULT_REPS[profile] if args.reps is None else int(args.reps)
    perms = _DEFAULT_PERMS[profile] if args.perms is None else int(args.perms)
    if profile == "simulation":
        designs = psc.SIMULATION_DESIGNS if args.designs is None else [d for d in args.designs.split(",") if d]
        return simulation_settings(replicates=reps, permutations=perms, designs=designs), None
    if profile == "weak_size":
        return weak_size_settings(replicates=reps, permutations=perms), None
    if profile == "dispersion":
        return dispersion_settings(replicates=reps, permutations=perms), None
    if args.export_dir is None:
        raise SystemExit("The plasmode profile needs --export-dir.")
    try:
        export = pc.load_export(args.export_dir)
    except pc.ExportContractError as exc:
        raise pc.contract_exit(exc) from exc
    return plasmode_settings(export, replicates=reps, fresh_replicates=int(args.fresh_reps), permutations=perms), export


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    if args.command == "preflight":
        settings, export = settings_from_args(args, "plasmode")
        result = preflight(args.out_dir, settings, export)
        print(json.dumps({"passed": result["passed"], "units": [
            {key: unit.get(key) for key in ("unit", "passed", "forced_rank", "statistics_differing",
                                            "p_values_differing", "error")}
            for unit in result["units"]
        ]}, indent=1))
        return 0 if result["passed"] else 1
    if args.command == "run":
        index, count = args.shard
        settings, export = settings_from_args(args, args.profile)
        evaluated = run_shard(args.out_dir, settings, index, count, export=export, limit=args.limit,
                              workers=args.workers)
        print(f"evaluated {evaluated} unit(s)", file=sys.stderr)
        return 0
    if args.command == "merge":
        print(json.dumps(merge(args.out_dir, args.profile), indent=1))
        return 0
    result = decide(args.out_dir)
    print(json.dumps({
        "verdict": result["verdict"],
        "null": {"verdict": result["null"]["verdict"], "failing": result["null"]["failing_criteria"]},
        "candidate": {"verdict": result["candidate"]["verdict"], "failing": result["candidate"]["failing_criteria"]},
        "gate": result["gate"],
        "complete": result["complete"],
    }, indent=1))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
