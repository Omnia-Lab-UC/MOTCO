#!/usr/bin/env python3
"""Predictive-subspace coordinates against the production measurement: a paired pilot (no R).

The driver behind ``openspec/changes/predictive-subspace-coordinates/`` (spec
``predictive-subspace-validation``, roadmap open item 1). Production uses the
cross-validated rank A of the stage-PLS fit both to regularize the fit and as
the number of dimensions ``delta``, ``angle`` and ``shape`` are measured in.
The stage labels are one-hot, so the fit's Y loadings have rank K − 1 whatever
A is, and the fitted stage labels move only along K − 1 directions of the
score space: the predictive subspace. The harness and ``motco.analysis`` offer
three coordinate systems on it (``latent_coordinates``, 0.24.0). This driver
reads them against production on matched units under a decision fixed before
any unit runs.

**Arms** (one production fit per unit; every arm at the unit's RRPP seed with
serial draws, net displacement, rank floor 2):

- ``pls_cv`` — production, in the fit's scores;
- ``predictive``, ``predictive_orthonormal``, ``predictive_fitted`` — the three
  coordinate systems on the predictive subspace of that same fit (the scores
  projected on it, the input coordinates on an orthonormal basis of it, the
  fitted stage indicators), built by ``pls_coordinates_from_pls`` without
  refitting;
- ``forced_k_minus_1`` — the fit forced to K − 1 components, in scores: the
  control that lowers the fit rank instead of the measured dimension;
- ``predictive_rank_plus_one`` — the ``predictive`` coordinates of the fit
  forced to A + 1 components (simulation only): how near the CV rank sits to
  over-aiming.

Each unit checks that the predictive arms measure the production fit: the Y
loadings have rank min(A, K − 1) (``rank``); the fitted stage indicators lie in
every predictive arm's column space to 1e-8 (``fitted_labels``); every predictive
arm's fitted model is production's, at its rank and with rotations and Y loadings
equal to 1e-8 (``same_fit``); their span-capture
shares are equal to 1e-8 (``one_subspace``, with span capture); and at
A ≤ K − 1 the ``predictive`` p-values equal ``pls_cv``'s (``identity_at_low_rank``).
A unit failing a check is recorded as ``inconsistent`` with the check and its
deviation and leaves every rate. A unit that raises is recorded as ``failed``
and the shard goes on; an arm that raises is recorded as failed on its unit.

**Profiles** (every profile's output lands under ``--out-dir``/``<profile>/``):

- ``simulation`` — a fresh seed family: base seed 2036 on the latent-space
  comparison's design generators, seed derivations, ``unit_dataset`` and
  ``cell_spec`` (``latent_space_comparison``, imported unchanged), under
  ``span_scale="planted"``, at 100 replicates × 199 permutations with span
  capture. Seven designs: ``n200`` and ``seaad`` (the SEA-AD table) at the
  generator's δ = 2.0; ``n200_delta1``, ``n200_delta0.5``, ``n200_delta0.25``
  and ``seaad_delta0.5`` (the stage-effect band); ``seaad_weak`` at the
  comparison's committed weak δ (its calibration bound by sha256). Eleven
  cells each: ``anchor``, ``translation`` (e = 1.0), magnitude at
  e ∈ {0.05, 0.1, 0.25}, orientation at e ∈ {0.1, 0.25, 0.5} with φ = 0 and
  shape at e ∈ {0.1, 0.25, 0.5}. One descriptive cell, ``surgery_orient_e1``,
  applies the default feature-set orientation surgery at e = 1.0 on its own
  design entry ``seaad_surgery`` (the SEA-AD table at δ = 2.0 without the
  stage-span construction). 7,800 units; rejection is p < α. Every cell of a
  (design, replicate) shares its null dataset; no dataset seed equals the
  comparison's (base seed 2035) or the probe's (4242).
- ``weak_size`` — the comparison's ``seaad_weak`` anchor and translation at
  replicates 0–499 × 999 permutations (base seed 2035), built by
  ``span_split_check.simulation_settings`` / ``simulation_units`` (imported
  unchanged) and checked against the span-split run's committed settings, so
  every unit is the span-split check's simulation unit. 1,000 units; p < α.
- ``plasmode`` — the SEA-AD astrocyte ``cohort_30`` export read through the
  plasmode's input contract (``plasmode_check``; the file hashes must equal
  the committed plasmode settings'). Labels, seeds and the real-noise span
  plant come from ``production_latent_space_check`` (``build_plasmode_unit``,
  ``plasmode_arm_params``, imported unchanged), so ``pls_cv`` is the
  production selection's arm. The anchor runs the ``paper`` stream's
  replicates 0–499, the plant at e ∈ {0.25, 0.5, 1.0} replicates 0–99, all at
  999 permutations; 800 units. Each arm is one ``analyze_trajectories`` call
  with the case-study parameters plus the arm's setting: ``pls_cv`` runs the CV;
  each predictive arm runs at ``forced_components`` = the unit's ``pls_cv``
  rank with its ``latent_coordinates`` (the amendment of 2026-11-12, before any
  array unit: the fitted model is the same ``fit_plsda_model`` call, so the CV
  is not repeated, and the arm's record carries ``pls_cv``'s AUROC); and
  ``forced_k_minus_1`` at ``forced_components = 2``. ``same_fit`` checks that
  every predictive arm fitted production's model. Rejection is p ≤ α. No record
  holds a donor identifier or a label vector.

**Parents** (sha256 in the settings; ``run``, ``merge`` and ``decide`` refuse a
changed file by name): the comparison's ``settings.json`` and calibration; the
span-split run's ``settings.json`` and ``production_table.csv``; the plasmode's
``settings.json`` and ``plasmode_summary.csv``; the production selection's
plasmode ``settings.json`` and ``rate_table.csv``.

Subcommands:

- ``run --profile simulation|weak_size|plasmode`` — evaluate (a shard of) the
  units into ``records/shard_<i>_of_<N>.jsonl``, one line per unit with every
  arm nested. Shards are round-robin; resume skips complete units;
  ``settings.json`` is bound by its hash and a resume under other settings
  exits naming the keys.
- ``merge --profile …`` — ``merged.jsonl`` (keep off the repository),
  ``rate_table.csv`` (per design, cell, arm and statistic: available and
  rejected counts, rate, Monte Carlo SE, failures), ``paired_table.csv``
  (every non-production arm against ``pls_cv``, and the two control arms
  against ``predictive``: both rates, the discordant counts b (only the
  reference rejects) and c (only the arm), the exact one-sided loss
  P(Bin(b + c, ½) ≥ b), gain P(Bin(b + c, ½) ≥ c) and two-sided p-values, 1
  when b + c = 0), ``cell_table.csv`` (per cell and arm: rank quantiles,
  AUROC, measured dimension, capture shares, in-sample stage R², median
  observed and ``population_latent`` statistics, median ``angle`` null 95th
  percentile), ``accounting.json`` (units, failures and inconsistent units by
  cell, duplicates, recorded core-hours, per-arm timing) and
  ``reproduction.json``: ``weak_size`` ``pls_cv`` against the span-split
  ``production_table.csv``; ``plasmode`` ``pls_cv`` anchor against the
  plasmode's ``plasmode_summary.csv`` (500 replicates), its plant rows against
  the production selection's ``rate_table.csv``, and the fidelity values
  (rank 5, AUROC 0.635556) on every anchor unit; all exact. The simulation
  profile is a fresh family and reproduces nothing.
- ``decide`` — ``predictive_subspace_decision.json`` (rule version 1) from the
  three merged profiles. The gate ``reproduction_failed`` (a failed
  reproduction or fidelity check, or any inconsistent unit; S1–G are still
  recorded), then per candidate arm:

  - **S1 size** — the anchor and translation ``delta``/``angle``/``shape`` rates
    in three families under the stage-span paper grade's count-plus-ceiling
    rule (``stage_span_paper_grade_decision.judge_rates``, imported
    unchanged): the simulation's seven designs (m = 42), ``weak_size``
    (m = 6) and the plasmode anchor (m = 3); all three must hold.
  - **S2 protected power** — the 54 target cells of the six designs at
    δ ≥ 0.25 and the plasmode's three plant ``delta`` cells (m = 57): at most
    the 95th percentile of Binomial(57, 0.025) one-sided exact paired loss
    p-values below 0.025, none below 0.05/57.
  - **S3 specificity** — the target cells' off-target rates: at δ = 2.0
    (m = 36) the bound rule of S1; in the band (m = 72) the one-sided exact
    paired excess p-value against ``pls_cv`` under the count-plus-ceiling rule
    (Binomial(72, 0.025), 0.05/72). Both must hold. The off-target readings
    above the S1 bound are counted over both families.
  - **G gain** — at least one one-sided exact paired gain p-value at most 0.05/18
    among the 18 target cells of the δ = 2.0 designs. The band's 36, the weak
    design's 9 and the plasmode's 3 target cells are recorded, never decisive.

  An unavailable reading fails its criterion. Per arm the first failure in the
  order ``size_excess``, ``power_loss``, ``specificity_excess``, ``no_gain``,
  else ``qualifies``; overall ``propose_paper_grade`` naming the qualifying arm
  with the fewest off-target readings above the bound (ties: ``predictive``,
  ``predictive_orthonormal``, ``predictive_fitted``), else ``keep_production``.
  Descriptive tables beside it never change it: the control arm, the
  rank-plus-one arm, observed against population statistics, the surgery
  cell, the weak design's target cells, production's band off-target rates,
  and rank and AUROC per design.

With ``--workers`` above 1, pin BLAS to one thread per process
(``OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1``).

Examples:

    export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
    EXPORT=../sead_sexdiff/results/motco_input/cohort_30
    # smoke: 2 replicates per cell
    python scripts/predictive_subspace_check.py run --profile simulation --reps 2 --out-dir /tmp/psc --workers 7
    python scripts/predictive_subspace_check.py run --profile weak_size --reps 2 --out-dir /tmp/psc
    python scripts/predictive_subspace_check.py run --profile plasmode --export-dir $EXPORT --reps 2 \\
        --plant-reps 2 --out-dir /tmp/psc --workers 7
    python scripts/predictive_subspace_check.py merge --profile simulation --out-dir /tmp/psc
    python scripts/predictive_subspace_check.py merge --profile weak_size --out-dir /tmp/psc
    python scripts/predictive_subspace_check.py merge --profile plasmode --out-dir /tmp/psc
    python scripts/predictive_subspace_check.py decide --out-dir /tmp/psc

    # one SLURM array task per shard (scripts/predictive_subspace_array.sbatch), then merge and decide
"""

from __future__ import annotations

import argparse
import dataclasses
import functools
import json
import math
import sys
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path
from time import perf_counter
from typing import Any, Iterable, Mapping, Sequence

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import latent_space_comparison as lsc  # noqa: E402  (designs, cells, seeds, settings files)
import plasmode_check as pc  # noqa: E402  (export contract, case-study analysis, fidelity)
import production_latent_space_check as plsc  # noqa: E402  (plasmode units, span plant, paired tests)
import span_split_check as ssc  # noqa: E402  (the weak-design control units)
import stage_span_paper_grade_decision as grade  # noqa: E402  (the count-plus-ceiling rule)

sys.path.pop(0)

DRIVER_VERSION = 1
RULE_VERSION = 1
PROFILES: tuple[str, ...] = ("simulation", "weak_size", "plasmode")
REFERENCE_ARM = "pls_cv"
#: The candidates, in the order of preference among qualifying arms with equal off-target counts.
CANDIDATE_ARMS: tuple[str, ...] = ("predictive", "predictive_orthonormal", "predictive_fitted")
CONTROL_ARM = "forced_k_minus_1"
RANK_PLUS_ONE_ARM = "predictive_rank_plus_one"
ARMS: tuple[str, ...] = (REFERENCE_ARM, *CANDIDATE_ARMS, CONTROL_ARM, RANK_PLUS_ONE_ARM)
PROFILE_ARMS: dict[str, tuple[str, ...]] = {
    "simulation": ARMS,
    "weak_size": ARMS[:-1],
    "plasmode": ARMS[:-1],
}
#: Paired comparisons written by ``merge``: (arm, reference arm).
PAIRS: tuple[tuple[str, str], ...] = (
    *((arm, REFERENCE_ARM) for arm in (*CANDIDATE_ARMS, CONTROL_ARM, RANK_PLUS_ONE_ARM)),
    (CONTROL_ARM, "predictive"),
    (RANK_PLUS_ONE_ARM, "predictive"),
)
STATISTICS: tuple[str, ...] = ("delta", "angle", "shape")
TARGETS: dict[str, str] = {"magnitude": "delta", "orientation": "angle", "shape": "shape"}
ALPHA = 0.05
REJECTION_RULES: dict[str, str] = {"simulation": "p < alpha", "weak_size": "p < alpha", "plasmode": "p <= alpha"}
CONSISTENCY_TOLERANCE = 1e-8
CONSISTENCY_CHECKS: tuple[str, ...] = ("rank", "fitted_labels", "same_fit", "one_subspace", "identity_at_low_rank")
#: Relative singular-value tolerance of the Y loadings' rank in the ``rank`` check.
RANK_TOLERANCE = 1e-10

REPO = Path(__file__).resolve().parents[1]
COMPARISON_DIR = REPO / "results" / "latent-space-comparison-2026-11-04"
SPAN_SPLIT_DIR = REPO / "results" / "span-split-2026-11-05" / "simulation"
PLASMODE_DIR = REPO / "results" / "plasmode-2026-11-02"
PRODUCTION_PLASMODE_DIR = REPO / "results" / "production-latent-space-2026-11-10" / "plasmode"

# Simulation profile (design D6).
BASE_SEED = 2036
SIMULATION_REPLICATES = 100
SIMULATION_PERMUTATIONS = 199
SPAN_SCALE = "planted"
STRONG_DESIGNS: tuple[str, ...] = ("n200", "seaad")
BAND_DESIGNS: tuple[str, ...] = ("n200_delta1", "n200_delta0.5", "n200_delta0.25", "seaad_delta0.5")
WEAK_DESIGN = "seaad_weak"
SURGERY_DESIGN = "seaad_surgery"
#: The designs read by S1 (every stage-span design).
STAGE_SPAN_DESIGNS: tuple[str, ...] = (*STRONG_DESIGNS, *BAND_DESIGNS, WEAK_DESIGN)
SIMULATION_DESIGNS: tuple[str, ...] = (*STAGE_SPAN_DESIGNS, SURGERY_DESIGN)
#: design -> (the comparison design it is built from, its stage effect δ; ``None`` is the committed weak δ).
DESIGN_SOURCES: dict[str, tuple[str, float | None]] = {
    "n200": ("n200", 2.0),
    "seaad": ("seaad", 2.0),
    "n200_delta1": ("n200", 1.0),
    "n200_delta0.5": ("n200", 0.5),
    "n200_delta0.25": ("n200", 0.25),
    "seaad_delta0.5": ("seaad", 0.5),
    "seaad_weak": ("seaad_weak", None),
    "seaad_surgery": ("seaad", 2.0),
}
GENERATOR_DELTA = 2.0
STAGE_EFFECT_KEYS: tuple[str, ...] = ("delta_methyl", "delta_expr", "delta_protein")
#: (cell, trajectory mode, effect, φ or None, support or None), the comparison's cell layout.
CELLS: tuple[tuple[str, str, float, float | None, int | None], ...] = (
    ("anchor", "none", 0.0, None, None),
    ("translation", "translation", 1.0, None, None),
    *((f"magnitude_e{effect:g}", "magnitude", effect, None, None) for effect in (0.05, 0.1, 0.25)),
    *((f"orient_e{effect:g}_phi0", "orientation", effect, 0.0, None) for effect in (0.1, 0.25, 0.5)),
    *((f"shape_e{effect:g}", "shape", effect, None, None) for effect in (0.1, 0.25, 0.5)),
)
SURGERY_CELL: tuple[str, str, float, float | None, int | None] = ("surgery_orient_e1", "orientation", 1.0, None, None)
CONTROL_CELLS: tuple[str, ...] = ("anchor", "translation")
TARGET_CELLS: tuple[str, ...] = tuple(cell[0] for cell in CELLS if cell[1] in TARGETS)

# Weak-size profile.
WEAK_REPLICATES = 500
WEAK_PERMUTATIONS = 999
#: The span-split run's unit settings the weak profile must share.
_SHARED_UNIT_KEYS: tuple[str, ...] = (
    "cells", "base_seed", "pls_cv_params", "min_components", "orientation_estimator", "weak_delta",
)

# Plasmode profile.
PLASMODE_DESIGN = plsc.PLASMODE_DESIGN
ANCHOR_CELL = plsc.ANCHOR_CELL
PLASMODE_ANCHOR_REPLICATES = 500
PLASMODE_PLANT_REPLICATES = 100
PLASMODE_PERMUTATIONS = 999
PLANT_CELLS: tuple[str, ...] = tuple(str(name) for name, _ in plsc.plasmode_cells() if name != ANCHOR_CELL)
#: The production selection's arm the plasmode's ``pls_cv`` is (no weighting, no sparsity).
_PRODUCTION_ARM_SETTINGS = {"block_weighting": None, "sparsity": None}
#: How the plasmode's predictive arms are fitted (the amendment of 2026-11-12, before any array unit).
PLASMODE_PREDICTIVE_FIT = "forced_components = the unit's pls_cv rank (the same final fit; no repeated CV)"

# Decision (design D7).
LOSS_ALPHA = plsc.LOSS_ALPHA
FAMILY_ALPHA = plsc.FAMILY_ALPHA
VERDICTS: tuple[str, ...] = ("reproduction_failed", "propose_paper_grade", "keep_production")
ARM_VERDICTS: tuple[str, ...] = ("size_excess", "power_loss", "specificity_excess", "no_gain", "qualifies")
_VERDICT_OF = {"S1": "size_excess", "S2": "power_loss", "S3": "specificity_excess", "G": "no_gain"}
CONSEQUENCES: dict[str, str] = {
    "reproduction_failed": (
        "Nothing is read: the production arm does not reproduce a committed parent, the plasmode's label-free "
        "latent space, or a unit's predictive arms do not measure its production fit. Find the drift before any "
        "rerun."
    ),
    "propose_paper_grade": (
        "A separate change runs the named arm at 500 x 999 on the stage-span paper grade's columns and the band; "
        "only its verdict may flip the harness and motco analyze defaults, with pins for every committed "
        "configuration, the CLI flag, attribution and the split test on the measured subspace, and the SEA-AD "
        "re-issue."
    ),
    "keep_production": (
        "Production keeps measuring in the fitted span's scores; the three predictive coordinate systems stay "
        "opt-in diagnostics, and the measurement-visibility page records what the measured dimension costs and "
        "why it was kept."
    ),
}
DECISION_NAME = "predictive_subspace_decision.json"


# ── Shared helpers ───────────────────────────────────────────────────────────


def _parent(path: Path) -> dict[str, str]:
    return plsc._parent(Path(path))


def profile_dir(out_dir: Path, profile: str) -> Path:
    return Path(out_dir) / profile


def _rejects(p_value: Any, alpha: float, rule: str) -> bool | None:
    return plsc._rejects(p_value, alpha, rule)


def _quantile(values: Iterable[Any], q: float) -> float | None:
    return plsc._quantile(values, q)


def _number(value: Any) -> float | None:
    return plsc._number(value)


def _canonical(value: Any) -> Any:
    return json.loads(json.dumps(value, sort_keys=True))


@functools.lru_cache(maxsize=1)
def _reference() -> Any:
    from motco.simulations.reference import load_reference

    return load_reference()


def unit_key(design: str, cell: str, replicate: int) -> str:
    return json.dumps([str(design), str(cell), int(replicate)])


def record_key(record: Mapping[str, Any]) -> str:
    return unit_key(record["design"], record["cell"], record["replicate"])


# ── Simulation profile: settings and units ───────────────────────────────────


def design_generator(name: str, comparison_designs: Mapping[str, Any]) -> dict[str, Any]:
    """The design's generator: its comparison design with the stage effect, the planted scale or the surgery."""

    source, delta = DESIGN_SOURCES[name]
    generator = dict(comparison_designs[source]["generator"])
    if delta is not None and float(delta) != GENERATOR_DELTA:
        generator.update({key: float(delta) for key in STAGE_EFFECT_KEYS})
    if name == SURGERY_DESIGN:
        # The default feature-set surgery: neither the stage-span construction nor its blocks.
        generator.pop("construction", None)
        generator.pop("span_layers", None)
    else:
        generator["span_scale"] = SPAN_SCALE
    return generator


def design_cells(design: str) -> list[str]:
    return [SURGERY_CELL[0]] if design == SURGERY_DESIGN else [cell[0] for cell in CELLS]


def simulation_settings(
    *,
    parent_dir: Path = COMPARISON_DIR,
    replicates: int = SIMULATION_REPLICATES,
    permutations: int = SIMULATION_PERMUTATIONS,
    designs: Sequence[str] = SIMULATION_DESIGNS,
    base_seed: int = BASE_SEED,
) -> dict[str, Any]:
    """Seven stage-span designs and the surgery design on the comparison's generators, at a fresh base seed."""

    unknown = sorted(set(designs) - set(SIMULATION_DESIGNS))
    if unknown:
        raise SystemExit(f"Unknown design(s) {unknown}; expected a subset of {list(SIMULATION_DESIGNS)}.")
    parent_dir = Path(parent_dir)
    parent = plsc.load_parent(parent_dir)
    parent_settings = parent["settings"]
    calibration = parent_dir / "calibration" / "calibration.json"
    if not calibration.exists():
        raise SystemExit(f"{calibration} does not exist; the comparison's calibration is required.")
    digest = plsc.sha256_file(calibration)
    expected = (parent_settings.get("weak_delta") or {}).get("calibration_sha256")
    if digest != expected:
        raise SystemExit(
            f"{calibration} has sha256 {digest}, not the comparison's calibration_sha256 {expected}; refusing a "
            "weak δ the comparison did not use."
        )
    comparison = lsc.build_settings(calibration=calibration, replicates=replicates, permutations=permutations)
    shared = {
        "designs": parent_settings["designs"],
        **{key: parent_settings[key] for key in ("pls_cv_params", "min_components", "orientation_estimator",
                                                 "weak_delta")},
    }
    differing = [key for key, value in shared.items() if _canonical(comparison[key]) != _canonical(value)]
    if differing:
        raise SystemExit(
            f"The comparison's generators differ from {parent_dir}/settings.json in: {', '.join(differing)}; the "
            "designs would not be the comparison's."
        )
    ordered = [name for name in SIMULATION_DESIGNS if name in designs]
    units = {
        "designs": {
            name: {"generator": design_generator(name, comparison["designs"]), "layers": list(lsc.LAYERS)}
            for name in ordered
        },
        "cells": [list(cell) for cell in CELLS] + [list(SURGERY_CELL)],
        "arms": [lsc.REFERENCE_ARM],
        "pls_cv_params": comparison["pls_cv_params"],
        "min_components": comparison["min_components"],
        "orientation_estimator": comparison["orientation_estimator"],
        "weak_delta": comparison["weak_delta"],
        "replicates": int(replicates),
        "permutations": int(permutations),
        "alpha": ALPHA,
        "base_seed": int(base_seed),
    }
    weak_value = float(comparison["weak_delta"]["value"])
    stage_effects = {}
    for name in ordered:
        delta = DESIGN_SOURCES[name][1]
        stage_effects[name] = weak_value if delta is None else float(delta)
    return {
        "driver_version": DRIVER_VERSION,
        "profile": "simulation",
        "parents": {
            "comparison_settings": _parent(parent_dir / "settings.json"),
            "comparison_calibration": _parent(calibration),
        },
        "units": units,
        "design_cells": {name: design_cells(name) for name in ordered},
        "stage_effects": stage_effects,
        "span_scale": SPAN_SCALE,
        "arms": list(PROFILE_ARMS["simulation"]),
        "replicates": int(replicates),
        "permutations": int(permutations),
        "alpha": ALPHA,
        "rejection": REJECTION_RULES["simulation"],
        "consistency_tolerance": CONSISTENCY_TOLERANCE,
    }


def simulation_units(settings: Mapping[str, Any]) -> list[lsc.Unit]:
    return [
        lsc.Unit(design, cell, replicate)
        for design, cells in settings["design_cells"].items()
        for cell in cells
        for replicate in range(int(settings["replicates"]))
    ]


# ── Weak-size profile: settings and units ────────────────────────────────────


def weak_size_settings(
    *,
    replicates: int = WEAK_REPLICATES,
    permutations: int = WEAK_PERMUTATIONS,
    parent_dir: Path = SPAN_SPLIT_DIR,
    comparison_dir: Path = COMPARISON_DIR,
) -> dict[str, Any]:
    """The span-split check's weak-design control units, checked against its committed settings."""

    units = ssc.simulation_settings(
        parent_dir=comparison_dir,
        replicates=0,
        control_replicates=replicates,
        permutations=permutations,
        designs=[WEAK_DESIGN],
    )
    split = ssc.load_parent(Path(parent_dir))
    theirs, ours = split["settings"]["comparison"], units["comparison"]
    differing = [key for key in _SHARED_UNIT_KEYS if _canonical(ours[key]) != _canonical(theirs[key])]
    if _canonical(ours["designs"][WEAK_DESIGN]) != _canonical(theirs["designs"].get(WEAK_DESIGN)):
        differing.append(f"designs.{WEAK_DESIGN}")
    if list(units["control_cells"]) != list(CONTROL_CELLS):
        differing.append("control_cells")
    if differing:
        raise SystemExit(
            f"The weak design's units differ from {parent_dir}/settings.json in: {', '.join(differing)}; they would "
            "not be the span-split check's units."
        )
    return {
        "driver_version": DRIVER_VERSION,
        "profile": "weak_size",
        "parents": {
            "span_split_settings": _parent(Path(parent_dir) / "settings.json"),
            "span_split_production": _parent(Path(parent_dir) / "production_table.csv"),
            "comparison_calibration": _parent(Path(comparison_dir) / "calibration" / "calibration.json"),
        },
        "units": units,
        "design": WEAK_DESIGN,
        "cells": list(CONTROL_CELLS),
        "arms": list(PROFILE_ARMS["weak_size"]),
        "replicates": int(replicates),
        "permutations": int(permutations),
        "alpha": ALPHA,
        "rejection": REJECTION_RULES["weak_size"],
        "consistency_tolerance": CONSISTENCY_TOLERANCE,
    }


def weak_size_units(settings: Mapping[str, Any]) -> list[lsc.Unit]:
    return ssc.simulation_units(settings["units"])


def _harness_units(settings: Mapping[str, Any]) -> dict[str, Any]:
    """The comparison-shaped settings ``lsc.unit_dataset`` / ``cell_spec`` / ``arm_params`` read."""

    return settings["units"] if settings["profile"] == "simulation" else settings["units"]["comparison"]


# ── Plasmode profile: settings and units ─────────────────────────────────────


def plasmode_settings(
    export: Any,
    *,
    parent_dir: Path = PLASMODE_DIR,
    production_dir: Path = PRODUCTION_PLASMODE_DIR,
    replicates: int = PLASMODE_ANCHOR_REPLICATES,
    plant_replicates: int = PLASMODE_PLANT_REPLICATES,
    permutations: int = PLASMODE_PERMUTATIONS,
) -> dict[str, Any]:
    """The plasmode's configuration and the production selection's plant, refused on another export."""

    parent_dir, production_dir = Path(parent_dir), Path(production_dir)
    parent = plsc.load_parent(parent_dir)
    parent_settings = parent["settings"]
    plsc.verify_export(export, {"export": parent_settings["export"]})
    production = plsc.load_parent(production_dir)["settings"]
    cells = plsc.plasmode_cells()
    differing = []
    if _canonical(production.get("plasmode")) != _canonical(parent_settings):
        differing.append("plasmode")
    if _canonical(production.get("cells")) != _canonical(cells):
        differing.append("cells")
    if _canonical(production.get("plant")) != _canonical(_plant_description()):
        differing.append("plant")
    if differing:
        raise SystemExit(
            f"The production selection's plasmode settings ({production_dir}/settings.json) differ from this "
            f"profile's in: {', '.join(differing)}; its committed plant rows would not be these units."
        )
    summary = parent_dir / "plasmode_summary.csv"
    return {
        "driver_version": DRIVER_VERSION,
        "profile": "plasmode",
        "parents": {
            "plasmode_settings": _parent(parent_dir / "settings.json"),
            "plasmode_summary": _parent(summary),
            "production_plasmode_settings": _parent(production_dir / "settings.json"),
            "production_plasmode_rates": _parent(production_dir / "rate_table.csv"),
        },
        "export": {"files": dict(export.files), "cohort": export.manifest.get("cohort")},
        "plasmode": parent_settings,
        "design": PLASMODE_DESIGN,
        "cells": cells,
        "replicates": {"anchor": int(replicates), "plant": int(plant_replicates)},
        "permutations": int(permutations),
        "fidelity": dict(parent_settings["fidelity"]),
        "attribution": False,
        "plant": _plant_description(),
        "predictive_arm_fit": PLASMODE_PREDICTIVE_FIT,
        "arms": list(PROFILE_ARMS["plasmode"]),
        "alpha": ALPHA,
        "rejection": REJECTION_RULES["plasmode"],
        "consistency_tolerance": CONSISTENCY_TOLERANCE,
    }


def _plant_description() -> dict[str, Any]:
    """The production selection's description of its span plant (``plsc.plasmode_settings``' ``plant``)."""

    return {
        "group": pc.GROUPS[1],
        "shift": "e * (pooled stage-s mean - mean of the pooled stage means), standardized, mapped back",
        "scaler": "the export's pooled per-feature mean and SD (ddof 0), the analysis's sample scaler",
        "geometry": "the planted group's centred stage configuration gains e times the pooled one",
        "labels": "the plasmode's within-stage permutation of its paper stream, never the real labels",
    }


def plasmode_units(settings: Mapping[str, Any]) -> list[lsc.Unit]:
    return plsc.plasmode_units(settings)


def plasmode_production_params(settings: Mapping[str, Any], seeds: Mapping[str, int]) -> Any:
    """The production selection's ``pls_cv`` analysis parameters for the unit's seeds."""

    view = {
        "plasmode": settings["plasmode"],
        "permutations": settings["permutations"],
        "arms": {plsc.REFERENCE_ARM: dict(_PRODUCTION_ARM_SETTINGS)},
    }
    return plsc.plasmode_arm_params(plsc.REFERENCE_ARM, view, seeds)


def plasmode_arm_params(arm: str, production: Any, *, n_stages: int, production_rank: int | None = None) -> Any:
    """The case-study parameters plus the arm's setting.

    A predictive arm runs at ``forced_components`` = ``production_rank`` (the
    unit's ``pls_cv`` rank): its final fit is the CV's final fit, so the CV is
    not repeated.
    """

    if arm == REFERENCE_ARM:
        return production
    if arm in CANDIDATE_ARMS:
        if production_rank is None:
            raise ValueError(f"The plasmode arm {arm!r} needs the unit's pls_cv rank.")
        return dataclasses.replace(production, latent_coordinates=arm, forced_components=int(production_rank))
    if arm == CONTROL_ARM:
        return dataclasses.replace(production, forced_components=int(n_stages) - 1)
    raise ValueError(f"The plasmode profile has no arm {arm!r}.")


def enumerate_units(settings: Mapping[str, Any]) -> list[lsc.Unit]:
    profile = settings["profile"]
    if profile == "simulation":
        return simulation_units(settings)
    if profile == "weak_size":
        return weak_size_units(settings)
    return plasmode_units(settings)


# ── One unit: arms, records and the consistency checks ───────────────────────


def _forced(params: Any, rank: int) -> Any:
    return dataclasses.replace(params, integration_params={**params.integration_params, "forced_components": int(rank)})


def harness_arm_record(result: Any, seconds: float) -> dict[str, Any]:
    """One arm's rank, dimension, AUROC, statistics, p-values, null width, capture and population statistics."""

    from motco.simulations.diagnostics import POPULATION_LATENT_CHECKPOINT

    metadata = result.latent_matrix_metadata
    geometry = result.realized_geometry
    return {
        "status": "completed",
        "rank": int(metadata["selected_lv"]),
        "component_selection": metadata.get("component_selection"),
        "measured_dimension": int(metadata.get("measured_dimension", metadata["n_features"])),
        "cv_mean_auroc": metadata.get("cv_mean_auroc"),
        "pair_statistics": lsc._finite(result.pair_statistics),
        "p_values": lsc._finite(result.p_values),
        "angle_null_q95": (result.null_summary.get("angle") or {}).get("q95"),
        "span_capture": geometry.get("span_capture"),
        "population_latent": (geometry.get("checkpoints", {}).get(POPULATION_LATENT_CHECKPOINT) or {}).get("joint"),
        "seconds": seconds,
    }


def analysis_arm_record(result: Any, seconds: float) -> dict[str, Any]:
    """One plasmode arm's record from the analysis result: statistics only, never a donor value."""

    record = result.record
    latent = record["latent"]
    coordinates = record["measurement"].get("latent_coordinates") or {}
    return {
        "status": "completed",
        "rank": int(latent["selected_lv"]),
        "component_selection": latent.get("component_selection"),
        "measured_dimension": int(coordinates.get("measured_dimension", latent["selected_lv"])),
        "cv_mean_auroc": latent.get("cv_mean_auroc"),
        "pair_statistics": lsc._finite(record["statistics"]),
        "p_values": lsc._finite(record["p_values"]),
        "angle_null_q95": (record.get("null_summary", {}).get("angle") or {}).get("q95"),
        "span_capture": None,
        "population_latent": None,
        "seconds": seconds,
    }


def _failed_arm(exc: Exception, seconds: float) -> dict[str, Any]:
    return {"status": "failed", "error_type": type(exc).__name__, "error_message": str(exc), "seconds": seconds}


def loadings_rank(model: Any, *, tolerance: float = RANK_TOLERANCE) -> int:
    """The numerical rank of the model's Y loadings (singular values above ``tolerance`` of the largest)."""

    singular = np.linalg.svd(np.asarray(model.y_loadings_, dtype=float), compute_uv=False)
    return int(np.sum(singular > tolerance * singular[0])) if singular.size and singular[0] > 0 else 0


def stage_r2(model: Any, features: np.ndarray, labels: Sequence[str]) -> float:
    """The in-sample share of the one-hot stage labels the fit explains (categories sorted as strings)."""

    names = np.asarray([str(label) for label in labels])
    categories = sorted(set(names.tolist()))
    onehot = (names[:, None] == np.asarray(categories)[None, :]).astype(float)
    fitted = np.asarray(model.predict(np.asarray(features, dtype=float)), dtype=float)
    centred = onehot - onehot.mean(axis=0)
    return float(1.0 - np.sum((onehot - fitted) ** 2) / np.sum(centred**2))


def _check(deviation: float | None, tolerance: float, *, applicable: bool = True) -> dict[str, Any]:
    if not applicable:
        return {"applicable": False, "passed": True, "deviation": None, "tolerance": tolerance}
    passed = deviation is not None and math.isfinite(float(deviation)) and float(deviation) <= tolerance
    return {"applicable": True, "passed": bool(passed), "deviation": deviation, "tolerance": tolerance}


def _share_deviation(captures: Sequence[Mapping[str, Any] | None]) -> float | None:
    """The largest difference of a span-capture share across arms (``None`` shares must agree as ``None``)."""

    if any(capture is None for capture in captures):
        return None
    worst = 0.0
    for key in ("group_difference", "pooled_stage", "rank"):
        values = [capture.get(key) for capture in captures if capture is not None]
        if all(value is None for value in values):
            continue
        if any(value is None for value in values):
            return math.inf
        worst = max(worst, float(max(values)) - float(min(values)))
    return worst


def _p_deviation(left: Mapping[str, Any], right: Mapping[str, Any]) -> float:
    worst = 0.0
    for statistic in STATISTICS:
        a, b = left.get(statistic), right.get(statistic)
        if a is None and b is None:
            continue
        if a is None or b is None:
            return math.inf
        worst = max(worst, abs(float(a) - float(b)))
    return worst


def consistency_checks(
    *,
    model: Any,
    scores: np.ndarray,
    n_stages: int,
    candidates: Mapping[str, np.ndarray],
    arms: Mapping[str, Mapping[str, Any]],
    with_capture: bool,
    models: Mapping[str, Any] | None = None,
    tolerance: float = CONSISTENCY_TOLERANCE,
) -> dict[str, dict[str, Any]]:
    """The five checks that the predictive arms measure the production fit; each with its deviation.

    ``candidates`` holds each completed predictive arm's latent matrix, ``models``
    its fitted model (default: ``model``, the production fit itself) and
    ``arms`` every arm's record.
    """

    rank = int(np.asarray(model.x_rotations_).shape[1])
    expected = min(rank, int(n_stages) - 1)
    fitted = np.asarray(scores, dtype=float) @ np.asarray(model.y_loadings_, dtype=float).T
    label_deviation = 0.0
    for matrix in candidates.values():
        matrix = np.asarray(matrix, dtype=float)
        mapping, *_ = np.linalg.lstsq(matrix, fitted, rcond=None)
        label_deviation = max(label_deviation, float(np.max(np.abs(matrix @ mapping - fitted))))
    production = arms.get(REFERENCE_ARM) or {}
    completed = {arm: arms[arm] for arm in CANDIDATE_ARMS if (arms.get(arm) or {}).get("status") == "completed"}
    fit_deviation = 0.0
    for arm, record in completed.items():
        fitted_model = (models or {}).get(arm, model)
        fit_deviation = max(fit_deviation, _model_deviation(model, fitted_model))
        if record.get("rank") != production.get("rank"):
            fit_deviation = math.inf
    captures = [record.get("span_capture") for record in completed.values()]
    low_rank = rank <= int(n_stages) - 1 and "predictive" in completed and production.get("status") == "completed"
    return {
        "rank": _check(float(abs(loadings_rank(model) - expected)), 0.0),
        "fitted_labels": _check(label_deviation if candidates else None, tolerance, applicable=bool(candidates)),
        "same_fit": _check(fit_deviation if completed else None, tolerance, applicable=bool(completed)),
        "one_subspace": _check(
            _share_deviation(captures) if len(captures) > 1 else None,
            tolerance,
            applicable=with_capture and len(captures) > 1,
        ),
        "identity_at_low_rank": _check(
            _p_deviation(completed["predictive"]["p_values"], production["p_values"]) if low_rank else None,
            0.0,
            applicable=low_rank,
        ),
    }


def _model_deviation(reference: Any, other: Any) -> float:
    """The largest difference of the rotations and Y loadings of two fitted models (inf on another shape)."""

    if other is reference:
        return 0.0
    worst = 0.0
    for name in ("x_rotations_", "y_loadings_"):
        left, right = np.asarray(getattr(reference, name), dtype=float), np.asarray(getattr(other, name), dtype=float)
        if left.shape != right.shape:
            return math.inf
        worst = max(worst, float(np.max(np.abs(left - right))) if left.size else 0.0)
    return worst


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
        "arms": arms,
        "runtime_seconds": timing,
    }


def _failed_unit(base: Mapping[str, Any], exc: Exception, start: float, timing: dict[str, float]) -> dict[str, Any]:
    timing["unit"] = perf_counter() - start
    return {**base, "status": "failed", "error_type": type(exc).__name__, "error_message": str(exc),
            "runtime_seconds": timing}


def evaluate_harness_unit(unit: lsc.Unit, settings: Mapping[str, Any]) -> dict[str, Any]:
    """One dataset, one production fit, every arm at the unit's RRPP seed; one record (simulation, weak_size)."""

    from motco.simulations.evaluation import (
        evaluate_semisynthetic_trajectory,
        integrate_semisynthetic_dataset,
        pls_coordinates_from_pls,
    )

    start = perf_counter()
    units = _harness_units(settings)
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
        dataset = lsc.unit_dataset(unit, units, reference=_reference())
        timing["generate"] = perf_counter() - mark
        params = lsc.arm_params(lsc.REFERENCE_ARM, units, seed=seed, n_stages=n_stages)
        mark = perf_counter()
        scores = integrate_semisynthetic_dataset(dataset, params)
        timing["fit"] = perf_counter() - mark
        model = scores.artifacts.model
        production = {
            "rank": int(scores.metadata["selected_lv"]),
            "cv_mean_auroc": scores.metadata.get("cv_mean_auroc"),
            "loadings_rank": loadings_rank(model),
            "stage_r2": stage_r2(
                model, scores.artifacts.features.to_numpy(dtype=float), dataset.metadata[params.stage_col].tolist()
            ),
        }
        mark = perf_counter()
        arms: dict[str, Any] = {REFERENCE_ARM: harness_arm_record(
            evaluate_semisynthetic_trajectory(dataset, params, latent=scores), 0.0
        )}
        arms[REFERENCE_ARM]["seconds"] = perf_counter() - mark
    except Exception as exc:  # recorded, never aborts the shard
        return _failed_unit(base, exc, start, timing)

    candidates: dict[str, np.ndarray] = {}
    for arm in settings["arms"]:
        if arm == REFERENCE_ARM:
            continue
        mark = perf_counter()
        try:
            if arm in CANDIDATE_ARMS:
                arm_params = dataclasses.replace(params, latent_coordinates=arm)
                latent = pls_coordinates_from_pls(dataset, scores, arm_params)
                candidates[arm] = latent.matrix.to_numpy(dtype=float)
                result = evaluate_semisynthetic_trajectory(dataset, arm_params, latent=latent)
            elif arm == CONTROL_ARM:
                result = evaluate_semisynthetic_trajectory(dataset, _forced(params, n_stages - 1))
            else:
                plus_one = dataclasses.replace(_forced(params, production["rank"] + 1), latent_coordinates="predictive")
                result = evaluate_semisynthetic_trajectory(dataset, plus_one)
            arms[arm] = harness_arm_record(result, perf_counter() - mark)
        except Exception as exc:  # recorded on the arm, the unit goes on
            arms[arm] = _failed_arm(exc, perf_counter() - mark)
    checks = consistency_checks(
        model=model,
        scores=scores.matrix.to_numpy(dtype=float),
        n_stages=n_stages,
        candidates=candidates,
        arms=arms,
        with_capture=bool(params.span_capture),
    )
    return _finish(base, arms, checks, production, start, timing)


def evaluate_plasmode_unit(unit: lsc.Unit, settings: Mapping[str, Any], export: Any) -> dict[str, Any]:
    """Build the unit once and analyse it under every arm serially; one record of seeds and statistics."""

    start = perf_counter()
    parent = settings["plasmode"]
    seeds = pc.unit_seeds(int(parent["base_seed"]), str(parent["stream"]), unit.replicate)
    effect = plsc.plasmode_effect(settings, unit.cell)
    n_stages = len(export.stage_order)
    base = {
        "profile": "plasmode",
        "design": unit.design,
        "cell": unit.cell,
        "trajectory_mode": "none" if unit.cell == ANCHOR_CELL else plsc.PLANT_MODE,
        "effect_size": effect,
        "replicate": unit.replicate,
        "stream": parent["stream"],
        "seeds": seeds,
        "rrpp_seed": seeds["rrpp"],
        "settings_hash": lsc.settings_hash(settings),
    }
    timing: dict[str, float] = {}
    try:
        mark = perf_counter()
        built = plsc.build_plasmode_unit(unit, settings, export)
        base["label_redraws"] = built["label_redraws"]
        timing["generate"] = perf_counter() - mark
        production_params = plasmode_production_params(settings, built["seeds"])
        mark = perf_counter()
        _, result = pc.analyse(built["blocks"], built["metadata"], export, production_params)
        arms: dict[str, Any] = {REFERENCE_ARM: analysis_arm_record(result, perf_counter() - mark)}
        artifacts = result.latent_artifacts
        model = artifacts.model
        scores = result.latent_scores.to_numpy(dtype=float)
        production = {
            "rank": int(result.record["latent"]["selected_lv"]),
            "cv_mean_auroc": result.record["latent"].get("cv_mean_auroc"),
            "loadings_rank": loadings_rank(model),
            "stage_r2": stage_r2(
                model, artifacts.features.to_numpy(dtype=float), built["metadata"][export.stage_col].tolist()
            ),
        }
        if unit.cell == ANCHOR_CELL:
            arms[REFERENCE_ARM]["fidelity"] = plsc.fidelity_outcome(result, settings["fidelity"])
    except Exception as exc:  # recorded, never aborts the shard
        return lsc._jsonable(_failed_unit(base, exc, start, timing))

    candidates: dict[str, np.ndarray] = {}
    models: dict[str, Any] = {}
    for arm in settings["arms"]:
        if arm == REFERENCE_ARM:
            continue
        mark = perf_counter()
        try:
            params = plasmode_arm_params(arm, production_params, n_stages=n_stages, production_rank=production["rank"])
            _, arm_result = pc.analyse(built["blocks"], built["metadata"], export, params)
            record = analysis_arm_record(arm_result, perf_counter() - mark)
            if arm in CANDIDATE_ARMS:
                candidates[arm] = arm_result.latent_scores.to_numpy(dtype=float)
                models[arm] = arm_result.latent_artifacts.model
                # Fitted at the pls_cv rank: the CV that chose it is pls_cv's.
                record["fit"] = "forced_at_pls_cv_rank"
                record["cv_mean_auroc"] = production["cv_mean_auroc"]
            arms[arm] = record
        except Exception as exc:  # recorded on the arm, the unit goes on
            arms[arm] = _failed_arm(exc, perf_counter() - mark)
    checks = consistency_checks(
        model=model, scores=scores, n_stages=n_stages, candidates=candidates, arms=arms, with_capture=False,
        models=models,
    )
    return lsc._jsonable(_finish(base, arms, checks, production, start, timing))


# ── Shards, resume and the run ───────────────────────────────────────────────


def read_records(paths: Iterable[Path]) -> list[dict[str, Any]]:
    """Every parsable unit line; a line torn by an interrupted write is dropped, so its unit is re-run."""

    records: list[dict[str, Any]] = []
    for path in paths:
        with Path(path).open(encoding="utf-8") as handle:
            for line in handle:
                try:
                    record = json.loads(line)
                except json.JSONDecodeError:
                    continue
                if record.get("status") in ("completed", "inconsistent", "failed"):
                    records.append(record)
    return records


_WORKER_EXPORT: Any = None


def _init_plasmode_worker(directory: str, files: Mapping[str, str]) -> None:
    global _WORKER_EXPORT
    export = pc.load_export(Path(directory))
    if export.files != dict(files):
        raise SystemExit(f"{directory} changed after the run started (file hashes differ).")
    _WORKER_EXPORT = export


def _evaluate_plasmode_pending(unit: lsc.Unit, settings: Mapping[str, Any]) -> dict[str, Any]:
    return evaluate_plasmode_unit(unit, settings, _WORKER_EXPORT)


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
        plsc.verify_parents(lsc.load_settings(directory))
    plsc.verify_parents(settings)
    lsc.ensure_settings(directory, settings)
    if profile == "plasmode":
        if export is None:
            raise SystemExit("The plasmode profile needs the export.")
        plsc.verify_export(export, settings)
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

    if profile in ("simulation", "weak_size"):
        if workers <= 1:
            for position, unit in enumerate(pending, start=1):
                write(position, unit, evaluate_harness_unit(unit, settings))
        else:
            with ProcessPoolExecutor(max_workers=workers) as pool:
                futures = {pool.submit(evaluate_harness_unit, unit, settings): unit for unit in pending}
                for position, future in enumerate(as_completed(futures), start=1):
                    write(position, futures[future], future.result())
    elif workers <= 1:
        for position, unit in enumerate(pending, start=1):
            write(position, unit, evaluate_plasmode_unit(unit, settings, export))
    else:
        with ProcessPoolExecutor(
            max_workers=workers, initializer=_init_plasmode_worker, initargs=(str(export.directory), dict(export.files))
        ) as pool:
            futures = {pool.submit(_evaluate_plasmode_pending, unit, settings): unit for unit in pending}
            for position, future in enumerate(as_completed(futures), start=1):
                write(position, futures[future], future.result())
    return len(pending)


# ── Merge: rate, paired and cell tables, accounting, reproduction ────────────


RATE_COLUMNS = (
    "design", "cell", "trajectory_mode", "effect_size", "arm", "statistic", "rate", "se", "count", "rejections",
    "failures", "inconsistent",
)
PAIRED_COLUMNS = (
    "design", "cell", "trajectory_mode", "effect_size", "statistic", "arm", "reference_arm", "n_paired",
    "rate_reference", "rate_arm", "difference", "b", "c", "p_loss", "p_gain", "p_two_sided",
)
_CELL_QUANTITIES: tuple[str, ...] = ("rank", "measured_dimension", "cv_mean_auroc", "pooled_stage")
CELL_COLUMNS = (
    "design", "cell", "trajectory_mode", "effect_size", "arm", "n",
    *(f"{quantity}_{suffix}" for quantity in _CELL_QUANTITIES for suffix in ("median", "q10", "q90")),
    "group_difference_median", "stage_r2_median",
    *(f"observed_{statistic}_median" for statistic in STATISTICS),
    *(f"population_{statistic}_median" for statistic in STATISTICS),
    "angle_null_q95_median",
)


def _cell_fields(record: Mapping[str, Any]) -> dict[str, Any]:
    return {
        "design": str(record["design"]),
        "cell": str(record["cell"]),
        "trajectory_mode": str(record["trajectory_mode"]),
        "effect_size": float(record["effect_size"]),
    }


def _by_cell(records: Iterable[Mapping[str, Any]]) -> dict[tuple[str, str], list[Mapping[str, Any]]]:
    grouped: dict[tuple[str, str], list[Mapping[str, Any]]] = {}
    for record in records:
        grouped.setdefault((str(record["design"]), str(record["cell"])), []).append(record)
    return grouped


def _arm(record: Mapping[str, Any], arm: str) -> Mapping[str, Any] | None:
    """The arm's completed record on a consistent, completed unit, else ``None``."""

    if record.get("status") != "completed":
        return None
    entry = (record.get("arms") or {}).get(arm)
    return entry if entry is not None and entry.get("status") == "completed" else None


def _p(record: Mapping[str, Any], arm: str, statistic: str) -> Any:
    entry = _arm(record, arm)
    return None if entry is None else (entry.get("p_values") or {}).get(statistic)


def rate_table(records: Sequence[Mapping[str, Any]], settings: Mapping[str, Any]) -> list[dict[str, Any]]:
    """Rejections under the profile's rule per design × cell × arm × statistic; inconsistent units excluded."""

    alpha, rule = float(settings["alpha"]), str(settings["rejection"])
    rows = []
    for _, group in sorted(_by_cell(records).items()):
        fields = _cell_fields(group[0])
        inconsistent = sum(1 for record in group if record.get("status") == "inconsistent")
        for arm in settings["arms"]:
            failures = sum(
                1
                for record in group
                if record.get("status") == "failed"
                or (record.get("status") == "completed" and _arm(record, arm) is None)
            )
            for statistic in STATISTICS:
                decisions = [_rejects(_p(record, arm, statistic), alpha, rule) for record in group]
                available = [bool(value) for value in decisions if value is not None]
                count, hits = len(available), sum(available)
                rate = hits / count if count else None
                rows.append({
                    **fields, "arm": arm, "statistic": statistic, "rate": rate,
                    "se": None if rate is None else math.sqrt(rate * (1.0 - rate) / count),
                    "count": count, "rejections": hits, "failures": failures, "inconsistent": inconsistent,
                })
    return rows


def paired_table(records: Sequence[Mapping[str, Any]], settings: Mapping[str, Any]) -> list[dict[str, Any]]:
    """Every pair of :data:`PAIRS` present in the profile, per cell and statistic, over units available in both."""

    alpha, rule = float(settings["alpha"]), str(settings["rejection"])
    arms = set(settings["arms"])
    rows = []
    for _, group in sorted(_by_cell(records).items()):
        fields = _cell_fields(group[0])
        for statistic in STATISTICS:
            for arm, reference in PAIRS:
                if arm not in arms or reference not in arms:
                    continue
                pairs = [
                    (_rejects(_p(record, reference, statistic), alpha, rule),
                     _rejects(_p(record, arm, statistic), alpha, rule))
                    for record in sorted(group, key=lambda r: int(r["replicate"]))
                ]
                pairs = [(ref, cand) for ref, cand in pairs if ref is not None and cand is not None]
                n = len(pairs)
                b = sum(1 for ref, cand in pairs if ref and not cand)
                c = sum(1 for ref, cand in pairs if cand and not ref)
                rows.append({
                    **fields, "statistic": statistic, "arm": arm, "reference_arm": reference, "n_paired": n,
                    "rate_reference": sum(ref for ref, _ in pairs) / n if n else None,
                    "rate_arm": sum(cand for _, cand in pairs) / n if n else None,
                    "difference": (c - b) / n if n else None,
                    "b": b, "c": c, **plsc.paired_p_values(b, c),
                })
    return rows


def cell_table(records: Sequence[Mapping[str, Any]], settings: Mapping[str, Any]) -> list[dict[str, Any]]:
    """Per cell and arm: rank, dimension, AUROC and capture quantiles, stage R², observed and population medians."""

    rows = []
    for _, group in sorted(_by_cell(records).items()):
        fields = _cell_fields(group[0])
        units = [record for record in group if record.get("status") == "completed"]
        for arm in settings["arms"]:
            done = [entry for entry in (_arm(record, arm) for record in units) if entry is not None]
            values = {
                "rank": [entry.get("rank") for entry in done],
                "measured_dimension": [entry.get("measured_dimension") for entry in done],
                "cv_mean_auroc": [entry.get("cv_mean_auroc") for entry in done],
                "pooled_stage": [(entry.get("span_capture") or {}).get("pooled_stage") for entry in done],
            }
            row: dict[str, Any] = {**fields, "arm": arm, "n": len(done)}
            for quantity, entries in values.items():
                for suffix, q in (("median", 0.5), ("q10", 0.1), ("q90", 0.9)):
                    row[f"{quantity}_{suffix}"] = _quantile(entries, q)
            row["group_difference_median"] = _quantile(
                ((entry.get("span_capture") or {}).get("group_difference") for entry in done), 0.5
            )
            row["stage_r2_median"] = _quantile(
                ((record.get("production") or {}).get("stage_r2") for record in units if _arm(record, arm)), 0.5
            )
            for statistic in STATISTICS:
                row[f"observed_{statistic}_median"] = _quantile(
                    ((entry.get("pair_statistics") or {}).get(statistic) for entry in done), 0.5
                )
                row[f"population_{statistic}_median"] = _quantile(
                    ((entry.get("population_latent") or {}).get(statistic) for entry in done), 0.5
                )
            row["angle_null_q95_median"] = _quantile((entry.get("angle_null_q95") for entry in done), 0.5)
            rows.append(row)
    return rows


def accounting(
    records: Sequence[Mapping[str, Any]], settings: Mapping[str, Any], *, expected: int, duplicates: int
) -> dict[str, Any]:
    """Units, failures and inconsistent units by cell, arm failures, duplicates, core-hours and per-arm timing."""

    failures: dict[str, int] = {}
    inconsistent: dict[str, int] = {}
    inconsistent_units: list[dict[str, Any]] = []
    arm_failures: dict[str, dict[str, int]] = {}
    arm_seconds: dict[str, list[float]] = {}
    unit_seconds: list[float] = []
    fit_seconds: list[float] = []
    for record in records:
        label = f"{record['design']}/{record['cell']}"
        timing = record.get("runtime_seconds") or {}
        if timing.get("unit") is not None:
            unit_seconds.append(float(timing["unit"]))
        if timing.get("fit") is not None:
            fit_seconds.append(float(timing["fit"]))
        status = record.get("status")
        if status == "failed":
            failures[label] = failures.get(label, 0) + 1
            continue
        if status == "inconsistent":
            inconsistent[label] = inconsistent.get(label, 0) + 1
            inconsistent_units.append({
                "unit": record_key(record),
                "checks": {
                    name: record["consistency"][name].get("deviation") for name in record["inconsistent_checks"]
                },
            })
        for arm, entry in (record.get("arms") or {}).items():
            if entry.get("status") != "completed":
                by_cell = arm_failures.setdefault(str(arm), {})
                by_cell[label] = by_cell.get(label, 0) + 1
            elif entry.get("seconds") is not None:
                arm_seconds.setdefault(str(arm), []).append(float(entry["seconds"]))
    return {
        "profile": settings["profile"],
        "settings_hash": lsc.settings_hash(settings),
        "expected_units": expected,
        "complete_units": len(records),
        "missing_units": expected - len(records),
        "failed_units": sum(failures.values()),
        "failed_units_by_cell": failures,
        "inconsistent_units": sum(inconsistent.values()),
        "inconsistent_units_by_cell": inconsistent,
        "inconsistent_unit_checks": inconsistent_units,
        "arm_failures_by_cell": arm_failures,
        "duplicate_records": duplicates,
        "recorded_core_hours": sum(unit_seconds) / 3600.0,
        "unit_seconds": {
            "median": _quantile(unit_seconds, 0.5),
            "q90": _quantile(unit_seconds, 0.9),
            "max": max(unit_seconds) if unit_seconds else None,
        },
        "fit_seconds_median": _quantile(fit_seconds, 0.5),
        "arm_seconds_median": {arm: _quantile(arm_seconds.get(arm, []), 0.5) for arm in settings["arms"]},
        "arm_seconds_total_hours": {arm: sum(arm_seconds.get(arm, [])) / 3600.0 for arm in settings["arms"]},
    }


def _compare(
    label: str, keys: Sequence[tuple[str, ...]], committed: Mapping[tuple[str, ...], Mapping[str, int]],
    reproduced: Mapping[tuple[str, ...], Mapping[str, int]], *, source: Mapping[str, Any], rule: str,
) -> dict[str, Any]:
    rows, differing = [], []
    for key in keys:
        theirs, mine = committed.get(key), reproduced.get(key)
        entry = {"key": list(key), "committed": theirs, "reproduced": mine}
        entry["matches"] = theirs is not None and theirs == mine
        rows.append(entry)
        if not entry["matches"]:
            differing.append(entry)
    return {
        "check": label,
        "source": dict(source),
        "rule": rule,
        "rows_compared": len(rows),
        "rows_matching": len(rows) - len(differing),
        "passed": bool(rows) and not differing,
        "rows": rows,
        "differing": differing,
    }


def _ours(rates: Sequence[Mapping[str, Any]]) -> dict[tuple[str, ...], dict[str, int]]:
    return {
        (str(row["design"]), str(row["cell"]), str(row["statistic"])): {
            "count": int(row["count"]), "rejections": int(row["rejections"])
        }
        for row in rates
        if row["arm"] == REFERENCE_ARM
    }


def weak_size_reproduction(
    rates: Sequence[Mapping[str, Any]], settings: Mapping[str, Any], committed_rows: Sequence[Mapping[str, Any]]
) -> dict[str, Any]:
    """``pls_cv`` against the span-split run's committed ``production_table.csv``: exact count and rejections."""

    committed = {
        (str(row["design"]), str(row["cell"]), str(row["statistic"])): {
            "count": int(row["count"]), "rejections": int(row["rejections"])
        }
        for row in committed_rows
    }
    keys = [(WEAK_DESIGN, cell, statistic) for cell in settings["cells"] for statistic in STATISTICS]
    check = _compare(
        "span_split_production", keys, committed, _ours(rates),
        source=settings["parents"]["span_split_production"],
        rule="pls_cv count and rejections (p < alpha) equal the span-split production_table.csv",
    )
    return {"profile": "weak_size", "passed": check["passed"], "checks": {"span_split_production": check}}


def plasmode_summary_rejections(rows: Sequence[Mapping[str, Any]]) -> dict[tuple[str, ...], dict[str, int]]:
    """The plasmode's committed anchor rejections: ``reject_<statistic>`` rows as count n and round(mean·n)."""

    out = {}
    for row in rows:
        if str(row["cell"]) != ANCHOR_CELL or str(row["arm"]) != "analysis" or str(row["subset"]) != "all":
            continue
        quantity = str(row["quantity"])
        if not quantity.startswith("reject_"):
            continue
        statistic = quantity.removeprefix("reject_")
        n = int(float(row["n"]))
        out[(PLASMODE_DESIGN, ANCHOR_CELL, statistic)] = {"count": n, "rejections": int(round(float(row["mean"]) * n))}
    return out


def plasmode_reproduction(
    records: Sequence[Mapping[str, Any]],
    rates: Sequence[Mapping[str, Any]],
    settings: Mapping[str, Any],
    summary_rows: Sequence[Mapping[str, Any]],
    production_rows: Sequence[Mapping[str, Any]],
) -> dict[str, Any]:
    """The anchor against the plasmode, the plant against the production selection, and the fidelity values."""

    ours = _ours(rates)
    anchor = _compare(
        "plasmode_anchor",
        [(PLASMODE_DESIGN, ANCHOR_CELL, statistic) for statistic in STATISTICS],
        plasmode_summary_rejections(summary_rows),
        ours,
        source=settings["parents"]["plasmode_summary"],
        rule="pls_cv anchor count and rejections (p <= alpha) equal plasmode_summary.csv's reject_* rows",
    )
    committed = {
        (str(row["design"]), str(row["cell"]), str(row["statistic"])): {
            "count": int(row["count"]), "rejections": int(row["rejections"])
        }
        for row in production_rows
        if str(row["arm"]) == plsc.REFERENCE_ARM
    }
    plant = _compare(
        "production_plant",
        [(PLASMODE_DESIGN, cell, statistic) for cell in PLANT_CELLS for statistic in STATISTICS],
        committed,
        ours,
        source=settings["parents"]["production_plasmode_rates"],
        rule="pls_cv plant count and rejections (p <= alpha) equal the production selection's rate_table.csv",
    )
    anchors = [record for record in records if record["cell"] == ANCHOR_CELL]
    unavailable, failing = [], []
    for record in anchors:
        fidelity = ((record.get("arms") or {}).get(REFERENCE_ARM) or {}).get("fidelity")
        if record.get("status") == "failed" or fidelity is None:
            unavailable.append(record_key(record))
        elif not fidelity.get("passed"):
            failing.append({"unit": record_key(record), "rank": fidelity.get("rank"),
                            "cv_mean_auroc_rounded": fidelity.get("cv_mean_auroc_rounded")})
    checked = len(anchors) - len(unavailable)
    fidelity_check = {
        "check": "fidelity",
        "source": settings["parents"]["plasmode_settings"],
        "rule": "every anchor unit's pls_cv rank and held-out AUROC (6 decimals) equal the plasmode's fidelity values",
        "expected": dict(settings["fidelity"]),
        "anchor_units": len(anchors),
        "units_checked": checked,
        "units_passed": checked - len(failing),
        "unavailable": unavailable,
        "passed": checked > 0 and not failing and not unavailable,
        "differing": failing,
    }
    checks = {"plasmode_anchor": anchor, "production_plant": plant, "fidelity": fidelity_check}
    return {"profile": "plasmode", "passed": all(check["passed"] for check in checks.values()), "checks": checks}


def merge(out_dir: Path, profile: str) -> dict[str, Any]:
    """Merge one profile's shards and write its tables, accounting and reproduction."""

    directory = profile_dir(out_dir, profile)
    settings = lsc.load_settings(directory)
    plsc.verify_parents(settings)
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
    if profile == "weak_size":
        source = plsc._resolve(settings["parents"]["span_split_production"]["path"])
        reproduction = weak_size_reproduction(rates, settings, lsc._read_csv(source))
    elif profile == "plasmode":
        reproduction = plasmode_reproduction(
            records,
            rates,
            settings,
            lsc._read_csv(plsc._resolve(settings["parents"]["plasmode_summary"]["path"])),
            lsc._read_csv(plsc._resolve(settings["parents"]["production_plasmode_rates"]["path"])),
        )
    else:
        reproduction = {
            "profile": "simulation",
            "passed": True,
            "checks": {},
            "note": "a fresh seed family (base seed 2036): there is no committed run to reproduce",
        }
    lsc._write_json(directory / "accounting.json", summary)
    lsc._write_json(directory / "reproduction.json", reproduction)
    return {
        "profile": profile,
        **{key: summary[key] for key in ("complete_units", "missing_units", "failed_units", "inconsistent_units",
                                         "duplicate_records", "recorded_core_hours")},
        "reproduction_passed": reproduction["passed"],
        "reproduction_differing": {
            name: len(check.get("differing", [])) for name, check in reproduction["checks"].items()
        },
    }


# ── Decision (pure functions) ────────────────────────────────────────────────

#: A reading: (profile, design, cell, statistic).
Reading = tuple[str, str, str, str]


def _cell_mode(cell: str) -> str:
    for name, mode, *_ in (*CELLS, SURGERY_CELL):
        if name == cell:
            return str(mode)
    if cell in PLANT_CELLS:
        return plsc.PLANT_MODE
    return "none"


def criteria_plan() -> dict[str, list[Reading]]:
    """The readings of every criterion, fixed before any unit.

    S1: anchor and translation at the seven stage-span designs (m = 42),
    ``weak_size`` (m = 6) and the plasmode anchor (m = 3). S2: the target
    cells of the six designs at δ ≥ 0.25 (m = 54) and the plasmode plant's
    ``delta`` (m = 3). S3: the target cells' off-target statistics at δ = 2.0
    (m = 36) and in the band (m = 72). G: the target cells at δ = 2.0
    (m = 18); recorded only: the band (36), the weak design (9) and the
    plasmode plant (3).
    """

    def targets(designs: Sequence[str]) -> list[Reading]:
        return [("simulation", design, cell, TARGETS[_cell_mode(cell)]) for design in designs for cell in TARGET_CELLS]

    def off_target(designs: Sequence[str]) -> list[Reading]:
        return [
            ("simulation", design, cell, statistic)
            for design in designs
            for cell in TARGET_CELLS
            for statistic in STATISTICS
            if statistic != TARGETS[_cell_mode(cell)]
        ]

    plant = [("plasmode", PLASMODE_DESIGN, cell, TARGETS[plsc.PLANT_MODE]) for cell in PLANT_CELLS]
    return {
        "S1_simulation": [
            ("simulation", design, cell, statistic)
            for design in STAGE_SPAN_DESIGNS
            for cell in CONTROL_CELLS
            for statistic in STATISTICS
        ],
        "S1_weak_size": [
            ("weak_size", WEAK_DESIGN, cell, statistic) for cell in CONTROL_CELLS for statistic in STATISTICS
        ],
        "S1_plasmode": [("plasmode", PLASMODE_DESIGN, ANCHOR_CELL, statistic) for statistic in STATISTICS],
        "S2": targets((*STRONG_DESIGNS, *BAND_DESIGNS)) + plant,
        "S3_strong": off_target(STRONG_DESIGNS),
        "S3_band": off_target(BAND_DESIGNS),
        "G": targets(STRONG_DESIGNS),
        "G_band": targets(BAND_DESIGNS),
        "G_weak": targets((WEAK_DESIGN,)),
        "G_plasmode": plant,
    }


def _index_rates(tables: Mapping[str, Mapping[str, Sequence[Mapping[str, Any]]]]) -> dict[tuple[str, ...], Any]:
    return {
        (profile, str(row["design"]), str(row["cell"]), str(row["arm"]), str(row["statistic"])): row
        for profile, table in tables.items()
        for row in table["rates"]
    }


def _index_paired(tables: Mapping[str, Mapping[str, Sequence[Mapping[str, Any]]]]) -> dict[tuple[str, ...], Any]:
    return {
        (profile, str(row["design"]), str(row["cell"]), str(row["statistic"]), str(row["arm"]),
         str(row["reference_arm"])): row
        for profile, table in tables.items()
        for row in table["paired"]
    }


def _observation(rates: Mapping[tuple[str, ...], Any], reading: Reading, arm: str) -> dict[str, Any]:
    profile, design, cell, statistic = reading
    entry = {"profile": profile, "design": design, "cell": cell, "arm": arm, "statistic": statistic}
    row = rates.get((profile, design, cell, arm, statistic))
    if row is None or int(row["count"]) <= 0:
        return {**entry, "status": "unavailable"}
    return {**entry, "n_available": int(row["count"]), "n_rejected": int(row["rejections"])}


def _paired(
    paired: Mapping[tuple[str, ...], Any], reading: Reading, arm: str, reference: str = REFERENCE_ARM
) -> dict[str, Any]:
    profile, design, cell, statistic = reading
    entry = {"profile": profile, "design": design, "cell": cell, "arm": arm, "reference_arm": reference,
             "statistic": statistic}
    row = paired.get((profile, design, cell, statistic, arm, reference))
    if row is None or int(row["n_paired"]) <= 0:
        return {**entry, "status": "unavailable"}
    b, c = int(row["b"]), int(row["c"])
    # Recomputed from the discordant counts: the table's p-values are rounded to six significant digits.
    return {
        **entry,
        "n_paired": int(row["n_paired"]),
        "rate_reference": _number(row["rate_reference"]),
        "rate_arm": _number(row["rate_arm"]),
        "difference": _number(row["difference"]),
        "b": b,
        "c": c,
        **plsc.paired_p_values(b, c),
    }


def allowed_count(m: int) -> int:
    """The 95th percentile of Binomial(m, 0.025): the readings a null run has below 0.025."""

    return plsc.loss_allowed(m)


def judge_paired(cells: Sequence[Mapping[str, Any]], *, p_field: str, rule: str) -> dict[str, Any]:
    """The count-plus-ceiling rule on one-sided exact paired p-values (``p_field``); unavailable fails."""

    m = len(cells)
    allowed = allowed_count(m) if m else None
    ceiling = FAMILY_ALPHA / m if m else None
    rows, below, below_ceiling, unavailable = [], 0, 0, 0
    for cell in cells:
        entry = dict(cell)
        if cell.get("status") == "unavailable":
            unavailable += 1
            rows.append(entry)
            continue
        p = float(cell[p_field])
        under = p < LOSS_ALPHA
        under_ceiling = ceiling is not None and p < ceiling
        below += under
        below_ceiling += under_ceiling
        entry["status"] = "below_ceiling" if under_ceiling else ("below_alpha" if under else "holds")
        rows.append(entry)
    holds = unavailable == 0 and allowed is not None and below <= allowed and below_ceiling == 0
    return {
        "rule": rule,
        "p_field": p_field,
        "holds": holds,
        "m": m,
        "alpha": LOSS_ALPHA,
        "allowed": allowed,
        "observed": below,
        "ceiling": ceiling,
        "n_below_ceiling": below_ceiling,
        "n_unavailable": unavailable,
        "observations": rows,
    }


def judge_gain(cells: Sequence[Mapping[str, Any]]) -> dict[str, Any]:
    """G: at least one exact paired gain p-value at most 0.05/m, and no unavailable reading."""

    judged = plsc.judge_gain(cells)
    judged["holds"] = bool(judged["holds"]) and judged["n_unavailable"] == 0
    return judged


def _above_bound(observations: Sequence[Mapping[str, Any]]) -> list[dict[str, Any]]:
    """The readings whose rejected count exceeds the S1 bound α + 2·SE at their own count."""

    out = []
    for row in observations:
        if row.get("status") == "unavailable":
            continue
        n, hits = int(row["n_available"]), int(row["n_rejected"])
        if hits > grade.exceedance_threshold(n):
            out.append({key: row[key] for key in ("profile", "design", "cell", "statistic", "n_available",
                                                  "n_rejected")} | {"rate": hits / n, "bound": grade.bound(n)})
    return out


def judge_arm(
    arm: str,
    plan: Mapping[str, Sequence[Reading]],
    rates: Mapping[tuple[str, ...], Any],
    paired: Mapping[tuple[str, ...], Any],
) -> dict[str, Any]:
    """S1, S2, S3 and G for one candidate arm; every criterion is evaluated whatever fails first."""

    s1 = {
        family: grade.judge_rates([_observation(rates, reading, arm) for reading in plan[f"S1_{family}"]])
        for family in ("simulation", "weak_size", "plasmode")
    }
    strong = grade.judge_rates([_observation(rates, reading, arm) for reading in plan["S3_strong"]])
    band_rows = []
    for reading in plan["S3_band"]:
        row = _paired(paired, reading, arm)
        production = _observation(rates, reading, REFERENCE_ARM)
        row["production_rate"] = (
            None if production.get("status") == "unavailable"
            else production["n_rejected"] / production["n_available"]
        )
        band_rows.append(row)
    band = judge_paired(band_rows, p_field="p_gain", rule="count_and_ceiling_on_paired_excess")
    off_target = [_observation(rates, reading, arm) for reading in (*plan["S3_strong"], *plan["S3_band"])]
    above = _above_bound(off_target)
    criteria = {
        "S1": {"holds": all(family["holds"] for family in s1.values()), "families": s1},
        "S2": judge_paired(
            [_paired(paired, reading, arm) for reading in plan["S2"]],
            p_field="p_loss",
            rule="count_and_ceiling_on_paired_loss",
        ),
        "S3": {
            "holds": strong["holds"] and band["holds"],
            "families": {"strong": strong, "band": band},
            "off_target_above_bound": len(above),
            "off_target_above_bound_readings": above,
        },
        "G": judge_gain([_paired(paired, reading, arm) for reading in plan["G"]]),
    }
    recorded_gains = {
        family: judge_gain([_paired(paired, reading, arm) for reading in plan[key]])
        for family, key in (("band", "G_band"), ("weak", "G_weak"), ("plasmode", "G_plasmode"))
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
        "verdict": _VERDICT_OF[failing[0]] if failing else "qualifies",
        "failing_criteria": failing,
        "failing_readings": {name: rows for name, rows in failing_readings.items() if name in failing},
        "criteria": criteria,
        "recorded_gains": recorded_gains,
    }


def _gate(
    reproduction: Mapping[str, Mapping[str, Any]], accounting_summaries: Mapping[str, Mapping[str, Any]]
) -> dict[str, Any]:
    failures = []
    for profile in PROFILES:
        entry = reproduction.get(profile) or {}
        for name, check in (entry.get("checks") or {}).items():
            if not check.get("passed"):
                failures.append({"profile": profile, "reason": f"{name} failed"})
        if not entry.get("passed", False) and not any(f["profile"] == profile for f in failures):
            failures.append({"profile": profile, "reason": "reproduction failed"})
        inconsistent = int((accounting_summaries.get(profile) or {}).get("inconsistent_units", 0))
        if inconsistent:
            failures.append({"profile": profile, "reason": f"{inconsistent} inconsistent unit(s)"})
    return {"passed": not failures, "failures": failures}


def decision(
    tables: Mapping[str, Mapping[str, Sequence[Mapping[str, Any]]]],
    *,
    reproduction: Mapping[str, Mapping[str, Any]],
    accounting_summaries: Mapping[str, Mapping[str, Any]],
) -> dict[str, Any]:
    """The predeclared verdict (rule version 1) from the three profiles' rate, paired and cell tables."""

    rates, paired = _index_rates(tables), _index_paired(tables)
    plan = criteria_plan()
    gate = _gate(reproduction, accounting_summaries)
    arms = {arm: judge_arm(arm, plan, rates, paired) for arm in CANDIDATE_ARMS}
    qualifying = [arm for arm in CANDIDATE_ARMS if arms[arm]["verdict"] == "qualifies"]
    named, reason = None, None
    if not gate["passed"]:
        verdict = "reproduction_failed"
    elif qualifying:
        verdict = "propose_paper_grade"
        counts = {arm: int(arms[arm]["criteria"]["S3"]["off_target_above_bound"]) for arm in qualifying}
        named = min(qualifying, key=lambda arm: (counts[arm], CANDIDATE_ARMS.index(arm)))
        tied = [arm for arm in qualifying if counts[arm] == counts[named]]
        reason = (
            f"{named} has the fewest off-target readings above the bound among the qualifying arms ({counts})"
            + (f"; tied with {[arm for arm in tied if arm != named]}, preferred as nearer production" if len(tied) > 1
               else "")
            + "."
        )
    else:
        verdict = "keep_production"
    return {
        "rule_version": RULE_VERSION,
        "verdict": verdict,
        "named_arm": named,
        "preference_reason": reason,
        "qualifying_arms": qualifying,
        "consequence": CONSEQUENCES[verdict],
        "consequences": dict(CONSEQUENCES),
        "gate": gate,
        "reference_arm": REFERENCE_ARM,
        "preference": list(CANDIDATE_ARMS),
        "alpha": ALPHA,
        "rejection": dict(REJECTION_RULES),
        "plan": {name: [list(reading) for reading in readings] for name, readings in plan.items()},
        "arms": arms,
        "descriptive": descriptive_tables(tables, plan),
        "notes": [
            "S1 and the delta = 2.0 half of S3 hold each candidate to its own rates under the stage-span paper "
            "grade's count-plus-ceiling rule; S2, G and the band half of S3 compare it with pls_cv on the same "
            "units, fit and draws (exact paired tests).",
            "The band's, the weak design's and the plasmode's gains are recorded with the same test and do not "
            "decide G.",
            "predictive_fitted measures in the model's response units and predictive_orthonormal in the input "
            "metric; delta is compared across arms only through p-values and rejections.",
        ],
    }


def _rate(rates: Mapping[tuple[str, ...], Any], reading: Reading, arm: str) -> float | None:
    row = rates.get((reading[0], reading[1], reading[2], arm, reading[3]))
    return None if row is None else _number(row.get("rate"))


def _cell_row(cells: Mapping[tuple[str, ...], Any], profile: str, design: str, cell: str, arm: str) -> Any:
    return cells.get((profile, design, cell, arm)) or {}


def _num(row: Mapping[str, Any], key: str) -> float | None:
    """A cell-table value as a number (the tables are read back as text)."""

    return _number(row.get(key))


def descriptive_tables(
    tables: Mapping[str, Mapping[str, Sequence[Mapping[str, Any]]]], plan: Mapping[str, Sequence[Reading]]
) -> dict[str, Any]:
    """The control and rank-plus-one arms, estimation, the surgery cell, the weak design, band rates, rank."""

    rates, paired = _index_rates(tables), _index_paired(tables)
    cells = {
        (profile, str(row["design"]), str(row["cell"]), str(row["arm"])): row
        for profile, table in tables.items()
        for row in table["cells"]
    }
    sim_targets = [
        ("simulation", design, cell, TARGETS[_cell_mode(cell)])
        for design in STAGE_SPAN_DESIGNS
        for cell in TARGET_CELLS
    ]
    control = []
    for reading in [*sim_targets, *plan["G_plasmode"]]:
        profile, design, cell, statistic = reading
        control.append({
            "profile": profile, "design": design, "cell": cell, "statistic": statistic,
            **{f"rate_{arm}": _rate(rates, reading, arm) for arm in (REFERENCE_ARM, "predictive", CONTROL_ARM)},
            **{f"pooled_stage_{arm}": _num(_cell_row(cells, profile, design, cell, arm), "pooled_stage_median")
               for arm in (REFERENCE_ARM, "predictive", CONTROL_ARM)},
            "against_pls_cv": _paired(paired, reading, CONTROL_ARM),
            "against_predictive": _paired(paired, reading, CONTROL_ARM, "predictive"),
        })
    plus_one = []
    for design in STAGE_SPAN_DESIGNS:
        for cell in TARGET_CELLS:
            for statistic in STATISTICS:
                reading = ("simulation", design, cell, statistic)
                plus_one.append({
                    "design": design, "cell": cell, "statistic": statistic,
                    "role": "target" if statistic == TARGETS[_cell_mode(cell)] else "off_target",
                    "rate_predictive": _rate(rates, reading, "predictive"),
                    "rate_predictive_rank_plus_one": _rate(rates, reading, RANK_PLUS_ONE_ARM),
                    "p_two_sided": _paired(paired, reading, RANK_PLUS_ONE_ARM, "predictive").get("p_two_sided"),
                })
    rank_sensitivity = [
        row for row in plus_one
        if _cell_mode(row["cell"]) == "orientation" and row["statistic"] == "shape"
        and row["design"] in (*STRONG_DESIGNS, *BAND_DESIGNS)
    ]
    estimation = []
    for profile, table in tables.items():
        for row in table["cells"]:
            if profile == "plasmode":
                continue
            estimation.append({
                "profile": profile,
                **{key: row.get(key) for key in ("design", "cell", "trajectory_mode", "arm")},
                "effect_size": _num(row, "effect_size"),
                "n": _num(row, "n"),
                **{f"{kind}_{statistic}_median": _num(row, f"{kind}_{statistic}_median")
                   for kind in ("observed", "population") for statistic in ("delta", "angle")},
            })
    orientation_estimation = [row for row in estimation if row["trajectory_mode"] == "orientation"]
    surgery = []
    for arm in ARMS:
        cell_row = _cell_row(cells, "simulation", SURGERY_DESIGN, SURGERY_CELL[0], arm)
        for statistic in STATISTICS:
            reading = ("simulation", SURGERY_DESIGN, SURGERY_CELL[0], statistic)
            surgery.append({
                "arm": arm, "statistic": statistic, "rate": _rate(rates, reading, arm),
                "group_difference_capture": _num(cell_row, "group_difference_median"),
                "rank_median": _num(cell_row, "rank_median"),
            })
    weak = [
        {"design": WEAK_DESIGN, "cell": reading[2], "statistic": reading[3],
         **{f"rate_{arm}": _rate(rates, reading, arm) for arm in ARMS}}
        for reading in plan["G_weak"]
    ]
    band_rates = [
        {"design": reading[1], "cell": reading[2], "statistic": reading[3],
         **{f"rate_{arm}": _rate(rates, reading, arm) for arm in (REFERENCE_ARM, *CANDIDATE_ARMS)}}
        for reading in plan["S3_band"]
    ]
    rank_auroc = []
    for profile, table in tables.items():
        designs = sorted({str(row["design"]) for row in table["cells"]})
        for design in designs:
            rows = [row for row in table["cells"] if str(row["design"]) == design and row["arm"] == REFERENCE_ARM]
            anchor = next((row for row in rows if row["cell"] in (ANCHOR_CELL, SURGERY_CELL[0])), rows[0] if rows
                          else {})
            medians = [_number(row.get("rank_median")) for row in rows]
            medians = [value for value in medians if value is not None]
            rank_auroc.append({
                "profile": profile, "design": design, "reference_cell": anchor.get("cell"),
                **{key: _num(anchor, key) for key in ("rank_median", "rank_q10", "rank_q90", "cv_mean_auroc_median",
                                                      "stage_r2_median", "pooled_stage_median")},
                "rank_median_range_over_cells": [min(medians), max(medians)] if medians else None,
            })
    return {
        "control_arm": control,
        "rank_plus_one": plus_one,
        "rank_sensitivity": rank_sensitivity,
        "estimation": estimation,
        "orientation_estimation": orientation_estimation,
        "surgery_cell": surgery,
        "weak_targets": weak,
        "band_off_target_rates": band_rates,
        "rank_and_auroc": rank_auroc,
    }


_TABLE_FILES: dict[str, str] = {"rates": "rate_table.csv", "paired": "paired_table.csv", "cells": "cell_table.csv"}


def _read_json(path: Path) -> dict[str, Any]:
    if not path.exists():
        raise SystemExit(f"{path} does not exist; run `merge` first.")
    return json.loads(path.read_text(encoding="utf-8"))


def decide(out_dir: Path) -> dict[str, Any]:
    """Read the three merged profiles and write ``predictive_subspace_decision.json``."""

    tables, reproduction, summaries, complete, inputs = {}, {}, {}, {}, {}
    for profile in PROFILES:
        directory = profile_dir(out_dir, profile)
        settings = lsc.load_settings(directory)
        plsc.verify_parents(settings)
        summary = _read_json(directory / "accounting.json")
        if summary.get("settings_hash") != lsc.settings_hash(settings):
            raise SystemExit(f"{directory}/accounting.json was merged under other settings; run `merge` again.")
        loaded = {}
        for name, filename in _TABLE_FILES.items():
            path = directory / filename
            if not path.exists():
                raise SystemExit(f"{path} does not exist; run `merge` first.")
            loaded[name] = lsc._read_csv(path)
        tables[profile] = loaded
        reproduction[profile] = _read_json(directory / "reproduction.json")
        summaries[profile] = summary
        complete[profile] = {
            key: summary[key]
            for key in ("expected_units", "complete_units", "missing_units", "failed_units", "inconsistent_units")
        }
        for name in ("settings.json", *_TABLE_FILES.values(), "accounting.json", "reproduction.json"):
            inputs[f"{profile}/{name}"] = plsc.sha256_file(directory / name)
    result = decision(tables, reproduction=reproduction, accounting_summaries=summaries)
    result["complete"] = complete
    result["reproduction"] = {
        profile: {
            "passed": reproduction[profile].get("passed"),
            **{name: {key: check.get(key) for key in ("passed", "rows_compared", "rows_matching", "units_checked",
                                                      "units_passed") if key in check}
               for name, check in (reproduction[profile].get("checks") or {}).items()},
        }
        for profile in PROFILES
    }
    result["inputs"] = inputs
    lsc._write_json(Path(out_dir) / DECISION_NAME, result)
    return result


# ── CLI ──────────────────────────────────────────────────────────────────────


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = parser.add_subparsers(dest="command", required=True)
    run = sub.add_parser("run", help="Evaluate (a shard of) one profile's units")
    run.add_argument("--profile", choices=PROFILES, required=True)
    run.add_argument("--out-dir", type=Path, required=True)
    run.add_argument("--shard", type=lsc.parse_shard, default=(0, 1), help="i/N, round-robin over the units")
    run.add_argument("--limit", type=int, default=None, help="Evaluate at most this many pending units")
    run.add_argument("--workers", type=int, default=1, help="Units evaluated in parallel (not part of the settings)")
    run.add_argument("--perms", type=int, default=None,
                     help=f"Permutations (simulation {SIMULATION_PERMUTATIONS}, weak_size {WEAK_PERMUTATIONS}, "
                          f"plasmode {PLASMODE_PERMUTATIONS})")
    run.add_argument("--reps", type=int, default=None,
                     help=f"Replicates per cell (simulation {SIMULATION_REPLICATES}, weak_size {WEAK_REPLICATES}, "
                          f"plasmode anchor {PLASMODE_ANCHOR_REPLICATES})")
    run.add_argument("--plant-reps", type=int, default=None,
                     help=f"Plasmode: replicates of the plant cells (default {PLASMODE_PLANT_REPLICATES})")
    run.add_argument("--designs", type=str, default=",".join(SIMULATION_DESIGNS), help="Simulation: comma-separated")
    run.add_argument("--export-dir", type=Path, default=None, help="Plasmode: the SEA-AD cohort_30 export")
    merge_parser = sub.add_parser("merge", help="Merge one profile's shards")
    merge_parser.add_argument("--profile", choices=PROFILES, required=True)
    merge_parser.add_argument("--out-dir", type=Path, required=True)
    decide_parser = sub.add_parser("decide", help=f"Write {DECISION_NAME} from the three profiles")
    decide_parser.add_argument("--out-dir", type=Path, required=True)
    return parser


def settings_from_args(args: argparse.Namespace) -> tuple[dict[str, Any], Any]:
    if args.profile == "simulation":
        return simulation_settings(
            replicates=SIMULATION_REPLICATES if args.reps is None else args.reps,
            permutations=SIMULATION_PERMUTATIONS if args.perms is None else args.perms,
            designs=[name for name in args.designs.split(",") if name],
        ), None
    if args.profile == "weak_size":
        return weak_size_settings(
            replicates=WEAK_REPLICATES if args.reps is None else args.reps,
            permutations=WEAK_PERMUTATIONS if args.perms is None else args.perms,
        ), None
    if args.export_dir is None:
        raise SystemExit("The plasmode profile needs --export-dir.")
    try:
        export = pc.load_export(args.export_dir)
    except pc.ExportContractError as exc:
        raise pc.contract_exit(exc) from exc
    settings = plasmode_settings(
        export,
        replicates=PLASMODE_ANCHOR_REPLICATES if args.reps is None else args.reps,
        plant_replicates=PLASMODE_PLANT_REPLICATES if args.plant_reps is None else args.plant_reps,
        permutations=PLASMODE_PERMUTATIONS if args.perms is None else args.perms,
    )
    return settings, export


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    if args.command == "run":
        index, count = args.shard
        settings, export = settings_from_args(args)
        evaluated = run_shard(args.out_dir, settings, index, count, export=export, limit=args.limit,
                              workers=args.workers)
        print(f"evaluated {evaluated} {args.profile} unit(s)", file=sys.stderr)
        return 0
    if args.command == "merge":
        print(json.dumps(merge(args.out_dir, args.profile), indent=1))
        return 0
    result = decide(args.out_dir)
    print(json.dumps({"verdict": result["verdict"], "named_arm": result["named_arm"],
                      "arms": {arm: entry["verdict"] for arm, entry in result["arms"].items()},
                      "gate": result["gate"], "complete": result["complete"]}, indent=1))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
