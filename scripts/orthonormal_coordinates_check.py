#!/usr/bin/env python3
"""Orthonormal span coordinates against the PLS scores: size, specificity and power on two parent runs (no R).

The driver behind ``openspec/changes/archive/2026-11-09-orthonormal-span-coordinates/`` (spec
``orthonormal-coordinates-validation``, roadmap item 15). The production
measurement computes ``delta``, ``angle`` and ``shape`` in PLS score
coordinates ``X_in·R``, a sheared image of the fitted span; the harness's
opt-in ``latent_coordinates="orthonormal"`` measures ``X_in·Q`` instead,
``Q = qr(R)``. This driver measures both on the same units before the
production measurement may change.

**Units.** Every unit generates one dataset, fits the production stage-PLS
space once (net displacement, rank floor 2, span capture) and evaluates that
fit twice at the unit's RRPP seed, serially: in score coordinates (``scores``,
the production measurement) and on the orthonormal basis of the same fit
(``orthonormal``, built by ``pls_orthonormal_from_pls`` without refitting).
The arms therefore share the data, the fit, the rank and every permutation
draw. Each unit checks, to 1e-8, that the orthonormal basis is orthonormal and
that the score matrix equals the orthonormal matrix times ``Qᵀ·R`` (same
span); a unit failing a check is recorded as inconsistent with the check and
its deviation. A unit that raises is recorded as failed and the shard goes on.

**Profiles.**

- ``paper_grade`` — the stage-span paper grade's units
  (``results/stage-span-paper-grade-2026-11-04/``): the five columns (``n200``,
  ``n80``, ``seaad``, ``n200_three_blocks``, ``n200_four_stages``) × 13 cells
  (the anchor; translation, orientation and shape at e ∈ {0.1, 0.25, 0.5} from
  ``stage_span_paper_grade_span_geometry.json``; magnitude at
  e ∈ {0.02, 0.05, 0.1} from ``stage_span_paper_grade_span_magnitude.json``) ×
  replicates 0–499 × 999 permutations, 32,500 units, under the planted scale.
  Unit parameters come from ``stage_span_leak_check.column_cells`` /
  ``unit_params`` (imported unchanged): the profiles' matched seeds, so every
  generator seed and the RRPP seed (the cell's ``evaluation.seed``, else the
  replicate seed, as the study engine does) equal the paper grade's.
- ``weak`` — the latent-space comparison's ``seaad_weak`` design (δ = 0.11, the
  committed calibration; held-out AUROC 0.617): its anchor and translation
  (e = 1.0) at replicates 0–499 × 999 permutations, 1,000 units. Settings come
  from ``span_split_check.simulation_settings`` and seeds and parameters from
  ``latent_space_comparison`` (both imported unchanged), so every unit is the
  span-split check's simulation unit; its settings must match the span-split
  run's committed ones.

**Parents.** The settings record the sha256 of every parent file (the two
paper-grade profiles and operating tables; the span-split run's settings and
production table, and the comparison's calibration). ``run``, ``merge`` and
``decide`` refuse to proceed when a recorded file has changed, naming it.

Subcommands (every profile's output lands under ``--out-dir``/``<profile>/``):

- ``run --profile paper_grade|weak`` — evaluate (a shard of) the units into
  ``records/shard_<i>_of_<N>.jsonl``, one line per unit. Shards are
  round-robin. The settings and their hash are written to ``settings.json``;
  a resume under other settings exits naming the differing keys. Resume skips
  complete units.
- ``merge --profile …`` — ``merged.jsonl`` (per unit; keep off the
  repository), ``rate_table.csv`` (per column or design, cell, statistic and
  arm: available and rejected counts at p < α, rate, Monte Carlo SE),
  ``paired_table.csv`` (per cell and statistic over units available in both
  arms: both rates, their difference, the discordant counts b (only
  ``scores`` rejects) and c (only ``orthonormal`` rejects), the exact
  two-sided paired p-value and the one-sided loss p-value), ``cell_table.csv``
  (per cell and arm: median rank, AUROC, capture, population B − A path-length
  difference as a share of A's path, pooled eigengap and ``angle`` null 95th
  percentile), ``accounting.json`` (units, failures by cell, inconsistent
  units, duplicates, recorded core-hours) and ``reproduction.json``. The
  ``paper_grade`` reproduction requires the ``scores`` arm's available and
  rejected counts to equal the committed ``design_point_operating.csv``
  ``n_available`` / ``n_rejected`` for every column, cell and statistic (design
  phases; the anchor from ``span_geometry``); the ``weak`` reproduction requires
  the committed span-split ``production_table.csv`` count and rejections for
  the weak design's anchor and translation.
- ``decide`` — ``orthonormal_coordinates_decision.json`` (rule version 1) from
  both merged profiles. The gate ``reproduction_failed`` (either reproduction
  failed or any unit is inconsistent; O1–O3 are still recorded), then:

  - **O1** size: the ``orthonormal`` arm's anchor and translation rates at the
    five columns (m = 60) and the weak design's anchor and translation
    (m = 6), each family under the stage-span paper grade's count-plus-ceiling
    rule (``stage_span_paper_grade_decision.judge_rates``, imported unchanged);
    both must hold.
  - **O2** specificity: the ``orthonormal`` arm's off-target rates of the
    orientation (``delta``, ``shape``), shape (``delta``, ``angle``) and
    magnitude (``angle``, ``shape``) cells (m = 90), same rule.
  - **O3** power non-inferiority: per target cell (m = 45), the one-sided exact
    paired p-value P(Binomial(b + c, ½) ≥ b) for a loss of ``orthonormal``
    against ``scores`` (1 when b + c = 0); it holds when at most the 95th
    percentile of Binomial(45, 0.025) cells are below 0.025 and none is below
    0.05/45.

  The verdict is the first failure in the order ``size_excess`` →
  ``specificity_excess`` → ``power_loss``, else ``adopt_orthonormal``; an
  unavailable reading fails its criterion. Descriptive tables beside it never
  change it: the watched cells (the 45° rotation's ``shape`` at ``n200`` and
  ``n200_four_stages``, both arms with the exact two-sided paired p-value),
  the paired target and off-target differences, the population residual per
  cell and arm (shape at e = 0.5 named per column), and the eigengap and
  ``angle`` null-width changes.

With ``--workers`` above 1, pin BLAS to one thread per process
(``OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1``).

Examples:

    export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
    # smoke
    python scripts/orthonormal_coordinates_check.py run --profile paper_grade --reps 2 --perms 19 \\
        --out-dir /tmp/occ --workers 7
    python scripts/orthonormal_coordinates_check.py run --profile weak --reps 2 --perms 19 --out-dir /tmp/occ
    python scripts/orthonormal_coordinates_check.py merge --profile paper_grade --out-dir /tmp/occ
    python scripts/orthonormal_coordinates_check.py merge --profile weak --out-dir /tmp/occ
    python scripts/orthonormal_coordinates_check.py decide --out-dir /tmp/occ

    # one SLURM array task per shard (scripts/orthonormal_coordinates_array.sbatch), then merge and decide
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
from scipy.stats import binom

sys.path.insert(0, str(Path(__file__).resolve().parent))
import latent_space_comparison as lsc  # noqa: E402  (the weak design's cells, seeds and parameters)
import span_split_check as ssc  # noqa: E402  (the weak design's settings and units)
import stage_span_leak_check as leak  # noqa: E402  (the paper grade's column cells and unit parameters)
import stage_span_paper_grade_decision as grade  # noqa: E402  (the count-plus-ceiling rule, operating tables)

sys.path.pop(0)

DRIVER_VERSION = 1
RULE_VERSION = 1
PROFILES: tuple[str, ...] = ("paper_grade", "weak")
ARMS: tuple[str, ...] = ("scores", "orthonormal")
STATISTICS: tuple[str, ...] = ("delta", "angle", "shape")
TARGETS: dict[str, str] = {"magnitude": "delta", "orientation": "angle", "shape": "shape"}
ALPHA = 0.05
REJECTION = "p < alpha"
CONSISTENCY_TOLERANCE = 1e-8
CONSISTENCY_CHECKS: tuple[str, ...] = ("orthonormal_basis", "same_span")

DEFAULT_REPLICATES = 500
DEFAULT_PERMUTATIONS = 999

REPO = Path(__file__).resolve().parents[1]
PROFILE_DIR = REPO / "examples" / "trajectory_power_study"
PAPER_GRADE_DIR = REPO / "results" / "stage-span-paper-grade-2026-11-04"
PAPER_GRADE_PROFILES: dict[str, Path] = {
    "span_geometry": PROFILE_DIR / "stage_span_paper_grade_span_geometry.json",
    "span_magnitude": PROFILE_DIR / "stage_span_paper_grade_span_magnitude.json",
}
SPAN_SPLIT_DIR = REPO / "results" / "span-split-2026-11-05" / "simulation"
PAPER_GRADE_COLUMNS: tuple[str, ...] = grade.COLUMNS
GEOMETRY_EFFECTS: tuple[float, ...] = grade.SPAN_GEOMETRY_EFFECTS
MAGNITUDE_EFFECTS: tuple[float, ...] = grade.SPAN_MAGNITUDE_EFFECTS
#: (mode, effects, source profile) — the paper grade's 13 cells per column.
PAPER_GRADE_CELLS: tuple[tuple[str, tuple[float, ...], str], ...] = (
    ("none", (0.0,), "span_geometry"),
    ("translation", GEOMETRY_EFFECTS, "span_geometry"),
    ("orientation", GEOMETRY_EFFECTS, "span_geometry"),
    ("shape", GEOMETRY_EFFECTS, "span_geometry"),
    ("magnitude", MAGNITUDE_EFFECTS, "span_magnitude"),
)
WEAK_DESIGN = "seaad_weak"
WEAK_CELLS: tuple[str, ...] = ssc.CONTROL_CELLS
#: The weak design's controls as (design, mode, effect): the comparison's anchor and translation (e = 1.0).
WEAK_CONTROLS: tuple[tuple[str, str, float], ...] = tuple(
    (WEAK_DESIGN, mode, float(effect)) for name, mode, effect, *_ in lsc.CELLS if name in WEAK_CELLS
)

#: O3 (design D7): one-sided loss level per cell, count quantile and family level.
LOSS_ALPHA = 0.025
LOSS_COUNT_QUANTILE = 0.95
LOSS_FAMILY_ALPHA = 0.05
#: The watched cells: the 45° rotation's ``shape`` (stage-span orientation at e = 0.5).
WATCHED: tuple[tuple[str, str, float, str], ...] = (
    ("n200", "orientation", 0.5, "shape"),
    ("n200_four_stages", "orientation", 0.5, "shape"),
)
RESIDUAL_MODE, RESIDUAL_EFFECT = "shape", 0.5

VERDICTS: tuple[str, ...] = ("reproduction_failed", "size_excess", "specificity_excess", "power_loss",
                             "adopt_orthonormal")
_VERDICT_OF = {"O1": "size_excess", "O2": "specificity_excess", "O3": "power_loss"}
CONSEQUENCES: dict[str, str] = {
    "reproduction_failed": (
        "Nothing is read: the scores arm does not reproduce its parent runs, or a unit's two arms do not measure "
        "one span. Find the drift before any rerun."
    ),
    "size_excess": (
        "Production stays on scores: the orthonormal coordinates exceed the size bound at the failing controls. "
        "The setting stays an opt-in diagnostic and the score shear a documented limitation."
    ),
    "specificity_excess": (
        "Production stays on scores: the orthonormal coordinates exceed the off-target bound at the failing cells. "
        "The setting stays an opt-in diagnostic and the score shear a documented limitation."
    ),
    "power_loss": (
        "Production stays on scores: the orthonormal coordinates lose target power against the scores at the "
        "failing cells. The report documents the trade-off and the setting stays an opt-in diagnostic."
    ),
    "adopt_orthonormal": (
        "A separate change may flip the harness default to orthonormal coordinates, pinning every committed "
        "configuration and diagnostic to scores, and open the motco analyze adoption with its own real-noise size "
        "check, attribution under orthonormal coordinates and the SEA-AD re-issue."
    ),
}
DECISION_NAME = "orthonormal_coordinates_decision.json"

#: Unit key: (column or design, trajectory mode, effect, replicate).
UnitKey = tuple[str, str, float, int]


# ── Shared helpers ───────────────────────────────────────────────────────────


def sha256_file(path: Path) -> str:
    return ssc.sha256_file(Path(path))


def _display_path(path: Path) -> str:
    return ssc._display_path(path)


def _resolve(stored: str) -> Path:
    return ssc._resolve(stored)


def _canonical(value: Any) -> Any:
    return ssc._canonical(value)


def _parent(path: Path) -> dict[str, str]:
    if not Path(path).exists():
        raise SystemExit(f"{path} does not exist; the parent file is required.")
    return {"path": _display_path(path), "sha256": sha256_file(path)}


def verify_parents(settings: Mapping[str, Any]) -> None:
    """Refuse when a parent file recorded in the settings is missing or has another sha256, naming it."""

    for role, entry in sorted(settings.get("parents", {}).items()):
        path = _resolve(entry["path"])
        if not path.exists():
            raise SystemExit(f"Parent {role} ({entry['path']}) does not exist; refusing.")
        digest = sha256_file(path)
        if digest != entry["sha256"]:
            raise SystemExit(
                f"Parent {role} ({entry['path']}) changed: sha256 {digest}, the settings record {entry['sha256']}; "
                "refusing."
            )


def profile_dir(out_dir: Path, profile: str) -> Path:
    return Path(out_dir) / profile


def _unit_key(key: UnitKey) -> str:
    column, mode, effect, replicate = key
    return json.dumps([column, mode, float(effect), int(replicate)])


def _record_key(record: Mapping[str, Any]) -> str:
    return json.dumps(
        [str(record["column"]), str(record["trajectory_mode"]), float(record["effect_size"]), int(record["replicate"])]
    )


def _rejects(p_value: Any, alpha: float = ALPHA) -> bool | None:
    return ssc._rejects(p_value, alpha, REJECTION)


def _median(values: Iterable[Any]) -> float | None:
    return ssc._median(values)


@functools.lru_cache(maxsize=1)
def _reference() -> Any:
    from motco.simulations.reference import load_reference

    return load_reference()


def generate_dataset(generator: Any) -> Any:
    from motco.simulations.semisynthetic import generate_semisynthetic_trajectory

    return generate_semisynthetic_trajectory(generator, reference=_reference())


# ── Settings and units ───────────────────────────────────────────────────────


@functools.lru_cache(maxsize=4)
def _profile_columns(path: str) -> dict[str, dict[str, Any]]:
    return leak.column_cells(Path(path))


def paper_grade_settings(
    *,
    replicates: int = DEFAULT_REPLICATES,
    permutations: int = DEFAULT_PERMUTATIONS,
    profiles: Mapping[str, Path] = PAPER_GRADE_PROFILES,
    parent_dir: Path = PAPER_GRADE_DIR,
) -> dict[str, Any]:
    """The stage-span paper grade's columns and cells, its profiles and operating tables bound by sha256."""

    parents: dict[str, dict[str, str]] = {}
    for name in sorted(profiles):
        parents[f"profile:{name}"] = _parent(profiles[name])
        parents[f"operating:{name}"] = _parent(Path(parent_dir) / name / "report" / "design_point_operating.csv")
    for mode, _, name in PAPER_GRADE_CELLS:
        columns = _profile_columns(str(profiles[name]))
        missing = [column for column in PAPER_GRADE_COLUMNS if mode not in columns.get(column, {})]
        if missing:
            raise SystemExit(f"Profile {name} has no {mode!r} cell at column(s) {missing}.")
    return {
        "driver_version": DRIVER_VERSION,
        "profile": "paper_grade",
        "parents": parents,
        "columns": list(PAPER_GRADE_COLUMNS),
        "cells": [[mode, float(effect), name] for mode, effects, name in PAPER_GRADE_CELLS for effect in effects],
        "span_scale": "planted",
        "span_capture": True,
        "replicates": int(replicates),
        "permutations": int(permutations),
        "alpha": ALPHA,
        "rejection": REJECTION,
        "arms": list(ARMS),
        "consistency_tolerance": CONSISTENCY_TOLERANCE,
    }


#: The span-split run's unit settings the weak profile must share.
_SHARED_UNIT_KEYS: tuple[str, ...] = (
    "cells", "base_seed", "pls_cv_params", "min_components", "orientation_estimator", "weak_delta",
)


def weak_settings(
    *,
    replicates: int = DEFAULT_REPLICATES,
    permutations: int = DEFAULT_PERMUTATIONS,
    parent_dir: Path = SPAN_SPLIT_DIR,
    comparison_dir: Path = ssc.COMPARISON_DIR,
) -> dict[str, Any]:
    """The span-split check's weak-design control units, checked against its committed settings."""

    units = ssc.simulation_settings(
        parent_dir=comparison_dir,
        replicates=0,
        control_replicates=replicates,
        permutations=permutations,
        designs=[WEAK_DESIGN],
    )
    split = ssc.load_parent(parent_dir)
    theirs, ours = split["settings"]["comparison"], units["comparison"]
    differing = [key for key in _SHARED_UNIT_KEYS if _canonical(ours[key]) != _canonical(theirs[key])]
    if _canonical(ours["designs"][WEAK_DESIGN]) != _canonical(theirs["designs"].get(WEAK_DESIGN)):
        differing.append(f"designs.{WEAK_DESIGN}")
    if differing:
        raise SystemExit(
            f"The weak design's units differ from {parent_dir}/settings.json in: {', '.join(differing)}; they would "
            "not be the span-split check's units."
        )
    return {
        "driver_version": DRIVER_VERSION,
        "profile": "weak",
        "parents": {
            "span_split_settings": _parent(Path(parent_dir) / "settings.json"),
            "span_split_production": _parent(Path(parent_dir) / "production_table.csv"),
            "comparison_calibration": _parent(Path(comparison_dir) / "calibration" / "calibration.json"),
        },
        "design": WEAK_DESIGN,
        "units": units,
        "cells": list(WEAK_CELLS),
        "replicates": int(replicates),
        "permutations": int(permutations),
        "alpha": ALPHA,
        "rejection": REJECTION,
        "arms": list(ARMS),
        "consistency_tolerance": CONSISTENCY_TOLERANCE,
    }


def _weak_cell(settings: Mapping[str, Any], cell: str) -> dict[str, Any]:
    return lsc.cell_spec(settings["units"]["comparison"], cell)


def enumerate_units(settings: Mapping[str, Any]) -> list[UnitKey]:
    replicates = range(int(settings["replicates"]))
    if settings["profile"] == "paper_grade":
        return [
            (column, mode, float(effect), replicate)
            for column in settings["columns"]
            for mode, effect, _ in settings["cells"]
            for replicate in replicates
        ]
    keys: list[UnitKey] = []
    for unit in ssc.simulation_units(settings["units"]):
        spec = _weak_cell(settings, unit.cell)
        keys.append((unit.design, spec["mode"], float(spec["effect"]), unit.replicate))
    return keys


def _cell_label(settings: Mapping[str, Any], mode: str, effect: float) -> str:
    if settings["profile"] == "weak":
        for cell in settings["cells"]:
            spec = _weak_cell(settings, cell)
            if spec["mode"] == mode and float(spec["effect"]) == float(effect):
                return cell
        raise KeyError(f"The weak design has no {mode!r} cell at {effect}.")
    return "anchor" if mode == "none" else f"{mode}_e{effect:g}"


def _source_profile(settings: Mapping[str, Any], mode: str) -> str:
    for cell_mode, _, name in settings["cells"]:
        if cell_mode == mode:
            return str(name)
    raise KeyError(f"The paper grade has no {mode!r} cell.")


def unit_params(key: UnitKey, settings: Mapping[str, Any]) -> tuple[Any, Any]:
    """The unit's generator and score-coordinate evaluation parameters (the RRPP seed is ``evaluation.seed``)."""

    column, mode, effect, replicate = key
    if settings["profile"] == "paper_grade":
        path = _resolve(settings["parents"][f"profile:{_source_profile(settings, mode)}"]["path"])
        generator, evaluation = leak.unit_params(
            _profile_columns(str(path)),
            (column, mode, float(effect), int(replicate)),
            span_scale=settings["span_scale"],
            permutations=int(settings["permutations"]),
        )
        # The study engine's rule: the cell's evaluation seed, else the replicate seed.
        seed = evaluation.seed if evaluation.seed is not None else generator.seed
        return generator, dataclasses.replace(evaluation, seed=seed, n_jobs=1)
    comparison = settings["units"]["comparison"]
    cell = _cell_label(settings, mode, effect)
    spec = lsc.cell_spec(comparison, cell)
    base_seed = comparison["base_seed"]
    generator = lsc._generator_params(
        comparison["designs"][column]["generator"],
        seed=lsc.generator_seed(base_seed, column, replicate),
        trajectory_mode=spec["mode"],
        group_effect_size=float(spec["effect"]),
        span_alignment=spec["phi"],
        span_support=spec["support"],
    )
    evaluation = lsc.arm_params(
        lsc.REFERENCE_ARM,
        comparison,
        seed=lsc.rrpp_seed(base_seed, column, cell, replicate),
        n_stages=int(comparison["designs"][column]["generator"]["n_stages"]),
    )
    return generator, evaluation


# ── One unit ─────────────────────────────────────────────────────────────────


def _finite(values: Mapping[str, Any]) -> dict[str, float | None]:
    return lsc._finite(values)


def _number(value: Any) -> float | None:
    if value is None or value == "":
        return None
    number = float(value)
    return number if math.isfinite(number) else None


def consistency_checks(scores: Any, orthonormal: Any, *, tolerance: float = CONSISTENCY_TOLERANCE) -> dict[str, Any]:
    """Largest deviation of ``QᵀQ = I`` and of ``T = Z·(Qᵀ·R)``; inconsistent above ``tolerance``."""

    basis = np.asarray(orthonormal.linear_map.basis, dtype=float)
    rotations = np.asarray(scores.artifacts.model.x_rotations_, dtype=float)
    through_basis = orthonormal.matrix.to_numpy(dtype=float) @ (basis.T @ rotations)
    deviations = {
        "orthonormal_basis": float(np.max(np.abs(basis.T @ basis - np.eye(basis.shape[1])))),
        "same_span": float(np.max(np.abs(scores.matrix.to_numpy(dtype=float) - through_basis))),
    }
    failing = [{"check": name, "deviation": value} for name, value in deviations.items() if not value <= tolerance]
    return {"tolerance": tolerance, "deviations": deviations, "consistent": not failing, "failing": failing}


def arm_record(result: Any) -> dict[str, Any]:
    """One arm's statistics, p-values, spectrum, null width, capture and population checkpoint."""

    from motco.simulations.diagnostics import POPULATION_LATENT_CHECKPOINT

    metadata = result.latent_matrix_metadata
    geometry = result.realized_geometry
    return {
        "rank": int(metadata["selected_lv"]),
        "cv_mean_auroc": metadata.get("cv_mean_auroc"),
        "pair_statistics": _finite(result.pair_statistics),
        "p_values": _finite(result.p_values),
        "pooled_eigengap": (result.config_spectrum.get("pooled") or {}).get("relative_eigengap"),
        "angle_null_q95": (result.null_summary.get("angle") or {}).get("q95"),
        "span_capture": geometry.get("span_capture"),
        "group_levels": list(result.group_levels),
        "population_latent": (geometry.get("checkpoints", {}).get(POPULATION_LATENT_CHECKPOINT) or {}).get("joint"),
    }


def measure_unit(dataset: Any, params: Any, *, tolerance: float = CONSISTENCY_TOLERANCE) -> dict[str, Any]:
    """Fit ``pls`` once and evaluate it in both coordinate systems at ``params.seed``."""

    from motco.simulations.evaluation import (
        evaluate_semisynthetic_trajectory,
        integrate_semisynthetic_dataset,
        pls_orthonormal_from_pls,
    )

    if params.integration_method != "pls" or params.latent_coordinates is not None:
        raise ValueError("measure_unit takes the score-coordinate pls parameters.")
    timing: dict[str, float] = {}
    mark = perf_counter()
    scores_latent = integrate_semisynthetic_dataset(dataset, params)
    timing["fit"] = perf_counter() - mark
    orthonormal_params = dataclasses.replace(params, latent_coordinates="orthonormal")
    orthonormal_latent = pls_orthonormal_from_pls(dataset, scores_latent, orthonormal_params)
    consistency = consistency_checks(scores_latent, orthonormal_latent, tolerance=tolerance)
    arms = {}
    for arm, arm_params, latent in (
        ("scores", params, scores_latent),
        ("orthonormal", orthonormal_params, orthonormal_latent),
    ):
        mark = perf_counter()
        arms[arm] = arm_record(evaluate_semisynthetic_trajectory(dataset, arm_params, latent=latent))
        timing[arm] = perf_counter() - mark
    return {"consistency": consistency, "arms": arms, "timing": timing}


def evaluate_unit(key: UnitKey, settings: Mapping[str, Any]) -> dict[str, Any]:
    """One unit's record; a failure is recorded with its error, never raised."""

    start = perf_counter()
    column, mode, effect, replicate = key
    base: dict[str, Any] = {
        "profile": settings["profile"],
        "column": column,
        "cell": _cell_label(settings, mode, effect),
        "trajectory_mode": mode,
        "effect_size": float(effect),
        "replicate": int(replicate),
        "settings_hash": lsc.settings_hash(settings),
    }
    timing: dict[str, float] = {}
    try:
        generator, params = unit_params(key, settings)
        base.update({"generator_seed": int(generator.seed), "rrpp_seed": int(params.seed)})
        mark = perf_counter()
        dataset = generate_dataset(generator)
        timing["generate"] = perf_counter() - mark
        measured = measure_unit(dataset, params, tolerance=float(settings["consistency_tolerance"]))
        timing.update(measured["timing"])
    except Exception as exc:  # recorded, never aborts the shard
        timing["unit"] = perf_counter() - start
        return {**base, "status": "failed", "error_type": type(exc).__name__, "error_message": str(exc),
                "runtime_seconds": timing}
    timing["unit"] = perf_counter() - start
    return {**base, "status": "completed", "consistency": measured["consistency"], "arms": measured["arms"],
            "runtime_seconds": timing}


# ── Run ──────────────────────────────────────────────────────────────────────


def read_records(paths: Iterable[Path]) -> dict[str, dict[str, Any]]:
    return ssc.read_records(paths, _record_key)


def run_shard(
    out_dir: Path, settings: Mapping[str, Any], index: int, count: int, *, limit: int | None = None, workers: int = 1
) -> int:
    """Evaluate this shard's missing units into its file; returns how many were evaluated."""

    profile = settings["profile"]
    directory = profile_dir(out_dir, profile)
    if (directory / "settings.json").exists():
        verify_parents(lsc.load_settings(directory))
    verify_parents(settings)
    lsc.ensure_settings(directory, settings)
    path = lsc.shard_path(directory, index, count)
    path.parent.mkdir(parents=True, exist_ok=True)
    done = read_records([path]) if path.exists() else {}
    # Rewritten with complete lines only, so an interrupted unit is re-run, not duplicated.
    with path.open("w", encoding="utf-8") as handle:
        for key in sorted(done):
            handle.write(lsc._line(done[key]))
    pending = [
        key for position, key in enumerate(enumerate_units(settings))
        if position % count == index and _unit_key(key) not in done
    ]
    if limit is not None:
        pending = pending[:limit]

    def append(position: int, record: Mapping[str, Any]) -> None:
        with path.open("a", encoding="utf-8") as handle:
            handle.write(lsc._line(record))
            handle.flush()
        print(f"[{profile} {index}/{count}] {position}/{len(pending)} {_record_key(record)} {record['status']}",
              file=sys.stderr, flush=True)

    if workers <= 1:
        for position, key in enumerate(pending, start=1):
            append(position, evaluate_unit(key, settings))
    else:
        with ProcessPoolExecutor(max_workers=workers) as pool:
            futures = [pool.submit(evaluate_unit, key, settings) for key in pending]
            for position, future in enumerate(as_completed(futures), start=1):
                append(position, future.result())
    return len(pending)


# ── Merge ────────────────────────────────────────────────────────────────────


RATE_COLUMNS = (
    "profile", "column", "cell", "trajectory_mode", "effect_size", "statistic", "arm", "n_units", "n_failed",
    "n_available", "n_rejected", "rate", "se",
)
PAIRED_COLUMNS = (
    "profile", "column", "cell", "trajectory_mode", "effect_size", "statistic", "n_paired", "rate_scores",
    "rate_orthonormal", "difference", "b", "c", "p_two_sided", "p_loss",
)
CELL_COLUMNS = (
    "profile", "column", "cell", "trajectory_mode", "effect_size", "arm", "n", "median_rank", "median_cv_auroc",
    "median_capture", "median_population_share", "median_population_difference", "median_pooled_eigengap",
    "median_angle_null_q95",
)


def paired_p_values(b: int, c: int) -> tuple[float, float]:
    """Exact paired (McNemar) p-values: two-sided, and one-sided for a loss of ``orthonormal`` (b ≥ c)."""

    n = int(b) + int(c)
    if n == 0:
        return 1.0, 1.0
    two_sided = min(1.0, 2.0 * float(binom.cdf(min(b, c), n, 0.5)))
    loss = float(binom.sf(int(b) - 1, n, 0.5))
    return two_sided, loss


def _cells(records: Sequence[Mapping[str, Any]]) -> dict[tuple[str, str, str, float], list[Mapping[str, Any]]]:
    grouped: dict[tuple[str, str, str, float], list[Mapping[str, Any]]] = {}
    for record in records:
        key = (str(record["column"]), str(record["cell"]), str(record["trajectory_mode"]), float(record["effect_size"]))
        grouped.setdefault(key, []).append(record)
    return grouped


def _cell_fields(profile: str, key: tuple[str, str, str, float]) -> dict[str, Any]:
    column, cell, mode, effect = key
    return {"profile": profile, "column": column, "cell": cell, "trajectory_mode": mode, "effect_size": effect}


def _arm_p(record: Mapping[str, Any], arm: str, statistic: str) -> Any:
    if record.get("status") != "completed":
        return None
    return ((record.get("arms") or {}).get(arm) or {}).get("p_values", {}).get(statistic)


def _population_difference(arm: Mapping[str, Any]) -> tuple[float | None, float | None]:
    """B − A population path-length difference and its share of A's path (group levels in contrast order)."""

    population, levels = arm.get("population_latent"), arm.get("group_levels") or []
    if not population or len(levels) < 2:
        return None, None
    lengths = population.get("path_lengths") or {}
    a, b = lengths.get(str(levels[0])), lengths.get(str(levels[1]))
    if a is None or b is None:
        return None, None
    return float(b) - float(a), (float(b) - float(a)) / float(a) if float(a) > 0 else None


def tables(records: Sequence[Mapping[str, Any]], settings: Mapping[str, Any]) -> dict[str, list[dict[str, Any]]]:
    """The rate, paired and cell tables of one profile's records."""

    profile = settings["profile"]
    alpha = float(settings["alpha"])
    rates, paired, cells = [], [], []
    for key, group in sorted(_cells(records).items()):
        fields = _cell_fields(profile, key)
        failed = sum(1 for record in group if record.get("status") != "completed")
        for statistic in STATISTICS:
            for arm in ARMS:
                decisions = [_rejects(_arm_p(record, arm, statistic), alpha) for record in group]
                available = [bool(value) for value in decisions if value is not None]
                count, hits = len(available), sum(available)
                rate = hits / count if count else None
                rates.append({
                    **fields, "statistic": statistic, "arm": arm, "n_units": len(group), "n_failed": failed,
                    "n_available": count, "n_rejected": hits, "rate": rate,
                    "se": None if rate is None else math.sqrt(rate * (1.0 - rate) / count),
                })
            pairs = [
                (_rejects(_arm_p(record, "scores", statistic), alpha),
                 _rejects(_arm_p(record, "orthonormal", statistic), alpha))
                for record in group
            ]
            pairs = [(s, o) for s, o in pairs if s is not None and o is not None]
            n = len(pairs)
            b = sum(1 for s, o in pairs if s and not o)
            c = sum(1 for s, o in pairs if o and not s)
            two_sided, loss = paired_p_values(b, c)
            rate_scores = sum(s for s, _ in pairs) / n if n else None
            rate_orthonormal = sum(o for _, o in pairs) / n if n else None
            paired.append({
                **fields, "statistic": statistic, "n_paired": n, "rate_scores": rate_scores,
                "rate_orthonormal": rate_orthonormal,
                "difference": None if n == 0 else (c - b) / n,
                "b": b, "c": c, "p_two_sided": two_sided, "p_loss": loss,
            })
        completed = [record for record in group if record.get("status") == "completed"]
        for arm in ARMS:
            entries = [record["arms"][arm] for record in completed]
            population = [_population_difference(entry) for entry in entries]
            cells.append({
                **fields, "arm": arm, "n": len(entries),
                "median_rank": _median(entry["rank"] for entry in entries),
                "median_cv_auroc": _median(entry.get("cv_mean_auroc") for entry in entries),
                "median_capture": _median((entry.get("span_capture") or {}).get("group_difference")
                                          for entry in entries),
                "median_population_share": _median(share for _, share in population),
                "median_population_difference": _median(difference for difference, _ in population),
                "median_pooled_eigengap": _median(entry.get("pooled_eigengap") for entry in entries),
                "median_angle_null_q95": _median(entry.get("angle_null_q95") for entry in entries),
            })
    return {"rates": rates, "paired": paired, "cells": cells}


def accounting(
    records: Sequence[Mapping[str, Any]], settings: Mapping[str, Any], *, expected: int, duplicates: int
) -> dict[str, Any]:
    failures: dict[str, int] = {}
    inconsistent = []
    unit_seconds: list[float] = []
    stages: dict[str, list[float]] = {}
    max_deviation = {name: 0.0 for name in CONSISTENCY_CHECKS}
    for record in records:
        label = f"{record['column']}/{record['cell']}"
        if record.get("status") != "completed":
            failures[label] = failures.get(label, 0) + 1
        else:
            consistency = record.get("consistency") or {}
            for name, value in (consistency.get("deviations") or {}).items():
                max_deviation[name] = max(max_deviation.get(name, 0.0), float(value))
            if not consistency.get("consistent", False):
                inconsistent.append({"unit": _record_key(record), "failing": consistency.get("failing")})
        runtime = record.get("runtime_seconds") or {}
        if runtime.get("unit") is not None:
            unit_seconds.append(float(runtime["unit"]))
        for stage in ("generate", "fit", "scores", "orthonormal"):
            if runtime.get(stage) is not None:
                stages.setdefault(stage, []).append(float(runtime[stage]))
    return {
        "profile": settings["profile"],
        "settings_hash": lsc.settings_hash(settings),
        "expected_units": expected,
        "complete_units": len(records),
        "missing_units": expected - len(records),
        "failed_units": sum(failures.values()),
        "failures_by_cell": failures,
        "inconsistent_units": len(inconsistent),
        "inconsistent": inconsistent,
        "max_consistency_deviation": max_deviation,
        "duplicate_records": duplicates,
        "recorded_core_hours": sum(unit_seconds) / 3600.0,
        "unit_seconds": {
            "median": _median(unit_seconds),
            "q90": float(np.quantile(unit_seconds, 0.9)) if unit_seconds else None,
            "max": max(unit_seconds) if unit_seconds else None,
        },
        "stage_seconds_median": {stage: _median(values) for stage, values in sorted(stages.items())},
    }


def committed_operating(settings: Mapping[str, Any]) -> dict[tuple[str, str, float], dict[str, dict[str, Any]]]:
    """The paper grade's committed counts per (column, mode, effect), design phases, the anchor from span_geometry."""

    loaded = {
        name: grade.load_rates(_resolve(settings["parents"][f"operating:{name}"]["path"]))
        for name in sorted({str(name) for _, _, name in settings["cells"]})
    }
    reference: dict[tuple[str, str, float], dict[str, dict[str, Any]]] = {}
    for column in settings["columns"]:
        for mode, effect, name in settings["cells"]:
            key = (column, mode, float(effect))
            if key in loaded[name]:
                reference[key] = loaded[name][key]
    return reference


def paper_grade_reproduction(rate_rows: Sequence[Mapping[str, Any]], settings: Mapping[str, Any]) -> dict[str, Any]:
    """The ``scores`` arm against the committed operating tables, exact on available and rejected counts."""

    committed = committed_operating(settings)
    ours = {
        (row["column"], row["trajectory_mode"], float(row["effect_size"]), row["statistic"]): row
        for row in rate_rows
        if row["arm"] == "scores"
    }
    differing, compared = [], 0
    for column in settings["columns"]:
        for mode, effect, _ in settings["cells"]:
            for statistic in STATISTICS:
                compared += 1
                theirs = (committed.get((column, mode, float(effect))) or {}).get(statistic)
                mine = ours.get((column, mode, float(effect), statistic))
                expected = None if theirs is None else (int(theirs["n_available"]), int(theirs["n_rejected"]))
                got = None if mine is None else (int(mine["n_available"]), int(mine["n_rejected"]))
                if expected is None or got != expected:
                    differing.append({
                        "column": column, "trajectory_mode": mode, "effect_size": float(effect), "statistic": statistic,
                        "committed": None if expected is None else dict(zip(("n_available", "n_rejected"), expected)),
                        "reproduced": None if got is None else dict(zip(("n_available", "n_rejected"), got)),
                    })
    return {
        "source": {key: settings["parents"][key] for key in sorted(settings["parents"]) if key.startswith("operating")},
        "rule": "scores n_available and n_rejected (p < alpha) equal the committed design_point_operating.csv",
        "rows_compared": compared,
        "passed": not differing,
        "differing": differing,
    }


def weak_reproduction(
    rate_rows: Sequence[Mapping[str, Any]], settings: Mapping[str, Any], committed_rows: Sequence[Mapping[str, Any]]
) -> dict[str, Any]:
    """The ``scores`` arm against the span-split check's committed production table, exact."""

    committed = {(row["design"], row["cell"], row["statistic"]): row for row in committed_rows}
    ours = {(row["column"], row["cell"], row["statistic"]): row for row in rate_rows if row["arm"] == "scores"}
    rows, differing = [], []
    for cell in settings["cells"]:
        for statistic in STATISTICS:
            key = (settings["design"], cell, statistic)
            theirs, mine = committed.get(key), ours.get(key)
            entry = {
                "design": settings["design"], "cell": cell, "statistic": statistic,
                "committed": None if theirs is None else {"count": int(theirs["count"]),
                                                          "rejections": int(theirs["rejections"])},
                "reproduced": None if mine is None else {"count": int(mine["n_available"]),
                                                         "rejections": int(mine["n_rejected"])},
            }
            entry["matches"] = (
                entry["committed"] is not None and entry["reproduced"] is not None
                and entry["committed"]["count"] == entry["reproduced"]["count"]
                and entry["committed"]["rejections"] == entry["reproduced"]["rejections"]
            )
            rows.append(entry)
            if not entry["matches"]:
                differing.append(entry)
    return {
        "source": settings["parents"]["span_split_production"],
        "rule": "scores count and rejections (p < alpha) equal the span-split check's production_table.csv",
        "rows_compared": len(rows),
        "rows": rows,
        "passed": not differing,
        "differing": differing,
    }


def merge(out_dir: Path, profile: str) -> dict[str, Any]:
    """Merge one profile's shards and write its tables, accounting and reproduction check."""

    directory = profile_dir(out_dir, profile)
    settings = lsc.load_settings(directory)
    verify_parents(settings)
    expected_keys = [_unit_key(key) for key in enumerate_units(settings)]
    raw: list[dict[str, Any]] = []
    for path in sorted((directory / "records").glob("shard_*.jsonl")):
        raw.extend(read_records([path]).values())
    stale = sorted({_record_key(r) for r in raw if r.get("settings_hash") != lsc.settings_hash(settings)})
    if stale:
        raise SystemExit(f"{len(stale)} record(s) were written under other settings, e.g. {stale[:3]}; refusing.")
    seen: dict[str, int] = {}
    for record in raw:
        seen[_record_key(record)] = seen.get(_record_key(record), 0) + 1
    by_key = {_record_key(r): r for r in raw}
    records = [by_key[key] for key in expected_keys if key in by_key]
    with (directory / "merged.jsonl").open("w", encoding="utf-8") as handle:
        handle.write("".join(lsc._line(record) for record in records))
    summary = accounting(
        records, settings, expected=len(expected_keys), duplicates=sum(count - 1 for count in seen.values())
    )
    merged = tables(records, settings)
    lsc._write_csv(directory / "rate_table.csv", merged["rates"], RATE_COLUMNS)
    lsc._write_csv(directory / "paired_table.csv", merged["paired"], PAIRED_COLUMNS)
    lsc._write_csv(directory / "cell_table.csv", merged["cells"], CELL_COLUMNS)
    if profile == "paper_grade":
        reproduction = paper_grade_reproduction(merged["rates"], settings)
    else:
        source = _resolve(settings["parents"]["span_split_production"]["path"])
        reproduction = weak_reproduction(merged["rates"], settings, lsc._read_csv(source))
    lsc._write_json(directory / "accounting.json", summary)
    lsc._write_json(directory / "reproduction.json", reproduction)
    return {
        "profile": profile,
        **{key: summary[key] for key in ("complete_units", "missing_units", "failed_units", "inconsistent_units",
                                         "recorded_core_hours")},
        "reproduction_passed": reproduction["passed"],
        "reproduction_differing": len(reproduction["differing"]),
    }


# ── Decision (pure functions) ────────────────────────────────────────────────


def _index(rows: Iterable[Mapping[str, Any]], *fields: str) -> dict[tuple, Mapping[str, Any]]:
    out: dict[tuple, Mapping[str, Any]] = {}
    for row in rows:
        out[tuple(float(row[f]) if f == "effect_size" else str(row[f]) for f in fields)] = row
    return out


def _observation(rates: Mapping[tuple, Mapping[str, Any]], column: str, mode: str, effect: float, statistic: str,
                 profile: str) -> dict[str, Any]:
    cell = {"profile": profile, "column": column, "trajectory_mode": mode, "effect_size": float(effect),
            "statistic": statistic}
    row = rates.get((column, mode, float(effect), statistic, "orthonormal"))
    if row is None or int(row["n_available"]) <= 0:
        return {**cell, "status": "unavailable"}
    return {**cell, "n_available": int(row["n_available"]), "n_rejected": int(row["n_rejected"])}


def _paired(paired: Mapping[tuple, Mapping[str, Any]], column: str, mode: str, effect: float,
            statistic: str) -> dict[str, Any] | None:
    row = paired.get((column, mode, float(effect), statistic))
    if row is None:
        return None
    return {
        "column": column, "trajectory_mode": mode, "effect_size": float(effect), "statistic": statistic,
        "n_paired": int(row["n_paired"]),
        "rate_scores": _number(row["rate_scores"]), "rate_orthonormal": _number(row["rate_orthonormal"]),
        "difference": _number(row["difference"]), "b": int(row["b"]), "c": int(row["c"]),
        "p_two_sided": _number(row["p_two_sided"]), "p_loss": _number(row["p_loss"]),
    }


def loss_allowed(m: int) -> int:
    """The 95th percentile of Binomial(m, 0.025): the cells a null run has below 0.025."""

    return int(binom.ppf(LOSS_COUNT_QUANTILE, m, LOSS_ALPHA))


def judge_losses(cells: Sequence[Mapping[str, Any]]) -> dict[str, Any]:
    """O3: the count-plus-ceiling rule on the one-sided paired loss p-values."""

    m = len(cells)
    allowed = loss_allowed(m) if m else None
    ceiling = LOSS_FAMILY_ALPHA / m if m else None
    rows, below, below_ceiling, unavailable = [], 0, 0, 0
    for cell in cells:
        entry = dict(cell)
        if cell.get("status") == "unavailable":
            unavailable += 1
            rows.append(entry)
            continue
        p = float(cell["p_loss"])
        under = p < LOSS_ALPHA
        under_ceiling = ceiling is not None and p < ceiling
        below += under
        below_ceiling += under_ceiling
        entry["status"] = "below_ceiling" if under_ceiling else ("loss" if under else "holds")
        rows.append(entry)
    holds = unavailable == 0 and allowed is not None and below <= allowed and below_ceiling == 0
    return {
        "rule": "count_and_ceiling_on_paired_loss",
        "holds": holds,
        "m": m,
        "loss_alpha": LOSS_ALPHA,
        "allowed_losses": allowed,
        "observed_losses": below,
        "ceiling": ceiling,
        "n_below_ceiling": below_ceiling,
        "n_unavailable": unavailable,
        "observations": rows,
    }


def _target_cells() -> list[tuple[str, str, float]]:
    return [
        (column, mode, float(effect))
        for column in PAPER_GRADE_COLUMNS
        for mode, effects in (("orientation", GEOMETRY_EFFECTS), ("shape", GEOMETRY_EFFECTS),
                              ("magnitude", MAGNITUDE_EFFECTS))
        for effect in effects
    ]


def decision(
    paper: Mapping[str, Sequence[Mapping[str, Any]]],
    weak: Mapping[str, Sequence[Mapping[str, Any]]],
    *,
    reproduction: Mapping[str, Mapping[str, Any]],
    accounting_by_profile: Mapping[str, Mapping[str, Any]],
) -> dict[str, Any]:
    """The predeclared verdict from both profiles' rate, paired and cell tables (rule version 1)."""

    paper_rates = _index(paper["rates"], "column", "trajectory_mode", "effect_size", "statistic", "arm")
    weak_rates = _index(weak["rates"], "column", "trajectory_mode", "effect_size", "statistic", "arm")
    paper_paired = _index(paper["paired"], "column", "trajectory_mode", "effect_size", "statistic")

    gate_failures = []
    for profile in PROFILES:
        if not (reproduction.get(profile) or {}).get("passed", False):
            gate_failures.append({"profile": profile, "reason": "reproduction failed"})
        inconsistent = int((accounting_by_profile.get(profile) or {}).get("inconsistent_units", 0) or 0)
        if inconsistent:
            gate_failures.append({"profile": profile, "reason": f"{inconsistent} inconsistent unit(s)"})
    gate = {"passed": not gate_failures, "failures": gate_failures}

    size_paper = [
        _observation(paper_rates, column, "none", 0.0, statistic, "paper_grade")
        for column in PAPER_GRADE_COLUMNS
        for statistic in STATISTICS
    ] + [
        _observation(paper_rates, column, "translation", effect, statistic, "paper_grade")
        for column in PAPER_GRADE_COLUMNS
        for effect in GEOMETRY_EFFECTS
        for statistic in STATISTICS
    ]
    size_weak = [
        _observation(weak_rates, column, mode, effect, statistic, "weak")
        for column, mode, effect in WEAK_CONTROLS
        for statistic in STATISTICS
    ]
    o1_paper, o1_weak = grade.judge_rates(size_paper), grade.judge_rates(size_weak)
    o1 = {"holds": o1_paper["holds"] and o1_weak["holds"], "families": {"paper_grade": o1_paper, "weak": o1_weak}}

    off_target = [
        _observation(paper_rates, column, mode, effect, statistic, "paper_grade")
        for column, mode, effect in _target_cells()
        for statistic in STATISTICS
        if statistic != TARGETS[mode]
    ]
    o2 = grade.judge_rates(off_target)

    losses = []
    for column, mode, effect in _target_cells():
        row = _paired(paper_paired, column, mode, effect, TARGETS[mode])
        if row is None or row["n_paired"] <= 0:
            losses.append({"column": column, "trajectory_mode": mode, "effect_size": effect,
                           "statistic": TARGETS[mode], "status": "unavailable"})
        else:
            losses.append(row)
    o3 = judge_losses(losses)

    criteria = {"O1": o1, "O2": o2, "O3": o3}
    failures = [
        {"criterion": "O1", "family": family, **row}
        for family, judged in o1["families"].items()
        for row in judged["observations"]
        if row["status"] != "holds"
    ]
    failures += [{"criterion": "O2", **row} for row in o2["observations"] if row["status"] != "holds"]
    failures += [{"criterion": "O3", **row} for row in o3["observations"] if row["status"] != "holds"]
    failing = [name for name in ("O1", "O2", "O3") if not criteria[name]["holds"]]
    if not gate["passed"]:
        verdict = "reproduction_failed"
    else:
        verdict = _VERDICT_OF[failing[0]] if failing else "adopt_orthonormal"
    return {
        "rule_version": RULE_VERSION,
        "verdict": verdict,
        "consequence": CONSEQUENCES[verdict],
        "gate": gate,
        "failing_criteria": failing,
        "criteria": criteria,
        "failures": failures,
        "alpha": ALPHA,
        "rejection": REJECTION,
        "descriptive": descriptive_tables(paper),
        "notes": [
            "O1 and O2 judge the orthonormal arm under the stage-span paper grade's count-plus-ceiling rule "
            "(imported); O3 compares it with the scores on the same units, data, fit and draws.",
            "Rates within a column share matched-seed noise, so exceedances can cluster; every failing reading is "
            "listed.",
        ],
    }


def descriptive_tables(paper: Mapping[str, Sequence[Mapping[str, Any]]]) -> dict[str, Any]:
    """Watched cells, paired differences, population residuals and spectrum changes; never the verdict."""

    paired = _index(paper["paired"], "column", "trajectory_mode", "effect_size", "statistic")
    cells = _index(paper["cells"], "column", "trajectory_mode", "effect_size", "arm")
    watched = [
        _paired(paired, column, mode, effect, statistic) or {
            "column": column, "trajectory_mode": mode, "effect_size": effect, "statistic": statistic,
            "status": "unavailable",
        }
        for column, mode, effect, statistic in WATCHED
    ]
    differences = []
    for column, mode, effect in _target_cells():
        for statistic in STATISTICS:
            row = _paired(paired, column, mode, effect, statistic)
            if row is not None:
                differences.append({**row, "role": "target" if statistic == TARGETS[mode] else "off_target"})
    keys = sorted({(column, mode, effect) for column, mode, effect, _ in cells})
    residuals, spectrum = [], []
    for column, mode, effect in keys:
        scores, orthonormal = cells.get((column, mode, effect, "scores")), cells.get((column, mode, effect,
                                                                                       "orthonormal"))
        if scores is None or orthonormal is None:
            continue
        residuals.append({
            "column": column, "trajectory_mode": mode, "effect_size": effect,
            "share_scores": _number(scores["median_population_share"]),
            "share_orthonormal": _number(orthonormal["median_population_share"]),
        })
        entry = {"column": column, "trajectory_mode": mode, "effect_size": effect}
        for field, name in (("median_pooled_eigengap", "eigengap"), ("median_angle_null_q95", "angle_null_q95")):
            before, after = _number(scores[field]), _number(orthonormal[field])
            entry[f"{name}_scores"] = before
            entry[f"{name}_orthonormal"] = after
            entry[f"{name}_change"] = None if before is None or after is None else after - before
        spectrum.append(entry)
    named = [row for row in residuals if row["trajectory_mode"] == RESIDUAL_MODE
             and row["effect_size"] == RESIDUAL_EFFECT]
    return {
        "watched_cells": watched,
        "paired_differences": differences,
        "population_residual_shape_e0.5": named,
        "population_residuals": residuals,
        "spectrum_changes": spectrum,
    }


_TABLE_FILES: dict[str, str] = {"rates": "rate_table.csv", "paired": "paired_table.csv", "cells": "cell_table.csv"}


def _profile_tables(directory: Path) -> dict[str, list[dict[str, Any]]]:
    return {name: ssc._read_csv(directory / filename) for name, filename in _TABLE_FILES.items()}


def _read_json(path: Path) -> dict[str, Any]:
    if not path.exists():
        raise SystemExit(f"{path} does not exist; run `merge` first.")
    return json.loads(path.read_text(encoding="utf-8"))


def decide(out_dir: Path) -> dict[str, Any]:
    """Read both merged profiles and write ``orthonormal_coordinates_decision.json``."""

    loaded, reproduction, summaries, inputs, complete = {}, {}, {}, {}, {}
    for profile in PROFILES:
        directory = profile_dir(out_dir, profile)
        settings = lsc.load_settings(directory)
        verify_parents(settings)
        summary = _read_json(directory / "accounting.json")
        if summary.get("settings_hash") != lsc.settings_hash(settings):
            raise SystemExit(f"{directory}/accounting.json was merged under other settings; run `merge` again.")
        loaded[profile] = _profile_tables(directory)
        reproduction[profile] = _read_json(directory / "reproduction.json")
        summaries[profile] = summary
        complete[profile] = {key: summary[key] for key in ("expected_units", "complete_units", "missing_units",
                                                           "failed_units", "inconsistent_units")}
        for name in ("rate_table.csv", "paired_table.csv", "cell_table.csv", "accounting.json", "reproduction.json",
                     "settings.json"):
            inputs[f"{profile}/{name}"] = sha256_file(directory / name)
    result = decision(loaded["paper_grade"], loaded["weak"], reproduction=reproduction,
                      accounting_by_profile=summaries)
    result["complete"] = complete
    result["reproduction"] = {profile: {key: reproduction[profile][key] for key in ("passed", "rows_compared")}
                              for profile in PROFILES}
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
    run.add_argument("--shard", type=str, default="0/1", help="i/N, round-robin over the units")
    run.add_argument("--limit", type=int, default=None, help="Evaluate at most this many pending units")
    run.add_argument("--workers", type=int, default=1)
    run.add_argument("--perms", type=int, default=DEFAULT_PERMUTATIONS)
    run.add_argument("--reps", type=int, default=DEFAULT_REPLICATES, help="Replicates per cell")
    merge_parser = sub.add_parser("merge", help="Merge one profile's shards")
    merge_parser.add_argument("--profile", choices=PROFILES, required=True)
    merge_parser.add_argument("--out-dir", type=Path, required=True)
    decide_parser = sub.add_parser("decide", help=f"Write {DECISION_NAME} from both profiles")
    decide_parser.add_argument("--out-dir", type=Path, required=True)
    return parser


def settings_for(profile: str, *, replicates: int, permutations: int) -> dict[str, Any]:
    if profile == "paper_grade":
        return paper_grade_settings(replicates=replicates, permutations=permutations)
    return weak_settings(replicates=replicates, permutations=permutations)


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    if args.command == "run":
        index, count = lsc.parse_shard(args.shard)
        settings = settings_for(args.profile, replicates=args.reps, permutations=args.perms)
        evaluated = run_shard(args.out_dir, settings, index, count, limit=args.limit, workers=args.workers)
        print(f"evaluated {evaluated} {args.profile} unit(s)", file=sys.stderr)
        return 0
    if args.command == "merge":
        print(json.dumps(merge(args.out_dir, args.profile), indent=1))
        return 0
    result = decide(args.out_dir)
    print(json.dumps({"verdict": result["verdict"], "failing": result["failing_criteria"],
                      "complete": result["complete"]}, indent=1))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
