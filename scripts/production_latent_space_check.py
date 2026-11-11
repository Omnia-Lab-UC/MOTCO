#!/usr/bin/env python3
"""Production latent-space selection: stage-PLS against its block-weighted and sparse variants (no R).

The driver behind ``openspec/changes/production-latent-space-selection/``
(spec ``production-latent-space-selection``). Two opt-in settings change the
production stage-PLS fit's input or regularization without changing its stage
objective: block weighting (``block_weighting``, 0.21.0) and sparse stage-PLS
(``sparsity``, 0.22.0). This driver measures both, and both together, against
production on matched units before any default moves.

**Arms** (the weighting × sparsity factorial; each arm is its own
cross-validated stage-PLS fit with net displacement and rank floor 2, at the
unit's RRPP seed with serial draws):

- ``pls_cv`` — production;
- ``weighted`` — ``block_weighting = "equal_variance"`` (1/√p_b per block);
- ``sparse_block`` — sparse stage-PLS-DA with the ``block`` scope over the
  default ladder 1.0/0.5/0.2/0.1/0.05 (``--ladder`` overrides it, recorded in
  the settings);
- ``weighted_sparse_block`` — both.

**Profiles** (every profile's output lands under ``--out-dir``/``<profile>/``):

- ``simulation`` — the latent-space comparison's three designs (``n200`` and
  ``seaad`` at δ = 2.0, ``seaad_weak`` at the committed δ = 0.11), base seed
  2035 and seed functions (``latent_space_comparison``, imported unchanged),
  at 100 replicates × 199 permutations, span capture on. Cells: the
  comparison's seven in-span cells by name (``anchor``, ``translation``,
  ``magnitude`` at e = 0.05, ``orient_e0.25_phi0``, ``orient_e0.5_phi0``,
  ``shape_e0.25``, ``shape_e0.5``) — their datasets and draws are the
  comparison's, so ``pls_cv`` must reproduce its ``rate_table.csv`` rows — and
  three larger plants (``magnitude_e0.25``, ``orient_e1.0_phi0``,
  ``shape_e1.0``) with their own RRPP seeds; 10 cells × 3 designs × 100 =
  3,000 units. Rejection is p < α, the comparison's rule.
- ``plasmode`` — the SEA-AD astrocyte ``cohort_30`` export read through the
  plasmode's input contract (``plasmode_check``, imported unchanged; the five
  files' sha256 must equal the committed plasmode settings', naming any that
  differ). Each replicate (0–99 of the plasmode's ``paper`` stream, base seed
  2034) permutes the sex labels within stage exactly as the plasmode does,
  and every arm analyses the unit with ``analyze_trajectories`` at the
  case-study configuration (six covariates, floor 2, net displacement, ``df``,
  attribution off, the plasmode's CV knobs) at 999 permutations. Cells: the
  ``anchor`` and the real-noise span plant at e ∈ {0.25, 0.5, 1.0}
  (``span_e0.25``, ``span_e0.5``, ``span_e1``); 4 cells × 100 = 400 units.
  Rejection is p ≤ α, the plasmode's rule. The real labels are never
  analysed and no record holds a donor identifier or a label vector.

**The real-noise span plant** (:func:`span_plant`). In the analysis's
standardized units (the export's pooled per-feature mean and SD, the scaler
``analyze_trajectories`` fits), let m̄_s be the pooled stage-s mean over both
groups and m̄ the mean of the m̄_s. Every donor labelled ``Male`` (under the
permuted labels) at stage s gains e · (m̄_s − m̄), mapped back through the
pooled SD, i.e. e · (x̄_s − x̄) in the block's units. The planted group's
centred stage configuration gains exactly e times the pooled one, so B is A
stretched by 1 + e about the centroid up to the permutation's sampling noise
— the stage-span magnitude edit on real noise, in the stage span by
construction. At e = 0 the export's frames are returned untouched.

**Identity checks.** The simulation's ``pls_cv`` count, rejections and failures
must equal the comparison's committed ``rate_table.csv`` for the seven shared
cells at every design (63 rows). Every plasmode anchor unit's ``pls_cv`` arm
must select the plasmode's label-free fidelity values (rank 5, held-out AUROC
0.635556 to six decimals); a mismatch is recorded on the unit, never raised.

**Parents.** The settings record the sha256 of the comparison's
``settings.json``, ``calibration/calibration.json`` and ``rate_table.csv`` and
of the plasmode's ``settings.json`` (with the export's five file hashes);
``run``, ``merge`` and ``decide`` refuse to proceed when one has changed,
naming it.

Subcommands:

- ``run --profile simulation|plasmode`` — evaluate (a shard of) the units
  into ``records/shard_<i>_of_<N>.jsonl``, one line per (unit, arm). Shards
  are round-robin; resume skips complete units; ``settings.json`` is bound by
  its hash and a resume under other settings exits naming the keys. A failing
  arm is recorded as failed for that unit and the shard goes on.
- ``merge --profile …`` — ``merged.jsonl`` (per (unit, arm); keep off the
  repository), ``rate_table.csv`` (per design, cell, arm and statistic:
  count, rejections, failures, rate, Monte Carlo SE), ``paired_table.csv``
  (per cell, statistic and candidate against ``pls_cv`` over units available
  in both: both rates, the difference, the discordant counts b (only
  ``pls_cv`` rejects) and c (only the candidate), exact one-sided loss
  P(Bin(b + c, ½) ≥ b) and gain P(Bin(b + c, ½) ≥ c) p-values and the
  two-sided one), ``cell_table.csv`` (per cell and arm: median and q10/q90 of
  rank, keep fraction, held-out AUROC and ``pooled_stage`` capture, median
  ``group_difference`` capture and population target, median support per
  block), ``accounting.json`` (units, failures by arm and cell, fidelity
  failures, duplicates, core-hours, per-arm timing) and ``reproduction.json``.
- ``decide`` — ``production_latent_space_decision.json`` (rule version 1) from
  both merged profiles: the gate ``reproduction_failed`` (S1–G still
  recorded), then per candidate arm

  - **S1 size** — the anchor and translation ``delta``/``angle``/``shape``
    rates at the three designs (m = 18) and the plasmode anchor (m = 3), each
    family under the stage-span paper grade's count-plus-ceiling rule on
    α + 2·√(α(1 − α)/R) (``stage_span_paper_grade_decision.judge_rates``,
    imported unchanged); both families must hold;
  - **S2 protected power** — the 16 target cells at the δ = 2.0 designs
    (magnitude ``delta`` at e 0.05/0.25, orientation ``angle`` and shape
    ``shape`` at e 0.25/0.5/1.0): at most the 95th percentile of
    Binomial(16, 0.025) exact paired loss p-values below 0.025, none below
    0.05/16;
  - **S3 specificity** — the off-target rates of those cells at e ≤ 0.5
    (m = 24) under the S1 rule; the e = 1.0 off-target rates are descriptive
    (production's own construction leak, the stage-span pilot 2026-11-04);
  - **G gain** — the weak design's 8 target cells and the plasmode's 3 span
    ``delta`` cells: one exact paired gain p-value at most 0.05/11.

  Per arm the first failure in the order ``size_excess``, ``power_loss``,
  ``specificity_excess``, ``no_gain``, else ``qualifies``; overall
  ``propose_paper_grade`` naming the first qualifying arm in the order
  ``weighted``, ``sparse_block``, ``weighted_sparse_block``, else
  ``keep_production``. Descriptive tables beside it never change it.

With ``--workers`` above 1, pin BLAS to one thread per process
(``OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1``).

Examples:

    export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
    EXPORT=../sead_sexdiff/results/motco_input/cohort_30
    # smoke: 2 replicates per cell
    python scripts/production_latent_space_check.py run --profile simulation --reps 2 --out-dir /tmp/plsc \\
        --workers 7
    python scripts/production_latent_space_check.py run --profile plasmode --export-dir $EXPORT --reps 2 \\
        --out-dir /tmp/plsc --workers 7
    python scripts/production_latent_space_check.py merge --profile simulation --out-dir /tmp/plsc
    python scripts/production_latent_space_check.py merge --profile plasmode --out-dir /tmp/plsc
    python scripts/production_latent_space_check.py decide --out-dir /tmp/plsc

    # one SLURM array task per shard (scripts/production_latent_space_array.sbatch), then merge and decide
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
from typing import Any, Callable, Iterable, Mapping, Sequence

import numpy as np
import pandas as pd
from scipy.stats import binom

sys.path.insert(0, str(Path(__file__).resolve().parent))
import latent_space_comparison as lsc  # noqa: E402  (designs, cells, seeds, settings files)
import plasmode_check as pc  # noqa: E402  (export contract, label permutation, case-study parameters)
import stage_span_paper_grade_decision as grade  # noqa: E402  (the count-plus-ceiling rule)

sys.path.pop(0)

DRIVER_VERSION = 1
RULE_VERSION = 1
PROFILES: tuple[str, ...] = ("simulation", "plasmode")
ARMS: tuple[str, ...] = ("pls_cv", "weighted", "sparse_block", "weighted_sparse_block")
REFERENCE_ARM = "pls_cv"
#: Candidates in the order of preference among qualifying arms (the simplest production rule first).
CANDIDATE_ARMS: tuple[str, ...] = ("weighted", "sparse_block", "weighted_sparse_block")
WEIGHTING = "equal_variance"
SPARSITY_SCOPE = "block"
#: (weighted, sparse) per arm.
ARM_FACTORS: dict[str, tuple[bool, bool]] = {
    "pls_cv": (False, False),
    "weighted": (True, False),
    "sparse_block": (False, True),
    "weighted_sparse_block": (True, True),
}
STATISTICS: tuple[str, ...] = ("delta", "angle", "shape")
TARGETS: dict[str, str] = {"magnitude": "delta", "orientation": "angle", "shape": "shape"}
ALPHA = 0.05
REJECTION_RULES: dict[str, str] = {"simulation": "p < alpha", "plasmode": "p <= alpha"}

REPO = Path(__file__).resolve().parents[1]
COMPARISON_DIR = REPO / "results" / "latent-space-comparison-2026-11-04"
PLASMODE_DIR = REPO / "results" / "plasmode-2026-11-02"

# Simulation profile (design D2).
SIMULATION_REPLICATES = 100
SIMULATION_PERMUTATIONS = 199
SHARED_CELLS: tuple[str, ...] = (
    "anchor", "translation", "magnitude", "orient_e0.25_phi0", "orient_e0.5_phi0", "shape_e0.25", "shape_e0.5",
)
#: (cell, trajectory mode, effect, φ or None, support or None), the comparison's cell layout.
ADDED_CELLS: tuple[tuple[str, str, float, float | None, int | None], ...] = (
    ("magnitude_e0.25", "magnitude", 0.25, None, None),
    ("orient_e1.0_phi0", "orientation", 1.0, 0.0, None),
    ("shape_e1.0", "shape", 1.0, None, None),
)
CONTROL_CELLS: tuple[str, ...] = ("anchor", "translation")
STRONG_DESIGNS: tuple[str, ...] = lsc.STRONG_DESIGNS
WEAK_DESIGN = "seaad_weak"
SIMULATION_BLOCKS: tuple[str, ...] = tuple(lsc.LAYERS)

# Plasmode profile (design D3).
PLASMODE_DESIGN = "plasmode"
PLASMODE_REPLICATES = 100
PLASMODE_PERMUTATIONS = 999
PLANT_EFFECTS: tuple[float, ...] = (0.25, 0.5, 1.0)
ANCHOR_CELL = "anchor"
PLANT_MODE = "magnitude"
PLASMODE_BLOCKS: tuple[str, ...] = pc.BLOCKS

# Decision (design D4).
SE_MULTIPLIER = 2.0
LOSS_ALPHA = 0.025
LOSS_COUNT_QUANTILE = 0.95
FAMILY_ALPHA = 0.05
#: S3 reads the target cells' off-target rates up to this effect; above it they are descriptive.
SPECIFICITY_MAX_EFFECT = 0.5
VERDICTS: tuple[str, ...] = ("reproduction_failed", "propose_paper_grade", "keep_production")
ARM_VERDICTS: tuple[str, ...] = ("size_excess", "power_loss", "specificity_excess", "no_gain", "qualifies")
_VERDICT_OF = {"S1": "size_excess", "S2": "power_loss", "S3": "specificity_excess", "G": "no_gain"}
CONSEQUENCES: dict[str, str] = {
    "reproduction_failed": (
        "Nothing is read: the production arm does not reproduce the comparison's committed rates or the "
        "plasmode's label-free latent space. Find the drift before any rerun."
    ),
    "propose_paper_grade": (
        "A separate change runs the named arm at paper grade (500 x 999 on the stage-span paper grade's cells and "
        "the plasmode's grade); only its verdict may flip the harness and motco analyze defaults, with pins for "
        "every committed configuration and the SEA-AD re-issue."
    ),
    "keep_production": (
        "Production keeps the stage-PLS space unchanged; block weighting and sparse stage-PLS stay opt-in and are "
        "documented as measured."
    ),
}
DELTA_METRIC_NOTE = (
    "Under block weighting the PLS space is fitted in the weighted metric (scale=False), so the weighted arms' "
    "delta is in another metric than production's; only p-values and rejections are compared across arms."
)
DECISION_NAME = "production_latent_space_decision.json"


# ── Shared helpers ───────────────────────────────────────────────────────────


def sha256_file(path: Path) -> str:
    return pc.sha256_file(Path(path))


def _display_path(path: Path) -> str:
    resolved = Path(path).resolve()
    try:
        return str(resolved.relative_to(REPO))
    except ValueError:
        return str(resolved)


def _resolve(stored: str) -> Path:
    path = Path(stored)
    return path if path.is_absolute() else REPO / path


def _canonical(value: Any) -> Any:
    return json.loads(json.dumps(value, sort_keys=True))


def _parent(path: Path) -> dict[str, str]:
    if not Path(path).exists():
        raise SystemExit(f"{path} does not exist; the parent file is required.")
    return {"path": _display_path(path), "sha256": sha256_file(path)}


def load_parent(directory: Path) -> dict[str, Any]:
    """A committed parent run's ``settings.json``, refused when its hash does not match its settings."""

    path = Path(directory) / "settings.json"
    if not path.exists():
        raise SystemExit(f"{path} does not exist; the parent run is required.")
    stored = json.loads(path.read_text(encoding="utf-8"))
    if stored.get("hash") != lsc.settings_hash(stored.get("settings", {})):
        raise SystemExit(f"{path} fails its hash: the parent's settings were edited.")
    return stored


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


def verify_export(export: Any, settings: Mapping[str, Any]) -> None:
    """Refuse an export whose file hashes differ from the settings', naming every differing file."""

    recorded = settings["export"]["files"]
    differing = sorted(
        name for name in set(export.files) | set(recorded) if export.files.get(name) != recorded.get(name)
    )
    if differing:
        raise SystemExit(
            f"The export's file hashes differ from the settings' for: {', '.join(differing)}; refusing."
        )


def profile_dir(out_dir: Path, profile: str) -> Path:
    return Path(out_dir) / profile


def _rejects(p_value: Any, alpha: float, rule: str) -> bool | None:
    """``None`` for an unavailable p-value."""

    if p_value is None or p_value == "":
        return None
    value = float(p_value)
    if not math.isfinite(value):
        return None
    return value < alpha if rule == "p < alpha" else value <= alpha


def _finite_list(values: Iterable[Any]) -> list[float]:
    return [float(value) for value in values if value is not None and math.isfinite(float(value))]


def _quantile(values: Iterable[Any], q: float) -> float | None:
    finite = _finite_list(values)
    return float(np.quantile(finite, q)) if finite else None


def _number(value: Any) -> float | None:
    if value is None or value == "":
        return None
    number = float(value)
    return number if math.isfinite(number) else None


def default_ladder() -> tuple[float, ...]:
    from motco.stats.sparse_pls import DEFAULT_KEEP_FRACTIONS

    return tuple(float(value) for value in DEFAULT_KEEP_FRACTIONS)


def arm_settings(ladder: Sequence[float] | None = None) -> dict[str, dict[str, Any]]:
    """Each arm's two settings; the arms differ in nothing else."""

    fractions = [float(value) for value in (default_ladder() if ladder is None else ladder)]
    return {
        arm: {
            "block_weighting": WEIGHTING if weighted else None,
            "sparsity": {"scope": SPARSITY_SCOPE, "keep_fractions": list(fractions)} if sparse else None,
        }
        for arm, (weighted, sparse) in ARM_FACTORS.items()
    }


@functools.lru_cache(maxsize=1)
def _reference() -> Any:
    from motco.simulations.reference import load_reference

    return load_reference()


# ── Simulation profile: settings, units and arms ─────────────────────────────


def simulation_settings(
    *,
    parent_dir: Path = COMPARISON_DIR,
    replicates: int = SIMULATION_REPLICATES,
    permutations: int = SIMULATION_PERMUTATIONS,
    designs: Sequence[str] = lsc.DESIGNS,
    ladder: Sequence[float] | None = None,
) -> dict[str, Any]:
    """The comparison's designs, seeds, weak δ and shared cells, plus the added cells; parents bound by sha256."""

    parent_dir = Path(parent_dir)
    parent = load_parent(parent_dir)
    parent_settings = parent["settings"]
    calibration = parent_dir / "calibration" / "calibration.json"
    if not calibration.exists():
        raise SystemExit(f"{calibration} does not exist; the comparison's calibration is required.")
    digest = sha256_file(calibration)
    expected = (parent_settings.get("weak_delta") or {}).get("calibration_sha256")
    if digest != expected:
        raise SystemExit(
            f"{calibration} has sha256 {digest}, not the comparison's calibration_sha256 {expected}; refusing a "
            "weak δ the comparison did not use."
        )
    comparison = lsc.build_settings(
        calibration=calibration, replicates=replicates, permutations=permutations, designs=designs
    )
    parent_cells = {str(cell[0]): list(cell) for cell in parent_settings["cells"]}
    missing = [name for name in SHARED_CELLS if name not in parent_cells]
    if missing:
        raise SystemExit(f"{parent_dir}/settings.json has no cell(s) {missing}; they cannot be shared.")
    comparison["cells"] = [parent_cells[name] for name in SHARED_CELLS] + [list(cell) for cell in ADDED_CELLS]
    comparison["arms"] = [lsc.REFERENCE_ARM]
    shared = {
        "designs": {name: parent_settings["designs"][name] for name in comparison["designs"]},
        **{
            key: parent_settings[key]
            for key in ("base_seed", "pls_cv_params", "min_components", "orientation_estimator", "weak_delta")
        },
    }
    differing = [key for key, value in shared.items() if _canonical(comparison[key]) != _canonical(value)]
    code_cells = {str(cell[0]): list(cell) for cell in lsc.CELLS}
    differing += [
        f"cells.{name}" for name in SHARED_CELLS if _canonical(parent_cells[name]) != _canonical(code_cells.get(name))
    ]
    if differing:
        raise SystemExit(
            f"The comparison's units differ from {parent_dir}/settings.json in: {', '.join(differing)}; "
            "the shared units would not be its datasets."
        )
    return {
        "driver_version": DRIVER_VERSION,
        "profile": "simulation",
        "parents": {
            "comparison_settings": _parent(parent_dir / "settings.json"),
            "comparison_calibration": _parent(calibration),
            "comparison_rate_table": _parent(parent_dir / "rate_table.csv"),
        },
        "comparison": comparison,
        "shared_cells": list(SHARED_CELLS),
        "added_cells": [cell[0] for cell in ADDED_CELLS],
        "arms": arm_settings(ladder),
        "span_capture": True,
        "replicates": int(replicates),
        "permutations": int(permutations),
        "alpha": ALPHA,
        "rejection": REJECTION_RULES["simulation"],
    }


def simulation_units(settings: Mapping[str, Any]) -> list[lsc.Unit]:
    """Design × cell × replicate, every design and cell of the settings."""

    return lsc.enumerate_units(settings["comparison"])


def simulation_arm_params(arm: str, settings: Mapping[str, Any], *, seed: int, n_stages: int) -> Any:
    """``pls_cv``'s evaluation parameters (the comparison's) with the arm's weighting and sparsity."""

    params = lsc.arm_params(lsc.REFERENCE_ARM, settings["comparison"], seed=seed, n_stages=n_stages)
    spec = settings["arms"][arm]
    integration = dict(params.integration_params)
    if spec["block_weighting"] is not None:
        integration["block_weighting"] = spec["block_weighting"]
    if spec["sparsity"] is not None:
        integration["sparsity"] = {
            "scope": spec["sparsity"]["scope"],
            "keep_fractions": [float(value) for value in spec["sparsity"]["keep_fractions"]],
        }
    return dataclasses.replace(params, integration_params=integration)


def simulation_arm_record(result: Any) -> dict[str, Any]:
    """One arm's statistics, p-values, rank, keep fraction, AUROC, capture, population target and support."""

    from motco.simulations.diagnostics import POPULATION_LATENT_CHECKPOINT

    metadata = result.latent_matrix_metadata
    sparsity = metadata.get("sparsity") or {}
    geometry = result.realized_geometry
    return {
        "rank": int(metadata["selected_lv"]),
        "keep_fraction": sparsity.get("selected_keep_fraction"),
        "support": sparsity.get("support"),
        "block_weighting": (metadata.get("integration_params") or {}).get("block_weighting"),
        "cv_mean_auroc": metadata.get("cv_mean_auroc"),
        "pair_statistics": lsc._finite(result.pair_statistics),
        "p_values": lsc._finite(result.p_values),
        "span_capture": geometry.get("span_capture"),
        "population_latent": (geometry.get("checkpoints", {}).get(POPULATION_LATENT_CHECKPOINT) or {}).get("joint"),
    }


def _failed(base: Mapping[str, Any], arm: str, exc: Exception, seconds: float) -> dict[str, Any]:
    return {
        **base,
        "arm": arm,
        "status": "failed",
        "error_type": type(exc).__name__,
        "error_message": str(exc),
        "arm_seconds": seconds,
    }


def evaluate_simulation_unit(unit: lsc.Unit, settings: Mapping[str, Any]) -> list[dict[str, Any]]:
    """Generate the unit's dataset once and evaluate every arm on it serially; one record per arm."""

    from motco.simulations.evaluation import evaluate_semisynthetic_trajectory

    start = perf_counter()
    comparison = settings["comparison"]
    cell = lsc.cell_spec(comparison, unit.cell)
    seed = lsc.rrpp_seed(comparison["base_seed"], unit.design, unit.cell, unit.replicate)
    n_stages = int(comparison["designs"][unit.design]["generator"]["n_stages"])
    base = {
        "profile": "simulation",
        "design": unit.design,
        "cell": unit.cell,
        "trajectory_mode": cell["mode"],
        "effect_size": float(cell["effect"]),
        "replicate": unit.replicate,
        "generator_seed": lsc.generator_seed(comparison["base_seed"], unit.design, unit.replicate),
        "rrpp_seed": seed,
        "settings_hash": lsc.settings_hash(settings),
    }
    records: list[dict[str, Any]] = []
    mark = perf_counter()
    try:
        dataset = lsc.unit_dataset(unit, comparison, reference=_reference())
    except Exception as exc:  # every arm is recorded as failed; the shard goes on
        records = [_failed(base, arm, exc, 0.0) for arm in ARMS]
        dataset = None
    generate_seconds = perf_counter() - mark
    if dataset is not None:
        for arm in ARMS:
            mark = perf_counter()
            try:
                params = simulation_arm_params(arm, settings, seed=seed, n_stages=n_stages)
                result = evaluate_semisynthetic_trajectory(dataset, params)
                records.append({**base, "arm": arm, "status": "completed", **simulation_arm_record(result),
                                "arm_seconds": perf_counter() - mark})
            except Exception as exc:  # recorded per arm, never aborts the unit
                records.append(_failed(base, arm, exc, perf_counter() - mark))
    unit_seconds = perf_counter() - start
    for record in records:
        record["generate_seconds"] = generate_seconds
        record["unit_runtime_seconds"] = unit_seconds
    return records


# ── Plasmode profile: settings, the span plant, units and arms ───────────────


def plasmode_cells(effects: Sequence[float] = PLANT_EFFECTS) -> list[list[Any]]:
    """``[name, effect]``: the anchor, then one span plant per effect."""

    return [[ANCHOR_CELL, 0.0], *([f"span_e{float(effect):g}", float(effect)] for effect in effects)]


def plasmode_settings(
    export: Any,
    *,
    parent_dir: Path = PLASMODE_DIR,
    replicates: int = PLASMODE_REPLICATES,
    plant_replicates: int | None = None,
    permutations: int = PLASMODE_PERMUTATIONS,
    effects: Sequence[float] = PLANT_EFFECTS,
    ladder: Sequence[float] | None = None,
) -> dict[str, Any]:
    """The plasmode's committed configuration and fidelity values, refused on another export."""

    parent_dir = Path(parent_dir)
    parent = load_parent(parent_dir)
    parent_settings = parent["settings"]
    verify_export(export, {"export": parent_settings["export"]})
    if any(float(effect) <= 0.0 for effect in effects):
        raise SystemExit(f"Plant effects must be positive, got {list(effects)}.")
    return {
        "driver_version": DRIVER_VERSION,
        "profile": "plasmode",
        "parents": {"plasmode_settings": _parent(parent_dir / "settings.json")},
        "export": {"files": dict(export.files), "cohort": export.manifest.get("cohort")},
        "plasmode": parent_settings,
        "design": PLASMODE_DESIGN,
        "cells": plasmode_cells(effects),
        "replicates": {
            "anchor": int(replicates),
            "plant": int(replicates if plant_replicates is None else plant_replicates),
        },
        "permutations": int(permutations),
        "fidelity": dict(parent_settings["fidelity"]),
        "attribution": False,
        "arms": arm_settings(ladder),
        "plant": {
            "group": pc.GROUPS[1],
            "shift": "e * (pooled stage-s mean - mean of the pooled stage means), standardized, mapped back",
            "scaler": "the export's pooled per-feature mean and SD (ddof 0), the analysis's sample scaler",
            "geometry": "the planted group's centred stage configuration gains e times the pooled one",
            "labels": "the plasmode's within-stage permutation of its paper stream, never the real labels",
        },
        "alpha": ALPHA,
        "rejection": REJECTION_RULES["plasmode"],
    }


def plasmode_units(settings: Mapping[str, Any]) -> list[lsc.Unit]:
    """Replicate-major, so a shard holds the matched cells of its replicates."""

    counts = settings["replicates"]
    cells = [str(name) for name, _ in settings["cells"]]
    return [
        lsc.Unit(PLASMODE_DESIGN, cell, replicate)
        for replicate in range(max(int(counts["anchor"]), int(counts["plant"])))
        for cell in cells
        if replicate < int(counts["anchor" if cell == ANCHOR_CELL else "plant"])
    ]


def plasmode_effect(settings: Mapping[str, Any], cell: str) -> float:
    for name, effect in settings["cells"]:
        if name == cell:
            return float(effect)
    raise KeyError(f"Unknown plasmode cell {cell!r}.")


def span_plant(
    blocks: Mapping[str, pd.DataFrame],
    stages: np.ndarray,
    stage_order: Sequence[str],
    labels: np.ndarray,
    effect: float,
    *,
    group: str,
) -> dict[str, pd.DataFrame]:
    """The blocks with the real-noise span plant at ``effect`` on the rows labelled ``group``.

    In each block's standardized units (the pooled per-feature mean and SD the
    analysis fits), the rows of ``group`` at stage s gain ``effect · (m̄_s −
    m̄)``, m̄_s the pooled stage-s mean over every row and m̄ the mean of the
    m̄_s; the shift is mapped back through the pooled SD. At ``effect == 0``
    the frames are returned untouched.
    """

    from motco.simulations.preprocessing import fit_block_scalers

    if float(effect) == 0.0:
        return dict(blocks)
    stages = np.asarray(stages).astype(str)
    planted_rows = np.asarray(labels).astype(str) == str(group)
    scalers = fit_block_scalers(dict(blocks), {})
    planted: dict[str, pd.DataFrame] = {}
    for name, frame in blocks.items():
        values = frame.to_numpy(dtype=float, copy=True)
        scaler = scalers[name]
        standardized = (values - scaler.mean) / scaler.scale
        means = np.vstack([standardized[stages == str(stage)].mean(axis=0) for stage in stage_order])
        shift = float(effect) * (means - means.mean(axis=0))
        for index, stage in enumerate(stage_order):
            rows = planted_rows & (stages == str(stage))
            values[rows] = values[rows] + shift[index] * scaler.scale
        planted[name] = pd.DataFrame(values, index=frame.index, columns=frame.columns)
    return planted


def build_plasmode_unit(unit: lsc.Unit, settings: Mapping[str, Any], export: Any) -> dict[str, Any]:
    """The unit's permuted labels (the plasmode's draw for the replicate), planted blocks and metadata."""

    parent = settings["plasmode"]
    base_seed, stream = int(parent["base_seed"]), str(parent["stream"])
    seeds = pc.unit_seeds(base_seed, stream, unit.replicate)
    labels, redraws = pc.draw_labels(
        export.real_labels, export.stages, export.stage_order, base_seed, stream, unit.replicate
    )
    effect = plasmode_effect(settings, unit.cell)
    metadata = export.metadata.copy()
    metadata[export.group_col] = labels.astype(str)
    return {
        "effect": effect,
        "seeds": seeds,
        "label_redraws": redraws,
        "blocks": span_plant(export.blocks, export.stages, export.stage_order, labels, effect, group=export.groups[1]),
        "metadata": metadata,
    }


def plasmode_arm_params(arm: str, settings: Mapping[str, Any], seeds: Mapping[str, int]) -> Any:
    """The plasmode's case-study parameters (attribution off) with the arm's weighting and sparsity."""

    from motco.analysis import SparsitySettings

    parent = {**settings["plasmode"], "attribution": None}
    params = pc.analysis_params(parent, seeds, permutations=int(settings["permutations"]))
    spec = settings["arms"][arm]
    sparsity = (
        None
        if spec["sparsity"] is None
        else SparsitySettings(spec["sparsity"]["scope"], tuple(float(v) for v in spec["sparsity"]["keep_fractions"]))
    )
    return dataclasses.replace(params, block_weighting=spec["block_weighting"], sparsity=sparsity)


def plasmode_arm_record(result: Any) -> dict[str, Any]:
    """One arm's statistics, p-values, rank, keep fraction, AUROC, support and block weights; no donor value."""

    record = result.record
    latent = record["latent"]
    sparsity = record.get("sparsity") or {}
    weighting = record.get("block_weighting")
    return {
        "rank": int(latent["selected_lv"]),
        "keep_fraction": sparsity.get("selected_keep_fraction"),
        "support": sparsity.get("support"),
        "block_weighting": None if weighting is None else weighting["kind"],
        "block_weights": None if weighting is None else {b["name"]: b["weight"] for b in weighting["blocks"]},
        "cv_mean_auroc": latent.get("cv_mean_auroc"),
        "pair_statistics": lsc._finite(record["statistics"]),
        "p_values": lsc._finite(record["p_values"]),
        "span_capture": None,
        "population_latent": None,
    }


def fidelity_outcome(result: Any, fidelity: Mapping[str, Any]) -> dict[str, Any]:
    """The anchor's production arm against the plasmode's label-free latent space; recorded, never raised."""

    latent = result.record["latent"]
    auroc = latent.get("cv_mean_auroc")
    decimals = int(pc.FIDELITY_GATE["decimals"])
    outcome: dict[str, Any] = {
        "rank": int(latent["selected_lv"]),
        "cv_mean_auroc_rounded": None if auroc is None else round(float(auroc), decimals),
        "expected": dict(fidelity),
    }
    try:
        pc.check_anchor_latent(result, fidelity)
    except pc.ConsistencyError as exc:
        return {**outcome, "passed": False, "message": str(exc)}
    return {**outcome, "passed": True}


def evaluate_plasmode_unit(unit: lsc.Unit, settings: Mapping[str, Any], export: Any) -> list[dict[str, Any]]:
    """Build the unit once and analyse it under every arm serially; one record per arm, seeds and statistics."""

    start = perf_counter()
    effect = plasmode_effect(settings, unit.cell)
    parent = settings["plasmode"]
    seeds = pc.unit_seeds(int(parent["base_seed"]), str(parent["stream"]), unit.replicate)
    base = {
        "profile": "plasmode",
        "design": unit.design,
        "cell": unit.cell,
        "trajectory_mode": "none" if unit.cell == ANCHOR_CELL else PLANT_MODE,
        "effect_size": effect,
        "replicate": unit.replicate,
        "stream": parent["stream"],
        "seeds": seeds,
        "rrpp_seed": seeds["rrpp"],
        "settings_hash": lsc.settings_hash(settings),
    }
    records: list[dict[str, Any]] = []
    mark = perf_counter()
    try:
        built = build_plasmode_unit(unit, settings, export)
        base["label_redraws"] = built["label_redraws"]
    except Exception as exc:  # every arm is recorded as failed; the shard goes on
        records = [_failed(base, arm, exc, 0.0) for arm in ARMS]
        built = None
    generate_seconds = perf_counter() - mark
    if built is not None:
        for arm in ARMS:
            mark = perf_counter()
            try:
                params = plasmode_arm_params(arm, settings, built["seeds"])
                _, result = pc.analyse(built["blocks"], built["metadata"], export, params)
                entry = {**base, "arm": arm, "status": "completed", **plasmode_arm_record(result)}
                if arm == REFERENCE_ARM and unit.cell == ANCHOR_CELL:
                    entry["fidelity"] = fidelity_outcome(result, settings["fidelity"])
                entry["arm_seconds"] = perf_counter() - mark
                records.append(entry)
            except Exception as exc:  # recorded per arm, never aborts the unit
                records.append(_failed(base, arm, exc, perf_counter() - mark))
    unit_seconds = perf_counter() - start
    for record in records:
        record["generate_seconds"] = generate_seconds
        record["unit_runtime_seconds"] = unit_seconds
    return [lsc._jsonable(record) for record in records]


# ── Shards, resume and the run ───────────────────────────────────────────────


def unit_key(design: str, cell: str, replicate: int) -> str:
    return json.dumps([str(design), str(cell), int(replicate)])


def record_key(record: Mapping[str, Any]) -> str:
    return unit_key(record["design"], record["cell"], record["replicate"])


def enumerate_units(settings: Mapping[str, Any]) -> list[lsc.Unit]:
    return simulation_units(settings) if settings["profile"] == "simulation" else plasmode_units(settings)


def read_records(paths: Iterable[Path]) -> list[dict[str, Any]]:
    """Every parsable line; a line torn by an interrupted write is dropped, so its unit is re-run."""

    records: list[dict[str, Any]] = []
    for path in paths:
        with Path(path).open(encoding="utf-8") as handle:
            for line in handle:
                try:
                    record = json.loads(line)
                except json.JSONDecodeError:
                    continue
                if record.get("status") in ("completed", "failed"):
                    records.append(record)
    return records


def complete_units(records: Iterable[Mapping[str, Any]]) -> dict[str, list[Mapping[str, Any]]]:
    """Records grouped by unit in arm order, keeping only units with every arm present (completed or failed)."""

    grouped: dict[str, dict[str, Mapping[str, Any]]] = {}
    for record in records:
        grouped.setdefault(record_key(record), {})[str(record["arm"])] = record
    return {key: [by_arm[arm] for arm in ARMS] for key, by_arm in grouped.items() if set(by_arm) == set(ARMS)}


_WORKER_EXPORT: Any = None


def _init_plasmode_worker(directory: str, files: Mapping[str, str]) -> None:
    global _WORKER_EXPORT
    export = pc.load_export(Path(directory))
    if export.files != dict(files):
        raise SystemExit(f"{directory} changed after the run started (file hashes differ).")
    _WORKER_EXPORT = export


def _evaluate_plasmode_pending(unit: lsc.Unit, settings: Mapping[str, Any]) -> list[dict[str, Any]]:
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
        verify_parents(lsc.load_settings(directory))
    verify_parents(settings)
    lsc.ensure_settings(directory, settings)
    if profile == "plasmode":
        if export is None:
            raise SystemExit("The plasmode profile needs the export.")
        verify_export(export, settings)
    path = lsc.shard_path(directory, index, count)
    path.parent.mkdir(parents=True, exist_ok=True)
    done = complete_units(read_records([path]) if path.exists() else [])
    # Rewritten with complete units only, so an interrupted unit is re-run, not duplicated.
    with path.open("w", encoding="utf-8") as handle:
        for key in sorted(done):
            handle.write("".join(lsc._line(record) for record in done[key]))
    pending = [
        unit
        for position, unit in enumerate(enumerate_units(settings))
        if position % count == index and unit_key(*unit.key) not in done
    ]
    if limit is not None:
        pending = pending[:limit]

    def write(position: int, unit: lsc.Unit, records: Sequence[Mapping[str, Any]]) -> None:
        with path.open("a", encoding="utf-8") as handle:
            handle.write("".join(lsc._line(record) for record in records))
            handle.flush()
        statuses = ",".join(str(record["status"])[0] for record in records)
        print(f"[{profile} {index}/{count}] {position}/{len(pending)} {unit.key} {statuses}", file=sys.stderr,
              flush=True)

    if profile == "simulation":
        if workers <= 1:
            for position, unit in enumerate(pending, start=1):
                write(position, unit, evaluate_simulation_unit(unit, settings))
        else:
            with ProcessPoolExecutor(max_workers=workers) as pool:
                futures = {pool.submit(evaluate_simulation_unit, unit, settings): unit for unit in pending}
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
    "failures",
)
PAIRED_COLUMNS = (
    "design", "cell", "trajectory_mode", "effect_size", "statistic", "arm", "n_paired", "rate_reference",
    "rate_candidate", "difference", "b", "c", "p_loss", "p_gain", "p_two_sided",
)
_QUANTITIES: tuple[str, ...] = ("rank", "keep_fraction", "cv_mean_auroc", "pooled_stage")


def cell_columns(blocks: Sequence[str]) -> tuple[str, ...]:
    return (
        "design", "cell", "trajectory_mode", "effect_size", "arm", "n",
        *(f"{quantity}_{suffix}" for quantity in _QUANTITIES for suffix in ("median", "q10", "q90")),
        "group_difference_median", "target_statistic", "population_target_median",
        *(f"support_{block}_median" for block in blocks),
    )


def _blocks_of(settings: Mapping[str, Any]) -> tuple[str, ...]:
    return SIMULATION_BLOCKS if settings["profile"] == "simulation" else PLASMODE_BLOCKS


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


def _p(record: Mapping[str, Any] | None, statistic: str) -> Any:
    if record is None or record.get("status") != "completed":
        return None
    return (record.get("p_values") or {}).get(statistic)


def paired_p_values(b: int, c: int) -> dict[str, float]:
    """Exact paired p-values from the discordant counts: one-sided loss (b), gain (c) and two-sided; 1 at b = c = 0."""

    n = int(b) + int(c)
    if n == 0:
        return {"p_loss": 1.0, "p_gain": 1.0, "p_two_sided": 1.0}
    return {
        "p_loss": float(binom.sf(int(b) - 1, n, 0.5)),
        "p_gain": float(binom.sf(int(c) - 1, n, 0.5)),
        "p_two_sided": min(1.0, 2.0 * float(binom.cdf(min(int(b), int(c)), n, 0.5))),
    }


def rate_table(records: Sequence[Mapping[str, Any]], settings: Mapping[str, Any]) -> list[dict[str, Any]]:
    """Rejections under the profile's rule, binomial SE, count and failures per design × cell × arm × statistic."""

    alpha, rule = float(settings["alpha"]), str(settings["rejection"])
    rows = []
    for _, group in sorted(_by_cell(records).items()):
        fields = _cell_fields(group[0])
        for arm in ARMS:
            arm_records = [record for record in group if record["arm"] == arm]
            failures = sum(1 for record in arm_records if record.get("status") != "completed")
            for statistic in STATISTICS:
                decisions = [_rejects(_p(record, statistic), alpha, rule) for record in arm_records]
                available = [bool(value) for value in decisions if value is not None]
                count, hits = len(available), sum(available)
                rate = hits / count if count else None
                rows.append({
                    **fields, "arm": arm, "statistic": statistic, "rate": rate,
                    "se": None if rate is None else math.sqrt(rate * (1.0 - rate) / count),
                    "count": count, "rejections": hits, "failures": failures,
                })
    return rows


def paired_table(records: Sequence[Mapping[str, Any]], settings: Mapping[str, Any]) -> list[dict[str, Any]]:
    """Each candidate against ``pls_cv`` per cell and statistic over the units available in both arms."""

    alpha, rule = float(settings["alpha"]), str(settings["rejection"])
    rows = []
    for _, group in sorted(_by_cell(records).items()):
        fields = _cell_fields(group[0])
        units: dict[int, dict[str, Mapping[str, Any]]] = {}
        for record in group:
            units.setdefault(int(record["replicate"]), {})[str(record["arm"])] = record
        for statistic in STATISTICS:
            for arm in CANDIDATE_ARMS:
                pairs = [
                    (_rejects(_p(arms.get(REFERENCE_ARM), statistic), alpha, rule),
                     _rejects(_p(arms.get(arm), statistic), alpha, rule))
                    for _, arms in sorted(units.items())
                ]
                pairs = [(ref, cand) for ref, cand in pairs if ref is not None and cand is not None]
                n = len(pairs)
                b = sum(1 for ref, cand in pairs if ref and not cand)
                c = sum(1 for ref, cand in pairs if cand and not ref)
                rows.append({
                    **fields, "statistic": statistic, "arm": arm, "n_paired": n,
                    "rate_reference": sum(ref for ref, _ in pairs) / n if n else None,
                    "rate_candidate": sum(cand for _, cand in pairs) / n if n else None,
                    "difference": (c - b) / n if n else None,
                    "b": b, "c": c, **paired_p_values(b, c),
                })
    return rows


def _mean_support(record: Mapping[str, Any], block: str) -> float | None:
    support = record.get("support")
    if not support:
        return None
    values = [float(component[block]) for component in support if block in component]
    return float(np.mean(values)) if values else None


def cell_table(records: Sequence[Mapping[str, Any]], settings: Mapping[str, Any]) -> list[dict[str, Any]]:
    """Per cell and arm: rank, keep fraction, AUROC and capture quantiles, population target and support."""

    blocks = _blocks_of(settings)
    rows = []
    for _, group in sorted(_by_cell(records).items()):
        fields = _cell_fields(group[0])
        target = TARGETS.get(fields["trajectory_mode"])
        for arm in ARMS:
            done = [r for r in group if r["arm"] == arm and r.get("status") == "completed"]
            values = {
                "rank": [r.get("rank") for r in done],
                "keep_fraction": [r.get("keep_fraction") for r in done],
                "cv_mean_auroc": [r.get("cv_mean_auroc") for r in done],
                "pooled_stage": [(r.get("span_capture") or {}).get("pooled_stage") for r in done],
            }
            row: dict[str, Any] = {**fields, "arm": arm, "n": len(done)}
            for quantity, entries in values.items():
                for suffix, q in (("median", 0.5), ("q10", 0.1), ("q90", 0.9)):
                    row[f"{quantity}_{suffix}"] = _quantile(entries, q)
            row["group_difference_median"] = _quantile(
                ((r.get("span_capture") or {}).get("group_difference") for r in done), 0.5
            )
            row["target_statistic"] = target
            row["population_target_median"] = (
                None if target is None else _quantile(((r.get("population_latent") or {}).get(target) for r in done),
                                                      0.5)
            )
            for block in blocks:
                row[f"support_{block}_median"] = _quantile((_mean_support(r, block) for r in done), 0.5)
            rows.append(row)
    return rows


def accounting(
    records: Sequence[Mapping[str, Any]], settings: Mapping[str, Any], *, expected: int, duplicates: int
) -> dict[str, Any]:
    """Units, failures by arm and cell, fidelity failures, duplicates, core-hours and per-arm timing."""

    units: dict[str, float] = {}
    generate: list[float] = []
    failures: dict[str, dict[str, int]] = {}
    arm_seconds: dict[str, list[float]] = {}
    fidelity_failures = []
    for record in records:
        key = record_key(record)
        if key not in units:
            units[key] = float(record.get("unit_runtime_seconds") or 0.0)
            if record.get("generate_seconds") is not None:
                generate.append(float(record["generate_seconds"]))
        label = f"{record['design']}/{record['cell']}"
        if record.get("status") != "completed":
            by_cell = failures.setdefault(str(record["arm"]), {})
            by_cell[label] = by_cell.get(label, 0) + 1
        elif record.get("arm_seconds") is not None:
            arm_seconds.setdefault(str(record["arm"]), []).append(float(record["arm_seconds"]))
        fidelity = record.get("fidelity")
        if fidelity is not None and not fidelity.get("passed"):
            fidelity_failures.append({"unit": key, "rank": fidelity.get("rank"),
                                      "cv_mean_auroc_rounded": fidelity.get("cv_mean_auroc_rounded")})
    unit_seconds = list(units.values())
    return {
        "profile": settings["profile"],
        "settings_hash": lsc.settings_hash(settings),
        "expected_units": expected,
        "complete_units": len(units),
        "missing_units": expected - len(units),
        "records": len(records),
        "failed_records": sum(sum(cells.values()) for cells in failures.values()),
        "failures_by_arm_and_cell": failures,
        "fidelity_failures": len(fidelity_failures),
        "fidelity_failed_units": fidelity_failures,
        "duplicate_records": duplicates,
        "recorded_core_hours": sum(unit_seconds) / 3600.0,
        "unit_seconds": {
            "median": _quantile(unit_seconds, 0.5),
            "q90": _quantile(unit_seconds, 0.9),
            "max": max(unit_seconds) if unit_seconds else None,
        },
        "generate_seconds_median": _quantile(generate, 0.5),
        "arm_seconds_median": {arm: _quantile(arm_seconds.get(arm, []), 0.5) for arm in ARMS},
        "arm_seconds_total_hours": {arm: sum(arm_seconds.get(arm, [])) / 3600.0 for arm in ARMS},
    }


def simulation_reproduction(
    rates: Sequence[Mapping[str, Any]], settings: Mapping[str, Any], committed_rows: Sequence[Mapping[str, Any]]
) -> dict[str, Any]:
    """``pls_cv`` against the comparison's committed rate table: exact count, rejections and failures."""

    committed = {
        (str(row["design"]), str(row["cell"]), str(row["statistic"])): row
        for row in committed_rows
        if row["arm"] == lsc.REFERENCE_ARM
    }
    ours = {
        (str(row["design"]), str(row["cell"]), str(row["statistic"])): row
        for row in rates
        if row["arm"] == REFERENCE_ARM
    }
    fields = ("count", "rejections", "failures")
    rows, differing = [], []
    for design in settings["comparison"]["designs"]:
        for cell in settings["shared_cells"]:
            for statistic in STATISTICS:
                key = (design, cell, statistic)
                theirs, mine = committed.get(key), ours.get(key)
                entry = {
                    "design": design, "cell": cell, "statistic": statistic,
                    "committed": None if theirs is None else {field: int(theirs[field]) for field in fields},
                    "reproduced": None if mine is None else {field: int(mine[field]) for field in fields},
                }
                entry["matches"] = entry["committed"] is not None and entry["committed"] == entry["reproduced"]
                rows.append(entry)
                if not entry["matches"]:
                    differing.append(entry)
    return {
        "profile": "simulation",
        "source": settings["parents"]["comparison_rate_table"],
        "rule": "pls_cv count, rejections (p < alpha) and failures equal the comparison's rate_table.csv",
        "rows_compared": len(rows),
        "rows_matching": len(rows) - len(differing),
        "passed": bool(rows) and not differing,
        "rows": rows,
        "differing": differing,
    }


def plasmode_reproduction(records: Sequence[Mapping[str, Any]], settings: Mapping[str, Any]) -> dict[str, Any]:
    """Every anchor unit's ``pls_cv`` arm against the plasmode's label-free fidelity values."""

    anchor = [r for r in records if r["cell"] == ANCHOR_CELL and r["arm"] == REFERENCE_ARM]
    unavailable = [record_key(r) for r in anchor if r.get("status") != "completed" or r.get("fidelity") is None]
    failing = [
        {"unit": record_key(r), "rank": r["fidelity"].get("rank"),
         "cv_mean_auroc_rounded": r["fidelity"].get("cv_mean_auroc_rounded")}
        for r in anchor
        if r.get("status") == "completed" and r.get("fidelity") is not None and not r["fidelity"].get("passed")
    ]
    checked = len(anchor) - len(unavailable)
    return {
        "profile": "plasmode",
        "source": settings["parents"]["plasmode_settings"],
        "rule": "every anchor unit's pls_cv rank and held-out AUROC (6 decimals) equal the plasmode's fidelity values",
        "expected": dict(settings["fidelity"]),
        "anchor_units": len(anchor),
        "units_checked": checked,
        "units_passed": checked - len(failing),
        "unavailable": unavailable,
        "passed": checked > 0 and not failing and not unavailable,
        "differing": failing,
    }


def merge(out_dir: Path, profile: str) -> dict[str, Any]:
    """Merge one profile's shards and write its tables, accounting and reproduction."""

    directory = profile_dir(out_dir, profile)
    settings = lsc.load_settings(directory)
    verify_parents(settings)
    raw = read_records(sorted((directory / "records").glob("shard_*.jsonl")))
    stale = sorted({record_key(r) for r in raw if r.get("settings_hash") != lsc.settings_hash(settings)})
    if stale:
        raise SystemExit(f"{len(stale)} unit(s) were written under other settings, e.g. {stale[:3]}; refusing.")
    seen: dict[tuple[str, str], int] = {}
    for record in raw:
        key = (record_key(record), str(record["arm"]))
        seen[key] = seen.get(key, 0) + 1
    units = complete_units(raw)
    expected = [unit_key(*unit.key) for unit in enumerate_units(settings)]
    records = [record for key in expected if key in units for record in units[key]]
    with (directory / "merged.jsonl").open("w", encoding="utf-8") as handle:
        handle.write("".join(lsc._line(record) for record in records))
    rates = rate_table(records, settings)
    lsc._write_csv(directory / "rate_table.csv", rates, RATE_COLUMNS)
    lsc._write_csv(directory / "paired_table.csv", paired_table(records, settings), PAIRED_COLUMNS)
    lsc._write_csv(directory / "cell_table.csv", cell_table(records, settings), cell_columns(_blocks_of(settings)))
    summary = accounting(records, settings, expected=len(expected),
                         duplicates=sum(count - 1 for count in seen.values()))
    if profile == "simulation":
        source = _resolve(settings["parents"]["comparison_rate_table"]["path"])
        reproduction = simulation_reproduction(rates, settings, lsc._read_csv(source))
    else:
        reproduction = plasmode_reproduction(records, settings)
    lsc._write_json(directory / "accounting.json", summary)
    lsc._write_json(directory / "reproduction.json", reproduction)
    return {
        "profile": profile,
        **{key: summary[key] for key in ("complete_units", "missing_units", "failed_records", "fidelity_failures",
                                         "duplicate_records", "recorded_core_hours")},
        "reproduction_passed": reproduction["passed"],
        "reproduction_differing": len(reproduction["differing"]),
    }


# ── Decision (pure functions) ────────────────────────────────────────────────


def _index(rows: Iterable[Mapping[str, Any]], *fields: str) -> dict[tuple[str, ...], Mapping[str, Any]]:
    return {tuple(str(row[field]) for field in fields): row for row in rows}


def _observation(
    rates: Mapping[tuple[str, ...], Mapping[str, Any]], design: str, cell: str, arm: str, statistic: str
) -> dict[str, Any]:
    entry = {"design": design, "cell": cell, "arm": arm, "statistic": statistic}
    row = rates.get((design, cell, arm, statistic))
    if row is None or int(row["count"]) <= 0:
        return {**entry, "status": "unavailable"}
    return {**entry, "n_available": int(row["count"]), "n_rejected": int(row["rejections"])}


def _paired_row(
    paired: Mapping[tuple[str, ...], Mapping[str, Any]], design: str, cell: str, arm: str, statistic: str
) -> dict[str, Any]:
    entry = {"design": design, "cell": cell, "arm": arm, "statistic": statistic}
    row = paired.get((design, cell, statistic, arm))
    if row is None or int(row["n_paired"]) <= 0:
        return {**entry, "status": "unavailable"}
    return {
        **entry,
        "trajectory_mode": row.get("trajectory_mode"),
        "effect_size": _number(row.get("effect_size")),
        "n_paired": int(row["n_paired"]),
        "rate_reference": _number(row["rate_reference"]),
        "rate_candidate": _number(row["rate_candidate"]),
        "difference": _number(row["difference"]),
        "b": int(row["b"]),
        "c": int(row["c"]),
        "p_loss": float(row["p_loss"]),
        "p_gain": float(row["p_gain"]),
        "p_two_sided": float(row["p_two_sided"]),
    }


def loss_allowed(m: int) -> int:
    """The 95th percentile of Binomial(m, 0.025): the cells a null run has below 0.025."""

    return int(binom.ppf(LOSS_COUNT_QUANTILE, m, LOSS_ALPHA))


def judge_losses(cells: Sequence[Mapping[str, Any]]) -> dict[str, Any]:
    """S2: the count-plus-ceiling rule on the one-sided exact paired loss p-values; unavailable fails."""

    m = len(cells)
    allowed = loss_allowed(m) if m else None
    ceiling = FAMILY_ALPHA / m if m else None
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


def judge_gain(cells: Sequence[Mapping[str, Any]]) -> dict[str, Any]:
    """G: at least one exact paired gain p-value at most 0.05/m; the pooled gain is reported, not judged."""

    m = len(cells)
    threshold = FAMILY_ALPHA / m if m else None
    rows, gains = [], []
    total_b = total_c = total_n = 0
    differences = []
    for cell in cells:
        entry = dict(cell)
        if cell.get("status") == "unavailable":
            rows.append(entry)
            continue
        gain = threshold is not None and float(cell["p_gain"]) <= threshold
        entry["status"] = "gain" if gain else "no_gain"
        rows.append(entry)
        if gain:
            gains.append(entry)
        total_b += int(cell["b"])
        total_c += int(cell["c"])
        total_n += int(cell["n_paired"])
        if cell.get("difference") is not None:
            differences.append(float(cell["difference"]))
    return {
        "rule": "any_exact_paired_gain_at_bonferroni",
        "holds": bool(gains),
        "m": m,
        "threshold": threshold,
        "gains": gains,
        "n_unavailable": sum(1 for cell in cells if cell.get("status") == "unavailable"),
        "pooled": {
            "b": total_b,
            "c": total_c,
            "n_paired": total_n,
            "difference": (total_c - total_b) / total_n if total_n else None,
            "mean_cell_difference": float(np.mean(differences)) if differences else None,
        },
        "observations": rows,
    }


def simulation_cell_specs() -> dict[str, tuple[str, float]]:
    """``cell -> (mode, effect)`` of the simulation profile's ten cells (the comparison's layout)."""

    shared = {str(cell[0]): (str(cell[1]), float(cell[2])) for cell in lsc.CELLS if cell[0] in SHARED_CELLS}
    return {name: shared[name] for name in SHARED_CELLS} | {
        str(cell[0]): (str(cell[1]), float(cell[2])) for cell in ADDED_CELLS
    }


def criteria_plan() -> dict[str, list[tuple[str, str, str]]]:
    """The ``(design, cell, statistic)`` readings of every criterion, fixed before any unit.

    S1: anchor and translation at the three designs (m = 18) and the plasmode
    anchor (m = 3); S2: the δ = 2.0 designs' 16 target cells; S3: their
    off-target rates at e ≤ 0.5 (m = 24), the e = 1.0 ones descriptive; G: the
    weak design's 8 target cells and the plasmode's 3 span ``delta`` cells.
    """

    specs = simulation_cell_specs()
    targets = [(cell, mode, effect) for cell, (mode, effect) in specs.items() if mode in TARGETS]
    plant_cells = [str(name) for name, _ in plasmode_cells() if name != ANCHOR_CELL]

    def off_target(keep: Callable[[float], bool]) -> list[tuple[str, str, str]]:
        return [
            (design, cell, statistic)
            for design in STRONG_DESIGNS
            for cell, mode, effect in targets
            if keep(effect)
            for statistic in STATISTICS
            if statistic != TARGETS[mode]
        ]

    return {
        "S1_simulation": [
            (design, cell, statistic) for design in lsc.DESIGNS for cell in CONTROL_CELLS for statistic in STATISTICS
        ],
        "S1_plasmode": [(PLASMODE_DESIGN, ANCHOR_CELL, statistic) for statistic in STATISTICS],
        "S2": [(design, cell, TARGETS[mode]) for design in STRONG_DESIGNS for cell, mode, _ in targets],
        "S3": off_target(lambda effect: effect <= SPECIFICITY_MAX_EFFECT),
        "S3_descriptive": off_target(lambda effect: effect > SPECIFICITY_MAX_EFFECT),
        "G": [(WEAK_DESIGN, cell, TARGETS[mode]) for cell, mode, _ in targets]
        + [(PLASMODE_DESIGN, cell, TARGETS[PLANT_MODE]) for cell in plant_cells],
    }


def judge_arm(
    arm: str,
    plan: Mapping[str, Sequence[tuple[str, str, str]]],
    rates: Mapping[tuple[str, ...], Mapping[str, Any]],
    paired: Mapping[tuple[str, ...], Mapping[str, Any]],
) -> dict[str, Any]:
    """S1, S2, S3 and G for one candidate arm; every criterion is evaluated whatever fails first."""

    s1_families = {
        family: grade.judge_rates([_observation(rates, *spec[:2], arm, spec[2]) for spec in plan[key]])
        for family, key in (("simulation", "S1_simulation"), ("plasmode", "S1_plasmode"))
    }
    criteria = {
        "S1": {"holds": all(family["holds"] for family in s1_families.values()), "families": s1_families},
        "S2": judge_losses([_paired_row(paired, *spec[:2], arm, spec[2]) for spec in plan["S2"]]),
        "S3": grade.judge_rates([_observation(rates, *spec[:2], arm, spec[2]) for spec in plan["S3"]]),
        "G": judge_gain([_paired_row(paired, *spec[:2], arm, spec[2]) for spec in plan["G"]]),
    }
    failing_readings = {
        "S1": [
            {"family": family, **row}
            for family, judged in s1_families.items()
            for row in judged["observations"]
            if row.get("status") != "holds"
        ],
        "S2": [row for row in criteria["S2"]["observations"] if row.get("status") != "holds"],
        "S3": [row for row in criteria["S3"]["observations"] if row.get("status") != "holds"],
        "G": [] if criteria["G"]["holds"] else criteria["G"]["observations"],
    }
    failing = [name for name in ("S1", "S2", "S3", "G") if not criteria[name]["holds"]]
    return {
        "verdict": _VERDICT_OF[failing[0]] if failing else "qualifies",
        "failing_criteria": failing,
        "failing_readings": {name: rows for name, rows in failing_readings.items() if name in failing},
        "criteria": criteria,
    }


def decision(
    simulation: Mapping[str, Sequence[Mapping[str, Any]]],
    plasmode: Mapping[str, Sequence[Mapping[str, Any]]],
    *,
    reproduction: Mapping[str, Mapping[str, Any]],
) -> dict[str, Any]:
    """The predeclared verdict (rule version 1) from both profiles' rate, paired and cell tables."""

    rates = _index([*simulation["rates"], *plasmode["rates"]], "design", "cell", "arm", "statistic")
    paired = _index([*simulation["paired"], *plasmode["paired"]], "design", "cell", "statistic", "arm")
    plan = criteria_plan()
    gate_failures = [
        {"profile": profile, "reason": "reproduction failed"}
        for profile in PROFILES
        if not (reproduction.get(profile) or {}).get("passed", False)
    ]
    gate = {"passed": not gate_failures, "failures": gate_failures}
    arms = {arm: judge_arm(arm, plan, rates, paired) for arm in CANDIDATE_ARMS}
    qualifying = [arm for arm in CANDIDATE_ARMS if arms[arm]["verdict"] == "qualifies"]
    if not gate["passed"]:
        verdict, named = "reproduction_failed", None
    elif qualifying:
        verdict, named = "propose_paper_grade", qualifying[0]
    else:
        verdict, named = "keep_production", None
    return {
        "rule_version": RULE_VERSION,
        "verdict": verdict,
        "named_arm": named,
        "qualifying_arms": qualifying,
        "consequence": CONSEQUENCES[verdict],
        "gate": gate,
        "reference_arm": REFERENCE_ARM,
        "preference": list(CANDIDATE_ARMS),
        "alpha": ALPHA,
        "rejection": dict(REJECTION_RULES),
        "plan": {name: [list(spec) for spec in specs] for name, specs in plan.items()},
        "arms": arms,
        "descriptive": descriptive_tables(simulation, plasmode, plan),
        "notes": [
            DELTA_METRIC_NOTE,
            "S1 and S3 hold each candidate to its own rates under the stage-span paper grade's count-plus-ceiling "
            "rule; S2 and G compare it with pls_cv on the same units, data and draws (exact paired tests).",
            "The e = 1.0 cells' off-target rates are descriptive: the stage-span pilot (2026-11-04) measured "
            "production's own construction leak there.",
        ],
    }


def descriptive_tables(
    simulation: Mapping[str, Sequence[Mapping[str, Any]]],
    plasmode: Mapping[str, Sequence[Mapping[str, Any]]],
    plan: Mapping[str, Sequence[tuple[str, str, str]]],
) -> dict[str, Any]:
    """Paired differences, rank and fraction, AUROC and capture, the e = 1.0 off-target rates; never the verdict."""

    rates = _index([*simulation["rates"], *plasmode["rates"]], "design", "cell", "arm", "statistic")
    roles: dict[tuple[str, ...], str] = {}
    for name, specs in plan.items():
        for spec in specs:
            roles.setdefault(tuple(spec), name)
    differences = [
        {
            **{key: row.get(key) for key in PAIRED_COLUMNS},
            "role": roles.get((str(row["design"]), str(row["cell"]), str(row["statistic"]))),
        }
        for row in [*simulation["paired"], *plasmode["paired"]]
    ]
    distributions = [
        {
            key: row.get(key)
            for key in ("design", "cell", "arm", "n", *(f"{q}_{s}" for q in _QUANTITIES for s in ("median", "q10",
                                                                                                  "q90")),
                        "group_difference_median", "population_target_median")
        }
        for row in [*simulation["cells"], *plasmode["cells"]]
    ]
    large = []
    for design, cell, statistic in plan["S3_descriptive"]:
        for arm in ARMS:
            row = rates.get((design, cell, arm, statistic))
            if row is None:
                continue
            count = int(row["count"])
            large.append({
                "design": design, "cell": cell, "arm": arm, "statistic": statistic, "count": count,
                "rejections": int(row["rejections"]), "rate": _number(row["rate"]),
                "bound": None if count <= 0 else ALPHA + SE_MULTIPLIER * math.sqrt(ALPHA * (1 - ALPHA) / count),
            })
    support = [
        {key: row.get(key) for key in ("cell", "arm", "n", "keep_fraction_median",
                                        *(f"support_{block}_median" for block in PLASMODE_BLOCKS))}
        for row in plasmode["cells"]
    ]
    return {
        "paired_differences": differences,
        "rank_fraction_auroc_capture": distributions,
        "large_plant_off_target": large,
        "plasmode_support": support,
        "delta_metric_note": DELTA_METRIC_NOTE,
    }


_TABLE_FILES: dict[str, str] = {"rates": "rate_table.csv", "paired": "paired_table.csv", "cells": "cell_table.csv"}


def _read_json(path: Path) -> dict[str, Any]:
    if not path.exists():
        raise SystemExit(f"{path} does not exist; run `merge` first.")
    return json.loads(path.read_text(encoding="utf-8"))


def _profile_tables(directory: Path) -> dict[str, list[dict[str, Any]]]:
    tables = {}
    for name, filename in _TABLE_FILES.items():
        path = directory / filename
        if not path.exists():
            raise SystemExit(f"{path} does not exist; run `merge` first.")
        tables[name] = lsc._read_csv(path)
    return tables


def decide(out_dir: Path) -> dict[str, Any]:
    """Read both merged profiles and write ``production_latent_space_decision.json``."""

    loaded, reproduction, complete, inputs = {}, {}, {}, {}
    for profile in PROFILES:
        directory = profile_dir(out_dir, profile)
        settings = lsc.load_settings(directory)
        verify_parents(settings)
        summary = _read_json(directory / "accounting.json")
        if summary.get("settings_hash") != lsc.settings_hash(settings):
            raise SystemExit(f"{directory}/accounting.json was merged under other settings; run `merge` again.")
        loaded[profile] = _profile_tables(directory)
        reproduction[profile] = _read_json(directory / "reproduction.json")
        complete[profile] = {
            key: summary[key]
            for key in ("expected_units", "complete_units", "missing_units", "failed_records", "fidelity_failures")
        }
        for name in ("settings.json", *(_TABLE_FILES.values()), "accounting.json", "reproduction.json"):
            inputs[f"{profile}/{name}"] = sha256_file(directory / name)
    result = decision(loaded["simulation"], loaded["plasmode"], reproduction=reproduction)
    result["complete"] = complete
    result["reproduction"] = {
        profile: {key: reproduction[profile].get(key) for key in ("passed", "rows_compared", "rows_matching",
                                                                   "units_checked", "units_passed")}
        for profile in PROFILES
    }
    result["inputs"] = inputs
    lsc._write_json(Path(out_dir) / DECISION_NAME, result)
    return result


# ── CLI ──────────────────────────────────────────────────────────────────────


def _ladder(text: str) -> tuple[float, ...]:
    try:
        return tuple(float(value) for value in text.split(",") if value.strip())
    except ValueError as exc:
        raise argparse.ArgumentTypeError(f"--ladder must be comma-separated fractions, got {text!r}.") from exc


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
                     help=f"Permutations (simulation {SIMULATION_PERMUTATIONS}, plasmode {PLASMODE_PERMUTATIONS})")
    run.add_argument("--reps", type=int, default=None,
                     help=f"Replicates per cell (simulation {SIMULATION_REPLICATES}, plasmode {PLASMODE_REPLICATES})")
    run.add_argument("--plant-reps", type=int, default=None, help="Plasmode: replicates of the plant cells")
    run.add_argument("--ladder", type=_ladder, default=None, help="Sparse arms' keep fractions (default the package's)")
    run.add_argument("--designs", type=str, default=",".join(lsc.DESIGNS), help="Simulation: comma-separated")
    run.add_argument("--parent-dir", type=Path, default=None, help="The committed parent run (default: the profile's)")
    run.add_argument("--export-dir", type=Path, default=None, help="Plasmode: the SEA-AD cohort_30 export")
    merge_parser = sub.add_parser("merge", help="Merge one profile's shards")
    merge_parser.add_argument("--profile", choices=PROFILES, required=True)
    merge_parser.add_argument("--out-dir", type=Path, required=True)
    decide_parser = sub.add_parser("decide", help=f"Write {DECISION_NAME} from both profiles")
    decide_parser.add_argument("--out-dir", type=Path, required=True)
    return parser


def settings_from_args(args: argparse.Namespace) -> tuple[dict[str, Any], Any]:
    if args.profile == "simulation":
        settings = simulation_settings(
            parent_dir=args.parent_dir or COMPARISON_DIR,
            replicates=SIMULATION_REPLICATES if args.reps is None else args.reps,
            permutations=SIMULATION_PERMUTATIONS if args.perms is None else args.perms,
            designs=[name for name in args.designs.split(",") if name],
            ladder=args.ladder,
        )
        return settings, None
    if args.export_dir is None:
        raise SystemExit("The plasmode profile needs --export-dir.")
    try:
        export = pc.load_export(args.export_dir)
    except pc.ExportContractError as exc:
        raise pc.contract_exit(exc) from exc
    settings = plasmode_settings(
        export,
        parent_dir=args.parent_dir or PLASMODE_DIR,
        replicates=PLASMODE_REPLICATES if args.reps is None else args.reps,
        plant_replicates=args.plant_reps,
        permutations=PLASMODE_PERMUTATIONS if args.perms is None else args.perms,
        ladder=args.ladder,
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
                      "complete": result["complete"]}, indent=1))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
