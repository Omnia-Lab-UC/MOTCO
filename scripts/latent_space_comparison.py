#!/usr/bin/env python3
"""Latent-space comparison: five measurement spaces across alignment with the stage span (no R).

The driver behind ``openspec/changes/archive/2026-11-04-latent-space-comparison/``. The stage-span
pilot (2026-11-04) found that the production stage-supervised PLS space sees a
group difference only through its component inside the fitted span. Ranking
spaces on one construction is circular (a plant inside the span favours
stage-PLS, one outside favours a label-free space), so this comparison varies
the unknown: the **alignment** φ of a stage-span orientation with the stage
span (``span_alignment``), at a fixed intended geometry, and reads every space's
visibility, size and protected cost on identical data.

**Arms** (every unit generates one dataset; all arms share one RRPP seed):

- ``pls_cv`` — production: stage PLS-DA, double-CV rank with floor 2;
- ``pca_matched`` — label-free PCA at that unit's ``pls_cv`` rank;
- ``pca_fixed`` — PCA at rank 2(K − 1), the orientation probe's arm;
- ``pls_residual_pca`` — the ``pls_cv`` scores plus K − 1 stage-residual
  components orthogonal to the PLS span, built on the same fit;
- ``concat`` — the standardized features, the reference ceiling.

The PLS model is fitted once per unit and shared by ``pls_cv``,
``pca_matched``'s rank and ``pls_residual_pca``. Every arm runs with net
displacement and span capture on.

**Designs** (three stages, methylation + expression, ``p_dmp`` 0.1, ρ = 0,
joint magnitude, span layers = the measured blocks, ``construction =
stage_span``): ``n200`` (n = 200, group ratio 0.5, δ = 2.0), ``seaad`` (the
SEA-AD table ``[[11, 10, 28], [9, 10, 12]]``, δ = 2.0) and ``seaad_weak`` (the
table at the calibrated δ). **Cells** (14 per design): the anchor, translation
e = 1.0, magnitude e = 0.05, orientation e ∈ {0.25, 0.5} × φ ∈ {0, 30, 60, 90}
(dense support), shape e ∈ {0.25, 0.5}, and orientation e = 0.5 at φ = 90 on a
20-feature support.

**Seeds** come from base seed 2035: the generator seed hashes (design,
replicate), so every cell at a replicate shares its null dataset; the RRPP
seed hashes (design, cell, replicate); calibration draws from its own
``calibration`` stream.

Subcommands:

- ``calibrate`` — the weak design's δ (design D9): the ladder
  {0.1, 0.11, 0.12, 0.13, 0.14, 0.15, 0.2, 0.25, 0.3, 0.4, 0.5} (0.11–0.14
  added on 2026-11-04, before any run unit, when the band fell between 0.1
  and 0.15), all three per-omic δ together, 40
  anchor replicates of the SEA-AD table per step, ``pls_cv`` fits only. The
  first step whose median held-out stage AUROC lies in [0.60, 0.67] is
  adopted, else the step closest to 0.636, marked ``band_missed``. Writes
  ``calibration.json`` (every step, the adoption) and
  ``calibration_units.csv``.
- ``run`` — evaluate (a shard of) the units into
  ``records/shard_<i>_of_<N>.jsonl``, one line per (unit, arm). It refuses to
  start without ``--calibration`` or an explicit ``--weak-delta``. The
  settings and their hash are written to ``settings.json``; a resume or merge
  into a directory written with different settings exits naming the keys.
  Resume skips complete units.
- ``merge`` — ``merged.jsonl`` (per-unit records; keep off the repository),
  ``rate_table.csv`` (rate and Monte Carlo SE per arm × design × cell ×
  statistic), ``capture_table.csv`` (span-capture shares, retained rank and
  the ``population_latent`` target per arm × design × cell) and
  ``accounting.json`` (units, failures by arm and cell, duplicates,
  recorded core-hours, per-arm timing).
- ``decide`` — ``latent_space_comparison_reading.json``, the predeclared
  advisory reading (:func:`reading`): ``pca_matched``, ``pca_fixed`` and
  ``pls_residual_pca`` against ``pls_cv`` on size (C1), out-of-span gain (C2)
  and protected power (C3); ``propose_paper_grade`` when an arm qualifies,
  else ``document_tradeoff``. ``concat`` is never eligible.

Examples:

    # calibration (local, pinned BLAS)
    python scripts/latent_space_comparison.py calibrate --out-dir /tmp/lsc/calibration --workers 7

    # smoke with an explicit δ
    python scripts/latent_space_comparison.py run --weak-delta 0.15 --reps 2 --perms 19 \\
        --out-dir /tmp/lsc/smoke --workers 7
    python scripts/latent_space_comparison.py merge --out-dir /tmp/lsc/smoke

    # one SLURM array task (scripts/latent_space_comparison_array.sbatch), then
    python scripts/latent_space_comparison.py merge --out-dir /abs/run
    python scripts/latent_space_comparison.py decide --out-dir /abs/run
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import sys
from concurrent.futures import ProcessPoolExecutor, as_completed
from dataclasses import dataclass
from pathlib import Path
from statistics import median
from time import perf_counter
from typing import Any, Iterable, Mapping, Sequence

import numpy as np

DRIVER_VERSION = 1
RULE_VERSION = 1
BASE_SEED = 2035
DEFAULT_REPLICATES = 100
DEFAULT_PERMUTATIONS = 199
DEFAULT_ALPHA = 0.05

ARMS: tuple[str, ...] = ("pls_cv", "pca_matched", "pca_fixed", "pls_residual_pca", "concat")
REFERENCE_ARM = "pls_cv"
CEILING_ARM = "concat"
ELIGIBLE_ARMS: tuple[str, ...] = ("pca_matched", "pca_fixed", "pls_residual_pca")
STATISTICS: tuple[str, ...] = ("delta", "angle", "shape")

LAYERS: list[str] = ["methylation", "expression"]
SEAAD_TABLE: list[list[int]] = [[11, 10, 28], [9, 10, 12]]
DESIGNS: tuple[str, ...] = ("n200", "seaad", "seaad_weak")
#: Designs at the generator's stage effect, where C2 and C3 apply.
STRONG_DESIGNS: tuple[str, ...] = ("n200", "seaad")

#: Stage-PLS cross-validation knobs of the study configs (production rank rule).
PLS_CV_PARAMS: dict[str, Any] = {
    "cv1_splits": 3,
    "cv2_splits": 4,
    "n_repeats": 5,
    "max_components": 20,
    "random_state": 1203,
}
MIN_COMPONENTS = 2
ORIENTATION_ESTIMATOR = "net_displacement"

#: (cell, trajectory mode, effect, φ or None, support or None) — design D10.
CELLS: tuple[tuple[str, str, float, float | None, int | None], ...] = (
    ("anchor", "none", 0.0, None, None),
    ("translation", "translation", 1.0, None, None),
    ("magnitude", "magnitude", 0.05, None, None),
    *(
        (f"orient_e{effect:g}_phi{phi}", "orientation", effect, float(phi), None)
        for effect in (0.25, 0.5)
        for phi in (0, 30, 60, 90)
    ),
    ("shape_e0.25", "shape", 0.25, None, None),
    ("shape_e0.5", "shape", 0.5, None, None),
    ("orient_e0.5_phi90_k20", "orientation", 0.5, 90.0, 20),
)

#: The statistic each mode targets.
TARGET_STATISTIC: dict[str, str | None] = {
    "none": None,
    "translation": None,
    "magnitude": "delta",
    "orientation": "angle",
    "shape": "shape",
}

#: Calibration (design D9).
#: Amended 2026-11-04 before any run unit: the original ladder (0.1, 0.15, 0.2, …)
#: put the band between its first two steps (job 933435), so 0.11–0.14 were added.
CALIBRATION_LADDER: tuple[float, ...] = (0.1, 0.11, 0.12, 0.13, 0.14, 0.15, 0.2, 0.25, 0.3, 0.4, 0.5)
CALIBRATION_REPLICATES = 40
CALIBRATION_BAND: tuple[float, float] = (0.60, 0.67)
CALIBRATION_TARGET = 0.636  # SEA-AD astrocytes cohort_30, held-out stage AUROC


# ── Designs, settings and seeds ──────────────────────────────────────────────


def design_generators(weak_delta: float) -> dict[str, dict[str, Any]]:
    """Generator parameters per design (mode, effect, alignment and seed are set per unit)."""

    common: dict[str, Any] = {
        "n_stages": 3,
        "p_dmp": 0.1,
        "baseline_continuity": 0.0,
        "magnitude_kind": "joint",
        "construction": "stage_span",
        "span_layers": list(LAYERS),
    }
    weak = {"delta_methyl": weak_delta, "delta_expr": weak_delta, "delta_protein": weak_delta}
    return {
        "n200": {**common, "n_samples": 200, "group_ratio": 0.5},
        "seaad": {**common, "group_stage_sizes": SEAAD_TABLE},
        "seaad_weak": {**common, "group_stage_sizes": SEAAD_TABLE, **weak},
    }


def build_settings(
    *,
    weak_delta: float | None = None,
    calibration: Path | None = None,
    replicates: int = DEFAULT_REPLICATES,
    permutations: int = DEFAULT_PERMUTATIONS,
    alpha: float = DEFAULT_ALPHA,
    base_seed: int = BASE_SEED,
    designs: Sequence[str] = DESIGNS,
) -> dict[str, Any]:
    """The run's settings; the weak design's δ comes from a calibration or is explicit."""

    if (weak_delta is None) == (calibration is None):
        raise SystemExit(
            "run needs the weak design's δ: pass --calibration <calibration.json> from `calibrate`, "
            "or an explicit --weak-delta (smoke and timing)."
            if weak_delta is None
            else "Pass --calibration or --weak-delta, not both."
        )
    if calibration is not None:
        adopted = load_calibration(calibration)
        weak = {
            "value": float(adopted["adopted"]["delta"]),
            "source": "calibration",
            "rule": adopted["adopted"]["rule"],
            "auroc_median": adopted["adopted"]["auroc_median"],
            "calibration_sha256": hashlib.sha256(calibration.read_bytes()).hexdigest(),
        }
    else:
        assert weak_delta is not None
        weak = {"value": float(weak_delta), "source": "explicit"}
    unknown = sorted(set(designs) - set(DESIGNS))
    if unknown:
        raise SystemExit(f"Unknown design(s) {unknown}; expected a subset of {list(DESIGNS)}.")
    generators = design_generators(weak["value"])
    return {
        "driver_version": DRIVER_VERSION,
        "designs": {
            name: {"generator": generators[name], "layers": list(LAYERS)} for name in DESIGNS if name in designs
        },
        "weak_delta": weak,
        "cells": [list(cell) for cell in CELLS],
        "arms": list(ARMS),
        "pls_cv_params": dict(PLS_CV_PARAMS),
        "min_components": MIN_COMPONENTS,
        "orientation_estimator": ORIENTATION_ESTIMATOR,
        "replicates": int(replicates),
        "permutations": int(permutations),
        "alpha": float(alpha),
        "base_seed": int(base_seed),
    }


def settings_hash(settings: Mapping[str, Any]) -> str:
    encoded = json.dumps(settings, sort_keys=True, separators=(",", ":")).encode("utf-8")
    return hashlib.sha256(encoded).hexdigest()


def _flatten(value: Any, prefix: str = "") -> dict[str, Any]:
    if isinstance(value, Mapping):
        out: dict[str, Any] = {}
        for key in value:
            out.update(_flatten(value[key], f"{prefix}.{key}" if prefix else str(key)))
        return out
    return {prefix: value}


def settings_difference(stored: Mapping[str, Any], current: Mapping[str, Any]) -> list[str]:
    left, right = _flatten(stored), _flatten(current)
    return sorted(key for key in set(left) | set(right) if left.get(key, "<absent>") != right.get(key, "<absent>"))


def ensure_settings(out_dir: Path, settings: Mapping[str, Any]) -> None:
    """Write the settings header, or refuse a directory written with different settings."""

    path = out_dir / "settings.json"
    if path.exists():
        stored = json.loads(path.read_text(encoding="utf-8"))
        if stored.get("hash") != settings_hash(settings):
            differing = settings_difference(stored.get("settings", {}), settings)
            raise SystemExit(
                f"{path} was written with different settings; differing: {', '.join(differing)}. "
                "Use a fresh --out-dir."
            )
        return
    out_dir.mkdir(parents=True, exist_ok=True)
    payload = {"hash": settings_hash(settings), "settings": settings}
    path.write_text(json.dumps(payload, indent=1, sort_keys=True) + "\n", encoding="utf-8")


def load_settings(out_dir: Path) -> dict[str, Any]:
    path = out_dir / "settings.json"
    if not path.exists():
        raise SystemExit(f"{path} does not exist; run `run` first.")
    stored = json.loads(path.read_text(encoding="utf-8"))
    settings = stored["settings"]
    if stored.get("hash") != settings_hash(settings):
        raise SystemExit(f"{path} is corrupt: its hash does not match its settings.")
    return settings


def _seed(*parts: Any) -> int:
    encoded = json.dumps(["latent-space-comparison", *parts], separators=(",", ":")).encode("utf-8")
    return int(hashlib.sha256(encoded).hexdigest()[:8], 16) & 0x7FFFFFFF


def generator_seed(base_seed: int, design: str, replicate: int) -> int:
    """Matched across cells: every cell of a (design, replicate) shares its null dataset."""

    return _seed(base_seed, design, "generator", replicate)


def rrpp_seed(base_seed: int, design: str, cell: str, replicate: int) -> int:
    """One RRPP seed per unit, shared by every arm."""

    return _seed(base_seed, design, cell, "rrpp", replicate)


def calibration_seed(base_seed: int, replicate: int) -> int:
    """The calibration stream: the same anchor dataset noise at every δ step."""

    return _seed(base_seed, "calibration", "generator", replicate)


@dataclass(frozen=True)
class Unit:
    design: str
    cell: str
    replicate: int

    @property
    def key(self) -> tuple[str, str, int]:
        return (self.design, self.cell, self.replicate)


def enumerate_units(settings: Mapping[str, Any]) -> list[Unit]:
    return [
        Unit(design, cell[0], replicate)
        for design in settings["designs"]
        for cell in settings["cells"]
        for replicate in range(int(settings["replicates"]))
    ]


def shard_units(units: Sequence[Unit], index: int, count: int) -> list[Unit]:
    """Round-robin, so every shard mixes designs and cells."""

    return [unit for position, unit in enumerate(units) if position % count == index]


# ── Unit evaluation ──────────────────────────────────────────────────────────


def _generator_params(generator: Mapping[str, Any], **overrides: Any) -> Any:
    from motco.simulations.semisynthetic import SemiSyntheticTrajectoryParams

    values = dict(generator)
    if values.get("group_stage_sizes") is not None:
        values["group_stage_sizes"] = tuple(tuple(int(n) for n in row) for row in values["group_stage_sizes"])
    if values.get("span_layers") is not None:
        values["span_layers"] = tuple(values["span_layers"])
    return SemiSyntheticTrajectoryParams(**{**values, **overrides})


def unit_dataset(unit: Unit, settings: Mapping[str, Any], reference: Any = None) -> Any:
    """The unit's dataset: the (design, replicate) null dataset plus the cell's plant."""

    from motco.simulations.reference import load_reference
    from motco.simulations.semisynthetic import generate_semisynthetic_trajectory

    cell = cell_spec(settings, unit.cell)
    params = _generator_params(
        settings["designs"][unit.design]["generator"],
        seed=generator_seed(settings["base_seed"], unit.design, unit.replicate),
        trajectory_mode=cell["mode"],
        group_effect_size=float(cell["effect"]),
        span_alignment=cell["phi"],
        span_support=cell["support"],
    )
    return generate_semisynthetic_trajectory(params, reference=reference if reference is not None else load_reference())


def cell_spec(settings: Mapping[str, Any], name: str) -> dict[str, Any]:
    for cell, mode, effect, phi, support in settings["cells"]:
        if cell == name:
            return {"cell": cell, "mode": mode, "effect": effect, "phi": phi, "support": support}
    raise KeyError(f"Unknown cell {name!r}.")


def arm_params(arm: str, settings: Mapping[str, Any], *, seed: int, n_stages: int, pls_rank: int | None = None) -> Any:
    """The ``SimulationEvaluationParams`` of one arm."""

    from motco.simulations.evaluation import SimulationEvaluationParams

    layers = {"layers": list(LAYERS)}
    if arm == "pls_cv":
        method, integration = "pls", {**settings["pls_cv_params"], **layers}
    elif arm == "pls_residual_pca":
        method, integration = "pls_residual_pca", {**settings["pls_cv_params"], **layers}
    elif arm == "pca_matched":
        if pls_rank is None:
            raise ValueError("pca_matched needs the unit's pls_cv rank.")
        method, integration = "pca", {"n_components": int(pls_rank), **layers}
    elif arm == "pca_fixed":
        method, integration = "pca", {"n_components": 2 * (n_stages - 1), **layers}
    elif arm == "concat":
        method, integration = "concat", dict(layers)
    else:
        raise ValueError(f"Unknown arm {arm!r}.")
    return SimulationEvaluationParams(
        integration_method=method,
        integration_params=integration,
        permutations=int(settings["permutations"]),
        n_jobs=1,
        seed=seed,
        orientation_estimator=settings["orientation_estimator"],
        min_components=int(settings["min_components"]),
        span_capture=True,
    )


def _finite(values: Mapping[str, Any]) -> dict[str, float | None]:
    out: dict[str, float | None] = {}
    for key, value in values.items():
        out[key] = float(value) if value is not None and math.isfinite(float(value)) else None
    return out


def _retained(metadata: Mapping[str, Any]) -> int:
    for key in ("selected_lv", "retained_dimensionality", "n_features"):
        if metadata.get(key) is not None:
            return int(metadata[key])
    raise KeyError("integration metadata records no dimensionality")


def evaluate_unit(unit: Unit, settings: Mapping[str, Any]) -> list[dict[str, Any]]:
    """Generate one dataset, fit stage-PLS once, and evaluate every arm on it; one record per arm."""

    from motco.simulations.diagnostics import POPULATION_LATENT_CHECKPOINT
    from motco.simulations.evaluation import (
        evaluate_semisynthetic_trajectory,
        integrate_semisynthetic_dataset,
        pls_residual_pca_from_pls,
    )

    unit_start = perf_counter()
    cell = cell_spec(settings, unit.cell)
    dataset = unit_dataset(unit, settings)
    n_stages = int(settings["designs"][unit.design]["generator"]["n_stages"])
    perm_seed = rrpp_seed(settings["base_seed"], unit.design, unit.cell, unit.replicate)
    stage_span = dataset.truth.get("stage_span") or {}
    base = {
        "design": unit.design,
        "cell": unit.cell,
        "trajectory_mode": cell["mode"],
        "effect_size": float(cell["effect"]),
        "span_alignment": cell["phi"],
        "span_support": cell["support"],
        "replicate": unit.replicate,
        "generator_seed": generator_seed(settings["base_seed"], unit.design, unit.replicate),
        "rrpp_seed": perm_seed,
        "settings_hash": settings_hash(settings),
        "intended_geometry": stage_span.get("intended_geometry"),
        "plant_norm": stage_span.get("plant_norm"),
        "out_of_span": stage_span.get("out_of_span"),
    }

    fit_start = perf_counter()
    pls_latent = None
    pls_error: Exception | None = None
    try:
        pls_latent = integrate_semisynthetic_dataset(dataset, arm_params("pls_cv", settings, seed=perm_seed,
                                                                          n_stages=n_stages))
    except Exception as exc:  # recorded per dependent arm, never aborts the unit
        pls_error = exc
    pls_fit_seconds = perf_counter() - fit_start
    pls_rank = None if pls_latent is None else int(pls_latent.metadata["selected_lv"])

    records = []
    for arm in settings["arms"]:
        arm_start = perf_counter()
        if arm in ("pls_cv", "pca_matched", "pls_residual_pca") and pls_latent is None:
            records.append(_failed(base, arm, pls_error or RuntimeError("no pls fit")))
            continue
        try:
            params = arm_params(arm, settings, seed=perm_seed, n_stages=n_stages, pls_rank=pls_rank)
            latent = None
            if arm == "pls_cv":
                latent = pls_latent
            elif arm == "pls_residual_pca":
                latent = pls_residual_pca_from_pls(dataset, pls_latent, params)
            result = evaluate_semisynthetic_trajectory(dataset, params, latent=latent)
        except Exception as exc:
            records.append(_failed(base, arm, exc))
            continue
        geometry = result.realized_geometry
        metadata = result.latent_matrix_metadata
        records.append(
            {
                **base,
                "arm": arm,
                "status": "completed",
                "integration_method": params.integration_method,
                "retained_dimensionality": _retained(metadata),
                "pls_rank": pls_rank,
                "cv_mean_auroc": metadata.get("cv_mean_auroc"),
                "pair_statistics": _finite(result.pair_statistics),
                "p_values": _finite(result.p_values),
                "span_capture": geometry.get("span_capture"),
                "population_latent": (geometry["checkpoints"].get(POPULATION_LATENT_CHECKPOINT) or {}).get("joint"),
                "runtime_seconds": {
                    "pls_fit": pls_fit_seconds if arm in ("pls_cv", "pls_residual_pca") else 0.0,
                    "arm": perf_counter() - arm_start,
                },
            }
        )
    unit_seconds = perf_counter() - unit_start
    for record in records:
        record["unit_runtime_seconds"] = unit_seconds
    return records


def _failed(base: Mapping[str, Any], arm: str, exc: Exception) -> dict[str, Any]:
    return {
        **base,
        "arm": arm,
        "status": "failed",
        "error_type": type(exc).__name__,
        "error_message": str(exc),
    }


# ── Persistence, resume and the run ──────────────────────────────────────────


def _jsonable(value: Any) -> Any:
    if isinstance(value, Mapping):
        return {str(key): _jsonable(item) for key, item in value.items()}
    if isinstance(value, list | tuple):
        return [_jsonable(item) for item in value]
    if isinstance(value, np.generic):
        return _jsonable(value.item())
    if isinstance(value, float) and not math.isfinite(value):
        return None
    return value


def _line(record: Mapping[str, Any]) -> str:
    return json.dumps(_jsonable(record), sort_keys=True, allow_nan=False) + "\n"


def shard_path(out_dir: Path, index: int, count: int) -> Path:
    return out_dir / "records" / f"shard_{index:03d}_of_{count:03d}.jsonl"


def read_records(paths: Iterable[Path]) -> list[dict[str, Any]]:
    records = []
    for path in paths:
        with path.open(encoding="utf-8") as handle:
            for line in handle:
                try:
                    records.append(json.loads(line))
                except json.JSONDecodeError:
                    continue  # a line torn by an interrupted write; its unit is re-run
    return records


def _unit_key(record: Mapping[str, Any]) -> tuple[str, str, int]:
    return (str(record["design"]), str(record["cell"]), int(record["replicate"]))


def complete_units(records: Iterable[Mapping[str, Any]], arms: Sequence[str]) -> dict[tuple[str, str, int], list]:
    """Records grouped by unit, keeping only units with every arm present (completed or failed)."""

    grouped: dict[tuple[str, str, int], dict[str, Mapping[str, Any]]] = {}
    for record in records:
        grouped.setdefault(_unit_key(record), {})[str(record["arm"])] = record
    return {key: [by_arm[a] for a in arms] for key, by_arm in grouped.items() if set(by_arm) == set(arms)}


def run_shard(
    out_dir: Path,
    settings: Mapping[str, Any],
    index: int,
    count: int,
    *,
    limit: int | None = None,
    workers: int = 1,
) -> int:
    """Evaluate this shard's missing units; returns how many were evaluated."""

    path = shard_path(out_dir, index, count)
    path.parent.mkdir(parents=True, exist_ok=True)
    done = complete_units(read_records([path]) if path.exists() else [], settings["arms"])
    # Rewrite the shard with complete units only, so a torn or partial unit
    # left by an interruption is re-run instead of duplicated.
    with path.open("w", encoding="utf-8") as handle:
        for key in sorted(done):
            handle.write("".join(_line(record) for record in done[key]))
    pending = [unit for unit in shard_units(enumerate_units(settings), index, count) if unit.key not in done]
    if limit is not None:
        pending = pending[:limit]

    def write(position: int, unit: Unit, records: Sequence[Mapping[str, Any]]) -> None:
        with path.open("a", encoding="utf-8") as handle:
            handle.write("".join(_line(record) for record in records))
            handle.flush()
        print(f"[shard {index}/{count}] {position}/{len(pending)} {unit.key}", file=sys.stderr, flush=True)

    if workers <= 1:
        for position, unit in enumerate(pending, start=1):
            write(position, unit, evaluate_unit(unit, settings))
    else:
        with ProcessPoolExecutor(max_workers=workers) as pool:
            futures = {pool.submit(evaluate_unit, unit, settings): unit for unit in pending}
            for position, future in enumerate(as_completed(futures), start=1):
                write(position, futures[future], future.result())
    return len(pending)


# ── Merge: rate table, capture table, accounting ─────────────────────────────


RATE_COLUMNS = ("design", "cell", "arm", "statistic", "rate", "se", "count", "rejections", "failures")


def rate_table(records: Iterable[Mapping[str, Any]], *, alpha: float) -> list[dict[str, Any]]:
    """Rejection rate (p < alpha), binomial SE and count per design × cell × arm × statistic."""

    tallies: dict[tuple[str, str, str, str], list[int]] = {}
    failures: dict[tuple[str, str, str], int] = {}
    for record in records:
        key3 = (str(record["design"]), str(record["cell"]), str(record["arm"]))
        if record.get("status") != "completed":
            failures[key3] = failures.get(key3, 0) + 1
            for statistic in STATISTICS:
                tallies.setdefault((*key3, statistic), [0, 0])
            continue
        for statistic in STATISTICS:
            tally = tallies.setdefault((*key3, statistic), [0, 0])
            p_value = (record.get("p_values") or {}).get(statistic)
            if p_value is None or not math.isfinite(float(p_value)):
                continue
            tally[0] += 1
            tally[1] += int(float(p_value) < alpha)
    rows = []
    for (design, cell, arm, statistic), (count, rejections) in sorted(tallies.items()):
        rate = rejections / count if count else None
        rows.append(
            {
                "design": design,
                "cell": cell,
                "arm": arm,
                "statistic": statistic,
                "rate": rate,
                "se": None if rate is None else math.sqrt(rate * (1.0 - rate) / count),
                "count": count,
                "rejections": rejections,
                "failures": failures.get((design, cell, arm), 0),
            }
        )
    return rows


CAPTURE_COLUMNS = (
    "design",
    "cell",
    "arm",
    "trajectory_mode",
    "effect_size",
    "span_alignment",
    "span_support",
    "n",
    "group_difference_mean",
    "group_difference_sd",
    "group_difference_q10",
    "group_difference_median",
    "group_difference_q90",
    "pooled_stage_median",
    "retained_rank_median",
    "random_share_median",
    "target_statistic",
    "population_latent_target_median",
    "intended_target",
    "population_latent_delta_median",
    "population_latent_angle_median",
    "population_latent_shape_median",
)


def _quantile(values: Sequence[float], q: float) -> float | None:
    return float(np.quantile(values, q)) if values else None


def capture_table(records: Iterable[Mapping[str, Any]]) -> list[dict[str, Any]]:
    """Span capture, retained rank and the realized target per design × cell × arm (completed records)."""

    grouped: dict[tuple[str, str, str], list[Mapping[str, Any]]] = {}
    for record in records:
        if record.get("status") == "completed":
            grouped.setdefault((str(record["design"]), str(record["cell"]), str(record["arm"])), []).append(record)
    rows = []
    for (design, cell, arm), group in sorted(grouped.items()):
        first = group[0]
        target = TARGET_STATISTIC.get(str(first["trajectory_mode"]))
        captures = [r.get("span_capture") or {} for r in group]
        difference = [float(c["group_difference"]) for c in captures if c.get("group_difference") is not None]
        pooled = [float(c["pooled_stage"]) for c in captures if c.get("pooled_stage") is not None]
        random_share = [float(c["random_share"]) for c in captures if c.get("random_share") is not None]
        realized = [
            float(value)
            for r in group
            if target is not None
            and (value := (r.get("population_latent") or {}).get(target)) is not None
        ]
        intended = [
            float(value)
            for r in group
            if target is not None and (value := (r.get("intended_geometry") or {}).get(target)) is not None
        ]
        # The population means through each replicate's own scaler and basis, every statistic.
        population = {
            statistic: [
                float(value) for r in group if (value := (r.get("population_latent") or {}).get(statistic)) is not None
            ]
            for statistic in STATISTICS
        }
        rows.append(
            {
                "design": design,
                "cell": cell,
                "arm": arm,
                "trajectory_mode": first["trajectory_mode"],
                "effect_size": first["effect_size"],
                "span_alignment": first.get("span_alignment"),
                "span_support": first.get("span_support"),
                "n": len(group),
                "group_difference_mean": float(np.mean(difference)) if difference else None,
                "group_difference_sd": float(np.std(difference, ddof=1)) if len(difference) > 1 else None,
                "group_difference_q10": _quantile(difference, 0.1),
                "group_difference_median": _quantile(difference, 0.5),
                "group_difference_q90": _quantile(difference, 0.9),
                "pooled_stage_median": _quantile(pooled, 0.5),
                "retained_rank_median": float(median(int(r["retained_dimensionality"]) for r in group)),
                "random_share_median": _quantile(random_share, 0.5),
                "target_statistic": target,
                "population_latent_target_median": _quantile(realized, 0.5),
                "intended_target": _quantile(intended, 0.5),
                **{
                    f"population_latent_{statistic}_median": _quantile(values, 0.5)
                    for statistic, values in population.items()
                },
            }
        )
    return rows


def accounting(
    records: Sequence[Mapping[str, Any]],
    settings: Mapping[str, Any],
    *,
    duplicates: int,
    missing: Sequence[tuple[str, str, int]],
) -> dict[str, Any]:
    """Units, failures by arm and cell, duplicates, recorded core-hours and per-arm timing."""

    units: dict[tuple[str, str, int], float] = {}
    failures: dict[str, dict[str, int]] = {}
    arm_seconds: dict[str, list[float]] = {}
    for record in records:
        units[_unit_key(record)] = float(record.get("unit_runtime_seconds") or 0.0)
        if record.get("status") != "completed":
            by_cell = failures.setdefault(str(record["arm"]), {})
            cell = f"{record['design']}/{record['cell']}"
            by_cell[cell] = by_cell.get(cell, 0) + 1
        else:
            arm_seconds.setdefault(str(record["arm"]), []).append(float(record["runtime_seconds"]["arm"]))
    unit_seconds = list(units.values())
    return {
        "settings_hash": settings_hash(settings),
        "expected_units": len(enumerate_units(settings)),
        "complete_units": len(units),
        "missing_units": len(missing),
        "records": len(records),
        "failed_records": sum(sum(cells.values()) for cells in failures.values()),
        "failures_by_arm_and_cell": failures,
        "duplicate_records": duplicates,
        "recorded_core_hours": sum(unit_seconds) / 3600.0,
        "unit_seconds": {
            "median": _quantile(unit_seconds, 0.5),
            "q90": _quantile(unit_seconds, 0.9),
            "max": max(unit_seconds) if unit_seconds else None,
        },
        "arm_seconds_median": {arm: _quantile(values, 0.5) for arm, values in sorted(arm_seconds.items())},
    }


def _write_csv(path: Path, rows: Sequence[Mapping[str, Any]], columns: Sequence[str]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(columns), lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow({key: _csv_value(row.get(key)) for key in columns})


def _csv_value(value: Any) -> Any:
    if isinstance(value, float):
        return f"{value:.6g}"
    return "" if value is None else value


def _write_json(path: Path, payload: Any) -> None:
    path.write_text(json.dumps(_jsonable(payload), indent=1, sort_keys=True, allow_nan=False) + "\n", encoding="utf-8")


def merge(out_dir: Path, settings: Mapping[str, Any]) -> dict[str, Any]:
    """Merge every shard, then write the rate table, the capture table and the accounting."""

    ensure_settings(out_dir, settings)
    paths = sorted((out_dir / "records").glob("shard_*.jsonl"))
    raw = read_records(paths)
    for record in raw:
        if record.get("settings_hash") != settings_hash(settings):
            raise SystemExit(
                f"A record of unit {_unit_key(record)} was written with settings hash "
                f"{record.get('settings_hash')!r}, not {settings_hash(settings)!r}; refusing to merge."
            )
    seen: dict[tuple[tuple[str, str, int], str], int] = {}
    for record in raw:
        key = (_unit_key(record), str(record["arm"]))
        seen[key] = seen.get(key, 0) + 1
    duplicates = sum(count - 1 for count in seen.values())
    units = complete_units(raw, settings["arms"])
    expected = {unit.key for unit in enumerate_units(settings)}
    missing = sorted(expected - set(units))
    records = [record for key in sorted(units) if key in expected for record in units[key]]
    with (out_dir / "merged.jsonl").open("w", encoding="utf-8") as handle:
        handle.write("".join(_line(record) for record in records))
    rates = rate_table(records, alpha=float(settings["alpha"]))
    _write_csv(out_dir / "rate_table.csv", rates, RATE_COLUMNS)
    _write_csv(out_dir / "capture_table.csv", capture_table(records), CAPTURE_COLUMNS)
    summary = accounting(records, settings, duplicates=duplicates, missing=missing)
    _write_json(out_dir / "accounting.json", summary)
    return {
        key: summary[key]
        for key in ("complete_units", "missing_units", "failed_records", "duplicate_records", "recorded_core_hours")
    }


# ── The advisory reading (pure functions) ────────────────────────────────────


def size_bound(alpha: float, replicates: int) -> float:
    return alpha + 2.0 * math.sqrt(alpha * (1.0 - alpha) / replicates)


def _number(value: Any) -> float | None:
    if value is None or value == "":
        return None
    return float(value)


def _index(rows: Iterable[Mapping[str, Any]], *keys: str) -> dict[tuple[str, ...], Mapping[str, Any]]:
    return {tuple(str(row[key]) for key in keys): row for row in rows}


#: The reading's criteria: (criterion, design set, cell, statistic, kind).
def _criteria_plan(designs: Sequence[str]) -> dict[str, list[tuple[str, str, str, str]]]:
    strong = [design for design in STRONG_DESIGNS if design in designs]
    return {
        "size": [
            (design, cell, statistic, "size")
            for design in designs
            for cell in ("anchor", "translation")
            for statistic in STATISTICS
        ],
        "gain": [(design, "orient_e0.5_phi90", "angle", "gain") for design in strong],
        "protected": [
            (design, cell, statistic, "no_loss")
            for design in strong
            for cell, statistic in (("orient_e0.5_phi0", "angle"), ("shape_e0.5", "shape"), ("magnitude", "delta"))
        ],
    }


def reading(
    rates: Sequence[Mapping[str, Any]],
    captures: Sequence[Mapping[str, Any]],
    *,
    alpha: float,
    designs: Sequence[str] = DESIGNS,
) -> dict[str, Any]:
    """Apply the predeclared advisory reading (design D11) to the merged tables.

    For each eligible arm against ``pls_cv``:

    - **C1 size** — at every design, the anchor and translation rates of
      ``delta``, ``angle`` and ``shape`` at or below α + 2·√(α(1 − α)/R);
    - **C2 gain** — at both δ = 2.0 designs, orientation e = 0.5 at φ = 90
      (dense) ``angle`` above ``pls_cv``'s by more than 2 pooled SE;
    - **C3 protected** — at both δ = 2.0 designs, orientation e = 0.5 at φ = 0
      ``angle``, shape e = 0.5 ``shape`` and magnitude e = 0.05 ``delta`` at
      least ``pls_cv``'s minus 2 pooled SE.

    Verdicts: ``size_excess`` (C1 fails), else ``no_gain`` (C2 fails), else
    ``trade_off`` (C3 fails), else ``qualifies``. A missing rate fails its
    check. The pooled SE of a difference is √(se_arm² + se_ref²).
    """

    index = _index(rates, "design", "cell", "arm", "statistic")
    plan = _criteria_plan(designs)

    def rate_of(design: str, cell: str, arm: str, statistic: str) -> tuple[float, float, int] | None:
        row = index.get((design, cell, arm, statistic))
        if row is None or _number(row.get("rate")) is None:
            return None
        return float(row["rate"]), float(row["se"]), int(row["count"])

    def check(arm: str, design: str, cell: str, statistic: str, kind: str) -> dict[str, Any]:
        entry: dict[str, Any] = {"design": design, "cell": cell, "statistic": statistic, "kind": kind}
        mine = rate_of(design, cell, arm, statistic)
        if kind == "size":
            if mine is None:
                return {**entry, "value": None, "threshold": None, "pass": False, "missing": True}
            bound = size_bound(alpha, mine[2])
            return {**entry, "value": mine[0], "count": mine[2], "threshold": bound, "pass": mine[0] <= bound}
        reference = rate_of(design, cell, REFERENCE_ARM, statistic)
        if mine is None or reference is None:
            return {**entry, "value": None, "reference": None, "threshold": None, "pass": False, "missing": True}
        margin = 2.0 * math.sqrt(mine[1] ** 2 + reference[1] ** 2)
        difference = mine[0] - reference[0]
        passed = difference > margin if kind == "gain" else difference >= -margin
        return {
            **entry,
            "value": mine[0],
            "reference": reference[0],
            "difference": difference,
            "threshold": margin,
            "pass": passed,
        }

    arms: dict[str, Any] = {}
    for arm in ELIGIBLE_ARMS:
        criteria = {}
        for name, checks in plan.items():
            results = [check(arm, *spec) for spec in checks]
            criteria[name] = {
                "pass": all(result["pass"] for result in results),
                "checks": results,
                "failing": [
                    {key: result[key] for key in ("design", "cell", "statistic")}
                    for result in results
                    if not result["pass"]
                ],
            }
        if not criteria["size"]["pass"]:
            verdict = "size_excess"
        elif not criteria["gain"]["pass"]:
            verdict = "no_gain"
        elif not criteria["protected"]["pass"]:
            verdict = "trade_off"
        else:
            verdict = "qualifies"
        arms[arm] = {
            "verdict": verdict,
            "failing": {name: block["failing"] for name, block in criteria.items() if block["failing"]},
            "criteria": criteria,
        }
    qualifying = sorted(arm for arm, entry in arms.items() if entry["verdict"] == "qualifies")
    return {
        "rule_version": RULE_VERSION,
        "advisory": True,
        "alpha": alpha,
        "reference_arm": REFERENCE_ARM,
        "eligible_arms": list(ELIGIBLE_ARMS),
        "ineligible_arms": {CEILING_ARM: "the reference ceiling; never eligible"},
        "verdict": "propose_paper_grade" if qualifying else "document_tradeoff",
        "qualifying_arms": qualifying,
        "arms": arms,
        "descriptive": descriptive_tables(rates, captures, alpha=alpha, designs=designs),
    }


def descriptive_tables(
    rates: Sequence[Mapping[str, Any]],
    captures: Sequence[Mapping[str, Any]],
    *,
    alpha: float,
    designs: Sequence[str] = DESIGNS,
) -> dict[str, Any]:
    """Visibility against φ, ``pls_cv``'s controls, the weak design and the sparse-minus-dense contrast."""

    rate_index = _index(rates, "design", "cell", "arm", "statistic")
    capture_index = _index(captures, "design", "cell", "arm")

    def rate(design: str, cell: str, arm: str, statistic: str) -> dict[str, Any] | None:
        row = rate_index.get((design, cell, arm, statistic))
        if row is None or _number(row.get("rate")) is None:
            return None
        return {"rate": float(row["rate"]), "se": float(row["se"]), "count": int(row["count"])}

    def capture(design: str, cell: str, arm: str, column: str = "group_difference_median") -> float | None:
        row = capture_index.get((design, cell, arm))
        return None if row is None else _number(row.get(column))

    visibility = []
    for design in designs:
        for arm in ARMS:
            for effect in (0.25, 0.5):
                for phi in (0, 30, 60, 90):
                    cell = f"orient_e{effect:g}_phi{phi}"
                    visibility.append(
                        {
                            "design": design,
                            "arm": arm,
                            "effect": effect,
                            "phi": phi,
                            "capture": capture(design, cell, arm),
                            "angle": rate(design, cell, arm, "angle"),
                        }
                    )
    controls = []
    for design in designs:
        for cell in ("anchor", "translation"):
            for statistic in STATISTICS:
                entry = rate(design, cell, REFERENCE_ARM, statistic)
                bound = None if entry is None else size_bound(alpha, entry["count"])
                controls.append(
                    {
                        "design": design,
                        "cell": cell,
                        "statistic": statistic,
                        "rate": entry,
                        "bound": bound,
                        "holds": None if entry is None or bound is None else entry["rate"] <= bound,
                    }
                )
    weak = []
    if "seaad_weak" in designs:
        for arm in ARMS:
            for cell, mode, *_ in CELLS:
                target = TARGET_STATISTIC[mode]
                weak.append(
                    {
                        "arm": arm,
                        "cell": cell,
                        "target_statistic": target,
                        "target": None if target is None else rate("seaad_weak", cell, arm, target),
                        "capture": capture("seaad_weak", cell, arm),
                        "rank": capture("seaad_weak", cell, arm, "retained_rank_median"),
                    }
                )
    contrast = []
    for design in designs:
        for arm in ARMS:
            sparse_cell, dense_cell = "orient_e0.5_phi90_k20", "orient_e0.5_phi90"
            power = [rate(design, sparse_cell, arm, "angle"), rate(design, dense_cell, arm, "angle")]
            sparse_row = capture_index.get((design, sparse_cell, arm))
            dense_row = capture_index.get((design, dense_cell, arm))
            capture_difference: float | None = None
            capture_se: float | None = None
            if sparse_row is not None and dense_row is not None:
                sparse_mean, dense_mean = (_number(row.get("group_difference_mean")) for row in (sparse_row, dense_row))
                sparse_sd, dense_sd = (_number(row.get("group_difference_sd")) for row in (sparse_row, dense_row))
                if sparse_mean is not None and dense_mean is not None:
                    capture_difference = sparse_mean - dense_mean
                if sparse_sd is not None and dense_sd is not None:
                    capture_se = math.sqrt(sparse_sd**2 / int(sparse_row["n"]) + dense_sd**2 / int(dense_row["n"]))
            power_difference = power_se = None
            if power[0] is not None and power[1] is not None:
                power_difference = power[0]["rate"] - power[1]["rate"]
                power_se = math.sqrt(power[0]["se"] ** 2 + power[1]["se"] ** 2)
            contrast.append(
                {
                    "design": design,
                    "arm": arm,
                    "capture_difference": capture_difference,
                    "capture_pooled_se": capture_se,
                    "angle_difference": power_difference,
                    "angle_pooled_se": power_se,
                }
            )
    return {
        "visibility_against_phi": visibility,
        "pls_cv_controls": controls,
        "weak_design": weak,
        "sparse_minus_dense": contrast,
    }


def _read_csv(path: Path) -> list[dict[str, Any]]:
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle))


def decide(out_dir: Path, settings: Mapping[str, Any]) -> dict[str, Any]:
    """Read the merged tables and write ``latent_space_comparison_reading.json``."""

    ensure_settings(out_dir, settings)
    accounting_path = out_dir / "accounting.json"
    if not accounting_path.exists():
        raise SystemExit(f"{accounting_path} does not exist; run `merge` first.")
    summary = json.loads(accounting_path.read_text(encoding="utf-8"))
    if summary.get("settings_hash") != settings_hash(settings):
        raise SystemExit("accounting.json was merged under different settings; run `merge` again.")
    result = reading(
        _read_csv(out_dir / "rate_table.csv"),
        _read_csv(out_dir / "capture_table.csv"),
        alpha=float(settings["alpha"]),
        designs=list(settings["designs"]),
    )
    result.update(
        {
            "settings_hash": settings_hash(settings),
            "replicates": settings["replicates"],
            "permutations": settings["permutations"],
            "weak_delta": settings["weak_delta"],
            "complete": summary["missing_units"] == 0,
            "missing_units": summary["missing_units"],
            "failed_records": summary["failed_records"],
        }
    )
    _write_json(out_dir / "latent_space_comparison_reading.json", result)
    return result


# ── Calibration of the weak design (design D9) ───────────────────────────────


def calibration_settings(*, replicates: int = CALIBRATION_REPLICATES, base_seed: int = BASE_SEED) -> dict[str, Any]:
    generator = design_generators(CALIBRATION_LADDER[0])["seaad"]
    return {
        "driver_version": DRIVER_VERSION,
        "design": "seaad",
        "generator": generator,
        "ladder": list(CALIBRATION_LADDER),
        "replicates": int(replicates),
        "band": list(CALIBRATION_BAND),
        "target": CALIBRATION_TARGET,
        "pls_cv_params": dict(PLS_CV_PARAMS),
        "min_components": MIN_COMPONENTS,
        "base_seed": int(base_seed),
    }


def calibration_unit(delta: float, replicate: int, settings: Mapping[str, Any]) -> dict[str, Any]:
    """One anchor replicate of the SEA-AD table at ``delta``: the ``pls_cv`` fit only, no permutations."""

    from motco.simulations.evaluation import integrate_semisynthetic_dataset
    from motco.simulations.reference import load_reference
    from motco.simulations.semisynthetic import generate_semisynthetic_trajectory

    start = perf_counter()
    seed = calibration_seed(settings["base_seed"], replicate)
    entry: dict[str, Any] = {"delta": delta, "replicate": replicate, "generator_seed": seed}
    try:
        params = _generator_params(
            settings["generator"],
            seed=seed,
            delta_methyl=delta,
            delta_expr=delta,
            delta_protein=delta,
        )
        dataset = generate_semisynthetic_trajectory(params, reference=load_reference())
        fit_settings = {**settings, "permutations": 0, "orientation_estimator": ORIENTATION_ESTIMATOR}
        evaluation = arm_params("pls_cv", fit_settings, seed=0, n_stages=int(settings["generator"]["n_stages"]))
        latent = integrate_semisynthetic_dataset(dataset, evaluation)
        entry.update(
            {
                "status": "completed",
                "auroc": float(latent.metadata["cv_mean_auroc"]),
                "rank": int(latent.metadata["selected_lv"]),
            }
        )
    except Exception as exc:  # a failed fit is recorded, never fatal
        entry.update({"status": "failed", "error_type": type(exc).__name__, "error_message": str(exc)})
    entry["seconds"] = perf_counter() - start
    return entry


def summarize_step(delta: float, units: Sequence[Mapping[str, Any]]) -> dict[str, Any]:
    completed = [u for u in units if u.get("status") == "completed"]
    aurocs = [float(u["auroc"]) for u in completed]
    ranks = [int(u["rank"]) for u in completed]
    return {
        "delta": delta,
        "replicates": len(units),
        "completed": len(completed),
        "auroc_median": _quantile(aurocs, 0.5),
        "auroc_q25": _quantile(aurocs, 0.25),
        "auroc_q75": _quantile(aurocs, 0.75),
        "rank_median": float(median(ranks)) if ranks else None,
    }


def adopt_calibration(
    steps: Sequence[Mapping[str, Any]],
    *,
    band: tuple[float, float] = CALIBRATION_BAND,
    target: float = CALIBRATION_TARGET,
) -> dict[str, Any]:
    """The first ascending step whose median AUROC lies in ``band``, else the closest to ``target``.

    The closest step is marked ``band_missed``; a tie goes to the smaller δ.
    """

    ordered = sorted((step for step in steps if step.get("auroc_median") is not None), key=lambda s: s["delta"])
    if not ordered:
        raise SystemExit("No calibration step produced an AUROC; nothing to adopt.")
    for step in ordered:
        if band[0] <= float(step["auroc_median"]) <= band[1]:
            return {"delta": step["delta"], "rule": "first_in_band", "auroc_median": step["auroc_median"]}
    closest = min(ordered, key=lambda s: abs(float(s["auroc_median"]) - target))
    return {"delta": closest["delta"], "rule": "band_missed", "auroc_median": closest["auroc_median"]}


CALIBRATION_UNIT_COLUMNS = ("delta", "replicate", "generator_seed", "status", "auroc", "rank", "seconds")


def calibrate(out_dir: Path, settings: Mapping[str, Any], *, workers: int = 1) -> dict[str, Any]:
    """Run every ladder step and write ``calibration.json`` and ``calibration_units.csv``."""

    out_dir.mkdir(parents=True, exist_ok=True)
    jobs = [(float(delta), replicate) for delta in settings["ladder"] for replicate in range(settings["replicates"])]
    if workers <= 1:
        units = [calibration_unit(delta, replicate, settings) for delta, replicate in jobs]
    else:
        with ProcessPoolExecutor(max_workers=workers) as pool:
            units = list(pool.map(calibration_unit, *zip(*jobs), [settings] * len(jobs)))
    steps = [
        summarize_step(float(delta), [u for u in units if u["delta"] == float(delta)]) for delta in settings["ladder"]
    ]
    adopted = adopt_calibration(steps, band=tuple(settings["band"]), target=float(settings["target"]))
    payload = {
        "settings": settings,
        "hash": settings_hash(settings),
        "steps": steps,
        "adopted": adopted,
        "failures": sum(1 for u in units if u.get("status") != "completed"),
        "core_hours": sum(float(u["seconds"]) for u in units) / 3600.0,
    }
    _write_json(out_dir / "calibration.json", payload)
    _write_csv(out_dir / "calibration_units.csv", units, CALIBRATION_UNIT_COLUMNS)
    return payload


def load_calibration(path: Path) -> dict[str, Any]:
    if not path.exists():
        raise SystemExit(f"{path} does not exist; run `calibrate` first or pass --weak-delta.")
    payload = json.loads(path.read_text(encoding="utf-8"))
    if payload.get("hash") != settings_hash(payload.get("settings", {})):
        raise SystemExit(f"{path} is corrupt: its hash does not match its settings.")
    if "adopted" not in payload:
        raise SystemExit(f"{path} adopts no δ.")
    return payload


# ── CLI ──────────────────────────────────────────────────────────────────────


def parse_shard(text: str) -> tuple[int, int]:
    try:
        index_text, count_text = text.split("/")
        index, count = int(index_text), int(count_text)
    except ValueError as exc:
        raise argparse.ArgumentTypeError(f"--shard must be i/N, got {text!r}.") from exc
    if count < 1 or not 0 <= index < count:
        raise argparse.ArgumentTypeError(f"--shard index must be in [0, N), got {text!r}.")
    return index, count


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Latent-space comparison across alignment (no R).")
    sub = parser.add_subparsers(dest="command", required=True)

    calibrate_command = sub.add_parser("calibrate", help="Calibrate the weak design's δ to SEA-AD's stage AUROC")
    calibrate_command.add_argument("--out-dir", type=Path, required=True, help="Calibration output directory")
    calibrate_command.add_argument(
        "--reps", type=int, default=CALIBRATION_REPLICATES, help=f"Replicates per step ({CALIBRATION_REPLICATES})"
    )
    calibrate_command.add_argument("--base-seed", type=int, default=BASE_SEED, help=f"Base seed ({BASE_SEED})")
    calibrate_command.add_argument("--workers", type=int, default=1, help="Fits run in parallel")

    run = sub.add_parser("run", help="Evaluate (a shard of) the units")
    run.add_argument("--out-dir", type=Path, required=True, help="Output directory")
    run.add_argument("--calibration", type=Path, default=None, help="calibration.json whose adopted δ is used")
    run.add_argument("--weak-delta", type=float, default=None, help="Explicit weak-design δ (smoke and timing)")
    run.add_argument("--reps", type=int, default=DEFAULT_REPLICATES, help=f"Replicates ({DEFAULT_REPLICATES})")
    run.add_argument("--perms", type=int, default=DEFAULT_PERMUTATIONS, help=f"Permutations ({DEFAULT_PERMUTATIONS})")
    run.add_argument("--alpha", type=float, default=DEFAULT_ALPHA, help=f"Reject at p < alpha ({DEFAULT_ALPHA})")
    run.add_argument("--base-seed", type=int, default=BASE_SEED, help=f"Base seed ({BASE_SEED})")
    run.add_argument("--designs", nargs="+", default=list(DESIGNS), choices=list(DESIGNS), help="Designs to run")
    run.add_argument("--shard", type=parse_shard, default=(0, 1), help="Run shard i of N (i/N; default 0/1)")
    run.add_argument("--workers", type=int, default=1, help="Units evaluated in parallel (not part of the settings)")
    run.add_argument("--limit", type=int, default=None, help="Evaluate at most this many pending units")

    for name, text in (("merge", "Merge every shard into the tables"), ("decide", "Write the advisory reading")):
        command = sub.add_parser(name, help=f"{text} (reads settings.json)")
        command.add_argument("--out-dir", type=Path, required=True, help="Output directory")
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    if args.command == "calibrate":
        payload = calibrate(
            args.out_dir, calibration_settings(replicates=args.reps, base_seed=args.base_seed), workers=args.workers
        )
        print(json.dumps({"adopted": payload["adopted"], "failures": payload["failures"]}, indent=1))
        return 0
    if args.command == "run":
        settings = build_settings(
            weak_delta=args.weak_delta,
            calibration=args.calibration,
            replicates=args.reps,
            permutations=args.perms,
            alpha=args.alpha,
            base_seed=args.base_seed,
            designs=args.designs,
        )
        ensure_settings(args.out_dir, settings)
        index, count = args.shard
        evaluated = run_shard(args.out_dir, settings, index, count, limit=args.limit, workers=args.workers)
        print(f"evaluated {evaluated} unit(s) into {shard_path(args.out_dir, index, count)}", file=sys.stderr)
        return 0
    settings = load_settings(args.out_dir)
    if args.command == "merge":
        print(json.dumps(merge(args.out_dir, settings), indent=1))
        return 0
    result = decide(args.out_dir, settings)
    print(json.dumps({key: result[key] for key in ("verdict", "qualifying_arms", "complete")}, indent=1))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
