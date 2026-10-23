#!/usr/bin/env python3
"""Orientation measurement probe: latent space × orientation estimator on matched data (no R).

The probe behind ``openspec/changes/archive/2026-10-12-orientation-measurement-probe/``. At the
Phase 6 SEA-AD design, orientation/``angle`` power is 0.300 at the strongest
construction, and the loss sits in the measurement: the stage-supervised PLS
space collapses a rotation the observed space preserves, and the signed PC1
``angle`` has a heavy null tail. The probe asks which measurement design
recovers orientation, without giving up Type I control, the other statistics'
power, or specificity.

An **arm** is a measurement space × an orientation estimator:

- spaces: ``pls_cv`` (production stage-PLS, CV rank — the reference),
  ``pls_forced`` (stage-PLS at rank ``2(K − 1)``), ``pls_within_group``
  (``K − 1`` components per group), ``pca`` (rank ``2(K − 1)``), and ``concat``
  (ceiling reference only; never advances);
- estimators: ``pc1`` and ``net_displacement``.

Each (design, cell, replicate) **unit** generates its dataset once; every space
is fitted once on it and evaluated under both estimators with one shared RRPP
seed, so arms differ only in measurement. The generator seed depends on
(design, replicate) only — matched across cells, as in the study engine — and
the RRPP seed on the unit.

Outputs (``--out-dir``):

- ``settings.json``: the probe settings and their hash. A resume or merge into a
  directory written with different settings exits with the differing keys.
- ``records/shard_<i>_of_<N>.jsonl``: one line per (unit, arm). Resume skips
  completed units.
- after ``--merge`` (automatic for an unsharded run): ``merged.jsonl``,
  ``rate_table.csv``, ``checkpoint_table.csv``, ``advancement_decision.json`` /
  ``.csv`` and ``summary.md``.

The advancement rule is predeclared (see :func:`advancement_decision`) and
advisory: it picks the arms for the paper-grade follow-up and changes no
production default.

Examples:

    # smoke
    python scripts/orientation_measurement_probe.py --reps 2 --perms 19 \
        --phase5-n-samples 160 --out-dir /tmp/omp

    # one SLURM array task, then the merge
    python scripts/orientation_measurement_probe.py --shard 3/40 --out-dir results/omp
    python scripts/orientation_measurement_probe.py --merge --out-dir results/omp
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import sys
from dataclasses import dataclass, replace
from pathlib import Path
from statistics import median
from time import perf_counter
from typing import Any, Iterable, Mapping, Sequence

import numpy as np

PROBE_VERSION = 1
RULE_VERSION = 1

SPACES: tuple[str, ...] = ("pls_cv", "pls_forced", "pls_within_group", "pca", "concat")
ESTIMATORS: tuple[str, ...] = ("pc1", "net_displacement")
REFERENCE_ARM = "pls_cv__pc1"
CEILING_SPACE = "concat"
STATISTICS: tuple[str, ...] = ("delta", "angle", "shape")

#: (cell name, trajectory mode, effect size).
CELLS: tuple[tuple[str, str, float], ...] = (
    ("anchor", "none", 0.0),
    ("orientation_0.5", "orientation", 0.5),
    ("orientation_1.0", "orientation", 1.0),
    ("magnitude_1.0", "magnitude", 1.0),
    ("shape_1.0", "shape", 1.0),
)

#: Stage-PLS cross-validation knobs of the Phase 5 / Phase 6 study configs.
PLS_CV_PARAMS: dict[str, Any] = {
    "cv1_splits": 3,
    "cv2_splits": 4,
    "n_repeats": 5,
    "max_components": 20,
    "random_state": 1203,
}

PHASE5_DEVIATION_NOTE = (
    "The Phase 5 design here uses the joint magnitude construction (magnitude_kind='joint'), not the "
    "methylation-only construction ('all') pinned in the Phase 5 paper-grade config "
    "(examples/trajectory_power_study/phase5_power_study.json). Its reference-arm rates are recomputed in "
    "this probe, and every comparison is against them."
)


def design_generators(phase5_n_samples: int = 1200) -> dict[str, dict[str, Any]]:
    """Generator parameters per design (mode, effect and seed are set per unit)."""

    return {
        "phase6": {
            "n_stages": 3,
            "group_stage_sizes": [[11, 10, 28], [9, 10, 12]],
            "p_dmp": 0.1,
            "baseline_continuity": 0.0,
            "magnitude_kind": "joint",
            "surgery_censoring": "error",
        },
        "phase5": {
            "n_samples": phase5_n_samples,
            "n_stages": 4,
            "group_ratio": 0.5,
            "p_dmp": 0.1,
            "baseline_continuity": 0.0,
            "magnitude_kind": "joint",
            "surgery_censoring": "error",
        },
    }


#: Measured blocks per design: Phase 6 measures methylation + expression.
DESIGN_LAYERS: dict[str, list[str] | None] = {"phase6": ["methylation", "expression"], "phase5": None}


def space_params(space: str, n_stages: int, layers: list[str] | None) -> tuple[str, dict[str, Any]]:
    """(integration_method, integration_params) of one measurement space."""

    selection: dict[str, Any] = {} if layers is None else {"layers": list(layers)}
    rank = 2 * (n_stages - 1)
    if space == "pls_cv":
        return "pls", {**PLS_CV_PARAMS, **selection}
    if space == "pls_forced":
        return "pls", {"forced_components": rank, **selection}
    if space == "pls_within_group":
        return "pls_within_group", {"components_per_group": n_stages - 1, **selection}
    if space == "pca":
        return "pca", {"n_components": rank, **selection}
    if space == "concat":
        return "concat", dict(selection)
    raise ValueError(f"Unknown space {space!r}.")


def arm_name(space: str, estimator: str) -> str:
    return f"{space}__{estimator}"


ARMS: tuple[str, ...] = tuple(arm_name(space, estimator) for space in SPACES for estimator in ESTIMATORS)


# ── Settings, seeds and units ────────────────────────────────────────────────


def build_settings(args: argparse.Namespace) -> dict[str, Any]:
    generators = design_generators(args.phase5_n_samples)
    return {
        "probe_version": PROBE_VERSION,
        "designs": {name: {"generator": generators[name], "layers": DESIGN_LAYERS[name]} for name in args.designs},
        "cells": [list(cell) for cell in CELLS],
        "spaces": list(SPACES),
        "estimators": list(ESTIMATORS),
        "replicates": args.reps,
        "permutations": args.perms,
        "base_seed": args.base_seed,
        "alpha": args.alpha,
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
                f"{path} was written with different probe settings; differing: {', '.join(differing)}. "
                "Use a fresh --out-dir."
            )
        return
    out_dir.mkdir(parents=True, exist_ok=True)
    payload = {"hash": settings_hash(settings), "settings": settings}
    path.write_text(json.dumps(payload, indent=1, sort_keys=True) + "\n", encoding="utf-8")


def _seed(*parts: Any) -> int:
    encoded = json.dumps(["orientation-measurement-probe", *parts], separators=(",", ":")).encode("utf-8")
    return int(hashlib.sha256(encoded).hexdigest()[:8], 16) & 0x7FFFFFFF


def generator_seed(base_seed: int, design: str, replicate: int) -> int:
    """Matched across cells: the same (design, replicate) gives the same seed."""

    return _seed(base_seed, design, "generator", replicate)


def rrpp_seed(base_seed: int, design: str, cell: str, replicate: int) -> int:
    """One RRPP seed per unit, shared by every arm."""

    return _seed(base_seed, design, cell, "rrpp", replicate)


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
    """Round-robin, so every shard mixes cheap and expensive designs."""

    return [unit for position, unit in enumerate(units) if position % count == index]


def parse_shard(text: str) -> tuple[int, int]:
    try:
        index_text, count_text = text.split("/")
        index, count = int(index_text), int(count_text)
    except ValueError as exc:
        raise argparse.ArgumentTypeError(f"--shard must be i/N, got {text!r}.") from exc
    if count < 1 or not 0 <= index < count:
        raise argparse.ArgumentTypeError(f"--shard index must be in [0, N), got {text!r}.")
    return index, count


# ── Unit evaluation ──────────────────────────────────────────────────────────


def _checkpoint_angle(geometry: Mapping[str, Any], checkpoint: str, scope: str) -> float | None:
    entry = (geometry.get("checkpoints") or {}).get(checkpoint, {}).get(scope)
    return None if entry is None else entry.get("angle")


def _finite(values: Mapping[str, float]) -> dict[str, float | None]:
    return {key: (float(value) if math.isfinite(float(value)) else None) for key, value in values.items()}


def _retained_dimensionality(metadata: Mapping[str, Any]) -> int:
    for key in ("selected_lv", "retained_dimensionality", "n_features"):
        if metadata.get(key) is not None:
            return int(metadata[key])
    raise KeyError("integration metadata records no dimensionality")


def evaluate_unit(unit: Unit, settings: Mapping[str, Any]) -> list[dict[str, Any]]:
    """Generate one dataset and evaluate every arm on it; one record per arm."""

    from motco.simulations.evaluation import (
        SimulationEvaluationParams,
        evaluate_semisynthetic_trajectory,
        integrate_semisynthetic_dataset,
    )
    from motco.simulations.reference import load_reference
    from motco.simulations.semisynthetic import (
        SemiSyntheticTrajectoryParams,
        generate_semisynthetic_trajectory,
    )

    design = settings["designs"][unit.design]
    generator = dict(design["generator"])
    if generator.get("group_stage_sizes") is not None:
        generator["group_stage_sizes"] = tuple(tuple(int(n) for n in row) for row in generator["group_stage_sizes"])
    cell = {name: (mode, effect) for name, mode, effect in settings["cells"]}[unit.cell]
    gen_seed = generator_seed(settings["base_seed"], unit.design, unit.replicate)
    perm_seed = rrpp_seed(settings["base_seed"], unit.design, unit.cell, unit.replicate)
    params = SemiSyntheticTrajectoryParams(
        seed=gen_seed, trajectory_mode=cell[0], group_effect_size=float(cell[1]), **generator
    )
    unit_start = perf_counter()
    dataset = generate_semisynthetic_trajectory(params, reference=load_reference())
    n_stages = int(generator["n_stages"])
    base = {
        "design": unit.design,
        "cell": unit.cell,
        "trajectory_mode": cell[0],
        "effect_size": float(cell[1]),
        "replicate": unit.replicate,
        "generator_seed": gen_seed,
        "rrpp_seed": perm_seed,
        "settings_hash": settings_hash(settings),
    }
    records = []
    for space in settings["spaces"]:
        method, integration_params = space_params(space, n_stages, design["layers"])
        evaluation = SimulationEvaluationParams(
            integration_method=method,  # type: ignore[arg-type]
            integration_params=integration_params,
            permutations=int(settings["permutations"]),
            n_jobs=1,
            seed=perm_seed,
            min_components=1,  # the 2026-10-12 probe ran the unfloored rank rule
        )
        space_start = perf_counter()
        try:
            latent = integrate_semisynthetic_dataset(dataset, evaluation)
        except Exception as exc:  # recorded per arm, never aborts the unit
            for estimator in settings["estimators"]:
                records.append(_failed(base, space, estimator, exc))
            continue
        fit_seconds = perf_counter() - space_start
        for estimator in settings["estimators"]:
            arm_start = perf_counter()
            try:
                result = evaluate_semisynthetic_trajectory(
                    dataset, replace(evaluation, orientation_estimator=estimator), latent=latent
                )
            except Exception as exc:
                records.append(_failed(base, space, estimator, exc))
                continue
            geometry = result.realized_geometry
            latent_checkpoint = next((name for name in geometry["checkpoints"] if name.endswith("_latent")), None)
            records.append(
                {
                    **base,
                    "arm": arm_name(space, estimator),
                    "space": space,
                    "estimator": estimator,
                    "status": "completed",
                    "integration_method": method,
                    "retained_dimensionality": _retained_dimensionality(result.latent_matrix_metadata),
                    "pair_statistics": _finite(result.pair_statistics),
                    "p_values": _finite(result.p_values),
                    "realized_angle": {
                        "population": _checkpoint_angle(geometry, "population_native", "methylation"),
                        "population_standardized": _checkpoint_angle(geometry, "population_standardized", "joint"),
                        "observed_standardized": _checkpoint_angle(geometry, "observed_standardized", "joint"),
                        "latent": None
                        if latent_checkpoint is None
                        else _checkpoint_angle(geometry, latent_checkpoint, "joint"),
                    },
                    "pooled_eigengap": (result.config_spectrum.get("pooled") or {}).get("relative_eigengap"),
                    "integration_metadata": result.latent_matrix_metadata,
                    "runtime_seconds": {"space_fit": fit_seconds, "arm": perf_counter() - arm_start},
                }
            )
    unit_seconds = perf_counter() - unit_start
    for record in records:
        record["unit_runtime_seconds"] = unit_seconds
    return records


def _failed(base: Mapping[str, Any], space: str, estimator: str, exc: Exception) -> dict[str, Any]:
    return {
        **base,
        "arm": arm_name(space, estimator),
        "space": space,
        "estimator": estimator,
        "status": "failed",
        "error_type": type(exc).__name__,
        "error_message": str(exc),
    }


# ── Persistence and resume ───────────────────────────────────────────────────


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


def complete_units(records: Iterable[Mapping[str, Any]], arms: Sequence[str]) -> dict[tuple[str, str, int], list]:
    """Records grouped by unit, keeping only units with every arm present exactly once."""

    grouped: dict[tuple[str, str, int], dict[str, Mapping[str, Any]]] = {}
    for record in records:
        key = (str(record["design"]), str(record["cell"]), int(record["replicate"]))
        grouped.setdefault(key, {})[str(record["arm"])] = record
    return {key: [arms_[a] for a in arms] for key, arms_ in grouped.items() if set(arms_) == set(arms)}


def run_shard(out_dir: Path, settings: Mapping[str, Any], index: int, count: int, *, limit: int | None = None) -> int:
    """Evaluate this shard's missing units; returns how many were evaluated."""

    path = shard_path(out_dir, index, count)
    path.parent.mkdir(parents=True, exist_ok=True)
    done = complete_units(read_records([path]) if path.exists() else [], arm_list(settings))
    # Rewrite the shard with complete units only, so a torn or partial unit
    # left by an interruption is re-run instead of duplicated.
    with path.open("w", encoding="utf-8") as handle:
        for key in sorted(done):
            for record in done[key]:
                handle.write(json.dumps(_jsonable(record), sort_keys=True, allow_nan=False) + "\n")
    pending = [unit for unit in shard_units(enumerate_units(settings), index, count) if unit.key not in done]
    if limit is not None:
        pending = pending[:limit]
    for position, unit in enumerate(pending, start=1):
        records = evaluate_unit(unit, settings)
        lines = "".join(json.dumps(_jsonable(r), sort_keys=True, allow_nan=False) + "\n" for r in records)
        with path.open("a", encoding="utf-8") as handle:
            handle.write(lines)
            handle.flush()
        print(f"[shard {index}/{count}] {position}/{len(pending)} {unit.key}", file=sys.stderr, flush=True)
    return len(pending)


def arm_list(settings: Mapping[str, Any]) -> list[str]:
    return [arm_name(space, estimator) for space in settings["spaces"] for estimator in settings["estimators"]]


# ── Rate table and advancement rule (pure functions) ─────────────────────────


RATE_COLUMNS = ("design", "cell", "arm", "statistic", "rate", "se", "count", "rejections")


def rate_table(records: Iterable[Mapping[str, Any]], *, alpha: float) -> list[dict[str, Any]]:
    """Rejection rate (p < alpha), binomial SE and count per design × cell × arm × statistic."""

    tallies: dict[tuple[str, str, str, str], list[int]] = {}
    for record in records:
        if record.get("status") != "completed":
            continue
        for statistic in STATISTICS:
            p_value = (record.get("p_values") or {}).get(statistic)
            if p_value is None or not math.isfinite(float(p_value)):
                continue
            key = (record["design"], record["cell"], record["arm"], statistic)
            tally = tallies.setdefault(key, [0, 0])
            tally[0] += 1
            tally[1] += int(float(p_value) < alpha)
    rows = []
    for (design, cell, arm, statistic), (count, rejections) in sorted(tallies.items()):
        rate = rejections / count
        rows.append(
            {
                "design": design,
                "cell": cell,
                "arm": arm,
                "statistic": statistic,
                "rate": rate,
                "se": math.sqrt(rate * (1.0 - rate) / count),
                "count": count,
                "rejections": rejections,
            }
        )
    return rows


def anchor_bound(alpha: float, replicates: int) -> float:
    return alpha + 2.0 * math.sqrt(alpha * (1.0 - alpha) / replicates)


def _pooled_se(left: Mapping[str, Any], right: Mapping[str, Any]) -> float:
    return math.sqrt(left["se"] ** 2 + right["se"] ** 2)


def advancement_decision(rows: Sequence[Mapping[str, Any]], *, alpha: float) -> dict[str, Any]:
    """Apply the predeclared advisory rule to a rate table.

    A non-reference, non-``concat`` arm advances only if all five criteria
    hold (``openspec/specs/orientation-measurement-probe``):

    1. anchor — at both designs, anchor ``delta``/``angle``/``shape`` rates at
       or below ``alpha + 2·sqrt(alpha(1 − alpha)/R)``;
    2. gain — Phase 6 orientation e = 1.0 ``angle`` above the reference by
       more than 2 pooled SE;
    3. no loss — Phase 5 orientation e = 1.0 ``angle`` not below the reference
       by more than 2 pooled SE;
    4. protected — at both designs, magnitude/``delta`` and shape/``shape`` at
       e = 1.0 not below the reference by more than 2 pooled SE;
    5. specificity — at both designs, magnitude/``angle`` and shape/``angle``
       at e = 1.0 not above the reference by more than 2 pooled SE.

    A criterion whose rates are missing fails. The pooled SE of a difference is
    ``sqrt(se_arm² + se_ref²)`` from the two arms' binomial rates.
    """

    index = {(r["design"], r["cell"], r["arm"], r["statistic"]): r for r in rows}
    arms = sorted({r["arm"] for r in rows})

    def compare(arm: str, design: str, cell: str, statistic: str, kind: str) -> dict[str, Any]:
        check: dict[str, Any] = {"design": design, "cell": cell, "statistic": statistic, "kind": kind}
        mine = index.get((design, cell, arm, statistic))
        if kind == "anchor":
            if mine is None:
                return {**check, "value": None, "threshold": None, "pass": False, "missing": True}
            threshold = anchor_bound(alpha, int(mine["count"]))
            return {**check, "value": mine["rate"], "threshold": threshold, "pass": mine["rate"] <= threshold}
        reference = index.get((design, cell, REFERENCE_ARM, statistic))
        if mine is None or reference is None:
            return {**check, "value": None, "reference": None, "threshold": None, "pass": False, "missing": True}
        margin = 2.0 * _pooled_se(mine, reference)
        difference = mine["rate"] - reference["rate"]
        passed = {"gain": difference > margin, "no_loss": difference >= -margin, "no_excess": difference <= margin}[
            kind
        ]
        return {
            **check,
            "value": mine["rate"],
            "reference": reference["rate"],
            "difference": difference,
            "threshold": margin,
            "pass": passed,
        }

    designs = ("phase6", "phase5")
    plan: dict[str, list[tuple[str, str, str, str]]] = {
        "anchor": [(d, "anchor", s, "anchor") for d in designs for s in STATISTICS],
        "gain": [("phase6", "orientation_1.0", "angle", "gain")],
        "no_loss": [("phase5", "orientation_1.0", "angle", "no_loss")],
        "protected": [
            (d, cell, statistic, "no_loss")
            for d in designs
            for cell, statistic in (("magnitude_1.0", "delta"), ("shape_1.0", "shape"))
        ],
        "specificity": [
            (d, cell, "angle", "no_excess") for d in designs for cell in ("magnitude_1.0", "shape_1.0")
        ],
    }
    decided: dict[str, Any] = {}
    for arm in arms:
        criteria = {}
        for name, checks in plan.items():
            results = [compare(arm, *check) for check in checks]
            criteria[name] = {"pass": all(result["pass"] for result in results), "checks": results}
        if arm == REFERENCE_ARM:
            role, advances, reason = "reference", False, "reference arm"
        elif arm.startswith(f"{CEILING_SPACE}__"):
            role, advances, reason = "ceiling_reference", False, "concat is reference-only and never advances"
        else:
            role = "candidate"
            failed = [name for name, criterion in criteria.items() if not criterion["pass"]]
            advances = not failed
            reason = "all criteria pass" if advances else f"failed: {', '.join(failed)}"
        decided[arm] = {"role": role, "advances": advances, "reason": reason, "criteria": criteria}
    return {
        "rule_version": RULE_VERSION,
        "alpha": alpha,
        "reference_arm": REFERENCE_ARM,
        "advisory": True,
        "advancing": sorted(arm for arm, entry in decided.items() if entry["advances"]),
        "arms": decided,
    }


def decision_rows(decision: Mapping[str, Any]) -> list[dict[str, Any]]:
    """CSV mirror of the decision: one row per arm × criterion check."""

    rows = []
    for arm, entry in decision["arms"].items():
        for criterion, block in entry["criteria"].items():
            for check in block["checks"]:
                rows.append(
                    {
                        "arm": arm,
                        "role": entry["role"],
                        "advances": entry["advances"],
                        "criterion": criterion,
                        "criterion_pass": block["pass"],
                        "design": check["design"],
                        "cell": check["cell"],
                        "statistic": check["statistic"],
                        "value": check.get("value"),
                        "reference": check.get("reference"),
                        "difference": check.get("difference"),
                        "threshold": check.get("threshold"),
                        "pass": check["pass"],
                    }
                )
    return rows


CHECKPOINT_COLUMNS = (
    "design",
    "cell",
    "arm",
    "n",
    "population",
    "population_standardized",
    "observed_standardized",
    "latent",
    "retained_dimensionality",
    "angle_rate",
)


def checkpoint_table(records: Iterable[Mapping[str, Any]], *, alpha: float) -> list[dict[str, Any]]:
    """Median realized ``angle`` per checkpoint and median retained rank, per design × cell × arm."""

    grouped: dict[tuple[str, str, str], list[Mapping[str, Any]]] = {}
    for record in records:
        if record.get("status") == "completed":
            grouped.setdefault((record["design"], record["cell"], record["arm"]), []).append(record)

    def med(values: Iterable[Any]) -> float | None:
        finite = [float(v) for v in values if v is not None and math.isfinite(float(v))]
        return float(median(finite)) if finite else None

    rows = []
    for (design, cell, arm), group in sorted(grouped.items()):
        p_values = [r["p_values"].get("angle") for r in group]
        decisions = [float(p) < alpha for p in p_values if p is not None]
        row: dict[str, Any] = {"design": design, "cell": cell, "arm": arm, "n": len(group)}
        for name in ("population", "population_standardized", "observed_standardized", "latent"):
            row[name] = med(r["realized_angle"].get(name) for r in group)
        row["retained_dimensionality"] = med(r["retained_dimensionality"] for r in group)
        row["angle_rate"] = float(np.mean(decisions)) if decisions else None
        rows.append(row)
    return rows


# ── Merge and render ─────────────────────────────────────────────────────────


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


def merge(out_dir: Path, settings: Mapping[str, Any]) -> dict[str, Any]:
    """Merge every shard, then write the rate table, decision, checkpoint table and summary."""

    paths = sorted((out_dir / "records").glob("shard_*.jsonl"))
    units = complete_units(read_records(paths), arm_list(settings))
    expected = {unit.key for unit in enumerate_units(settings)}
    missing = sorted(expected - set(units))
    records = [record for key in sorted(units) if key in expected for record in units[key]]
    with (out_dir / "merged.jsonl").open("w", encoding="utf-8") as handle:
        for record in records:
            handle.write(json.dumps(_jsonable(record), sort_keys=True, allow_nan=False) + "\n")
    alpha = float(settings["alpha"])
    rates = rate_table(records, alpha=alpha)
    _write_csv(out_dir / "rate_table.csv", rates, RATE_COLUMNS)
    checkpoints = checkpoint_table(records, alpha=alpha)
    _write_csv(out_dir / "checkpoint_table.csv", checkpoints, CHECKPOINT_COLUMNS)
    decision = advancement_decision(rates, alpha=alpha)
    decision["complete"] = not missing
    decision["missing_units"] = len(missing)
    (out_dir / "advancement_decision.json").write_text(
        json.dumps(_jsonable(decision), indent=1, sort_keys=True, allow_nan=False) + "\n", encoding="utf-8"
    )
    decision_columns = (
        "arm", "role", "advances", "criterion", "criterion_pass", "design", "cell", "statistic",
        "value", "reference", "difference", "threshold", "pass",
    )  # fmt: skip
    _write_csv(out_dir / "advancement_decision.csv", decision_rows(decision), decision_columns)
    failures = sum(1 for record in records if record.get("status") != "completed")
    summary = render_summary(settings, rates, checkpoints, decision, n_units=len(units), failures=failures)
    (out_dir / "summary.md").write_text(summary, encoding="utf-8")
    return {"units": len(units), "missing": len(missing), "records": len(records), "failures": failures}


def _rate_text(row: Mapping[str, Any] | None) -> str:
    return "—" if row is None else f"{row['rate']:.3f} ± {row['se']:.3f}"


#: Rate-table columns of the summary: (cell, statistic).
SUMMARY_RATE_COLUMNS: tuple[tuple[str, str], ...] = (
    ("anchor", "delta"),
    ("anchor", "angle"),
    ("anchor", "shape"),
    ("orientation_0.5", "angle"),
    ("orientation_1.0", "angle"),
    ("magnitude_1.0", "delta"),
    ("magnitude_1.0", "angle"),
    ("shape_1.0", "shape"),
    ("shape_1.0", "angle"),
)


def render_summary(
    settings: Mapping[str, Any],
    rates: Sequence[Mapping[str, Any]],
    checkpoints: Sequence[Mapping[str, Any]],
    decision: Mapping[str, Any],
    *,
    n_units: int,
    failures: int,
) -> str:
    index = {(r["design"], r["cell"], r["arm"], r["statistic"]): r for r in rates}
    arms = arm_list(settings)
    lines = [
        "# Orientation measurement probe",
        "",
        "Generated by `scripts/orientation_measurement_probe.py --merge`. Exploratory and advisory: it ranks",
        "measurement arms for the paper-grade follow-up and certifies none.",
        "",
        "## Settings",
        "",
        f"- settings hash: `{settings_hash(settings)}`",
        f"- replicates per cell: {settings['replicates']}; RRPP permutations: {settings['permutations']}; "
        f"alpha: {settings['alpha']} (reject at p < alpha); base seed: {settings['base_seed']}",
        f"- complete units: {n_units}; failed arm records: {failures}; "
        f"missing units: {decision.get('missing_units', 0)}",
        f"- reference arm: `{REFERENCE_ARM}`; `concat` is a ceiling reference and never advances",
        "",
    ]
    for design, block in settings["designs"].items():
        lines.append(f"- `{design}` generator: `{json.dumps(block['generator'], sort_keys=True)}`; "
                     f"layers: `{block['layers'] or 'all'}`")
    lines += ["", "## Phase 5 deviation", "", PHASE5_DEVIATION_NOTE, "", "## Rejection rates (rate ± MC SE)", ""]
    header = "| arm | " + " | ".join(f"{cell} `{statistic}`" for cell, statistic in SUMMARY_RATE_COLUMNS) + " |"
    for design in settings["designs"]:
        lines += [f"### {design}", "", header, "|" + "---|" * (len(SUMMARY_RATE_COLUMNS) + 1)]
        for arm in arms:
            cells = [_rate_text(index.get((design, cell, arm, statistic))) for cell, statistic in SUMMARY_RATE_COLUMNS]
            lines.append(f"| `{arm}` | " + " | ".join(cells) + " |")
        lines.append("")
    lines += [
        "## Checkpoint table (median realized `angle`, degrees; median retained rank)",
        "",
        "| design | cell | arm | n | population | observed std. | latent | rank | `angle` rate |",
        "|---|---|---|---|---|---|---|---|---|",
    ]

    def fmt(value: Any, digits: int = 1) -> str:
        return "—" if value is None else f"{value:.{digits}f}"

    for row in checkpoints:
        if row["cell"] == "anchor" or row["cell"].startswith("orientation"):
            lines.append(
                f"| {row['design']} | {row['cell']} | `{row['arm']}` | {row['n']} | {fmt(row['population'])} | "
                f"{fmt(row['observed_standardized'])} | {fmt(row['latent'])} | "
                f"{fmt(row['retained_dimensionality'])} | {fmt(row['angle_rate'], 3)} |"
            )
    lines += [
        "",
        "## Advancement decision",
        "",
        f"Advancing arms: {', '.join(f'`{a}`' for a in decision['advancing']) or 'none'}.",
        "",
        "| arm | role | anchor | gain | no_loss | protected | specificity | advances | reason |",
        "|---|---|---|---|---|---|---|---|---|",
    ]
    for arm, entry in decision["arms"].items():
        marks = " | ".join("pass" if entry["criteria"][name]["pass"] else "fail" for name in entry["criteria"])
        lines.append(f"| `{arm}` | {entry['role']} | {marks} | {entry['advances']} | {entry['reason']} |")
    lines += ["", "Per-check values and thresholds: `advancement_decision.json` / `.csv`.", ""]
    return "\n".join(lines)


# ── CLI ──────────────────────────────────────────────────────────────────────


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(description="Orientation measurement probe (no R).")
    p.add_argument("--out-dir", type=Path, required=True, help="Output directory")
    p.add_argument("--reps", type=int, default=100, help="Replicates per cell (default: 100)")
    p.add_argument("--perms", type=int, default=199, help="RRPP permutations (default: 199)")
    p.add_argument("--base-seed", type=int, default=900, help="Base seed (default: 900)")
    p.add_argument("--alpha", type=float, default=0.05, help="Rejection threshold, p < alpha (default: 0.05)")
    p.add_argument(
        "--phase5-n-samples", type=int, default=1200, help="Phase 5 design sample size (default: 1200)"
    )
    p.add_argument(
        "--designs", nargs="+", default=["phase6", "phase5"], choices=["phase6", "phase5"], help="Designs to run"
    )
    p.add_argument("--shard", type=parse_shard, default=(0, 1), help="Run shard i of N (i/N; default 0/1)")
    p.add_argument(
        "--merge", action="store_true", help="Merge every shard and render the outputs (reads settings.json)"
    )
    p.add_argument("--limit", type=int, default=None, help="Evaluate at most this many pending units (timing)")
    return p


def load_settings(out_dir: Path) -> dict[str, Any]:
    path = out_dir / "settings.json"
    if not path.exists():
        raise SystemExit(f"{path} does not exist; nothing to merge.")
    stored = json.loads(path.read_text(encoding="utf-8"))
    settings = stored["settings"]
    if stored.get("hash") != settings_hash(settings):
        raise SystemExit(f"{path} is corrupt: its hash does not match its settings.")
    return settings


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    if args.merge:
        # The merge reads the settings the run wrote, so it needs no run flags.
        settings = load_settings(args.out_dir)
    else:
        settings = build_settings(args)
        ensure_settings(args.out_dir, settings)
    if not args.merge:
        index, count = args.shard
        evaluated = run_shard(args.out_dir, settings, index, count, limit=args.limit)
        print(f"evaluated {evaluated} unit(s) into {shard_path(args.out_dir, index, count)}", file=sys.stderr)
        if count > 1 or args.limit is not None:
            return 0
    outcome = merge(args.out_dir, settings)
    print(json.dumps(outcome), file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

