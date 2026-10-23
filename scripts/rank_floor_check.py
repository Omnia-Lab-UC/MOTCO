#!/usr/bin/env python3
"""Rank-floor check at a Phase 6 construction calibrated to select rank 1 (no R).

The gate behind ``openspec/changes/minimum-latent-rank/`` (design D7, spec
``rank-floor-validation``): the size and power of the rank-2 floor of the PLS
latent space (``TrajectoryAnalysisParams.min_components=2``, MOTCO 0.10.0)
against the unfloored 0.9.0 rule, measured where the two rules differ.

Each **unit** is one (cell, replicate):

- **cells**: the zero-effect anchor, the ``translation`` negative control at its
  effect, and ``magnitude``, ``orientation`` and ``shape`` at the Phase 6 power
  profile's effect levels (``group_effect_size``);
- **data**: ``generate_semisynthetic_trajectory`` at the Phase 6 small-n design
  (three stages, ``group_stage_sizes = [[11, 10, 28], [9, 10, 12]]``,
  methylation + expression measured, joint magnitude), with the **stage signal
  weakened** through the per-omic ``delta_*`` and ``p_dmp`` until the unfloored
  rule selects rank 1 in at least half of the null replicates. The generator
  and RRPP seeds depend on the replicate index only, so every cell and both
  arms of replicate ``r`` share the baseline draw and the permutation draws.

Each unit is analysed with :func:`motco.analysis.analyze_trajectories` twice,
with the same CV settings (the harness defaults), permutation count and RRPP
seed, no covariates and the net-displacement estimator:

- ``floor_1``: ``min_components=1`` (the 0.9.0 rule; measured, not judged);
- ``floor_2``: ``min_components=2`` (the candidate default).

Each arm records its p-values, selected rank, mean CV AUROC, null summary and
whether ``angle`` was measured in one dimension (rank 1: the two orientations
are parallel or antiparallel by construction and the ``angle`` null is a
two-point distribution).

**Calibration.** Before any paper-grade run, ``calibrate`` runs the anchor cell
in the ``floor_1`` arm (no permutations: the CV rank does not depend on RRPP)
and reports the fraction of replicates that selected rank 1, the rank
distribution, the mean CV AUROC and the pool headroom of every power cell at
the current strengths; every invocation appends to ``calibration.json`` so the
iterations are recorded. The paper-grade run is admissible only once that
fraction is at or above ``calibration_threshold`` (0.5). The decision reports
the fraction realized in the run itself.

The **decision** (pure, :func:`rank_floor_decision`): for the anchor and the
translation control, the ``floor_2`` arm's rejection rate (``p ≤ alpha``) over
its R available replicates must be at or below ``alpha + k·√(alpha(1 − alpha)/R)``
for ``delta``, ``angle`` and ``shape``. Every statistic holding at both controls
gives ``adopt_floor``, otherwise ``diagnose_floor_excess`` naming the exceeding
cell and statistic. The ``floor_1`` arm's rates are reported and never enter the
decision. Power rows are descriptive: both arms' rates, their difference with
its pooled Monte Carlo SE, and the selected-rank distribution per arm.

Outputs (``--out-dir``):

- ``settings.json``: settings and their hash; a resume into a directory
  written with different settings is refused;
- ``calibration.json``: one entry per ``calibrate`` invocation;
- ``records/shard_<i>_of_<N>.jsonl``: one line per unit (both arms); resume
  skips recorded units;
- after ``merge`` (automatic for an unsharded run): ``merged.jsonl``,
  ``rank_floor_pvalues.csv``, ``rank_floor_replicates.csv`` (per unit × arm:
  rank, CV AUROC, one-dimensional ``angle`` marker, ``angle`` null width),
  ``rank_floor_rates.csv`` and ``rank_floor_decision.json``.

Examples:

    # calibration (anchor only, no permutations; appends to calibration.json)
    python scripts/rank_floor_check.py calibrate --reps 20 --delta-methyl 0.5 --delta-expr 0.5 \\
        --out-dir /tmp/rank-floor-smoke

    # smoke profile at the settled strengths
    python scripts/rank_floor_check.py run --reps 20 --perms 99 --delta-methyl 0.5 --delta-expr 0.5 \\
        --out-dir /tmp/rank-floor-smoke

    # one SLURM array task, then the merge
    python scripts/rank_floor_check.py run --shard 3/50 --out-dir results/rank-floor
    python scripts/rank_floor_check.py merge --out-dir results/rank-floor
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import sys
from concurrent.futures import ProcessPoolExecutor
from dataclasses import asdict
from pathlib import Path
from time import perf_counter
from typing import Any, Iterable, Mapping, Sequence

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import covariate_type_i_check as d10  # noqa: E402  (Phase 6 generator settings, sharding, CSV helpers)

sys.path.pop(0)

SCRIPT_VERSION = 1
CHECK = "rank_floor"
STATISTICS: tuple[str, ...] = d10.STATISTICS
ARMS: tuple[str, ...] = ("floor_1", "floor_2")
ARM_FLOOR: dict[str, int] = {"floor_1": 1, "floor_2": 2}
#: The arm the decision judges; the other is measured only.
JUDGED_ARM = "floor_2"
UNFLOORED_ARM = "floor_1"
MODES: tuple[str, ...] = ("magnitude", "orientation", "shape")
ANCHOR = "anchor"
CONTROL_MODES: tuple[str, ...] = (ANCHOR, "translation")
LAYERS: tuple[str, ...] = d10.LAYERS
CALIBRATION_STATISTIC = "fraction of unfloored (floor_1) anchor replicates whose CV selected rank 1"

#: Phase 6 effect levels per mode (``phase6_small_n_magnitude.json`` below its
#: saturation, ``phase6_small_n_study.json`` for orientation and shape).
DEFAULT_EFFECTS: dict[str, list[float]] = {
    "magnitude": [0.02, 0.05, 0.1, 0.25],
    "orientation": [0.25, 0.5, 0.75, 1.0],
    "shape": [0.25, 0.5, 0.75, 1.0],
}
DEFAULT_TRANSLATION_EFFECT = 1.0
#: The generator's per-omic effect sizes and ``p_dmp`` at the Phase 6 profile.
PHASE6_STRENGTHS: dict[str, float] = {"delta_methyl": 2.0, "delta_expr": 2.0, "delta_protein": 2.0, "p_dmp": 0.1}
STRENGTH_FIELDS: tuple[str, ...] = tuple(PHASE6_STRENGTHS)


# ── Settings ─────────────────────────────────────────────────────────────────


def _floats(text: str) -> list[float]:
    return [float(value) for value in text.split(",") if value.strip()]


def build_settings(args: argparse.Namespace) -> dict[str, Any]:
    generator = dict(d10.PHASE6_GENERATOR)
    for field in STRENGTH_FIELDS:
        generator[field] = float(getattr(args, field))
    return {
        "script_version": SCRIPT_VERSION,
        "check": CHECK,
        "replicates": int(args.reps),
        "permutations": int(args.perms),
        "base_seed": int(args.base_seed),
        "alpha": float(args.alpha),
        "k": float(args.k),
        "generator": generator,
        "layers": list(LAYERS),
        "methylation_units": "mvalue",
        "cv": dict(d10.PLS_CV_PARAMS),
        "orientation_estimator": "net_displacement",
        "arms": {arm: ARM_FLOOR[arm] for arm in ARMS},
        "effects": {mode: _floats(getattr(args, f"{mode}_effects")) for mode in MODES},
        "controls": {ANCHOR: 0.0, "translation": float(args.translation_effect)},
        "calibration": {"statistic": CALIBRATION_STATISTIC, "threshold": float(args.calibration_threshold)},
        "rule": {
            "rate": "share of available replicates with p <= alpha",
            "judged_arm": JUDGED_ARM,
            "controls": "anchor and translation, every statistic, rate <= alpha + k*sqrt(alpha(1-alpha)/R)",
            "decision": "adopt_floor when the judged arm holds every control statistic, else diagnose_floor_excess",
            "power": "descriptive: both arms' rates, difference, pooled SE, rank distributions",
        },
    }


settings_hash = d10.settings_hash
ensure_settings = d10.ensure_settings
load_settings = d10.load_settings
parse_shard = d10.parse_shard
shard_path = d10.shard_path


def _seed(*parts: Any) -> int:
    encoded = json.dumps([CHECK, *parts], separators=(",", ":")).encode("utf-8")
    return int(hashlib.sha256(encoded).hexdigest()[:8], 16) & 0x7FFFFFFF


def replicate_seeds(base_seed: int, replicate: int) -> dict[str, int]:
    """Generator and RRPP seeds depend on the replicate only: matched across cells and arms."""

    return {"generator": _seed(base_seed, "generator", replicate), "rrpp": _seed(base_seed, "rrpp", replicate)}


def cells(settings: Mapping[str, Any]) -> list[tuple[str, float]]:
    """(mode, effect) per cell: anchor, translation, then every mode's effects."""

    out = [(ANCHOR, 0.0), ("translation", float(settings["controls"]["translation"]))]
    for mode in MODES:
        out.extend((mode, float(effect)) for effect in settings["effects"][mode])
    return out


def units(settings: Mapping[str, Any]) -> list[tuple[str, float, int]]:
    """Every (mode, effect, replicate), replicate-major; shards take every N-th."""

    return [
        (mode, effect, replicate)
        for replicate in range(int(settings["replicates"]))
        for mode, effect in cells(settings)
    ]


def unit_key(mode: str, effect: float, replicate: int) -> str:
    return json.dumps([mode, float(effect), int(replicate)])


# ── Construction ─────────────────────────────────────────────────────────────


def generator_params(mode: str, effect: float, replicate: int, settings: Mapping[str, Any]):
    from motco.simulations.semisynthetic import SemiSyntheticTrajectoryParams

    generator = dict(settings["generator"])
    generator["trajectory_mode"] = "magnitude" if mode == ANCHOR else mode
    generator["group_effect_size"] = float(effect)
    generator["group_stage_sizes"] = tuple(tuple(row) for row in generator["group_stage_sizes"])
    seeds = replicate_seeds(int(settings["base_seed"]), replicate)
    return SemiSyntheticTrajectoryParams(seed=seeds["generator"], **generator), seeds


def headroom_report(settings: Mapping[str, Any], reference: Any = None) -> dict[str, Any]:
    """Expected pool headroom of every pool-limited power cell at the settings' strengths."""

    from motco.simulations.semisynthetic import expected_surgery_headroom

    report: dict[str, Any] = {}
    for mode, effect in cells(settings):
        params, _ = generator_params(mode, effect, 0, settings)
        headroom = expected_surgery_headroom(params, reference=reference)
        if headroom is not None:
            report[unit_key(mode, effect, 0)[:-4] + "]"] = asdict(headroom)
    return report


def build_unit(
    mode: str, effect: float, replicate: int, settings: Mapping[str, Any], reference: Any = None
) -> dict[str, Any]:
    """Generate the unit's dataset: methylation (beta) and expression blocks plus metadata."""

    from motco.simulations.reference import load_reference
    from motco.simulations.semisynthetic import generate_semisynthetic_trajectory

    params, seeds = generator_params(mode, effect, replicate, settings)
    if reference is None:
        reference = load_reference()
    dataset = generate_semisynthetic_trajectory(params, reference=reference)
    metadata = dataset.metadata.reset_index(drop=True)
    return {
        "seeds": seeds,
        "metadata": metadata,
        "blocks": {"methylation": dataset.methylation, "expression": dataset.expression},
        "group_labels": sorted(set(metadata["group"].astype(str))),
    }


# ── Unit evaluation ──────────────────────────────────────────────────────────


def _trajectory_dataset(built: Mapping[str, Any], settings: Mapping[str, Any]):
    from motco.analysis import OmicsBlock, TrajectoryDataset

    metadata = built["metadata"]
    transforms = {"methylation": "beta_to_mvalue", "expression": "identity"}
    return TrajectoryDataset(
        blocks=[OmicsBlock(name, built["blocks"][name], transforms[name]) for name in settings["layers"]],
        metadata=metadata,
        sample_col="sample_id",
        group_col="group",
        stage_col="stage",
        stage_order=[str(stage) for stage in sorted(set(metadata["stage"].tolist()))],
    )


def _arm_record(result: Any, floor: int) -> dict[str, Any]:
    record = result.record
    latent = record["latent"]
    return {
        "min_components": floor,
        "statistics": record["statistics"],
        "p_values": record["p_values"],
        "null_summary": record["null_summary"],
        "selected_lv": latent["selected_lv"],
        "component_selection": latent["component_selection"],
        "cv_mean_auroc": latent.get("cv_mean_auroc"),
        "angle_one_dimensional": int(latent["selected_lv"]) == 1,
        "pooled_spectrum": record["config_spectrum"].get("pooled"),
        "error": None,
    }


def _failed_arm(floor: int, exc: Exception) -> dict[str, Any]:
    return {
        "min_components": floor,
        "statistics": {},
        "p_values": {},
        "null_summary": {},
        "selected_lv": None,
        "component_selection": None,
        "cv_mean_auroc": None,
        "angle_one_dimensional": None,
        "pooled_spectrum": None,
        "error": str(exc),
    }


def analyze_arms(
    built: Mapping[str, Any], settings: Mapping[str, Any], *, n_jobs: int = 1, permutations: int | None = None
) -> dict[str, Any]:
    """Analyse one built dataset in both arms: same blocks, CV settings and RRPP seed, different floor."""

    from motco.analysis import TrajectoryAnalysisError, TrajectoryAnalysisParams, analyze_trajectories

    dataset = _trajectory_dataset(built, settings)
    arms: dict[str, Any] = {}
    for arm in ARMS:
        floor = ARM_FLOOR[arm]
        params = TrajectoryAnalysisParams(
            permutations=int(settings["permutations"] if permutations is None else permutations),
            seed=built["seeds"]["rrpp"],
            n_jobs=n_jobs,
            orientation_estimator=settings["orientation_estimator"],
            min_components=floor,
            **settings["cv"],
        )
        try:
            arms[arm] = _arm_record(analyze_trajectories(dataset, params), floor)
        except TrajectoryAnalysisError as exc:
            arms[arm] = _failed_arm(floor, exc)
    return arms


def evaluate_unit(
    mode: str, effect: float, replicate: int, settings: Mapping[str, Any], *, n_jobs: int = 1, reference: Any = None
) -> dict[str, Any]:
    start = perf_counter()
    built = build_unit(mode, effect, replicate, settings, reference=reference)
    arms = analyze_arms(built, settings, n_jobs=n_jobs)
    return d10._jsonable(
        {
            "key": unit_key(mode, effect, replicate),
            "mode": mode,
            "effect": float(effect),
            "replicate": replicate,
            "settings_hash": settings_hash(settings),
            "seeds": built["seeds"],
            "group_labels": built["group_labels"],
            "arms": arms,
            "runtime_seconds": perf_counter() - start,
        }
    )


# ── Calibration ──────────────────────────────────────────────────────────────


def _calibration_replicate(replicate: int, settings: Mapping[str, Any]) -> dict[str, Any]:
    """The anchor cell, ``floor_1`` arm, no permutations: the CV-selected rank and mean AUROC."""

    from motco.analysis import TrajectoryAnalysisParams, analyze_trajectories

    built = build_unit(ANCHOR, 0.0, replicate, settings)
    params = TrajectoryAnalysisParams(
        permutations=0,
        orientation_estimator=settings["orientation_estimator"],
        min_components=ARM_FLOOR[UNFLOORED_ARM],
        **settings["cv"],
    )
    latent = analyze_trajectories(_trajectory_dataset(built, settings), params).record["latent"]
    return {"replicate": replicate, "selected_lv": int(latent["selected_lv"]), "cv_mean_auroc": latent["cv_mean_auroc"]}


def calibration_summary(entries: Sequence[Mapping[str, Any]], threshold: float) -> dict[str, Any]:
    ranks = [int(entry["selected_lv"]) for entry in entries]
    aurocs = [float(entry["cv_mean_auroc"]) for entry in entries if entry.get("cv_mean_auroc") is not None]
    fraction = (sum(1 for rank in ranks if rank == 1) / len(ranks)) if ranks else None
    return {
        "statistic": CALIBRATION_STATISTIC,
        "threshold": float(threshold),
        "replicates": len(ranks),
        "rank1_fraction": fraction,
        "rank_distribution": rank_distribution(ranks),
        "mean_cv_auroc": float(np.mean(aurocs)) if aurocs else None,
        "met": fraction is not None and fraction >= threshold,
    }


def calibrate(out_dir: Path, settings: Mapping[str, Any], *, workers: int = 1) -> dict[str, Any]:
    """Run the calibration replicates, append the entry to ``calibration.json`` and return it."""

    replicates = list(range(int(settings["replicates"])))
    if workers > 1:
        with ProcessPoolExecutor(max_workers=workers) as pool:
            entries = list(pool.map(_calibration_replicate, replicates, [settings] * len(replicates)))
    else:
        entries = [_calibration_replicate(replicate, settings) for replicate in replicates]
    summary = calibration_summary(entries, float(settings["calibration"]["threshold"]))
    entry = {
        "strengths": {field: settings["generator"][field] for field in STRENGTH_FIELDS},
        "base_seed": settings["base_seed"],
        "cv": settings["cv"],
        "headroom": headroom_report(settings),
        **summary,
        "per_replicate": entries,
    }
    out_dir.mkdir(parents=True, exist_ok=True)
    path = out_dir / "calibration.json"
    history = json.loads(path.read_text(encoding="utf-8")) if path.exists() else []
    history.append(d10._jsonable(entry))
    path.write_text(json.dumps(history, indent=1, allow_nan=False) + "\n", encoding="utf-8")
    return entry


# ── Shards, resume, merge ────────────────────────────────────────────────────


def read_records(paths: Iterable[Path]) -> dict[str, dict[str, Any]]:
    """Records keyed by unit; torn lines (an interrupted write) are skipped and re-run."""

    records: dict[str, dict[str, Any]] = {}
    for path in paths:
        with path.open(encoding="utf-8") as handle:
            for line in handle:
                try:
                    record = json.loads(line)
                except json.JSONDecodeError:
                    continue
                if set(record.get("arms", {})) == set(ARMS):
                    records[str(record["key"])] = record
    return records


def check_headroom(settings: Mapping[str, Any], reference: Any = None) -> None:
    """Refuse settings whose pool-limited cells exceed the generator's headroom."""

    over = [key for key, headroom in headroom_report(settings, reference).items() if not headroom["fits"]]
    if over:
        raise SystemExit(
            f"These cells exceed the pool headroom at the settings' strengths: {over}. "
            "Lower the effects or raise p_dmp."
        )


def run_shard(
    out_dir: Path, settings: Mapping[str, Any], index: int, count: int, *, n_jobs: int = 1, limit: int | None = None
) -> int:
    """Evaluate this shard's missing units; returns how many were evaluated."""

    path = shard_path(out_dir, index, count)
    path.parent.mkdir(parents=True, exist_ok=True)
    own = read_records([path]) if path.exists() else {}
    with path.open("w", encoding="utf-8") as handle:
        for key in sorted(own):
            handle.write(json.dumps(own[key], sort_keys=True, allow_nan=False) + "\n")
    others = [other for other in sorted(path.parent.glob("shard_*.jsonl")) if other != path]
    done = {**read_records(others), **own}
    pending = [
        unit
        for position, unit in enumerate(units(settings))
        if position % count == index and unit_key(*unit) not in done
    ]
    if limit is not None:
        pending = pending[:limit]
    from motco.simulations.reference import load_reference

    reference = load_reference()
    check_headroom(settings, reference)
    for position, unit in enumerate(pending, start=1):
        record = evaluate_unit(*unit, settings, n_jobs=n_jobs, reference=reference)
        with path.open("a", encoding="utf-8") as handle:
            handle.write(json.dumps(record, sort_keys=True, allow_nan=False) + "\n")
            handle.flush()
        print(f"[shard {index}/{count}] {position}/{len(pending)} {record['key']}", file=sys.stderr, flush=True)
    return len(pending)


def p_value_rows(records: Mapping[str, Mapping[str, Any]]) -> list[dict[str, Any]]:
    rows = []
    for key in sorted(records):
        record = records[key]
        for arm in ARMS:
            p_values = record["arms"][arm]["p_values"]
            for statistic in STATISTICS:
                rows.append(
                    {
                        "mode": record["mode"],
                        "effect": record["effect"],
                        "replicate": record["replicate"],
                        "arm": arm,
                        "statistic": statistic,
                        "p_value": p_values.get(statistic),
                    }
                )
    return rows


REPLICATE_COLUMNS = (
    "mode",
    "effect",
    "replicate",
    "arm",
    "min_components",
    "selected_lv",
    "cv_mean_auroc",
    "angle_one_dimensional",
    "angle_null_q95",
    "angle_null_sd",
    "error",
)


def replicate_rows(records: Mapping[str, Mapping[str, Any]]) -> list[dict[str, Any]]:
    rows = []
    for key in sorted(records):
        record = records[key]
        for arm in ARMS:
            entry = record["arms"][arm]
            angle_null = (entry.get("null_summary") or {}).get("angle") or {}
            rows.append(
                {
                    "mode": record["mode"],
                    "effect": record["effect"],
                    "replicate": record["replicate"],
                    "arm": arm,
                    "min_components": entry.get("min_components"),
                    "selected_lv": entry.get("selected_lv"),
                    "cv_mean_auroc": entry.get("cv_mean_auroc"),
                    "angle_one_dimensional": entry.get("angle_one_dimensional"),
                    "angle_null_q95": angle_null.get("q95"),
                    "angle_null_sd": angle_null.get("sd"),
                    "error": entry.get("error"),
                }
            )
    return rows


def rank_rows(records: Mapping[str, Mapping[str, Any]]) -> list[dict[str, Any]]:
    """Per unit × arm: the selected rank (for the rank distributions and the calibration)."""

    return [
        {
            "mode": record["mode"],
            "effect": record["effect"],
            "replicate": record["replicate"],
            "arm": arm,
            "selected_lv": record["arms"][arm].get("selected_lv"),
        }
        for key in sorted(records)
        for record in [records[key]]
        for arm in ARMS
    ]


# ── Rates and decision (pure) ────────────────────────────────────────────────


anchor_bound = d10.anchor_bound


def rank_distribution(values: Iterable[Any]) -> dict[str, int]:
    counts: dict[int, int] = {}
    for value in values:
        if value is not None:
            counts[int(value)] = counts.get(int(value), 0) + 1
    return {str(rank): counts[rank] for rank in sorted(counts)}


def _cell_label(mode: str, effect: float) -> tuple[str, float]:
    return str(mode), float(effect)


def rejection_tallies(rows: Iterable[Mapping[str, Any]], *, alpha: float) -> dict[tuple, list[int]]:
    """(mode, effect, arm, statistic) → [available replicates, rejections]."""

    tallies: dict[tuple, list[int]] = {}
    for row in rows:
        p_value = row.get("p_value")
        if p_value is None or not math.isfinite(float(p_value)):
            continue
        key = (*_cell_label(row["mode"], row["effect"]), str(row["arm"]), str(row["statistic"]))
        tally = tallies.setdefault(key, [0, 0])
        tally[0] += 1
        tally[1] += int(float(p_value) <= alpha)
    return tallies


def _rate(tally: Sequence[int] | None) -> tuple[int, int, float | None, float | None]:
    replicates, rejections = tally if tally else (0, 0)
    if not replicates:
        return 0, 0, None, None
    rate = rejections / replicates
    return replicates, rejections, rate, math.sqrt(rate * (1.0 - rate) / replicates)


RATE_COLUMNS = (
    "mode",
    "effect",
    "statistic",
    "replicates_floor_1",
    "rejections_floor_1",
    "rate_floor_1",
    "se_floor_1",
    "replicates_floor_2",
    "rejections_floor_2",
    "rate_floor_2",
    "se_floor_2",
    "difference",
    "pooled_se",
    "bound",
    "verdict",
    "rank_distribution_floor_1",
    "rank_distribution_floor_2",
    "rank1_fraction_floor_1",
    "rank1_fraction_floor_2",
)


def rate_rows(
    rows: Iterable[Mapping[str, Any]],
    *,
    alpha: float,
    k: float,
    ranks: Iterable[Mapping[str, Any]] = (),
    cell_order: Sequence[tuple[str, float]] | None = None,
) -> list[dict[str, Any]]:
    """Per cell × statistic: both arms' rates, their difference and pooled SE, rank distributions.

    The ``floor_2`` arm carries a verdict at the control cells (anchor,
    translation); every other row is descriptive.
    """

    rows = list(rows)
    tallies = rejection_tallies(rows, alpha=alpha)
    rank_by_cell_arm: dict[tuple, list[Any]] = {}
    for entry in ranks:
        rank_by_cell_arm.setdefault((*_cell_label(entry["mode"], entry["effect"]), str(entry["arm"])), []).append(
            entry.get("selected_lv")
        )
    seen = {(*_cell_label(row["mode"], row["effect"]),) for row in rows}
    seen |= {(*_cell_label(entry["mode"], entry["effect"]),) for entry in ranks}
    order = list(cell_order) if cell_order is not None else sorted(seen)
    out = []
    for mode, effect in order:
        if (mode, effect) not in seen:
            continue
        for statistic in STATISTICS:
            r1, j1, rate1, se1 = _rate(tallies.get((mode, effect, "floor_1", statistic)))
            r2, j2, rate2, se2 = _rate(tallies.get((mode, effect, "floor_2", statistic)))
            difference = None if rate1 is None or rate2 is None else rate2 - rate1
            pooled = None if se1 is None or se2 is None else math.sqrt(se1**2 + se2**2)
            bound = anchor_bound(alpha, r2, k) if mode in CONTROL_MODES and r2 else None
            if mode not in CONTROL_MODES:
                verdict = "descriptive"
            elif rate2 is None:
                verdict = "unavailable"
            else:
                verdict = "holds" if rate2 <= bound else "exceeds"
            ranks1 = [v for v in rank_by_cell_arm.get((mode, effect, "floor_1"), []) if v is not None]
            ranks2 = [v for v in rank_by_cell_arm.get((mode, effect, "floor_2"), []) if v is not None]
            out.append(
                {
                    "mode": mode,
                    "effect": effect,
                    "statistic": statistic,
                    "replicates_floor_1": r1,
                    "rejections_floor_1": j1,
                    "rate_floor_1": rate1,
                    "se_floor_1": se1,
                    "replicates_floor_2": r2,
                    "rejections_floor_2": j2,
                    "rate_floor_2": rate2,
                    "se_floor_2": se2,
                    "difference": difference,
                    "pooled_se": pooled,
                    "bound": bound,
                    "verdict": verdict,
                    "rank_distribution_floor_1": json.dumps(rank_distribution(ranks1)),
                    "rank_distribution_floor_2": json.dumps(rank_distribution(ranks2)),
                    "rank1_fraction_floor_1": (sum(1 for v in ranks1 if int(v) == 1) / len(ranks1)) if ranks1 else None,
                    "rank1_fraction_floor_2": (sum(1 for v in ranks2 if int(v) == 1) / len(ranks2)) if ranks2 else None,
                }
            )
    return out


def rank_floor_decision(
    rows: Iterable[Mapping[str, Any]],
    *,
    alpha: float,
    k: float,
    ranks: Iterable[Mapping[str, Any]] = (),
    calibration_threshold: float = 0.5,
) -> dict[str, Any]:
    """``adopt_floor`` when the ``floor_2`` arm holds the one-sided anchor bound at both controls.

    ``rows`` are per unit × arm × statistic p-values (``mode``, ``effect``,
    ``replicate``, ``arm``, ``statistic``, ``p_value``); ``ranks`` are per
    unit × arm selected ranks, which feed the calibration evidence (the
    ``floor_1`` anchor's rank-1 fraction) and the rank distributions. The
    ``floor_1`` arm is reported and never enters the decision.
    """

    rows = list(rows)
    ranks = list(ranks)
    tallies = rejection_tallies(rows, alpha=alpha)
    controls: dict[str, dict[str, dict[str, Any]]] = {}
    verdicts: dict[str, dict[str, str]] = {}
    exceeding: list[dict[str, Any]] = []
    for mode in CONTROL_MODES:
        effects = sorted({float(row["effect"]) for row in rows if str(row["mode"]) == mode})
        effects += sorted({float(e["effect"]) for e in ranks if str(e["mode"]) == mode} - set(effects))
        effect = effects[0] if effects else (0.0 if mode == ANCHOR else None)
        controls[mode] = {"effect": effect}
        verdicts[mode] = {}
        for arm in ARMS:
            controls[mode][arm] = {}
            for statistic in STATISTICS:
                tally = tallies.get((mode, effect, arm, statistic)) if effect is not None else None
                replicates, rejections, rate, se = _rate(tally)
                bound = anchor_bound(alpha, replicates, k) if replicates else None
                if arm != JUDGED_ARM:
                    verdict = "descriptive"
                elif rate is None:
                    verdict = "unavailable"
                else:
                    verdict = "holds" if rate <= bound else "exceeds"
                controls[mode][arm][statistic] = {
                    "replicates": replicates,
                    "rejections": rejections,
                    "rate": rate,
                    "monte_carlo_se": se,
                    "bound": bound,
                    "verdict": verdict,
                }
                if arm == JUDGED_ARM:
                    verdicts[mode][statistic] = verdict
                    if verdict == "exceeds":
                        exceeding.append({"mode": mode, "statistic": statistic, "rate": rate, "bound": bound})
    admissible = all(v == "holds" for per_mode in verdicts.values() for v in per_mode.values())
    anchor_ranks = [
        e.get("selected_lv") for e in ranks if str(e["mode"]) == ANCHOR and str(e["arm"]) == UNFLOORED_ARM
    ]
    anchor_ranks = [int(v) for v in anchor_ranks if v is not None]
    fraction = (sum(1 for v in anchor_ranks if v == 1) / len(anchor_ranks)) if anchor_ranks else None
    return {
        "rule": (
            f"{JUDGED_ARM}-arm rate (p <= alpha) <= alpha + k * sqrt(alpha * (1 - alpha) / R) for delta, angle "
            "and shape at the zero-effect anchor and the translation control; adopt_floor when every "
            "statistic holds, otherwise diagnose_floor_excess"
        ),
        "alpha": alpha,
        "k": k,
        "judged_arm": JUDGED_ARM,
        "controls": controls,
        "verdicts": verdicts,
        "exceeding": exceeding,
        "decision": "adopt_floor" if admissible else "diagnose_floor_excess",
        "floor_1_role": (
            "descriptive (the unfloored 0.9.0 rule; its angle at rank 1 is a two-point null); "
            "not part of the decision"
        ),
        "calibration": {
            "statistic": CALIBRATION_STATISTIC,
            "threshold": calibration_threshold,
            "replicates": len(anchor_ranks),
            "rank1_fraction": fraction,
            "rank_distribution": rank_distribution(anchor_ranks),
            "met": fraction is not None and fraction >= calibration_threshold,
        },
    }


def write_decision(
    out_dir: Path,
    rows: Sequence[Mapping[str, Any]],
    settings: Mapping[str, Any],
    records: Mapping[str, Mapping[str, Any]],
) -> dict[str, Any]:
    alpha, k = float(settings["alpha"]), float(settings["k"])
    ranks = rank_rows(records)
    decision = rank_floor_decision(
        rows, alpha=alpha, k=k, ranks=ranks, calibration_threshold=float(settings["calibration"]["threshold"])
    )
    rates = rate_rows(rows, alpha=alpha, k=k, ranks=ranks, cell_order=cells(settings))
    decision["power"] = [
        {key: row[key] for key in RATE_COLUMNS}
        for row in rates
        if row["mode"] not in CONTROL_MODES
    ]
    decision["selected_rank"] = {
        arm: {
            unit_key(mode, effect, 0)[:-4] + "]": rank_distribution(
                r["arms"][arm].get("selected_lv")
                for r in records.values()
                if r["mode"] == mode and float(r["effect"]) == float(effect)
            )
            for mode, effect in cells(settings)
        }
        for arm in ARMS
    }
    decision["one_dimensional_angle"] = {
        arm: sum(1 for r in records.values() if r["arms"][arm].get("angle_one_dimensional")) for arm in ARMS
    }
    decision["failures"] = {
        arm: sum(1 for r in records.values() if r["arms"][arm].get("error") is not None) for arm in ARMS
    }
    decision["settings_hash"] = settings_hash(settings)
    decision["construction"] = {
        "generator": settings["generator"],
        "strengths": {field: settings["generator"][field] for field in STRENGTH_FIELDS},
        "phase6_strengths": dict(PHASE6_STRENGTHS),
        "layers": settings["layers"],
        "cv": settings["cv"],
        "arms": settings["arms"],
        "effects": settings["effects"],
        "controls": settings["controls"],
        "permutations": settings["permutations"],
        "base_seed": settings["base_seed"],
    }
    (out_dir / "rank_floor_decision.json").write_text(
        json.dumps(d10._jsonable(decision), indent=2, allow_nan=False) + "\n", encoding="utf-8"
    )
    d10._write_csv(out_dir / "rank_floor_rates.csv", rates, RATE_COLUMNS)
    return decision


def merge(out_dir: Path, settings: Mapping[str, Any]) -> dict[str, Any]:
    """Merge every shard, then write the tables, the rates and the decision."""

    records = read_records(sorted((out_dir / "records").glob("shard_*.jsonl")))
    expected = [unit_key(*unit) for unit in units(settings)]
    with (out_dir / "merged.jsonl").open("w", encoding="utf-8") as handle:
        for key in sorted(records):
            handle.write(json.dumps(records[key], sort_keys=True, allow_nan=False) + "\n")
    rows = p_value_rows(records)
    d10._write_csv(
        out_dir / "rank_floor_pvalues.csv", rows, ("mode", "effect", "replicate", "arm", "statistic", "p_value")
    )
    d10._write_csv(out_dir / "rank_floor_replicates.csv", replicate_rows(records), REPLICATE_COLUMNS)
    decision = write_decision(out_dir, rows, settings, records)
    return {
        "units": len(records),
        "expected": len(expected),
        "missing": [key for key in expected if key not in records][:20],
        "decision": decision["decision"],
        "verdicts": decision["verdicts"],
        "exceeding": decision["exceeding"],
        "controls": {
            mode: {arm: {s: decision["controls"][mode][arm][s]["rate"] for s in STATISTICS} for arm in ARMS}
            for mode in CONTROL_MODES
        },
        "bound": {mode: decision["controls"][mode][JUDGED_ARM]["delta"]["bound"] for mode in CONTROL_MODES},
        "calibration": {
            "rank1_fraction": decision["calibration"]["rank1_fraction"],
            "met": decision["calibration"]["met"],
        },
        "one_dimensional_angle": decision["one_dimensional_angle"],
        "failures": decision["failures"],
    }


# ── CLI ──────────────────────────────────────────────────────────────────────


def add_construction_arguments(parser: argparse.ArgumentParser) -> None:
    parser.add_argument("--out-dir", type=Path, required=True, help="Output directory")
    parser.add_argument("--reps", type=int, default=500, help="Replicates (default: 500)")
    parser.add_argument("--perms", type=int, default=999, help="RRPP permutations (default: 999)")
    parser.add_argument("--base-seed", type=int, default=2026, help="Parent seed (default: 2026)")
    parser.add_argument("--alpha", type=float, default=0.05, help="Rejection threshold, p <= alpha (default: 0.05)")
    parser.add_argument("--k", type=float, default=2.0, help="Bound multiplier k (default: 2)")
    for field, default in PHASE6_STRENGTHS.items():
        parser.add_argument(
            f"--{field.replace('_', '-')}",
            dest=field,
            type=float,
            default=default,
            help=f"Generator {field} (Phase 6: {default}); weaken to reach the rank-1 calibration",
        )
    parser.add_argument(
        "--magnitude-effects", type=str, default=",".join(map(str, DEFAULT_EFFECTS["magnitude"])),
        help="Comma-separated magnitude effects (default: Phase 6 levels below saturation)",
    )
    parser.add_argument(
        "--orientation-effects", type=str, default=",".join(map(str, DEFAULT_EFFECTS["orientation"])),
        help="Comma-separated orientation effects (default: Phase 6 levels)",
    )
    parser.add_argument(
        "--shape-effects", type=str, default=",".join(map(str, DEFAULT_EFFECTS["shape"])),
        help="Comma-separated shape effects (default: Phase 6 levels)",
    )
    parser.add_argument(
        "--translation-effect", type=float, default=DEFAULT_TRANSLATION_EFFECT,
        help="Translation control effect (default: 1.0)",
    )
    parser.add_argument(
        "--calibration-threshold", type=float, default=0.5,
        help="Rank-1 fraction of the unfloored anchor arm required before a paper-grade run (default: 0.5)",
    )


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Rank-floor check at a rank-1-selecting Phase 6 construction (no R).")
    sub = parser.add_subparsers(dest="command", required=True)
    run = sub.add_parser("run", help="Evaluate (a shard of) the units; an unsharded run also merges")
    add_construction_arguments(run)
    run.add_argument("--shard", type=parse_shard, default=(0, 1), help="Run shard i of N (i/N; default 0/1)")
    run.add_argument("--n-jobs", type=int, default=1, help="Jobs for CV and RRPP (not part of the settings)")
    run.add_argument("--limit", type=int, default=None, help="Evaluate at most this many pending units")
    cal = sub.add_parser("calibrate", help="Anchor cell, floor_1 arm, no permutations: the rank-1 fraction")
    add_construction_arguments(cal)
    cal.add_argument("--workers", type=int, default=1, help="Replicates evaluated in parallel")
    merge_parser = sub.add_parser("merge", help="Merge every shard and write the decision (reads settings.json)")
    merge_parser.add_argument("--out-dir", type=Path, required=True, help="Output directory")
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    if args.command == "merge":
        settings = load_settings(args.out_dir)
    elif args.command == "calibrate":
        settings = build_settings(args)
        entry = calibrate(args.out_dir, settings, workers=args.workers)
        print(json.dumps({key: value for key, value in entry.items() if key != "per_replicate"}, indent=1))
        return 0
    else:
        settings = build_settings(args)
        ensure_settings(args.out_dir, settings)
        index, count = args.shard
        evaluated = run_shard(args.out_dir, settings, index, count, n_jobs=args.n_jobs, limit=args.limit)
        print(f"evaluated {evaluated} unit(s) into {shard_path(args.out_dir, index, count)}", file=sys.stderr)
        if count > 1 or args.limit is not None:
            return 0
    print(json.dumps(merge(args.out_dir, settings), indent=1))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
