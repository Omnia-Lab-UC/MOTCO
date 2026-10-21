#!/usr/bin/env python3
"""Latent-adjustment Type I check at the Phase 6 SEA-AD design (no R).

The gate behind ``openspec/changes/covariate-aware-latent-construction/`` (design
D8, spec ``latent-adjustment-validation``). The case study may report
latent-adjusted ``motco analyze`` results only after this check names the
conditioning it uses as admissible.

Each **replicate** reuses the covariate-adjustment Type I construction (D10,
:mod:`covariate_type_i_check`, imported, not copied): a null semi-synthetic
dataset at the Phase 6 small-n design, a numeric covariate ``z`` and a
three-level categorical covariate confounded with group × stage, and their
effect injected into the measured features. The seeds are D10's, so at equal
strengths a replicate's blocks and covariates equal the D10 replicate's. On top
of that it declares ``extra_covariates`` pure-noise numeric covariates
(``x1``, ``x2``, …; i.i.d. N(0, 1), never injected), so the residualization
runs at the case study's 8–10 covariates and the inner CV folds at its
degrees of freedom.

Each replicate is analysed with :func:`motco.analysis.analyze_trajectories`
three times, with every covariate declared, the same CV ``random_state`` and the
same RRPP seed:

- ``baseline``: covariates in the trajectory model only (the 0.8.0 basis);
- ``stage``: ``covariate_latent_adjustment="stage"``;
- ``group_stage``: ``covariate_latent_adjustment="group_stage"``.

Each arm records its p-values, selected rank, the Spearman correlation of its
first latent component with ``z`` and, for the residualized arms, the record's
``latent_adjustment`` block. An arm the analysis rejects (for example a CV
training fold missing a factor level) is recorded with its error and is
unavailable for that replicate.

**Calibration.** The injection must make the covariate-blind basis spend its
first component on the nuisance axis: the baseline arm's median ``|ρ|`` must be
at least ``calibration_threshold`` (0.8) on the smoke profile before any
paper-grade run; otherwise ``--s-z`` is raised and recorded. At D10's
``s_z = 3`` the median was 0.64 over the smoke replicates (2026-10-20); it
rises with ``s_z`` and saturates (0.74 at 6, 0.78 at 10, 0.80 at 20, 0.81 at
30, 0.815 at 40), because the covariate-blind fit re-standardizes every feature
and only the injected fraction of features can carry ``z``. The default is
therefore ``s_z = 40``; ``s_c`` stays at D10's 2.

The **decision** (pure, :func:`latent_adjustment_type_i_decision`): per
residualized arm and statistic, the rejection rate (``p ≤ alpha``) over its R
available replicates must be at or below ``alpha + k·√(alpha(1 − alpha)/R)``.
An arm holding every statistic is admissible. ``adopt_stage`` when ``stage`` is
admissible, ``adopt_group_stage`` when only ``group_stage`` is, otherwise
``diagnose_latent_excess``. The baseline arm is reported with the calibration
evidence and never enters the decision.

Outputs (``--out-dir``):

- ``settings.json``: settings and their hash; a resume into a directory
  written with different settings is refused;
- ``records/shard_<i>_of_<N>.jsonl``: one line per replicate (all arms);
  resume skips recorded replicates;
- after ``merge`` (automatic for an unsharded run): ``merged.jsonl``,
  ``latent_adjustment_type_i_pvalues.csv``,
  ``latent_adjustment_type_i_replicates.csv`` (per replicate × arm: rank,
  first-component correlation, residual df, error),
  ``latent_adjustment_type_i_rates.csv`` and
  ``latent_adjustment_type_i_decision.json``.

Examples:

    # smoke profile (calibration)
    python scripts/latent_adjustment_type_i_check.py run --reps 20 --perms 99 --out-dir /tmp/latent-smoke

    # one SLURM array task, then the merge
    python scripts/latent_adjustment_type_i_check.py run --shard 3/50 --out-dir results/latent
    python scripts/latent_adjustment_type_i_check.py merge --out-dir results/latent
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import sys
from pathlib import Path
from time import perf_counter
from typing import Any, Iterable, Mapping, Sequence

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import covariate_type_i_check as d10  # noqa: E402  (the D10 construction, sharding and settings helpers)

sys.path.pop(0)

SCRIPT_VERSION = 1
CHECK = "latent_adjustment_type_i"
STATISTICS: tuple[str, ...] = d10.STATISTICS
ARMS: tuple[str, ...] = ("baseline", "stage", "group_stage")
#: Residualized arms in decision preference order.
RESIDUALIZED_ARMS: tuple[str, ...] = ("stage", "group_stage")
ARM_SETTING: dict[str, str | None] = {"baseline": None, "stage": "stage", "group_stage": "group_stage"}
CALIBRATION_STATISTIC = "Spearman correlation of the first latent component's scores with the numeric covariate z"
_QUANTILES: tuple[tuple[str, float], ...] = (
    ("min", 0.0),
    ("q05", 0.05),
    ("q25", 0.25),
    ("q50", 0.5),
    ("q75", 0.75),
    ("q95", 0.95),
    ("max", 1.0),
)


# ── Settings ─────────────────────────────────────────────────────────────────


def build_settings(args: argparse.Namespace) -> dict[str, Any]:
    settings = d10.build_settings(args)
    settings["script_version"] = SCRIPT_VERSION
    settings["check"] = CHECK
    settings["arms"] = {arm: ARM_SETTING[arm] for arm in ARMS}
    settings["covariates"]["extra_numeric"] = {
        "count": int(args.extra_covariates),
        "distribution": "iid N(0, 1), not injected",
    }
    settings["calibration"] = {"statistic": CALIBRATION_STATISTIC, "threshold": float(args.calibration_threshold)}
    return settings


settings_hash = d10.settings_hash
ensure_settings = d10.ensure_settings
load_settings = d10.load_settings
parse_shard = d10.parse_shard
shard_path = d10.shard_path


def _seed(*parts: Any) -> int:
    encoded = json.dumps([CHECK, *parts], separators=(",", ":")).encode("utf-8")
    return int(hashlib.sha256(encoded).hexdigest()[:8], 16) & 0x7FFFFFFF


def extra_covariate_names(settings: Mapping[str, Any]) -> list[str]:
    return [f"x{index + 1}" for index in range(int(settings["covariates"]["extra_numeric"]["count"]))]


def declared_covariates(settings: Mapping[str, Any]) -> list[str]:
    """Every covariate the arms declare, in model order."""

    return ["z", "category", *extra_covariate_names(settings)]


# ── Replicate evaluation ─────────────────────────────────────────────────────


def build_replicate(replicate: int, settings: Mapping[str, Any], reference: Any = None) -> dict[str, Any]:
    """The D10 replicate plus the pure-noise covariates (drawn from their own seed)."""

    built = d10.build_replicate(replicate, settings, reference=reference)
    names = extra_covariate_names(settings)
    if names:
        seed = _seed(int(settings["base_seed"]), "extra_covariates", replicate)
        noise = np.random.default_rng(seed).normal(size=(built["metadata"].shape[0], len(names)))
        built["metadata"] = built["metadata"].assign(**{name: noise[:, j] for j, name in enumerate(names)})
        built["seeds"] = {**built["seeds"], "extra_covariates": seed}
    return built


def _spearman(a: np.ndarray, b: np.ndarray) -> float:
    from scipy.stats import spearmanr

    return float(spearmanr(a, b).statistic)


def _arm_record(result: Any, z: np.ndarray) -> dict[str, Any]:
    record = result.record
    scores = result.latent_scores.to_numpy(dtype=float)
    return {
        "statistics": record["statistics"],
        "p_values": record["p_values"],
        "selected_lv": record["latent"]["selected_lv"],
        "component_selection": record["latent"]["component_selection"],
        "first_component_spearman_z": _spearman(scores[:, 0], z),
        "latent_adjustment": record.get("latent_adjustment"),
        "error": None,
    }


def analyze_arms(built: Mapping[str, Any], settings: Mapping[str, Any], *, n_jobs: int = 1) -> dict[str, Any]:
    """Analyse one built dataset in every arm: same covariates, CV settings and RRPP seed."""

    from motco.analysis import (
        OmicsBlock,
        TrajectoryAnalysisError,
        TrajectoryAnalysisParams,
        TrajectoryDataset,
        analyze_trajectories,
    )

    metadata = built["metadata"]
    dataset = TrajectoryDataset(
        blocks=[OmicsBlock(name, built["blocks"][name], "identity") for name in settings["layers"]],
        metadata=metadata,
        sample_col="sample_id",
        group_col="group",
        stage_col="stage",
        stage_order=[str(stage) for stage in sorted(set(metadata["stage"].tolist()))],
        covariate_cols=declared_covariates(settings),
        categorical_covariates=["category"],
    )
    z = metadata["z"].to_numpy(dtype=float)
    arms: dict[str, Any] = {}
    for arm in ARMS:
        params = TrajectoryAnalysisParams(
            permutations=int(settings["permutations"]),
            seed=built["seeds"]["rrpp"],
            n_jobs=n_jobs,
            orientation_estimator=settings["orientation_estimator"],
            covariate_latent_adjustment=ARM_SETTING[arm],
            **settings["cv"],
        )
        try:
            arms[arm] = _arm_record(analyze_trajectories(dataset, params), z)
        except TrajectoryAnalysisError as exc:
            arms[arm] = {
                "statistics": {},
                "p_values": {},
                "selected_lv": None,
                "component_selection": None,
                "first_component_spearman_z": None,
                "latent_adjustment": None,
                "error": str(exc),
            }
    return arms


def evaluate_replicate(
    replicate: int, settings: Mapping[str, Any], *, n_jobs: int = 1, reference: Any = None
) -> dict[str, Any]:
    start = perf_counter()
    built = build_replicate(replicate, settings, reference=reference)
    arms = analyze_arms(built, settings, n_jobs=n_jobs)
    return d10._jsonable(
        {
            "replicate": replicate,
            "settings_hash": settings_hash(settings),
            "seeds": built["seeds"],
            "group_labels": built["group_labels"],
            "covariate_cells": d10.covariate_cell_summary(built["metadata"], built["group_labels"]),
            "arms": arms,
            "runtime_seconds": perf_counter() - start,
        }
    )


# ── Shards, resume, merge ────────────────────────────────────────────────────


def read_records(paths: Iterable[Path]) -> dict[int, dict[str, Any]]:
    """Records keyed by replicate; torn lines (an interrupted write) are skipped and re-run."""

    records: dict[int, dict[str, Any]] = {}
    for path in paths:
        with path.open(encoding="utf-8") as handle:
            for line in handle:
                try:
                    record = json.loads(line)
                except json.JSONDecodeError:
                    continue
                if set(record.get("arms", {})) == set(ARMS):
                    records[int(record["replicate"])] = record
    return records


def run_shard(
    out_dir: Path, settings: Mapping[str, Any], index: int, count: int, *, n_jobs: int = 1, limit: int | None = None
) -> int:
    """Evaluate this shard's missing replicates; returns how many were evaluated.

    A replicate recorded in any shard file under ``records/`` counts as done;
    ``merge`` keeps one record per replicate.
    """

    path = shard_path(out_dir, index, count)
    path.parent.mkdir(parents=True, exist_ok=True)
    own = read_records([path]) if path.exists() else {}
    with path.open("w", encoding="utf-8") as handle:
        for replicate in sorted(own):
            handle.write(json.dumps(own[replicate], sort_keys=True, allow_nan=False) + "\n")
    others = [other for other in sorted(path.parent.glob("shard_*.jsonl")) if other != path]
    done = {**read_records(others), **own}
    pending = [r for r in range(int(settings["replicates"])) if r % count == index and r not in done]
    if limit is not None:
        pending = pending[:limit]
    from motco.simulations.reference import load_reference

    reference = load_reference()
    for position, replicate in enumerate(pending, start=1):
        record = evaluate_replicate(replicate, settings, n_jobs=n_jobs, reference=reference)
        with path.open("a", encoding="utf-8") as handle:
            handle.write(json.dumps(record, sort_keys=True, allow_nan=False) + "\n")
            handle.flush()
        print(f"[shard {index}/{count}] {position}/{len(pending)} replicate {replicate}", file=sys.stderr, flush=True)
    return len(pending)


def p_value_rows(records: Mapping[int, Mapping[str, Any]]) -> list[dict[str, Any]]:
    rows = []
    for replicate in sorted(records):
        for arm in ARMS:
            p_values = records[replicate]["arms"][arm]["p_values"]
            for statistic in STATISTICS:
                rows.append(
                    {"replicate": replicate, "arm": arm, "statistic": statistic, "p_value": p_values.get(statistic)}
                )
    return rows


REPLICATE_COLUMNS = (
    "replicate",
    "arm",
    "selected_lv",
    "first_component_spearman_z",
    "residual_df",
    "cv_min_training_df",
    "error",
)


def replicate_rows(records: Mapping[int, Mapping[str, Any]]) -> list[dict[str, Any]]:
    rows = []
    for replicate in sorted(records):
        for arm in ARMS:
            entry = records[replicate]["arms"][arm]
            block = entry.get("latent_adjustment") or {}
            rows.append(
                {
                    "replicate": replicate,
                    "arm": arm,
                    "selected_lv": entry.get("selected_lv"),
                    "first_component_spearman_z": entry.get("first_component_spearman_z"),
                    "residual_df": block.get("residual_df"),
                    "cv_min_training_df": block.get("cv_min_training_df"),
                    "error": entry.get("error"),
                }
            )
    return rows


# ── Decision (pure) ──────────────────────────────────────────────────────────


RATE_COLUMNS = ("arm", "statistic", "replicates", "rejections", "rate", "bound", "verdict")
anchor_bound = d10.anchor_bound


def rejection_rates(rows: Iterable[Mapping[str, Any]], *, alpha: float, k: float) -> list[dict[str, Any]]:
    """Per arm × statistic: available replicates R, rejections (``p ≤ alpha``), rate and bound.

    Only the residualized arms carry a verdict; the baseline arm is descriptive.
    """

    tallies: dict[tuple[str, str], list[int]] = {}
    for row in rows:
        p_value = row.get("p_value")
        if p_value is None or not math.isfinite(float(p_value)):
            continue
        tally = tallies.setdefault((str(row["arm"]), str(row["statistic"])), [0, 0])
        tally[0] += 1
        tally[1] += int(float(p_value) <= alpha)
    out = []
    for arm in ARMS:
        for statistic in STATISTICS:
            replicates, rejections = tallies.get((arm, statistic), [0, 0])
            rate = rejections / replicates if replicates else None
            bound = anchor_bound(alpha, replicates, k) if replicates else None
            if arm not in RESIDUALIZED_ARMS:
                verdict = "descriptive"
            elif rate is None or bound is None:
                verdict = "unavailable"
            else:
                verdict = "holds" if rate <= bound else "exceeds"
            out.append(
                {
                    "arm": arm,
                    "statistic": statistic,
                    "replicates": replicates,
                    "rejections": rejections,
                    "rate": rate,
                    "bound": bound,
                    "verdict": verdict,
                }
            )
    return out


def correlation_summary(values: Iterable[float | None]) -> dict[str, Any]:
    """Distribution of ``|ρ|`` over the replicates that report one."""

    magnitudes = np.array([abs(float(v)) for v in values if v is not None and math.isfinite(float(v))])
    if magnitudes.size == 0:
        return {"replicates": 0, "median_abs": None, "quantiles_abs": {}}
    return {
        "replicates": int(magnitudes.size),
        "median_abs": float(np.median(magnitudes)),
        "quantiles_abs": {name: float(np.quantile(magnitudes, q)) for name, q in _QUANTILES},
    }


def latent_adjustment_type_i_decision(
    rows: Iterable[Mapping[str, Any]],
    *,
    alpha: float,
    k: float,
    correlations: Mapping[str, Sequence[float | None]] | None = None,
    calibration_threshold: float = 0.8,
) -> dict[str, Any]:
    """Admissible conditionings and the adoption decision from the residualized arms' rates.

    ``correlations`` (per arm, the per-replicate first-component Spearman
    correlations with ``z``) only feeds the reported calibration evidence.
    """

    rates = rejection_rates(rows, alpha=alpha, k=k)
    by_arm: dict[str, dict[str, Mapping[str, Any]]] = {arm: {} for arm in ARMS}
    for row in rates:
        by_arm[row["arm"]][row["statistic"]] = row
    verdicts = {arm: {s: by_arm[arm][s]["verdict"] for s in STATISTICS} for arm in RESIDUALIZED_ARMS}
    admissible = [arm for arm in RESIDUALIZED_ARMS if all(v == "holds" for v in verdicts[arm].values())]
    if "stage" in admissible:
        decision = "adopt_stage"
    elif "group_stage" in admissible:
        decision = "adopt_group_stage"
    else:
        decision = "diagnose_latent_excess"
    correlations = correlations or {}
    evidence = {arm: correlation_summary(correlations.get(arm, ())) for arm in ARMS}
    baseline_median = evidence["baseline"]["median_abs"]
    return {
        "rule": (
            "residualized-arm rate (p <= alpha) <= alpha + k * sqrt(alpha * (1 - alpha) / R) for every statistic; "
            "adopt_stage if stage is admissible, adopt_group_stage if only group_stage is, "
            "otherwise diagnose_latent_excess"
        ),
        "alpha": alpha,
        "k": k,
        "replicates": {arm: {s: by_arm[arm][s]["replicates"] for s in STATISTICS} for arm in ARMS},
        "bound": {arm: {s: by_arm[arm][s]["bound"] for s in STATISTICS} for arm in ARMS},
        "rates": {arm: {s: by_arm[arm][s]["rate"] for s in STATISTICS} for arm in ARMS},
        "verdicts": verdicts,
        "admissible": admissible,
        "decision": decision,
        "baseline_role": "descriptive (covariates in the trajectory model only, covariate-blind basis); "
        "not part of the decision",
        "calibration": {
            "statistic": CALIBRATION_STATISTIC,
            "threshold": calibration_threshold,
            "baseline_median_abs": baseline_median,
            "met": baseline_median is not None and baseline_median >= calibration_threshold,
            "by_arm": evidence,
        },
    }


def write_decision(
    out_dir: Path,
    rows: Sequence[Mapping[str, Any]],
    settings: Mapping[str, Any],
    records: Mapping[int, Mapping[str, Any]],
) -> dict[str, Any]:
    alpha, k = float(settings["alpha"]), float(settings["k"])
    correlations = {
        arm: [records[r]["arms"][arm].get("first_component_spearman_z") for r in sorted(records)] for arm in ARMS
    }
    decision = latent_adjustment_type_i_decision(
        rows,
        alpha=alpha,
        k=k,
        correlations=correlations,
        calibration_threshold=float(settings["calibration"]["threshold"]),
    )
    decision["selected_rank"] = {
        arm: _rank_distribution(records[r]["arms"][arm].get("selected_lv") for r in records) for arm in ARMS
    }
    decision["failures"] = {
        arm: sum(1 for r in records if records[r]["arms"][arm].get("error") is not None) for arm in ARMS
    }
    decision["settings_hash"] = settings_hash(settings)
    decision["construction"] = {
        "generator": settings["generator"],
        "layers": settings["layers"],
        "covariates": settings["covariates"],
        "declared_covariates": declared_covariates(settings),
        "injection": settings["injection"],
        "permutations": settings["permutations"],
        "base_seed": settings["base_seed"],
        "cv": settings["cv"],
    }
    (out_dir / "latent_adjustment_type_i_decision.json").write_text(
        json.dumps(d10._jsonable(decision), indent=2, allow_nan=False) + "\n", encoding="utf-8"
    )
    d10._write_csv(
        out_dir / "latent_adjustment_type_i_rates.csv", rejection_rates(rows, alpha=alpha, k=k), RATE_COLUMNS
    )
    return decision


def _rank_distribution(values: Iterable[Any]) -> dict[str, int]:
    counts: dict[int, int] = {}
    for value in values:
        if value is not None:
            counts[int(value)] = counts.get(int(value), 0) + 1
    return {str(rank): counts[rank] for rank in sorted(counts)}


def merge(out_dir: Path, settings: Mapping[str, Any]) -> dict[str, Any]:
    """Merge every shard, then write the tables, the rates and the decision."""

    records = read_records(sorted((out_dir / "records").glob("shard_*.jsonl")))
    expected = int(settings["replicates"])
    with (out_dir / "merged.jsonl").open("w", encoding="utf-8") as handle:
        for replicate in sorted(records):
            handle.write(json.dumps(records[replicate], sort_keys=True, allow_nan=False) + "\n")
    rows = p_value_rows(records)
    d10._write_csv(
        out_dir / "latent_adjustment_type_i_pvalues.csv", rows, ("replicate", "arm", "statistic", "p_value")
    )
    d10._write_csv(out_dir / "latent_adjustment_type_i_replicates.csv", replicate_rows(records), REPLICATE_COLUMNS)
    decision = write_decision(out_dir, rows, settings, records)
    return {
        "replicates": len(records),
        "expected": expected,
        "missing": sorted(set(range(expected)) - set(records))[:20],
        "decision": decision["decision"],
        "admissible": decision["admissible"],
        "rates": decision["rates"],
        "bound": decision["bound"]["stage"],
        "calibration": {
            "baseline_median_abs": decision["calibration"]["baseline_median_abs"],
            "met": decision["calibration"]["met"],
        },
        "selected_rank": decision["selected_rank"],
        "failures": decision["failures"],
    }


# ── CLI ──────────────────────────────────────────────────────────────────────


def add_construction_arguments(run: argparse.ArgumentParser) -> None:
    """The run flags shared with the power comparison: construction, strengths, sharding."""

    run.add_argument("--out-dir", type=Path, required=True, help="Output directory")
    run.add_argument("--reps", type=int, default=500, help="Replicates (default: 500)")
    run.add_argument("--perms", type=int, default=999, help="RRPP permutations (default: 999)")
    run.add_argument("--base-seed", type=int, default=2026, help="Parent seed, shared with D10 (default: 2026)")
    run.add_argument("--alpha", type=float, default=0.05, help="Rejection threshold, p <= alpha (default: 0.05)")
    run.add_argument("--k", type=float, default=2.0, help="Bound multiplier k (default: 2)")
    run.add_argument("--a-g", type=float, default=0.5, help="Covariate group effect a_g (default: 0.5)")
    run.add_argument("--a-s", type=float, default=1.0, help="Covariate stage effect a_s (default: 1.0)")
    run.add_argument("--a-gs", type=float, default=2.0, help="Covariate group x stage effect a_gs (default: 2.0)")
    run.add_argument(
        "--categorical-logit-scale", type=float, default=0.5, help="Categorical logit scale c (default: 0.5)"
    )
    run.add_argument(
        "--s-z", type=float, default=40.0, help="Numeric covariate injection strength (default: 40, calibrated)"
    )
    run.add_argument("--s-c", type=float, default=2.0, help="Categorical covariate injection strength (default: 2.0)")
    run.add_argument(
        "--feature-fraction", type=float, default=0.2, help="Fraction of features per direction (default: 0.2)"
    )
    run.add_argument(
        "--extra-covariates",
        type=int,
        default=6,
        help="Pure-noise numeric covariates declared beside z and the factor (default: 6, i.e. 8 covariates)",
    )
    run.add_argument("--shard", type=parse_shard, default=(0, 1), help="Run shard i of N (i/N; default 0/1)")
    run.add_argument("--n-jobs", type=int, default=1, help="Jobs for CV and RRPP (not part of the settings)")
    run.add_argument("--limit", type=int, default=None, help="Evaluate at most this many pending units")


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Latent-adjustment Type I check at the Phase 6 design (no R).")
    sub = parser.add_subparsers(dest="command", required=True)
    run = sub.add_parser("run", help="Evaluate (a shard of) the replicates; an unsharded run also merges")
    add_construction_arguments(run)
    run.add_argument(
        "--calibration-threshold",
        type=float,
        default=0.8,
        help="Baseline median |Spearman| of the first component with z required (default: 0.8)",
    )
    merge_parser = sub.add_parser("merge", help="Merge every shard and write the decision (reads settings.json)")
    merge_parser.add_argument("--out-dir", type=Path, required=True, help="Output directory")
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    if args.command == "merge":
        settings = load_settings(args.out_dir)
    else:
        settings = build_settings(args)
        ensure_settings(args.out_dir, settings)
        index, count = args.shard
        evaluated = run_shard(args.out_dir, settings, index, count, n_jobs=args.n_jobs, limit=args.limit)
        print(f"evaluated {evaluated} replicate(s) into {shard_path(args.out_dir, index, count)}", file=sys.stderr)
        if count > 1 or args.limit is not None:
            return 0
    print(json.dumps(merge(args.out_dir, settings), indent=1))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
