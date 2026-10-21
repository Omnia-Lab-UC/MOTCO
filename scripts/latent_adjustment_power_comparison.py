#!/usr/bin/env python3
"""Latent-adjustment power comparison at the Phase 6 SEA-AD design (no R).

The second validation of ``openspec/changes/covariate-aware-latent-construction/``
(design D9, spec ``latent-adjustment-validation``): does the covariate-residualized
basis recover power the covariate-blind (0.8.0) basis loses to a covariate-aligned
nuisance axis, without costing power when there is no such axis?

Each **unit** is one (column, cell, replicate):

- **cells**: the zero-effect anchor, the translation control at its effect, and
  magnitude, orientation and shape at the Phase 6 effect levels
  (``group_effect_size``);
- **columns**: ``nuisance`` (the calibrated injection of the latent-adjustment
  Type I check, :mod:`latent_adjustment_type_i_check`) and ``clean``
  (``s_z = s_c = 0``: the same covariates are declared, nothing is injected).

A unit's dataset is the Type I construction with the cell's trajectory mode and
effect: the generator, covariate and RRPP seeds depend on the replicate index
only (matched seeds across cells, columns and arms), so the two columns of one
replicate share the generated blocks and covariates and differ only in the
injection, and the nuisance column's anchor at replicate ``r`` is the Type I
check's replicate ``r``. Each unit is analysed in the three arms of the Type I
check (``baseline``, ``stage``, ``group_stage``) with every covariate declared,
the same CV settings and the same RRPP seed, recording p-values, the CV-selected
rank and the first component's correlation with ``z``.

The **decision** (pure, :func:`latent_adjustment_power_decision`; advisory),
per residualized arm against the baseline arm, following the study engine's
``column_decision`` convention (rate = share of available replicates with
``p ≤ alpha``, Monte Carlo SE ``√(r(1 − r)/R)``, pooled SE of two rates
``√(SE_a² + SE_b²)``):

- **gain** (``nuisance`` column, top effect of each gain target, by default
  orientation/``angle``): ``rate_arm − rate_baseline > gain_se_multiplier ·
  pooled SE``;
- **protected** (``nuisance`` column, top effect of each protected pair, by
  default magnitude/``delta`` and shape/``shape``): ``rate_baseline − rate_arm
  ≤ loss_se_multiplier · pooled SE``;
- **no loss** (``clean`` column, every statistic at every effect of every
  mode): ``rate_baseline − rate_arm ≤ loss_se_multiplier · pooled SE``;
- **controls** (both columns, anchor and translation, every statistic):
  ``rate_arm ≤ alpha + k·√(alpha(1 − alpha)/R)``.

Verdicts are ``candidate_qualifies``, ``candidate_fails`` or ``unavailable``
(a quantity the rule needs has no replicate). The baseline arm is measured,
not judged. The rule and its multipliers are part of the hashed settings, so
they are fixed before any run.

Outputs (``--out-dir``): ``settings.json``; ``records/shard_<i>_of_<N>.jsonl``
(one line per unit, all arms; resume skips recorded units); after ``merge``:
``merged.jsonl``, ``latent_adjustment_power_pvalues.csv``,
``latent_adjustment_power_rates.csv`` (per column × mode × effect × arm ×
statistic: R, rejections, rate, SE, selected-rank distribution, rank-1 fraction)
and ``latent_adjustment_power_decision.json``.

Examples:

    # smoke profile
    python scripts/latent_adjustment_power_comparison.py run --reps 4 --perms 99 --out-dir /tmp/latent-power-smoke

    # one SLURM array task, then the merge
    python scripts/latent_adjustment_power_comparison.py run --shard 3/100 --out-dir results/latent-power
    python scripts/latent_adjustment_power_comparison.py merge --out-dir results/latent-power
"""

from __future__ import annotations

import argparse
import json
import math
import sys
from pathlib import Path
from time import perf_counter
from typing import Any, Iterable, Mapping, Sequence

sys.path.insert(0, str(Path(__file__).resolve().parent))
import covariate_type_i_check as d10  # noqa: E402  (shared settings, CSV and seed helpers)
import latent_adjustment_type_i_check as type_i  # noqa: E402  (construction, arms, arm records)

sys.path.pop(0)

SCRIPT_VERSION = 1
CHECK = "latent_adjustment_power"
STATISTICS: tuple[str, ...] = type_i.STATISTICS
ARMS: tuple[str, ...] = type_i.ARMS
RESIDUALIZED_ARMS: tuple[str, ...] = type_i.RESIDUALIZED_ARMS
COLUMNS: tuple[str, ...] = ("nuisance", "clean")
MODES: tuple[str, ...] = ("magnitude", "orientation", "shape")
ANCHOR = "anchor"
CONTROL_MODES: tuple[str, ...] = (ANCHOR, "translation")

#: Phase 6 effect levels per mode (``phase6_small_n_magnitude.json`` below its
#: saturation, ``phase6_small_n_study.json`` for orientation and shape).
DEFAULT_EFFECTS: dict[str, list[float]] = {
    "magnitude": [0.02, 0.05, 0.1, 0.25],
    "orientation": [0.25, 0.5, 0.75, 1.0],
    "shape": [0.25, 0.5, 0.75, 1.0],
}
DEFAULT_TRANSLATION_EFFECT = 1.0


# ── Settings ─────────────────────────────────────────────────────────────────


def _floats(text: str) -> list[float]:
    return [float(value) for value in text.split(",") if value.strip()]


def build_settings(args: argparse.Namespace) -> dict[str, Any]:
    settings = type_i.build_settings(args)
    settings["script_version"] = SCRIPT_VERSION
    settings["check"] = CHECK
    settings.pop("calibration")
    injection = settings["injection"]
    settings["columns"] = {
        "nuisance": {"s_z": injection["s_z"], "s_c": injection["s_c"]},
        "clean": {"s_z": 0.0, "s_c": 0.0},
    }
    settings["effects"] = {mode: _floats(getattr(args, f"{mode}_effects")) for mode in MODES}
    settings["controls"] = {ANCHOR: 0.0, "translation": float(args.translation_effect)}
    settings["rule"] = {
        "rate": "share of available replicates with p <= alpha",
        "gain_targets": [["orientation", "angle"]],
        "protected": [["magnitude", "delta"], ["shape", "shape"]],
        "no_loss": "clean column, every statistic at every effect of every mode",
        "controls": "both columns, anchor and translation, every statistic, rate <= alpha + k*sqrt(alpha(1-alpha)/R)",
        "gain_se_multiplier": float(args.gain_se_multiplier),
        "loss_se_multiplier": float(args.loss_se_multiplier),
    }
    return settings


settings_hash = d10.settings_hash
ensure_settings = d10.ensure_settings
load_settings = d10.load_settings
parse_shard = d10.parse_shard
shard_path = d10.shard_path


def cells(settings: Mapping[str, Any]) -> list[tuple[str, float]]:
    """(mode, effect) per cell: anchor, translation, then every mode's effects."""

    out = [(ANCHOR, 0.0), ("translation", float(settings["controls"]["translation"]))]
    for mode in MODES:
        out.extend((mode, float(effect)) for effect in settings["effects"][mode])
    return out


def units(settings: Mapping[str, Any]) -> list[tuple[str, str, float, int]]:
    """Every (column, mode, effect, replicate) in a fixed order; shards take every N-th."""

    return [
        (column, mode, effect, replicate)
        for replicate in range(int(settings["replicates"]))
        for column in COLUMNS
        for mode, effect in cells(settings)
    ]


def unit_key(column: str, mode: str, effect: float, replicate: int) -> str:
    return json.dumps([column, mode, float(effect), int(replicate)])


def cell_settings(settings: Mapping[str, Any], column: str, mode: str, effect: float) -> dict[str, Any]:
    """The Type I settings with this cell's generator mode/effect and this column's injection."""

    generator = dict(settings["generator"])
    generator["trajectory_mode"] = "magnitude" if mode == ANCHOR else mode
    generator["group_effect_size"] = float(effect)
    injection = {**settings["injection"], **settings["columns"][column]}
    return {**settings, "generator": generator, "injection": injection}


# ── Unit evaluation ──────────────────────────────────────────────────────────


def build_unit(
    column: str, mode: str, effect: float, replicate: int, settings: Mapping[str, Any], reference: Any = None
) -> dict[str, Any]:
    return type_i.build_replicate(replicate, cell_settings(settings, column, mode, effect), reference=reference)


def evaluate_unit(
    column: str,
    mode: str,
    effect: float,
    replicate: int,
    settings: Mapping[str, Any],
    *,
    n_jobs: int = 1,
    reference: Any = None,
) -> dict[str, Any]:
    start = perf_counter()
    built = build_unit(column, mode, effect, replicate, settings, reference=reference)
    arms = type_i.analyze_arms(built, settings, n_jobs=n_jobs)
    return d10._jsonable(
        {
            "key": unit_key(column, mode, effect, replicate),
            "column": column,
            "mode": mode,
            "effect": float(effect),
            "replicate": replicate,
            "settings_hash": settings_hash(settings),
            "seeds": built["seeds"],
            "arms": arms,
            "runtime_seconds": perf_counter() - start,
        }
    )


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
            entry = record["arms"][arm]
            for statistic in STATISTICS:
                rows.append(
                    {
                        "column": record["column"],
                        "mode": record["mode"],
                        "effect": float(record["effect"]),
                        "replicate": int(record["replicate"]),
                        "arm": arm,
                        "statistic": statistic,
                        "p_value": entry["p_values"].get(statistic),
                        "selected_lv": entry.get("selected_lv"),
                    }
                )
    return rows


# ── Rates and decision (pure) ────────────────────────────────────────────────


RATE_COLUMNS = (
    "column",
    "mode",
    "effect",
    "arm",
    "statistic",
    "replicates",
    "rejections",
    "rate",
    "se",
    "rank_distribution",
    "rank1_fraction",
)


def rate_table(rows: Iterable[Mapping[str, Any]], *, alpha: float) -> list[dict[str, Any]]:
    """Per column × mode × effect × arm × statistic: R, rejections, rate, SE and the selected ranks.

    The rank distribution counts each replicate once (from its p-value rows),
    including replicates whose statistic is unavailable.
    """

    tallies: dict[tuple, list[int]] = {}
    ranks: dict[tuple, dict[int, int]] = {}
    seen: set[tuple] = set()
    for row in rows:
        cell = (str(row["column"]), str(row["mode"]), float(row["effect"]), str(row["arm"]))
        tally = tallies.setdefault((*cell, str(row["statistic"])), [0, 0])
        p_value = row.get("p_value")
        if p_value is not None and math.isfinite(float(p_value)):
            tally[0] += 1
            tally[1] += int(float(p_value) <= alpha)
        replicate = (*cell, int(row.get("replicate", -1)))
        if replicate not in seen:
            seen.add(replicate)
            rank = row.get("selected_lv")
            if rank is not None:
                counts = ranks.setdefault(cell, {})
                counts[int(rank)] = counts.get(int(rank), 0) + 1
    out = []
    def order(key: tuple) -> tuple:
        column, mode, effect, arm, statistic = key
        return (
            COLUMNS.index(column) if column in COLUMNS else len(COLUMNS),
            mode,
            effect,
            ARMS.index(arm) if arm in ARMS else len(ARMS),
            statistic,
        )

    for key in sorted(tallies, key=order):
        replicates, rejections = tallies[key]
        rate = rejections / replicates if replicates else None
        counts = ranks.get(key[:4], {})
        total = sum(counts.values())
        out.append(
            {
                "column": key[0],
                "mode": key[1],
                "effect": key[2],
                "arm": key[3],
                "statistic": key[4],
                "replicates": replicates,
                "rejections": rejections,
                "rate": rate,
                "se": math.sqrt(rate * (1.0 - rate) / replicates) if rate is not None else None,
                "rank_distribution": ";".join(f"{rank}:{counts[rank]}" for rank in sorted(counts)),
                "rank1_fraction": counts.get(1, 0) / total if total else None,
            }
        )
    return out


def _lookup(table: Sequence[Mapping[str, Any]]) -> dict[tuple, Mapping[str, Any]]:
    return {(r["column"], r["mode"], float(r["effect"]), r["arm"], r["statistic"]): r for r in table}


def _pooled_se(a: Mapping[str, Any], b: Mapping[str, Any]) -> float:
    return math.sqrt(float(a["se"] or 0.0) ** 2 + float(b["se"] or 0.0) ** 2)


def _available(*rows: Mapping[str, Any] | None) -> bool:
    return all(row is not None and row.get("rate") is not None for row in rows)


def latent_adjustment_power_decision(
    rows: Iterable[Mapping[str, Any]],
    *,
    alpha: float,
    k: float,
    effects: Mapping[str, Sequence[float]],
    controls: Mapping[str, float],
    gain_targets: Sequence[Sequence[str]] = (("orientation", "angle"),),
    protected: Sequence[Sequence[str]] = (("magnitude", "delta"), ("shape", "shape")),
    gain_se_multiplier: float = 2.0,
    loss_se_multiplier: float = 2.0,
) -> dict[str, Any]:
    """Per residualized arm against the baseline: gain, protected, no-loss and control reads with a verdict."""

    table = rate_table(rows, alpha=alpha)
    lookup = _lookup(table)
    top = {mode: max(float(e) for e in values) for mode, values in effects.items()}
    arms: dict[str, Any] = {}
    for arm in RESIDUALIZED_ARMS:
        reads: dict[str, list[dict[str, Any]]] = {"gain": [], "protected": [], "no_loss": [], "controls": []}

        def compare(column: str, mode: str, effect: float, statistic: str, kind: str) -> None:
            candidate = lookup.get((column, mode, effect, arm, statistic))
            baseline = lookup.get((column, mode, effect, "baseline", statistic))
            entry: dict[str, Any] = {"column": column, "mode": mode, "effect": effect, "statistic": statistic}
            if not _available(candidate, baseline):
                entry["passed"] = None
            else:
                assert candidate is not None and baseline is not None
                pooled = _pooled_se(candidate, baseline)
                difference = float(candidate["rate"]) - float(baseline["rate"])
                multiplier = gain_se_multiplier if kind == "gain" else loss_se_multiplier
                threshold = multiplier * pooled
                entry.update(
                    {
                        "rate": candidate["rate"],
                        "baseline_rate": baseline["rate"],
                        "difference": difference,
                        "pooled_se": pooled,
                        "threshold": threshold,
                        "passed": difference > threshold if kind == "gain" else -difference <= threshold,
                    }
                )
            reads[kind].append(entry)

        for mode, statistic in gain_targets:
            compare("nuisance", mode, top[mode], statistic, "gain")
        for mode, statistic in protected:
            compare("nuisance", mode, top[mode], statistic, "protected")
        for mode in effects:
            for effect in effects[mode]:
                for statistic in STATISTICS:
                    compare("clean", mode, float(effect), statistic, "no_loss")
        for column in COLUMNS:
            for mode, effect in controls.items():
                for statistic in STATISTICS:
                    row = lookup.get((column, mode, float(effect), arm, statistic))
                    entry = {"column": column, "mode": mode, "effect": float(effect), "statistic": statistic}
                    if row is None or row.get("rate") is None:
                        entry["passed"] = None
                    else:
                        bound = d10.anchor_bound(alpha, int(row["replicates"]), k)
                        entry.update({"rate": row["rate"], "bound": bound, "passed": float(row["rate"]) <= bound})
                    reads["controls"].append(entry)
        outcomes = [entry["passed"] for entries in reads.values() for entry in entries]
        if any(outcome is None for outcome in outcomes):
            verdict = "unavailable"
        elif all(outcomes):
            verdict = "candidate_qualifies"
        else:
            verdict = "candidate_fails"
        failing = {kind: [e for e in entries if e["passed"] is False] for kind, entries in reads.items()}
        arms[arm] = {"verdict": verdict, "failing": {k_: v for k_, v in failing.items() if v}, "reads": reads}
    return {
        "rule": {
            "rate": "share of available replicates with p <= alpha",
            "gain": "nuisance column, top effect: rate_arm - rate_baseline > gain_se_multiplier * pooled SE",
            "protected": "nuisance column, top effect: rate_baseline - rate_arm <= loss_se_multiplier * pooled SE",
            "no_loss": "clean column, every effect: rate_baseline - rate_arm <= loss_se_multiplier * pooled SE",
            "controls": "both columns: rate_arm <= alpha + k * sqrt(alpha * (1 - alpha) / R)",
            "gain_targets": [list(pair) for pair in gain_targets],
            "protected_pairs": [list(pair) for pair in protected],
            "gain_se_multiplier": gain_se_multiplier,
            "loss_se_multiplier": loss_se_multiplier,
            "alpha": alpha,
            "k": k,
        },
        "top_effects": top,
        "arms": arms,
        "verdicts": {arm: arms[arm]["verdict"] for arm in RESIDUALIZED_ARMS},
        "baseline_role": "measured, not judged (covariates in the trajectory model only, covariate-blind basis)",
    }


def write_decision(out_dir: Path, rows: Sequence[Mapping[str, Any]], settings: Mapping[str, Any]) -> dict[str, Any]:
    alpha, k = float(settings["alpha"]), float(settings["k"])
    rule = settings["rule"]
    decision = latent_adjustment_power_decision(
        rows,
        alpha=alpha,
        k=k,
        effects=settings["effects"],
        controls=settings["controls"],
        gain_targets=rule["gain_targets"],
        protected=rule["protected"],
        gain_se_multiplier=float(rule["gain_se_multiplier"]),
        loss_se_multiplier=float(rule["loss_se_multiplier"]),
    )
    decision["settings_hash"] = settings_hash(settings)
    decision["construction"] = {
        "generator": settings["generator"],
        "layers": settings["layers"],
        "covariates": settings["covariates"],
        "declared_covariates": type_i.declared_covariates(settings),
        "columns": settings["columns"],
        "effects": settings["effects"],
        "controls": settings["controls"],
        "permutations": settings["permutations"],
        "replicates": settings["replicates"],
        "base_seed": settings["base_seed"],
        "cv": settings["cv"],
    }
    (out_dir / "latent_adjustment_power_decision.json").write_text(
        json.dumps(d10._jsonable(decision), indent=2, allow_nan=False) + "\n", encoding="utf-8"
    )
    d10._write_csv(out_dir / "latent_adjustment_power_rates.csv", rate_table(rows, alpha=alpha), RATE_COLUMNS)
    return decision


def merge(out_dir: Path, settings: Mapping[str, Any]) -> dict[str, Any]:
    """Merge every shard, then write the p-value table, the rates and the decision."""

    records = read_records(sorted((out_dir / "records").glob("shard_*.jsonl")))
    expected = len(units(settings))
    with (out_dir / "merged.jsonl").open("w", encoding="utf-8") as handle:
        for key in sorted(records):
            handle.write(json.dumps(records[key], sort_keys=True, allow_nan=False) + "\n")
    rows = p_value_rows(records)
    d10._write_csv(
        out_dir / "latent_adjustment_power_pvalues.csv",
        rows,
        ("column", "mode", "effect", "replicate", "arm", "statistic", "p_value", "selected_lv"),
    )
    decision = write_decision(out_dir, rows, settings)
    failures = {arm: sum(1 for r in records.values() if r["arms"][arm].get("error") is not None) for arm in ARMS}
    return {"units": len(records), "expected": expected, "verdicts": decision["verdicts"], "failures": failures}


# ── CLI ──────────────────────────────────────────────────────────────────────


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Latent-adjustment power comparison at the Phase 6 design (no R).")
    sub = parser.add_subparsers(dest="command", required=True)
    run = sub.add_parser("run", help="Evaluate (a shard of) the units; an unsharded run also merges")
    type_i.add_construction_arguments(run)
    run.add_argument(
        "--magnitude-effects", default=",".join(map(str, DEFAULT_EFFECTS["magnitude"])), help="Magnitude effects"
    )
    run.add_argument(
        "--orientation-effects", default=",".join(map(str, DEFAULT_EFFECTS["orientation"])), help="Orientation effects"
    )
    run.add_argument("--shape-effects", default=",".join(map(str, DEFAULT_EFFECTS["shape"])), help="Shape effects")
    run.add_argument(
        "--translation-effect", type=float, default=DEFAULT_TRANSLATION_EFFECT, help="Translation control effect"
    )
    run.add_argument("--gain-se-multiplier", type=float, default=2.0, help="Gain threshold in pooled SEs (default: 2)")
    run.add_argument("--loss-se-multiplier", type=float, default=2.0, help="Loss tolerance in pooled SEs (default: 2)")
    merge_parser = sub.add_parser("merge", help="Merge every shard and write the decision (reads settings.json)")
    merge_parser.add_argument("--out-dir", type=Path, required=True, help="Output directory")
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    if args.command == "merge":
        settings = load_settings(args.out_dir)
    else:
        args.calibration_threshold = 0.8  # read by the shared Type I settings builder, then dropped
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
