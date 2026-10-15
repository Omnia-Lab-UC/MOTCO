#!/usr/bin/env python3
"""Record-level consistency checks for a net-displacement power comparison run.

Two checks the study engine does not compute (change
net-displacement-power-comparison, design D5), both from merged JSONL:

(a) **Shared statistics are identical.** The orientation estimator enters only
    the ``angle`` statistic, so for every (mode, effect, replicate) the
    ``delta`` and ``shape`` p-values of every non-baseline design column must
    equal the baseline column's, and the observed statistics must agree to
    within ``OBSERVED_RELATIVE_TOLERANCE`` (the columns run in separate
    processes, so BLAS may round differently in the last digits). Any larger
    difference is an estimator-plumbing bug, not a finding.
(b) **Baseline reproduces the parents.** The baseline column (``power_primary``)
    shares the parent studies' generator seeds, evaluation parameters and RRPP
    seed, so wherever a parent ``power_primary`` record of the same mode,
    effect and replicate exists, the generator seed and every p-value must be
    equal. Observed statistics are reported as a bit-identical count and the
    largest relative difference (floating point across code revisions, never
    enough to change a p-value). A pair is compared only when both records
    share the magnitude construction (``truth_metadata.transform.magnitude_kind``,
    recorded only by magnitude cells): a parent run under a different
    ``magnitude_kind`` generated different magnitude datasets from the same
    seeds, so those pairs are counted as excluded, not as mismatches.

Both checks report counts, not just booleans. Exit status is 0 when there is no seed
or p-value mismatch and no shared-statistic mismatch, 1 otherwise.

Example:

    python scripts/net_displacement_power_checks.py \
        --merged results/net-displacement-power-<date>/phase6/merged.jsonl \
        --parent results/phase6-small-n-2026-10-09/merged.jsonl \
        --parent results/phase6-small-n-magnitude-2026-10-09/merged.jsonl \
        --out results/net-displacement-power-<date>/phase6/consistency_checks.json
"""

from __future__ import annotations

import argparse
import json
import math
import sys
from collections.abc import Iterable, Mapping, Sequence
from pathlib import Path
from typing import Any

from motco.simulations.grid import SimulationReplicateResult, read_replicate_results
from motco.simulations.study.spectrum import record_design_point

#: Version of these checks; written to the output.
CHECKS_VERSION = 2

#: Largest relative difference in a shared observed statistic still treated as rounding (check a).
OBSERVED_RELATIVE_TOLERANCE = 1e-9

STATISTICS = ("delta", "angle", "shape")
SHARED_STATISTICS = ("delta", "shape")

#: At most this many mismatch examples are listed per check; the counts are complete.
MAX_EXAMPLES = 20

Key = tuple[str, float, int]


def _key(record: SimulationReplicateResult) -> Key:
    meta = record.cell_metadata or {}
    effect = meta.get("effect_size")
    return (
        str(meta.get("trajectory_mode", "none")),
        math.nan if effect is None else float(effect),
        int(record.replicate_index),
    )


def _point_label(point: Mapping[str, Any]) -> str:
    return json.dumps(point, sort_keys=True, default=str)


def _same(left: Any, right: Any) -> bool:
    if left is None or right is None:
        return left is right
    left, right = float(left), float(right)
    if math.isnan(left) or math.isnan(right):
        return math.isnan(left) and math.isnan(right)
    return left == right


def _relative_difference(left: Any, right: Any) -> float | None:
    if left is None or right is None:
        return None
    left, right = float(left), float(right)
    if not (math.isfinite(left) and math.isfinite(right)):
        return None
    scale = max(abs(left), abs(right))
    return 0.0 if scale == 0 else abs(left - right) / scale


def _magnitude_construction(record: SimulationReplicateResult) -> Any:
    return ((record.truth_metadata or {}).get("transform") or {}).get("magnitude_kind")


def _completed(records: Iterable[SimulationReplicateResult]) -> list[SimulationReplicateResult]:
    return [record for record in records if record.status == "completed"]


def check_shared_statistics(
    records: Sequence[SimulationReplicateResult], *, statistics: Sequence[str] = SHARED_STATISTICS
) -> dict[str, Any]:
    """Check (a): every non-baseline column's shared statistics equal the baseline's."""

    columns: dict[str, dict[Key, SimulationReplicateResult]] = {}
    baseline_label: str | None = None
    for record in _completed(records):
        point = record_design_point(record)
        if point is None:
            continue
        label = _point_label(point)
        columns.setdefault(label, {})[_key(record)] = record
        if record.phase == "power_primary":
            baseline_label = label
    if baseline_label is None:
        return {"status": "no_baseline", "statistics": list(statistics), "columns": []}
    baseline = columns[baseline_label]
    results = []
    total_mismatches = 0
    for label, column in sorted(columns.items()):
        if label == baseline_label:
            continue
        compared = 0
        p_mismatch = 0
        stat_mismatch = 0
        bit_identical = 0
        max_relative: float = 0.0
        examples: list[dict[str, Any]] = []
        unpaired = sorted(set(column) ^ set(baseline), key=str)
        for key in sorted(set(column) & set(baseline), key=str):
            mine, ref = column[key], baseline[key]
            compared += 1
            for statistic in statistics:
                p_equal = _same(mine.p_values.get(statistic), ref.p_values.get(statistic))
                left, right = mine.pair_statistics.get(statistic), ref.pair_statistics.get(statistic)
                s_equal = _same(left, right)
                bit_identical += s_equal
                if not s_equal:
                    relative = _relative_difference(left, right)
                    if relative is not None:
                        max_relative = max(max_relative, relative)
                    s_equal = relative is not None and relative <= OBSERVED_RELATIVE_TOLERANCE
                p_mismatch += not p_equal
                stat_mismatch += not s_equal
                if not (p_equal and s_equal) and len(examples) < MAX_EXAMPLES:
                    examples.append(
                        {
                            "trajectory_mode": key[0],
                            "effect_size": key[1],
                            "replicate_index": key[2],
                            "statistic": statistic,
                            "p_value": mine.p_values.get(statistic),
                            "baseline_p_value": ref.p_values.get(statistic),
                            "observed": mine.pair_statistics.get(statistic),
                            "baseline_observed": ref.pair_statistics.get(statistic),
                        }
                    )
        total_mismatches += p_mismatch + stat_mismatch
        results.append(
            {
                "design_point": json.loads(label),
                "pairs_compared": compared,
                "unpaired_records": len(unpaired),
                "comparisons": compared * len(statistics),
                "p_value_mismatches": p_mismatch,
                "observed_statistic_mismatches": stat_mismatch,
                "observed_statistics_bit_identical": bit_identical,
                "observed_statistics_max_relative_difference": max_relative,
                "examples": examples,
            }
        )
    return {
        "status": "identical" if total_mismatches == 0 else "mismatch",
        "statistics": list(statistics),
        "observed_relative_tolerance": OBSERVED_RELATIVE_TOLERANCE,
        "baseline": json.loads(baseline_label),
        "columns": results,
    }


def check_parent_reproduction(
    records: Sequence[SimulationReplicateResult],
    parents: Mapping[str, Sequence[SimulationReplicateResult]],
    *,
    statistics: Sequence[str] = STATISTICS,
) -> dict[str, Any]:
    """Check (b): baseline records equal each parent's primary records where one exists."""

    baseline = {_key(r): r for r in _completed(records) if r.phase == "power_primary"}
    covered: set[Key] = set()
    per_parent = []
    any_mismatch = False
    for name, parent_records in parents.items():
        index = {_key(r): r for r in _completed(parent_records) if r.phase == "power_primary"}
        paired = sorted(set(baseline) & set(index), key=str)
        shared = [k for k in paired if _magnitude_construction(baseline[k]) == _magnitude_construction(index[k])]
        excluded = len(paired) - len(shared)
        covered.update(shared)
        seed_mismatch = 0
        p_mismatch = 0
        bit_identical = 0
        max_relative: float = 0.0
        examples: list[dict[str, Any]] = []
        cells: set[tuple[str, float]] = set()
        for key in shared:
            mine, ref = baseline[key], index[key]
            cells.add((key[0], key[1]))
            seed_equal = mine.generator_seed == ref.generator_seed
            seed_mismatch += not seed_equal
            bad = [s for s in statistics if not _same(mine.p_values.get(s), ref.p_values.get(s))]
            p_mismatch += len(bad)
            identical = True
            for statistic in statistics:
                left, right = mine.pair_statistics.get(statistic), ref.pair_statistics.get(statistic)
                if not _same(left, right):
                    identical = False
                    relative = _relative_difference(left, right)
                    if relative is not None:
                        max_relative = max(max_relative, relative)
            bit_identical += identical
            if (bad or not seed_equal) and len(examples) < MAX_EXAMPLES:
                examples.append(
                    {
                        "trajectory_mode": key[0],
                        "effect_size": key[1],
                        "replicate_index": key[2],
                        "generator_seed": mine.generator_seed,
                        "parent_generator_seed": ref.generator_seed,
                        "p_value_mismatch_statistics": bad,
                    }
                )
        any_mismatch = any_mismatch or seed_mismatch > 0 or p_mismatch > 0
        per_parent.append(
            {
                "parent": name,
                "parent_records": len(index),
                "records_compared": len(shared),
                "records_excluded_magnitude_construction": excluded,
                "cells_compared": [{"trajectory_mode": m, "effect_size": e} for m, e in sorted(cells)],
                "generator_seed_mismatches": seed_mismatch,
                "p_value_mismatches": p_mismatch,
                "observed_statistics_bit_identical": bit_identical,
                "observed_statistics_max_relative_difference": max_relative,
                "examples": examples,
            }
        )
    uncovered_cells = sorted({(key[0], key[1]) for key in set(baseline) - covered})
    return {
        "status": "reproduced" if not any_mismatch else "mismatch",
        "statistics": list(statistics),
        "baseline_records": len(baseline),
        "baseline_records_with_parent": len(covered),
        "baseline_cells_with_uncovered_records": [{"trajectory_mode": m, "effect_size": e} for m, e in uncovered_cells],
        "parents": per_parent,
    }


def run(merged: Path, parents: Sequence[Path], out: Path | None) -> dict[str, Any]:
    # read_replicate_results returns [] for a missing file; a typo must not pass as "nothing to compare".
    missing = [str(path) for path in (merged, *parents) if not Path(path).is_file()]
    if missing:
        raise FileNotFoundError(f"merged JSONL not found: {missing}")
    records = read_replicate_results(merged)
    parent_records = {str(path): read_replicate_results(path) for path in parents}
    payload: dict[str, Any] = {
        "checks_version": CHECKS_VERSION,
        "merged": str(merged),
        "records": len(records),
        "completed_records": len(_completed(records)),
        "shared_statistics": check_shared_statistics(records),
        "parent_reproduction": check_parent_reproduction(records, parent_records) if parents else None,
    }
    if out is not None:
        out.parent.mkdir(parents=True, exist_ok=True)
        out.write_text(json.dumps(payload, indent=2, default=_json_default) + "\n", encoding="utf-8")
    return payload


def _json_default(value: Any) -> Any:
    if isinstance(value, float) and math.isnan(value):
        return None
    raise TypeError(f"Object of type {type(value).__name__} is not JSON serializable")


def failed(payload: Mapping[str, Any]) -> bool:
    shared = payload["shared_statistics"]["status"] == "mismatch"
    parent = payload["parent_reproduction"]
    return shared or (parent is not None and parent["status"] == "mismatch")


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(prog="net_displacement_power_checks", description=__doc__.split("\n\n")[0])
    parser.add_argument("--merged", type=Path, required=True, help="The comparison run's merged.jsonl.")
    parser.add_argument(
        "--parent", type=Path, action="append", default=[], help="A parent run's merged.jsonl (repeatable)."
    )
    parser.add_argument("--out", type=Path, default=None, help="Where to write the checks JSON.")
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    try:
        payload = run(args.merged, args.parent, args.out)
    except FileNotFoundError as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 2
    shared = payload["shared_statistics"]
    print(f"(a) shared statistics {shared['statistics']}: {shared['status']}")
    for column in shared["columns"]:
        print(
            f"    {column['design_point']}: {column['pairs_compared']} pairs, "
            f"{column['p_value_mismatches']} p-value / {column['observed_statistic_mismatches']} observed "
            f"mismatches (tolerance {shared['observed_relative_tolerance']:.0e}), observed bit-identical "
            f"{column['observed_statistics_bit_identical']}, max relative difference "
            f"{column['observed_statistics_max_relative_difference']:.3g}, {column['unpaired_records']} unpaired"
        )
    parent = payload["parent_reproduction"]
    if parent is not None:
        print(
            f"(b) parent reproduction: {parent['status']} — {parent['baseline_records_with_parent']} of "
            f"{parent['baseline_records']} baseline records have a parent record"
        )
        for entry in parent["parents"]:
            print(
                f"    {entry['parent']}: {entry['records_compared']} compared "
                f"({entry['records_excluded_magnitude_construction']} excluded, other magnitude construction), "
                f"{entry['generator_seed_mismatches']} seed / {entry['p_value_mismatches']} p-value mismatches, "
                f"observed bit-identical {entry['observed_statistics_bit_identical']}, "
                f"max relative difference {entry['observed_statistics_max_relative_difference']:.3g}"
            )
    if args.out is not None:
        print(f"  wrote {args.out}")
    return 1 if failed(payload) else 0


if __name__ == "__main__":
    sys.exit(main())
