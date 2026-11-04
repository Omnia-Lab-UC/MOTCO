#!/usr/bin/env python3
"""Compare a study's shared zero-effect anchor against a reference run's anchor.

Two configs that share ``base_seed``, the primary matched-seed family, and the
generator parameters (other than parameters no mode reads at zero effect, such
as ``magnitude_kind``) generate the same anchor dataset at every replicate index
and evaluate it with the same RRPP draws. Their anchor records therefore
reproduce each other — the observed ``delta``/``angle``/``shape``, the RRPP
p-values, the selected PLS rank, and the generator seed — even though their
``parameter_signature`` differs (the dataclass and the resolved-mode set do).

This script joins the two anchors on ``replicate_index`` and writes one row per
replicate with both sides' values and their absolute differences, plus a summary
line on stdout. Statistics are compared at a relative tolerance (BLAS last-bit
differences between hosts are expected); seeds must match exactly, and p-values
must match exactly when both runs used the same permutation count (the RRPP
draws are then identical). When the counts differ (a 199-permutation pilot
against a 999-permutation reference) the p-values are reported but excluded
from the match verdict, since they cannot agree beyond Monte Carlo resolution.

``--by-column`` compares design-grid studies: every zero-effect anchor of the
``power_primary``/``power_design`` phases is keyed by its design column (the
record's ``cell_metadata.design_point``; ``power_primary`` is the baseline
column; the ``type_i_baseline`` anchor is a different draw and is excluded),
the runs are joined on (design column, replicate index), the CSV gains the
column coordinates, and the summary prints one line per column, listing a
column the other run lacks with zero pairs. The verdict covers the paired
columns only. Columns are keyed by their design point and never by
``generator_seed``, which every column shares at a replicate index within a
seed family. Without the flag the single-anchor comparison is unchanged.

``--by-cell`` pairs every cell, not only the anchors, across runs whose design
grids name different axes (the stage-span pilot's surgery siblings against the
design-guidance pilot). Each ``power_primary``/``power_design`` record is
filtered by ``--candidate-where`` / ``--reference-where`` (``AXIS=JSON``, the
record's design point must hold that value), its design point is projected onto
the ``--project`` axes, and the runs are joined on (projected design column,
trajectory mode, effect size, replicate index). After filtering and projection a
key may name one cell per run. The summary prints one line per (column, mode,
effect), lists a cell the other run lacks with zero pairs, and the verdict
covers the paired cells only.

Example:

    python scripts/anchor_reproduction.py \\
        --candidate results/phase5-magnitude-2026-10-08/merged.jsonl \\
        --reference results/phase5-2026-09-10/merged.jsonl \\
        --out results/phase5-magnitude-2026-10-08/report/anchor_reproduction.csv

    python scripts/anchor_reproduction.py --by-column \\
        --candidate results/design-guidance-slice-<date>/core_magnitude/merged.jsonl \\
        --reference results/design-guidance-slice-<date>/core_geometry/merged.jsonl \\
        --out results/design-guidance-slice-<date>/reproduction/siblings.csv

    python scripts/anchor_reproduction.py --by-cell --project generator.n_samples \\
        --candidate-where 'generator.group_stage_sizes=null' \\
        --reference-where 'generator.n_stages=3' \\
        --reference-where 'evaluation.integration_params.layers=["methylation", "expression"]' \\
        --candidate results/stage-span-pilot-<date>/surgery_geometry/merged.jsonl \\
        --reference results/design-guidance-pilot-2026-10-28/core_geometry/merged.jsonl \\
        --out results/stage-span-pilot-<date>/reproduction/surgery_geometry_vs_design_guidance.csv
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd

STATISTICS = ("delta", "angle", "shape")
#: Phases whose zero-effect anchor sits in a design column (``power_primary`` is the baseline column).
DESIGN_PHASES = ("power_primary", "power_design")
#: The CSV column keying a design column: its design point as sorted JSON.
COLUMN_KEY = "design_column"
#: The join keys of the cell-level pairing beyond the design column.
CELL_KEYS = ("trajectory_mode", "effect_size", "replicate_index")


def _record_row(record: dict) -> dict:
    row = {
        "cell_id": record["cell_id"],
        "replicate_index": int(record["replicate_index"]),
        "status": record.get("status"),
        "generator_seed": record.get("generator_seed"),
        "selected_lv": (record.get("integration_metadata") or {}).get("selected_lv"),
        "permutations": (record.get("runtime_metadata") or {}).get("permutations"),
    }
    for statistic in STATISTICS:
        row[f"observed_{statistic}"] = (record.get("pair_statistics") or {}).get(statistic)
        row[f"p_{statistic}"] = (record.get("p_values") or {}).get(statistic)
    return row


def _anchor_rows(path: Path, cell_id: str | None) -> pd.DataFrame:
    rows = []
    with path.open(encoding="utf-8") as handle:
        for line in handle:
            if not line.strip():
                continue
            record = json.loads(line)
            if cell_id is not None:
                if record.get("cell_id") != cell_id:
                    continue
            elif not (record.get("cell_metadata") or {}).get("zero_effect_anchor"):
                continue
            rows.append(_record_row(record))
    if not rows:
        raise SystemExit(f"No anchor records found in {path}")
    frame = pd.DataFrame(rows)
    if frame["cell_id"].nunique() != 1:
        raise SystemExit(f"{path}: anchor records span several cells: {sorted(frame['cell_id'].unique())}")
    return frame.sort_values("replicate_index").reset_index(drop=True)


def column_key(point: dict) -> str:
    return json.dumps(point, sort_keys=True)


def _column_anchor_rows(path: Path) -> pd.DataFrame:
    """Every design column's zero-effect anchor, keyed by the column's design point.

    ``power_primary`` records carry the baseline design point (``{}`` in a study
    without a design grid); ``power_design`` records carry their own;
    ``type_i_baseline`` anchors are a different draw and are left out.
    """

    rows = []
    with path.open(encoding="utf-8") as handle:
        for line in handle:
            if not line.strip():
                continue
            record = json.loads(line)
            meta = record.get("cell_metadata") or {}
            if not meta.get("zero_effect_anchor") or record.get("phase") not in DESIGN_PHASES:
                continue
            row = _record_row(record)
            row[COLUMN_KEY] = column_key(meta.get("design_point") or {})
            row["is_baseline"] = record.get("phase") == "power_primary"
            rows.append(row)
    if not rows:
        raise SystemExit(f"No design-column anchor records found in {path}")
    frame = pd.DataFrame(rows)
    cells = frame.groupby(COLUMN_KEY)["cell_id"].nunique()
    if (cells > 1).any():
        crowded = sorted(cells[cells > 1].index)
        raise SystemExit(f"{path}: several anchor cells at one design column: {crowded}")
    return frame.sort_values([COLUMN_KEY, "replicate_index"]).reset_index(drop=True)


def parse_where(clauses: list[str] | None) -> dict[str, object]:
    """``AXIS=JSON`` clauses as ``{axis: value}``; a value that is not JSON is taken as a string."""

    where: dict[str, object] = {}
    for clause in clauses or []:
        axis, sep, raw = clause.partition("=")
        if not sep or not axis:
            raise SystemExit(f"A where clause must read AXIS=VALUE, got {clause!r}.")
        try:
            where[axis] = json.loads(raw)
        except json.JSONDecodeError:
            where[axis] = raw
        # The comparison is on the JSON form, so a tuple and a list name the same value.
    return where


def _cell_rows(path: Path, *, project: list[str], where: dict[str, object]) -> pd.DataFrame:
    """Every design-phase record kept by ``where``, keyed by its projected column, mode and effect."""

    rows = []
    with path.open(encoding="utf-8") as handle:
        for line in handle:
            if not line.strip():
                continue
            record = json.loads(line)
            if record.get("phase") not in DESIGN_PHASES:
                continue
            meta = record.get("cell_metadata") or {}
            point = meta.get("design_point") or {}
            if any(axis not in point or json.dumps(point[axis]) != json.dumps(value) for axis, value in where.items()):
                continue
            missing = [axis for axis in project if axis not in point]
            if missing:
                raise SystemExit(f"{path}: cell {record.get('cell_id')} has no design axis {missing} to project on.")
            row = _record_row(record)
            row[COLUMN_KEY] = column_key({axis: point[axis] for axis in project})
            row["trajectory_mode"] = str(meta.get("trajectory_mode"))
            row["effect_size"] = float(meta.get("effect_size", 0.0))
            rows.append(row)
    if not rows:
        raise SystemExit(f"No design-phase records left in {path} after filtering.")
    frame = pd.DataFrame(rows)
    cells = frame.groupby([COLUMN_KEY, "trajectory_mode", "effect_size"])["cell_id"].nunique()
    if (cells > 1).any():
        crowded = [list(key) for key in cells[cells > 1].index]
        raise SystemExit(
            f"{path}: several cells share a (projected column, mode, effect): {crowded}; "
            "filter with --candidate-where/--reference-where or project on more axes."
        )
    return frame.sort_values([COLUMN_KEY, *CELL_KEYS]).reset_index(drop=True)


def compare(
    candidate: pd.DataFrame,
    reference: pd.DataFrame,
    rtol: float,
    *,
    by_column: bool = False,
    by_cell: bool = False,
) -> pd.DataFrame:
    if by_cell:
        keys: str | list[str] = [COLUMN_KEY, *CELL_KEYS]
    else:
        keys = [COLUMN_KEY, "replicate_index"] if by_column else "replicate_index"
    merged = candidate.merge(reference, on=keys, how="inner", suffixes=("_candidate", "_reference"))
    merged["seed_match"] = merged["generator_seed_candidate"] == merged["generator_seed_reference"]
    merged["lv_match"] = merged["selected_lv_candidate"] == merged["selected_lv_reference"]
    for statistic in STATISTICS:
        c = merged[f"observed_{statistic}_candidate"].to_numpy(dtype=float)
        r = merged[f"observed_{statistic}_reference"].to_numpy(dtype=float)
        merged[f"abs_diff_{statistic}"] = np.abs(c - r)
        merged[f"{statistic}_match"] = np.isclose(c, r, rtol=rtol, atol=0.0)
        merged[f"p_{statistic}_match"] = (
            merged[f"p_{statistic}_candidate"] == merged[f"p_{statistic}_reference"]
        )
    merged["p_values_comparable"] = merged["permutations_candidate"] == merged["permutations_reference"]
    columns = ["seed_match", "lv_match"] + [f"{s}_match" for s in STATISTICS]
    merged["all_match"] = merged[columns].all(axis=1) & (
        ~merged["p_values_comparable"] | merged[[f"p_{s}_match" for s in STATISTICS]].all(axis=1)
    )
    if by_column or by_cell:
        merged = _with_column_coordinates(merged, keys=list(CELL_KEYS) if by_cell else ["replicate_index"])
    return merged


def _with_column_coordinates(merged: pd.DataFrame, *, keys: list[str] | None = None) -> pd.DataFrame:
    """One column per design axis, read back from the column key, ahead of the remaining join keys."""

    keys = keys or ["replicate_index"]
    points = [json.loads(key) for key in merged[COLUMN_KEY]]
    axes = sorted({axis for point in points for axis in point})
    coordinates = pd.DataFrame(
        {axis: [json.dumps(point.get(axis)) if axis in point else None for point in points] for axis in axes},
        index=merged.index,
    )
    rest = [column for column in merged.columns if column not in (COLUMN_KEY, *keys)]
    return pd.concat([merged[[COLUMN_KEY]], coordinates, merged[[*keys, *rest]]], axis=1)


def _column_order(candidate: pd.DataFrame, reference: pd.DataFrame) -> list[tuple[str, bool]]:
    """Every column either run has, the baseline first, with its baseline flag."""

    baseline = {
        key: bool(flag)
        for frame in (reference, candidate)
        for key, flag in zip(frame[COLUMN_KEY], frame["is_baseline"], strict=True)
    }
    return sorted(baseline.items(), key=lambda item: (not item[1], item[0]))


def _print_column_summary(
    candidate: pd.DataFrame, reference: pd.DataFrame, merged: pd.DataFrame, rtol: float
) -> None:
    n_pairs = len(merged)
    paired_columns = merged[COLUMN_KEY].nunique()
    print(
        f"candidate {len(candidate)} anchor records over {candidate[COLUMN_KEY].nunique()} design column(s) vs "
        f"reference {len(reference)} anchor records over {reference[COLUMN_KEY].nunique()}; "
        f"{n_pairs} (design column, replicate index) pairs over {paired_columns} paired column(s)"
    )
    for key, is_baseline in _column_order(candidate, reference):
        pairs = merged[merged[COLUMN_KEY] == key]
        n_candidate = int((candidate[COLUMN_KEY] == key).sum())
        n_reference = int((reference[COLUMN_KEY] == key).sum())
        label = f"  column {key}{' [baseline]' if is_baseline else ''}"
        if pairs.empty:
            print(f"{label}: 0 pairs (candidate {n_candidate} records, reference {n_reference} records)")
            continue
        n = len(pairs)
        parts = [
            f"{n} pairs (candidate {n_candidate}, reference {n_reference})",
            f"seeds identical {int(pairs['seed_match'].sum())}/{n}",
            f"rank identical {int(pairs['lv_match'].sum())}/{n}",
        ]
        for statistic in STATISTICS:
            parts.append(
                f"{statistic} within rtol {rtol:g} {int(pairs[f'{statistic}_match'].sum())}/{n} "
                f"(max abs diff {pairs[f'abs_diff_{statistic}'].max():.3g})"
            )
        if bool(pairs["p_values_comparable"].all()):
            identical = int(pairs[[f"p_{s}_match" for s in STATISTICS]].all(axis=1).sum())
            parts.append(f"p-values identical {identical}/{n}")
        else:
            candidate_perms = sorted(int(v) for v in pairs["permutations_candidate"].unique())
            reference_perms = sorted(int(v) for v in pairs["permutations_reference"].unique())
            parts.append(f"p-values not comparable ({candidate_perms} vs {reference_perms})")
        parts.append(f"all comparable fields match {int(pairs['all_match'].sum())}/{n}")
        print(f"{label}: " + "; ".join(parts))
    if n_pairs and not bool(merged["p_values_comparable"].all()):
        print("  p-values not comparable on some pairs (permutation counts differ); excluded from the verdict")
    print(f"  all comparable fields match: {int(merged['all_match'].sum())}/{n_pairs} (paired columns only)")


def _print_cell_summary(
    candidate: pd.DataFrame, reference: pd.DataFrame, merged: pd.DataFrame, rtol: float
) -> None:
    cell = [COLUMN_KEY, "trajectory_mode", "effect_size"]
    n_pairs = len(merged)
    print(
        f"candidate {len(candidate)} records over {len(candidate.groupby(cell))} cell(s) vs "
        f"reference {len(reference)} records over {len(reference.groupby(cell))}; "
        f"{n_pairs} (column, mode, effect, replicate index) pairs over "
        f"{len(merged.groupby(cell)) if n_pairs else 0} paired cell(s)"
    )
    keys = sorted(
        {tuple(row) for frame in (candidate, reference) for row in frame[cell].itertuples(index=False)},
        key=lambda key: (key[0], key[1], key[2]),
    )
    for column, mode, effect in keys:
        def select(frame: pd.DataFrame) -> pd.DataFrame:
            return frame[
                (frame[COLUMN_KEY] == column) & (frame["trajectory_mode"] == mode) & (frame["effect_size"] == effect)
            ]

        pairs = select(merged)
        n_candidate, n_reference = len(select(candidate)), len(select(reference))
        label = f"  cell {column} {mode}@{effect:g}"
        if pairs.empty:
            print(f"{label}: 0 pairs (candidate {n_candidate} records, reference {n_reference} records)")
            continue
        n = len(pairs)
        parts = [
            f"{n} pairs (candidate {n_candidate}, reference {n_reference})",
            f"seeds identical {int(pairs['seed_match'].sum())}/{n}",
            f"rank identical {int(pairs['lv_match'].sum())}/{n}",
        ]
        for statistic in STATISTICS:
            parts.append(f"{statistic} within rtol {rtol:g} {int(pairs[f'{statistic}_match'].sum())}/{n}")
        if bool(pairs["p_values_comparable"].all()):
            identical = int(pairs[[f"p_{s}_match" for s in STATISTICS]].all(axis=1).sum())
            parts.append(f"p-values identical {identical}/{n}")
        else:
            parts.append("p-values not comparable")
        parts.append(f"all comparable fields match {int(pairs['all_match'].sum())}/{n}")
        print(f"{label}: " + "; ".join(parts))
    if n_pairs and not bool(merged["p_values_comparable"].all()):
        print("  p-values not comparable on some pairs (permutation counts differ); excluded from the verdict")
    print(f"  all comparable fields match: {int(merged['all_match'].sum())}/{n_pairs} (paired cells only)")


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--candidate", type=Path, required=True, help="Merged JSONL of the run under check.")
    parser.add_argument(
        "--reference", type=Path, required=True, help="Merged JSONL (or anchor extract) of the reference run."
    )
    parser.add_argument(
        "--candidate-cell", default=None, help="Anchor cell id in the candidate (default: the zero_effect_anchor cell)."
    )
    parser.add_argument(
        "--reference-cell", default=None, help="Anchor cell id in the reference (default: the zero_effect_anchor cell)."
    )
    parser.add_argument("--rtol", type=float, default=1e-6, help="Relative tolerance for the observed statistics.")
    parser.add_argument("--out", type=Path, required=True, help="CSV path for the per-replicate comparison.")
    parser.add_argument(
        "--by-column",
        action="store_true",
        help="Pair every design column's anchor on (design column, replicate index) instead of one anchor cell.",
    )
    parser.add_argument(
        "--by-cell",
        action="store_true",
        help="Pair every design-phase cell on (projected design column, mode, effect, replicate index).",
    )
    parser.add_argument(
        "--project", action="append", default=None, help="--by-cell: a design axis the column key keeps (repeatable)."
    )
    parser.add_argument(
        "--candidate-where", action="append", default=None, help="--by-cell: AXIS=JSON filter on candidate records."
    )
    parser.add_argument(
        "--reference-where", action="append", default=None, help="--by-cell: AXIS=JSON filter on reference records."
    )
    args = parser.parse_args(argv)

    cell_only = (args.project, args.candidate_where, args.reference_where)
    if not args.by_cell and any(value is not None for value in cell_only):
        parser.error("--project/--candidate-where/--reference-where apply to --by-cell only.")
    if args.by_cell:
        if args.by_column:
            parser.error("--by-cell and --by-column are exclusive.")
        if args.candidate_cell is not None or args.reference_cell is not None:
            parser.error("--by-cell reads every design-phase cell; --candidate-cell/--reference-cell do not apply.")
        if not args.project:
            parser.error("--by-cell needs at least one --project axis.")
        candidate = _cell_rows(args.candidate, project=args.project, where=parse_where(args.candidate_where))
        reference = _cell_rows(args.reference, project=args.project, where=parse_where(args.reference_where))
        merged = compare(candidate, reference, args.rtol, by_cell=True)
        args.out.parent.mkdir(parents=True, exist_ok=True)
        merged.to_csv(args.out, index=False)
        _print_cell_summary(candidate, reference, merged, args.rtol)
        print(f"Wrote {args.out}")
        # Nothing paired is a failed reproduction, not a vacuous match.
        return 0 if len(merged) and bool(merged["all_match"].all()) else 1

    if args.by_column:
        if args.candidate_cell is not None or args.reference_cell is not None:
            parser.error(
                "--by-column reads every design column's anchor; --candidate-cell/--reference-cell do not apply."
            )
        candidate = _column_anchor_rows(args.candidate)
        reference = _column_anchor_rows(args.reference)
        merged = compare(candidate, reference, args.rtol, by_column=True)
        args.out.parent.mkdir(parents=True, exist_ok=True)
        merged.to_csv(args.out, index=False)
        _print_column_summary(candidate, reference, merged, args.rtol)
        print(f"Wrote {args.out}")
        return 0 if bool(merged["all_match"].all()) else 1

    candidate = _anchor_rows(args.candidate, args.candidate_cell)
    reference = _anchor_rows(args.reference, args.reference_cell)
    merged = compare(candidate, reference, args.rtol)
    args.out.parent.mkdir(parents=True, exist_ok=True)
    merged.to_csv(args.out, index=False)

    n = len(merged)
    print(
        f"candidate {candidate['cell_id'].iat[0]} ({len(candidate)} records) vs "
        f"reference {reference['cell_id'].iat[0]} ({len(reference)} records); "
        f"{n} replicate indices in common"
    )
    print(
        f"  seeds identical: {int(merged['seed_match'].sum())}/{n}; "
        f"selected rank identical: {int(merged['lv_match'].sum())}/{n}"
    )
    for statistic in STATISTICS:
        within = int(merged[f"{statistic}_match"].sum())
        same_p = int(merged[f"p_{statistic}_match"].sum())
        print(
            f"  {statistic}: observed within rtol {args.rtol:g} on {within}/{n} "
            f"(max abs diff {merged[f'abs_diff_{statistic}'].max():.3g}); p-values identical on {same_p}/{n}"
        )
    if not bool(merged["p_values_comparable"].all()):
        candidate_perms = sorted(int(v) for v in merged["permutations_candidate"].unique())
        reference_perms = sorted(int(v) for v in merged["permutations_reference"].unique())
        print(
            f"  p-values not comparable (permutation counts {candidate_perms} vs {reference_perms}); "
            "excluded from the verdict"
        )
    print(f"  all comparable fields match: {int(merged['all_match'].sum())}/{n}")
    print(f"Wrote {args.out}")
    return 0 if bool(merged["all_match"].all()) else 1


if __name__ == "__main__":
    raise SystemExit(main())
