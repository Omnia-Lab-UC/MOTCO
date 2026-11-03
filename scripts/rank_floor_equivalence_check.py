#!/usr/bin/env python3
"""Does the rank floor change the selected rank at a pre-floor design? (no R, no RRPP)

The check behind ``openspec/changes/simulation-presentation-figures/`` (design D5).
A study run before the rank floor (``min_components=1``) may be presented as the
current measurement (``min_components=2``) only when the floor does not bind at
its design. The pooled PLS is refit once at the final rank, so equal final ranks
under the two floors mean identical latent scores, statistics and p-values; this
script compares them on rebuilt datasets.

For every baseline cell of ``--config`` (phases ``power_primary`` and
``type_i_baseline``) and every replicate index ``0 … --replicates − 1``, it
rebuilds the dataset from the config's matched seeds
(:func:`~motco.simulations.grid.derive_replicate_seed`), fits the PLS latent
space with the cell's evaluation parameters at ``min_components=1`` and ``2``
(:func:`~motco.simulations.evaluation.integrate_semisynthetic_dataset`) and
records the selected rank and the mean CV AUROC under each. Equal mean AUROC
means no outer fold picked rank 1 at floor 1 either; a differing AUROC with an
equal final rank changes nothing downstream and is reported, not judged.

Outputs under ``--out-dir``:

- ``rank_floor_equivalence.csv``: one row per dataset (cell, phase, mode,
  effect, replicate, rank and mean CV AUROC at floor 1 and floor 2);
- ``summary.json``: the dataset count, the rank mismatches, the AUROC
  differences, the one-sided 95% Clopper–Pearson upper bound on the mismatch
  rate and the verdict (``ranks_agree`` | ``ranks_differ``).

``--from-jsonl`` writes the same two files from the JSON lines of an earlier run
of this logic (keys ``cell_id``, ``phase``, ``mode``, ``effect``, ``rep``,
``lv1``, ``lv2``, ``auroc1``, ``auroc2``) instead of fitting anything.

Pin BLAS to one thread before a ``--workers N > 1`` run
(``OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1``).

Example:

    python scripts/rank_floor_equivalence_check.py \\
        --config examples/trajectory_power_study/net_displacement_power_phase5.json \\
        --replicates 8 --workers 8 --out-dir results/rank-floor-equivalence-<date>
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import sys
from concurrent.futures import ProcessPoolExecutor
from dataclasses import replace
from pathlib import Path
from typing import Any, Iterable, Sequence

from scipy.stats import beta

BASELINE_PHASES: tuple[str, ...] = ("power_primary", "type_i_baseline")
FLOORS: tuple[int, int] = (1, 2)
CSV_NAME = "rank_floor_equivalence.csv"
SUMMARY_NAME = "summary.json"
CSV_COLUMNS: tuple[str, ...] = (
    "cell_id",
    "phase",
    "mode",
    "effect",
    "replicate",
    "rank_floor_1",
    "rank_floor_2",
    "mean_cv_auroc_floor_1",
    "mean_cv_auroc_floor_2",
)
#: Confidence of the one-sided upper bound on the mismatch rate.
CONFIDENCE = 0.95


class EquivalenceCheckError(ValueError):
    """Raised when the check's inputs are unusable."""


def baseline_cells(config_path: str | Path) -> dict[str, Any]:
    """The config's baseline cells by ``cell_id``."""

    from motco.simulations.study.config import load_study_config
    from motco.simulations.study.enumerate import enumerate_study

    grid = enumerate_study(load_study_config(config_path))
    return {cell.cell_id: cell for cell in grid.cells if cell.phase in BASELINE_PHASES}


def ordered_tasks(cells: dict[str, Any], replicates: Sequence[int]) -> list[tuple[Any, int]]:
    """(cell, replicate) pairs, cells ordered by phase then metadata, replicates ascending."""

    ordered = sorted(cells.values(), key=lambda cell: (cell.phase, str(cell.metadata)))
    return [(cell, int(replicate)) for cell in ordered for replicate in replicates]


def evaluate_task(task: tuple[Any, int]) -> dict[str, Any]:
    """Rebuild one dataset and fit its latent space at both floors."""

    from motco.simulations.evaluation import integrate_semisynthetic_dataset
    from motco.simulations.grid import derive_replicate_seed
    from motco.simulations.semisynthetic import generate_semisynthetic_trajectory

    cell, replicate = task
    seed = derive_replicate_seed(cell, replicate)
    dataset = generate_semisynthetic_trajectory(replace(cell.generator_params, seed=seed))
    row: dict[str, Any] = {
        "cell_id": cell.cell_id,
        "phase": cell.phase,
        "mode": cell.metadata.get("trajectory_mode"),
        "effect": cell.metadata.get("effect_size"),
        "replicate": replicate,
    }
    for floor in FLOORS:
        latent = integrate_semisynthetic_dataset(dataset, replace(cell.evaluation_params, min_components=floor))
        row[f"rank_floor_{floor}"] = int(latent.metadata["selected_lv"])
        row[f"mean_cv_auroc_floor_{floor}"] = latent.metadata.get("cv_mean_auroc")
    return row


def rows_from_jsonl(path: str | Path) -> list[dict[str, Any]]:
    """Rows of the CSV from an earlier run's JSON lines; a line lacking a key is refused."""

    keys = ("cell_id", "phase", "mode", "effect", "rep", "lv1", "lv2", "auroc1", "auroc2")
    rows: list[dict[str, Any]] = []
    with open(path, encoding="utf-8") as handle:
        for number, line in enumerate(handle, start=1):
            if not line.strip():
                continue
            record = json.loads(line)
            missing = [key for key in keys if key not in record]
            if missing:
                raise EquivalenceCheckError(f"{path}:{number} lacks {missing}.")
            rows.append(
                {
                    "cell_id": record["cell_id"],
                    "phase": record["phase"],
                    "mode": record["mode"],
                    "effect": record["effect"],
                    "replicate": int(record["rep"]),
                    "rank_floor_1": int(record["lv1"]),
                    "rank_floor_2": int(record["lv2"]),
                    "mean_cv_auroc_floor_1": record["auroc1"],
                    "mean_cv_auroc_floor_2": record["auroc2"],
                }
            )
    return rows


def upper_bound(mismatches: int, count: int, confidence: float = CONFIDENCE) -> float:
    """One-sided Clopper–Pearson upper bound on a binomial rate (``1 − (1 − c)^(1/n)`` at zero)."""

    if count <= 0:
        raise EquivalenceCheckError("the upper bound needs at least one dataset.")
    if mismatches >= count:
        return 1.0
    return float(beta.ppf(confidence, mismatches + 1, count - mismatches))


def _auroc_differs(a: Any, b: Any) -> bool:
    if a is None or b is None:
        return (a is None) != (b is None)
    a, b = float(a), float(b)
    if math.isnan(a) or math.isnan(b):
        return math.isnan(a) != math.isnan(b)
    return a != b


def summarize(rows: Iterable[dict[str, Any]]) -> dict[str, Any]:
    """Count, rank mismatches, AUROC differences, the 95% upper bound and the verdict."""

    rows = list(rows)
    if not rows:
        raise EquivalenceCheckError("no datasets to summarize.")
    keys = [(row["cell_id"], int(row["replicate"])) for row in rows]
    if len(set(keys)) != len(keys):
        raise EquivalenceCheckError("a (cell, replicate) pair appears more than once.")
    mismatched = [row for row in rows if int(row["rank_floor_1"]) != int(row["rank_floor_2"])]
    auroc = [row for row in rows if _auroc_differs(row["mean_cv_auroc_floor_1"], row["mean_cv_auroc_floor_2"])]
    ranks_1 = sorted({int(row["rank_floor_1"]) for row in rows})
    ranks_2 = sorted({int(row["rank_floor_2"]) for row in rows})
    return {
        "count": len(rows),
        "cells": len({row["cell_id"] for row in rows}),
        "replicates": sorted({int(row["replicate"]) for row in rows}),
        "mismatches": len(mismatched),
        "mismatched": [{"cell_id": row["cell_id"], "replicate": int(row["replicate"])} for row in mismatched],
        "auroc_differences": len(auroc),
        "ranks_floor_1": ranks_1,
        "ranks_floor_2": ranks_2,
        "upper_bound_95": upper_bound(len(mismatched), len(rows)),
        "upper_bound_method": "one-sided 95% Clopper-Pearson on the rank-mismatch rate",
        "verdict": "ranks_agree" if not mismatched else "ranks_differ",
    }


def write_outputs(
    rows: Sequence[dict[str, Any]], out_dir: str | Path, *, provenance: dict[str, Any]
) -> dict[str, Path]:
    """Write the CSV and the summary (with the run's arguments under ``provenance``)."""

    summary = summarize(rows)
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    csv_path = out_dir / CSV_NAME
    with open(csv_path, "w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(CSV_COLUMNS))
        writer.writeheader()
        for row in rows:
            writer.writerow({column: row[column] for column in CSV_COLUMNS})
    summary_path = out_dir / SUMMARY_NAME
    summary_path.write_text(json.dumps({**summary, "provenance": provenance}, indent=2) + "\n", encoding="utf-8")
    return {"csv": csv_path, "summary": summary_path}


def run_check(config: str | Path, replicates: int, workers: int) -> list[dict[str, Any]]:
    """Fit every (baseline cell, replicate) at both floors, in task order."""

    if replicates < 1:
        raise EquivalenceCheckError("--replicates must be at least 1.")
    if workers < 1:
        raise EquivalenceCheckError("--workers must be at least 1.")
    tasks = ordered_tasks(baseline_cells(config), range(replicates))
    print(f"{len(tasks)} datasets", file=sys.stderr, flush=True)
    if workers == 1:
        return [evaluate_task(task) for task in tasks]
    with ProcessPoolExecutor(max_workers=workers) as pool:
        return list(pool.map(evaluate_task, tasks))


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(prog="rank_floor_equivalence_check", description=(__doc__ or "").split("\n\n")[0])
    parser.add_argument("--config", type=Path, help="Study config whose baseline cells are rebuilt.")
    parser.add_argument("--replicates", type=int, default=8, help="Replicate indices 0 … N − 1 per cell (default 8).")
    parser.add_argument("--workers", type=int, default=1, help="Worker processes (default 1).")
    parser.add_argument("--out-dir", type=Path, required=True, help="Where to write the CSV and the summary.")
    parser.add_argument("--from-jsonl", type=Path, help="Summarize an earlier run's JSON lines instead of fitting.")
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    try:
        if args.from_jsonl is not None:
            rows = rows_from_jsonl(args.from_jsonl)
            provenance: dict[str, Any] = {"from_jsonl": str(args.from_jsonl)}
        else:
            if args.config is None:
                raise EquivalenceCheckError("--config is required unless --from-jsonl is given.")
            rows = run_check(args.config, args.replicates, args.workers)
            provenance = {"config": str(args.config), "replicates": args.replicates, "workers": args.workers}
        paths = write_outputs(rows, args.out_dir, provenance=provenance)
    except EquivalenceCheckError as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 1
    summary = json.loads(paths["summary"].read_text(encoding="utf-8"))
    print(
        f"{summary['count']} datasets, {summary['mismatches']} rank mismatches "
        f"(95% upper bound {summary['upper_bound_95']:.4f}), {summary['auroc_differences']} AUROC differences: "
        f"{summary['verdict']}"
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
