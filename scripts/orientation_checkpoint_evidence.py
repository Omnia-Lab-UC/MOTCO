#!/usr/bin/env python3
"""Reproduce the orientation measurement probe's motivating evidence from a study ``merged.jsonl``.

The orientation measurement probe (``openspec/changes/archive/2026-10-12-orientation-measurement-probe/``)
is motivated by two observations on the Phase 6 small-n records:

- **Checkpoint table.** The orientation surgery's rotation is realized exactly
  in the population (arccos(1 − e)) and survives into the observed standardized
  space, but the stage-supervised PLS latent space collapses it — unless CV
  happens to keep more than ``n_stages`` components.
- **Flip-tail table.** The signed PC1 ``angle`` has a heavy null tail at the
  zero-effect anchor (q95 near 137°), so on an across-replicate oracle
  threshold the signed angle barely separates e = 1 from the anchor, while its
  folded counterpart ``min(a, 180 − a)`` does.

Both tables are read from committed-format records only:
``realized_geometry.checkpoints``, ``integration_metadata.selected_lv`` and
``p_values.angle``. One design column is selected by
``cell_metadata.design_point``: ``baseline`` is the column the ``power_primary``
records carry, or pass the design point as a JSON object.

Example:

    python scripts/orientation_checkpoint_evidence.py \
        --merged results/phase6-small-n-2026-10-09/merged.jsonl --column baseline
"""

from __future__ import annotations

import argparse
import json
import math
import sys
from collections import defaultdict
from pathlib import Path
from typing import Any, Iterable, Mapping, Sequence

import numpy as np

#: Checkpoint columns of the table. The population column reads the native
#: methylation block: the orientation surgery relocates methylation sites, so
#: that block carries its exact rotation, arccos(1 - e).
CHECKPOINT_COLUMNS: tuple[tuple[str, str, str], ...] = (
    ("population", "population_native", "methylation"),
    ("population_standardized", "population_standardized", "joint"),
    ("observed_standardized", "observed_standardized", "joint"),
)


def _design_key(record: Mapping[str, Any]) -> str:
    return json.dumps(record.get("cell_metadata", {}).get("design_point"), sort_keys=True)


def load_column(path: Path, column: str) -> list[dict[str, Any]]:
    """Completed records of the selected design column, in file order."""

    records = []
    with path.open(encoding="utf-8") as handle:
        for line in handle:
            if line.strip():
                record = json.loads(line)
                if record.get("status", "completed") == "completed":
                    records.append(record)
    if column == "baseline":
        keys = {_design_key(r) for r in records if r.get("phase") == "power_primary"}
        if len(keys) != 1:
            raise SystemExit(f"Expected one baseline design point among power_primary records, found {len(keys)}.")
        target = keys.pop()
    else:
        target = json.dumps(json.loads(column), sort_keys=True)
    selected = [
        r for r in records if r.get("phase", "").startswith("power") and _design_key(r) == target
    ]
    if not selected:
        raise SystemExit(f"No power records carry design point {target}.")
    return selected


def _angle(record: Mapping[str, Any], checkpoint: str, scope: str) -> float:
    entry = record.get("realized_geometry", {}).get("checkpoints", {}).get(checkpoint, {}).get(scope)
    value = None if entry is None else entry.get("angle")
    return math.nan if value is None else float(value)


def _median(values: Iterable[float]) -> float:
    array = np.asarray([v for v in values if math.isfinite(v)], dtype=float)
    return float(np.median(array)) if array.size else math.nan


def _rejects(record: Mapping[str, Any], alpha: float) -> bool | None:
    p_value = (record.get("p_values") or {}).get("angle")
    if p_value is None or not math.isfinite(float(p_value)):
        return None
    return float(p_value) < alpha


def _rate(records: Sequence[Mapping[str, Any]], alpha: float) -> float:
    decisions = [d for d in (_rejects(r, alpha) for r in records) if d is not None]
    return float(np.mean(decisions)) if decisions else math.nan


def _by_cell(records: list[dict[str, Any]]) -> dict[tuple[str, float], list[dict[str, Any]]]:
    cells: dict[tuple[str, float], list[dict[str, Any]]] = defaultdict(list)
    for record in records:
        metadata = record["cell_metadata"]
        cells[(str(metadata.get("trajectory_mode")), float(metadata.get("effect_size") or 0.0))].append(record)
    return cells


def _n_stages(records: Sequence[Mapping[str, Any]]) -> int:
    values = {int(r["truth_metadata"]["n_stages"]) for r in records}
    if len(values) != 1:
        raise SystemExit(f"Records mix stage counts {sorted(values)}.")
    return values.pop()


def checkpoint_table(
    records: list[dict[str, Any]],
    *,
    latent_checkpoint: str = "pls_latent",
    rank_split: int | None = None,
    alpha: float = 0.05,
) -> list[dict[str, Any]]:
    """Median realized ``angle`` per checkpoint for the anchor and every orientation cell.

    Rows are split by selected latent rank (``all``, ``<= rank_split``,
    ``> rank_split``; ``rank_split`` defaults to ``n_stages``), each with its
    ``angle`` rejection rate at ``p < alpha``.
    """

    split = _n_stages(records) if rank_split is None else rank_split
    rows = []
    for (mode, effect), cell in sorted(_by_cell(records).items()):
        if mode not in {"orientation", "none"}:
            continue
        subsets = {
            "all": cell,
            f"<={split}": [r for r in cell if int(r["integration_metadata"]["selected_lv"]) <= split],
            f">{split}": [r for r in cell if int(r["integration_metadata"]["selected_lv"]) > split],
        }
        for label, subset in subsets.items():
            if not subset:
                continue
            row: dict[str, Any] = {"mode": mode, "effect_size": effect, "rank": label, "n": len(subset)}
            for name, checkpoint, scope in CHECKPOINT_COLUMNS:
                row[name] = _median(_angle(r, checkpoint, scope) for r in subset)
            row["latent"] = _median(_angle(r, latent_checkpoint, "joint") for r in subset)
            row["angle_rejection_rate"] = _rate(subset, alpha)
            rows.append(row)
    return rows


def _fold(angles: np.ndarray) -> np.ndarray:
    return np.minimum(angles, 180.0 - angles)


def flip_tail_table(
    records: list[dict[str, Any]],
    *,
    latent_checkpoint: str = "pls_latent",
) -> dict[str, Any]:
    """Anchor ``angle`` quantiles and oracle separability of signed vs folded angles.

    Separability is the fraction of an orientation cell's replicates whose angle
    exceeds the anchor's across-replicate 95th percentile of the same angle.
    """

    cells = _by_cell(records)
    anchors = [key for key in cells if key[0] == "none"]
    if len(anchors) != 1:
        raise SystemExit(f"Expected one zero-effect anchor cell, found {len(anchors)}.")
    anchor = cells[anchors[0]]
    spaces = {"observed": ("observed_standardized", "joint"), "latent": (latent_checkpoint, "joint")}
    quantiles: dict[str, dict[str, float]] = {}
    thresholds: dict[tuple[str, str], float] = {}
    for space, (checkpoint, scope) in spaces.items():
        values = np.asarray([_angle(r, checkpoint, scope) for r in anchor], dtype=float)
        values = values[np.isfinite(values)]
        quantiles[space] = {f"q{int(q * 100)}": float(np.quantile(values, q)) for q in (0.5, 0.9, 0.95, 0.99)}
        thresholds[(space, "signed")] = float(np.quantile(values, 0.95))
        thresholds[(space, "folded")] = float(np.quantile(_fold(values), 0.95))
    separability = []
    for (mode, effect), cell in sorted(cells.items()):
        if mode != "orientation":
            continue
        row: dict[str, Any] = {"effect_size": effect, "n": len(cell)}
        for space, (checkpoint, scope) in spaces.items():
            values = np.asarray([_angle(r, checkpoint, scope) for r in cell], dtype=float)
            values = values[np.isfinite(values)]
            row[f"{space}_signed"] = float(np.mean(values > thresholds[(space, "signed")]))
            row[f"{space}_folded"] = float(np.mean(_fold(values) > thresholds[(space, "folded")]))
        separability.append(row)
    return {"anchor_quantiles": quantiles, "separability": separability}


def _fmt(value: Any, digits: int = 1) -> str:
    if isinstance(value, float):
        return "—" if not math.isfinite(value) else f"{value:.{digits}f}"
    return str(value)


def render_markdown(checkpoints: list[dict[str, Any]], flip: dict[str, Any], *, source: str) -> str:
    lines = [
        "# Orientation checkpoint evidence",
        "",
        f"Source: `{source}`. Generated by `scripts/orientation_checkpoint_evidence.py`.",
        "",
        "## Checkpoint table (median realized `angle`, degrees)",
        "",
        "| mode | e | rank | n | population | population std. | observed std. | latent | `angle` rate |",
        "|---|---|---|---|---|---|---|---|---|",
    ]
    for row in checkpoints:
        lines.append(
            f"| {row['mode']} | {row['effect_size']:.2f} | {row['rank']} | {row['n']} | "
            f"{_fmt(row['population'])} | {_fmt(row['population_standardized'])} | "
            f"{_fmt(row['observed_standardized'])} | {_fmt(row['latent'])} | "
            f"{_fmt(row['angle_rejection_rate'], 3)} |"
        )
    lines += [
        "",
        "## Flip-tail table",
        "",
        "Anchor `angle` quantiles (signed, degrees):",
        "",
        "| space | q50 | q90 | q95 | q99 |",
        "|---|---|---|---|---|",
    ]
    for space, quantiles in flip["anchor_quantiles"].items():
        lines.append(
            f"| {space} | {_fmt(quantiles['q50'])} | {_fmt(quantiles['q90'])} | "
            f"{_fmt(quantiles['q95'])} | {_fmt(quantiles['q99'])} |"
        )
    lines += [
        "",
        "Fraction of replicates above the anchor's 95th percentile (oracle separability):",
        "",
        "| e | n | observed signed | observed folded | latent signed | latent folded |",
        "|---|---|---|---|---|---|",
    ]
    for row in flip["separability"]:
        lines.append(
            f"| {row['effect_size']:.2f} | {row['n']} | {_fmt(row['observed_signed'], 3)} | "
            f"{_fmt(row['observed_folded'], 3)} | {_fmt(row['latent_signed'], 3)} | "
            f"{_fmt(row['latent_folded'], 3)} |"
        )
    return "\n".join(lines) + "\n"


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(description="Reproduce the orientation checkpoint and flip-tail tables.")
    p.add_argument("--merged", type=Path, required=True, help="Study merged.jsonl")
    p.add_argument(
        "--column",
        default="baseline",
        help="Design column: 'baseline' (the power_primary design point) or a JSON design point",
    )
    p.add_argument("--latent-checkpoint", default="pls_latent", help="Latent checkpoint name (default: pls_latent)")
    p.add_argument("--rank-split", type=int, default=None, help="Selected-rank split (default: n_stages)")
    p.add_argument("--alpha", type=float, default=0.05, help="Rejection threshold, p < alpha (default: 0.05)")
    p.add_argument("--out", type=Path, default=None, help="Write the Markdown here instead of stdout")
    p.add_argument("--json", type=Path, default=None, help="Also write both tables as JSON")
    return p


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    records = load_column(args.merged, args.column)
    checkpoints = checkpoint_table(
        records, latent_checkpoint=args.latent_checkpoint, rank_split=args.rank_split, alpha=args.alpha
    )
    flip = flip_tail_table(records, latent_checkpoint=args.latent_checkpoint)
    markdown = render_markdown(checkpoints, flip, source=str(args.merged))
    if args.out is None:
        sys.stdout.write(markdown)
    else:
        args.out.write_text(markdown, encoding="utf-8")
    if args.json is not None:
        payload = {"checkpoints": checkpoints, "flip_tail": flip}
        args.json.write_text(json.dumps(payload, indent=1, sort_keys=True, allow_nan=False) + "\n", encoding="utf-8")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
