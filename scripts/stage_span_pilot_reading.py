#!/usr/bin/env python3
"""Advisory reading of the stage-span pilot: propose the paper-grade follow-up or revise the construction.

Reads the four pilot reports (``openspec/changes/archive/2026-11-04-stage-span-constructions/``, design
D11) — ``{surgery,span}_{magnitude,geometry}/report`` under ``--results-dir`` —
and writes ``stage_span_pilot_reading.json``. From each report it reads
``design_point_operating.csv`` (rejection rates per design column, mode, effect
and statistic) and ``span_capture.csv`` (the capture shares per cell). The
bound is ``α + 2·√(α(1 − α)/R)``, with ``R`` the row's available replicates.

- **R1** — the stage-span ``translation`` control holds size: at every column
  and effect its ``delta``, ``angle`` and ``shape`` rates are at or below the bound.
- **R2** — the stage-span plant is captured: in every stage-span ``magnitude``,
  ``orientation`` and ``shape`` cell the median ``group_difference`` share is at
  least 0.90.
- **R3** — the stage-span modes are specific: in those cells every statistic
  other than the mode's target (``magnitude``→``delta``, ``orientation``→``angle``,
  ``shape``→``shape``) is at or below the bound.

The verdict is ``propose_paper_grade`` when all three hold and
``revise_construction`` otherwise, naming every failing criterion, cell and
statistic. A missing row or share counts as a failure (``unavailable``). Target
power is reported, not judged. Descriptive tables: both constructions' target
power and median capture per (column, mode, effect), the stage-span − surgery
``magnitude`` ``delta`` power difference with its pooled SE
(``√(SE₁² + SE₂²)``), and each profile's anchor size as read, unjudged.

Example:

    python scripts/stage_span_pilot_reading.py --results-dir results/stage-span-pilot-<date>
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
from collections.abc import Mapping, Sequence
from pathlib import Path
from typing import Any

import pandas as pd

#: Version of this reading rule; written to the output.
READING_RULE_VERSION = 1

ALPHA = 0.05
SE_MULTIPLIER = 2.0
CAPTURE_FLOOR = 0.90
STATISTICS: tuple[str, ...] = ("delta", "angle", "shape")
TARGETS: dict[str, str] = {"magnitude": "delta", "orientation": "angle", "shape": "shape"}
CONSTRUCTIONS: tuple[str, ...] = ("surgery", "span")
FAMILIES: tuple[str, ...] = ("magnitude", "geometry")
OUT_NAME = "stage_span_pilot_reading.json"
VERDICTS: tuple[str, ...] = ("propose_paper_grade", "revise_construction")

#: Columns that follow the design axes in the report's design-point tables.
_OPERATING_FIRST = "phase"
_CAPTURE_FIRST = "phase"


class ReadingError(ValueError):
    """Raised when a report cannot support the reading."""


def bound(n_available: int, *, alpha: float = ALPHA, k: float = SE_MULTIPLIER) -> float | None:
    """``α + k·√(α(1 − α)/R)``; ``None`` when no replicate is available."""

    if n_available <= 0:
        return None
    return alpha + k * math.sqrt(alpha * (1.0 - alpha) / n_available)


def _axes(frame: pd.DataFrame, first: str) -> list[str]:
    columns = list(frame.columns)
    return columns[: columns.index(first)] if first in columns else []


def _label(row: Mapping[str, Any], axes: Sequence[str]) -> str:
    def text(value: Any) -> str:
        if value is None or (isinstance(value, float) and math.isnan(value)):
            return "null"
        if isinstance(value, float) and value.is_integer():
            return str(int(value))
        return str(value)

    return "; ".join(f"{axis}={text(row[axis])}" for axis in axes) or "baseline"


def _number(value: Any) -> float | None:
    if value is None:
        return None
    number = float(value)
    return number if math.isfinite(number) else None


def load_report(report_dir: Path) -> dict[str, Any]:
    """The operating rows and capture rows of one report, keyed by (column, mode, effect)."""

    operating_path = report_dir / "design_point_operating.csv"
    capture_path = report_dir / "span_capture.csv"
    for path in (operating_path, capture_path):
        if not path.exists():
            raise ReadingError(f"Missing {path}; the reading needs the design-point and span-capture tables.")
    operating = pd.read_csv(operating_path)
    capture = pd.read_csv(capture_path)
    operating_axes = _axes(operating, _OPERATING_FIRST)
    capture_axes = _axes(capture, _CAPTURE_FIRST)
    if operating_axes != capture_axes:
        raise ReadingError(f"{report_dir}: design axes differ between the tables ({operating_axes} vs {capture_axes}).")

    rates: dict[tuple[str, str, float], dict[str, dict[str, Any]]] = {}
    for row in operating.to_dict("records"):
        key = (_label(row, operating_axes), str(row["trajectory_mode"]), float(row["effect_size"]))
        rates.setdefault(key, {})[str(row["statistic"])] = {
            "rate": _number(row["rejection_rate"]),
            "se": _number(row["monte_carlo_se"]),
            "n_available": int(row["n_available"]),
        }
    shares: dict[tuple[str, str, float], dict[str, Any]] = {}
    for row in capture.to_dict("records"):
        if row["phase"] not in ("power_primary", "power_design"):
            continue
        key = (_label(row, capture_axes), str(row["trajectory_mode"]), float(row["effect_size"]))
        shares[key] = {
            "construction": str(row["construction"]),
            "n_group_difference_available": int(row["n_group_difference_available"]),
            "median_group_difference": _number(row["median_group_difference"]),
            "median_pooled_stage": _number(row["median_pooled_stage"]),
            "median_rank": _number(row["median_rank"]),
            "median_population_latent_target": _number(row["median_population_latent_target"]),
            "median_intended_target": _number(row["median_intended_target"]),
        }
    return {"rates": rates, "shares": shares, "files": {str(p): _sha256(p) for p in (operating_path, capture_path)}}


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _within(entry: Mapping[str, Any] | None) -> tuple[bool, dict[str, Any]]:
    """Whether a rate holds the bound, with the observation."""

    if entry is None or entry.get("rate") is None or entry["n_available"] <= 0:
        return False, {"status": "unavailable"}
    limit = bound(entry["n_available"])
    assert limit is not None
    holds = entry["rate"] <= limit
    return holds, {
        "status": "holds" if holds else "exceeds",
        "rate": entry["rate"],
        "n_available": entry["n_available"],
        "bound": limit,
    }


def read(reports: Mapping[tuple[str, str], Mapping[str, Any]]) -> dict[str, Any]:
    """The reading over the four loaded reports, keyed by (construction, family)."""

    missing = [f"{c}_{f}" for c in CONSTRUCTIONS for f in FAMILIES if (c, f) not in reports]
    if missing:
        raise ReadingError(f"Missing report(s): {missing}.")

    failures: list[dict[str, Any]] = []
    criteria: dict[str, dict[str, Any]] = {name: {"holds": True, "observations": []} for name in ("R1", "R2", "R3")}

    def record(criterion: str, cell: Mapping[str, Any], holds: bool, observation: Mapping[str, Any]) -> None:
        entry = {**cell, **observation}
        criteria[criterion]["observations"].append(entry)
        if not holds:
            criteria[criterion]["holds"] = False
            failures.append({"criterion": criterion, **entry})

    span_geometry = reports[("span", "geometry")]
    translation = sorted(key for key in span_geometry["rates"] if key[1] == "translation")
    if not translation:
        record("R1", {"profile": "span_geometry", "trajectory_mode": "translation"}, False, {"status": "unavailable"})
    for column, mode, effect in translation:
        for statistic in STATISTICS:
            holds, observation = _within(span_geometry["rates"][(column, mode, effect)].get(statistic))
            cell = {"profile": "span_geometry", "column": column, "trajectory_mode": mode, "effect_size": effect}
            record("R1", {**cell, "statistic": statistic}, holds, observation)

    for family in FAMILIES:
        report = reports[("span", family)]
        profile = f"span_{family}"
        cells = sorted(key for key in report["rates"] if key[1] in TARGETS and key[2] != 0.0)
        if not cells:
            for criterion in ("R2", "R3"):
                record(criterion, {"profile": profile}, False, {"status": "unavailable"})
        for column, mode, effect in cells:
            cell = {"profile": profile, "column": column, "trajectory_mode": mode, "effect_size": effect}
            share = report["shares"].get((column, mode, effect)) or {}
            median = share.get("median_group_difference")
            if median is None:
                record("R2", cell, False, {"status": "unavailable"})
            else:
                holds = median >= CAPTURE_FLOOR
                record(
                    "R2",
                    cell,
                    holds,
                    {
                        "status": "holds" if holds else "below_floor",
                        "median_group_difference": median,
                        "floor": CAPTURE_FLOOR,
                        "n_available": share.get("n_group_difference_available"),
                    },
                )
            for statistic in STATISTICS:
                if statistic == TARGETS[mode]:
                    continue
                holds, observation = _within(report["rates"][(column, mode, effect)].get(statistic))
                record("R3", {**cell, "statistic": statistic}, holds, observation)

    verdict = "propose_paper_grade" if all(entry["holds"] for entry in criteria.values()) else "revise_construction"
    return {
        "rule_version": READING_RULE_VERSION,
        "verdict": verdict,
        "rationale": _rationale(verdict, failures),
        "alpha": ALPHA,
        "se_multiplier": SE_MULTIPLIER,
        "capture_floor": CAPTURE_FLOOR,
        "criteria": criteria,
        "failures": failures,
        "target_power_and_capture": _power_and_capture(reports),
        "magnitude_delta_difference": _magnitude_difference(reports),
        "anchor_size": _anchor_size(reports),
    }


def _rationale(verdict: str, failures: Sequence[Mapping[str, Any]]) -> str:
    if verdict == "propose_paper_grade":
        return (
            "R1, R2 and R3 hold: the stage-span translation holds size, every stage-span plant is captured at "
            f"a median of at least {CAPTURE_FLOOR:g}, and every off-target rate holds the bound."
        )
    named = sorted({str(failure["criterion"]) for failure in failures})
    return f"{', '.join(named)} fail on {len(failures)} observation(s); see failures."


def _power_and_capture(reports: Mapping[tuple[str, str], Mapping[str, Any]]) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    for family in FAMILIES:
        keys = sorted(
            {key for construction in CONSTRUCTIONS for key in reports[(construction, family)]["rates"]},
            key=lambda key: (key[0], key[1], key[2]),
        )
        for column, mode, effect in keys:
            if mode == "none":
                continue
            row: dict[str, Any] = {
                "family": family,
                "column": column,
                "trajectory_mode": mode,
                "effect_size": effect,
                "target_statistic": TARGETS.get(mode),
            }
            for construction in CONSTRUCTIONS:
                report = reports[(construction, family)]
                rates = report["rates"].get((column, mode, effect)) or {}
                share = report["shares"].get((column, mode, effect)) or {}
                target = TARGETS.get(mode)
                if target is not None:
                    entry = rates.get(target) or {}
                    row[f"{construction}_target_rate"] = entry.get("rate")
                    row[f"{construction}_target_se"] = entry.get("se")
                    row[f"{construction}_target_n"] = entry.get("n_available")
                else:
                    for statistic in STATISTICS:
                        row[f"{construction}_{statistic}_rate"] = (rates.get(statistic) or {}).get("rate")
                row[f"{construction}_median_group_difference"] = share.get("median_group_difference")
                row[f"{construction}_median_pooled_stage"] = share.get("median_pooled_stage")
                row[f"{construction}_median_rank"] = share.get("median_rank")
                row[f"{construction}_median_population_latent_target"] = share.get("median_population_latent_target")
                row[f"{construction}_median_intended_target"] = share.get("median_intended_target")
            rows.append(row)
    return rows


def _magnitude_difference(reports: Mapping[tuple[str, str], Mapping[str, Any]]) -> list[dict[str, Any]]:
    span, surgery = reports[("span", "magnitude")], reports[("surgery", "magnitude")]
    rows = []
    for key in sorted(k for k in span["rates"] if k[1] == "magnitude"):
        a = span["rates"][key].get("delta") or {}
        b = (surgery["rates"].get(key) or {}).get("delta") or {}
        if a.get("rate") is None or b.get("rate") is None:
            rows.append({"column": key[0], "effect_size": key[2], "status": "unavailable"})
            continue
        pooled = math.sqrt(float(a.get("se") or 0.0) ** 2 + float(b.get("se") or 0.0) ** 2)
        rows.append(
            {
                "column": key[0],
                "effect_size": key[2],
                "span_rate": a["rate"],
                "surgery_rate": b["rate"],
                "difference": a["rate"] - b["rate"],
                "pooled_se": pooled,
            }
        )
    return rows


def _anchor_size(reports: Mapping[tuple[str, str], Mapping[str, Any]]) -> list[dict[str, Any]]:
    rows = []
    for construction in CONSTRUCTIONS:
        for family in FAMILIES:
            for (column, mode, effect), rates in sorted(reports[(construction, family)]["rates"].items()):
                if mode != "none" or effect != 0.0:
                    continue
                row: dict[str, Any] = {"profile": f"{construction}_{family}", "column": column}
                for statistic in STATISTICS:
                    entry = rates.get(statistic) or {}
                    row[f"{statistic}_rate"] = entry.get("rate")
                    row[f"{statistic}_n"] = entry.get("n_available")
                    row[f"{statistic}_bound"] = bound(int(entry.get("n_available") or 0))
                rows.append(row)
    return rows


def run(results_dir: Path, out: Path | None = None) -> dict[str, Any]:
    reports = {
        (construction, family): load_report(results_dir / f"{construction}_{family}" / "report")
        for construction in CONSTRUCTIONS
        for family in FAMILIES
    }
    reading = read(reports)
    reading["inputs"] = {name: digest for report in reports.values() for name, digest in report["files"].items()}
    target = out or results_dir / OUT_NAME
    target.parent.mkdir(parents=True, exist_ok=True)
    target.write_text(json.dumps(reading, indent=1, sort_keys=True, allow_nan=False, default=_json_default) + "\n")
    return reading


def _json_default(value: Any) -> Any:
    raise TypeError(f"Not JSON-serializable: {value!r}")


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument(
        "--results-dir",
        type=Path,
        required=True,
        help="Directory holding {surgery,span}_{magnitude,geometry}/report.",
    )
    parser.add_argument("--out", type=Path, default=None, help=f"Output path (default: <results-dir>/{OUT_NAME}).")
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    try:
        reading = run(args.results_dir, args.out)
    except ReadingError as exc:
        raise SystemExit(str(exc)) from exc
    print(f"Stage-span pilot reading: {reading['verdict']} — {reading['rationale']}")
    for failure in reading["failures"]:
        print(f"  {failure}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
