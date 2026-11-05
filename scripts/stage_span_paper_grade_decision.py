#!/usr/bin/env python3
"""Predeclared decision of the stage-span paper grade (``stage-span-paper-grade``, design D5, D6 and D9).

Reads the three paper-grade reports under ``--results-dir``
(``span_geometry/report``, ``span_magnitude/report``, ``surgery_geometry/report``:
``design_point_operating.csv`` and ``span_capture.csv``) and writes
``stage_span_paper_grade_decision.json`` (rule version 1). Rows are keyed by a
canonical column name read from the design axes: ``n200`` and ``n80`` (three
stages, methylation + expression), ``seaad`` (the SEA-AD table),
``n200_three_blocks`` and ``n200_four_stages``.

**Bound rule (D5).** A criterion over m rejection rates holds only if

- the number of rates above ``α + 2·√(α(1 − α)/R)`` (R the row's available
  replicates) is at most the 95th percentile of Binomial(m, p₂), where p₂ is
  the exact probability that a Binomial(R, α)/R rate exceeds that bound
  (computed at the median R of the criterion's rows), and
- no rate is above the Bonferroni ceiling: the smallest rate a Binomial(R, α)
  rate exceeds with probability at most 0.05/m (computed per row).

Comparisons use the rejected counts, so they are exact. Rates within a column
share noise (matched seeds), so the count is approximate when exceedances
cluster; the JSON names every exceeding row.

**Criteria (D6).**

- **P1** size: the anchor of every column and the stage-span ``translation``
  cells, ``delta``/``angle``/``shape`` (m = 5·3 + 5·3·3 = 60).
- **P2** capture: the median ``group_difference`` share is at least 0.90 in
  every stage-span ``magnitude``, ``orientation`` and ``shape`` cell (no
  multiplicity).
- **P3** specificity: in those cells, every statistic other than the mode's
  target (``magnitude``→``delta``, ``orientation``→``angle``, ``shape``→``shape``),
  under the bound rule (m = 15·2 + 30·2 = 90).

A missing row, rate or share is ``unavailable`` and fails its criterion. The
verdict is the first failure in the order ``size_excess`` → ``capture_shortfall``
→ ``specificity_excess``, else ``sensitivity_established``; every failure is
listed.

**Descriptive (D9), not judged.** Target power ± SE per stage-span cell,
``informative_from`` (the smallest effect with target power ≥ 0.80) per column
and mode, capture quantiles, median rank, the median ``population_latent``
target beside the intended one, cells with failed units, and the visibility
contrast: the surgeries' target power and capture at the SEA-AD and four-stage
columns from this run, and at the other three columns the design-guidance
slice's surgery power (``--slice-report``) and the stage-span pilot's surgery
capture (``--pilot-capture``), each labelled with its source and grade, or
unavailable.

Example:

    python scripts/stage_span_paper_grade_decision.py --results-dir results/stage-span-paper-grade-<date>
"""

from __future__ import annotations

import argparse
import ast
import hashlib
import json
import math
import statistics
import sys
from collections.abc import Iterable, Mapping, Sequence
from pathlib import Path
from typing import Any

import pandas as pd
from scipy.stats import binom

sys.path.insert(0, str(Path(__file__).resolve().parent))
try:
    from stage_span_pilot_reading import _number, bound
finally:
    sys.path.pop(0)

#: Version of this decision rule; written to the output.
DECISION_RULE_VERSION = 1

ALPHA = 0.05
SE_MULTIPLIER = 2.0
FAMILY_ALPHA = 0.05
COUNT_QUANTILE = 0.95
CAPTURE_FLOOR = 0.90
POWER_FLOOR = 0.80
STATISTICS: tuple[str, ...] = ("delta", "angle", "shape")
TARGETS: dict[str, str] = {"magnitude": "delta", "orientation": "angle", "shape": "shape"}
COLUMNS: tuple[str, ...] = ("n200", "n80", "seaad", "n200_three_blocks", "n200_four_stages")
SURGERY_COLUMNS: tuple[str, ...] = ("seaad", "n200_four_stages")
SPAN_GEOMETRY_EFFECTS: tuple[float, ...] = (0.1, 0.25, 0.5)
SPAN_MAGNITUDE_EFFECTS: tuple[float, ...] = (0.02, 0.05, 0.1)
SURGERY_EFFECTS: tuple[float, ...] = (0.25, 0.5, 1.0)
PROFILES: tuple[str, ...] = ("span_geometry", "span_magnitude", "surgery_geometry")
OUT_NAME = "stage_span_paper_grade_decision.json"
VERDICTS: tuple[str, ...] = ("size_excess", "capture_shortfall", "specificity_excess", "sensitivity_established")
_VERDICT_OF = {"P1": "size_excess", "P2": "capture_shortfall", "P3": "specificity_excess"}

SLICE_REPORT = Path("results/design-guidance-slice-2026-10-29/core_geometry/report")
PILOT_CAPTURE = Path("results/stage-span-pilot-2026-11-04/surgery_geometry/report/span_capture.csv")
SLICE_SOURCE = {"source": "design-guidance slice 2026-10-29, core_geometry", "grade": "paper (500 x 999)"}
PILOT_SOURCE = {"source": "stage-span pilot 2026-11-04, surgery_geometry", "grade": "pilot (100 x 199)"}
RUN_SOURCE = {"source": "this run, surgery_geometry", "grade": "paper (500 x 999)"}

SEAAD_TABLE = ((11, 10, 28), (9, 10, 12))
_N, _SIZES, _STAGES, _LAYERS = (
    "generator.n_samples",
    "generator.group_stage_sizes",
    "generator.n_stages",
    "evaluation.integration_params.layers",
)
_DESIGN_PHASES = ("power_primary", "power_design")

#: Row key: (column, trajectory mode, effect size).
Key = tuple[str, str, float]


class DecisionError(ValueError):
    """Raised when a required report cannot support the decision."""


# --------------------------------------------------------------------------- #
# Bound rule (D5)
# --------------------------------------------------------------------------- #


def exceedance_threshold(n_available: int, *, alpha: float = ALPHA, k: float = SE_MULTIPLIER) -> int:
    """The largest rejected count whose rate is at or below the bound ``α + k·SE``."""

    limit = bound(n_available, alpha=alpha, k=k)
    if limit is None:
        raise DecisionError("The bound needs at least one available replicate.")
    return math.floor(limit * n_available)


def exceedance_probability(n_available: int, *, alpha: float = ALPHA, k: float = SE_MULTIPLIER) -> float:
    """p₂: the exact probability that a Binomial(R, α)/R rate exceeds ``α + k·SE``."""

    return float(binom.sf(exceedance_threshold(n_available, alpha=alpha, k=k), n_available, alpha))


def allowed_exceedances(m: int, p_exceed: float, *, quantile: float = COUNT_QUANTILE) -> int:
    """The ``quantile`` point of Binomial(m, p₂): the exceedance count a null run stays within."""

    return int(binom.ppf(quantile, m, p_exceed))


def ceiling_count(n_available: int, m: int, *, alpha: float = ALPHA, family_alpha: float = FAMILY_ALPHA) -> int:
    """The smallest rejected count a Binomial(R, α) count exceeds with probability at most ``family_alpha/m``."""

    return int(binom.ppf(1.0 - family_alpha / m, n_available, alpha))


def bonferroni_ceiling(n_available: int, m: int, *, alpha: float = ALPHA, family_alpha: float = FAMILY_ALPHA) -> float:
    """The Bonferroni ceiling as a rate."""

    return ceiling_count(n_available, m, alpha=alpha, family_alpha=family_alpha) / n_available


def judge_rates(observations: Sequence[Mapping[str, Any]]) -> dict[str, Any]:
    """Apply the bound rule to rows carrying ``n_rejected`` and ``n_available`` (or ``status: unavailable``)."""

    m = len(observations)
    rows: list[dict[str, Any]] = []
    available = [int(row["n_available"]) for row in observations if row.get("status") != "unavailable"]
    median_r = int(statistics.median_low(available)) if available else None
    p_exceed = exceedance_probability(median_r) if median_r else None
    allowed = allowed_exceedances(m, p_exceed) if p_exceed is not None and m else None
    exceed_count = 0
    above_ceiling = 0
    unavailable = 0
    for row in observations:
        entry = dict(row)
        if row.get("status") == "unavailable":
            unavailable += 1
            rows.append(entry)
            continue
        n_available, n_rejected = int(row["n_available"]), int(row["n_rejected"])
        threshold = exceedance_threshold(n_available)
        ceiling = ceiling_count(n_available, m)
        exceeds = n_rejected > threshold
        over = n_rejected > ceiling
        exceed_count += exceeds
        above_ceiling += over
        entry.update(
            {
                "rate": n_rejected / n_available,
                "bound": bound(n_available),
                "ceiling": ceiling / n_available,
                "status": "above_ceiling" if over else ("exceeds" if exceeds else "holds"),
            }
        )
        rows.append(entry)
    holds = unavailable == 0 and allowed is not None and exceed_count <= allowed and above_ceiling == 0
    return {
        "rule": "count_and_ceiling",
        "holds": holds,
        "m": m,
        "median_n_available": median_r,
        "p_exceed": p_exceed,
        "allowed_exceedances": allowed,
        "observed_exceedances": exceed_count,
        "ceiling_at_median": bonferroni_ceiling(median_r, m) if median_r and m else None,
        "n_above_ceiling": above_ceiling,
        "n_unavailable": unavailable,
        "observations": rows,
    }


# --------------------------------------------------------------------------- #
# Loading
# --------------------------------------------------------------------------- #


def _is_null(value: Any) -> bool:
    if value is None:
        return True
    if isinstance(value, float) and math.isnan(value):
        return True
    return isinstance(value, str) and value.strip() in ("", "None", "null", "nan")


def _parse(value: Any) -> Any:
    if isinstance(value, str):
        try:
            return ast.literal_eval(value)
        except (ValueError, SyntaxError):
            return value
    return value


def column_name(row: Mapping[str, Any]) -> str:
    """The canonical column of a report row, from whichever design axes the table carries."""

    sizes = None if _is_null(row.get(_SIZES)) else _parse(row.get(_SIZES))
    if sizes is not None:
        table = tuple(tuple(int(count) for count in group) for group in sizes)
        base = "seaad" if table == SEAAD_TABLE else f"table{table}"
    else:
        base = f"n{int(float(row[_N]))}"
    if _LAYERS in row and _is_null(row[_LAYERS]):
        base += "_three_blocks"
    if _STAGES in row and not _is_null(row[_STAGES]) and int(float(row[_STAGES])) == 4:
        base += "_four_stages"
    return base


def _axes(frame: pd.DataFrame) -> list[str]:
    columns = list(frame.columns)
    return columns[: columns.index("phase")] if "phase" in columns else []


def _design_rows(frame: pd.DataFrame) -> Iterable[dict[str, Any]]:
    axes = _axes(frame)
    for row in frame.to_dict("records"):
        if row.get("phase") not in _DESIGN_PHASES:
            continue
        yield {**row, "_column": column_name({axis: row[axis] for axis in axes})}


def load_rates(path: Path) -> dict[Key, dict[str, dict[str, Any]]]:
    """Rejection counts and rates per (column, mode, effect) and statistic from ``design_point_operating.csv``."""

    rates: dict[Key, dict[str, dict[str, Any]]] = {}
    for row in _design_rows(pd.read_csv(path)):
        key = (row["_column"], str(row["trajectory_mode"]), float(row["effect_size"]))
        rates.setdefault(key, {})[str(row["statistic"])] = {
            "n_rejected": int(row["n_rejected"]),
            "n_available": int(row["n_available"]),
            "n_replicates": int(row["n_replicates"]),
            "rate": _number(row["rejection_rate"]),
            "se": _number(row["monte_carlo_se"]),
            "median_selected_lv": _number(row.get("median_selected_lv")),
        }
    return rates


def load_capture(path: Path) -> dict[Key, dict[str, Any]]:
    """Capture shares and realized geometry per (column, mode, effect) from ``span_capture.csv``."""

    shares: dict[Key, dict[str, Any]] = {}
    for row in _design_rows(pd.read_csv(path)):
        key = (row["_column"], str(row["trajectory_mode"]), float(row["effect_size"]))
        shares[key] = {
            "construction": str(row["construction"]),
            "n_group_difference_available": int(row["n_group_difference_available"]),
            "q10_group_difference": _number(row["q10_group_difference"]),
            "median_group_difference": _number(row["median_group_difference"]),
            "q90_group_difference": _number(row["q90_group_difference"]),
            "median_pooled_stage": _number(row["median_pooled_stage"]),
            "median_rank": _number(row["median_rank"]),
            "median_population_latent_target": _number(row["median_population_latent_target"]),
            "median_intended_target": _number(row["median_intended_target"]),
        }
    return shares


def load_report(report_dir: Path) -> dict[str, Any]:
    operating, capture = report_dir / "design_point_operating.csv", report_dir / "span_capture.csv"
    for path in (operating, capture):
        if not path.exists():
            raise DecisionError(f"Missing {path}; the decision needs the design-point and span-capture tables.")
    return {
        "rates": load_rates(operating),
        "shares": load_capture(capture),
        "files": {str(path): _sha256(path) for path in (operating, capture)},
    }


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


# --------------------------------------------------------------------------- #
# Criteria
# --------------------------------------------------------------------------- #


def _rate_observation(rates: Mapping[Key, Mapping[str, Any]], key: Key, statistic: str, profile: str) -> dict[str, Any]:
    cell = {
        "profile": profile,
        "column": key[0],
        "trajectory_mode": key[1],
        "effect_size": key[2],
        "statistic": statistic,
    }
    entry = (rates.get(key) or {}).get(statistic)
    if entry is None or entry["n_available"] <= 0:
        return {**cell, "status": "unavailable"}
    return {**cell, "n_rejected": entry["n_rejected"], "n_available": entry["n_available"]}


def _stage_span_cells() -> list[tuple[str, Key]]:
    cells = [
        ("span_magnitude", (column, "magnitude", effect)) for column in COLUMNS for effect in SPAN_MAGNITUDE_EFFECTS
    ]
    cells += [
        ("span_geometry", (column, mode, effect))
        for column in COLUMNS
        for mode in ("orientation", "shape")
        for effect in SPAN_GEOMETRY_EFFECTS
    ]
    return cells


def decide(reports: Mapping[str, Mapping[str, Any]]) -> dict[str, Any]:
    missing = [profile for profile in PROFILES if profile not in reports]
    if missing:
        raise DecisionError(f"Missing report(s): {missing}.")
    geometry = reports["span_geometry"]

    size = [
        _rate_observation(geometry["rates"], (column, "none", 0.0), statistic, "span_geometry")
        for column in COLUMNS
        for statistic in STATISTICS
    ]
    size += [
        _rate_observation(geometry["rates"], (column, "translation", effect), statistic, "span_geometry")
        for column in COLUMNS
        for effect in SPAN_GEOMETRY_EFFECTS
        for statistic in STATISTICS
    ]
    p1 = judge_rates(size)

    capture_rows = []
    for profile, key in _stage_span_cells():
        cell = {"profile": profile, "column": key[0], "trajectory_mode": key[1], "effect_size": key[2]}
        share = reports[profile]["shares"].get(key) or {}
        median = share.get("median_group_difference")
        if median is None:
            capture_rows.append({**cell, "status": "unavailable"})
        else:
            capture_rows.append(
                {
                    **cell,
                    "status": "holds" if median >= CAPTURE_FLOOR else "below_floor",
                    "median_group_difference": median,
                    "n_available": share.get("n_group_difference_available"),
                }
            )
    p2 = {
        "rule": "every_cell_at_floor",
        "floor": CAPTURE_FLOOR,
        "holds": all(row["status"] == "holds" for row in capture_rows),
        "m": len(capture_rows),
        "observations": capture_rows,
    }

    specificity = [
        _rate_observation(reports[profile]["rates"], key, statistic, profile)
        for profile, key in _stage_span_cells()
        for statistic in STATISTICS
        if statistic != TARGETS[key[1]]
    ]
    p3 = judge_rates(specificity)

    criteria = {"P1": p1, "P2": p2, "P3": p3}
    failures = [
        {"criterion": name, **row}
        for name, criterion in criteria.items()
        for row in criterion["observations"]
        if row["status"] != "holds"
    ]
    failing = [name for name in ("P1", "P2", "P3") if not criteria[name]["holds"]]
    verdict = _VERDICT_OF[failing[0]] if failing else "sensitivity_established"
    return {
        "rule_version": DECISION_RULE_VERSION,
        "verdict": verdict,
        "failing_criteria": failing,
        "rationale": _rationale(verdict, criteria),
        "alpha": ALPHA,
        "se_multiplier": SE_MULTIPLIER,
        "family_alpha": FAMILY_ALPHA,
        "count_quantile": COUNT_QUANTILE,
        "capture_floor": CAPTURE_FLOOR,
        "criteria": criteria,
        "failures": failures,
        "notes": [
            "Rejections are p < alpha (the study's rule). The bound rule's count assumes independent rates; rates "
            "within a column share matched-seed noise, so exceedances can cluster. Every exceeding row is listed.",
            "A failure row with status 'exceeds' does not by itself fail a bound-rule criterion; the criterion fails "
            "when the count passes the allowed number, a rate passes the ceiling, or a row is unavailable.",
        ],
    }


def _rationale(verdict: str, criteria: Mapping[str, Mapping[str, Any]]) -> str:
    parts = []
    for name in ("P1", "P3"):
        c = criteria[name]
        parts.append(
            f"{name}: {c['observed_exceedances']} of {c['m']} rates above the bound (allowed "
            f"{c['allowed_exceedances']}), {c['n_above_ceiling']} above the ceiling, {c['n_unavailable']} unavailable"
        )
    below = sum(row["status"] != "holds" for row in criteria["P2"]["observations"])
    parts.insert(
        1, f"P2: {below} of {criteria['P2']['m']} cells below the {CAPTURE_FLOOR:g} capture floor or unavailable"
    )
    return f"{verdict}. " + "; ".join(parts) + "."


# --------------------------------------------------------------------------- #
# Descriptive outputs (D9)
# --------------------------------------------------------------------------- #


def target_power(reports: Mapping[str, Mapping[str, Any]]) -> list[dict[str, Any]]:
    rows = []
    for profile, key in _stage_span_cells():
        target = TARGETS[key[1]]
        entry = (reports[profile]["rates"].get(key) or {}).get(target) or {}
        share = reports[profile]["shares"].get(key) or {}
        rows.append(
            {
                "profile": profile,
                "column": key[0],
                "trajectory_mode": key[1],
                "effect_size": key[2],
                "target_statistic": target,
                "rate": entry.get("rate"),
                "se": entry.get("se"),
                "n_available": entry.get("n_available"),
                "q10_group_difference": share.get("q10_group_difference"),
                "median_group_difference": share.get("median_group_difference"),
                "q90_group_difference": share.get("q90_group_difference"),
                "median_rank": share.get("median_rank"),
                "median_population_latent_target": share.get("median_population_latent_target"),
                "median_intended_target": share.get("median_intended_target"),
            }
        )
    return rows


def informative_from(rows: Sequence[Mapping[str, Any]]) -> list[dict[str, Any]]:
    """Per column and mode, the smallest effect whose target rate is at least 0.80, or ``None``."""

    groups: dict[tuple[str, str], list[Mapping[str, Any]]] = {}
    for row in rows:
        groups.setdefault((row["column"], row["trajectory_mode"]), []).append(row)
    out = []
    for (column, mode), group in sorted(groups.items()):
        reached = sorted(
            row["effect_size"] for row in group if row.get("rate") is not None and row["rate"] >= POWER_FLOOR
        )
        out.append(
            {
                "column": column,
                "trajectory_mode": mode,
                "power_floor": POWER_FLOOR,
                "informative_from": reached[0] if reached else None,
                "effects": sorted(row["effect_size"] for row in group),
            }
        )
    return out


def incomplete_cells(reports: Mapping[str, Mapping[str, Any]]) -> list[dict[str, Any]]:
    """Cells where some replicates have no reading (failed units), per profile."""

    rows = []
    for profile in PROFILES:
        for key, entries in sorted(reports[profile]["rates"].items()):
            worst = min(entries.values(), key=lambda entry: entry["n_available"])
            if worst["n_available"] < worst["n_replicates"]:
                rows.append(
                    {
                        "profile": profile,
                        "column": key[0],
                        "trajectory_mode": key[1],
                        "effect_size": key[2],
                        "n_replicates": worst["n_replicates"],
                        "n_available": worst["n_available"],
                        "n_missing": worst["n_replicates"] - worst["n_available"],
                    }
                )
    return rows


def visibility_contrast(
    surgery: Mapping[str, Any],
    slice_rates: Mapping[Key, Mapping[str, Any]] | None,
    pilot_shares: Mapping[Key, Mapping[str, Any]] | None,
    *,
    slice_path: Path,
    pilot_path: Path,
) -> list[dict[str, Any]]:
    """The surgeries' target power and capture per column, labelled with source and grade."""

    rows = []
    for column in COLUMNS:
        for mode in ("orientation", "shape"):
            for effect in SURGERY_EFFECTS:
                key = (column, mode, effect)
                target = TARGETS[mode]
                row: dict[str, Any] = {
                    "column": column,
                    "trajectory_mode": mode,
                    "effect_size": effect,
                    "target_statistic": target,
                }
                if column in SURGERY_COLUMNS:
                    power_source, power_entry = RUN_SOURCE, (surgery["rates"].get(key) or {}).get(target)
                    capture_source, share = RUN_SOURCE, surgery["shares"].get(key)
                else:
                    power_source = {**SLICE_SOURCE, "path": str(slice_path)}
                    power_entry = None if slice_rates is None else (slice_rates.get(key) or {}).get(target)
                    capture_source = {**PILOT_SOURCE, "path": str(pilot_path)}
                    share = None if pilot_shares is None else pilot_shares.get(key)
                if power_entry is None:
                    row["power"] = {**power_source, "status": "unavailable"}
                else:
                    row["power"] = {
                        **power_source,
                        "status": "available",
                        "rate": power_entry.get("rate"),
                        "se": power_entry.get("se"),
                        "n_available": power_entry.get("n_available"),
                    }
                median = None if share is None else share.get("median_group_difference")
                if median is None:
                    row["capture"] = {**capture_source, "status": "unavailable"}
                else:
                    row["capture"] = {
                        **capture_source,
                        "status": "available",
                        "median_group_difference": median,
                        "q10_group_difference": share.get("q10_group_difference"),
                        "q90_group_difference": share.get("q90_group_difference"),
                    }
                rows.append(row)
    return rows


# --------------------------------------------------------------------------- #
# Driver
# --------------------------------------------------------------------------- #


def run(
    results_dir: Path,
    *,
    slice_report: Path = SLICE_REPORT,
    pilot_capture: Path = PILOT_CAPTURE,
    out: Path | None = None,
) -> dict[str, Any]:
    reports = {profile: load_report(results_dir / profile / "report") for profile in PROFILES}
    decision = decide(reports)
    power = target_power(reports)
    slice_operating = slice_report / "design_point_operating.csv"
    slice_rates = load_rates(slice_operating) if slice_operating.exists() else None
    pilot_shares = load_capture(pilot_capture) if pilot_capture.exists() else None
    decision["descriptive"] = {
        "target_power_and_capture": power,
        "informative_from": informative_from(power),
        "incomplete_cells": incomplete_cells(reports),
        "visibility_contrast": visibility_contrast(
            reports["surgery_geometry"], slice_rates, pilot_shares, slice_path=slice_operating, pilot_path=pilot_capture
        ),
    }
    inputs = {name: digest for report in reports.values() for name, digest in report["files"].items()}
    for path in (slice_operating, pilot_capture):
        inputs[str(path)] = _sha256(path) if path.exists() else "unavailable"
    decision["inputs"] = inputs
    target = out or results_dir / OUT_NAME
    target.parent.mkdir(parents=True, exist_ok=True)
    target.write_text(json.dumps(decision, indent=1, sort_keys=True, allow_nan=False) + "\n")
    return decision


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument(
        "--results-dir",
        type=Path,
        required=True,
        help="Directory holding {span_geometry,span_magnitude,surgery_geometry}/report.",
    )
    parser.add_argument("--slice-report", type=Path, default=SLICE_REPORT, help="The slice's core_geometry report.")
    parser.add_argument(
        "--pilot-capture", type=Path, default=PILOT_CAPTURE, help="The pilot's surgery span_capture.csv."
    )
    parser.add_argument("--out", type=Path, default=None, help=f"Output path (default: <results-dir>/{OUT_NAME}).")
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    try:
        decision = run(args.results_dir, slice_report=args.slice_report, pilot_capture=args.pilot_capture, out=args.out)
    except DecisionError as exc:
        raise SystemExit(str(exc)) from exc
    print(f"Stage-span paper-grade decision: {decision['rationale']}")
    for failure in decision["failures"]:
        print(f"  {failure}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
