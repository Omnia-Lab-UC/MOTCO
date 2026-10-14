#!/usr/bin/env python3
"""Select the follow-up change from the two net-displacement Type I anchor decisions.

Reads the ``report/anchor_decision.json`` written for each committed
net-displacement Type I profile (Phase 6 and Phase 5) and writes
``follow_up_decision.json``:

- ``paper_grade_power_comparison`` when the lead arm — cross-validated stage
  PLS with the net-displacement estimator — has verdict ``size_holds`` at both
  designs: compare it against ``pls_cv × pc1`` at paper grade next;
- ``diagnose_anchor_excess`` otherwise: diagnose the anchor excess first
  (endpoint-stage heteroscedasticity under full-residual permutation is the
  first hypothesis).

Per design it also records a descriptive localization over every
net-displacement column: ``every_column`` when all exceed, ``some_columns`` when
some do, ``none`` when none does. The selection is advisory; it states the rule
version and the inputs it read (path and sha256).

Example:

    python scripts/net_displacement_follow_up.py \
        --phase6 results/net-displacement-type-i-<date>/phase6/report/anchor_decision.json \
        --phase5 results/net-displacement-type-i-<date>/phase5/report/anchor_decision.json \
        --out results/net-displacement-type-i-<date>/follow_up_decision.json
"""

from __future__ import annotations

import argparse
import hashlib
import json
import sys
from collections.abc import Mapping
from pathlib import Path
from typing import Any

#: Version of this selection rule; written to the output.
FOLLOW_UP_RULE_VERSION = 1

METHOD = "evaluation.integration_method"
FORCED = "evaluation.integration_params.forced_components"
NCOMP = "evaluation.integration_params.n_components"
ESTIMATOR = "evaluation.orientation_estimator"
NET_DISPLACEMENT = "net_displacement"

#: Coordinates of the lead arm: cross-validated stage PLS × net displacement.
LEAD_COLUMN: dict[str, Any] = {METHOD: "pls", FORCED: None, NCOMP: None, ESTIMATOR: NET_DISPLACEMENT}

DESIGNS = ("phase6", "phase5")


class FollowUpError(ValueError):
    """Raised when an anchor decision cannot support the selection."""


def _matches(point: Mapping[str, Any], coordinates: Mapping[str, Any]) -> bool:
    return all(point.get(axis) == value for axis, value in coordinates.items())


def _design_summary(design: str, decision: Mapping[str, Any]) -> dict[str, Any]:
    points = list(decision.get("points") or [])
    leads = [point for point in points if _matches(point.get("design_point") or {}, LEAD_COLUMN)]
    if len(leads) != 1:
        raise FollowUpError(
            f"{design}: expected exactly one lead column {LEAD_COLUMN} in the anchor decision, found {len(leads)}."
        )
    (lead,) = leads
    if lead.get("verdict") not in ("size_holds", "excess"):
        raise FollowUpError(
            f"{design}: the lead column's anchor verdict is {lead.get('verdict')!r}; the selection needs "
            "'size_holds' or 'excess' (an unavailable anchor is a failed run, not a Type I result)."
        )
    net = [point for point in points if (point.get("design_point") or {}).get(ESTIMATOR) == NET_DISPLACEMENT]
    exceeding = [point for point in net if point.get("verdict") == "excess"]
    if not exceeding:
        localization = "none"
    elif len(exceeding) == len(net):
        localization = "every_column"
    else:
        localization = "some_columns"
    return {
        "lead_verdict": lead["verdict"],
        "lead_statistics": {
            entry["statistic"]: {
                key: entry.get(key) for key in ("rate", "monte_carlo_se", "n_available", "bound", "classification")
            }
            for entry in lead.get("statistics") or []
        },
        "localization": localization,
        "net_displacement_columns": [
            {"design_point": point.get("design_point"), "verdict": point.get("verdict")} for point in net
        ],
    }


def select_follow_up(decisions: Mapping[str, Mapping[str, Any]]) -> dict[str, Any]:
    """Pure selection over the parsed anchor decisions, keyed by design name."""

    missing = [design for design in DESIGNS if design not in decisions]
    if missing:
        raise FollowUpError(f"anchor decisions missing for design(s): {missing}.")
    designs = {design: _design_summary(design, decisions[design]) for design in DESIGNS}
    holds = all(summary["lead_verdict"] == "size_holds" for summary in designs.values())
    verdict = "paper_grade_power_comparison" if holds else "diagnose_anchor_excess"
    if holds:
        rationale = (
            "The lead arm (pls_cv x net_displacement) holds size on every named statistic at both designs: "
            "compare it against pls_cv x pc1 at paper grade."
        )
    else:
        exceeded = [design for design, summary in designs.items() if summary["lead_verdict"] == "excess"]
        rationale = (
            f"The lead arm (pls_cv x net_displacement) exceeds the anchor bound at {exceeded}: diagnose the "
            "anchor excess before any power comparison (first hypothesis: endpoint-stage heteroscedasticity "
            "under full-residual permutation)."
        )
    return {
        "rule_version": FOLLOW_UP_RULE_VERSION,
        "lead_column": LEAD_COLUMN,
        "verdict": verdict,
        "rationale": rationale,
        "designs": designs,
        "advisory": True,
    }


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def run(paths: Mapping[str, Path], out: Path) -> dict[str, Any]:
    decisions = {design: json.loads(Path(path).read_text(encoding="utf-8")) for design, path in paths.items()}
    payload = select_follow_up(decisions)
    payload["inputs"] = {design: {"path": str(path), "sha256": _sha256(Path(path))} for design, path in paths.items()}
    out = Path(out)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    return payload


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(prog="net_displacement_follow_up", description=__doc__.split("\n\n")[0])
    parser.add_argument("--phase6", type=Path, required=True, help="Phase 6 profile's anchor_decision.json.")
    parser.add_argument("--phase5", type=Path, required=True, help="Phase 5 profile's anchor_decision.json.")
    parser.add_argument("--out", type=Path, required=True, help="Where to write follow_up_decision.json.")
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    try:
        payload = run({"phase6": args.phase6, "phase5": args.phase5}, args.out)
    except FollowUpError as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 1
    print(f"Follow-up: {payload['verdict']} — {payload['rationale']}")
    for design, summary in payload["designs"].items():
        print(f"  {design}: lead {summary['lead_verdict']}, localization {summary['localization']}")
    print(f"  wrote {args.out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
