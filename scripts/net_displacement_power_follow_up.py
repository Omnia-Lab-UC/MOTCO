#!/usr/bin/env python3
"""Select the orientation estimator from the two net-displacement power column decisions.

Reads the ``report/column_decision.json`` written for each committed
net-displacement power comparison profile (Phase 6, ``gain`` criterion; Phase 5,
``no_loss`` criterion) and writes ``comparison_decision.json``:

- ``adopt_net_displacement`` when the candidate — cross-validated stage PLS
  with the net-displacement estimator — has verdict ``candidate_qualifies`` at
  both designs: the follow-up is the production-default change;
- ``keep_pc1`` otherwise, naming per design the candidate's verdict and the
  failing criteria and pairs.

The selection is advisory: it does not change the production default. It states
the rule version and the inputs it read (path and sha256).

Example:

    python scripts/net_displacement_power_follow_up.py \
        --phase6 results/net-displacement-power-<date>/phase6/report/column_decision.json \
        --phase5 results/net-displacement-power-<date>/phase5/report/column_decision.json \
        --out results/net-displacement-power-<date>/comparison_decision.json
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
COMPARISON_RULE_VERSION = 1

METHOD = "evaluation.integration_method"
ESTIMATOR = "evaluation.orientation_estimator"

#: Coordinates of the candidate: cross-validated stage PLS × net displacement.
CANDIDATE: dict[str, Any] = {METHOD: "pls", ESTIMATOR: "net_displacement"}

DESIGNS = ("phase6", "phase5")

#: Column-decision verdicts the selection can read; anything else is a failed run.
JUDGED_VERDICTS = ("candidate_qualifies", "candidate_fails", "unavailable")

FOLLOW_UP_CHANGE = (
    "Flip SimulationEvaluationParams.orientation_estimator's production default to 'net_displacement' "
    "in a separate change that re-pins the committed signatures, as default-joint-magnitude-construction did."
)


class ComparisonError(ValueError):
    """Raised when a column decision cannot support the selection."""


def _design_summary(design: str, decision: Mapping[str, Any]) -> dict[str, Any]:
    if not isinstance(decision, Mapping) or "verdict" not in decision or "rule" not in decision:
        raise ComparisonError(f"{design}: input is not a column decision (missing 'verdict' or 'rule').")
    candidate = (decision.get("rule") or {}).get("candidate")
    if candidate != CANDIDATE:
        raise ComparisonError(
            f"{design}: the column decision's candidate is {candidate}; expected {CANDIDATE}."
        )
    verdict = decision["verdict"]
    if verdict not in JUDGED_VERDICTS:
        raise ComparisonError(
            f"{design}: the column decision's verdict is {verdict!r}; the selection needs one of "
            f"{list(JUDGED_VERDICTS)} (no candidate or no design grid is a failed run, not a comparison)."
        )
    columns = [column for column in decision.get("columns") or [] if column.get("role") == "candidate"]
    if len(columns) != 1:
        raise ComparisonError(f"{design}: expected exactly one candidate column, found {len(columns)}.")
    (column,) = columns
    criteria = {criterion["name"]: criterion for criterion in column.get("criteria") or []}
    target = (criteria.get("target") or {}).get("observations") or {}
    return {
        "verdict": verdict,
        "criterion": decision["rule"].get("criterion"),
        "failing": dict(decision.get("failing") or {}),
        "unavailable": [name for name, criterion in criteria.items() if criterion.get("passed") is None],
        "target": {
            key: target.get(key)
            for key in ("target", "rate", "se", "baseline_rate", "baseline_se", "difference", "pooled_se", "threshold")
        },
    }


def select_estimator(decisions: Mapping[str, Mapping[str, Any]]) -> dict[str, Any]:
    """Pure selection over the parsed column decisions, keyed by design name."""

    missing = [design for design in DESIGNS if design not in decisions]
    if missing:
        raise ComparisonError(f"column decisions missing for design(s): {missing}.")
    designs = {design: _design_summary(design, decisions[design]) for design in DESIGNS}
    adopt = all(summary["verdict"] == "candidate_qualifies" for summary in designs.values())
    if adopt:
        verdict = "adopt_net_displacement"
        rationale = (
            "The candidate (pls_cv x net_displacement) qualifies against pls_cv x pc1 at both designs: "
            f"target, anchor, every mandatory control and every protected pair. Follow-up: {FOLLOW_UP_CHANGE}"
        )
    else:
        verdict = "keep_pc1"
        reasons = []
        for design, summary in designs.items():
            if summary["verdict"] == "candidate_qualifies":
                continue
            if summary["failing"]:
                detail = "; ".join(f"{name}: {', '.join(pairs)}" for name, pairs in summary["failing"].items())
            else:
                detail = f"unevaluable criteria {summary['unavailable']}"
            reasons.append(f"{design} {summary['verdict']} ({detail})")
        rationale = (
            "The candidate (pls_cv x net_displacement) does not qualify at both designs — "
            + "; ".join(reasons)
            + ". The production orientation estimator stays PC1."
        )
    return {
        "rule_version": COMPARISON_RULE_VERSION,
        "candidate": CANDIDATE,
        "verdict": verdict,
        "rationale": rationale,
        "follow_up_change": FOLLOW_UP_CHANGE if adopt else None,
        "designs": designs,
        "advisory": True,
    }


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def run(paths: Mapping[str, Path], out: Path) -> dict[str, Any]:
    decisions: dict[str, Any] = {}
    for design, path in paths.items():
        try:
            decisions[design] = json.loads(Path(path).read_text(encoding="utf-8"))
        except (OSError, json.JSONDecodeError) as exc:
            raise ComparisonError(f"{design}: cannot read column decision {path}: {exc}") from exc
    payload = select_estimator(decisions)
    payload["inputs"] = {design: {"path": str(path), "sha256": _sha256(Path(path))} for design, path in paths.items()}
    out = Path(out)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    return payload


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(prog="net_displacement_power_follow_up", description=__doc__.split("\n\n")[0])
    parser.add_argument("--phase6", type=Path, required=True, help="Phase 6 profile's column_decision.json.")
    parser.add_argument("--phase5", type=Path, required=True, help="Phase 5 profile's column_decision.json.")
    parser.add_argument("--out", type=Path, required=True, help="Where to write comparison_decision.json.")
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    try:
        payload = run({"phase6": args.phase6, "phase5": args.phase5}, args.out)
    except ComparisonError as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 1
    print(f"Comparison: {payload['verdict']} — {payload['rationale']}")
    for design, summary in payload["designs"].items():
        print(f"  {design}: {summary['verdict']} ({summary['criterion']})")
    print(f"  wrote {args.out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
