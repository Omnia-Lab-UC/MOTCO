#!/usr/bin/env python3
"""Post-process MOTCO trajectory power study shards: merge + report.

Subcommands:

    merge   Combine all shard_*.jsonl files in a directory into one merged.jsonl.
    report  Build summaries, specificity matrix, power curves, Type I table,
            figures, and acceptance-target report from merged JSONL.

Typical interactive use after a cluster array completes:

    python scripts/motco_study.py merge --out-dir results/
    python scripts/motco_study.py report --config study.yaml --out-dir results/
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

from motco.simulations.grid import read_replicate_results
from motco.simulations.study.config import load_study_config
from motco.simulations.study.enumerate import enumerate_study
from motco.simulations.study.merge import discover_shard_paths, merge_shards
from motco.simulations.study.phase4 import summarize_attribution
from motco.simulations.study.report import (
    build_phase4_frames,
    build_report_frames,
    render_design_point_power,
    render_phase4_figures,
    render_power_curves,
    render_rank_ladder,
    render_specificity_matrix,
    render_type_i_plot,
    write_driver_report,
    write_phase4_report,
    write_report_contract,
    write_report_csvs,
)
from motco.simulations.study.summary import (
    summarize_combined_rule,
    summarize_study,
)
from motco.simulations.study.targets import (
    evaluate_anchor_decision,
    evaluate_design_point_decision,
    evaluate_rank_decision,
    evaluate_targets,
    write_anchor_decision,
    write_design_point_decision,
    write_rank_decision,
    write_target_report,
)


def _cmd_merge(args: argparse.Namespace) -> int:
    out_dir = Path(args.out_dir)
    shards = discover_shard_paths(out_dir)
    if not shards:
        print(f"No shard files found under {out_dir}.", file=sys.stderr)
        return 1
    merged_path = out_dir / "merged.jsonl"
    records = merge_shards(shards, out_path=merged_path)
    print(f"Merged {len(shards)} shards → {merged_path} ({len(records)} records).")
    return 0


def _cmd_report(args: argparse.Namespace) -> int:
    out_dir = Path(args.out_dir)
    merged_path = Path(args.merged) if args.merged else out_dir / "merged.jsonl"
    if not merged_path.exists():
        print(f"Merged JSONL not found: {merged_path}. Run `merge` first.", file=sys.stderr)
        return 1
    records = read_replicate_results(merged_path)
    config = load_study_config(args.config)
    per_stat = summarize_study(records, alpha=config.alpha)
    combined = summarize_combined_rule(records, alpha=config.alpha)
    frames = build_report_frames(per_stat, combined, records, alpha=config.alpha)
    report_dir = out_dir / "report"
    report_dir.mkdir(parents=True, exist_ok=True)
    csv_paths = write_report_csvs(frames, report_dir)
    figure_paths = {
        "specificity_matrix": render_specificity_matrix(frames.specificity_matrix, report_dir / "specificity_matrix.png"),
        "power_curves": render_power_curves(frames.power_curves, report_dir / "power_curves.png"),
        "type_i_plot": render_type_i_plot(frames.type_i_table, report_dir / "type_i.png", alpha=config.alpha),
    }
    evaluations = evaluate_targets(config.acceptance, per_stat, records)
    target_paths = write_target_report(evaluations, report_dir)

    design_paths: dict[str, Path] = {}
    if not frames.design_point_operating.empty:
        design_paths["design_point_power"] = render_design_point_power(
            frames.design_point_operating, report_dir / "design_point_power.png"
        )
    if config.acceptance.design_point is not None:
        decision = evaluate_design_point_decision(
            records, config.acceptance.design_point, alpha=config.alpha
        )
        design_paths.update(write_design_point_decision(decision, report_dir))
        print(f"Design-point decision: {decision.verdict.upper()} — {decision.rationale}")
    # Written only when the design grid declares the retained-rank axis.
    ladder = render_rank_ladder(frames.design_point_operating, report_dir / "rank_ladder.png")
    if ladder is not None:
        design_paths["rank_ladder"] = ladder
    if config.acceptance.rank_decision is not None:
        rank_decision = evaluate_rank_decision(records, config.acceptance.rank_decision, alpha=config.alpha)
        design_paths.update(write_rank_decision(rank_decision, report_dir))
        print(f"Retained-rank decision: {rank_decision.verdict.upper()} — {rank_decision.rationale}")
    # Written only when the config declares the anchor Type I rule.
    if config.acceptance.anchor_decision is not None:
        anchor_decision = evaluate_anchor_decision(records, config.acceptance.anchor_decision, alpha=config.alpha)
        design_paths.update(write_anchor_decision(anchor_decision, report_dir))
        print(f"Anchor Type I decision: {anchor_decision.verdict.upper()} — {anchor_decision.rationale}")

    contract = config.report_contract
    phase4_paths: dict[str, Path] = {}
    attribution_frame = None
    if config.acceptance.gate.enabled:
        expected_units = sum(cell.n_replicates for cell in enumerate_study(config).cells)
        frames4 = build_phase4_frames(
            config.acceptance.gate, per_stat, records, expected_units=expected_units
        )
        attribution_frame = frames4.attribution
        phase4_paths = {
            **write_phase4_report(frames4, report_dir),
            **render_phase4_figures(frames4, report_dir, contract=contract),
        }
        print(f"Phase 4 gate decision: {frames4.decision.decision.upper()} — {frames4.decision.rationale}")
        for run in frames4.decision.confirmation_runs:
            print(f"  confirmation re-run required: {run['action']}")

    # Written only when the config declares a report contract; configs without
    # one produce exactly the files above.
    contract_paths: dict[str, Path] = {}
    if contract is not None:
        if attribution_frame is None:
            attribution_frame = summarize_attribution(records)
        contract_paths = {
            "driver_report": write_driver_report(attribution_frame, contract, report_dir),
            "report_contract": write_report_contract(config, records, report_dir),
        }

    for key, path in {
        **csv_paths,
        **figure_paths,
        **target_paths,
        **design_paths,
        **phase4_paths,
        **contract_paths,
    }.items():
        print(f"  {key}: {path}")
    return 0


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(prog="motco_study", description=__doc__)
    sub = parser.add_subparsers(dest="command", required=True)

    p_merge = sub.add_parser("merge", help="Merge shard_*.jsonl into merged.jsonl")
    p_merge.add_argument("--out-dir", type=Path, required=True)
    p_merge.set_defaults(func=_cmd_merge)

    p_report = sub.add_parser("report", help="Build summaries, reports, and figures from merged JSONL")
    p_report.add_argument("--config", type=Path, required=True)
    p_report.add_argument("--out-dir", type=Path, required=True)
    p_report.add_argument("--merged", type=Path, default=None, help="Path to merged JSONL (defaults to <out-dir>/merged.jsonl).")
    p_report.set_defaults(func=_cmd_report)

    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    return int(args.func(args))


if __name__ == "__main__":
    sys.exit(main())
