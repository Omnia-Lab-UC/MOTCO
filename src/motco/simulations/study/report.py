"""Build paper-ready tables and figures from study summaries."""

from __future__ import annotations

import json
from collections.abc import Iterable, Sequence
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

from motco.simulations.grid import (
    SEED_FAMILY_KEY,
    SimulationReplicateResult,
    SimulationSummaryResult,
    summarize_realized_surgery,
)
from motco.simulations.study.config import Phase4GateConfig, ReportContract, StudyConfig
from motco.simulations.study.phase4 import (
    CHECKPOINT_ORDER as _CHECKPOINT_ORDER,
)
from motco.simulations.study.phase4 import (
    MEASUREMENT_SPACES as _MEASUREMENT_SPACES,
)
from motco.simulations.study.phase4 import (
    Phase4GateDecision,
    build_driver_report,
    build_operating_frame,
    evaluate_phase4_gate,
    localize_off_diagonal,
    summarize_attribution,
    summarize_pls_selection,
    summarize_realized_geometry,
)
from motco.simulations.study.spectrum import (
    RANK_AXIS,
    record_component_selection,
    record_forced_rank,
    resolve_operating_by_design_point,
    resolve_orientation_by_continuity,
    resolve_stage_signal_by_design_point,
    stratify_power_by_eigengap,
    summarize_config_spectrum,
)
from motco.simulations.study.summary import CombinedRuleSummary


class StudyReportError(ValueError):
    """Raised when summaries cannot be turned into a report."""


@dataclass(frozen=True)
class ReportFrames:
    """Paper-ready tables built from study summaries.

    ``config_spectrum`` and ``eigengap_stratified_power`` read the recorded
    latent configuration spectra; they default to empty frames so a caller
    reporting a record set written before that field existed still builds every
    pre-existing table unchanged. ``realized_surgery`` does the same for the
    requested-vs-realized surgery recorded in truth metadata.
    ``continuity_resolved_orientation`` is empty — and its CSV is not written —
    unless the record set actually spans more than one baseline continuity
    value, so a study that holds the axis fixed reports exactly what it did
    before the axis existed. ``design_point_operating`` likewise is empty — and
    unwritten — unless the record set contains design-grid cells, and so is
    ``design_point_stage_signal``, the per-column held-out CV stage AUROC that
    sits beside it.
    """

    specificity_matrix: pd.DataFrame
    power_curves: pd.DataFrame
    type_i_table: pd.DataFrame
    config_spectrum: pd.DataFrame = field(default_factory=pd.DataFrame)
    eigengap_stratified_power: pd.DataFrame = field(default_factory=pd.DataFrame)
    realized_surgery: pd.DataFrame = field(default_factory=pd.DataFrame)
    continuity_resolved_orientation: pd.DataFrame = field(default_factory=pd.DataFrame)
    design_point_operating: pd.DataFrame = field(default_factory=pd.DataFrame)
    design_point_stage_signal: pd.DataFrame = field(default_factory=pd.DataFrame)


@dataclass(frozen=True)
class Phase4ReportFrames:
    """Phase 4 structured outputs, all traceable to the merged records.

    Every claim a findings report makes should be readable out of one of these
    tables; they are the primary output and the prose is secondary.
    """

    operating: pd.DataFrame
    geometry: pd.DataFrame
    pls_selection: pd.DataFrame
    attribution: pd.DataFrame
    localization: pd.DataFrame
    gate: pd.DataFrame
    decision: Phase4GateDecision


def build_phase4_frames(
    gate: Phase4GateConfig,
    summaries: Sequence[SimulationSummaryResult],
    records: Sequence[SimulationReplicateResult],
    *,
    expected_units: int | None = None,
) -> Phase4ReportFrames:
    """Build every Phase 4 structured output plus the gate decision."""

    decision = evaluate_phase4_gate(gate, summaries, records, expected_units=expected_units)
    return Phase4ReportFrames(
        operating=build_operating_frame(summaries, records),
        geometry=summarize_realized_geometry(records),
        pls_selection=summarize_pls_selection(records),
        attribution=summarize_attribution(records),
        localization=localize_off_diagonal(records),
        gate=decision.to_frame(),
        decision=decision,
    )


def write_phase4_report(frames: Phase4ReportFrames, out_dir: Path) -> dict[str, Path]:
    """Write the Phase 4 tables as CSV and the gate decision as JSON."""

    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    paths = {
        "phase4_operating": out_dir / "phase4_operating.csv",
        "phase4_geometry": out_dir / "phase4_geometry.csv",
        "phase4_pls_selection": out_dir / "phase4_pls_selection.csv",
        "phase4_attribution": out_dir / "phase4_attribution.csv",
        "phase4_localization": out_dir / "phase4_localization.csv",
        "phase4_gate": out_dir / "phase4_gate.csv",
    }
    frames.operating.to_csv(paths["phase4_operating"], index=False)
    frames.geometry.to_csv(paths["phase4_geometry"], index=False)
    frames.pls_selection.to_csv(paths["phase4_pls_selection"], index=False)
    frames.attribution.to_csv(paths["phase4_attribution"], index=False)
    frames.localization.to_csv(paths["phase4_localization"], index=False)
    frames.gate.to_csv(paths["phase4_gate"], index=False)

    decision_path = out_dir / "phase4_gate_decision.json"
    decision_path.write_text(
        json.dumps(
            {
                "decision": frames.decision.decision,
                "rationale": frames.decision.rationale,
                "confirmation_runs": list(frames.decision.confirmation_runs),
                "observations": [
                    {
                        "rule": observation.rule,
                        "kind": observation.kind,
                        "trajectory_mode": observation.trajectory_mode,
                        "statistic": observation.statistic,
                        "cell_id": observation.cell_id,
                        "effect_size": observation.effect_size,
                        "met": observation.met,
                        "detail": observation.detail,
                        "observations": observation.observations,
                    }
                    for observation in frames.decision.observations
                ],
            },
            indent=2,
            default=_json_default,
        )
        + "\n",
        encoding="utf-8",
    )
    paths["phase4_gate_decision"] = decision_path
    return paths


def write_driver_report(attribution: pd.DataFrame, contract: ReportContract, out_dir: Path) -> Path:
    """Write ``driver_report.csv`` — the declared component only (see D2)."""

    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    path = out_dir / "driver_report.csv"
    build_driver_report(attribution, contract).to_csv(path, index=False)
    return path


def build_report_contract_echo(
    config: StudyConfig, records: Sequence[SimulationReplicateResult]
) -> dict[str, Any]:
    """Resolve the declared contract against the merged records.

    Beside the three declared fields the echo carries facts read from the
    records: the worker count they were produced with (``n_jobs``, asserted
    uniform — a mixed set means two signatures were merged and is refused), the
    configured worker count, and the shared zero-effect anchor's cell id, the
    modes whose ``0.00`` power point it resolves, and ``counted_as: 1``. One
    plain-language statement per contract item accompanies the values so the
    findings report can cite the echo rather than prose.
    """

    contract = config.report_contract
    if contract is None:
        raise StudyReportError("build_report_contract_echo requires a config with a report_contract block.")

    # The evaluator records ``params.n_jobs`` verbatim; ``None`` means serial and
    # is the same measurement as 1, so both normalise to 1. A record without the
    # key predates the field (or came from an evaluator that does not write it).
    n_jobs_values = sorted(
        {
            1 if record.runtime_metadata["n_jobs"] is None else int(record.runtime_metadata["n_jobs"])
            for record in records
            if record.status == "completed" and "n_jobs" in record.runtime_metadata
        }
    )
    if len(n_jobs_values) > 1:
        raise StudyReportError(
            f"merged records carry more than one n_jobs value: {n_jobs_values}; the contract echo "
            "cannot name one worker count for the run."
        )
    n_jobs = n_jobs_values[0] if n_jobs_values else None
    # ``evaluation.n_jobs`` is Optional in the evaluator's dataclass (None = serial);
    # RRPP records the resolved integer, so normalise the same way here.
    config_n_jobs = 1 if config.evaluation.n_jobs is None else int(config.evaluation.n_jobs)

    anchors: dict[str, list[str]] = {}
    for record in records:
        meta = record.cell_metadata
        if record.phase == "power_primary" and meta.get("zero_effect_anchor"):
            anchors.setdefault(record.cell_id, [str(mode) for mode in meta.get("resolves_modes", [])])
    if len(anchors) > 1:
        raise StudyReportError(
            f"records contain more than one primary shared zero-effect anchor: {sorted(anchors)}."
        )
    anchor: dict[str, Any] | None = None
    if anchors:
        (cell_id, modes), = anchors.items()
        anchor = {"cell_id": cell_id, "resolves_modes": modes, "counted_as": 1}

    component = contract.driver_component
    statements = {
        "driver_component": (
            f"Driver tables (driver_report.csv) and the attribution figure present the {component!r} "
            "attribution component only. phase4_attribution.csv keeps every component."
        ),
        "cross_replicate_driver_agreement": (
            "Cross-replicate top-k Jaccard and sign agreement are descriptive: they remain in "
            "phase4_attribution.csv, are absent from driver_report.csv and the attribution figure, "
            "and are not a driver-stability claim. The driver set is re-drawn every replicate, so "
            "no cross-replicate stability claim is made."
        ),
        "n_jobs": (
            "Records carry no n_jobs value (written before the field existed or by an evaluator that "
            f"does not record it); the configured value is {config_n_jobs}."
            if n_jobs is None
            else (
                f"Every record was produced with n_jobs = {n_jobs}"
                + (
                    ", the configured value."
                    if n_jobs == config_n_jobs
                    else f", which differs from the configured value {config_n_jobs}."
                )
                + f" The run-time override policy is {contract.n_jobs_override!r}."
            )
        ),
        "zero_effect_anchor": (
            "No shared zero-effect anchor is present in the records."
            if anchor is None
            else (
                f"Cell {anchor['cell_id']} is the shared zero-effect anchor; every mode in "
                f"{anchor['resolves_modes']} reads its 0.00 power point from it, and it is counted as "
                "one measurement (one cell), not once per mode."
            )
        ),
    }
    return {
        "driver_component": component,
        "cross_replicate_driver_agreement": contract.cross_replicate_driver_agreement,
        "n_jobs_override": contract.n_jobs_override,
        "n_jobs": n_jobs,
        "config_n_jobs": config_n_jobs,
        "zero_effect_anchor": anchor,
        "statements": statements,
    }


def write_report_contract(
    config: StudyConfig, records: Sequence[SimulationReplicateResult], out_dir: Path
) -> Path:
    """Write ``report_contract.json`` (see :func:`build_report_contract_echo`)."""

    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    path = out_dir / "report_contract.json"
    payload = build_report_contract_echo(config, records)
    path.write_text(json.dumps(payload, indent=2, default=_json_default) + "\n", encoding="utf-8")
    return path


def _json_default(value: object) -> object:
    if isinstance(value, float) and not np.isfinite(value):
        return None
    if isinstance(value, np.generic):
        item = value.item()
        return None if isinstance(item, float) and not np.isfinite(item) else item
    raise TypeError(f"Object of type {type(value).__name__} is not JSON serializable")


def build_specificity_matrix(
    summaries: Sequence[SimulationSummaryResult],
    records: Sequence[SimulationReplicateResult],
) -> pd.DataFrame:
    """Mode × statistic matrix of rejection rates (rate ± SE) at the top effect size.

    Rows are trajectory modes (from cell metadata), columns are statistics. For modes
    with multiple effect sizes (power cells), the entry for that mode uses the
    *largest* effect size; for negative-control / null modes the entry uses null cells.
    """

    cell_meta = _cell_metadata_index(records)
    rows: list[dict] = []
    for summary in summaries:
        meta = cell_meta.get(summary.cell_id)
        if meta is None:
            continue
        mode = _resolve_mode(meta)
        effect_size = float(meta.get("effect_size", 0.0) or 0.0)
        varied_axis = meta.get("varied_axis")
        if varied_axis is not None:
            continue
        rows.append(
            {
                "trajectory_mode": mode,
                "statistic": summary.statistic,
                "effect_size": effect_size,
                "cell_id": summary.cell_id,
                "phase": summary.phase,
                "rejection_rate": summary.rejection_rate,
                "monte_carlo_se": summary.monte_carlo_se,
                "available_replicates": summary.available_replicates,
                "completed_replicates": summary.completed_replicates,
            }
        )
    if not rows:
        return pd.DataFrame(
            columns=[
                "trajectory_mode",
                "statistic",
                "effect_size",
                "rejection_rate",
                "monte_carlo_se",
                "available_replicates",
            ]
        )
    frame = pd.DataFrame(rows)
    # for each (mode, statistic), keep the largest effect_size row
    idx = frame.groupby(["trajectory_mode", "statistic"])["effect_size"].idxmax()
    frame = frame.loc[idx].reset_index(drop=True)
    return frame.sort_values(["trajectory_mode", "statistic"]).reset_index(drop=True)


def build_power_curves(
    summaries: Sequence[SimulationSummaryResult],
    records: Sequence[SimulationReplicateResult],
) -> pd.DataFrame:
    """One row per (mode, statistic, effect_size) over baseline (non-OFAT) cells.

    Each row carries its cell's ``censored_fraction`` and
    ``duplicate_construction`` flag, so a curve point built on a censored — or
    duplicated — construction can be annotated rather than plotted as if the
    requested effect were the realized one.
    """

    cell_meta = _cell_metadata_index(records)
    surgery = {
        summary.cell_id: summary for summary in summarize_realized_surgery(records)
    }
    duplicates = find_duplicate_constructions(records)
    rows: list[dict] = []
    for summary in summaries:
        meta = cell_meta.get(summary.cell_id)
        if meta is None:
            continue
        if meta.get("varied_axis") is not None:
            continue
        mode = _resolve_mode(meta)
        effect_size = meta.get("effect_size")
        if effect_size is None:
            continue
        rows.append(
            {
                "trajectory_mode": mode,
                "statistic": summary.statistic,
                "effect_size": float(effect_size),
                "cell_id": summary.cell_id,
                "phase": summary.phase,
                "rejection_rate": summary.rejection_rate,
                "monte_carlo_se": summary.monte_carlo_se,
                "available_replicates": summary.available_replicates,
                "completed_replicates": summary.completed_replicates,
                "censored_fraction": (
                    surgery[summary.cell_id].censored_fraction
                    if summary.cell_id in surgery
                    else None
                ),
                "duplicate_construction": summary.cell_id in duplicates,
            }
        )
    frame = pd.DataFrame(
        rows,
        columns=[
            "trajectory_mode",
            "statistic",
            "effect_size",
            "cell_id",
            "phase",
            "rejection_rate",
            "monte_carlo_se",
            "available_replicates",
            "completed_replicates",
            "censored_fraction",
            "duplicate_construction",
        ],
    )
    if frame.empty:
        return frame
    return frame.sort_values(
        ["trajectory_mode", "statistic", "effect_size"]
    ).reset_index(drop=True)


#: A same-family cell pair sharing more than this fraction of identical realized
#: constructions is reported as a duplicated construction, not two measurements.
DUPLICATE_CONSTRUCTION_THRESHOLD = 0.05


def build_realized_surgery(records: Sequence[SimulationReplicateResult]) -> pd.DataFrame:
    """Per-cell requested-vs-realized surgery, with duplicate-construction flags.

    One row per cell: the nominal (requested) surgery size, the realized size
    (mean/min/max), and the fraction of replicates whose surgery was censored.
    Cells whose mode performs no pool-limited surgery carry ``NA`` in those
    columns rather than zeros.

    ``duplicate_of`` names the same-family cell this one shares its realized
    construction with, when more than
    ``DUPLICATE_CONSTRUCTION_THRESHOLD`` of replicate indices censored to the
    same realized size — the condition under which matched-seed cells generate
    byte-identical datasets and must not be read as independent measurements.
    """

    rows = [
        {
            "cell_id": summary.cell_id,
            "phase": summary.phase,
            "trajectory_mode": summary.trajectory_mode,
            "effect_size": summary.effect_size,
            "completed_replicates": summary.completed_replicates,
            "nominal_size": summary.nominal_size,
            "realized_mean": summary.realized_mean,
            "realized_min": summary.realized_min,
            "realized_max": summary.realized_max,
            "censored_replicates": summary.censored_replicates,
            "censored_fraction": summary.censored_fraction,
        }
        for summary in summarize_realized_surgery(records)
    ]
    frame = pd.DataFrame(
        rows,
        columns=[
            "cell_id",
            "phase",
            "trajectory_mode",
            "effect_size",
            "completed_replicates",
            "nominal_size",
            "realized_mean",
            "realized_min",
            "realized_max",
            "censored_replicates",
            "censored_fraction",
        ],
    )
    duplicates = find_duplicate_constructions(records)
    frame["duplicate_of"] = frame["cell_id"].map(
        lambda cell_id: ", ".join(sorted(duplicates.get(cell_id, ()))) or None
    )
    frame["duplicate_construction"] = frame["duplicate_of"].notna()
    if frame.empty:
        return frame
    return frame.sort_values(["trajectory_mode", "effect_size", "cell_id"]).reset_index(drop=True)


def find_duplicate_constructions(
    records: Sequence[SimulationReplicateResult],
    *,
    threshold: float = DUPLICATE_CONSTRUCTION_THRESHOLD,
) -> dict[str, set[str]]:
    """Cell pairs whose realized constructions coincide too often to be independent.

    Two cells in one matched-seed family draw the *same* generator seed at a
    given replicate index, and the generator's draw sequence depends only on that
    seed and the realized surgery size. So "both replicates censored to the same
    realized size" is exactly the condition for byte-identical datasets — it can
    be read straight off the records, with no regeneration.

    Returns a symmetric ``cell_id -> {cell_id, ...}`` map of the pairs whose
    coinciding fraction exceeds ``threshold``.
    """

    # (seed family, replicate index) -> [(cell_id, realized size), ...]
    by_slot: dict[tuple[str, int], list[tuple[str, int]]] = {}
    replicate_counts: dict[str, int] = {}
    for record in records:
        if record.status != "completed":
            continue
        transform = dict((record.truth_metadata or {}).get("transform") or {})
        if not transform.get("censored", False):
            continue
        realized = _realized_size(transform)
        if realized is None:
            continue
        family = str(record.cell_metadata.get(SEED_FAMILY_KEY) or record.cell_id)
        by_slot.setdefault((family, record.replicate_index), []).append(
            (record.cell_id, realized)
        )
        replicate_counts[record.cell_id] = replicate_counts.get(record.cell_id, 0) + 1

    coincidences: dict[tuple[str, str], int] = {}
    for entries in by_slot.values():
        for i, (cell_a, size_a) in enumerate(entries):
            for cell_b, size_b in entries[i + 1 :]:
                if cell_a == cell_b or size_a != size_b:
                    continue
                key = (cell_a, cell_b) if cell_a < cell_b else (cell_b, cell_a)
                coincidences[key] = coincidences.get(key, 0) + 1

    duplicates: dict[str, set[str]] = {}
    for (cell_a, cell_b), shared in coincidences.items():
        # denominator: the replicates the smaller cell could have shared
        available = min(replicate_counts[cell_a], replicate_counts[cell_b])
        if available == 0 or shared / available <= threshold:
            continue
        duplicates.setdefault(cell_a, set()).add(cell_b)
        duplicates.setdefault(cell_b, set()).add(cell_a)
    return duplicates


def _realized_size(transform: dict) -> int | None:
    for key in ("orientation_relocated", "translation_set_size", "shape_relocated"):
        if key in transform:
            return int(transform[key])
    return None


def build_type_i_table(
    summaries: Sequence[SimulationSummaryResult],
    combined: Sequence[CombinedRuleSummary],
    records: Sequence[SimulationReplicateResult],
) -> pd.DataFrame:
    """Null-cell table: per-statistic + combined-rule rejection rates."""

    cell_meta = _cell_metadata_index(records)
    by_cell: dict[str, dict] = {}
    for summary in summaries:
        meta = cell_meta.get(summary.cell_id)
        if meta is None:
            continue
        if not summary.phase.startswith("type_i_"):
            continue
        entry = by_cell.setdefault(
            summary.cell_id,
            {
                "cell_id": summary.cell_id,
                "phase": summary.phase,
                "trajectory_mode": _resolve_mode(meta),
                "varied_axis": meta.get("varied_axis"),
                "varied_value": meta.get("varied_value"),
                "completed_replicates": summary.completed_replicates,
            },
        )
        entry[f"{summary.statistic}_rate"] = summary.rejection_rate
        entry[f"{summary.statistic}_se"] = summary.monte_carlo_se
        entry[f"{summary.statistic}_available"] = summary.available_replicates

    combined_by_cell = {row.cell_id: row for row in combined}
    for cell_id, entry in by_cell.items():
        crow = combined_by_cell.get(cell_id)
        entry["combined_rate"] = crow.rejection_rate if crow else None
        entry["combined_se"] = crow.monte_carlo_se if crow else None
        entry["combined_available"] = crow.available_replicates if crow else 0

    frame = pd.DataFrame(list(by_cell.values()))
    if frame.empty:
        return frame
    return frame.sort_values(["trajectory_mode", "varied_axis", "varied_value", "cell_id"]).reset_index(drop=True)


def assert_production_component_selection(
    records: Sequence[SimulationReplicateResult],
) -> None:
    """Refuse to report records whose latent rank was forced without being declared.

    ``forced_components`` is a rank diagnostic (see ``scripts/latent_rank_probe.py``);
    a study report built over such records would present a hand-picked latent
    dimensionality as the study's operating point. The parameter signature
    already keeps forced-rank cells out of production shards, so this is the
    cheap second line: any record that says so explicitly is rejected here.

    The one legitimate exception is a study that *declares* the retained rank
    as a design-grid axis (``evaluation.integration_params.forced_components``):
    its forced columns are predeclared measurements of the same datasets at
    other ranks, sit at a design point carrying that rank, and are invisible to
    every baseline reader. Such a record passes only when its design point's
    rank equals the rank it recorded. Records written before
    ``component_selection`` existed carry no marker and pass unchanged.
    """

    offenders: set[str] = set()
    for record in records:
        mode = record_component_selection(record)
        if mode is None or mode == "cv":
            continue
        declared = record_forced_rank(record)
        recorded = (record.integration_metadata or {}).get("selected_lv")
        if declared is not None and recorded is not None and int(recorded) == declared:
            continue
        offenders.add(str(record.cell_id))
    if offenders:
        names = sorted(offenders)
        raise StudyReportError(
            "Study report requires cross-validated PLS component selection unless the retained rank "
            f"is a declared design axis; {len(names)} cell(s) recorded an undeclared forced latent rank: "
            f"{', '.join(names[:5])}" + (" ..." if len(names) > 5 else "")
        )


def build_report_frames(
    summaries: Sequence[SimulationSummaryResult],
    combined: Sequence[CombinedRuleSummary],
    records: Sequence[SimulationReplicateResult],
    *,
    alpha: float = 0.05,
) -> ReportFrames:
    """Build every always-on report frame from summaries and merged records."""

    assert_production_component_selection(records)
    return ReportFrames(
        specificity_matrix=build_specificity_matrix(summaries, records),
        power_curves=build_power_curves(summaries, records),
        type_i_table=build_type_i_table(summaries, combined, records),
        config_spectrum=summarize_config_spectrum(records),
        eigengap_stratified_power=stratify_power_by_eigengap(records, alpha=alpha),
        realized_surgery=build_realized_surgery(records),
        continuity_resolved_orientation=resolve_orientation_by_continuity(records, alpha=alpha),
        design_point_operating=resolve_operating_by_design_point(records, alpha=alpha),
        design_point_stage_signal=resolve_stage_signal_by_design_point(records),
    )


def write_report_csvs(
    frames: ReportFrames,
    out_dir: Path,
) -> dict[str, Path]:
    """Write the report frames as CSVs under ``out_dir``."""

    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    paths = {
        "specificity_matrix": out_dir / "specificity_matrix.csv",
        "power_curves": out_dir / "power_curves.csv",
        "type_i_table": out_dir / "type_i_table.csv",
        "config_spectrum": out_dir / "config_spectrum.csv",
        "eigengap_stratified_power": out_dir / "eigengap_stratified_power.csv",
        "realized_surgery": out_dir / "realized_surgery.csv",
    }
    frames.specificity_matrix.to_csv(paths["specificity_matrix"], index=False)
    frames.power_curves.to_csv(paths["power_curves"], index=False)
    frames.type_i_table.to_csv(paths["type_i_table"], index=False)
    frames.config_spectrum.to_csv(paths["config_spectrum"], index=False)
    frames.eigengap_stratified_power.to_csv(paths["eigengap_stratified_power"], index=False)
    frames.realized_surgery.to_csv(paths["realized_surgery"], index=False)
    # Written only when the study actually swept the continuity axis; its
    # absence is the report's statement that the axis was held fixed.
    if not frames.continuity_resolved_orientation.empty:
        paths["continuity_resolved_orientation"] = out_dir / "continuity_resolved_orientation.csv"
        frames.continuity_resolved_orientation.to_csv(
            paths["continuity_resolved_orientation"], index=False
        )
    # Written only when the study enumerated a design grid.
    if not frames.design_point_operating.empty:
        paths["design_point_operating"] = out_dir / "design_point_operating.csv"
        frames.design_point_operating.to_csv(paths["design_point_operating"], index=False)
    # The stage-signal table sits beside the operating table under the same
    # condition; the operating table itself is unchanged by it.
    if not frames.design_point_stage_signal.empty:
        paths["design_point_stage_signal"] = out_dir / "design_point_stage_signal.csv"
        frames.design_point_stage_signal.to_csv(paths["design_point_stage_signal"], index=False)
    return paths


def render_specificity_matrix(frame: pd.DataFrame, out_path: Path) -> Path:
    """Render the specificity matrix as a heatmap PNG."""

    import matplotlib.pyplot as plt

    out_path = Path(out_path)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    if frame.empty:
        fig, ax = plt.subplots(figsize=(4, 3))
        ax.set_title("Specificity matrix (empty)")
        ax.axis("off")
        fig.savefig(out_path, dpi=150, bbox_inches="tight")
        plt.close(fig)
        return out_path
    pivot = frame.pivot(index="trajectory_mode", columns="statistic", values="rejection_rate")
    fig, ax = plt.subplots(figsize=(1.2 * pivot.shape[1] + 2, 0.6 * pivot.shape[0] + 1.5))
    im = ax.imshow(pivot.to_numpy(dtype=float), cmap="viridis", vmin=0.0, vmax=1.0, aspect="auto")
    ax.set_xticks(range(pivot.shape[1]))
    ax.set_xticklabels(list(pivot.columns), rotation=45, ha="right")
    ax.set_yticks(range(pivot.shape[0]))
    ax.set_yticklabels(list(pivot.index))
    ax.set_xlabel("statistic")
    ax.set_ylabel("trajectory_mode")
    ax.set_title("Rejection rate by mode × statistic")
    values = pivot.to_numpy(dtype=float)
    for i in range(pivot.shape[0]):
        for j in range(pivot.shape[1]):
            value = values[i, j]
            if np.isfinite(value):
                ax.text(j, i, f"{value:.2f}", ha="center", va="center", color="white", fontsize=8)
    fig.colorbar(im, ax=ax, label="rejection rate")
    fig.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    return out_path


def render_type_i_plot(
    frame: pd.DataFrame,
    out_path: Path,
    *,
    alpha: float = 0.05,
) -> Path:
    """Plot Type I rates (per statistic + combined) vs alpha."""

    import matplotlib.pyplot as plt

    out_path = Path(out_path)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    if frame.empty:
        fig, ax = plt.subplots(figsize=(4, 3))
        ax.set_title("Type I table (empty)")
        ax.axis("off")
        fig.savefig(out_path, dpi=150, bbox_inches="tight")
        plt.close(fig)
        return out_path
    statistics = [c[:-5] for c in frame.columns if c.endswith("_rate")]
    fig, ax = plt.subplots(figsize=(6, 4))
    x = np.arange(len(frame))
    width = 0.8 / max(len(statistics), 1)
    for i, stat in enumerate(statistics):
        rate_col = f"{stat}_rate"
        se_col = f"{stat}_se"
        rates = frame[rate_col].to_numpy(dtype=float)
        ses = frame[se_col].to_numpy(dtype=float) if se_col in frame.columns else np.zeros_like(rates)
        ax.bar(x + i * width - 0.4 + width / 2, rates, width=width, yerr=ses, label=stat, capsize=2)
    ax.axhline(alpha, color="black", linestyle="--", linewidth=1, label=f"alpha={alpha}")
    labels = [
        f"{str(row.trajectory_mode)}\n({str(row.varied_axis) if row.varied_axis else 'baseline'})"
        for row in frame.itertuples(index=False)
    ]
    ax.set_xticks(x)
    ax.set_xticklabels(labels, rotation=45, ha="right", fontsize=8)
    ax.set_ylabel("rejection rate")
    ax.set_title("Type I rejection rates on null cells")
    ax.legend(loc="upper right", fontsize=8)
    fig.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    return out_path


def render_power_curves(
    frame: pd.DataFrame,
    out_path: Path,
) -> Path:
    """Render a panel grid of rejection-rate curves vs effect size."""

    import matplotlib.pyplot as plt

    out_path = Path(out_path)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    if frame.empty:
        fig, ax = plt.subplots(figsize=(4, 3))
        ax.set_title("Power curves (empty)")
        ax.axis("off")
        fig.savefig(out_path, dpi=150, bbox_inches="tight")
        plt.close(fig)
        return out_path

    modes = sorted(frame["trajectory_mode"].dropna().unique().tolist())
    statistics = sorted(frame["statistic"].dropna().unique().tolist())
    n_modes = len(modes)
    fig, axes = plt.subplots(
        1, n_modes, figsize=(3.4 * n_modes, 3.2), sharey=True, squeeze=False
    )
    for j, mode in enumerate(modes):
        ax = axes[0, j]
        sub = frame[frame["trajectory_mode"] == mode]
        for stat in statistics:
            stat_sub = sub[sub["statistic"] == stat].sort_values("effect_size")
            if stat_sub.empty:
                continue
            ax.errorbar(
                stat_sub["effect_size"].to_numpy(dtype=float),
                stat_sub["rejection_rate"].to_numpy(dtype=float),
                yerr=stat_sub["monte_carlo_se"].to_numpy(dtype=float),
                marker="o",
                label=stat,
                capsize=2,
            )
        ax.set_title(mode)
        ax.set_xlabel("effect size")
        ax.set_ylim(-0.02, 1.02)
        ax.grid(True, alpha=0.3)
        if j == 0:
            ax.set_ylabel("rejection rate")
        ax.legend(fontsize=8)
    fig.suptitle("Power curves by trajectory mode")
    fig.tight_layout(rect=(0, 0, 1, 0.95))
    fig.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    return out_path


def render_design_point_power(
    frame: pd.DataFrame,
    out_path: Path,
    *,
    trajectory_mode: str = "orientation",
    statistic: str = "angle",
    x_axis: str = "generator.n_samples",
    line_axis: str | None = "generator.baseline_continuity",
) -> Path:
    """Plot one statistic's power across the design grid, at the top effect per point.

    ``x_axis`` is drawn along the horizontal axis and ``line_axis`` selects one
    line per value; when either is not a column of ``frame`` the first (and
    second) design axes present are used instead. Each point is annotated with
    the median recorded eigengap so the geometry that governs the power is read
    beside it. Returns an empty placeholder figure when ``frame`` is empty.
    """

    import matplotlib.pyplot as plt

    out_path = Path(out_path)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    if frame.empty:
        return _empty_figure(out_path, "Design-point power (no design grid)")

    known = set(_DESIGN_POINT_NON_AXIS_COLUMNS)
    axes = [column for column in frame.columns if column not in known]
    if x_axis not in axes:
        x_axis = axes[0] if axes else ""
    if line_axis not in axes:
        remaining = [axis for axis in axes if axis != x_axis]
        line_axis = remaining[0] if remaining else None
    sub = frame[
        (frame["trajectory_mode"] == trajectory_mode)
        & (frame["statistic"] == statistic)
        & (frame["effect_size"].astype(float) > 0.0)
    ]
    if sub.empty or not x_axis:
        return _empty_figure(out_path, f"Design-point power ({trajectory_mode}/{statistic}: no rows)")
    # Axis values are plotted as numbers; a coordinate that is not numeric (the
    # ``null`` of a nested evaluation axis, say) has no position and is dropped
    # rather than crashing the figure. A rank-only grid thus draws its forced
    # columns and simply omits the cross-validated one — the rank ladder is the
    # figure that reads that column.
    sub = sub.assign(**{x_axis: pd.to_numeric(sub[x_axis], errors="coerce")}).dropna(subset=[x_axis])
    if line_axis:
        sub = sub.assign(**{line_axis: pd.to_numeric(sub[line_axis], errors="coerce")}).dropna(
            subset=[line_axis]
        )
    if sub.empty:
        return _empty_figure(out_path, f"Design-point power ({trajectory_mode}/{statistic}: no numeric axis)")
    group_keys = [x_axis] + ([line_axis] if line_axis else [])
    top = sub.loc[sub.groupby(group_keys)["effect_size"].idxmax()]

    fig, ax = plt.subplots(figsize=(5.2, 3.6))
    line_values = sorted(top[line_axis].dropna().unique().tolist()) if line_axis else [None]
    for line_index, value in enumerate(line_values):
        block = top if value is None else top[top[line_axis] == value]
        block = block.sort_values(x_axis)
        ax.errorbar(
            block[x_axis].to_numpy(dtype=float),
            block["rejection_rate"].to_numpy(dtype=float),
            yerr=block["monte_carlo_se"].to_numpy(dtype=float),
            marker="o",
            capsize=2,
            label=None if value is None or line_axis is None else f"{line_axis.split('.')[-1]} = {value:g}",
        )
        for _, row in block.iterrows():
            gap = row.get("median_eigengap")
            if gap is not None and np.isfinite(float(gap)):
                ax.annotate(
                    f"gap {float(gap):.3f}",
                    (float(row[x_axis]), float(row["rejection_rate"])),
                    textcoords="offset points",
                    xytext=(4, 4 + 9 * line_index),  # stagger lines so labels do not overlap
                    fontsize=7,
                )
    ax.set_xlabel(x_axis.split(".")[-1])
    ax.set_ylabel(f"{trajectory_mode}/{statistic} rejection rate (top effect)")
    ax.set_ylim(-0.02, 1.02)
    ax.grid(True, alpha=0.3)
    if line_axis:
        ax.legend(fontsize=8)
    ax.set_title("Power across the design grid")
    fig.tight_layout()
    fig.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    return out_path


def render_rank_ladder(
    frame: pd.DataFrame,
    out_path: Path,
    *,
    axis: str = RANK_AXIS,
    statistics: Sequence[str] = ("delta", "angle", "shape"),
) -> Path | None:
    """Plot each statistic's rejection rate against the retained rank, per mode.

    Built from the design-point operating frame alone. One panel per trajectory
    mode (the zero-effect anchor first, as ``none``); within a panel each
    statistic's rate at the mode's largest enumerated effect is drawn against
    the forced rank with Monte Carlo error bars. The cross-validated column
    (``axis`` value ``null``) is placed at its recorded median selected rank
    with a distinct marker and a ``CV`` label so it is never read as a forced
    rung. Returns ``None`` — and writes nothing — when ``axis`` is not a column
    of ``frame``, so a study without the rank axis produces no such figure.
    """

    import matplotlib.pyplot as plt

    if frame.empty or axis not in frame.columns:
        return None
    out_path = Path(out_path)
    out_path.parent.mkdir(parents=True, exist_ok=True)

    modes = sorted(frame["trajectory_mode"].dropna().unique().tolist(), key=lambda m: (m != "none", m))
    fig, axes = plt.subplots(1, len(modes), figsize=(3.6 * len(modes), 3.4), sharey=True, squeeze=False)
    forced_ranks = sorted(pd.to_numeric(frame[axis], errors="coerce").dropna().unique().tolist())
    for j, mode in enumerate(modes):
        ax = axes[0, j]
        block = frame[frame["trajectory_mode"] == mode]
        effects = block["effect_size"].astype(float)
        top = 0.0 if mode == "none" else float(effects.max())
        block = block[effects == top]
        for stat in statistics:
            rows = block[block["statistic"] == stat]
            if rows.empty:
                continue
            is_cv = rows[axis].isna()
            forced = rows[~is_cv].assign(rank=pd.to_numeric(rows.loc[~is_cv, axis], errors="coerce"))
            forced = forced.sort_values("rank")
            line = None
            if not forced.empty:
                line = ax.errorbar(
                    forced["rank"].to_numpy(dtype=float),
                    forced["rejection_rate"].to_numpy(dtype=float),
                    yerr=forced["monte_carlo_se"].to_numpy(dtype=float),
                    marker="o",
                    capsize=2,
                    label=stat,
                )
            cv = rows[is_cv]
            if not cv.empty:
                color = line[0].get_color() if line is not None else None
                x = pd.to_numeric(cv["median_selected_lv"], errors="coerce").to_numpy(dtype=float)
                y = cv["rejection_rate"].to_numpy(dtype=float)
                ax.errorbar(
                    x,
                    y,
                    yerr=cv["monte_carlo_se"].to_numpy(dtype=float),
                    fmt="*",
                    markersize=11,
                    capsize=2,
                    color=color,
                    markeredgecolor="black",
                    label=f"{stat} (CV)" if line is None else None,
                )
                for xi, yi in zip(x, y):
                    if np.isfinite(xi) and np.isfinite(yi):
                        ax.annotate("CV", (xi, yi), textcoords="offset points", xytext=(5, -10), fontsize=7)
        title = "zero-effect anchor (none)" if mode == "none" else f"{mode} @ e = {top:g}"
        ax.set_title(title, fontsize=9)
        ax.set_xlabel("retained rank")
        if forced_ranks:
            ax.set_xticks(forced_ranks)
        ax.set_ylim(-0.02, 1.02)
        ax.grid(True, alpha=0.3)
        if j == 0:
            ax.set_ylabel("rejection rate")
        ax.legend(fontsize=7)
    fig.suptitle("Retained-rank ladder (★ = cross-validated column at its median selected rank)")
    fig.tight_layout(rect=(0, 0, 1, 0.93))
    fig.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    return out_path


#: Every non-axis column of ``resolve_operating_by_design_point``'s frame.
_DESIGN_POINT_NON_AXIS_COLUMNS: tuple[str, ...] = (
    "phase",
    "is_baseline",
    "trajectory_mode",
    "effect_size",
    "statistic",
    "n_cells",
    "n_replicates",
    "n_available",
    "n_rejected",
    "rejection_rate",
    "monte_carlo_se",
    "n_eigengap_available",
    "mean_eigengap",
    "q33_eigengap",
    "median_eigengap",
    "q67_eigengap",
    "n_angle_null_available",
    "median_angle_null_q95",
    "iqr_angle_null_q95",
    "sd_angle_null_q95",
    "n_selected_lv_available",
    "median_selected_lv",
    "min_selected_lv",
    "max_selected_lv",
    "component_selection",
)


def render_geometry_checkpoints(frame: pd.DataFrame, out_path: Path) -> Path:
    """Plot each statistic's realized geometry across checkpoints, by mode.

    Every checkpoint is drawn on its own axis with its measurement space in the
    label, because a standardized-feature distance and a PLS-latent distance are
    not on one scale and must not be read off a shared axis.
    """

    import matplotlib.pyplot as plt

    out_path = Path(out_path)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    joint = (
        frame[(frame["scope"] == "joint") & (frame["statistic"].isin(("delta", "angle", "shape")))]
        if not frame.empty
        else frame
    )
    if joint.empty:
        return _empty_figure(out_path, "Realized geometry (empty)")

    checkpoints = [c for c in _CHECKPOINT_ORDER if c in set(joint["checkpoint"])]
    statistics = [s for s in ("delta", "angle", "shape") if s in set(joint["statistic"])]
    fig, axes = plt.subplots(
        len(statistics),
        len(checkpoints),
        figsize=(3.2 * len(checkpoints), 2.8 * len(statistics)),
        squeeze=False,
    )
    for row, statistic in enumerate(statistics):
        for column, checkpoint in enumerate(checkpoints):
            ax = axes[row][column]
            block = joint[(joint["statistic"] == statistic) & (joint["checkpoint"] == checkpoint)]
            for mode in sorted(block["trajectory_mode"].dropna().unique().tolist()):
                sub = block[block["trajectory_mode"] == mode].sort_values("effect_size")
                if sub.empty:
                    continue
                ax.plot(
                    sub["effect_size"].to_numpy(dtype=float),
                    sub["mean"].to_numpy(dtype=float),
                    marker="o",
                    label=mode,
                )
            space = _MEASUREMENT_SPACES.get(checkpoint, "unknown")
            if row == 0:
                ax.set_title(f"{checkpoint}\n({space})", fontsize=8)
            if column == 0:
                ax.set_ylabel(statistic)
            ax.set_xlabel("effect size", fontsize=8)
            ax.grid(True, alpha=0.3)
    axes[0][-1].legend(fontsize=7)
    fig.suptitle("Realized geometry by checkpoint (axes are per-space, not comparable)")
    fig.tight_layout(rect=(0, 0, 1, 0.94))
    fig.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    return out_path


def render_selected_components(frame: pd.DataFrame, out_path: Path) -> Path:
    """Plot the distribution of selected PLS component counts per cell."""

    import matplotlib.pyplot as plt

    out_path = Path(out_path)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    block = frame[frame["n_selected_lv"] > 0] if not frame.empty else frame
    if block.empty:
        return _empty_figure(out_path, "Selected PLS components (empty)")

    block = block.sort_values(["trajectory_mode", "effect_size"])
    labels = [
        f"{row.trajectory_mode}@{row.effect_size}" for row in block.itertuples(index=False)
    ]
    x = np.arange(len(block))
    means = block["selected_lv_mean"].to_numpy(dtype=float)
    lows = means - block["selected_lv_min"].to_numpy(dtype=float)
    highs = block["selected_lv_max"].to_numpy(dtype=float) - means
    fig, ax = plt.subplots(figsize=(max(6.0, 0.5 * len(block)), 3.6))
    ax.errorbar(x, means, yerr=np.vstack([lows, highs]), fmt="o", capsize=3)
    ax.set_xticks(x)
    ax.set_xticklabels(labels, rotation=45, ha="right", fontsize=7)
    ax.set_ylabel("selected latent variables")
    ax.set_title("Selected PLS components (mean with observed min/max)")
    ax.grid(True, axis="y", alpha=0.3)
    fig.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    return out_path


def render_attribution_stability(
    frame: pd.DataFrame, out_path: Path, *, contract: ReportContract | None = None
) -> Path:
    """Plot attribution stability for one component.

    Without a contract (the historical figure) the ``observed`` component's
    cross-replicate top-k Jaccard and sign agreement are drawn beside bootstrap
    sign stability under the title "Attribution stability". Under a contract
    only the *within-replicate* bootstrap series are drawn, for the declared
    component, and the title says so — cross-replicate agreement is descriptive
    and never plotted as stability.
    """

    import matplotlib.pyplot as plt

    out_path = Path(out_path)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    component = "observed" if contract is None else contract.driver_component
    series: tuple[tuple[str, str], ...]
    if contract is None:
        title = "Attribution stability"
        series = (
            ("top_k_jaccard", "cross-replicate top-k Jaccard"),
            ("sign_agreement", "cross-replicate sign agreement"),
            ("bootstrap_sign_stability_mean", "bootstrap sign stability"),
        )
    else:
        title = "Within-replicate bootstrap stability"
        series = (
            ("bootstrap_sign_stability_mean", "bootstrap sign stability"),
            ("bootstrap_top_k_frequency_mean", "bootstrap top-k selection frequency"),
        )
    block = frame[frame["component"] == component] if not frame.empty else frame
    if block.empty:
        return _empty_figure(out_path, f"{title} (empty)")

    block = block.sort_values(["trajectory_mode", "effect_size", "transition_id"])
    labels = [
        f"{row.trajectory_mode}@{row.effect_size}\n{row.transition_id}"
        for row in block.itertuples(index=False)
    ]
    x = np.arange(len(block))
    fig, ax = plt.subplots(figsize=(max(6.0, 0.7 * len(block)), 3.6))
    for column, label in series:
        if column not in block.columns:
            continue
        ax.plot(x, block[column].to_numpy(dtype=float), marker="o", label=label)
    ax.set_xticks(x)
    ax.set_xticklabels(labels, rotation=45, ha="right", fontsize=7)
    ax.set_ylim(-0.02, 1.02)
    ax.set_ylabel("agreement")
    ax.set_title(f"{title} ({component} component)")
    ax.grid(True, alpha=0.3)
    ax.legend(fontsize=7)
    fig.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    return out_path


def render_phase4_figures(
    frames: Phase4ReportFrames, out_dir: Path, *, contract: ReportContract | None = None
) -> dict[str, Path]:
    """Render every Phase 4 figure alongside the existing study outputs.

    ``contract`` is forwarded to :func:`render_attribution_stability`; without it
    every figure is byte-for-byte the pre-contract output.
    """

    out_dir = Path(out_dir)
    return {
        "phase4_geometry_plot": render_geometry_checkpoints(
            frames.geometry, out_dir / "phase4_geometry.png"
        ),
        "phase4_selected_components_plot": render_selected_components(
            frames.pls_selection, out_dir / "phase4_selected_components.png"
        ),
        "phase4_attribution_stability_plot": render_attribution_stability(
            frames.attribution, out_dir / "phase4_attribution_stability.png", contract=contract
        ),
    }


def _empty_figure(out_path: Path, title: str) -> Path:
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(figsize=(4, 3))
    ax.set_title(title)
    ax.axis("off")
    fig.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    return out_path


def _cell_metadata_index(records: Iterable[SimulationReplicateResult]) -> dict[str, dict]:
    index: dict[str, dict] = {}
    for record in records:
        if record.cell_id in index:
            continue
        index[record.cell_id] = dict(record.cell_metadata)
    return index


def _resolve_mode(meta: dict) -> str:
    mode = meta.get("trajectory_mode")
    if mode is not None:
        return str(mode)
    # baseline Type I cells from enumerate_type_i_grid carry no trajectory_mode metadata;
    # they are always the `none` mode (group_effect_size = 0).
    return "none"


__all__ = [
    "DUPLICATE_CONSTRUCTION_THRESHOLD",
    "Phase4ReportFrames",
    "ReportFrames",
    "StudyReportError",
    "assert_production_component_selection",
    "build_phase4_frames",
    "build_power_curves",
    "build_realized_surgery",
    "build_report_contract_echo",
    "build_report_frames",
    "build_specificity_matrix",
    "build_type_i_table",
    "find_duplicate_constructions",
    "render_power_curves",
    "render_specificity_matrix",
    "render_attribution_stability",
    "render_design_point_power",
    "render_geometry_checkpoints",
    "render_phase4_figures",
    "render_rank_ladder",
    "render_selected_components",
    "render_type_i_plot",
    "write_driver_report",
    "write_phase4_report",
    "write_report_contract",
    "write_report_csvs",
]
