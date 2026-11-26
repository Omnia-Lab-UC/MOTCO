#!/usr/bin/env python3
"""Unbalance at equal dispersion: whether the group-dispersion run's one failing reading replicates and what it follows.

The driver behind ``openspec/changes/unbalance-equal-dispersion-diagnosis/`` (spec
``unbalance-diagnosis-validation``, roadmap open item 16). The group-dispersion run
(``group_dispersion_check``, ``results/group-dispersion-2026-11-24/``) read
``unbalance_excess`` on one reading of 24: at equal dispersion, ``delta`` at n = 80 with a
first-group share of 0.8 (22 / 5 samples per stage) rejected 44 times of 500 against a
ceiling of 40. That design is both the most unbalanced and the smallest, so the committed
readings cannot separate small cells from unbalance. Every unit here is a null unit at equal
group dispersion.

**The designs** (base seed 2043, three stages, three blocks, δ = 2.0, the anchor cell only,
999 permutations; samples per stage in the larger / smaller group):

- ``n80_r0.8`` — the rerun gate: the group-dispersion run's design under its own name and
  seeds, 500 replicates, whose counts must reproduce the committed ``columns`` tables exactly.
- Seven new designs at 1,000 replicates: ``n80_r0.8_b`` (22 / 5, a fresh replicate of the
  failing design), ``n200_r0.925`` (62 / 5), ``n42_r0.5`` (7 / 7, the smallest balanced cell
  the production fit can run), ``n80_r0.74`` (20 / 7), ``n200_r0.895`` (60 / 7), ``n80_r0.65``
  (18 / 9) and ``n80_r0.5_b`` (14 / 13, a fresh balanced reference). 7,500 units.

Every generator is the committed ``n80_r0.8`` generator, read from the group-dispersion run's
``columns/settings.json``, with only ``n_samples`` and ``group_ratio`` changed. A new name gives
new dataset seeds (``latent_space_comparison.generator_seed``); the rerun's name gives the
committed datasets.

**The unit** is the group-dispersion run's (``group_dispersion_check.evaluate_unit``, imported
unchanged): the univariate baseline's unit (one dataset, one cross-validated production fit,
``delta``, ``angle`` and ``shape``, the split test, the interaction scan, the false-discovery
reading and ``pca_matched``) and then the within-group arm on a refit at the recorded rank.

**The corrected arms** (recorded, never judged). After the parent's unit the dataset is
regenerated, production is refitted at the recorded rank with no cross-validation, and
``delta``, ``angle`` and ``shape`` are read under ``RRPP(..., residual_correction="df")`` and
``"leverage"`` with unrestricted draws at the unit's RRPP seed and the harness's p-value rule.
Each correction rescales the reduced residuals once before the draws, so its draws use the
uncorrected draws' permutation indices. ``same_statistics`` checks that the refit's observed
statistics equal the production record's exactly; a unit where they differ is inconsistent for
these arms only and leaves their rates.

**Parents.** ``columns/settings.json``, ``columns/rate_table.csv``, ``columns/cell_table.csv``,
``columns/within_group_table.csv`` and ``group_dispersion_decision.json`` under
``results/group-dispersion-2026-11-24/`` are bound by sha256; a changed file is refused by name
in ``run``, ``merge`` and ``decide``.

Subcommands (outputs under ``--out-dir``):

- ``run`` — (a shard of) the units into ``records/shard_<i>_of_<N>.jsonl``; ``settings.json``
  holds the settings and their hash, a resume or merge under other settings exits naming the
  differing keys, and a resume skips complete units.
- ``merge`` — ``rate_table.csv`` (the parent's rows, the within-group rows and the corrected rows
  ``df_<statistic>`` and ``leverage_<statistic>`` under ``p < alpha``), ``cell_table.csv``,
  ``within_group_table.csv``, ``corrected_table.csv`` (per design: units, parent completions,
  refit completions and failures, inconsistent units, each correction's failures),
  ``accounting.json`` and ``reproduction.json`` (the rerun design against the committed
  ``columns/rate_table.csv`` and ``columns/within_group_table.csv``, exact on available and
  rejected counts for every arm, every differing row with both counts).
- ``decide`` — ``unbalance_diagnosis_decision.json`` (rule version 1), from the merged
  aggregates alone: the gate ``reproduction_failed``; size of ``delta``, ``angle`` and ``shape``
  at p < 0.05 on every new design, in three families under the stage-span paper grade's
  count-plus-ceiling rule (``judge_rates``, imported unchanged): ``balanced_small`` (``n42_r0.5``,
  m = 3), ``unbalanced_small`` (``n80_r0.8_b``, ``n200_r0.925``, ``n80_r0.74``, ``n200_r0.895``,
  m = 12) and ``reference`` (``n80_r0.65``, ``n80_r0.5_b``, m = 6); a reading with fewer than 900
  available units is unavailable and fails its family; the rerun design enters none. The
  verdict, the first that applies: ``reproduction_failed`` | ``small_cells`` |
  ``small_unbalanced_cells`` | ``moderate_excess`` | ``not_replicated``, with its consequence.
  Beside it, never changing it: the pooled ``delta`` count of ``n80_r0.8`` and ``n80_r0.8_b``,
  every design's counts by statistic, a label per family for every recorded arm (``production``,
  ``scan``, ``fdr``, ``total``, ``span``, ``complement``, ``pca_matched``, the within-group and the
  corrected statistics), the selected rank and held-out AUROC by design, and every arm's rate.

**The re-read** of the group-dispersion families, fixed before any unit and written into this
run's decision from the committed ``group_dispersion_decision.json`` (whose files are never
written to). Under ``small_cells`` or ``small_unbalanced_cells`` the five group-dispersion
designs with a 5-per-stage group (``n80_r0.8``, ``n80_r0.8_s1.5``, ``n80_r0.8_s2``,
``n80_r0.2_s1.5``, ``n80_r0.2_s2``) are set aside and the four families (21, 42, 54 and 54
readings) are judged again by the same rule; under ``not_replicated`` the three ``n80_r0.8``
readings of ``unbalanced_equal`` are replaced by ``n80_r0.8_b``'s and the families are judged in
full; under ``moderate_excess`` or ``reproduction_failed`` nothing is re-read and
``unbalance_excess`` stands. The re-read verdict follows the group-dispersion order, with that
verdict's consequence and the scope of the readings it covers.

With ``--workers`` above 1, pin BLAS to one thread per process
(``OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1``).

Examples:

    export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1

    # smoke
    python scripts/unbalance_diagnosis_check.py run --reps 2 --out-dir /tmp/ud --workers 7
    python scripts/unbalance_diagnosis_check.py merge --out-dir /tmp/ud
    python scripts/unbalance_diagnosis_check.py decide --out-dir /tmp/ud

    # one SLURM array task (scripts/unbalance_diagnosis_array.sbatch), then merge and decide
"""

from __future__ import annotations

import argparse
import json
import sys
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path
from time import perf_counter
from typing import Any, Iterable, Mapping, Sequence

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import group_dispersion_check as gd  # noqa: E402  (the unit, the within-group arm, merge helpers, the families)
import interaction_scan_check as isc  # noqa: E402  (the size observation of one rate row)
import latent_space_comparison as lsc  # noqa: E402  (seeds, settings binding, CSV)
import span_split_check as ssc  # noqa: E402  (parent loading, records, paths)
import stage_span_paper_grade_decision as grade  # noqa: E402  (the count-plus-ceiling rule)
import univariate_baseline_check as ubc  # noqa: E402  (the unit's fields and inputs, its tables)

sys.path.pop(0)

DRIVER = "unbalance_diagnosis_check"
DRIVER_VERSION = 1
RULE_VERSION = 1
PROFILE = "diagnosis"

BASE_SEED = gd.BASE_SEED
DEFAULT_PERMUTATIONS = gd.DEFAULT_PERMUTATIONS
RERUN_REPLICATES = 500
NEW_REPLICATES = 1000
#: A judged reading with fewer available units is unavailable and fails its family.
MIN_AVAILABLE = 900

REPO = Path(__file__).resolve().parents[1]
PARENT_DIR = REPO / "results" / "group-dispersion-2026-11-24"
PARENT_PROFILE = "columns"
PARENT_FILE_NAMES: tuple[str, ...] = (
    "columns/settings.json", "columns/rate_table.csv", "columns/cell_table.csv", "columns/within_group_table.csv",
    gd.DECISION_NAME,
)
PARENT_SOURCE = "group dispersion 2026-11-24, columns profile and decision"

#: The group-dispersion design rerun as the gate, under its own name and seeds.
RERUN_DESIGN = "n80_r0.8"
CELLS = gd.COLUMN_CELLS
CELL = CELLS[0][0]

#: design -> (n_samples, first-group share, larger and smaller group's samples per stage, family).
NEW_DESIGNS: dict[str, tuple[int, float, int, int, str]] = {
    "n80_r0.8_b": (80, 0.8, 22, 5, "unbalanced_small"),
    "n200_r0.925": (200, 0.925, 62, 5, "unbalanced_small"),
    "n42_r0.5": (42, 0.5, 7, 7, "balanced_small"),
    "n80_r0.74": (80, 0.74, 20, 7, "unbalanced_small"),
    "n200_r0.895": (200, 0.895, 60, 7, "unbalanced_small"),
    "n80_r0.65": (80, 0.65, 18, 9, "reference"),
    "n80_r0.5_b": (80, 0.5, 14, 13, "reference"),
}
RERUN_CELLS: tuple[int, int] = (22, 5)
FRESH_REPLICATE = "n80_r0.8_b"

STATISTICS: tuple[str, ...] = gd.STATISTICS
FAMILIES: tuple[str, ...] = ("balanced_small", "unbalanced_small", "reference")
FAMILY_M_EXPECTED: dict[str, int] = {"balanced_small": 3, "unbalanced_small": 12, "reference": 6}

CORRECTIONS: tuple[str, ...] = ("df", "leverage")
CORRECTED_RULE = ubc.RULE_STRICT
CORRECTED_ARMS: tuple[str, ...] = tuple(f"{c}_{statistic}" for c in CORRECTIONS for statistic in STATISTICS)
RATE_ARMS: tuple[str, ...] = (*gd.RATE_ARMS, *CORRECTED_ARMS)
RECORDED_ARMS: tuple[str, ...] = (*gd.RECORDED_ARMS, *CORRECTED_ARMS)
CORRECTED_COLUMNS: tuple[str, ...] = (
    *ubc.CELL_FIELDS, "units", "parent_completed", "refit_completed", "refit_failed", "inconsistent",
    *(f"{c}_failed" for c in CORRECTIONS),
)
WITHIN_COUNT_COLUMNS: tuple[str, ...] = ("units", "parent_completed", "within_completed", "within_failed",
                                         "inconsistent")

VERDICTS: tuple[str, ...] = (
    "reproduction_failed", "small_cells", "small_unbalanced_cells", "moderate_excess", "not_replicated",
)
_LIMIT = (
    " The readings hold for simulated Gaussian noise at equal group dispersion, 658 features on three blocks, "
    "δ = 2.0, three stages and no covariates; no balanced design below 7 samples per stage was read, and the "
    "re-read rule was written after the group-dispersion families' rates were known. A cohort's own table is read "
    "at its own design. No verdict changes what `motco analyze` runs."
)
_DOCUMENTATION = (
    " The documentation states the smallest cell size from which `delta`, `angle` and `shape` held and the designs "
    "and statistics that failed, and applies the re-read verdict's consequence with its scope; roadmap open item 16 "
    "moves to Done."
)
CONSEQUENCES: dict[str, str] = {
    "reproduction_failed": (
        "Nothing is read: the rerun of n80_r0.8 does not reproduce the group-dispersion run's committed counts, so "
        "the generator, the fit or the environment drifted. Find the drift before any rerun; a mismatch is "
        "diagnosed, never waived. Nothing is re-read and the group-dispersion verdict unbalance_excess stands."
    ),
    "small_cells": (
        "Small cells lose size with balanced groups, so unbalance is not needed: the balanced design at 7 samples "
        "per stage fails." + _DOCUMENTATION + _LIMIT
    ),
    "small_unbalanced_cells": (
        "A group of 5 or 7 samples per stage in an unbalanced design loses size, and balanced groups at 7 samples "
        "per stage hold. Whether balanced groups would protect 5 samples per stage was not read." + _DOCUMENTATION
        + _LIMIT
    ),
    "moderate_excess": (
        "The excess is not confined to small cells: n = 80 loses size at 9 or 13 samples per stage. Nothing is "
        "re-read: unbalance_excess stands, the roadmap keeps a diagnosis at n = 80 open and the dispersion families "
        "stay uninterpreted." + _LIMIT
    ),
    "not_replicated": (
        "The committed reading is not reproduced on fresh datasets. The documentation states that it did not "
        "reproduce, gives the pooled `delta` count of n80_r0.8 and n80_r0.8_b, and applies the re-read verdict's "
        "consequence; roadmap open item 16 moves to Done." + _LIMIT
    ),
}
DECISION_NAME = "unbalance_diagnosis_decision.json"

#: The re-read of the group-dispersion families under each verdict (design D7).
SET_ASIDE_VERDICTS: tuple[str, ...] = ("small_cells", "small_unbalanced_cells")
REPLACE_VERDICTS: tuple[str, ...] = ("not_replicated",)
SET_ASIDE_DESIGNS: tuple[str, ...] = ("n80_r0.8", "n80_r0.8_s1.5", "n80_r0.8_s2", "n80_r0.2_s1.5", "n80_r0.2_s2")
REREAD_M_SET_ASIDE: dict[str, int] = {
    "unbalanced_equal": 21, "balanced_unequal": 42, "wider_smaller": 54, "wider_larger": 54,
}
#: The fields a committed observation carries before ``judge_rates`` adds its own.
OBSERVATION_FIELDS: tuple[str, ...] = (
    "profile", "design", "column", "share", "scale", "cell", "arm", "n_rejected", "n_available", "status",
)
REREAD_SCOPES: dict[str, str] = {
    "set_aside": (
        "Readings at 13 samples per stage or more in every cell: the five group-dispersion designs at n = 80 with "
        "shares 0.8 and 0.2 (a 5-per-stage group) are set aside."
    ),
    "replaced": (
        "Every group-dispersion design, with the three readings of n80_r0.8 replaced by those of its fresh "
        "replicate n80_r0.8_b."
    ),
}


# ── Designs ──────────────────────────────────────────────────────────────────


def design_entries() -> list[dict[str, Any]]:
    """The rerun design, then the seven new designs in their declared order."""

    rerun = {"design": RERUN_DESIGN, "role": "rerun", "n_samples": 80, "share": 0.8,
             "larger_per_stage": RERUN_CELLS[0], "smaller_per_stage": RERUN_CELLS[1], "family": None,
             "replicates": RERUN_REPLICATES}
    new = [{"design": name, "role": "new", "n_samples": n, "share": share, "larger_per_stage": larger,
            "smaller_per_stage": smaller, "family": family, "replicates": NEW_REPLICATES}
           for name, (n, share, larger, smaller, family) in NEW_DESIGNS.items()]
    return [rerun, *new]


def design_generator(committed: Mapping[str, Any], n_samples: int, share: float) -> dict[str, Any]:
    """The committed ``n80_r0.8`` generator with only ``n_samples`` and ``group_ratio`` changed."""

    generator = dict(committed)
    generator["n_samples"] = int(n_samples)
    generator["group_ratio"] = float(share)
    return generator


# ── Parents ──────────────────────────────────────────────────────────────────


def parent_binding(parent_dir: Path | None = None) -> tuple[dict[str, Any], dict[str, Any]]:
    """The bound parent files and the committed ``columns`` settings (hash-checked)."""

    directory = Path(parent_dir) if parent_dir is not None else PARENT_DIR
    stored = ssc.load_parent(directory / PARENT_PROFILE)
    files = ubc._bind(directory / name for name in PARENT_FILE_NAMES)
    binding = {
        "dir": ssc._display_path(directory),
        "files": files,
        "columns_settings_hash": stored["hash"],
        "source": PARENT_SOURCE,
    }
    return binding, stored["settings"]


def committed_generator(committed: Mapping[str, Any], directory: Path) -> dict[str, Any]:
    """The committed generator of ``n80_r0.8``, checked to be an equal-dispersion three-block design."""

    designs = committed["unit_settings"]["units"]["designs"]
    if RERUN_DESIGN not in designs:
        raise SystemExit(f"{directory}/settings.json has no design {RERUN_DESIGN!r}.")
    generator = dict(designs[RERUN_DESIGN]["generator"])
    if (designs[RERUN_DESIGN].get("layers") is not None or generator.get("group_dispersion_scale") is not None
            or float(generator["group_ratio"]) != 0.8 or int(generator["n_samples"]) != 80):
        raise SystemExit(f"{directory}/settings.json: design {RERUN_DESIGN!r} is not the equal-dispersion "
                         "three-block design at n = 80 and share 0.8.")
    return generator


def _check_parent_units(committed: Mapping[str, Any], directory: Path) -> None:
    """The committed units are the ones ``group_dispersion_check.unit_settings`` writes today, at their grade.

    Rebuilding the committed unit settings from its own designs and cells must give them back exactly, so the
    units this driver builds with the same function differ from the parent's only in their designs and grade.
    """

    units = committed["unit_settings"]
    replicates = {int(value) for value in units["replicates"].values()}
    if len(replicates) != 1:
        raise SystemExit(f"{directory}/settings.json: the parent's units run unequal replicates {sorted(replicates)}.")
    generators = {name: entry["generator"] for name, entry in units["units"]["designs"].items()}
    rebuilt = gd.unit_settings(generators, units["units"]["cells"], replicates=replicates.pop(),
                               permutations=int(units["permutations"]))
    rebuilt["designs"] = list(units["designs"])
    if ubc._canonical(rebuilt) != ubc._canonical(units):
        differing = lsc.settings_difference(ubc._canonical(units), ubc._canonical(rebuilt))
        shown = ", ".join(differing[:8]) + (f" and {len(differing) - 8} more" if len(differing) > 8 else "")
        raise SystemExit(f"{directory}/settings.json differs from the units group_dispersion_check writes in "
                         f"{shown}; its designs would not be read as the parent read them.")
    if [list(cell) for cell in units["units"]["cells"]] != [list(cell) for cell in CELLS]:
        raise SystemExit(f"{directory}/settings.json: the parent's cells are not {[cell[0] for cell in CELLS]}.")


# ── Settings ─────────────────────────────────────────────────────────────────


def _select(entries: Sequence[dict[str, Any]], designs: Sequence[str] | None) -> list[dict[str, Any]]:
    if designs is None:
        return list(entries)
    known = [entry["design"] for entry in entries]
    unknown = sorted(set(designs) - set(known))
    if unknown:
        raise SystemExit(f"Unknown design(s) {unknown}; expected a subset of {known}.")
    return [entry for entry in entries if entry["design"] in designs]


def build_settings(
    *,
    replicates: int | None = None,
    permutations: int = DEFAULT_PERMUTATIONS,
    designs: Sequence[str] | None = None,
    parent_dir: Path | None = None,
) -> dict[str, Any]:
    """The run's settings: the rerun and the new designs on the committed generator, the parent bound.

    ``replicates`` (a smoke) overrides every design's count; by default the rerun runs 500 and each new design
    1,000. The unit settings are ``group_dispersion_check.unit_settings`` on these designs and the anchor cell, so
    ``evaluate_unit`` treats every unit as the parent's; their ``replicates`` entry is the new designs' count and
    each design's own count is in ``replicates`` here.
    """

    directory = Path(parent_dir) if parent_dir is not None else PARENT_DIR
    binding, committed = parent_binding(directory)
    _check_parent_units(committed, directory / PARENT_PROFILE)
    base = committed_generator(committed, directory / PARENT_PROFILE)
    entries = _select(design_entries(), designs)
    if replicates is not None:
        entries = [{**entry, "replicates": int(replicates)} for entry in entries]
    generators = {entry["design"]: base if entry["design"] == RERUN_DESIGN
                  else design_generator(base, entry["n_samples"], entry["share"]) for entry in entries}
    unit_replicates = NEW_REPLICATES if replicates is None else int(replicates)
    units = gd.unit_settings(generators, CELLS, replicates=unit_replicates, permutations=permutations)
    return {
        "driver": DRIVER,
        "driver_version": DRIVER_VERSION,
        "profile": PROFILE,
        "parent": binding,
        "parent_files": dict(binding["files"]),
        "rerun_design": RERUN_DESIGN if RERUN_DESIGN in generators else None,
        "replicates": {entry["design"]: int(entry["replicates"]) for entry in entries},
        "design_index": {entry["design"]: {key: entry[key] for key in entry if key not in ("design", "replicates")}
                         for entry in entries},
        "judged_statistics": list(STATISTICS),
        "families": {family: [entry["design"] for entry in entries if entry["family"] == family]
                     for family in FAMILIES},
        "min_available": MIN_AVAILABLE,
        "unit": "group_dispersion_check.evaluate_unit on unit_settings (the parent's stage_span settings with the "
                "designs and cells replaced), pca_matched the one recorded space, then the within-group arm",
        "corrected": {
            "corrections": list(CORRECTIONS),
            "rule": CORRECTED_RULE,
            "definition": "the dataset regenerated, production refitted at the recorded rank with no "
                          "cross-validation (forced_components), delta, angle and shape under "
                          "motco.stats.permutation.RRPP(..., residual_correction=<correction>) with unrestricted "
                          "draws at the unit's RRPP seed, p-values by the harness's empirical rule; recorded, never "
                          "judged",
        },
        "unit_settings": units,
    }


def check_bound_files(settings: Mapping[str, Any]) -> None:
    """Refuse, by name, a changed parent file."""

    ubc.check_parent_files(settings)


# ── The corrected arms ───────────────────────────────────────────────────────


def corrected_reading(
    dataset: Any, matrix: Any, params: Any, *, residual_correction: str | None, return_nulls: bool = False
) -> dict[str, Any]:
    """``delta``, ``angle`` and ``shape`` and their p-values under unrestricted draws with ``residual_correction``.

    The harness's statistics and p-values (``evaluate_semisynthetic_trajectory``) with only ``residual_correction``
    added to the RRPP call: ``None`` reproduces its p-values exactly. ``scale`` records the correction's
    multipliers (one value for ``df``, the range for ``leverage``).
    """

    from motco.simulations.evaluation import (
        _empirical_p_value,
        _extract_null_distributions,
        _extract_pair_statistics,
        build_simulation_trajectory_design,
    )
    from motco.stats.permutation import RRPP, residual_correction_scale
    from motco.stats.trajectory import estimate_difference

    design = build_simulation_trajectory_design(dataset.metadata, group_col=params.group_col,
                                                stage_col=params.stage_col)
    estimator = params.resolved_orientation_estimator
    deltas, angles, shapes, _ = estimate_difference(
        Y=matrix, model_matrix=design.model_full, LS_means=design.ls_means, contrast=design.contrast,
        return_spectra=True, orientation_estimator=estimator,
    )
    shape_available = len(design.stage_levels) >= 3
    observed = _extract_pair_statistics(deltas, angles, shapes, shape_available=shape_available)
    dist_delta, dist_angle, dist_shape, _ = RRPP(
        matrix, design.model_full, design.model_reduced, design.ls_means, design.contrast,
        permutations=int(params.permutations), n_jobs=params.n_jobs, progress=params.progress, seed=params.seed,
        return_eigengaps=True, orientation_estimator=estimator, residual_correction=residual_correction,
    )
    nulls = _extract_null_distributions(dist_delta, dist_angle, dist_shape, shape_available=shape_available)
    p_values = {statistic: _empirical_p_value(nulls[statistic], value) for statistic, value in observed.items()
                if statistic in nulls and np.isfinite(value)}
    out: dict[str, Any] = {"pair_statistics": lsc._finite(observed), "p_values": lsc._finite(p_values)}
    if residual_correction is not None:
        scale = residual_correction_scale(design.model_reduced, residual_correction)
        out["scale"] = {"min": float(scale.min()), "max": float(scale.max())}
    if return_nulls:
        out["nulls"] = {statistic: [float(value) for value in values] for statistic, values in nulls.items()}
    return out


def corrected_arms(unit: lsc.Unit, units: Mapping[str, Any], record: Mapping[str, Any]) -> dict[str, Any]:
    """The corrected arms of one completed unit, on one refit; a failure is recorded, never raised."""

    from motco.simulations.evaluation import integrate_semisynthetic_dataset

    try:
        fields = ubc.unit_fields(unit, units)
        dataset, params = ubc.unit_inputs(unit, units, fields)
        rank = int(record["production"]["rank"])
        latent = integrate_semisynthetic_dataset(dataset, gd.forced_rank_params(params, rank))
    except Exception as exc:  # recorded for these arms only
        return {"status": "failed", "error_type": type(exc).__name__, "error_message": str(exc)}
    out: dict[str, Any] = {"status": "completed", "rank": rank}
    statistics: Mapping[str, Any] | None = None
    for correction in CORRECTIONS:
        try:
            reading = corrected_reading(dataset, latent.matrix, params, residual_correction=correction)
        except Exception as exc:  # recorded for this correction only
            out[correction] = {"status": "failed", "error_type": type(exc).__name__, "error_message": str(exc)}
            continue
        if statistics is None:
            statistics = reading["pair_statistics"]
        out[correction] = {
            "status": "completed",
            "scale": reading["scale"],
            "p_values": reading["p_values"],
            "rejects": {statistic: ubc.rejects(reading["p_values"].get(statistic), CORRECTED_RULE)
                        for statistic in STATISTICS},
        }
    # With both corrections failed there is nothing to compare: neither consistent nor inconsistent.
    out["pair_statistics"] = statistics
    out["same_statistics"] = (None if statistics is None
                              else gd.same_statistics(statistics, record["production"]["pair_statistics"]))
    return out


# ── Units and the run ────────────────────────────────────────────────────────


def run_units(settings: Mapping[str, Any]) -> list[lsc.Unit]:
    """Design-major, then cell, then replicate (the parent's order), at each design's own count."""

    units = settings["unit_settings"]
    return [lsc.Unit(design, cell, replicate)
            for design in units["designs"]
            for cell in ubc._cell_names(units)
            for replicate in range(int(settings["replicates"][design]))]


def evaluate_unit(unit: lsc.Unit, units: Mapping[str, Any]) -> dict[str, Any]:
    """The group-dispersion run's unit, unchanged, then the corrected arms on the same fit, statistics and seed."""

    start = perf_counter()
    record = dict(gd.evaluate_unit(unit, units))
    mark = perf_counter()
    if record.get("status") == "completed":
        record["corrected"] = corrected_arms(unit, units, record)
    else:
        record["corrected"] = {"status": "skipped", "reason": "the parent's unit failed"}
    runtime = dict(record.get("runtime_seconds") or {})
    runtime["group_dispersion_unit"] = runtime.get("unit")
    runtime["corrected"] = perf_counter() - mark
    runtime["unit"] = perf_counter() - start
    record["runtime_seconds"] = runtime
    return lsc._jsonable(record)


def run_shard(
    out_dir: Path,
    settings: Mapping[str, Any],
    index: int,
    count: int,
    *,
    limit: int | None = None,
    workers: int = 1,
) -> int:
    """Evaluate this shard's missing units into its file; returns how many were evaluated."""

    check_bound_files(settings)
    directory = Path(out_dir)
    lsc.ensure_settings(directory, settings)
    units_settings = settings["unit_settings"]
    path = lsc.shard_path(directory, index, count)
    path.parent.mkdir(parents=True, exist_ok=True)
    done = ssc.read_records([path], ubc.record_key) if path.exists() else {}
    # Rewritten with complete lines only, so an interrupted unit is re-run, not duplicated.
    with path.open("w", encoding="utf-8") as handle:
        for key in sorted(done):
            handle.write(lsc._line(done[key]))
    units = run_units(settings)
    pending = [unit for position, unit in enumerate(units)
               if position % count == index and ubc.unit_key(unit) not in done]
    if limit is not None:
        pending = pending[:limit]

    def append(position: int, record: Mapping[str, Any]) -> None:
        with path.open("a", encoding="utf-8") as handle:
            handle.write(lsc._line(record))
            handle.flush()
        print(f"[{index}/{count}] {position}/{len(pending)} {ubc.record_key(record)} {record['status']}",
              file=sys.stderr, flush=True)

    if workers <= 1:
        for position, unit in enumerate(pending, start=1):
            append(position, evaluate_unit(unit, units_settings))
    else:
        with ProcessPoolExecutor(max_workers=workers) as pool:
            futures = [pool.submit(evaluate_unit, unit, units_settings) for unit in pending]
            for position, future in enumerate(as_completed(futures), start=1):
                append(position, future.result())
    return len(pending)


# ── Merge ────────────────────────────────────────────────────────────────────


def _corrected(record: Mapping[str, Any]) -> Mapping[str, Any]:
    return record.get("corrected") or {}


def corrected_rows(
    records: Sequence[Mapping[str, Any]], units: Mapping[str, Any]
) -> tuple[list[dict[str, Any]], list[dict[str, Any]]]:
    """The corrected arms' rate rows and their per-design accounting.

    A rate counts the units whose refit completed and is consistent and whose correction completed; ``failures``
    counts every other unit of the design (a failed parent unit, a failed refit or correction, an inconsistent
    refit).
    """

    rates, table = [], []
    for fields, group in ubc._cell_groups(records, units):
        completed = [r for r in group if r.get("status") == "completed"]
        refit = [r for r in completed if _corrected(r).get("status") == "completed"]
        consistent = [r for r in refit if _corrected(r).get("same_statistics") is True]
        for correction in CORRECTIONS:
            read = [r for r in consistent if (_corrected(r).get(correction) or {}).get("status") == "completed"]
            for statistic in STATISTICS:
                rates.append({
                    **fields, "arm": f"{correction}_{statistic}", "rule": CORRECTED_RULE,
                    **ssc._rate_row([_corrected(r)[correction]["rejects"].get(statistic) for r in read]),
                    "failures": len(group) - len(read),
                })
        table.append({
            **fields, "units": len(group), "parent_completed": len(completed), "refit_completed": len(refit),
            "refit_failed": len(completed) - len(refit),
            "inconsistent": sum(1 for r in refit if _corrected(r).get("same_statistics") is False),
            **{f"{c}_failed": sum(1 for r in refit if (_corrected(r).get(c) or {}).get("status") != "completed")
               for c in CORRECTIONS},
        })
    return rates, table


def rate_rows(records: Sequence[Mapping[str, Any]], units: Mapping[str, Any]) -> list[dict[str, Any]]:
    """The group-dispersion run's rate rows with the corrected rows after each cell's within-group rows."""

    parent = gd.rate_rows(records, units)
    corrected, _ = corrected_rows(records, units)
    designs = {design: position for position, design in enumerate(units["designs"])}
    cells = {cell: position for position, cell in enumerate(ubc._cell_names(units))}
    arms = {arm: position for position, arm in enumerate(RATE_ARMS)}
    return sorted([*parent, *corrected],
                  key=lambda row: (designs[str(row["design"])], cells[str(row["cell"])], arms[str(row["arm"])]))


def corrected_accounting(table: Sequence[Mapping[str, Any]]) -> dict[str, Any]:
    def by_design(key: str) -> dict[str, int]:
        return {f"{row['design']}/{row['cell']}": int(row[key]) for row in table if int(row[key])}

    return {
        "completed": sum(int(row["refit_completed"]) for row in table),
        "refit_failed": sum(int(row["refit_failed"]) for row in table),
        "inconsistent": sum(int(row["inconsistent"]) for row in table),
        "refit_failed_by_design": by_design("refit_failed"),
        "inconsistent_by_design": by_design("inconsistent"),
        **{f"{c}_failed_by_design": by_design(f"{c}_failed") for c in CORRECTIONS},
    }


def _parent_table(settings: Mapping[str, Any], name: str) -> list[dict[str, Any]]:
    return lsc._read_csv(ssc._resolve(settings["parent"]["dir"]) / PARENT_PROFILE / name)


def rerun_reproduction(
    rates: Sequence[Mapping[str, Any]],
    within: Sequence[Mapping[str, Any]],
    committed_rates: Sequence[Mapping[str, Any]],
    committed_within: Sequence[Mapping[str, Any]],
    settings: Mapping[str, Any],
) -> dict[str, Any]:
    """The rerun design against the committed ``columns`` tables, exact on every count, every arm.

    Each rate row is compared on its available and rejected counts for every arm the committed table holds
    (the parent arms and the within-group arms); the within-group table's row on every count column. The corrected
    arms are new and have no committed row.
    """

    base = {"design": RERUN_DESIGN, "source": settings["parent"]["source"], "files": dict(settings["parent_files"])}
    if settings.get("rerun_design") != RERUN_DESIGN:
        return {**base, "passed": False, "rows_compared": 0, "differing": [], "unresolved_rows": [],
                "reason": f"the rerun design {RERUN_DESIGN!r} is not in this run"}
    compared_arms = set(gd.RATE_ARMS)
    mine: dict[tuple[str, str, str], Any] = {
        ("rate_table", str(r["cell"]), str(r["arm"])): ubc._counts(r) for r in rates
        if str(r["design"]) == RERUN_DESIGN and str(r["arm"]) in compared_arms}
    theirs: dict[tuple[str, str, str], Any] = {
        ("rate_table", str(r["cell"]), str(r["arm"])): ubc._counts(r) for r in committed_rates
        if str(r["design"]) == RERUN_DESIGN}

    def within_counts(row: Mapping[str, Any]) -> dict[str, int]:
        return {column: int(row[column]) for column in WITHIN_COUNT_COLUMNS}

    mine.update({("within_group_table", str(r["cell"]), "within_group"): within_counts(r) for r in within
                 if str(r["design"]) == RERUN_DESIGN})
    theirs.update({("within_group_table", str(r["cell"]), "within_group"): within_counts(r) for r in committed_within
                   if str(r["design"]) == RERUN_DESIGN})
    differing, unresolved = [], []
    compared = 0
    for key in [*theirs, *(key for key in mine if key not in theirs)]:
        table, cell, arm = key
        if key not in mine or key not in theirs:
            unresolved.append({"table": table, "cell": cell, "arm": arm,
                               "missing_from": "this run" if key not in mine else "the committed table"})
            continue
        compared += 1
        wanted, got = theirs[key], mine[key]
        if got != wanted:
            if table == "rate_table":
                assert isinstance(wanted, tuple) and isinstance(got, tuple)
                differing.append({"table": table, "cell": cell, "arm": arm,
                                  "committed": {"count": wanted[0], "rejections": wanted[1]},
                                  "reproduced": {"count": got[0], "rejections": got[1]}})
            else:
                differing.append({"table": table, "cell": cell, "arm": arm, "committed": wanted, "reproduced": got})
    return {**base, "passed": compared > 0 and not differing and not unresolved, "rows_compared": compared,
            "differing": differing, "unresolved_rows": unresolved}


def merge(out_dir: Path) -> dict[str, Any]:
    """Merge the shards and write the tables, the accounting and the reproduction check."""

    directory = Path(out_dir)
    settings = lsc.load_settings(directory)
    check_bound_files(settings)
    units = settings["unit_settings"]
    unit_hash = lsc.settings_hash(units)
    expected_keys = [ubc.unit_key(unit) for unit in run_units(settings)]
    raw: list[dict[str, Any]] = []
    for path in sorted((directory / "records").glob("shard_*.jsonl")):
        raw.extend(ssc.read_records([path], ubc.record_key).values())
    stale = sorted({ubc.record_key(r) for r in raw if r.get("settings_hash") != unit_hash})
    if stale:
        raise SystemExit(f"{len(stale)} record(s) were written under other settings, e.g. {stale[:3]}; refusing.")
    seen: dict[str, int] = {}
    for record in raw:
        seen[ubc.record_key(record)] = seen.get(ubc.record_key(record), 0) + 1
    by_key = {ubc.record_key(r): r for r in raw}
    records = [by_key[key] for key in expected_keys if key in by_key]
    with (directory / "merged.jsonl").open("w", encoding="utf-8") as handle:
        handle.write("".join(lsc._line(record) for record in records))
    _, within_table = gd.within_group_rows(records, units)
    _, corrected_table = corrected_rows(records, units)
    summary = {
        **ubc.accounting(records, units, expected=len(expected_keys),
                         duplicates=sum(count - 1 for count in seen.values())),
        "profile": PROFILE,
        "settings_hash": lsc.settings_hash(settings),
        "unit_settings_hash": unit_hash,
        "expected_by_design": dict(settings["replicates"]),
        "within_group": gd.within_accounting(within_table),
        "corrected": corrected_accounting(corrected_table),
    }
    rates = rate_rows(records, units)
    lsc._write_csv(directory / "rate_table.csv", rates, ubc.RATE_COLUMNS)
    lsc._write_csv(directory / "cell_table.csv", ubc.cell_rows(records, units), ubc.CELL_COLUMNS)
    lsc._write_csv(directory / "within_group_table.csv", within_table, gd.WITHIN_COLUMNS)
    lsc._write_csv(directory / "corrected_table.csv", corrected_table, CORRECTED_COLUMNS)
    lsc._write_json(directory / "accounting.json", summary)
    result = rerun_reproduction(rates, within_table, _parent_table(settings, "rate_table.csv"),
                                _parent_table(settings, "within_group_table.csv"), settings)
    lsc._write_json(directory / "reproduction.json", result)
    return {
        **{key: summary[key] for key in ("complete_units", "missing_units", "failed_units", "recorded_core_hours")},
        "within_group_inconsistent": summary["within_group"]["inconsistent"],
        "corrected_inconsistent": summary["corrected"]["inconsistent"],
        "reproduction_passed": result["passed"],
        "reproduction_differing": len(result["differing"]),
    }


# ── Decision (pure functions) ────────────────────────────────────────────────


def size_observation(
    row: Mapping[str, Any] | None, *, min_available: int = MIN_AVAILABLE, **cell: Any
) -> dict[str, Any]:
    """One judged reading; fewer than ``min_available`` available units makes it unavailable."""

    observation = isc._size_observation(row, **cell)
    if observation.get("status") != "unavailable" and int(observation["n_available"]) < min_available:
        return {**cell, "status": "unavailable", "reason": f"fewer than {min_available} available units",
                "n_available_read": observation["n_available"], "n_rejected_read": observation["n_rejected"]}
    return observation


def _observation(block: Mapping[str, Any], rates: Mapping[tuple, Mapping[str, Any]], design: str, arm: str,
                 *, profile: str = PROFILE) -> dict[str, Any]:
    meta = block["settings"]["design_index"][design]
    return size_observation(
        rates.get((design, CELL, arm)), min_available=int(block["settings"]["min_available"]),
        profile=profile, design=design, n_samples=meta["n_samples"], share=meta["share"],
        larger_per_stage=meta["larger_per_stage"], smaller_per_stage=meta["smaller_per_stage"], cell=CELL, arm=arm,
    )


def family_readings(block: Mapping[str, Any], arms: Sequence[str] = STATISTICS) -> dict[str, list[dict[str, Any]]]:
    """Every family's observations: each member design's anchor, one per arm; the rerun design enters none."""

    rates = ubc._index(block["rates"], "design", "cell", "arm")
    members = block["settings"]["families"]
    return {family: [_observation(block, rates, design, arm) for design in members.get(family, []) for arm in arms]
            for family in FAMILIES}


def size_reading(block: Mapping[str, Any]) -> dict[str, Any]:
    """``delta``, ``angle`` and ``shape`` judged in the three families."""

    readings = family_readings(block)
    return {family: {**grade.judge_rates(readings[family]), "m_expected": FAMILY_M_EXPECTED[family]}
            for family in FAMILIES}


def verdict_of(reproduction: Mapping[str, Any], holds: Mapping[str, bool]) -> str:
    """The first verdict that applies, in the fixed order."""

    if not reproduction.get("passed"):
        return "reproduction_failed"
    if not holds["balanced_small"]:
        return "small_cells"
    if not holds["unbalanced_small"]:
        return "small_unbalanced_cells"
    if not holds["reference"]:
        return "moderate_excess"
    return "not_replicated"


def pooled_delta(block: Mapping[str, Any]) -> dict[str, Any]:
    """The pooled ``delta`` count of the committed design's rerun and its fresh replicate; never judged."""

    rates = ubc._index(block["rates"], "design", "cell", "arm")
    parts = {design: ubc._counts(rates.get((design, CELL, "delta"))) for design in (RERUN_DESIGN, FRESH_REPLICATE)}
    available = sum(counts[0] for counts in parts.values() if counts is not None)
    rejected = sum(counts[1] for counts in parts.values() if counts is not None)
    rate = rejected / available if available else None
    return {
        "designs": list(parts),
        "by_design": {design: None if counts is None else {"n_available": counts[0], "n_rejected": counts[1]}
                      for design, counts in parts.items()},
        "n_available": available,
        "n_rejected": rejected,
        "rate": rate,
        "se": None if rate is None else (rate * (1.0 - rate) / available) ** 0.5,
    }


# ── The re-read of the group-dispersion families ────────────────────────────


def committed_observations(parent: Mapping[str, Any]) -> dict[str, list[dict[str, Any]]]:
    """The committed group-dispersion families' observations, stripped of what ``judge_rates`` added."""

    size = parent["criteria"]["size"]
    return {family: [{key: row[key] for key in OBSERVATION_FIELDS if key in row and
                      (key != "status" or row[key] == "unavailable")} for row in size[family]["observations"]]
            for family in gd.FAMILIES}


def group_dispersion_verdict(reproduction_passed: bool, holds: Mapping[str, bool]) -> str:
    """The group-dispersion run's verdict order (``group_dispersion_check.decision``) on judged families."""

    if not reproduction_passed:
        return "reproduction_failed"
    if not holds["unbalanced_equal"]:
        return "unbalance_excess"
    if not holds["balanced_unequal"]:
        return "dispersion_excess"
    if not (holds["wider_smaller"] and holds["wider_larger"]):
        return "unbalanced_dispersion_excess"
    return "size_holds"


def reread_families(
    parent: Mapping[str, Any],
    *,
    set_aside: Sequence[str] = (),
    replacement: Sequence[Mapping[str, Any]] | None = None,
) -> dict[str, Any]:
    """The group-dispersion families judged again on the committed counts.

    ``set_aside`` removes every reading of the named designs; ``replacement`` takes the place of the three
    ``n80_r0.8`` readings of ``unbalanced_equal`` (at their position). With neither, the committed families and
    verdict come back.
    """

    observations = committed_observations(parent)
    out: dict[str, list[dict[str, Any]]] = {}
    for family in gd.FAMILIES:
        rows = [row for row in observations[family] if row["design"] not in set_aside]
        if replacement is not None and family == "unbalanced_equal":
            position = next(index for index, row in enumerate(rows) if row["design"] == RERUN_DESIGN)
            rows = [row for row in rows if row["design"] != RERUN_DESIGN]
            rows[position:position] = [dict(row) for row in replacement]
        out[family] = rows
    size = {family: grade.judge_rates(out[family]) for family in gd.FAMILIES}
    holds = {family: size[family]["holds"] for family in gd.FAMILIES}
    passed = bool((parent.get("reproduction") or {}).get("passed"))
    verdict = group_dispersion_verdict(passed, holds)
    return {"size": size, "holds": holds, "verdict": verdict}


def _family_summary(judged: Mapping[str, Any], m_expected: int) -> dict[str, Any]:
    return {**{key: judged[key] for key in (
        "holds", "m", "median_n_available", "allowed_exceedances", "observed_exceedances", "ceiling_at_median",
        "n_above_ceiling", "n_unavailable")}, "m_expected": m_expected}


def reread(parent: Mapping[str, Any], verdict: str, block: Mapping[str, Any]) -> dict[str, Any]:
    """The re-read this run's verdict calls for (design D7); the group-dispersion files are only read."""

    committed = {"verdict": parent["verdict"], "families": {family: {key: entry[key] for key in (
        "holds", "m", "observed_exceedances", "n_above_ceiling")} for family, entry in parent["families"].items()}}
    if verdict in SET_ASIDE_VERDICTS:
        branch, scope = "set_aside", REREAD_SCOPES["set_aside"]
        judged = reread_families(parent, set_aside=SET_ASIDE_DESIGNS)
        m_expected = REREAD_M_SET_ASIDE
        changed: dict[str, Any] = {"set_aside": list(SET_ASIDE_DESIGNS), "replaced": None}
    elif verdict in REPLACE_VERDICTS:
        branch, scope = "replaced", REREAD_SCOPES["replaced"]
        rates = ubc._index(block["rates"], "design", "cell", "arm")
        replacement = [_observation(block, rates, FRESH_REPLICATE, arm) for arm in STATISTICS]
        judged = reread_families(parent, replacement=replacement)
        m_expected = dict(gd.FAMILY_M_EXPECTED)
        changed = {"set_aside": [], "replaced": {"design": RERUN_DESIGN, "by": FRESH_REPLICATE,
                                                 "family": "unbalanced_equal", "readings": replacement}}
    else:
        return {
            "performed": False,
            "branch": None,
            "reason": f"under {verdict} nothing is re-read; the group-dispersion verdict {parent['verdict']} stands",
            "committed": committed,
            "verdict": parent["verdict"],
            "consequence": None,
            "scope": None,
            "set_aside": [],
            "replaced": None,
        }
    size = judged["size"]
    return {
        "performed": True,
        "branch": branch,
        "rule": "stage_span_paper_grade_decision.judge_rates (count plus Bonferroni ceiling) per family, the "
                "group-dispersion verdict order",
        "committed": committed,
        **changed,
        "families": {family: _family_summary(size[family], m_expected[family]) for family in gd.FAMILIES},
        "failing": {family: ssc._failing_observations(size[family]) for family in gd.FAMILIES
                    if not size[family]["holds"]},
        "criteria": {"size": size},
        "verdict": judged["verdict"],
        "consequence": gd.CONSEQUENCES[judged["verdict"]],
        "scope": scope,
        "note": "The re-read rule was fixed before this run's units and after the group-dispersion families' rates "
                "were known.",
    }


# ── What is reported beside the verdict ──────────────────────────────────────


def recorded_arm_labels(block: Mapping[str, Any]) -> dict[str, Any]:
    """Each recorded arm's size label per family, under the same rule; never part of the verdict."""

    out: dict[str, Any] = {}
    for arm in RECORDED_ARMS:
        readings = family_readings(block, arms=(arm,))
        labels = {}
        for family in FAMILIES:
            judged = grade.judge_rates(readings[family])
            labels[family] = {
                "label": "holds" if judged["holds"] else "exceeds",
                "m": judged["m"],
                "allowed_exceedances": judged["allowed_exceedances"],
                "observed_exceedances": judged["observed_exceedances"],
                "n_above_ceiling": judged["n_above_ceiling"],
                "n_unavailable": judged["n_unavailable"],
                "above_ceiling": [{key: row.get(key) for key in ("design", "rate")}
                                  for row in judged["observations"] if row.get("status") == "above_ceiling"],
            }
        rule = ubc.RULES.get(arm, CORRECTED_RULE if arm in CORRECTED_ARMS else gd.WITHIN_RULE)
        out[arm] = {"rule": rule, "families": labels}
    return out


def design_counts(block: Mapping[str, Any]) -> list[dict[str, Any]]:
    """Every design's counts by statistic, the rerun design included."""

    rates = ubc._index(block["rates"], "design", "cell", "arm")
    out = []
    for design, meta in block["settings"]["design_index"].items():
        counts = {}
        for statistic in STATISTICS:
            row = rates.get((design, CELL, statistic))
            pair = ubc._counts(row)
            counts[statistic] = None if pair is None else {
                "n_available": pair[0], "n_rejected": pair[1], "rate": ubc._finite(row.get("rate") if row else None)}
        out.append({"design": design, **meta, "counts": counts})
    return _ordered(block, out)


def _ordered(block: Mapping[str, Any], rows: list[dict[str, Any]]) -> list[dict[str, Any]]:
    """Rows in the run's design order (stored settings sort their keys)."""

    order = {design: position for position, design in enumerate(block["settings"]["unit_settings"]["designs"])}
    return sorted(rows, key=lambda row: order[row["design"]])


def fit_by_design(block: Mapping[str, Any]) -> list[dict[str, Any]]:
    """The selected rank and held-out AUROC per design."""

    cells = ubc._index(block["cells"], "design", "cell")
    out = []
    for design, meta in block["settings"]["design_index"].items():
        row = cells.get((design, CELL)) or {}
        out.append({"design": design, **meta, **{key: ubc._finite(row.get(key)) for key in (
            "completed", "rank_q10", "rank_median", "rank_q90", "auroc_median")}})
    return _ordered(block, out)


def arm_rates(block: Mapping[str, Any]) -> list[dict[str, Any]]:
    """Every arm's rate per design."""

    rates = ubc._index(block["rates"], "design", "cell", "arm")
    out = [{"design": design, **meta, "rates": {arm: ubc._rate_of(rates.get((design, CELL, arm)))
                                                for arm in RATE_ARMS}}
           for design, meta in block["settings"]["design_index"].items()]
    return _ordered(block, out)


def decision(block: Mapping[str, Any], parent: Mapping[str, Any]) -> dict[str, Any]:
    """Apply rule version 1.

    ``block`` holds ``settings``, ``accounting``, ``rates``, ``cells`` (rows) and ``reproduction``; ``parent`` is
    the committed group-dispersion decision. Every family is recorded whatever the verdict.
    """

    size = size_reading(block)
    reproduction = block.get("reproduction") or {}
    holds = {family: size[family]["holds"] for family in FAMILIES}
    verdict = verdict_of(reproduction, holds)
    failing: dict[str, Any] = {}
    if verdict == "reproduction_failed":
        failing["reproduction"] = {key: reproduction.get(key) for key in ("differing", "unresolved_rows", "reason")
                                   if reproduction.get(key)}
    for family in FAMILIES:
        if not holds[family]:
            failing[family] = ssc._failing_observations(size[family])
    accounting = block["accounting"]
    return {
        "rule_version": RULE_VERSION,
        "alpha": ubc.ALPHA,
        "judged_statistics": list(STATISTICS),
        "rule": "stage_span_paper_grade_decision.judge_rates (count plus Bonferroni ceiling) per family; a reading "
                f"with fewer than {block['settings']['min_available']} available units is unavailable",
        "settings_hash": lsc.settings_hash(block["settings"]),
        "parent": block["settings"].get("parent"),
        "complete": accounting["missing_units"] == 0,
        "failed_units": accounting["failed_units"],
        "failures_by_cell": accounting.get("failures_by_cell") or {},
        "within_group": accounting.get("within_group") or {},
        "corrected": accounting.get("corrected") or {},
        "reproduction": {key: reproduction.get(key) for key in ("passed", "rows_compared")},
        "verdict": verdict,
        "consequence": CONSEQUENCES[verdict],
        "consequences": dict(CONSEQUENCES),
        "families": {family: {key: size[family][key] for key in (
            "holds", "m", "m_expected", "median_n_available", "allowed_exceedances", "observed_exceedances",
            "ceiling_at_median", "n_above_ceiling", "n_unavailable")} for family in FAMILIES},
        "failing": failing,
        "criteria": {"size": size},
        "pooled_delta": pooled_delta(block),
        "reread": reread(parent, verdict, block),
        "counts": design_counts(block),
        "recorded_arms": recorded_arm_labels(block),
        "fit": fit_by_design(block),
        "rates": arm_rates(block),
    }


def _tables(out_dir: Path) -> dict[str, Any]:
    directory = Path(out_dir)
    settings = lsc.load_settings(directory)
    check_bound_files(settings)
    path = directory / "accounting.json"
    if not path.exists():
        raise SystemExit(f"{path} does not exist; run `merge` first.")
    summary = json.loads(path.read_text(encoding="utf-8"))
    if summary.get("settings_hash") != lsc.settings_hash(settings):
        raise SystemExit(f"{path} was merged under other settings; run `merge` again.")
    return {
        "settings": settings,
        "accounting": summary,
        "rates": ubc._read_csv(directory / "rate_table.csv"),
        "cells": ubc._read_csv(directory / "cell_table.csv"),
        "reproduction": json.loads((directory / "reproduction.json").read_text(encoding="utf-8")),
    }


def decide(out_dir: Path) -> dict[str, Any]:
    """Read the merged aggregates and the committed group-dispersion decision; write the decision."""

    block = _tables(out_dir)
    parent_path = ssc._resolve(block["settings"]["parent"]["dir"]) / gd.DECISION_NAME
    parent = json.loads(parent_path.read_text(encoding="utf-8"))
    result = decision(block, parent)
    lsc._write_json(Path(out_dir) / DECISION_NAME, result)
    return result


# ── CLI ──────────────────────────────────────────────────────────────────────


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = parser.add_subparsers(dest="command", required=True)
    run = sub.add_parser("run", help="Evaluate (a shard of) the units")
    run.add_argument("--out-dir", type=Path, required=True)
    run.add_argument("--shard", type=str, default="0/1", help="i/N, round-robin over the units")
    run.add_argument("--limit", type=int, default=None, help="Evaluate at most this many pending units")
    run.add_argument("--workers", type=int, default=1)
    run.add_argument("--perms", type=int, default=DEFAULT_PERMUTATIONS)
    run.add_argument("--reps", type=int, default=None,
                     help=f"Replicates of every design (a smoke); default {RERUN_REPLICATES} for the rerun design "
                          f"and {NEW_REPLICATES} for each new design")
    run.add_argument("--designs", type=str, default=None, help="Comma-separated designs; default: all")
    merge_parser = sub.add_parser("merge", help="Merge the shards")
    merge_parser.add_argument("--out-dir", type=Path, required=True)
    decide_parser = sub.add_parser("decide", help=f"Write {DECISION_NAME}")
    decide_parser.add_argument("--out-dir", type=Path, required=True)
    return parser


def settings_from_args(args: argparse.Namespace) -> dict[str, Any]:
    chosen = None if args.designs is None else [name for name in args.designs.split(",") if name]
    return build_settings(replicates=args.reps, permutations=args.perms, designs=chosen)


def _summary(result: Mapping[str, Any]) -> dict[str, Any]:
    return {
        "verdict": result["verdict"],
        "reproduction": result["reproduction"],
        "families": {family: {key: entry[key] for key in ("holds", "m", "observed_exceedances", "n_above_ceiling",
                                                           "n_unavailable")}
                     for family, entry in result["families"].items()},
        "pooled_delta": {key: result["pooled_delta"][key] for key in ("n_rejected", "n_available", "rate")},
        "reread": {key: result["reread"][key] for key in ("performed", "branch", "verdict")},
        "complete": result["complete"],
    }


def main(argv: Iterable[str] | None = None) -> int:
    args = build_parser().parse_args(None if argv is None else list(argv))
    if args.command == "run":
        index, count = lsc.parse_shard(args.shard)
        evaluated = run_shard(args.out_dir, settings_from_args(args), index, count, limit=args.limit,
                              workers=args.workers)
        print(f"evaluated {evaluated} unit(s)", file=sys.stderr)
        return 0
    if args.command == "merge":
        print(json.dumps(merge(args.out_dir), indent=1))
        return 0
    print(json.dumps(_summary(decide(args.out_dir)), indent=1))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
