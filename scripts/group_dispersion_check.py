#!/usr/bin/env python3
"""Group-specific dispersion: whether the trajectory test holds its size when the groups' noise dispersions differ.

The driver behind ``openspec/changes/group-specific-dispersion/`` (spec
``group-dispersion-validation``, roadmap open item 8). RRPP exchanges the
reduced model's residuals between the two groups, so its null assumes one noise
dispersion for both. Every unit here is a null unit: the groups share their
trajectory and differ in the noise standard deviation of the second group
(``SemiSyntheticTrajectoryParams.group_dispersion_scale`` = s, 0.26.0) and in
their share of the samples (``group_ratio``, the first group's share). At a
share above 0.5 the wider group is the smaller one, below 0.5 the larger one.

**The unit** is the univariate baseline's (``univariate_baseline_check.evaluate_unit``,
imported unchanged), handed a parent-shaped ``stage_span`` settings block
(``unit_settings``: the parent's settings with only the designs and cells
replaced): one dataset, one cross-validated production fit (stage PLS on
methylation, expression and proteomics, 658 features, rank floor 2, net
displacement), ``delta``, ``angle`` and ``shape`` with their null distributions,
then the split test, the interaction scan, the false-discovery reading and
``pca_matched`` at the unit's RRPP seed. Its record equals the parent driver's
for the same generator parameters, cell and seeds.

**The within-group arm** (recorded, never judged). After the parent's unit, the
dataset is regenerated, production is refitted at the recorded rank with no
cross-validation (``forced_components``), and ``delta``, ``angle`` and ``shape``
are read under ``RRPP(..., strata=<group positions>)`` at the unit's RRPP seed
with the harness's extraction and p-value helpers: each draw permutes the
reduced residuals within one group, so it keeps each group's own dispersion.
``same_statistics`` checks that the refit's observed statistics equal the
production record's exactly; a unit where they differ is inconsistent for this
arm only and leaves its rates.

**Profiles** (base seed 2043, the parent's; three stages, three blocks; every
cell at 500 replicates and 999 permutations):

- ``core`` — n = 200, δ = 2.0 (the univariate baseline's ``n200``). The rerun
  design ``n200`` (the parent's name and seeds: the reproduction gate) and 17
  new designs: equal dispersion at shares 0.65 and 0.8, and s ∈ {1.25, 1.5, 2}
  at shares 0.5, 0.65, 0.8, 0.35 and 0.2. Cells: the anchor and the translation
  at e = 1.0, both under the planted scale. 18,000 units.
- ``columns`` — the univariate baseline's ``n80`` and ``n500`` (δ = 2.0) and
  ``n200_delta0.5``, and the stage-effect scaling run's ``band_n320`` (δ = 0.05
  at n = 320, read from its committed settings). At each: equal dispersion at
  share 0.8, and s ∈ {1.5, 2} at shares 0.5, 0.8 and 0.2; the anchor only.
  28 designs, 14,000 units. Each column's balanced equal-dispersion reading is
  the committed one.

New designs are named ``<column>_r<share>`` (equal dispersion) and
``<column>_r<share>_s<scale>``, so their dataset seeds
(``latent_space_comparison.generator_seed``) are new at base seed 2043, and the
two cells of one design and replicate share their null dataset.

**Parents.** ``settings.json``, ``rate_table.csv`` and ``cell_table.csv`` of the
univariate baseline's ``stage_span`` profile
(``results/univariate-baseline-2026-11-23/stage_span/``) and of the
stage-effect scaling run's ``band`` profile
(``results/stage-effect-scaling-2026-11-24/band/``) are bound by sha256 in both
profiles' settings; a changed file is refused by name in ``run``, ``merge`` and
``decide``.

Subcommands (outputs under ``--out-dir``):

- ``run --profile core|columns`` — (a shard of) the units into
  ``<profile>/records/shard_<i>_of_<N>.jsonl``; ``<profile>/settings.json`` holds
  the settings and their hash, a resume or merge under other settings exits
  naming the differing keys, and a resume skips complete units.
- ``merge --profile …`` — the parent's ``rate_table.csv`` (every parent arm, plus
  the rows ``within_group_<statistic>`` under ``p < alpha``) and
  ``cell_table.csv``, ``within_group_table.csv`` (per design and cell: units,
  parent completions, within-group completions, failures and inconsistent
  units), ``accounting.json`` and, for ``core``, ``reproduction.json`` (the rerun
  design against the committed rate table, exact on available and rejected
  counts for both cells and every parent arm, every differing row with both
  counts).
- ``decide`` — ``group_dispersion_decision.json`` (rule version 1): the gate
  ``reproduction_failed``; size of ``delta``, ``angle`` and ``shape`` at p < 0.05
  on every control cell of every new design, in four families over both
  profiles under the stage-span paper grade's count-plus-ceiling rule
  (``judge_rates``, imported unchanged): ``unbalanced_equal`` (m = 24),
  ``balanced_unequal`` (m = 42), ``wider_smaller`` (m = 60) and ``wider_larger``
  (m = 60); the verdict, the first that applies: ``reproduction_failed`` |
  ``unbalance_excess`` | ``dispersion_excess`` | ``unbalanced_dispersion_excess``
  | ``size_holds``, with its consequence. Beside it, never changing it: the map
  (per statistic, column and share, the smallest scale whose reading is above
  its family's ceiling, or none up to the largest scale read), every
  conservative reading (below 0.05 − 2 SE), a label per family for
  ``production``, ``scan``, ``fdr``, ``total``, ``span``, ``complement``,
  ``pca_matched`` and the three within-group statistics, the selected rank and
  held-out AUROC by design and cell, and every arm's rate per design and cell.

With ``--workers`` above 1, pin BLAS to one thread per process
(``OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1``).

Examples:

    export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1

    # smoke
    python scripts/group_dispersion_check.py run --profile core --reps 2 --out-dir /tmp/gd --workers 7
    python scripts/group_dispersion_check.py run --profile columns --reps 2 --out-dir /tmp/gd --workers 7
    python scripts/group_dispersion_check.py merge --profile core --out-dir /tmp/gd
    python scripts/group_dispersion_check.py merge --profile columns --out-dir /tmp/gd
    python scripts/group_dispersion_check.py decide --out-dir /tmp/gd

    # one SLURM array task (scripts/group_dispersion_array.sbatch), then merge and decide
"""

from __future__ import annotations

import argparse
import dataclasses
import json
import math
import sys
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path
from time import perf_counter
from typing import Any, Iterable, Mapping, Sequence

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import interaction_scan_check as isc  # noqa: E402  (the size observation of one rate row)
import latent_space_comparison as lsc  # noqa: E402  (seeds, settings binding, CSV)
import span_split_check as ssc  # noqa: E402  (parent loading, records, paths)
import stage_span_paper_grade_decision as grade  # noqa: E402  (the count-plus-ceiling rule)
import univariate_baseline_check as ubc  # noqa: E402  (the unit, its arms, its tables)

sys.path.pop(0)

DRIVER = "group_dispersion_check"
DRIVER_VERSION = 1
RULE_VERSION = 1
PROFILES: tuple[str, ...] = ("core", "columns")

BASE_SEED = ubc.STAGE_SPAN_BASE_SEED
DEFAULT_PERMUTATIONS = ubc.DEFAULT_PERMUTATIONS
DEFAULT_REPLICATES = 500

REPO = Path(__file__).resolve().parents[1]
UNIVARIATE_DIR = REPO / "results" / "univariate-baseline-2026-11-23" / "stage_span"
BAND_DIR = REPO / "results" / "stage-effect-scaling-2026-11-24" / "band"
PARENT_FILE_NAMES: tuple[str, ...] = ("settings.json", "rate_table.csv", "cell_table.csv")
#: parent -> (committed directory, source label).
PARENTS: dict[str, tuple[Path, str]] = {
    "univariate_baseline": (UNIVARIATE_DIR, "univariate baseline 2026-11-23, stage_span profile"),
    "stage_effect_band": (BAND_DIR, "stage-effect scaling 2026-11-24, band profile"),
}

#: The parent's design rerun as the reproduction gate, under its own name and seeds.
RERUN_DESIGN = "n200"
CORE_COLUMN = "n200"
CORE_SCALES: tuple[float, ...] = (1.25, 1.5, 2.0)
CORE_SHARES: tuple[float, ...] = (0.5, 0.65, 0.8, 0.35, 0.2)
CORE_EQUAL_SHARES: tuple[float, ...] = (0.65, 0.8)
#: column -> the parent it is read from.
COLUMN_SOURCES: dict[str, str] = {
    "n80": "univariate_baseline",
    "n500": "univariate_baseline",
    "n200_delta0.5": "univariate_baseline",
    "band_n320": "stage_effect_band",
}
COLUMN_SCALES: tuple[float, ...] = (1.5, 2.0)
COLUMN_SHARES: tuple[float, ...] = (0.5, 0.8, 0.2)
COLUMN_EQUAL_SHARES: tuple[float, ...] = (0.8,)
#: The band column's committed stage effect, checked against the committed settings.
BAND_DELTA = 0.05

CONTROL_CELLS: tuple[str, ...] = ubc.CONTROL_CELLS
CORE_CELLS = tuple(cell for cell in ubc.STAGE_SPAN_CELLS if cell[0] in CONTROL_CELLS)
COLUMN_CELLS = tuple(cell for cell in ubc.STAGE_SPAN_CELLS if cell[0] == "anchor")

STATISTICS: tuple[str, ...] = ubc.STATISTICS
FAMILIES: tuple[str, ...] = ("unbalanced_equal", "balanced_unequal", "wider_smaller", "wider_larger")
FAMILY_M_EXPECTED: dict[str, int] = {
    "unbalanced_equal": 24, "balanced_unequal": 42, "wider_smaller": 60, "wider_larger": 60,
}
WITHIN_RULE = ubc.RULE_STRICT
WITHIN_ARMS: tuple[str, ...] = tuple(f"within_group_{statistic}" for statistic in STATISTICS)
RATE_ARMS: tuple[str, ...] = (*ubc.ARMS, *WITHIN_ARMS)
RECORDED_ARMS: tuple[str, ...] = (
    "production", "scan", "fdr", "total", "span", "complement", "pca_matched", *WITHIN_ARMS,
)
WITHIN_COLUMNS: tuple[str, ...] = (
    *ubc.CELL_FIELDS, "units", "parent_completed", "within_completed", "within_failed", "inconsistent",
)

VERDICTS: tuple[str, ...] = (
    "reproduction_failed", "unbalance_excess", "dispersion_excess", "unbalanced_dispersion_excess", "size_holds",
)
_LIMIT = (
    " The readings hold for simulated Gaussian noise at 658 features, one scale applied to every feature of one "
    "group, three stages and no covariates; a difference in correlation structure was not simulated and no real "
    "cohort's dispersion was read. No verdict changes what `motco analyze` runs. Under any verdict the recorded "
    "arms say whether the per-feature tests share the reading, and the roadmap's \"Not yet established\" entries "
    "on unbalance and on noise are updated."
)
CONSEQUENCES: dict[str, str] = {
    "reproduction_failed": (
        "Nothing is read: the rerun of n200 does not reproduce its committed counts, so either the generator is "
        "not unchanged while the scale is absent or the environment drifted. Find the drift before any rerun; a "
        "mismatch is diagnosed, never waived."
    ),
    "unbalance_excess": (
        "Unequal groups lose size before any difference in dispersion: unbalance at equal dispersion is diagnosed "
        "first, as its own change, and the dispersion families are reported and not interpreted." + _LIMIT
    ),
    "dispersion_excess": (
        "Size needs equal group dispersion even with equal groups. The documentation states that a rejection is "
        "evidence of a difference in trajectory or in dispersion, for the statistics and scales the map gives, in "
        "balanced designs too, and the interpretability statement of each failing statistic gains that clause. A "
        "case study that rejects reports the two groups' residual dispersion in the scores beside the rejection. "
        "The roadmap opens a diagnosis of the excess the within-group draws leave, and after it a change for a "
        "null or a fit that holds size, and notes that open item 3's bias correction needs a per-group trace."
        + _LIMIT
    ),
    "unbalanced_dispersion_excess": (
        "Size needs equal group dispersion when the groups are unequal. The documentation states that a rejection "
        "in an unbalanced design is evidence of a difference in trajectory or in dispersion, for the shares, "
        "statistics and scales the map gives, and that balanced designs held up to the largest scale read; the "
        "interpretability statement of each failing statistic gains that clause. A case study that rejects reports "
        "the two groups' residual dispersion in the scores beside the rejection. The roadmap opens a diagnosis of "
        "the excess the within-group draws leave, and notes that open item 3's bias correction needs a per-group "
        "trace." + _LIMIT
    ),
    "size_holds": (
        "The documentation states the range of scales, shares, sample sizes and stage effects over which size "
        "held." + _LIMIT
    ),
}
DECISION_NAME = "group_dispersion_decision.json"


# ── Designs ──────────────────────────────────────────────────────────────────


def design_name(column: str, share: float, scale: float | None) -> str:
    """``<column>_r<share>`` at equal dispersion, ``<column>_r<share>_s<scale>`` otherwise."""

    name = f"{column}_r{share:g}"
    return name if scale is None else f"{name}_s{scale:g}"


def family_of(share: float, scale: float | None) -> str | None:
    """The size family of a design; ``None`` for a balanced design at equal dispersion (no family)."""

    if scale is None:
        return None if share == 0.5 else "unbalanced_equal"
    if share == 0.5:
        return "balanced_unequal"
    # The scale is the second group's, and the share the first group's: above 0.5 the wider group is the smaller.
    return "wider_smaller" if share > 0.5 else "wider_larger"


def _entry(column: str, share: float, scale: float | None, *, design: str | None = None) -> dict[str, Any]:
    return {
        "design": design or design_name(column, share, scale),
        "column": column,
        "share": float(share),
        "scale": None if scale is None else float(scale),
        "covariance_ratio": None if scale is None else float(scale) ** 2,
        "family": None if design is not None else family_of(share, scale),
    }


def _crossed(column: str, shares: Sequence[float], scales: Sequence[float], equal: Sequence[float]) -> list[dict]:
    entries = []
    for share in shares:
        if share in equal:
            entries.append(_entry(column, share, None))
        entries.extend(_entry(column, share, scale) for scale in scales)
    return entries


def core_entries() -> list[dict[str, Any]]:
    """The rerun design, then the 17 new designs, share-major."""

    rerun = _entry(CORE_COLUMN, 0.5, None, design=RERUN_DESIGN)
    return [rerun, *_crossed(CORE_COLUMN, CORE_SHARES, CORE_SCALES, CORE_EQUAL_SHARES)]


def column_entries() -> list[dict[str, Any]]:
    """Seven designs per column, column-major."""

    return [entry for column in COLUMN_SOURCES
            for entry in _crossed(column, COLUMN_SHARES, COLUMN_SCALES, COLUMN_EQUAL_SHARES)]


def design_generator(committed: Mapping[str, Any], share: float, scale: float | None) -> dict[str, Any]:
    """The committed design's generator with only ``group_ratio`` and ``group_dispersion_scale`` changed."""

    generator = {key: value for key, value in committed.items() if key != "group_dispersion_scale"}
    generator["group_ratio"] = float(share)
    if scale is not None:
        generator["group_dispersion_scale"] = float(scale)
    return generator


# ── Parents ──────────────────────────────────────────────────────────────────


def _parent_dirs(parent_dirs: Mapping[str, Path] | None) -> dict[str, Path]:
    dirs = {name: directory for name, (directory, _) in PARENTS.items()}
    dirs.update({name: Path(path) for name, path in (parent_dirs or {}).items()})
    return dirs


def parent_binding(parent_dirs: Mapping[str, Path] | None = None) -> tuple[dict[str, Any], dict[str, Any]]:
    """Each parent's bound files and settings hash, and its committed settings (hash-checked)."""

    binding: dict[str, Any] = {}
    committed: dict[str, Any] = {}
    for name, directory in _parent_dirs(parent_dirs).items():
        stored = ssc.load_parent(directory)
        binding[name] = {
            "dir": ssc._display_path(directory),
            "files": ubc._bind(directory / file_name for file_name in PARENT_FILE_NAMES),
            "settings_hash": stored["hash"],
            "source": PARENTS[name][1],
        }
        committed[name] = stored["settings"]
    return binding, committed


def _parent_units(settings: Mapping[str, Any]) -> Mapping[str, Any]:
    """A parent's parent-shaped unit settings: the baseline's own, the scaling run's ``unit_settings``."""

    return settings["unit_settings"] if "unit_settings" in settings else settings


def committed_generator(column: str, committed: Mapping[str, Any], directory: Path) -> dict[str, Any]:
    """The committed generator of ``column`` from its parent's settings, checked against what it must be."""

    units = _parent_units(committed)["units"]
    designs = units["designs"]
    if column not in designs:
        raise SystemExit(f"{directory}/settings.json has no design {column!r}.")
    generator = dict(designs[column]["generator"])
    if column in ubc.STAGE_SPAN_DESIGNS:
        if generator != ubc.stage_span_generator(*ubc.STAGE_SPAN_DESIGNS[column]):
            raise SystemExit(f"{directory}/settings.json: design {column!r} is not the univariate baseline's.")
    elif column == "band_n320":
        deltas = {float(generator[key]) for key in ubc.STAGE_EFFECT_KEYS}
        if deltas != {BAND_DELTA}:
            raise SystemExit(f"{directory}/settings.json: band_n320 has stage effect {sorted(deltas)}, not the "
                             f"committed {BAND_DELTA}.")
    if designs[column].get("layers") is not None or float(generator["group_ratio"]) != 0.5:
        raise SystemExit(f"{directory}/settings.json: design {column!r} is not a balanced three-block design.")
    return generator


def _check_parent_fit(committed: Mapping[str, Any], reference: Mapping[str, Any], directory: Path) -> None:
    """The parent's units were fitted, seeded and planted as this driver's are (its grade may differ)."""

    units = _parent_units(committed)
    keys = ("pls_cv_params", "min_components", "orientation_estimator", "base_seed")
    differing = [f"units.{key}" for key in keys if units["units"][key] != reference["units"][key]]
    cells = {cell[0]: list(cell) for cell in units["units"]["cells"]}
    differing += [f"cells.{cell[0]}" for cell in CORE_CELLS if cells.get(cell[0]) != list(cell)]
    differing += [key for key in ("recorded_spaces", "rules", "population_null_tolerance")
                  if units[key] != reference[key]]
    if differing:
        raise SystemExit(f"{directory}/settings.json differs from this driver's unit in {differing}; its designs "
                         "would not be read as the parent read them.")


# ── Settings ─────────────────────────────────────────────────────────────────


def unit_settings(
    designs: Mapping[str, Mapping[str, Any]],
    cells: Sequence[Sequence[Any]],
    *,
    replicates: int,
    permutations: int,
) -> dict[str, Any]:
    """The parent's ``stage_span`` settings with only the designs and cells replaced.

    Every other key (rules, arms, CV knobs, rank floor, estimator, recorded space, base seed and seed functions) is
    the one ``univariate_baseline_check.stage_span_settings`` writes, so ``evaluate_unit`` treats the units as the
    parent's. Every cell is a control cell, so every cell runs ``replicates`` units.
    """

    settings = ubc.stage_span_settings(
        designs=[RERUN_DESIGN], replicates=replicates, control_replicates=replicates,
        permutations=permutations, base_seed=BASE_SEED, recorded_spaces=ubc.DEFAULT_RECORDED_SPACES,
    )
    settings["units"]["designs"] = {name: {"generator": dict(generator), "layers": None}
                                    for name, generator in designs.items()}
    settings["units"]["cells"] = [list(cell) for cell in cells]
    settings["designs"] = list(designs)
    return settings


def _select(entries: Sequence[dict[str, Any]], designs: Sequence[str] | None, profile: str) -> list[dict[str, Any]]:
    if designs is None:
        return list(entries)
    known = [entry["design"] for entry in entries]
    unknown = sorted(set(designs) - set(known))
    if unknown:
        raise SystemExit(f"Unknown {profile} design(s) {unknown}; expected a subset of {known}.")
    return [entry for entry in entries if entry["design"] in designs]


def profile_settings(
    profile: str,
    *,
    replicates: int = DEFAULT_REPLICATES,
    permutations: int = DEFAULT_PERMUTATIONS,
    designs: Sequence[str] | None = None,
    parent_dirs: Mapping[str, Path] | None = None,
) -> dict[str, Any]:
    """One profile's settings: its designs on the committed generators, both parents bound."""

    if profile not in PROFILES:
        raise SystemExit(f"Unknown profile {profile!r}; expected one of {list(PROFILES)}.")
    dirs = _parent_dirs(parent_dirs)
    binding, committed = parent_binding(parent_dirs)
    entries = _select(core_entries() if profile == "core" else column_entries(), designs, profile)
    cells = CORE_CELLS if profile == "core" else COLUMN_CELLS
    generators: dict[str, dict[str, Any]] = {}
    for entry in entries:
        source = COLUMN_SOURCES.get(entry["column"], "univariate_baseline")
        base = committed_generator(entry["column"], committed[source], dirs[source])
        if entry["design"] == RERUN_DESIGN:
            generators[entry["design"]] = base
        else:
            generators[entry["design"]] = design_generator(base, entry["share"], entry["scale"])
    units = unit_settings(generators, cells, replicates=replicates, permutations=permutations)
    for name in PARENTS:
        _check_parent_fit(committed[name], units, dirs[name])
    files: dict[str, str] = {}
    for block in binding.values():
        files.update(block["files"])
    return {
        "driver": DRIVER,
        "driver_version": DRIVER_VERSION,
        "profile": profile,
        "parents": binding,
        "parent_files": files,
        "rerun_design": RERUN_DESIGN if profile == "core" and RERUN_DESIGN in generators else None,
        "design_index": {entry["design"]: {key: entry[key] for key in entry if key != "design"} for entry in entries},
        "judged_statistics": list(STATISTICS),
        "families": list(FAMILIES),
        "unit": "univariate_baseline_check.evaluate_unit on unit_settings (the parent's stage_span settings with "
                "the designs and cells replaced), pca_matched the one recorded space",
        "within_group": {
            "strata": "group",
            "rule": WITHIN_RULE,
            "definition": "the dataset regenerated, production refitted at the recorded rank with no "
                          "cross-validation (forced_components), delta, angle and shape under "
                          "motco.stats.permutation.RRPP(..., strata=<group positions>) at the unit's RRPP seed, "
                          "p-values by the harness's empirical rule; recorded, never judged",
        },
        "unit_settings": units,
    }


def check_bound_files(settings: Mapping[str, Any]) -> None:
    """Refuse, by name, a changed parent file."""

    ubc.check_parent_files(settings)


def profile_dir(out_dir: Path, profile: str) -> Path:
    return Path(out_dir) / profile


# ── The within-group arm ─────────────────────────────────────────────────────


def group_strata(metadata: Any, group_col: str, stage_col: str) -> np.ndarray:
    """Each sample's group position in the design's group levels."""

    from motco.simulations.evaluation import build_simulation_trajectory_design

    design = build_simulation_trajectory_design(metadata, group_col=group_col, stage_col=stage_col)
    position = {str(level): index for index, level in enumerate(design.group_levels)}
    return np.array([position[str(label)] for label in metadata[group_col].tolist()], dtype=int)


def strata_reading(dataset: Any, matrix: Any, params: Any, *, strata: np.ndarray) -> dict[str, Any]:
    """``delta``, ``angle`` and ``shape`` and their p-values under draws restricted to ``strata``.

    The harness's statistics and p-values (``evaluate_semisynthetic_trajectory``) with only ``strata`` added to
    the RRPP call: one stratum holding every sample reproduces its p-values exactly.
    """

    from motco.simulations.evaluation import (
        _empirical_p_value,
        _extract_null_distributions,
        _extract_pair_statistics,
        build_simulation_trajectory_design,
    )
    from motco.stats.permutation import RRPP
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
        return_eigengaps=True, orientation_estimator=estimator, strata=strata,
    )
    nulls = _extract_null_distributions(dist_delta, dist_angle, dist_shape, shape_available=shape_available)
    p_values = {statistic: _empirical_p_value(nulls[statistic], value) for statistic, value in observed.items()
                if statistic in nulls and np.isfinite(value)}
    return {"pair_statistics": lsc._finite(observed), "p_values": lsc._finite(p_values)}


def forced_rank_params(params: Any, rank: int) -> Any:
    """The unit's production parameters at a forced rank: the same final fit, no cross-validation."""

    return dataclasses.replace(params, integration_params={**dict(params.integration_params or {}),
                                                           "forced_components": int(rank)})


def same_statistics(left: Mapping[str, Any], right: Mapping[str, Any]) -> bool:
    """Exact equality of ``delta``, ``angle`` and ``shape``; two unavailable values agree."""

    return all(left.get(statistic) == right.get(statistic) for statistic in STATISTICS)


def within_group_reading(unit: lsc.Unit, units: Mapping[str, Any], record: Mapping[str, Any]) -> dict[str, Any]:
    """The within-group arm of one completed unit; a failure is recorded, never raised."""

    from motco.simulations.evaluation import integrate_semisynthetic_dataset

    try:
        fields = ubc.unit_fields(unit, units)
        dataset, params = ubc.unit_inputs(unit, units, fields)
        rank = int(record["production"]["rank"])
        latent = integrate_semisynthetic_dataset(dataset, forced_rank_params(params, rank))
        strata = group_strata(dataset.metadata, params.group_col, params.stage_col)
        reading = strata_reading(dataset, latent.matrix, params, strata=strata)
    except Exception as exc:  # recorded for this arm only
        return {"status": "failed", "error_type": type(exc).__name__, "error_message": str(exc)}
    return {
        "status": "completed",
        "strata": "group",
        "rank": rank,
        "same_statistics": same_statistics(reading["pair_statistics"], record["production"]["pair_statistics"]),
        **reading,
        "rejects": {statistic: ubc.rejects(reading["p_values"].get(statistic), WITHIN_RULE)
                    for statistic in STATISTICS},
    }


# ── Units and the run ────────────────────────────────────────────────────────


def profile_units(settings: Mapping[str, Any]) -> list[lsc.Unit]:
    """Design-major, then cell, then replicate (the parent's order)."""

    return ubc.profile_units(settings["unit_settings"])


def evaluate_unit(unit: lsc.Unit, units: Mapping[str, Any]) -> dict[str, Any]:
    """The parent's unit, unchanged, then the within-group arm on the same fit, statistics and seed."""

    start = perf_counter()
    record = dict(ubc.evaluate_unit(unit, units))
    mark = perf_counter()
    if record.get("status") == "completed":
        record["within_group"] = within_group_reading(unit, units, record)
    else:
        record["within_group"] = {"status": "skipped", "reason": "the parent's unit failed"}
    runtime = dict(record.get("runtime_seconds") or {})
    runtime["parent_unit"] = runtime.get("unit")
    runtime["within_group"] = perf_counter() - mark
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

    profile = settings["profile"]
    check_bound_files(settings)
    directory = profile_dir(out_dir, profile)
    lsc.ensure_settings(directory, settings)
    units_settings = settings["unit_settings"]
    path = lsc.shard_path(directory, index, count)
    path.parent.mkdir(parents=True, exist_ok=True)
    done = ssc.read_records([path], ubc.record_key) if path.exists() else {}
    # Rewritten with complete lines only, so an interrupted unit is re-run, not duplicated.
    with path.open("w", encoding="utf-8") as handle:
        for key in sorted(done):
            handle.write(lsc._line(done[key]))
    units = profile_units(settings)
    pending = [unit for position, unit in enumerate(units)
               if position % count == index and ubc.unit_key(unit) not in done]
    if limit is not None:
        pending = pending[:limit]

    def append(position: int, record: Mapping[str, Any]) -> None:
        with path.open("a", encoding="utf-8") as handle:
            handle.write(lsc._line(record))
            handle.flush()
        print(f"[{profile} {index}/{count}] {position}/{len(pending)} {ubc.record_key(record)} {record['status']}",
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


def _within(record: Mapping[str, Any]) -> Mapping[str, Any]:
    return record.get("within_group") or {}


def within_group_rows(
    records: Sequence[Mapping[str, Any]], units: Mapping[str, Any]
) -> tuple[list[dict[str, Any]], list[dict[str, Any]]]:
    """The within-group arm's rate rows and its per-cell accounting.

    A rate counts the completed, consistent units; ``failures`` counts every other unit of the cell (a failed
    parent unit, a failed arm, an inconsistent refit).
    """

    rates, table = [], []
    for fields, group in ubc._cell_groups(records, units):
        completed = [r for r in group if r.get("status") == "completed"]
        read = [r for r in completed if _within(r).get("status") == "completed"]
        consistent = [r for r in read if _within(r).get("same_statistics") is True]
        for statistic, arm in zip(STATISTICS, WITHIN_ARMS, strict=True):
            rates.append({
                **fields, "arm": arm, "rule": WITHIN_RULE,
                **ssc._rate_row([(_within(r).get("rejects") or {}).get(statistic) for r in consistent]),
                "failures": len(group) - len(consistent),
            })
        table.append({
            **fields, "units": len(group), "parent_completed": len(completed), "within_completed": len(read),
            "within_failed": len(completed) - len(read), "inconsistent": len(read) - len(consistent),
        })
    return rates, table


def rate_rows(records: Sequence[Mapping[str, Any]], units: Mapping[str, Any]) -> list[dict[str, Any]]:
    """The parent's rate rows with the within-group rows after each cell's parent arms."""

    parent = ubc.rate_rows(records, units)
    within, _ = within_group_rows(records, units)
    designs = {design: position for position, design in enumerate(units["designs"])}
    cells = {cell: position for position, cell in enumerate(ubc._cell_names(units))}
    arms = {arm: position for position, arm in enumerate(RATE_ARMS)}
    return sorted([*parent, *within],
                  key=lambda row: (designs[str(row["design"])], cells[str(row["cell"])], arms[str(row["arm"])]))


def within_accounting(table: Sequence[Mapping[str, Any]]) -> dict[str, Any]:
    def by_cell(key: str) -> dict[str, int]:
        return {f"{row['design']}/{row['cell']}": int(row[key]) for row in table if int(row[key])}

    return {
        "completed": sum(int(row["within_completed"]) for row in table),
        "failed": sum(int(row["within_failed"]) for row in table),
        "inconsistent": sum(int(row["inconsistent"]) for row in table),
        "failed_by_cell": by_cell("within_failed"),
        "inconsistent_by_cell": by_cell("inconsistent"),
    }


def rerun_reproduction(
    rates: Sequence[Mapping[str, Any]], committed: Sequence[Mapping[str, Any]], settings: Mapping[str, Any]
) -> dict[str, Any]:
    """The rerun design's available and rejected counts against the committed table, both cells, every parent arm."""

    parent = settings["parents"]["univariate_baseline"]
    base = {"design": RERUN_DESIGN, "source": parent["source"], "files": dict(parent["files"])}
    if settings.get("rerun_design") != RERUN_DESIGN:
        return {**base, "passed": False, "rows_compared": 0, "differing": [], "unresolved_rows": [],
                "reason": f"the rerun design {RERUN_DESIGN!r} is not in this run"}
    wanted_arms = set(ubc.ARMS)
    mine = {(str(r["cell"]), str(r["arm"])): r for r in rates
            if str(r["design"]) == RERUN_DESIGN and str(r["arm"]) in wanted_arms}
    theirs = {(str(r["cell"]), str(r["arm"])): r for r in committed
              if str(r["design"]) == RERUN_DESIGN and str(r["cell"]) in CONTROL_CELLS}
    differing, unresolved = [], []
    compared = 0
    for key in [*mine, *(key for key in theirs if key not in mine)]:
        if key not in mine or key not in theirs:
            unresolved.append({"cell": key[0], "arm": key[1],
                               "missing_from": "this run" if key not in mine else "the committed table"})
            continue
        compared += 1
        wanted, got = ubc._counts(theirs[key]), ubc._counts(mine[key])
        if got != wanted:
            assert wanted is not None
            differing.append({"cell": key[0], "arm": key[1],
                              "committed": {"count": wanted[0], "rejections": wanted[1]},
                              "reproduced": None if got is None else {"count": got[0], "rejections": got[1]}})
    return {**base, "passed": compared > 0 and not differing and not unresolved, "rows_compared": compared,
            "differing": differing, "unresolved_rows": unresolved}


def merge(out_dir: Path, profile: str) -> dict[str, Any]:
    """Merge one profile's shards and write its tables, accounting and (``core``) reproduction check."""

    directory = profile_dir(out_dir, profile)
    settings = lsc.load_settings(directory)
    check_bound_files(settings)
    units = settings["unit_settings"]
    unit_hash = lsc.settings_hash(units)
    expected_keys = [ubc.unit_key(unit) for unit in profile_units(settings)]
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
    _, within_table = within_group_rows(records, units)
    summary = {
        **ubc.accounting(records, units, expected=len(expected_keys),
                         duplicates=sum(count - 1 for count in seen.values())),
        "profile": profile,
        "settings_hash": lsc.settings_hash(settings),
        "unit_settings_hash": unit_hash,
        "within_group": within_accounting(within_table),
    }
    rates = rate_rows(records, units)
    lsc._write_csv(directory / "rate_table.csv", rates, ubc.RATE_COLUMNS)
    lsc._write_csv(directory / "cell_table.csv", ubc.cell_rows(records, units), ubc.CELL_COLUMNS)
    lsc._write_csv(directory / "within_group_table.csv", within_table, WITHIN_COLUMNS)
    lsc._write_json(directory / "accounting.json", summary)
    result = None
    if profile == "core":
        parent = settings["parents"]["univariate_baseline"]
        committed = lsc._read_csv(ssc._resolve(parent["dir"]) / "rate_table.csv")
        result = rerun_reproduction(rates, committed, settings)
        lsc._write_json(directory / "reproduction.json", result)
    return {
        "profile": profile,
        **{key: summary[key] for key in ("complete_units", "missing_units", "failed_units", "recorded_core_hours")},
        "within_group_inconsistent": summary["within_group"]["inconsistent"],
        "reproduction_passed": None if result is None else result["passed"],
        "reproduction_differing": None if result is None else len(result["differing"]),
    }


# ── Decision (pure functions) ────────────────────────────────────────────────


def _blocks(core: Mapping[str, Any], columns: Mapping[str, Any]) -> list[tuple[str, Mapping[str, Any]]]:
    return [("core", core), ("columns", columns)]


def family_readings(
    core: Mapping[str, Any], columns: Mapping[str, Any], arms: Sequence[str] = STATISTICS
) -> dict[str, list[dict[str, Any]]]:
    """Every family's observations: each new design's control cells, one per arm; the rerun design enters none."""

    out: dict[str, list[dict[str, Any]]] = {family: [] for family in FAMILIES}
    for profile, block in _blocks(core, columns):
        rates = ubc._index(block["rates"], "design", "cell", "arm")
        units = block["settings"]["unit_settings"]
        for design in units["designs"]:
            meta = block["settings"]["design_index"][design]
            if meta["family"] is None:
                continue
            for cell in ubc._cell_names(units):
                for arm in arms:
                    out[meta["family"]].append(isc._size_observation(
                        rates.get((design, cell, arm)), profile=profile, design=design, column=meta["column"],
                        share=meta["share"], scale=meta["scale"], cell=cell, arm=arm,
                    ))
    return out


def size_reading(core: Mapping[str, Any], columns: Mapping[str, Any]) -> dict[str, Any]:
    """``delta``, ``angle`` and ``shape`` judged in the four families."""

    readings = family_readings(core, columns)
    return {family: {**grade.judge_rates(readings[family]), "m_expected": FAMILY_M_EXPECTED[family]}
            for family in FAMILIES}


def conservative_limit(n_available: int) -> float:
    """The rate below which a reading is conservative: ``α − k·SE`` at ``n_available`` (the bound's mirror)."""

    alpha, k = grade.ALPHA, grade.SE_MULTIPLIER
    return alpha - k * math.sqrt(alpha * (1.0 - alpha) / n_available)


def conservative_readings(size: Mapping[str, Any]) -> list[dict[str, Any]]:
    """Every judged reading below ``α − 2·SE``, flagged; its family's judgement is unchanged."""

    out = []
    for family in FAMILIES:
        for row in size[family]["observations"]:
            if row.get("status") == "unavailable":
                continue
            limit = conservative_limit(int(row["n_available"]))
            if float(row["rate"]) < limit:
                out.append({"family": family, **{key: row[key] for key in (
                    "profile", "design", "column", "share", "scale", "cell", "arm", "n_rejected", "n_available",
                    "rate")}, "limit": limit})
    return out


def _column_order(core: Mapping[str, Any], columns: Mapping[str, Any]) -> list[str]:
    """The core column, then the columns in their declared order (stored settings sort their keys)."""

    present = {meta["column"] for _, block in _blocks(core, columns)
               for meta in block["settings"]["design_index"].values()}
    return [column for column in (CORE_COLUMN, *COLUMN_SOURCES) if column in present]


def failing_scale_map(size: Mapping[str, Any], column_order: Sequence[str]) -> list[dict[str, Any]]:
    """Per statistic, column and share: the smallest scale whose reading is above its family's ceiling.

    A design's reading is above the ceiling when any of its control cells is; ``by_cell`` gives each cell's own
    smallest failing scale. ``none_failed`` is true when no scaled reading there is above its ceiling, up to
    ``largest_scale_read``. The equal-dispersion reading at that share, where one was run, is given beside it.
    """

    grouped: dict[tuple[str, str, float], list[Mapping[str, Any]]] = {}
    for family in FAMILIES:
        for row in size[family]["observations"]:
            grouped.setdefault((str(row["arm"]), str(row["column"]), float(row["share"])), []).append(row)
    shares = sorted({key[2] for key in grouped}, key=lambda s: (CORE_SHARES + COLUMN_SHARES).index(s))
    out = []
    for statistic in STATISTICS:
        for column in column_order:
            for share in shares:
                rows = grouped.get((statistic, column, share))
                if not rows:
                    continue
                scaled = [row for row in rows if row["scale"] is not None]
                if not scaled:
                    continue
                scales = sorted({float(row["scale"]) for row in scaled})
                failing = sorted({float(row["scale"]) for row in scaled if row.get("status") == "above_ceiling"})
                cells = list(dict.fromkeys(str(row["cell"]) for row in scaled))
                by_cell = {}
                for cell in cells:
                    cell_failing = sorted({float(row["scale"]) for row in scaled
                                           if str(row["cell"]) == cell and row.get("status") == "above_ceiling"})
                    by_cell[cell] = cell_failing[0] if cell_failing else None
                equal = {str(row["cell"]): row.get("status") for row in rows if row["scale"] is None}
                out.append({
                    "statistic": statistic, "column": column, "share": share,
                    "family": family_of(share, scales[0]),
                    "scales_read": scales,
                    "largest_scale_read": scales[-1],
                    "first_failing_scale": failing[0] if failing else None,
                    "none_failed": not failing,
                    "unavailable_scales": sorted({float(row["scale"]) for row in scaled
                                                  if row.get("status") == "unavailable"}),
                    "by_cell": by_cell,
                    "rates": [{"scale": float(row["scale"]), "cell": row["cell"], "status": row.get("status"),
                               "rate": row.get("rate")} for row in sorted(scaled, key=lambda r: float(r["scale"]))],
                    "equal_dispersion": equal or None,
                })
    return out


def recorded_arm_labels(core: Mapping[str, Any], columns: Mapping[str, Any]) -> dict[str, Any]:
    """Each recorded arm's size label per family, under the same rule; never part of the verdict."""

    out: dict[str, Any] = {}
    for arm in RECORDED_ARMS:
        readings = family_readings(core, columns, arms=(arm,))
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
                "above_ceiling": [{key: row.get(key) for key in ("design", "cell", "rate")}
                                  for row in judged["observations"] if row.get("status") == "above_ceiling"],
            }
        out[arm] = {"rule": ubc.RULES.get(arm, WITHIN_RULE), "families": labels}
    return out


def fit_by_design(core: Mapping[str, Any], columns: Mapping[str, Any]) -> list[dict[str, Any]]:
    """The selected rank and held-out AUROC per design and cell."""

    out = []
    for profile, block in _blocks(core, columns):
        cells = ubc._index(block["cells"], "design", "cell")
        units = block["settings"]["unit_settings"]
        for design in units["designs"]:
            meta = block["settings"]["design_index"][design]
            for cell in ubc._cell_names(units):
                row = cells.get((design, cell)) or {}
                out.append({"profile": profile, "design": design, "cell": cell, **meta,
                            **{key: ubc._finite(row.get(key)) for key in (
                                "completed", "rank_q10", "rank_median", "rank_q90", "auroc_median")}})
    return out


def arm_rates(core: Mapping[str, Any], columns: Mapping[str, Any]) -> list[dict[str, Any]]:
    """Every arm's rate per design and cell."""

    out = []
    for profile, block in _blocks(core, columns):
        rates = ubc._index(block["rates"], "design", "cell", "arm")
        units = block["settings"]["unit_settings"]
        for design in units["designs"]:
            meta = block["settings"]["design_index"][design]
            for cell in ubc._cell_names(units):
                out.append({"profile": profile, "design": design, "cell": cell, **meta,
                            "rates": {arm: ubc._rate_of(rates.get((design, cell, arm))) for arm in RATE_ARMS}})
    return out


def decision(core: Mapping[str, Any], columns: Mapping[str, Any]) -> dict[str, Any]:
    """Apply rule version 1.

    ``core`` and ``columns`` hold ``settings``, ``accounting``, ``rates`` and ``cells`` (rows); ``core`` also
    ``reproduction``. Every family is recorded whatever the verdict.
    """

    size = size_reading(core, columns)
    reproduction = core.get("reproduction") or {}
    failing: dict[str, Any] = {}
    holds = {family: size[family]["holds"] for family in FAMILIES}
    if not reproduction.get("passed"):
        verdict = "reproduction_failed"
        failing["reproduction"] = {key: reproduction.get(key) for key in ("differing", "unresolved_rows", "reason")
                                   if reproduction.get(key)}
    elif not holds["unbalanced_equal"]:
        verdict = "unbalance_excess"
    elif not holds["balanced_unequal"]:
        verdict = "dispersion_excess"
    elif not (holds["wider_smaller"] and holds["wider_larger"]):
        verdict = "unbalanced_dispersion_excess"
    else:
        verdict = "size_holds"
    for family in FAMILIES:
        if not holds[family]:
            failing[family] = ssc._failing_observations(size[family])
    blocks = dict(_blocks(core, columns))
    return {
        "rule_version": RULE_VERSION,
        "alpha": ubc.ALPHA,
        "judged_statistics": list(STATISTICS),
        "rule": "stage_span_paper_grade_decision.judge_rates (count plus Bonferroni ceiling) per family",
        "settings_hashes": {name: lsc.settings_hash(block["settings"]) for name, block in blocks.items()},
        "parents": core["settings"].get("parents"),
        "complete": all(block["accounting"]["missing_units"] == 0 for block in blocks.values()),
        "failed_units": {name: block["accounting"]["failed_units"] for name, block in blocks.items()},
        "within_group": {name: (block["accounting"].get("within_group") or {}) for name, block in blocks.items()},
        "reproduction": {key: reproduction.get(key) for key in ("passed", "rows_compared")},
        "verdict": verdict,
        "consequence": CONSEQUENCES[verdict],
        "consequences": dict(CONSEQUENCES),
        "families": {family: {key: size[family][key] for key in (
            "holds", "m", "m_expected", "median_n_available", "allowed_exceedances", "observed_exceedances",
            "ceiling_at_median", "n_above_ceiling", "n_unavailable")} for family in FAMILIES},
        "failing": failing,
        "criteria": {"size": size},
        "map": failing_scale_map(size, _column_order(core, columns)),
        "conservative": conservative_readings(size),
        "recorded_arms": recorded_arm_labels(core, columns),
        "fit": fit_by_design(core, columns),
        "rates": arm_rates(core, columns),
    }


def _profile_tables(out_dir: Path, profile: str) -> dict[str, Any]:
    directory = profile_dir(out_dir, profile)
    settings = lsc.load_settings(directory)
    check_bound_files(settings)
    path = directory / "accounting.json"
    if not path.exists():
        raise SystemExit(f"{path} does not exist; run `merge --profile {profile}` first.")
    summary = json.loads(path.read_text(encoding="utf-8"))
    if summary.get("settings_hash") != lsc.settings_hash(settings):
        raise SystemExit(f"{path} was merged under other settings; run `merge` again.")
    block: dict[str, Any] = {
        "settings": settings,
        "accounting": summary,
        "rates": ubc._read_csv(directory / "rate_table.csv"),
        "cells": ubc._read_csv(directory / "cell_table.csv"),
    }
    if profile == "core":
        block["reproduction"] = json.loads((directory / "reproduction.json").read_text(encoding="utf-8"))
    return block


def decide(out_dir: Path) -> dict[str, Any]:
    """Read both merged profiles; write the decision."""

    core = _profile_tables(out_dir, "core")
    columns = _profile_tables(out_dir, "columns")
    if core["settings"]["parent_files"] != columns["settings"]["parent_files"]:
        raise SystemExit("The two profiles were written against different parent files; refusing.")
    result = decision(core, columns)
    lsc._write_json(Path(out_dir) / DECISION_NAME, result)
    return result


# ── CLI ──────────────────────────────────────────────────────────────────────


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = parser.add_subparsers(dest="command", required=True)
    run = sub.add_parser("run", help="Evaluate (a shard of) one profile's units")
    run.add_argument("--profile", choices=PROFILES, required=True)
    run.add_argument("--out-dir", type=Path, required=True)
    run.add_argument("--shard", type=str, default="0/1", help="i/N, round-robin over the units")
    run.add_argument("--limit", type=int, default=None, help="Evaluate at most this many pending units")
    run.add_argument("--workers", type=int, default=1)
    run.add_argument("--perms", type=int, default=DEFAULT_PERMUTATIONS)
    run.add_argument("--reps", type=int, default=DEFAULT_REPLICATES,
                     help=f"Replicates of every cell (default {DEFAULT_REPLICATES})")
    run.add_argument("--designs", type=str, default=None, help="Comma-separated designs; default: all")
    merge_parser = sub.add_parser("merge", help="Merge one profile's shards")
    merge_parser.add_argument("--profile", choices=PROFILES, required=True)
    merge_parser.add_argument("--out-dir", type=Path, required=True)
    decide_parser = sub.add_parser("decide", help=f"Write {DECISION_NAME} from both profiles")
    decide_parser.add_argument("--out-dir", type=Path, required=True)
    return parser


def settings_from_args(args: argparse.Namespace) -> dict[str, Any]:
    chosen = None if args.designs is None else [name for name in args.designs.split(",") if name]
    return profile_settings(args.profile, replicates=args.reps, permutations=args.perms, designs=chosen)


def _summary(result: Mapping[str, Any]) -> dict[str, Any]:
    return {
        "verdict": result["verdict"],
        "reproduction": result["reproduction"],
        "families": {family: {key: entry[key] for key in ("holds", "m", "observed_exceedances", "n_above_ceiling",
                                                           "n_unavailable")}
                     for family, entry in result["families"].items()},
        "first_failing": {f"{e['statistic']}/{e['column']}/r{e['share']:g}": e["first_failing_scale"]
                          for e in result["map"]},
        "complete": result["complete"],
    }


def main(argv: Iterable[str] | None = None) -> int:
    args = build_parser().parse_args(None if argv is None else list(argv))
    if args.command == "run":
        index, count = lsc.parse_shard(args.shard)
        evaluated = run_shard(args.out_dir, settings_from_args(args), index, count, limit=args.limit,
                              workers=args.workers)
        print(f"evaluated {evaluated} {args.profile} unit(s)", file=sys.stderr)
        return 0
    if args.command == "merge":
        print(json.dumps(merge(args.out_dir, args.profile), indent=1))
        return 0
    print(json.dumps(_summary(decide(args.out_dir)), indent=1))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
