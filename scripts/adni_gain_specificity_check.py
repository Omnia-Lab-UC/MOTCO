#!/usr/bin/env python3
"""Gain, protected power and specificity at ADNI's design: plants on ADNI's blocks and stage-span simulation (no R).

The driver behind ``openspec/changes/adni-gain-and-specificity/`` (spec
``cohort-design-power-validation``, roadmap open items 1 and 2). It reads
``predictive_orthonormal`` against production's scores (``pls_cv``) at ADNI's
design in two profiles, under one decision fixed before any unit runs, and
states what production itself detects on ADNI's noise. No package code
changes; every arm is one ``analyze_trajectories`` call on existing options.

**Profile ``plants``** (designs D2–D7): ADNI's blocks under within-stage-permuted
sex labels with a known shift added to the second declared group's rows, stage
by stage. The export is read in place through ``adni_anchor_check``'s contract
and units, imported unchanged (base seed 2039, the anchor's seed family), so a
cohort anchor unit is the anchor run's unit. Points: ``n750`` (750 participants
subsampled within cells; the other 751 are its held-out reference) and
``cohort`` (all 1,501; no held-out reference). A plant edits a
covariate-adjusted stage configuration — the stage coefficients of
``features ~ stage indicators + centred covariate columns`` in the analysis
sample's pooled-standardized coordinates, the covariates coded as
``analyze_trajectories`` codes them — through the generator's
``edit_stage_configuration``, and maps the edit back through the pooled SD.
References: ``sample`` (the unit's own analysis sample: an exact off-target
truth) and ``heldout`` (the participants outside it: a direction no fit was
built on; ``n750`` only). Cells ``<reference>_<mode>_e<effect>``:

- controls (500 replicates): ``anchor``; ``sample_translation_e1``;
- plants (200): ``sample_`` magnitude 0.25 / 0.5 / 1, orientation 0.25 / 0.5 /
  1, shape 0.5 / 1 at both points; ``heldout_`` magnitude 0.5 / 1 / 2,
  orientation 0.5 / 1, shape 1 at ``n750``;
- recorded, never judged (200): ``sample_out_of_span_e0.25`` and ``e0.5``, a
  rotation toward a seeded unit direction orthogonal to the configuration's
  centred stage means (``alignment_degrees`` 90).

4,200 units at ``n750`` and 3,000 at ``cohort``. Own seeds (sha256 of
``[driver name, 2040, point, replicate, purpose]``): ``direction`` (the
out-of-span direction, shared by both effects of a replicate) and
``reference_split`` (the held-out halves of the reliability reading).

**Profile ``simulation``** (design D9): ``adni_design_check``'s construction,
imported unchanged at its calibrated strengths and ADNI's table, with
``construction: stage_span``, ``span_layers`` the measured blocks and
``span_scale: planted``; base seed 2040. Cells: ``anchor`` and
``population_translation_e1`` (500), ``population_`` magnitude 0.5 / 1 / 2,
orientation 0.25 / 0.5 / 1 and shape 1 / 2 (200): 2,600 units. It reads no ADNI
file and runs on the cluster (``scripts/adni_gain_specificity_array.sbatch``).

**Arms** (design D8), on one fit at the case study's configuration with the
``df`` correction and the unit's RRPP seed: ``pls_cv`` (cross-validated on the
unit's data; a ``cohort`` anchor unit is forced at the pre-flight's rank 8),
``predictive_orthonormal`` forced at its rank, and both under within-stage
draws (``__within_stage``, recorded). One cross-validation per unit.

**Checks** (design D11): the anchor driver's ``rank``, ``same_fit``,
``measured_dimension``, ``same_statistics`` and ``cell_table``, imported, plus
``plant`` (plants: the second group's raw stage means moved by the declared
shift and the first group's did not, in pooled SDs to 1e-8) or
``construction`` (simulation: the truth records the stage-span construction at
the cell's mode and effect). A unit failing one is ``inconsistent`` and leaves
every rate.

**Readings** (design D16, never decisive): per plant its norm per stage, the
largest relative change of a pooled SD, and the shares of the plant and of its
reference configuration inside the fitted span and the predictive subspace
(``LinearLatentMap``); on anchor units the shares of each configuration and, at
``n750``, the held-out configuration's reliability between two halves split
within stage.

**Parents** (sha256 in both profiles' settings; every subcommand refuses a
changed file by name): the ADNI anchor run's ``rate_table.csv`` (the
reproduction) and ``adni_anchor_decision.json`` (the size scope ``reopened``
cites), and the ADNI-design run's ``settings.json`` and
``adni_design_pvalues.csv`` (the simulation pre-flight). The export's five
files are bound in the plants settings.

Subcommands (``--profile plants | simulation`` on ``preflight``, ``run`` and
``merge``; ``--export-dir`` required wherever the plants profile is read, no
default; each profile writes under ``<out>/<profile>/``):

- ``preflight`` — ``settings.json`` and ``preflight.json``. Plants: the anchor
  driver's pre-flight units (cohort replicates 0 and 1 select rank 8 at AUROC
  0.599514 and equal the forced call). Simulation: anchor replicates 0 and 1 of
  the ADNI-design run, through its own evaluator under its committed settings,
  give its committed p-values exactly. ``run`` refuses to start without a pass.
- ``run`` — ``records/shard_<i>_of_<N>.jsonl``; round-robin shards, resume of
  complete units, refusal of other settings naming the keys.
- ``merge`` — ``merged.jsonl`` (keep off the repository), ``rate_table.csv``,
  ``paired_table.csv``, ``cell_table.csv``, ``accounting.json`` and, for
  plants, ``reproduction.json`` (the cohort anchor's counts against the anchor
  run's committed counts).
- ``decide`` — ``adni_gain_specificity_decision.json`` (rule version 1) in
  ``<out>``: the verdict ``gate_failed`` | ``production_excess`` |
  ``size_excess`` | ``power_loss`` | ``specificity_excess`` | ``no_gain`` |
  ``reopened`` (designs D12, D13), what it licenses (D14), production's
  statement (D15) and the descriptive readings (D16).

Records hold seeds and statistics only: no participant identifier, label
vector, subsample membership, fold assignment, configuration or shift vector.
The export is read in place on the machine that holds it, never copied; the
sbatch array runs both profiles, the plants profile with ``GAIN_EXPORT``. With
``--workers`` above 1, pin BLAS in the shell.

Examples (the smoke: four replicates per cell):

    export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
    EXPORT=../adni_sexdiff/results/motco_input
    OUT=$HOME/motco-runs/adni-gain-smoke
    python scripts/adni_gain_specificity_check.py preflight --profile plants --export-dir $EXPORT --out-dir $OUT --control-reps 4 --plant-reps 4
    python scripts/adni_gain_specificity_check.py run --profile plants --export-dir $EXPORT --out-dir $OUT --control-reps 4 --plant-reps 4 --workers 7
    python scripts/adni_gain_specificity_check.py merge --profile plants --export-dir $EXPORT --out-dir $OUT
    python scripts/adni_gain_specificity_check.py preflight --profile simulation --out-dir $OUT --control-reps 4 --plant-reps 4
    python scripts/adni_gain_specificity_check.py run --profile simulation --out-dir $OUT --control-reps 4 --plant-reps 4
    python scripts/adni_gain_specificity_check.py merge --profile simulation --out-dir $OUT
    python scripts/adni_gain_specificity_check.py decide --export-dir $EXPORT --out-dir $OUT
"""  # noqa: E501

from __future__ import annotations

import argparse
import hashlib
import json
import math
import sys
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path
from time import perf_counter
from typing import Any, Callable, Iterable, Mapping, NamedTuple, Sequence

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
import adni_anchor_check as aac  # noqa: E402  (contract, units, seeds, arm parameters, checks)
import adni_design_check as adc  # noqa: E402  (the ADNI-design construction, its seeds and unit builder)
import latent_space_comparison as lsc  # noqa: E402  (settings hashing and files, JSON and CSV helpers)
import predictive_subspace_check as psc  # noqa: E402  (arm records, the paired judgements)
import production_latent_space_check as plsc  # noqa: E402  (parents, the exact paired tests, the gain rule)
import stage_span_paper_grade_decision as grade  # noqa: E402  (the count-plus-ceiling rule)

sys.path.pop(0)

DRIVER_NAME = "adni-gain-specificity"
DRIVER_VERSION = 1
RULE_VERSION = 1
#: This driver's own seed family (``direction``, ``reference_split``) and the simulation's base seed.
BASE_SEED = 2040
OWN_PURPOSES: tuple[str, ...] = ("direction", "reference_split")
PROFILES: tuple[str, ...] = ("plants", "simulation")
PLANTS, SIMULATION = PROFILES

#: Subsampled design points of the plants profile (participants); the full cohort is the last point.
POINT_SIZES: tuple[int, ...] = (750,)
COHORT_POINT = aac.COHORT_POINT
CONTROL_REPLICATES = 500
PLANT_REPLICATES = 200
PERMUTATIONS = 999
PREFLIGHT_REPLICATES: tuple[int, ...] = aac.PREFLIGHT_REPLICATES
CV: dict[str, Any] = dict(aac.CV)
#: The simulated table (``adni_design_check.TABLES``) and an optional CV override (``None``: the parent's).
SIMULATION_TABLE = "adni"
SIMULATION_CV: dict[str, Any] | None = None

STATISTICS: tuple[str, ...] = aac.STATISTICS
ALPHA = aac.ALPHA
REJECTION = aac.REJECTION
CONSISTENCY_TOLERANCE = aac.CONSISTENCY_TOLERANCE
#: The ``plant`` check's tolerance, in pooled SDs of the unplanted sample.
PLANT_TOLERANCE = 1e-8
POWER_FLOOR = 0.80

SCORES, PREDICTIVE = aac.SCORES, aac.PREDICTIVE
WITHIN_STAGE = aac.WITHIN_STAGE
REFERENCE_ARM, CANDIDATE_ARM = aac.REFERENCE_ARM, aac.CANDIDATE_ARM
ARMS: tuple[str, ...] = aac.OPERATING_ARMS
WITHIN_STAGE_ARMS: tuple[str, ...] = aac.WITHIN_STAGE_OPERATING_ARMS
#: (arm, reference arm, comparison) of the paired table.
PAIRS: tuple[tuple[str, str, str], ...] = (
    (CANDIDATE_ARM, REFERENCE_ARM, "coordinates"),
    (f"{CANDIDATE_ARM}{WITHIN_STAGE}", f"{REFERENCE_ARM}{WITHIN_STAGE}", "coordinates"),
    (f"{REFERENCE_ARM}{WITHIN_STAGE}", REFERENCE_ARM, "draws"),
    (f"{CANDIDATE_ARM}{WITHIN_STAGE}", CANDIDATE_ARM, "draws"),
)
CHECKS: dict[str, tuple[str, ...]] = {
    PLANTS: (*aac.CONSISTENCY_CHECKS, "plant"),
    SIMULATION: (*aac.CONSISTENCY_CHECKS, "construction"),
}

# ── Cells (designs D6, D7, D9) ───────────────────────────────────────────────

ANCHOR = "anchor"
SAMPLE, HELDOUT, POPULATION = "sample", "heldout", "population"
MAGNITUDE, ORIENTATION, SHAPE, TRANSLATION = "magnitude", "orientation", "shape", "translation"
OUT_OF_SPAN = "out_of_span"
#: The out-of-span rotation's tilt: the plane leaves the stage span entirely.
OUT_OF_SPAN_ALIGNMENT = 90.0
TARGETS: dict[str, str] = {MAGNITUDE: "delta", ORIENTATION: "angle", SHAPE: "shape", OUT_OF_SPAN: "angle"}
TRANSLATION_EFFECT = 1.0
SAMPLE_EFFECTS: dict[str, tuple[float, ...]] = {MAGNITUDE: (0.25, 0.5, 1.0), ORIENTATION: (0.25, 0.5, 1.0),
                                                SHAPE: (0.5, 1.0)}
HELDOUT_EFFECTS: dict[str, tuple[float, ...]] = {MAGNITUDE: (0.5, 1.0, 2.0), ORIENTATION: (0.5, 1.0), SHAPE: (1.0,)}
OUT_OF_SPAN_EFFECTS: tuple[float, ...] = (0.25, 0.5)
SIMULATION_EFFECTS: dict[str, tuple[float, ...]] = {MAGNITUDE: (0.5, 1.0, 2.0), ORIENTATION: (0.25, 0.5, 1.0),
                                                    SHAPE: (1.0, 2.0)}
CONTROL, PLANT, RECORDED = "control", "plant", "out_of_span"
NOUNS: dict[str, str] = {MAGNITUDE: "stretch", ORIENTATION: "rotation", SHAPE: "bend"}

REPO = Path(__file__).resolve().parents[1]
#: The ADNI anchor run (item 21) and the ADNI-design run (item 9), bound by sha256.
ANCHOR_DIR = REPO / "results" / "adni-real-noise-anchor-2026-11-17"
DESIGN_DIR = REPO / "results" / "adni-design-2026-11-02" / "run"
PARENT_FILES: dict[str, tuple[str, str]] = {
    "anchor_rates": ("anchor", "rate_table.csv"),
    "anchor_decision": ("anchor", aac.DECISION_NAME),
    "design_settings": ("design", "settings.json"),
    "design_pvalues": ("design", "adni_design_pvalues.csv"),
}

PREFLIGHT_NAME = "preflight.json"
REPRODUCTION_NAME = "reproduction.json"
DECISION_NAME = "adni_gain_specificity_decision.json"


def cell_name(reference: str, mode: str, effect: float) -> str:
    return f"{reference}_{mode}_e{float(effect):g}"


def _cell(name: str, kind: str, points: Sequence[str], *, reference: str | None = None, mode: str | None = None,
          effect: float | None = None) -> dict[str, Any]:
    return {"name": name, "kind": kind, "reference": reference, "mode": mode,
            "effect": None if effect is None else float(effect),
            "target": TARGETS.get(str(mode)) if kind != CONTROL else None, "points": list(points)}


def plants_cells(points: Sequence[Mapping[str, Any]]) -> list[dict[str, Any]]:
    """Design D6's cells, each with the points it runs at, in enumeration order."""

    every = [str(point["name"]) for point in points]
    subsampled = [str(point["name"]) for point in points if point["subsampled"]]
    cells = [_cell(ANCHOR, CONTROL, every),
             _cell(cell_name(SAMPLE, TRANSLATION, TRANSLATION_EFFECT), CONTROL, every, reference=SAMPLE,
                   mode=TRANSLATION, effect=TRANSLATION_EFFECT)]
    for mode, effects in SAMPLE_EFFECTS.items():
        cells += [_cell(cell_name(SAMPLE, mode, e), PLANT, every, reference=SAMPLE, mode=mode, effect=e)
                  for e in effects]
    for mode, effects in HELDOUT_EFFECTS.items():
        cells += [_cell(cell_name(HELDOUT, mode, e), PLANT, subsampled, reference=HELDOUT, mode=mode, effect=e)
                  for e in effects]
    cells += [_cell(cell_name(SAMPLE, OUT_OF_SPAN, e), RECORDED, every, reference=SAMPLE, mode=OUT_OF_SPAN, effect=e)
              for e in OUT_OF_SPAN_EFFECTS]
    return cells


def simulation_cells(point: str) -> list[dict[str, Any]]:
    """Design D9's cells at the simulated table."""

    cells = [_cell(ANCHOR, CONTROL, [point]),
             _cell(cell_name(POPULATION, TRANSLATION, TRANSLATION_EFFECT), CONTROL, [point], reference=POPULATION,
                   mode=TRANSLATION, effect=TRANSLATION_EFFECT)]
    for mode, effects in SIMULATION_EFFECTS.items():
        cells += [_cell(cell_name(POPULATION, mode, e), PLANT, [point], reference=POPULATION, mode=mode, effect=e)
                  for e in effects]
    return cells


def cell_spec(settings: Mapping[str, Any], name: str) -> dict[str, Any]:
    for cell in settings["cells"]:
        if cell["name"] == name:
            return dict(cell)
    raise KeyError(f"Unknown cell {name!r}.")


def cell_replicates(settings: Mapping[str, Any], cell: Mapping[str, Any]) -> int:
    return int(settings["replicates"][CONTROL if cell["kind"] == CONTROL else PLANT])


def point_names(settings: Mapping[str, Any]) -> list[str]:
    return aac.point_names(settings)


def point_spec(settings: Mapping[str, Any], name: str) -> dict[str, Any]:
    return aac.point_spec(settings, name)


def subsampled_points(settings: Mapping[str, Any]) -> list[str]:
    return [str(point["name"]) for point in settings["points"] if point["subsampled"]]


def point_cells(settings: Mapping[str, Any], point: str) -> list[dict[str, Any]]:
    return [cell for cell in settings["cells"] if point in cell["points"]]


# ── Arms and settings ────────────────────────────────────────────────────────


def arm_specs() -> dict[str, dict[str, Any]]:
    """The four operating arms of the anchor driver (design D8), in evaluation order."""

    every = aac.arm_specs(())
    return {arm: dict(every[arm]) for arm in ARMS}


def _seed(*parts: Any) -> int:
    encoded = json.dumps([DRIVER_NAME, *parts], separators=(",", ":")).encode("utf-8")
    return int(hashlib.sha256(encoded).hexdigest()[:8], 16) & 0x7FFFFFFF


def own_seed(settings: Mapping[str, Any], point: str, replicate: int, purpose: str) -> int:
    if purpose not in OWN_PURPOSES:
        raise ValueError(f"Unknown purpose {purpose!r}; this driver's own purposes are {list(OWN_PURPOSES)}.")
    return _seed(int(settings["own_base_seed"]), point, int(replicate), purpose)


def parents(anchor_dir: Path | None = None, design_dir: Path | None = None) -> dict[str, dict[str, str]]:
    directories = {"anchor": Path(ANCHOR_DIR if anchor_dir is None else anchor_dir),
                   "design": Path(DESIGN_DIR if design_dir is None else design_dir)}
    return {role: plsc._parent(directories[source] / name) for role, (source, name) in PARENT_FILES.items()}


def _replicates(control: int | None, plant: int | None) -> dict[str, int]:
    control = CONTROL_REPLICATES if control is None else int(control)
    plant = PLANT_REPLICATES if plant is None else int(plant)
    if control < max(PREFLIGHT_REPLICATES) + 1 or plant < 1:
        raise SystemExit(f"--control-reps must be at least {max(PREFLIGHT_REPLICATES) + 1} and --plant-reps at "
                         "least 1.")
    return {CONTROL: control, PLANT: plant}


def _analysis(cv: Mapping[str, Any], min_components: int, orientation_estimator: str) -> dict[str, Any]:
    return {
        "cv": dict(cv),
        "min_components": int(min_components),
        "orientation_estimator": str(orientation_estimator),
        "residual_correction": aac.RESIDUAL_CORRECTION,
        "missing_covariates": aac.MISSING_COVARIATES,
        "permutation_strata": aac.STRATA,
        "latent_coordinates": PREDICTIVE,
        "n_jobs": 1,
    }


def plants_settings(
    cohort: aac.Cohort,
    *,
    control_replicates: int | None = None,
    plant_replicates: int | None = None,
    permutations: int | None = None,
    point_sizes: Sequence[int] | None = None,
    cv: Mapping[str, Any] | None = None,
    anchor_dir: Path | None = None,
    design_dir: Path | None = None,
) -> dict[str, Any]:
    """The plants profile's settings; unset arguments take the module's values at call time."""

    contract = cohort.contract
    if contract.transform != "identity":
        raise SystemExit(f"The plant maps a shift back through the pooled SD; the contract's transform is "
                         f"{contract.transform!r}, not 'identity'.")
    permutations = PERMUTATIONS if permutations is None else int(permutations)
    if permutations < 1:
        raise SystemExit("--perms must be at least 1.")
    points = aac.design_points(cohort.table(), POINT_SIZES if point_sizes is None else tuple(point_sizes))
    specs = arm_specs()
    return {
        "driver": DRIVER_NAME,
        "driver_version": DRIVER_VERSION,
        "rule_version": RULE_VERSION,
        "profile": PLANTS,
        "contract": contract.to_dict(),
        "export": {"files": dict(cohort.files), "participants": int(len(cohort.metadata)),
                   "features": {name: int(frame.shape[1]) for name, frame in cohort.blocks.items()}},
        "parents": parents(anchor_dir, design_dir),
        "base_seed": aac.BASE_SEED,
        "seed_derivation": f"the anchor driver's: sha256 of [{aac.DRIVER_NAME!r}, {aac.BASE_SEED}, point, replicate, "
                           "purpose(, attempt)]; purposes subsample, labels, rrpp, folds; shared by every cell of a "
                           "point at a replicate",
        "own_base_seed": BASE_SEED,
        "own_seed_derivation": f"sha256 of [{DRIVER_NAME!r}, {BASE_SEED}, point, replicate, purpose]; purposes "
                               "direction, reference_split",
        "max_attempts": aac.MAX_ATTEMPTS,
        "points": points,
        "cells": plants_cells(points),
        "replicates": _replicates(control_replicates, plant_replicates),
        "permutations": permutations,
        "cohort_rank": int(contract.fidelity_rank),
        "arms": list(specs),
        "arm_specs": specs,
        "analysis": _analysis(CV if cv is None else cv, aac.MIN_COMPONENTS, aac.ORIENTATION_ESTIMATOR),
        "plant": {
            "planted_group": contract.groups[1],
            "configuration": "stage coefficients of features ~ stage indicators + covariate columns (coded as "
                             "analyze_trajectories codes them, centred on the reference's means) in the analysis "
                             "sample's pooled-standardized coordinates; the group label is not in the model",
            "edit": "motco.simulations.semisynthetic.edit_stage_configuration(z, mode, effect); the shift is "
                    "edited - z, mapped back through the unplanted sample's pooled SD and added to every row of the "
                    "planted group at each stage under the unit's permuted labels",
            "out_of_span": {"alignment_degrees": OUT_OF_SPAN_ALIGNMENT,
                            "direction": "standard normal over the joint features (purpose direction), projected off "
                                         "the sample configuration's centred stage means and normalized"},
            "tolerance": PLANT_TOLERANCE,
        },
        "fit": "pls_cv cross-validates on the unit's data and every other arm is forced at its rank; a cohort anchor "
               "unit forces every arm at cohort_rank",
        "alpha": ALPHA,
        "rejection": REJECTION,
        "consistency_tolerance": CONSISTENCY_TOLERANCE,
        "preflight_replicates": list(PREFLIGHT_REPLICATES),
    }


def simulation_point(table: str) -> str:
    return f"table_{table}"


def construction_settings(
    *, table: str, replicates: int, permutations: int, cv: Mapping[str, Any] | None = None
) -> dict[str, Any]:
    """``adni_design_check``'s settings at its calibrated strengths and the table, under the stage-span plant."""

    args = adc.build_parser().parse_args(["run", "--out-dir", "unused", "--table", str(table), "--reps",
                                          str(int(replicates)), "--perms", str(int(permutations)), "--base-seed",
                                          str(BASE_SEED)])
    construction = adc.build_settings(args)
    if cv is not None:
        construction["cv"] = dict(cv)
    construction["generator"].update({"construction": "stage_span", "span_layers": list(construction["layers"]),
                                      "span_scale": "planted"})
    return construction


def simulation_settings(
    *,
    control_replicates: int | None = None,
    plant_replicates: int | None = None,
    permutations: int | None = None,
    table: str | None = None,
    cv: Mapping[str, Any] | None = None,
    anchor_dir: Path | None = None,
    design_dir: Path | None = None,
) -> dict[str, Any]:
    """The simulation profile's settings; unset arguments take the module's values at call time."""

    replicates = _replicates(control_replicates, plant_replicates)
    permutations = PERMUTATIONS if permutations is None else int(permutations)
    if permutations < 1:
        raise SystemExit("--perms must be at least 1.")
    table = SIMULATION_TABLE if table is None else str(table)
    construction = construction_settings(table=table, replicates=max(replicates.values()), permutations=permutations,
                                         cv=SIMULATION_CV if cv is None else cv)
    point = simulation_point(table)
    sizes = construction["generator"]["group_stage_sizes"]
    specs = arm_specs()
    return {
        "driver": DRIVER_NAME,
        "driver_version": DRIVER_VERSION,
        "rule_version": RULE_VERSION,
        "profile": SIMULATION,
        "parents": parents(anchor_dir, design_dir),
        "base_seed": BASE_SEED,
        "seed_derivation": "adni_design_check's, at this base seed: every seed derives from the base seed and the "
                           "replicate, so every cell of a replicate shares the baseline, the covariates, the phase "
                           "draw and the RRPP seed",
        "construction": construction,
        "points": [{"name": point, "participants": int(sum(map(sum, sizes))),
                    "table": [list(map(int, row)) for row in sizes], "subsampled": False}],
        "cells": simulation_cells(point),
        "replicates": replicates,
        "permutations": permutations,
        "arms": list(specs),
        "arm_specs": specs,
        "analysis": _analysis(construction["cv"], int(construction["min_components"]),
                              str(construction["orientation_estimator"])),
        "fit": "pls_cv cross-validates on the unit's data and every other arm is forced at its rank",
        "alpha": ALPHA,
        "rejection": REJECTION,
        "consistency_tolerance": CONSISTENCY_TOLERANCE,
        "preflight_replicates": list(PREFLIGHT_REPLICATES),
    }


class Unit(NamedTuple):
    point: str
    cell: str
    replicate: int


def unit_key(point: str, cell: str, replicate: int) -> str:
    return f"{point}/{cell}/{int(replicate)}"


def record_key(record: Mapping[str, Any]) -> str:
    return unit_key(record["point"], record["cell"], record["replicate"])


def enumerate_units(settings: Mapping[str, Any]) -> list[Unit]:
    """Replicate-major, so a partial run covers every point and cell."""

    most = max(int(value) for value in settings["replicates"].values())
    return [
        Unit(point, cell["name"], replicate)
        for replicate in range(most)
        for point in point_names(settings)
        for cell in point_cells(settings, point)
        if replicate < cell_replicates(settings, cell)
    ]


def expected_units(settings: Mapping[str, Any]) -> dict[tuple[str, str], int]:
    return {(point, cell["name"]): cell_replicates(settings, cell) for point in point_names(settings)
            for cell in point_cells(settings, point)}


# ── The configuration and the plant (designs D3, D4, D5, D7) ────────────────


def _transforms(contract: aac.CohortContract) -> dict[str, str]:
    return {name: contract.transform for name in contract.block_names}


def sample_scalers(blocks: Mapping[str, pd.DataFrame], contract: aac.CohortContract) -> dict[str, Any]:
    """The analysis's pooled scalers of these blocks (``fit_block_scalers``, unweighted)."""

    from motco.simulations.preprocessing import fit_block_scalers

    return fit_block_scalers({name: blocks[name] for name in contract.block_names}, _transforms(contract))


def standardized(blocks: Mapping[str, pd.DataFrame], scalers: Mapping[str, Any],
                 contract: aac.CohortContract) -> np.ndarray:
    """Rows in the scalers' pooled-standardized joint coordinates, blocks in declared order."""

    from motco.simulations.preprocessing import concatenate_blocks, transform_blocks

    names = contract.block_names
    joined = transform_blocks({name: blocks[name] for name in names}, scalers, _transforms(contract))
    return concatenate_blocks(joined, names).to_numpy(dtype=float)


def joint_scale(scalers: Mapping[str, Any], contract: aac.CohortContract) -> np.ndarray:
    """Raw units per standardized unit, per joint feature."""

    return np.concatenate([np.asarray(scalers[name].scale, dtype=float) / float(scalers[name].weight)
                           for name in contract.block_names])


def stage_positions(metadata: pd.DataFrame, contract: aac.CohortContract) -> np.ndarray:
    """Each row's stage position in the declared order (never sorted)."""

    order = {stage: index for index, stage in enumerate(contract.stage_order)}
    return np.asarray([order[str(stage)] for stage in metadata[contract.stage_col]], dtype=int)


def covariate_columns(metadata: pd.DataFrame, contract: aac.CohortContract) -> np.ndarray:
    """The covariate design columns as ``analyze_trajectories`` codes them (``build_covariate_matrix``)."""

    from motco.stats.design import build_covariate_matrix

    return build_covariate_matrix(metadata.reset_index(drop=True), contract.covariate_names,
                                  contract.categorical_covariates).matrix


def adjusted_configuration(
    features: np.ndarray, stages: np.ndarray, covariates: np.ndarray, n_stages: int, *,
    reference_row: np.ndarray | None = None,
) -> np.ndarray:
    """The stage coefficients (stages × features) of ``features ~ stage indicators + (covariates − reference_row)``.

    ``reference_row`` defaults to the covariates' mean. A design without full
    column rank raises ``ValueError``.
    """

    features = np.asarray(features, dtype=float)
    covariates = np.asarray(covariates, dtype=float).reshape(features.shape[0], -1)
    stages = np.asarray(stages, dtype=int)
    indicators = (stages[:, None] == np.arange(int(n_stages))[None, :]).astype(float)
    row = covariates.mean(axis=0) if reference_row is None else np.asarray(reference_row, dtype=float)
    design = np.hstack([indicators, covariates - row])
    rank = int(np.linalg.matrix_rank(design))
    if rank < design.shape[1]:
        raise ValueError(f"The covariate-adjusted stage fit is rank-deficient: rank {rank} of {design.shape[1]} "
                         f"columns on {design.shape[0]} participants.")
    coefficients, *_ = np.linalg.lstsq(design, features, rcond=None)
    return coefficients[: int(n_stages)]


def heldout_rows(cohort: aac.Cohort, rows: np.ndarray) -> np.ndarray:
    """The cohort's rows outside the sample, in the cohort's row order."""

    return np.setdiff1d(np.arange(len(cohort.metadata)), np.asarray(rows, dtype=int))


def reference_configuration(
    cohort: aac.Cohort, built: Mapping[str, Any], scalers: Mapping[str, Any], reference: str
) -> np.ndarray:
    """The adjusted configuration of ``sample`` or ``heldout`` in the sample's pooled-standardized coordinates."""

    contract = cohort.contract
    if reference == SAMPLE:
        blocks, metadata = built["blocks"], built["metadata"]
    elif reference == HELDOUT:
        rows = heldout_rows(cohort, built["rows"])
        if rows.size == 0:
            raise ValueError("The held-out reference is empty: the point holds every participant.")
        blocks = {name: frame.iloc[rows] for name, frame in cohort.blocks.items()}
        metadata = cohort.metadata.iloc[rows].reset_index(drop=True)
    else:
        raise ValueError(f"Unknown reference {reference!r}.")
    return adjusted_configuration(standardized(blocks, scalers, contract), stage_positions(metadata, contract),
                                  covariate_columns(metadata, contract), len(contract.stage_order))


def out_of_span_vector(z: np.ndarray, seed: int) -> np.ndarray:
    """A seeded unit direction orthogonal to every centred stage mean of ``z`` (stages × features)."""

    z = np.asarray(z, dtype=float)
    centred = z - z.mean(axis=0)
    left, singular, right = np.linalg.svd(centred, full_matrices=False)
    keep = singular > 1e-12 * max(float(singular[0]) if singular.size else 0.0, np.finfo(float).tiny)
    basis = right[keep].T
    vector = np.random.default_rng(int(seed)).standard_normal(z.shape[1])
    for _ in range(2):  # project twice: orthogonal to round-off
        vector = vector - basis @ (basis.T @ vector)
    return vector / np.linalg.norm(vector)


def plant_shift(
    z: np.ndarray, mode: str, effect: float, *, direction: np.ndarray | None = None
) -> tuple[np.ndarray, dict[str, Any]]:
    """The edit of ``z`` minus ``z`` (stages × features) and the edit's parameters."""

    from motco.simulations.semisynthetic import edit_stage_configuration

    z = np.asarray(z, dtype=float)
    if mode == OUT_OF_SPAN:
        if direction is None:
            raise ValueError("The out-of-span rotation needs its direction.")
        edited, details = edit_stage_configuration(z, ORIENTATION, float(effect),
                                                   alignment_degrees=OUT_OF_SPAN_ALIGNMENT, out_of_span=direction)
    else:
        edited, details = edit_stage_configuration(z, mode, float(effect))
    return edited - z, details


def apply_plant(
    blocks: Mapping[str, pd.DataFrame], metadata: pd.DataFrame, contract: aac.CohortContract, shift_raw: np.ndarray
) -> dict[str, pd.DataFrame]:
    """The blocks with ``shift_raw[s]`` (raw units) added to every planted-group row at stage ``s``.

    A zero shift returns the frames untouched (the same objects).
    """

    shift_raw = np.asarray(shift_raw, dtype=float)
    if not np.any(shift_raw):
        return dict(blocks)
    planted_rows = metadata[contract.group_col].astype(str).to_numpy() == str(contract.groups[1])
    stages = stage_positions(metadata, contract)
    out: dict[str, pd.DataFrame] = {}
    start = 0
    for name in contract.block_names:
        frame = blocks[name]
        width = frame.shape[1]
        values = frame.to_numpy(dtype=float).copy()
        for stage in range(shift_raw.shape[0]):
            values[planted_rows & (stages == stage)] += shift_raw[stage, start:start + width]
        out[name] = pd.DataFrame(values, index=frame.index, columns=frame.columns)
        start += width
    return out


def check_plant(
    original: Mapping[str, pd.DataFrame], planted: Mapping[str, pd.DataFrame], metadata: pd.DataFrame,
    contract: aac.CohortContract, shift_raw: np.ndarray | None, scale: np.ndarray,
) -> dict[str, Any]:
    """The second group's raw stage means moved by the declared shift and the first group's did not (pooled SDs)."""

    if shift_raw is None:
        return {**psc._check(None, PLANT_TOLERANCE, applicable=False), "failing": {}}
    labels = metadata[contract.group_col].astype(str).to_numpy()
    stages = stage_positions(metadata, contract)
    before = np.hstack([original[name].to_numpy(dtype=float) for name in contract.block_names])
    after = np.hstack([planted[name].to_numpy(dtype=float) for name in contract.block_names])
    deviations = {"first_group": 0.0, "second_group": 0.0}
    for stage in range(np.asarray(shift_raw).shape[0]):
        for key, group, expected in (("first_group", contract.groups[0], np.zeros(before.shape[1])),
                                     ("second_group", contract.groups[1], np.asarray(shift_raw)[stage])):
            rows = (labels == str(group)) & (stages == stage)
            if not rows.any():
                deviations[key] = math.inf
                continue
            moved = after[rows].mean(axis=0) - before[rows].mean(axis=0)
            deviations[key] = max(deviations[key], float(np.max(np.abs((moved - expected) / scale))))
    worst = max(deviations.values())
    check = psc._check(worst, PLANT_TOLERANCE)
    return {**check, "failing": {key: value for key, value in deviations.items()
                                 if not (math.isfinite(value) and value <= PLANT_TOLERANCE)}}


def model_inputs(displacements_raw: np.ndarray, scale: np.ndarray, model: Any) -> np.ndarray:
    """Raw displacements in a fitted model's input coordinates (analysis scaler, then the model's own scaling)."""

    return np.asarray(displacements_raw, dtype=float) / scale / np.asarray(model._x_std, dtype=float)


def capture_shares(vectors: np.ndarray, model: Any) -> dict[str, float | None]:
    """Shares of ``Σ‖v‖²`` inside the fitted span and inside the predictive subspace (model input coordinates)."""

    from motco.simulations.diagnostics import LinearLatentMap, span_share

    vectors = np.atleast_2d(np.asarray(vectors, dtype=float))
    if not float(np.sum(vectors**2)) > 0:
        return {"span": None, "predictive": None}
    return {"span": span_share(vectors, LinearLatentMap.from_pls(model).basis),
            "predictive": span_share(vectors, LinearLatentMap.from_pls_predictive(model, PREDICTIVE).basis)}


def split_within_stage(stages: np.ndarray, seed: int, n_stages: int) -> tuple[np.ndarray, np.ndarray]:
    """Two halves of the rows, each stage split by one seeded permutation (the first half takes the floor)."""

    rng = np.random.default_rng(int(seed))
    first, second = [], []
    for stage in range(int(n_stages)):
        rows = np.flatnonzero(np.asarray(stages) == stage)
        order = rng.permutation(rows.size)
        half = rows.size // 2
        first.append(rows[order[:half]])
        second.append(rows[order[half:]])
    return np.sort(np.concatenate(first)), np.sort(np.concatenate(second))


def reliability(
    features: np.ndarray, stages: np.ndarray, metadata: pd.DataFrame, contract: aac.CohortContract, seed: int
) -> dict[str, Any]:
    """The centred adjusted configuration's split-half agreement: correlation and the cross-product's eigenvalues.

    The cross-product is taken in stage space, ``(A·Bᵀ + B·Aᵀ) / 2`` (stages ×
    stages) for the halves' centred configurations ``A`` and ``B``: its
    expectation is the true configuration's Gram matrix, whose nonzero
    eigenvalues are the configuration's squared singular values. (The
    feature-space product ``AᵀB`` has the same expectation but carries ± noise
    eigenvalue pairs of the order of ‖A‖·‖B‖.)
    """

    n_stages = len(contract.stage_order)
    halves = split_within_stage(stages, seed, n_stages)
    centred = []
    for rows in halves:
        part = metadata.iloc[rows].reset_index(drop=True)
        z = adjusted_configuration(features[rows], stages[rows], covariate_columns(part, contract), n_stages)
        centred.append(z - z.mean(axis=0))
    a, b = centred
    norm = float(np.linalg.norm(a) * np.linalg.norm(b))
    cross = (a @ b.T + b @ a.T) / 2.0
    eigenvalues = np.sort(np.linalg.eigvalsh(cross))[::-1]
    first, second = float(eigenvalues[0]), float(eigenvalues[1])
    correlation = float(np.sum(a * b) / norm) if norm > 0 else None
    return {
        "seed": int(seed),
        "correlation": correlation,
        # Spearman–Brown: the full reference's reliability from the split-half correlation.
        "full_reference": (2.0 * correlation / (1.0 + correlation)
                           if correlation is not None and correlation > -1.0 else None),
        "eigenvalues": [first, second],
        "eigenvalue_ratio": max(second, 0.0) / first if first > 0 else None,
    }


def build_plant(
    cohort: aac.Cohort, settings: Mapping[str, Any], built: Mapping[str, Any], cell: Mapping[str, Any]
) -> dict[str, Any]:
    """The unit's planted blocks and, in memory only, the shift and the configurations the readings need."""

    contract = cohort.contract
    scalers = sample_scalers(built["blocks"], contract)
    scale = joint_scale(scalers, contract)
    configurations: dict[str, np.ndarray] = {}
    if cell["kind"] == CONTROL and cell["mode"] is None:
        references = [SAMPLE, *([HELDOUT] if point_spec(settings, built["point"])["subsampled"] else [])]
        for reference in references:
            configurations[reference] = reference_configuration(cohort, built, scalers, reference)
        return {"blocks": dict(built["blocks"]), "shift_std": None, "shift_raw": None, "scale": scale,
                "configurations": configurations, "edit": None, "direction_seed": None}
    reference = str(cell["reference"])
    z = reference_configuration(cohort, built, scalers, reference)
    configurations[reference] = z
    direction_seed = None
    direction = None
    if cell["mode"] == OUT_OF_SPAN:
        direction_seed = own_seed(settings, built["point"], built["replicate"], "direction")
        direction = out_of_span_vector(z, direction_seed)
    shift_std, details = plant_shift(z, str(cell["mode"]), float(cell["effect"]), direction=direction)
    shift_raw = shift_std * scale
    return {"blocks": apply_plant(built["blocks"], built["metadata"], contract, shift_raw), "shift_std": shift_std,
            "shift_raw": shift_raw, "scale": scale, "configurations": configurations, "edit": details,
            "direction_seed": direction_seed}


def plant_summary(planted: Mapping[str, Any], contract: aac.CohortContract) -> dict[str, Any] | None:
    """The plant's record-safe scalars: norm per stage (pooled SDs), largest entry, the scale drift."""

    shift = planted["shift_std"]
    if shift is None:
        return None
    after = joint_scale(sample_scalers(planted["blocks"], contract), contract)
    return {
        "norm_per_stage": [float(value) for value in np.linalg.norm(shift, axis=1)],
        "max_abs_entry": float(np.max(np.abs(shift))),
        "scale_drift": float(np.max(np.abs(after / planted["scale"] - 1.0))),
        "edit": {key: value for key, value in (planted["edit"] or {}).items()},
    }


# ── Arms on one fit (design D8) ──────────────────────────────────────────────


def _analyse(dataset: Any, params: Any) -> Any:
    return aac._analyse(dataset, params)


def arm_params(arm: str, settings: Mapping[str, Any], seed: int, *, rank: int | None) -> Any:
    """``rank`` None cross-validates (``pls_cv`` only); otherwise the arm is forced at ``rank``."""

    spec = settings["arm_specs"][arm]
    return aac.analysis_params(settings, seed, forced=rank, coordinates=spec["coordinates"], draws=spec["draws"],
                               correction=spec["residual_correction"])


def evaluate_arms(
    dataset: Any, settings: Mapping[str, Any], seed: int, *, cross_validate: bool, forced_rank: int | None,
    timing: dict[str, float],
) -> tuple[Any, dict[str, Any], dict[str, Any], dict[str, int]]:
    """Production first (raises to the caller), then every other arm forced at its rank (a failure is recorded)."""

    specs = settings["arm_specs"]
    mark = perf_counter()
    production = _analyse(dataset, arm_params(REFERENCE_ARM, settings, seed,
                                              rank=None if cross_validate else int(forced_rank or 0)))
    timing["production"] = perf_counter() - mark
    rank = int(production.record["latent"]["selected_lv"])
    fit = "cross_validated" if cross_validate else "forced_at_cohort_rank"
    arms = {REFERENCE_ARM: aac.arm_record(REFERENCE_ARM, production, timing["production"], specs[REFERENCE_ARM],
                                          fit)}
    results = {REFERENCE_ARM: production}
    forced: dict[str, int] = {} if cross_validate else {REFERENCE_ARM: rank}
    for arm in settings["arms"]:
        if arm == REFERENCE_ARM:
            continue
        forced[arm] = rank
        mark = perf_counter()
        try:
            result = _analyse(dataset, arm_params(arm, settings, seed, rank=rank))
            arms[arm] = aac.arm_record(arm, result, perf_counter() - mark, specs[arm],
                                       "forced_at_cv_rank" if cross_validate else "forced_at_cohort_rank")
            results[arm] = result
        except Exception as exc:  # recorded on the arm, the unit goes on
            arms[arm] = aac.failed_arm(arm, exc, perf_counter() - mark, specs[arm])
    return production, arms, results, forced


def recorded_table(result: Any, groups: Sequence[str], stages: Sequence[str]) -> list[list[int]]:
    sizes = result.record["design"]["cell_sizes"]
    return [[int(sizes[str(group)][str(stage)]) for stage in stages] for group in groups]


def common_checks(
    *, arms: Mapping[str, Any], results: Mapping[str, Any], forced: Mapping[str, int], settings: Mapping[str, Any],
    n_stages: int, expected_table: Sequence[Sequence[int]], realized: Mapping[str, Sequence[Sequence[int]]],
) -> dict[str, dict[str, Any]]:
    """The anchor driver's five checks, imported."""

    order = list(settings["arms"])
    return {
        "rank": aac.check_rank(arms, forced),
        "same_fit": aac.check_same_fit({arm: result.latent_artifacts.model for arm, result in results.items()}, order),
        "measured_dimension": aac.check_measured_dimension(arms, settings["arm_specs"], n_stages),
        "same_statistics": aac.check_same_statistics(arms),
        "cell_table": aac.check_cell_table(expected_table, realized),
    }


def _finish(base: dict[str, Any], checks: Mapping[str, Mapping[str, Any]], arms: Mapping[str, Any],
            names: Sequence[str], start: float, timing: dict[str, float], **extra: Any) -> dict[str, Any]:
    failing = [name for name in names if not checks[name]["passed"]]
    timing["unit"] = perf_counter() - start
    return lsc._jsonable({
        **base,
        **extra,
        "status": "inconsistent" if failing else "completed",
        "inconsistent_checks": failing,
        "consistency": dict(checks),
        "arms": dict(arms),
        "runtime_seconds": timing,
    })


def _failed_unit(base: Mapping[str, Any], exc: Exception, start: float, timing: dict[str, float]) -> dict[str, Any]:
    return lsc._jsonable(psc._failed_unit(base, exc, start, timing))


def _production_entry(production: Any, fit: str) -> dict[str, Any]:
    latent = production.record["latent"]
    return {"rank": int(latent["selected_lv"]), "component_selection": latent.get("component_selection"),
            "cv_mean_auroc": latent.get("cv_mean_auroc"), "fit": fit}


def _cell_fields(cell: Mapping[str, Any]) -> dict[str, Any]:
    return {key: cell[key] for key in ("kind", "reference", "mode", "effect", "target")}


# ── One plants unit ──────────────────────────────────────────────────────────


def evaluate_plants_unit(unit: Unit, settings: Mapping[str, Any], cohort: aac.Cohort) -> dict[str, Any]:
    """Build the unit and its plant, fit production once, evaluate every arm; one record of seeds and statistics."""

    start = perf_counter()
    contract = cohort.contract
    point = point_spec(settings, unit.point)
    cell = cell_spec(settings, unit.cell)
    seeds = aac.unit_seeds(settings, unit.point, unit.replicate)
    seed = int(seeds["rrpp"])
    n_stages = len(contract.stage_order)
    cohort_anchor = not point["subsampled"] and cell["name"] == ANCHOR
    base: dict[str, Any] = {
        "profile": PLANTS, "point": unit.point, "cell": unit.cell, "replicate": int(unit.replicate),
        **_cell_fields(cell), "participants": int(point["participants"]), "seeds": dict(seeds), "rrpp_seed": seed,
        "settings_hash": lsc.settings_hash(settings),
    }
    timing: dict[str, float] = {}
    try:
        mark = perf_counter()
        built = aac.build_unit(cohort, settings, unit.point, unit.replicate)
        base.update({"subsample_redraws": built["subsample_redraws"], "label_redraws": built["label_redraws"],
                     "cell_table": built["table"]})
        planted = build_plant(cohort, settings, built, cell)
        if planted["direction_seed"] is not None:
            base["seeds"]["direction"] = planted["direction_seed"]
        dataset = aac.make_dataset(planted["blocks"], built["metadata"], contract)
        timing["build"] = perf_counter() - mark
        production, arms, results, forced = evaluate_arms(
            dataset, settings, seed, cross_validate=not cohort_anchor, forced_rank=int(settings["cohort_rank"]),
            timing=timing)
    except Exception as exc:  # recorded, never aborts the shard
        return _failed_unit(base, exc, start, timing)

    base["production"] = _production_entry(production, "forced_at_cohort_rank" if cohort_anchor else
                                           "cross_validated")
    realized = {"labels": built["table"]}
    realized.update({arm: recorded_table(result, contract.groups, contract.stage_order)
                     for arm, result in results.items()})
    checks = common_checks(arms=arms, results=results, forced=forced, settings=settings, n_stages=n_stages,
                           expected_table=point["table"], realized=realized)
    checks["plant"] = check_plant(built["blocks"], planted["blocks"], built["metadata"], contract,
                                  planted["shift_raw"], planted["scale"])

    mark = perf_counter()
    extra: dict[str, Any] = {}
    try:
        extra["plant"] = plant_summary(planted, contract)
        extra.update(plant_readings(planted, production, contract, cell))
    except Exception as exc:  # readings, never decisive
        extra["readings_error"] = f"{type(exc).__name__}: {exc}"
    if cell["name"] == ANCHOR and point["subsampled"]:
        split_seed = own_seed(settings, unit.point, unit.replicate, "reference_split")
        base["seeds"]["reference_split"] = split_seed
        try:
            extra["reliability"] = heldout_reliability(cohort, built, split_seed)
        except Exception as exc:  # a reading, never decisive
            extra["reliability"] = {"seed": split_seed, "error": f"{type(exc).__name__}: {exc}"}
    timing["readings"] = perf_counter() - mark
    return _finish(base, checks, arms, CHECKS[PLANTS], start, timing, **extra)


def plant_readings(planted: Mapping[str, Any], production: Any, contract: aac.CohortContract,
                   cell: Mapping[str, Any]) -> dict[str, Any]:
    """Capture shares in the fitted span and the predictive subspace: the plant's, and each configuration's."""

    model = production.latent_artifacts.model
    after = joint_scale(sample_scalers(planted["blocks"], contract), contract)
    captures: dict[str, Any] = {}
    if planted["shift_raw"] is not None:
        captures["plant"] = capture_shares(model_inputs(planted["shift_raw"], after, model), model)
    for reference, z in planted["configurations"].items():
        centred = (z - z.mean(axis=0)) * planted["scale"]
        key = "reference" if cell["reference"] == reference else reference
        captures[key] = capture_shares(model_inputs(centred, after, model), model)
    return {"capture": captures}


def heldout_reliability(cohort: aac.Cohort, built: Mapping[str, Any], seed: int) -> dict[str, Any]:
    """The held-out configuration's split-half reliability, in the sample's pooled-standardized coordinates."""

    contract = cohort.contract
    rows = heldout_rows(cohort, built["rows"])
    blocks = {name: frame.iloc[rows] for name, frame in cohort.blocks.items()}
    metadata = cohort.metadata.iloc[rows].reset_index(drop=True)
    features = standardized(blocks, sample_scalers(built["blocks"], contract), contract)
    return reliability(features, stage_positions(metadata, contract), metadata, contract, seed)


# ── One simulated unit ───────────────────────────────────────────────────────


_REFERENCE: Any = None


def _reference() -> Any:
    global _REFERENCE
    if _REFERENCE is None:
        from motco.simulations.reference import load_reference

        _REFERENCE = load_reference()
    return _REFERENCE


def construction_mode(cell: Mapping[str, Any]) -> tuple[str, float]:
    """The parent's (mode, effect) of a cell: the anchor is its own mode at zero effect."""

    if cell["name"] == ANCHOR:
        return adc.ANCHOR, 0.0
    return str(cell["mode"]), float(cell["effect"])


def check_construction(truth: Mapping[str, Any], cell: Mapping[str, Any]) -> dict[str, Any]:
    """The dataset's truth records the stage-span construction under the planted scale at the cell's mode and effect."""

    mode, effect = construction_mode(cell)
    expected = {"construction": "stage_span", "span_scale": "planted",
                "trajectory_mode": "magnitude" if mode == adc.ANCHOR else mode, "group_effect_size": float(effect)}
    recorded = {key: truth.get(key) for key in expected}
    differing = sorted(key for key in expected if recorded[key] != expected[key])
    return {**psc._check(math.inf if differing else 0.0, 0.0), "expected": expected, "recorded": recorded,
            "differing": differing}


def evaluate_simulation_unit(unit: Unit, settings: Mapping[str, Any]) -> dict[str, Any]:
    """Generate the replicate at the cell's stage-span plant, fit production once, evaluate every arm."""

    start = perf_counter()
    construction = settings["construction"]
    point = point_spec(settings, unit.point)
    cell = cell_spec(settings, unit.cell)
    mode, effect = construction_mode(cell)
    base: dict[str, Any] = {
        "profile": SIMULATION, "point": unit.point, "cell": unit.cell, "replicate": int(unit.replicate),
        **_cell_fields(cell), "participants": int(point["participants"]), "settings_hash": lsc.settings_hash(settings),
    }
    timing: dict[str, float] = {}
    try:
        mark = perf_counter()
        built = adc.build_unit(mode, effect, unit.replicate, construction, reference=_reference())
        seed = int(built["seeds"]["rrpp"])
        base.update({"seeds": dict(built["seeds"]), "rrpp_seed": seed})
        dataset = adc.trajectory_dataset(built, construction, covariates=True)
        timing["build"] = perf_counter() - mark
        production, arms, results, forced = evaluate_arms(dataset, settings, seed, cross_validate=True,
                                                          forced_rank=None, timing=timing)
    except Exception as exc:  # recorded, never aborts the shard
        return _failed_unit(base, exc, start, timing)

    base["production"] = _production_entry(production, "cross_validated")
    groups = list(built["group_labels"])
    stages = list(dataset.stages)
    metadata = built["metadata"]
    labels = [[int(((metadata["group"].astype(str) == group) & (metadata["stage"].astype(str) == stage)).sum())
               for stage in stages] for group in groups]
    base["cell_table"] = labels
    realized = {"labels": labels}
    realized.update({arm: recorded_table(result, groups, stages) for arm, result in results.items()})
    checks = common_checks(arms=arms, results=results, forced=forced, settings=settings, n_stages=len(stages),
                           expected_table=point["table"], realized=realized)
    truth = built["dataset"].truth
    checks["construction"] = check_construction(truth, cell)
    span = truth.get("stage_span") or {}
    extra = {"truth": {
        "intended_geometry": {key: (span.get("intended_geometry") or {}).get(key) for key in STATISTICS},
        "path_length": span.get("path_length"), "mean_step": span.get("mean_step"),
        "plant_norm": span.get("plant_norm"), "scale_iterations": span.get("scale_iterations"),
    }}
    return _finish(base, checks, arms, CHECKS[SIMULATION], start, timing, **extra)


# ── Settings files, verification and pre-flights (design D10) ───────────────


def profile_dir(out_dir: Path, profile: str) -> Path:
    if profile not in PROFILES:
        raise SystemExit(f"--profile must be one of {list(PROFILES)}, got {profile!r}.")
    return Path(out_dir) / profile


def verify(settings: Mapping[str, Any], cohort: aac.Cohort | None = None) -> None:
    """Refuse a changed parent file, or an export whose files differ from the settings', by name."""

    plsc.verify_parents(settings)
    if cohort is not None:
        plsc.verify_export(cohort, settings)


def _prepare(directory: Path, settings: Mapping[str, Any], cohort: aac.Cohort | None) -> None:
    if (directory / "settings.json").exists():
        verify(lsc.load_settings(directory))
    verify(settings, cohort)
    lsc.ensure_settings(directory, settings)


def load_preflight(directory: Path, settings: Mapping[str, Any]) -> dict[str, Any] | None:
    path = Path(directory) / PREFLIGHT_NAME
    if not path.exists():
        return None
    record = json.loads(path.read_text(encoding="utf-8"))
    return record if record.get("settings_hash") == lsc.settings_hash(settings) else None


def plants_preflight(out_dir: Path, settings: Mapping[str, Any], cohort: aac.Cohort) -> dict[str, Any]:
    """The anchor driver's pre-flight units under this driver's settings; writes ``plants/preflight.json``."""

    directory = profile_dir(out_dir, PLANTS)
    _prepare(directory, settings, cohort)
    units = []
    for replicate in settings["preflight_replicates"]:
        try:
            units.append(aac.preflight_unit(int(replicate), settings, cohort))
        except Exception as exc:  # a pre-flight that raises fails
            units.append({"unit": aac.unit_key(COHORT_POINT, int(replicate)), "passed": False,
                          "error": f"{type(exc).__name__}: {exc}"})
    result = {
        "profile": PLANTS,
        "settings_hash": lsc.settings_hash(settings),
        "rule": "the cross-validated pls_cv on the cohort anchor's units selects the contract's rank and AUROC "
                "(rounded), and its statistics and p-values equal the call forced at that rank exactly",
        "expected": dict(settings["contract"]["fidelity"]),
        "units": units,
        "passed": bool(units) and all(unit["passed"] for unit in units),
    }
    lsc._write_json(directory / PREFLIGHT_NAME, result)
    return result


def committed_pvalues(path: Path, replicate: int) -> dict[str, str]:
    """The ADNI-design run's committed anchor p-values (as written) at one replicate."""

    found: dict[str, str] = {}
    for row in lsc._read_csv(Path(path)):
        if row["mode"] == adc.ANCHOR and float(row["effect"]) == 0.0 and int(row["replicate"]) == int(replicate):
            found[str(row["statistic"])] = str(row["p_value"])
    return found


def simulation_preflight_unit(replicate: int, parent_settings: Mapping[str, Any], pvalues: Path) -> dict[str, Any]:
    """One anchor unit of the ADNI-design run, by its own evaluator under its committed settings."""

    mark = perf_counter()
    record = adc.evaluate_unit(adc.ANCHOR, 0.0, int(replicate), parent_settings, reference=_reference())
    seconds = perf_counter() - mark
    observed = (record["arms"][adc.ADJUSTED] or {}).get("p_values") or {}
    expected = committed_pvalues(pvalues, replicate)
    differing = []
    for statistic in STATISTICS:
        value, committed = observed.get(statistic), expected.get(statistic)
        if value is None or committed in (None, "") or float(value) != float(committed):
            differing.append({"statistic": statistic, "committed": committed, "observed": value})
    return {
        "unit": unit_key("adni_design", adc.ANCHOR, replicate),
        "replicate": int(replicate),
        "committed": expected,
        "observed": {statistic: observed.get(statistic) for statistic in STATISTICS},
        "differing": differing,
        "seconds": seconds,
        "passed": not differing,
    }


def simulation_preflight(out_dir: Path, settings: Mapping[str, Any]) -> dict[str, Any]:
    """Anchor replicates 0 and 1 of the ADNI-design run against its committed p-values (``preflight.json``)."""

    directory = profile_dir(out_dir, SIMULATION)
    _prepare(directory, settings, None)
    stored = plsc.load_parent(plsc._resolve(settings["parents"]["design_settings"]["path"]).parent)
    pvalues = plsc._resolve(settings["parents"]["design_pvalues"]["path"])
    units = []
    for replicate in settings["preflight_replicates"]:
        try:
            units.append(simulation_preflight_unit(int(replicate), stored["settings"], pvalues))
        except Exception as exc:  # a pre-flight that raises fails
            units.append({"unit": unit_key("adni_design", adc.ANCHOR, int(replicate)), "passed": False,
                          "error": f"{type(exc).__name__}: {exc}"})
    result = {
        "profile": SIMULATION,
        "settings_hash": lsc.settings_hash(settings),
        "rule": "anchor replicates of the ADNI-design run, analysed by its own unit evaluator under its committed "
                "settings, give its committed p-values exactly",
        "parent": {key: settings["parents"][key] for key in ("design_settings", "design_pvalues")},
        "units": units,
        "passed": bool(units) and all(unit["passed"] for unit in units),
    }
    lsc._write_json(directory / PREFLIGHT_NAME, result)
    return result


# ── Shards, resume and the run ───────────────────────────────────────────────


def read_records(paths: Iterable[Path]) -> list[dict[str, Any]]:
    """Every parsable unit line; a line torn by an interrupted write is dropped, so its unit is re-run."""

    return psc.read_records(paths)


_WORKER_COHORT: aac.Cohort | None = None


def _init_worker(directory: str | None, contract: aac.CohortContract | None, files: Mapping[str, str] | None) -> None:
    global _WORKER_COHORT
    if directory is None or contract is None:
        return
    cohort = aac.load_cohort(Path(directory), contract)
    if cohort.files != dict(files or {}):
        raise SystemExit(f"{directory} changed after the run started (file hashes differ).")
    _WORKER_COHORT = cohort


def evaluate_unit(unit: Unit, settings: Mapping[str, Any], cohort: aac.Cohort | None = None) -> dict[str, Any]:
    if settings["profile"] == PLANTS:
        if cohort is None:
            raise ValueError("The plants profile needs the cohort.")
        return evaluate_plants_unit(unit, settings, cohort)
    return evaluate_simulation_unit(unit, settings)


def _evaluate_pending(unit: Unit, settings: Mapping[str, Any]) -> dict[str, Any]:
    return evaluate_unit(unit, settings, _WORKER_COHORT)


def run_shard(
    out_dir: Path,
    settings: Mapping[str, Any],
    index: int,
    count: int,
    *,
    cohort: aac.Cohort | None = None,
    limit: int | None = None,
    workers: int = 1,
) -> int:
    """Evaluate this shard's missing units of the settings' profile; returns how many were evaluated."""

    profile = str(settings["profile"])
    directory = profile_dir(out_dir, profile)
    if profile == PLANTS and cohort is None:
        raise SystemExit("The plants profile reads the export: --export-dir is required.")
    _prepare(directory, settings, cohort if profile == PLANTS else None)
    checked = load_preflight(directory, settings)
    if checked is None:
        raise SystemExit(f"{directory / PREFLIGHT_NAME} is missing or was written under other settings; run "
                         f"`preflight --profile {profile}` first.")
    if not checked.get("passed"):
        raise SystemExit(f"The pre-flight failed ({directory / PREFLIGHT_NAME}); the run stops.")
    path = lsc.shard_path(directory, index, count)
    path.parent.mkdir(parents=True, exist_ok=True)
    done = {record_key(record): record for record in (read_records([path]) if path.exists() else [])}
    # Rewritten with complete units only, so a torn line is re-run, not duplicated.
    with path.open("w", encoding="utf-8") as handle:
        handle.write("".join(lsc._line(done[key]) for key in sorted(done)))
    pending = [unit for position, unit in enumerate(enumerate_units(settings))
               if position % count == index and unit_key(*unit) not in done]
    if limit is not None:
        pending = pending[:limit]

    def write(position: int, unit: Unit, record: Mapping[str, Any]) -> None:
        with path.open("a", encoding="utf-8") as handle:
            handle.write(lsc._line(record))
            handle.flush()
        print(f"[{profile} {index}/{count}] {position}/{len(pending)} {unit_key(*unit)} {record['status']}",
              file=sys.stderr, flush=True)

    if workers <= 1:
        for position, unit in enumerate(pending, start=1):
            write(position, unit, evaluate_unit(unit, settings, cohort))
        return len(pending)
    initargs = ((str(cohort.directory), cohort.contract, dict(cohort.files)) if cohort is not None and
                profile == PLANTS else (None, None, None))
    with ProcessPoolExecutor(max_workers=workers, initializer=_init_worker, initargs=initargs) as pool:
        futures = {pool.submit(_evaluate_pending, unit, settings): unit for unit in pending}
        for position, future in enumerate(as_completed(futures), start=1):
            write(position, futures[future], future.result())
    return len(pending)


# ── Merge: rate, paired and cell tables, accounting, reproduction ───────────


RATE_COLUMNS: tuple[str, ...] = (
    "profile", "point", "participants", "cell", "kind", "reference", "mode", "effect", "role", "arm", "coordinates",
    "draws", "statistic", "rate", "se", "count", "rejections", "failures", "inconsistent",
)
PAIRED_COLUMNS: tuple[str, ...] = (
    "profile", "point", "cell", "kind", "reference", "mode", "effect", "role", "statistic", "arm", "reference_arm",
    "comparison", "n_paired", "rate_reference", "rate_arm", "difference", "b", "c", "p_loss", "p_gain",
    "p_two_sided",
)
CELL_COLUMNS: tuple[str, ...] = ("profile", "point", "cell", "quantity", "arm", "statistic", "key", "value")


def role(cell: Mapping[str, Any], statistic: str) -> str:
    if cell["kind"] == CONTROL:
        return CONTROL
    return "target" if statistic == cell["target"] else "off_target"


def _rejects(p_value: Any, settings: Mapping[str, Any]) -> bool | None:
    return plsc._rejects(p_value, float(settings["alpha"]), str(settings["rejection"]))


def _p(record: Mapping[str, Any], arm: str, statistic: str) -> Any:
    return aac._p(record, arm, statistic)


def _group(records: Iterable[Mapping[str, Any]]) -> dict[tuple[str, str], list[Mapping[str, Any]]]:
    grouped: dict[tuple[str, str], list[Mapping[str, Any]]] = {}
    for record in records:
        grouped.setdefault((str(record["point"]), str(record["cell"])), []).append(record)
    return grouped


def _cell_columns(settings: Mapping[str, Any], cell: Mapping[str, Any]) -> dict[str, Any]:
    return {"kind": cell["kind"], "reference": cell["reference"] or "", "mode": cell["mode"] or "",
            "effect": cell["effect"]}


def rate_table(records: Sequence[Mapping[str, Any]], settings: Mapping[str, Any]) -> list[dict[str, Any]]:
    """Rejections (p ≤ α) per point × cell × arm × statistic; failed and inconsistent units excluded."""

    grouped = _group(records)
    specs = settings["arm_specs"]
    rows = []
    for point in point_names(settings):
        participants = point_spec(settings, point)["participants"]
        for cell in point_cells(settings, point):
            group = grouped.get((point, cell["name"]), [])
            inconsistent = sum(1 for record in group if record.get("status") == "inconsistent")
            for arm in settings["arms"]:
                failures = sum(1 for record in group if record.get("status") == "failed"
                               or (record.get("status") == "completed" and aac._arm(record, arm) is None))
                spec = specs[arm]
                for statistic in STATISTICS:
                    decisions = [_rejects(_p(record, arm, statistic), settings) for record in group]
                    available = [bool(value) for value in decisions if value is not None]
                    count, hits = len(available), sum(available)
                    rate = hits / count if count else None
                    rows.append({
                        "profile": settings["profile"], "point": point, "participants": participants,
                        "cell": cell["name"], **_cell_columns(settings, cell), "role": role(cell, statistic),
                        "arm": arm, "coordinates": spec["coordinates"], "draws": spec["draws"],
                        "statistic": statistic, "rate": rate,
                        "se": None if rate is None else math.sqrt(rate * (1.0 - rate) / count),
                        "count": count, "rejections": hits, "failures": failures, "inconsistent": inconsistent,
                    })
    return rows


def paired_table(records: Sequence[Mapping[str, Any]], settings: Mapping[str, Any]) -> list[dict[str, Any]]:
    """Each pair of :data:`PAIRS` on the same units, by the exact tests on the discordant ones."""

    grouped = _group(records)
    rows = []
    for point in point_names(settings):
        for cell in point_cells(settings, point):
            group = sorted(grouped.get((point, cell["name"]), []), key=lambda record: int(record["replicate"]))
            for statistic in STATISTICS:
                for arm, reference, comparison in PAIRS:
                    pairs = [(_rejects(_p(record, reference, statistic), settings),
                              _rejects(_p(record, arm, statistic), settings)) for record in group]
                    both = [(ref, cand) for ref, cand in pairs if ref is not None and cand is not None]
                    n = len(both)
                    b = sum(1 for ref, cand in both if ref and not cand)
                    c = sum(1 for ref, cand in both if cand and not ref)
                    rows.append({
                        "profile": settings["profile"], "point": point, "cell": cell["name"],
                        **_cell_columns(settings, cell), "role": role(cell, statistic), "statistic": statistic,
                        "arm": arm, "reference_arm": reference, "comparison": comparison, "n_paired": n,
                        "rate_reference": sum(ref for ref, _ in both) / n if n else None,
                        "rate_arm": sum(cand for _, cand in both) / n if n else None,
                        "difference": (c - b) / n if n else None,
                        "b": b, "c": c, **plsc.paired_p_values(b, c),
                    })
    return rows


def _median(values: Iterable[Any]) -> float | None:
    return plsc._quantile(values, 0.5)


def cell_rows(records: Sequence[Mapping[str, Any]], settings: Mapping[str, Any]) -> list[dict[str, Any]]:
    """The cell table in long form: one quantity per row; blank keys are not applicable."""

    grouped = _group(records)
    rows: list[dict[str, Any]] = []
    profile = settings["profile"]
    stage_labels = (settings["contract"]["stage_order"] if profile == PLANTS
                    else [str(index) for index in range(int(settings["construction"]["generator"]["n_stages"]))])
    floor_rank = len(stage_labels) - 1

    def add(point: str, cell: str, quantity: str, value: Any, **keys: Any) -> None:
        rows.append({"profile": profile, "point": point, "cell": cell, "quantity": quantity, **keys, "value": value})

    for point in point_names(settings):
        for cell in point_cells(settings, point):
            name = cell["name"]
            everything = grouped.get((point, name), [])
            units = [record for record in everything if record.get("status") == "completed"]
            add(point, name, "units_completed", len(units))
            ranks = [int(record["production"]["rank"]) for record in units]
            for rank in sorted(set(ranks)):
                add(point, name, "cv_rank_units", ranks.count(rank), key=str(rank))
            add(point, name, "cv_rank_median", _median(ranks))
            add(point, name, "cv_auroc_median", _median(record["production"].get("cv_mean_auroc") for record in units))
            for arm in settings["arms"]:
                add(point, name, "angle_null_q95_median",
                    _median((aac._arm(record, arm) or {}).get("angle_null_q95") for record in units), arm=arm)
                for statistic in STATISTICS:
                    add(point, name, "statistic_median",
                        _median(((aac._arm(record, arm) or {}).get("pair_statistics") or {}).get(statistic)
                                for record in units), arm=arm, statistic=statistic)
            classes = (("k_minus_1", [r for r in units if int(r["production"]["rank"]) == floor_rank]),
                       ("above_k_minus_1", [r for r in units if int(r["production"]["rank"]) > floor_rank]))
            for arm in settings["arms"]:
                for label, members in classes:
                    for statistic in STATISTICS:
                        decisions = [_rejects(_p(record, arm, statistic), settings) for record in members]
                        available = [bool(value) for value in decisions if value is not None]
                        add(point, name, "rank_class_count", len(available), arm=arm, statistic=statistic, key=label)
                        add(point, name, "rank_class_rejections", sum(available), arm=arm, statistic=statistic,
                            key=label)
            if profile == PLANTS:
                _plants_cell_rows(add, point, name, units, stage_labels)
                add(point, name, "subsample_redraws",
                    sum(int(record.get("subsample_redraws") or 0) for record in everything))
                add(point, name, "label_redraws", sum(int(record.get("label_redraws") or 0) for record in everything))
            else:
                truths = [record.get("truth") or {} for record in units]
                for statistic in STATISTICS:
                    add(point, name, "intended_median",
                        _median((truth.get("intended_geometry") or {}).get(statistic) for truth in truths),
                        statistic=statistic)
                add(point, name, "path_length_median", _median(truth.get("path_length") for truth in truths))
                for index, label in enumerate(stage_labels):
                    add(point, name, "plant_norm_median",
                        _median((truth.get("plant_norm") or [None] * len(stage_labels))[index] for truth in truths),
                        key=label)
    return rows


def _plants_cell_rows(add: Callable[..., None], point: str, name: str, units: Sequence[Mapping[str, Any]],
                      stage_labels: Sequence[str]) -> None:
    plants = [record.get("plant") for record in units if record.get("plant")]
    if plants:
        for index, label in enumerate(stage_labels):
            add(point, name, "plant_norm_median", _median(plant["norm_per_stage"][index] for plant in plants),
                key=label)
        add(point, name, "plant_max_abs_entry_median", _median(plant["max_abs_entry"] for plant in plants))
        add(point, name, "scale_drift_median", _median(plant["scale_drift"] for plant in plants))
        add(point, name, "scale_drift_max", max(float(plant["scale_drift"]) for plant in plants))
    captures = [record.get("capture") or {} for record in units]
    keys = sorted({key for capture in captures for key in capture})
    for key in keys:
        for space in ("span", "predictive"):
            add(point, name, "capture_median", _median((capture.get(key) or {}).get(space) for capture in captures),
                key=f"{key}:{space}")
    readings = [record.get("reliability") for record in units if record.get("reliability")]
    if readings:
        add(point, name, "reliability_median", _median(entry.get("correlation") for entry in readings))
        add(point, name, "reliability_full_reference_median",
            _median(entry.get("full_reference") for entry in readings))
        add(point, name, "reliability_eigenvalue_ratio_median",
            _median(entry.get("eigenvalue_ratio") for entry in readings))
        for index in range(2):
            add(point, name, "reliability_eigenvalue_median",
                _median((entry.get("eigenvalues") or [None, None])[index] for entry in readings), key=str(index + 1))
        add(point, name, "reliability_errors", sum(1 for entry in readings if entry.get("error")))


def accounting(records: Sequence[Mapping[str, Any]], settings: Mapping[str, Any], *, duplicates: int) -> dict:
    """Units, failures and inconsistent units by point and cell, arm failures, redraws, core-hours, timing."""

    expected = expected_units(settings)
    by_cell: dict[tuple[str, str], dict[str, int]] = {}
    inconsistent_units: list[dict[str, Any]] = []
    failed_units: list[dict[str, Any]] = []
    arm_failures: dict[str, int] = {}
    arm_seconds: dict[str, list[float]] = {}
    unit_seconds: dict[tuple[str, str], list[float]] = {}
    for record in records:
        key = (str(record["point"]), str(record["cell"]))
        counts = by_cell.setdefault(key, {"complete": 0, "failed": 0, "inconsistent": 0})
        counts["complete"] += 1
        timing = record.get("runtime_seconds") or {}
        if timing.get("unit") is not None:
            unit_seconds.setdefault(key, []).append(float(timing["unit"]))
        status = record.get("status")
        if status == "failed":
            counts["failed"] += 1
            failed_units.append({"unit": record_key(record), "error_type": record.get("error_type"),
                                 "error_message": record.get("error_message")})
            continue
        if status == "inconsistent":
            counts["inconsistent"] += 1
            inconsistent_units.append({
                "unit": record_key(record),
                "checks": {name: {key: value for key, value in record["consistency"][name].items()
                                  if key in ("deviation", "failing_arms", "failing", "differing")}
                           for name in record["inconsistent_checks"]},
            })
        for arm, entry in (record.get("arms") or {}).items():
            if entry.get("status") != "completed":
                arm_failures[str(arm)] = arm_failures.get(str(arm), 0) + 1
            elif entry.get("seconds") is not None:
                arm_seconds.setdefault(str(arm), []).append(float(entry["seconds"]))
    total_expected = sum(expected.values())
    complete = sum(counts["complete"] for counts in by_cell.values())
    all_seconds = [value for values in unit_seconds.values() for value in values]
    cells = []
    for (point, cell), wanted in expected.items():
        counts = by_cell.get((point, cell), {"complete": 0, "failed": 0, "inconsistent": 0})
        cells.append({"point": point, "cell": cell, "expected": int(wanted), **counts,
                      "usable": counts["complete"] - counts["failed"] - counts["inconsistent"],
                      "unit_seconds_median": _median(unit_seconds.get((point, cell), []))})
    return {
        "profile": settings["profile"],
        "settings_hash": lsc.settings_hash(settings),
        "expected_units": total_expected,
        "complete_units": complete,
        "missing_units": total_expected - complete,
        "failed_units": sum(counts["failed"] for counts in by_cell.values()),
        "inconsistent_units": sum(counts["inconsistent"] for counts in by_cell.values()),
        "by_cell": cells,
        "failed_unit_errors": failed_units,
        "inconsistent_unit_checks": inconsistent_units,
        "arm_failures": arm_failures,
        "duplicate_records": int(duplicates),
        "subsample_redraws": sum(int(record.get("subsample_redraws") or 0) for record in records),
        "label_redraws": sum(int(record.get("label_redraws") or 0) for record in records),
        "recorded_core_hours": sum(all_seconds) / 3600.0,
        "arm_seconds_median": {arm: _median(arm_seconds.get(arm, [])) for arm in settings["arms"]},
    }


def anchor_reproduction(rates_rows: Sequence[Mapping[str, Any]], committed_rows: Sequence[Mapping[str, Any]]) -> dict:
    """The cohort anchor's counts for the four arms against the anchor run's committed rate table."""

    ours = {(str(row["arm"]), str(row["statistic"])): row for row in rates_rows
            if str(row["point"]) == COHORT_POINT and str(row["cell"]) == ANCHOR}
    theirs = {(str(row["arm"]), str(row["statistic"])): row for row in committed_rows
              if str(row["point"]) == COHORT_POINT}
    compared, differing = [], []
    for arm in ARMS:
        for statistic in STATISTICS:
            committed, observed = theirs.get((arm, statistic)), ours.get((arm, statistic))
            entry = {
                "arm": arm, "statistic": statistic,
                "committed": None if committed is None else {"count": int(float(committed["count"])),
                                                             "rejections": int(float(committed["rejections"]))},
                "observed": None if observed is None else {"count": int(float(observed["count"])),
                                                           "rejections": int(float(observed["rejections"]))},
            }
            entry["equal"] = entry["committed"] is not None and entry["committed"] == entry["observed"]
            compared.append(entry)
            if not entry["equal"]:
                differing.append(entry)
    return {"rule": "the cohort anchor's available and rejected counts for the four arms equal the anchor run's "
                    "committed rate table", "rows": compared, "differing": differing, "passed": not differing}


def merge(out_dir: Path, profile: str, export_dir: Path | None = None) -> dict[str, Any]:
    """Merge a profile's shards and write its tables, the accounting and (plants) the reproduction."""

    directory = profile_dir(out_dir, profile)
    settings = lsc.load_settings(directory)
    verify(settings)
    if profile == PLANTS:
        if export_dir is None:
            raise SystemExit("The plants profile reads the export: --export-dir is required.")
        aac.verify_export_files(export_dir, settings)
    raw = read_records(sorted((directory / "records").glob("shard_*.jsonl")))
    stale = sorted({record_key(r) for r in raw if r.get("settings_hash") != lsc.settings_hash(settings)})
    if stale:
        raise SystemExit(f"{len(stale)} unit(s) were written under other settings, e.g. {stale[:3]}; refusing.")
    seen: dict[str, int] = {}
    by_key: dict[str, dict[str, Any]] = {}
    for record in raw:
        key = record_key(record)
        seen[key] = seen.get(key, 0) + 1
        by_key[key] = record
    records = [by_key[unit_key(*unit)] for unit in enumerate_units(settings) if unit_key(*unit) in by_key]
    with (directory / "merged.jsonl").open("w", encoding="utf-8") as handle:
        handle.write("".join(lsc._line(record) for record in records))
    rates = rate_table(records, settings)
    lsc._write_csv(directory / "rate_table.csv", rates, RATE_COLUMNS)
    lsc._write_csv(directory / "paired_table.csv", paired_table(records, settings), PAIRED_COLUMNS)
    lsc._write_csv(directory / "cell_table.csv", cell_rows(records, settings), CELL_COLUMNS)
    summary = accounting(records, settings, duplicates=sum(n - 1 for n in seen.values()))
    lsc._write_json(directory / "accounting.json", summary)
    out = {key: summary[key] for key in ("expected_units", "complete_units", "missing_units", "failed_units",
                                         "inconsistent_units", "duplicate_records", "recorded_core_hours")}
    if profile == PLANTS:
        committed = lsc._read_csv(plsc._resolve(settings["parents"]["anchor_rates"]["path"]))
        reproduced = anchor_reproduction(lsc._read_csv(directory / "rate_table.csv"), committed)
        reproduced["settings_hash"] = lsc.settings_hash(settings)
        lsc._write_json(directory / REPRODUCTION_NAME, reproduced)
        out["reproduction_passed"] = reproduced["passed"]
    return out


# ── Decision (pure functions; designs D12–D16) ───────────────────────────────


VERDICTS: tuple[str, ...] = (
    "gate_failed", "production_excess", "size_excess", "power_loss", "specificity_excess", "no_gain", "reopened",
)
CONSEQUENCES: dict[str, str] = {
    "gate_failed": (
        "No reading is interpreted. The cause is fixed and the affected units are rerun under the same rule."
    ),
    "production_excess": (
        "Production's unrestricted null with df exceeds size on a control at the named readings. No reading of "
        "predictive_orthonormal is interpreted. A diagnosis of that control comes before any other reading on ADNI. "
        "If the control is a translation on ADNI's blocks, the roadmap states that a group main effect along the "
        "stage program reaches the interaction statistics at ADNI's design."
    ),
    "size_excess": (
        "predictive_orthonormal exceeds size on a control this change adds, where production holds. Item 21's hold "
        "on the global null does not extend to that control. Open item 1's second reopening condition is not met "
        "for ADNI, and it stays an opt-in diagnostic."
    ),
    "power_loss": (
        "predictive_orthonormal loses target power against production at the named plants. The condition is not "
        "met for ADNI, and it stays an opt-in diagnostic. The measurement-visibility page records which "
        "differences its two measured dimensions give up at ADNI's design."
    ),
    "specificity_excess": (
        "predictive_orthonormal adds off-target rejections over production at the named cells, where the kind of "
        "the planted difference is known exactly. The condition is not met for ADNI, and it stays an opt-in "
        "diagnostic."
    ),
    "no_gain": (
        "No simulated target cell at ADNI's table gains at the family threshold, so the measured dimension buys "
        "nothing at ADNI's design in simulation. The condition is not met, and it stays an opt-in diagnostic. A "
        "gain recorded on ADNI's noise without a simulated one is a lead on the simulation of ADNI's design, not a "
        "qualification."
    ),
    "reopened": (
        "At ADNI's design predictive_orthonormal holds size on the added controls, loses no target power, adds no "
        "off-target rejection where the truth is exact, and gains over production in simulation. With item 21's "
        "size reading, open item 1's second reopening condition is met for ADNI's export at 750 participants and "
        "the full cohort. It stays an opt-in diagnostic: a separate change may propose its paper-grade run with a "
        "rank rule for a fixed measured dimension, and only that run can change a default."
    ),
}
STANDING = (
    "In every verdict the SEA-AD anchor's reading stands for SEA-AD's design, keep_production stands, and the "
    "out-of-span cells are reported beside the verdict and were not judged."
)
LICENSES: dict[str, str] = {
    "L1": "The reading speaks for ADNI's export, the two design points and the simulated table. Nothing is extended "
          "to 80, 200 or 500 participants, where only size was read.",
    "L2": "At the full cohort, power and specificity were read on in-sample plants and in simulation. No held-out "
          "plant exists there. The held-out readings are those of 750 participants.",
    "L3": "reopened permits a proposal. It changes no default and adds no option. The later change carries the rank "
          "rule open item 1 asks for.",
    "L4": "Cohorts do not vouch for one another. The SEA-AD veto stands for SEA-AD's design, and ROSMAP has no "
          "reading.",
    "L5": "The out-of-span cells license nothing. Outside the stage span the split test's complement is the "
          "detector, and its reading on ADNI's noise is not this change's.",
}
CRITERIA_ORDER: tuple[tuple[str, str], ...] = (
    ("S1_production", "production_excess"), ("S1_candidate", "size_excess"), ("S2", "power_loss"),
    ("S3", "specificity_excess"), ("G", "no_gain"),
)

#: A reading: (profile, point, cell, statistic).
Reading = tuple[str, str, str, str]


def default_plan_settings() -> dict[str, dict[str, Any]]:
    """Each profile's points and cells at the module's values, for a decision with a profile missing."""

    plants_points = [{"name": aac.point_name(size), "participants": size, "subsampled": True}
                     for size in POINT_SIZES] + [{"name": COHORT_POINT, "participants": None, "subsampled": False}]
    point = simulation_point(SIMULATION_TABLE)
    return {
        PLANTS: {"profile": PLANTS, "points": plants_points, "cells": plants_cells(plants_points),
                 "replicates": {CONTROL: 0, PLANT: 0}, "permutations": 0},
        SIMULATION: {"profile": SIMULATION, "points": [{"name": point, "participants": None, "subsampled": False}],
                     "cells": simulation_cells(point), "replicates": {CONTROL: 0, PLANT: 0}, "permutations": 0},
    }


def _readings(settings: Mapping[str, Any], *, kinds: Sequence[str], references: Sequence[str] | None = None,
              role_: str, points: Sequence[str] | None = None, names: Sequence[str] | None = None) -> list[Reading]:
    """The (profile, point, cell, statistic) readings of the cells matching the filters, in enumeration order."""

    out = []
    for point in point_names(settings) if points is None else points:
        for cell in point_cells(settings, point):
            if cell["kind"] not in kinds or (names is not None and cell["name"] not in names):
                continue
            if references is not None and cell["reference"] not in references:
                continue
            for statistic in STATISTICS:
                if role(cell, statistic) == role_:
                    out.append((str(settings["profile"]), point, cell["name"], statistic))
    return out


def criteria_plan(plants: Mapping[str, Any], simulation: Mapping[str, Any]) -> dict[str, list[Reading]]:
    """Every family of design D12 and the recorded families beside it, as readings."""

    subsampled = subsampled_points(plants)
    translation = [cell["name"] for cell in plants["cells"] if cell["kind"] == CONTROL and cell["mode"]]
    size_plants = (_readings(plants, kinds=(CONTROL,), role_=CONTROL, points=subsampled, names=[ANCHOR])
                   + _readings(plants, kinds=(CONTROL,), role_=CONTROL, names=translation))
    return {
        "S1_plants": size_plants,
        "S1_simulation": _readings(simulation, kinds=(CONTROL,), role_=CONTROL),
        "S2": (_readings(simulation, kinds=(PLANT,), role_="target")
               + _readings(plants, kinds=(PLANT,), references=(HELDOUT,), role_="target")
               + _readings(plants, kinds=(PLANT,), references=(SAMPLE,), role_="target")),
        "S3_simulation": _readings(simulation, kinds=(PLANT,), role_="off_target"),
        "S3_in_sample": _readings(plants, kinds=(PLANT,), references=(SAMPLE,), role_="off_target"),
        "G": _readings(simulation, kinds=(PLANT,), role_="target"),
        "gain_heldout": _readings(plants, kinds=(PLANT,), references=(HELDOUT,), role_="target"),
        "gain_in_sample": _readings(plants, kinds=(PLANT,), references=(SAMPLE,), role_="target"),
        "heldout_off_target": _readings(plants, kinds=(PLANT,), references=(HELDOUT,), role_="off_target"),
        "out_of_span": (_readings(plants, kinds=(RECORDED,), role_="target")
                        + _readings(plants, kinds=(RECORDED,), role_="off_target")),
    }


def index_rates(rows: Iterable[Mapping[str, Any]]) -> dict[tuple[str, ...], Mapping[str, Any]]:
    return {(str(row["profile"]), str(row["point"]), str(row["cell"]), str(row["arm"]), str(row["statistic"])): row
            for row in rows}


def index_paired(rows: Iterable[Mapping[str, Any]]) -> dict[tuple[str, ...], Mapping[str, Any]]:
    return {(str(row["profile"]), str(row["point"]), str(row["cell"]), str(row["statistic"]), str(row["arm"]),
             str(row["reference_arm"])): row for row in rows}


def observation(rates: Mapping[tuple[str, ...], Mapping[str, Any]], reading: Reading, arm: str) -> dict[str, Any]:
    profile, point, cell, statistic = reading
    entry = {"profile": profile, "point": point, "cell": cell, "arm": arm, "statistic": statistic}
    row = rates.get((profile, point, cell, arm, statistic))
    if row is None or row.get("count") in (None, "") or int(float(row["count"])) <= 0:
        return {**entry, "status": "unavailable"}
    return {**entry, "n_available": int(float(row["count"])), "n_rejected": int(float(row["rejections"]))}


def paired_observation(paired: Mapping[tuple[str, ...], Mapping[str, Any]], reading: Reading, arm: str,
                       reference: str) -> dict[str, Any]:
    profile, point, cell, statistic = reading
    entry = {"profile": profile, "point": point, "cell": cell, "statistic": statistic, "arm": arm,
             "reference_arm": reference}
    row = paired.get((profile, point, cell, statistic, arm, reference))
    if row is None or row.get("n_paired") in (None, "") or int(float(row["n_paired"])) <= 0:
        return {**entry, "status": "unavailable"}
    return {
        **entry,
        "n_paired": int(float(row["n_paired"])),
        "rate_reference": plsc._number(row.get("rate_reference")),
        "rate_candidate": plsc._number(row.get("rate_arm")),
        "difference": plsc._number(row.get("difference")),
        "b": int(float(row["b"])),
        "c": int(float(row["c"])),
        "p_loss": float(row["p_loss"]),
        "p_gain": float(row["p_gain"]),
        "p_two_sided": float(row["p_two_sided"]),
    }


def judge_size(rates: Mapping[tuple[str, ...], Mapping[str, Any]], readings: Sequence[Reading], arm: str) -> dict:
    judged = grade.judge_rates([observation(rates, reading, arm) for reading in readings])
    judged["failing_readings"] = [row for row in judged["observations"] if row.get("status") != "holds"]
    return judged


def _paired_family(paired: Mapping[tuple[str, ...], Mapping[str, Any]], readings: Sequence[Reading], *,
                   arm: str = CANDIDATE_ARM, reference: str = REFERENCE_ARM) -> list[dict[str, Any]]:
    return [paired_observation(paired, reading, arm, reference) for reading in readings]


def judge_family(cells: Sequence[Mapping[str, Any]], *, p_field: str, rule: str) -> dict[str, Any]:
    judged = psc.judge_paired(cells, p_field=p_field, rule=rule)
    judged["failing_readings"] = [row for row in judged["observations"] if row.get("status") != "holds"]
    return judged


def judge_gain(cells: Sequence[Mapping[str, Any]]) -> dict[str, Any]:
    """``production_latent_space_check.judge_gain`` with an unavailable reading failing (the pilot's wrapper)."""

    return psc.judge_gain(cells)


def gate(
    preflights: Mapping[str, Mapping[str, Any] | None], summaries: Mapping[str, Mapping[str, Any] | None],
    reproduction: Mapping[str, Any] | None,
) -> dict[str, Any]:
    failures = []
    for profile in PROFILES:
        record = preflights.get(profile)
        if record is None:
            failures.append(f"the {profile} pre-flight record is missing or was written under other settings")
        elif not record.get("passed"):
            failures.append(f"the {profile} pre-flight did not pass")
        summary = summaries.get(profile)
        if summary is None:
            failures.append(f"the {profile} profile is not merged")
        else:
            inconsistent = int(summary.get("inconsistent_units", 0) or 0)
            if inconsistent:
                failures.append(f"{inconsistent} inconsistent unit(s) in the {profile} profile")
    if reproduction is None:
        failures.append("the cohort anchor's reproduction is missing")
    elif not reproduction.get("passed"):
        failures.append(f"the cohort anchor does not reproduce the anchor run's counts "
                        f"({len(reproduction.get('differing') or [])} differing)")
    return {"passed": not failures, "failures": failures}


def evaluate_criteria(plan: Mapping[str, Sequence[Reading]], rates: Mapping[tuple[str, ...], Mapping[str, Any]],
                      paired: Mapping[tuple[str, ...], Mapping[str, Any]]) -> dict[str, Any]:
    """S1 for both arms, S2, S3 and G, every one evaluated; the recorded gains beside them."""

    size = {
        arm: {"plants": judge_size(rates, plan["S1_plants"], arm),
              "simulation": judge_size(rates, plan["S1_simulation"], arm)}
        for arm in (REFERENCE_ARM, CANDIDATE_ARM)
    }
    for entry in size.values():
        entry["holds"] = entry["plants"]["holds"] and entry["simulation"]["holds"]
    s3 = {
        "simulation": judge_family(_paired_family(paired, plan["S3_simulation"]), p_field="p_gain",
                                   rule="count_and_ceiling_on_paired_excess"),
        "in_sample": judge_family(_paired_family(paired, plan["S3_in_sample"]), p_field="p_gain",
                                  rule="count_and_ceiling_on_paired_excess"),
    }
    return {
        "S1_production": size[REFERENCE_ARM],
        "S1_candidate": size[CANDIDATE_ARM],
        "S2": judge_family(_paired_family(paired, plan["S2"]), p_field="p_loss",
                           rule="count_and_ceiling_on_paired_loss"),
        "S3": {**s3, "holds": s3["simulation"]["holds"] and s3["in_sample"]["holds"]},
        "G": judge_gain(_paired_family(paired, plan["G"])),
        "recorded_gains": {
            "heldout": judge_gain(_paired_family(paired, plan["gain_heldout"])),
            "in_sample": judge_gain(_paired_family(paired, plan["gain_in_sample"])),
            "note": "judged by G's rule and recorded; they never qualify the candidate",
        },
    }


def named_readings(criterion: str, judged: Mapping[str, Any]) -> list[dict[str, Any]]:
    if criterion in ("S1_production", "S1_candidate"):
        return [*judged["plants"]["failing_readings"], *judged["simulation"]["failing_readings"]]
    if criterion == "S3":
        return [*judged["simulation"]["failing_readings"], *judged["in_sample"]["failing_readings"]]
    if criterion == "G":
        return [row for row in judged["observations"] if row.get("status") in ("no_gain", "unavailable")]
    return list(judged["failing_readings"])


def verdict_of(gated: Mapping[str, Any], criteria: Mapping[str, Any]) -> tuple[str, list[str]]:
    """The first verdict that applies, and every failing criterion in order."""

    failing = [name for name, _ in CRITERIA_ORDER if not criteria[name]["holds"]]
    if not gated["passed"]:
        return "gate_failed", failing
    for name, verdict in CRITERIA_ORDER:
        if name in failing:
            return verdict, failing
    return "reopened", failing


def indicative(settings_by_profile: Mapping[str, Mapping[str, Any] | None],
               summaries: Mapping[str, Mapping[str, Any] | None]) -> dict[str, Any]:
    reasons = []
    for profile in PROFILES:
        settings = settings_by_profile.get(profile)
        if settings is None:
            reasons.append(f"the {profile} profile has no settings")
            continue
        replicates = settings["replicates"]
        if int(replicates[CONTROL]) < CONTROL_REPLICATES:
            reasons.append(f"{profile}: {replicates[CONTROL]} replicates per control (paper grade "
                           f"{CONTROL_REPLICATES})")
        if int(replicates[PLANT]) < PLANT_REPLICATES:
            reasons.append(f"{profile}: {replicates[PLANT]} replicates per plant (paper grade {PLANT_REPLICATES})")
        if int(settings["permutations"]) < PERMUTATIONS:
            reasons.append(f"{profile}: {settings['permutations']} permutations (paper grade {PERMUTATIONS})")
        for entry in (summaries.get(profile) or {}).get("by_cell", []):
            cell = cell_spec(settings, entry["cell"])
            grade_ = CONTROL_REPLICATES if cell["kind"] == CONTROL else PLANT_REPLICATES
            if int(entry.get("usable", 0)) < grade_:
                reasons.append(f"{profile}: {entry.get('usable', 0)} usable units at {entry['point']}/{entry['cell']}")
    return {"indicative": bool(reasons), "reasons": reasons}


def _bound_phrase(reference: str) -> str:
    return ("a held-out plant bounds the detectable effect from above" if reference == HELDOUT
            else "an in-sample plant bounds it from below")


def power_statement(plants: Mapping[str, Any], rates: Mapping[tuple[str, ...], Mapping[str, Any]]) -> list[dict]:
    """Design D15: per point, reference and target statistic, ``informative_from`` on ``pls_cv`` and a sentence."""

    out = []
    for point in point_names(plants):
        participants = point_spec(plants, point).get("participants")
        cells = point_cells(plants, point)
        for reference in (SAMPLE, HELDOUT):
            for mode in (MAGNITUDE, ORIENTATION, SHAPE):
                ladder = sorted((cell for cell in cells if cell["kind"] == PLANT and cell["reference"] == reference
                                 and cell["mode"] == mode), key=lambda cell: float(cell["effect"]))
                if not ladder:
                    continue
                statistic = TARGETS[mode]
                effects = []
                for cell in ladder:
                    seen = observation(rates, (PLANTS, point, cell["name"], statistic), REFERENCE_ARM)
                    rate = (None if seen.get("status") == "unavailable"
                            else seen["n_rejected"] / seen["n_available"])
                    effects.append({"cell": cell["name"], "effect": float(cell["effect"]), "rate": rate,
                                    "count": seen.get("n_available"), "rejections": seen.get("n_rejected")})
                available = [entry for entry in effects if entry["rate"] is not None]
                informative = [entry for entry in available if entry["rate"] >= POWER_FLOOR]
                first = min(informative, key=lambda entry: entry["effect"]) if informative else None
                top = max(available, key=lambda entry: (entry["rate"], -entry["effect"])) if available else None
                kind = "a held-out" if reference == HELDOUT else "an in-sample"
                where = f"On ADNI's blocks at {point}" + (f" ({participants} participants)" if participants else "")
                if first is not None:
                    sentence = (f"{where}, production's `{statistic}` rejects {kind} {NOUNS[mode]} of e = "
                                f"{first['effect']:g} in {first['rejections']} of {first['count']} units "
                                f"({first['rate']:.3f}, floor {POWER_FLOOR:.2f}); {_bound_phrase(reference)}.")
                elif top is not None:
                    sentence = (f"{where}, production's `{statistic}` stays below the {POWER_FLOOR:.2f} floor on "
                                f"{kind} {NOUNS[mode]} up to e = {ladder[-1]['effect']:g}; its largest rate is "
                                f"{top['rate']:.3f} at e = {top['effect']:g} ({top['rejections']} of {top['count']}); "
                                f"{_bound_phrase(reference)}.")
                else:
                    sentence = f"{where}, production's `{statistic}` has no available {kind} {NOUNS[mode]} reading."
                out.append({
                    "point": point, "reference": reference, "mode": mode, "statistic": statistic,
                    "floor": POWER_FLOOR, "effects": effects,
                    "informative_from": first["effect"] if first else None,
                    "informative_rate": first["rate"] if first else None,
                    "max_rate": top["rate"] if top else None, "effect_at_max": top["effect"] if top else None,
                    "statement": sentence,
                })
    return out


def production_statement(plants: Mapping[str, Any], plan: Mapping[str, Sequence[Reading]],
                         rates: Mapping[tuple[str, ...], Mapping[str, Any]]) -> dict[str, Any]:
    """Design D15: production's power ladder and its own off-target family on the in-sample plants."""

    off_target = judge_size(rates, plan["S3_in_sample"], REFERENCE_ARM)
    return {
        "arm": REFERENCE_ARM,
        "power": power_statement(plants, rates),
        "off_target": {
            "verdict": "specific" if off_target["holds"] else "leaks",
            "named_readings": off_target["failing_readings"],
            "family": off_target,
            "rule": "production's off-target readings on the in-sample plants as one family under the "
                    "count-plus-ceiling rule",
        },
        "note": "Read on pls_cv alone; it never changes the verdict.",
    }


def _cell_value(cells_rows: Sequence[Mapping[str, Any]], profile: str, point: str, cell: str, quantity: str,
                **keys: str) -> float | None:
    for row in cells_rows:
        where = (str(row["profile"]), str(row["point"]), str(row["cell"]), str(row["quantity"]))
        if where == (profile, point, cell, quantity) and all(str(row.get(key) or "") == value
                                                             for key, value in keys.items()):
            return plsc._number(row["value"])
    return None


def descriptive_tables(
    plants: Mapping[str, Any], simulation: Mapping[str, Any], plan: Mapping[str, Sequence[Reading]],
    rates_rows: Sequence[Mapping[str, Any]], paired_rows: Sequence[Mapping[str, Any]],
    cells_rows: Sequence[Mapping[str, Any]],
) -> dict[str, Any]:
    """Design D16: none of these changes the verdict or production's statement."""

    rates = index_rates(rates_rows)
    paired = index_paired(paired_rows)
    bracket = []
    for point in subsampled_points(plants):
        cells = {(cell["reference"], cell["mode"], cell["effect"]): cell["name"] for cell in point_cells(plants, point)
                 if cell["kind"] == PLANT}
        for (reference, mode, effect), name in sorted(cells.items(), key=lambda item: str(item[0])):
            if reference != SAMPLE or (HELDOUT, mode, effect) not in cells:
                continue
            statistic = TARGETS[str(mode)]
            for arm in (REFERENCE_ARM, CANDIDATE_ARM):
                bracket.append({
                    "point": point, "mode": mode, "effect": effect, "statistic": statistic, "arm": arm,
                    "sample": observation(rates, (PLANTS, point, name, statistic), arm),
                    "heldout": observation(rates, (PLANTS, point, cells[(HELDOUT, mode, effect)], statistic), arm),
                })
    within = {
        arm: {"plants": judge_size(rates, plan["S1_plants"], arm),
              "simulation": judge_size(rates, plan["S1_simulation"], arm)}
        for arm in WITHIN_STAGE_ARMS
    }
    out_of_span = []
    for reading in plan["out_of_span"]:
        profile, point, cell, statistic = reading
        out_of_span.append({
            "point": point, "cell": cell, "statistic": statistic,
            "rates": {arm: observation(rates, reading, arm) for arm in ARMS},
            "paired": paired_observation(paired, reading, CANDIDATE_ARM, REFERENCE_ARM),
            "cv_rank_median": _cell_value(cells_rows, profile, point, cell, "cv_rank_median"),
            "capture": {key: _cell_value(cells_rows, profile, point, cell, "capture_median", key=key)
                        for key in ("plant:span", "plant:predictive", "reference:span", "reference:predictive")},
        })
    heldout = []
    for reading in plan["heldout_off_target"]:
        profile, point, cell, statistic = reading
        heldout.append({"point": point, "cell": cell, "statistic": statistic,
                        "rates": {arm: observation(rates, reading, arm) for arm in (REFERENCE_ARM, CANDIDATE_ARM)},
                        "reliability_median": _cell_value(cells_rows, profile, point, ANCHOR, "reliability_median")})
    geometry = []
    for cell in simulation.get("cells", []):
        for point in cell["points"]:
            geometry.append({
                "point": point, "cell": cell["name"],
                "intended": {statistic: _cell_value(cells_rows, SIMULATION, point, cell["name"], "intended_median",
                                                    statistic=statistic) for statistic in STATISTICS},
                "observed_median": {arm: {statistic: _cell_value(cells_rows, SIMULATION, point, cell["name"],
                                                                 "statistic_median", arm=arm, statistic=statistic)
                                          for statistic in STATISTICS} for arm in (REFERENCE_ARM, CANDIDATE_ARM)},
            })
    return {
        "note": "Every table here is a reading beside the verdict and production's statement and never changes them.",
        "bracket": bracket,
        "within_stage": {
            "size": within,
            "paired": [{key: row.get(key) for key in PAIRED_COLUMNS} for row in paired_rows
                       if str(row["comparison"]) == "draws" and str(row["kind"]) != CONTROL],
        },
        "out_of_span": out_of_span,
        "heldout_off_target": heldout,
        "simulation_geometry": geometry,
        "rank_classes": [{key: row.get(key) for key in CELL_COLUMNS} for row in cells_rows
                         if str(row["quantity"]) in ("rank_class_count", "rank_class_rejections")],
        "cells": [{key: row.get(key) for key in CELL_COLUMNS} for row in cells_rows
                  if str(row["quantity"]) not in ("rank_class_count", "rank_class_rejections")],
        "rates": [{key: row.get(key) for key in RATE_COLUMNS} for row in rates_rows],
        "paired": [{key: row.get(key) for key in PAIRED_COLUMNS} for row in paired_rows],
    }


def decision(
    settings_by_profile: Mapping[str, Mapping[str, Any] | None],
    tables: Mapping[str, Mapping[str, Sequence[Mapping[str, Any]]]],
    summaries: Mapping[str, Mapping[str, Any] | None],
    preflights: Mapping[str, Mapping[str, Any] | None],
    reproduction: Mapping[str, Any] | None,
) -> dict[str, Any]:
    """The predeclared verdict (rule version 1), what it licenses, production's statement, the descriptive readings."""

    defaults = default_plan_settings()
    plants = settings_by_profile.get(PLANTS) or defaults[PLANTS]
    simulation = settings_by_profile.get(SIMULATION) or defaults[SIMULATION]
    plan = criteria_plan(plants, simulation)
    rates_rows = [row for profile in PROFILES for row in (tables.get(profile) or {}).get("rates", [])]
    paired_rows = [row for profile in PROFILES for row in (tables.get(profile) or {}).get("paired", [])]
    cells_rows = [row for profile in PROFILES for row in (tables.get(profile) or {}).get("cells", [])]
    rates, paired = index_rates(rates_rows), index_paired(paired_rows)
    gated = gate(preflights, summaries, reproduction)
    criteria = evaluate_criteria(plan, rates, paired)
    verdict, failing = verdict_of(gated, criteria)
    named: list[dict[str, Any]] = []
    for name, outcome in CRITERIA_ORDER:
        if outcome == verdict:
            named = named_readings(name, criteria[name])
    return {
        "rule_version": RULE_VERSION,
        "driver_version": DRIVER_VERSION,
        "verdict": verdict,
        "consequence": f"{CONSEQUENCES[verdict]} {STANDING}",
        "standing": STANDING,
        "failing_criteria": failing,
        "named_readings": named,
        "consequences": dict(CONSEQUENCES),
        "order": [verdict_name for _, verdict_name in CRITERIA_ORDER],
        "gate": gated,
        "indicative": indicative(settings_by_profile, summaries),
        "criteria": criteria,
        "licenses": dict(LICENSES),
        "scoping_rule": dict(aac.SCOPING_RULE),
        "production_statement": None if verdict == "gate_failed" else production_statement(plants, plan, rates),
        "plan": {name: [list(reading) for reading in readings] for name, readings in plan.items()},
        "alpha": ALPHA,
        "rejection": REJECTION,
        "reproduction": reproduction,
        "descriptive": descriptive_tables(plants, simulation, plan, rates_rows, paired_rows, cells_rows),
    }


_TABLE_FILES: dict[str, str] = {"rates": "rate_table.csv", "paired": "paired_table.csv", "cells": "cell_table.csv"}


def _profile_inputs(directory: Path) -> tuple[Any, Any, Any, Any, dict[str, str]]:
    """A profile's settings, tables, accounting and pre-flight; ``None`` where a file is absent."""

    if not (directory / "settings.json").exists():
        return None, None, None, None, {}
    settings = lsc.load_settings(directory)
    verify(settings)
    checked = load_preflight(directory, settings)
    summary_path = directory / "accounting.json"
    if not summary_path.exists():
        return settings, None, None, checked, {}
    summary = json.loads(summary_path.read_text(encoding="utf-8"))
    if summary.get("settings_hash") != lsc.settings_hash(settings):
        raise SystemExit(f"{summary_path} was merged under other settings; run `merge` again.")
    tables = {}
    for name, filename in _TABLE_FILES.items():
        path = directory / filename
        if not path.exists():
            raise SystemExit(f"{path} does not exist; run `merge` first.")
        tables[name] = lsc._read_csv(path)
    hashes = {f"{directory.name}/{filename}": plsc.sha256_file(directory / filename)
              for filename in ("settings.json", *_TABLE_FILES.values(), "accounting.json")}
    if (directory / PREFLIGHT_NAME).exists():
        hashes[f"{directory.name}/{PREFLIGHT_NAME}"] = plsc.sha256_file(directory / PREFLIGHT_NAME)
    return settings, tables, summary, checked, hashes


def decide(out_dir: Path, export_dir: Path) -> dict[str, Any]:
    """Read both merged profiles (a missing one fails the gate) and write the decision."""

    out_dir = Path(out_dir)
    settings_by_profile: dict[str, Any] = {}
    tables: dict[str, Any] = {}
    summaries: dict[str, Any] = {}
    preflights: dict[str, Any] = {}
    inputs: dict[str, str] = {}
    for profile in PROFILES:
        directory = profile_dir(out_dir, profile)
        settings, table, summary, checked, hashes = _profile_inputs(directory)
        if profile == PLANTS and settings is not None:
            aac.verify_export_files(export_dir, settings)
        settings_by_profile[profile], tables[profile], summaries[profile], preflights[profile] = (
            settings, table, summary, checked)
        inputs.update(hashes)
    reproduction = None
    path = profile_dir(out_dir, PLANTS) / REPRODUCTION_NAME
    if path.exists() and settings_by_profile.get(PLANTS) is not None:
        stored = json.loads(path.read_text(encoding="utf-8"))
        if stored.get("settings_hash") == lsc.settings_hash(settings_by_profile[PLANTS]):
            reproduction = stored
            inputs[f"{PLANTS}/{REPRODUCTION_NAME}"] = plsc.sha256_file(path)
    result = decision(settings_by_profile, tables, summaries, preflights, reproduction)
    result["complete"], result["preflight"] = {}, {}
    for profile in PROFILES:
        summary, checked = summaries.get(profile), preflights.get(profile)
        result["complete"][profile] = None if summary is None else {
            key: summary[key] for key in ("expected_units", "complete_units", "missing_units", "failed_units",
                                          "inconsistent_units")}
        result["preflight"][profile] = None if checked is None else {"passed": checked.get("passed"),
                                                                     "units": checked.get("units")}
    result["inputs"] = inputs
    lsc._write_json(out_dir / DECISION_NAME, result)
    return result


# ── CLI ──────────────────────────────────────────────────────────────────────


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = parser.add_subparsers(dest="command", required=True)

    def common(command: argparse.ArgumentParser, *, profile: bool = True) -> None:
        if profile:
            command.add_argument("--profile", choices=PROFILES, required=True)
        command.add_argument("--export-dir", type=Path, default=None,
                             help="The ADNI export, read in place (adni_sexdiff/results/motco_input); required "
                                  "wherever the plants profile is read")
        command.add_argument("--out-dir", type=Path, required=True)

    def units(command: argparse.ArgumentParser) -> None:
        command.add_argument("--control-reps", type=int, default=None,
                             help=f"Replicates per control cell (default {CONTROL_REPLICATES})")
        command.add_argument("--plant-reps", type=int, default=None,
                             help=f"Replicates per plant cell (default {PLANT_REPLICATES})")
        command.add_argument("--perms", type=int, default=None, help=f"Permutations (default {PERMUTATIONS})")

    pre = sub.add_parser("preflight", help="Write a profile's settings and run its pre-flight")
    common(pre)
    units(pre)
    run = sub.add_parser("run", help="Evaluate (a shard of) a profile's units")
    common(run)
    units(run)
    run.add_argument("--shard", type=lsc.parse_shard, default=(0, 1), help="i/N, round-robin over the units")
    run.add_argument("--limit", type=int, default=None, help="Evaluate at most this many pending units")
    run.add_argument("--workers", type=int, default=1, help="Units evaluated in parallel (not part of the settings)")
    merge_parser = sub.add_parser("merge", help="Merge a profile's shards and write its tables")
    common(merge_parser)
    decide_parser = sub.add_parser("decide", help=f"Write {DECISION_NAME}")
    common(decide_parser, profile=False)
    return parser


def _load(directory: Path | None) -> aac.Cohort:
    if directory is None:
        raise SystemExit("The plants profile reads the export: --export-dir is required.")
    try:
        return aac.load_cohort(directory, aac.CONTRACTS["adni"])
    except aac.pc.ExportContractError as exc:
        raise aac.contract_exit(exc) from exc


def settings_from_args(args: argparse.Namespace) -> tuple[dict[str, Any], aac.Cohort | None]:
    counts = {"control_replicates": args.control_reps, "plant_replicates": args.plant_reps,
              "permutations": args.perms}
    if args.profile == PLANTS:
        cohort = _load(args.export_dir)
        return plants_settings(cohort, **counts), cohort
    return simulation_settings(**counts), None


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    if args.command in ("preflight", "run"):
        settings, cohort = settings_from_args(args)
        if args.command == "preflight":
            if args.profile == PLANTS:
                assert cohort is not None
                result = plants_preflight(args.out_dir, settings, cohort)
            else:
                result = simulation_preflight(args.out_dir, settings)
            print(json.dumps({"profile": args.profile, "passed": result["passed"], "units": [
                {key: unit.get(key) for key in ("unit", "passed", "fidelity", "forced_rank", "statistics_differing",
                                                "p_values_differing", "differing", "error") if key in unit}
                for unit in result["units"]
            ]}, indent=1))
            return 0 if result["passed"] else 1
        index, count = args.shard
        evaluated = run_shard(args.out_dir, settings, index, count, cohort=cohort, limit=args.limit,
                              workers=args.workers)
        print(f"evaluated {evaluated} unit(s)", file=sys.stderr)
        return 0
    if args.command == "merge":
        print(json.dumps(merge(args.out_dir, args.profile, args.export_dir), indent=1))
        return 0
    if args.export_dir is None:
        raise SystemExit("decide reads the plants profile: --export-dir is required.")
    result = decide(args.out_dir, args.export_dir)
    print(json.dumps({
        "verdict": result["verdict"],
        "failing_criteria": result["failing_criteria"],
        "gate": result["gate"],
        "indicative": result["indicative"],
        "production_off_target": (result["production_statement"] or {}).get("off_target", {}).get("verdict"),
        "complete": result["complete"],
    }, indent=1))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
