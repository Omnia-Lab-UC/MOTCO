#!/usr/bin/env python3
"""Split interaction test validation: size, out-of-span power and localization on two parent runs (no R).

The driver behind ``openspec/changes/archive/2026-11-05-span-split-interaction-test/`` (spec
``span-split-validation``, roadmap item 16). The production measurement sees a
group difference only through its component inside the fitted stage-PLS span;
:func:`motco.stats.interaction.interaction_split_test` tests the group × stage
interaction energy of the standardized features inside that span (``span``), in
its orthogonal complement (``complement``) and in all of them (``total``), each
by a trace pseudo-F under Freedman–Lane permutation. This driver measures the
test before a case study may use it, on units that reproduce two committed runs.

**Profiles.**

- ``simulation`` — the latent-space comparison's units
  (``results/latent-space-comparison-2026-11-04/``): its three designs
  (``n200``, ``seaad``, ``seaad_weak``), fourteen cells, base seed 2035 and
  seed derivations, imported from ``scripts/latent_space_comparison.py``
  unchanged. The weak δ is read from the comparison's committed
  ``calibration/calibration.json``; its sha256 must equal the comparison's
  settings, and the designs, cells and seeds must equal them too. The anchor
  and translation run replicates 0–499, every other cell 0–99 (the
  comparison's units), all at 999 permutations. Each unit generates the
  comparison's dataset, fits ``pls_cv`` once, evaluates it (net displacement,
  rank floor 2, span capture, null distributions kept) and runs the split
  test on the fit's standardized features in the PLS input coordinates
  against ``qr(x_rotations_)``, with the trajectory models, at the unit's RRPP
  seed: once uncorrected (the pseudo-F and the interaction-SS p-values) and
  once under ``df`` (the corrected SS p-value). Rejection is p < α, the
  comparison's rule.
- ``plasmode`` — the SEA-AD plasmode's units (``results/plasmode-2026-11-02/``):
  its export contract, within-stage label permutation, 20-feature plant and
  seeds (base seed 2034, ``paper`` stream), read through
  ``scripts/plasmode_check.py`` unchanged. The export's file hashes must equal
  the plasmode's committed settings (a mismatch exits naming the files), and
  ``--fidelity`` must be a passed gate on this export. Cells: the anchor,
  e = 2.0 and e = 3.0, 500 replicates × 999 permutations. Each unit is one
  ``analyze_trajectories`` at the plasmode's case-study parameters (six
  covariates, rank floor 2, net displacement, ``df``, harness CV) with
  attribution off and ``span_split_test`` on, then the plasmode's latent
  checks (every anchor unit must reproduce the gate's rank and AUROC). The
  split test is rerun uncorrected on the same inputs for the SS diagnostic.
  Rejection is p ≤ α, the plasmode's rule. A record holds seeds and
  statistics only, never a donor identifier, label vector or donor-level value.

A unit that raises is recorded as failed with its error; the shard continues.

**Seeds** are the parents': every simulation unit at replicates 0–99 is the
comparison's dataset and RRPP seed, every plasmode unit the plasmode's labels,
plant and RRPP seed. The split test draws RRPP's serial sequence at that seed,
so its draws are the production statistics' draws.

Subcommands (every output lands under ``--out-dir``/``<profile>/``):

- ``run --profile simulation|plasmode`` — evaluate (a shard of) the units into
  ``records/shard_<i>_of_<N>.jsonl``, one line per unit. The settings and
  their hash are written to ``settings.json``; a resume or merge under other
  settings exits naming the differing keys. Resume skips complete units.
- ``merge --profile …`` — ``merged.jsonl`` (per unit; keep off the
  repository), ``split_table.csv`` (rejection rate, Monte Carlo SE and count
  per design/cell × part × statistic: ``f`` the pseudo-F, ``ss`` the
  uncorrected interaction SS, ``ss_df`` the interaction SS under ``df``; the
  median observed in-span share), the production rates
  (``production_table.csv`` with ``capture_table.csv`` for the simulation,
  ``trajectory_table.csv`` for the plasmode), ``accounting.json`` and
  ``reproduction.json``. The simulation reproduction recomputes ``pls_cv``'s
  rejections from the first 199 draws of each shared unit (replicates 0–99)
  and requires the committed ``rate_table.csv`` count and rejections for every
  design, cell and statistic; the plasmode reproduction requires the
  committed ``plasmode_summary.csv`` rates of ``delta``, ``angle`` and
  ``shape`` exactly per cell.
- ``decide`` — ``span_split_decision.json`` (rule version 1): the gate
  ``reproduction_failed`` (no criterion read), then C1 simulation size of the
  three parts' pseudo-F at the anchor and translation of every design (m = 18)
  and C2 real-noise size at the plasmode anchor (m = 3), both under the stage-span
  paper grade's count-plus-ceiling rule (imported); C3 the complement's
  pseudo-F rate for orientation e = 0.5 at φ = 90 (dense) above ``pls_cv``'s
  ``angle`` on the same units by more than 2 pooled SE at ``n200`` and
  ``seaad``. Verdict: ``size_excess`` | ``no_out_of_span_gain`` |
  ``split_validated``, with its consequence and descriptive tables
  (localization, share against capture, span part against production,
  complement against the comparison's arms, plasmode power, the SS
  diagnostic) that never change it.

With ``--workers`` above 1, pin BLAS to one thread per process
(``OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1``).

Examples:

    export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
    EXPORT=../sead_sexdiff/results/motco_input/cohort_30

    # smoke
    python scripts/span_split_check.py run --profile simulation --reps 2 --control-reps 2 --perms 19 \\
        --out-dir /tmp/ssc --workers 7
    python scripts/plasmode_check.py fidelity --export-dir $EXPORT --out-dir /tmp/ssc/fidelity
    python scripts/span_split_check.py run --profile plasmode --export-dir $EXPORT \\
        --fidelity /tmp/ssc/fidelity/fidelity.json --reps 2 --perms 19 --out-dir /tmp/ssc --workers 6
    python scripts/span_split_check.py merge --profile simulation --out-dir /tmp/ssc
    python scripts/span_split_check.py merge --profile plasmode --out-dir /tmp/ssc
    python scripts/span_split_check.py decide --out-dir /tmp/ssc

    # one SLURM array task (scripts/span_split_check_array.sbatch), then merge and decide
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
from typing import Any, Callable, Iterable, Mapping, Sequence

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import latent_space_comparison as lsc  # noqa: E402  (designs, cells, seeds, settings files)
import plasmode_check as pc  # noqa: E402  (export contract, units, case-study parameters)
import stage_span_paper_grade_decision as grade  # noqa: E402  (the count-plus-ceiling rule)

sys.path.pop(0)

DRIVER_VERSION = 1
RULE_VERSION = 1
PROFILES: tuple[str, ...] = ("simulation", "plasmode")
PARTS: tuple[str, ...] = ("total", "span", "complement")
#: ``f`` the pseudo-F (the test), ``ss`` the uncorrected interaction SS, ``ss_df`` the interaction SS under ``df``.
SPLIT_STATISTICS: tuple[str, ...] = ("f", "ss", "ss_df")
STATISTICS: tuple[str, ...] = lsc.STATISTICS
ALPHA = 0.05
SE_MULTIPLIER = 2.0

DEFAULT_PERMUTATIONS = 999
REPRODUCTION_PERMUTATIONS = 199
CONTROL_CELLS: tuple[str, ...] = ("anchor", "translation")
CONTROL_REPLICATES = 500
SIMULATION_REPLICATES = 100
PLASMODE_CELLS: tuple[float, ...] = (0.0, 2.0, 3.0)
PLASMODE_REPLICATES = 500
GAIN_DESIGNS: tuple[str, ...] = lsc.STRONG_DESIGNS
GAIN_CELL = "orient_e0.5_phi90"
POWER_EFFECTS: tuple[float, ...] = (2.0, 3.0)

#: The rejection rule of each profile: its parent's, so its rates pair with the parent's.
REJECTION_RULES: dict[str, str] = {"simulation": "p < alpha", "plasmode": "p <= alpha"}

REPO = Path(__file__).resolve().parents[1]
COMPARISON_DIR = REPO / "results" / "latent-space-comparison-2026-11-04"
PLASMODE_DIR = REPO / "results" / "plasmode-2026-11-02"
COMPARISON_SOURCE = {"source": "latent-space comparison 2026-11-04", "grade": "100 × 199"}

VERDICTS: tuple[str, ...] = ("reproduction_failed", "size_excess", "no_out_of_span_gain", "split_validated")
CONSEQUENCES: dict[str, str] = {
    "reproduction_failed": (
        "Nothing is read: the units do not reproduce their parent runs, so no rate pairs with the comparison or "
        "the plasmode. Find the drift before any rerun."
    ),
    "size_excess": (
        "A case study may not report the split interaction test: its size is not established at the failing "
        "parts and designs."
    ),
    "no_out_of_span_gain": (
        "A case study may not report the split interaction test as seeing what the production statistics miss: "
        "its complement part does not reject an out-of-span rotation more often than angle does."
    ),
    "split_validated": (
        "A case study may report the split interaction test, reading its parts as the measurement-visibility "
        "documentation states: span = the fitted span, a complement rejection = energy outside it, a non-rejection "
        "uninformative without power at the design."
    ),
}

DECISION_NAME = "span_split_decision.json"


# ── Shared helpers ───────────────────────────────────────────────────────────


def sha256_file(path: Path) -> str:
    return pc.sha256_file(Path(path))


def _rejects(p_value: Any, alpha: float, rule: str) -> bool | None:
    """``None`` for an unavailable p-value."""

    if p_value is None or p_value == "":
        return None
    value = float(p_value)
    if not math.isfinite(value):
        return None
    return value < alpha if rule == "p < alpha" else value <= alpha


def _rate_row(rejections: Sequence[bool | None]) -> dict[str, Any]:
    available = [bool(value) for value in rejections if value is not None]
    count = len(available)
    hits = sum(available)
    rate = hits / count if count else None
    return {
        "rate": rate,
        "se": None if rate is None else math.sqrt(rate * (1.0 - rate) / count),
        "count": count,
        "rejections": hits,
    }


def _median(values: Iterable[Any]) -> float | None:
    finite = [float(value) for value in values if value is not None and math.isfinite(float(value))]
    return float(np.median(finite)) if finite else None


def _display_path(path: Path) -> str:
    resolved = Path(path).resolve()
    try:
        return str(resolved.relative_to(REPO))
    except ValueError:
        return str(resolved)


def _resolve(stored: str) -> Path:
    path = Path(stored)
    return path if path.is_absolute() else REPO / path


def load_parent(directory: Path) -> dict[str, Any]:
    """A committed parent run's ``settings.json``, refused when its hash does not match its settings."""

    path = Path(directory) / "settings.json"
    if not path.exists():
        raise SystemExit(f"{path} does not exist; the parent run is required.")
    stored = json.loads(path.read_text(encoding="utf-8"))
    if stored.get("hash") != lsc.settings_hash(stored.get("settings", {})):
        raise SystemExit(f"{path} fails its hash: the parent's settings were edited.")
    return stored


def _canonical(value: Any) -> Any:
    return json.loads(json.dumps(value, sort_keys=True))


def profile_dir(out_dir: Path, profile: str) -> Path:
    return Path(out_dir) / profile


def _write_csv(path: Path, rows: Sequence[Mapping[str, Any]], columns: Sequence[str]) -> None:
    lsc._write_csv(path, rows, columns)


def _write_json(path: Path, payload: Any) -> None:
    lsc._write_json(path, payload)


def _read_csv(path: Path) -> list[dict[str, Any]]:
    if not path.exists():
        raise SystemExit(f"{path} does not exist; run `merge` first.")
    return lsc._read_csv(path)


# ── Simulation profile: settings and units ───────────────────────────────────


def simulation_settings(
    *,
    parent_dir: Path = COMPARISON_DIR,
    calibration: Path | None = None,
    replicates: int = SIMULATION_REPLICATES,
    control_replicates: int = CONTROL_REPLICATES,
    permutations: int = DEFAULT_PERMUTATIONS,
    designs: Sequence[str] = lsc.DESIGNS,
) -> dict[str, Any]:
    """The comparison's designs, cells, seeds and weak δ, checked against its committed settings."""

    parent = load_parent(parent_dir)
    parent_settings = parent["settings"]
    if calibration is None:
        calibration = Path(parent_dir) / "calibration" / "calibration.json"
    calibration = Path(calibration)
    if not calibration.exists():
        raise SystemExit(f"{calibration} does not exist; the comparison's calibration is required.")
    digest = sha256_file(calibration)
    expected = (parent_settings.get("weak_delta") or {}).get("calibration_sha256")
    if digest != expected:
        raise SystemExit(
            f"{calibration} has sha256 {digest}, not the comparison's calibration_sha256 {expected}; refusing a "
            "weak δ the comparison did not use."
        )
    comparison = lsc.build_settings(
        calibration=calibration, replicates=int(parent_settings["replicates"]), permutations=permutations,
        designs=designs,
    )
    comparison["arms"] = [lsc.REFERENCE_ARM]
    shared = {
        "designs": {name: parent_settings["designs"][name] for name in comparison["designs"]},
        **{
            key: parent_settings[key]
            for key in ("cells", "base_seed", "pls_cv_params", "min_components", "orientation_estimator", "weak_delta")
        },
    }
    differing = [key for key, value in shared.items() if _canonical(comparison[key]) != _canonical(value)]
    if differing:
        raise SystemExit(
            f"The comparison's units differ from {parent_dir}/settings.json in: {', '.join(differing)}; "
            "the shared units would not be its datasets."
        )
    return {
        "driver_version": DRIVER_VERSION,
        "profile": "simulation",
        "parent": {
            "dir": _display_path(parent_dir),
            "settings_hash": parent["hash"],
            "rate_table_sha256": sha256_file(Path(parent_dir) / "rate_table.csv"),
            "shared_replicates": int(parent_settings["replicates"]),
            "permutations": int(parent_settings["permutations"]),
        },
        "comparison": comparison,
        "control_cells": list(CONTROL_CELLS),
        "replicates": {"control": int(control_replicates), "other": int(replicates)},
        "permutations": int(permutations),
        "reproduction_permutations": min(REPRODUCTION_PERMUTATIONS, int(permutations)),
        "alpha": ALPHA,
        "rejection": REJECTION_RULES["simulation"],
        "split": {
            "parts": list(PARTS),
            "coordinates": "standardized features in the PLS input coordinates (x - _x_mean) / _x_std",
            "basis": "qr(x_rotations_) of the unit's pls_cv fit",
            "models": "full group x stage and reduced group + stage trajectory models",
            "statistics": {"f": "pseudo-F, uncorrected", "ss": "interaction SS, uncorrected", "ss_df": "SS, df"},
        },
    }


def simulation_units(settings: Mapping[str, Any]) -> list[lsc.Unit]:
    control = set(settings["control_cells"])
    counts = settings["replicates"]
    return [
        lsc.Unit(design, cell[0], replicate)
        for design in settings["comparison"]["designs"]
        for cell in settings["comparison"]["cells"]
        for replicate in range(int(counts["control" if cell[0] in control else "other"]))
    ]


def _simulation_key(record: Mapping[str, Any]) -> str:
    return json.dumps([str(record["design"]), str(record["cell"]), int(record["replicate"])])


def _unit_key(unit: lsc.Unit) -> str:
    return json.dumps([unit.design, unit.cell, unit.replicate])


# ── Simulation profile: one unit ─────────────────────────────────────────────


def production_params(settings: Mapping[str, Any], *, seed: int, n_stages: int) -> Any:
    """``pls_cv``'s evaluation parameters at the profile's permutations, null distributions kept."""

    params = lsc.arm_params(lsc.REFERENCE_ARM, settings["comparison"], seed=seed, n_stages=n_stages)
    return dataclasses.replace(params, include_null_distributions=True)


def split_inputs(dataset: Any, latent: Any, params: Any) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """The fit's standardized features in its input coordinates, its basis and the trajectory models."""

    from motco.simulations.evaluation import build_simulation_trajectory_design

    linear_map = latent.linear_map
    features = linear_map.inputs(latent.artifacts.features.to_numpy(dtype=float))
    design = build_simulation_trajectory_design(
        dataset.metadata, group_col=params.group_col, stage_col=params.stage_col
    )
    return (
        features,
        np.asarray(linear_map.basis, dtype=float),
        np.asarray(design.model_full, dtype=float),
        np.asarray(design.model_reduced, dtype=float),
    )


def split_record(
    features: np.ndarray, basis: np.ndarray, full: np.ndarray, reduced: np.ndarray, *, permutations: int, seed: int
) -> dict[str, Any]:
    """The split test uncorrected, plus the interaction-SS p-value under ``df`` from a second run."""

    from motco.stats.interaction import interaction_split_test

    plain = interaction_split_test(features, full, reduced, basis, permutations=permutations, seed=seed)
    corrected = interaction_split_test(
        features, full, reduced, basis, permutations=permutations, seed=seed, residual_correction="df"
    )
    parts = {}
    for name in PARTS:
        entry = plain.parts[name].to_dict()
        entry["p_value_ss_df"] = corrected.parts[name].p_value_ss
        parts[name] = entry
    return {"parts": parts, "span_share": plain.span_share}


def production_record(result: Any, reproduction_permutations: int) -> dict[str, Any]:
    """``pls_cv``'s rank, AUROC, statistics, p-values (all draws and the first ones) and capture."""

    from motco.simulations.diagnostics import POPULATION_LATENT_CHECKPOINT
    from motco.simulations.evaluation import _empirical_p_value

    nulls = result.null_distributions or {}
    early = {
        statistic: _empirical_p_value(nulls[statistic][:reproduction_permutations], observed)
        for statistic, observed in result.pair_statistics.items()
        if statistic in nulls and np.isfinite(observed) and reproduction_permutations > 0
    }
    geometry = result.realized_geometry
    metadata = result.latent_matrix_metadata
    return {
        "rank": int(metadata["selected_lv"]),
        "cv_mean_auroc": metadata.get("cv_mean_auroc"),
        "pair_statistics": lsc._finite(result.pair_statistics),
        "p_values": lsc._finite(result.p_values),
        "p_values_reproduction": lsc._finite(early),
        "span_capture": geometry.get("span_capture"),
        "population_latent": (geometry.get("checkpoints", {}).get(POPULATION_LATENT_CHECKPOINT) or {}).get("joint"),
    }


def evaluate_simulation_unit(unit: lsc.Unit, settings: Mapping[str, Any]) -> dict[str, Any]:
    """One dataset, one ``pls_cv`` fit and evaluation, and the split test at the unit's RRPP seed."""

    from motco.simulations.evaluation import evaluate_semisynthetic_trajectory, integrate_semisynthetic_dataset

    start = perf_counter()
    comparison = settings["comparison"]
    cell = lsc.cell_spec(comparison, unit.cell)
    seed = lsc.rrpp_seed(comparison["base_seed"], unit.design, unit.cell, unit.replicate)
    base = {
        "profile": "simulation",
        "design": unit.design,
        "cell": unit.cell,
        "trajectory_mode": cell["mode"],
        "effect_size": float(cell["effect"]),
        "span_alignment": cell["phi"],
        "span_support": cell["support"],
        "replicate": unit.replicate,
        "generator_seed": lsc.generator_seed(comparison["base_seed"], unit.design, unit.replicate),
        "rrpp_seed": seed,
        "settings_hash": lsc.settings_hash(settings),
    }
    timing: dict[str, float] = {}
    try:
        dataset = lsc.unit_dataset(unit, comparison)
        n_stages = int(comparison["designs"][unit.design]["generator"]["n_stages"])
        params = production_params(settings, seed=seed, n_stages=n_stages)
        mark = perf_counter()
        latent = integrate_semisynthetic_dataset(dataset, params)
        timing["fit"] = perf_counter() - mark
        mark = perf_counter()
        result = evaluate_semisynthetic_trajectory(dataset, params, latent=latent)
        timing["evaluation"] = perf_counter() - mark
        production = production_record(result, int(settings["reproduction_permutations"]))
        mark = perf_counter()
        split = split_record(
            *split_inputs(dataset, latent, params), permutations=int(settings["permutations"]), seed=seed
        )
        timing["split"] = perf_counter() - mark
    except Exception as exc:  # recorded, never aborts the shard
        timing["unit"] = perf_counter() - start
        return {**base, "status": "failed", "error_type": type(exc).__name__, "error_message": str(exc),
                "runtime_seconds": timing}
    timing["unit"] = perf_counter() - start
    return {**base, "status": "completed", "production": production, "split": split, "runtime_seconds": timing}


# ── Plasmode profile: settings, units and one unit ───────────────────────────


def plasmode_settings(
    export: Any,
    fidelity: Mapping[str, Any],
    *,
    parent_dir: Path = PLASMODE_DIR,
    replicates: int = PLASMODE_REPLICATES,
    permutations: int = DEFAULT_PERMUTATIONS,
    cells: Sequence[float] = PLASMODE_CELLS,
) -> dict[str, Any]:
    """The plasmode's committed settings, refused on another export, gate or cell."""

    parent = load_parent(parent_dir)
    parent_settings = parent["settings"]
    committed = parent_settings["export"]["files"]
    differing = sorted(
        name for name in set(export.files) | set(committed) if export.files.get(name) != committed.get(name)
    )
    if differing:
        raise SystemExit(
            f"The export's file hashes differ from {parent_dir}/settings.json for: {', '.join(differing)}; "
            "refusing another export."
        )
    if _canonical(dict(fidelity)) != _canonical(parent_settings["fidelity"]):
        raise SystemExit(
            f"The fidelity gate {dict(fidelity)} differs from the plasmode's {parent_settings['fidelity']}."
        )
    unknown = [float(effect) for effect in cells if float(effect) not in [float(e) for e in parent_settings["cells"]]]
    if unknown:
        raise SystemExit(f"Cell(s) {unknown} are not the plasmode's ({parent_settings['cells']}).")
    return {
        "driver_version": DRIVER_VERSION,
        "profile": "plasmode",
        "parent": {
            "dir": _display_path(parent_dir),
            "settings_hash": parent["hash"],
            "summary_sha256": sha256_file(Path(parent_dir) / "plasmode_summary.csv"),
        },
        "plasmode": parent_settings,
        "cells": [float(effect) for effect in cells],
        "replicates": int(replicates),
        "permutations": int(permutations),
        "fidelity": dict(fidelity),
        "attribution": False,
        "span_split_test": True,
        "alpha": ALPHA,
        "rejection": REJECTION_RULES["plasmode"],
    }


def plasmode_units(settings: Mapping[str, Any]) -> list[tuple[float, int]]:
    """Replicate-major, so a shard holds the matched cells of its replicates."""

    return [
        (float(effect), replicate) for replicate in range(int(settings["replicates"])) for effect in settings["cells"]
    ]


def _plasmode_key(record: Mapping[str, Any]) -> str:
    return pc.unit_key(float(record["effect"]), int(record["replicate"]))


def plasmode_params(settings: Mapping[str, Any], seeds: Mapping[str, int]) -> Any:
    """The plasmode's case-study parameters with attribution off and the split test on."""

    parent = {**settings["plasmode"], "attribution": None}
    params = pc.analysis_params(parent, seeds, permutations=int(settings["permutations"]))
    return dataclasses.replace(params, span_split_test=True)


def uncorrected_split(dataset: Any, result: Any, params: Any) -> Any:
    """The split test on the analysis's own inputs without the residual correction (the SS diagnostic)."""

    from motco import analysis
    from motco.simulations.evaluation import build_simulation_trajectory_design
    from motco.stats.interaction import interaction_split_test

    features, basis = analysis._span_split_inputs(result.latent_artifacts)
    design = build_simulation_trajectory_design(
        analysis._working_frame(dataset), group_col=analysis._GROUP, stage_col=analysis._STAGE
    )
    full, reduced, _ = analysis._covariate_adjusted_design(
        design.model_full, design.model_reduced, design.ls_means, dataset.covariate_design
    )
    plain = interaction_split_test(features, full, reduced, basis, permutations=params.permutations, seed=params.seed)
    block = result.record["span_split"]["parts"]
    for name in PARTS:
        observed = plain.parts[name]
        if (observed.ss_interaction, observed.rss) != (block[name]["ss_interaction"], block[name]["rss"]):
            raise pc.ConsistencyError(f"split_inputs: the rebuilt inputs do not reproduce the {name} part's sums.")
    return plain


def evaluate_plasmode_unit(effect: float, replicate: int, settings: Mapping[str, Any], export: Any) -> dict[str, Any]:
    """One case-study analysis with the split test; seeds and statistics only."""

    start = perf_counter()
    parent = settings["plasmode"]
    base = {
        "profile": "plasmode",
        "key": pc.unit_key(effect, replicate),
        "effect": float(effect),
        "cell": pc.cell_label(effect),
        "replicate": int(replicate),
        "stream": parent["stream"],
        "seeds": pc.unit_seeds(int(parent["base_seed"]), str(parent["stream"]), int(replicate)),
        "settings_hash": lsc.settings_hash(settings),
    }
    try:
        built = pc.build_unit(effect, replicate, parent, export)
        params = plasmode_params(settings, built["seeds"])
        dataset, result = pc.analyse(built["blocks"], built["metadata"], export, params)
        pc.check_consistency(result)
        if float(effect) == pc.ANCHOR:
            pc.check_anchor_latent(result, settings["fidelity"])
        plain = uncorrected_split(dataset, result, params)
    except Exception as exc:  # recorded, never aborts the shard
        return {**base, "status": "failed", "error_type": type(exc).__name__, "error_message": str(exc),
                "runtime_seconds": perf_counter() - start}
    return lsc._jsonable(
        {
            **base,
            "status": "completed",
            "seeds": built["seeds"],
            "label_redraws": built["label_redraws"],
            "analysis": pc.analysis_summary(result),
            "span_split": result.record["span_split"],
            "ss_uncorrected": {name: plain.parts[name].p_value_ss for name in PARTS},
            "runtime_seconds": perf_counter() - start,
        }
    )


# ── Shards, resume and the run ───────────────────────────────────────────────


def read_records(paths: Iterable[Path], key: Callable[[Mapping[str, Any]], str]) -> dict[str, dict[str, Any]]:
    """Records by unit; a line torn by an interrupted write is dropped, so its unit is re-run."""

    out: dict[str, dict[str, Any]] = {}
    for path in paths:
        with Path(path).open(encoding="utf-8") as handle:
            for line in handle:
                try:
                    record = json.loads(line)
                except json.JSONDecodeError:
                    continue
                if record.get("status") in ("completed", "failed"):
                    out[key(record)] = record
    return out


_WORKER_EXPORT: Any = None


def _init_plasmode_worker(directory: str, files: Mapping[str, str]) -> None:
    global _WORKER_EXPORT
    export = pc.load_export(Path(directory))
    if export.files != dict(files):
        raise SystemExit(f"{directory} changed after the run started (file hashes differ).")
    _WORKER_EXPORT = export


def _evaluate_plasmode_pending(unit: tuple[float, int], settings: Mapping[str, Any]) -> dict[str, Any]:
    return evaluate_plasmode_unit(*unit, settings, _WORKER_EXPORT)


def run_shard(
    out_dir: Path,
    settings: Mapping[str, Any],
    index: int,
    count: int,
    *,
    export: Any = None,
    limit: int | None = None,
    workers: int = 1,
) -> int:
    """Evaluate this shard's missing units into its file; returns how many were evaluated."""

    profile = settings["profile"]
    directory = profile_dir(out_dir, profile)
    lsc.ensure_settings(directory, settings)
    path = lsc.shard_path(directory, index, count)
    path.parent.mkdir(parents=True, exist_ok=True)
    if profile == "simulation":
        keyed: list[tuple[str, Any]] = [(_unit_key(unit), unit) for unit in simulation_units(settings)]
        record_key: Callable[[Mapping[str, Any]], str] = _simulation_key
    else:
        if export is None or export.files != settings["plasmode"]["export"]["files"]:
            raise SystemExit("The export's file hashes differ from the settings'; refusing to append.")
        keyed = [(pc.unit_key(*unit), unit) for unit in plasmode_units(settings)]
        record_key = _plasmode_key
    done = read_records([path], record_key) if path.exists() else {}
    # Rewritten with complete lines only, so an interrupted unit is re-run, not duplicated.
    with path.open("w", encoding="utf-8") as handle:
        for key in sorted(done):
            handle.write(lsc._line(done[key]))
    pending = [unit for position, (key, unit) in enumerate(keyed) if position % count == index and key not in done]
    if limit is not None:
        pending = pending[:limit]

    def append(position: int, record: Mapping[str, Any]) -> None:
        with path.open("a", encoding="utf-8") as handle:
            handle.write(lsc._line(record))
            handle.flush()
        print(f"[{profile} {index}/{count}] {position}/{len(pending)} {record_key(record)} {record['status']}",
              file=sys.stderr, flush=True)

    if profile == "simulation":
        if workers <= 1:
            for position, unit in enumerate(pending, start=1):
                append(position, evaluate_simulation_unit(unit, settings))
        else:
            with ProcessPoolExecutor(max_workers=workers) as pool:
                futures = [pool.submit(evaluate_simulation_unit, unit, settings) for unit in pending]
                for position, future in enumerate(as_completed(futures), start=1):
                    append(position, future.result())
    elif workers <= 1:
        for position, unit in enumerate(pending, start=1):
            append(position, evaluate_plasmode_unit(*unit, settings, export))
    else:
        with ProcessPoolExecutor(
            max_workers=workers, initializer=_init_plasmode_worker, initargs=(str(export.directory), dict(export.files))
        ) as pool:
            for position, record in enumerate(
                pool.map(_evaluate_plasmode_pending, pending, [settings] * len(pending)), start=1
            ):
                append(position, record)
    return len(pending)


# ── Merge ────────────────────────────────────────────────────────────────────


SIMULATION_SPLIT_COLUMNS = (
    "design", "cell", "trajectory_mode", "effect_size", "span_alignment", "span_support", "part", "statistic",
    "rate", "se", "count", "rejections", "failures", "span_share_median",
)
PLASMODE_SPLIT_COLUMNS = (
    "cell", "effect", "part", "statistic", "rate", "se", "count", "rejections", "failures", "span_share_median",
)
PRODUCTION_COLUMNS = ("design", "cell", "statistic", "rate", "se", "count", "rejections", "failures")
TRAJECTORY_COLUMNS = ("cell", "effect", "statistic", "rate", "se", "count", "rejections", "failures")
CAPTURE_COLUMNS = (
    "design", "cell", "n", "rank_median", "cv_mean_auroc_median", "group_difference_median", "pooled_stage_median",
    "span_share_median",
)

_P_FIELD = {"f": "p_value_f", "ss": "p_value_ss", "ss_df": "p_value_ss_df"}


def _simulation_split_p(record: Mapping[str, Any], part: str, statistic: str) -> Any:
    return ((record.get("split") or {}).get("parts", {}).get(part) or {}).get(_P_FIELD[statistic])


def _plasmode_split_p(record: Mapping[str, Any], part: str, statistic: str) -> Any:
    if statistic == "ss":
        return (record.get("ss_uncorrected") or {}).get(part)
    # The analysis applies ``df``: its SS p-value is the corrected one.
    field = "p_value_f" if statistic == "f" else "p_value_ss"
    return ((record.get("span_split") or {}).get("parts", {}).get(part) or {}).get(field)


def _span_share(record: Mapping[str, Any]) -> Any:
    block = record.get("split") if record.get("profile") == "simulation" else record.get("span_split")
    return (block or {}).get("span_share")


def _grouped(records: Iterable[Mapping[str, Any]], key: Callable[[Mapping[str, Any]], tuple]) -> dict[tuple, list]:
    groups: dict[tuple, list] = {}
    for record in records:
        groups.setdefault(key(record), []).append(record)
    return groups


def simulation_tables(records: Sequence[Mapping[str, Any]], settings: Mapping[str, Any]) -> dict[str, list]:
    """The split, production and capture tables; failed units are counted by cell, not dropped."""

    alpha, rule = float(settings["alpha"]), str(settings["rejection"])
    split_rows, production_rows, capture_rows = [], [], []
    for (design, cell), group in sorted(_grouped(records, lambda r: (str(r["design"]), str(r["cell"]))).items()):
        completed = [r for r in group if r.get("status") == "completed"]
        failures = len(group) - len(completed)
        first = group[0]
        share = _median(_span_share(r) for r in completed)
        for part in PARTS:
            for statistic in SPLIT_STATISTICS:
                split_rows.append(
                    {
                        "design": design,
                        "cell": cell,
                        "trajectory_mode": first.get("trajectory_mode"),
                        "effect_size": first.get("effect_size"),
                        "span_alignment": first.get("span_alignment"),
                        "span_support": first.get("span_support"),
                        "part": part,
                        "statistic": statistic,
                        **_rate_row(
                            [_rejects(_simulation_split_p(r, part, statistic), alpha, rule) for r in completed]
                        ),
                        "failures": failures,
                        "span_share_median": share,
                    }
                )
        for statistic in STATISTICS:
            production_rows.append(
                {
                    "design": design,
                    "cell": cell,
                    "statistic": statistic,
                    **_rate_row(
                        [_rejects(r["production"]["p_values"].get(statistic), alpha, rule) for r in completed]
                    ),
                    "failures": failures,
                }
            )
        captures = [r["production"].get("span_capture") or {} for r in completed]
        capture_rows.append(
            {
                "design": design,
                "cell": cell,
                "n": len(completed),
                "rank_median": _median(r["production"]["rank"] for r in completed),
                "cv_mean_auroc_median": _median(r["production"].get("cv_mean_auroc") for r in completed),
                "group_difference_median": _median(c.get("group_difference") for c in captures),
                "pooled_stage_median": _median(c.get("pooled_stage") for c in captures),
                "span_share_median": share,
            }
        )
    return {"split": split_rows, "production": production_rows, "capture": capture_rows}


def plasmode_tables(records: Sequence[Mapping[str, Any]], settings: Mapping[str, Any]) -> dict[str, list]:
    alpha, rule = float(settings["alpha"]), str(settings["rejection"])
    split_rows, trajectory_rows = [], []
    groups = _grouped(records, lambda r: (float(r["effect"]),))
    for effect in settings["cells"]:
        group = groups.get((float(effect),), [])
        completed = [r for r in group if r.get("status") == "completed"]
        failures = len(group) - len(completed)
        share = _median(_span_share(r) for r in completed)
        label = pc.cell_label(float(effect))
        for part in PARTS:
            for statistic in SPLIT_STATISTICS:
                split_rows.append(
                    {
                        "cell": label,
                        "effect": float(effect),
                        "part": part,
                        "statistic": statistic,
                        **_rate_row([_rejects(_plasmode_split_p(r, part, statistic), alpha, rule) for r in completed]),
                        "failures": failures,
                        "span_share_median": share,
                    }
                )
        for statistic in STATISTICS:
            trajectory_rows.append(
                {
                    "cell": label,
                    "effect": float(effect),
                    "statistic": statistic,
                    **_rate_row(
                        [_rejects(r["analysis"]["p_values"].get(statistic), alpha, rule) for r in completed]
                    ),
                    "failures": failures,
                }
            )
    return {"split": split_rows, "trajectory": trajectory_rows}


def simulation_reproduction(
    records: Sequence[Mapping[str, Any]], settings: Mapping[str, Any], committed: Sequence[Mapping[str, Any]]
) -> dict[str, Any]:
    """``pls_cv``'s rejections from the first draws of the shared units against the comparison's rate table."""

    shared = int(settings["parent"]["shared_replicates"])
    pseudo = [
        {
            "design": r["design"],
            "cell": r["cell"],
            "arm": lsc.REFERENCE_ARM,
            "status": r.get("status"),
            "p_values": (r.get("production") or {}).get("p_values_reproduction") or {},
        }
        for r in records
        if int(r["replicate"]) < shared
    ]
    reproduced = {
        (row["design"], row["cell"], row["statistic"]): row
        for row in lsc.rate_table(pseudo, alpha=float(settings["alpha"]))
    }
    designs = set(settings["comparison"]["designs"])
    differing = []
    compared = 0
    for row in committed:
        if row["arm"] != lsc.REFERENCE_ARM or row["design"] not in designs:
            continue
        compared += 1
        mine = reproduced.get((row["design"], row["cell"], row["statistic"]))
        expected = (int(row["count"]), int(row["rejections"]))
        observed = None if mine is None else (int(mine["count"]), int(mine["rejections"]))
        if observed != expected:
            differing.append(
                {
                    "design": row["design"],
                    "cell": row["cell"],
                    "statistic": row["statistic"],
                    "committed": {"count": expected[0], "rejections": expected[1]},
                    "reproduced": None if observed is None else {"count": observed[0], "rejections": observed[1]},
                }
            )
    return {
        "profile": "simulation",
        "passed": compared > 0 and not differing,
        "rows_compared": compared,
        "differing": differing,
        "permutations": int(settings["reproduction_permutations"]),
        "shared_replicates": shared,
        "source": settings["parent"]["dir"] + "/rate_table.csv",
        "source_sha256": settings["parent"]["rate_table_sha256"],
    }


def plasmode_reproduction(
    records: Sequence[Mapping[str, Any]], settings: Mapping[str, Any], committed: Sequence[Mapping[str, Any]]
) -> dict[str, Any]:
    """The trajectory rejection rates per cell against the plasmode's committed summary, exactly."""

    alpha = float(settings["alpha"])
    index = {
        (float(row["effect"]), row["quantity"]): row
        for row in committed
        if row["arm"] == "analysis" and row["subset"] == "all" and row["quantity"].startswith("reject_")
    }
    rows, compared = [], 0
    for effect in settings["cells"]:
        completed = [r for r in records if float(r["effect"]) == float(effect) and r.get("status") == "completed"]
        for statistic in STATISTICS:
            values = [
                float(pc._rejects(r["analysis"]["p_values"].get(statistic), alpha))
                for r in completed
                if r["analysis"]["p_values"].get(statistic) is not None
            ]
            n, mean, _ = pc._mean_se(values)
            row = index.get((float(effect), f"reject_{statistic}"))
            compared += 1
            committed_rate = None if row is None else float(row["mean"])
            committed_n = None if row is None else int(row["n"])
            rows.append(
                {
                    "cell": pc.cell_label(float(effect)),
                    "statistic": statistic,
                    "committed": {"rate": committed_rate, "n": committed_n},
                    "reproduced": {"rate": mean, "n": n},
                    "matches": row is not None and committed_n == n and committed_rate == mean,
                }
            )
    return {
        "profile": "plasmode",
        "passed": compared > 0 and all(row["matches"] for row in rows),
        "rows_compared": compared,
        "rows": rows,
        "differing": [row for row in rows if not row["matches"]],
        "source": settings["parent"]["dir"] + "/plasmode_summary.csv",
        "source_sha256": settings["parent"]["summary_sha256"],
    }


def accounting(
    records: Sequence[Mapping[str, Any]], settings: Mapping[str, Any], *, expected: int, duplicates: int
) -> dict[str, Any]:
    profile = settings["profile"]
    failures: dict[str, int] = {}
    unit_seconds, stages = [], {}
    for record in records:
        label = f"{record['design']}/{record['cell']}" if profile == "simulation" else str(record["cell"])
        if record.get("status") != "completed":
            failures[label] = failures.get(label, 0) + 1
        runtime = record.get("runtime_seconds")
        if isinstance(runtime, Mapping):
            unit_seconds.append(float(runtime.get("unit") or 0.0))
            for stage in ("fit", "evaluation", "split"):
                if runtime.get(stage) is not None:
                    stages.setdefault(stage, []).append(float(runtime[stage]))
        elif runtime is not None:
            unit_seconds.append(float(runtime))
    return {
        "profile": profile,
        "settings_hash": lsc.settings_hash(settings),
        "expected_units": expected,
        "complete_units": len(records),
        "missing_units": expected - len(records),
        "failed_units": sum(failures.values()),
        "failures_by_cell": failures,
        "duplicate_records": duplicates,
        "recorded_core_hours": sum(unit_seconds) / 3600.0,
        "unit_seconds": {
            "median": _median(unit_seconds),
            "q90": float(np.quantile(unit_seconds, 0.9)) if unit_seconds else None,
            "max": max(unit_seconds) if unit_seconds else None,
        },
        "stage_seconds_median": {stage: _median(values) for stage, values in sorted(stages.items())},
    }


def merge(out_dir: Path, profile: str) -> dict[str, Any]:
    """Merge one profile's shards and write its tables, accounting and reproduction check."""

    directory = profile_dir(out_dir, profile)
    settings = lsc.load_settings(directory)
    if profile == "simulation":
        expected_keys = [_unit_key(unit) for unit in simulation_units(settings)]
        record_key: Callable[[Mapping[str, Any]], str] = _simulation_key
    else:
        expected_keys = [pc.unit_key(*unit) for unit in plasmode_units(settings)]
        record_key = _plasmode_key
    paths = sorted((directory / "records").glob("shard_*.jsonl"))
    raw: list[dict[str, Any]] = []
    for path in paths:
        raw.extend(read_records([path], record_key).values())
    stale = sorted({record_key(r) for r in raw if r.get("settings_hash") != lsc.settings_hash(settings)})
    if stale:
        raise SystemExit(f"{len(stale)} record(s) were written under other settings, e.g. {stale[:3]}; refusing.")
    seen: dict[str, int] = {}
    for record in raw:
        seen[record_key(record)] = seen.get(record_key(record), 0) + 1
    by_key = {record_key(r): r for r in raw}
    records = [by_key[key] for key in expected_keys if key in by_key]
    with (directory / "merged.jsonl").open("w", encoding="utf-8") as handle:
        handle.write("".join(lsc._line(record) for record in records))
    summary = accounting(
        records, settings, expected=len(expected_keys), duplicates=sum(count - 1 for count in seen.values())
    )
    parent = _resolve(settings["parent"]["dir"])
    if profile == "simulation":
        tables = simulation_tables(records, settings)
        _write_csv(directory / "split_table.csv", tables["split"], SIMULATION_SPLIT_COLUMNS)
        _write_csv(directory / "production_table.csv", tables["production"], PRODUCTION_COLUMNS)
        _write_csv(directory / "capture_table.csv", tables["capture"], CAPTURE_COLUMNS)
        source = parent / "rate_table.csv"
        if sha256_file(source) != settings["parent"]["rate_table_sha256"]:
            raise SystemExit(f"{source} changed since the settings were written (sha256 differs).")
        reproduction = simulation_reproduction(records, settings, lsc._read_csv(source))
    else:
        tables = plasmode_tables(records, settings)
        _write_csv(directory / "split_table.csv", tables["split"], PLASMODE_SPLIT_COLUMNS)
        _write_csv(directory / "trajectory_table.csv", tables["trajectory"], TRAJECTORY_COLUMNS)
        source = parent / "plasmode_summary.csv"
        if sha256_file(source) != settings["parent"]["summary_sha256"]:
            raise SystemExit(f"{source} changed since the settings were written (sha256 differs).")
        reproduction = plasmode_reproduction(records, settings, lsc._read_csv(source))
    _write_json(directory / "accounting.json", summary)
    _write_json(directory / "reproduction.json", reproduction)
    return {
        "profile": profile,
        **{key: summary[key] for key in ("complete_units", "missing_units", "failed_units", "recorded_core_hours")},
        "reproduction_passed": reproduction["passed"],
        "reproduction_differing": len(reproduction["differing"]),
    }


# ── Decision (pure functions) ────────────────────────────────────────────────


def _number(value: Any) -> float | None:
    if value is None or value == "":
        return None
    return float(value)


def _rate(row: Mapping[str, Any] | None) -> dict[str, Any] | None:
    if row is None or _number(row.get("rate")) is None:
        return None
    return {"rate": float(row["rate"]), "se": float(row["se"]), "count": int(row["count"])}


def _size_observation(row: Mapping[str, Any] | None, **cell: Any) -> dict[str, Any]:
    if row is None or int(row.get("count") or 0) <= 0:
        return {**cell, "status": "unavailable"}
    return {**cell, "n_rejected": int(row["rejections"]), "n_available": int(row["count"])}


def _failing_observations(judged: Mapping[str, Any]) -> list[dict[str, Any]]:
    return [
        {key: value for key, value in row.items() if key not in ("bound", "ceiling")}
        for row in judged["observations"]
        if row.get("status") != "holds"
    ]


def _grade_label(settings: Mapping[str, Any], count: int | None) -> str:
    return f"{count if count is not None else '?'} × {settings['permutations']}"


def decision(
    simulation: Mapping[str, Any],
    plasmode: Mapping[str, Any],
    *,
    comparison_rates: Sequence[Mapping[str, Any]] = (),
) -> dict[str, Any]:
    """Apply rule version 1 to the merged tables of both profiles.

    ``simulation`` and ``plasmode`` hold ``settings``, ``split`` (rows),
    ``reproduction`` and ``accounting``; the simulation also ``production`` and
    ``capture``, the plasmode ``trajectory``. ``comparison_rates`` is the
    comparison's committed rate table, read only by a descriptive table.
    """

    failed = [
        {"profile": name, "differing": block["reproduction"].get("differing", [])}
        for name, block in (("simulation", simulation), ("plasmode", plasmode))
        if not block["reproduction"].get("passed")
    ]
    header = {
        "rule_version": RULE_VERSION,
        "alpha": ALPHA,
        "se_multiplier": SE_MULTIPLIER,
        "rejection": {name: block["settings"]["rejection"] for name, block in (("simulation", simulation),
                                                                              ("plasmode", plasmode))},
        "settings_hashes": {
            "simulation": lsc.settings_hash(simulation["settings"]),
            "plasmode": lsc.settings_hash(plasmode["settings"]),
        },
        "complete": all(block["accounting"]["missing_units"] == 0 for block in (simulation, plasmode)),
        "failed_units": {
            "simulation": simulation["accounting"]["failed_units"],
            "plasmode": plasmode["accounting"]["failed_units"],
        },
        "reproduction": {
            "simulation": {k: simulation["reproduction"].get(k) for k in ("passed", "rows_compared", "source")},
            "plasmode": {k: plasmode["reproduction"].get(k) for k in ("passed", "rows_compared", "source")},
        },
    }
    if failed:
        return {**header, "verdict": "reproduction_failed", "consequence": CONSEQUENCES["reproduction_failed"],
                "failing": {"reproduction": failed}}

    sim_split = {(r["design"], r["cell"], r["part"], r["statistic"]): r for r in simulation["split"]}
    production = {(r["design"], r["cell"], r["statistic"]): r for r in simulation["production"]}
    pm_split = {(r["cell"], r["part"], r["statistic"]): r for r in plasmode["split"]}
    designs = list(simulation["settings"]["comparison"]["designs"])

    c1 = grade.judge_rates(
        [
            _size_observation(sim_split.get((design, cell, part, "f")), profile="simulation", design=design,
                              cell=cell, part=part)
            for design in designs
            for cell in CONTROL_CELLS
            for part in PARTS
        ]
    )
    c2 = grade.judge_rates(
        [
            _size_observation(pm_split.get(("anchor", part, "f")), profile="plasmode", cell="anchor", part=part)
            for part in PARTS
        ]
    )
    gain_checks = []
    for design in GAIN_DESIGNS:
        mine = _rate(sim_split.get((design, GAIN_CELL, "complement", "f")))
        reference = _rate(production.get((design, GAIN_CELL, "angle")))
        entry: dict[str, Any] = {"design": design, "cell": GAIN_CELL}
        if mine is None or reference is None:
            gain_checks.append({**entry, "complement": mine, "production_angle": reference, "pass": False,
                                "missing": True})
            continue
        margin = SE_MULTIPLIER * math.sqrt(mine["se"] ** 2 + reference["se"] ** 2)
        difference = mine["rate"] - reference["rate"]
        gain_checks.append(
            {**entry, "complement": mine, "production_angle": reference, "difference": difference,
             "pooled_se": margin / SE_MULTIPLIER, "threshold": margin, "pass": difference > margin}
        )
    c3 = {"pass": all(check["pass"] for check in gain_checks), "checks": gain_checks,
          "failing": [check["design"] for check in gain_checks if not check["pass"]]}
    criteria = {
        "C1": {"name": "simulation size", **c1},
        "C2": {"name": "real-noise size", **c2},
        "C3": {"name": "out-of-span gain", **c3},
    }
    failing: dict[str, Any] = {}
    if not c1["holds"]:
        failing["C1"] = _failing_observations(c1)
    if not c2["holds"]:
        failing["C2"] = _failing_observations(c2)
    if not c3["pass"]:
        failing["C3"] = c3["failing"]
    if "C1" in failing or "C2" in failing:
        verdict = "size_excess"
    elif "C3" in failing:
        verdict = "no_out_of_span_gain"
    else:
        verdict = "split_validated"
    return {
        **header,
        "verdict": verdict,
        "consequence": CONSEQUENCES[verdict],
        "failing": failing,
        "criteria": criteria,
        "descriptive": descriptive_tables(simulation, plasmode, comparison_rates=comparison_rates),
    }


def descriptive_tables(
    simulation: Mapping[str, Any],
    plasmode: Mapping[str, Any],
    *,
    comparison_rates: Sequence[Mapping[str, Any]] = (),
) -> dict[str, Any]:
    """Localization, share against capture, span against production, complement against the comparison,
    plasmode power and the SS diagnostic; each row names its source and grade."""

    sim_settings, pm_settings = simulation["settings"], plasmode["settings"]
    sim_split = {(r["design"], r["cell"], r["part"], r["statistic"]): r for r in simulation["split"]}
    production = {(r["design"], r["cell"], r["statistic"]): r for r in simulation["production"]}
    capture = {(r["design"], r["cell"]): r for r in simulation["capture"]}
    pm_split = {(r["cell"], r["part"], r["statistic"]): r for r in plasmode["split"]}
    trajectory = {(r["cell"], r["statistic"]): r for r in plasmode["trajectory"]}
    comparison = {(r["design"], r["cell"], r["arm"], r["statistic"]): r for r in comparison_rates}
    designs = list(sim_settings["comparison"]["designs"])
    cells = [cell[0] for cell in sim_settings["comparison"]["cells"]]

    def source(settings: Mapping[str, Any], row: Mapping[str, Any] | None) -> dict[str, str]:
        count = None if row is None else int(row.get("count") or 0)
        return {"source": f"this run, {settings['profile']} profile", "grade": _grade_label(settings, count)}

    localization, share, span_production, against = [], [], [], []
    for design in designs:
        for cell in cells:
            parts = {part: _rate(sim_split.get((design, cell, part, "f"))) for part in PARTS}
            anchor_row = sim_split.get((design, cell, "total", "f"))
            localization.append({"design": design, "cell": cell, **parts, **source(sim_settings, anchor_row)})
            captured = capture.get((design, cell)) or {}
            share.append(
                {
                    "design": design,
                    "cell": cell,
                    "span_share_median": _number(captured.get("span_share_median")),
                    "group_difference_median": _number(captured.get("group_difference_median")),
                    "pooled_stage_median": _number(captured.get("pooled_stage_median")),
                    "rank_median": _number(captured.get("rank_median")),
                    **source(sim_settings, anchor_row),
                }
            )
            span_production.append(
                {
                    "design": design,
                    "cell": cell,
                    "span": parts["span"],
                    **{statistic: _rate(production.get((design, cell, statistic))) for statistic in STATISTICS},
                    **source(sim_settings, anchor_row),
                }
            )
            against.append(
                {
                    "design": design,
                    "cell": cell,
                    "complement": {**(parts["complement"] or {}), **source(sim_settings, anchor_row)}
                    if parts["complement"] is not None
                    else None,
                    "comparison_concat_angle": _with_source(comparison.get((design, cell, lsc.CEILING_ARM, "angle"))),
                    "comparison_pls_cv_angle": _with_source(
                        comparison.get((design, cell, lsc.REFERENCE_ARM, "angle"))
                    ),
                }
            )
    power = []
    for effect in POWER_EFFECTS:
        label = pc.cell_label(effect)
        row = pm_split.get((label, "total", "f"))
        power.append(
            {
                "cell": label,
                **{part: _rate(pm_split.get((label, part, "f"))) for part in PARTS},
                **{statistic: _rate(trajectory.get((label, statistic))) for statistic in STATISTICS},
                **source(pm_settings, row),
            }
        )
    diagnostic = []
    for design in designs:
        for cell in CONTROL_CELLS:
            for part in PARTS:
                row = sim_split.get((design, cell, part, "f"))
                diagnostic.append(
                    {
                        "profile": "simulation",
                        "design": design,
                        "cell": cell,
                        "part": part,
                        **{s: _rate(sim_split.get((design, cell, part, s))) for s in SPLIT_STATISTICS},
                        "bound": None if row is None else lsc.size_bound(ALPHA, max(int(row["count"]), 1)),
                        **source(sim_settings, row),
                    }
                )
    for part in PARTS:
        row = pm_split.get(("anchor", part, "f"))
        diagnostic.append(
            {
                "profile": "plasmode",
                "design": "seaad_export",
                "cell": "anchor",
                "part": part,
                **{s: _rate(pm_split.get(("anchor", part, s))) for s in SPLIT_STATISTICS},
                "bound": None if row is None else lsc.size_bound(ALPHA, max(int(row["count"]), 1)),
                **source(pm_settings, row),
            }
        )
    return {
        "localization": localization,
        "share_against_capture": share,
        "span_against_production": span_production,
        "complement_against_comparison": against,
        "plasmode_power": power,
        "ss_diagnostic": diagnostic,
    }


def _with_source(row: Mapping[str, Any] | None) -> dict[str, Any] | None:
    rate = _rate(row)
    return None if rate is None else {**rate, **COMPARISON_SOURCE}


def _profile_tables(directory: Path, profile: str) -> dict[str, Any]:
    settings = lsc.load_settings(directory)
    summary = json.loads((directory / "accounting.json").read_text(encoding="utf-8")) if (
        directory / "accounting.json"
    ).exists() else None
    if summary is None:
        raise SystemExit(f"{directory}/accounting.json does not exist; run `merge --profile {profile}` first.")
    if summary.get("settings_hash") != lsc.settings_hash(settings):
        raise SystemExit(f"{directory}/accounting.json was merged under other settings; run `merge` again.")
    block: dict[str, Any] = {
        "settings": settings,
        "accounting": summary,
        "reproduction": json.loads((directory / "reproduction.json").read_text(encoding="utf-8")),
        "split": _read_csv(directory / "split_table.csv"),
    }
    if profile == "simulation":
        block["production"] = _read_csv(directory / "production_table.csv")
        block["capture"] = _read_csv(directory / "capture_table.csv")
    else:
        block["trajectory"] = _read_csv(directory / "trajectory_table.csv")
    return block


def decide(out_dir: Path) -> dict[str, Any]:
    """Read both profiles' merged tables and write ``span_split_decision.json``."""

    simulation = _profile_tables(profile_dir(out_dir, "simulation"), "simulation")
    plasmode = _profile_tables(profile_dir(out_dir, "plasmode"), "plasmode")
    comparison = _resolve(simulation["settings"]["parent"]["dir"]) / "rate_table.csv"
    result = decision(simulation, plasmode, comparison_rates=lsc._read_csv(comparison))
    _write_json(Path(out_dir) / DECISION_NAME, result)
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
    run.add_argument("--reps", type=int, default=None,
                     help=f"Replicates per cell (simulation: non-control cells, default {SIMULATION_REPLICATES}; "
                          f"plasmode: default {PLASMODE_REPLICATES})")
    run.add_argument("--control-reps", type=int, default=CONTROL_REPLICATES,
                     help="Simulation: replicates of the anchor and translation")
    run.add_argument("--designs", type=str, default=",".join(lsc.DESIGNS), help="Simulation: comma-separated")
    run.add_argument("--parent-dir", type=Path, default=None, help="The committed parent run (default: the profile's)")
    run.add_argument("--calibration", type=Path, default=None,
                     help="Simulation: the comparison's calibration.json (default: the parent's)")
    run.add_argument("--export-dir", type=Path, default=None, help="Plasmode: the SEA-AD cohort_30 export")
    run.add_argument("--fidelity", type=Path, default=None, help="Plasmode: fidelity.json of a passed gate")
    merge_parser = sub.add_parser("merge", help="Merge one profile's shards")
    merge_parser.add_argument("--profile", choices=PROFILES, required=True)
    merge_parser.add_argument("--out-dir", type=Path, required=True)
    decide_parser = sub.add_parser("decide", help="Write span_split_decision.json from both profiles")
    decide_parser.add_argument("--out-dir", type=Path, required=True)
    return parser


def _settings_from_args(args: argparse.Namespace) -> tuple[dict[str, Any], Any]:
    if args.profile == "simulation":
        settings = simulation_settings(
            parent_dir=args.parent_dir or COMPARISON_DIR,
            calibration=args.calibration,
            replicates=SIMULATION_REPLICATES if args.reps is None else args.reps,
            control_replicates=args.control_reps,
            permutations=args.perms,
            designs=[name for name in args.designs.split(",") if name],
        )
        return settings, None
    if args.export_dir is None or args.fidelity is None:
        raise SystemExit("The plasmode profile needs --export-dir and --fidelity.")
    try:
        export = pc.load_export(args.export_dir)
    except pc.ExportContractError as exc:
        raise pc.contract_exit(exc) from exc
    settings = plasmode_settings(
        export,
        pc.check_fidelity(args.fidelity, export),
        parent_dir=args.parent_dir or PLASMODE_DIR,
        replicates=PLASMODE_REPLICATES if args.reps is None else args.reps,
        permutations=args.perms,
    )
    return settings, export


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    if args.command == "run":
        index, count = lsc.parse_shard(args.shard)
        settings, export = _settings_from_args(args)
        evaluated = run_shard(args.out_dir, settings, index, count, export=export, limit=args.limit,
                              workers=args.workers)
        print(f"evaluated {evaluated} {args.profile} unit(s)", file=sys.stderr)
        return 0
    if args.command == "merge":
        print(json.dumps(merge(args.out_dir, args.profile), indent=1))
        return 0
    result = decide(args.out_dir)
    failing = {name: len(entries) for name, entries in result["failing"].items()}
    print(json.dumps({"verdict": result["verdict"], "failing": failing, "complete": result["complete"]}, indent=1))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
