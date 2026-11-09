#!/usr/bin/env python3
"""Per-feature interaction scan validation: size, familywise control and sparse power on the split-test units.

The driver behind ``openspec/changes/interaction-feature-scan/`` (spec
``interaction-scan-validation``, roadmap item 16, thread 3). Every test MOTCO
had pooled the group × stage difference over many dimensions, so none saw a
difference on a few features: the plasmode's 20-feature plant was rejected at
the null rate by ``delta``/``angle``/``shape`` and by all three parts of the
split interaction test. :func:`motco.stats.interaction.interaction_feature_scan`
tests each feature's interaction F on its own and combines the features by
their maximum (Westfall–Young step-down max-T). This driver measures it before
a case study may use it, on the split-test validation's units
(``results/span-split-2026-11-05/``), which it must reproduce.

**Profiles.**

- ``simulation`` — the split-test validation's simulation units: the
  latent-space comparison's three designs (``n200``, ``seaad``, ``seaad_weak``,
  498 features), fourteen cells, base seed 2035 and seed derivations, rebuilt
  through ``span_split_check.simulation_settings``. At the parent's grade the
  rebuilt settings must equal the committed ``simulation/settings.json``; the
  sha256 of that file and of the parent's ``split_table.csv`` and
  ``production_table.csv`` are recorded, and a changed file is refused by name
  in ``run``, ``merge`` and ``decide``. The anchor and translation run
  replicates 0–499, every other cell 0–99, at 999 permutations. Each unit runs
  the parent unit unchanged (the comparison's dataset, one ``pls_cv`` fit and
  evaluation, the split test uncorrected and under ``df``) and then the scan on
  the split test's own features (the fit's standardized features in the PLS
  input coordinates, a per-feature affine map the scan ignores) and models, at
  the unit's RRPP seed. It records the scan's counts, maximum, global p-value,
  hits per block, the per-block share of null maxima and, from the dataset's
  population stage means, the hits on population-null features (whose
  population B − A difference is constant across stages within 1e-9 of their
  largest absolute population mean).
- ``plasmode`` — the split-test validation's plasmode units: the SEA-AD
  astrocyte ``cohort_30`` export under the plasmode's contract (its file
  hashes must equal the plasmode's committed settings; ``--fidelity`` must be
  a passed gate on it), base seed 2034, ``paper`` stream, the anchor, e = 2.0
  and e = 3.0 × 500 replicates × 999 permutations, rebuilt through
  ``span_split_check.plasmode_settings`` and checked against the committed
  ``plasmode/settings.json`` the same way. Each unit is one
  ``analyze_trajectories`` at the plasmode's case-study parameters (six
  covariates, rank floor 2, net displacement, ``df``, harness CV) with the
  split test on, the scan on and attribution off, then the plasmode's
  consistency check and, at the anchor, its latent check. The planted set's
  counts (hits inside and outside it, the planted features among the 20
  largest F, in total and per block) are computed from the unit's plant and
  recorded as counts. A record holds seeds, statistics and counts only, never
  a donor identifier, label vector, feature list or donor-level value.

A unit that raises is recorded as failed with its error; the shard continues.
The scan's rejection rule is adjusted p ≤ 0.05 in both profiles (exact at
B = 999); the reproduction uses each parent's own rule (p < α in the
simulation, p ≤ α in the plasmode).

Subcommands (every output lands under ``--out-dir``/``<profile>/``):

- ``run --profile simulation|plasmode`` — evaluate (a shard of) the units into
  ``records/shard_<i>_of_<N>.jsonl``. The settings and their hash go to
  ``settings.json``; a resume or merge under other settings exits naming the
  differing keys, and a resume skips complete units.
- ``merge --profile …`` — ``merged.jsonl`` (per unit; keep off the
  repository), ``scan_table.csv`` (the scan's global rejection rate, Monte
  Carlo SE, count and median hit count per design/cell), ``block_balance.csv``
  (per block the mean share of null maxima beside its share of features),
  ``comparator_table.csv`` (the scan and the six existing tests on the same
  units), ``false_hit_table.csv`` (simulation: units with a hit on a
  population-null feature) or ``recovery_table.csv`` (plasmode planted cells:
  units with a hit outside the planted set, mean recall, hit precision and
  top-20 precision), ``accounting.json`` and ``reproduction.json`` (the split
  test's pseudo-F rows and the trajectory statistics' rows against the
  parent's committed tables, exact in count and rejections; differing rows
  listed).
- ``decide`` — ``interaction_scan_decision.json`` (rule version 1): the gate
  ``reproduction_failed`` (no criterion read); S1 simulation size of the
  scan's global test at the anchor and translation of every design (m = 6),
  S2 real-noise size at the plasmode anchor (m = 1) and S3 strong familywise
  control (units with a hit outside the planted set at e = 2.0 and 3.0, m = 2),
  all under the stage-span paper grade's count-plus-ceiling rule (imported);
  G the scan's global rate at e = 3.0 above the largest of ``delta``,
  ``angle``, ``shape`` and the split test's three parts on the same units by
  more than 2 pooled SE. Verdict: ``size_excess`` | ``fwer_excess`` |
  ``no_sparse_gain`` | ``scan_validated``, with its consequence and
  descriptive tables (plasmode power and recovery beside attribution's
  committed planted-set precision, block balance, simulation localization,
  simulation false hits) that never change it.

With ``--workers`` above 1, pin BLAS to one thread per process
(``OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1``).

Examples:

    export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
    EXPORT=../sead_sexdiff/results/motco_input/cohort_30

    # smoke
    python scripts/interaction_scan_check.py run --profile simulation --reps 2 --control-reps 2 --perms 19 \\
        --out-dir /tmp/isc --workers 7
    python scripts/plasmode_check.py fidelity --export-dir $EXPORT --out-dir /tmp/isc/fidelity
    python scripts/interaction_scan_check.py run --profile plasmode --export-dir $EXPORT \\
        --fidelity /tmp/isc/fidelity/fidelity.json --reps 2 --perms 19 --out-dir /tmp/isc --workers 6
    python scripts/interaction_scan_check.py merge --profile simulation --out-dir /tmp/isc
    python scripts/interaction_scan_check.py merge --profile plasmode --out-dir /tmp/isc
    python scripts/interaction_scan_check.py decide --out-dir /tmp/isc

    # one SLURM array task (scripts/interaction_scan_check_array.sbatch), then merge and decide
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
import latent_space_comparison as lsc  # noqa: E402  (units, seeds, settings binding, CSV helpers)
import plasmode_check as pc  # noqa: E402  (export contract, units, planted set)
import span_split_check as ssc  # noqa: E402  (the parent run's settings and unit steps)
import stage_span_paper_grade_decision as grade  # noqa: E402  (the count-plus-ceiling rule)

sys.path.pop(0)

DRIVER_VERSION = 1
RULE_VERSION = 1
PROFILES: tuple[str, ...] = ("simulation", "plasmode")
PARTS: tuple[str, ...] = ssc.PARTS
STATISTICS: tuple[str, ...] = ssc.STATISTICS
ALPHA = 0.05
SE_MULTIPLIER = 2.0
SCAN_RULE = "p <= alpha"
TOP_K = 20
POPULATION_NULL_TOLERANCE = 1e-9

CONTROL_CELLS: tuple[str, ...] = ssc.CONTROL_CELLS
PLANTED_EFFECTS: tuple[float, ...] = (2.0, 3.0)
GAIN_EFFECT = 3.0
#: The six existing tests G compares the scan with, in the order a tie is broken.
EXISTING_TESTS: tuple[str, ...] = ("delta", "angle", "shape", "split_total", "split_span", "split_complement")
COMPARATOR_TESTS: tuple[str, ...] = ("scan", *EXISTING_TESTS)

REPO = Path(__file__).resolve().parents[1]
SPAN_SPLIT_DIR = REPO / "results" / "span-split-2026-11-05"
#: The parent's committed files each profile binds by sha256.
PARENT_FILES: dict[str, tuple[str, ...]] = {
    "simulation": ("settings.json", "split_table.csv", "production_table.csv"),
    "plasmode": ("settings.json", "split_table.csv", "trajectory_table.csv"),
}
PARENT_SOURCE = "split-test validation 2026-11-05"
ATTRIBUTION_SOURCE = "plasmode 2026-11-02"

VERDICTS: tuple[str, ...] = ("reproduction_failed", "size_excess", "fwer_excess", "no_sparse_gain", "scan_validated")
_PHASE6 = (
    " No verdict reopens Phase 6 or changes its attribution leg: the plasmode's attribution_leg_failed stands, and "
    "the scan is a new test with its own claim."
)
CONSEQUENCES: dict[str, str] = {
    "reproduction_failed": (
        "Nothing is read: the units do not reproduce the split-test validation, so no rate pairs with the split "
        "test or the trajectory statistics. Find the drift before any rerun." + _PHASE6
    ),
    "size_excess": (
        "A case study may not report the per-feature interaction scan: the size of its global test is not "
        "established at the failing cells." + _PHASE6
    ),
    "fwer_excess": (
        "A case study may report at most the scan's global test, not the features it names: features outside the "
        "planted set were hits more often than familywise control allows under the plasmode's partial null."
        + _PHASE6
    ),
    "no_sparse_gain": (
        "A case study may not report the scan as seeing what the other tests miss: at the plasmode's e = 3.0 it "
        "does not reject more often than the best existing test by more than 2 pooled SE." + _PHASE6
    ),
    "scan_validated": (
        "A case study may report the per-feature interaction scan, read as the measurement-visibility "
        "documentation states: a hit names a feature whose covariate-adjusted group difference changes across "
        "stages, with the familywise error rate over the scanned features controlled at 0.05; a scan without hits "
        "is uninformative without power at the design and effect." + _PHASE6
    ),
}

DECISION_NAME = "interaction_scan_decision.json"


# ── Shared helpers ───────────────────────────────────────────────────────────


def _canonical(value: Any) -> Any:
    return json.loads(json.dumps(value, sort_keys=True))


def profile_dir(out_dir: Path, profile: str) -> Path:
    return Path(out_dir) / profile


def _rate_row(rejections: Sequence[bool | None]) -> dict[str, Any]:
    return ssc._rate_row(rejections)


def _rejects(p_value: Any, rule: str = SCAN_RULE) -> bool | None:
    return ssc._rejects(p_value, ALPHA, rule)


def _median(values: Iterable[Any]) -> float | None:
    return ssc._median(values)


def _mean_se(values: Iterable[Any]) -> tuple[int, float | None, float | None]:
    finite = [float(value) for value in values if value is not None and math.isfinite(float(value))]
    if not finite:
        return 0, None, None
    mean = float(np.mean(finite))
    se = float(np.std(finite, ddof=1) / math.sqrt(len(finite))) if len(finite) > 1 else None
    return len(finite), mean, se


def _number(value: Any) -> float | None:
    return ssc._number(value)


def _write_csv(path: Path, rows: Sequence[Mapping[str, Any]], columns: Sequence[str]) -> None:
    lsc._write_csv(path, rows, columns)


def _write_json(path: Path, payload: Any) -> None:
    lsc._write_json(path, payload)


def _read_csv(path: Path) -> list[dict[str, Any]]:
    return ssc._read_csv(path)


def _grouped(records: Iterable[Mapping[str, Any]], key: Callable[[Mapping[str, Any]], tuple]) -> dict[tuple, list]:
    return ssc._grouped(records, key)


# ── Parent binding ───────────────────────────────────────────────────────────


def _committed(parent_dir: Path) -> dict[str, Any]:
    return ssc.load_parent(parent_dir)


def _require_equal(rebuilt: Mapping[str, Any], committed: Mapping[str, Any], path: Path) -> None:
    if _canonical(dict(rebuilt)) != _canonical(dict(committed)):
        differing = lsc.settings_difference(_canonical(dict(committed)), _canonical(dict(rebuilt)))
        raise SystemExit(
            f"{path} differs from the settings span_split_check rebuilds at its grade; differing: "
            f"{', '.join(differing)}. The units would not be the split-test validation's."
        )


def _parent_block(parent_dir: Path, profile: str, committed: Mapping[str, Any]) -> dict[str, Any]:
    return {
        "dir": ssc._display_path(parent_dir),
        "settings_hash": committed["hash"],
        "files": {name: ssc.sha256_file(Path(parent_dir) / name) for name in PARENT_FILES[profile]},
    }


def check_parent_files(settings: Mapping[str, Any]) -> None:
    """Refuse, by name, a committed parent file whose sha256 differs from the one the settings record."""

    directory = ssc._resolve(settings["parent"]["dir"])
    for name, digest in settings["parent"]["files"].items():
        path = directory / name
        if not path.exists():
            raise SystemExit(f"{path} does not exist; the parent run's committed file is required.")
        if ssc.sha256_file(path) != digest:
            raise SystemExit(f"{path} has changed since the settings were written (sha256 differs); refusing.")


def _scan_description() -> dict[str, Any]:
    return {
        "function": "motco.stats.interaction.interaction_feature_scan",
        "statistic": "per-feature interaction F (SS_int / df_int) / (RSS / df_res)",
        "combination": "maximum F; Westfall-Young step-down max-T adjusted p-values",
        "hit_alpha": ALPHA,
        "rejection": SCAN_RULE,
        "top_k": TOP_K,
    }


# ── Simulation profile: settings and units ───────────────────────────────────


def simulation_settings(
    *,
    parent_dir: Path = SPAN_SPLIT_DIR / "simulation",
    comparison_dir: Path = ssc.COMPARISON_DIR,
    replicates: int = ssc.SIMULATION_REPLICATES,
    control_replicates: int = ssc.CONTROL_REPLICATES,
    permutations: int = ssc.DEFAULT_PERMUTATIONS,
    designs: Sequence[str] = lsc.DESIGNS,
) -> dict[str, Any]:
    """The split-test validation's simulation settings at this run's grade, its committed ones checked."""

    committed = _committed(parent_dir)
    stored = committed["settings"]
    rebuilt = ssc.simulation_settings(
        parent_dir=comparison_dir,
        replicates=int(stored["replicates"]["other"]),
        control_replicates=int(stored["replicates"]["control"]),
        permutations=int(stored["permutations"]),
        designs=list(stored["comparison"]["designs"]),
    )
    _require_equal(rebuilt, stored, Path(parent_dir) / "settings.json")
    span_split = ssc.simulation_settings(
        parent_dir=comparison_dir,
        replicates=replicates,
        control_replicates=control_replicates,
        permutations=permutations,
        designs=designs,
    )
    return {
        "driver_version": DRIVER_VERSION,
        "profile": "simulation",
        "parent": _parent_block(parent_dir, "simulation", committed),
        "span_split": span_split,
        "replicates": dict(span_split["replicates"]),
        "permutations": int(permutations),
        "alpha": ALPHA,
        "rejection": {"scan": SCAN_RULE, "parent": ssc.REJECTION_RULES["simulation"]},
        "scan": {
            **_scan_description(),
            "features": "the split test's features: standardized features in the PLS input coordinates",
            "models": "full group x stage and reduced group + stage trajectory models",
            "seed": "the unit's RRPP seed",
            "residual_correction": None,
        },
        "population_null_tolerance": POPULATION_NULL_TOLERANCE,
    }


def simulation_units(settings: Mapping[str, Any]) -> list[lsc.Unit]:
    return ssc.simulation_units(settings["span_split"])


# ── Simulation profile: one unit ─────────────────────────────────────────────


def _block_labels(columns: Sequence[str]) -> np.ndarray:
    return np.asarray([str(column).split("__", 1)[0] for column in columns], dtype=object)


def scan_summary(scan: Any, blocks: np.ndarray) -> dict[str, Any]:
    """The scan's counts, maximum, global p-value, hits and per-block null-maximum shares; no feature names."""

    names = list(dict.fromkeys(blocks.tolist()))
    permuted = scan.p_value_adjusted is not None
    hits = scan.hits(ALPHA) if permuted else np.array([], dtype=int)
    per_block: dict[str, dict[str, Any]] = {}
    for name in names:
        member = blocks == name
        entry: dict[str, Any] = {
            "n_features": int(np.count_nonzero(member)),
            "n_degenerate": int(np.count_nonzero(scan.degenerate[member])),
        }
        if permuted:
            entry["hits"] = int(np.count_nonzero(member[hits]))
            entry["null_max_share"] = float(np.mean(blocks[scan.null_max_column] == name))
        per_block[name] = entry
    return {
        "n_features": scan.n_features,
        "n_degenerate": scan.n_degenerate,
        "df_interaction": scan.df_interaction,
        "df_residual": scan.df_residual,
        "max_f": scan.max_f,
        "max_block": str(blocks[scan.max_column]),
        "p_value_global": scan.p_value_global,
        "hits": int(len(hits)),
        "blocks": per_block,
    }


def population_null_mask(dataset: Any, columns: Sequence[str], group_col: str, tolerance: float) -> np.ndarray | None:
    """Features whose population B − A difference is constant across stages, or ``None`` without population means.

    The difference is taken in integration units from ``dataset.population_trajectories`` (the second sorted
    group minus the first) and is constant when its range over stages is at most ``tolerance`` times the
    feature's largest absolute population mean.
    """

    population = getattr(dataset, "population_trajectories", None)
    if population is None:
        return None
    groups = sorted(str(value) for value in dataset.metadata[group_col].astype(str).unique())
    split = [str(column).split("__", 1) for column in columns]
    mask = np.zeros(len(columns), dtype=bool)
    for layer in dict.fromkeys(parts[0] for parts in split):
        positions = [index for index, parts in enumerate(split) if parts[0] == layer]
        frame = population.layers[layer][[split[index][1] for index in positions]]
        first = frame.xs(groups[0], level=0)
        second = frame.xs(groups[1], level=0).loc[first.index]
        difference = second.to_numpy(dtype=float) - first.to_numpy(dtype=float)
        scale = np.abs(frame.to_numpy(dtype=float)).max(axis=0)
        mask[positions] = np.ptp(difference, axis=0) <= tolerance * scale
    return mask


def evaluate_simulation_unit(unit: lsc.Unit, settings: Mapping[str, Any]) -> dict[str, Any]:
    """The split-test validation's unit unchanged, then the scan on the split test's inputs at its seed."""

    from motco.simulations.evaluation import evaluate_semisynthetic_trajectory, integrate_semisynthetic_dataset
    from motco.stats.interaction import interaction_feature_scan

    start = perf_counter()
    span_split = settings["span_split"]
    comparison = span_split["comparison"]
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
    permutations = int(settings["permutations"])
    try:
        dataset = lsc.unit_dataset(unit, comparison)
        n_stages = int(comparison["designs"][unit.design]["generator"]["n_stages"])
        params = ssc.production_params(span_split, seed=seed, n_stages=n_stages)
        mark = perf_counter()
        latent = integrate_semisynthetic_dataset(dataset, params)
        timing["fit"] = perf_counter() - mark
        mark = perf_counter()
        result = evaluate_semisynthetic_trajectory(dataset, params, latent=latent)
        timing["evaluation"] = perf_counter() - mark
        production = ssc.production_record(result, int(span_split["reproduction_permutations"]))
        features, basis, full, reduced = ssc.split_inputs(dataset, latent, params)
        mark = perf_counter()
        split = ssc.split_record(features, basis, full, reduced, permutations=permutations, seed=seed)
        timing["split"] = perf_counter() - mark
        mark = perf_counter()
        scan = interaction_feature_scan(features, full, reduced, permutations=permutations, seed=seed)
        timing["scan"] = perf_counter() - mark
        columns = [str(column) for column in latent.artifacts.features.columns]
        summary = scan_summary(scan, _block_labels(columns))
        null_mask = population_null_mask(
            dataset, columns, params.group_col, float(settings["population_null_tolerance"])
        )
        hits = scan.hits(ALPHA) if scan.p_value_adjusted is not None else np.array([], dtype=int)
        summary["n_population_null"] = None if null_mask is None else int(null_mask.sum())
        summary["hits_population_null"] = None if null_mask is None else int(null_mask[hits].sum())
    except Exception as exc:  # recorded, never aborts the shard
        timing["unit"] = perf_counter() - start
        return {**base, "status": "failed", "error_type": type(exc).__name__, "error_message": str(exc),
                "runtime_seconds": timing}
    timing["unit"] = perf_counter() - start
    return {**base, "status": "completed", "production": production, "split": split, "scan": summary,
            "runtime_seconds": timing}


# ── Plasmode profile: settings, units and one unit ───────────────────────────


def plasmode_settings(
    export: Any,
    fidelity: Mapping[str, Any],
    *,
    parent_dir: Path = SPAN_SPLIT_DIR / "plasmode",
    plasmode_dir: Path = ssc.PLASMODE_DIR,
    replicates: int = ssc.PLASMODE_REPLICATES,
    permutations: int = ssc.DEFAULT_PERMUTATIONS,
    cells: Sequence[float] = ssc.PLASMODE_CELLS,
) -> dict[str, Any]:
    """The split-test validation's plasmode settings at this run's grade, its committed ones checked."""

    committed = _committed(parent_dir)
    stored = committed["settings"]
    rebuilt = ssc.plasmode_settings(
        export,
        fidelity,
        parent_dir=plasmode_dir,
        replicates=int(stored["replicates"]),
        permutations=int(stored["permutations"]),
        cells=[float(effect) for effect in stored["cells"]],
    )
    _require_equal(rebuilt, stored, Path(parent_dir) / "settings.json")
    span_split = ssc.plasmode_settings(
        export, fidelity, parent_dir=plasmode_dir, replicates=replicates, permutations=permutations, cells=cells
    )
    return {
        "driver_version": DRIVER_VERSION,
        "profile": "plasmode",
        "parent": _parent_block(parent_dir, "plasmode", committed),
        "span_split": span_split,
        "cells": [float(effect) for effect in cells],
        "replicates": int(replicates),
        "permutations": int(permutations),
        "interaction_scan": True,
        "alpha": ALPHA,
        "rejection": {"scan": SCAN_RULE, "parent": ssc.REJECTION_RULES["plasmode"]},
        "scan": {
            **_scan_description(),
            "features": "the analysis's pooled-standardized joint features",
            "models": "the covariate-adjusted full and reduced trajectory models",
            "seed": "the unit's RRPP seed",
            "residual_correction": "df",
        },
    }


def plasmode_units(settings: Mapping[str, Any]) -> list[tuple[float, int]]:
    return ssc.plasmode_units(settings["span_split"])


def plasmode_params(settings: Mapping[str, Any], seeds: Mapping[str, int]) -> Any:
    """The parent's analysis parameters (split test on, attribution off) plus the scan."""

    return dataclasses.replace(ssc.plasmode_params(settings["span_split"], seeds), interaction_scan=True)


def planted_counts(table: Any, plant: Mapping[str, Mapping[str, Any]], *, top_k: int = TOP_K) -> dict[str, Any]:
    """Hits inside and outside the planted set and the planted features among the ``top_k`` largest F.

    Counts only: no feature name leaves this function.
    """

    joint = (table["block"].astype(str) + "__" + table["feature"].astype(str)).to_numpy()
    planted_names = {f"{name}__{feature}" for name in plant for feature in plant[name]["features"]}
    in_planted = np.isin(joint, list(planted_names))
    if int(in_planted.sum()) != len(planted_names):
        raise pc.ConsistencyError("planted: the planted features are not all in the scan table.")
    degenerate = table["degenerate"].to_numpy(dtype=bool)
    f_values = table["f"].to_numpy(dtype=float)
    adjusted = table["p_value_adjusted"].to_numpy(dtype=float)
    hit = ~degenerate & (adjusted <= ALPHA)
    valid = np.flatnonzero(~degenerate)
    top = valid[np.lexsort((valid, -f_values[valid]))][:top_k]
    in_top = np.zeros(len(joint), dtype=bool)
    in_top[top] = True
    blocks = table["block"].astype(str).to_numpy()

    def counts(member: np.ndarray) -> dict[str, int]:
        return {
            "hits": int(np.count_nonzero(hit & member)),
            "hits_in_planted": int(np.count_nonzero(hit & in_planted & member)),
            "hits_outside_planted": int(np.count_nonzero(hit & ~in_planted & member)),
            "top_k_in_planted": int(np.count_nonzero(in_top & in_planted & member)),
            "planted": int(np.count_nonzero(in_planted & member)),
        }

    return {
        "top_k": int(top_k),
        **counts(np.ones(len(joint), dtype=bool)),
        "blocks": {name: counts(blocks == name) for name in dict.fromkeys(blocks.tolist())},
    }


def _plasmode_scan_record(block: Mapping[str, Any], table: Any) -> dict[str, Any]:
    """The analysis's ``interaction_scan`` block without feature names, checked against the scan table."""

    hits = block.get("hits") or []
    if "p_value_adjusted" in table.columns:
        adjusted = table["p_value_adjusted"].to_numpy(dtype=float)
        table_hits = int(np.count_nonzero(~table["degenerate"].to_numpy(dtype=bool) & (adjusted <= ALPHA)))
        if table_hits != len(hits):
            raise pc.ConsistencyError(f"scan: the table has {table_hits} hits, the record block {len(hits)}.")
    return {
        "n_features": block["n_features"],
        "n_degenerate": block["n_degenerate"],
        "df_interaction": block["df_interaction"],
        "df_residual": block["df_residual"],
        "max_f": block["max_f"],
        "max_block": block["max_block"],
        "p_value_global": block.get("p_value_global"),
        "hits": len(hits),
        "blocks": block["blocks"],
        "permutations": block["permutations"],
        "seed": block["seed"],
        "residual_correction": block["residual_correction"],
        "draws_match_rrpp": block.get("draws_match_rrpp"),
    }


def evaluate_plasmode_unit(effect: float, replicate: int, settings: Mapping[str, Any], export: Any) -> dict[str, Any]:
    """One case-study analysis with the split test and the scan; seeds, statistics and counts only."""

    start = perf_counter()
    parent = settings["span_split"]["plasmode"]
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
        _, result = pc.analyse(built["blocks"], built["metadata"], export, params)
        pc.check_consistency(result)
        if float(effect) == pc.ANCHOR:
            pc.check_anchor_latent(result, settings["span_split"]["fidelity"])
        scan = _plasmode_scan_record(result.record["interaction_scan"], result.interaction_scan)
        planted = planted_counts(result.interaction_scan, built["plant"])
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
            "scan": scan,
            "planted": planted,
            "runtime_seconds": perf_counter() - start,
        }
    )


# ── Shards, resume and the run ───────────────────────────────────────────────


_WORKER_EXPORT: Any = None


def _init_plasmode_worker(directory: str, files: Mapping[str, str]) -> None:
    global _WORKER_EXPORT
    export = pc.load_export(Path(directory))
    if export.files != dict(files):
        raise SystemExit(f"{directory} changed after the run started (file hashes differ).")
    _WORKER_EXPORT = export


def _evaluate_plasmode_pending(unit: tuple[float, int], settings: Mapping[str, Any]) -> dict[str, Any]:
    return evaluate_plasmode_unit(*unit, settings, _WORKER_EXPORT)


def _keys(settings: Mapping[str, Any]) -> tuple[list[tuple[str, Any]], Callable[[Mapping[str, Any]], str]]:
    if settings["profile"] == "simulation":
        return [(ssc._unit_key(unit), unit) for unit in simulation_units(settings)], ssc._simulation_key
    return [(pc.unit_key(*unit), unit) for unit in plasmode_units(settings)], ssc._plasmode_key


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
    check_parent_files(settings)
    directory = profile_dir(out_dir, profile)
    lsc.ensure_settings(directory, settings)
    path = lsc.shard_path(directory, index, count)
    path.parent.mkdir(parents=True, exist_ok=True)
    if profile == "plasmode" and (
        export is None or export.files != settings["span_split"]["plasmode"]["export"]["files"]
    ):
        raise SystemExit("The export's file hashes differ from the settings'; refusing to append.")
    keyed, record_key = _keys(settings)
    done = ssc.read_records([path], record_key) if path.exists() else {}
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


SIMULATION_CELL = ("design", "cell", "trajectory_mode", "effect_size", "span_alignment", "span_support")
PLASMODE_CELL = ("cell", "effect")
RATE = ("rate", "se", "count", "rejections", "failures")
SCAN_COLUMNS = {
    "simulation": (*SIMULATION_CELL, *RATE, "hits_median", "max_f_median"),
    "plasmode": (*PLASMODE_CELL, *RATE, "hits_median", "max_f_median"),
}
BLOCK_COLUMNS = {
    "simulation": ("design", "cell", "block", "n", "feature_share", "null_max_share_mean", "null_max_share_se",
                   "hits_mean"),
    "plasmode": ("cell", "effect", "block", "n", "feature_share", "null_max_share_mean", "null_max_share_se",
                 "hits_mean"),
}
COMPARATOR_COLUMNS = {
    "simulation": ("design", "cell", "test", "rule", *RATE),
    "plasmode": ("cell", "effect", "test", "rule", *RATE),
}
FALSE_HIT_COLUMNS = ("design", "cell", "n_population_null_median", *RATE)
RECOVERY_COLUMNS = (
    "cell", "effect", "false_hit_rate", "false_hit_se", "false_hit_count", "false_hit_rejections", "failures",
    "recall_mean", "recall_se", "hit_precision_mean", "hit_precision_se", "hit_precision_n", "top_k_precision_mean",
    "top_k_precision_se", "planted_size",
)


def _scan_p(record: Mapping[str, Any]) -> Any:
    return (record.get("scan") or {}).get("p_value_global")


def _test_p(record: Mapping[str, Any], test: str) -> Any:
    """The p-value of one comparator test in a completed record."""

    if test == "scan":
        return _scan_p(record)
    if test.startswith("split_"):
        part = test.split("_", 1)[1]
        block = record.get("split") if record.get("profile") == "simulation" else record.get("span_split")
        return ((block or {}).get("parts", {}).get(part) or {}).get("p_value_f")
    values = (record.get("production") or {}).get("p_values") if record.get("profile") == "simulation" else (
        (record.get("analysis") or {}).get("p_values")
    )
    return (values or {}).get(test)


def _test_rule(settings: Mapping[str, Any], test: str) -> str:
    return SCAN_RULE if test == "scan" else str(settings["rejection"]["parent"])


def _cell_groups(records: Sequence[Mapping[str, Any]], settings: Mapping[str, Any]) -> list[tuple[dict, list]]:
    """(cell fields, records) per design/cell or per plasmode cell, in a stable order."""

    if settings["profile"] == "simulation":
        out = []
        groups = _grouped(records, lambda r: (str(r["design"]), str(r["cell"])))
        for (design, cell), group in sorted(groups.items()):
            first = group[0]
            out.append(({"design": design, "cell": cell, **{key: first.get(key) for key in SIMULATION_CELL[2:]}},
                        group))
        return out
    groups = _grouped(records, lambda r: (float(r["effect"]),))
    return [
        ({"cell": pc.cell_label(float(effect)), "effect": float(effect)}, groups.get((float(effect),), []))
        for effect in settings["cells"]
    ]


def scan_rows(records: Sequence[Mapping[str, Any]], settings: Mapping[str, Any]) -> list[dict[str, Any]]:
    rows = []
    for fields, group in _cell_groups(records, settings):
        completed = [r for r in group if r.get("status") == "completed"]
        rows.append(
            {
                **fields,
                **_rate_row([_rejects(_scan_p(r)) for r in completed]),
                "failures": len(group) - len(completed),
                "hits_median": _median(r["scan"].get("hits") for r in completed),
                "max_f_median": _median(r["scan"].get("max_f") for r in completed),
            }
        )
    return rows


def block_rows(records: Sequence[Mapping[str, Any]], settings: Mapping[str, Any]) -> list[dict[str, Any]]:
    rows = []
    for fields, group in _cell_groups(records, settings):
        completed = [r for r in group if r.get("status") == "completed"]
        names: list[str] = []
        for record in completed:
            names.extend(name for name in record["scan"]["blocks"] if name not in names)
        for name in names:
            entries = [r["scan"]["blocks"].get(name) or {} for r in completed]
            shares = [r["scan"]["blocks"][name]["n_features"] / r["scan"]["n_features"] for r in completed
                      if name in r["scan"]["blocks"]]
            n, mean, se = _mean_se(entry.get("null_max_share") for entry in entries)
            rows.append(
                {
                    **{key: fields[key] for key in BLOCK_COLUMNS[settings["profile"]][:-6]},
                    "block": name,
                    "n": n,
                    "feature_share": _median(shares),
                    "null_max_share_mean": mean,
                    "null_max_share_se": se,
                    "hits_mean": _mean_se(entry.get("hits") for entry in entries)[1],
                }
            )
    return rows


def comparator_rows(records: Sequence[Mapping[str, Any]], settings: Mapping[str, Any]) -> list[dict[str, Any]]:
    """The scan and the six existing tests on the same units, each under its rule."""

    rows = []
    for fields, group in _cell_groups(records, settings):
        completed = [r for r in group if r.get("status") == "completed"]
        for test in COMPARATOR_TESTS:
            rule = _test_rule(settings, test)
            rows.append(
                {
                    **{key: fields[key] for key in COMPARATOR_COLUMNS[settings["profile"]][:-7]},
                    "test": test,
                    "rule": rule,
                    **_rate_row([_rejects(_test_p(r, test), rule) for r in completed]),
                    "failures": len(group) - len(completed),
                }
            )
    return rows


def false_hit_rows(records: Sequence[Mapping[str, Any]], settings: Mapping[str, Any]) -> list[dict[str, Any]]:
    """Simulation: units with a hit on a population-null feature, among units that have such features."""

    rows = []
    for fields, group in _cell_groups(records, settings):
        completed = [r for r in group if r.get("status") == "completed"]
        with_null = [r for r in completed if (r["scan"].get("n_population_null") or 0) > 0]
        rows.append(
            {
                "design": fields["design"],
                "cell": fields["cell"],
                "n_population_null_median": _median(r["scan"].get("n_population_null") for r in completed),
                **_rate_row([int(r["scan"].get("hits_population_null") or 0) > 0 for r in with_null]),
                "failures": len(group) - len(completed),
            }
        )
    return rows


def recovery_rows(records: Sequence[Mapping[str, Any]], settings: Mapping[str, Any]) -> list[dict[str, Any]]:
    """Plasmode planted cells: false hits outside the planted set, recall, hit and top-k precision."""

    rows = []
    for fields, group in _cell_groups(records, settings):
        if float(fields["effect"]) == pc.ANCHOR:
            continue
        completed = [r for r in group if r.get("status") == "completed"]
        planted = [r["planted"] for r in completed]
        false_hits = _rate_row([int(p["hits_outside_planted"]) > 0 for p in planted])
        size = _median(p["planted"] for p in planted)
        _, recall, recall_se = _mean_se(p["hits_in_planted"] / p["planted"] for p in planted if p["planted"])
        precision_n, precision, precision_se = _mean_se(
            p["hits_in_planted"] / p["hits"] for p in planted if p["hits"]
        )
        _, top, top_se = _mean_se(p["top_k_in_planted"] / min(p["top_k"], p["planted"]) for p in planted
                                  if p["planted"])
        rows.append(
            {
                **fields,
                "false_hit_rate": false_hits["rate"],
                "false_hit_se": false_hits["se"],
                "false_hit_count": false_hits["count"],
                "false_hit_rejections": false_hits["rejections"],
                "failures": len(group) - len(completed),
                "recall_mean": recall,
                "recall_se": recall_se,
                "hit_precision_mean": precision,
                "hit_precision_se": precision_se,
                "hit_precision_n": precision_n,
                "top_k_precision_mean": top,
                "top_k_precision_se": top_se,
                "planted_size": size,
            }
        )
    return rows


def reproduction(
    comparators: Sequence[Mapping[str, Any]],
    settings: Mapping[str, Any],
    committed_split: Sequence[Mapping[str, Any]],
    committed_trajectory: Sequence[Mapping[str, Any]],
) -> dict[str, Any]:
    """The split test's pseudo-F rows and the trajectory rows against the parent's committed tables, exactly."""

    simulation = settings["profile"] == "simulation"
    if simulation:
        in_run = set(settings["span_split"]["comparison"]["designs"])
    else:
        in_run = {pc.cell_label(float(effect)) for effect in settings["cells"]}

    def keep(row: Mapping[str, Any]) -> bool:
        return (row["design"] if simulation else row["cell"]) in in_run

    def place(row: Mapping[str, Any]) -> tuple:
        return (row["design"], row["cell"]) if simulation else (row["cell"],)

    mine = {(*place(row), row["test"]): row for row in comparators}
    expected: list[tuple[tuple, Mapping[str, Any]]] = []
    for row in committed_split:
        if row["statistic"] == "f" and keep(row):
            expected.append(((*place(row), f"split_{row['part']}"), row))
    for row in committed_trajectory:
        if keep(row):
            expected.append(((*place(row), row["statistic"]), row))
    differing = []
    for key, row in expected:
        reproduced = mine.get(key)
        wanted = (int(row["count"]), int(row["rejections"]))
        got = None if reproduced is None else (int(reproduced["count"]), int(reproduced["rejections"]))
        if got != wanted:
            differing.append(
                {
                    "cell": list(key[:-1]),
                    "test": key[-1],
                    "committed": {"count": wanted[0], "rejections": wanted[1]},
                    "reproduced": None if got is None else {"count": got[0], "rejections": got[1]},
                }
            )
    names = PARENT_FILES[settings["profile"]]
    return {
        "profile": settings["profile"],
        "passed": bool(expected) and not differing,
        "rows_compared": len(expected),
        "differing": differing,
        "rule": settings["rejection"]["parent"],
        "sources": {name: settings["parent"]["files"][name] for name in names[1:]},
        "source_dir": settings["parent"]["dir"],
    }


def accounting(
    records: Sequence[Mapping[str, Any]], settings: Mapping[str, Any], *, expected: int, duplicates: int
) -> dict[str, Any]:
    profile = settings["profile"]
    failures: dict[str, int] = {}
    unit_seconds: list[float] = []
    stages: dict[str, list[float]] = {}
    for record in records:
        label = f"{record['design']}/{record['cell']}" if profile == "simulation" else str(record["cell"])
        if record.get("status") != "completed":
            failures[label] = failures.get(label, 0) + 1
        runtime = record.get("runtime_seconds")
        if isinstance(runtime, Mapping):
            unit_seconds.append(float(runtime.get("unit") or 0.0))
            for stage in ("fit", "evaluation", "split", "scan"):
                if runtime.get(stage) is not None:
                    stages.setdefault(stage, []).append(float(runtime[stage]))
        elif runtime is not None:
            unit_seconds.append(float(runtime))
    scan_share = None
    if stages.get("scan") and unit_seconds:
        scan_share = float(sum(stages["scan"]) / sum(unit_seconds)) if sum(unit_seconds) > 0 else None
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
        "scan_share_of_unit_seconds": scan_share,
    }


def merge(out_dir: Path, profile: str) -> dict[str, Any]:
    """Merge one profile's shards and write its tables, accounting and reproduction check."""

    directory = profile_dir(out_dir, profile)
    settings = lsc.load_settings(directory)
    check_parent_files(settings)
    keyed, record_key = _keys(settings)
    expected_keys = [key for key, _ in keyed]
    raw: list[dict[str, Any]] = []
    for path in sorted((directory / "records").glob("shard_*.jsonl")):
        raw.extend(ssc.read_records([path], record_key).values())
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
    comparators = comparator_rows(records, settings)
    _write_csv(directory / "scan_table.csv", scan_rows(records, settings), SCAN_COLUMNS[profile])
    _write_csv(directory / "block_balance.csv", block_rows(records, settings), BLOCK_COLUMNS[profile])
    _write_csv(directory / "comparator_table.csv", comparators, COMPARATOR_COLUMNS[profile])
    if profile == "simulation":
        _write_csv(directory / "false_hit_table.csv", false_hit_rows(records, settings), FALSE_HIT_COLUMNS)
    else:
        _write_csv(directory / "recovery_table.csv", recovery_rows(records, settings), RECOVERY_COLUMNS)
    parent = ssc._resolve(settings["parent"]["dir"])
    names = PARENT_FILES[profile]
    result = reproduction(comparators, settings, lsc._read_csv(parent / names[1]), lsc._read_csv(parent / names[2]))
    _write_json(directory / "accounting.json", summary)
    _write_json(directory / "reproduction.json", result)
    return {
        "profile": profile,
        **{key: summary[key] for key in ("complete_units", "missing_units", "failed_units", "recorded_core_hours")},
        "reproduction_passed": result["passed"],
        "reproduction_differing": len(result["differing"]),
    }


# ── Decision (pure functions) ────────────────────────────────────────────────


def _rate(row: Mapping[str, Any] | None, prefix: str = "") -> dict[str, Any] | None:
    if row is None or _number(row.get(f"{prefix}rate")) is None:
        return None
    return {"rate": float(row[f"{prefix}rate"]), "se": float(row[f"{prefix}se"]), "count": int(row[f"{prefix}count"])}


def _size_observation(row: Mapping[str, Any] | None, *, prefix: str = "", **cell: Any) -> dict[str, Any]:
    if row is None or int(row.get(f"{prefix}count") or 0) <= 0:
        return {**cell, "status": "unavailable"}
    return {**cell, "n_rejected": int(row[f"{prefix}rejections"]), "n_available": int(row[f"{prefix}count"])}


def _failing_observations(judged: Mapping[str, Any]) -> list[dict[str, Any]]:
    return ssc._failing_observations(judged)


def _source(settings: Mapping[str, Any], count: Any) -> dict[str, str]:
    shown = "?" if count in (None, "") else int(count)
    return {"source": f"this run, {settings['profile']} profile", "grade": f"{shown} × {settings['permutations']}"}


def decision(
    simulation: Mapping[str, Any],
    plasmode: Mapping[str, Any],
    *,
    attribution_precision: Sequence[Mapping[str, Any]] = (),
) -> dict[str, Any]:
    """Apply rule version 1 to both profiles' merged tables.

    ``simulation`` and ``plasmode`` hold ``settings``, ``accounting``,
    ``reproduction``, ``scan``, ``blocks`` and ``comparators`` (rows); the
    simulation also ``false_hits``, the plasmode ``recovery``.
    ``attribution_precision`` holds the plasmode's committed planted-set
    precision rows, read only by a descriptive table.
    """

    profiles = (("simulation", simulation), ("plasmode", plasmode))
    header = {
        "rule_version": RULE_VERSION,
        "alpha": ALPHA,
        "se_multiplier": SE_MULTIPLIER,
        "scan_rejection": SCAN_RULE,
        "settings_hashes": {name: lsc.settings_hash(block["settings"]) for name, block in profiles},
        "complete": all(block["accounting"]["missing_units"] == 0 for _, block in profiles),
        "failed_units": {name: block["accounting"]["failed_units"] for name, block in profiles},
        "reproduction": {
            name: {key: block["reproduction"].get(key) for key in ("passed", "rows_compared", "rule", "source_dir")}
            for name, block in profiles
        },
    }
    failed = [
        {"profile": name, "differing": block["reproduction"].get("differing", [])}
        for name, block in profiles
        if not block["reproduction"].get("passed")
    ]
    if failed:
        return {**header, "verdict": "reproduction_failed", "consequence": CONSEQUENCES["reproduction_failed"],
                "failing": {"reproduction": failed}}

    sim_scan = {(r["design"], r["cell"]): r for r in simulation["scan"]}
    pm_scan = {r["cell"]: r for r in plasmode["scan"]}
    recovery = {r["cell"]: r for r in plasmode["recovery"]}
    comparators = {(r["cell"], r["test"]): r for r in plasmode["comparators"]}
    designs = list(simulation["settings"]["span_split"]["comparison"]["designs"])

    s1 = grade.judge_rates(
        [
            _size_observation(sim_scan.get((design, cell)), profile="simulation", design=design, cell=cell)
            for design in designs
            for cell in CONTROL_CELLS
        ]
    )
    s2 = grade.judge_rates([_size_observation(pm_scan.get("anchor"), profile="plasmode", cell="anchor")])
    s3 = grade.judge_rates(
        [
            _size_observation(recovery.get(pc.cell_label(effect)), prefix="false_hit_", profile="plasmode",
                              cell=pc.cell_label(effect), quantity="units with a hit outside the planted set")
            for effect in PLANTED_EFFECTS
        ]
    )
    label = pc.cell_label(GAIN_EFFECT)
    scan_rate = _rate(comparators.get((label, "scan")))
    existing = {test: _rate(comparators.get((label, test))) for test in EXISTING_TESTS}
    available = {test: rate for test, rate in existing.items() if rate is not None}
    g: dict[str, Any] = {"cell": label, "scan": scan_rate, "existing": existing}
    if scan_rate is None or len(available) < len(EXISTING_TESTS):
        g.update({"pass": False, "missing": [test for test in EXISTING_TESTS if existing[test] is None]
                  + ([] if scan_rate is not None else ["scan"])})
    else:
        comparator = max(EXISTING_TESTS, key=lambda test: (available[test]["rate"], -EXISTING_TESTS.index(test)))
        best = available[comparator]
        pooled = math.sqrt(scan_rate["se"] ** 2 + best["se"] ** 2)
        difference = scan_rate["rate"] - best["rate"]
        g.update({"comparator": comparator, "comparator_rate": best, "difference": difference, "pooled_se": pooled,
                  "threshold": SE_MULTIPLIER * pooled, "pass": difference > SE_MULTIPLIER * pooled})
    criteria = {
        "S1": {"name": "simulation size", **s1},
        "S2": {"name": "real-noise size", **s2},
        "S3": {"name": "strong familywise control", **s3},
        "G": {"name": "sparse gain", **g},
    }
    failing: dict[str, Any] = {}
    for name, judged in (("S1", s1), ("S2", s2), ("S3", s3)):
        if not judged["holds"]:
            failing[name] = _failing_observations(judged)
    if not g["pass"]:
        failing["G"] = {key: g.get(key) for key in ("cell", "comparator", "difference", "pooled_se", "missing")
                        if g.get(key) is not None}
    if "S1" in failing or "S2" in failing:
        verdict = "size_excess"
    elif "S3" in failing:
        verdict = "fwer_excess"
    elif "G" in failing:
        verdict = "no_sparse_gain"
    else:
        verdict = "scan_validated"
    return {
        **header,
        "verdict": verdict,
        "consequence": CONSEQUENCES[verdict],
        "failing": failing,
        "criteria": criteria,
        "descriptive": descriptive_tables(simulation, plasmode, attribution_precision=attribution_precision),
    }


def descriptive_tables(
    simulation: Mapping[str, Any],
    plasmode: Mapping[str, Any],
    *,
    attribution_precision: Sequence[Mapping[str, Any]] = (),
) -> dict[str, Any]:
    """Plasmode power and recovery, block balance, simulation localization and false hits; each row sourced."""

    sim_settings, pm_settings = simulation["settings"], plasmode["settings"]
    pm_comparators = {(r["cell"], r["test"]): r for r in plasmode["comparators"]}
    sim_comparators = {(r["design"], r["cell"], r["test"]): r for r in simulation["comparators"]}
    recovery = {r["cell"]: r for r in plasmode["recovery"]}
    attribution: dict[str, list[dict[str, Any]]] = {}
    for row in attribution_precision:
        attribution.setdefault(str(row["cell"]), []).append(
            {
                "arm": row["arm"],
                "component": row["component"],
                "planted_precision": _number(row["mean"]),
                "se": _number(row.get("se")),
                "n": None if row.get("n") in (None, "") else int(row["n"]),
                "source": ATTRIBUTION_SOURCE,
                "grade": "500 × 999, B = 200",
            }
        )
    power = []
    for effect in PLANTED_EFFECTS:
        label = pc.cell_label(effect)
        scan_row = pm_comparators.get((label, "scan"))
        rec = recovery.get(label) or {}
        power.append(
            {
                "cell": label,
                "scan": _rate(scan_row),
                "recall": _number(rec.get("recall_mean")),
                "hit_precision": _number(rec.get("hit_precision_mean")),
                "top_k_precision": _number(rec.get("top_k_precision_mean")),
                "false_hit_rate": _number(rec.get("false_hit_rate")),
                **{test: _rate(pm_comparators.get((label, test))) for test in EXISTING_TESTS},
                "attribution_planted_precision": attribution.get(label, []),
                **_source(pm_settings, None if scan_row is None else scan_row.get("count")),
            }
        )
    balance = []
    for name, block in (("simulation", simulation), ("plasmode", plasmode)):
        for row in block["blocks"]:
            balance.append(
                {
                    "profile": name,
                    **({"design": row["design"]} if name == "simulation" else {}),
                    "cell": row["cell"],
                    "block": row["block"],
                    "feature_share": _number(row.get("feature_share")),
                    "null_max_share": _number(row.get("null_max_share_mean")),
                    "null_max_share_se": _number(row.get("null_max_share_se")),
                    **_source(block["settings"], row.get("n")),
                }
            )
    localization = []
    for row in simulation["scan"]:
        design, cell = row["design"], row["cell"]
        localization.append(
            {
                "design": design,
                "cell": cell,
                **{test: _rate(sim_comparators.get((design, cell, test))) for test in COMPARATOR_TESTS},
                "rules": {"scan": SCAN_RULE, "existing": sim_settings["rejection"]["parent"]},
                **_source(sim_settings, row.get("count")),
            }
        )
    false_hits = [
        {
            "design": row["design"],
            "cell": row["cell"],
            "n_population_null_median": _number(row.get("n_population_null_median")),
            "rate": _number(row.get("rate")),
            "se": _number(row.get("se")),
            "count": None if row.get("count") in (None, "") else int(row["count"]),
            **_source(sim_settings, row.get("count")),
        }
        for row in simulation["false_hits"]
    ]
    return {
        "plasmode_power_and_recovery": power,
        "block_balance": balance,
        "simulation_localization": localization,
        "simulation_false_hits": false_hits,
    }


def _profile_tables(directory: Path, profile: str) -> dict[str, Any]:
    settings = lsc.load_settings(directory)
    check_parent_files(settings)
    path = directory / "accounting.json"
    if not path.exists():
        raise SystemExit(f"{path} does not exist; run `merge --profile {profile}` first.")
    summary = json.loads(path.read_text(encoding="utf-8"))
    if summary.get("settings_hash") != lsc.settings_hash(settings):
        raise SystemExit(f"{path} was merged under other settings; run `merge` again.")
    block: dict[str, Any] = {
        "settings": settings,
        "accounting": summary,
        "reproduction": json.loads((directory / "reproduction.json").read_text(encoding="utf-8")),
        "scan": _read_csv(directory / "scan_table.csv"),
        "blocks": _read_csv(directory / "block_balance.csv"),
        "comparators": _read_csv(directory / "comparator_table.csv"),
    }
    if profile == "simulation":
        block["false_hits"] = _read_csv(directory / "false_hit_table.csv")
    else:
        block["recovery"] = _read_csv(directory / "recovery_table.csv")
    return block


def attribution_rows(settings: Mapping[str, Any]) -> list[dict[str, Any]]:
    """The plasmode's committed planted-set precision of attribution (adjusted arm, principal block, all units)."""

    parent = settings["span_split"]["parent"]
    path = ssc._resolve(parent["dir"]) / "plasmode_summary.csv"
    if ssc.sha256_file(path) != parent["summary_sha256"]:
        raise SystemExit(f"{path} has changed since the settings were written (sha256 differs); refusing.")
    labels = {pc.cell_label(effect) for effect in PLANTED_EFFECTS}
    return [
        row
        for row in lsc._read_csv(path)
        if row["quantity"] == "planted_precision" and row["arm"] == pc.PRIMARY_ARM and row["block"] == pc.PRINCIPAL
        and row["subset"] == "all" and row["cell"] in labels
    ]


def decide(out_dir: Path) -> dict[str, Any]:
    """Read both profiles' merged tables and write ``interaction_scan_decision.json``."""

    simulation = _profile_tables(profile_dir(out_dir, "simulation"), "simulation")
    plasmode = _profile_tables(profile_dir(out_dir, "plasmode"), "plasmode")
    result = decision(simulation, plasmode, attribution_precision=attribution_rows(plasmode["settings"]))
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
    run.add_argument("--perms", type=int, default=ssc.DEFAULT_PERMUTATIONS)
    run.add_argument("--reps", type=int, default=None,
                     help=f"Replicates per cell (simulation: non-control cells, default {ssc.SIMULATION_REPLICATES}; "
                          f"plasmode: default {ssc.PLASMODE_REPLICATES})")
    run.add_argument("--control-reps", type=int, default=ssc.CONTROL_REPLICATES,
                     help="Simulation: replicates of the anchor and translation")
    run.add_argument("--designs", type=str, default=",".join(lsc.DESIGNS), help="Simulation: comma-separated")
    run.add_argument("--parent-dir", type=Path, default=SPAN_SPLIT_DIR,
                     help="The committed split-test validation run (default: results/span-split-2026-11-05)")
    run.add_argument("--export-dir", type=Path, default=None, help="Plasmode: the SEA-AD cohort_30 export")
    run.add_argument("--fidelity", type=Path, default=None, help="Plasmode: fidelity.json of a passed gate")
    merge_parser = sub.add_parser("merge", help="Merge one profile's shards")
    merge_parser.add_argument("--profile", choices=PROFILES, required=True)
    merge_parser.add_argument("--out-dir", type=Path, required=True)
    decide_parser = sub.add_parser("decide", help=f"Write {DECISION_NAME} from both profiles")
    decide_parser.add_argument("--out-dir", type=Path, required=True)
    return parser


def _settings_from_args(args: argparse.Namespace) -> tuple[dict[str, Any], Any]:
    if args.profile == "simulation":
        settings = simulation_settings(
            parent_dir=args.parent_dir / "simulation",
            replicates=ssc.SIMULATION_REPLICATES if args.reps is None else args.reps,
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
        parent_dir=args.parent_dir / "plasmode",
        replicates=ssc.PLASMODE_REPLICATES if args.reps is None else args.reps,
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
    failing = {name: (len(entries) if isinstance(entries, list) else 1) for name, entries in result["failing"].items()}
    print(json.dumps({"verdict": result["verdict"], "failing": failing, "complete": result["complete"]}, indent=1))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
