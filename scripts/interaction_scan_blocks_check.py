#!/usr/bin/env python3
"""Per-feature interaction scan by block: the block-wise combination against the joint scan on paired units.

The driver behind ``openspec/changes/interaction-scan-block-maxima/`` (spec
``interaction-scan-block-validation``, roadmap item 16, thread 3). The joint
scan (:func:`motco.stats.interaction.interaction_feature_scan`, MOTCO 0.20.0,
``scan_validated``) takes one maximum F over every feature, so its familywise
bar is set by whichever block most often gives the largest null F: on the
SEA-AD export ATAC gave 0.618 of the null maxima with half the features. The
by-block combination (``blocks=``, ``block_bonferroni``, MOTCO 0.23.0) runs the
step-down within each block on the same draws and tests each of the K blocks
at α/K. This driver measures it against the joint scan before a case study may
use it.

**Arms.** Every unit runs the scan twice on the same features, models, seed
and residual correction: ``joint`` (no labels, the 0.20.0 scan) and
``by_block`` (each feature's block as its label). Both arms share every F;
they differ only in their p-values and hits. A scan rejects globally, and a
feature is a hit, at adjusted p ≤ 0.05. No unit fits a latent space: the scan
never reads one.

**Profiles** (every output lands under ``--out-dir``/``<profile>/``):

- ``simulation`` — the interaction-scan validation's simulation units (the
  latent-space comparison's three designs, fourteen cells, base seed 2035, 498
  features; anchor and translation × 500, other cells × 100, 999
  permutations). The settings are rebuilt through
  ``interaction_scan_check.simulation_settings`` and must equal the committed
  ``results/interaction-scan-2026-11-09/simulation/settings.json`` (differing
  keys named). A unit generates the comparison's dataset, takes the harness's
  pooled-standardized features on the measured blocks (a per-feature affine
  image of the parent's PLS input coordinates, which every F ignores) and the
  trajectory models from the design, and runs both arms at the unit's RRPP
  seed without a correction, as the parent did. Population-null features
  (population B − A constant across stages) are the parent's mask.
- ``plasmode`` — the interaction-scan validation's SEA-AD plasmode units (the
  anchor, e = 2.0 and 3.0 × 500 at 999 permutations, base seed 2034, ``paper``
  stream). The settings come through ``interaction_scan_check.plasmode_settings``
  with the parent's committed fidelity reading and must equal its committed
  settings; the export's file hashes must equal the plasmode's. A unit is
  ``plasmode_check.build_unit`` and the analysis's own scan inputs
  (``motco.analysis._interaction_scan_inputs``) at the parent's parameters,
  both arms with ``df``. The fidelity gate and the anchor latent check are not
  repeated: they check the latent space, which the scan does not use.
- ``plasmode_unequal`` — the same export and units with the RNA block reduced
  to 500 of its 5,000 features: the first 500 of a permutation of its columns
  (seed by sha256 of the driver's name, base seed 2034 and ``unequal_subset``),
  kept in export order and recorded by the sha256 of the feature list. The 10
  planted RNA features of a unit are drawn from the 500 by the plasmode's own
  unit builder. It has no parent run to reproduce.

A record holds seeds, statistics and counts only: per arm the global p-value,
the hits per block and the planted hits and hits outside the planted set per
block (plasmode profiles) or the hits on population-null features per block
(simulation); per unit the joint arm's per-block null-maximum shares, the 20
largest F's planted count, the Spearman correlation of the two block maxima
across draws, and the global p-value a minP combination between blocks would
give. Never a feature name, donor identifier, label vector or donor-level
value. A unit that raises is recorded as failed; the shard continues.

Subcommands:

- ``run --profile …`` — evaluate (a shard of) the units into
  ``records/shard_<i>_of_<N>.jsonl``; the settings and their hash go to
  ``settings.json``, a resume or merge under other settings exits naming the
  differing keys, and a resume skips complete units. Every parent file's
  sha256 is recorded and a changed file is refused by name in ``run``,
  ``merge`` and ``decide``.
- ``merge --profile …`` — ``merged.jsonl`` (per unit; keep off the
  repository), ``scan_table.csv`` (per cell and arm: the global rejection
  rate, Monte Carlo SE, count and median hits), ``block_table.csv`` (per cell,
  arm and block: units with a hit in the block, mean hits, the joint arm's
  null-maximum share, and in the plasmode profiles recall, hit precision and
  units with a hit outside the planted set), ``false_hit_table.csv`` (per cell
  and arm: units with a hit on a population-null feature, or outside the
  planted set), ``paired_table.csv`` (per cell: the arms' discordant global
  rejections and exact loss and gain p-values), ``minp_table.csv`` (per cell:
  the joint, Bonferroni and minP global rates and the median correlation of
  block maxima), ``allocation_table.csv`` (plasmode profiles: per planted cell
  and block, the per-unit difference in planted hits, by block minus joint,
  and its exact sign test), ``accounting.json`` and, for ``simulation`` and
  ``plasmode``, ``reproduction.json`` (the joint arm against the parent's
  committed tables, exact at their written precision).
- ``decide`` — ``interaction_scan_blocks_decision.json`` (rule version 1):
  ``reproduction_failed`` gates every criterion; S1 the by-block global rate
  at the anchor and translation of every simulation design (m = 6), S2 at both
  exports' anchors (m = 2), S3 units with a hit outside the planted set at
  e = 2.0 and 3.0 on both exports (m = 4), all under the stage-span paper
  grade's count-plus-ceiling rule; P no exact one-sided paired loss p-value of
  the by-block global rejections against the joint ones below 0.05/4 at
  e = 2.0 and 3.0 on both exports; G the exact one-sided sign test of the
  per-unit gain in planted reduced-block hits on the unequal export at
  p ≤ 0.025 at e = 2.0 or 3.0. Verdict: ``size_excess`` | ``fwer_excess`` |
  ``power_loss`` | ``no_allocation_gain`` | ``by_block_validated``, with its
  consequence (the joint scan stays the default under every verdict) and
  descriptive tables that never change it.

With ``--workers`` above 1, pin BLAS to one thread per process
(``OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1``).

Examples:

    export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
    EXPORT=../sead_sexdiff/results/motco_input/cohort_30

    # smoke
    python scripts/interaction_scan_blocks_check.py run --profile simulation --reps 2 --control-reps 2 \\
        --perms 99 --out-dir /tmp/isb --workers 7
    python scripts/interaction_scan_blocks_check.py run --profile plasmode --export-dir $EXPORT --reps 2 \\
        --perms 99 --out-dir /tmp/isb --workers 6
    python scripts/interaction_scan_blocks_check.py run --profile plasmode_unequal --export-dir $EXPORT \\
        --reps 2 --perms 99 --out-dir /tmp/isb --workers 6
    python scripts/interaction_scan_blocks_check.py merge --profile simulation --out-dir /tmp/isb   # etc.
    python scripts/interaction_scan_blocks_check.py decide --out-dir /tmp/isb

    # one SLURM array task (scripts/interaction_scan_blocks_check_array.sbatch), then merge and decide
"""

from __future__ import annotations

import argparse
import dataclasses
import hashlib
import json
import math
import sys
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path
from time import perf_counter
from typing import Any, Callable, Iterable, Mapping, Sequence

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import interaction_scan_check as isc  # noqa: E402  (the parent run's settings, units, summaries and tables)
import latent_space_comparison as lsc  # noqa: E402  (units, seeds, settings binding, CSV helpers)
import plasmode_check as pc  # noqa: E402  (export contract, units, planted set)
import span_split_check as ssc  # noqa: E402  (parent loading, record reading, rate rows)
import stage_span_paper_grade_decision as grade  # noqa: E402  (the count-plus-ceiling rule)

sys.path.pop(0)

DRIVER = "interaction_scan_blocks_check"
DRIVER_VERSION = 1
RULE_VERSION = 1
PROFILES: tuple[str, ...] = ("simulation", "plasmode", "plasmode_unequal")
PLASMODE_PROFILES: tuple[str, ...] = ("plasmode", "plasmode_unequal")
#: The profiles whose joint arm must reproduce the parent run.
REPRODUCED_PROFILES: tuple[str, ...] = ("simulation", "plasmode")
ARMS: tuple[str, ...] = ("joint", "by_block")
ALPHA = 0.05
SCAN_RULE = isc.SCAN_RULE
TOP_K = isc.TOP_K
CONTROL_CELLS: tuple[str, ...] = isc.CONTROL_CELLS
PLANTED_EFFECTS: tuple[float, ...] = isc.PLANTED_EFFECTS
POWER_M = 4
POWER_LEVEL = ALPHA / POWER_M
GAIN_LEVEL = 0.025

UNEQUAL_BLOCK = "rna"
UNEQUAL_SIZE = 500
UNEQUAL_PURPOSE = "unequal_subset"

REPO = Path(__file__).resolve().parents[1]
PARENT_DIR = REPO / "results" / "interaction-scan-2026-11-09"
#: The parent's committed files each profile binds by sha256 (the unequal profile binds the plasmode's settings).
PARENT_FILES: dict[str, tuple[str, ...]] = {
    "simulation": ("settings.json", "scan_table.csv", "block_balance.csv", "false_hit_table.csv"),
    "plasmode": ("settings.json", "scan_table.csv", "block_balance.csv", "recovery_table.csv"),
    "plasmode_unequal": ("settings.json",),
}
PARENT_PROFILE: dict[str, str] = {"simulation": "simulation", "plasmode": "plasmode", "plasmode_unequal": "plasmode"}
PARENT_SOURCE = "interaction-scan validation 2026-11-09"
#: The parent's tables and the columns its joint arm must reproduce exactly (at their written precision).
REPRODUCED_COLUMNS: dict[str, dict[str, tuple[str, ...]]] = {
    "simulation": {
        "scan_table.csv": ("count", "rejections", "hits_median"),
        "block_balance.csv": ("n", "hits_mean", "null_max_share_mean"),
        "false_hit_table.csv": ("count", "rejections"),
    },
    "plasmode": {
        "scan_table.csv": ("count", "rejections", "hits_median"),
        "block_balance.csv": ("n", "hits_mean", "null_max_share_mean"),
        "recovery_table.csv": (
            "false_hit_count", "false_hit_rejections", "recall_mean", "hit_precision_mean", "top_k_precision_mean",
        ),
    },
}

VERDICTS: tuple[str, ...] = (
    "reproduction_failed", "size_excess", "fwer_excess", "power_loss", "no_allocation_gain", "by_block_validated",
)
_DEFAULT = " Under every verdict the joint scan stays the default; the by-block combination stays opt-in."
CONSEQUENCES: dict[str, str] = {
    "reproduction_failed": (
        "Nothing is read: the joint arm does not reproduce the interaction-scan validation, so no rate pairs with "
        "it. Find the drift before any rerun." + _DEFAULT
    ),
    "size_excess": (
        "The by-block scan is not validated: the size of its global test is not established at the failing cells; "
        "it is documented with the failing criterion and cells." + _DEFAULT
    ),
    "fwer_excess": (
        "The by-block scan is not validated for naming features: features outside the planted set were hits more "
        "often than familywise control allows under the plasmode's partial null; it is documented with the failing "
        "cells." + _DEFAULT
    ),
    "power_loss": (
        "The by-block scan is not validated: its global test loses power against the joint scan at the failing "
        "cells; it is documented with the failing cells and counts." + _DEFAULT
    ),
    "no_allocation_gain": (
        "The by-block scan is not validated: on the unequal export it does not name more planted features of the "
        "small block than the joint scan; it is documented with the sign tests." + _DEFAULT
    ),
    "by_block_validated": (
        "The by-block scan is validated for exports whose blocks differ in size or noise: it controls the "
        "familywise error over all scanned features at 0.05, split equally across the blocks, keeps the joint "
        "scan's global power and names more features of the small block. Only two blocks were measured." + _DEFAULT
    ),
}

DECISION_NAME = "interaction_scan_blocks_decision.json"


# ── Shared helpers ───────────────────────────────────────────────────────────


def _canonical(value: Any) -> Any:
    return json.loads(json.dumps(value, sort_keys=True))


def profile_dir(out_dir: Path, profile: str) -> Path:
    return Path(out_dir) / profile


def _rejects(p_value: Any) -> bool | None:
    return ssc._rejects(p_value, ALPHA, SCAN_RULE)


def _median(values: Iterable[Any]) -> float | None:
    return ssc._median(values)


def _mean_se(values: Iterable[Any]) -> tuple[int, float | None, float | None]:
    return isc._mean_se(values)


def _number(value: Any) -> float | None:
    return ssc._number(value)


def _integer(value: Any) -> int | None:
    return None if value in (None, "") else int(float(value))


def _write_csv(path: Path, rows: Sequence[Mapping[str, Any]], columns: Sequence[str]) -> None:
    lsc._write_csv(path, rows, columns)


def _write_json(path: Path, payload: Any) -> None:
    lsc._write_json(path, payload)


def _read_csv(path: Path) -> list[dict[str, Any]]:
    return ssc._read_csv(path)


def _written(value: Any) -> str:
    """A value as the CSV writer writes it (floats to 6 significant digits, ``None`` empty)."""

    return str(lsc._csv_value(value))


def _upper_binomial(successes: int, trials: int) -> float:
    """P(Binomial(trials, ½) ≥ successes), exactly; 1 when there are no trials."""

    if trials <= 0 or successes <= 0:
        return 1.0
    return float(sum(math.comb(trials, k) for k in range(successes, trials + 1)) / 2.0**trials)


# ── Parent binding ───────────────────────────────────────────────────────────


def _require_equal(rebuilt: Mapping[str, Any], committed: Mapping[str, Any], path: Path) -> None:
    if _canonical(dict(rebuilt)) != _canonical(dict(committed)):
        differing = lsc.settings_difference(_canonical(dict(committed)), _canonical(dict(rebuilt)))
        raise SystemExit(
            f"{path} differs from the settings interaction_scan_check rebuilds at its grade; differing: "
            f"{', '.join(differing)}. The units would not be the interaction-scan validation's."
        )


def _parent_block(parent_dir: Path, profile: str, committed: Mapping[str, Any]) -> dict[str, Any]:
    return {
        "dir": ssc._display_path(parent_dir),
        "source": PARENT_SOURCE,
        "settings_hash": committed["hash"],
        "files": {name: ssc.sha256_file(Path(parent_dir) / name) for name in PARENT_FILES[profile]},
    }


def check_parent_files(settings: Mapping[str, Any]) -> None:
    """Refuse, by name, a parent file (this run's or the parent's own) whose sha256 differs from the recorded one."""

    directory = ssc._resolve(settings["parent"]["dir"])
    for name, digest in settings["parent"]["files"].items():
        path = directory / name
        if not path.exists():
            raise SystemExit(f"{path} does not exist; the parent run's committed file is required.")
        if ssc.sha256_file(path) != digest:
            raise SystemExit(f"{path} has changed since the settings were written (sha256 differs); refusing.")
    isc.check_parent_files(settings["interaction_scan"])


def _scan_description(**extra: Any) -> dict[str, Any]:
    return {
        "function": "motco.stats.interaction.interaction_feature_scan",
        "statistic": "per-feature interaction F (SS_int / df_int) / (RSS / df_res)",
        "arms": {
            "joint": "no block labels: maximum F, Westfall-Young step-down max-T adjusted p-values",
            "by_block": "the feature's block as its label: step-down within each block, Bonferroni over the K "
            "blocks with a non-degenerate feature (adjusted p = min(1, K * within-block p))",
        },
        "hit_alpha": ALPHA,
        "rejection": SCAN_RULE,
        "top_k": TOP_K,
        "minp": "descriptive: each block maximum's tail probability under the block's own null maxima (the "
        "observed configuration counted as one of B + 1), the smallest compared with every draw's smallest",
        **extra,
    }


# ── Settings ─────────────────────────────────────────────────────────────────


def simulation_settings(
    *,
    parent_dir: Path = PARENT_DIR / "simulation",
    span_split_dir: Path = isc.SPAN_SPLIT_DIR / "simulation",
    replicates: int = ssc.SIMULATION_REPLICATES,
    control_replicates: int = ssc.CONTROL_REPLICATES,
    permutations: int = ssc.DEFAULT_PERMUTATIONS,
    designs: Sequence[str] = lsc.DESIGNS,
) -> dict[str, Any]:
    """The interaction-scan validation's simulation settings at this run's grade, its committed ones checked."""

    committed = ssc.load_parent(parent_dir)
    stored = committed["settings"]
    rebuilt = isc.simulation_settings(
        parent_dir=span_split_dir,
        replicates=int(stored["replicates"]["other"]),
        control_replicates=int(stored["replicates"]["control"]),
        permutations=int(stored["permutations"]),
        designs=list(stored["span_split"]["comparison"]["designs"]),
    )
    _require_equal(rebuilt, stored, Path(parent_dir) / "settings.json")
    current = isc.simulation_settings(
        parent_dir=span_split_dir,
        replicates=replicates,
        control_replicates=control_replicates,
        permutations=permutations,
        designs=designs,
    )
    return {
        "driver_version": DRIVER_VERSION,
        "profile": "simulation",
        "parent": _parent_block(parent_dir, "simulation", committed),
        "interaction_scan": current,
        "replicates": dict(current["replicates"]),
        "permutations": int(permutations),
        "alpha": ALPHA,
        "arms": list(ARMS),
        "rejection": SCAN_RULE,
        "scan": _scan_description(
            features="the harness's pooled-standardized features on the measured blocks (no latent fit)",
            models="full group x stage and reduced group + stage trajectory models",
            seed="the unit's RRPP seed",
            residual_correction=None,
        ),
        "population_null_tolerance": current["population_null_tolerance"],
    }


def _features_sha256(features: Sequence[str]) -> str:
    return hashlib.sha256(json.dumps(list(features), separators=(",", ":")).encode("utf-8")).hexdigest()


def unequal_seed(base_seed: int) -> int:
    """The subset's seed, derived as the plasmode derives its seeds, from this driver's name."""

    encoded = json.dumps([DRIVER, int(base_seed), UNEQUAL_PURPOSE], separators=(",", ":")).encode("utf-8")
    return int(hashlib.sha256(encoded).hexdigest()[:8], 16) & 0x7FFFFFFF


def _subset_features(export: Any, block: str, size: int, seed: int) -> list[str]:
    names = [str(name) for name in export.blocks[block].columns]
    if not 0 < size < len(names):
        raise SystemExit(f"The unequal subset needs 0 < size < {len(names)} {block} features, got {size}.")
    chosen = np.sort(np.random.default_rng(seed).permutation(len(names))[:size])
    return [names[index] for index in chosen]


def unequal_subset(export: Any, *, base_seed: int, block: str = UNEQUAL_BLOCK, size: int = UNEQUAL_SIZE) -> dict:
    """The fixed reduction of one block: a seeded subset in export order, recorded by its feature list's sha256."""

    seed = unequal_seed(base_seed)
    features = _subset_features(export, block, int(size), seed)
    return {
        "block": block,
        "size": int(size),
        "of": int(export.blocks[block].shape[1]),
        "seed": seed,
        "derivation": f"sha256 of [{DRIVER!r}, base seed, {UNEQUAL_PURPOSE!r}]",
        "selection": "the first `size` positions of a permutation of the block's columns, kept in export order",
        "features_sha256": _features_sha256(features),
    }


def reduced_export(export: Any, subset: Mapping[str, Any]) -> Any:
    """The export with ``subset``'s block reduced to its features, refused if the feature list's sha256 differs."""

    features = _subset_features(export, str(subset["block"]), int(subset["size"]), int(subset["seed"]))
    if _features_sha256(features) != subset["features_sha256"]:
        raise SystemExit("The export's reduced feature list differs from the settings' (sha256 differs); refusing.")
    blocks = dict(export.blocks)
    blocks[str(subset["block"])] = blocks[str(subset["block"])][features]
    return dataclasses.replace(export, blocks=blocks, _cache={})


def plasmode_settings(
    export: Any,
    *,
    profile: str = "plasmode",
    parent_dir: Path = PARENT_DIR / "plasmode",
    span_split_dir: Path = isc.SPAN_SPLIT_DIR / "plasmode",
    plasmode_dir: Path = ssc.PLASMODE_DIR,
    replicates: int = ssc.PLASMODE_REPLICATES,
    permutations: int = ssc.DEFAULT_PERMUTATIONS,
    cells: Sequence[float] = ssc.PLASMODE_CELLS,
    unequal_size: int = UNEQUAL_SIZE,
) -> dict[str, Any]:
    """A plasmode profile's settings: the parent's at this run's grade, its committed ones and the export checked."""

    if profile not in PLASMODE_PROFILES:
        raise ValueError(f"Unknown plasmode profile {profile!r}; allowed: {list(PLASMODE_PROFILES)}.")
    committed = ssc.load_parent(parent_dir)
    stored = committed["settings"]
    # The parent's fidelity reading: the gate checks the latent space, which the scan never reads.
    fidelity = dict(stored["span_split"]["fidelity"])
    rebuilt = isc.plasmode_settings(
        export,
        fidelity,
        parent_dir=span_split_dir,
        plasmode_dir=plasmode_dir,
        replicates=int(stored["replicates"]),
        permutations=int(stored["permutations"]),
        cells=[float(effect) for effect in stored["cells"]],
    )
    _require_equal(rebuilt, stored, Path(parent_dir) / "settings.json")
    current = isc.plasmode_settings(
        export,
        fidelity,
        parent_dir=span_split_dir,
        plasmode_dir=plasmode_dir,
        replicates=replicates,
        permutations=permutations,
        cells=cells,
    )
    settings: dict[str, Any] = {
        "driver_version": DRIVER_VERSION,
        "profile": profile,
        "parent": _parent_block(parent_dir, profile, committed),
        "interaction_scan": current,
        "cells": [float(effect) for effect in cells],
        "replicates": int(replicates),
        "permutations": int(permutations),
        "alpha": ALPHA,
        "arms": list(ARMS),
        "rejection": SCAN_RULE,
        "scan": _scan_description(
            features="the analysis's unweighted pooled-standardized joint features (no latent fit)",
            models="the covariate-adjusted full and reduced trajectory models",
            seed="the unit's RRPP seed",
            residual_correction="df",
        ),
    }
    if profile == "plasmode_unequal":
        base_seed = int(current["span_split"]["plasmode"]["base_seed"])
        settings["subset"] = unequal_subset(export, base_seed=base_seed, size=unequal_size)
    return settings


def simulation_units(settings: Mapping[str, Any]) -> list[lsc.Unit]:
    return isc.simulation_units(settings["interaction_scan"])


def plasmode_units(settings: Mapping[str, Any]) -> list[tuple[float, int]]:
    return isc.plasmode_units(settings["interaction_scan"])


# ── One unit: the two arms ───────────────────────────────────────────────────


def block_labels(columns: Sequence[Any]) -> np.ndarray:
    return isc._block_labels([str(column) for column in columns])


def run_arms(
    features: np.ndarray,
    full: np.ndarray,
    reduced: np.ndarray,
    labels: np.ndarray,
    *,
    permutations: int,
    seed: int | None,
    residual_correction: str | None,
) -> dict[str, Any]:
    """Both arms on the same inputs and draws; they must share every F."""

    from motco.stats.interaction import interaction_feature_scan

    scans = {
        arm: interaction_feature_scan(
            features,
            full,
            reduced,
            permutations=permutations,
            seed=seed,
            residual_correction=residual_correction,
            blocks=None if arm == "joint" else list(labels),
        )
        for arm in ARMS
    }
    if not np.array_equal(scans["joint"].f, scans["by_block"].f, equal_nan=True):
        raise pc.ConsistencyError("arms: the joint and by-block scans do not share every F.")
    return scans


def _hits(scan: Any) -> np.ndarray:
    return scan.hits(ALPHA) if scan.p_value_adjusted is not None else np.array([], dtype=int)


def block_maxima(scan: Any) -> dict[str, Any]:
    """The by-block scan's null block maxima: their Spearman correlation (two blocks) and the minP global p-value."""

    from scipy.stats import spearmanr

    tested = [block for block in scan.blocks if block.tested]
    null = scan.null_block_max_f
    if null is None or not tested:
        return {"n_blocks_tested": len(tested), "spearman": None, "minp_p_value_global": None}
    draws = null.shape[0]
    observed_p = np.empty(len(tested))
    null_q = np.empty_like(null)
    for slot, block in enumerate(tested):
        values = null[:, slot]
        ordered = np.sort(values)
        at_least = draws - np.searchsorted(ordered, values, side="left")  # #{draws ≥ each draw's value}
        observed_p[slot] = (1.0 + draws - np.searchsorted(ordered, block.max_f, side="left")) / (1.0 + draws)
        null_q[:, slot] = (at_least + (block.max_f >= values)) / (1.0 + draws)
    smallest = float(observed_p.min())
    minp = float((1.0 + np.sum(null_q.min(axis=1) <= smallest)) / (1.0 + draws))
    rho = None
    if len(tested) == 2 and np.ptp(null[:, 0]) > 0 and np.ptp(null[:, 1]) > 0:
        rho = float(spearmanr(null[:, 0], null[:, 1]).statistic)
    return {"n_blocks_tested": len(tested), "spearman": rho, "minp_p_value_global": minp}


def scan_record(scans: Mapping[str, Any], labels: np.ndarray) -> dict[str, Any]:
    """The arm-independent quantities: counts, df, the maxima and the joint arm's per-block null-maximum shares."""

    joint, by_block = scans["joint"], scans["by_block"]
    names = list(dict.fromkeys(labels.tolist()))
    entries = {block.name: block for block in by_block.blocks}
    blocks: dict[str, dict[str, Any]] = {}
    for name in names:
        member = labels == name
        entry: dict[str, Any] = {
            "n_features": int(np.count_nonzero(member)),
            "n_degenerate": int(np.count_nonzero(joint.degenerate[member])),
            "max_f": entries[name].max_f,
        }
        if joint.null_max_column is not None:
            entry["null_max_share"] = float(np.mean(labels[joint.null_max_column] == name))
        blocks[name] = entry
    return {
        "n_features": joint.n_features,
        "n_degenerate": joint.n_degenerate,
        "df_interaction": joint.df_interaction,
        "df_residual": joint.df_residual,
        "max_f": joint.max_f,
        "max_block": str(labels[joint.max_column]),
        "blocks": blocks,
    }


def arm_record(scan: Any, labels: np.ndarray, *, masks: Mapping[str, np.ndarray] | None = None) -> dict[str, Any]:
    """One arm's global p-value and hits, in total and per block, and its hits on each named feature mask."""

    masks = masks or {}
    hits = _hits(scan)
    hit = np.zeros(labels.size, dtype=bool)
    hit[hits] = True
    names = list(dict.fromkeys(labels.tolist()))
    record: dict[str, Any] = {"p_value_global": scan.p_value_global, "hits": int(hit.sum())}
    for key, mask in masks.items():
        record[key] = int(np.count_nonzero(hit & mask))
    blocks: dict[str, dict[str, Any]] = {}
    by_name = {} if scan.blocks is None else {block.name: block for block in scan.blocks}
    for name in names:
        member = labels == name
        entry: dict[str, Any] = {"hits": int(np.count_nonzero(hit & member))}
        for key, mask in masks.items():
            entry[key] = int(np.count_nonzero(hit & member & mask))
        if name in by_name and scan.p_value_adjusted is not None:
            entry["p_value"] = by_name[name].p_value
        blocks[name] = entry
    record["blocks"] = blocks
    if scan.n_blocks_tested is not None:
        record["n_blocks_tested"] = scan.n_blocks_tested
    return record


def _base_simulation(unit: lsc.Unit, settings: Mapping[str, Any]) -> dict[str, Any]:
    comparison = settings["interaction_scan"]["span_split"]["comparison"]
    cell = lsc.cell_spec(comparison, unit.cell)
    return {
        "profile": "simulation",
        "design": unit.design,
        "cell": unit.cell,
        "trajectory_mode": cell["mode"],
        "effect_size": float(cell["effect"]),
        "span_alignment": cell["phi"],
        "span_support": cell["support"],
        "replicate": unit.replicate,
        "generator_seed": lsc.generator_seed(comparison["base_seed"], unit.design, unit.replicate),
        "rrpp_seed": lsc.rrpp_seed(comparison["base_seed"], unit.design, unit.cell, unit.replicate),
        "settings_hash": lsc.settings_hash(settings),
    }


def simulation_inputs(dataset: Any, params: Any) -> tuple[Any, np.ndarray, np.ndarray]:
    """The harness's pooled-standardized features on the measured blocks and the trajectory models; no fit."""

    from motco.simulations.evaluation import build_simulation_trajectory_design
    from motco.simulations.preprocessing import concatenate_blocks, fit_omics_preprocessor, selected_layers

    preprocessor = fit_omics_preprocessor(dataset, selected_layers(params.integration_params))
    features = concatenate_blocks(preprocessor.transform_dataset(dataset)).reset_index(drop=True)
    design = build_simulation_trajectory_design(
        dataset.metadata, group_col=params.group_col, stage_col=params.stage_col
    )
    return features, np.asarray(design.model_full, dtype=float), np.asarray(design.model_reduced, dtype=float)


def evaluate_simulation_unit(unit: lsc.Unit, settings: Mapping[str, Any]) -> dict[str, Any]:
    """Both arms on one simulation unit's standardized features at its RRPP seed, uncorrected."""

    start = perf_counter()
    span_split = settings["interaction_scan"]["span_split"]
    comparison = span_split["comparison"]
    base = _base_simulation(unit, settings)
    timing: dict[str, float] = {}
    try:
        dataset = lsc.unit_dataset(unit, comparison)
        n_stages = int(comparison["designs"][unit.design]["generator"]["n_stages"])
        params = ssc.production_params(span_split, seed=base["rrpp_seed"], n_stages=n_stages)
        features, full, reduced = simulation_inputs(dataset, params)
        timing["inputs"] = perf_counter() - start
        columns = [str(column) for column in features.columns]
        labels = block_labels(columns)
        mark = perf_counter()
        scans = run_arms(
            features.to_numpy(dtype=float), full, reduced, labels,
            permutations=int(settings["permutations"]), seed=base["rrpp_seed"], residual_correction=None,
        )
        timing["scans"] = perf_counter() - mark
        null_mask = isc.population_null_mask(
            dataset, columns, params.group_col, float(settings["population_null_tolerance"])
        )
        summary = scan_record(scans, labels)
        masks = {} if null_mask is None else {"hits_population_null": null_mask}
        summary["n_population_null"] = None if null_mask is None else int(null_mask.sum())
        for name, entry in summary["blocks"].items():
            entry["n_population_null"] = None if null_mask is None else int(null_mask[labels == name].sum())
        arms = {arm: arm_record(scans[arm], labels, masks=masks) for arm in ARMS}
        maxima = block_maxima(scans["by_block"])
    except Exception as exc:  # recorded, never aborts the shard
        timing["unit"] = perf_counter() - start
        return {**base, "status": "failed", "error_type": type(exc).__name__, "error_message": str(exc),
                "runtime_seconds": timing}
    timing["unit"] = perf_counter() - start
    return lsc._jsonable({**base, "status": "completed", "scan": summary, "arms": arms, "block_maxima": maxima,
                          "runtime_seconds": timing})


def plasmode_dataset(blocks: Mapping[str, Any], metadata: Any, export: Any) -> Any:
    """The case-study dataset of one assignment, as ``plasmode_check.analyse`` builds it."""

    from motco.analysis import OmicsBlock, TrajectoryDataset

    return TrajectoryDataset(
        blocks=[OmicsBlock(name, blocks[name], "identity") for name in pc.BLOCKS],
        metadata=metadata,
        sample_col=export.sample_col,
        group_col=export.group_col,
        stage_col=export.stage_col,
        stage_order=list(export.stage_order),
        covariate_cols=list(pc.COVARIATES),
        missing_covariates=pc.MISSING_COVARIATES,
    )


def _scan_table(scan: Any, columns: Sequence[str]) -> Any:
    """The scan-table columns ``interaction_scan_check.planted_counts`` reads."""

    import pandas as pd

    split = [str(column).split("__", 1) for column in columns]
    return pd.DataFrame(
        {
            "block": [parts[0] for parts in split],
            "feature": [parts[1] for parts in split],
            "f": scan.f,
            "p_value_adjusted": scan.p_value_adjusted,
            "degenerate": scan.degenerate,
        }
    )


_SHARED_PLANTED = ("top_k_in_planted", "planted")
_ARM_PLANTED = ("hits", "hits_in_planted", "hits_outside_planted")


def planted_record(scans: Mapping[str, Any], columns: Sequence[str], plant: Mapping[str, Any]) -> dict[str, Any]:
    """Per arm, the hits inside and outside the planted set (total and per block); once, the top-k count."""

    counts = {arm: isc.planted_counts(_scan_table(scans[arm], columns), plant, top_k=TOP_K) for arm in ARMS}
    joint = counts["joint"]
    shared = {
        "top_k": joint["top_k"],
        **{key: joint[key] for key in _SHARED_PLANTED},
        "blocks": {name: {key: entry[key] for key in _SHARED_PLANTED} for name, entry in joint["blocks"].items()},
    }
    per_arm = {
        arm: {
            **{key: counts[arm][key] for key in _ARM_PLANTED},
            "blocks": {
                name: {key: entry[key] for key in _ARM_PLANTED} for name, entry in counts[arm]["blocks"].items()
            },
        }
        for arm in ARMS
    }
    return {"shared": shared, "arms": per_arm}


def evaluate_plasmode_unit(effect: float, replicate: int, settings: Mapping[str, Any], export: Any) -> dict[str, Any]:
    """Both arms on one plasmode unit's analysis scan inputs, with ``df``; seeds, statistics and counts only.

    ``export`` is the export the units are built from (reduced for ``plasmode_unequal``).
    """

    from motco import analysis

    start = perf_counter()
    scan_settings = settings["interaction_scan"]
    parent = scan_settings["span_split"]["plasmode"]
    base = {
        "profile": settings["profile"],
        "key": pc.unit_key(effect, replicate),
        "effect": float(effect),
        "cell": pc.cell_label(effect),
        "replicate": int(replicate),
        "stream": parent["stream"],
        "seeds": pc.unit_seeds(int(parent["base_seed"]), str(parent["stream"]), int(replicate)),
        "settings_hash": lsc.settings_hash(settings),
    }
    timing: dict[str, float] = {}
    try:
        built = pc.build_unit(effect, replicate, parent, export)
        params = isc.plasmode_params(scan_settings, built["seeds"])
        inputs = analysis._interaction_scan_inputs(plasmode_dataset(built["blocks"], built["metadata"], export), params)
        timing["inputs"] = perf_counter() - start
        columns = [str(column) for column in inputs.features.columns]
        labels = block_labels(columns)
        mark = perf_counter()
        scans = run_arms(
            inputs.features.to_numpy(dtype=float), inputs.model_full, inputs.model_reduced, labels,
            permutations=params.permutations, seed=params.seed,
            residual_correction=params.residual_correction if params.permutations > 0 else None,
        )
        timing["scans"] = perf_counter() - mark
        summary = scan_record(scans, labels)
        planted = planted_record(scans, columns, built["plant"])
        arms = {arm: arm_record(scans[arm], labels) for arm in ARMS}
        for arm in ARMS:
            arms[arm]["planted"] = planted["arms"][arm]
        summary["planted"] = planted["shared"]
        maxima = block_maxima(scans["by_block"])
    except Exception as exc:  # recorded, never aborts the shard
        timing["unit"] = perf_counter() - start
        return {**base, "status": "failed", "error_type": type(exc).__name__, "error_message": str(exc),
                "runtime_seconds": timing}
    timing["unit"] = perf_counter() - start
    return lsc._jsonable(
        {
            **base,
            "status": "completed",
            "seeds": built["seeds"],
            "label_redraws": built["label_redraws"],
            "scan": summary,
            "arms": arms,
            "block_maxima": maxima,
            "runtime_seconds": timing,
        }
    )


# ── Shards, resume and the run ───────────────────────────────────────────────


_WORKER_EXPORT: Any = None


def unit_export(export: Any, settings: Mapping[str, Any]) -> Any:
    """The export the profile's units are built from: reduced for ``plasmode_unequal``."""

    return reduced_export(export, settings["subset"]) if "subset" in settings else export


def _init_plasmode_worker(directory: str, files: Mapping[str, str], subset: Mapping[str, Any] | None) -> None:
    global _WORKER_EXPORT
    export = pc.load_export(Path(directory))
    if export.files != dict(files):
        raise SystemExit(f"{directory} changed after the run started (file hashes differ).")
    _WORKER_EXPORT = export if subset is None else reduced_export(export, subset)


def _evaluate_plasmode_pending(unit: tuple[float, int], settings: Mapping[str, Any]) -> dict[str, Any]:
    return evaluate_plasmode_unit(*unit, settings, _WORKER_EXPORT)


def _keys(settings: Mapping[str, Any]) -> tuple[list[tuple[str, Any]], Callable[[Mapping[str, Any]], str]]:
    if settings["profile"] == "simulation":
        return [(ssc._unit_key(unit), unit) for unit in simulation_units(settings)], ssc._simulation_key
    return [(pc.unit_key(*unit), unit) for unit in plasmode_units(settings)], ssc._plasmode_key


def _export_files(settings: Mapping[str, Any]) -> Mapping[str, str]:
    return settings["interaction_scan"]["span_split"]["plasmode"]["export"]["files"]


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
    """Evaluate this shard's missing units into its file; returns how many were evaluated.

    ``export`` is the full export as loaded under the contract; the unequal profile reduces it here.
    """

    profile = settings["profile"]
    check_parent_files(settings)
    directory = profile_dir(out_dir, profile)
    lsc.ensure_settings(directory, settings)
    path = lsc.shard_path(directory, index, count)
    path.parent.mkdir(parents=True, exist_ok=True)
    if profile in PLASMODE_PROFILES and (export is None or export.files != _export_files(settings)):
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
        units_from = unit_export(export, settings)
        for position, unit in enumerate(pending, start=1):
            append(position, evaluate_plasmode_unit(*unit, settings, units_from))
    else:
        initargs = (str(export.directory), dict(export.files), settings.get("subset"))
        with ProcessPoolExecutor(max_workers=workers, initializer=_init_plasmode_worker, initargs=initargs) as pool:
            for position, record in enumerate(
                pool.map(_evaluate_plasmode_pending, pending, [settings] * len(pending)), start=1
            ):
                append(position, record)
    return len(pending)


# ── Merge ────────────────────────────────────────────────────────────────────


SIMULATION_CELL = isc.SIMULATION_CELL
PLASMODE_CELL = isc.PLASMODE_CELL
RATE = ("rate", "se", "count", "rejections", "failures")


def _cell_columns(profile: str) -> tuple[str, ...]:
    return SIMULATION_CELL if profile == "simulation" else PLASMODE_CELL


def _short_cell(profile: str) -> tuple[str, ...]:
    return ("design", "cell") if profile == "simulation" else ("cell", "effect")


def table_columns(profile: str) -> dict[str, tuple[str, ...]]:
    cell, short = _cell_columns(profile), _short_cell(profile)
    block = (*short, "arm", "block", "n", "feature_share", "null_max_share_mean", "hit_unit_rate", "hit_unit_se",
             "hits_mean")
    if profile == "simulation":
        block = (*block, "population_null_unit_rate")
        false_hits = (*short, "arm", "n_population_null_median", *RATE)
    else:
        block = (*block, "outside_unit_rate", "recall_mean", "recall_se", "hit_precision_mean", "hit_precision_n")
        false_hits = (*short, "arm", *RATE)
    columns = {
        "scan_table.csv": (*cell, "arm", *RATE, "hits_median"),
        "block_table.csv": block,
        "false_hit_table.csv": false_hits,
        "paired_table.csv": (*short, "n", "both", "joint_only", "by_block_only", "neither", "p_loss", "p_gain"),
        "minp_table.csv": (*short, "count", "joint_rate", "bonferroni_rate", "minp_rate", "minp_rejections",
                           "spearman_median"),
    }
    if profile != "simulation":
        columns["allocation_table.csv"] = (*short, "block", "n", "gains", "losses", "ties", "mean_difference",
                                           "p_sign")
    return columns


def _cell_groups(records: Sequence[Mapping[str, Any]], settings: Mapping[str, Any]) -> list[tuple[dict, list]]:
    return isc._cell_groups(records, settings["interaction_scan"])


def _completed(group: Sequence[Mapping[str, Any]]) -> list[Mapping[str, Any]]:
    return [record for record in group if record.get("status") == "completed"]


def _planted_cell(fields: Mapping[str, Any]) -> bool:
    return "effect" in fields and float(fields["effect"]) != pc.ANCHOR


def scan_rows(records: Sequence[Mapping[str, Any]], settings: Mapping[str, Any]) -> list[dict[str, Any]]:
    rows = []
    for fields, group in _cell_groups(records, settings):
        completed = _completed(group)
        for arm in ARMS:
            rows.append(
                {
                    **fields,
                    "arm": arm,
                    **ssc._rate_row([_rejects(r["arms"][arm]["p_value_global"]) for r in completed]),
                    "failures": len(group) - len(completed),
                    "hits_median": _median(r["arms"][arm]["hits"] for r in completed),
                }
            )
    return rows


def _block_names(completed: Sequence[Mapping[str, Any]]) -> list[str]:
    names: list[str] = []
    for record in completed:
        names.extend(name for name in record["scan"]["blocks"] if name not in names)
    return names


def block_rows(records: Sequence[Mapping[str, Any]], settings: Mapping[str, Any]) -> list[dict[str, Any]]:
    profile = settings["profile"]
    rows = []
    for fields, group in _cell_groups(records, settings):
        completed = _completed(group)
        short = {key: fields[key] for key in _short_cell(profile)}
        for arm in ARMS:
            for name in _block_names(completed):
                present = [r for r in completed if name in r["scan"]["blocks"]]
                entries = [r["arms"][arm]["blocks"][name] for r in present]
                hit_units = ssc._rate_row([entry["hits"] > 0 for entry in entries])
                row: dict[str, Any] = {
                    **short,
                    "arm": arm,
                    "block": name,
                    "n": len(present),
                    "feature_share": _median(r["scan"]["blocks"][name]["n_features"] / r["scan"]["n_features"]
                                             for r in present),
                    "null_max_share_mean": _mean_se(r["scan"]["blocks"][name].get("null_max_share")
                                                    for r in present)[1],
                    "hit_unit_rate": hit_units["rate"],
                    "hit_unit_se": hit_units["se"],
                    "hits_mean": _mean_se(entry["hits"] for entry in entries)[1],
                }
                if profile == "simulation":
                    with_null = [(r, e) for r, e in zip(present, entries)
                                 if (r["scan"]["blocks"][name].get("n_population_null") or 0) > 0]
                    row["population_null_unit_rate"] = ssc._rate_row(
                        [int(e.get("hits_population_null") or 0) > 0 for _, e in with_null]
                    )["rate"]
                elif _planted_cell(fields):
                    planted = [(r["arms"][arm]["planted"]["blocks"][name], r["scan"]["planted"]["blocks"][name])
                               for r in present]
                    outside = ssc._rate_row([p["hits_outside_planted"] > 0 for p, _ in planted])
                    row["outside_unit_rate"] = outside["rate"]
                    _, row["recall_mean"], row["recall_se"] = _mean_se(
                        p["hits_in_planted"] / s["planted"] for p, s in planted if s["planted"]
                    )
                    row["hit_precision_n"], row["hit_precision_mean"], _ = _mean_se(
                        p["hits_in_planted"] / p["hits"] for p, _ in planted if p["hits"]
                    )
                rows.append(row)
    return rows


def false_hit_rows(records: Sequence[Mapping[str, Any]], settings: Mapping[str, Any]) -> list[dict[str, Any]]:
    """Units with a hit on a population-null feature (simulation) or outside the planted set (plasmode cells)."""

    profile = settings["profile"]
    rows = []
    for fields, group in _cell_groups(records, settings):
        if profile != "simulation" and not _planted_cell(fields):
            continue
        completed = _completed(group)
        short = {key: fields[key] for key in _short_cell(profile)}
        for arm in ARMS:
            if profile == "simulation":
                with_null = [r for r in completed if (r["scan"].get("n_population_null") or 0) > 0]
                flags = [int(r["arms"][arm].get("hits_population_null") or 0) > 0 for r in with_null]
                extra = {"n_population_null_median": _median(r["scan"].get("n_population_null") for r in completed)}
            else:
                flags = [int(r["arms"][arm]["planted"]["hits_outside_planted"]) > 0 for r in completed]
                extra = {}
            rows.append({**short, "arm": arm, **extra, **ssc._rate_row(flags),
                         "failures": len(group) - len(completed)})
    return rows


def paired_rows(records: Sequence[Mapping[str, Any]], settings: Mapping[str, Any]) -> list[dict[str, Any]]:
    """Per cell, the arms' discordant global rejections and the exact one-sided loss and gain p-values."""

    rows = []
    for fields, group in _cell_groups(records, settings):
        pairs = []
        for record in _completed(group):
            joint = _rejects(record["arms"]["joint"]["p_value_global"])
            by_block = _rejects(record["arms"]["by_block"]["p_value_global"])
            if joint is not None and by_block is not None:
                pairs.append((joint, by_block))
        joint_only = sum(1 for joint, by_block in pairs if joint and not by_block)
        by_block_only = sum(1 for joint, by_block in pairs if by_block and not joint)
        discordant = joint_only + by_block_only
        rows.append(
            {
                **{key: fields[key] for key in _short_cell(settings["profile"])},
                "n": len(pairs),
                "both": sum(1 for joint, by_block in pairs if joint and by_block),
                "joint_only": joint_only,
                "by_block_only": by_block_only,
                "neither": sum(1 for joint, by_block in pairs if not joint and not by_block),
                "p_loss": _upper_binomial(joint_only, discordant),
                "p_gain": _upper_binomial(by_block_only, discordant),
            }
        )
    return rows


def minp_rows(records: Sequence[Mapping[str, Any]], settings: Mapping[str, Any]) -> list[dict[str, Any]]:
    rows = []
    for fields, group in _cell_groups(records, settings):
        completed = [r for r in _completed(group) if r["block_maxima"].get("minp_p_value_global") is not None]
        minp = ssc._rate_row([_rejects(r["block_maxima"]["minp_p_value_global"]) for r in completed])
        rows.append(
            {
                **{key: fields[key] for key in _short_cell(settings["profile"])},
                "count": len(completed),
                **{
                    f"{'joint' if arm == 'joint' else 'bonferroni'}_rate": ssc._rate_row(
                        [_rejects(r["arms"][arm]["p_value_global"]) for r in completed]
                    )["rate"]
                    for arm in ARMS
                },
                "minp_rate": minp["rate"],
                "minp_rejections": minp["rejections"],
                "spearman_median": _median(r["block_maxima"].get("spearman") for r in completed),
            }
        )
    return rows


def allocation_rows(records: Sequence[Mapping[str, Any]], settings: Mapping[str, Any]) -> list[dict[str, Any]]:
    """Plasmode planted cells: per block, the per-unit gain in planted hits (by block − joint) and its sign test."""

    rows = []
    for fields, group in _cell_groups(records, settings):
        if not _planted_cell(fields):
            continue
        completed = _completed(group)
        for name in _block_names(completed):
            differences = [
                r["arms"]["by_block"]["planted"]["blocks"][name]["hits_in_planted"]
                - r["arms"]["joint"]["planted"]["blocks"][name]["hits_in_planted"]
                for r in completed
                if name in r["scan"]["blocks"]
            ]
            gains = sum(1 for d in differences if d > 0)
            losses = sum(1 for d in differences if d < 0)
            rows.append(
                {
                    **{key: fields[key] for key in _short_cell(settings["profile"])},
                    "block": name,
                    "n": len(differences),
                    "gains": gains,
                    "losses": losses,
                    "ties": len(differences) - gains - losses,
                    "mean_difference": float(np.mean(differences)) if differences else None,
                    "p_sign": _upper_binomial(gains, gains + losses),
                }
            )
    return rows


def parent_view(record: Mapping[str, Any]) -> dict[str, Any]:
    """The joint arm in the parent driver's record shape, for its table functions."""

    if record.get("status") != "completed":
        return dict(record)
    joint = record["arms"]["joint"]
    scan = record["scan"]
    view_scan: dict[str, Any] = {
        "n_features": scan["n_features"],
        "max_f": scan["max_f"],
        "p_value_global": joint["p_value_global"],
        "hits": joint["hits"],
        "blocks": {
            name: {
                "n_features": entry["n_features"],
                "n_degenerate": entry["n_degenerate"],
                "hits": joint["blocks"][name]["hits"],
                **({"null_max_share": entry["null_max_share"]} if "null_max_share" in entry else {}),
            }
            for name, entry in scan["blocks"].items()
        },
    }
    view = {key: value for key, value in record.items() if key not in ("scan", "arms", "block_maxima")}
    if record.get("profile") == "simulation":
        view_scan["n_population_null"] = scan.get("n_population_null")
        view_scan["hits_population_null"] = joint.get("hits_population_null")
    else:
        planted = scan["planted"]
        view["planted"] = {"top_k": planted["top_k"], "planted": planted["planted"],
                           "top_k_in_planted": planted["top_k_in_planted"], **joint["planted"]}
    view["scan"] = view_scan
    view["profile"] = PARENT_PROFILE[str(record["profile"])]
    return view


def _parent_tables(records: Sequence[Mapping[str, Any]], settings: Mapping[str, Any]) -> dict[str, list[dict]]:
    """The joint arm's rows of the parent's tables, by the parent's own functions."""

    views = [parent_view(record) for record in records]
    parent = settings["interaction_scan"]
    tables = {"scan_table.csv": isc.scan_rows(views, parent), "block_balance.csv": isc.block_rows(views, parent)}
    if settings["profile"] == "simulation":
        tables["false_hit_table.csv"] = isc.false_hit_rows(views, parent)
    else:
        tables["recovery_table.csv"] = isc.recovery_rows(views, parent)
    return tables


def _row_key(profile: str, table: str, row: Mapping[str, Any]) -> tuple[str, ...]:
    place = (str(row["design"]), str(row["cell"])) if profile == "simulation" else (str(row["cell"]),)
    return (*place, str(row["block"])) if table == "block_balance.csv" else place


def reproduction(
    records: Sequence[Mapping[str, Any]],
    settings: Mapping[str, Any],
    committed: Mapping[str, Sequence[Mapping[str, Any]]],
) -> dict[str, Any]:
    """The joint arm against the parent's committed tables, exact at their written precision."""

    profile = settings["profile"]
    parent = settings["interaction_scan"]
    in_run = (
        set(parent["span_split"]["comparison"]["designs"])
        if profile == "simulation"
        else {pc.cell_label(float(effect)) for effect in parent["cells"]}
    )
    mine = _parent_tables(records, settings)
    differing = []
    compared = 0
    for table, columns in REPRODUCED_COLUMNS[profile].items():
        reproduced = {_row_key(profile, table, row): row for row in mine[table]}
        for row in committed[table]:
            if (row["design"] if profile == "simulation" else row["cell"]) not in in_run:
                continue
            key = _row_key(profile, table, row)
            compared += 1
            got = reproduced.get(key)
            wanted = {column: str(row[column]) for column in columns}
            found = None if got is None else {column: _written(got.get(column)) for column in columns}
            if found != wanted:
                differing.append({"table": table, "cell": list(key), "committed": wanted, "reproduced": found})
    return {
        "profile": profile,
        "passed": compared > 0 and not differing,
        "rows_compared": compared,
        "differing": differing,
        "rule": SCAN_RULE,
        "arm": "joint",
        "sources": {name: settings["parent"]["files"][name] for name in REPRODUCED_COLUMNS[profile]},
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
        runtime = record.get("runtime_seconds") or {}
        unit_seconds.append(float(runtime.get("unit") or 0.0))
        for stage in ("inputs", "scans"):
            if runtime.get(stage) is not None:
                stages.setdefault(stage, []).append(float(runtime[stage]))
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
    """Merge one profile's shards and write its tables, accounting and (for two profiles) reproduction check."""

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
    builders: dict[str, Callable[..., list[dict[str, Any]]]] = {
        "scan_table.csv": scan_rows,
        "block_table.csv": block_rows,
        "false_hit_table.csv": false_hit_rows,
        "paired_table.csv": paired_rows,
        "minp_table.csv": minp_rows,
        "allocation_table.csv": allocation_rows,
    }
    for name, columns in table_columns(profile).items():
        _write_csv(directory / name, builders[name](records, settings), columns)
    _write_json(directory / "accounting.json", summary)
    out: dict[str, Any] = {
        "profile": profile,
        **{key: summary[key] for key in ("complete_units", "missing_units", "failed_units", "recorded_core_hours")},
    }
    if profile in REPRODUCED_PROFILES:
        parent = ssc._resolve(settings["parent"]["dir"])
        committed = {name: lsc._read_csv(parent / name) for name in REPRODUCED_COLUMNS[profile]}
        result = reproduction(records, settings, committed)
        _write_json(directory / "reproduction.json", result)
        out.update({"reproduction_passed": result["passed"], "reproduction_differing": len(result["differing"])})
    return out


# ── Decision (pure functions) ────────────────────────────────────────────────


def _index(rows: Iterable[Mapping[str, Any]], *keys: str) -> dict[tuple, Mapping[str, Any]]:
    return {tuple(str(row[key]) for key in keys): row for row in rows}


def _size_observation(row: Mapping[str, Any] | None, **cell: Any) -> dict[str, Any]:
    if row is None or int(row.get("count") or 0) <= 0:
        return {**cell, "status": "unavailable"}
    return {**cell, "n_rejected": int(row["rejections"]), "n_available": int(row["count"])}


def _source(settings: Mapping[str, Any], count: Any) -> dict[str, str]:
    shown = "?" if count in (None, "") else int(float(count))
    return {"source": f"this run, {settings['profile']} profile", "grade": f"{shown} × {settings['permutations']}"}


def decision(
    simulation: Mapping[str, Any],
    plasmode: Mapping[str, Any],
    unequal: Mapping[str, Any],
    *,
    parent_recovery: Sequence[Mapping[str, Any]] = (),
) -> dict[str, Any]:
    """Apply rule version 1 to the three profiles' merged tables.

    Each profile holds ``settings``, ``accounting`` and its tables by file name
    (``scan_table.csv`` …); ``simulation`` and ``plasmode`` also hold
    ``reproduction``. ``parent_recovery`` holds the parent's committed plasmode
    recovery rows, read only by a descriptive table.
    """

    profiles = (("simulation", simulation), ("plasmode", plasmode), ("plasmode_unequal", unequal))
    header = {
        "rule_version": RULE_VERSION,
        "alpha": ALPHA,
        "scan_rejection": SCAN_RULE,
        "settings_hashes": {name: lsc.settings_hash(block["settings"]) for name, block in profiles},
        "complete": all(block["accounting"]["missing_units"] == 0 for _, block in profiles),
        "failed_units": {name: block["accounting"]["failed_units"] for name, block in profiles},
        "reproduction": {
            name: {key: block["reproduction"].get(key) for key in ("passed", "rows_compared", "rule", "source_dir")}
            for name, block in profiles[:2]
        },
    }
    failed = [
        {"profile": name, "differing": block["reproduction"].get("differing", [])}
        for name, block in profiles[:2]
        if not block["reproduction"].get("passed")
    ]
    if failed:
        return {**header, "verdict": "reproduction_failed", "consequence": CONSEQUENCES["reproduction_failed"],
                "failing": {"reproduction": failed}}

    sim_scan = _index(simulation["scan_table.csv"], "design", "cell", "arm")
    designs = list(simulation["settings"]["interaction_scan"]["span_split"]["comparison"]["designs"])
    exports = (("plasmode", plasmode), ("plasmode_unequal", unequal))
    s1 = grade.judge_rates(
        [
            _size_observation(sim_scan.get((design, cell, "by_block")), profile="simulation", design=design, cell=cell)
            for design in designs
            for cell in CONTROL_CELLS
        ]
    )
    s2 = grade.judge_rates(
        [
            _size_observation(_index(block["scan_table.csv"], "cell", "arm").get(("anchor", "by_block")),
                              profile=name, cell="anchor")
            for name, block in exports
        ]
    )
    s3 = grade.judge_rates(
        [
            _size_observation(
                _index(block["false_hit_table.csv"], "cell", "arm").get((pc.cell_label(effect), "by_block")),
                profile=name, cell=pc.cell_label(effect), quantity="units with a hit outside the planted set",
            )
            for name, block in exports
            for effect in PLANTED_EFFECTS
        ]
    )
    p_rows = []
    for name, block in exports:
        paired = _index(block["paired_table.csv"], "cell")
        for effect in PLANTED_EFFECTS:
            row = paired.get((pc.cell_label(effect),))
            entry: dict[str, Any] = {"profile": name, "cell": pc.cell_label(effect)}
            if row is None or int(row.get("n") or 0) <= 0:
                entry.update({"status": "unavailable", "pass": False})
            else:
                p_loss = float(row["p_loss"])
                entry.update({"n": int(row["n"]), "joint_only": int(row["joint_only"]),
                              "by_block_only": int(row["by_block_only"]), "p_loss": p_loss,
                              "pass": p_loss >= POWER_LEVEL})
            p_rows.append(entry)
    p = {"level": POWER_LEVEL, "m": POWER_M, "holds": all(row["pass"] for row in p_rows), "observations": p_rows}
    allocation = _index(unequal.get("allocation_table.csv", []), "cell", "block")
    reduced = str(unequal["settings"]["subset"]["block"])
    g_rows = []
    for effect in PLANTED_EFFECTS:
        row = allocation.get((pc.cell_label(effect), reduced))
        entry = {"profile": "plasmode_unequal", "cell": pc.cell_label(effect), "block": reduced}
        if row is None or int(row.get("n") or 0) <= 0:
            entry.update({"status": "unavailable", "pass": False})
        else:
            p_sign = float(row["p_sign"])
            entry.update({"n": int(row["n"]), "gains": int(row["gains"]), "losses": int(row["losses"]),
                          "ties": int(row["ties"]), "p_sign": p_sign, "pass": p_sign <= GAIN_LEVEL})
        g_rows.append(entry)
    g = {"level": GAIN_LEVEL, "holds": any(row["pass"] for row in g_rows), "observations": g_rows}
    criteria = {
        "S1": {"name": "simulation size", **s1},
        "S2": {"name": "real-noise size", **s2},
        "S3": {"name": "strong familywise control", **s3},
        "P": {"name": "protected power", **p},
        "G": {"name": "allocation gain", **g},
    }
    failing: dict[str, Any] = {}
    for name, judged in (("S1", s1), ("S2", s2), ("S3", s3)):
        if not judged["holds"]:
            failing[name] = isc._failing_observations(judged)
    if not p["holds"]:
        failing["P"] = [row for row in p_rows if not row["pass"]]
    if not g["holds"]:
        failing["G"] = g_rows
    if "S1" in failing or "S2" in failing:
        verdict = "size_excess"
    elif "S3" in failing:
        verdict = "fwer_excess"
    elif "P" in failing:
        verdict = "power_loss"
    elif "G" in failing:
        verdict = "no_allocation_gain"
    else:
        verdict = "by_block_validated"
    return {
        **header,
        "verdict": verdict,
        "consequence": CONSEQUENCES[verdict],
        "failing": failing,
        "criteria": criteria,
        "descriptive": descriptive_tables(simulation, plasmode, unequal, parent_recovery=parent_recovery),
    }


def descriptive_tables(
    simulation: Mapping[str, Any],
    plasmode: Mapping[str, Any],
    unequal: Mapping[str, Any],
    *,
    parent_recovery: Sequence[Mapping[str, Any]] = (),
) -> dict[str, Any]:
    """Block balance, recall per block, the minP comparison, localization and false hits; each row sourced."""

    profiles = (("simulation", simulation), ("plasmode", plasmode), ("plasmode_unequal", unequal))
    balance = []
    for name, block in profiles:
        for row in block["block_table.csv"]:
            if name != "simulation" and _planted_cell(row):
                continue
            balance.append(
                {
                    "profile": name,
                    **({"design": row["design"]} if name == "simulation" else {}),
                    "cell": row["cell"],
                    "arm": row["arm"],
                    "block": row["block"],
                    "feature_share": _number(row.get("feature_share")),
                    "null_max_share": _number(row.get("null_max_share_mean")),
                    "hit_unit_rate": _number(row.get("hit_unit_rate")),
                    "hits_mean": _number(row.get("hits_mean")),
                    **_source(block["settings"], row.get("n")),
                }
            )
    committed = {str(row["cell"]): row for row in parent_recovery}
    recall = []
    for name, block in profiles[1:]:
        for row in block["block_table.csv"]:
            if not _planted_cell(row):
                continue
            parent = committed.get(str(row["cell"])) if name == "plasmode" else None
            recall.append(
                {
                    "profile": name,
                    "cell": row["cell"],
                    "arm": row["arm"],
                    "block": row["block"],
                    "recall": _number(row.get("recall_mean")),
                    "hit_precision": _number(row.get("hit_precision_mean")),
                    "outside_unit_rate": _number(row.get("outside_unit_rate")),
                    **_source(block["settings"], row.get("n")),
                    "parent": None if parent is None else {
                        "recall": _number(parent.get("recall_mean")),
                        "hit_precision": _number(parent.get("hit_precision_mean")),
                        "top_k_precision": _number(parent.get("top_k_precision_mean")),
                        "scope": "both blocks, joint scan",
                        "source": PARENT_SOURCE,
                        "grade": "500 × 999",
                    },
                }
            )
    minp = [
        {
            "profile": name,
            **({"design": row["design"]} if name == "simulation" else {}),
            "cell": row["cell"],
            "joint_rate": _number(row.get("joint_rate")),
            "bonferroni_rate": _number(row.get("bonferroni_rate")),
            "minp_rate": _number(row.get("minp_rate")),
            "spearman_median": _number(row.get("spearman_median")),
            **_source(block["settings"], row.get("count")),
        }
        for name, block in profiles
        for row in block["minp_table.csv"]
    ]
    localization = [
        {
            "design": row["design"],
            "cell": row["cell"],
            "arm": row["arm"],
            "rate": _number(row.get("rate")),
            "se": _number(row.get("se")),
            "hits_median": _number(row.get("hits_median")),
            **_source(simulation["settings"], row.get("count")),
        }
        for row in simulation["scan_table.csv"]
    ]
    false_hits = [
        {
            "design": row["design"],
            "cell": row["cell"],
            "arm": row["arm"],
            "rate": _number(row.get("rate")),
            "count": _integer(row.get("count")),
            **_source(simulation["settings"], row.get("count")),
        }
        for row in simulation["false_hit_table.csv"]
    ]
    return {
        "block_balance": balance,
        "recall_by_block": recall,
        "minp_comparison": minp,
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
    block: dict[str, Any] = {"settings": settings, "accounting": summary}
    for name in table_columns(profile):
        block[name] = _read_csv(directory / name)
    if profile in REPRODUCED_PROFILES:
        block["reproduction"] = json.loads((directory / "reproduction.json").read_text(encoding="utf-8"))
    return block


def decide(out_dir: Path) -> dict[str, Any]:
    """Read the three profiles' merged tables and write ``interaction_scan_blocks_decision.json``."""

    tables = {profile: _profile_tables(profile_dir(out_dir, profile), profile) for profile in PROFILES}
    parent = ssc._resolve(tables["plasmode"]["settings"]["parent"]["dir"])
    result = decision(
        tables["simulation"], tables["plasmode"], tables["plasmode_unequal"],
        parent_recovery=lsc._read_csv(parent / "recovery_table.csv"),
    )
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
                          f"plasmode profiles: default {ssc.PLASMODE_REPLICATES})")
    run.add_argument("--control-reps", type=int, default=ssc.CONTROL_REPLICATES,
                     help="Simulation: replicates of the anchor and translation")
    run.add_argument("--designs", type=str, default=",".join(lsc.DESIGNS), help="Simulation: comma-separated")
    run.add_argument("--parent-dir", type=Path, default=PARENT_DIR,
                     help="The committed interaction-scan validation run "
                          "(default: results/interaction-scan-2026-11-09)")
    run.add_argument("--export-dir", type=Path, default=None, help="Plasmode profiles: the SEA-AD cohort_30 export")
    merge_parser = sub.add_parser("merge", help="Merge one profile's shards")
    merge_parser.add_argument("--profile", choices=PROFILES, required=True)
    merge_parser.add_argument("--out-dir", type=Path, required=True)
    decide_parser = sub.add_parser("decide", help=f"Write {DECISION_NAME} from the three profiles")
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
    if args.export_dir is None:
        raise SystemExit(f"The {args.profile} profile needs --export-dir.")
    try:
        export = pc.load_export(args.export_dir)
    except pc.ExportContractError as exc:
        raise pc.contract_exit(exc) from exc
    settings = plasmode_settings(
        export,
        profile=args.profile,
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
