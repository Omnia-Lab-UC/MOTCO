#!/usr/bin/env python3
"""Univariate interaction baseline: the production trajectory test against the per-feature interaction test (no R).

The driver behind ``openspec/changes/univariate-interaction-baseline/`` (spec
``univariate-baseline-validation``, roadmap open item 13). The first
comparison a reader makes is against the default analysis with no
integration: a group × stage interaction test on each feature, corrected for
multiple testing. This driver judges the production trajectory test against
it on matched simulated units, under a rule fixed before any unit.

**Arms** (every unit generates one dataset, fits production once and reads
every arm at the unit's RRPP seed, so all arms see the same permutations):

- ``production`` — the omnibus of the production test: rejects when the
  smallest of the ``delta``, ``angle`` and ``shape`` p-values is at or below
  0.05 / 3. ``delta``, ``angle`` and ``shape`` are also recorded singly
  (p < 0.05), and ``production_minp`` (the smallest of the three calibrated on
  the unit's joint permutation draws, p ≤ 0.05) is recorded and decides nothing.
- Baseline arms, each one test at 0.05 on the per-feature interaction F of
  the fit's standardized features (the split test's features, in the PLS
  input coordinates; the F is invariant to that per-feature affine map) under
  the trajectory test's full and reduced models: ``scan``, the interaction
  scan's maximum (global p ≤ 0.05); ``fdr``, Benjamini–Hochberg at 0.05 over
  the upper-tail F p-values of the non-degenerate features, rejecting when
  any feature is a discovery; ``total``, the split test's ``total`` part
  (p < 0.05).
- Recorded arms that decide nothing: the split test's ``span`` and
  ``complement`` (p < 0.05), and the trajectory statistics on ``concat`` and on
  PCA at the unit's production rank (``pca_matched``), each with the unit's
  production parameters changed only in the integration method, singly
  (p < 0.05) and as the omnibus (≤ 0.05 / 3). ``--recorded-spaces`` (bound in
  the settings) chooses the spaces; a space left out records its arms as
  unavailable. The default is ``pca_matched`` alone: ``concat`` costs about
  80 % of a unit, because its ``shape`` aligns the stage means in the ambient
  dimension on every draw, so it was dropped at the approval gate (the design's
  first compute lever) and stays available as an option.

**Profiles.**

- ``parent`` — the interaction-scan validation's simulation units for the
  ``n200`` and ``seaad`` designs (two blocks, base seed 2035, fourteen cells;
  the anchor and translation at replicates 0–499, every other cell 0–99,
  999 permutations; 4,400 units), rebuilt through
  ``interaction_scan_check.simulation_settings``. The scan validation's
  settings are rebuilt at its grade and must equal its committed
  ``simulation/settings.json``. Its counts of ``scan``, the three split parts
  and the three trajectory statistics must equal the committed
  ``comparator_table.csv`` exactly; its other readings are descriptive.
- ``stage_span`` — a fresh seed family (base seed 2043, the latent-space
  comparison's ``generator_seed`` and ``rrpp_seed``) on methylation,
  expression and proteomics (658 features; ``layers`` and ``span_layers``
  absent), three stages. Designs: balanced n = 80, 200 and 500 at δ = 2.0 and
  n = 200 with every per-omic δ at 0.5. Thirteen cells each: the anchor,
  translation e = 1.0, magnitude e ∈ {0.05, 0.1, 0.25}, in-span rotation and
  bend e ∈ {0.1, 0.25, 0.5} (these and the controls under the planted scale),
  and the out-of-span rotation e = 0.5 at φ = 90 and its 20-feature form under
  the null scale. Controls at 500 replicates, other cells at 200, 999
  permutations; 12,800 units.
- ``surgery`` — the design-guidance paper-grade slice's three-block columns
  at n = 80, 200 and 500 (``design_guidance_core_{magnitude,geometry}_study.json``
  through ``stage_span_leak_check.column_cells`` and ``unit_params``; the
  replicate seeds, and the evaluation seed the study engine resolves: the
  cell's own when set). Eleven cells each: the anchor (from the geometry
  profile), translation e = 1.0, magnitude e ∈ {0.02, 0.05, 0.1}, and
  orientation and shape e ∈ {0.25, 0.5, 1.0}; 500 replicates × 999
  permutations; 16,500 units. The ``pls_cv`` counts of ``delta``, ``angle``
  and ``shape`` must equal the slice's committed ``design_point_operating.csv``
  tables exactly.

The parents' committed files are bound by sha256 in ``settings.json`` and a
changed file is refused by name in ``run``, ``merge`` and ``decide``. A unit
that raises is recorded as failed with its error; the shard continues. A
record holds seeds, statistics and counts only, never a per-feature array.

Subcommands (every output lands under ``--out-dir``/``<profile>/``, the
decision under ``--out-dir``):

- ``run --profile parent|stage_span|surgery`` — evaluate (a shard of) the
  units into ``records/shard_<i>_of_<N>.jsonl``. The settings and their hash
  go to ``settings.json``; a resume or merge under other settings exits naming
  the differing keys, and a resume skips complete units.
- ``merge --profile …`` — ``merged.jsonl`` (per unit; keep off the
  repository), ``rate_table.csv`` (per design, cell and arm: rule, rate, Monte
  Carlo SE, available and rejected counts, failures), ``paired_table.csv``
  (per design, cell and baseline arm against ``production`` on the units
  available in both: both rates, the discordant counts and the exact paired
  p-values, imported from ``production_latent_space_check``),
  ``cell_table.csv`` (rank quantiles, held-out AUROC, capture shares, median
  scan hits and ``fdr`` discoveries, the realized false discovery proportion on
  population-null features), ``accounting.json`` and, for ``parent`` and
  ``surgery``, ``reproduction.json`` (every differing row with both counts).
- ``decide`` — ``univariate_baseline_decision.json`` (rule version 1) from
  the three merged profiles: the gate ``reproduction_failed``; size of
  ``production``, ``scan``, ``fdr`` and ``total`` at the anchors and
  translations of ``stage_span`` (m = 8) and ``surgery`` (m = 6) under the
  stage-span paper grade's count-plus-ceiling rule (``judge_rates``,
  imported unchanged), an arm admissible when both families hold; each of the
  71 signal cells labelled by exact one-sided paired tests of ``production``
  against each admissible baseline arm (``production_leads`` at 0.025 / 71
  against every arm with a gain of at least 0.10, ``baseline_leads`` at
  0.025 / (71·A) against any arm with a loss of at least 0.10, else
  ``unresolved``); the verdict ``reproduction_failed`` |
  ``production_size_excess`` | ``no_admissible_baseline`` |
  ``production_leads`` | ``baseline_leads`` | ``complementary`` |
  ``unresolved``; the characterization grade ``names_kind`` | ``detects_only``
  | ``no_detected_cell`` (the off-target statistics of the in-span
  ``stage_span`` cells that ``production`` detects in at least 0.80 of units,
  as one ``judge_rates`` family), beside the verdict and never changing it;
  and descriptive tables (the label map by plant family and design, every
  arm's rate per cell, discoveries, the fit and the parent's comparison).

With ``--workers`` above 1, pin BLAS to one thread per process
(``OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1``).

Examples:

    export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1

    # smoke
    python scripts/univariate_baseline_check.py run --profile parent --reps 2 --control-reps 2 \\
        --out-dir /tmp/ubc --workers 7
    python scripts/univariate_baseline_check.py run --profile stage_span --reps 2 --control-reps 2 \\
        --out-dir /tmp/ubc --workers 7
    python scripts/univariate_baseline_check.py run --profile surgery --reps 2 --control-reps 2 \\
        --out-dir /tmp/ubc --workers 7
    python scripts/univariate_baseline_check.py merge --profile parent --out-dir /tmp/ubc   # and the others
    python scripts/univariate_baseline_check.py decide --out-dir /tmp/ubc

    # one SLURM array task (scripts/univariate_baseline_array.sbatch), then merge and decide
"""

from __future__ import annotations

import argparse
import dataclasses
import functools
import json
import math
import sys
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path
from time import perf_counter
from typing import Any, Iterable, Mapping, Sequence

import numpy as np
from numpy.typing import ArrayLike

sys.path.insert(0, str(Path(__file__).resolve().parent))
import interaction_scan_check as isc  # noqa: E402  (the parent run's settings, units, population-null mask)
import latent_space_comparison as lsc  # noqa: E402  (seeds, datasets, settings binding, CSV helpers)
import production_latent_space_check as plsc  # noqa: E402  (the exact paired p-values)
import span_split_check as ssc  # noqa: E402  (parent loading, split inputs, production parameters)
import stage_span_leak_check as slc  # noqa: E402  (study-engine columns and unit parameters)
import stage_span_paper_grade_decision as grade  # noqa: E402  (the count-plus-ceiling rule)

sys.path.pop(0)

DRIVER = "univariate_baseline_check"
DRIVER_VERSION = 1
RULE_VERSION = 1
PROFILES: tuple[str, ...] = ("parent", "stage_span", "surgery")
DECIDED_PROFILES: tuple[str, ...] = ("stage_span", "surgery")
REPRODUCED_PROFILES: tuple[str, ...] = ("parent", "surgery")

ALPHA = 0.05
OMNIBUS_DIVISOR = 3
LABEL_ALPHA = 0.025
MARGIN = 0.10
DETECTION_FLOOR = 0.80
SIGNAL_CELLS_EXPECTED = 71
GRADED_CELLS_EXPECTED = 36
DEFAULT_PERMUTATIONS = 999
POPULATION_NULL_TOLERANCE = isc.POPULATION_NULL_TOLERANCE
#: Tolerance on a rate difference compared with the margin, so 0.6 − 0.5 counts as 0.10.
_MARGIN_TOLERANCE = 1e-9

STATISTICS: tuple[str, ...] = ("delta", "angle", "shape")
TARGETS: dict[str, str] = {"magnitude": "delta", "orientation": "angle", "shape": "shape"}
OFF_TARGET: dict[str, tuple[str, ...]] = {
    mode: tuple(statistic for statistic in STATISTICS if statistic != target) for mode, target in TARGETS.items()
}
CONTROL_CELLS: tuple[str, ...] = ("anchor", "translation")

#: The judged arms: production and the three baseline combinations.
BASELINE_ARMS: tuple[str, ...] = ("scan", "fdr", "total")
JUDGED_ARMS: tuple[str, ...] = ("production", *BASELINE_ARMS)
RECORDED_SPACES: tuple[str, ...] = ("concat", "pca_matched")
#: The spaces a run records unless told otherwise; ``concat`` was dropped at the approval gate (2026-11-23).
DEFAULT_RECORDED_SPACES: tuple[str, ...] = ("pca_matched",)
#: Every arm of the rate table, each under the rule it was validated under (design D4).
ARMS: tuple[str, ...] = (
    "production", *STATISTICS, "production_minp", *BASELINE_ARMS, "span", "complement",
    *(name for space in RECORDED_SPACES for name in (space, *(f"{space}_{s}" for s in STATISTICS))),
)
RULE_STRICT = "p < alpha"
RULE_INCLUSIVE = "p <= alpha"
RULE_OMNIBUS = "min p <= alpha / 3"
RULE_FDR = "min BH-adjusted p <= alpha"
RULES: dict[str, str] = {
    "production": RULE_OMNIBUS,
    **{statistic: RULE_STRICT for statistic in STATISTICS},
    "production_minp": RULE_INCLUSIVE,
    "scan": RULE_INCLUSIVE,
    "fdr": RULE_FDR,
    "total": RULE_STRICT,
    "span": RULE_STRICT,
    "complement": RULE_STRICT,
    **{space: RULE_OMNIBUS for space in RECORDED_SPACES},
    **{f"{space}_{s}": RULE_STRICT for space in RECORDED_SPACES for s in STATISTICS},
}
_LEVELS: dict[str, tuple[float, bool]] = {
    RULE_STRICT: (ALPHA, False),
    RULE_INCLUSIVE: (ALPHA, True),
    RULE_OMNIBUS: (ALPHA / OMNIBUS_DIVISOR, True),
    RULE_FDR: (ALPHA, True),
}
#: The parent's committed comparator tests and the arm each one is.
PARENT_TESTS: dict[str, str] = {
    "scan": "scan", "delta": "delta", "angle": "angle", "shape": "shape",
    "split_total": "total", "split_span": "span", "split_complement": "complement",
}

REPO = Path(__file__).resolve().parents[1]
SCAN_DIR = REPO / "results" / "interaction-scan-2026-11-09" / "simulation"
PARENT_FILE_NAMES: tuple[str, ...] = ("settings.json", "comparator_table.csv")
PARENT_DESIGNS: tuple[str, ...] = ("n200", "seaad")
PARENT_SOURCE = "interaction-scan validation 2026-11-09, simulation profile"

STAGE_SPAN_BASE_SEED = 2043
#: design -> (sample size, every per-omic stage effect δ); balanced groups, three stages.
STAGE_SPAN_DESIGNS: dict[str, tuple[int, float]] = {
    "n80": (80, 2.0), "n200": (200, 2.0), "n500": (500, 2.0), "n200_delta0.5": (200, 0.5),
}
STAGE_EFFECT_KEYS: tuple[str, ...] = ("delta_methyl", "delta_expr", "delta_protein")
PLANTED = "planted"
#: (cell, trajectory mode, effect, φ or None, support or None, span scale or None for the null scale).
STAGE_SPAN_CELLS: tuple[tuple[str, str, float, float | None, int | None, str | None], ...] = (
    ("anchor", "none", 0.0, None, None, PLANTED),
    ("translation", "translation", 1.0, None, None, PLANTED),
    *((f"magnitude_e{effect:g}", "magnitude", effect, None, None, PLANTED) for effect in (0.05, 0.1, 0.25)),
    *((f"orient_e{effect:g}_phi0", "orientation", effect, 0.0, None, PLANTED) for effect in (0.1, 0.25, 0.5)),
    *((f"shape_e{effect:g}", "shape", effect, None, None, PLANTED) for effect in (0.1, 0.25, 0.5)),
    ("orient_e0.5_phi90", "orientation", 0.5, 90.0, None, None),
    ("orient_e0.5_phi90_k20", "orientation", 0.5, 90.0, 20, None),
)

SLICE_DIR = REPO / "results" / "design-guidance-slice-2026-10-29"
CONFIG_DIR = REPO / "examples" / "trajectory_power_study"
#: The slice's two profiles: their study config and committed operating table.
SURGERY_SOURCES: dict[str, dict[str, Path]] = {
    "magnitude": {
        "config": CONFIG_DIR / "design_guidance_core_magnitude_study.json",
        "operating": SLICE_DIR / "core_magnitude" / "report" / "design_point_operating.csv",
    },
    "geometry": {
        "config": CONFIG_DIR / "design_guidance_core_geometry_study.json",
        "operating": SLICE_DIR / "core_geometry" / "report" / "design_point_operating.csv",
    },
}
SURGERY_COLUMNS: tuple[str, ...] = ("n80_three_blocks", "n200_three_blocks", "n500_three_blocks")
#: (cell, trajectory mode, effect, the slice profile it is read from).
SURGERY_CELLS: tuple[tuple[str, str, float, str], ...] = (
    ("anchor", "none", 0.0, "geometry"),
    ("translation", "translation", 1.0, "geometry"),
    *((f"magnitude_e{effect:g}", "magnitude", effect, "magnitude") for effect in (0.02, 0.05, 0.1)),
    *((f"orientation_e{effect:g}", "orientation", effect, "geometry") for effect in (0.25, 0.5, 1.0)),
    *((f"shape_e{effect:g}", "shape", effect, "geometry") for effect in (0.25, 0.5, 1.0)),
)
SURGERY_SOURCE = "design-guidance slice 2026-10-29, core profiles"

DEFAULT_REPLICATES: dict[str, dict[str, int]] = {
    "parent": {"control": ssc.CONTROL_REPLICATES, "other": ssc.SIMULATION_REPLICATES},
    "stage_span": {"control": 500, "other": 200},
    "surgery": {"control": 500, "other": 500},
}

VERDICTS: tuple[str, ...] = (
    "reproduction_failed", "production_size_excess", "no_admissible_baseline",
    "production_leads", "baseline_leads", "complementary", "unresolved",
)
LABELS: tuple[str, ...] = ("production_leads", "baseline_leads", "unresolved")
GRADES: tuple[str, ...] = ("names_kind", "detects_only", "no_detected_cell")
_LIMIT = (
    " The reading holds at 658 simulated features with Gaussian noise; no real-noise reading exists on "
    "methylation, expression and proteomics (roadmap open item 2). No verdict changes what `motco analyze` runs."
)
CONSEQUENCES: dict[str, str] = {
    "reproduction_failed": (
        "Nothing is read: the units do not reproduce the scan validation or the design-guidance slice, so no rate "
        "pairs with the committed ones. Find the drift before any rerun."
    ),
    "production_size_excess": (
        "Nothing is claimed against the default analysis: the production omnibus does not hold size at the "
        "failing controls, so a detection advantage could not be told from excess rejections." + _LIMIT
    ),
    "no_admissible_baseline": (
        "Nothing is claimed against the default analysis: no per-feature combination holds size at these "
        "designs, which is itself the finding the report states." + _LIMIT
    ),
    "production_leads": (
        "The documentation may state that the trajectory test detects differences the default per-feature "
        "analysis does not, at these designs and this dimension, and that the per-feature analysis led in no "
        "cell." + _LIMIT
    ),
    "baseline_leads": (
        "The documentation states that integration gives no detection advantage at these designs: a case study "
        "reads the per-feature test first, and the trajectory statistics are presented as a characterization "
        "only under the grade names_kind." + _LIMIT
    ),
    "complementary": (
        "The documentation states the map of which side leads by plant family and design, and a case study "
        "reports both analyses under a multiplicity policy fixed in advance." + _LIMIT
    ),
    "unresolved": "Nothing is claimed either way: no cell separates the two analyses by the margin." + _LIMIT,
}
DECISION_NAME = "univariate_baseline_decision.json"


# ── Shared helpers ───────────────────────────────────────────────────────────


def _canonical(value: Any) -> Any:
    return json.loads(json.dumps(value, sort_keys=True))


def profile_dir(out_dir: Path, profile: str) -> Path:
    return Path(out_dir) / profile


def _finite(value: Any) -> float | None:
    if value is None or value == "":
        return None
    number = float(value)
    return number if math.isfinite(number) else None


def rejects(p_value: Any, rule: str) -> bool | None:
    """The rejection of one p-value under ``rule``; ``None`` for an unavailable p-value."""

    value = _finite(p_value)
    if value is None:
        return None
    level, inclusive = _LEVELS[rule]
    return value <= level if inclusive else value < level


def _median(values: Iterable[Any]) -> float | None:
    return ssc._median(values)


def _quantile(values: Iterable[Any], q: float) -> float | None:
    finite = [float(value) for value in values if _finite(value) is not None]
    return float(np.quantile(finite, q)) if finite else None


def _mean(values: Iterable[Any]) -> float | None:
    finite = [float(value) for value in values if _finite(value) is not None]
    return float(np.mean(finite)) if finite else None


def _read_csv(path: Path) -> list[dict[str, Any]]:
    return ssc._read_csv(path)


def unit_key(unit: lsc.Unit) -> str:
    return ssc._unit_key(unit)


def record_key(record: Mapping[str, Any]) -> str:
    return ssc._simulation_key(record)


# ── Parent binding ───────────────────────────────────────────────────────────


def _bind(paths: Iterable[Path]) -> dict[str, str]:
    """The sha256 of each committed file, keyed by its repository-relative path."""

    out: dict[str, str] = {}
    for path in paths:
        if not Path(path).exists():
            raise SystemExit(f"{path} does not exist; the committed parent file is required.")
        out[ssc._display_path(path)] = ssc.sha256_file(path)
    return out


def check_parent_files(settings: Mapping[str, Any]) -> None:
    """Refuse, by name, a committed parent file whose sha256 differs from the one the settings record."""

    for stored, digest in settings.get("parent_files", {}).items():
        path = ssc._resolve(stored)
        if not path.exists():
            raise SystemExit(f"{path} does not exist; the committed parent file is required.")
        if ssc.sha256_file(path) != digest:
            raise SystemExit(f"{path} has changed since the settings were written (sha256 differs); refusing.")
    if settings["profile"] == "parent":
        isc.check_parent_files(settings["scan_validation"])


def _require_equal(rebuilt: Mapping[str, Any], committed: Mapping[str, Any], path: Path) -> None:
    if _canonical(dict(rebuilt)) != _canonical(dict(committed)):
        differing = lsc.settings_difference(_canonical(dict(committed)), _canonical(dict(rebuilt)))
        raise SystemExit(
            f"{path} differs from the settings interaction_scan_check rebuilds at its grade; differing: "
            f"{', '.join(differing)}. The units would not be the scan validation's."
        )


def _spaces(recorded_spaces: Sequence[str]) -> list[str]:
    unknown = sorted(set(recorded_spaces) - set(RECORDED_SPACES))
    if unknown:
        raise SystemExit(f"Unknown recorded space(s) {unknown}; expected a subset of {list(RECORDED_SPACES)}.")
    return [space for space in RECORDED_SPACES if space in recorded_spaces]


def _common_settings(
    profile: str, replicates: Mapping[str, int], permutations: int, recorded_spaces: Sequence[str]
) -> dict[str, Any]:
    return {
        "driver": DRIVER,
        "driver_version": DRIVER_VERSION,
        "profile": profile,
        "replicates": {"control": int(replicates["control"]), "other": int(replicates["other"])},
        "control_cells": list(CONTROL_CELLS),
        "permutations": int(permutations),
        "alpha": ALPHA,
        "rules": dict(RULES),
        "arms": list(ARMS),
        "baseline_arms": list(BASELINE_ARMS),
        "baseline": {
            "statistic": "per-feature group x stage interaction F (SS_int / df_int) / (RSS / df_res)",
            "features": "the production fit's standardized features in the PLS input coordinates",
            "models": "the trajectory test's full (group x stage) and reduced (group + stage) models",
            "scan": "motco.stats.interaction.interaction_feature_scan at the unit's RRPP seed, uncorrected",
            "fdr": "scipy.stats.f.sf at the scan's degrees of freedom, scipy.stats.false_discovery_control "
                   "(method 'bh') over the non-degenerate features",
            "split": "motco.stats.interaction.interaction_split_test on qr(x_rotations_) at the unit's RRPP seed, "
                     "uncorrected",
        },
        "recorded_spaces": _spaces(recorded_spaces),
        "recorded_space_definitions": {
            "concat": "the unit's production parameters with integration_method 'concat' on the same blocks",
            "pca_matched": "the unit's production parameters with integration_method 'pca' at the unit's "
                           "production rank on the same blocks",
        },
        "population_null_tolerance": POPULATION_NULL_TOLERANCE,
    }


# ── Parent profile ───────────────────────────────────────────────────────────


def parent_settings(
    *,
    scan_dir: Path = SCAN_DIR,
    replicates: int = DEFAULT_REPLICATES["parent"]["other"],
    control_replicates: int = DEFAULT_REPLICATES["parent"]["control"],
    permutations: int = DEFAULT_PERMUTATIONS,
    designs: Sequence[str] = PARENT_DESIGNS,
    recorded_spaces: Sequence[str] = DEFAULT_RECORDED_SPACES,
) -> dict[str, Any]:
    """The scan validation's simulation units for ``designs``; its committed settings rebuilt and checked."""

    unknown = sorted(set(designs) - set(PARENT_DESIGNS))
    if unknown:
        raise SystemExit(f"Unknown parent design(s) {unknown}; expected a subset of {list(PARENT_DESIGNS)}.")
    scan_dir = Path(scan_dir)
    committed = ssc.load_parent(scan_dir)
    stored = committed["settings"]
    span_split_dir = ssc._resolve(stored["parent"]["dir"])
    rebuilt = isc.simulation_settings(
        parent_dir=span_split_dir,
        replicates=int(stored["replicates"]["other"]),
        control_replicates=int(stored["replicates"]["control"]),
        permutations=int(stored["permutations"]),
        designs=list(stored["span_split"]["comparison"]["designs"]),
    )
    _require_equal(rebuilt, stored, scan_dir / "settings.json")
    scan_validation = isc.simulation_settings(
        parent_dir=span_split_dir,
        replicates=replicates,
        control_replicates=control_replicates,
        permutations=permutations,
        designs=[name for name in PARENT_DESIGNS if name in designs],
    )
    return {
        **_common_settings("parent", scan_validation["replicates"], permutations, recorded_spaces),
        "parent_files": _bind(scan_dir / name for name in PARENT_FILE_NAMES),
        "parent_settings_hash": committed["hash"],
        "parent_source": PARENT_SOURCE,
        "scan_validation": scan_validation,
        "designs": list(scan_validation["span_split"]["comparison"]["designs"]),
        "production": "the split-test validation's pls_cv parameters (two blocks), null distributions kept",
    }


# ── Stage-span profile ───────────────────────────────────────────────────────


def stage_span_generator(n_samples: int, delta: float) -> dict[str, Any]:
    """A balanced three-stage stage-span design on all three blocks (``span_layers`` absent)."""

    return {
        "n_stages": 3,
        "n_samples": int(n_samples),
        "group_ratio": 0.5,
        "p_dmp": 0.1,
        "baseline_continuity": 0.0,
        "magnitude_kind": "joint",
        "construction": "stage_span",
        **{key: float(delta) for key in STAGE_EFFECT_KEYS},
    }


def stage_span_settings(
    *,
    replicates: int = DEFAULT_REPLICATES["stage_span"]["other"],
    control_replicates: int = DEFAULT_REPLICATES["stage_span"]["control"],
    permutations: int = DEFAULT_PERMUTATIONS,
    designs: Sequence[str] = tuple(STAGE_SPAN_DESIGNS),
    base_seed: int = STAGE_SPAN_BASE_SEED,
    recorded_spaces: Sequence[str] = DEFAULT_RECORDED_SPACES,
) -> dict[str, Any]:
    """Four three-block stage-span designs at a fresh base seed, the comparison's seed derivations."""

    unknown = sorted(set(designs) - set(STAGE_SPAN_DESIGNS))
    if unknown:
        raise SystemExit(f"Unknown stage-span design(s) {unknown}; expected a subset of {list(STAGE_SPAN_DESIGNS)}.")
    ordered = [name for name in STAGE_SPAN_DESIGNS if name in designs]
    return {
        **_common_settings("stage_span", {"control": control_replicates, "other": replicates}, permutations,
                           recorded_spaces),
        "parent_files": {},
        "units": {
            "designs": {
                name: {"generator": stage_span_generator(*STAGE_SPAN_DESIGNS[name]), "layers": None}
                for name in ordered
            },
            "cells": [list(cell) for cell in STAGE_SPAN_CELLS],
            "base_seed": int(base_seed),
            "seed_functions": "latent_space_comparison.generator_seed / rrpp_seed",
            "pls_cv_params": dict(lsc.PLS_CV_PARAMS),
            "min_components": lsc.MIN_COMPONENTS,
            "orientation_estimator": lsc.ORIENTATION_ESTIMATOR,
        },
        "designs": ordered,
        "production": "pls with the comparison's CV knobs on all three blocks, rank floor 2, net displacement, "
                      "span capture, null distributions kept",
    }


def stage_span_cell(settings: Mapping[str, Any], name: str) -> dict[str, Any]:
    for cell, mode, effect, phi, support, scale in settings["units"]["cells"]:
        if cell == name:
            return {"cell": cell, "mode": mode, "effect": effect, "phi": phi, "support": support, "span_scale": scale}
    raise KeyError(f"Unknown stage-span cell {name!r}.")


def stage_span_generator_params(unit: lsc.Unit, settings: Mapping[str, Any]) -> Any:
    units = settings["units"]
    cell = stage_span_cell(settings, unit.cell)
    return lsc._generator_params(
        units["designs"][unit.design]["generator"],
        seed=lsc.generator_seed(units["base_seed"], unit.design, unit.replicate),
        trajectory_mode=cell["mode"],
        group_effect_size=float(cell["effect"]),
        span_alignment=cell["phi"],
        span_support=cell["support"],
        span_scale=cell["span_scale"],
    )


def stage_span_params(settings: Mapping[str, Any], *, seed: int) -> Any:
    from motco.simulations.evaluation import SimulationEvaluationParams

    units = settings["units"]
    return SimulationEvaluationParams(
        integration_method="pls",
        integration_params=dict(units["pls_cv_params"]),
        permutations=int(settings["permutations"]),
        n_jobs=1,
        seed=int(seed),
        include_null_distributions=True,
        orientation_estimator=units["orientation_estimator"],
        min_components=int(units["min_components"]),
        span_capture=True,
    )


# ── Surgery profile ──────────────────────────────────────────────────────────


@functools.lru_cache(maxsize=None)
def _study_columns(config: str) -> dict[str, dict[str, Any]]:
    """The study engine's columns of one committed config (``stage_span_leak_check.column_cells``), cached."""

    return slc.column_cells(ssc._resolve(config))


def surgery_cell(settings: Mapping[str, Any], name: str) -> dict[str, Any]:
    for cell, mode, effect, source in settings["cells"]:
        if cell == name:
            return {"cell": cell, "mode": mode, "effect": effect, "source": source}
    raise KeyError(f"Unknown surgery cell {name!r}.")


def surgery_params(unit: lsc.Unit, settings: Mapping[str, Any]) -> tuple[Any, Any]:
    """The study engine's generator and evaluation parameters of one slice unit, null distributions kept.

    ``stage_span_leak_check.unit_params`` sets the replicate seed and the unit's effect; the evaluation seed is
    resolved as the study engine resolves it (the cell's own when set, else the replicate seed).
    """

    cell = surgery_cell(settings, unit.cell)
    columns = _study_columns(settings["sources"][cell["source"]]["config"])
    generator, evaluation = slc.unit_params(
        columns,
        (unit.design, cell["mode"], float(cell["effect"]), int(unit.replicate)),
        span_scale=None,
        permutations=int(settings["permutations"]),
    )
    seed = evaluation.seed if evaluation.seed is not None else generator.seed
    return generator, dataclasses.replace(evaluation, seed=int(seed), include_null_distributions=True)


def surgery_settings(
    *,
    replicates: int = DEFAULT_REPLICATES["surgery"]["other"],
    control_replicates: int = DEFAULT_REPLICATES["surgery"]["control"],
    permutations: int = DEFAULT_PERMUTATIONS,
    columns: Sequence[str] = SURGERY_COLUMNS,
    sources: Mapping[str, Mapping[str, Path]] = SURGERY_SOURCES,
    recorded_spaces: Sequence[str] = DEFAULT_RECORDED_SPACES,
) -> dict[str, Any]:
    """The slice's three-block columns from its committed configs, every column and cell checked present."""

    unknown = sorted(set(columns) - set(SURGERY_COLUMNS))
    if unknown:
        raise SystemExit(f"Unknown surgery column(s) {unknown}; expected a subset of {list(SURGERY_COLUMNS)}.")
    ordered = [name for name in SURGERY_COLUMNS if name in columns]
    bound = {
        name: {"config": ssc._display_path(paths["config"]), "operating": ssc._display_path(paths["operating"])}
        for name, paths in sources.items()
    }
    sizes: dict[str, int] = {}
    for cell, mode, effect, source in SURGERY_CELLS:
        study = _study_columns(bound[source]["config"])
        for column in ordered:
            if column not in study or mode not in study[column]:
                raise SystemExit(f"{bound[source]['config']} has no {mode!r} cell in column {column!r}.")
            sizes[column] = int(study[column][mode].generator_params.n_samples)
    return {
        **_common_settings("surgery", {"control": control_replicates, "other": replicates}, permutations,
                           recorded_spaces),
        "parent_files": _bind(path for paths in sources.values() for path in (paths["config"], paths["operating"])),
        "parent_source": SURGERY_SOURCE,
        "sources": bound,
        "designs": ordered,
        "n_samples": sizes,
        "cells": [list(cell) for cell in SURGERY_CELLS],
        "production": "the study cell's evaluation parameters (pls on all three blocks) with span capture and "
                      "null distributions kept; the evaluation seed the study engine resolves",
    }


# ── Units ────────────────────────────────────────────────────────────────────


def _cell_names(settings: Mapping[str, Any]) -> list[str]:
    if settings["profile"] == "parent":
        return [cell[0] for cell in settings["scan_validation"]["span_split"]["comparison"]["cells"]]
    if settings["profile"] == "stage_span":
        return [cell[0] for cell in settings["units"]["cells"]]
    return [cell[0] for cell in settings["cells"]]


def profile_units(settings: Mapping[str, Any]) -> list[lsc.Unit]:
    """Design-major, then cell, then replicate; the controls at the control count."""

    if settings["profile"] == "parent":
        return isc.simulation_units(settings["scan_validation"])
    counts = settings["replicates"]
    control = set(settings["control_cells"])
    return [
        lsc.Unit(design, cell, replicate)
        for design in settings["designs"]
        for cell in _cell_names(settings)
        for replicate in range(int(counts["control" if cell in control else "other"]))
    ]


def signal_cells(settings: Mapping[str, Any]) -> list[tuple[str, str]]:
    control = set(settings["control_cells"])
    return [(design, cell) for design in settings["designs"] for cell in _cell_names(settings) if cell not in control]


def unit_fields(unit: lsc.Unit, settings: Mapping[str, Any]) -> dict[str, Any]:
    """The unit's cell description and both seeds, without generating anything."""

    profile = settings["profile"]
    fields: dict[str, Any] = {"profile": profile, "design": unit.design, "cell": unit.cell,
                              "replicate": int(unit.replicate)}
    if profile == "parent":
        comparison = settings["scan_validation"]["span_split"]["comparison"]
        cell = lsc.cell_spec(comparison, unit.cell)
        base = int(comparison["base_seed"])
        fields.update(
            trajectory_mode=cell["mode"], effect_size=float(cell["effect"]), span_alignment=cell["phi"],
            span_support=cell["support"], span_scale=None,
            generator_seed=lsc.generator_seed(base, unit.design, unit.replicate),
            rrpp_seed=lsc.rrpp_seed(base, unit.design, unit.cell, unit.replicate),
        )
    elif profile == "stage_span":
        cell = stage_span_cell(settings, unit.cell)
        base = int(settings["units"]["base_seed"])
        fields.update(
            trajectory_mode=cell["mode"], effect_size=float(cell["effect"]), span_alignment=cell["phi"],
            span_support=cell["support"], span_scale=cell["span_scale"],
            generator_seed=lsc.generator_seed(base, unit.design, unit.replicate),
            rrpp_seed=lsc.rrpp_seed(base, unit.design, unit.cell, unit.replicate),
        )
    else:
        cell = surgery_cell(settings, unit.cell)
        generator, evaluation = surgery_params(unit, settings)
        fields.update(
            trajectory_mode=cell["mode"], effect_size=float(cell["effect"]), span_alignment=None,
            span_support=None, span_scale=None, generator_seed=int(generator.seed), rrpp_seed=int(evaluation.seed),
        )
    return fields


def unit_inputs(unit: lsc.Unit, settings: Mapping[str, Any], fields: Mapping[str, Any]) -> tuple[Any, Any]:
    """The unit's dataset and its production evaluation parameters."""

    from motco.simulations.reference import load_reference
    from motco.simulations.semisynthetic import generate_semisynthetic_trajectory

    profile = settings["profile"]
    if profile == "parent":
        span_split = settings["scan_validation"]["span_split"]
        comparison = span_split["comparison"]
        n_stages = int(comparison["designs"][unit.design]["generator"]["n_stages"])
        dataset = lsc.unit_dataset(unit, comparison)
        return dataset, ssc.production_params(span_split, seed=int(fields["rrpp_seed"]), n_stages=n_stages)
    if profile == "stage_span":
        dataset = generate_semisynthetic_trajectory(stage_span_generator_params(unit, settings),
                                                    reference=load_reference())
        return dataset, stage_span_params(settings, seed=int(fields["rrpp_seed"]))
    generator, evaluation = surgery_params(unit, settings)
    return generate_semisynthetic_trajectory(generator, reference=load_reference()), evaluation


# ── The arms of one unit ─────────────────────────────────────────────────────


def omnibus_p(p_values: Mapping[str, Any]) -> float | None:
    """The smallest available trajectory p-value, compared with alpha / 3 (Bonferroni)."""

    finite = [value for statistic in STATISTICS if (value := _finite(p_values.get(statistic))) is not None]
    return min(finite) if finite else None


def minp_p_value(observed: Mapping[str, Any], nulls: Mapping[str, Sequence[float]]) -> float | None:
    """The smallest trajectory p-value calibrated on the joint permutation draws (Westfall–Young min-P).

    Each statistic's observed value and its B draws form one set of B + 1 values; each value's p-value is the
    share of the set at or above it (the observed one's is the evaluator's empirical p-value, a non-finite draw
    counting as below every value). The result is the share of the B + 1 rows whose smallest p-value is at or
    below the observed row's, so it is never below 1 / (B + 1). Statistics with a non-finite observed value are
    left out; ``None`` when none remains.
    """

    names = [s for s in STATISTICS if s in nulls and _finite(observed.get(s)) is not None]
    if not names:
        return None
    lengths = {len(nulls[s]) for s in names}
    if len(lengths) != 1:
        raise ValueError(f"The null draws differ in length across statistics: {sorted(lengths)}.")
    rows = []
    for statistic in names:
        values = np.concatenate([[float(observed[statistic])], np.asarray(nulls[statistic], dtype=float)])
        values = np.where(np.isfinite(values), values, -np.inf)
        at_or_above = values.size - np.searchsorted(np.sort(values), values, side="left")
        rows.append(at_or_above / values.size)
    smallest = np.min(np.vstack(rows), axis=0)
    return float(np.count_nonzero(smallest <= smallest[0]) / smallest.size)


def bh_adjusted(
    f_values: ArrayLike, degenerate: ArrayLike, df_interaction: int, df_residual: int
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Upper-tail F p-values and Benjamini–Hochberg adjusted p-values over the non-degenerate features.

    Returns ``(tested, p, adjusted)``; ``p`` and ``adjusted`` are ``NaN`` outside the tested family.
    """

    from scipy.stats import f as f_distribution
    from scipy.stats import false_discovery_control

    f = np.asarray(f_values, dtype=float)
    tested = ~np.asarray(degenerate, dtype=bool) & np.isfinite(f)
    p = np.full(f.shape, np.nan)
    adjusted = np.full(f.shape, np.nan)
    if tested.any():
        p[tested] = f_distribution.sf(f[tested], int(df_interaction), int(df_residual))
        adjusted[tested] = false_discovery_control(p[tested], method="bh")
    return tested, p, adjusted


def fdr_reading(
    f_values: ArrayLike,
    degenerate: ArrayLike,
    df_interaction: int,
    df_residual: int,
    *,
    null_mask: np.ndarray | None = None,
) -> dict[str, Any]:
    """The ``fdr`` arm: the smallest adjusted p-value, the discoveries and those on population-null features."""

    tested, p, adjusted = bh_adjusted(f_values, degenerate, df_interaction, df_residual)
    discovered = tested & (np.nan_to_num(adjusted, nan=np.inf) <= ALPHA)
    return {
        "n_tested": int(tested.sum()),
        "n_degenerate": int((~tested).sum()),
        "min_p": float(np.nanmin(p)) if tested.any() else None,
        "min_adjusted": float(np.nanmin(adjusted)) if tested.any() else None,
        "discoveries": int(discovered.sum()),
        "discoveries_population_null": None if null_mask is None else int((discovered & null_mask).sum()),
    }


def recorded_space_params(params: Any, space: str, *, rank: int) -> Any:
    """The unit's production parameters with only the integration method changed, on the same blocks."""

    layers = (params.integration_params or {}).get("layers")
    blocks = {} if layers is None else {"layers": list(layers)}
    if space == "concat":
        return dataclasses.replace(params, integration_method="concat", integration_params=blocks)
    if space == "pca_matched":
        return dataclasses.replace(params, integration_method="pca",
                                   integration_params={"n_components": int(rank), **blocks})
    raise ValueError(f"Unknown recorded space {space!r}.")


def _capture(result: Any) -> dict[str, Any]:
    capture = result.realized_geometry.get("span_capture") or {}
    return {key: capture.get(key) for key in ("group_difference", "pooled_stage", "random_share")}


def arm_readings(
    production: Mapping[str, Any],
    split: Mapping[str, Any],
    scan: Mapping[str, Any],
    fdr: Mapping[str, Any],
    spaces: Mapping[str, Mapping[str, Any]],
) -> dict[str, dict[str, Any]]:
    """Every arm's p-value and its rejection under the arm's rule."""

    p_values: dict[str, Any] = {
        "production": production.get("p_min"),
        **{s: (production.get("p_values") or {}).get(s) for s in STATISTICS},
        "production_minp": production.get("p_minp"),
        "scan": scan.get("p_value_global"),
        "fdr": fdr.get("min_adjusted"),
        **{part: (split.get(part) or {}).get("p_value_f") for part in ssc.PARTS},
    }
    for space, block in spaces.items():
        p_values[space] = block.get("p_min")
        for s in STATISTICS:
            p_values[f"{space}_{s}"] = (block.get("p_values") or {}).get(s)
    return {arm: {"p_value": _finite(p_values.get(arm)), "rejects": rejects(p_values.get(arm), RULES[arm])}
            for arm in ARMS}


def evaluate_unit(unit: lsc.Unit, settings: Mapping[str, Any]) -> dict[str, Any]:
    """One dataset, one production fit and evaluation, then every arm at the unit's RRPP seed."""

    from motco.simulations.evaluation import evaluate_semisynthetic_trajectory, integrate_semisynthetic_dataset
    from motco.stats.interaction import interaction_feature_scan, interaction_split_test

    start = perf_counter()
    timing: dict[str, float] = {}
    base: dict[str, Any] = {}
    try:
        base = {**unit_fields(unit, settings), "settings_hash": lsc.settings_hash(settings)}
        seed = int(base["rrpp_seed"])
        permutations = int(settings["permutations"])
        dataset, params = unit_inputs(unit, settings, base)
        mark = perf_counter()
        latent = integrate_semisynthetic_dataset(dataset, params)
        timing["fit"] = perf_counter() - mark
        mark = perf_counter()
        result = evaluate_semisynthetic_trajectory(dataset, params, latent=latent)
        timing["evaluation"] = perf_counter() - mark
        metadata = result.latent_matrix_metadata
        rank = int(metadata["selected_lv"])
        production = {
            "rank": rank,
            "cv_mean_auroc": metadata.get("cv_mean_auroc"),
            "pair_statistics": lsc._finite(result.pair_statistics),
            "p_values": lsc._finite(result.p_values),
            "p_min": omnibus_p(result.p_values),
            "p_minp": minp_p_value(result.pair_statistics, result.null_distributions or {}),
            "span_capture": _capture(result),
        }
        features, basis, full, reduced = ssc.split_inputs(dataset, latent, params)
        if latent.artifacts is None:
            raise RuntimeError("The production fit carries no artifacts; the baseline needs its features.")
        columns = [str(column) for column in latent.artifacts.features.columns]
        blocks = list(dict.fromkeys(column.split("__", 1)[0] for column in columns))
        mark = perf_counter()
        split_result = interaction_split_test(features, full, reduced, basis, permutations=permutations, seed=seed)
        timing["split"] = perf_counter() - mark
        split: dict[str, Any] = {
            name: {key: getattr(split_result.parts[name], key) for key in ("dimension", "f", "p_value_f", "p_value_ss")}
            for name in ssc.PARTS
        }
        split["span_share"] = split_result.span_share
        mark = perf_counter()
        scan_result = interaction_feature_scan(features, full, reduced, permutations=permutations, seed=seed)
        timing["scan"] = perf_counter() - mark
        null_mask = isc.population_null_mask(
            dataset, columns, params.group_col, float(settings["population_null_tolerance"])
        )
        hits = scan_result.hits(ALPHA) if scan_result.p_value_adjusted is not None else np.array([], dtype=int)
        scan = {
            "n_features": scan_result.n_features,
            "n_degenerate": scan_result.n_degenerate,
            "df_interaction": scan_result.df_interaction,
            "df_residual": scan_result.df_residual,
            "max_f": scan_result.max_f,
            "p_value_global": scan_result.p_value_global,
            "hits": int(len(hits)),
            "hits_population_null": None if null_mask is None else int(null_mask[hits].sum()),
        }
        fdr = fdr_reading(scan_result.f, scan_result.degenerate, scan_result.df_interaction,
                          scan_result.df_residual, null_mask=null_mask)
        spaces: dict[str, dict[str, Any]] = {}
        for space in settings["recorded_spaces"]:
            mark = perf_counter()
            space_result = evaluate_semisynthetic_trajectory(
                dataset, recorded_space_params(params, space, rank=rank)
            )
            timing[space] = perf_counter() - mark
            spaces[space] = {
                "dimension": lsc._retained(space_result.latent_matrix_metadata),
                "p_values": lsc._finite(space_result.p_values),
                "p_min": omnibus_p(space_result.p_values),
            }
    except Exception as exc:  # recorded, never aborts the shard
        timing["unit"] = perf_counter() - start
        fallback = {"profile": settings["profile"], "design": unit.design, "cell": unit.cell,
                    "replicate": int(unit.replicate)}
        return {**fallback, **base, "status": "failed", "error_type": type(exc).__name__,
                "error_message": str(exc), "runtime_seconds": timing}
    timing["unit"] = perf_counter() - start
    return lsc._jsonable({
        **base,
        "status": "completed",
        "blocks": blocks,
        "n_features": len(columns),
        "n_population_null": None if null_mask is None else int(null_mask.sum()),
        "production": production,
        "split": split,
        "scan": scan,
        "fdr": fdr,
        "spaces": spaces,
        "arms": arm_readings(production, split, scan, fdr, spaces),
        "runtime_seconds": timing,
    })


# ── Shards, resume and the run ───────────────────────────────────────────────


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
    check_parent_files(settings)
    directory = profile_dir(out_dir, profile)
    lsc.ensure_settings(directory, settings)
    path = lsc.shard_path(directory, index, count)
    path.parent.mkdir(parents=True, exist_ok=True)
    done = ssc.read_records([path], record_key) if path.exists() else {}
    # Rewritten with complete lines only, so an interrupted unit is re-run, not duplicated.
    with path.open("w", encoding="utf-8") as handle:
        for key in sorted(done):
            handle.write(lsc._line(done[key]))
    units = profile_units(settings)
    pending = [unit for position, unit in enumerate(units) if position % count == index and unit_key(unit) not in done]
    if limit is not None:
        pending = pending[:limit]

    def append(position: int, record: Mapping[str, Any]) -> None:
        with path.open("a", encoding="utf-8") as handle:
            handle.write(lsc._line(record))
            handle.flush()
        print(f"[{profile} {index}/{count}] {position}/{len(pending)} {record_key(record)} {record['status']}",
              file=sys.stderr, flush=True)

    if workers <= 1:
        for position, unit in enumerate(pending, start=1):
            append(position, evaluate_unit(unit, settings))
    else:
        with ProcessPoolExecutor(max_workers=workers) as pool:
            futures = [pool.submit(evaluate_unit, unit, settings) for unit in pending]
            for position, future in enumerate(as_completed(futures), start=1):
                append(position, future.result())
    return len(pending)


# ── Merge ────────────────────────────────────────────────────────────────────


CELL_FIELDS: tuple[str, ...] = ("design", "cell", "trajectory_mode", "effect_size", "span_alignment", "span_support")
RATE_COLUMNS: tuple[str, ...] = (*CELL_FIELDS, "arm", "rule", "rate", "se", "count", "rejections", "failures")
PAIRED_COLUMNS: tuple[str, ...] = (
    *CELL_FIELDS, "arm", "n_paired", "rate_production", "rate_arm", "difference", "b", "c",
    "p_gain", "p_loss", "p_two_sided",
)
CELL_COLUMNS: tuple[str, ...] = (
    *CELL_FIELDS, "completed", "failures", "rank_q10", "rank_median", "rank_q90", "auroc_median",
    "capture_group_difference_median", "capture_pooled_stage_median", "pca_matched_rank_equal", "concat_dimension",
    "n_features", "scan_hits_median", "fdr_discoveries_median", "n_population_null_median",
    "scan_fdp_mean", "scan_fdp_units", "fdr_fdp_mean", "fdr_fdp_units",
)


def _cell_groups(records: Sequence[Mapping[str, Any]], settings: Mapping[str, Any]) -> list[tuple[dict, list]]:
    """(cell fields, records) per design and cell, in the settings' order; cells without records included."""

    grouped = ssc._grouped(records, lambda r: (str(r["design"]), str(r["cell"])))
    out = []
    for design in settings["designs"]:
        for cell in _cell_names(settings):
            fields = {key: value for key, value in unit_fields(lsc.Unit(design, cell, 0), settings).items()
                      if key in CELL_FIELDS}
            out.append((fields, grouped.get((design, cell), [])))
    return out


def _arm_rejection(record: Mapping[str, Any], arm: str) -> bool | None:
    return ((record.get("arms") or {}).get(arm) or {}).get("rejects")


def rate_rows(records: Sequence[Mapping[str, Any]], settings: Mapping[str, Any]) -> list[dict[str, Any]]:
    """Per design, cell and arm: the rule, rate, SE, available and rejected counts; failed units left out."""

    rows = []
    for fields, group in _cell_groups(records, settings):
        completed = [r for r in group if r.get("status") == "completed"]
        for arm in ARMS:
            rows.append({
                **fields, "arm": arm, "rule": RULES[arm],
                **ssc._rate_row([_arm_rejection(r, arm) for r in completed]),
                "failures": len(group) - len(completed),
            })
    return rows


def paired_row(pairs: Sequence[tuple[bool | None, bool | None]]) -> dict[str, Any]:
    """``production`` against one arm over the units where both are available.

    ``b`` counts units where only the arm rejects (production's loss), ``c`` where only production rejects (its
    gain); the exact one-sided p-values are ``production_latent_space_check.paired_p_values``.
    """

    usable = [(bool(production), bool(arm)) for production, arm in pairs if production is not None and arm is not None]
    n = len(usable)
    b = sum(1 for production, arm in usable if arm and not production)
    c = sum(1 for production, arm in usable if production and not arm)
    return {
        "n_paired": n,
        "rate_production": sum(production for production, _ in usable) / n if n else None,
        "rate_arm": sum(arm for _, arm in usable) / n if n else None,
        "difference": (c - b) / n if n else None,
        "b": b,
        "c": c,
        **plsc.paired_p_values(b, c),
    }


def paired_rows(records: Sequence[Mapping[str, Any]], settings: Mapping[str, Any]) -> list[dict[str, Any]]:
    rows = []
    for fields, group in _cell_groups(records, settings):
        completed = [r for r in group if r.get("status") == "completed"]
        for arm in BASELINE_ARMS:
            pairs = [(_arm_rejection(r, "production"), _arm_rejection(r, arm)) for r in completed]
            row = paired_row(pairs)
            rows.append({**fields, "arm": arm, **{key: row[key] for key in PAIRED_COLUMNS[len(CELL_FIELDS) + 1:]}})
    return rows


def _fdp(discoveries: Any, on_null: Any) -> float | None:
    """The realized false discovery proportion V / max(R, 1); ``None`` without the population-null set."""

    if on_null is None or discoveries is None:
        return None
    return int(on_null) / max(int(discoveries), 1)


def cell_rows(records: Sequence[Mapping[str, Any]], settings: Mapping[str, Any]) -> list[dict[str, Any]]:
    """Per design and cell: the fit, the capture shares, discoveries and realized false discovery proportions."""

    rows = []
    for fields, group in _cell_groups(records, settings):
        completed = [r for r in group if r.get("status") == "completed"]
        ranks = [r["production"]["rank"] for r in completed]
        scan_fdp = [_fdp(r["scan"].get("hits"), r["scan"].get("hits_population_null")) for r in completed]
        fdr_fdp = [_fdp(r["fdr"].get("discoveries"), r["fdr"].get("discoveries_population_null")) for r in completed]
        matched = [int(r["spaces"]["pca_matched"]["dimension"]) == int(r["production"]["rank"])
                   for r in completed if "pca_matched" in r["spaces"]]
        rows.append({
            **fields,
            "completed": len(completed),
            "failures": len(group) - len(completed),
            "rank_q10": _quantile(ranks, 0.1),
            "rank_median": _quantile(ranks, 0.5),
            "rank_q90": _quantile(ranks, 0.9),
            "auroc_median": _median(r["production"].get("cv_mean_auroc") for r in completed),
            "capture_group_difference_median": _median(
                (r["production"].get("span_capture") or {}).get("group_difference") for r in completed),
            "capture_pooled_stage_median": _median(
                (r["production"].get("span_capture") or {}).get("pooled_stage") for r in completed),
            "pca_matched_rank_equal": (sum(matched) / len(matched)) if matched else None,
            "concat_dimension": _median((r["spaces"].get("concat") or {}).get("dimension") for r in completed),
            "n_features": _median(r.get("n_features") for r in completed),
            "scan_hits_median": _median(r["scan"].get("hits") for r in completed),
            "fdr_discoveries_median": _median(r["fdr"].get("discoveries") for r in completed),
            "n_population_null_median": _median(r.get("n_population_null") for r in completed),
            "scan_fdp_mean": _mean(scan_fdp),
            "scan_fdp_units": sum(1 for value in scan_fdp if value is not None),
            "fdr_fdp_mean": _mean(fdr_fdp),
            "fdr_fdp_units": sum(1 for value in fdr_fdp if value is not None),
        })
    return rows


def accounting(
    records: Sequence[Mapping[str, Any]], settings: Mapping[str, Any], *, expected: int, duplicates: int
) -> dict[str, Any]:
    failures: dict[str, int] = {}
    errors: dict[str, int] = {}
    unit_seconds: list[float] = []
    stages: dict[str, list[float]] = {}
    for record in records:
        if record.get("status") != "completed":
            label = f"{record['design']}/{record['cell']}"
            failures[label] = failures.get(label, 0) + 1
            error = str(record.get("error_type"))
            errors[error] = errors.get(error, 0) + 1
        runtime = record.get("runtime_seconds") or {}
        unit_seconds.append(float(runtime.get("unit") or 0.0))
        for stage, value in runtime.items():
            if stage != "unit" and value is not None:
                stages.setdefault(stage, []).append(float(value))
    return {
        "profile": settings["profile"],
        "settings_hash": lsc.settings_hash(settings),
        "expected_units": expected,
        "complete_units": len(records),
        "missing_units": expected - len(records),
        "failed_units": sum(failures.values()),
        "failures_by_cell": failures,
        "failures_by_error": errors,
        "duplicate_records": duplicates,
        "recorded_core_hours": sum(unit_seconds) / 3600.0,
        "unit_seconds": {
            "median": _median(unit_seconds),
            "q90": _quantile(unit_seconds, 0.9),
            "max": max(unit_seconds) if unit_seconds else None,
        },
        "stage_seconds_median": {stage: _median(values) for stage, values in sorted(stages.items())},
    }


def _differing(key: Sequence[Any], test: str, wanted: tuple[int, int], got: tuple[int, int] | None) -> dict:
    return {
        "cell": list(key),
        "test": test,
        "committed": {"count": wanted[0], "rejections": wanted[1]},
        "reproduced": None if got is None else {"count": got[0], "rejections": got[1]},
    }


def _counts(row: Mapping[str, Any] | None) -> tuple[int, int] | None:
    return None if row is None else (int(row["count"]), int(row["rejections"]))


def parent_reproduction(
    rates: Sequence[Mapping[str, Any]], settings: Mapping[str, Any], committed: Sequence[Mapping[str, Any]]
) -> dict[str, Any]:
    """``scan``, the split parts and the trajectory statistics against the committed comparator table, exactly."""

    mine = {(r["design"], r["cell"], r["arm"]): r for r in rates}
    designs = set(settings["designs"])
    expected = [row for row in committed if row["design"] in designs and row["test"] in PARENT_TESTS]
    differing = []
    for row in expected:
        wanted = (int(row["count"]), int(row["rejections"]))
        got = _counts(mine.get((row["design"], row["cell"], PARENT_TESTS[row["test"]])))
        if got != wanted:
            differing.append(_differing((row["design"], row["cell"]), row["test"], wanted, got))
    return {
        "profile": "parent",
        "passed": bool(expected) and not differing,
        "rows_compared": len(expected),
        "differing": differing,
        "source": PARENT_SOURCE,
        "files": dict(settings["parent_files"]),
    }


def surgery_reproduction(
    rates: Sequence[Mapping[str, Any]], settings: Mapping[str, Any], operating: Mapping[str, Sequence[Mapping]]
) -> dict[str, Any]:
    """``pls_cv``'s ``delta``, ``angle`` and ``shape`` against the slice's operating tables, exactly.

    ``operating`` maps each slice profile to its committed rows; a three-block row has an empty layers column.
    """

    mine = {(r["design"], r["cell"], r["arm"]): r for r in rates}
    indexed: dict[tuple, Mapping[str, Any]] = {}
    for source, rows in operating.items():
        for row in rows:
            if row["evaluation.integration_params.layers"] != "":
                continue
            indexed[(source, int(row["generator.n_samples"]), str(row["trajectory_mode"]),
                     float(row["effect_size"]), str(row["statistic"]))] = row
    compared = 0
    differing = []
    missing = []
    for column in settings["designs"]:
        n_samples = int(settings["n_samples"][column])
        for cell, mode, effect, source in settings["cells"]:
            for statistic in STATISTICS:
                row = indexed.get((source, n_samples, mode, float(effect), statistic))
                if row is None:
                    missing.append({"cell": [column, cell], "test": statistic, "source": source})
                    continue
                compared += 1
                wanted = (int(row["n_available"]), int(row["n_rejected"]))
                got = _counts(mine.get((column, cell, statistic)))
                if got != wanted:
                    differing.append(_differing((column, cell), statistic, wanted, got))
    return {
        "profile": "surgery",
        "passed": compared > 0 and not differing and not missing,
        "rows_compared": compared,
        "differing": differing,
        "unresolved_rows": missing,
        "source": SURGERY_SOURCE,
        "files": dict(settings["parent_files"]),
    }


def reproduction(rates: Sequence[Mapping[str, Any]], settings: Mapping[str, Any]) -> dict[str, Any] | None:
    if settings["profile"] == "parent":
        committed = lsc._read_csv(_parent_path(settings, "comparator_table.csv"))
        return parent_reproduction(rates, settings, committed)
    if settings["profile"] == "surgery":
        operating = {
            source: lsc._read_csv(ssc._resolve(paths["operating"])) for source, paths in settings["sources"].items()
        }
        return surgery_reproduction(rates, settings, operating)
    return None


def _parent_path(settings: Mapping[str, Any], name: str) -> Path:
    for stored in settings["parent_files"]:
        if Path(stored).name == name:
            return ssc._resolve(stored)
    raise KeyError(f"No bound parent file named {name!r}.")


def merge(out_dir: Path, profile: str) -> dict[str, Any]:
    """Merge one profile's shards and write its tables, accounting and reproduction check."""

    directory = profile_dir(out_dir, profile)
    settings = lsc.load_settings(directory)
    check_parent_files(settings)
    expected_keys = [unit_key(unit) for unit in profile_units(settings)]
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
    summary = accounting(records, settings, expected=len(expected_keys),
                         duplicates=sum(count - 1 for count in seen.values()))
    rates = rate_rows(records, settings)
    lsc._write_csv(directory / "rate_table.csv", rates, RATE_COLUMNS)
    lsc._write_csv(directory / "paired_table.csv", paired_rows(records, settings), PAIRED_COLUMNS)
    lsc._write_csv(directory / "cell_table.csv", cell_rows(records, settings), CELL_COLUMNS)
    lsc._write_json(directory / "accounting.json", summary)
    result = reproduction(rates, settings)
    if result is not None:
        lsc._write_json(directory / "reproduction.json", result)
    return {
        "profile": profile,
        **{key: summary[key] for key in ("complete_units", "missing_units", "failed_units", "recorded_core_hours")},
        "reproduction_passed": None if result is None else result["passed"],
        "reproduction_differing": None if result is None else len(result["differing"]),
    }


# ── Decision (pure functions) ────────────────────────────────────────────────


def _index(rows: Sequence[Mapping[str, Any]], *keys: str) -> dict[tuple, Mapping[str, Any]]:
    return {tuple(str(row[key]) for key in keys): row for row in rows}


def _rate_of(row: Mapping[str, Any] | None) -> dict[str, Any] | None:
    if row is None or int(row.get("count") or 0) <= 0:
        return None
    return {"rate": int(row["rejections"]) / int(row["count"]), "count": int(row["count"]),
            "rejections": int(row["rejections"])}


def _signal_cells(block: Mapping[str, Any]) -> list[tuple[str, str]]:
    return signal_cells(block["settings"])


def _cell_meta(block: Mapping[str, Any], design: str, cell: str) -> dict[str, Any]:
    fields = unit_fields(lsc.Unit(design, cell, 0), block["settings"])
    return {key: fields[key] for key in CELL_FIELDS if key not in ("design", "cell")}


def size_reading(blocks: Mapping[str, Mapping[str, Any]]) -> dict[str, Any]:
    """Each judged arm's controls in the two families under ``judge_rates``; admissible when both hold."""

    out: dict[str, Any] = {}
    for arm in JUDGED_ARMS:
        families: dict[str, Any] = {}
        for profile in DECIDED_PROFILES:
            block = blocks[profile]
            rates = _index(block["rates"], "design", "cell", "arm")
            families[profile] = grade.judge_rates([
                isc._size_observation(rates.get((design, cell, arm)), profile=profile, design=design, cell=cell,
                                      arm=arm)
                for design in block["settings"]["designs"]
                for cell in CONTROL_CELLS
            ])
        out[arm] = {"admissible": all(family["holds"] for family in families.values()), "families": families}
    return out


def label_cells(
    blocks: Mapping[str, Mapping[str, Any]], admissible: Sequence[str], *, m: int
) -> list[dict[str, Any]]:
    """Each signal cell's label from the paired exact tests against every admissible baseline arm."""

    a = len(admissible)
    gain_level = LABEL_ALPHA / m if m else None
    loss_level = LABEL_ALPHA / (m * a) if m and a else None
    labels = []
    for profile in DECIDED_PROFILES:
        block = blocks[profile]
        paired = _index(block["paired"], "design", "cell", "arm")
        for design, cell in _signal_cells(block):
            entry: dict[str, Any] = {"profile": profile, "design": design, "cell": cell,
                                     **_cell_meta(block, design, cell)}
            against: dict[str, Any] = {}
            unavailable = []
            for arm in admissible:
                row = paired.get((design, cell, arm))
                n = 0 if row is None else int(row["n_paired"])
                if n == 0:
                    unavailable.append(arm)
                    continue
                b, c = int(row["b"]), int(row["c"])
                p = plsc.paired_p_values(b, c)
                difference = (c - b) / n
                against[arm] = {
                    "n_paired": n, "b": b, "c": c, "difference": difference,
                    "rate_production": float(row["rate_production"]), "rate_arm": float(row["rate_arm"]),
                    "p_gain": p["p_gain"], "p_loss": p["p_loss"],
                    "gain": bool(p["p_gain"] <= gain_level and difference >= MARGIN - _MARGIN_TOLERANCE),
                    "loss": bool(p["p_loss"] <= loss_level and -difference >= MARGIN - _MARGIN_TOLERANCE),
                }
            leading = [arm for arm, reading in against.items() if reading["loss"]]
            if not admissible:
                label, reason = "unresolved", "no admissible baseline arm"
            elif leading:
                label, reason = "baseline_leads", None
            elif not unavailable and all(reading["gain"] for reading in against.values()):
                label, reason = "production_leads", None
            else:
                label, reason = "unresolved", ("unavailable" if unavailable else None)
            entry.update({"label": label, "against": against})
            if leading:
                entry["leading_arms"] = leading
            if unavailable:
                entry["unavailable"] = unavailable
            if reason is not None:
                entry["reason"] = reason
            labels.append(entry)
    return labels


def plant_family(profile: str, mode: str, alignment: Any) -> str:
    if profile == "stage_span":
        out_of_span = alignment is not None and float(alignment) > 0
        return "stage_span_out_of_span" if out_of_span else f"stage_span_{mode}"
    return f"surgery_{mode}"


def characterization_grade(stage_span: Mapping[str, Any]) -> dict[str, Any]:
    """Off-target statistics of the in-span stage-span cells ``production`` detects in at least 0.80 of units."""

    rates = _index(stage_span["rates"], "design", "cell", "arm")
    detected = []
    considered = 0
    for design, cell in _signal_cells(stage_span):
        meta = _cell_meta(stage_span, design, cell)
        mode, alignment = meta["trajectory_mode"], meta["span_alignment"]
        if mode not in TARGETS or (alignment is not None and float(alignment) > 0):
            continue
        considered += 1
        production = _rate_of(rates.get((design, cell, "production")))
        if production is None or production["rejections"] < DETECTION_FLOOR * production["count"] - 1e-9:
            continue
        detected.append({"design": design, "cell": cell, "trajectory_mode": mode, "production": production})
    observations = [
        isc._size_observation(rates.get((entry["design"], entry["cell"], statistic)), design=entry["design"],
                              cell=entry["cell"], statistic=statistic)
        for entry in detected
        for statistic in OFF_TARGET[entry["trajectory_mode"]]
    ]
    if not detected:
        return {"grade": "no_detected_cell", "considered_cells": considered, "detected_cells": [], "family": None}
    judged = grade.judge_rates(observations)
    result = {
        "grade": "names_kind" if judged["holds"] else "detects_only",
        "considered_cells": considered,
        "detection_floor": DETECTION_FLOOR,
        "detected_cells": detected,
        "family": judged,
    }
    if not judged["holds"]:
        result["failing"] = ssc._failing_observations(judged)
    return result


def off_target_reported(blocks: Mapping[str, Mapping[str, Any]]) -> list[dict[str, Any]]:
    """Off-target rates of the cells the grade leaves out: every surgery cell and the out-of-span cells."""

    rows = []
    for profile in DECIDED_PROFILES:
        block = blocks[profile]
        rates = _index(block["rates"], "design", "cell", "arm")
        for design, cell in _signal_cells(block):
            meta = _cell_meta(block, design, cell)
            mode, alignment = meta["trajectory_mode"], meta["span_alignment"]
            if mode not in TARGETS:
                continue
            if profile == "stage_span" and (alignment is None or float(alignment) == 0):
                continue
            rows.append({
                "profile": profile, "design": design, "cell": cell, "trajectory_mode": mode,
                "production": _rate_of(rates.get((design, cell, "production"))),
                **{statistic: _rate_of(rates.get((design, cell, statistic))) for statistic in OFF_TARGET[mode]},
            })
    return rows


def descriptive_tables(
    blocks: Mapping[str, Mapping[str, Any]], labels: Sequence[Mapping[str, Any]], parent: Mapping[str, Any] | None
) -> dict[str, Any]:
    """The label map, every arm's rate per cell, discoveries, the fit and the parent's comparison."""

    label_of = {(entry["profile"], entry["design"], entry["cell"]): entry["label"] for entry in labels}
    tally: dict[tuple[str, str], dict[str, int]] = {}
    for entry in labels:
        family = plant_family(entry["profile"], entry["trajectory_mode"], entry["span_alignment"])
        counts = tally.setdefault((family, entry["design"]), {label: 0 for label in LABELS})
        counts[entry["label"]] += 1
    label_map = [{"family": family, "design": design, **counts} for (family, design), counts in tally.items()]
    cell_rates, discoveries, fit = [], [], []
    for profile in DECIDED_PROFILES:
        block = blocks[profile]
        rates = _index(block["rates"], "design", "cell", "arm")
        cells = _index(block["cells"], "design", "cell")
        source = {"source": f"this run, {profile} profile"}
        for design in block["settings"]["designs"]:
            for cell in _cell_names(block["settings"]):
                meta = _cell_meta(block, design, cell)
                cell_rates.append({
                    "profile": profile, "design": design, "cell": cell, **meta,
                    "label": label_of.get((profile, design, cell)),
                    "rates": {arm: _rate_of(rates.get((design, cell, arm))) for arm in ARMS}, **source,
                })
                row = cells.get((design, cell)) or {}
                discoveries.append({
                    "profile": profile, "design": design, "cell": cell,
                    **{key: _finite(row.get(key)) for key in ("scan_hits_median", "fdr_discoveries_median",
                                                              "n_population_null_median", "scan_fdp_mean",
                                                              "fdr_fdp_mean")},
                    **source,
                })
                fit.append({
                    "profile": profile, "design": design, "cell": cell,
                    **{key: _finite(row.get(key)) for key in ("rank_q10", "rank_median", "rank_q90", "auroc_median",
                                                              "capture_group_difference_median",
                                                              "capture_pooled_stage_median",
                                                              "pca_matched_rank_equal")},
                    **source,
                })
    parent_rows = []
    if parent is not None:
        rates = _index(parent["rates"], "design", "cell", "arm")
        for design in parent["settings"]["designs"]:
            for cell in _cell_names(parent["settings"]):
                parent_rows.append({
                    "design": design, "cell": cell,
                    **{arm: _rate_of(rates.get((design, cell, arm)))
                       for arm in ("production", "fdr", *PARENT_TESTS.values())},
                    "source": "this run, parent profile",
                })
    return {"map": label_map, "cell_rates": cell_rates, "discoveries": discoveries, "fit": fit,
            "parent": parent_rows}


def decision(
    parent: Mapping[str, Any], stage_span: Mapping[str, Any], surgery: Mapping[str, Any]
) -> dict[str, Any]:
    """Apply rule version 1 to the three merged profiles.

    Each block holds ``settings``, ``accounting``, ``rates``, ``paired`` and ``cells`` (rows); ``parent`` and
    ``surgery`` also ``reproduction``. Every criterion is recorded whatever the verdict.
    """

    blocks = {"parent": parent, "stage_span": stage_span, "surgery": surgery}
    decided = {name: blocks[name] for name in DECIDED_PROFILES}
    m = sum(len(_signal_cells(block)) for block in decided.values())
    size = size_reading(decided)
    admissible = [arm for arm in BASELINE_ARMS if size[arm]["admissible"]]
    labels = label_cells(decided, admissible, m=m)
    characterization = characterization_grade(stage_span)
    header = {
        "rule_version": RULE_VERSION,
        "alpha": ALPHA,
        "omnibus_level": ALPHA / OMNIBUS_DIVISOR,
        "margin": MARGIN,
        "m_signal_cells": m,
        "m_signal_cells_expected": SIGNAL_CELLS_EXPECTED,
        "admissible_baseline_arms": admissible,
        "levels": {
            "production_leads": LABEL_ALPHA / m if m else None,
            "baseline_leads": LABEL_ALPHA / (m * len(admissible)) if m and admissible else None,
        },
        "settings_hashes": {name: lsc.settings_hash(block["settings"]) for name, block in blocks.items()},
        "complete": all(block["accounting"]["missing_units"] == 0 for block in blocks.values()),
        "failed_units": {name: block["accounting"]["failed_units"] for name, block in blocks.items()},
        "reproduction": {
            name: {key: (blocks[name].get("reproduction") or {}).get(key) for key in ("passed", "rows_compared")}
            for name in REPRODUCED_PROFILES
        },
    }
    failing: dict[str, Any] = {}
    failed = [{"profile": name, "differing": (blocks[name].get("reproduction") or {}).get("differing", [])}
              for name in REPRODUCED_PROFILES if not (blocks[name].get("reproduction") or {}).get("passed")]
    leads = [entry for entry in labels if entry["label"] == "production_leads"]
    loses = [entry for entry in labels if entry["label"] == "baseline_leads"]
    if failed:
        verdict = "reproduction_failed"
        failing["reproduction"] = failed
    elif not size["production"]["admissible"]:
        verdict = "production_size_excess"
        failing["production_size"] = {
            family: ssc._failing_observations(judged)
            for family, judged in size["production"]["families"].items() if not judged["holds"]
        }
    elif not admissible:
        verdict = "no_admissible_baseline"
    elif leads and loses:
        verdict = "complementary"
    elif leads:
        verdict = "production_leads"
    elif loses:
        verdict = "baseline_leads"
    else:
        verdict = "unresolved"
    inadmissible = {
        arm: {family: ssc._failing_observations(judged) for family, judged in size[arm]["families"].items()
              if not judged["holds"]}
        for arm in BASELINE_ARMS if not size[arm]["admissible"]
    }
    return {
        **header,
        "verdict": verdict,
        "consequence": CONSEQUENCES[verdict],
        "failing": failing,
        "inadmissible_baseline_arms": inadmissible,
        "criteria": {
            "size": size,
            "labels": labels,
            "label_cells": {
                label: [[entry["profile"], entry["design"], entry["cell"]] for entry in labels
                        if entry["label"] == label]
                for label in LABELS
            },
        },
        "characterization": {**characterization, "off_target_reported": off_target_reported(decided)},
        "descriptive": descriptive_tables(decided, labels, parent),
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
        "rates": _read_csv(directory / "rate_table.csv"),
        "paired": _read_csv(directory / "paired_table.csv"),
        "cells": _read_csv(directory / "cell_table.csv"),
    }
    if profile in REPRODUCED_PROFILES:
        block["reproduction"] = json.loads((directory / "reproduction.json").read_text(encoding="utf-8"))
    return block


def decide(out_dir: Path) -> dict[str, Any]:
    """Read the three merged profiles and write ``univariate_baseline_decision.json``."""

    blocks = {profile: _profile_tables(profile_dir(out_dir, profile), profile) for profile in PROFILES}
    result = decision(blocks["parent"], blocks["stage_span"], blocks["surgery"])
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
    run.add_argument("--reps", type=int, default=None, help="Replicates of every non-control cell "
                                                            "(default: parent 100, stage_span 200, surgery 500)")
    run.add_argument("--control-reps", type=int, default=None, help="Replicates of the anchor and translation "
                                                                    "(default 500)")
    run.add_argument("--recorded-spaces", type=str, default=",".join(DEFAULT_RECORDED_SPACES),
                     help="Comma-separated recorded spaces to evaluate (concat, pca_matched; default pca_matched), "
                          "or 'none'; the arms of a space left out are recorded as unavailable")
    run.add_argument("--designs", type=str, default=None,
                     help="Comma-separated designs (parent, stage_span) or columns (surgery); default: all")
    merge_parser = sub.add_parser("merge", help="Merge one profile's shards")
    merge_parser.add_argument("--profile", choices=PROFILES, required=True)
    merge_parser.add_argument("--out-dir", type=Path, required=True)
    decide_parser = sub.add_parser("decide", help=f"Write {DECISION_NAME} from the three profiles")
    decide_parser.add_argument("--out-dir", type=Path, required=True)
    return parser


def settings_from_args(args: argparse.Namespace) -> dict[str, Any]:
    defaults = DEFAULT_REPLICATES[args.profile]
    replicates = defaults["other"] if args.reps is None else args.reps
    control = defaults["control"] if args.control_reps is None else args.control_reps
    chosen = None if args.designs is None else [name for name in args.designs.split(",") if name]
    spaces = [] if args.recorded_spaces == "none" else [name for name in args.recorded_spaces.split(",") if name]
    common = {"replicates": replicates, "control_replicates": control, "permutations": args.perms,
              "recorded_spaces": spaces}
    if args.profile == "parent":
        return parent_settings(designs=chosen or PARENT_DESIGNS, **common)
    if args.profile == "stage_span":
        return stage_span_settings(designs=chosen or tuple(STAGE_SPAN_DESIGNS), **common)
    return surgery_settings(columns=chosen or SURGERY_COLUMNS, **common)


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    if args.command == "run":
        index, count = lsc.parse_shard(args.shard)
        evaluated = run_shard(args.out_dir, settings_from_args(args), index, count, limit=args.limit,
                              workers=args.workers)
        print(f"evaluated {evaluated} {args.profile} unit(s)", file=sys.stderr)
        return 0
    if args.command == "merge":
        print(json.dumps(merge(args.out_dir, args.profile), indent=1))
        return 0
    result = decide(args.out_dir)
    print(json.dumps({
        "verdict": result["verdict"],
        "admissible_baseline_arms": result["admissible_baseline_arms"],
        "labels": {label: len(cells) for label, cells in result["criteria"]["label_cells"].items()},
        "grade": result["characterization"]["grade"],
        "complete": result["complete"],
    }, indent=1))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
