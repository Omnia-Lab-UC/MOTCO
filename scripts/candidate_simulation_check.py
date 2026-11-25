#!/usr/bin/env python3
"""Simulation-only candidate comparison: the measurement candidates against production, by design family (no R).

The driver behind ``openspec/changes/simulation-only-candidate-comparison/``
(spec ``candidate-simulation-validation``). Three measurement candidates were
kept out of production by real-noise readings on convenience data. On
2026-11-25 the user decided to re-read them in simulation alone, across more
conditions, with those readings kept as recorded limits and not as vetoes. No
real export is read here.

**Arms** (one dataset per unit; every arm at the unit's RRPP seed with serial
draws, net displacement and rank floor 2):

- ``pls_cv`` — production: stage-PLS on the unit's blocks, rank by double CV;
- ``predictive_orthonormal`` — the K − 1 stage-predictive directions of that
  same fit (``pls_coordinates_from_pls`` on the production fit; no second fit
  and no cross-validation); its record carries the fit's rank and AUROC;
- ``weighted`` — its own cross-validated fit with ``block_weighting =
  "equal_variance"`` (``integration_params`` alone changed);
- ``sparse_block`` — its own cross-validated sparse fit, ``block`` scope over
  the default ladder 1.0/0.5/0.2/0.1/0.05, on ``n200``, ``n200_delta0.5`` and
  ``band_n320`` only; elsewhere it is not run.

On the ``adni_design`` profile each arm is one ``analyze_trajectories`` call at
the ADNI gain run's configuration (covariates, the ``df`` correction):
``pls_cv`` cross-validates, ``predictive_orthonormal`` is forced at its rank
(``adni_gain_specificity_check.evaluate_arms``, imported unchanged, its two
operating arms), and ``weighted`` adds ``block_weighting`` to the CV call.

**Profiles** (base seed 2043 unless stated; controls are the anchor and the
translation; 999 permutations):

- ``committed`` — the univariate baseline's ``n80``, ``n200``, ``n500`` and
  ``n200_delta0.5`` and the stage-effect scaling run's ``n320_delta1`` (its
  thirteen stage-span cells) and ``band_n320`` (eleven cells, δ_c = 0.05 from its
  committed calibration), at their committed seeds; controls at 500, other
  cells at 200. 18,800 units.
- ``conditions`` — five new three-block designs with the baseline's thirteen
  cells: ``n200_stages4`` and ``n200_delta0.5_stages4``;
  ``seaad_table`` and ``seaad_table_delta0.5`` (F 11 / 10 / 28, M 9 / 10 / 12
  through ``group_stage_sizes``); ``n200_share0.8``. 16,000 units.
- ``dispersion`` — null only: ``stage_dispersion_scale`` (0.27.0) at four
  patterns on ``seaad_table`` and ``seaad_table_delta0.5`` and the two strong
  ones on ``n200``; anchor and translation at 500. 10,000 units. Recorded arms,
  never decisive: ``pls_cv`` and ``predictive_orthonormal`` under within-stage
  draws on the same latent matrix and seed, and every arm's residual dispersion
  by stage in its measured coordinates.
- ``surgery`` — the design-guidance slice's three-block columns at n = 80 and
  200 (``univariate_baseline_check.surgery_settings``), eleven cells at 500.
  11,000 units.
- ``adni_design`` — the ADNI gain run's simulated design
  (``adni_gain_specificity_check.simulation_settings``, base seed 2040): ten
  cells, controls at 500, plants at 200. 2,600 units.

New design names derive new dataset seeds through
``latent_space_comparison.generator_seed``; every cell of a design and replicate
shares its null dataset.

**Parents** (sha256 in the settings; ``run``, ``merge`` and ``decide`` refuse a
changed file by name): ``settings.json`` and ``rate_table.csv`` of
``results/univariate-baseline-2026-11-23/{stage_span,surgery}/``,
``results/stage-effect-scaling-2026-11-24/{saturated,band}/`` and
``results/adni-gain-specificity-2026-11-17/simulation/``, and
``results/stage-effect-scaling-2026-11-24/calibration/calibration.json``. Each
parent driver's settings builder, called at its committed grade, must rebuild
its committed ``settings.json`` exactly, so the units here are its units.

Subcommands (outputs under ``--out-dir``/``<profile>/``, the decision under
``--out-dir``):

- ``run --profile … --arms …`` — (a shard of) the units into
  ``records/shard_<i>_of_<N>.jsonl``, one line per (unit, arm). ``settings.json``
  holds the settings and their hash; a resume or merge under other settings
  exits naming the keys. ``--arms`` (``all`` or a comma list) selects the arms;
  each selection is recorded in ``settings.json`` outside the hash, and a resume
  skips complete (unit, arm) keys, so a later ``--arms sparse_block`` appends to
  the same records. Use the same shard count when appending.
- ``merge --profile …`` — ``merged.jsonl`` (keep off the repository),
  ``rate_table.csv`` (per design, cell, arm and statistic: rule, status
  ``complete`` | ``partial`` | ``not_run``, rate, SE, counts, failures),
  ``paired_table.csv`` (each candidate against ``pls_cv`` on shared units, and on
  ``dispersion`` the within-stage draws against the unrestricted ones),
  ``cell_table.csv`` (per design, cell and arm: rank, AUROC, captures, keep
  fraction, support, observed and population statistics, residual dispersion),
  ``rank_table.csv`` (counts by ``pls_cv``'s rank), ``accounting.json`` and, for
  ``committed``, ``surgery`` and ``adni_design``, ``reproduction.json``.
- ``decide`` — ``candidate_simulation_decision.json`` (rule version 1): the
  gate ``reproduction_failed``; per candidate and family (``core``,
  ``weak_stage``, ``four_stages``, ``unbalanced``, ``adni_design``) S1 size, S2
  protected power, S3 specificity (the bound at δ = 2.0, the paired excess over
  ``pls_cv`` below it) and G gain, labelled ``size_excess`` | ``power_loss`` |
  ``specificity_excess`` | ``no_gain`` | ``qualifies`` | ``not_run`` (``partial``
  marked); D, size under stage dispersion (with the fallback to a paired excess
  test and ``production_dispersion_excess`` when ``pls_cv`` fails); I, the
  incidental-route labels (never decisive); each arm's verdict
  ``qualifies_everywhere`` | ``scoped_candidate`` | ``keep_production``; the run's
  ``reproduction_failed`` | ``propose_default_change`` | ``scoped_candidates`` |
  ``keep_production``; and descriptive tables.

Records hold seeds and statistics only. With ``--workers`` above 1, pin BLAS to
one thread per process (``OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
MKL_NUM_THREADS=1``).

Examples (the smoke: two replicates per cell, every arm):

    export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
    for P in committed conditions dispersion surgery adni_design; do
      python scripts/candidate_simulation_check.py run --profile $P --arms all --reps 2 --control-reps 2 \\
          --out-dir /tmp/csc --workers 7
      python scripts/candidate_simulation_check.py merge --profile $P --out-dir /tmp/csc
    done
    python scripts/candidate_simulation_check.py decide --out-dir /tmp/csc

    # SLURM: scripts/candidate_simulation_array.sbatch (a profile, an arms selection), then merge and decide
"""

from __future__ import annotations

import argparse
import dataclasses
import json
import math
import os
import sys
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path
from time import perf_counter
from typing import Any, Iterable, Mapping, Sequence

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import adni_anchor_check as aac  # noqa: E402  (the ADNI arms' parameters and records)
import adni_design_check as adc  # noqa: E402  (the ADNI-design unit builder)
import adni_gain_specificity_check as ags  # noqa: E402  (the ADNI-design settings, arms and checks)
import latent_space_comparison as lsc  # noqa: E402  (seeds, settings binding, JSON and CSV helpers)
import predictive_size_excess_check as pse  # noqa: E402  (residual dispersion by stage)
import predictive_subspace_check as psc  # noqa: E402  (the paired count-plus-ceiling rule, the gain rule)
import production_latent_space_check as plsc  # noqa: E402  (arm settings, the exact paired p-values)
import span_split_check as ssc  # noqa: E402  (parent loading, paths, rate rows)
import stage_effect_scaling_check as ses  # noqa: E402  (the scaling run's designs and cells)
import stage_span_paper_grade_decision as grade  # noqa: E402  (the count-plus-ceiling rule)
import univariate_baseline_check as ubc  # noqa: E402  (stage-span and surgery units, cells, labelling rule)

sys.path.pop(0)

DRIVER = "candidate_simulation_check"
DRIVER_VERSION = 1
RULE_VERSION = 1
PROFILES: tuple[str, ...] = ("committed", "conditions", "dispersion", "surgery", "adni_design")
REPRODUCED_PROFILES: tuple[str, ...] = ("committed", "surgery", "adni_design")
#: The run's tiers, each launched on its own approval; a tier not launched leaves its families not_run.
TIERS: dict[int, dict[str, Any]] = {
    1: {"profiles": ["committed", "conditions", "dispersion"],
        "arms": ["pls_cv", "predictive_orthonormal", "weighted"]},
    2: {"profiles": ["surgery", "adni_design"], "arms": ["pls_cv", "predictive_orthonormal", "weighted"]},
    3: {"profiles": ["committed"], "arms": ["sparse_block"]},
}

REFERENCE_ARM = "pls_cv"
PREDICTIVE_ARM = "predictive_orthonormal"
WEIGHTED_ARM = "weighted"
SPARSE_ARM = "sparse_block"
ARMS: tuple[str, ...] = (REFERENCE_ARM, PREDICTIVE_ARM, WEIGHTED_ARM, SPARSE_ARM)
CANDIDATE_ARMS: tuple[str, ...] = (PREDICTIVE_ARM, WEIGHTED_ARM, SPARSE_ARM)
WITHIN_STAGE = "__within_stage"
WITHIN_STAGE_ARMS: tuple[str, ...] = (f"{REFERENCE_ARM}{WITHIN_STAGE}", f"{PREDICTIVE_ARM}{WITHIN_STAGE}")
SPARSE_DESIGNS: tuple[str, ...] = ("n200", "n200_delta0.5", "band_n320")
WEIGHTING = plsc.WEIGHTING
SPARSITY_SCOPE = plsc.SPARSITY_SCOPE
#: (arm, reference arm, comparison) of the paired table; the last three on the dispersion profile only.
PAIRS: tuple[tuple[str, str, str], ...] = (
    (PREDICTIVE_ARM, REFERENCE_ARM, "coordinates"),
    (WEIGHTED_ARM, REFERENCE_ARM, "fit"),
    (SPARSE_ARM, REFERENCE_ARM, "fit"),
)
DISPERSION_PAIRS: tuple[tuple[str, str, str], ...] = (
    (WITHIN_STAGE_ARMS[0], REFERENCE_ARM, "draws"),
    (WITHIN_STAGE_ARMS[1], PREDICTIVE_ARM, "draws"),
    (WITHIN_STAGE_ARMS[1], WITHIN_STAGE_ARMS[0], "coordinates"),
)

STATISTICS: tuple[str, ...] = ("delta", "angle", "shape")
TARGETS: dict[str, str] = {"magnitude": "delta", "orientation": "angle", "shape": "shape"}
ALPHA = 0.05
STRICT, INCLUSIVE = "p < alpha", "p <= alpha"
#: The rule each profile's parent used: the stage-span and surgery parents p < α, the ADNI gain run p ≤ α.
REJECTION_RULES: dict[str, str] = {
    "committed": STRICT, "conditions": STRICT, "dispersion": STRICT, "surgery": STRICT, "adni_design": INCLUSIVE,
}
BASE_SEED = ubc.STAGE_SPAN_BASE_SEED
DEFAULT_PERMUTATIONS = 999
#: The generator's stage effect: designs at it are judged by the bound in S3, designs below by paired excess.
CORE_DELTA = 2.0
CONTROL_CELLS: tuple[str, ...] = ubc.CONTROL_CELLS

DEFAULT_REPLICATES: dict[str, dict[str, int]] = {
    "committed": {"control": 500, "other": 200},
    "conditions": {"control": 500, "other": 200},
    "dispersion": {"control": 500, "other": 500},
    "surgery": {"control": 500, "other": 500},
    "adni_design": {"control": ags.CONTROL_REPLICATES, "other": ags.PLANT_REPLICATES},
}
EXPECTED_UNITS: dict[str, int] = {
    "committed": 18_800, "conditions": 16_000, "dispersion": 10_000, "surgery": 11_000, "adni_design": 2_600,
}

# ── Designs ──────────────────────────────────────────────────────────────────

#: committed design -> (its parent, sample size, stage effect).
COMMITTED_DESIGNS: dict[str, tuple[str, int, float]] = {
    **{name: ("univariate_stage_span", *ubc.STAGE_SPAN_DESIGNS[name]) for name in ("n80", "n200", "n500",
                                                                                    "n200_delta0.5")},
    "n320_delta1": ("scaling_saturated", *ses.SATURATED_DESIGNS["n320_delta1"]),
}
BAND_DESIGN = "band_n320"
SURGERY_COLUMNS: tuple[str, ...] = ("n80_three_blocks", "n200_three_blocks")
SEAAD_TABLE: tuple[tuple[int, ...], tuple[int, ...]] = ((11, 10, 28), (9, 10, 12))
#: condition design -> what differs from the baseline's balanced three-stage n = 200 design at δ = 2.0.
CONDITION_DESIGNS: dict[str, dict[str, Any]] = {
    "n200_stages4": {"n_samples": 200, "delta": 2.0, "n_stages": 4},
    "n200_delta0.5_stages4": {"n_samples": 200, "delta": 0.5, "n_stages": 4},
    "seaad_table": {"table": SEAAD_TABLE, "delta": 2.0},
    "seaad_table_delta0.5": {"table": SEAAD_TABLE, "delta": 0.5},
    "n200_share0.8": {"n_samples": 200, "delta": 2.0, "share": 0.8},
}
#: Stage SD multipliers in stage order (design D4): the within-stage run's measured variance patterns.
DISPERSION_PATTERNS: dict[str, tuple[float, ...]] = {
    "last_narrow": (1.058, 1.130, 0.896),
    "last_wide": (0.897, 0.960, 1.066),
    "last_narrow_strong": (1.120, 1.278, 0.802),
    "last_wide_strong": (0.805, 0.921, 1.137),
}
STRONG_PATTERNS: tuple[str, ...] = ("last_narrow_strong", "last_wide_strong")
#: equal-dispersion design -> the patterns read on it.
DISPERSION_REFERENCES: dict[str, tuple[str, ...]] = {
    "seaad_table": tuple(DISPERSION_PATTERNS),
    "seaad_table_delta0.5": tuple(DISPERSION_PATTERNS),
    "n200": STRONG_PATTERNS,
}
ADNI_DESIGN = ags.simulation_point(ags.SIMULATION_TABLE)

FAMILIES: dict[str, tuple[str, ...]] = {
    "core": ("n80", "n200", "n500"),
    "weak_stage": ("n200_delta0.5", "n320_delta1", BAND_DESIGN),
    "four_stages": ("n200_stages4", "n200_delta0.5_stages4"),
    "unbalanced": ("seaad_table", "seaad_table_delta0.5", "n200_share0.8"),
    "adni_design": (ADNI_DESIGN,),
}
FAMILY_OF: dict[str, str] = {design: family for family, designs in FAMILIES.items() for design in designs}

# ── Parents ──────────────────────────────────────────────────────────────────

REPO = Path(__file__).resolve().parents[1]
RESULTS = REPO / "results"
UNIVARIATE = "univariate-baseline-2026-11-23"
SCALING = "stage-effect-scaling-2026-11-24"
ADNI_GAIN = "adni-gain-specificity-2026-11-17"
#: parent -> (committed directory under results/, bound files, source label).
PARENTS: dict[str, tuple[str, tuple[str, ...], str]] = {
    "univariate_stage_span": (f"{UNIVARIATE}/stage_span", ("settings.json", "rate_table.csv"),
                              "univariate baseline 2026-11-23, stage_span profile"),
    "univariate_surgery": (f"{UNIVARIATE}/surgery", ("settings.json", "rate_table.csv"),
                           "univariate baseline 2026-11-23, surgery profile"),
    "scaling_saturated": (f"{SCALING}/saturated", ("settings.json", "rate_table.csv"),
                          "stage-effect scaling 2026-11-24, saturated profile"),
    "scaling_band": (f"{SCALING}/band", ("settings.json", "rate_table.csv"),
                     "stage-effect scaling 2026-11-24, band profile"),
    "scaling_calibration": (f"{SCALING}/calibration", ("calibration.json",),
                            "stage-effect scaling 2026-11-24, calibration"),
    "adni_gain_simulation": (f"{ADNI_GAIN}/simulation", ("settings.json", "rate_table.csv"),
                             "ADNI gain and specificity 2026-11-17, simulation profile"),
}
PROFILE_PARENTS: dict[str, tuple[str, ...]] = {
    "committed": ("univariate_stage_span", "scaling_saturated", "scaling_band", "scaling_calibration"),
    "conditions": ("univariate_stage_span",),
    "dispersion": ("univariate_stage_span",),
    "surgery": ("univariate_surgery",),
    "adni_design": ("adni_gain_simulation",),
}

VERDICTS: tuple[str, ...] = ("reproduction_failed", "propose_default_change", "scoped_candidates", "keep_production")
ARM_VERDICTS: tuple[str, ...] = ("qualifies_everywhere", "scoped_candidate", "keep_production")
FAMILY_LABELS: tuple[str, ...] = ("size_excess", "power_loss", "specificity_excess", "no_gain", "qualifies", "not_run")
_LABEL_OF = {"S1": "size_excess", "S2": "power_loss", "S3": "specificity_excess", "G": "no_gain"}
INCIDENTAL_LABELS: tuple[str, ...] = ("candidate_loses", "candidate_gains", "unresolved")
_LIMIT = (
    " Every reading is simulated: 658 Gaussian features on three blocks at 1 to 8 features per sample; the SEA-AD "
    "and ADNI real-noise readings did not enter the decision and stay recorded as limits, and size at many features "
    "per sample is unread (roadmap open item 15). No verdict changes what the harness or `motco analyze` runs."
)
CONSEQUENCES: dict[str, str] = {
    "reproduction_failed": (
        "Nothing is read: a rerun of committed designs does not reproduce its committed counts, so either the "
        "generator or the fit is not unchanged or the environment drifted. Find the drift before any rerun; a "
        "mismatch is diagnosed, never waived. The criteria are recorded for the diagnosis only."
    ),
    "propose_default_change": (
        "A separate change flips the named arm's default in the harness and `motco analyze`, with pins for every "
        "committed configuration and a re-issue of the case studies." + _LIMIT
    ),
    "scoped_candidates": (
        "The arms stay opt-in. The measurement-visibility page carries the map of where each candidate is better "
        "and worse than production; a rule that recommends an arm from a cohort's design is a separate change."
        + _LIMIT
    ),
    "keep_production": "The arms stay opt-in diagnostics, now read at paper grade in simulation." + _LIMIT,
}
DECISION_NAME = "candidate_simulation_decision.json"


# ── Small helpers ────────────────────────────────────────────────────────────


def _finite(value: Any) -> float | None:
    if value is None or value == "":
        return None
    number = float(value)
    return number if math.isfinite(number) else None


def _quantile(values: Iterable[Any], q: float) -> float | None:
    finite = [float(v) for v in values if _finite(v) is not None]
    return float(np.quantile(finite, q)) if finite else None


def _median(values: Iterable[Any]) -> float | None:
    return _quantile(values, 0.5)


def _share(values: Iterable[Any]) -> float | None:
    known = [bool(v) for v in values if v is not None]
    return sum(known) / len(known) if known else None


def rejects(p_value: Any, rule: str) -> bool | None:
    """The rejection of one p-value under ``rule``; ``None`` for an unavailable p-value."""

    value = _finite(p_value)
    if value is None:
        return None
    return value < ALPHA if rule == STRICT else value <= ALPHA


def profile_dir(out_dir: Path, profile: str) -> Path:
    if profile not in PROFILES:
        raise SystemExit(f"--profile must be one of {list(PROFILES)}, got {profile!r}.")
    return Path(out_dir) / profile


def record_key(record: Mapping[str, Any]) -> str:
    return json.dumps([str(record["design"]), str(record["cell"]), int(record["replicate"]), str(record["arm"])])


def unit_arm_key(unit: lsc.Unit, arm: str) -> str:
    return json.dumps([unit.design, unit.cell, int(unit.replicate), arm])


def _canonical(value: Any) -> Any:
    return json.loads(json.dumps(value, sort_keys=True))


# ── Parents ──────────────────────────────────────────────────────────────────


def parent_dir(name: str) -> Path:
    return RESULTS / PARENTS[name][0]


def bind_parents(names: Sequence[str]) -> dict[str, Any]:
    """Each parent's directory, bound files (sha256) and source label."""

    out = {}
    for name in names:
        directory, files, source = PARENTS[name]
        out[name] = {"dir": f"results/{directory}", "files": ubc._bind(parent_dir(name) / f for f in files),
                     "source": source}
    return out


def committed_settings(name: str) -> dict[str, Any]:
    """A parent's committed ``settings.json`` settings, refused when its hash fails."""

    return ssc.load_parent(parent_dir(name))["settings"]


def require_rebuilt(name: str, rebuilt: Mapping[str, Any]) -> None:
    """The parent driver's builder at its committed grade must rebuild the committed settings exactly."""

    committed = _canonical(committed_settings(name))
    if _canonical(dict(rebuilt)) != committed:
        differing = lsc.settings_difference(committed, _canonical(dict(rebuilt)))
        raise SystemExit(
            f"{parent_dir(name)}/settings.json differs from the settings its driver rebuilds at its committed grade; "
            f"differing: {', '.join(differing[:12])}. The units would not be the parent's."
        )


def check_bound_files(settings: Mapping[str, Any]) -> None:
    """Refuse, by name, a bound parent file that is missing or whose sha256 differs."""

    for stored, digest in sorted(settings.get("parent_files", {}).items()):
        path = ssc._resolve(stored)
        if not path.exists():
            raise SystemExit(f"{path} does not exist; the committed parent file is required.")
        if ssc.sha256_file(path) != digest:
            raise SystemExit(f"{path} has changed since the settings were written (sha256 differs); refusing.")


# ── Groups: parent-shaped unit settings ──────────────────────────────────────


def _harness_settings(
    generators: Mapping[str, Mapping[str, Any]],
    cells: Sequence[Sequence[Any]],
    *,
    replicates: int,
    control_replicates: int,
    permutations: int,
) -> dict[str, Any]:
    """The univariate baseline's ``stage_span`` settings with only the designs and cells replaced.

    Every other key (CV knobs, rank floor, estimator, base seed and seed functions) is the one
    ``univariate_baseline_check.stage_span_settings`` writes, so its ``unit_fields`` and ``unit_inputs`` build the
    parent's datasets and production parameters.
    """

    settings = ubc.stage_span_settings(
        designs=["n200"], replicates=replicates, control_replicates=control_replicates, permutations=permutations,
        base_seed=BASE_SEED, recorded_spaces=ubc.DEFAULT_RECORDED_SPACES,
    )
    settings["units"]["designs"] = {name: {"generator": dict(generator), "layers": None}
                                    for name, generator in generators.items()}
    settings["units"]["cells"] = [list(cell) for cell in cells]
    settings["designs"] = list(generators)
    return settings


def _group(name: str, kind: str, settings: Mapping[str, Any]) -> dict[str, Any]:
    return {"name": name, "kind": kind, "settings": dict(settings)}


def _check_committed_unit(committed: Mapping[str, Any], group: Mapping[str, Any], name: str, design: str) -> None:
    """A committed design's generator, cells and unit knobs equal the parent's, so its seeds and units do."""

    units = committed["unit_settings"]["units"] if "unit_settings" in committed else committed["units"]
    mine = group["units"]
    differing = []
    if units["designs"].get(design) != mine["designs"][design]:
        differing.append(f"designs.{design}")
    if [list(cell) for cell in units["cells"]] != [list(cell) for cell in mine["cells"]]:
        differing.append("cells")
    differing += [key for key in ("base_seed", "pls_cv_params", "min_components", "orientation_estimator")
                  if units[key] != mine[key]]
    if differing:
        raise SystemExit(f"{parent_dir(name)}/settings.json differs from this driver's {design!r} in {differing}; "
                         "its units would not be rerun.")


def committed_groups(*, replicates: int, control_replicates: int, permutations: int,
                     designs: Sequence[str] | None = None) -> list[dict[str, Any]]:
    """The thirteen-cell committed designs and the band design, each checked against its parent."""

    require_rebuilt("univariate_stage_span", ubc.stage_span_settings())
    require_rebuilt("scaling_saturated", ses.saturated_settings())
    require_rebuilt("scaling_band", ses.band_settings(RESULTS / SCALING))
    chosen = list(COMMITTED_DESIGNS) + [BAND_DESIGN] if designs is None else list(designs)
    groups = []
    stage_span = {name: COMMITTED_DESIGNS[name][1:] for name in COMMITTED_DESIGNS if name in chosen}
    if stage_span:
        group = ses.unit_settings(stage_span, ubc.STAGE_SPAN_CELLS, replicates=replicates,
                                  control_replicates=control_replicates, permutations=permutations)
        for design in stage_span:
            parent = COMMITTED_DESIGNS[design][0]
            _check_committed_unit(committed_settings(parent), group, parent, design)
        groups.append(_group("stage_span", "stage_span", group))
    if BAND_DESIGN in chosen:
        band = committed_settings("scaling_band")
        delta = float(band["calibration"]["delta"])
        size, multiple = ses.BAND_DESIGNS[BAND_DESIGN]
        group = ses.unit_settings({BAND_DESIGN: (size, multiple * delta)}, ses.BAND_CELLS, replicates=replicates,
                                  control_replicates=control_replicates, permutations=permutations)
        _check_committed_unit(band, group, "scaling_band", BAND_DESIGN)
        groups.append(_group("band", "stage_span", group))
    return groups


def condition_generator(name: str) -> dict[str, Any]:
    """A condition design's generator: the baseline's three-block design with the condition's change."""

    spec = CONDITION_DESIGNS[name]
    generator = ubc.stage_span_generator(int(spec.get("n_samples", 200)), float(spec["delta"]))
    if "n_stages" in spec:
        generator["n_stages"] = int(spec["n_stages"])
    if "share" in spec:
        generator["group_ratio"] = float(spec["share"])
    if "table" in spec:
        # The table fixes every cell; the proportional sample size would conflict with it.
        generator.pop("n_samples")
        generator["group_stage_sizes"] = [list(row) for row in spec["table"]]
    return generator


def conditions_groups(*, replicates: int, control_replicates: int, permutations: int,
                      designs: Sequence[str] | None = None) -> list[dict[str, Any]]:
    chosen = [name for name in CONDITION_DESIGNS if designs is None or name in designs]
    generators = {name: condition_generator(name) for name in chosen}
    group = _harness_settings(generators, ubc.STAGE_SPAN_CELLS, replicates=replicates,
                              control_replicates=control_replicates, permutations=permutations)
    return [_group("stage_span", "stage_span", group)]


def dispersion_design_name(reference: str, pattern: str) -> str:
    return f"{reference}_{pattern}"


def dispersion_entries() -> list[dict[str, Any]]:
    """The ten dispersion designs: (name, its equal-dispersion reference, the pattern)."""

    return [{"design": dispersion_design_name(reference, pattern), "reference_design": reference, "pattern": pattern}
            for reference, patterns in DISPERSION_REFERENCES.items() for pattern in patterns]


def reference_generator(reference: str) -> dict[str, Any]:
    """The equal-dispersion design's generator: a condition design's, or the committed ``n200``'s."""

    if reference in CONDITION_DESIGNS:
        return condition_generator(reference)
    units = committed_settings("univariate_stage_span")["units"]
    return dict(units["designs"][reference]["generator"])


def dispersion_generator(reference: str, pattern: str) -> dict[str, Any]:
    return {**reference_generator(reference), "stage_dispersion_scale": list(DISPERSION_PATTERNS[pattern])}


def dispersion_groups(*, replicates: int, control_replicates: int, permutations: int,
                      designs: Sequence[str] | None = None) -> list[dict[str, Any]]:
    entries = [entry for entry in dispersion_entries() if designs is None or entry["design"] in designs]
    generators = {entry["design"]: dispersion_generator(entry["reference_design"], entry["pattern"])
                  for entry in entries}
    cells = [cell for cell in ubc.STAGE_SPAN_CELLS if cell[0] in CONTROL_CELLS]
    group = _harness_settings(generators, cells, replicates=replicates, control_replicates=control_replicates,
                              permutations=permutations)
    return [_group("stage_span", "stage_span", group)]


def surgery_groups(*, replicates: int, control_replicates: int, permutations: int,
                   designs: Sequence[str] | None = None) -> list[dict[str, Any]]:
    require_rebuilt("univariate_surgery", ubc.surgery_settings())
    columns = [column for column in SURGERY_COLUMNS if designs is None or column in designs]
    group = ubc.surgery_settings(columns=columns, replicates=replicates, control_replicates=control_replicates,
                                 permutations=permutations, recorded_spaces=ubc.DEFAULT_RECORDED_SPACES)
    return [_group("surgery", "surgery", group)]


def adni_groups(*, replicates: int, control_replicates: int, permutations: int,
                designs: Sequence[str] | None = None) -> list[dict[str, Any]]:
    """The gain run's simulation settings at this grade, with its two operating arms only."""

    require_rebuilt("adni_gain_simulation", ags.simulation_settings())
    if designs is not None and ADNI_DESIGN not in designs:
        raise SystemExit(f"The adni_design profile has one design, {ADNI_DESIGN!r}.")
    gain = ags.simulation_settings(control_replicates=control_replicates, plant_replicates=replicates,
                                   permutations=permutations)
    gain["arms"] = [ags.REFERENCE_ARM, ags.CANDIDATE_ARM]
    gain["arm_specs"] = {arm: gain["arm_specs"][arm] for arm in gain["arms"]}
    return [_group("adni", "adni", gain)]


_GROUP_BUILDERS = {
    "committed": committed_groups, "conditions": conditions_groups, "dispersion": dispersion_groups,
    "surgery": surgery_groups, "adni_design": adni_groups,
}


def _known_designs(profile: str) -> list[str]:
    if profile == "committed":
        return [*COMMITTED_DESIGNS, BAND_DESIGN]
    if profile == "conditions":
        return list(CONDITION_DESIGNS)
    if profile == "dispersion":
        return [entry["design"] for entry in dispersion_entries()]
    if profile == "surgery":
        return list(SURGERY_COLUMNS)
    return [ADNI_DESIGN]


# ── Profile settings ─────────────────────────────────────────────────────────


def _group_designs(group: Mapping[str, Any]) -> list[str]:
    settings = group["settings"]
    if group["kind"] == "adni":
        return ags.point_names(settings)
    return list(settings["designs"])


def _design_meta(profile: str, group: Mapping[str, Any], design: str) -> dict[str, Any]:
    settings = group["settings"]
    if group["kind"] == "adni":
        generator = settings["construction"]["generator"]
        table = [list(row) for row in ags.point_spec(settings, design)["table"]]
    elif group["kind"] == "surgery":
        study = ubc._study_columns(settings["sources"]["geometry"]["config"])
        params = next(iter(study[design].values())).generator_params
        generator = {key: getattr(params, key) for key in ("n_stages", "delta_methyl", "group_ratio")}
        generator["n_samples"] = int(settings["n_samples"][design])
        table = None
    else:
        generator = settings["units"]["designs"][design]["generator"]
        table = generator.get("group_stage_sizes")
    dispersion = next((entry for entry in dispersion_entries() if entry["design"] == design), None)
    return {
        "profile": profile,
        "group": group["name"],
        "family": FAMILY_OF.get(design),
        "delta": float(generator["delta_methyl"]),
        "n_stages": int(generator.get("n_stages", 3)),
        "n_samples": int(sum(map(sum, table))) if table is not None else int(generator["n_samples"]),
        "group_stage_sizes": None if table is None else [list(map(int, row)) for row in table],
        "share": None if table is not None else float(generator.get("group_ratio", 0.5)),
        "stage_dispersion_scale": generator.get("stage_dispersion_scale"),
        "pattern": None if dispersion is None else dispersion["pattern"],
        "reference_design": None if dispersion is None else dispersion["reference_design"],
        "sparse": design in SPARSE_DESIGNS,
    }


def _reproduction_plan(profile: str, groups: Sequence[Mapping[str, Any]]) -> list[dict[str, Any]]:
    """Which designs are compared with which committed rate table, and for which arms."""

    if profile not in REPRODUCED_PROFILES:
        return []
    plan = []
    for group in groups:
        for design in _group_designs(group):
            if profile == "committed":
                parent = "scaling_band" if design == BAND_DESIGN else COMMITTED_DESIGNS[design][0]
                plan.append({"design": design, "parent": parent, "arms": [REFERENCE_ARM], "table": "stage_span"})
            elif profile == "surgery":
                plan.append({"design": design, "parent": "univariate_surgery", "arms": [REFERENCE_ARM],
                             "table": "stage_span"})
            else:
                plan.append({"design": design, "parent": "adni_gain_simulation",
                             "arms": [REFERENCE_ARM, PREDICTIVE_ARM], "table": "adni"})
    return plan


def arm_catalogue() -> dict[str, Any]:
    """What each arm is; the arms differ in nothing else."""

    ladder = list(plsc.default_ladder())
    return {
        REFERENCE_ARM: "production: pls on the unit's blocks, rank by double CV (floor 2), net displacement, scores",
        PREDICTIVE_ARM: "the production fit's predictive subspace in the input metric (latent_coordinates="
                        "'predictive_orthonormal'), built by pls_coordinates_from_pls with no second fit; on "
                        "adni_design the gain run's forced call at pls_cv's rank",
        WEIGHTED_ARM: {"integration_params": {"block_weighting": WEIGHTING},
                       "fit": "its own cross-validated fit; delta in the weighted metric"},
        SPARSE_ARM: {"integration_params": {"sparsity": {"scope": SPARSITY_SCOPE, "keep_fractions": ladder}},
                     "fit": "its own cross-validated fit over rank x keep fraction", "designs": list(SPARSE_DESIGNS)},
    }


def profile_settings(
    profile: str,
    *,
    replicates: int | None = None,
    control_replicates: int | None = None,
    permutations: int = DEFAULT_PERMUTATIONS,
    designs: Sequence[str] | None = None,
) -> dict[str, Any]:
    """One profile's settings: its groups of parent-shaped unit settings, its designs, its parents."""

    if profile not in PROFILES:
        raise SystemExit(f"Unknown profile {profile!r}; expected one of {list(PROFILES)}.")
    if designs is not None:
        unknown = sorted(set(designs) - set(_known_designs(profile)))
        if unknown:
            raise SystemExit(f"Unknown {profile} design(s) {unknown}; expected a subset of {_known_designs(profile)}.")
    defaults = DEFAULT_REPLICATES[profile]
    reps = defaults["other"] if replicates is None else int(replicates)
    control = defaults["control"] if control_replicates is None else int(control_replicates)
    if profile == "dispersion":
        reps = control  # null only: every cell is a control
    groups = _GROUP_BUILDERS[profile](replicates=reps, control_replicates=control, permutations=int(permutations),
                                      designs=designs)
    parents = bind_parents(PROFILE_PARENTS[profile])
    files: dict[str, str] = {}
    for block in parents.values():
        files.update(block["files"])
    for group in groups:
        settings = group["settings"]
        files.update(settings.get("parent_files") or {})
        for entry in (settings.get("parents") or {}).values():
            files[entry["path"]] = entry["sha256"]
    design_index = {design: _design_meta(profile, group, design)
                    for group in groups for design in _group_designs(group)}
    return {
        "driver": DRIVER,
        "driver_version": DRIVER_VERSION,
        "rule_version": RULE_VERSION,
        "profile": profile,
        "alpha": ALPHA,
        "rejection": REJECTION_RULES[profile],
        "permutations": int(permutations),
        "replicates": {"control": control, "other": reps},
        "arms": list(ARMS),
        "arm_definitions": arm_catalogue(),
        "sparse_designs": list(SPARSE_DESIGNS),
        "recorded": ({"within_stage": list(WITHIN_STAGE_ARMS), "residual_dispersion": list(ARMS[:3])}
                     if profile == "dispersion" else None),
        "parents": parents,
        "parent_files": dict(sorted(files.items())),
        "groups": groups,
        "design_index": design_index,
        "reproduction": _reproduction_plan(profile, groups),
    }


# ── Units ────────────────────────────────────────────────────────────────────


def group_units(group: Mapping[str, Any]) -> list[lsc.Unit]:
    """Design-major, then cell, then replicate; controls at the control count."""

    settings = group["settings"]
    if group["kind"] == "adni":
        return [lsc.Unit(point, cell["name"], replicate)
                for point in ags.point_names(settings)
                for cell in ags.point_cells(settings, point)
                for replicate in range(ags.cell_replicates(settings, cell))]
    return ubc.profile_units(settings)


def profile_units(settings: Mapping[str, Any]) -> list[lsc.Unit]:
    return [unit for group in settings["groups"] for unit in group_units(group)]


def group_of(settings: Mapping[str, Any], design: str) -> dict[str, Any]:
    for group in settings["groups"]:
        if design in _group_designs(group):
            return group
    raise KeyError(f"Design {design!r} is not in the {settings['profile']} profile.")


def arm_applies(settings: Mapping[str, Any], design: str, arm: str) -> bool:
    if arm == SPARSE_ARM:
        return design in settings["sparse_designs"]
    return arm in ARMS


def _adni_cell_fields(cell: Mapping[str, Any]) -> dict[str, Any]:
    return {"trajectory_mode": cell["mode"] or "none", "effect_size": float(cell["effect"] or 0.0),
            "span_alignment": None, "span_support": None, "span_scale": "planted"}


def unit_fields(unit: lsc.Unit, settings: Mapping[str, Any]) -> dict[str, Any]:
    """The unit's cell description and both seeds (``None`` on ``adni_design`` until built)."""

    group = group_of(settings, unit.design)
    fields: dict[str, Any] = {"profile": settings["profile"], "design": unit.design, "cell": unit.cell,
                              "replicate": int(unit.replicate)}
    if group["kind"] == "adni":
        fields.update(_adni_cell_fields(ags.cell_spec(group["settings"], unit.cell)),
                      generator_seed=None, rrpp_seed=None)
        return fields
    parent = ubc.unit_fields(unit, group["settings"])
    fields.update({key: parent[key] for key in ("trajectory_mode", "effect_size", "span_alignment", "span_support",
                                                "span_scale", "generator_seed", "rrpp_seed")})
    return fields


def cell_kind(settings: Mapping[str, Any], design: str, cell: str) -> dict[str, Any]:
    """The cell's mode, effect, kind (control, in_span, out_of_span, surgery) and target statistic."""

    group = group_of(settings, design)
    if group["kind"] == "adni":
        spec = ags.cell_spec(group["settings"], cell)
        fields = _adni_cell_fields(spec)
        kind = "control" if spec["kind"] == ags.CONTROL else "in_span"
    else:
        fields = unit_fields(lsc.Unit(design, cell, 0), settings)
        mode, alignment = fields["trajectory_mode"], fields["span_alignment"]
        if cell in CONTROL_CELLS:
            kind = "control"
        elif group["kind"] == "surgery":
            kind = "surgery"
        elif alignment is not None and float(alignment) > 0:
            kind = "out_of_span"
        else:
            kind = "in_span"
    mode = fields["trajectory_mode"]
    return {"trajectory_mode": mode, "effect_size": fields["effect_size"], "span_alignment": fields["span_alignment"],
            "span_support": fields["span_support"], "kind": kind,
            "target": TARGETS.get(mode) if kind != "control" else None}


def design_cells(settings: Mapping[str, Any], design: str) -> list[str]:
    group = group_of(settings, design)
    if group["kind"] == "adni":
        return [cell["name"] for cell in ags.point_cells(group["settings"], design)]
    return ubc._cell_names(group["settings"])


def expected_counts(settings: Mapping[str, Any]) -> dict[tuple[str, str], int]:
    counts: dict[tuple[str, str], int] = {}
    for unit in profile_units(settings):
        counts[(unit.design, unit.cell)] = counts.get((unit.design, unit.cell), 0) + 1
    return counts


# ── One unit: the arms ───────────────────────────────────────────────────────


def _failed(base: Mapping[str, Any], arm: str, exc: BaseException | str, seconds: float = 0.0) -> dict[str, Any]:
    if isinstance(exc, BaseException):
        error = {"error_type": type(exc).__name__, "error_message": str(exc)}
    else:
        error = {"error_type": "Skipped", "error_message": exc}
    return {**base, "arm": arm, "status": "failed", **error, "seconds": {"arm": seconds}}


def harness_arm_record(result: Any, seconds: float) -> dict[str, Any]:
    """One harness arm's rank, dimension, AUROC, statistics, p-values, capture, population statistics, sparsity."""

    from motco.simulations.diagnostics import POPULATION_LATENT_CHECKPOINT

    metadata = result.latent_matrix_metadata
    geometry = result.realized_geometry or {}
    capture = geometry.get("span_capture") or {}
    population = ((geometry.get("checkpoints") or {}).get(POPULATION_LATENT_CHECKPOINT) or {}).get("joint") or {}
    sparsity = metadata.get("sparsity") or {}
    return {
        "status": "completed",
        "rank": int(metadata["selected_lv"]),
        "component_selection": metadata.get("component_selection"),
        "measured_dimension": int(metadata.get("measured_dimension", metadata.get("n_features", 0))),
        "cv_mean_auroc": metadata.get("cv_mean_auroc"),
        "pair_statistics": lsc._finite(result.pair_statistics),
        "p_values": lsc._finite(result.p_values),
        "span_capture": {key: capture.get(key) for key in ("pooled_stage", "group_difference")},
        "population": {statistic: population.get(statistic) for statistic in STATISTICS},
        "keep_fraction": sparsity.get("selected_keep_fraction"),
        "support": sparsity.get("support"),
        "seconds": {"arm": seconds},
    }


def within_stage_reading(dataset: Any, params: Any, latent: Any, unrestricted: Mapping[str, Any]) -> dict[str, Any]:
    """The arm's statistics and p-values under within-stage draws on the same latent matrix and seed."""

    from motco.simulations.evaluation import evaluate_semisynthetic_trajectory

    mark = perf_counter()
    result = evaluate_semisynthetic_trajectory(dataset, dataclasses.replace(params, permutation_strata="stage"),
                                               latent=latent)
    statistics = lsc._finite(result.pair_statistics)
    return {
        "pair_statistics": statistics,
        "p_values": lsc._finite(result.p_values),
        "same_statistics": all(statistics.get(s) == unrestricted.get(s) for s in STATISTICS),
        "seconds": perf_counter() - mark,
    }


def residual_dispersion(dataset: Any, params: Any, latent: Any) -> dict[str, Any]:
    """The reduced model's residual dispersion by stage in the arm's measured coordinates, over the pooled one."""

    from motco.simulations.evaluation import build_simulation_trajectory_design

    design = build_simulation_trajectory_design(dataset.metadata, group_col=params.group_col,
                                                stage_col=params.stage_col)
    return pse.stage_dispersion(latent.matrix.to_numpy(dtype=float), design.model_reduced,
                                dataset.metadata[params.stage_col], design.stage_levels)


def _own_fit_params(params: Any, arm: str, settings: Mapping[str, Any]) -> Any:
    """Production's parameters with only the arm's ``integration_params`` (from the settings) added."""

    changes = settings["arm_definitions"][arm]["integration_params"]
    return dataclasses.replace(params, integration_params={**dict(params.integration_params or {}), **changes})


def _evaluate_harness(unit: lsc.Unit, group: Mapping[str, Any], settings: Mapping[str, Any],
                      arms: Sequence[str]) -> list[dict[str, Any]]:
    """Generate once; fit production once for ``pls_cv`` and the predictive arm; fit each own-fit arm."""

    from motco.simulations.evaluation import (
        evaluate_semisynthetic_trajectory,
        integrate_semisynthetic_dataset,
        pls_coordinates_from_pls,
    )

    start = perf_counter()
    base = {**unit_fields(unit, settings), "settings_hash": lsc.settings_hash(settings)}
    recorded = settings.get("recorded") is not None
    records: dict[str, dict[str, Any]] = {}
    shared: dict[str, float] = {}
    try:
        mark = perf_counter()
        dataset, params = ubc.unit_inputs(unit, group["settings"], base)
        shared["generate"] = perf_counter() - mark
    except Exception as exc:  # every arm fails with the dataset; the shard goes on
        return _finish_call([_failed(base, arm, exc) for arm in arms], shared, start)

    if REFERENCE_ARM in arms or PREDICTIVE_ARM in arms:
        latent = None
        try:
            mark = perf_counter()
            latent = integrate_semisynthetic_dataset(dataset, params)
            shared["fit"] = perf_counter() - mark
        except Exception as exc:
            for arm in (REFERENCE_ARM, PREDICTIVE_ARM):
                if arm in arms:
                    records[arm] = _failed(base, arm, exc)
        production: dict[str, Any] | None = None
        if latent is not None and REFERENCE_ARM in arms:
            mark = perf_counter()
            try:
                result = evaluate_semisynthetic_trajectory(dataset, params, latent=latent)
                production = harness_arm_record(result, perf_counter() - mark)
                if recorded:
                    production["within_stage"] = within_stage_reading(dataset, params, latent,
                                                                      production["pair_statistics"])
                    production["stage_dispersion"] = residual_dispersion(dataset, params, latent)
                records[REFERENCE_ARM] = {**base, "arm": REFERENCE_ARM, **production}
            except Exception as exc:
                records[REFERENCE_ARM] = _failed(base, REFERENCE_ARM, exc, perf_counter() - mark)
        if latent is not None and PREDICTIVE_ARM in arms:
            mark = perf_counter()
            try:
                predictive_params = dataclasses.replace(params, latent_coordinates=PREDICTIVE_ARM)
                predictive = pls_coordinates_from_pls(dataset, latent, predictive_params)
                result = evaluate_semisynthetic_trajectory(dataset, predictive_params, latent=predictive)
                entry = harness_arm_record(result, perf_counter() - mark)
                entry["consistency"] = predictive_consistency(entry, latent, production,
                                                              n_stages=int(dataset.metadata[params.stage_col].nunique()))
                if recorded:
                    entry["within_stage"] = within_stage_reading(dataset, predictive_params, predictive,
                                                                 entry["pair_statistics"])
                    entry["stage_dispersion"] = residual_dispersion(dataset, predictive_params, predictive)
                records[PREDICTIVE_ARM] = {**base, "arm": PREDICTIVE_ARM, **entry}
            except Exception as exc:
                records[PREDICTIVE_ARM] = _failed(base, PREDICTIVE_ARM, exc, perf_counter() - mark)

    for arm in (WEIGHTED_ARM, SPARSE_ARM):
        if arm not in arms:
            continue
        if not arm_applies(settings, unit.design, arm):
            continue  # not run here; the merged tables mark it
        mark = perf_counter()
        try:
            own = _own_fit_params(params, arm, settings)
            latent_own = integrate_semisynthetic_dataset(dataset, own)
            result = evaluate_semisynthetic_trajectory(dataset, own, latent=latent_own)
            entry = harness_arm_record(result, perf_counter() - mark)
            if recorded:
                entry["stage_dispersion"] = residual_dispersion(dataset, own, latent_own)
            records[arm] = {**base, "arm": arm, **entry}
        except Exception as exc:
            records[arm] = _failed(base, arm, exc, perf_counter() - mark)
    return _finish_call([records[arm] for arm in arms if arm in records], shared, start)


def predictive_consistency(entry: Mapping[str, Any], latent: Any, production: Mapping[str, Any] | None, *,
                           n_stages: int) -> dict[str, Any]:
    """Recorded, never required: rank and AUROC equal the production fit's; p-values at a fit rank ≤ K − 1."""

    metadata = latent.metadata
    rank, auroc = int(metadata["selected_lv"]), metadata.get("cv_mean_auroc")
    out: dict[str, Any] = {
        "rank_equal": int(entry["rank"]) == rank,
        "auroc_equal": entry.get("cv_mean_auroc") == auroc,
        "low_rank": rank <= n_stages - 1,
        "identity_at_low_rank": None,
    }
    if out["low_rank"] and production is not None:
        mine, theirs = entry["p_values"], production["p_values"]
        out["identity_at_low_rank"] = all(mine.get(s) == theirs.get(s) for s in STATISTICS)
    return out


def _finish_call(records: list[dict[str, Any]], shared: Mapping[str, float], start: float) -> list[dict[str, Any]]:
    """Each record's timing: its arm's, the call's shared generation and fit, and its share of the call."""

    total = perf_counter() - start
    for record in records:
        seconds = dict(record.get("seconds") or {})
        seconds.update({key: value for key, value in shared.items()})
        seconds["call"] = total
        seconds["call_share"] = total / len(records)
        seconds["call_arms"] = len(records)
        record["seconds"] = seconds
    return [lsc._jsonable(record) for record in records]


def _evaluate_adni(unit: lsc.Unit, group: Mapping[str, Any], settings: Mapping[str, Any],
                   arms: Sequence[str]) -> list[dict[str, Any]]:
    """The gain run's simulated unit: its two operating arms on one CV fit, and the weighted CV call."""

    start = perf_counter()
    gain = group["settings"]
    construction = gain["construction"]
    point = ags.point_spec(gain, unit.design)
    cell = ags.cell_spec(gain, unit.cell)
    mode, effect = ags.construction_mode(cell)
    base = {**unit_fields(unit, settings), "settings_hash": lsc.settings_hash(settings)}
    shared: dict[str, float] = {}
    records: dict[str, dict[str, Any]] = {}
    try:
        mark = perf_counter()
        built = adc.build_unit(mode, effect, unit.replicate, construction, reference=ags._reference())
        seed = int(built["seeds"]["rrpp"])
        base.update(generator_seed=int(built["seeds"].get("generator", seed)), rrpp_seed=seed,
                    seeds=dict(built["seeds"]))
        dataset = adc.trajectory_dataset(built, construction, covariates=True)
        shared["generate"] = perf_counter() - mark
    except Exception as exc:
        return _finish_call([_failed(base, arm, exc) for arm in arms], shared, start)

    timing: dict[str, float] = {}
    results: dict[str, Any] = {}
    entries: dict[str, Any] = {}
    forced: dict[str, int] = {}
    if REFERENCE_ARM in arms or PREDICTIVE_ARM in arms:
        try:
            _, entries, results, forced = ags.evaluate_arms(dataset, gain, seed, cross_validate=True,
                                                            forced_rank=None, timing=timing)
        except Exception as exc:
            entries = {arm: {**aac.failed_arm(arm, exc, 0.0, gain["arm_specs"][arm])} for arm in gain["arms"]}
    if WEIGHTED_ARM in arms:
        mark = perf_counter()
        try:
            weighting = settings["arm_definitions"][WEIGHTED_ARM]["integration_params"]["block_weighting"]
            params = dataclasses.replace(aac.analysis_params(gain, seed), block_weighting=weighting)
            result = aac._analyse(dataset, params)
            entries[WEIGHTED_ARM] = {**psc.analysis_arm_record(result, perf_counter() - mark),
                                     "fit": "cross_validated", "block_weighting": weighting}
        except Exception as exc:
            entries[WEIGHTED_ARM] = psc._failed_arm(exc, perf_counter() - mark)

    try:
        failing = adni_checks(built, dataset, cell, point, gain, entries, results, forced)
    except Exception as exc:  # a check that cannot be computed fails every arm of the unit; the shard goes on
        return _finish_call([_failed(base, arm, exc) for arm in arms], shared, start)
    for arm in arms:
        entry = entries.get(arm)
        if entry is None:
            continue
        record = {**base, "arm": arm, **{k: v for k, v in entry.items() if k != "seconds"},
                  "seconds": {"arm": entry.get("seconds") or timing.get("production", 0.0)}}
        if record.get("status") == "completed" and failing:
            record["status"] = "inconsistent"
        record["inconsistent_checks"] = failing
        records[arm] = record
    shared.update({f"adni_{key}": value for key, value in timing.items()})
    return _finish_call([records[arm] for arm in arms if arm in records], shared, start)


def adni_checks(built: Mapping[str, Any], dataset: Any, cell: Mapping[str, Any], point: Mapping[str, Any],
                gain: Mapping[str, Any], entries: Mapping[str, Any], results: Mapping[str, Any],
                forced: Mapping[str, int]) -> list[str]:
    """The gain run's checks on its two arms and the construction; the failing ones' names.

    A failing check marks every arm of the unit ``inconsistent``, as the gain run left every rate.
    """

    stages = list(dataset.stages)
    groups = list(built["group_labels"])
    metadata = built["metadata"]
    labels = [[int(((metadata["group"].astype(str) == g) & (metadata["stage"].astype(str) == s)).sum())
               for s in stages] for g in groups]
    checks: dict[str, Any] = {"construction": ags.check_construction(built["dataset"].truth, cell)}
    if results:
        own = {arm: entries[arm] for arm in gain["arms"] if arm in entries}
        realized = {"labels": labels, **{arm: ags.recorded_table(result, groups, stages) for arm, result in
                                         results.items()}}
        checks.update(ags.common_checks(arms=own, results=results, forced=forced, settings=gain,
                                        n_stages=len(stages), expected_table=point["table"], realized=realized))
    else:
        checks["cell_table"] = aac.check_cell_table(point["table"], {"labels": labels})
    return sorted(name for name, check in checks.items() if not check["passed"])


def evaluate_unit(unit: lsc.Unit, settings: Mapping[str, Any], arms: Sequence[str]) -> list[dict[str, Any]]:
    """One unit's records for ``arms`` (the applicable ones); a failure is recorded, never raised."""

    group = group_of(settings, unit.design)
    wanted = [arm for arm in ARMS if arm in arms and arm_applies(settings, unit.design, arm)]
    if not wanted:
        return []
    if group["kind"] == "adni":
        return _evaluate_adni(unit, group, settings, wanted)
    return _evaluate_harness(unit, group, settings, wanted)


# ── Settings files and the run ───────────────────────────────────────────────


def parse_arms(text: str) -> tuple[str, ...]:
    """``all`` or a comma list of arms, in the arms' order."""

    if text.strip() == "all":
        return ARMS
    names = [name.strip() for name in text.split(",") if name.strip()]
    unknown = sorted(set(names) - set(ARMS))
    if unknown or not names:
        raise SystemExit(f"--arms must be 'all' or a comma list of {list(ARMS)}; got {text!r}.")
    return tuple(arm for arm in ARMS if arm in names)


def _write_atomic(path: Path, payload: Mapping[str, Any]) -> None:
    temporary = path.with_name(f".{path.name}.{os.getpid()}.tmp")
    temporary.write_text(json.dumps(payload, indent=1, sort_keys=True) + "\n", encoding="utf-8")
    os.replace(temporary, path)


def ensure_profile_settings(directory: Path, settings: Mapping[str, Any], arms: Sequence[str]) -> None:
    """Write the settings, or refuse other settings naming the keys; record the arm selection outside the hash."""

    path = Path(directory) / "settings.json"
    selection = list(arms)
    if path.exists():
        stored = json.loads(path.read_text(encoding="utf-8"))
        if stored.get("hash") != lsc.settings_hash(settings):
            differing = lsc.settings_difference(stored.get("settings", {}), _canonical(dict(settings)))
            raise SystemExit(f"{path} was written with different settings; differing: {', '.join(differing)}. "
                             "Use a fresh --out-dir.")
        selections = stored.get("arm_selections") or []
        if selection not in selections:
            _write_atomic(path, {**stored, "arm_selections": [*selections, selection]})
        return
    path.parent.mkdir(parents=True, exist_ok=True)
    _write_atomic(path, {"hash": lsc.settings_hash(settings), "settings": settings, "arm_selections": [selection]})


def arm_selections(directory: Path) -> list[list[str]]:
    path = Path(directory) / "settings.json"
    return json.loads(path.read_text(encoding="utf-8")).get("arm_selections") or [] if path.exists() else []


def read_records(paths: Iterable[Path]) -> list[dict[str, Any]]:
    """Every complete line (a line torn by an interrupted write is dropped, so its (unit, arm) is re-run)."""

    out = []
    for path in paths:
        with Path(path).open(encoding="utf-8") as handle:
            for line in handle:
                try:
                    record = json.loads(line)
                except json.JSONDecodeError:
                    continue
                if record.get("status") in ("completed", "failed", "inconsistent"):
                    out.append(record)
    return out


def run_shard(
    out_dir: Path,
    settings: Mapping[str, Any],
    index: int,
    count: int,
    *,
    arms: Sequence[str],
    limit: int | None = None,
    workers: int = 1,
) -> int:
    """Evaluate this shard's missing (unit, arm) keys into its file; returns how many units were evaluated."""

    profile = settings["profile"]
    check_bound_files(settings)
    directory = profile_dir(out_dir, profile)
    ensure_profile_settings(directory, settings, arms)
    path = lsc.shard_path(directory, index, count)
    path.parent.mkdir(parents=True, exist_ok=True)
    done = {record_key(record): record for record in read_records([path])} if path.exists() else {}
    # Rewritten with complete lines only, so an interrupted unit is re-run, not duplicated.
    with path.open("w", encoding="utf-8") as handle:
        for key in sorted(done):
            handle.write(lsc._line(done[key]))
    pending: list[tuple[lsc.Unit, tuple[str, ...]]] = []
    for position, unit in enumerate(profile_units(settings)):
        if position % count != index:
            continue
        wanted = tuple(arm for arm in arms if arm_applies(settings, unit.design, arm)
                       and unit_arm_key(unit, arm) not in done)
        if wanted:
            pending.append((unit, wanted))
    if limit is not None:
        pending = pending[:limit]

    def append(position: int, unit: lsc.Unit, records: Sequence[Mapping[str, Any]]) -> None:
        with path.open("a", encoding="utf-8") as handle:
            for record in records:
                handle.write(lsc._line(record))
            handle.flush()
        statuses = ",".join(f"{r['arm']}={r['status']}" for r in records)
        print(f"[{profile} {index}/{count}] {position}/{len(pending)} {unit.design}/{unit.cell}/{unit.replicate} "
              f"{statuses}", file=sys.stderr, flush=True)

    if workers <= 1:
        for position, (unit, wanted) in enumerate(pending, start=1):
            append(position, unit, evaluate_unit(unit, settings, wanted))
    else:
        with ProcessPoolExecutor(max_workers=workers) as pool:
            futures = {pool.submit(evaluate_unit, unit, settings, wanted): unit for unit, wanted in pending}
            for position, future in enumerate(as_completed(futures), start=1):
                append(position, futures[future], future.result())
    return len(pending)


# ── Merge ────────────────────────────────────────────────────────────────────


CELL_FIELDS: tuple[str, ...] = (
    "design", "cell", "trajectory_mode", "effect_size", "span_alignment", "span_support", "kind",
)
RATE_COLUMNS: tuple[str, ...] = (
    *CELL_FIELDS, "arm", "statistic", "role", "rule", "status", "expected", "records", "rate", "se", "count",
    "rejections", "failures", "inconsistent",
)
PAIRED_COLUMNS: tuple[str, ...] = (
    *CELL_FIELDS, "statistic", "role", "arm", "reference_arm", "comparison", "status", "n_paired", "rate_reference",
    "rate_arm", "difference", "b", "c", "p_loss", "p_gain", "p_two_sided",
)
SUPPORT_BLOCKS: tuple[str, ...] = ("methylation", "expression", "proteomics")
CELL_COLUMNS: tuple[str, ...] = (
    *CELL_FIELDS, "arm", "completed", "failures", "inconsistent", "rank_q10", "rank_median", "rank_q90",
    "auroc_median", "measured_dimension_median", "capture_pooled_stage_median", "capture_group_difference_median",
    "keep_fraction_q10", "keep_fraction_median", "keep_fraction_q90",
    *(f"support_{block}_median" for block in SUPPORT_BLOCKS),
    *(f"observed_{s}_median" for s in STATISTICS), *(f"population_{s}_median" for s in STATISTICS),
    *(f"dispersion_stage_{k}_median" for k in range(4)), "dispersion_max_min_median",
    "within_stage_same_statistics", "rank_equal", "auroc_equal", "low_rank_units", "identity_at_low_rank",
    "seconds_median",
)
RANK_COLUMNS: tuple[str, ...] = ("design", "cell", "rank", "arm", "statistic", "count", "rejections")


def _base_arm(arm: str) -> str:
    return arm.removesuffix(WITHIN_STAGE)


def _p(record: Mapping[str, Any] | None, arm: str, statistic: str) -> Any:
    """The p-value of ``arm`` (a within-stage arm reads its base arm's nested reading) on a completed record."""

    if record is None or record.get("status") != "completed":
        return None
    if arm.endswith(WITHIN_STAGE):
        return ((record.get("within_stage") or {}).get("p_values") or {}).get(statistic)
    return (record.get("p_values") or {}).get(statistic)


def rate_arms(settings: Mapping[str, Any]) -> tuple[str, ...]:
    return (*ARMS, *WITHIN_STAGE_ARMS) if settings.get("recorded") is not None else ARMS


def pair_specs(settings: Mapping[str, Any]) -> tuple[tuple[str, str, str], ...]:
    return (*PAIRS, *DISPERSION_PAIRS) if settings.get("recorded") is not None else PAIRS


def _role(kind: Mapping[str, Any], statistic: str) -> str:
    if kind["kind"] == "control":
        return "control"
    return "target" if statistic == kind["target"] else "off_target"


def _grouped(records: Iterable[Mapping[str, Any]]) -> dict[tuple[str, str, str], list[Mapping[str, Any]]]:
    """Records by (design, cell, base arm)."""

    out: dict[tuple[str, str, str], list[Mapping[str, Any]]] = {}
    for record in records:
        out.setdefault((str(record["design"]), str(record["cell"]), str(record["arm"])), []).append(record)
    return out


def _cells(settings: Mapping[str, Any]) -> list[tuple[str, str, dict[str, Any]]]:
    out = []
    for design in settings["design_index"]:
        for cell in design_cells(settings, design):
            kind = cell_kind(settings, design, cell)
            out.append((design, cell, {"design": design, "cell": cell, **{key: kind[key] for key in CELL_FIELDS[2:]},
                                       "target": kind["target"]}))
    return out


def _status(n_records: int, expected: int) -> str:
    if n_records == 0:
        return "not_run"
    return "complete" if n_records >= expected else "partial"


def rate_rows(records: Sequence[Mapping[str, Any]], settings: Mapping[str, Any]) -> list[dict[str, Any]]:
    """Per design, cell, arm and statistic: counts over completed records; failed and inconsistent left out."""

    grouped = _grouped(records)
    expected = expected_counts(settings)
    rule = settings["rejection"]
    rows = []
    for design, cell, fields in _cells(settings):
        for arm in rate_arms(settings):
            group = grouped.get((design, cell, _base_arm(arm)), [])
            applies = arm_applies(settings, design, _base_arm(arm))
            n_expected = expected[(design, cell)] if applies else 0
            status = _status(len(group), n_expected) if applies else "not_run"
            for statistic in STATISTICS:
                decisions = [rejects(_p(record, arm, statistic), rule) for record in group]
                rows.append({
                    **{key: fields[key] for key in CELL_FIELDS}, "arm": arm, "statistic": statistic,
                    "role": _role(fields, statistic), "rule": rule, "status": status, "expected": n_expected,
                    "records": len(group), **ssc._rate_row(decisions),
                    "failures": sum(1 for r in group if r.get("status") == "failed"),
                    "inconsistent": sum(1 for r in group if r.get("status") == "inconsistent"),
                })
    return rows


def paired_row(pairs: Sequence[tuple[bool | None, bool | None]]) -> dict[str, Any]:
    """An arm against its reference on the units where both are available.

    ``b`` counts units where only the reference rejects (the arm's loss), ``c`` where only the arm rejects (its gain);
    the exact one-sided p-values are ``production_latent_space_check.paired_p_values``.
    """

    usable = [(bool(ref), bool(arm)) for ref, arm in pairs if ref is not None and arm is not None]
    n = len(usable)
    b = sum(1 for ref, arm in usable if ref and not arm)
    c = sum(1 for ref, arm in usable if arm and not ref)
    return {
        "n_paired": n,
        "rate_reference": sum(ref for ref, _ in usable) / n if n else None,
        "rate_arm": sum(arm for _, arm in usable) / n if n else None,
        "difference": (c - b) / n if n else None,
        "b": b,
        "c": c,
        **plsc.paired_p_values(b, c),
    }


def paired_rows(records: Sequence[Mapping[str, Any]], settings: Mapping[str, Any],
                rates: Sequence[Mapping[str, Any]]) -> list[dict[str, Any]]:
    by_unit: dict[tuple[str, str, int, str], Mapping[str, Any]] = {
        (str(r["design"]), str(r["cell"]), int(r["replicate"]), str(r["arm"])): r for r in records
    }
    status = {(str(r["design"]), str(r["cell"]), str(r["arm"])): str(r["status"]) for r in rates}
    expected = expected_counts(settings)
    rule = settings["rejection"]
    rows = []
    for design, cell, fields in _cells(settings):
        replicates = range(expected[(design, cell)])
        for statistic in STATISTICS:
            for arm, reference, comparison in pair_specs(settings):
                pairs = [(rejects(_p(by_unit.get((design, cell, r, _base_arm(reference))), reference, statistic), rule),
                          rejects(_p(by_unit.get((design, cell, r, _base_arm(arm))), arm, statistic), rule))
                         for r in replicates]
                rows.append({**{key: fields[key] for key in CELL_FIELDS}, "statistic": statistic,
                             "role": _role(fields, statistic), "arm": arm, "reference_arm": reference,
                             "comparison": comparison, "status": status.get((design, cell, arm), "not_run"),
                             **paired_row(pairs)})
    return rows


def _support(record: Mapping[str, Any], block: str) -> float | None:
    support = record.get("support")
    if not support:
        return None
    values = [float(component[block]) for component in support if block in component]
    return float(np.mean(values)) if values else None


def cell_rows(records: Sequence[Mapping[str, Any]], settings: Mapping[str, Any]) -> list[dict[str, Any]]:
    """Per design, cell and arm: the fit, captures, sparsity, statistics, residual dispersion and checks."""

    grouped = _grouped(records)
    rows = []
    for design, cell, fields in _cells(settings):
        for arm in ARMS:
            group = grouped.get((design, cell, arm), [])
            done = [r for r in group if r.get("status") == "completed"]
            consistency = [r.get("consistency") or {} for r in done]
            low = [c for c in consistency if c.get("low_rank")]
            dispersion = [r.get("stage_dispersion") or {} for r in done]
            rows.append({
                **{key: fields[key] for key in CELL_FIELDS}, "arm": arm, "completed": len(done),
                "failures": sum(1 for r in group if r.get("status") == "failed"),
                "inconsistent": sum(1 for r in group if r.get("status") == "inconsistent"),
                "rank_q10": _quantile((r.get("rank") for r in done), 0.1),
                "rank_median": _median(r.get("rank") for r in done),
                "rank_q90": _quantile((r.get("rank") for r in done), 0.9),
                "auroc_median": _median(r.get("cv_mean_auroc") for r in done),
                "measured_dimension_median": _median(r.get("measured_dimension") for r in done),
                "capture_pooled_stage_median": _median((r.get("span_capture") or {}).get("pooled_stage")
                                                       for r in done),
                "capture_group_difference_median": _median((r.get("span_capture") or {}).get("group_difference")
                                                           for r in done),
                "keep_fraction_q10": _quantile((r.get("keep_fraction") for r in done), 0.1),
                "keep_fraction_median": _median(r.get("keep_fraction") for r in done),
                "keep_fraction_q90": _quantile((r.get("keep_fraction") for r in done), 0.9),
                **{f"support_{block}_median": _median(_support(r, block) for r in done) for block in SUPPORT_BLOCKS},
                **{f"observed_{s}_median": _median((r.get("pair_statistics") or {}).get(s) for r in done)
                   for s in STATISTICS},
                **{f"population_{s}_median": _median((r.get("population") or {}).get(s) for r in done)
                   for s in STATISTICS},
                **{f"dispersion_stage_{k}_median": _median(d.get(f"stage_{k}") for d in dispersion) for k in range(4)},
                "dispersion_max_min_median": _median(d.get("max_min_ratio") for d in dispersion),
                "within_stage_same_statistics": _share((r.get("within_stage") or {}).get("same_statistics")
                                                       for r in done if r.get("within_stage")),
                "rank_equal": _share(c.get("rank_equal") for c in consistency),
                "auroc_equal": _share(c.get("auroc_equal") for c in consistency),
                "low_rank_units": len(low),
                "identity_at_low_rank": _share(c.get("identity_at_low_rank") for c in low),
                "seconds_median": _median((r.get("seconds") or {}).get("arm") for r in done),
            })
    return rows


def rank_rows(records: Sequence[Mapping[str, Any]], settings: Mapping[str, Any]) -> list[dict[str, Any]]:
    """Per design, cell and ``pls_cv`` rank: each arm's counts on the units of that stratum."""

    by_unit: dict[tuple[str, str, int], dict[str, Mapping[str, Any]]] = {}
    for record in records:
        by_unit.setdefault((str(record["design"]), str(record["cell"]), int(record["replicate"])), {})[
            str(record["arm"])] = record
    rule = settings["rejection"]
    strata: dict[tuple[str, str, int, str, str], list[bool | None]] = {}
    for (design, cell, _), arms in by_unit.items():
        production = arms.get(REFERENCE_ARM)
        if production is None or production.get("status") != "completed":
            continue
        rank = int(production["rank"])
        for arm in rate_arms(settings):
            record = arms.get(_base_arm(arm))
            for statistic in STATISTICS:
                strata.setdefault((design, cell, rank, arm, statistic), []).append(
                    rejects(_p(record, arm, statistic), rule))
    order = {design: position for position, design in enumerate(settings["design_index"])}
    rows = []
    for (design, cell, rank, arm, statistic), decisions in sorted(
            strata.items(), key=lambda item: (order[item[0][0]], item[0][1], item[0][2], item[0][3], item[0][4])):
        rate = ssc._rate_row(decisions)
        if rate["count"]:
            rows.append({"design": design, "cell": cell, "rank": rank, "arm": arm, "statistic": statistic,
                         "count": rate["count"], "rejections": rate["rejections"]})
    return rows


def accounting(records: Sequence[Mapping[str, Any]], settings: Mapping[str, Any], *, duplicates: int,
               selections: Sequence[Sequence[str]]) -> dict[str, Any]:
    expected = expected_counts(settings)
    selected = sorted({arm for selection in selections for arm in selection}, key=ARMS.index)
    per_arm: dict[str, dict[str, int]] = {}
    for arm in ARMS:
        per_arm[arm] = {
            "expected": sum(n for (design, _), n in expected.items() if arm_applies(settings, design, arm)),
            "records": 0, "completed": 0, "failed": 0, "inconsistent": 0,
        }
    failures: dict[str, int] = {}
    errors: dict[str, int] = {}
    seconds = 0.0
    arm_seconds: dict[tuple[str, str], list[float]] = {}
    design_seconds: dict[str, dict[str, list[float]]] = {}
    for record in records:
        arm = str(record["arm"])
        status = str(record.get("status"))
        per_arm[arm]["records"] += 1
        per_arm[arm][status if status in ("completed", "failed", "inconsistent") else "failed"] += 1
        if status != "completed":
            label = f"{record['design']}/{record['cell']}/{arm}"
            failures[label] = failures.get(label, 0) + 1
            key = str(record.get("error_type") or status)
            errors[key] = errors.get(key, 0) + 1
        timing = record.get("seconds") or {}
        seconds += float(timing.get("call_share") or 0.0)
        if timing.get("arm") is not None:
            arm_seconds.setdefault((str(record["design"]), arm), []).append(float(timing["arm"]))
        shares = design_seconds.setdefault(str(record["design"]), {"generate": [], "fit": [], "call": []})
        for stage in ("generate", "fit", "call"):
            if timing.get(stage) is not None:
                shares[stage].append(float(timing[stage]))
    for arm in ARMS:
        per_arm[arm]["missing"] = per_arm[arm]["expected"] - per_arm[arm]["records"] if arm in selected else 0
    return {
        "profile": settings["profile"],
        "settings_hash": lsc.settings_hash(settings),
        "expected_units": sum(expected.values()),
        "arm_selections": [list(selection) for selection in selections],
        "selected_arms": selected,
        "arms": per_arm,
        "missing_records": sum(entry["missing"] for entry in per_arm.values()),
        "failed_records": sum(entry["failed"] for entry in per_arm.values()),
        "inconsistent_records": sum(entry["inconsistent"] for entry in per_arm.values()),
        "failures_by_cell": failures,
        "failures_by_error": errors,
        "duplicate_records": duplicates,
        "recorded_core_hours": seconds / 3600.0,
        "arm_seconds_median": {f"{design}/{arm}": _median(values) for (design, arm), values in
                               sorted(arm_seconds.items())},
        "design_seconds_median": {design: {stage: _median(values) for stage, values in stages.items()}
                                  for design, stages in design_seconds.items()},
    }


def _committed_rows(parent: str) -> list[dict[str, Any]]:
    return ssc._read_csv(parent_dir(parent) / "rate_table.csv")


def reproduction(rates: Sequence[Mapping[str, Any]], settings: Mapping[str, Any],
                 committed: Mapping[str, Sequence[Mapping[str, Any]]] | None = None) -> dict[str, Any] | None:
    """Each planned design's arms against its committed rate table, exact on available and rejected counts.

    ``committed`` maps a parent to its rows (default: read from the bound files).
    """

    plan = settings.get("reproduction") or []
    if not plan:
        return None
    mine = {(str(r["design"]), str(r["cell"]), str(r["arm"]), str(r["statistic"])): r for r in rates}
    tables: dict[str, dict[tuple[str, str, str, str], Mapping[str, Any]]] = {}
    for entry in plan:
        parent = entry["parent"]
        if parent in tables:
            continue
        rows = (committed or {}).get(parent)
        rows = _committed_rows(parent) if rows is None else rows
        if entry["table"] == "adni":
            tables[parent] = {(str(r["point"]), str(r["cell"]), str(r["arm"]), str(r["statistic"])): r for r in rows}
        else:
            # The stage-span and surgery parents name each statistic's arm by the statistic itself.
            tables[parent] = {(str(r["design"]), str(r["cell"]), REFERENCE_ARM, str(r["arm"])): r for r in rows
                              if str(r["arm"]) in STATISTICS}
    differing, unresolved = [], []
    compared = 0
    for entry in plan:
        theirs = tables[entry["parent"]]
        for cell in design_cells(settings, entry["design"]):
            for arm in entry["arms"]:
                for statistic in STATISTICS:
                    key = (entry["design"], cell, arm, statistic)
                    row = theirs.get(key)
                    if row is None:
                        unresolved.append({"design": key[0], "cell": cell, "arm": arm, "statistic": statistic,
                                           "parent": entry["parent"]})
                        continue
                    compared += 1
                    wanted = (int(row["count"]), int(row["rejections"]))
                    got_row = mine.get(key)
                    got = None if got_row is None else (int(got_row["count"]), int(got_row["rejections"]))
                    if got != wanted:
                        differing.append({"design": key[0], "cell": cell, "arm": arm, "statistic": statistic,
                                          "committed": {"count": wanted[0], "rejections": wanted[1]},
                                          "reproduced": None if got is None else {"count": got[0],
                                                                                  "rejections": got[1]}})
    return {
        "profile": settings["profile"],
        "passed": compared > 0 and not differing and not unresolved,
        "rows_compared": compared,
        "differing": differing,
        "unresolved_rows": unresolved,
        "sources": {name: settings["parents"][name] for name in sorted({e["parent"] for e in plan})
                    if name in settings["parents"]},
    }


def merge(out_dir: Path, profile: str) -> dict[str, Any]:
    """Merge one profile's shards; write its tables, accounting and (reproduced profiles) reproduction."""

    directory = profile_dir(out_dir, profile)
    settings = lsc.load_settings(directory)
    check_bound_files(settings)
    digest = lsc.settings_hash(settings)
    raw: list[dict[str, Any]] = []
    for path in sorted((directory / "records").glob("shard_*.jsonl")):
        raw.extend(read_records([path]))
    stale = sorted({record_key(r) for r in raw if r.get("settings_hash") != digest})
    if stale:
        raise SystemExit(f"{len(stale)} record(s) were written under other settings, e.g. {stale[:3]}; refusing.")
    seen: dict[str, int] = {}
    by_key: dict[str, dict[str, Any]] = {}
    for record in raw:
        key = record_key(record)
        seen[key] = seen.get(key, 0) + 1
        by_key[key] = record
    order = [unit_arm_key(unit, arm) for unit in profile_units(settings) for arm in ARMS]
    records = [by_key[key] for key in order if key in by_key]
    with (directory / "merged.jsonl").open("w", encoding="utf-8") as handle:
        handle.write("".join(lsc._line(record) for record in records))
    summary = accounting(records, settings, duplicates=sum(count - 1 for count in seen.values()),
                         selections=arm_selections(directory))
    rates = rate_rows(records, settings)
    lsc._write_csv(directory / "rate_table.csv", rates, RATE_COLUMNS)
    lsc._write_csv(directory / "paired_table.csv", paired_rows(records, settings, rates), PAIRED_COLUMNS)
    lsc._write_csv(directory / "cell_table.csv", cell_rows(records, settings), CELL_COLUMNS)
    lsc._write_csv(directory / "rank_table.csv", rank_rows(records, settings), RANK_COLUMNS)
    lsc._write_json(directory / "accounting.json", summary)
    result = reproduction(rates, settings)
    if result is not None:
        lsc._write_json(directory / "reproduction.json", result)
    return {
        "profile": profile,
        **{key: summary[key] for key in ("missing_records", "failed_records", "inconsistent_records",
                                         "duplicate_records", "recorded_core_hours")},
        "reproduction_passed": None if result is None else result["passed"],
        "reproduction_differing": None if result is None else len(result["differing"]),
    }


# ── Decision (pure functions over the merged tables) ─────────────────────────


class Index:
    """The merged profiles' tables, indexed for the criteria."""

    def __init__(self, blocks: Mapping[str, Mapping[str, Any] | None]) -> None:
        self.blocks = {name: block for name, block in blocks.items() if block is not None}
        self.rates: dict[tuple[str, str, str, str], Mapping[str, Any]] = {}
        self.paired: dict[tuple[str, str, str, str, str], Mapping[str, Any]] = {}
        self.designs: dict[str, dict[str, Any]] = {}
        self.cells: dict[str, list[dict[str, Any]]] = {}
        for block in self.blocks.values():
            for row in block["rates"]:
                self.rates[(str(row["design"]), str(row["cell"]), str(row["arm"]), str(row["statistic"]))] = row
            for row in block["paired"]:
                self.paired[(str(row["design"]), str(row["cell"]), str(row["statistic"]), str(row["arm"]),
                             str(row["reference_arm"]))] = row
            settings = block["settings"]
            for design, meta in settings["design_index"].items():
                self.designs[design] = dict(meta)
                cells: dict[str, dict[str, Any]] = {}
                for row in block["rates"]:
                    if str(row["design"]) == design and str(row["cell"]) not in cells:
                        cells[str(row["cell"])] = {
                            "cell": str(row["cell"]), "kind": str(row["kind"]),
                            "trajectory_mode": str(row["trajectory_mode"]),
                            "effect_size": _finite(row.get("effect_size")),
                            "span_alignment": _finite(row.get("span_alignment")),
                        }
                self.cells[design] = list(cells.values())

    def design_status(self, design: str, arm: str) -> str:
        """``complete`` | ``partial`` | ``not_run`` for an arm at a design (``not_run`` when its profile is absent)."""

        if design not in self.designs:
            return "not_run"
        statuses = {str(self.rates[(design, cell["cell"], arm, s)]["status"]) for cell in self.cells[design]
                    for s in STATISTICS if (design, cell["cell"], arm, s) in self.rates}
        if not statuses or statuses == {"not_run"}:
            return "not_run"
        return "complete" if statuses == {"complete"} else "partial"

    def observation(self, design: str, cell: str, arm: str, statistic: str) -> dict[str, Any]:
        entry = {"profile": self.designs.get(design, {}).get("profile"), "design": design, "cell": cell, "arm": arm,
                 "statistic": statistic}
        row = self.rates.get((design, cell, arm, statistic))
        if row is None or int(row["count"] or 0) <= 0:
            return {**entry, "status": "unavailable"}
        return {**entry, "n_available": int(row["count"]), "n_rejected": int(row["rejections"])}

    def pair(self, design: str, cell: str, statistic: str, arm: str, reference: str = REFERENCE_ARM) -> dict[str, Any]:
        entry = {"profile": self.designs.get(design, {}).get("profile"), "design": design, "cell": cell, "arm": arm,
                 "reference_arm": reference, "statistic": statistic}
        row = self.paired.get((design, cell, statistic, arm, reference))
        if row is None or int(row["n_paired"] or 0) <= 0:
            return {**entry, "status": "unavailable"}
        b, c = int(row["b"]), int(row["c"])
        # Recomputed from the discordant counts: the table's p-values are rounded to six significant digits.
        return {**entry, "n_paired": int(row["n_paired"]), "rate_reference": _finite(row["rate_reference"]),
                "rate_arm": _finite(row["rate_arm"]), "difference": _finite(row["difference"]), "b": b, "c": c,
                **plsc.paired_p_values(b, c)}

    def cells_of(self, design: str, kind: str) -> list[dict[str, Any]]:
        return [cell for cell in self.cells.get(design, []) if cell["kind"] == kind]


def _judge_rates(observations: Sequence[Mapping[str, Any]]) -> dict[str, Any] | None:
    return grade.judge_rates(observations) if observations else None


def _judge_paired(cells: Sequence[Mapping[str, Any]], *, p_field: str, rule: str) -> dict[str, Any] | None:
    return psc.judge_paired(cells, p_field=p_field, rule=rule) if cells else None


def _holds(judged: Mapping[str, Any] | None) -> bool:
    return judged is None or bool(judged["holds"])


def family_criteria(index: Index, arm: str, family: str) -> dict[str, Any]:
    """S1, S2, S3 and G of one candidate in one family, every criterion recorded; the label, the first failure."""

    designs = FAMILIES[family]
    statuses = {design: index.design_status(design, arm) for design in designs}
    run = [design for design in designs if statuses[design] != "not_run"]
    if not run:
        return {"label": "not_run", "partial": False, "designs": statuses, "criteria": None}
    partial = any(status != "complete" for status in statuses.values())
    controls = [(d, cell["cell"]) for d in run for cell in index.cells_of(d, "control")]
    targets = [(d, cell["cell"], TARGETS[cell["trajectory_mode"]]) for d in run
               for cell in index.cells_of(d, "in_span")]
    s1 = grade.judge_rates([index.observation(d, c, arm, s) for d, c in controls for s in STATISTICS])
    s2 = psc.judge_paired([index.pair(d, c, t, arm) for d, c, t in targets], p_field="p_loss",
                          rule="count_and_ceiling_on_paired_loss")
    strong_designs = [d for d in run if float(index.designs[d]["delta"]) >= CORE_DELTA]
    off_target = [(d, c, s) for d, c, t in targets for s in STATISTICS if s != t]
    strong = _judge_rates([index.observation(d, c, arm, s) for d, c, s in off_target if d in strong_designs])
    band_rows = []
    for d, c, s in off_target:
        if d in strong_designs:
            continue
        row = index.pair(d, c, s, arm)
        production = index.observation(d, c, REFERENCE_ARM, s)
        row["production_rate"] = (None if production.get("status") == "unavailable"
                                  else production["n_rejected"] / production["n_available"])
        band_rows.append(row)
    band = _judge_paired(band_rows, p_field="p_gain", rule="count_and_ceiling_on_paired_excess")
    above = psc._above_bound([index.observation(d, c, arm, s) for d, c, s in off_target])
    s3 = {"holds": _holds(strong) and _holds(band), "families": {"strong": strong, "band": band},
          "strong_designs": strong_designs, "off_target_above_bound": len(above),
          "off_target_above_bound_readings": above}
    g = psc.judge_gain([index.pair(d, c, t, arm) for d, c, t in targets])
    criteria = {"S1": s1, "S2": s2, "S3": s3, "G": g}
    failing = [name for name in ("S1", "S2", "S3", "G") if not criteria[name]["holds"]]
    return {"label": _LABEL_OF[failing[0]] if failing else "qualifies", "partial": partial, "designs": statuses,
            "failing_criteria": failing, "criteria": criteria}


def production_size(index: Index) -> dict[str, Any]:
    """``pls_cv``'s own S1 per family (recorded; the candidates are judged against their own rates)."""

    out = {}
    for family, designs in FAMILIES.items():
        run = [d for d in designs if index.design_status(d, REFERENCE_ARM) != "not_run"]
        observations = [index.observation(d, cell["cell"], REFERENCE_ARM, s) for d in run
                        for cell in index.cells_of(d, "control") for s in STATISTICS]
        out[family] = None if not observations else grade.judge_rates(observations)
    return out


def dispersion_designs(index: Index) -> list[str]:
    block = index.blocks.get("dispersion")
    return [] if block is None else list(block["settings"]["design_index"])


def dispersion_criterion(index: Index, arms: Sequence[str] = ARMS) -> dict[str, Any]:
    """D: each arm's controls on the dispersion designs; a paired excess test against ``pls_cv`` if it fails."""

    designs = dispersion_designs(index)
    readings = {arm: [(d, cell["cell"], s) for d in designs if index.design_status(d, arm) != "not_run"
                      for cell in index.cells_of(d, "control") for s in STATISTICS] for arm in arms}
    own = {arm: _judge_rates([index.observation(d, c, arm, s) for d, c, s in readings[arm]]) for arm in arms}
    production = own.get(REFERENCE_ARM)
    production_holds = production is not None and bool(production["holds"])
    finding = None
    if production is not None and not production_holds:
        finding = {"finding": "production_dispersion_excess",
                   "failing": [row for row in production["observations"] if row.get("status") != "holds"]}
    labels: dict[str, Any] = {}
    for arm in arms:
        if not readings[arm]:
            labels[arm] = {"label": "not_run", "rule": None, "judged": None}
            continue
        if arm == REFERENCE_ARM or production_holds or production is None:
            judged = own[arm]
            rule = "count_and_ceiling"
        else:
            judged = psc.judge_paired([index.pair(d, c, s, arm) for d, c, s in readings[arm]], p_field="p_gain",
                                      rule="count_and_ceiling_on_paired_excess")
            rule = "count_and_ceiling_on_paired_excess"
        holds = bool(judged and judged["holds"])
        label = "holds" if holds else ("production_dispersion_excess" if arm == REFERENCE_ARM else "fails")
        labels[arm] = {"label": label, "rule": rule, "judged": judged, "own_rates": own[arm]}
    recorded = {arm: _judge_rates([index.observation(d, cell["cell"], arm, s) for d in designs
                                   for cell in index.cells_of(d, "control") for s in STATISTICS
                                   if index.design_status(d, arm) != "not_run"])
                for arm in WITHIN_STAGE_ARMS}
    return {"designs": designs, "m_expected": 60, "production_holds": production_holds if production else None,
            "production_finding": finding, "arms": labels,
            "recorded_within_stage": {arm: None if judged is None else {
                "label": "holds" if judged["holds"] else "exceeds", "m": judged["m"],
                "observed_exceedances": judged["observed_exceedances"], "n_above_ceiling": judged["n_above_ceiling"]}
                for arm, judged in recorded.items()}}


def incidental_cells(index: Index) -> list[tuple[str, str, str]]:
    """(design, cell, target) of every surgery target cell and every out-of-span rotation cell."""

    out = []
    for design in index.designs:
        for kind in ("surgery", "out_of_span"):
            out += [(design, cell["cell"], TARGETS[cell["trajectory_mode"]]) for cell in index.cells_of(design, kind)]
    return out


def incidental_labels(index: Index, arm: str) -> dict[str, Any]:
    """I: each incidental cell labelled against ``pls_cv`` by the univariate baseline's rule (0.025 / m, margin)."""

    cells = [(d, c, t) for d, c, t in incidental_cells(index) if index.design_status(d, arm) != "not_run"]
    m = len(cells)
    level = ubc.LABEL_ALPHA / m if m else None
    labels = []
    for design, cell, target in cells:
        row = index.pair(design, cell, target, arm)
        entry = {"design": design, "cell": cell, "statistic": target, "kind": next(
            c["kind"] for c in index.cells[design] if c["cell"] == cell)}
        if row.get("status") == "unavailable" or level is None:
            labels.append({**entry, "label": "unresolved", "reason": "unavailable"})
            continue
        difference = float(row["difference"])
        loses = row["p_loss"] <= level and -difference >= ubc.MARGIN - ubc._MARGIN_TOLERANCE
        gains = row["p_gain"] <= level and difference >= ubc.MARGIN - ubc._MARGIN_TOLERANCE
        label = "candidate_loses" if loses else ("candidate_gains" if gains else "unresolved")
        labels.append({**entry, "label": label, "rate_reference": row["rate_reference"], "rate_arm": row["rate_arm"],
                       "difference": difference, "b": row["b"], "c": row["c"], "p_loss": row["p_loss"],
                       "p_gain": row["p_gain"]})
    return {"m": m, "level": level, "margin": ubc.MARGIN, "labels": labels,
            "counts": {label: sum(1 for e in labels if e["label"] == label) for label in INCIDENTAL_LABELS}}


def arm_verdict(families: Mapping[str, Mapping[str, Any]], dispersion: Mapping[str, Any]) -> str:
    everywhere = all(entry["label"] == "qualifies" and not entry["partial"] for entry in families.values())
    if everywhere and dispersion.get("label") == "holds":
        return "qualifies_everywhere"
    if any(entry["label"] == "qualifies" for entry in families.values()):
        return "scoped_candidate"
    return "keep_production"


def arm_consequence(arm: str, verdict: str, families: Mapping[str, Mapping[str, Any]],
                    dispersion: Mapping[str, Any], incidental: Mapping[str, Any]) -> str:
    qualifying = [f for f, e in families.items() if e["label"] == "qualifies"]
    other = {f: e["label"] + (" (partial)" if e["partial"] else "") for f, e in families.items()
             if e["label"] != "qualifies"}
    loses = [f"{e['design']}/{e['cell']}" for e in incidental["labels"] if e["label"] == "candidate_loses"]
    gains = [f"{e['design']}/{e['cell']}" for e in incidental["labels"] if e["label"] == "candidate_gains"]
    text = {
        "qualifies_everywhere": f"{arm} qualifies in every family and holds size under stage dispersion.",
        "scoped_candidate": f"{arm} qualifies in {', '.join(qualifying)}" + (
            "" if not other else f"; elsewhere: {', '.join(f'{k} {v}' for k, v in other.items())}") + ".",
        "keep_production": f"{arm} qualifies in no family" + (
            "" if not other else f" ({', '.join(f'{k} {v}' for k, v in other.items())})") + ".",
    }[verdict]
    text += f" Stage dispersion: {dispersion.get('label')}."
    text += (f" Incidental routes (never decisive): loses at {', '.join(loses)}." if loses
             else " Incidental routes (never decisive): no cell lost by the margin.")
    if gains:
        text += f" Gains at {', '.join(gains)}."
    return text


def _gate(index: Index) -> dict[str, Any]:
    entries = {}
    failures = []
    for profile in REPRODUCED_PROFILES:
        block = index.blocks.get(profile)
        if block is None:
            entries[profile] = {"status": "not_run"}
            continue
        result = block.get("reproduction") or {}
        entries[profile] = {"status": "passed" if result.get("passed") else "failed",
                            "rows_compared": result.get("rows_compared")}
        if not result.get("passed"):
            failures.append({"profile": profile, "differing": result.get("differing", []),
                             "unresolved_rows": result.get("unresolved_rows", [])})
    return {"passed": not failures, "profiles": entries, "failures": failures}


def decision(blocks: Mapping[str, Mapping[str, Any] | None]) -> dict[str, Any]:
    """Apply rule version 1 to the merged profiles (``None`` for a profile not run)."""

    index = Index(blocks)
    gate = _gate(index)
    dispersion = dispersion_criterion(index)
    arms: dict[str, Any] = {}
    for arm in CANDIDATE_ARMS:
        families = {family: family_criteria(index, arm, family) for family in FAMILIES}
        incidental = incidental_labels(index, arm)
        verdict = arm_verdict(families, dispersion["arms"][arm])
        above = sum(int(entry["criteria"]["S3"]["off_target_above_bound"]) for entry in families.values()
                    if entry["criteria"] is not None)
        arms[arm] = {
            "verdict": verdict,
            "families": {family: {key: entry[key] for key in ("label", "partial", "designs")}
                         for family, entry in families.items()},
            "dispersion": dispersion["arms"][arm]["label"],
            "incidental": incidental["counts"],
            "off_target_above_bound": above,
            "consequence": arm_consequence(arm, verdict, families, dispersion["arms"][arm], incidental),
            "criteria": {family: entry.get("criteria") for family, entry in families.items()},
            "failing_criteria": {family: entry.get("failing_criteria") for family, entry in families.items()
                                 if entry.get("failing_criteria")},
            "incidental_labels": incidental,
        }
    everywhere = [arm for arm in CANDIDATE_ARMS if arms[arm]["verdict"] == "qualifies_everywhere"]
    scoped = [arm for arm in CANDIDATE_ARMS if arms[arm]["verdict"] == "scoped_candidate"]
    named, reason = None, None
    if not gate["passed"]:
        verdict = "reproduction_failed"
    elif everywhere:
        verdict = "propose_default_change"
        counts = {arm: arms[arm]["off_target_above_bound"] for arm in everywhere}
        named = min(everywhere, key=lambda arm: (counts[arm], CANDIDATE_ARMS.index(arm)))
        tied = [arm for arm in everywhere if counts[arm] == counts[named] and arm != named]
        reason = (f"{named} has the fewest off-target readings above the bound among the arms that qualify "
                  f"everywhere ({counts})" + (f"; tied with {tied}, preferred by the order {list(CANDIDATE_ARMS)}"
                                              if tied else "") + ".")
    elif scoped:
        verdict = "scoped_candidates"
    else:
        verdict = "keep_production"
    return {
        "rule_version": RULE_VERSION,
        "verdict": verdict,
        "named_arm": named,
        "preference_reason": reason,
        "scoped_arms": scoped if verdict == "scoped_candidates" else [],
        "consequence": CONSEQUENCES[verdict],
        "consequences": dict(CONSEQUENCES),
        "gate": gate,
        "profiles": {profile: ("merged" if profile in index.blocks else "not_run") for profile in PROFILES},
        "complete": {profile: block["accounting"]["missing_records"] == 0 for profile, block in index.blocks.items()},
        "settings_hashes": {profile: lsc.settings_hash(block["settings"]) for profile, block in index.blocks.items()},
        "families": {family: list(designs) for family, designs in FAMILIES.items()},
        "reference_arm": REFERENCE_ARM,
        "preference": list(CANDIDATE_ARMS),
        "rejection": dict(REJECTION_RULES),
        "arms": arms,
        "production": {"size": {family: None if judged is None else {
            key: judged[key] for key in ("holds", "m", "observed_exceedances", "allowed_exceedances",
                                         "n_above_ceiling", "n_unavailable")} for family, judged in
            production_size(index).items()}},
        "dispersion": dispersion,
        "descriptive": descriptive_tables(index, arms),
        "notes": [
            "S1 and the delta = 2.0 half of S3 hold each candidate to its own rates under the stage-span paper "
            "grade's count-plus-ceiling rule; S2, G and the half of S3 below delta = 2.0 compare it with pls_cv on "
            "the same units and draws (exact paired tests).",
            "weighted measures delta in the weighted metric; the criteria compare rejections, not statistics.",
            "The incidental-route labels and the within-stage readings never change a family label or a verdict.",
        ],
    }


# ── Descriptive tables (design D7) ───────────────────────────────────────────


def descriptive_tables(index: Index, arms: Mapping[str, Any]) -> dict[str, Any]:
    """The fit, the statistics, production's off-target rates below δ = 2.0, the dispersion readings, the strata."""

    fit, statistics, rates_out, dispersion_rows = [], [], [], []
    for profile, block in index.blocks.items():
        for row in block["cells"]:
            design = str(row["design"])
            meta = index.designs.get(design, {})
            entry = {"profile": profile, "design": design, "family": meta.get("family"), "cell": row["cell"],
                     "arm": row["arm"], "completed": int(row["completed"] or 0)}
            fit.append({**entry, **{key: _finite(row.get(key)) for key in (
                "rank_q10", "rank_median", "rank_q90", "auroc_median", "measured_dimension_median",
                "capture_pooled_stage_median", "capture_group_difference_median", "keep_fraction_median",
                *(f"support_{block_name}_median" for block_name in SUPPORT_BLOCKS))}})
            statistics.append({**entry, **{f"{kind}_{s}": _finite(row.get(f"{kind}_{s}_median"))
                                           for kind in ("observed", "population") for s in STATISTICS}})
            if profile == "dispersion":
                dispersion_rows.append({**entry, "pattern": meta.get("pattern"),
                                        "reference_design": meta.get("reference_design"),
                                        **{key: _finite(row.get(key)) for key in (
                                            *(f"dispersion_stage_{k}_median" for k in range(4)),
                                            "dispersion_max_min_median", "within_stage_same_statistics")}})
        for row in block["rates"]:
            if str(row["status"]) == "not_run":
                continue
            rates_out.append({"profile": profile, **{key: row[key] for key in (
                "design", "cell", "arm", "statistic", "role", "status", "count", "rejections")},
                "rate": _finite(row.get("rate"))})
    below = []
    for family, designs in FAMILIES.items():
        for design in designs:
            if design not in index.designs or float(index.designs[design]["delta"]) >= CORE_DELTA:
                continue
            for cell in index.cells_of(design, "in_span"):
                target = TARGETS[cell["trajectory_mode"]]
                for statistic in STATISTICS:
                    if statistic == target:
                        continue
                    observed = index.observation(design, cell["cell"], REFERENCE_ARM, statistic)
                    if observed.get("status") == "unavailable":
                        continue
                    below.append({"family": family, "design": design, "cell": cell["cell"], "statistic": statistic,
                                  "n_available": observed["n_available"], "n_rejected": observed["n_rejected"],
                                  "rate": observed["n_rejected"] / observed["n_available"]})
    return {"fit": fit, "statistics": statistics, "rates": rates_out, "production_off_target_below_core": below,
            "dispersion": dispersion_rows, "rank_strata": rank_strata(index, arms)}


def rank_strata(index: Index, arms: Mapping[str, Any]) -> list[dict[str, Any]]:
    """Counts by ``pls_cv``'s rank at every cell where an arm loses: an S2 loss or an incidental ``candidate_loses``."""

    losses: list[tuple[str, str, str, str]] = []
    for arm, entry in arms.items():
        for criteria in entry["criteria"].values():
            if criteria is None:
                continue
            for row in criteria["S2"]["observations"]:
                if row.get("status") in ("below_alpha", "below_ceiling"):
                    losses.append((arm, row["design"], row["cell"], row["statistic"]))
        for row in entry["incidental_labels"]["labels"]:
            if row["label"] == "candidate_loses":
                losses.append((arm, row["design"], row["cell"], row["statistic"]))
    rows = []
    for arm, design, cell, statistic in losses:
        profile = index.designs.get(design, {}).get("profile")
        block = index.blocks.get(str(profile))
        if block is None:
            continue
        strata: dict[int, dict[str, Any]] = {}
        for row in block["ranks"]:
            if (str(row["design"]), str(row["cell"]), str(row["statistic"])) != (design, cell, statistic):
                continue
            if str(row["arm"]) not in (arm, REFERENCE_ARM):
                continue
            stratum = strata.setdefault(int(row["rank"]), {})
            stratum[str(row["arm"])] = {"count": int(row["count"]), "rejections": int(row["rejections"])}
        rows.append({"arm": arm, "design": design, "cell": cell, "statistic": statistic,
                     "strata": [{"rank": rank, **values} for rank, values in sorted(strata.items())]})
    return rows


# ── decide ───────────────────────────────────────────────────────────────────


def load_profile_tables(out_dir: Path, profile: str) -> dict[str, Any] | None:
    """A merged profile's tables, ``None`` when it was never merged; refused when merged under other settings."""

    directory = profile_dir(out_dir, profile)
    path = directory / "accounting.json"
    if not (directory / "settings.json").exists() or not path.exists():
        return None
    settings = lsc.load_settings(directory)
    check_bound_files(settings)
    summary = json.loads(path.read_text(encoding="utf-8"))
    if summary.get("settings_hash") != lsc.settings_hash(settings):
        raise SystemExit(f"{path} was merged under other settings; run `merge --profile {profile}` again.")
    block: dict[str, Any] = {
        "settings": settings,
        "accounting": summary,
        "rates": ssc._read_csv(directory / "rate_table.csv"),
        "paired": ssc._read_csv(directory / "paired_table.csv"),
        "cells": ssc._read_csv(directory / "cell_table.csv"),
        "ranks": ssc._read_csv(directory / "rank_table.csv"),
    }
    if profile in REPRODUCED_PROFILES:
        block["reproduction"] = json.loads((directory / "reproduction.json").read_text(encoding="utf-8"))
    return block


def decide(out_dir: Path) -> dict[str, Any]:
    blocks = {profile: load_profile_tables(out_dir, profile) for profile in PROFILES}
    if not any(blocks.values()):
        raise SystemExit(f"No merged profile under {out_dir}; run `merge` first.")
    result = decision(blocks)
    lsc._write_json(Path(out_dir) / DECISION_NAME, result)
    return result


# ── CLI ──────────────────────────────────────────────────────────────────────


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = parser.add_subparsers(dest="command", required=True)
    run = sub.add_parser("run", help="Evaluate (a shard of) one profile's (unit, arm) keys")
    run.add_argument("--profile", choices=PROFILES, required=True)
    run.add_argument("--out-dir", type=Path, required=True)
    run.add_argument("--arms", type=str, required=True, help=f"'all' or a comma list of {list(ARMS)}")
    run.add_argument("--shard", type=str, default="0/1", help="i/N, round-robin over the units")
    run.add_argument("--limit", type=int, default=None, help="Evaluate at most this many pending units")
    run.add_argument("--workers", type=int, default=1)
    run.add_argument("--perms", type=int, default=DEFAULT_PERMUTATIONS)
    run.add_argument("--reps", type=int, default=None, help="Replicates of every non-control cell (profile default)")
    run.add_argument("--control-reps", type=int, default=None, help="Replicates of the controls (default 500)")
    run.add_argument("--designs", type=str, default=None, help="Comma-separated designs; default: all")
    merge_parser = sub.add_parser("merge", help="Merge one profile's shards")
    merge_parser.add_argument("--profile", choices=PROFILES, required=True)
    merge_parser.add_argument("--out-dir", type=Path, required=True)
    decide_parser = sub.add_parser("decide", help=f"Write {DECISION_NAME} from the merged profiles")
    decide_parser.add_argument("--out-dir", type=Path, required=True)
    return parser


def settings_from_args(args: argparse.Namespace) -> dict[str, Any]:
    chosen = None if args.designs is None else [name for name in args.designs.split(",") if name]
    return profile_settings(args.profile, replicates=args.reps, control_replicates=args.control_reps,
                            permutations=args.perms, designs=chosen)


def _summary(result: Mapping[str, Any]) -> dict[str, Any]:
    return {
        "verdict": result["verdict"],
        "named_arm": result["named_arm"],
        "gate": result["gate"]["profiles"],
        "arms": {arm: {"verdict": entry["verdict"], "dispersion": entry["dispersion"],
                       "families": {f: e["label"] + (" (partial)" if e["partial"] else "")
                                    for f, e in entry["families"].items()},
                       "incidental": entry["incidental"]} for arm, entry in result["arms"].items()},
        "production_dispersion": result["dispersion"]["production_holds"],
        "profiles": result["profiles"],
    }


def main(argv: Iterable[str] | None = None) -> int:
    args = build_parser().parse_args(None if argv is None else list(argv))
    if args.command == "run":
        index, count = lsc.parse_shard(args.shard)
        evaluated = run_shard(args.out_dir, settings_from_args(args), index, count, arms=parse_arms(args.arms),
                              limit=args.limit, workers=args.workers)
        print(f"evaluated {evaluated} {args.profile} unit(s)", file=sys.stderr)
        return 0
    if args.command == "merge":
        print(json.dumps(merge(args.out_dir, args.profile), indent=1))
        return 0
    print(json.dumps(_summary(decide(args.out_dir)), indent=1))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
