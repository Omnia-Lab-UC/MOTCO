#!/usr/bin/env python3
"""Operating characteristics at the ADNI case study's design (no R).

The check behind ``openspec/changes/adni-design-operating-characteristics/``
(spec ``adni-design-validation``). The ADNI serum-metabolomics case study
(``adni_sexdiff``) ran its primary ``motco analyze`` at n = 1,501 with six
covariate design columns, the rank floor of 2 and a held-out CV stage AUROC of
0.600, and rejected nothing; its reading rule 3 calls the ``angle`` and
``shape`` non-rejections uninformative because no power statement exists at
that design. This check measures size and power there.

**Construction** (by import of :mod:`case_study_type_i_check`, itself built on
the D10 construction of :mod:`covariate_type_i_check`): the semi-synthetic
generator at ADNI's sex × stage table ``[[237, 310, 125], [223, 448, 158]]``,
three stages, methylation + expression (367 + 131 features), joint magnitude,
one weakened per-omic stage effect ``--delta`` with ``p_dmp`` 0.1. Six
covariate design columns of ADNI's kinds, in model order:

- ``z``: D10's numeric covariate confounded with group × stage, injected along
  a fixed direction at ``--s-z`` (ADNI's ``age``);
- ``category``: D10's three-level covariate on the same linear predictor,
  **declared** as a factor (drop-first ``category[c1]``, ``category[c2]``) and
  not injected (ADNI's ``apoe4``);
- ``phase``: a binary covariate drawn per group × stage cell at the table's
  shares (ADNI's ADNI2GO proportions) and injected into the **expression block
  only** as a dense offset ``X_f += s_p·sd_f·d_f·phase`` (``d_f`` i.i.d.
  N(0, 1) drawn once from the base seed, ``sd_f`` the block's pooled SD as
  received), so the covariate-blind basis follows it as ADNI's follows the
  p180 run offset (ADNI's ``phase[ADNI2GO]``); entered as a numeric 0/1 column;
- ``x1``, ``x2``: i.i.d. N(0, 1), never injected (ADNI's ``education`` and
  ``statin``).

Each unit records the realized phase offset per block the way ADNI's export
did (median |d| with the pooled within-phase SD; |corr| of the standardized
block's PC1 with ``phase``).

Every seed (generator, ``z`` and ``category``, the injection directions, the
noise covariates, the phase draw and direction, RRPP) derives from base seed
**2032** and the replicate index, a family no committed check uses; every cell
of a replicate shares them and differs only in the surgery.

**Cells**: the zero-effect anchor, ``translation`` at e = 1.0, magnitude at
0.1, 0.25, 0.5 and 1.0, orientation and shape at 0.25, 0.5, 0.75 and 1.0
(14 cells). Each paper-grade **unit** (cell, replicate) is analysed once
through :func:`motco.analysis.analyze_trajectories` in the ``adjusted`` arm:
the six covariates in the trajectory model, ``min_components = 2``, net
displacement, the harness CV knobs, the settings' permutations and the
replicate's RRPP seed, no latent adjustment and no residual correction (the
case study's primary as issued). A ``TrajectoryAnalysisError`` is recorded as
the unit's error; the unit is then unavailable.

**Calibration** (``calibrate``; anchor units only) analyses each replicate in
the ``adjusted`` and the ``unadjusted`` arm (no covariates; same data, CV
settings and RRPP seed; the two latent spaces must be identical) and marks
three criteria: (1) the ``adjusted`` arm's mean held-out CV AUROC inside
``--auroc-band`` (default 0.55–0.65; ADNI 0.600); (2) its median first-
component |Spearman ρ| with ``phase`` at or above ``--rho-threshold``
(default 0.4; ADNI 0.542); (3) the ``unadjusted`` arm rejecting above α on at
least one statistic (the construction confounds). The knobs are ``--delta``,
``--s-z`` and ``--s-p``; every invocation appends to ``calibration.json``.

**Decision** (pure, :func:`adni_design_decision`; ``decide``):

- size: ``delta``, ``angle`` and ``shape`` at the anchor and the translation
  control, each against α + k·√(α(1 − α)/R) (k = 2, R the available
  replicates; 0.0695 at R = 500);
- power: per target pair (magnitude/``delta``, orientation/``angle``,
  shape/``shape``) the rate and SE at every effect beside the selected-rank
  distribution, and ``informative_from``: the smallest effect at which the
  rate reaches the 0.80 floor, or null. Off-target and rank-stratified rates
  (ranks 2–3, 4–7, 8+) are descriptive;
- verdict: ``diagnose_design_excess`` when a size reading exceeds (no power
  statement), else ``power_statement`` when a target statistic has an
  informative effect, else ``non_rejections_uninformative``; per statistic a
  measured ``statement`` for the case study; a run below 999 permutations (or
  with an unavailable size reading) is indicative.

Outputs (``--out-dir``):

- ``calibration.json`` (``calibrate``): one entry per invocation;
- ``settings.json`` (``run``): settings and their hash; a resume under other
  settings is refused;
- ``records/shard_<i>_of_<N>.jsonl``: one line per unit; resume skips units
  recorded in any shard file and re-runs torn lines;
- after ``merge`` (automatic for an unsharded run): ``adni_design_pvalues.csv``,
  ``adni_design_replicates.csv``, ``adni_design_rates.csv`` and
  ``adni_design_rank_strata.csv``;
- ``decide``: ``adni_design_decision.json`` and ``size_readings.csv``.

Examples:

    # calibration (anchor units, two arms; appends to calibration.json)
    python scripts/adni_design_check.py calibrate --reps 20 --perms 99 --delta 0.03 --s-z 3 --s-p 1.5 \\
        --out-dir results/adni-design/calibration --workers 8

    # smoke profile at the calibrated strengths, the defaults (99 permutations: indicative)
    python scripts/adni_design_check.py run --reps 20 --perms 99 --out-dir results/adni-design/smoke
    python scripts/adni_design_check.py decide --out-dir results/adni-design/smoke

    # SLURM: settings.json first, then the array (scripts/adni_design_array.sbatch), then merge and decide
    python scripts/adni_design_check.py run --shard 0/200 --limit 0 --out-dir results/adni-design/run
    sbatch --array=0-199 --export=ALL,ADNI_OUT=$PWD/results/adni-design/run,N_SHARDS=200 \\
        scripts/adni_design_array.sbatch
    python scripts/adni_design_check.py merge --out-dir results/adni-design/run
    python scripts/adni_design_check.py decide --out-dir results/adni-design/run
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import sys
from concurrent.futures import ProcessPoolExecutor
from dataclasses import asdict
from pathlib import Path
from time import perf_counter
from typing import Any, Iterable, Mapping, Sequence

import numpy as np
import pandas as pd

from motco.analysis import (
    OmicsBlock,
    TrajectoryAnalysisError,
    TrajectoryAnalysisParams,
    TrajectoryDataset,
    analyze_trajectories,
)

sys.path.insert(0, str(Path(__file__).resolve().parent))
import case_study_type_i_check as cs  # noqa: E402  (construction, settings, sharding, summaries)
import covariate_type_i_check as d10  # noqa: E402  (the D10 construction and its helpers)

sys.path.pop(0)

SCRIPT_VERSION = 1
CHECK = "adni_design"
STATISTICS: tuple[str, ...] = cs.STATISTICS
ANCHOR = cs.ANCHOR
CONTROL_MODES: tuple[str, ...] = cs.CONTROL_MODES
TRANSLATION_EFFECT = 1.0
MODES: tuple[str, ...] = ("magnitude", "orientation", "shape")
EFFECTS: dict[str, list[float]] = {
    "magnitude": [0.1, 0.25, 0.5, 1.0],
    "orientation": [0.25, 0.5, 0.75, 1.0],
    "shape": [0.25, 0.5, 0.75, 1.0],
}
#: What each mode's surgery changes, as the statements name it.
DIFFERENCE_NOUNS: dict[str, str] = {"magnitude": "size", "orientation": "orientation", "shape": "shape"}
#: Target pairs (mode, statistic) of the power rule.
TARGETS: tuple[tuple[str, str], ...] = (("magnitude", "delta"), ("orientation", "angle"), ("shape", "shape"))
ADJUSTED = "adjusted"
UNADJUSTED = "unadjusted"
ARM_COVARIATES: dict[str, bool] = {ADJUSTED: True, UNADJUSTED: False}
#: The paper-grade unit runs the primary only; calibration adds the descriptive unadjusted arm.
PAPER_ARMS: tuple[str, ...] = (ADJUSTED,)
CALIBRATION_ARMS: tuple[str, ...] = (ADJUSTED, UNADJUSTED)
BASE_SEED = 2032
DEFAULT_REPLICATES = 500
DEFAULT_PERMUTATIONS = 999
#: The calibrated strengths (calibration.json, 2026-11-02; the starting points were 0.03 / 3 / 1.5).
DEFAULT_DELTA = 0.02
DEFAULT_S_Z = 0.1
DEFAULT_S_P = 1.75
NOISE_COVARIATES = 2
CATEGORICAL = "category"
PHASE = "phase"
#: The block the phase offset lands on (131 features, standing for p180's 133).
PHASE_BLOCK = "expression"
DEFAULT_AUROC_BAND = (0.55, 0.65)
DEFAULT_RHO_THRESHOLD = 0.4
ALPHA = 0.05
K = 2.0
POWER_FLOOR = 0.80
#: Selected-rank buckets of the descriptive strata: (label, largest rank in the bucket; None = open).
RANK_STRATA: tuple[tuple[str, int | None], ...] = (("2-3", 3), ("4-7", 7), ("8+", None))
VERDICTS: tuple[str, ...] = ("power_statement", "non_rejections_uninformative", "diagnose_design_excess")
CONSEQUENCES: dict[str, str] = {
    "power_statement": (
        "adni_sexdiff replaces its reading rule 3 with the measured statements in its own change; MOTCO documents "
        "them beside the SEA-AD statement; the Phase 6 exit review (roadmap item 10) reads this verdict"
    ),
    "non_rejections_uninformative": (
        "adni_sexdiff keeps its reading rule 3, now citing the measured rates, in its own change; the Phase 6 exit "
        "review (roadmap item 10) reads it as a design without power"
    ),
    "diagnose_design_excess": (
        "no power statement is issued; the roadmap names the exceeding readings and a diagnosis of them as the next "
        "step before any statement is adopted"
    ),
}
_ADNI_SOURCE = "adni_sexdiff results/motco_analysis/result.json (MOTCO 0.10.0, 2026-10-27 to 31)"
#: Cell tables with the phase-like covariate's per-cell shares (rows = the generator's group order A, B).
TABLES: dict[str, dict[str, Any]] = {
    "adni": {
        "group_stage_sizes": [[237, 310, 125], [223, 448, 158]],
        "phase_counts": [[143, 193, 53], [125, 227, 76]],
        "phase_shares": [[143 / 237, 193 / 310, 53 / 125], [125 / 223, 227 / 448, 76 / 158]],
        "groups": ["Female", "Male"],
        "stages": ["CN", "MCI", "AD"],
        "source": {
            "table": _ADNI_SOURCE,
            "phase_shares": "adni_sexdiff results/motco_input/metadata.csv + covariates.csv (ADNI2GO per sex x stage)",
        },
        "sha256": {
            "result.json": "d4a5cbf90dda74763dea91b9f9be8df47439da0e5fd4c475a8546fe2a6f581f0",
            "metadata.csv": "9c445b7de5a50ee14c0771f2ebf57416120e46f6daa3703e48ff90308cecf221",
            "covariates.csv": "06999e827d4f6eb6987bfa9c63bb24f9dfd4c07dfac96ce0b55500e8cd3fccfd",
        },
        #: ADNI's primary, read 2026-11-02 (descriptive; nothing is judged against these).
        "reference": {
            "selected_lv": 8,
            "cv_mean_auroc": 0.5995,
            "cv_mean_auroc_range_cv_selected": [0.580, 0.616],
            "first_component_spearman_phase": 0.542,
            "pooled_relative_eigengap": 0.744,
            "angle_null_q95": 35.44,
            "phase_offset": {
                "p180": {"median_abs_d": 1.12, "abs_d_gt_1": 80, "n_features": 133, "abs_corr_pc1_phase": 0.73},
                "nmr": {"median_abs_d": 0.21, "abs_d_gt_1": 0, "n_features": 165, "abs_corr_pc1_phase": 0.13},
            },
            "design_columns": ["age", "education", "apoe4[1]", "apoe4[2]", "phase[ADNI2GO]", "statin[1]"],
            "phase_pooled_share": 0.5443,
            "p_values": {"delta": 0.1116, "angle": 0.6879, "shape": 0.7105},
        },
    },
}


# ── Settings ─────────────────────────────────────────────────────────────────


def build_settings(args: argparse.Namespace) -> dict[str, Any]:
    """The case-study construction at the table, the six ADNI-kind columns, the phase nuisance and the rules."""

    low, high = args.auroc_band
    parent_args = cs.build_parser().parse_args(
        [
            "run",
            "--out-dir",
            "unused",
            "--reps",
            str(int(args.reps)),
            "--perms",
            str(int(args.perms)),
            "--base-seed",
            str(int(args.base_seed)),
            "--alpha",
            repr(ALPHA),
            "--k",
            repr(K),
            "--delta",
            repr(float(args.delta)),
            "--s-z",
            repr(float(args.s_z)),
            "--noise-covariates",
            str(NOISE_COVARIATES),
            "--auroc-band",
            f"{float(low)!r},{float(high)!r}",
            "--rho-threshold",
            repr(float(args.rho_threshold)),
            "--translation-effect",
            repr(TRANSLATION_EFFECT),
        ]
    )
    settings = cs.build_settings(parent_args)
    settings["script_version"] = SCRIPT_VERSION
    settings["check"] = CHECK
    table = TABLES[args.table]
    settings["table"] = {"name": args.table, **table}
    settings["generator"]["group_stage_sizes"] = [list(row) for row in table["group_stage_sizes"]]
    covariates = settings["covariates"]
    covariates["categorical"]["role"] = (
        f"declared as a factor (drop-first {CATEGORICAL}[c1], {CATEGORICAL}[c2]); not injected (s_c = 0)"
    )
    covariates["phase"] = {
        "shares": [list(row) for row in table["phase_shares"]],
        "distribution": "Bernoulli per group x stage cell at the table's shares",
        "column": "numeric 0/1 (one design column; the LS-mean reference row sits at the pooled proportion)",
    }
    covariates["declared"] = declared_covariates(settings)
    covariates["categorical_declared"] = [CATEGORICAL]
    settings["injection"]["s_p"] = float(args.s_p)
    settings["injection"]["phase"] = {
        "block": PHASE_BLOCK,
        "direction": "dense, iid N(0, 1) per feature of the block, drawn once from the base seed",
        "effect": "X_f += s_p * sd_f * d_f * phase, after the z injection",
        "scale": "pooled feature SD of the block as received (after the z injection)",
    }
    settings["arms"] = {
        ADJUSTED: {"covariates": True, "role": "every unit (the case study's primary as issued)"},
        UNADJUSTED: {"covariates": False, "role": "calibration only (descriptive non-vacuity criterion)"},
    }
    settings["effects"] = {mode: list(EFFECTS[mode]) for mode in MODES}
    settings["calibration"] = {
        "auroc_band": [float(low), float(high)],
        "rho_threshold": float(args.rho_threshold),
        "alpha": ALPHA,
        "criteria": {
            "adjusted_auroc": "adjusted-arm mean held-out CV stage AUROC at the anchor inside auroc_band",
            "adjusted_rho_phase": (
                "adjusted-arm median |Spearman| of the first latent component with phase at the anchor "
                ">= rho_threshold"
            ),
            "unadjusted_rejects": (
                "unadjusted-arm anchor rejection rate (p <= alpha) above alpha on at least one statistic"
            ),
        },
    }
    settings["rule"] = {
        "alpha": ALPHA,
        "k": K,
        "size": (
            "anchor and translation, every statistic: rate (p <= alpha) <= alpha + k*sqrt(alpha(1-alpha)/R), "
            "R the available replicates; size holds when all six readings hold"
        ),
        "power_floor": POWER_FLOOR,
        "power": (
            "per target pair (magnitude/delta, orientation/angle, shape/shape): informative_from = the smallest "
            "effect whose rate >= power_floor, else null; off-target and rank-stratified rates are descriptive"
        ),
        "targets": [list(target) for target in TARGETS],
        "verdict": (
            "diagnose_design_excess when any size reading exceeds; else power_statement when a target statistic "
            "has an informative effect; else non_rejections_uninformative"
        ),
        "rank_strata": [label for label, _ in RANK_STRATA],
        "indicative": f"fewer than {DEFAULT_PERMUTATIONS} permutations, or an unavailable size reading",
    }
    return settings


settings_hash = cs.settings_hash
ensure_settings = cs.ensure_settings
load_settings = cs.load_settings
parse_shard = cs.parse_shard
shard_path = cs.shard_path
unit_key = cs.unit_key
anchor_bound = cs.anchor_bound


def _seed(*parts: Any) -> int:
    encoded = json.dumps([CHECK, *parts], separators=(",", ":")).encode("utf-8")
    return int(hashlib.sha256(encoded).hexdigest()[:8], 16) & 0x7FFFFFFF


def phase_seed(base_seed: int, replicate: int) -> int:
    return _seed(int(base_seed), "phase", int(replicate))


def declared_covariates(settings: Mapping[str, Any]) -> list[str]:
    """Every covariate the adjusted arm declares, in model order: ``z``, the factor, ``phase``, the noise."""

    return ["z", CATEGORICAL, PHASE, *cs.noise_covariate_names(settings)]


def strengths(settings: Mapping[str, Any]) -> dict[str, float]:
    return {**cs.strengths(settings), "s_p": float(settings["injection"]["s_p"])}


# ── Grid ─────────────────────────────────────────────────────────────────────


def cells(settings: Mapping[str, Any]) -> list[tuple[str, float]]:
    """(mode, effect) per cell: the anchor, the translation control, then each mode's effects."""

    out = [(ANCHOR, 0.0), ("translation", float(settings["controls"]["translation"]))]
    for mode in MODES:
        out.extend((mode, float(effect)) for effect in settings["effects"].get(mode, []))
    return out


def units(settings: Mapping[str, Any]) -> list[tuple[str, float, int]]:
    """Every (mode, effect, replicate), replicate-major; shards take every N-th."""

    return [
        (mode, effect, replicate)
        for replicate in range(int(settings["replicates"]))
        for mode, effect in cells(settings)
    ]


def headroom_report(settings: Mapping[str, Any], reference: Any = None) -> dict[str, Any]:
    """Expected pool headroom of every pool-limited cell at the table and δ."""

    from motco.simulations.semisynthetic import expected_surgery_headroom

    report: dict[str, Any] = {}
    for mode, effect in cells(settings):
        headroom = expected_surgery_headroom(cs.generator_params(mode, effect, 0, settings), reference=reference)
        if headroom is not None:
            report[json.dumps([mode, float(effect)])] = asdict(headroom)
    return report


def check_headroom(settings: Mapping[str, Any], reference: Any = None) -> None:
    """Refuse settings with a cell beyond the generator's surgery headroom."""

    over = [key for key, headroom in headroom_report(settings, reference).items() if not headroom["fits"]]
    if over:
        raise SystemExit(f"These cells exceed the pool headroom at the table and delta: {over}.")


# ── Construction ─────────────────────────────────────────────────────────────


def table_cells(group_stage_sizes: Sequence[Sequence[int]]) -> tuple[np.ndarray, np.ndarray]:
    """(group index, stage index) per sample of a table, rows then stages (the phase draw without a dataset)."""

    group, stage = [], []
    for g, row in enumerate(group_stage_sizes):
        for s, size in enumerate(row):
            group.extend([g] * int(size))
            stage.extend([s] * int(size))
    return np.asarray(group, dtype=int), np.asarray(stage, dtype=int)


def draw_phase(
    group_index: np.ndarray, stage_index: np.ndarray, shares: Sequence[Sequence[float]], rng: np.random.Generator
) -> np.ndarray:
    """``phase_i ~ Bernoulli(share[group_i][stage_i])`` as 0/1 floats."""

    probabilities = np.asarray(shares, dtype=float)[group_index, stage_index]
    return (rng.random(group_index.size) < probabilities).astype(float)


def phase_direction(n_features: int, settings: Mapping[str, Any]) -> np.ndarray:
    """The dense offset direction ``d`` (i.i.d. N(0, 1)); drawn from the base seed, shared by every replicate."""

    rng = np.random.default_rng(_seed(int(settings["base_seed"]), "phase_direction", int(n_features)))
    return rng.normal(size=int(n_features))


def phase_offset_summary(values: np.ndarray, phase: np.ndarray) -> dict[str, Any]:
    """ADNI's export measure: per-feature standardized phase difference and the block PC1's |corr| with phase."""

    a, b = values[phase == 1.0], values[phase == 0.0]
    d = (a.mean(axis=0) - b.mean(axis=0)) / np.sqrt((a.var(axis=0) + b.var(axis=0)) / 2)
    standardized = (values - values.mean(axis=0)) / values.std(axis=0)
    u = np.linalg.svd(standardized, full_matrices=False)[0]
    return {
        "median_abs_d": float(np.median(np.abs(d))),
        "abs_d_gt_1": int((np.abs(d) > 1).sum()),
        "n_features": int(values.shape[1]),
        "abs_corr_pc1_phase": float(abs(np.corrcoef(u[:, 0], phase)[0, 1])),
    }


def build_unit(
    mode: str, effect: float, replicate: int, settings: Mapping[str, Any], reference: Any = None
) -> dict[str, Any]:
    """The case-study unit at the cell's surgery, plus ``phase`` drawn per cell and injected into one block.

    ``cs.build_unit`` gives D10's replicate (``z`` and ``category`` drawn, ``z``
    injected) with the noise covariates; ``phase`` comes from this check's own
    replicate seed, and its offset lands on the expression block only.
    """

    built = cs.build_unit(mode, effect, replicate, settings, reference=reference)
    metadata = built["metadata"]
    labels = list(cs.generator_params(mode, effect, replicate, settings).group_labels)
    group_index = np.asarray([labels.index(str(label)) for label in metadata["group"]], dtype=int)
    stage_index = metadata["stage"].to_numpy(dtype=int)
    seed = phase_seed(int(settings["base_seed"]), replicate)
    phase = draw_phase(group_index, stage_index, settings["covariates"]["phase"]["shares"], np.random.default_rng(seed))
    frame = built["blocks"][PHASE_BLOCK]
    values = frame.to_numpy(dtype=float)
    sd = values.std(axis=0, ddof=1)
    offset = float(settings["injection"]["s_p"]) * sd * phase_direction(values.shape[1], settings)
    blocks = dict(built["blocks"])
    blocks[PHASE_BLOCK] = pd.DataFrame(values + np.outer(phase, offset), index=frame.index, columns=frame.columns)
    built.update(
        metadata=metadata.assign(**{PHASE: phase}),
        blocks=blocks,
        seeds={**built["seeds"], "phase": seed},
        phase_offset=offset,
        phase_realized={
            name: phase_offset_summary(blocks[name].to_numpy(dtype=float), phase) for name in settings["layers"]
        },
    )
    return built


def phase_cell_shares(metadata: pd.DataFrame, group_labels: Sequence[str]) -> list[list[float]]:
    """Realized ``phase`` share per group × stage cell (rows = ``group_labels``)."""

    stages = sorted(set(metadata["stage"].tolist()))
    return [
        [
            float(metadata.loc[(metadata["group"].astype(str) == label) & (metadata["stage"] == stage), PHASE].mean())
            for stage in stages
        ]
        for label in group_labels
    ]


# ── Unit evaluation ──────────────────────────────────────────────────────────


ARM_FIELDS: tuple[str, ...] = (
    "statistics",
    "p_values",
    "null_summary",
    "selected_lv",
    "component_selection",
    "cv_mean_auroc",
    "first_component_spearman_phase",
    "first_component_spearman_z",
    "pooled_relative_eigengap",
    "covariate_design_columns",
    "latent_digest",
    "error",
)


def analysis_params(settings: Mapping[str, Any], seed: int, *, n_jobs: int = 1) -> TrajectoryAnalysisParams:
    return TrajectoryAnalysisParams(
        permutations=int(settings["permutations"]),
        seed=seed,
        n_jobs=n_jobs,
        orientation_estimator=settings["orientation_estimator"],
        min_components=int(settings["min_components"]),
        **settings["cv"],
    )


def trajectory_dataset(built: Mapping[str, Any], settings: Mapping[str, Any], *, covariates: bool) -> TrajectoryDataset:
    metadata = built["metadata"]
    covariate_kwargs: dict[str, Any] = (
        {"covariate_cols": declared_covariates(settings), "categorical_covariates": [CATEGORICAL]} if covariates else {}
    )
    return TrajectoryDataset(
        blocks=[OmicsBlock(name, built["blocks"][name], "identity") for name in settings["layers"]],
        metadata=metadata,
        sample_col="sample_id",
        group_col="group",
        stage_col="stage",
        stage_order=[str(stage) for stage in sorted(set(metadata["stage"].tolist()))],
        **covariate_kwargs,
    )


def latent_digest(scores: np.ndarray) -> str:
    return hashlib.sha256(np.ascontiguousarray(scores, dtype=np.float64).tobytes()).hexdigest()


def _arm_record(result: Any, metadata: pd.DataFrame) -> dict[str, Any]:
    record = result.record
    latent = record["latent"]
    scores = result.latent_scores.to_numpy(dtype=float)
    return {
        "statistics": record["statistics"],
        "p_values": record["p_values"],
        "null_summary": record["null_summary"],
        "selected_lv": latent["selected_lv"],
        "component_selection": latent["component_selection"],
        "cv_mean_auroc": latent.get("cv_mean_auroc"),
        "first_component_spearman_phase": cs._spearman(scores[:, 0], metadata[PHASE].to_numpy(dtype=float)),
        "first_component_spearman_z": cs._spearman(scores[:, 0], metadata["z"].to_numpy(dtype=float)),
        "pooled_relative_eigengap": (record["config_spectrum"].get("pooled") or {}).get("relative_eigengap"),
        "covariate_design_columns": record["covariates"]["design_columns"] if "covariates" in record else None,
        "latent_digest": latent_digest(scores),
        "error": None,
    }


def _failed_arm(exc: Exception) -> dict[str, Any]:
    failed: dict[str, Any] = {field: None for field in ARM_FIELDS}
    failed.update({"statistics": {}, "p_values": {}, "null_summary": {}, "error": str(exc)})
    return failed


def analyze_arms(
    built: Mapping[str, Any], settings: Mapping[str, Any], arms: Sequence[str], *, n_jobs: int = 1
) -> dict[str, Any]:
    """Analyse one built unit in each arm: same data, CV settings and RRPP seed.

    The covariate-blind latent space is shared, so when both arms run their
    scores and ranks must be identical; a difference is a ``RuntimeError``. An
    arm the analysis rejects is recorded with its error.
    """

    params = analysis_params(settings, built["seeds"]["rrpp"], n_jobs=n_jobs)
    records: dict[str, Any] = {}
    results: dict[str, Any] = {}
    for arm in arms:
        try:
            result = analyze_trajectories(trajectory_dataset(built, settings, covariates=ARM_COVARIATES[arm]), params)
        except TrajectoryAnalysisError as exc:
            records[arm] = _failed_arm(exc)
            continue
        results[arm] = result
        records[arm] = _arm_record(result, built["metadata"])
    if len(results) > 1:
        first, *others = results
        if any(
            records[arm]["latent_digest"] != records[first]["latent_digest"]
            or records[arm]["selected_lv"] != records[first]["selected_lv"]
            for arm in others
        ):
            raise RuntimeError(f"the {list(results)} latent spaces differ; the unit is not recorded.")
    return records


def evaluate_unit(
    mode: str,
    effect: float,
    replicate: int,
    settings: Mapping[str, Any],
    *,
    arms: Sequence[str] = PAPER_ARMS,
    n_jobs: int = 1,
    reference: Any = None,
) -> dict[str, Any]:
    start = perf_counter()
    built = build_unit(mode, effect, replicate, settings, reference=reference)
    records = analyze_arms(built, settings, arms, n_jobs=n_jobs)
    return d10._jsonable(
        {
            "key": unit_key(mode, effect, replicate),
            "mode": mode,
            "effect": float(effect),
            "replicate": replicate,
            "settings_hash": settings_hash(settings),
            "seeds": built["seeds"],
            "group_labels": built["group_labels"],
            "phase_realized": built["phase_realized"],
            "phase_cell_shares": phase_cell_shares(built["metadata"], built["group_labels"]),
            "arms": records,
            "runtime_seconds": perf_counter() - start,
        }
    )


# ── Shards and resume ────────────────────────────────────────────────────────


def read_records(paths: Iterable[Path], arms: Sequence[str] = PAPER_ARMS) -> dict[str, dict[str, Any]]:
    """Records keyed by unit; torn lines (an interrupted write) and partial records are re-run."""

    records: dict[str, dict[str, Any]] = {}
    for path in paths:
        with path.open(encoding="utf-8") as handle:
            for line in handle:
                try:
                    record = json.loads(line)
                except json.JSONDecodeError:
                    continue
                if set(record.get("arms", {})) == set(arms):
                    records[str(record["key"])] = record
    return records


def run_shard(
    out_dir: Path, settings: Mapping[str, Any], index: int, count: int, *, n_jobs: int = 1, limit: int | None = None
) -> int:
    """Evaluate this shard's missing units; returns how many were evaluated.

    A unit recorded in any shard file under ``records/`` counts as done; this
    shard's file is first rewritten with its complete records, so a torn line
    is re-run instead of duplicated. The headroom of every cell is checked
    before the first unit.
    """

    path = shard_path(out_dir, index, count)
    path.parent.mkdir(parents=True, exist_ok=True)
    own = read_records([path]) if path.exists() else {}
    with path.open("w", encoding="utf-8") as handle:
        for key in sorted(own):
            handle.write(json.dumps(own[key], sort_keys=True, allow_nan=False) + "\n")
    others = [other for other in sorted(path.parent.glob("shard_*.jsonl")) if other != path]
    done = {**read_records(others), **own}
    pending = [
        unit
        for position, unit in enumerate(units(settings))
        if position % count == index and unit_key(*unit) not in done
    ]
    if limit is not None:
        pending = pending[:limit]
    from motco.simulations.reference import load_reference

    reference = load_reference()
    check_headroom(settings, reference)
    for position, unit in enumerate(pending, start=1):
        record = evaluate_unit(*unit, settings, n_jobs=n_jobs, reference=reference)
        with path.open("a", encoding="utf-8") as handle:
            handle.write(json.dumps(record, sort_keys=True, allow_nan=False) + "\n")
            handle.flush()
        print(
            f"[shard {index}/{count}] {position}/{len(pending)} {record['key']} {record['runtime_seconds']:.1f}s",
            file=sys.stderr,
            flush=True,
        )
    return len(pending)


# ── Tables ───────────────────────────────────────────────────────────────────


PVALUE_COLUMNS = ("mode", "effect", "replicate", "statistic", "p_value")
REPLICATE_COLUMNS = (
    "mode",
    "effect",
    "replicate",
    "selected_lv",
    "component_selection",
    "cv_mean_auroc",
    "first_component_spearman_phase",
    "first_component_spearman_z",
    "pooled_relative_eigengap",
    "angle_null_q95",
    "angle_null_sd",
    "phase_median_abs_d",
    "phase_abs_corr_pc1",
    "other_block_median_abs_d",
    "other_block_abs_corr_pc1",
    "runtime_seconds",
    "error",
)
RATE_COLUMNS = (
    "mode",
    "effect",
    "statistic",
    "role",
    "units",
    "failures",
    "replicates",
    "rejections",
    "rate",
    "monte_carlo_se",
    "selected_rank",
    "cv_mean_auroc_mean",
    "median_abs_rho_phase",
    "median_abs_rho_z",
)
STRATA_COLUMNS = (
    "mode",
    "effect",
    "statistic",
    "role",
    "stratum",
    "replicates",
    "rejections",
    "rate",
    "pooled_replicates",
    "pooled_rate",
)
SIZE_COLUMNS = ("mode", "effect", "statistic", "replicates", "rejections", "rate", "monte_carlo_se", "bound", "verdict")
_MODE_ORDER: dict[str, int] = {mode: position for position, mode in enumerate((*CONTROL_MODES, *MODES))}


def _cell_order(cell: tuple[str, float]) -> tuple[int, float]:
    return _MODE_ORDER.get(cell[0], len(_MODE_ORDER)), cell[1]


def _number(value: Any) -> float | None:
    """A finite float from a record or a CSV cell (``""`` and ``None`` are missing)."""

    if value is None or value == "":
        return None
    number = float(value)
    return number if math.isfinite(number) else None


def p_value_rows(records: Mapping[str, Mapping[str, Any]]) -> list[dict[str, Any]]:
    rows = []
    for key in sorted(records):
        record = records[key]
        p_values = record["arms"][ADJUSTED].get("p_values") or {}
        for statistic in STATISTICS:
            rows.append(
                {
                    "mode": record["mode"],
                    "effect": record["effect"],
                    "replicate": record["replicate"],
                    "statistic": statistic,
                    "p_value": p_values.get(statistic),
                }
            )
    return rows


def replicate_rows(records: Mapping[str, Mapping[str, Any]]) -> list[dict[str, Any]]:
    rows = []
    for key in sorted(records):
        record = records[key]
        entry = record["arms"][ADJUSTED]
        angle_null = (entry.get("null_summary") or {}).get("angle") or {}
        realized = record.get("phase_realized") or {}
        other = next((realized[name] for name in sorted(realized) if name != PHASE_BLOCK), {})
        rows.append(
            {
                "mode": record["mode"],
                "effect": record["effect"],
                "replicate": record["replicate"],
                "selected_lv": entry.get("selected_lv"),
                "component_selection": entry.get("component_selection"),
                "cv_mean_auroc": entry.get("cv_mean_auroc"),
                "first_component_spearman_phase": entry.get("first_component_spearman_phase"),
                "first_component_spearman_z": entry.get("first_component_spearman_z"),
                "pooled_relative_eigengap": entry.get("pooled_relative_eigengap"),
                "angle_null_q95": angle_null.get("q95"),
                "angle_null_sd": angle_null.get("sd"),
                "phase_median_abs_d": (realized.get(PHASE_BLOCK) or {}).get("median_abs_d"),
                "phase_abs_corr_pc1": (realized.get(PHASE_BLOCK) or {}).get("abs_corr_pc1_phase"),
                "other_block_median_abs_d": other.get("median_abs_d"),
                "other_block_abs_corr_pc1": other.get("abs_corr_pc1_phase"),
                "runtime_seconds": record.get("runtime_seconds"),
                "error": entry.get("error"),
            }
        )
    return rows


def rejection_tallies(rows: Iterable[Mapping[str, Any]], *, alpha: float) -> dict[tuple, list[int]]:
    """(mode, effect, statistic) → [available replicates, rejections (``p ≤ alpha``)]."""

    tallies: dict[tuple, list[int]] = {}
    for row in rows:
        p_value = _number(row.get("p_value"))
        if p_value is None:
            continue
        tally = tallies.setdefault((str(row["mode"]), float(row["effect"]), str(row["statistic"])), [0, 0])
        tally[0] += 1
        tally[1] += int(p_value <= alpha)
    return tallies


def role(mode: str, statistic: str) -> str:
    if mode in CONTROL_MODES:
        return "size"
    return "target" if (mode, statistic) in TARGETS else "off_target"


def _by_cell(replicate_rows_: Iterable[Mapping[str, Any]]) -> dict[tuple[str, float], list[Mapping[str, Any]]]:
    out: dict[tuple[str, float], list[Mapping[str, Any]]] = {}
    for row in replicate_rows_:
        out.setdefault((str(row["mode"]), float(row["effect"])), []).append(row)
    return out


def _ranks(rows: Iterable[Mapping[str, Any]]) -> list[int]:
    return [int(number) for row in rows if (number := _number(row.get("selected_lv"))) is not None]


def rate_rows(
    rows: Iterable[Mapping[str, Any]], replicate_rows_: Iterable[Mapping[str, Any]], *, alpha: float
) -> list[dict[str, Any]]:
    """Per cell × statistic: R, rejections, rate, Monte Carlo SE, the rank distribution and the latent summaries."""

    tallies = rejection_tallies(rows, alpha=alpha)
    by_cell = _by_cell(replicate_rows_)
    out = []
    for mode, effect in sorted(by_cell, key=_cell_order):
        cell_rows = by_cell[(mode, effect)]
        auroc = cs.auroc_summary(_number(row.get("cv_mean_auroc")) for row in cell_rows)
        rho_phase = cs.correlation_summary(_number(row.get("first_component_spearman_phase")) for row in cell_rows)
        rho_z = cs.correlation_summary(_number(row.get("first_component_spearman_z")) for row in cell_rows)
        ranks = json.dumps(cs.rank_distribution(_ranks(cell_rows)), separators=(",", ":"))
        for statistic in STATISTICS:
            replicates, rejections, rate, se = cs._rate(tallies.get((mode, effect, statistic)))
            out.append(
                {
                    "mode": mode,
                    "effect": effect,
                    "statistic": statistic,
                    "role": role(mode, statistic),
                    "units": len(cell_rows),
                    "failures": sum(1 for row in cell_rows if row.get("error") not in (None, "")),
                    "replicates": replicates,
                    "rejections": rejections,
                    "rate": rate,
                    "monte_carlo_se": se,
                    "selected_rank": ranks,
                    "cv_mean_auroc_mean": auroc["mean"],
                    "median_abs_rho_phase": rho_phase["median_abs"],
                    "median_abs_rho_z": rho_z["median_abs"],
                }
            )
    return out


def rank_stratum(rank: int | None) -> str:
    """The bucket of a selected rank (``unavailable`` when the unit reports none)."""

    if rank is None:
        return "unavailable"
    for label, upper in RANK_STRATA:
        if upper is None or rank <= upper:
            return label
    raise AssertionError("the last stratum is open")


def rank_strata_rows(
    rows: Iterable[Mapping[str, Any]], replicate_rows_: Iterable[Mapping[str, Any]], *, alpha: float
) -> list[dict[str, Any]]:
    """Per cell × statistic × selected-rank bucket: R, rejections and rate, beside the pooled rate.

    Every available p-value falls in exactly one bucket, so the buckets'
    replicates sum to the pooled R.
    """

    ranks: dict[tuple[str, float, int], int | None] = {}
    for row in replicate_rows_:
        number = _number(row.get("selected_lv"))
        ranks[(str(row["mode"]), float(row["effect"]), int(float(row["replicate"])))] = (
            int(number) if number is not None else None
        )
    tallies: dict[tuple[str, float, str], dict[str, list[int]]] = {}
    for row in rows:
        p_value = _number(row.get("p_value"))
        if p_value is None:
            continue
        mode, effect, replicate = str(row["mode"]), float(row["effect"]), int(float(row["replicate"]))
        stratum = rank_stratum(ranks.get((mode, effect, replicate)))
        tally = tallies.setdefault((mode, effect, str(row["statistic"])), {}).setdefault(stratum, [0, 0])
        tally[0] += 1
        tally[1] += int(p_value <= alpha)
    out = []
    labels = [label for label, _ in RANK_STRATA]
    for mode, effect, statistic in sorted(tallies, key=lambda key: (*_cell_order(key[:2]), STATISTICS.index(key[2]))):
        strata = tallies[(mode, effect, statistic)]
        pooled_replicates = sum(tally[0] for tally in strata.values())
        pooled_rejections = sum(tally[1] for tally in strata.values())
        for stratum in labels + sorted(set(strata) - set(labels)):
            replicates, rejections = strata.get(stratum, [0, 0])
            out.append(
                {
                    "mode": mode,
                    "effect": effect,
                    "statistic": statistic,
                    "role": role(mode, statistic),
                    "stratum": stratum,
                    "replicates": replicates,
                    "rejections": rejections,
                    "rate": rejections / replicates if replicates else None,
                    "pooled_replicates": pooled_replicates,
                    "pooled_rate": pooled_rejections / pooled_replicates if pooled_replicates else None,
                }
            )
    return out


def merge(out_dir: Path, settings: Mapping[str, Any]) -> dict[str, Any]:
    """Merge every shard into the p-value, replicate, rate and rank-strata tables."""

    out_dir = Path(out_dir)
    records = read_records(sorted((out_dir / "records").glob("shard_*.jsonl")))
    stale = sorted({str(record["settings_hash"]) for record in records.values()} - {settings_hash(settings)})
    if stale:
        raise SystemExit(f"{out_dir}: records written under other settings ({stale}); refusing to merge.")
    expected = [unit_key(*unit) for unit in units(settings)]
    alpha = float(settings["alpha"])
    rows = p_value_rows(records)
    replicates = replicate_rows(records)
    rates = rate_rows(rows, replicates, alpha=alpha)
    d10._write_csv(out_dir / "adni_design_pvalues.csv", rows, PVALUE_COLUMNS)
    d10._write_csv(out_dir / "adni_design_replicates.csv", replicates, REPLICATE_COLUMNS)
    d10._write_csv(out_dir / "adni_design_rates.csv", rates, RATE_COLUMNS)
    strata = rank_strata_rows(rows, replicates, alpha=alpha)
    d10._write_csv(out_dir / "adni_design_rank_strata.csv", strata, STRATA_COLUMNS)
    return {
        "units": len(records),
        "expected": len(expected),
        "missing": [key for key in expected if key not in records][:20],
        "failures": sum(1 for row in replicates if row["error"] is not None),
        "rates": {f"{row['mode']}@{row['effect']:g}/{row['statistic']}": row["rate"] for row in rates},
        "runtime_seconds": cs._quantiles(
            sorted(float(record["runtime_seconds"]) for record in records.values() if "runtime_seconds" in record)
        ),
    }


# ── Calibration ──────────────────────────────────────────────────────────────


def _anchor(records: Iterable[Mapping[str, Any]]) -> list[Mapping[str, Any]]:
    return [record for record in records if str(record["mode"]) == ANCHOR]


def _arm_values(records: Iterable[Mapping[str, Any]], arm: str, field: str) -> list[Any]:
    return [record["arms"][arm].get(field) for record in records]


def _median(values: Iterable[Any]) -> float | None:
    finite = [number for value in values if (number := _number(value)) is not None]
    return float(np.median(finite)) if finite else None


def calibration_criteria(records: Iterable[Mapping[str, Any]], settings: Mapping[str, Any]) -> dict[str, Any]:
    """The three calibration criteria on the anchor records; ``met`` is their conjunction.

    (1) the ``adjusted`` arm's mean held-out CV AUROC inside the band; (2) its
    median first-component ``|ρ|`` with ``phase`` at or above the threshold;
    (3) the ``unadjusted`` arm rejecting above α on at least one statistic.
    Each is evaluated on its own.
    """

    anchor = _anchor(records)
    calibration = settings["calibration"]
    low, high = (float(value) for value in calibration["auroc_band"])
    threshold = float(calibration["rho_threshold"])
    alpha = float(calibration["alpha"])
    auroc = cs.auroc_summary(_arm_values(anchor, ADJUSTED, "cv_mean_auroc"))
    rho = cs.correlation_summary(_arm_values(anchor, ADJUSTED, "first_component_spearman_phase"))
    rates: dict[str, float | None] = {}
    replicates: dict[str, int] = {}
    for statistic in STATISTICS:
        p_values = [_number((record["arms"][UNADJUSTED].get("p_values") or {}).get(statistic)) for record in anchor]
        finite = [p for p in p_values if p is not None]
        replicates[statistic] = len(finite)
        rates[statistic] = sum(p <= alpha for p in finite) / len(finite) if finite else None
    criteria = {
        "adjusted_auroc": {
            "band": [low, high],
            **auroc,
            "met": auroc["mean"] is not None and low <= auroc["mean"] <= high,
        },
        "adjusted_rho_phase": {
            "threshold": threshold,
            **rho,
            "met": rho["median_abs"] is not None and rho["median_abs"] >= threshold,
        },
        "unadjusted_rejects": {
            "alpha": alpha,
            "replicates": replicates,
            "rates": rates,
            "met": any(rate is not None and rate > alpha for rate in rates.values()),
        },
    }
    return {"criteria": criteria, "met": all(entry["met"] for entry in criteria.values())}


def calibration_descriptives(records: Iterable[Mapping[str, Any]]) -> dict[str, Any]:
    """Placement descriptives at the anchor (none judged): rank, ``|ρ|`` with ``z``, offset, eigengap, null width."""

    anchor = _anchor(records)
    arms = [arm for arm in CALIBRATION_ARMS if anchor and arm in anchor[0]["arms"]]
    angle_q95 = [((r["arms"][ADJUSTED].get("null_summary") or {}).get("angle") or {}).get("q95") for r in anchor]
    realized = [record.get("phase_realized") or {} for record in anchor]
    blocks = sorted({name for entry in realized for name in entry})
    return {
        "selected_rank": {arm: cs.rank_distribution(_arm_values(anchor, arm, "selected_lv")) for arm in arms},
        "first_component_spearman_z": cs.correlation_summary(
            _arm_values(anchor, ADJUSTED, "first_component_spearman_z")
        ),
        "pooled_relative_eigengap": {
            "median": _median(_arm_values(anchor, ADJUSTED, "pooled_relative_eigengap")),
            "quantiles": cs._quantiles(
                sorted(
                    n
                    for v in _arm_values(anchor, ADJUSTED, "pooled_relative_eigengap")
                    if (n := _number(v)) is not None
                )
            ),
        },
        "angle_null_q95": {"median": _median(angle_q95)},
        "phase_offset": {
            name: {
                "median_abs_d": _median(entry.get(name, {}).get("median_abs_d") for entry in realized),
                "abs_corr_pc1_phase": _median(entry.get(name, {}).get("abs_corr_pc1_phase") for entry in realized),
                "abs_d_gt_1": _median(entry.get(name, {}).get("abs_d_gt_1") for entry in realized),
            }
            for name in blocks
        },
    }


def calibration_entry(records: Sequence[Mapping[str, Any]], settings: Mapping[str, Any]) -> dict[str, Any]:
    """One ``calibration.json`` entry: strengths, R, permutations, the criteria, ``met`` and the descriptives."""

    anchor = _anchor(records)
    arms = [arm for arm in CALIBRATION_ARMS if anchor and arm in anchor[0]["arms"]]
    return {
        "check": CHECK,
        "table": settings["table"]["name"],
        "strengths": strengths(settings),
        "base_seed": settings["base_seed"],
        "replicates": len(anchor),
        "permutations": settings["permutations"],
        "cv": settings["cv"],
        "min_components": settings["min_components"],
        "declared_covariates": declared_covariates(settings),
        **calibration_criteria(anchor, settings),
        "descriptives": calibration_descriptives(anchor),
        "adni_reference": settings["table"].get("reference"),
        "failures": {arm: sum(1 for r in anchor if r["arms"][arm].get("error") is not None) for arm in arms},
        "runtime_seconds": cs._quantiles(sorted(float(r["runtime_seconds"]) for r in anchor)),
        "per_replicate": [
            {
                "replicate": record["replicate"],
                "phase_realized": record.get("phase_realized"),
                **{
                    arm: {
                        field: record["arms"][arm].get(field)
                        for field in (
                            "selected_lv",
                            "cv_mean_auroc",
                            "first_component_spearman_phase",
                            "first_component_spearman_z",
                            "pooled_relative_eigengap",
                            "p_values",
                            "error",
                        )
                    }
                    for arm in arms
                },
            }
            for record in sorted(anchor, key=lambda record: int(record["replicate"]))
        ],
    }


record_calibration = cs.record_calibration


def _calibration_unit(replicate: int, settings: Mapping[str, Any]) -> dict[str, Any]:
    return evaluate_unit(ANCHOR, 0.0, replicate, settings, arms=CALIBRATION_ARMS)


def calibrate(out_dir: Path, settings: Mapping[str, Any], *, workers: int = 1) -> dict[str, Any]:
    """Evaluate the anchor units of every replicate in both arms and append the criteria to ``calibration.json``."""

    replicates = list(range(int(settings["replicates"])))
    if workers > 1:
        with ProcessPoolExecutor(max_workers=workers) as pool:
            records = list(pool.map(_calibration_unit, replicates, [settings] * len(replicates)))
    else:
        records = [_calibration_unit(replicate, settings) for replicate in replicates]
    entry = calibration_entry(records, settings)
    record_calibration(Path(out_dir), entry)
    return entry


# ── Decision (pure) ──────────────────────────────────────────────────────────


def _verdict(rate: float | None, bound: float | None) -> str:
    if rate is None or bound is None:
        return "unavailable"
    return "holds" if rate <= bound else "exceeds"


def size_readings(rows: Iterable[Mapping[str, Any]], settings: Mapping[str, Any]) -> list[dict[str, Any]]:
    """The six readings: anchor and translation × ``delta``, ``angle``, ``shape`` against the one-sided bound."""

    rule = settings["rule"]
    alpha, k = float(rule["alpha"]), float(rule["k"])
    tallies = rejection_tallies(rows, alpha=alpha)
    out = []
    for mode in CONTROL_MODES:
        effect = float(settings["controls"][mode])
        for statistic in STATISTICS:
            replicates, rejections, rate, se = cs._rate(tallies.get((mode, effect, statistic)))
            bound = anchor_bound(alpha, replicates, k) if replicates else None
            out.append(
                {
                    "mode": mode,
                    "effect": effect,
                    "statistic": statistic,
                    "replicates": replicates,
                    "rejections": rejections,
                    "rate": rate,
                    "monte_carlo_se": se,
                    "bound": bound,
                    "verdict": _verdict(rate, bound),
                }
            )
    return out


def _rank_text(ranks: Sequence[int]) -> str:
    if not ranks:
        return "no CV-selected rank recorded"
    low, high, median = min(ranks), max(ranks), float(np.median(ranks))
    span = f"{low}" if low == high else f"{low}–{high}"
    return f"CV-selected rank {span} (median {median:g})"


def power_curves(
    rows: Iterable[Mapping[str, Any]], replicate_rows_: Iterable[Mapping[str, Any]], settings: Mapping[str, Any]
) -> dict[str, dict[str, Any]]:
    """Per target statistic: rate, SE and selected ranks at every effect, ``informative_from`` and the maximum."""

    rule = settings["rule"]
    alpha, floor = float(rule["alpha"]), float(rule["power_floor"])
    tallies = rejection_tallies(rows, alpha=alpha)
    by_cell = _by_cell(replicate_rows_)
    out: dict[str, dict[str, Any]] = {}
    for mode, statistic in TARGETS:
        effects = []
        for effect in settings["effects"].get(mode, []):
            replicates, rejections, rate, se = cs._rate(tallies.get((mode, float(effect), statistic)))
            ranks = _ranks(by_cell.get((mode, float(effect)), []))
            effects.append(
                {
                    "effect": float(effect),
                    "replicates": replicates,
                    "rejections": rejections,
                    "rate": rate,
                    "monte_carlo_se": se,
                    "selected_rank": cs.rank_distribution(ranks),
                    "selected_rank_median": float(np.median(ranks)) if ranks else None,
                    "_ranks": ranks,
                }
            )
        available = [entry for entry in effects if entry["rate"] is not None]
        informative = [entry for entry in available if entry["rate"] >= floor]
        first = min(informative, key=lambda entry: entry["effect"]) if informative else None
        top = max(available, key=lambda entry: (entry["rate"], -entry["effect"])) if available else None
        out[statistic] = {
            "mode": mode,
            "statistic": statistic,
            "floor": floor,
            "effects": effects,
            "informative_from": first["effect"] if first else None,
            "floor_held_from_informative": (
                all(entry["rate"] >= floor for entry in available if entry["effect"] >= first["effect"])
                if first
                else None
            ),
            "max_rate": top["rate"] if top else None,
            "effect_at_max": top["effect"] if top else None,
            "construction_maximum": max((entry["effect"] for entry in effects), default=None),
            "_first": first,
            "_top": top,
        }
    return out


def statement(curve: Mapping[str, Any]) -> str:
    """The measured sentence the case study can adopt for one target statistic."""

    statistic, mode = curve["statistic"], curve["mode"]
    noun = DIFFERENCE_NOUNS[mode]
    article = "an" if noun[0] in "aeiou" else "a"
    first, top = curve["_first"], curve["_top"]
    if first is not None:
        return (
            f"At the ADNI design `{statistic}` rejects {article} {noun} difference of e = {first['effect']:g} (the "
            f"{mode} construction) with power {first['rate']:.3f} (floor {curve['floor']:.2f}; SE "
            f"{first['monte_carlo_se']:.3f}, R = {first['replicates']}), at {_rank_text(first['_ranks'])}; a "
            f"non-rejection argues against {noun} differences of e ≥ {first['effect']:g} at this design."
        )
    if top is None:
        return f"At the ADNI design `{statistic}` has no available {mode} replicates; no statement is possible."
    return (
        f"At the ADNI design `{statistic}` reaches at most {top['rate']:.3f} (SE {top['monte_carlo_se']:.3f}, "
        f"R = {top['replicates']}) at e = {top['effect']:g} of the {mode} construction, at "
        f"{_rank_text(top['_ranks'])}; a non-rejection is uninformative up to the construction maximum "
        f"(e = {curve['construction_maximum']:g}; floor {curve['floor']:.2f})."
    )


def size_sentence(readings: Sequence[Mapping[str, Any]]) -> str:
    def reading(row: Mapping[str, Any]) -> str:
        rate = "unavailable" if row["rate"] is None else f"{row['rate']:.3f}"
        bound = "" if row["bound"] is None else f" (bound {row['bound']:.4f}, R = {row['replicates']})"
        return f"{row['mode']} `{row['statistic']}` {rate}{bound}"

    exceeding = [row for row in readings if row["verdict"] == "exceeds"]
    unavailable = [row for row in readings if row["verdict"] == "unavailable"]
    listed = "; ".join(reading(row) for row in readings)
    if exceeding:
        return (
            f"Size exceeds at the ADNI design: {'; '.join(reading(row) for row in exceeding)} above the one-sided "
            f"bound. All six readings: {listed}."
        )
    if unavailable:
        return f"Size is not established at the ADNI design: {len(unavailable)} reading(s) unavailable. {listed}."
    return f"Size holds at the ADNI design: all six readings at or below the one-sided bound ({listed})."


def run_descriptives(replicate_rows_: Iterable[Mapping[str, Any]], settings: Mapping[str, Any]) -> dict[str, Any]:
    """The calibration placement realized in the run's anchor units (descriptive; the unadjusted arm never runs)."""

    anchor = [row for row in replicate_rows_ if str(row["mode"]) == ANCHOR]
    calibration = settings["calibration"]
    low, high = (float(value) for value in calibration["auroc_band"])
    threshold = float(calibration["rho_threshold"])
    auroc = cs.auroc_summary(_number(row.get("cv_mean_auroc")) for row in anchor)
    rho = cs.correlation_summary(_number(row.get("first_component_spearman_phase")) for row in anchor)
    return {
        "role": "descriptive; the calibration criteria realized at the run's anchor units",
        "adjusted_auroc": {
            **auroc,
            "band": [low, high],
            "in_band": auroc["mean"] is not None and low <= auroc["mean"] <= high,
        },
        "adjusted_rho_phase": {
            **rho,
            "threshold": threshold,
            "met": rho["median_abs"] is not None and rho["median_abs"] >= threshold,
        },
        "unadjusted_rejects": "not run in the paper-grade profile (calibration only)",
        "selected_rank": cs.rank_distribution(_ranks(anchor)),
        "first_component_spearman_z": cs.correlation_summary(
            _number(row.get("first_component_spearman_z")) for row in anchor
        ),
        "pooled_relative_eigengap_median": _median(row.get("pooled_relative_eigengap") for row in anchor),
        "angle_null_q95_median": _median(row.get("angle_null_q95") for row in anchor),
        "phase_offset": {
            "median_abs_d": _median(row.get("phase_median_abs_d") for row in anchor),
            "abs_corr_pc1_phase": _median(row.get("phase_abs_corr_pc1") for row in anchor),
            "other_block_median_abs_d": _median(row.get("other_block_median_abs_d") for row in anchor),
            "other_block_abs_corr_pc1_phase": _median(row.get("other_block_abs_corr_pc1") for row in anchor),
        },
    }


def _public(curve: Mapping[str, Any]) -> dict[str, Any]:
    """The curve without its private fields (``_first``, ``_top``, ``_ranks``)."""

    public = {key: value for key, value in curve.items() if not key.startswith("_")}
    public["effects"] = [{k: v for k, v in entry.items() if not k.startswith("_")} for entry in curve["effects"]]
    return public


def adni_design_decision(
    rows: Iterable[Mapping[str, Any]], replicate_rows_: Iterable[Mapping[str, Any]], settings: Mapping[str, Any]
) -> dict[str, Any]:
    """The verdict from the merged p-value rows and replicate rows (records or CSV rows).

    ``diagnose_design_excess`` when any of the six size readings exceeds its
    bound (no statement is issued); otherwise ``power_statement`` when a target
    statistic reaches the power floor at some effect, else
    ``non_rejections_uninformative``. Off-target and rank-stratified rates are
    descriptive and never enter.
    """

    rows = list(rows)
    replicate_rows_ = list(replicate_rows_)
    rule = settings["rule"]
    alpha = float(rule["alpha"])
    readings = size_readings(rows, settings)
    exceeding = [row for row in readings if row["verdict"] == "exceeds"]
    unavailable = [row for row in readings if row["verdict"] == "unavailable"]
    curves = power_curves(rows, replicate_rows_, settings)
    if exceeding:
        verdict = "diagnose_design_excess"
    elif any(curve["informative_from"] is not None for curve in curves.values()):
        verdict = "power_statement"
    else:
        verdict = "non_rejections_uninformative"
    statements = {} if exceeding else {statistic: statement(curve) for statistic, curve in curves.items()}
    permutations = int(settings["permutations"])
    reasons = {
        "below_default_permutations": permutations if permutations < DEFAULT_PERMUTATIONS else None,
        "unavailable_size_readings": [f"{row['mode']}/{row['statistic']}" for row in unavailable],
    }
    return {
        "rule": rule,
        "verdict": verdict,
        "consequence": CONSEQUENCES[verdict],
        "indicative": bool(reasons["below_default_permutations"] is not None or unavailable),
        "indicative_reasons": reasons,
        "size": {
            "holds": not exceeding and not unavailable,
            "sentence": size_sentence(readings),
            "readings": readings,
            "exceeding": [
                {key: row[key] for key in ("mode", "effect", "statistic", "replicates", "rate", "bound")}
                for row in exceeding
            ],
        },
        "power": {statistic: _public(curve) for statistic, curve in curves.items()},
        "statements": statements,
        "descriptive": {
            "role": "never enters the verdict",
            "cells": rate_rows(rows, replicate_rows_, alpha=alpha),
            "rank_strata": rank_strata_rows(rows, replicate_rows_, alpha=alpha),
        },
        "calibration_realized": run_descriptives(replicate_rows_, settings),
        "adni_reference": settings["table"].get("reference"),
        "construction": {
            "table": {key: settings["table"][key] for key in ("name", "group_stage_sizes") if key in settings["table"]},
            "strengths": strengths(settings),
            "generator": settings["generator"],
            "layers": settings["layers"],
            "cv": settings["cv"],
            "min_components": settings["min_components"],
            "orientation_estimator": settings["orientation_estimator"],
            "declared_covariates": declared_covariates(settings),
            "covariates": settings["covariates"],
            "injection": settings["injection"],
            "controls": settings["controls"],
            "effects": settings["effects"],
            "replicates": settings["replicates"],
            "permutations": permutations,
            "base_seed": settings["base_seed"],
        },
        "settings_hash": settings_hash(settings),
    }


def _read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle))


def decide(run_dir: Path, out_dir: Path | None = None) -> dict[str, Any]:
    """Read one merged run and write ``adni_design_decision.json`` and ``size_readings.csv``."""

    run_dir = Path(run_dir)
    out_dir = Path(out_dir) if out_dir is not None else run_dir
    settings = load_settings(run_dir)
    names = ("adni_design_pvalues.csv", "adni_design_replicates.csv")
    missing = [name for name in names if not (run_dir / name).exists()]
    if missing:
        raise SystemExit(f"{run_dir}: {', '.join(missing)} missing; run merge first.")
    rows = _read_csv(run_dir / names[0])
    replicates = _read_csv(run_dir / names[1])
    decision = adni_design_decision(rows, replicates, settings)
    decision["units"] = {
        "recorded": len(replicates),
        "expected": len(units(settings)),
        "failures": sum(1 for row in replicates if row.get("error")),
    }
    decision["run_directory"] = str(run_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    d10._write_csv(out_dir / "size_readings.csv", decision["size"]["readings"], SIZE_COLUMNS)
    (out_dir / "adni_design_decision.json").write_text(
        json.dumps(d10._jsonable(decision), indent=2, allow_nan=False) + "\n", encoding="utf-8"
    )
    return decision


# ── CLI ──────────────────────────────────────────────────────────────────────


def add_construction_arguments(parser: argparse.ArgumentParser) -> None:
    parser.add_argument("--out-dir", type=Path, required=True, help="Output directory")
    parser.add_argument("--table", choices=sorted(TABLES), default="adni", help="Cell table (default: adni)")
    parser.add_argument(
        "--reps", type=int, default=DEFAULT_REPLICATES, help=f"Replicates (default: {DEFAULT_REPLICATES})"
    )
    parser.add_argument(
        "--perms", type=int, default=DEFAULT_PERMUTATIONS, help=f"RRPP permutations (default: {DEFAULT_PERMUTATIONS})"
    )
    parser.add_argument("--base-seed", type=int, default=BASE_SEED, help=f"Seed family (default: {BASE_SEED})")
    parser.add_argument(
        "--delta",
        type=float,
        default=DEFAULT_DELTA,
        help=f"Per-omic stage effect on delta_methyl, delta_expr, delta_protein (default: {DEFAULT_DELTA})",
    )
    parser.add_argument(
        "--s-z", type=float, default=DEFAULT_S_Z, help=f"Numeric covariate injection strength (default: {DEFAULT_S_Z})"
    )
    parser.add_argument(
        "--s-p",
        type=float,
        default=DEFAULT_S_P,
        help=f"Phase offset strength on the expression block (default: {DEFAULT_S_P})",
    )
    parser.add_argument(
        "--auroc-band",
        type=cs.parse_band,
        default=DEFAULT_AUROC_BAND,
        help="Calibration band low,high for the adjusted arm's mean held-out CV AUROC (default: 0.55,0.65)",
    )
    parser.add_argument(
        "--rho-threshold",
        type=float,
        default=DEFAULT_RHO_THRESHOLD,
        help="Calibration threshold on the median first-component |Spearman| with phase (default: 0.4)",
    )


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Operating characteristics at the ADNI case study's design (no R).")
    sub = parser.add_subparsers(dest="command", required=True)
    cal = sub.add_parser("calibrate", help="Anchor units in the adjusted and unadjusted arms: the three criteria")
    add_construction_arguments(cal)
    cal.add_argument("--workers", type=int, default=1, help="Replicates evaluated in parallel")
    run = sub.add_parser("run", help="Evaluate (a shard of) the 14 cells' units; an unsharded run also merges")
    add_construction_arguments(run)
    run.add_argument("--shard", type=parse_shard, default=(0, 1), help="Run shard i of N (i/N; default 0/1)")
    run.add_argument("--n-jobs", type=int, default=1, help="Jobs for CV and RRPP (not part of the settings)")
    run.add_argument("--limit", type=int, default=None, help="Evaluate at most this many pending units")
    merge_parser = sub.add_parser("merge", help="Merge every shard into the tables (reads settings.json)")
    merge_parser.add_argument("--out-dir", type=Path, required=True, help="Run directory")
    decide_parser = sub.add_parser("decide", help="The verdict over one merged run")
    decide_parser.add_argument("--out-dir", type=Path, required=True, help="Where to write the decision")
    decide_parser.add_argument("--run-dir", type=Path, default=None, help="Merged run directory (default: --out-dir)")
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    if args.command == "decide":
        decision = decide(args.run_dir if args.run_dir is not None else args.out_dir, args.out_dir)
        summary = {key: decision[key] for key in ("verdict", "indicative", "indicative_reasons", "statements")}
        summary["size"] = decision["size"]["sentence"]
        summary["informative_from"] = {s: c["informative_from"] for s, c in decision["power"].items()}
        print(json.dumps(summary, indent=1, ensure_ascii=False))
        return 0
    if args.command == "merge":
        print(json.dumps(merge(args.out_dir, load_settings(args.out_dir)), indent=1))
        return 0
    settings = build_settings(args)
    if args.command == "calibrate":
        entry = calibrate(args.out_dir, settings, workers=args.workers)
        print(json.dumps({key: value for key, value in entry.items() if key != "per_replicate"}, indent=1))
        return 0
    ensure_settings(args.out_dir, settings)
    index, count = args.shard
    evaluated = run_shard(args.out_dir, settings, index, count, n_jobs=args.n_jobs, limit=args.limit)
    print(f"evaluated {evaluated} unit(s) into {shard_path(args.out_dir, index, count)}", file=sys.stderr)
    if count > 1 or args.limit is not None:
        return 0
    print(json.dumps(merge(args.out_dir, settings), indent=1))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
