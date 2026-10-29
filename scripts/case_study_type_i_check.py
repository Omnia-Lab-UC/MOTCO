#!/usr/bin/env python3
"""Case-study configuration Type I check at the Phase 6 SEA-AD design (no R).

The gate behind ``openspec/changes/case-study-configuration-type-i/`` (spec
``case-study-configuration-validation``): the size of the real-data analysis
(:func:`motco.analysis.analyze_trajectories`) at the SEA-AD case study's
configuration — the rank floor of 2, six numeric covariate design columns, a
weak stage signal and a nuisance covariate the covariate-blind basis follows.
It composes two validated constructions by import:

- the covariate Type I construction (D10, :mod:`covariate_type_i_check`): a
  null Phase 6 replicate (three stages, ``group_stage_sizes = [[11, 10, 28],
  [9, 10, 12]]``, methylation + expression, joint magnitude), a numeric
  covariate ``z = a_g·g + a_s·s + a_gs·g·s + ε`` confounded with group, stage
  and their interaction, and its effect injected into the features
  (``X_f += s_z·sd_f·w_f·z``) along directions drawn once from the base seed.
  The categorical covariate D10 draws stays in the RNG stream but is neither
  declared nor injected (``s_c = 0``);
- the weakened stage signal of the rank-floor check (:mod:`rank_floor_check`):
  one ``--delta`` applied to the three per-omic ``delta_*``, with ``p_dmp`` at
  the Phase 6 value 0.1 so the translation surgery keeps its headroom.

Beside ``z`` the check declares five pure-noise numeric covariates (``x1`` …
``x5``; i.i.d. N(0, 1) from their own seed, never injected), so the trajectory
model carries six numeric design columns and no factor, as the case study's
does.

Each **unit** is one (cell, replicate): the zero-effect anchor or the
``translation`` negative control at e = 1.0. Seeds depend on the replicate
index only, so both cells and every arm of replicate ``r`` share the baseline
draw, the covariate values and the permutation draws.

Each unit is analysed three times at ``min_components=2`` with the
net-displacement estimator, the harness CV knobs, the same ``random_state``,
permutation count and RRPP seed:

- ``adjusted``: the six covariates in the trajectory model only (the case
  study's primary analysis; judged first);
- ``unadjusted``: no covariates (descriptive — evidence the construction
  confounds; never judged);
- ``latent_stage``: the six covariates with ``covariate_latent_adjustment =
  "stage"`` (the case study's sensitivity arm; judged second, the fallback).

The ``adjusted`` and ``unadjusted`` latent spaces must be identical (the
covariate-blind basis); the unit fails loudly otherwise. Each arm records its
statistics, p-values, selected rank, component selection, mean held-out CV
stage AUROC, the Spearman correlation of its first latent component with
``z``, its null summary, pooled spectrum, covariate design columns, the
``latent_adjustment`` block (``latent_stage``) and any analysis error (the arm
is then unavailable for that unit; the other arms proceed).

**Calibration.** Before any paper-grade run, ``calibrate`` evaluates the
anchor units of R replicates in all three arms and marks three criteria:
(1) the ``adjusted`` arm's mean held-out CV AUROC lies inside ``--auroc-band``
(default 0.60–0.75; the SEA-AD cohorts read 0.636–0.697); (2) the
``unadjusted`` arm's median first-component ``|ρ|`` with ``z`` is at least
``--rho-threshold`` (default 0.8; the real basis reads 0.84); (3) the
``unadjusted`` arm rejects above α on at least one statistic (the construction
confounds, so the gate is not vacuous). The knobs are ``--delta`` and
``--s-z``; every invocation appends to ``calibration.json`` so the iterations
are recorded, and ``met`` is the conjunction. The decision file reports the
same criteria realized in the run itself.

The **decision** (pure, :func:`case_study_type_i_decision`): for the anchor and
the translation control, a judged arm's rejection rate (``p ≤ alpha``) over its
R available replicates must be at or below ``alpha + k·√(alpha(1 − alpha)/R)``
for ``delta``, ``angle`` and ``shape``. ``primary_validated`` when the
``adjusted`` arm holds every statistic at both controls; otherwise
``adopt_latent_stage`` when the ``latent_stage`` arm does; otherwise
``diagnose_configuration_excess`` naming every exceeding arm, control and
statistic. Over several stage-signal strengths the overall decision is the
conjunction (:func:`combine_decisions`; the ``combine`` subcommand): an arm is
admissible overall only when it is admissible in every run.

Outputs (``--out-dir``):

- ``settings.json``: settings and their hash; a resume into a directory
  written with different settings is refused;
- ``calibration.json``: one entry per ``calibrate`` invocation;
- ``records/shard_<i>_of_<N>.jsonl``: one line per unit (all arms); resume
  skips units recorded in any shard file;
- after ``merge`` (automatic for an unsharded run): ``merged.jsonl``,
  ``case_study_type_i_pvalues.csv``, ``case_study_type_i_replicates.csv``
  (per unit × arm: rank, CV AUROC, first-component correlation, ``angle``
  null width, residual df, error), ``case_study_type_i_rates.csv`` and
  ``case_study_type_i_decision.json``.

Examples:

    # calibration (anchor only, three arms; appends to calibration.json)
    python scripts/case_study_type_i_check.py calibrate --reps 20 --perms 99 --delta 0.1 --s-z 3 \\
        --out-dir results/case-study-type-i/calibration --workers 4

    # smoke profile at the settled strengths
    python scripts/case_study_type_i_check.py run --reps 20 --perms 99 --delta 0.1 --s-z 3 \\
        --out-dir /tmp/case-study-smoke

    # one SLURM array task, then the merge
    python scripts/case_study_type_i_check.py run --shard 3/50 --delta 0.1 --s-z 3 --out-dir results/cs
    python scripts/case_study_type_i_check.py merge --out-dir results/cs

    # the conjunction over two strengths
    python scripts/case_study_type_i_check.py combine \\
        --decision results/cs-a/case_study_type_i_decision.json \\
        --decision results/cs-b/case_study_type_i_decision.json \\
        --out results/case_study_type_i_overall_decision.json
"""

from __future__ import annotations

import argparse
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

sys.path.insert(0, str(Path(__file__).resolve().parent))
import covariate_type_i_check as d10  # noqa: E402  (the D10 construction, sharding and settings helpers)

sys.path.pop(0)

SCRIPT_VERSION = 1
CHECK = "case_study_type_i"
STATISTICS: tuple[str, ...] = d10.STATISTICS
ARMS: tuple[str, ...] = ("adjusted", "unadjusted", "latent_stage")
#: The arms the decision judges, in preference order; ``unadjusted`` is descriptive.
JUDGED_ARMS: tuple[str, ...] = ("adjusted", "latent_stage")
DESCRIPTIVE_ARM = "unadjusted"
ARM_COVARIATES: dict[str, bool] = {"adjusted": True, "unadjusted": False, "latent_stage": True}
ARM_LATENT_ADJUSTMENT: dict[str, str | None] = {"adjusted": None, "unadjusted": None, "latent_stage": "stage"}
MIN_COMPONENTS = 2
ANCHOR = "anchor"
CONTROL_MODES: tuple[str, ...] = (ANCHOR, "translation")
DEFAULT_TRANSLATION_EFFECT = 1.0
LAYERS: tuple[str, ...] = d10.LAYERS
#: The per-omic stage effects ``--delta`` sets; ``p_dmp`` stays at the Phase 6 value.
STRENGTH_FIELDS: tuple[str, ...] = ("delta_methyl", "delta_expr", "delta_protein")
P_DMP = 0.1
DEFAULT_DELTA = 0.1
DEFAULT_S_Z = 3.0
DEFAULT_NOISE_COVARIATES = 5
DEFAULT_AUROC_BAND = (0.60, 0.75)
DEFAULT_RHO_THRESHOLD = 0.8
#: D10's categorical covariate is drawn (keeping its RNG stream) but neither declared nor injected.
CATEGORICAL_LOGIT_SCALE = 0.5
S_C = 0.0
DECISIONS: tuple[str, ...] = ("primary_validated", "adopt_latent_stage", "diagnose_configuration_excess")
_QUANTILES: tuple[tuple[str, float], ...] = (
    ("min", 0.0),
    ("q05", 0.05),
    ("q25", 0.25),
    ("q50", 0.5),
    ("q75", 0.75),
    ("q95", 0.95),
    ("max", 1.0),
)


# ── Settings ─────────────────────────────────────────────────────────────────


def parse_band(text: str) -> tuple[float, float]:
    try:
        low, high = (float(value) for value in text.split(","))
    except ValueError as exc:
        raise argparse.ArgumentTypeError(f"--auroc-band must be low,high, got {text!r}.") from exc
    if not 0.0 <= low < high <= 1.0:
        raise argparse.ArgumentTypeError(f"--auroc-band must satisfy 0 <= low < high <= 1, got {text!r}.")
    return low, high


def build_settings(args: argparse.Namespace) -> dict[str, Any]:
    """The D10 settings at the check's fixed choices, plus the strengths, arms, controls and criteria."""

    d10_args = argparse.Namespace(
        reps=args.reps,
        perms=args.perms,
        base_seed=args.base_seed,
        alpha=args.alpha,
        k=args.k,
        a_g=args.a_g,
        a_s=args.a_s,
        a_gs=args.a_gs,
        categorical_logit_scale=CATEGORICAL_LOGIT_SCALE,
        s_z=args.s_z,
        s_c=S_C,
        feature_fraction=args.feature_fraction,
    )
    settings = d10.build_settings(d10_args)
    settings["script_version"] = SCRIPT_VERSION
    settings["check"] = CHECK
    generator = dict(d10.PHASE6_GENERATOR)
    for field in STRENGTH_FIELDS:
        generator[field] = float(args.delta)
    generator["p_dmp"] = P_DMP
    settings["generator"] = generator
    settings["min_components"] = MIN_COMPONENTS
    settings["covariates"]["categorical"]["role"] = "drawn (the D10 RNG stream is kept); neither declared nor injected"
    settings["covariates"]["noise_numeric"] = {
        "count": int(args.noise_covariates),
        "distribution": "iid N(0, 1), not injected",
    }
    settings["covariates"]["declared"] = declared_covariates(settings)
    settings["arms"] = {
        arm: {"covariates": ARM_COVARIATES[arm], "latent_adjustment": ARM_LATENT_ADJUSTMENT[arm]} for arm in ARMS
    }
    settings["controls"] = {ANCHOR: 0.0, "translation": float(args.translation_effect)}
    low, high = args.auroc_band
    settings["calibration"] = {
        "auroc_band": [float(low), float(high)],
        "rho_threshold": float(args.rho_threshold),
        "alpha": float(args.alpha),
        "criteria": {
            "adjusted_auroc": "adjusted-arm mean held-out CV stage AUROC at the anchor inside auroc_band",
            "unadjusted_rho": (
                "unadjusted-arm median |Spearman| of the first latent component with z at the anchor "
                ">= rho_threshold"
            ),
            "unadjusted_rejects": (
                "unadjusted-arm anchor rejection rate (p <= alpha) above alpha on at least one statistic"
            ),
        },
    }
    settings["rule"] = {
        "rate": "share of available replicates with p <= alpha",
        "judged_arms": list(JUDGED_ARMS),
        "controls": "anchor and translation, every statistic, rate <= alpha + k*sqrt(alpha(1-alpha)/R)",
        "decision": (
            "primary_validated when adjusted holds every control statistic; else adopt_latent_stage when "
            "latent_stage does; else diagnose_configuration_excess"
        ),
        "unadjusted": "descriptive (evidence the construction confounds); never judged",
    }
    return settings


settings_hash = d10.settings_hash
ensure_settings = d10.ensure_settings
load_settings = d10.load_settings
parse_shard = d10.parse_shard
shard_path = d10.shard_path
anchor_bound = d10.anchor_bound


def _seed(*parts: Any) -> int:
    encoded = json.dumps([CHECK, *parts], separators=(",", ":")).encode("utf-8")
    return int(hashlib.sha256(encoded).hexdigest()[:8], 16) & 0x7FFFFFFF


def noise_covariate_names(settings: Mapping[str, Any]) -> list[str]:
    return [f"x{index + 1}" for index in range(int(settings["covariates"]["noise_numeric"]["count"]))]


def declared_covariates(settings: Mapping[str, Any]) -> list[str]:
    """Every covariate the covariate arms declare, in model order: ``z`` then the noise covariates."""

    return ["z", *noise_covariate_names(settings)]


def strengths(settings: Mapping[str, Any]) -> dict[str, float]:
    generator = settings["generator"]
    return {
        "delta": float(generator[STRENGTH_FIELDS[0]]),
        **{field: float(generator[field]) for field in STRENGTH_FIELDS},
        "p_dmp": float(generator["p_dmp"]),
        "s_z": float(settings["injection"]["s_z"]),
    }


# ── Grid ─────────────────────────────────────────────────────────────────────


def cells(settings: Mapping[str, Any]) -> list[tuple[str, float]]:
    """(mode, effect) per cell: the anchor, then the translation control."""

    return [(ANCHOR, 0.0), ("translation", float(settings["controls"]["translation"]))]


def units(settings: Mapping[str, Any]) -> list[tuple[str, float, int]]:
    """Every (mode, effect, replicate), replicate-major; shards take every N-th."""

    return [
        (mode, effect, replicate)
        for replicate in range(int(settings["replicates"]))
        for mode, effect in cells(settings)
    ]


def unit_key(mode: str, effect: float, replicate: int) -> str:
    return json.dumps([mode, float(effect), int(replicate)])


# ── Construction ─────────────────────────────────────────────────────────────


def cell_settings(mode: str, effect: float, settings: Mapping[str, Any]) -> dict[str, Any]:
    """The settings with the cell's surgery folded into ``generator`` (what D10's builder reads)."""

    generator = dict(settings["generator"])
    generator["trajectory_mode"] = "magnitude" if mode == ANCHOR else mode
    generator["group_effect_size"] = float(effect)
    return {**settings, "generator": generator}


def generator_params(mode: str, effect: float, replicate: int, settings: Mapping[str, Any]):
    from motco.simulations.semisynthetic import SemiSyntheticTrajectoryParams

    generator = dict(cell_settings(mode, effect, settings)["generator"])
    generator["group_stage_sizes"] = tuple(tuple(row) for row in generator["group_stage_sizes"])
    seeds = d10.replicate_seeds(int(settings["base_seed"]), replicate)
    return SemiSyntheticTrajectoryParams(seed=seeds["generator"], **generator)


def headroom_report(settings: Mapping[str, Any], reference: Any = None) -> dict[str, Any]:
    """Expected pool headroom of the pool-limited control at the settings' strengths."""

    from motco.simulations.semisynthetic import expected_surgery_headroom

    report: dict[str, Any] = {}
    for mode, effect in cells(settings):
        headroom = expected_surgery_headroom(generator_params(mode, effect, 0, settings), reference=reference)
        if headroom is not None:
            report[json.dumps([mode, float(effect)])] = asdict(headroom)
    return report


def check_headroom(settings: Mapping[str, Any], reference: Any = None) -> None:
    """Refuse settings whose pool-limited control exceeds the generator's headroom."""

    over = [key for key, headroom in headroom_report(settings, reference).items() if not headroom["fits"]]
    if over:
        raise SystemExit(f"These cells exceed the pool headroom at the settings' strengths: {over}.")


def build_unit(
    mode: str, effect: float, replicate: int, settings: Mapping[str, Any], reference: Any = None
) -> dict[str, Any]:
    """The D10 replicate at the cell's surgery, plus the noise covariates (drawn from their own seed).

    D10's seeds depend on the base seed and the replicate only, so the anchor
    and the translation unit of one replicate share the generator seed, ``z``,
    the injection and the RRPP seed, and differ only in the surgery.
    """

    built = d10.build_replicate(replicate, cell_settings(mode, effect, settings), reference=reference)
    names = noise_covariate_names(settings)
    if names:
        seed = _seed(int(settings["base_seed"]), "noise_covariates", replicate)
        noise = np.random.default_rng(seed).normal(size=(built["metadata"].shape[0], len(names)))
        built["metadata"] = built["metadata"].assign(**{name: noise[:, j] for j, name in enumerate(names)})
        built["seeds"] = {**built["seeds"], "noise_covariates": seed}
    return built


# ── Unit evaluation ──────────────────────────────────────────────────────────


ARM_FIELDS: tuple[str, ...] = (
    "statistics",
    "p_values",
    "null_summary",
    "selected_lv",
    "component_selection",
    "cv_mean_auroc",
    "first_component_spearman_z",
    "pooled_spectrum",
    "covariate_design_columns",
    "latent_adjustment",
    "error",
)


def _spearman(a: np.ndarray, b: np.ndarray) -> float:
    from scipy.stats import spearmanr

    return float(spearmanr(a, b).statistic)


def _arm_record(result: Any, z: np.ndarray) -> dict[str, Any]:
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
        "first_component_spearman_z": _spearman(scores[:, 0], z),
        "pooled_spectrum": record["config_spectrum"].get("pooled"),
        "covariate_design_columns": record["covariates"]["design_columns"] if "covariates" in record else None,
        "latent_adjustment": record.get("latent_adjustment"),
        "error": None,
    }


def _failed_arm(exc: Exception) -> dict[str, Any]:
    failed: dict[str, Any] = {field: None for field in ARM_FIELDS}
    failed.update({"statistics": {}, "p_values": {}, "null_summary": {}, "error": str(exc)})
    return failed


def analyze_arms(built: Mapping[str, Any], settings: Mapping[str, Any], *, n_jobs: int = 1) -> dict[str, Any]:
    """Analyse one built unit in every arm: same blocks, CV settings and RRPP seed.

    The ``adjusted`` and ``unadjusted`` arms must share the latent space (the
    covariate-blind basis); a difference is a ``RuntimeError``. An arm the
    analysis rejects is recorded with its error and the other arms proceed.
    """

    from motco.analysis import (
        OmicsBlock,
        TrajectoryAnalysisError,
        TrajectoryAnalysisParams,
        TrajectoryDataset,
        analyze_trajectories,
    )

    metadata = built["metadata"]
    common: dict[str, Any] = {
        "blocks": [OmicsBlock(name, built["blocks"][name], "identity") for name in settings["layers"]],
        "metadata": metadata,
        "sample_col": "sample_id",
        "group_col": "group",
        "stage_col": "stage",
        "stage_order": [str(stage) for stage in sorted(set(metadata["stage"].tolist()))],
    }
    datasets = {
        True: TrajectoryDataset(**common, covariate_cols=declared_covariates(settings)),
        False: TrajectoryDataset(**common),
    }
    z = metadata["z"].to_numpy(dtype=float)
    arms: dict[str, Any] = {}
    results: dict[str, Any] = {}
    for arm in ARMS:
        params = TrajectoryAnalysisParams(
            permutations=int(settings["permutations"]),
            seed=built["seeds"]["rrpp"],
            n_jobs=n_jobs,
            orientation_estimator=settings["orientation_estimator"],
            covariate_latent_adjustment=ARM_LATENT_ADJUSTMENT[arm],
            min_components=int(settings["min_components"]),
            **settings["cv"],
        )
        try:
            result = analyze_trajectories(datasets[ARM_COVARIATES[arm]], params)
        except TrajectoryAnalysisError as exc:
            arms[arm] = _failed_arm(exc)
            continue
        results[arm] = result
        arms[arm] = _arm_record(result, z)
    if "adjusted" in results and "unadjusted" in results:
        adjusted, unadjusted = results["adjusted"], results["unadjusted"]
        if not (
            np.array_equal(adjusted.latent_scores.to_numpy(), unadjusted.latent_scores.to_numpy())
            and arms["adjusted"]["selected_lv"] == arms["unadjusted"]["selected_lv"]
        ):
            raise RuntimeError("the adjusted and unadjusted latent spaces differ.")
    return arms


def evaluate_unit(
    mode: str, effect: float, replicate: int, settings: Mapping[str, Any], *, n_jobs: int = 1, reference: Any = None
) -> dict[str, Any]:
    start = perf_counter()
    built = build_unit(mode, effect, replicate, settings, reference=reference)
    arms = analyze_arms(built, settings, n_jobs=n_jobs)
    return d10._jsonable(
        {
            "key": unit_key(mode, effect, replicate),
            "mode": mode,
            "effect": float(effect),
            "replicate": replicate,
            "settings_hash": settings_hash(settings),
            "seeds": built["seeds"],
            "group_labels": built["group_labels"],
            "covariate_cells": d10.covariate_cell_summary(built["metadata"], built["group_labels"]),
            "arms": arms,
            "runtime_seconds": perf_counter() - start,
        }
    )


# ── Summaries ────────────────────────────────────────────────────────────────


def _finite(values: Iterable[Any]) -> list[float]:
    return [float(value) for value in values if value is not None and math.isfinite(float(value))]


def _quantiles(values: Sequence[float]) -> dict[str, float]:
    return {name: float(np.quantile(values, q)) for name, q in _QUANTILES} if values else {}


def rank_distribution(values: Iterable[Any]) -> dict[str, int]:
    counts: dict[int, int] = {}
    for value in values:
        if value is not None:
            counts[int(value)] = counts.get(int(value), 0) + 1
    return {str(rank): counts[rank] for rank in sorted(counts)}


def auroc_summary(values: Iterable[Any]) -> dict[str, Any]:
    """Mean and quantiles of the held-out CV AUROC over the units that report one."""

    aurocs = _finite(values)
    return {
        "replicates": len(aurocs),
        "mean": float(np.mean(aurocs)) if aurocs else None,
        "quantiles": _quantiles(aurocs),
    }


def correlation_summary(values: Iterable[Any]) -> dict[str, Any]:
    """Distribution of ``|ρ|`` over the units that report one."""

    magnitudes = [abs(value) for value in _finite(values)]
    return {
        "replicates": len(magnitudes),
        "median_abs": float(np.median(magnitudes)) if magnitudes else None,
        "quantiles_abs": _quantiles(magnitudes),
    }


def _arm_values(records: Iterable[Mapping[str, Any]], arm: str, field: str) -> list[Any]:
    return [record["arms"][arm].get(field) for record in records]


# ── Calibration ──────────────────────────────────────────────────────────────


def calibration_criteria(records: Iterable[Mapping[str, Any]], settings: Mapping[str, Any]) -> dict[str, Any]:
    """The three calibration criteria on the anchor records; ``met`` is their conjunction.

    (1) the ``adjusted`` arm's mean held-out CV AUROC inside the band; (2) the
    ``unadjusted`` arm's median first-component ``|ρ|`` with ``z`` at or above
    the threshold; (3) the ``unadjusted`` arm rejecting above α on at least one
    statistic. Each is evaluated on its own.
    """

    anchor = [record for record in records if str(record["mode"]) == ANCHOR]
    calibration = settings["calibration"]
    low, high = (float(value) for value in calibration["auroc_band"])
    threshold = float(calibration["rho_threshold"])
    alpha = float(calibration.get("alpha", settings["alpha"]))
    auroc = auroc_summary(_arm_values(anchor, "adjusted", "cv_mean_auroc"))
    rho = correlation_summary(_arm_values(anchor, DESCRIPTIVE_ARM, "first_component_spearman_z"))
    tallies = rejection_tallies(p_value_rows({record["key"]: record for record in anchor}), alpha=alpha)
    rates: dict[str, float | None] = {}
    replicates: dict[str, int] = {}
    for statistic in STATISTICS:
        available, _, rate, _ = _rate(tallies.get((ANCHOR, 0.0, DESCRIPTIVE_ARM, statistic)))
        rates[statistic] = rate
        replicates[statistic] = available
    criteria = {
        "adjusted_auroc": {
            "band": [low, high],
            **auroc,
            "met": auroc["mean"] is not None and low <= auroc["mean"] <= high,
        },
        "unadjusted_rho": {
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


def calibration_entry(records: Sequence[Mapping[str, Any]], settings: Mapping[str, Any]) -> dict[str, Any]:
    """One ``calibration.json`` entry: strengths, R, permutations, the criteria and ``met``."""

    anchor = [record for record in records if str(record["mode"]) == ANCHOR]
    return {
        "strengths": strengths(settings),
        "base_seed": settings["base_seed"],
        "replicates": len(anchor),
        "permutations": settings["permutations"],
        "cv": settings["cv"],
        "min_components": settings["min_components"],
        "declared_covariates": declared_covariates(settings),
        **calibration_criteria(anchor, settings),
        "failures": {arm: sum(1 for r in anchor if r["arms"][arm].get("error") is not None) for arm in ARMS},
        "selected_rank": {arm: rank_distribution(_arm_values(anchor, arm, "selected_lv")) for arm in ARMS},
        "per_replicate": [
            {
                "replicate": record["replicate"],
                **{
                    arm: {
                        "selected_lv": record["arms"][arm].get("selected_lv"),
                        "cv_mean_auroc": record["arms"][arm].get("cv_mean_auroc"),
                        "first_component_spearman_z": record["arms"][arm].get("first_component_spearman_z"),
                        "p_values": record["arms"][arm].get("p_values"),
                        "error": record["arms"][arm].get("error"),
                    }
                    for arm in ARMS
                },
            }
            for record in sorted(anchor, key=lambda record: int(record["replicate"]))
        ],
    }


def record_calibration(out_dir: Path, entry: Mapping[str, Any]) -> list[dict[str, Any]]:
    """Append the entry to ``calibration.json`` (never overwritten) and return the history."""

    out_dir.mkdir(parents=True, exist_ok=True)
    path = out_dir / "calibration.json"
    history = json.loads(path.read_text(encoding="utf-8")) if path.exists() else []
    history.append(d10._jsonable(entry))
    path.write_text(json.dumps(history, indent=1, allow_nan=False) + "\n", encoding="utf-8")
    return history


def _calibration_unit(replicate: int, settings: Mapping[str, Any]) -> dict[str, Any]:
    return evaluate_unit(ANCHOR, 0.0, replicate, settings)


def calibrate(out_dir: Path, settings: Mapping[str, Any], *, workers: int = 1) -> dict[str, Any]:
    """Evaluate the anchor units of every replicate in all three arms and record the criteria."""

    replicates = list(range(int(settings["replicates"])))
    if workers > 1:
        with ProcessPoolExecutor(max_workers=workers) as pool:
            records = list(pool.map(_calibration_unit, replicates, [settings] * len(replicates)))
    else:
        records = [_calibration_unit(replicate, settings) for replicate in replicates]
    entry = calibration_entry(records, settings)
    record_calibration(out_dir, entry)
    return entry


# ── Shards, resume, merge ────────────────────────────────────────────────────


def read_records(paths: Iterable[Path]) -> dict[str, dict[str, Any]]:
    """Records keyed by unit; torn lines (an interrupted write) and partial records are re-run."""

    records: dict[str, dict[str, Any]] = {}
    for path in paths:
        with path.open(encoding="utf-8") as handle:
            for line in handle:
                try:
                    record = json.loads(line)
                except json.JSONDecodeError:
                    continue
                if set(record.get("arms", {})) == set(ARMS):
                    records[str(record["key"])] = record
    return records


def run_shard(
    out_dir: Path, settings: Mapping[str, Any], index: int, count: int, *, n_jobs: int = 1, limit: int | None = None
) -> int:
    """Evaluate this shard's missing units; returns how many were evaluated.

    A unit recorded in any shard file under ``records/`` counts as done, so a
    re-run under a different shard count skips recorded work; this shard's file
    is rewritten with its complete records first, so a torn line is re-run
    instead of duplicated.
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
        print(f"[shard {index}/{count}] {position}/{len(pending)} {record['key']}", file=sys.stderr, flush=True)
    return len(pending)


def p_value_rows(records: Mapping[str, Mapping[str, Any]]) -> list[dict[str, Any]]:
    rows = []
    for key in sorted(records):
        record = records[key]
        for arm in ARMS:
            p_values = record["arms"][arm].get("p_values") or {}
            for statistic in STATISTICS:
                rows.append(
                    {
                        "mode": record["mode"],
                        "effect": record["effect"],
                        "replicate": record["replicate"],
                        "arm": arm,
                        "statistic": statistic,
                        "p_value": p_values.get(statistic),
                    }
                )
    return rows


REPLICATE_COLUMNS = (
    "mode",
    "effect",
    "replicate",
    "arm",
    "selected_lv",
    "component_selection",
    "cv_mean_auroc",
    "first_component_spearman_z",
    "angle_null_q95",
    "angle_null_sd",
    "residual_df",
    "cv_min_training_df",
    "error",
)


def replicate_rows(records: Mapping[str, Mapping[str, Any]]) -> list[dict[str, Any]]:
    rows = []
    for key in sorted(records):
        record = records[key]
        for arm in ARMS:
            entry = record["arms"][arm]
            angle_null = (entry.get("null_summary") or {}).get("angle") or {}
            block = entry.get("latent_adjustment") or {}
            rows.append(
                {
                    "mode": record["mode"],
                    "effect": record["effect"],
                    "replicate": record["replicate"],
                    "arm": arm,
                    "selected_lv": entry.get("selected_lv"),
                    "component_selection": entry.get("component_selection"),
                    "cv_mean_auroc": entry.get("cv_mean_auroc"),
                    "first_component_spearman_z": entry.get("first_component_spearman_z"),
                    "angle_null_q95": angle_null.get("q95"),
                    "angle_null_sd": angle_null.get("sd"),
                    "residual_df": block.get("residual_df"),
                    "cv_min_training_df": block.get("cv_min_training_df"),
                    "error": entry.get("error"),
                }
            )
    return rows


# ── Rates and decision (pure) ────────────────────────────────────────────────


def rejection_tallies(rows: Iterable[Mapping[str, Any]], *, alpha: float) -> dict[tuple, list[int]]:
    """(mode, effect, arm, statistic) → [available replicates, rejections]."""

    tallies: dict[tuple, list[int]] = {}
    for row in rows:
        p_value = row.get("p_value")
        if p_value is None or not math.isfinite(float(p_value)):
            continue
        key = (str(row["mode"]), float(row["effect"]), str(row["arm"]), str(row["statistic"]))
        tally = tallies.setdefault(key, [0, 0])
        tally[0] += 1
        tally[1] += int(float(p_value) <= alpha)
    return tallies


def _rate(tally: Sequence[int] | None) -> tuple[int, int, float | None, float | None]:
    replicates, rejections = tally if tally else (0, 0)
    if not replicates:
        return 0, 0, None, None
    rate = rejections / replicates
    return replicates, rejections, rate, math.sqrt(rate * (1.0 - rate) / replicates)


def control_effects(rows: Iterable[Mapping[str, Any]]) -> dict[str, float | None]:
    """The anchor at 0 and the translation control at its largest effect in the rows (``None`` if absent)."""

    translation = [float(row["effect"]) for row in rows if str(row["mode"]) == "translation"]
    return {ANCHOR: 0.0, "translation": max(translation) if translation else None}


RATE_COLUMNS = (
    "mode",
    "effect",
    "arm",
    "statistic",
    "replicates",
    "rejections",
    "rate",
    "monte_carlo_se",
    "bound",
    "verdict",
)


def rejection_rates(rows: Iterable[Mapping[str, Any]], *, alpha: float, k: float) -> list[dict[str, Any]]:
    """Per control × arm × statistic: R available replicates, rejections (``p ≤ alpha``), rate and bound.

    The judged arms carry ``holds`` | ``exceeds`` | ``unavailable``; the
    ``unadjusted`` arm is ``descriptive``.
    """

    rows = list(rows)
    tallies = rejection_tallies(rows, alpha=alpha)
    effects = control_effects(rows)
    out = []
    for mode in CONTROL_MODES:
        effect = effects[mode]
        for arm in ARMS:
            for statistic in STATISTICS:
                tally = tallies.get((mode, effect, arm, statistic)) if effect is not None else None
                replicates, rejections, rate, se = _rate(tally)
                bound = anchor_bound(alpha, replicates, k) if replicates else None
                if arm not in JUDGED_ARMS:
                    verdict = "descriptive"
                elif rate is None or bound is None:
                    verdict = "unavailable"
                else:
                    verdict = "holds" if rate <= bound else "exceeds"
                out.append(
                    {
                        "mode": mode,
                        "effect": effect,
                        "arm": arm,
                        "statistic": statistic,
                        "replicates": replicates,
                        "rejections": rejections,
                        "rate": rate,
                        "monte_carlo_se": se,
                        "bound": bound,
                        "verdict": verdict,
                    }
                )
    return out


def case_study_type_i_decision(rows: Iterable[Mapping[str, Any]], *, alpha: float, k: float) -> dict[str, Any]:
    """The decision from the p-value rows (``mode``, ``effect``, ``replicate``, ``arm``, ``statistic``, ``p_value``).

    ``primary_validated`` when the ``adjusted`` arm holds the one-sided anchor
    bound for every statistic at both controls; otherwise ``adopt_latent_stage``
    when the ``latent_stage`` arm does; otherwise
    ``diagnose_configuration_excess``. The ``unadjusted`` arm never enters.
    """

    rows = list(rows)
    rates = rejection_rates(rows, alpha=alpha, k=k)
    effects = control_effects(rows)
    controls: dict[str, dict[str, Any]] = {mode: {"effect": effects[mode]} for mode in CONTROL_MODES}
    verdicts: dict[str, dict[str, dict[str, str]]] = {arm: {mode: {} for mode in CONTROL_MODES} for arm in JUDGED_ARMS}
    exceeding: list[dict[str, Any]] = []
    unavailable: list[dict[str, Any]] = []
    for row in rates:
        mode, arm, statistic = row["mode"], row["arm"], row["statistic"]
        controls[mode].setdefault(arm, {})[statistic] = {
            key: row[key] for key in ("replicates", "rejections", "rate", "monte_carlo_se", "bound", "verdict")
        }
        if arm not in JUDGED_ARMS:
            continue
        verdicts[arm][mode][statistic] = row["verdict"]
        if row["verdict"] == "exceeds":
            exceeding.append(
                {
                    "arm": arm,
                    "mode": mode,
                    "effect": row["effect"],
                    "statistic": statistic,
                    "rate": row["rate"],
                    "bound": row["bound"],
                }
            )
        elif row["verdict"] == "unavailable":
            unavailable.append({"arm": arm, "mode": mode, "statistic": statistic})
    admissible = [
        arm
        for arm in JUDGED_ARMS
        if all(verdict == "holds" for per_mode in verdicts[arm].values() for verdict in per_mode.values())
    ]
    if "adjusted" in admissible:
        decision = "primary_validated"
    elif "latent_stage" in admissible:
        decision = "adopt_latent_stage"
    else:
        decision = "diagnose_configuration_excess"
    return {
        "rule": (
            "judged-arm rate (p <= alpha) <= alpha + k * sqrt(alpha * (1 - alpha) / R) for delta, angle and "
            "shape at the zero-effect anchor and the translation control; primary_validated when the adjusted "
            "arm holds every statistic at both controls, else adopt_latent_stage when the latent_stage arm "
            "does, else diagnose_configuration_excess"
        ),
        "alpha": alpha,
        "k": k,
        "judged_arms": list(JUDGED_ARMS),
        "controls": controls,
        "verdicts": verdicts,
        "admissible": admissible,
        "exceeding": exceeding,
        "unavailable": unavailable,
        "decision": decision,
        "unadjusted_role": "descriptive (no covariates; evidence the construction confounds); not part of the decision",
    }


def combine_decisions(decisions: Sequence[Mapping[str, Any]]) -> dict[str, Any]:
    """The conjunction over runs: an arm is admissible overall only when it is admissible in every run."""

    if not decisions:
        raise ValueError("combine_decisions needs at least one decision.")
    admissible = [arm for arm in JUDGED_ARMS if all(arm in list(d.get("admissible", [])) for d in decisions)]
    if "adjusted" in admissible:
        decision = "primary_validated"
    elif "latent_stage" in admissible:
        decision = "adopt_latent_stage"
    else:
        decision = "diagnose_configuration_excess"
    return {
        "rule": "conjunction over runs: an arm is admissible overall only when it is admissible in every run",
        "runs": [
            {
                "label": d.get("label"),
                "strengths": (d.get("construction") or {}).get("strengths"),
                "settings_hash": d.get("settings_hash"),
                "decision": d.get("decision"),
                "admissible": list(d.get("admissible", [])),
                "exceeding": list(d.get("exceeding", [])),
            }
            for d in decisions
        ],
        "admissible": admissible,
        "decision": decision,
    }


def write_decision(
    out_dir: Path,
    rows: Sequence[Mapping[str, Any]],
    settings: Mapping[str, Any],
    records: Mapping[str, Mapping[str, Any]],
) -> dict[str, Any]:
    alpha, k = float(settings["alpha"]), float(settings["k"])
    decision = case_study_type_i_decision(rows, alpha=alpha, k=k)
    decision["calibration"] = {
        "statistics": settings["calibration"]["criteria"],
        **calibration_criteria(list(records.values()), settings),
    }
    decision["latent"] = {
        arm: {
            json.dumps([mode, float(effect)]): {
                "selected_rank": rank_distribution(_arm_values(cell_records, arm, "selected_lv")),
                "cv_mean_auroc": auroc_summary(_arm_values(cell_records, arm, "cv_mean_auroc")),
                "first_component_spearman_z": correlation_summary(
                    _arm_values(cell_records, arm, "first_component_spearman_z")
                ),
            }
            for mode, effect in cells(settings)
            for cell_records in [
                [r for r in records.values() if r["mode"] == mode and float(r["effect"]) == float(effect)]
            ]
        }
        for arm in ARMS
    }
    decision["failures"] = {
        arm: sum(1 for r in records.values() if r["arms"][arm].get("error") is not None) for arm in ARMS
    }
    decision["units"] = len(records)
    decision["settings_hash"] = settings_hash(settings)
    decision["construction"] = {
        "generator": settings["generator"],
        "strengths": strengths(settings),
        "layers": settings["layers"],
        "cv": settings["cv"],
        "min_components": settings["min_components"],
        "orientation_estimator": settings["orientation_estimator"],
        "covariates": settings["covariates"],
        "declared_covariates": declared_covariates(settings),
        "injection": settings["injection"],
        "arms": settings["arms"],
        "controls": settings["controls"],
        "permutations": settings["permutations"],
        "replicates": settings["replicates"],
        "base_seed": settings["base_seed"],
    }
    (out_dir / "case_study_type_i_decision.json").write_text(
        json.dumps(d10._jsonable(decision), indent=2, allow_nan=False) + "\n", encoding="utf-8"
    )
    d10._write_csv(out_dir / "case_study_type_i_rates.csv", rejection_rates(rows, alpha=alpha, k=k), RATE_COLUMNS)
    return decision


def merge(out_dir: Path, settings: Mapping[str, Any]) -> dict[str, Any]:
    """Merge every shard, then write the tables, the rates and the decision."""

    records = read_records(sorted((out_dir / "records").glob("shard_*.jsonl")))
    expected = [unit_key(*unit) for unit in units(settings)]
    with (out_dir / "merged.jsonl").open("w", encoding="utf-8") as handle:
        for key in sorted(records):
            handle.write(json.dumps(records[key], sort_keys=True, allow_nan=False) + "\n")
    rows = p_value_rows(records)
    d10._write_csv(
        out_dir / "case_study_type_i_pvalues.csv",
        rows,
        ("mode", "effect", "replicate", "arm", "statistic", "p_value"),
    )
    d10._write_csv(out_dir / "case_study_type_i_replicates.csv", replicate_rows(records), REPLICATE_COLUMNS)
    decision = write_decision(out_dir, rows, settings, records)
    return {
        "units": len(records),
        "expected": len(expected),
        "missing": [key for key in expected if key not in records][:20],
        "decision": decision["decision"],
        "admissible": decision["admissible"],
        "verdicts": decision["verdicts"],
        "exceeding": decision["exceeding"],
        "unavailable": decision["unavailable"],
        "controls": {
            mode: {
                arm: {s: decision["controls"][mode][arm][s]["rate"] for s in STATISTICS} for arm in ARMS
            }
            for mode in CONTROL_MODES
        },
        "bound": {mode: decision["controls"][mode]["adjusted"]["delta"]["bound"] for mode in CONTROL_MODES},
        "calibration": {
            "met": decision["calibration"]["met"],
            "adjusted_auroc_mean": decision["calibration"]["criteria"]["adjusted_auroc"]["mean"],
            "unadjusted_median_abs_rho": decision["calibration"]["criteria"]["unadjusted_rho"]["median_abs"],
            "unadjusted_anchor_rates": decision["calibration"]["criteria"]["unadjusted_rejects"]["rates"],
        },
        "failures": decision["failures"],
    }


# ── CLI ──────────────────────────────────────────────────────────────────────


def add_construction_arguments(parser: argparse.ArgumentParser) -> None:
    parser.add_argument("--out-dir", type=Path, required=True, help="Output directory")
    parser.add_argument("--reps", type=int, default=500, help="Replicates (default: 500)")
    parser.add_argument("--perms", type=int, default=999, help="RRPP permutations (default: 999)")
    parser.add_argument("--base-seed", type=int, default=2026, help="Parent seed, shared with D10 (default: 2026)")
    parser.add_argument("--alpha", type=float, default=0.05, help="Rejection threshold, p <= alpha (default: 0.05)")
    parser.add_argument("--k", type=float, default=2.0, help="Bound multiplier k (default: 2)")
    parser.add_argument(
        "--delta",
        type=float,
        default=DEFAULT_DELTA,
        help=f"Per-omic stage effect applied to delta_methyl, delta_expr and delta_protein (default: {DEFAULT_DELTA}; "
        f"p_dmp stays {P_DMP})",
    )
    parser.add_argument(
        "--s-z", type=float, default=DEFAULT_S_Z, help=f"Numeric covariate injection strength (default: {DEFAULT_S_Z})"
    )
    parser.add_argument(
        "--noise-covariates",
        type=int,
        default=DEFAULT_NOISE_COVARIATES,
        help=f"Pure-noise numeric covariates declared beside z (default: {DEFAULT_NOISE_COVARIATES}, i.e. six columns)",
    )
    parser.add_argument(
        "--auroc-band",
        type=parse_band,
        default=DEFAULT_AUROC_BAND,
        help="Calibration band low,high for the adjusted arm's mean held-out CV AUROC (default: 0.60,0.75)",
    )
    parser.add_argument(
        "--rho-threshold",
        type=float,
        default=DEFAULT_RHO_THRESHOLD,
        help="Calibration threshold on the unadjusted arm's median first-component |Spearman| with z (default: 0.8)",
    )
    parser.add_argument("--a-g", type=float, default=0.5, help="Covariate group effect a_g (default: 0.5)")
    parser.add_argument("--a-s", type=float, default=1.0, help="Covariate stage effect a_s (default: 1.0)")
    parser.add_argument("--a-gs", type=float, default=2.0, help="Covariate group x stage effect a_gs (default: 2.0)")
    parser.add_argument(
        "--feature-fraction", type=float, default=0.2, help="Fraction of features per direction (default: 0.2)"
    )
    parser.add_argument(
        "--translation-effect",
        type=float,
        default=DEFAULT_TRANSLATION_EFFECT,
        help="Translation control effect (default: 1.0)",
    )


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Case-study configuration Type I check at the Phase 6 design (no R)."
    )
    sub = parser.add_subparsers(dest="command", required=True)
    run = sub.add_parser("run", help="Evaluate (a shard of) the units; an unsharded run also merges")
    add_construction_arguments(run)
    run.add_argument("--shard", type=parse_shard, default=(0, 1), help="Run shard i of N (i/N; default 0/1)")
    run.add_argument("--n-jobs", type=int, default=1, help="Jobs for CV and RRPP (not part of the settings)")
    run.add_argument("--limit", type=int, default=None, help="Evaluate at most this many pending units")
    cal = sub.add_parser("calibrate", help="Anchor units, all three arms: the three calibration criteria")
    add_construction_arguments(cal)
    cal.add_argument("--workers", type=int, default=1, help="Replicates evaluated in parallel")
    merge_parser = sub.add_parser("merge", help="Merge every shard and write the decision (reads settings.json)")
    merge_parser.add_argument("--out-dir", type=Path, required=True, help="Output directory")
    combine = sub.add_parser("combine", help="The conjunction of several runs' decisions")
    combine.add_argument(
        "--decision", type=Path, action="append", required=True, help="A run's case_study_type_i_decision.json"
    )
    combine.add_argument("--out", type=Path, required=True, help="Where to write the overall decision")
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    if args.command == "combine":
        decisions = []
        for path in args.decision:
            decision = json.loads(path.read_text(encoding="utf-8"))
            decision.setdefault("label", str(path.parent.name))
            decisions.append(decision)
        overall = combine_decisions(decisions)
        args.out.parent.mkdir(parents=True, exist_ok=True)
        args.out.write_text(json.dumps(d10._jsonable(overall), indent=2, allow_nan=False) + "\n", encoding="utf-8")
        print(json.dumps(overall, indent=1))
        return 0
    if args.command == "merge":
        settings = load_settings(args.out_dir)
    elif args.command == "calibrate":
        settings = build_settings(args)
        entry = calibrate(args.out_dir, settings, workers=args.workers)
        print(json.dumps({key: value for key, value in entry.items() if key != "per_replicate"}, indent=1))
        return 0
    else:
        settings = build_settings(args)
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
