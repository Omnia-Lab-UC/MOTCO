#!/usr/bin/env python3
"""Corrected-permutation check at the case study's configuration (no R).

The validation behind ``openspec/changes/residual-df-corrected-permutation/``
(spec ``corrected-permutation-validation``). The case-study diagnosis
(:mod:`case_study_diagnosis`, 2026-10-30) returned ``residual_df_shrinkage``:
RRPP permutes reduced-model residuals whose mean square is RSS/n, so the
covariate-adjusted ``delta`` rejected at about 0.066. MOTCO 0.13.0 offers the
opt-in ``residual_correction`` (``"df"``: κ = √(n / (n − r)); ``"leverage"``:
1 / √(1 − h_ii)). This check measures the size of both corrections on data the
diagnosis never saw, and their power cost, through the shipped analysis path.

**Construction** (by import of :mod:`case_study_type_i_check`, itself built on
the D10 construction of :mod:`covariate_type_i_check`): the Phase 6 replicate
(three stages, methylation + expression, joint magnitude), a numeric covariate
``z`` confounded with group × stage and injected into the features at
``s_z = 10``, five pure-noise covariates (six numeric design columns), the rank
floor of 2 and the net-displacement estimator. Two cohort tables, read from the
case study's ``cohort_30`` ``result.json`` files (MOTCO 0.10.0, six covariate
columns each; the files' sha256 are recorded in the settings):

- ``astrocyte`` ``[[11, 10, 28], [9, 10, 12]]`` (n = 80, κ = √(80/70));
- ``microglia`` ``[[11, 10, 24], [9, 11, 11]]`` (n = 76, κ = √(76/66)).

Every seed (generator, ``z``, the noise covariates, RRPP; the injection
directions too) derives from base seed **2031** and the replicate index, a
family no committed check uses, so this is a confirmatory replication of the
diagnosis rather than a re-reading of its replicates. Every cell of a
replicate shares that replicate's seeds.

**Profiles** (``--profile``):

- ``type_i``: the anchor and the translation control at e = 1.0, per table and
  stage effect δ (``--delta``; the four runs are both tables at 0.1 and 0.05);
- ``power``: the astrocyte table at the generator's full stage effect
  (δ = 2.0), with the anchor, the translation control, magnitude at 0.02, 0.05,
  0.1 and 0.25, and orientation and shape at 0.25, 0.5, 0.75 and 1.0
  (14 cells).

**Units.** Each (cell, replicate) is analysed three times through
:func:`motco.analysis.analyze_trajectories` with the six covariates, the same
CV settings, permutation count and RRPP seed, differing only in
``residual_correction`` ∈ {``None``, ``"df"``, ``"leverage"``} (arms ``none``,
``df``, ``leverage``). The three latent-score matrices, selected ranks and
statistics must be identical; any difference raises and the unit is not
recorded. An arm the analysis rejects (for example a leverage-1 sample) is
recorded with its error and is unavailable for that unit.

**Decision** (pure, :func:`corrected_permutation_decision`; ``decide`` over
the five merged runs). A reading is one run × control × statistic; R its
available replicates; α = 0.05. A correction holds size when both tiers hold:

1. pooled per cohort table × statistic, over the table's four Type I readings
   (2 strengths × 2 controls): rate ≤ α + 2·√(α(1 − α)/R_pooled);
2. every reading (the four Type I runs' and the power profile's controls):
   rate ≤ α + 3·√(α(1 − α)/R).

``adopt_df`` when only ``df`` holds, ``adopt_leverage`` when only ``leverage``
does, ``no_correction_holds`` when neither does. When both hold, ``leverage``
is adopted only with a paired power gain over ``df`` at one-sided p ≤ 0.05 on
at least one target pair (magnitude/``delta``, orientation/``angle``,
shape/``shape``; each pair pooled over its mode's four effects, replicate-
clustered exact sign flip) and no such loss on any; otherwise ``df``. The
uncorrected arm goes through the same rule descriptively
(``uncorrected_excess_reproduced``: its pooled ``delta`` exceeds in either
table). A run below 999 permutations, or a missing run, marks the verdict
indicative. Because a table's Type I readings share replicate seeds (the two
strengths of one replicate are strongly dependent), the decision also reports
the pooled tier with a replicate-clustered standard error
(:func:`clustered_pooled_rows`), descriptively: it never changes the verdict.

Outputs (one ``--out-dir`` per run):

- ``settings.json``: settings and their hash; a resume under other settings is
  refused;
- ``records/shard_<i>_of_<N>.jsonl``: one line per unit (all three arms);
  resume skips units recorded in any shard file and re-runs torn lines;
- after ``merge`` (automatic for an unsharded run):
  ``corrected_permutation_pvalues.csv``, ``corrected_permutation_rates.csv``
  (per cell × arm × statistic: R, rejections, rate, rank distribution, CV AUROC
  and ``|ρ|`` with ``z``) and, for ``power``,
  ``corrected_permutation_contrasts.csv`` (each correction against none and
  ``leverage`` against ``df``, per cell and pooled per mode);
- ``decide`` (into its own ``--out-dir``): ``corrected_permutation_decision.json``
  and ``size_readings.csv``.

Examples:

    # smoke profile (99 permutations: indicative)
    python scripts/corrected_permutation_check.py run --profile type_i --table astrocyte --delta 0.1 \\
        --reps 20 --perms 99 --out-dir results/corrected-permutation-smoke/type-i-astrocyte-delta0.1

    # one SLURM array task, then the merge
    python scripts/corrected_permutation_check.py run --profile power --shard 3/140 --out-dir results/cp/power
    python scripts/corrected_permutation_check.py merge --out-dir results/cp/power

    # the verdict over the five runs
    python scripts/corrected_permutation_check.py decide --runs results/cp/type-i-* results/cp/power \\
        --out-dir results/cp
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import sys
from dataclasses import asdict
from pathlib import Path
from time import perf_counter
from typing import Any, Iterable, Mapping, Sequence

import numpy as np

from motco.analysis import (
    OmicsBlock,
    TrajectoryAnalysisError,
    TrajectoryAnalysisParams,
    TrajectoryDataset,
    analyze_trajectories,
)

sys.path.insert(0, str(Path(__file__).resolve().parent))
import case_study_diagnosis as diagnosis  # noqa: E402  (the replicate-clustered sign-flip contrast)
import case_study_type_i_check as cs  # noqa: E402  (construction, settings, sharding)
import covariate_type_i_check as d10  # noqa: E402  (the D10 construction and its helpers)

sys.path.pop(0)

SCRIPT_VERSION = 1
CHECK = "corrected_permutation"
STATISTICS: tuple[str, ...] = cs.STATISTICS
ANCHOR = cs.ANCHOR
CONTROL_MODES: tuple[str, ...] = cs.CONTROL_MODES
UNCORRECTED = "none"
ARM_CORRECTION: dict[str, str | None] = {UNCORRECTED: None, "df": "df", "leverage": "leverage"}
ARMS: tuple[str, ...] = tuple(ARM_CORRECTION)
#: The judged arms; the uncorrected arm is descriptive.
CORRECTIONS: tuple[str, ...] = ("df", "leverage")
PROFILES: tuple[str, ...] = ("type_i", "power")
BASE_SEED = 2031
S_Z = 10.0
NOISE_COVARIATES = 5
DEFAULT_REPLICATES = 500
DEFAULT_PERMUTATIONS = 999
TYPE_I_DELTAS: tuple[float, ...] = (0.1, 0.05)
POWER_DELTA = 2.0
POWER_TABLE = "astrocyte"
MODES: tuple[str, ...] = ("magnitude", "orientation", "shape")
#: The latent-adjustment power grid (Phase 6 effect levels below saturation).
POWER_EFFECTS: dict[str, list[float]] = {
    "magnitude": [0.02, 0.05, 0.1, 0.25],
    "orientation": [0.25, 0.5, 0.75, 1.0],
    "shape": [0.25, 0.5, 0.75, 1.0],
}
#: Target pairs of the tie-break, (mode, statistic).
TARGETS: tuple[tuple[str, str], ...] = (("magnitude", "delta"), ("orientation", "angle"), ("shape", "shape"))
#: The case study's ``cohort_30`` cell tables (six covariate columns each), read on ``ing`` 2026-10-30.
TABLES: dict[str, dict[str, Any]] = {
    "astrocyte": {
        "group_stage_sizes": [[11, 10, 28], [9, 10, 12]],
        "source": "sead_sexdiff results/motco_analysis/cohort_30/result.json (MOTCO 0.10.0)",
        "sha256": "34a744c28a062a9d9a1b17b793cba353a753fcac76b0f13b883565dfa9d212d3",
    },
    "microglia": {
        "group_stage_sizes": [[11, 10, 24], [9, 11, 11]],
        "source": "sead_sexdiff results/microglia/motco_analysis/cohort_30/result.json (MOTCO 0.10.0)",
        "sha256": "0dba07a78b72ba437074bfdaf731df82b28c7567b39940c51abe4b3a7c5752de",
    },
}
ALPHA = 0.05
POOLED_K = 2.0
READING_K = 3.0
CONTRAST_ALPHA = 0.05
VERDICTS: tuple[str, ...] = ("adopt_df", "adopt_leverage", "no_correction_holds")
CONSEQUENCES: dict[str, str] = {
    "adopt_df": (
        "sead_sexdiff re-issues both cell types' primaries with --residual-correction df, pinned to the release "
        "that contains it (0.13.0), in its own change; MOTCO documents df as validated at this configuration and "
        "keeps it opt-in"
    ),
    "adopt_leverage": (
        "sead_sexdiff re-issues both cell types' primaries with --residual-correction leverage, pinned to the "
        "release that contains it (0.13.0), in its own change; MOTCO documents leverage as validated at this "
        "configuration and keeps it opt-in"
    ),
    "no_correction_holds": (
        "the correction ships but is documented as not validated at this configuration; the case study keeps "
        "reading covariate-adjusted delta against its realized size and the correction stays unvalidated; the "
        "roadmap records the failing readings and names a diagnosis as the next step"
    ),
}
_RUN_DIRECTORY = "type-i-{table}-delta{delta:g}"


# ── Settings ─────────────────────────────────────────────────────────────────


def build_settings(args: argparse.Namespace) -> dict[str, Any]:
    """The case-study construction at the profile's table, δ and cells, with the arms and the rule."""

    profile = args.profile
    delta = float(args.delta) if args.delta is not None else (POWER_DELTA if profile == "power" else TYPE_I_DELTAS[0])
    if profile == "power" and (args.table != POWER_TABLE or delta != POWER_DELTA):
        raise SystemExit(f"the power profile runs the {POWER_TABLE} table at delta {POWER_DELTA}.")
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
            "--delta",
            repr(delta),
            "--s-z",
            repr(S_Z),
            "--noise-covariates",
            str(NOISE_COVARIATES),
        ]
    )
    settings = cs.build_settings(parent_args)
    settings.pop("k")
    settings["script_version"] = SCRIPT_VERSION
    settings["check"] = CHECK
    settings["profile"] = profile
    table = TABLES[args.table]
    settings["table"] = {"name": args.table, **table}
    settings["generator"]["group_stage_sizes"] = [list(row) for row in table["group_stage_sizes"]]
    settings["arms"] = {arm: {"residual_correction": ARM_CORRECTION[arm]} for arm in ARMS}
    settings["effects"] = {mode: list(POWER_EFFECTS[mode]) for mode in MODES} if profile == "power" else {}
    low, high = settings["calibration"]["auroc_band"]
    settings["calibration"] = {
        "auroc_band": [low, high],
        "rho_threshold": settings["calibration"]["rho_threshold"],
        "role": "descriptive; read from the uncorrected arm at the anchor, never judged",
        "criteria": {
            "auroc": "uncorrected-arm mean held-out CV stage AUROC at the anchor inside auroc_band",
            "rho": (
                "uncorrected-arm median |Spearman| of the first latent component with z at the anchor "
                ">= rho_threshold"
            ),
        },
        "not_run": "the parent's unadjusted-rejects criterion (this check runs no unadjusted arm)",
    }
    settings["rule"] = {
        "alpha": ALPHA,
        "reading": (
            "every reading (Type I runs and the power controls): rate <= alpha + reading_k*sqrt(alpha(1-alpha)/R)"
        ),
        "reading_k": READING_K,
        "pooled": (
            "per cohort table x statistic over the table's Type I readings: "
            "rate <= alpha + pooled_k*sqrt(alpha(1-alpha)/R_pooled)"
        ),
        "pooled_k": POOLED_K,
        "holds": "a correction holds size when every pooled and every reading verdict holds",
        "tie_break": (
            "when both hold: leverage only with a clustered one-sided sign-flip gain over df at p <= contrast_alpha "
            "on at least one target pair and no such loss on any; otherwise df"
        ),
        "targets": [list(target) for target in TARGETS],
        "contrast_alpha": CONTRAST_ALPHA,
        "uncorrected": "descriptive; never enters the verdict",
        "indicative": f"fewer than {DEFAULT_PERMUTATIONS} permutations in any run, or a missing run",
    }
    return settings


settings_hash = cs.settings_hash
ensure_settings = cs.ensure_settings
load_settings = cs.load_settings
parse_shard = cs.parse_shard
shard_path = cs.shard_path
unit_key = cs.unit_key


def run_label(settings: Mapping[str, Any]) -> str:
    """``power`` or ``type-i-<table>-delta<δ>``: the run's directory name and its rows' ``run``."""

    if settings["profile"] == "power":
        return "power"
    return _RUN_DIRECTORY.format(table=settings["table"]["name"], delta=cs.strengths(settings)["delta"])


def expected_runs() -> list[str]:
    """The five runs the paper-grade decision reads."""

    return [
        _RUN_DIRECTORY.format(table=table, delta=delta) for table in TABLES for delta in TYPE_I_DELTAS
    ] + ["power"]


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
    """Expected pool headroom of every pool-limited cell at the run's table and δ."""

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
        raise SystemExit(f"These cells exceed the pool headroom at the run's table and delta: {over}.")


# ── Unit evaluation ──────────────────────────────────────────────────────────


ARM_FIELDS: tuple[str, ...] = (
    "statistics",
    "p_values",
    "null_summary",
    "selected_lv",
    "component_selection",
    "cv_mean_auroc",
    "first_component_spearman_z",
    "residual_correction",
    "latent_digest",
    "error",
)


def arm_params(
    settings: Mapping[str, Any], seed: int, correction: str | None, *, n_jobs: int = 1
) -> TrajectoryAnalysisParams:
    return TrajectoryAnalysisParams(
        permutations=int(settings["permutations"]),
        seed=seed,
        n_jobs=n_jobs,
        orientation_estimator=settings["orientation_estimator"],
        min_components=int(settings["min_components"]),
        residual_correction=correction,
        **settings["cv"],
    )


def latent_digest(scores: np.ndarray) -> str:
    return hashlib.sha256(np.ascontiguousarray(scores, dtype=np.float64).tobytes()).hexdigest()


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
        "first_component_spearman_z": cs._spearman(scores[:, 0], z),
        "residual_correction": record["measurement"].get("residual_correction"),
        "latent_digest": latent_digest(scores),
        "error": None,
    }


def _failed_arm(exc: Exception) -> dict[str, Any]:
    failed: dict[str, Any] = {field: None for field in ARM_FIELDS}
    failed.update({"statistics": {}, "p_values": {}, "null_summary": {}, "error": str(exc)})
    return failed


def check_arms_agree(results: Mapping[str, Any]) -> None:
    """Raise unless every analysed arm has the same latent scores, selected rank and statistics."""

    arms = list(results)
    if len(arms) < 2:
        return
    first = results[arms[0]]
    differing = [
        arm
        for arm in arms[1:]
        if not (
            np.array_equal(first.latent_scores.to_numpy(), results[arm].latent_scores.to_numpy())
            and first.record["latent"]["selected_lv"] == results[arm].record["latent"]["selected_lv"]
            and first.record["statistics"] == results[arm].record["statistics"]
        )
    ]
    if differing:
        raise RuntimeError(
            f"arm(s) {differing} differ from arm {arms[0]!r} in their latent scores, selected rank or statistics; "
            "the unit is not recorded."
        )


def analyze_arms(built: Mapping[str, Any], settings: Mapping[str, Any], *, n_jobs: int = 1) -> dict[str, Any]:
    """Analyse one built unit under every arm: same data, covariates, CV settings and RRPP seed."""

    metadata = built["metadata"]
    dataset = TrajectoryDataset(
        blocks=[OmicsBlock(name, built["blocks"][name], "identity") for name in settings["layers"]],
        metadata=metadata,
        sample_col="sample_id",
        group_col="group",
        stage_col="stage",
        stage_order=[str(stage) for stage in sorted(set(metadata["stage"].tolist()))],
        covariate_cols=cs.declared_covariates(settings),
    )
    z = metadata["z"].to_numpy(dtype=float)
    arms: dict[str, Any] = {}
    results: dict[str, Any] = {}
    for arm in ARMS:
        params = arm_params(settings, built["seeds"]["rrpp"], ARM_CORRECTION[arm], n_jobs=n_jobs)
        try:
            result = analyze_trajectories(dataset, params)
        except TrajectoryAnalysisError as exc:
            arms[arm] = _failed_arm(exc)
            continue
        results[arm] = result
        arms[arm] = _arm_record(result, z)
    check_arms_agree(results)
    return arms


def evaluate_unit(
    mode: str, effect: float, replicate: int, settings: Mapping[str, Any], *, n_jobs: int = 1, reference: Any = None
) -> dict[str, Any]:
    start = perf_counter()
    built = cs.build_unit(mode, effect, replicate, settings, reference=reference)
    arms = analyze_arms(built, settings, n_jobs=n_jobs)
    return d10._jsonable(
        {
            "key": unit_key(mode, effect, replicate),
            "run": run_label(settings),
            "mode": mode,
            "effect": float(effect),
            "replicate": replicate,
            "settings_hash": settings_hash(settings),
            "seeds": built["seeds"],
            "group_labels": built["group_labels"],
            "arms": arms,
            "runtime_seconds": perf_counter() - start,
        }
    )


# ── Shards and resume ────────────────────────────────────────────────────────


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

    A unit recorded in any shard file under ``records/`` counts as done; this
    shard's file is first rewritten with its complete records, so a torn line
    is re-run instead of duplicated.
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


# ── Merge ────────────────────────────────────────────────────────────────────


PVALUE_COLUMNS = ("run", "mode", "effect", "replicate", "arm", "statistic", "p_value")
RATE_COLUMNS = (
    "run",
    "mode",
    "effect",
    "arm",
    "statistic",
    "units",
    "failures",
    "replicates",
    "rejections",
    "rate",
    "monte_carlo_se",
    "selected_rank",
    "cv_mean_auroc_mean",
    "median_abs_rho_z",
)
CONTRAST_COLUMNS = (
    "run",
    "scope",
    "mode",
    "effect",
    "statistic",
    "arm_b",
    "arm_a",
    "role",
    "readings",
    "clusters",
    "b_rate",
    "a_rate",
    "b_rejections",
    "a_rejections",
    "b_only",
    "a_only",
    "net",
    "p_b_greater",
    "p_a_greater",
)
#: (arm B, arm A) pairs of the contrast table: each correction against none, then leverage against df.
CONTRAST_PAIRS: tuple[tuple[str, str], ...] = (("df", UNCORRECTED), ("leverage", UNCORRECTED), ("leverage", "df"))
_MODE_ORDER: dict[str, int] = {mode: position for position, mode in enumerate((*CONTROL_MODES, *MODES))}


def _cell_order(cell: tuple[str, float]) -> tuple[int, float]:
    return _MODE_ORDER.get(cell[0], len(_MODE_ORDER)), cell[1]


def p_value_rows(records: Mapping[str, Mapping[str, Any]], run: str) -> list[dict[str, Any]]:
    rows = []
    for key in sorted(records):
        record = records[key]
        for arm in ARMS:
            p_values = record["arms"][arm].get("p_values") or {}
            for statistic in STATISTICS:
                rows.append(
                    {
                        "run": run,
                        "mode": record["mode"],
                        "effect": record["effect"],
                        "replicate": record["replicate"],
                        "arm": arm,
                        "statistic": statistic,
                        "p_value": p_values.get(statistic),
                    }
                )
    return rows


def _finite(value: Any) -> bool:
    return value is not None and value != "" and math.isfinite(float(value))


def rejection_tallies(rows: Iterable[Mapping[str, Any]], *, alpha: float) -> dict[tuple, list[int]]:
    """(run, mode, effect, arm, statistic) → [available replicates, rejections (``p ≤ alpha``)]."""

    tallies: dict[tuple, list[int]] = {}
    for row in rows:
        if not _finite(row.get("p_value")):
            continue
        key = (str(row["run"]), str(row["mode"]), float(row["effect"]), str(row["arm"]), str(row["statistic"]))
        tally = tallies.setdefault(key, [0, 0])
        tally[0] += 1
        tally[1] += int(float(row["p_value"]) <= alpha)
    return tallies


def rate_rows(records: Mapping[str, Mapping[str, Any]], run: str, *, alpha: float) -> list[dict[str, Any]]:
    """Per cell × arm × statistic: R, rejections, rate, Monte Carlo SE and the arm's latent summaries."""

    tallies = rejection_tallies(p_value_rows(records, run), alpha=alpha)
    by_cell: dict[tuple[str, float], list[Mapping[str, Any]]] = {}
    for record in records.values():
        by_cell.setdefault((str(record["mode"]), float(record["effect"])), []).append(record)
    rows = []
    for mode, effect in sorted(by_cell, key=_cell_order):
        cell_records = by_cell[(mode, effect)]
        for arm in ARMS:
            entries = [record["arms"][arm] for record in cell_records]
            ranks = cs.rank_distribution(entry.get("selected_lv") for entry in entries)
            auroc = cs.auroc_summary(entry.get("cv_mean_auroc") for entry in entries)
            rho = cs.correlation_summary(entry.get("first_component_spearman_z") for entry in entries)
            for statistic in STATISTICS:
                replicates, rejections, rate, se = cs._rate(tallies.get((run, mode, effect, arm, statistic)))
                rows.append(
                    {
                        "run": run,
                        "mode": mode,
                        "effect": effect,
                        "arm": arm,
                        "statistic": statistic,
                        "units": len(entries),
                        "failures": sum(1 for entry in entries if entry.get("error") is not None),
                        "replicates": replicates,
                        "rejections": rejections,
                        "rate": rate,
                        "monte_carlo_se": se,
                        "selected_rank": json.dumps(ranks, separators=(",", ":")),
                        "cv_mean_auroc_mean": auroc["mean"],
                        "median_abs_rho_z": rho["median_abs"],
                    }
                )
    return rows


def contrast_rows(
    rows: Iterable[Mapping[str, Any]], effects: Mapping[str, Sequence[float]], *, alpha: float
) -> list[dict[str, Any]]:
    """Paired contrasts per cell and pooled per mode, for each correction vs none and leverage vs df.

    ``p_b_greater`` is the replicate-clustered exact one-sided sign-flip p-value
    that arm B rejects more often than arm A; ``p_a_greater`` the reverse. The
    mode-pooled ``leverage``/``df`` rows of the target pairs are the tie-break.
    """

    rows = list(rows)
    runs = sorted({str(row["run"]) for row in rows})
    out = []
    for statistic in STATISTICS:
        rejected = diagnosis._rejections(rows, statistic, alpha=alpha)
        for run in runs:
            cell_list = sorted(
                {(str(row["mode"]), float(row["effect"])) for row in rows if str(row["run"]) == run}, key=_cell_order
            )
            scopes: list[tuple[str, str, float | None, set[tuple]]] = [
                ("cell", mode, effect, {(run, mode, effect)}) for mode, effect in cell_list
            ]
            scopes += [
                ("mode", mode, None, {(run, mode, float(effect)) for effect in effects[mode]})
                for mode in MODES
                if effects.get(mode)
            ]
            for scope, mode, effect, readings in scopes:
                for arm_b, arm_a in CONTRAST_PAIRS:
                    forward = diagnosis.paired_contrast(rejected, arm_b, arm_a, readings=readings)
                    reverse = diagnosis.paired_contrast(rejected, arm_a, arm_b, readings=readings)
                    paired = forward["readings"]
                    role = (
                        "tie_break"
                        if scope == "mode" and (arm_b, arm_a) == ("leverage", "df") and (mode, statistic) in TARGETS
                        else "descriptive"
                    )
                    out.append(
                        {
                            "run": run,
                            "scope": scope,
                            "mode": mode,
                            "effect": effect,
                            "statistic": statistic,
                            "arm_b": arm_b,
                            "arm_a": arm_a,
                            "role": role,
                            "readings": paired,
                            "clusters": forward["clusters"],
                            "b_rate": forward["b_rejections"] / paired if paired else None,
                            "a_rate": forward["a_rejections"] / paired if paired else None,
                            "b_rejections": forward["b_rejections"],
                            "a_rejections": forward["a_rejections"],
                            "b_only": forward["b_only"],
                            "a_only": forward["a_only"],
                            "net": forward["net"],
                            "p_b_greater": forward["p_value"],
                            "p_a_greater": reverse["p_value"],
                        }
                    )
    return out


def merge(out_dir: Path, settings: Mapping[str, Any]) -> dict[str, Any]:
    """Merge every shard of one run into its p-value, rate and (power) contrast tables."""

    out_dir = Path(out_dir)
    records = read_records(sorted((out_dir / "records").glob("shard_*.jsonl")))
    stale = sorted({str(record["settings_hash"]) for record in records.values()} - {settings_hash(settings)})
    if stale:
        raise SystemExit(f"{out_dir}: records written under other settings ({stale}); refusing to merge.")
    expected = [unit_key(*unit) for unit in units(settings)]
    run = run_label(settings)
    alpha = float(settings["alpha"])
    rows = p_value_rows(records, run)
    d10._write_csv(out_dir / "corrected_permutation_pvalues.csv", rows, PVALUE_COLUMNS)
    rates = rate_rows(records, run, alpha=alpha)
    d10._write_csv(out_dir / "corrected_permutation_rates.csv", rates, RATE_COLUMNS)
    if settings["effects"]:
        contrasts = contrast_rows(rows, settings["effects"], alpha=alpha)
        d10._write_csv(out_dir / "corrected_permutation_contrasts.csv", contrasts, CONTRAST_COLUMNS)
    return {
        "run": run,
        "units": len(records),
        "expected": len(expected),
        "missing": [key for key in expected if key not in records][:20],
        "failures": {
            arm: sum(1 for record in records.values() if record["arms"][arm].get("error") is not None) for arm in ARMS
        },
        "control_rates": {
            f"{row['mode']}/{row['arm']}/{row['statistic']}": row["rate"]
            for row in rates
            if row["mode"] in CONTROL_MODES
        },
        "runtime_seconds": cs._quantiles(
            sorted(float(record["runtime_seconds"]) for record in records.values() if "runtime_seconds" in record)
        ),
    }


# ── Decision (pure) ──────────────────────────────────────────────────────────


SIZE_COLUMNS = (
    "tier",
    "table",
    "run",
    "mode",
    "effect",
    "statistic",
    "arm",
    "role",
    "replicates",
    "rejections",
    "rate",
    "bound",
    "verdict",
    "clusters",
    "design_effect",
    "effective_replicates",
)


def size_bound(alpha: float, replicates: int, k: float) -> float:
    return d10.anchor_bound(alpha, replicates, k)


def _verdict(rate: float | None, bound: float | None) -> str:
    if rate is None or bound is None:
        return "unavailable"
    return "holds" if rate <= bound else "exceeds"


def size_rows_by_tier(
    rows: Iterable[Mapping[str, Any]],
    settings: Mapping[str, Mapping[str, Any]],
    *,
    alpha: float,
    reading_k: float,
    pooled_k: float,
) -> list[dict[str, Any]]:
    """Every reading (run × control × arm × statistic) and every pooled rate (table × arm × statistic)."""

    tallies = rejection_tallies(rows, alpha=alpha)
    out: list[dict[str, Any]] = []
    pooled: dict[tuple[str, str, str], list[int]] = {}
    for run in sorted(settings):
        run_settings = settings[run]
        table = str(run_settings["table"]["name"])
        type_i = run_settings["profile"] == "type_i"
        for mode in CONTROL_MODES:
            effect = float(run_settings["controls"][mode])
            for arm in ARMS:
                for statistic in STATISTICS:
                    replicates, rejections = tallies.get((run, mode, effect, arm, statistic), (0, 0))
                    rate = rejections / replicates if replicates else None
                    bound = size_bound(alpha, replicates, reading_k) if replicates else None
                    out.append(
                        {
                            "tier": "reading",
                            "table": table,
                            "run": run,
                            "mode": mode,
                            "effect": effect,
                            "statistic": statistic,
                            "arm": arm,
                            "role": "judged" if arm in CORRECTIONS else "descriptive",
                            "replicates": replicates,
                            "rejections": rejections,
                            "rate": rate,
                            "bound": bound,
                            "verdict": _verdict(rate, bound),
                        }
                    )
                    if type_i:
                        tally = pooled.setdefault((table, arm, statistic), [0, 0])
                        tally[0] += replicates
                        tally[1] += rejections
    for (table, arm, statistic), (replicates, rejections) in sorted(pooled.items()):
        rate = rejections / replicates if replicates else None
        bound = size_bound(alpha, replicates, pooled_k) if replicates else None
        out.append(
            {
                "tier": "pooled",
                "table": table,
                "run": None,
                "mode": None,
                "effect": None,
                "statistic": statistic,
                "arm": arm,
                "role": "judged" if arm in CORRECTIONS else "descriptive",
                "replicates": replicates,
                "rejections": rejections,
                "rate": rate,
                "bound": bound,
                "verdict": _verdict(rate, bound),
            }
        )
    return out


def clustered_pooled_rows(
    rows: Iterable[Mapping[str, Any]],
    settings: Mapping[str, Mapping[str, Any]],
    *,
    alpha: float,
    pooled_k: float,
) -> list[dict[str, Any]]:
    """The pooled tier with a replicate-clustered standard error (descriptive sensitivity; never judged).

    A table's Type I readings share replicate seeds, so they are not
    R_pooled independent replicates. Per table × arm × statistic, with c_r the
    rejections and m_r the available readings of replicate r, the pooled rate
    is p = Σc / N (N = Σm), its ratio-estimator variance
    R_c/(R_c − 1)·Σ(c_r − p·m_r)²/N² over the R_c replicates, and its design
    effect that variance over p(1 − p)/N. The bound is
    alpha + pooled_k·√(alpha(1 − alpha)/R_eff) with R_eff = N / design effect.
    """

    tables = {run: str(s["table"]["name"]) for run, s in settings.items() if s["profile"] == "type_i"}
    controls = {run: {mode: float(settings[run]["controls"][mode]) for mode in CONTROL_MODES} for run in tables}
    clusters: dict[tuple[str, str, str], dict[int, list[int]]] = {}
    for row in rows:
        run = str(row["run"])
        mode = str(row["mode"])
        if run not in tables or mode not in CONTROL_MODES or float(row["effect"]) != controls[run][mode]:
            continue
        if not _finite(row.get("p_value")):
            continue
        tally = clusters.setdefault((tables[run], str(row["arm"]), str(row["statistic"])), {}).setdefault(
            int(row["replicate"]), [0, 0]
        )
        tally[0] += int(float(row["p_value"]) <= alpha)
        tally[1] += 1
    out = []
    for table in sorted(set(tables.values())):
        for arm in ARMS:
            for statistic in STATISTICS:
                per_replicate = list(clusters.get((table, arm, statistic), {}).values())
                count = len(per_replicate)
                total = sum(m for _, m in per_replicate)
                rejections = sum(c for c, _ in per_replicate)
                rate = rejections / total if total else None
                design_effect = effective = bound = None
                if rate is not None and count > 1 and 0.0 < rate < 1.0:
                    variance = count / (count - 1) * sum((c - rate * m) ** 2 for c, m in per_replicate) / total**2
                    design_effect = variance / (rate * (1.0 - rate) / total)
                    effective = total / design_effect
                    bound = size_bound(alpha, effective, pooled_k)
                out.append(
                    {
                        "tier": "pooled_clustered",
                        "table": table,
                        "run": None,
                        "mode": None,
                        "effect": None,
                        "statistic": statistic,
                        "arm": arm,
                        "role": "sensitivity",
                        "replicates": total,
                        "rejections": rejections,
                        "rate": rate,
                        "bound": bound,
                        "verdict": _verdict(rate, bound),
                        "clusters": count,
                        "design_effect": design_effect,
                        "effective_replicates": effective,
                    }
                )
    return out


def _failing(size: Sequence[Mapping[str, Any]], arm: str) -> list[dict[str, Any]]:
    keys = ("tier", "table", "run", "mode", "statistic", "replicates", "rejections", "rate", "bound", "verdict")
    return [{key: row[key] for key in keys} for row in size if row["arm"] == arm and row["verdict"] != "holds"]


def _mode_readings(run: str, effects: Mapping[str, Sequence[float]], mode: str) -> set[tuple]:
    return {(run, mode, float(effect)) for effect in effects.get(mode, [])}


def tie_break(
    power_rows: Sequence[Mapping[str, Any]],
    run: str | None,
    effects: Mapping[str, Sequence[float]],
    *,
    alpha: float,
    contrast_alpha: float,
) -> dict[str, Any]:
    """``leverage`` against ``df`` on each target pair, pooled over the mode's effects and clustered by replicate."""

    pairs = []
    for mode, statistic in TARGETS:
        rejected = diagnosis._rejections(power_rows, statistic, alpha=alpha)
        readings = _mode_readings(run, effects, mode) if run is not None else set()
        gain = diagnosis.paired_contrast(rejected, "leverage", "df", readings=readings)
        loss = diagnosis.paired_contrast(rejected, "df", "leverage", readings=readings)
        pairs.append(
            {
                "mode": mode,
                "statistic": statistic,
                "readings": gain["readings"],
                "clusters": gain["clusters"],
                "leverage_rejections": gain["b_rejections"],
                "df_rejections": gain["a_rejections"],
                "leverage_only": gain["b_only"],
                "df_only": gain["a_only"],
                "p_gain": gain["p_value"],
                "p_loss": loss["p_value"],
                "gain": gain["readings"] > 0 and gain["p_value"] <= contrast_alpha,
                "loss": gain["readings"] > 0 and loss["p_value"] <= contrast_alpha,
            }
        )
    available = all(pair["readings"] > 0 for pair in pairs)
    leverage = available and any(pair["gain"] for pair in pairs) and not any(pair["loss"] for pair in pairs)
    return {
        "status": "computed" if available else "unavailable",
        "pairs": pairs,
        "choice": "leverage" if leverage else "df",
    }


def power_cost(
    power_rows: Sequence[Mapping[str, Any]],
    run: str | None,
    effects: Mapping[str, Sequence[float]],
    *,
    alpha: float,
) -> list[dict[str, Any]]:
    """Each correction against none, per mode (pooled over its effects) and statistic: rates and paired p-values."""

    out = []
    if run is None:
        return out
    for statistic in STATISTICS:
        rejected = diagnosis._rejections(power_rows, statistic, alpha=alpha)
        for mode in MODES:
            readings = _mode_readings(run, effects, mode)
            if not readings:
                continue
            for correction in CORRECTIONS:
                loss = diagnosis.paired_contrast(rejected, UNCORRECTED, correction, readings=readings)
                gain = diagnosis.paired_contrast(rejected, correction, UNCORRECTED, readings=readings)
                paired = loss["readings"]
                out.append(
                    {
                        "correction": correction,
                        "mode": mode,
                        "statistic": statistic,
                        "readings": paired,
                        "uncorrected_rate": loss["b_rejections"] / paired if paired else None,
                        "corrected_rate": loss["a_rejections"] / paired if paired else None,
                        "uncorrected_only": loss["b_only"],
                        "corrected_only": loss["a_only"],
                        "p_loss": loss["p_value"],
                        "p_gain": gain["p_value"],
                    }
                )
    return out


def corrected_permutation_decision(
    size_rows: Iterable[Mapping[str, Any]],
    power_rows: Iterable[Mapping[str, Any]],
    settings: Mapping[str, Mapping[str, Any]],
) -> dict[str, Any]:
    """The verdict from the merged p-value rows and every run's settings (keyed by run label).

    ``size_rows`` are the control rows of every run (``run``, ``mode``,
    ``effect``, ``replicate``, ``arm``, ``statistic``, ``p_value``);
    ``power_rows`` the power profile's rows. The rule constants are read from
    the settings, which must agree on them.
    """

    rules = {json.dumps(run_settings["rule"], sort_keys=True) for run_settings in settings.values()}
    if len(rules) != 1:
        raise ValueError("the runs' settings disagree on the decision rule.")
    rule = next(iter(settings.values()))["rule"]
    alpha, reading_k, pooled_k = float(rule["alpha"]), float(rule["reading_k"]), float(rule["pooled_k"])
    contrast_alpha = float(rule["contrast_alpha"])
    size_rows = list(size_rows)
    size = size_rows_by_tier(size_rows, settings, alpha=alpha, reading_k=reading_k, pooled_k=pooled_k)
    clustered = clustered_pooled_rows(size_rows, settings, alpha=alpha, pooled_k=pooled_k)
    holds = {
        arm: bool(verdicts) and all(verdict == "holds" for verdict in verdicts)
        for arm in ARMS
        for verdicts in [[row["verdict"] for row in size if row["arm"] == arm]]
    }
    power_runs = [run for run, run_settings in settings.items() if run_settings["profile"] == "power"]
    if len(power_runs) > 1:
        raise ValueError(f"at most one power run; got {power_runs}.")
    power_run = power_runs[0] if power_runs else None
    effects = settings[power_run]["effects"] if power_run is not None else {}
    power_rows = [row for row in power_rows if power_run is not None and str(row["run"]) == power_run]
    choice: dict[str, Any] | None = None
    if holds["df"] and holds["leverage"]:
        choice = tie_break(power_rows, power_run, effects, alpha=alpha, contrast_alpha=contrast_alpha)
        verdict = f"adopt_{choice['choice']}"
    elif holds["df"]:
        verdict = "adopt_df"
    elif holds["leverage"]:
        verdict = "adopt_leverage"
    else:
        verdict = "no_correction_holds"
    uncorrected_pooled = [
        {key: row[key] for key in ("table", "statistic", "replicates", "rejections", "rate", "bound", "verdict")}
        for row in size
        if row["tier"] == "pooled" and row["arm"] == UNCORRECTED
    ]
    missing = [run for run in expected_runs() if run not in settings]
    permutations = {run: int(run_settings["permutations"]) for run, run_settings in settings.items()}
    below = sorted(run for run, count in permutations.items() if count < DEFAULT_PERMUTATIONS)
    return {
        "rule": rule,
        "verdict": verdict,
        "consequence": CONSEQUENCES[verdict],
        "indicative": bool(below or missing),
        "indicative_reasons": {"below_default_permutations": below, "missing_runs": missing},
        "holds": {arm: holds[arm] for arm in CORRECTIONS},
        "failing": {arm: _failing(size, arm) for arm in CORRECTIONS},
        "tie_break": choice,
        "uncorrected": {
            "role": "descriptive; never enters the verdict",
            "holds": holds[UNCORRECTED],
            "pooled": uncorrected_pooled,
            "uncorrected_excess_reproduced": any(
                row["statistic"] == "delta" and row["verdict"] == "exceeds" for row in uncorrected_pooled
            ),
            "failing": _failing(size, UNCORRECTED),
        },
        "power_cost": power_cost(power_rows, power_run, effects, alpha=alpha),
        "readings": [row for row in size if row["tier"] == "reading"],
        "pooled": [row for row in size if row["tier"] == "pooled"],
        "clustered_sensitivity": {
            "role": (
                "descriptive; never enters the verdict. The pooled tier with a replicate-clustered standard error, "
                "because a table's Type I readings share replicate seeds (decided 2026-10-30, before any "
                "paper-grade result)"
            ),
            "holds": {
                arm: all(row["verdict"] == "holds" for row in size if row["arm"] == arm and row["tier"] == "reading")
                and all(row["verdict"] == "holds" for row in clustered if row["arm"] == arm)
                for arm in ARMS
            },
            "pooled": clustered,
        },
        "runs": {
            run: {
                "profile": run_settings["profile"],
                "table": run_settings["table"]["name"],
                "delta": cs.strengths(run_settings)["delta"],
                "replicates": int(run_settings["replicates"]),
                "permutations": permutations[run],
                "base_seed": int(run_settings["base_seed"]),
                "settings_hash": settings_hash(run_settings),
            }
            for run, run_settings in sorted(settings.items())
        },
    }


def calibration_report(rate_rows_: Iterable[Mapping[str, Any]], settings: Mapping[str, Any]) -> dict[str, Any]:
    """The descriptive calibration criteria of one run, from its rates table (uncorrected arm, anchor)."""

    calibration = settings["calibration"]
    low, high = (float(value) for value in calibration["auroc_band"])
    threshold = float(calibration["rho_threshold"])
    anchor = [row for row in rate_rows_ if row["mode"] == ANCHOR and row["arm"] == UNCORRECTED]
    auroc = float(anchor[0]["cv_mean_auroc_mean"]) if anchor and _finite(anchor[0]["cv_mean_auroc_mean"]) else None
    rho = float(anchor[0]["median_abs_rho_z"]) if anchor and _finite(anchor[0]["median_abs_rho_z"]) else None
    return {
        "role": calibration["role"],
        "auroc": {"band": [low, high], "mean": auroc, "in_band": auroc is not None and low <= auroc <= high},
        "rho": {"threshold": threshold, "median_abs": rho, "met": rho is not None and rho >= threshold},
        "selected_rank": json.loads(anchor[0]["selected_rank"]) if anchor else {},
    }


def read_run(run_dir: Path) -> dict[str, Any]:
    run_dir = Path(run_dir)
    settings = load_settings(run_dir)
    missing = [
        name
        for name in ("corrected_permutation_pvalues.csv", "corrected_permutation_rates.csv")
        if not (run_dir / name).exists()
    ]
    if missing:
        raise SystemExit(f"{run_dir}: {', '.join(missing)} missing; run merge first.")
    return {
        "settings": settings,
        "rows": diagnosis._read_csv(run_dir / "corrected_permutation_pvalues.csv"),
        "rates": diagnosis._read_csv(run_dir / "corrected_permutation_rates.csv"),
    }


def decide(run_dirs: Sequence[Path], out_dir: Path) -> dict[str, Any]:
    """Read the merged runs and write ``corrected_permutation_decision.json`` and ``size_readings.csv``."""

    runs = {}
    for run_dir in run_dirs:
        run = read_run(run_dir)
        label = run_label(run["settings"])
        if label in runs:
            raise SystemExit(f"decide needs distinct runs; {label} given twice.")
        runs[label] = run
    seeds = {int(run["settings"]["base_seed"]) for run in runs.values()}
    if len(seeds) != 1:
        raise SystemExit(f"decide needs runs sharing one base seed; got {sorted(seeds)}.")
    settings = {label: run["settings"] for label, run in runs.items()}
    size_rows = [row for run in runs.values() for row in run["rows"] if row["mode"] in CONTROL_MODES]
    power_rows = [row for run in runs.values() if run["settings"]["profile"] == "power" for row in run["rows"]]
    decision = corrected_permutation_decision(size_rows, power_rows, settings)
    decision["calibration"] = {label: calibration_report(run["rates"], run["settings"]) for label, run in runs.items()}
    decision["run_directories"] = {
        run_label(run["settings"]): str(Path(run_dir)) for run_dir, run in zip(run_dirs, runs.values())
    }
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    d10._write_csv(
        out_dir / "size_readings.csv",
        [*decision["readings"], *decision["pooled"], *decision["clustered_sensitivity"]["pooled"]],
        SIZE_COLUMNS,
    )
    (out_dir / "corrected_permutation_decision.json").write_text(
        json.dumps(d10._jsonable(decision), indent=2, allow_nan=False) + "\n", encoding="utf-8"
    )
    return decision


# ── CLI ──────────────────────────────────────────────────────────────────────


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Corrected-permutation check at the case study's configuration (no R)."
    )
    sub = parser.add_subparsers(dest="command", required=True)
    run = sub.add_parser("run", help="Evaluate (a shard of) one run's units; an unsharded run also merges")
    run.add_argument("--profile", choices=PROFILES, required=True, help="type_i or power")
    run.add_argument(
        "--table", choices=sorted(TABLES), default=POWER_TABLE, help="Cohort cell table (default: astrocyte)"
    )
    run.add_argument(
        "--delta",
        type=float,
        default=None,
        help=f"Per-omic stage effect (default: {TYPE_I_DELTAS[0]} for type_i, {POWER_DELTA} for power)",
    )
    run.add_argument("--reps", type=int, default=DEFAULT_REPLICATES, help=f"Replicates (default: {DEFAULT_REPLICATES})")
    run.add_argument(
        "--perms", type=int, default=DEFAULT_PERMUTATIONS, help=f"RRPP permutations (default: {DEFAULT_PERMUTATIONS})"
    )
    run.add_argument("--base-seed", type=int, default=BASE_SEED, help=f"Seed family (default: {BASE_SEED})")
    run.add_argument("--out-dir", type=Path, required=True, help="Output directory (one per run)")
    run.add_argument("--shard", type=parse_shard, default=(0, 1), help="Run shard i of N (i/N; default 0/1)")
    run.add_argument("--n-jobs", type=int, default=1, help="Jobs for CV and RRPP (not part of the settings)")
    run.add_argument("--limit", type=int, default=None, help="Evaluate at most this many pending units")
    merge_parser = sub.add_parser("merge", help="Merge every shard of one run (reads settings.json)")
    merge_parser.add_argument("--out-dir", type=Path, required=True, help="Run directory")
    decide_parser = sub.add_parser("decide", help="The verdict over the merged runs")
    decide_parser.add_argument("--runs", type=Path, nargs="+", required=True, help="Merged run directories")
    decide_parser.add_argument("--out-dir", type=Path, required=True, help="Where to write the decision")
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    if args.command == "decide":
        decision = decide(args.runs, args.out_dir)
        summary = {key: decision[key] for key in ("verdict", "indicative", "indicative_reasons", "holds")}
        summary["uncorrected_excess_reproduced"] = decision["uncorrected"]["uncorrected_excess_reproduced"]
        summary["clustered_sensitivity_holds"] = decision["clustered_sensitivity"]["holds"]
        print(json.dumps(summary, indent=1))
        return 0
    if args.command == "merge":
        settings = load_settings(args.out_dir)
        print(json.dumps(merge(args.out_dir, settings), indent=1))
        return 0
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
