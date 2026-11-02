#!/usr/bin/env python3
"""Covariate-adjusted attribution check at the case study's configuration (no R).

The gate behind ``openspec/changes/covariate-adjusted-attribution/`` (spec
``adjusted-attribution-validation``): does covariate-adjusted orientation
attribution, as :func:`motco.analysis.analyze_trajectories` runs it, remove a
confounded covariate's contribution from the drivers, and does it recover and
stabilize planted orientation drivers as well as an oracle that knows the
injection?

**Construction** (by import of :mod:`case_study_type_i_check`, which composes
the D10 construction of :mod:`covariate_type_i_check`): the Phase 6 replicate
on the astrocyte ``cohort_30`` table ``[[11, 10, 28], [9, 10, 12]]`` (n = 80),
methylation + expression (498 joint features), joint magnitude, at the full
Phase 6 stage effect (δ = 2.0) so the planted orientation has a trajectory to
rotate. A covariate ``z = a_g·g + a_s·s + a_gs·g·s + ε`` confounded with group,
stage and their interaction is injected into a fixed 20 % feature subset
(``X_f += s_z·sd_f·w_f·z``), and five pure-noise covariates sit beside it (six
numeric design columns). Seeds derive from base seed 2033 (no committed check
uses it) and the replicate only, so the three cells of a replicate share the
baseline draw, ``z``, the injection, the RRPP seed and the bootstrap seed.

**Cells.** The zero-effect anchor and ``orientation`` at e = 0.5 and 1.0, 500
replicates each, 999 permutations, B = 200 bootstrap replicates, top-k 20.

**Arms**, on one measurement per unit:

- ``adjusted``: ``analyze_trajectories`` with the six covariates, attribution,
  rank floor 2, net displacement and ``residual_correction="df"`` — the shipped
  path end to end;
- ``oracle``: the attribution primitive on the analysed features minus the
  exact injection in analysed coordinates, ``X̃ − z ⊗ (v / σ)`` with ``v`` the
  raw-unit injection and ``σ`` the pooled analysed scale, so its cell means are
  the uninjected cell means up to one constant vector;
- ``unadjusted``: the primitive on the analysed features (descriptive).

All three use the analysis's exposed estimator (``result.latent_artifacts``)
and share the bootstrap seed, so they resample the same rows. Each unit first
checks that the estimator reproduces the latent scores and that the adjusted
means project to the analysis's stage means; a failed check raises and the unit
is not recorded.

**Metrics** per arm, block and component (transitions ``s->s+1`` and the
``principal`` block; observed, PLS-captured and residual): top-k truth
precision over min(k, |truth|), top-k contamination (the share of the top-k
that is injected and not a truth driver), the contrast's projection on the unit
injection direction, cosine and top-k overlap with the oracle, the mean
bootstrap sign stability and top-k frequency over the arm's own top-k, and the
valid and rank-deficient bootstrap counts.

**Decision** (pure, :func:`adjusted_attribution_decision`), read on the
observed principal-orientation contrast at 2 paired standard errors:

- gates (rule version 2, amended 2026-11-02): unadjusted − oracle
  contamination ≥ 0.25 at the anchor, |unadjusted − oracle projection| ≥ 0.2
  at every cell, and oracle truth precision ≥ 0.5 at orientation e = 1.0 —
  else ``construction_uninformative``. The contamination gate was declared at
  every cell; at the orientation cells the planted drivers crowd every
  injected feature out of the unadjusted top-k, so it moved to the anchor and
  the projection gate shows the confounding reaches the contrast there;
- criteria: adjusted − oracle contamination ≤ 0.05 and residual confounding
  ≤ 0.10 × the unadjusted arm's at every cell; truth precision and sign
  stability ≥ −0.05 and top-k frequency ≥ −0.10 against the oracle at both
  orientation cells — ``adjusted_attribution_validated`` when all hold, else
  ``diagnose_adjusted_attribution``.

The unadjusted arm, the captured and residual components, the transition
blocks and the readings conditional on ``angle`` rejection are reported
beside the verdict and never change it. A verdict below 999 permutations,
B = 200 or 500 replicates per cell is indicative.

**Calibration** (``calibrate``): R = 20 × 99 permutations × B = 50 at the
anchor and orientation e = 1.0, starting at ``--s-z`` (10) and stepping by 5
until the confounding criterion is met (the anchor's contamination gap and
both cells' projection gap reach their gate margins); a strength that meets it
but fails recoverability stops the check for the user's decision. Every
iteration is appended to ``calibration.json``.

Outputs (``--out-dir``): ``settings.json``; ``records/shard_<i>_of_<N>.jsonl``
(one line per unit; resume skips recorded units); after ``merge``:
``merged.jsonl``, ``adjusted_attribution_metrics.csv`` (per unit × arm ×
block × component) and ``adjusted_attribution_summary.csv`` (per cell × arm ×
block × component: mean and SE of every metric and of every paired difference
against the oracle); after ``decide``: ``adjusted_attribution_decision.json``.

With ``--workers`` above 1, pin BLAS to one thread per process
(``OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1``), as the
cluster arrays do; otherwise every worker spawns a full thread pool and the
machine is oversubscribed.

Examples:

    export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1

    # calibration (appends every iteration to calibration.json)
    python scripts/adjusted_attribution_check.py calibrate --out-dir results/aa/calibration --workers 4

    # smoke at the calibrated strength (indicative)
    python scripts/adjusted_attribution_check.py run --reps 20 --perms 99 --bootstrap 50 --s-z 10 \\
        --calibration results/aa/calibration/calibration.json --out-dir results/aa/smoke --workers 4

    # paper grade, sharded, then merge and decide
    python scripts/adjusted_attribution_check.py run --shard 3/30 --s-z 10 --out-dir results/aa/run
    python scripts/adjusted_attribution_check.py merge --out-dir results/aa/run
    python scripts/adjusted_attribution_check.py decide --out-dir results/aa/run
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
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
import case_study_type_i_check as cs  # noqa: E402  (construction, settings, sharding)
import covariate_type_i_check as d10  # noqa: E402  (the D10 construction and its helpers)

sys.path.pop(0)

SCRIPT_VERSION = 1
CHECK = "adjusted_attribution"
BASE_SEED = 2033
DELTA = 2.0
DEFAULT_S_Z = 10.0
S_Z_STEP = 5.0
DEFAULT_MAX_S_Z = 50.0
NOISE_COVARIATES = 5
DEFAULT_REPLICATES = 500
DEFAULT_PERMUTATIONS = 999
DEFAULT_BOOTSTRAP = 200
DEFAULT_TOP_K = 20
CALIBRATION_REPLICATES = 20
CALIBRATION_PERMUTATIONS = 99
CALIBRATION_BOOTSTRAP = 50
RESIDUAL_CORRECTION = "df"
ANCHOR = cs.ANCHOR
ORIENTATION_EFFECTS: tuple[float, ...] = (0.5, 1.0)
#: The astrocyte ``cohort_30`` table, item 12's cohort (``corrected_permutation_check.TABLES``).
TABLE: dict[str, Any] = {
    "name": "astrocyte",
    "group_stage_sizes": [[11, 10, 28], [9, 10, 12]],
    "source": "sead_sexdiff results/motco_analysis/cohort_30/result.json (MOTCO 0.10.0)",
    "sha256": "34a744c28a062a9d9a1b17b793cba353a753fcac76b0f13b883565dfa9d212d3",
}
ARMS: tuple[str, ...] = ("adjusted", "oracle", "unadjusted")
COMPARED_ARMS: tuple[str, ...] = ("adjusted", "unadjusted")
COMPONENTS: tuple[str, ...] = ("observed", "pls_captured", "residual")
PRINCIPAL = "principal"
ALPHA = 0.05
SE_MULTIPLIER = 2.0
#: The predeclared margins (design D10); changing one needs the user's approval recorded in design.md.
MARGINS: dict[str, float] = {
    "contamination": 0.05,
    "residual_confounding_ratio": 0.10,
    "truth_precision": -0.05,
    "sign_stability": -0.05,
    "top_k_frequency": -0.10,
    "contamination_gate": 0.25,
    "projection_gate": 0.2,
    "recoverability_gate": 0.5,
}
#: 2 = the gates amended 2026-11-02 (design D10): the contamination gate at the anchor only, plus a
#: projection gate at every cell, because the planted drivers crowd every injected feature out of the
#: unadjusted top-k at the orientation cells (calibration.json iterations 0-8, stop_exhausted).
RULE_VERSION = 2
METRICS: tuple[str, ...] = (
    "truth_precision",
    "contamination",
    "projection",
    "cosine_oracle",
    "overlap_oracle",
    "sign_stability",
    "top_k_frequency",
    "valid_replicates",
)
#: Metrics also summarized as paired differences against the oracle.
PAIRED_METRICS: tuple[str, ...] = (
    "truth_precision",
    "contamination",
    "projection",
    "sign_stability",
    "top_k_frequency",
)
VERDICTS: tuple[str, ...] = (
    "adjusted_attribution_validated",
    "diagnose_adjusted_attribution",
    "construction_uninformative",
)
CONSEQUENCES: dict[str, str] = {
    "adjusted_attribution_validated": (
        "covariate-adjusted attribution is validated at the case study's configuration; item 12 (the plasmode "
        "demonstration on SEA-AD astrocyte cohort_30) may start on the released version (0.14.0)"
    ),
    "diagnose_adjusted_attribution": (
        "the attribution leg stays withheld; a diagnosis of the failing criterion is proposed before item 12"
    ),
    "construction_uninformative": (
        "the attribution leg stays withheld; the construction is returned to the user, because a failed gate "
        "means the check could not judge adjustment"
    ),
}


# ── Settings ─────────────────────────────────────────────────────────────────


def build_settings(args: argparse.Namespace) -> dict[str, Any]:
    """The case-study construction on the astrocyte table at δ = 2.0, with the cells, arms and rule."""

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
            repr(DELTA),
            "--s-z",
            repr(float(args.s_z)),
            "--noise-covariates",
            str(NOISE_COVARIATES),
        ]
    )
    settings = cs.build_settings(parent_args)
    for key in ("k", "controls", "calibration", "rule"):
        settings.pop(key)
    settings["script_version"] = SCRIPT_VERSION
    settings["check"] = CHECK
    settings["table"] = dict(TABLE)
    settings["generator"]["group_stage_sizes"] = [list(row) for row in TABLE["group_stage_sizes"]]
    settings["cells"] = [[mode, effect] for mode, effect in cells_from_effects(ORIENTATION_EFFECTS)]
    settings["residual_correction"] = RESIDUAL_CORRECTION
    settings["attribution"] = {
        "bootstrap_replicates": int(args.bootstrap),
        "top_k": int(args.top_k),
        "bootstrap_seed": "derived per replicate from the base seed (shared by every arm and cell)",
        "zero_tolerance": 1e-12,
    }
    settings["arms"] = {
        "adjusted": "analyze_trajectories with the six covariates, attribution and residual_correction df",
        "oracle": "the primitive on X - z (x) (v / sigma): the exact injection removed in analysed coordinates",
        "unadjusted": "the primitive on the analysed features (descriptive)",
    }
    settings["rule"] = {
        "reading": "observed component of the principal-orientation block",
        "se_multiplier": SE_MULTIPLIER,
        "margins": dict(MARGINS),
        "version": RULE_VERSION,
        "gates": (
            "contamination: mean(unadjusted - oracle) - k*SE >= contamination_gate at the anchor; projection: "
            "|mean projection difference (unadjusted - oracle)| - k*SE >= projection_gate at every cell; "
            "recoverability: oracle mean truth precision at orientation e = 1.0 >= recoverability_gate"
        ),
        "criteria": (
            "contamination: mean(adjusted - oracle) + k*SE <= contamination at every cell; residual confounding: "
            "|mean projection difference (adjusted - oracle)| + k*SE <= residual_confounding_ratio * |mean "
            "projection difference (unadjusted - oracle)| at every cell; truth precision, sign stability and top-k "
            "frequency: mean(adjusted - oracle) - k*SE >= their margin at both orientation cells"
        ),
        "verdicts": list(VERDICTS),
        "indicative": (
            f"fewer than {DEFAULT_PERMUTATIONS} permutations, {DEFAULT_BOOTSTRAP} bootstrap replicates or "
            f"{DEFAULT_REPLICATES} available replicates in any cell"
        ),
        "descriptive": (
            "the unadjusted arm's criteria, the PLS-captured and residual components, the transition blocks and "
            "the readings conditional on angle rejection (adjusted arm, p <= alpha); none changes the verdict"
        ),
    }
    calibration = getattr(args, "calibration", None)
    settings["adopted_calibration"] = (
        adopted_calibration(Path(calibration), float(args.s_z)) if calibration is not None else None
    )
    return settings


settings_hash = cs.settings_hash
ensure_settings = cs.ensure_settings
load_settings = cs.load_settings
parse_shard = cs.parse_shard
shard_path = cs.shard_path
unit_key = cs.unit_key


def _seed(*parts: Any) -> int:
    encoded = json.dumps([CHECK, *parts], separators=(",", ":")).encode("utf-8")
    return int(hashlib.sha256(encoded).hexdigest()[:8], 16) & 0x7FFFFFFF


def bootstrap_seed(base_seed: int, replicate: int) -> int:
    return _seed(int(base_seed), "bootstrap", int(replicate))


def adopted_calibration(path: Path, s_z: float) -> dict[str, Any]:
    """The calibration iteration the settings adopt: the last one at ``s_z`` that met both criteria."""

    if not path.exists():
        raise SystemExit(f"{path} does not exist; run calibrate first or drop --calibration.")
    history = json.loads(path.read_text(encoding="utf-8"))
    for index in range(len(history) - 1, -1, -1):
        entry = history[index]
        if float(entry["s_z"]) == float(s_z) and entry["met"]:
            digest = hashlib.sha256(json.dumps(entry, sort_keys=True).encode("utf-8")).hexdigest()
            return {"iteration": index, "s_z": float(s_z), "entry_sha256": digest}
    raise SystemExit(f"{path} holds no iteration at s_z = {s_z:g} that met both criteria.")


# ── Grid ─────────────────────────────────────────────────────────────────────


def cells_from_effects(effects: Sequence[float]) -> list[tuple[str, float]]:
    return [(ANCHOR, 0.0), *(("orientation", float(effect)) for effect in effects)]


def cells(settings: Mapping[str, Any]) -> list[tuple[str, float]]:
    """(mode, effect) per cell: the anchor, then orientation at each effect."""

    return [(str(mode), float(effect)) for mode, effect in settings["cells"]]


def cell_label(mode: str, effect: float) -> str:
    return ANCHOR if mode == ANCHOR else f"{mode}@{float(effect):g}"


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
    """Refuse settings with a cell beyond the generator's surgery headroom, naming it."""

    over = [key for key, headroom in headroom_report(settings, reference).items() if not headroom["fits"]]
    if over:
        raise SystemExit(f"These cells exceed the pool headroom at the table and delta: {over}.")


# ── Unit evaluation ──────────────────────────────────────────────────────────


class ConsistencyError(RuntimeError):
    """A unit's exposed measurement does not reproduce the analysis; the unit is not recorded."""


def injection_vector(dataset: Any, settings: Mapping[str, Any]) -> np.ndarray:
    """The raw-unit injection ``v_f = s_z · sd_f · w_f`` that D10's builder added, in joint feature order."""

    from motco.simulations.preprocessing import apply_unit_transform

    raw = {
        "methylation": apply_unit_transform("methylation", dataset.methylation.to_numpy(dtype=float), "beta_to_mvalue"),
        "expression": dataset.expression.to_numpy(dtype=float),
    }
    n_features = sum(raw[name].shape[1] for name in settings["layers"])
    w, _ = d10.feature_directions(n_features, settings)
    sd = np.concatenate([raw[name].std(axis=0, ddof=1) for name in settings["layers"]])
    return float(settings["injection"]["s_z"]) * sd * w


def top_k_indices(contrast: np.ndarray, feature_names: Sequence[str], k: int) -> np.ndarray:
    """The k largest ``|contrast|``, ties by feature name (the compact record's rule)."""

    order = np.lexsort((np.asarray(feature_names, dtype=object).astype(str), -np.abs(contrast)))
    return order[:k]


def contrast_metrics(
    contrast: np.ndarray | None,
    oracle_contrast: np.ndarray | None,
    *,
    truth: np.ndarray,
    injected: np.ndarray,
    direction: np.ndarray,
    bootstrap: Any,
    feature_names: Sequence[str],
    top_k: int,
) -> dict[str, Any]:
    """Every metric of one contrast vector; ``None`` where the contrast (or the oracle's) is unavailable."""

    record: dict[str, Any] = {name: None for name in METRICS}
    record["available"] = contrast is not None
    record["truth_size"] = int(truth.sum())
    record["valid_replicates"] = int(bootstrap.valid_replicates)
    if contrast is None:
        return record
    selected = top_k_indices(contrast, feature_names, top_k)
    record["truth_precision"] = (
        float(truth[selected].sum() / min(top_k, int(truth.sum()))) if truth.any() else None
    )
    record["contamination"] = float((injected[selected] & ~truth[selected]).sum() / len(selected))
    record["projection"] = float(contrast @ direction)
    sign = np.asarray(bootstrap.sign_stability, dtype=float)[selected]
    frequency = np.asarray(bootstrap.top_k_selection_frequency, dtype=float)[selected]
    record["sign_stability"] = float(np.nanmean(sign)) if np.isfinite(sign).any() else None
    record["top_k_frequency"] = float(np.nanmean(frequency)) if np.isfinite(frequency).any() else None
    if oracle_contrast is not None:
        norms = float(np.linalg.norm(contrast) * np.linalg.norm(oracle_contrast))
        record["cosine_oracle"] = float(contrast @ oracle_contrast / norms) if norms > 0 else None
        oracle_selected = top_k_indices(oracle_contrast, feature_names, top_k)
        record["overlap_oracle"] = float(len(set(selected) & set(oracle_selected)) / len(selected))
    return record


def _blocks(result: Any) -> dict[str, Any]:
    blocks: dict[str, Any] = {record.transition_id: record for record in result.transitions}
    blocks[result.principal_orientation.component_id] = result.principal_orientation
    return blocks


def arm_metrics(
    result: Any,
    oracle: Any,
    *,
    truth: Any,
    injected: np.ndarray,
    direction: np.ndarray,
    top_k: int,
) -> dict[str, Any]:
    """Per block and component, the metrics of one arm's attribution result."""

    names = list(result.feature_names)
    bootstrap = {(record.transition_id, record.component): record for record in result.bootstrap_records}
    oracle_blocks = _blocks(oracle)
    out: dict[str, Any] = {}
    for block_id, block in _blocks(result).items():
        mask = truth.global_mask if block_id == PRINCIPAL else truth.per_transition[block_id]
        out[block_id] = {
            component: contrast_metrics(
                getattr(block, component).directional_contrast,
                getattr(oracle_blocks[block_id], component).directional_contrast,
                truth=np.asarray(mask, dtype=bool),
                injected=injected,
                direction=direction,
                bootstrap=bootstrap[(block_id, component)],
                feature_names=names,
                top_k=top_k,
            )
            for component in COMPONENTS
        }
    return out


def check_consistency(result: Any) -> None:
    """The exposed estimator reproduces the latent scores, and the adjusted means project to the stage means."""

    artifacts = result.latent_artifacts
    scores = result.latent_scores.to_numpy(dtype=float)
    transformed = np.asarray(artifacts.model.transform(artifacts.features.to_numpy(dtype=float)), dtype=float)
    if not np.allclose(transformed, scores, rtol=0.0, atol=1e-10 * max(1.0, float(np.abs(scores).max()))):
        raise ConsistencyError("latent_scores: the exposed estimator does not reproduce the latent scores.")
    means = result.attribution.means
    projected = np.asarray(model_transform(artifacts.model, means[list(artifacts.features.columns)]), dtype=float)
    stage_means = result.stage_means.drop(columns=["group", "stage"]).to_numpy(dtype=float)
    labels_match = list(means["group"]) == list(result.stage_means["group"]) and list(means["stage"]) == list(
        result.stage_means["stage"]
    )
    if not labels_match or not np.allclose(
        projected, stage_means, rtol=1e-8, atol=1e-8 * max(1.0, float(np.abs(stage_means).max()))
    ):
        raise ConsistencyError("stage_means: the adjusted means do not project to the analysis's stage means.")


def model_transform(model: Any, frame: pd.DataFrame) -> np.ndarray:
    return np.asarray(model.transform(frame.to_numpy(dtype=float)), dtype=float)


def analysis_params(settings: Mapping[str, Any], seeds: Mapping[str, int], *, n_jobs: int = 1) -> Any:
    from motco.analysis import TrajectoryAnalysisParams
    from motco.simulations.evaluation import AttributionDiagnosticSettings

    attribution = settings["attribution"]
    return TrajectoryAnalysisParams(
        permutations=int(settings["permutations"]),
        seed=seeds["rrpp"],
        n_jobs=n_jobs,
        orientation_estimator=settings["orientation_estimator"],
        min_components=int(settings["min_components"]),
        residual_correction=settings["residual_correction"],
        attribution=AttributionDiagnosticSettings(
            enabled=True,
            bootstrap_replicates=int(attribution["bootstrap_replicates"]),
            bootstrap_seed=seeds["bootstrap"],
            top_k=int(attribution["top_k"]),
            zero_tolerance=float(attribution["zero_tolerance"]),
        ),
        **settings["cv"],
    )


def analyse(built: Mapping[str, Any], settings: Mapping[str, Any], *, n_jobs: int = 1) -> tuple[Any, Any]:
    """The ``adjusted`` arm: the unit's dataset and its covariate-adjusted analysis, consistency-checked."""

    from motco.analysis import OmicsBlock, TrajectoryDataset, analyze_trajectories

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
    result = analyze_trajectories(dataset, analysis_params(settings, built["seeds"], n_jobs=n_jobs))
    check_consistency(result)
    return dataset, result


def oracle_inputs(
    built: Mapping[str, Any], dataset: Any, result: Any, settings: Mapping[str, Any]
) -> tuple[pd.DataFrame, np.ndarray]:
    """The oracle's features ``X̃ − z ⊗ (v / σ)`` and the injection ``v / σ`` in analysed coordinates."""

    artifacts = result.latent_artifacts
    injection = injection_vector(built["dataset"], settings) / np.asarray(artifacts.original_scales, dtype=float)
    z = dataset.metadata["z"].to_numpy(dtype=float)
    return artifacts.features - np.outer(z, injection), injection


def evaluate_built(built: Mapping[str, Any], settings: Mapping[str, Any], *, n_jobs: int = 1) -> dict[str, Any]:
    """The three arms on one built unit: one analysis, its frozen model and features, three attributions."""

    from motco.simulations.attribution_diagnostics import derive_truth_driver_features
    from motco.stats.attribution import analyze_orientation_attribution

    dataset, result = analyse(built, settings, n_jobs=n_jobs)
    artifacts = result.latent_artifacts
    features = artifacts.features
    groups, stages = list(dataset.groups), list(dataset.stages)
    labels = pd.DataFrame(
        {
            "group": dataset.metadata["group"].astype(str).to_numpy(),
            "stage": dataset.metadata["stage"].astype(str).to_numpy(),
        },
        index=features.index,
    )
    oracle_features, injection = oracle_inputs(built, dataset, result, settings)
    seeds = built["seeds"]
    attribution = settings["attribution"]
    primitive = {
        "groups": groups,
        "stages": stages,
        "bootstrap_replicates": int(attribution["bootstrap_replicates"]),
        "bootstrap_seed": seeds["bootstrap"],
        "top_k": int(attribution["top_k"]),
        "zero_tolerance": float(attribution["zero_tolerance"]),
    }
    results = {
        "adjusted": result.attribution,
        "oracle": analyze_orientation_attribution(oracle_features, labels, artifacts.model, **primitive),
        "unadjusted": analyze_orientation_attribution(features, labels, artifacts.model, **primitive),
    }
    truth = derive_truth_driver_features(built["dataset"], list(features.columns), stages, settings["layers"])
    injected = injection != 0.0
    direction = injection / np.linalg.norm(injection)
    arms = {
        arm: {
            "blocks": arm_metrics(
                results[arm],
                results["oracle"],
                truth=truth,
                injected=injected,
                direction=direction,
                top_k=int(attribution["top_k"]),
            ),
            "rank_deficient_replicates": results[arm].config.rank_deficient_replicates,
            "mean_source": results[arm].config.mean_source,
        }
        for arm in ARMS
    }
    record = result.record
    return {
        "analysis": {
            "statistics": record["statistics"],
            "p_values": record["p_values"],
            "selected_lv": record["latent"]["selected_lv"],
            "cv_mean_auroc": record["latent"].get("cv_mean_auroc"),
        },
        "truth": {
            "n_drivers": int(truth.global_mask.sum()),
            "n_drivers_by_transition": {key: int(mask.sum()) for key, mask in truth.per_transition.items()},
            "n_injected": int(injected.sum()),
            "n_injected_drivers": int((injected & truth.global_mask).sum()),
        },
        "arms": arms,
    }


def build_unit(
    mode: str, effect: float, replicate: int, settings: Mapping[str, Any], reference: Any = None
) -> dict[str, Any]:
    built = cs.build_unit(mode, effect, replicate, settings, reference=reference)
    built["seeds"] = {**built["seeds"], "bootstrap": bootstrap_seed(int(settings["base_seed"]), replicate)}
    return built


def evaluate_unit(
    mode: str, effect: float, replicate: int, settings: Mapping[str, Any], *, n_jobs: int = 1, reference: Any = None
) -> dict[str, Any]:
    start = perf_counter()
    built = build_unit(mode, effect, replicate, settings, reference=reference)
    evaluated = evaluate_built(built, settings, n_jobs=n_jobs)
    return d10._jsonable(
        {
            "key": unit_key(mode, effect, replicate),
            "mode": mode,
            "effect": float(effect),
            "replicate": replicate,
            "settings_hash": settings_hash(settings),
            "seeds": built["seeds"],
            "group_labels": built["group_labels"],
            **evaluated,
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


def pending_units(out_dir: Path, settings: Mapping[str, Any], index: int, count: int) -> list[tuple[str, float, int]]:
    """This shard's units not yet recorded in any shard file; rewrites this shard's file with its complete lines."""

    path = shard_path(out_dir, index, count)
    path.parent.mkdir(parents=True, exist_ok=True)
    own = read_records([path]) if path.exists() else {}
    with path.open("w", encoding="utf-8") as handle:
        for key in sorted(own):
            handle.write(json.dumps(own[key], sort_keys=True, allow_nan=False) + "\n")
    others = [other for other in sorted(path.parent.glob("shard_*.jsonl")) if other != path]
    done = {**read_records(others), **own}
    return [
        unit
        for position, unit in enumerate(units(settings))
        if position % count == index and unit_key(*unit) not in done
    ]


def _evaluate_pending(unit: tuple[str, float, int], settings: Mapping[str, Any]) -> dict[str, Any]:
    return evaluate_unit(*unit, settings)


def run_shard(
    out_dir: Path,
    settings: Mapping[str, Any],
    index: int,
    count: int,
    *,
    n_jobs: int = 1,
    workers: int = 1,
    limit: int | None = None,
) -> int:
    """Evaluate this shard's missing units, appending one line each; returns how many were evaluated."""

    pending = pending_units(out_dir, settings, index, count)
    if limit is not None:
        pending = pending[:limit]
    from motco.simulations.reference import load_reference

    reference = load_reference()
    check_headroom(settings, reference)
    path = shard_path(out_dir, index, count)

    def append(position: int, record: Mapping[str, Any]) -> None:
        with path.open("a", encoding="utf-8") as handle:
            handle.write(json.dumps(record, sort_keys=True, allow_nan=False) + "\n")
            handle.flush()
        print(
            f"[shard {index}/{count}] {position}/{len(pending)} {record['key']} {record['runtime_seconds']:.1f}s",
            file=sys.stderr,
            flush=True,
        )

    if workers > 1:
        with ProcessPoolExecutor(max_workers=workers) as pool:
            for position, record in enumerate(
                pool.map(_evaluate_pending, pending, [settings] * len(pending)), start=1
            ):
                append(position, record)
    else:
        for position, unit in enumerate(pending, start=1):
            append(position, evaluate_unit(*unit, settings, n_jobs=n_jobs, reference=reference))
    return len(pending)


# ── Merge ────────────────────────────────────────────────────────────────────


METRIC_COLUMNS: tuple[str, ...] = (
    "mode",
    "effect",
    "replicate",
    "arm",
    "block",
    "component",
    "available",
    "truth_size",
    *METRICS,
    "rank_deficient_replicates",
    "angle_p_value",
    "selected_lv",
)
SUMMARY_COLUMNS: tuple[str, ...] = ("mode", "effect", "arm", "block", "component", "quantity", "n", "mean", "se")


def metric_rows(records: Mapping[str, Mapping[str, Any]]) -> list[dict[str, Any]]:
    """One row per unit × arm × block × component, sorted by unit key."""

    rows: list[dict[str, Any]] = []
    for key in sorted(records):
        record = records[key]
        analysis = record["analysis"]
        for arm in ARMS:
            arm_record = record["arms"][arm]
            for block, components in arm_record["blocks"].items():
                for component in COMPONENTS:
                    values = components[component]
                    rows.append(
                        {
                            "mode": record["mode"],
                            "effect": float(record["effect"]),
                            "replicate": int(record["replicate"]),
                            "arm": arm,
                            "block": block,
                            "component": component,
                            "available": bool(values["available"]),
                            "truth_size": values["truth_size"],
                            **{name: values.get(name) for name in METRICS},
                            "rank_deficient_replicates": arm_record.get("rank_deficient_replicates"),
                            "angle_p_value": analysis["p_values"].get("angle"),
                            "selected_lv": analysis["selected_lv"],
                        }
                    )
    return rows


def _mean_se(values: Sequence[Any]) -> tuple[int, float | None, float | None]:
    finite = [float(value) for value in values if value is not None and math.isfinite(float(value))]
    if not finite:
        return 0, None, None
    mean = float(np.mean(finite))
    se = float(np.std(finite, ddof=1) / math.sqrt(len(finite))) if len(finite) > 1 else None
    return len(finite), mean, se


def summary_rows(rows: Iterable[Mapping[str, Any]]) -> list[dict[str, Any]]:
    """Per cell × arm × block × component, mean and SE of every metric and of every paired oracle difference.

    A paired difference uses the units where both the arm's and the oracle's
    value are available.
    """

    grouped: dict[tuple[str, float, str, str], dict[str, dict[int, Mapping[str, Any]]]] = {}
    for row in rows:
        cell = (str(row["mode"]), float(row["effect"]), str(row["block"]), str(row["component"]))
        grouped.setdefault(cell, {}).setdefault(str(row["arm"]), {})[int(row["replicate"])] = row
    out: list[dict[str, Any]] = []
    for (mode, effect, block, component), arms in sorted(grouped.items()):
        oracle = arms.get("oracle", {})
        for arm in ARMS:
            by_replicate = arms.get(arm, {})
            base = {"mode": mode, "effect": effect, "arm": arm, "block": block, "component": component}
            for name in METRICS:
                n, mean, se = _mean_se([row.get(name) for row in by_replicate.values()])
                out.append({**base, "quantity": name, "n": n, "mean": mean, "se": se})
            if arm == "oracle":
                continue
            for name in PAIRED_METRICS:
                differences = [
                    float(row[name]) - float(oracle[replicate][name])
                    for replicate, row in by_replicate.items()
                    if replicate in oracle and row.get(name) is not None and oracle[replicate].get(name) is not None
                ]
                n, mean, se = _mean_se(differences)
                out.append({**base, "quantity": f"{name}_minus_oracle", "n": n, "mean": mean, "se": se})
    return out


def merge(out_dir: Path, settings: Mapping[str, Any]) -> dict[str, Any]:
    """Merge every shard and write the per-unit metrics and the per-cell summary."""

    records = read_records(sorted((out_dir / "records").glob("shard_*.jsonl")))
    stale = sorted(key for key, record in records.items() if record.get("settings_hash") != settings_hash(settings))
    if stale:
        raise SystemExit(f"{len(stale)} record(s) were written under other settings, e.g. {stale[:3]}.")
    expected = [unit_key(*unit) for unit in units(settings)]
    with (out_dir / "merged.jsonl").open("w", encoding="utf-8") as handle:
        for key in sorted(records):
            handle.write(json.dumps(records[key], sort_keys=True, allow_nan=False) + "\n")
    rows = metric_rows(records)
    d10._write_csv(out_dir / "adjusted_attribution_metrics.csv", rows, METRIC_COLUMNS)
    d10._write_csv(out_dir / "adjusted_attribution_summary.csv", summary_rows(rows), SUMMARY_COLUMNS)
    return {
        "units": len(records),
        "expected": len(expected),
        "missing": [key for key in expected if key not in records][:20],
    }


# ── Decision ─────────────────────────────────────────────────────────────────


def _lookup(summary: Iterable[Mapping[str, Any]]) -> dict[tuple, Mapping[str, Any]]:
    return {
        (str(row["mode"]), float(row["effect"]), str(row["arm"]), str(row["block"]), str(row["component"]),
         str(row["quantity"])): row
        for row in summary
    }


def _entry(
    criterion: str,
    cell: tuple[str, float],
    row: Mapping[str, Any] | None,
    *,
    direction: str,
    bound: float | None,
    k: float,
) -> dict[str, Any]:
    """One bound reading: ``upper`` holds when mean + k·SE ≤ bound, ``lower`` when mean − k·SE ≥ bound."""

    n = int(row["n"]) if row is not None else 0
    mean = None if row is None or row["mean"] in (None, "") else float(row["mean"])
    se = None if row is None or row["se"] in (None, "") else float(row["se"])
    entry: dict[str, Any] = {
        "criterion": criterion,
        "cell": cell_label(*cell),
        "n": n,
        "mean": mean,
        "se": se,
        "bound": bound,
        "direction": direction,
    }
    if mean is None or se is None or bound is None:
        return {**entry, "value": None, "verdict": "unavailable"}
    value = mean + k * se if direction == "upper" else mean - k * se
    holds = value <= bound if direction == "upper" else value >= bound
    return {**entry, "value": value, "verdict": "holds" if holds else "fails"}


def criteria_readings(
    summary: Iterable[Mapping[str, Any]],
    settings: Mapping[str, Any],
    *,
    arm: str = "adjusted",
    block: str = PRINCIPAL,
    component: str = "observed",
) -> dict[str, Any]:
    """The gates and the four criteria of one arm, block and component against the oracle."""

    lookup = _lookup(summary)
    k = float(settings["rule"]["se_multiplier"])
    margins = settings["rule"]["margins"]
    all_cells = cells(settings)
    orientation_cells = [cell for cell in all_cells if cell[0] != ANCHOR]

    def row(cell: tuple[str, float], which: str, quantity: str) -> Mapping[str, Any] | None:
        return lookup.get((cell[0], float(cell[1]), which, block, component, quantity))

    gates = [
        _entry(
            "contamination_gate",
            cell,
            row(cell, "unadjusted", "contamination_minus_oracle"),
            direction="lower",
            bound=float(margins["contamination_gate"]),
            k=k,
        )
        for cell in all_cells
        if cell[0] == ANCHOR
    ]
    for cell in all_cells:
        gate = _entry(
            "projection_gate",
            cell,
            row(cell, "unadjusted", "projection_minus_oracle"),
            direction="lower",
            bound=float(margins["projection_gate"]),
            k=k,
        )
        if gate["mean"] is not None and gate["se"] is not None:
            # The confounding's sign follows the group order; its size is what the gate reads.
            gate["value"] = abs(gate["mean"]) - k * gate["se"]
            gate["verdict"] = "holds" if gate["value"] >= gate["bound"] else "fails"
        gates.append(gate)
    top = max(orientation_cells, key=lambda cell: cell[1])
    oracle_precision = row(top, "oracle", "truth_precision")
    recoverable = None if oracle_precision is None or oracle_precision["mean"] in (None, "") else float(
        oracle_precision["mean"]
    )
    gates.append(
        {
            "criterion": "recoverability_gate",
            "cell": cell_label(*top),
            "n": int(oracle_precision["n"]) if oracle_precision is not None else 0,
            "mean": recoverable,
            "se": None,
            "bound": float(margins["recoverability_gate"]),
            "direction": "lower",
            "value": recoverable,
            "verdict": (
                "unavailable"
                if recoverable is None
                else ("holds" if recoverable >= float(margins["recoverability_gate"]) else "fails")
            ),
        }
    )
    criteria = [
        _entry(
            "contamination",
            cell,
            row(cell, arm, "contamination_minus_oracle"),
            direction="upper",
            bound=float(margins["contamination"]),
            k=k,
        )
        for cell in all_cells
    ]
    for cell in all_cells:
        adjusted = row(cell, arm, "projection_minus_oracle")
        unadjusted = row(cell, "unadjusted", "projection_minus_oracle")
        reference = None if unadjusted is None or unadjusted["mean"] in (None, "") else abs(float(unadjusted["mean"]))
        entry = _entry(
            "residual_confounding",
            cell,
            adjusted,
            direction="upper",
            bound=None if reference is None else float(margins["residual_confounding_ratio"]) * reference,
            k=k,
        )
        if entry["mean"] is not None and entry["se"] is not None and entry["bound"] is not None:
            entry["value"] = abs(entry["mean"]) + k * entry["se"]
            entry["verdict"] = "holds" if entry["value"] <= entry["bound"] else "fails"
        entry["unadjusted_mean"] = None if unadjusted is None or reference is None else float(unadjusted["mean"])
        criteria.append(entry)
    for name, quantity in (
        ("truth_precision", "truth_precision_minus_oracle"),
        ("sign_stability", "sign_stability_minus_oracle"),
        ("top_k_frequency", "top_k_frequency_minus_oracle"),
    ):
        criteria.extend(
            _entry(name, cell, row(cell, arm, quantity), direction="lower", bound=float(margins[name]), k=k)
            for cell in orientation_cells
        )
    return {"arm": arm, "block": block, "component": component, "gates": gates, "criteria": criteria}


def _failing(entries: Iterable[Mapping[str, Any]]) -> list[dict[str, str]]:
    return [
        {"criterion": str(entry["criterion"]), "cell": str(entry["cell"]), "verdict": str(entry["verdict"])}
        for entry in entries
        if entry["verdict"] != "holds"
    ]


def verdict_from_summary(summary: Iterable[Mapping[str, Any]], settings: Mapping[str, Any]) -> dict[str, Any]:
    """The gates, the criteria and the verdict on the observed principal contrast.

    The gates come first: a failed or unavailable gate makes the verdict
    ``construction_uninformative``. Otherwise every criterion must hold for
    ``adjusted_attribution_validated``; anything else is
    ``diagnose_adjusted_attribution`` naming the failing criteria and cells.
    """

    primary = criteria_readings(summary, settings)
    gate_failures = _failing(primary["gates"])
    criterion_failures = _failing(primary["criteria"])
    if gate_failures:
        verdict, failing = "construction_uninformative", gate_failures
    elif criterion_failures:
        verdict, failing = "diagnose_adjusted_attribution", criterion_failures
    else:
        verdict, failing = "adjusted_attribution_validated", []
    return {
        "verdict": verdict,
        "consequence": CONSEQUENCES[verdict],
        "failing": failing,
        "gates": primary["gates"],
        "criteria": primary["criteria"],
    }


def adjusted_attribution_decision(
    rows: Sequence[Mapping[str, Any]],
    settings: Mapping[str, Any],
    *,
    alpha: float = ALPHA,
) -> dict[str, Any]:
    """The verdict (:func:`verdict_from_summary`) on the per-unit metric rows, with every reading beside it."""

    summary = summary_rows(rows)
    verdict = verdict_from_summary(summary, settings)
    lookup = _lookup(summary)
    replicates = {
        cell_label(*cell): int(
            (lookup.get((cell[0], float(cell[1]), "adjusted", PRINCIPAL, "observed", "valid_replicates")) or {"n": 0})[
                "n"
            ]
        )
        for cell in cells(settings)
    }
    reasons = []
    if int(settings["permutations"]) < DEFAULT_PERMUTATIONS:
        reasons.append(f"{settings['permutations']} permutations < {DEFAULT_PERMUTATIONS}")
    if int(settings["attribution"]["bootstrap_replicates"]) < DEFAULT_BOOTSTRAP:
        reasons.append(f"B = {settings['attribution']['bootstrap_replicates']} < {DEFAULT_BOOTSTRAP}")
    short = {label: count for label, count in replicates.items() if count < DEFAULT_REPLICATES}
    if short:
        reasons.append(f"fewer than {DEFAULT_REPLICATES} available replicates: {short}")
    blocks = sorted({str(row["block"]) for row in rows if str(row["block"]) != PRINCIPAL})
    rejected = [row for row in rows if row.get("angle_p_value") is not None and float(row["angle_p_value"]) <= alpha]
    descriptive = {
        "unadjusted_arm": criteria_readings(summary, settings, arm="unadjusted"),
        "components": {
            component: criteria_readings(summary, settings, component=component)
            for component in ("pls_captured", "residual")
        },
        "transitions": {block: criteria_readings(summary, settings, block=block) for block in blocks},
        "angle_rejected": {
            "alpha": alpha,
            "units": {
                cell_label(*cell): len(
                    {
                        int(row["replicate"])
                        for row in rejected
                        if str(row["mode"]) == cell[0] and float(row["effect"]) == float(cell[1])
                    }
                )
                for cell in cells(settings)
            },
            "readings": criteria_readings(summary_rows(rejected), settings) if rejected else None,
        },
        "rank_deficient_replicates": {
            cell_label(*cell): _mean_se(
                [
                    row["rank_deficient_replicates"]
                    for row in rows
                    if str(row["mode"]) == cell[0]
                    and float(row["effect"]) == float(cell[1])
                    and row["arm"] == "adjusted"
                    and row["block"] == PRINCIPAL
                    and row["component"] == "observed"
                ]
            )[1]
            for cell in cells(settings)
        },
    }
    return {
        "check": CHECK,
        "settings_hash": settings_hash(settings),
        "adopted_calibration": settings.get("adopted_calibration"),
        "verdict": verdict["verdict"],
        "consequence": verdict["consequence"],
        "failing": verdict["failing"],
        "indicative": bool(reasons),
        "indicative_reasons": reasons,
        "replicates": replicates,
        "permutations": int(settings["permutations"]),
        "bootstrap_replicates": int(settings["attribution"]["bootstrap_replicates"]),
        "margins": dict(settings["rule"]["margins"]),
        "se_multiplier": float(settings["rule"]["se_multiplier"]),
        "gates": verdict["gates"],
        "criteria": verdict["criteria"],
        "descriptive": descriptive,
    }


def _read_metric_rows(path: Path) -> list[dict[str, Any]]:
    records = pd.read_csv(path).to_dict(orient="records")
    return [{str(key): (None if pd.isna(value) else value) for key, value in row.items()} for row in records]


def decide(out_dir: Path, settings: Mapping[str, Any]) -> dict[str, Any]:
    path = out_dir / "adjusted_attribution_metrics.csv"
    if not path.exists():
        raise SystemExit(f"{path} does not exist; run merge first.")
    decision = adjusted_attribution_decision(_read_metric_rows(path), settings)
    target = out_dir / "adjusted_attribution_decision.json"
    target.write_text(json.dumps(d10._jsonable(decision), indent=2, allow_nan=False) + "\n", encoding="utf-8")
    return decision


# ── Calibration ──────────────────────────────────────────────────────────────


def calibration_criteria(rows: Sequence[Mapping[str, Any]], settings: Mapping[str, Any]) -> dict[str, Any]:
    """Confounding and recoverability on the calibration cells' means.

    Confounding (rule version 2) is met when the unadjusted − oracle
    contamination at the anchor reaches the contamination gate margin and the
    |unadjusted − oracle projection| reaches the projection gate margin in both
    cells; recoverability when the oracle's truth precision at the largest
    orientation reaches its margin.
    """

    summary = _lookup(summary_rows(rows))
    margins = settings["rule"]["margins"]

    def mean_of(mode: str, effect: float, quantity: str) -> tuple[float | None, int]:
        row = summary.get((mode, float(effect), "unadjusted", PRINCIPAL, "observed", quantity))
        return (None if row is None or row["mean"] is None else float(row["mean"])), (0 if row is None else row["n"])

    contamination = {}
    projection = {}
    for mode, effect in cells(settings):
        label = cell_label(mode, effect)
        if mode == ANCHOR:
            mean, n = mean_of(mode, effect, "contamination_minus_oracle")
            contamination[label] = {
                "mean": mean,
                "n": n,
                "met": mean is not None and mean >= float(margins["contamination_gate"]),
            }
        mean, n = mean_of(mode, effect, "projection_minus_oracle")
        projection[label] = {
            "mean": mean,
            "n": n,
            "met": mean is not None and abs(mean) >= float(margins["projection_gate"]),
        }
    top = max((cell for cell in cells(settings) if cell[0] != ANCHOR), key=lambda cell: cell[1])
    precision = summary.get((top[0], float(top[1]), "oracle", PRINCIPAL, "observed", "truth_precision"))
    precision_mean = None if precision is None else precision["mean"]
    met = all(entry["met"] for entry in (*contamination.values(), *projection.values()))
    return {
        "rule_version": RULE_VERSION,
        "confounding": {
            "contamination_bound": float(margins["contamination_gate"]),
            "projection_bound": float(margins["projection_gate"]),
            "contamination": contamination,
            "projection": projection,
            "met": met,
        },
        "recoverability": {
            "bound": float(margins["recoverability_gate"]),
            "cell": cell_label(*top),
            "oracle_truth_precision": precision_mean,
            "met": precision_mean is not None and float(precision_mean) >= float(margins["recoverability_gate"]),
        },
    }


def calibration_entry(records: Mapping[str, Mapping[str, Any]], settings: Mapping[str, Any]) -> dict[str, Any]:
    criteria = calibration_criteria(metric_rows(records), settings)
    return {
        "s_z": float(settings["injection"]["s_z"]),
        "base_seed": int(settings["base_seed"]),
        "replicates": int(settings["replicates"]),
        "permutations": int(settings["permutations"]),
        "bootstrap_replicates": int(settings["attribution"]["bootstrap_replicates"]),
        "cells": [cell_label(*cell) for cell in cells(settings)],
        "units": len(records),
        **criteria,
        "met": bool(criteria["confounding"]["met"] and criteria["recoverability"]["met"]),
        "selected_rank": cs.rank_distribution(record["analysis"]["selected_lv"] for record in records.values()),
    }


def record_calibration(out_dir: Path, entry: Mapping[str, Any]) -> list[dict[str, Any]]:
    """Append the entry to ``calibration.json`` (never overwritten) and return the history."""

    out_dir.mkdir(parents=True, exist_ok=True)
    path = out_dir / "calibration.json"
    history = json.loads(path.read_text(encoding="utf-8")) if path.exists() else []
    history.append(d10._jsonable(entry))
    path.write_text(json.dumps(history, indent=1, allow_nan=False) + "\n", encoding="utf-8")
    return history


def calibration_settings(args: argparse.Namespace, s_z: float) -> dict[str, Any]:
    namespace = argparse.Namespace(**{**vars(args), "s_z": s_z, "calibration": None})
    settings = build_settings(namespace)
    top = max(ORIENTATION_EFFECTS)
    settings["cells"] = [[ANCHOR, 0.0], ["orientation", top]]
    return settings


def _calibration_unit(unit: tuple[str, float, int], settings: Mapping[str, Any]) -> dict[str, Any]:
    return evaluate_unit(*unit, settings)


def calibrate_once(settings: Mapping[str, Any], *, workers: int = 1) -> dict[str, Any]:
    from motco.simulations.reference import load_reference

    check_headroom(settings, load_reference())
    pending = units(settings)
    if workers > 1:
        with ProcessPoolExecutor(max_workers=workers) as pool:
            evaluated = list(pool.map(_calibration_unit, pending, [settings] * len(pending)))
    else:
        evaluated = [_calibration_unit(unit, settings) for unit in pending]
    return calibration_entry({record["key"]: record for record in evaluated}, settings)


def calibration_step(entry: Mapping[str, Any], s_z: float, max_s_z: float) -> tuple[str, float | None]:
    """What follows an iteration: ``adopt``, ``stop_recoverability``, ``stop_exhausted`` or ``raise``."""

    if entry["confounding"]["met"]:
        return ("adopt", None) if entry["recoverability"]["met"] else ("stop_recoverability", None)
    following = s_z + S_Z_STEP
    return ("raise", following) if following <= max_s_z + 1e-9 else ("stop_exhausted", None)


def calibrate(args: argparse.Namespace, run_once: Any = calibrate_once) -> dict[str, Any]:
    """Iterate ``s_z`` upward by 5 until confounding is met; stop on a recoverability failure."""

    s_z = float(args.s_z)
    while True:
        entry = dict(run_once(calibration_settings(args, s_z), workers=int(args.workers)))
        outcome, following = calibration_step(entry, s_z, float(args.max_s_z))
        entry["outcome"] = outcome
        record_calibration(args.out_dir, entry)
        print(json.dumps({key: entry[key] for key in ("s_z", "confounding", "recoverability", "met", "outcome")}),
              file=sys.stderr, flush=True)
        if outcome != "raise":
            return entry
        assert following is not None
        s_z = following


# ── CLI ──────────────────────────────────────────────────────────────────────


def add_construction_arguments(parser: argparse.ArgumentParser, *, calibration: bool = False) -> None:
    parser.add_argument("--out-dir", type=Path, required=True, help="Output directory")
    reps, perms, bootstrap = (
        (CALIBRATION_REPLICATES, CALIBRATION_PERMUTATIONS, CALIBRATION_BOOTSTRAP)
        if calibration
        else (DEFAULT_REPLICATES, DEFAULT_PERMUTATIONS, DEFAULT_BOOTSTRAP)
    )
    parser.add_argument("--reps", type=int, default=reps, help=f"Replicates per cell (default: {reps})")
    parser.add_argument("--perms", type=int, default=perms, help=f"RRPP permutations (default: {perms})")
    parser.add_argument("--bootstrap", type=int, default=bootstrap, help=f"Bootstrap replicates (default: {bootstrap})")
    parser.add_argument("--top-k", type=int, default=DEFAULT_TOP_K, help=f"Top-k size (default: {DEFAULT_TOP_K})")
    parser.add_argument("--base-seed", type=int, default=BASE_SEED, help=f"Base seed (default: {BASE_SEED})")
    parser.add_argument(
        "--s-z", type=float, default=DEFAULT_S_Z, help=f"Injection strength of z (default: {DEFAULT_S_Z:g})"
    )


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Covariate-adjusted attribution check at the case study's configuration (no R)."
    )
    sub = parser.add_subparsers(dest="command", required=True)
    cal = sub.add_parser("calibrate", help="Raise s_z until the covariate contaminates unadjusted attribution")
    add_construction_arguments(cal, calibration=True)
    cal.add_argument("--max-s-z", type=float, default=DEFAULT_MAX_S_Z, help="Largest s_z tried (default: 50)")
    cal.add_argument("--workers", type=int, default=1, help="Units evaluated in parallel")
    run = sub.add_parser("run", help="Evaluate (a shard of) the units")
    add_construction_arguments(run)
    run.add_argument("--calibration", type=Path, default=None, help="calibration.json whose iteration is adopted")
    run.add_argument("--shard", type=parse_shard, default=(0, 1), help="Run shard i of N (i/N; default 0/1)")
    run.add_argument("--n-jobs", type=int, default=1, help="Jobs for CV and RRPP (not part of the settings)")
    run.add_argument("--workers", type=int, default=1, help="Units evaluated in parallel (not part of the settings)")
    run.add_argument("--limit", type=int, default=None, help="Evaluate at most this many pending units")
    for name, text in (("merge", "Merge every shard into the metric tables"), ("decide", "Write the decision")):
        command = sub.add_parser(name, help=f"{text} (reads settings.json)")
        command.add_argument("--out-dir", type=Path, required=True, help="Output directory")
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    if args.command == "calibrate":
        entry = calibrate(args)
        print(json.dumps(entry, indent=1))
        return 0 if entry["outcome"] == "adopt" else 3
    if args.command == "run":
        settings = build_settings(args)
        ensure_settings(args.out_dir, settings)
        index, count = args.shard
        evaluated = run_shard(
            args.out_dir, settings, index, count, n_jobs=args.n_jobs, workers=args.workers, limit=args.limit
        )
        print(f"evaluated {evaluated} unit(s) into {shard_path(args.out_dir, index, count)}", file=sys.stderr)
        return 0
    settings = load_settings(args.out_dir)
    if args.command == "merge":
        print(json.dumps(merge(args.out_dir, settings), indent=1))
        return 0
    decision = decide(args.out_dir, settings)
    print(json.dumps({key: decision[key] for key in ("verdict", "failing", "indicative", "replicates")}, indent=1))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
