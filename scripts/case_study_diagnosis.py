#!/usr/bin/env python3
"""Diagnosis of the case-study configuration's size excess (no R).

The diagnosis behind ``openspec/changes/case-study-configuration-diagnosis/``
(spec ``case-study-configuration-diagnosis``). The case-study configuration
Type I check (:mod:`case_study_type_i_check`, 2026-10-29) returned
``diagnose_configuration_excess``: the covariate-adjusted model (the case
study's primary) exceeded on translation ``delta`` at δ = 0.1, and the
``stage`` latent-adjusted fallback exceeded on anchor ``angle`` and ``shape``.
This script reruns that check's exact data under arms that separate the
candidate causes, and issues predeclared verdicts.

**Parent binding.** Every run is bound to one committed run directory of the
parent check (``--parent``). Its ``settings.json`` must pass its recorded hash,
and :func:`case_study_type_i_check.build_settings` must rebuild it from the
recorded arguments to the same hash (so construction code that drifted since
the parent run is refused). The diagnosis's own ``settings.json`` carries the
parent's settings, path and hash; a resume against another parent is refused.
``--reps`` takes a prefix of the parent's replicates; ``--perms`` below the
parent's marks the run as a **smoke** run (reproduction ``not_compared``,
verdicts indicative). Every unit (the anchor or the translation control at
e = 1.0, one replicate) regenerates the parent's dataset, covariates, CV
settings and RRPP seed.

**Arms** (all on the unit's data and RRPP seed):

- the **trajectory-model factorial** on the covariate-blind basis. One
  :func:`motco.analysis.analyze_trajectories` call with the six covariates
  (``z``, ``x1`` … ``x5``) fits the basis and gives the record; the script's
  inference path (the analysis's recoded design, ``TrajectoryDataset``
  covariate designs, ``_covariate_adjusted_design``, ``estimate_difference``
  and serial ``RRPP``) then runs ``z_noise{0,2,5}`` on the latent scores ``Y``
  (``z`` estimated) and ``oracle_noise{0,2,5}`` on ``Y0``, the same scores with
  ``z``'s injected effect removed exactly. The noise columns are prefixes of
  the parent's draw, so the arms are nested; ``p_red`` (reduced-model columns)
  is 4, 6, 9 for the oracle arms and 5, 7, 10 for the ``z`` arms. The script's
  ``z_noise5`` must equal the analysis record's statistics and p-values on every
  unit, or the unit fails;
- the **oracle**: ``X0`` is D10's build at ``s_z = 0`` (same seeds), ``X − X0``
  must be ``z`` times a fixed vector per block, the fitted covariate-blind map
  (block scalers, the PLS centring and ``ddof = 1`` scaling, the PLS rotations)
  is reconstructed and must reproduce the analysis's latent scores to 1e-9, and
  ``Y0`` is ``X0`` through that map (``Y − Y0`` is then ``z`` times a fixed
  latent vector). A failed check is a hard error;
- the **residual-df correction**, one per factorial arm (``<arm>_corrected``):
  RRPP on ``ŷ_red + κ·e_red`` with κ = √(n / (n − p_red)), the arm's models and
  the unit's seed (so the same permutation indices), the uncorrected observed
  statistics, and p-values by the analysis's rule;
- the **latent-adjustment arms**, each its own analysis with the six covariates
  and ``covariate_latent_adjustment = "stage"``: ``latent_stage_cv`` (CV rank at
  the floor of 2, the parent's arm), ``latent_stage_floor1`` (floor 1;
  descriptive) and ``latent_stage_rank2`` (forced rank 2, no CV). An analysis
  error makes that arm unavailable for the unit; the others proceed.

**Reproduction.** ``z_noise5`` and ``latent_stage_cv`` are compared with the
parent's committed per-unit p-values (as k/(B + 1) integers); any mismatch
gates every verdict (``reproduction_failed``).

**Contrasts.** Paired on the same units: per replicate ``r``, ``d_r`` is the
number of readings (run × control) arm B rejects and arm A does not, minus the
reverse, and the exact one-sided p-value is P(Σ ε_r d_r ≥ Σ d_r) over
independent equiprobable signs ε_r (each replicate's readings flip together).
Predeclared on ``delta``: C1 ``oracle_noise5`` vs ``oracle_noise0``, C2
``z_noise5`` vs ``z_noise0``, C3 ``z_noise0`` vs ``oracle_noise0``, C4
``z_noise5_corrected`` vs ``z_noise0_corrected``; the same on ``angle`` and
``shape``, and per reading, are descriptive.

**Verdicts** (pure, :func:`diagnosis_decision`), behind the reproduction gate:

- the covariate-adjusted ``delta``: ``residual_df_shrinkage`` when C1 or C2 is
  significant at one-sided 0.05, ``z_noise5_corrected`` holds the parent's
  one-sided bound (α + k·√(α(1 − α)/R)) on every statistic at both controls in
  both runs and C4 is not significant; else ``collinear_covariate`` when C3 is
  significant; else ``unresolved`` (next step: the R = 1,000 replication);
- the ``stage`` latent adjustment's ``angle``/``shape``: ``rank_selection``
  when ``latent_stage_rank2`` holds the bound on every statistic at both
  controls in both runs, else ``residualized_basis`` naming each exceeding run,
  control and statistic. The floor-1 arm's rates and rank-1 fraction are
  descriptive.

Each verdict carries its consequence for the case study and for MOTCO.

Outputs (one ``--out-dir`` per parent run):

- ``settings.json``: the diagnosis settings (with the parent's) and their hash;
- ``records/shard_<i>_of_<N>.jsonl``: one line per unit (every arm, the seeds,
  the reconstruction error and the rank-one residuals); resume skips units
  recorded in any shard file and re-runs torn lines;
- after ``merge`` (automatic for an unsharded run): ``merged.jsonl``,
  ``diagnosis_pvalues.csv``, ``diagnosis_rates.csv`` (per run × control × arm ×
  statistic: R, rejections, rate, bound, verdict), ``diagnosis_rank_table.csv``
  (anchor units by selected rank 1, 2, 3, ≥ 4 for ``z_noise5`` and the latent
  arms) and ``diagnosis_reproduction.json``;
- ``decide`` (across runs, into its own ``--out-dir``): ``diagnosis_contrasts.csv``
  and ``case_study_diagnosis_decision.json``.

Examples:

    # smoke profile against one parent (99 permutations: not compared, indicative)
    python scripts/case_study_diagnosis.py run --parent results/case-study-type-i-2026-10-29/delta0.1 \\
        --reps 20 --perms 99 --out-dir results/case-study-diagnosis/smoke-delta0.1

    # one SLURM array task, then the merge
    python scripts/case_study_diagnosis.py run --parent results/case-study-type-i-2026-10-29/delta0.1 \\
        --shard 3/50 --out-dir results/case-study-diagnosis/delta0.1
    python scripts/case_study_diagnosis.py merge --out-dir results/case-study-diagnosis/delta0.1

    # the verdicts over both runs
    python scripts/case_study_diagnosis.py decide --run results/case-study-diagnosis/delta0.1 \\
        --run results/case-study-diagnosis/delta0.05 --out-dir results/case-study-diagnosis
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import sys
from pathlib import Path
from time import perf_counter
from typing import Any, Iterable, Mapping, Sequence

import numpy as np
import pandas as pd

from motco import analysis
from motco.analysis import OmicsBlock, TrajectoryAnalysisError, TrajectoryAnalysisParams, TrajectoryDataset
from motco.simulations.evaluation import (
    _empirical_p_value,
    _extract_null_distributions,
    _extract_pair_statistics,
    _summarize_null_distributions,
    build_simulation_trajectory_design,
)
from motco.simulations.preprocessing import concatenate_blocks, fit_block_scalers, transform_blocks
from motco.stats.permutation import RRPP
from motco.stats.trajectory import estimate_betas, estimate_difference

sys.path.insert(0, str(Path(__file__).resolve().parent))
import case_study_type_i_check as cs  # noqa: E402  (the parent check: construction, settings, sharding)

sys.path.pop(0)

d10 = cs.d10
SCRIPT_VERSION = 1
CHECK = "case_study_diagnosis"
REPO_ROOT = Path(__file__).resolve().parents[1]
STATISTICS: tuple[str, ...] = cs.STATISTICS
CONTROL_MODES: tuple[str, ...] = cs.CONTROL_MODES
ANCHOR = cs.ANCHOR
#: The parent's declared covariates, in model order; ``z_noise5`` is exactly that model.
PARENT_COVARIATES: tuple[str, ...] = ("z", "x1", "x2", "x3", "x4", "x5")

#: Factorial arm → (``z``'s injection removed exactly from the latent scores, noise covariates in the model).
FACTORIAL: dict[str, tuple[bool, int]] = {
    "oracle_noise0": (True, 0),
    "oracle_noise2": (True, 2),
    "oracle_noise5": (True, 5),
    "z_noise0": (False, 0),
    "z_noise2": (False, 2),
    "z_noise5": (False, 5),
}
FACTORIAL_ARMS: tuple[str, ...] = tuple(FACTORIAL)
CORRECTED_SUFFIX = "_corrected"
CORRECTED_ARMS: tuple[str, ...] = tuple(f"{arm}{CORRECTED_SUFFIX}" for arm in FACTORIAL_ARMS)
#: Latent arm → its ``TrajectoryAnalysisParams`` rank settings (all with the six covariates and ``stage``).
LATENT_ARMS: dict[str, dict[str, int]] = {
    "latent_stage_cv": {"min_components": 2},
    "latent_stage_floor1": {"min_components": 1},
    "latent_stage_rank2": {"min_components": 2, "forced_components": 2},
}
LATENT_ADJUSTMENT = "stage"
ARMS: tuple[str, ...] = (*FACTORIAL_ARMS, *CORRECTED_ARMS, *LATENT_ARMS)
CASE_STUDY_ARM = "z_noise5"
#: Diagnosis arm → the parent arm whose committed p-values it must reproduce.
REPRODUCED: dict[str, str] = {"z_noise5": "adjusted", "latent_stage_cv": "latent_stage"}
RANK_TABLE_ARMS: tuple[str, ...] = (CASE_STUDY_ARM, *LATENT_ARMS)
RANK_BUCKETS: tuple[str, ...] = ("1", "2", "3", ">=4")
#: Predeclared contrasts (label, arm B, arm A): does B reject more often than A?
CONTRASTS: tuple[tuple[str, str, str], ...] = (
    ("C1", "oracle_noise5", "oracle_noise0"),
    ("C2", "z_noise5", "z_noise0"),
    ("C3", "z_noise0", "oracle_noise0"),
    ("C4", "z_noise5_corrected", "z_noise0_corrected"),
)
JUDGED_STATISTIC = "delta"
CONTRAST_ALPHA = 0.05
#: sklearn's ``PLSRegression(scale=True)`` divides the centred features by their ``ddof=1`` SD.
PLS_SCALE_DDOF = 1
RECONSTRUCTION_TOLERANCE = 1e-9
RANK_ONE_TOLERANCE = 1e-9
ADJUSTED_VERDICTS: tuple[str, ...] = ("residual_df_shrinkage", "collinear_covariate", "unresolved")
LATENT_VERDICTS: tuple[str, ...] = ("rank_selection", "residualized_basis")
DECISIONS: tuple[str, ...] = ("diagnosed", "reproduction_failed")
CONSEQUENCES: dict[str, str] = {
    "residual_df_shrinkage": (
        "a MOTCO change offering the residual-df corrected permutation (opt-in, byte-identical when absent), "
        "validated by its own Type I and power check before the case study may use it; then a sead_sexdiff "
        "re-issue"
    ),
    "collinear_covariate": (
        "the case study keeps reading delta against the realized size; the next diagnostic is the confounding "
        "strength (a_gs) at one covariate column"
    ),
    "unresolved": "the R = 1,000 replication of the case-study configuration Type I check is the next step",
    "rank_selection": (
        "the stage latent adjustment is admissible at a fixed rank of 2 at this configuration; sead_sexdiff may "
        "adopt it in its own change"
    ),
    "residualized_basis": "the stage latent adjustment stays inadmissible at this configuration",
}
NEXT_STEPS: dict[str, str] = {
    "residual_df_shrinkage": "MOTCO change: the opt-in residual-df corrected permutation with its own validation",
    "collinear_covariate": "the confounding-strength (a_gs) diagnostic at one covariate column",
    "unresolved": "the R = 1,000 replication",
}


# ── Parent binding and settings ──────────────────────────────────────────────


def _display_path(path: Path) -> str:
    """The path relative to the repository root when inside it (identical on every clone), else absolute."""

    resolved = Path(path).resolve()
    try:
        return resolved.relative_to(REPO_ROOT).as_posix()
    except ValueError:
        return str(resolved)


def _resolve_path(text: str) -> Path:
    path = Path(text)
    return path if path.is_absolute() else REPO_ROOT / path


def parent_arguments(parent: Mapping[str, Any]) -> argparse.Namespace:
    """The parent check's construction arguments, read back from its recorded settings."""

    return argparse.Namespace(
        reps=parent["replicates"],
        perms=parent["permutations"],
        base_seed=parent["base_seed"],
        alpha=parent["alpha"],
        k=parent["k"],
        a_g=parent["covariates"]["numeric"]["a_g"],
        a_s=parent["covariates"]["numeric"]["a_s"],
        a_gs=parent["covariates"]["numeric"]["a_gs"],
        s_z=parent["injection"]["s_z"],
        feature_fraction=parent["injection"]["feature_fraction"],
        delta=parent["generator"][cs.STRENGTH_FIELDS[0]],
        noise_covariates=parent["covariates"]["noise_numeric"]["count"],
        auroc_band=tuple(parent["calibration"]["auroc_band"]),
        rho_threshold=parent["calibration"]["rho_threshold"],
        translation_effect=parent["controls"]["translation"],
    )


def load_parent(parent_dir: Path) -> dict[str, Any]:
    """The parent run's settings and hash; refuses a corrupt or unreproducible parent.

    The settings file must pass its recorded hash, and the current construction
    code must rebuild the same settings (hence the same hash) from the recorded
    arguments.
    """

    path = Path(parent_dir) / "settings.json"
    if not path.exists():
        raise SystemExit(f"parent {parent_dir}: {path} does not exist.")
    stored = json.loads(path.read_text(encoding="utf-8"))
    settings, recorded = stored.get("settings"), stored.get("hash")
    if not isinstance(settings, dict) or cs.settings_hash(settings) != recorded:
        raise SystemExit(f"parent {parent_dir}: settings.json fails its recorded hash {recorded}; refusing.")
    if settings.get("check") != cs.CHECK:
        raise SystemExit(f"parent {parent_dir}: not a {cs.CHECK} run (check {settings.get('check')!r}); refusing.")
    try:
        rebuilt = cs.build_settings(parent_arguments(settings))
    except (KeyError, TypeError, ValueError) as exc:
        raise SystemExit(f"parent {parent_dir}: its settings cannot be rebuilt ({exc!r}); refusing.") from exc
    if cs.settings_hash(rebuilt) != recorded:
        differing = sorted(key for key in set(settings) | set(rebuilt) if settings.get(key) != rebuilt.get(key))
        raise SystemExit(
            f"parent {parent_dir}: the current construction code rebuilds its settings to hash "
            f"{cs.settings_hash(rebuilt)}, not the recorded {recorded} (differing: {', '.join(differing)}); "
            "refusing."
        )
    if tuple(cs.declared_covariates(settings)) != PARENT_COVARIATES:
        raise SystemExit(
            f"parent {parent_dir}: declares covariates {cs.declared_covariates(settings)}, not "
            f"{list(PARENT_COVARIATES)}; the factorial needs the six-column model."
        )
    return {"settings": settings, "hash": recorded}


def build_settings(
    parent_dir: Path, parent: Mapping[str, Any], *, reps: int | None = None, perms: int | None = None
) -> dict[str, Any]:
    """The diagnosis settings: the parent's (with its path and hash), the replicate prefix and permutations."""

    settings = parent["settings"]
    replicates = int(settings["replicates"]) if reps is None else int(reps)
    permutations = int(settings["permutations"]) if perms is None else int(perms)
    if not 1 <= replicates <= int(settings["replicates"]):
        raise SystemExit(f"--reps must be a prefix of the parent's {settings['replicates']} replicates.")
    if not 1 <= permutations <= int(settings["permutations"]):
        raise SystemExit(f"--perms must be between 1 and the parent's {settings['permutations']}.")
    return {
        "script_version": SCRIPT_VERSION,
        "check": CHECK,
        "parent": {
            "label": Path(parent_dir).resolve().name,
            "path": _display_path(Path(parent_dir)),
            "hash": parent["hash"],
            "settings": settings,
        },
        "replicates": replicates,
        "permutations": permutations,
        "smoke": permutations < int(settings["permutations"]),
        "alpha": float(settings["alpha"]),
        "k": float(settings["k"]),
        "arms": {
            "factorial": {
                arm: {
                    "latent_scores": "Y0 (z's injection removed exactly)" if oracle else "Y (covariate-blind basis)",
                    "covariates": arm_covariates(arm),
                }
                for arm, (oracle, _) in FACTORIAL.items()
            },
            "corrected": "each factorial arm with RRPP on yhat_red + kappa * e_red, kappa = sqrt(n / (n - p_red))",
            "latent": {
                arm: {"covariates": list(PARENT_COVARIATES), "latent_adjustment": LATENT_ADJUSTMENT, **kwargs}
                for arm, kwargs in LATENT_ARMS.items()
            },
            "reproduced": dict(REPRODUCED),
        },
        "oracle": {
            "features": "D10 build at s_z = 0 with the same seeds; X - X0 rank one in z per block",
            "projection": "analysis block scalers, PLS centring and ddof=1 scaling, pls_loadings as rotations",
            "reconstruction_tolerance": RECONSTRUCTION_TOLERANCE,
            "rank_one_tolerance": RANK_ONE_TOLERANCE,
        },
        "contrasts": {
            "predeclared": [{"label": label, "b": b, "a": a} for label, b, a in CONTRASTS],
            "statistic": JUDGED_STATISTIC,
            "test": "exact one-sided replicate-clustered sign flip on discordant rejections",
            "alpha": CONTRAST_ALPHA,
        },
    }


def ensure_settings(out_dir: Path, settings: Mapping[str, Any]) -> None:
    """Write the settings header, or refuse a directory written against another parent or other settings."""

    path = Path(out_dir) / "settings.json"
    if path.exists():
        stored = json.loads(path.read_text(encoding="utf-8"))
        written = (stored.get("settings") or {}).get("parent") or {}
        if written.get("hash") != settings["parent"]["hash"]:
            raise SystemExit(
                f"{out_dir} was written against parent {written.get('path')} (hash {written.get('hash')}); "
                f"refusing to resume it against {settings['parent']['path']} (hash {settings['parent']['hash']})."
            )
    d10.ensure_settings(Path(out_dir), settings)


def parent_view(settings: Mapping[str, Any]) -> dict[str, Any]:
    """The parent's settings at the diagnosis's replicates and permutations (what the construction reads)."""

    return {
        **settings["parent"]["settings"],
        "replicates": int(settings["replicates"]),
        "permutations": int(settings["permutations"]),
    }


def units(settings: Mapping[str, Any]) -> list[tuple[str, float, int]]:
    return cs.units(parent_view(settings))


def arm_covariates(arm: str) -> list[str]:
    """A factorial arm's covariates in model order: ``z`` unless oracle, then the first k noise covariates."""

    oracle, noise = FACTORIAL[arm.removesuffix(CORRECTED_SUFFIX)]
    return ([] if oracle else ["z"]) + list(PARENT_COVARIATES[1 : 1 + noise])


# ── Construction and the oracle ──────────────────────────────────────────────


def _without_injection(parent: Mapping[str, Any]) -> dict[str, Any]:
    return {**parent, "injection": {**parent["injection"], "s_z": 0.0}}


def build_unit(
    mode: str, effect: float, replicate: int, settings: Mapping[str, Any], reference: Any = None
) -> dict[str, Any]:
    """The parent's unit, plus the same unit built without the covariate injection (``X0``)."""

    parent = settings["parent"]["settings"]
    built = cs.build_unit(mode, effect, replicate, parent, reference=reference)
    clean = cs.build_unit(mode, effect, replicate, _without_injection(parent), reference=reference)
    if built["seeds"] != clean["seeds"] or not built["metadata"].equals(clean["metadata"]):
        raise RuntimeError(f"unit {cs.unit_key(mode, effect, replicate)}: the s_z = 0 build drew other data.")
    return {**built, "blocks_without_injection": clean["blocks"]}


def rank_one_residual(difference: np.ndarray, z: np.ndarray) -> float:
    """Max |D − z cᵀ| over max(1, max |D|), with c the least-squares fit of D on z (no intercept)."""

    difference = np.asarray(difference, dtype=float)
    coefficients = z @ difference / float(z @ z)
    residual = difference - np.outer(z, coefficients)
    return float(np.max(np.abs(residual)) / max(1.0, float(np.max(np.abs(difference)))))


def injection_residuals(built: Mapping[str, Any], layers: Sequence[str]) -> dict[str, float]:
    """Per block, the rank-one residual of ``X − X0`` in ``z``; above tolerance is a hard error."""

    z = built["metadata"]["z"].to_numpy(dtype=float)
    residuals = {
        name: rank_one_residual(
            built["blocks"][name].to_numpy(dtype=float) - built["blocks_without_injection"][name].to_numpy(dtype=float),
            z,
        )
        for name in layers
    }
    over = {name: value for name, value in residuals.items() if not value <= RANK_ONE_TOLERANCE}
    if over:
        raise RuntimeError(f"X - X0 is not rank one in z: residual {over} above {RANK_ONE_TOLERANCE}.")
    return residuals


class CovariateBlindProjection:
    """The analysis's fitted covariate-blind map from features to latent scores.

    Block scalers fitted on the unit's features (the analysis's own calls),
    then the PLS centring and ``ddof = 1`` scaling of the standardized matrix,
    then ``pls_loadings`` (``x_rotations_``) as the rotations.
    """

    def __init__(self, blocks: Mapping[str, pd.DataFrame], loadings: pd.DataFrame, layers: Sequence[str]) -> None:
        self.layers = list(layers)
        self.transforms = dict.fromkeys(self.layers, "identity")
        frames = {name: blocks[name] for name in self.layers}
        self.scalers = fit_block_scalers(frames, self.transforms)
        standardized = self._standardized(frames)
        features = [f"{block}__{feature}" for block, feature in zip(loadings["block"], loadings["feature"])]
        if features != list(standardized.columns):
            raise RuntimeError("pls_loadings rows do not follow the analysis's joint feature order.")
        values = standardized.to_numpy(dtype=float)
        self.x_mean = values.mean(axis=0)
        x_std = (values - self.x_mean).std(axis=0, ddof=PLS_SCALE_DDOF)
        x_std[x_std == 0.0] = 1.0
        self.x_std = x_std
        self.components = [column for column in loadings.columns if column not in ("block", "feature")]
        self.rotations = loadings[self.components].to_numpy(dtype=float)

    def _standardized(self, frames: Mapping[str, pd.DataFrame]) -> pd.DataFrame:
        return concatenate_blocks(transform_blocks(frames, self.scalers, self.transforms), self.layers)

    def __call__(self, blocks: Mapping[str, pd.DataFrame]) -> np.ndarray:
        values = self._standardized({name: blocks[name] for name in self.layers}).to_numpy(dtype=float)
        return ((values - self.x_mean) / self.x_std) @ self.rotations


def oracle_latent(
    built: Mapping[str, Any], result: Any, layers: Sequence[str]
) -> tuple[pd.DataFrame, dict[str, float]]:
    """``Y0``: the injection-free features through the reconstructed fitted map, with its checks.

    The reconstruction must reproduce the analysis's latent scores to
    ``RECONSTRUCTION_TOLERANCE`` and ``Y − Y0`` must be rank one in ``z``;
    either failing is a hard error.
    """

    scores = result.latent_scores
    project = CovariateBlindProjection(built["blocks"], result.pls_loadings, layers)
    error = float(np.max(np.abs(project(built["blocks"]) - scores.to_numpy(dtype=float))))
    if not error <= RECONSTRUCTION_TOLERANCE:
        raise RuntimeError(
            f"the reconstructed covariate-blind projection misses the analysis's latent scores by {error:.3g} "
            f"(tolerance {RECONSTRUCTION_TOLERANCE})."
        )
    clean = pd.DataFrame(project(built["blocks_without_injection"]), index=scores.index, columns=scores.columns)
    z = built["metadata"]["z"].to_numpy(dtype=float)
    shift = rank_one_residual(scores.to_numpy(dtype=float) - clean.to_numpy(dtype=float), z)
    if not shift <= RANK_ONE_TOLERANCE:
        raise RuntimeError(f"Y - Y0 is not rank one in z: residual {shift:.3g} above {RANK_ONE_TOLERANCE}.")
    return clean, {"reconstruction_max_abs_error": error, "latent_shift_rank_one_residual": shift}


# ── Inference path ───────────────────────────────────────────────────────────


def working_frame(dataset: TrajectoryDataset) -> pd.DataFrame:
    """Group and stage recoded to ordinal strings, exactly as ``analyze_trajectories`` recodes them."""

    group_codes = {label: str(index) for index, label in enumerate(dataset.groups)}
    stage_codes = {label: str(index) for index, label in enumerate(dataset.stages)}
    return pd.DataFrame(
        {
            "group": dataset.metadata[dataset.group_col].astype(str).map(group_codes).to_numpy(),
            "stage": dataset.metadata[dataset.stage_col].astype(str).map(stage_codes).to_numpy(),
        }
    )


def trajectory_design(dataset: TrajectoryDataset) -> Any:
    return build_simulation_trajectory_design(working_frame(dataset), group_col="group", stage_col="stage")


def rescaled_residuals(Y: Any, model_reduced: np.ndarray, kappa: float) -> np.ndarray:
    """``ŷ_red + κ·e_red``: RRPP refits it to ``ŷ_red`` and permutes ``κ·e_red``."""

    values = np.asarray(Y, dtype=float)
    fitted = np.asarray(model_reduced, dtype=float) @ np.asarray(estimate_betas(model_reduced, values), dtype=float)
    return fitted + kappa * (values - fitted)


def residual_df_kappa(n: int, p_red: int) -> float:
    return math.sqrt(n / (n - p_red))


def infer(
    Y: Any,
    design: Any,
    covariates: Any,
    *,
    permutations: int,
    seed: int,
    orientation_estimator: str,
    kappa: float | None = None,
) -> dict[str, Any]:
    """Statistics and RRPP p-values of one trajectory model, by the analysis's rule.

    With ``kappa`` the null permutes the reduced model's residuals scaled by
    ``kappa`` (the observed statistics stay those of ``Y``).
    """

    full, reduced, ls_means = analysis._covariate_adjusted_design(
        design.model_full, design.model_reduced, design.ls_means, covariates
    )
    shape_available = len(design.stage_levels) >= 3
    deltas, angles, shapes = estimate_difference(
        Y, full, ls_means, design.contrast, orientation_estimator=orientation_estimator
    )
    observed = _extract_pair_statistics(deltas, angles, shapes, shape_available=shape_available)
    null_Y = Y if kappa is None else rescaled_residuals(Y, reduced, kappa)
    dist_delta, dist_angle, dist_shape = RRPP(
        null_Y,
        full,
        reduced,
        ls_means,
        design.contrast,
        permutations=permutations,
        n_jobs=1,
        progress=False,
        seed=seed,
        orientation_estimator=orientation_estimator,
    )
    null = _extract_null_distributions(dist_delta, dist_angle, dist_shape, shape_available=shape_available)
    summary = _summarize_null_distributions(null)
    return {
        "statistics": observed,
        "p_values": {
            statistic: _empirical_p_value(null[statistic], value)
            for statistic, value in observed.items()
            if statistic in null and np.isfinite(value)
        },
        "null_q95": {statistic: entry.get("q95") for statistic, entry in summary.items()},
        "p_red": int(reduced.shape[1]),
        "kappa": kappa,
    }


# ── Unit evaluation ──────────────────────────────────────────────────────────


ARM_FIELDS: tuple[str, ...] = (
    "statistics",
    "p_values",
    "null_q95",
    "selected_lv",
    "component_selection",
    "cv_mean_auroc",
    "p_red",
    "kappa",
    "covariate_design_columns",
    "error",
)


def _failed_arm(error: str) -> dict[str, Any]:
    failed: dict[str, Any] = dict.fromkeys(ARM_FIELDS)
    failed.update({"statistics": {}, "p_values": {}, "null_q95": {}, "error": error})
    return failed


def _latent_fields(record: Mapping[str, Any]) -> dict[str, Any]:
    latent = record["latent"]
    return {
        "selected_lv": latent["selected_lv"],
        "component_selection": latent["component_selection"],
        "cv_mean_auroc": latent.get("cv_mean_auroc"),
    }


def check_identity(arm: Mapping[str, Any], record: Mapping[str, Any]) -> None:
    """The script's ``z_noise5`` must equal the analysis record exactly; a difference is a hard error."""

    statistics, p_values = d10._jsonable(arm["statistics"]), d10._jsonable(arm["p_values"])
    if statistics != record["statistics"] or p_values != record["p_values"]:
        raise RuntimeError(
            f"the script's {CASE_STUDY_ARM} differs from analyze_trajectories: statistics {statistics} vs "
            f"{record['statistics']}, p-values {p_values} vs {record['p_values']}."
        )


def analysis_params(settings: Mapping[str, Any], seed: int, **overrides: Any) -> TrajectoryAnalysisParams:
    parent = settings["parent"]["settings"]
    kwargs: dict[str, Any] = {
        "permutations": int(settings["permutations"]),
        "seed": seed,
        "n_jobs": 1,
        "orientation_estimator": parent["orientation_estimator"],
        "min_components": int(parent["min_components"]),
        **parent["cv"],
    }
    kwargs.update(overrides)
    return TrajectoryAnalysisParams(**kwargs)


def factorial_arms(
    Y: pd.DataFrame,
    Y0: pd.DataFrame,
    common: Mapping[str, Any],
    design: Any,
    *,
    permutations: int,
    seed: int,
    orientation_estimator: str,
    latent_fields: Mapping[str, Any],
) -> dict[str, dict[str, Any]]:
    """The six factorial arms and their corrected twins, on one basis and one RRPP seed."""

    arms: dict[str, dict[str, Any]] = {}
    for arm, (oracle, _) in FACTORIAL.items():
        columns = arm_covariates(arm)
        covariates = TrajectoryDataset(**common, covariate_cols=columns).covariate_design if columns else None
        latent = Y0 if oracle else Y
        kwargs = {"permutations": permutations, "seed": seed, "orientation_estimator": orientation_estimator}
        uncorrected = infer(latent, design, covariates, **kwargs)
        kappa = residual_df_kappa(latent.shape[0], uncorrected["p_red"])
        corrected = infer(latent, design, covariates, kappa=kappa, **kwargs)
        shared = {**latent_fields, "covariate_design_columns": columns, "error": None}
        arms[arm] = {**uncorrected, **shared}
        arms[f"{arm}{CORRECTED_SUFFIX}"] = {**corrected, **shared}
    return arms


def latent_arm(dataset: TrajectoryDataset, settings: Mapping[str, Any], seed: int, arm: str, p_red: int | None):
    params = analysis_params(settings, seed, covariate_latent_adjustment=LATENT_ADJUSTMENT, **LATENT_ARMS[arm])
    try:
        result = analysis.analyze_trajectories(dataset, params)
    except TrajectoryAnalysisError as exc:
        return _failed_arm(str(exc))
    record = result.record
    return {
        "statistics": record["statistics"],
        "p_values": record["p_values"],
        "null_q95": {statistic: entry.get("q95") for statistic, entry in record["null_summary"].items()},
        **_latent_fields(record),
        "p_red": p_red,
        "kappa": None,
        "covariate_design_columns": record["covariates"]["design_columns"],
        "error": None,
    }


def evaluate_unit(
    mode: str, effect: float, replicate: int, settings: Mapping[str, Any], *, reference: Any = None
) -> dict[str, Any]:
    """Every arm of one unit; the checks (reconstruction, rank one, identity) are hard errors."""

    start = perf_counter()
    parent = settings["parent"]["settings"]
    layers = list(parent["layers"])
    built = build_unit(mode, effect, replicate, settings, reference=reference)
    metadata = built["metadata"]
    seed = int(built["seeds"]["rrpp"])
    common: dict[str, Any] = {
        "blocks": [OmicsBlock(name, built["blocks"][name], "identity") for name in layers],
        "metadata": metadata,
        "sample_col": "sample_id",
        "group_col": "group",
        "stage_col": "stage",
        "stage_order": [str(stage) for stage in sorted(set(metadata["stage"].tolist()))],
    }
    six = TrajectoryDataset(**common, covariate_cols=list(PARENT_COVARIATES))
    checks: dict[str, Any] = {"injection_rank_one_residual": injection_residuals(built, layers)}
    arms: dict[str, dict[str, Any]] = {}
    p_red: int | None = None
    try:
        blind = analysis.analyze_trajectories(six, analysis_params(settings, seed))
    except TrajectoryAnalysisError as exc:
        arms.update({arm: _failed_arm(str(exc)) for arm in (*FACTORIAL_ARMS, *CORRECTED_ARMS)})
    else:
        Y = blind.latent_scores
        Y0, oracle_checks = oracle_latent(built, blind, layers)
        checks.update(oracle_checks)
        arms.update(
            factorial_arms(
                Y,
                Y0,
                common,
                trajectory_design(six),
                permutations=int(settings["permutations"]),
                seed=seed,
                orientation_estimator=parent["orientation_estimator"],
                latent_fields=_latent_fields(blind.record),
            )
        )
        check_identity(arms[CASE_STUDY_ARM], blind.record)
        p_red = arms[CASE_STUDY_ARM]["p_red"]
    for arm in LATENT_ARMS:
        arms[arm] = latent_arm(six, settings, seed, arm, p_red)
    return d10._jsonable(
        {
            "key": cs.unit_key(mode, effect, replicate),
            "mode": mode,
            "effect": float(effect),
            "replicate": replicate,
            "settings_hash": cs.settings_hash(settings),
            "parent_hash": settings["parent"]["hash"],
            "seeds": built["seeds"],
            "group_labels": built["group_labels"],
            "covariate_cells": d10.covariate_cell_summary(metadata, built["group_labels"]),
            "checks": checks,
            "arms": {arm: arms[arm] for arm in ARMS},
            "runtime_seconds": perf_counter() - start,
        }
    )


# ── Shards, resume ───────────────────────────────────────────────────────────


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
    out_dir: Path, settings: Mapping[str, Any], index: int, count: int, *, limit: int | None = None
) -> int:
    """Evaluate this shard's missing units (every N-th); returns how many were evaluated.

    A unit recorded in any shard file under ``records/`` counts as done; this
    shard's file is rewritten with its complete records first, so a torn line
    is re-run instead of duplicated.
    """

    path = cs.shard_path(out_dir, index, count)
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
        if position % count == index and cs.unit_key(*unit) not in done
    ]
    if limit is not None:
        pending = pending[:limit]
    from motco.simulations.reference import load_reference

    reference = load_reference()
    cs.check_headroom(parent_view(settings), reference)
    for position, unit in enumerate(pending, start=1):
        record = evaluate_unit(*unit, settings, reference=reference)
        with path.open("a", encoding="utf-8") as handle:
            handle.write(json.dumps(record, sort_keys=True, allow_nan=False) + "\n")
            handle.flush()
        print(
            f"[shard {index}/{count}] {position}/{len(pending)} {record['key']} {record['runtime_seconds']:.1f}s",
            file=sys.stderr,
            flush=True,
        )
    return len(pending)


# ── Tables (pure) ────────────────────────────────────────────────────────────


PVALUE_COLUMNS = ("run", "mode", "effect", "replicate", "arm", "statistic", "p_value")
RATE_COLUMNS = (
    "run",
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
RANK_COLUMNS = ("run", "arm", "rank_bucket", "statistic", "units", "replicates", "rejections", "rate")
CONTRAST_COLUMNS = (
    "contrast",
    "arm_b",
    "arm_a",
    "statistic",
    "role",
    "scope",
    "readings",
    "clusters",
    "b_rejections",
    "a_rejections",
    "b_only",
    "a_only",
    "net",
    "p_value",
)


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


def _runs(rows: Iterable[Mapping[str, Any]]) -> list[str]:
    return list(dict.fromkeys(str(row["run"]) for row in rows))


def rate_rows(rows: Iterable[Mapping[str, Any]], *, alpha: float, k: float) -> list[dict[str, Any]]:
    """Per run × control × arm × statistic: R available replicates, rejections (``p ≤ alpha``), rate, bound, verdict.

    Every arm carries ``holds`` | ``exceeds`` | ``unavailable`` against the
    parent check's one-sided bound α + k·√(α(1 − α)/R).
    """

    rows = list(rows)
    out = []
    for run in _runs(rows):
        run_rows = [row for row in rows if str(row["run"]) == run]
        tallies = cs.rejection_tallies(run_rows, alpha=alpha)
        effects = cs.control_effects(run_rows)
        for mode in CONTROL_MODES:
            effect = effects[mode]
            for arm in ARMS:
                for statistic in STATISTICS:
                    tally = tallies.get((mode, effect, arm, statistic)) if effect is not None else None
                    replicates, rejections, rate, se = cs._rate(tally)
                    bound = cs.anchor_bound(alpha, replicates, k) if replicates else None
                    if rate is None or bound is None:
                        verdict = "unavailable"
                    else:
                        verdict = "holds" if rate <= bound else "exceeds"
                    out.append(
                        {
                            "run": run,
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


def rank_bucket(selected_lv: Any) -> str | None:
    if selected_lv is None or selected_lv == "":
        return None
    rank = int(float(selected_lv))
    return ">=4" if rank >= 4 else str(rank)


def rank_rows(records: Mapping[str, Mapping[str, Any]], run: str, *, alpha: float) -> list[dict[str, Any]]:
    """Anchor units of ``z_noise5`` and the latent arms by selected rank (1, 2, 3, ≥ 4), with rejections per statistic.

    Every available anchor unit (the arm recorded no error) falls in exactly
    one bucket; ``replicates`` counts its finite p-values for the statistic.
    """

    anchor = [record for record in records.values() if str(record["mode"]) == ANCHOR]
    out = []
    for arm in RANK_TABLE_ARMS:
        available = [record["arms"][arm] for record in anchor if record["arms"][arm].get("error") is None]
        for bucket in RANK_BUCKETS:
            members = [entry for entry in available if rank_bucket(entry.get("selected_lv")) == bucket]
            for statistic in STATISTICS:
                p_values = [
                    float(entry["p_values"][statistic])
                    for entry in members
                    if _finite((entry.get("p_values") or {}).get(statistic))
                ]
                rejections = sum(1 for p_value in p_values if p_value <= alpha)
                out.append(
                    {
                        "run": run,
                        "arm": arm,
                        "rank_bucket": bucket,
                        "statistic": statistic,
                        "units": len(members),
                        "replicates": len(p_values),
                        "rejections": rejections,
                        "rate": rejections / len(p_values) if p_values else None,
                    }
                )
    return out


def _p_key(p_value: Any, scale: int) -> int | None:
    """A p-value as the integer k of k/scale (``None`` when unavailable)."""

    return int(round(float(p_value) * scale)) if _finite(p_value) else None


def read_parent_p_values(parent_dir: Path) -> dict[tuple[str, float, int, str, str], Any]:
    path = Path(parent_dir) / "case_study_type_i_pvalues.csv"
    if not path.exists():
        raise SystemExit(f"parent {parent_dir}: {path.name} does not exist; nothing to compare with.")
    with path.open(encoding="utf-8", newline="") as handle:
        return {
            (row["mode"], float(row["effect"]), int(row["replicate"]), row["arm"], row["statistic"]): (
                float(row["p_value"]) if row["p_value"] else None
            )
            for row in csv.DictReader(handle)
        }


def reproduction_report(
    records: Mapping[str, Mapping[str, Any]],
    settings: Mapping[str, Any],
    parent_p_values: Mapping[tuple[str, float, int, str, str], Any] | None,
) -> dict[str, Any]:
    """The judged arms' p-values against the parent's, as k/(B + 1) integers; ``not_compared`` for smoke runs."""

    parent = settings["parent"]
    header = {
        "parent": {key: parent[key] for key in ("label", "path", "hash")},
        "settings_hash": cs.settings_hash(settings),
        "permutations": int(settings["permutations"]),
        "parent_permutations": int(parent["settings"]["permutations"]),
        "units": len(records),
    }
    if settings["smoke"] or parent_p_values is None:
        return {
            **header,
            "status": "not_compared",
            "reason": "fewer permutations than the parent (smoke run)",
            "arms": {arm: {"parent_arm": parent_arm, "compared": 0, "mismatches": 0} for arm, parent_arm in
                     REPRODUCED.items()},
            "mismatches": [],
            "mismatching_units": [],
        }
    scale = int(parent["settings"]["permutations"]) + 1
    mismatches: list[dict[str, Any]] = []
    arms = {}
    for arm, parent_arm in REPRODUCED.items():
        compared = 0
        before = len(mismatches)
        for key in sorted(records):
            record = records[key]
            p_values = record["arms"][arm].get("p_values") or {}
            for statistic in STATISTICS:
                unit = (str(record["mode"]), float(record["effect"]), int(record["replicate"]), parent_arm, statistic)
                ours = _p_key(p_values.get(statistic), scale)
                theirs = _p_key(parent_p_values.get(unit), scale) if unit in parent_p_values else "absent"
                compared += 1
                if ours != theirs:
                    mismatches.append(
                        {"key": key, "arm": arm, "parent_arm": parent_arm, "statistic": statistic,
                         "diagnosis": ours, "parent": theirs}
                    )
        arms[arm] = {"parent_arm": parent_arm, "compared": compared, "mismatches": len(mismatches) - before}
    return {
        **header,
        "status": "failed" if mismatches else "holds",
        "scale": scale,
        "arms": arms,
        "mismatches": mismatches[:200],
        "mismatching_units": sorted({entry["key"] for entry in mismatches}),
    }


# ── Contrasts and decision (pure) ────────────────────────────────────────────


def sign_flip_p_value(differences: Iterable[int]) -> float:
    """Exact one-sided P(Σ ε_r d_r ≥ Σ d_r) over independent equiprobable signs ε_r (by convolution)."""

    values = [int(value) for value in differences]
    observed = sum(values)
    magnitudes = [abs(value) for value in values if value != 0]
    total = sum(magnitudes)
    distribution = np.zeros(2 * total + 1)
    distribution[total] = 1.0
    for magnitude in magnitudes:
        # The support never exceeds ±total, so the rolls never wrap mass around.
        distribution = 0.5 * (np.roll(distribution, magnitude) + np.roll(distribution, -magnitude))
    return float(min(1.0, distribution[total + observed :].sum()))


def _rejections(rows: Iterable[Mapping[str, Any]], statistic: str, *, alpha: float) -> dict[tuple, bool]:
    """(run, mode, effect, replicate, arm) → rejected, over the rows with a finite p-value for ``statistic``."""

    return {
        (str(row["run"]), str(row["mode"]), float(row["effect"]), int(row["replicate"]), str(row["arm"])): float(
            row["p_value"]
        )
        <= alpha
        for row in rows
        if str(row["statistic"]) == statistic and _finite(row.get("p_value"))
    }


def paired_contrast(
    rejected: Mapping[tuple, bool], arm_b: str, arm_a: str, *, readings: set[tuple] | None = None
) -> dict[str, Any]:
    """Discordant counts and the replicate-clustered sign-flip p-value that B rejects more often than A.

    A reading is (run, mode, effect) for one replicate; it enters only when
    both arms have a p-value. Each replicate contributes one ``d_r`` summed
    over its readings.
    """

    by_replicate: dict[int, int] = {}
    b_rejections = a_rejections = b_only = a_only = paired = 0
    for (run, mode, effect, replicate, arm), b in rejected.items():
        if arm != arm_b or (readings is not None and (run, mode, effect) not in readings):
            continue
        a = rejected.get((run, mode, effect, replicate, arm_a))
        if a is None:
            continue
        paired += 1
        b_rejections += int(b)
        a_rejections += int(a)
        b_only += int(b and not a)
        a_only += int(a and not b)
        by_replicate[replicate] = by_replicate.get(replicate, 0) + int(b and not a) - int(a and not b)
    return {
        "readings": paired,
        "clusters": len(by_replicate),
        "b_rejections": b_rejections,
        "a_rejections": a_rejections,
        "b_only": b_only,
        "a_only": a_only,
        "net": b_only - a_only,
        "p_value": sign_flip_p_value(by_replicate.values()),
    }


def contrast_rows(rows: Iterable[Mapping[str, Any]], *, alpha: float) -> list[dict[str, Any]]:
    """Every predeclared contrast on every statistic: clustered over both runs and controls, then per reading.

    Only the clustered ``delta`` rows are judged; the rest are descriptive.
    """

    rows = list(rows)
    readings = sorted({(str(r["run"]), str(r["mode"]), float(r["effect"])) for r in rows})
    out = []
    for statistic in STATISTICS:
        rejected = _rejections(rows, statistic, alpha=alpha)
        for label, arm_b, arm_a in CONTRASTS:
            head = {"contrast": label, "arm_b": arm_b, "arm_a": arm_a, "statistic": statistic}
            role = "judged" if statistic == JUDGED_STATISTIC else "descriptive"
            out.append({**head, "role": role, "scope": "clustered", **paired_contrast(rejected, arm_b, arm_a)})
            for reading in readings:
                scope = f"{reading[0]}/{reading[1]}"
                result = paired_contrast(rejected, arm_b, arm_a, readings={reading})
                out.append({**head, "role": "descriptive", "scope": scope, **result})
    return out


def _arm_verdicts(rates: Sequence[Mapping[str, Any]], arm: str) -> list[dict[str, Any]]:
    return [
        {key: row[key] for key in ("run", "mode", "statistic", "rate", "bound", "verdict")}
        for row in rates
        if row["arm"] == arm
    ]


def _nested_rates(rates: Sequence[Mapping[str, Any]]) -> dict[str, Any]:
    nested: dict[str, Any] = {}
    for row in rates:
        nested.setdefault(row["run"], {}).setdefault(row["mode"], {}).setdefault(row["arm"], {})[row["statistic"]] = {
            key: row[key] for key in ("replicates", "rejections", "rate", "bound", "verdict")
        }
    return nested


def floor1_summary(rank_table: Iterable[Mapping[str, Any]]) -> dict[str, Any]:
    """Per run: the floor-1 arm's anchor units by rank bucket and its rank-1 fraction (descriptive)."""

    summary: dict[str, Any] = {}
    for row in rank_table:
        if row["arm"] != "latent_stage_floor1" or row["statistic"] != JUDGED_STATISTIC:
            continue
        entry = summary.setdefault(str(row["run"]), {"units_by_rank": {}})
        entry["units_by_rank"][str(row["rank_bucket"])] = int(row["units"])
    for entry in summary.values():
        total = sum(entry["units_by_rank"].values())
        entry["units"] = total
        entry["rank1_fraction"] = entry["units_by_rank"].get("1", 0) / total if total else None
    return summary


def diagnosis_decision(
    rows: Iterable[Mapping[str, Any]],
    rank_table: Iterable[Mapping[str, Any]],
    reproductions: Sequence[Mapping[str, Any]],
    *,
    alpha: float,
    k: float,
) -> dict[str, Any]:
    """The verdicts from the merged p-value rows of every run, the rank rows and each run's reproduction.

    Any reproduction mismatch gives ``reproduction_failed`` and no verdict.
    Otherwise the covariate-adjusted ``delta`` verdict is
    ``residual_df_shrinkage`` (C1 or C2 significant, ``z_noise5_corrected``
    holding every statistic at both controls in every run, C4 not significant),
    else ``collinear_covariate`` (C3 significant), else ``unresolved``; the
    latent verdict is ``rank_selection`` when ``latent_stage_rank2`` holds
    everywhere, else ``residualized_basis``. A run that was not compared
    (fewer permutations than its parent) makes every verdict indicative.
    """

    rows = list(rows)
    rank_table = list(rank_table)
    rates = rate_rows(rows, alpha=alpha, k=k)
    contrasts = contrast_rows(rows, alpha=alpha)
    clustered = [row for row in contrasts if row["scope"] == "clustered"]
    judged = {row["contrast"]: row for row in clustered if row["statistic"] == JUDGED_STATISTIC}
    mismatching = sorted({key for report in reproductions for key in report.get("mismatching_units", [])})
    statuses = {str(report["parent"]["label"]): report["status"] for report in reproductions}
    failed = any(status == "failed" for status in statuses.values())
    indicative = any(status != "holds" for status in statuses.values())
    reproduction = {
        "status": "failed" if failed else ("not_compared" if indicative else "holds"),
        "runs": {
            str(report["parent"]["label"]): {
                "status": report["status"],
                "arms": report.get("arms"),
                "mismatching_units": report.get("mismatching_units", []),
            }
            for report in reproductions
        },
        "mismatching_units": mismatching,
    }
    decision: dict[str, Any] = {
        "rule": {
            "rate": "share of available replicates with p <= alpha, against alpha + k * sqrt(alpha(1 - alpha)/R)",
            "contrasts": (
                "exact one-sided sign flip on per-replicate discordant rejections (d_r over both runs and both "
                f"controls), significant at p <= {CONTRAST_ALPHA}; C1-C4 on delta judged"
            ),
            "adjusted_delta": (
                "residual_df_shrinkage when C1 or C2 is significant, z_noise5_corrected holds every statistic at "
                "both controls in every run and C4 is not significant; else collinear_covariate when C3 is "
                "significant; else unresolved"
            ),
            "latent_stage": (
                "rank_selection when latent_stage_rank2 holds every statistic at both controls in every run; else "
                "residualized_basis; latent_stage_floor1 descriptive"
            ),
            "gate": "any reproduction mismatch gives reproduction_failed and no verdict",
        },
        "alpha": alpha,
        "k": k,
        "contrast_alpha": CONTRAST_ALPHA,
        "runs": _runs(rows),
        "reproduction": reproduction,
        "indicative": indicative,
        "contrasts": clustered,
        "rates": _nested_rates(rates),
        "floor1": floor1_summary(rank_table),
    }
    if failed:
        return {**decision, "decision": "reproduction_failed", "verdicts": None}

    significant = {label: row["p_value"] <= CONTRAST_ALPHA for label, row in judged.items()}
    corrected = _arm_verdicts(rates, f"{CASE_STUDY_ARM}{CORRECTED_SUFFIX}")
    corrected_holds = bool(corrected) and all(row["verdict"] == "holds" for row in corrected)
    if (significant["C1"] or significant["C2"]) and corrected_holds and not significant["C4"]:
        adjusted = "residual_df_shrinkage"
    elif significant["C3"]:
        adjusted = "collinear_covariate"
    else:
        adjusted = "unresolved"
    rank2 = _arm_verdicts(rates, "latent_stage_rank2")
    rank2_failing = [row for row in rank2 if row["verdict"] != "holds"]
    latent = "rank_selection" if rank2 and not rank2_failing else "residualized_basis"
    decision["verdicts"] = {
        "adjusted_delta": {
            "verdict": adjusted,
            "consequence": CONSEQUENCES[adjusted],
            "next_step": NEXT_STEPS[adjusted],
            "indicative": indicative,
            "conditions": {
                "column_contrast_significant": {"C1": significant["C1"], "C2": significant["C2"]},
                "corrected_case_study_model_holds": corrected_holds,
                "corrected_not_holding": [row for row in corrected if row["verdict"] != "holds"],
                "corrected_column_contrast_significant": significant["C4"],
                "collinear_contrast_significant": significant["C3"],
                "p_values": {label: row["p_value"] for label, row in judged.items()},
            },
        },
        "latent_stage": {
            "verdict": latent,
            "consequence": CONSEQUENCES[latent],
            "indicative": indicative,
            "rank2_not_holding": rank2_failing,
            "cv_exceeding": [row for row in _arm_verdicts(rates, "latent_stage_cv") if row["verdict"] == "exceeds"],
            "floor1": "descriptive; does not enter the verdict",
        },
    }
    return {**decision, "decision": "diagnosed"}


# ── Merge and decide ─────────────────────────────────────────────────────────


def _write_csv(path: Path, rows: Sequence[Mapping[str, Any]], columns: Sequence[str]) -> None:
    d10._write_csv(path, rows, columns)


def _read_csv(path: Path) -> list[dict[str, Any]]:
    with path.open(encoding="utf-8", newline="") as handle:
        return [{key: (value if value != "" else None) for key, value in row.items()} for row in csv.DictReader(handle)]


def merge(out_dir: Path, settings: Mapping[str, Any], parent_dir: Path | None = None) -> dict[str, Any]:
    """Merge every shard, then write the p-value, rate and rank tables and the reproduction report."""

    out_dir = Path(out_dir)
    records = read_records(sorted((out_dir / "records").glob("shard_*.jsonl")))
    expected = [cs.unit_key(*unit) for unit in units(settings)]
    with (out_dir / "merged.jsonl").open("w", encoding="utf-8") as handle:
        for key in sorted(records):
            handle.write(json.dumps(records[key], sort_keys=True, allow_nan=False) + "\n")
    run = str(settings["parent"]["label"])
    alpha, k = float(settings["alpha"]), float(settings["k"])
    rows = p_value_rows(records, run)
    _write_csv(out_dir / "diagnosis_pvalues.csv", rows, PVALUE_COLUMNS)
    rates = rate_rows(rows, alpha=alpha, k=k)
    _write_csv(out_dir / "diagnosis_rates.csv", rates, RATE_COLUMNS)
    _write_csv(out_dir / "diagnosis_rank_table.csv", rank_rows(records, run, alpha=alpha), RANK_COLUMNS)
    parent_p_values = None
    if not settings["smoke"]:
        parent_dir = Path(parent_dir) if parent_dir is not None else _resolve_path(settings["parent"]["path"])
        bound = load_parent(parent_dir)
        if bound["hash"] != settings["parent"]["hash"]:
            raise SystemExit(f"parent {parent_dir} has hash {bound['hash']}, not the recorded parent's.")
        parent_p_values = read_parent_p_values(parent_dir)
    reproduction = reproduction_report(records, settings, parent_p_values)
    (out_dir / "diagnosis_reproduction.json").write_text(
        json.dumps(d10._jsonable(reproduction), indent=2, allow_nan=False) + "\n", encoding="utf-8"
    )
    checks = [record["checks"] for record in records.values()]
    return {
        "run": run,
        "units": len(records),
        "expected": len(expected),
        "missing": [key for key in expected if key not in records][:20],
        "smoke": settings["smoke"],
        "reproduction": {key: reproduction[key] for key in ("status", "arms")},
        "failures": {
            arm: sum(1 for record in records.values() if record["arms"][arm].get("error") is not None) for arm in ARMS
        },
        "max_reconstruction_error": max(
            (c["reconstruction_max_abs_error"] for c in checks if "reconstruction_max_abs_error" in c), default=None
        ),
        "max_rank_one_residual": max(
            (max(c["injection_rank_one_residual"].values()) for c in checks), default=None
        ),
        "exceeding": [
            {key: row[key] for key in ("mode", "arm", "statistic", "rate", "bound")}
            for row in rates
            if row["verdict"] == "exceeds"
        ],
    }


def read_run(run_dir: Path) -> dict[str, Any]:
    """A merged run: its settings, p-value rows, rank rows and reproduction report."""

    run_dir = Path(run_dir)
    settings = d10.load_settings(run_dir)
    missing = [
        name
        for name in ("diagnosis_pvalues.csv", "diagnosis_rank_table.csv", "diagnosis_reproduction.json")
        if not (run_dir / name).exists()
    ]
    if missing:
        raise SystemExit(f"{run_dir}: {', '.join(missing)} missing; run merge first.")
    return {
        "settings": settings,
        "rows": _read_csv(run_dir / "diagnosis_pvalues.csv"),
        "rank_rows": _read_csv(run_dir / "diagnosis_rank_table.csv"),
        "reproduction": json.loads((run_dir / "diagnosis_reproduction.json").read_text(encoding="utf-8")),
    }


def decide(run_dirs: Sequence[Path], out_dir: Path) -> dict[str, Any]:
    """Read the merged runs, write ``diagnosis_contrasts.csv`` and ``case_study_diagnosis_decision.json``."""

    runs = [read_run(run_dir) for run_dir in run_dirs]
    labels = [str(run["settings"]["parent"]["label"]) for run in runs]
    if len(set(labels)) != len(labels):
        raise SystemExit(f"decide needs runs against distinct parents; got {labels}.")
    shared = {
        (
            run["settings"]["alpha"],
            run["settings"]["k"],
            run["settings"]["parent"]["settings"]["base_seed"],
        )
        for run in runs
    }
    if len(shared) != 1:
        raise SystemExit("decide needs runs sharing alpha, k and the base seed (the replicate clusters).")
    alpha, k, _ = shared.pop()
    rows = [row for run in runs for row in run["rows"]]
    rank_table = [row for run in runs for row in run["rank_rows"]]
    decision = diagnosis_decision(rows, rank_table, [run["reproduction"] for run in runs], alpha=alpha, k=k)
    decision["run_settings"] = [
        {
            "label": label,
            "settings_hash": cs.settings_hash(run["settings"]),
            "parent_hash": run["settings"]["parent"]["hash"],
            "parent_path": run["settings"]["parent"]["path"],
            "strengths": cs.strengths(run["settings"]["parent"]["settings"]),
            "replicates": run["settings"]["replicates"],
            "permutations": run["settings"]["permutations"],
            "smoke": run["settings"]["smoke"],
        }
        for label, run in zip(labels, runs)
    ]
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    _write_csv(out_dir / "diagnosis_contrasts.csv", contrast_rows(rows, alpha=alpha), CONTRAST_COLUMNS)
    (out_dir / "case_study_diagnosis_decision.json").write_text(
        json.dumps(d10._jsonable(decision), indent=2, allow_nan=False) + "\n", encoding="utf-8"
    )
    return decision


# ── CLI ──────────────────────────────────────────────────────────────────────


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Diagnosis of the case-study configuration's size excess (no R).")
    sub = parser.add_subparsers(dest="command", required=True)
    run = sub.add_parser("run", help="Evaluate (a shard of) the units; an unsharded run also merges")
    run.add_argument("--parent", type=Path, required=True, help="Committed run directory of the parent check")
    run.add_argument("--out-dir", type=Path, required=True, help="Output directory (one per parent)")
    run.add_argument("--reps", type=int, default=None, help="Replicates, a prefix of the parent's (default: all)")
    run.add_argument(
        "--perms", type=int, default=None, help="RRPP permutations (default: the parent's; fewer marks a smoke run)"
    )
    run.add_argument("--shard", type=cs.parse_shard, default=(0, 1), help="Run shard i of N (i/N; default 0/1)")
    run.add_argument("--limit", type=int, default=None, help="Evaluate at most this many pending units")
    merge_parser = sub.add_parser("merge", help="Merge every shard and write the per-run tables")
    merge_parser.add_argument("--out-dir", type=Path, required=True, help="Output directory")
    merge_parser.add_argument(
        "--parent", type=Path, default=None, help="Parent run directory (default: the recorded path)"
    )
    decide_parser = sub.add_parser("decide", help="Contrasts and verdicts over the merged runs")
    decide_parser.add_argument("--run", type=Path, action="append", required=True, help="A merged run directory")
    decide_parser.add_argument("--out-dir", type=Path, required=True, help="Where to write the decision")
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    if args.command == "decide":
        decision = decide(args.run, args.out_dir)
        print(json.dumps({key: decision[key] for key in ("decision", "indicative", "verdicts")}, indent=1))
        return 0
    if args.command == "merge":
        settings = d10.load_settings(args.out_dir)
        print(json.dumps(merge(args.out_dir, settings, args.parent), indent=1))
        return 0
    parent = load_parent(args.parent)
    settings = build_settings(args.parent, parent, reps=args.reps, perms=args.perms)
    ensure_settings(args.out_dir, settings)
    index, count = args.shard
    evaluated = run_shard(args.out_dir, settings, index, count, limit=args.limit)
    print(f"evaluated {evaluated} unit(s) into {cs.shard_path(args.out_dir, index, count)}", file=sys.stderr)
    if count > 1 or args.limit is not None:
        return 0
    print(json.dumps(merge(args.out_dir, settings, args.parent), indent=1))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
