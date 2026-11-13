#!/usr/bin/env python3
"""Why the predictive-subspace coordinates exceed size on the SEA-AD plasmode anchor: a diagnosis (no R).

The driver behind ``openspec/changes/predictive-size-excess-diagnosis/`` (spec
``predictive-size-excess-diagnosis``, roadmap open item 2). The
predictive-subspace pilot (2026-11-12) read ``keep_production`` because on the
plasmode anchor its three predictive coordinate systems rejected ``delta`` in
65–68 of 500 replicates and ``shape`` in 38–39, where production rejected 21
and 13 on the same units. The guess: the two predictive directions of the
rank-5 fit hold the in-sample stage fit, and their residuals are not
exchangeable across stages, which RRPP's unrestricted draws assume. This
driver tests the guess twice, under a decision fixed before any unit runs:

- **Check W** keeps the coordinates and changes the null: RRPP's draws
  restricted within stage (the reduced model's residuals permuted among the
  donors of each stage, never across stages).
- **Check O** keeps the null and changes the coordinates: the fitted stage
  indicators computed out of fold (each donor's row from a stage-PLS model
  fitted without it).

No package code changes: both repairs live here, guarded by an identity with
the package on every unit.

**Units** (``predictive_subspace_check.plasmode_settings`` / ``plasmode_units``
and ``plasmode_check.unit_seeds``, imported unchanged): the pilot's 800
plasmode units, the anchor at the ``paper`` stream's replicates 0–499 and the
real-noise span plant at e ∈ {0.25, 0.5, 1.0} at replicates 0–99, 999
permutations, the case-study configuration (six covariates, ``df`` residual
correction, net displacement, rank floor 2). Each unit carries the pilot's
label draw, plant and RRPP seed; rejection is p ≤ α. The settings must agree
with the pilot's committed plasmode settings, and the export's files with the
plasmode's.

**One fit per unit.** An anchor unit fits at ``forced_components`` = the
plasmode's committed fidelity rank (5) without cross-validation: the anchor
does not move the label-free features, so every anchor unit's CV selects that
rank, and ``preflight`` checks that the forced and cross-validated calls agree
exactly. A plant unit runs the production CV once; every other arm is forced at
its rank. The pilot's ``same_fit`` check runs on every unit.

**Arms** (fourteen per unit, every one at the unit's RRPP seed, serial draws):

- reference — ``pls_cv`` (scores) and ``predictive``, ``predictive_orthonormal``,
  ``predictive_fitted`` (the in-sample predictive coordinates), each one
  ``analyze_trajectories`` call: these must reproduce the pilot's counts;
- decisive — ``predictive__within_stage``, ``predictive_orthonormal__within_stage``
  and ``predictive_fitted__within_stage`` (check W) and ``fitted_out_of_fold``
  (check O);
- descriptive — ``pls_cv__within_stage``, ``fitted_out_of_fold__within_stage``
  and the four reference coordinate systems analysed with no covariates in the
  trajectory model (``<arm>__no_covariates``, ``df`` kept).

The within-stage and out-of-fold arms run through the driver's inference path:
the covariate-adjusted models from ``motco.analysis._trajectory_inputs``,
``estimate_difference`` under net displacement, the null from serial ``RRPP``
with ``residual_correction="df"`` or from the within-stage loop (which mirrors
``RRPP``'s serial path draw for draw; with one stratum it is ``RRPP``), and
``_empirical_p_value``. One restricted index sequence per unit is shared by its
within-stage arms. The out-of-fold coordinates use four stage-stratified,
shuffled folds from a seed derived by sha256 from the driver's name, base seed
2037, the replicate and ``folds``; each held-out row is
``transform(X_held) · y_loadingsᵀ`` of ``fit_plsda_model`` fitted on the other
folds' rows of the analysis's PLS input features at the unit's rank and
scaling. The fold seed is recorded, the assignment is not.

**Checks on every unit.** The pilot's five consistency checks (``rank``,
``fitted_labels``, ``same_fit``, ``one_subspace``, ``identity_at_low_rank``);
``same_fit_covariate_free`` (each covariate-free arm fitted production's model);
and ``path_identity`` (the driver's path under unrestricted draws returns each
reference arm's statistics to 1e-10 and its p-values exactly). A failing unit is
``inconsistent``, naming the check, the arm and the deviation, and leaves every
rate. A unit that raises is ``failed``; an arm that raises is failed on its unit.
Each unit also records, per coordinate system, the reduced model's residual
dispersion by stage relative to the pooled one, and the in-sample stage R²
beside the out-of-fold explained share. Records hold seeds and statistics only:
no donor identifier, label vector or fold assignment.

**Parents** (sha256 in the settings; ``preflight``, ``run``, ``merge`` and
``decide`` refuse a changed file by name): the pilot's plasmode ``settings.json``
and ``rate_table.csv``, and through ``plasmode_settings`` the plasmode's
``settings.json`` and ``plasmode_summary.csv`` and the production selection's
plasmode ``settings.json`` and ``rate_table.csv``.

Subcommands:

- ``preflight`` — writes ``settings.json`` and ``preflight.json``: the
  cross-validated ``pls_cv`` on anchor replicates 0 and 1 must select the
  fidelity rank and AUROC and give the forced call's statistics and p-values
  exactly. ``run`` refuses to start without a passed pre-flight.
- ``run`` — evaluates (a shard of) the units into
  ``records/shard_<i>_of_<N>.jsonl``. Shards are round-robin; resume skips
  complete units; a resume under other settings exits naming the keys.
- ``merge`` — ``merged.jsonl`` (keep off the repository); ``rate_table.csv``
  (per cell, arm and statistic: available and rejected counts, rate, Monte
  Carlo SE, failures); ``paired_table.csv`` (each non-reference arm against its
  counterpart, and every arm against ``pls_cv``: both rates, the discordant
  counts b (only the reference rejects) and c (only the arm), the exact
  one-sided ``p_fewer`` = P(Bin(b + c, ½) ≥ b), ``p_more`` = P(Bin(b + c, ½) ≥ c)
  and two-sided p-values, 1 when b + c = 0); ``cell_table.csv`` (per cell and
  arm: rank, the dispersion medians by stage and their ratio, both stage-fit
  readings, the median ``angle`` null 95th percentile); ``accounting.json``; and
  ``reproduction.json`` (the four reference arms' counts in every cell against
  the pilot's committed ``rate_table.csv``, exact).
- ``decide`` — ``predictive_size_excess_decision.json`` (rule version 1). The
  gate ``reproduction_failed`` (a failed reproduction, an inconsistent unit, or
  a missing or failed pre-flight; the grades are still recorded). Per decisive
  arm on the anchor: ``restored`` when its ``delta``, ``angle`` and ``shape``
  rejections hold the stage-span paper grade's count-plus-ceiling rule
  (``judge_rates``, imported unchanged; m = 3); else ``reduced`` when the
  one-sided exact paired p-value for fewer rejections than its counterpart is at
  most 0.025 on ``delta`` or ``shape``; else ``unchanged`` (an unavailable
  reading is ``unchanged`` and named). Check W takes the weakest grade of its
  three arms, check O its arm's. Verdict ``both_restore`` | ``null_restores`` |
  ``out_of_fold_restores`` | ``partly_explained`` | ``unexplained``, each with
  its fixed consequence. Descriptive tables beside it never change it: the
  plant cells, the descriptive arms, the covariate-free arms, the dispersion
  and the null widths.

With ``--workers`` above 1, pin BLAS to one thread per process
(``OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1``).

Examples:

    export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
    EXPORT=../sead_sexdiff/results/motco_input/cohort_30
    # smoke: two units per cell
    python scripts/predictive_size_excess_check.py preflight --export-dir $EXPORT --reps 2 --plant-reps 2 \\
        --out-dir /tmp/pse
    python scripts/predictive_size_excess_check.py run --export-dir $EXPORT --reps 2 --plant-reps 2 \\
        --out-dir /tmp/pse --workers 7
    python scripts/predictive_size_excess_check.py merge --out-dir /tmp/pse
    python scripts/predictive_size_excess_check.py decide --out-dir /tmp/pse

    # one SLURM array task per shard (scripts/predictive_size_excess_array.sbatch), then merge and decide
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
from typing import Any, Iterable, Mapping, Sequence

import numpy as np
import pandas as pd

from motco import analysis
from motco.simulations.evaluation import (
    _empirical_p_value,
    _extract_null_distributions,
    _extract_pair_statistics,
    _summarize_null_distributions,
)
from motco.stats.permutation import RRPP, residual_correction_scale
from motco.stats.pls import fit_plsda_model
from motco.stats.trajectory import estimate_betas, estimate_difference

sys.path.insert(0, str(Path(__file__).resolve().parent))
import latent_space_comparison as lsc  # noqa: E402  (units, settings hashing, records, tables)
import plasmode_check as pc  # noqa: E402  (export contract, seeds, the case-study analysis)
import predictive_subspace_check as psc  # noqa: E402  (the pilot's plasmode settings, units, arms and checks)
import production_latent_space_check as plsc  # noqa: E402  (span plant, parents, exact paired tests)
import stage_span_paper_grade_decision as grade  # noqa: E402  (the count-plus-ceiling rule)

sys.path.pop(0)

DRIVER_NAME = "predictive-size-excess"
DRIVER_VERSION = 1
RULE_VERSION = 1

REFERENCE_ARM = psc.REFERENCE_ARM
PREDICTIVE_ARMS: tuple[str, ...] = psc.CANDIDATE_ARMS
REFERENCE_ARMS: tuple[str, ...] = (REFERENCE_ARM, *PREDICTIVE_ARMS)
WITHIN_STAGE = "__within_stage"
NO_COVARIATES = "__no_covariates"
OUT_OF_FOLD_ARM = "fitted_out_of_fold"
WITHIN_STAGE_ARMS: tuple[str, ...] = tuple(f"{arm}{WITHIN_STAGE}" for arm in PREDICTIVE_ARMS)
DECISIVE_ARMS: tuple[str, ...] = (*WITHIN_STAGE_ARMS, OUT_OF_FOLD_ARM)
COVARIATE_FREE_ARMS: tuple[str, ...] = tuple(f"{arm}{NO_COVARIATES}" for arm in REFERENCE_ARMS)
DESCRIPTIVE_ARMS: tuple[str, ...] = (
    f"{REFERENCE_ARM}{WITHIN_STAGE}",
    f"{OUT_OF_FOLD_ARM}{WITHIN_STAGE}",
    *COVARIATE_FREE_ARMS,
)
ARMS: tuple[str, ...] = (*REFERENCE_ARMS, *DECISIVE_ARMS, *DESCRIPTIVE_ARMS)

OUT_OF_FOLD = "out_of_fold"
SCORES = "scores"
COORDINATE_SYSTEMS: tuple[str, ...] = (SCORES, *PREDICTIVE_ARMS, OUT_OF_FOLD)


def _arm_spec(arm: str) -> dict[str, str]:
    """An arm's coordinate system, draws, covariates and role (design D7)."""

    base = arm.removesuffix(WITHIN_STAGE).removesuffix(NO_COVARIATES)
    coordinates = SCORES if base == REFERENCE_ARM else (OUT_OF_FOLD if base == OUT_OF_FOLD_ARM else base)
    if arm in REFERENCE_ARMS:
        role = "reference"
    elif arm in DECISIVE_ARMS:
        role = "check_W" if arm in WITHIN_STAGE_ARMS else "check_O"
    else:
        role = "descriptive"
    return {
        "coordinates": coordinates,
        "draws": "within_stage" if arm.endswith(WITHIN_STAGE) else "unrestricted",
        "covariates": "none" if arm.endswith(NO_COVARIATES) else "six",
        "role": role,
    }


ARM_SPECS: dict[str, dict[str, str]] = {arm: _arm_spec(arm) for arm in ARMS}

#: Each non-reference arm's counterpart: the same coordinate system under unrestricted draws and six
#: covariates; the out-of-fold arms pair with ``predictive_fitted`` under the same draws.
COUNTERPARTS: dict[str, str] = {
    **{f"{arm}{WITHIN_STAGE}": arm for arm in REFERENCE_ARMS},
    OUT_OF_FOLD_ARM: "predictive_fitted",
    f"{OUT_OF_FOLD_ARM}{WITHIN_STAGE}": f"predictive_fitted{WITHIN_STAGE}",
    **{f"{arm}{NO_COVARIATES}": arm for arm in REFERENCE_ARMS},
}

STATISTICS: tuple[str, ...] = ("delta", "angle", "shape")
#: The statistics the pilot's anchor exceeded on: the ``reduced`` grade reads these.
REDUCED_STATISTICS: tuple[str, ...] = ("delta", "shape")
ALPHA = 0.05
REJECTION = "p <= alpha"
#: One-sided level of the ``reduced`` grade's exact paired test (0.05 / 2 statistics).
REDUCED_ALPHA = 0.025
GRADES: tuple[str, ...] = ("unchanged", "reduced", "restored")  # weakest first
VERDICTS: tuple[str, ...] = (
    "reproduction_failed", "both_restore", "null_restores", "out_of_fold_restores", "partly_explained", "unexplained",
)
CHECKS: dict[str, tuple[str, ...]] = {"W": WITHIN_STAGE_ARMS, "O": (OUT_OF_FOLD_ARM,)}

PATH_IDENTITY_TOLERANCE = 1e-10
CONSISTENCY_TOLERANCE = psc.CONSISTENCY_TOLERANCE
CONSISTENCY_CHECKS: tuple[str, ...] = (*psc.CONSISTENCY_CHECKS, "same_fit_covariate_free", "path_identity")

RESIDUAL_CORRECTION = "df"
ORIENTATION_ESTIMATOR = "net_displacement"
FOLDS = 4
FOLD_BASE_SEED = 2037
FOLD_PURPOSE = "folds"
PREFLIGHT_REPLICATES: tuple[int, ...] = (0, 1)

REPO = Path(__file__).resolve().parents[1]
PILOT_PLASMODE_DIR = REPO / "results" / "predictive-subspace-2026-11-12" / "plasmode"
PLASMODE_DIR = psc.PLASMODE_DIR
PRODUCTION_PLASMODE_DIR = psc.PRODUCTION_PLASMODE_DIR
DESIGN = psc.PLASMODE_DESIGN
ANCHOR_CELL = psc.ANCHOR_CELL
PLANT_CELLS: tuple[str, ...] = psc.PLANT_CELLS
PRODUCTION_PLANT_CELL = "span_e0.5"
ANCHOR_REPLICATES = psc.PLASMODE_ANCHOR_REPLICATES
PLANT_REPLICATES = psc.PLASMODE_PLANT_REPLICATES
PERMUTATIONS = psc.PLASMODE_PERMUTATIONS
#: The pilot's plasmode settings this run must share (replicates may only be fewer).
_PILOT_SHARED_KEYS: tuple[str, ...] = (
    "driver_version", "profile", "design", "cells", "plasmode", "plant", "fidelity", "export", "attribution",
    "predictive_arm_fit", "arms", "alpha", "rejection", "consistency_tolerance", "permutations", "parents",
)

PREFLIGHT_NAME = "preflight.json"
DECISION_NAME = "predictive_size_excess_decision.json"

CONSEQUENCES: dict[str, str] = {
    "reproduction_failed": (
        "Nothing is read: a reference arm does not reproduce the pilot's committed counts, a unit's inference path "
        "or fit departs from the analysis, or the pre-flight did not pass. Find the drift before any rerun; a "
        "mismatch is diagnosed, never waived."
    ),
    "both_restore": (
        "The unrestricted draws are what fails in the predictive subspace on this export, and the out-of-fold "
        "reading confirms that the in-sample fit is where the excess lives. A follow-up change may propose "
        "within-stage draws as an opt-in option and must read their size and power on the pilot's simulated "
        "designs, the weak design and this anchor before predictive_orthonormal is re-read."
    ),
    "null_restores": (
        "The unrestricted draws are what fails in the predictive subspace on this export. A follow-up change may "
        "propose within-stage draws as an opt-in option and must read their size and power on the pilot's "
        "simulated designs, the weak design and this anchor before predictive_orthonormal is re-read."
    ),
    "out_of_fold_restores": (
        "The excess lives in the in-sample fit and a restricted null does not repair it. A repair needs a "
        "cross-fitted measurement, which belongs to open item 2's held-out half and open item 5's cross-fitted "
        "arm. The predictive values stay diagnostics at SEA-AD's design."
    ),
    "partly_explained": (
        "The guess is incomplete. The report names which check reduced and by how much; the remainder waits for "
        "open item 2's cohorts."
    ),
    "unexplained": (
        "The guess is rejected and the roadmap's prediction for ADNI is withdrawn. The covariate-free arms and the "
        "dispersion table are the leads."
    ),
}


# ── Settings, parents and units ──────────────────────────────────────────────


def _canonical(value: Any) -> Any:
    return json.loads(json.dumps(value, sort_keys=True))


def _seed(*parts: Any) -> int:
    encoded = json.dumps([DRIVER_NAME, *parts], separators=(",", ":")).encode("utf-8")
    return int(hashlib.sha256(encoded).hexdigest()[:8], 16) & 0x7FFFFFFF


def fold_seed(replicate: int) -> int:
    """The replicate's fold seed: sha256 of the driver's name, base seed 2037, the replicate and ``folds``."""

    return _seed(FOLD_BASE_SEED, int(replicate), FOLD_PURPOSE)


def check_pilot(pilot: Mapping[str, Any], committed: Mapping[str, Any], *, source: Path) -> None:
    """Refuse settings whose units would not be the pilot's committed plasmode units, naming the keys."""

    differing = [key for key in _PILOT_SHARED_KEYS if _canonical(pilot.get(key)) != _canonical(committed.get(key))]
    for kind in ("anchor", "plant"):
        ours, theirs = int(pilot["replicates"][kind]), int((committed.get("replicates") or {}).get(kind, -1))
        if ours > theirs:
            differing.append(f"replicates.{kind}")
    if differing:
        raise SystemExit(
            f"The plasmode units differ from the pilot's ({source}) in: {', '.join(differing)}; they would not be "
            "the units the excess was read on."
        )


def build_settings(
    export: Any,
    *,
    pilot_dir: Path = PILOT_PLASMODE_DIR,
    parent_dir: Path = PLASMODE_DIR,
    production_dir: Path = PRODUCTION_PLASMODE_DIR,
    replicates: int = ANCHOR_REPLICATES,
    plant_replicates: int = PLANT_REPLICATES,
    permutations: int = PERMUTATIONS,
) -> dict[str, Any]:
    """The pilot's plasmode profile (refused on another export or other units) plus this driver's arms."""

    pilot_dir = Path(pilot_dir)
    pilot = psc.plasmode_settings(
        export,
        parent_dir=parent_dir,
        production_dir=production_dir,
        replicates=replicates,
        plant_replicates=plant_replicates,
        permutations=permutations,
    )
    committed = plsc.load_parent(pilot_dir)["settings"]
    check_pilot(pilot, committed, source=pilot_dir / "settings.json")
    rates = pilot_dir / "rate_table.csv"
    return {
        "driver": DRIVER_NAME,
        "driver_version": DRIVER_VERSION,
        "pilot": pilot,
        "parents": {
            **pilot["parents"],
            "pilot_plasmode_settings": psc._parent(pilot_dir / "settings.json"),
            "pilot_plasmode_rates": psc._parent(rates),
        },
        "design": DESIGN,
        "anchor_rank": int(pilot["fidelity"]["selected_lv"]),
        "arms": list(ARMS),
        "arm_specs": {arm: dict(spec) for arm, spec in ARM_SPECS.items()},
        "counterparts": dict(COUNTERPARTS),
        "decisive_arms": list(DECISIVE_ARMS),
        "checks": {name: list(arms) for name, arms in CHECKS.items()},
        "fit": {
            "anchor": "every arm at forced_components = the plasmode's fidelity rank (no cross-validation)",
            "plant": "pls_cv cross-validates once; every other arm at forced_components = its rank",
        },
        "residual_correction": RESIDUAL_CORRECTION,
        "orientation_estimator": ORIENTATION_ESTIMATOR,
        "within_stage": (
            "per draw and per stage in declared order, idx[rows] = rows[rng.permutation(rows.size)] with "
            "rng = default_rng(the unit's RRPP seed); RRPP's reduced fit, df scale and statistics otherwise"
        ),
        "out_of_fold": {
            "folds": FOLDS,
            "stratified_by": "stage",
            "seed": f"sha256 of [{DRIVER_NAME!r}, {FOLD_BASE_SEED}, replicate, {FOLD_PURPOSE!r}]",
            "coordinates": "transform(X_held) . y_loadings^T of fit_plsda_model on the other folds' rows, at the "
                           "unit's rank and scaling",
        },
        "preflight_replicates": list(PREFLIGHT_REPLICATES),
        "alpha": ALPHA,
        "rejection": REJECTION,
        "reduced_alpha": REDUCED_ALPHA,
        "path_identity_tolerance": PATH_IDENTITY_TOLERANCE,
        "consistency_tolerance": CONSISTENCY_TOLERANCE,
    }


def enumerate_units(settings: Mapping[str, Any]) -> list[lsc.Unit]:
    """The pilot's plasmode units, replicate-major."""

    return psc.plasmode_units(settings["pilot"])


def unit_key(design: str, cell: str, replicate: int) -> str:
    return psc.unit_key(design, cell, replicate)


def record_key(record: Mapping[str, Any]) -> str:
    return unit_key(record["design"], record["cell"], record["replicate"])


def verify(settings: Mapping[str, Any], export: Any = None) -> None:
    """Refuse a changed parent file, or an export whose files differ from the settings', by name."""

    plsc.verify_parents(settings)
    if export is not None:
        plsc.verify_export(export, settings["pilot"])


# ── The inference path (design D4) and the within-stage draws (D5) ───────────


def reduced_fit(
    Y: Any, model_reduced: Any, residual_correction: str | None
) -> tuple[np.ndarray, np.ndarray]:
    """``RRPP``'s reduced fitted values and (corrected) residuals, computed as ``RRPP`` computes them."""

    reduced = np.asarray(model_reduced, dtype=float)
    scale = None if residual_correction is None else residual_correction_scale(reduced, residual_correction)
    frame = pd.DataFrame(Y)
    betas = estimate_betas(model_reduced, frame)
    fitted = np.matmul(model_reduced, betas)
    fitted.index = frame.index
    residuals = frame - fitted
    fitted_np, residuals_np = np.asarray(fitted), np.asarray(residuals)
    if scale is not None:
        residuals_np = residuals_np * scale[:, None]
    return fitted_np, residuals_np


def within_stage_indices(strata: Any, order: Sequence[str], permutations: int, seed: int) -> np.ndarray:
    """``permutations`` × n draws, each permuting the rows of every stratum among themselves.

    Per draw, the strata are visited in ``order`` and ``idx[rows] =
    rows[rng.permutation(rows.size)]``; with one stratum covering every row a
    draw is ``rng.permutation(n)``, ``RRPP``'s.
    """

    labels = np.asarray(strata).astype(str)
    groups = [np.flatnonzero(labels == str(stratum)) for stratum in order]
    covered = np.sort(np.concatenate(groups)) if groups else np.array([], dtype=int)
    if covered.size != labels.size or not np.array_equal(covered, np.arange(labels.size)):
        raise ValueError(f"Every row must belong to exactly one of the strata {list(order)}.")
    rng = np.random.default_rng(seed)
    draws = np.empty((int(permutations), labels.size), dtype=np.intp)
    for draw in range(int(permutations)):
        for rows in groups:
            draws[draw, rows] = rows[rng.permutation(rows.size)]
    return draws


def restricted_null(
    Y: Any,
    model_full: Any,
    model_reduced: Any,
    ls_means: Any,
    contrast: list[list[int]],
    indices: np.ndarray,
    *,
    orientation_estimator: str,
    residual_correction: str | None,
) -> tuple[list[np.ndarray], list[np.ndarray], list[np.ndarray]]:
    """``RRPP``'s serial null with the given draws: ``ŷ_red + e*[idx]`` and ``estimate_difference`` per draw."""

    fitted, residuals = reduced_fit(Y, model_reduced, residual_correction)
    deltas, angles, shapes = [], [], []
    for idx in np.asarray(indices):
        d, a, s = estimate_difference(
            fitted + residuals[idx, :], model_full, ls_means, contrast, orientation_estimator=orientation_estimator
        )
        deltas.append(d)
        angles.append(a)
        shapes.append(s)
    return deltas, angles, shapes


def infer(
    latent: Any,
    inputs: Any,
    *,
    permutations: int,
    seed: int,
    indices: np.ndarray | None = None,
    orientation_estimator: str = ORIENTATION_ESTIMATOR,
    residual_correction: str | None = RESIDUAL_CORRECTION,
) -> dict[str, Any]:
    """Statistics, p-values and the ``angle`` null 95th percentile of one latent matrix, by the analysis's rule.

    ``indices`` (``permutations`` × n) selects the within-stage null; without
    it the null is serial ``RRPP`` at ``seed``.
    """

    contrast = inputs.design.contrast
    shape_available = len(inputs.stage_levels) >= 3
    deltas, angles, shapes, _ = estimate_difference(
        Y=latent,
        model_matrix=inputs.model_full,
        LS_means=inputs.ls_means,
        contrast=contrast,
        return_spectra=True,
        orientation_estimator=orientation_estimator,
    )
    observed = _extract_pair_statistics(deltas, angles, shapes, shape_available=shape_available)
    if indices is None:
        dist = RRPP(
            latent, inputs.model_full, inputs.model_reduced, inputs.ls_means, contrast,
            permutations=int(permutations), n_jobs=1, progress=False, seed=int(seed),
            orientation_estimator=orientation_estimator, residual_correction=residual_correction,
        )
    else:
        if np.asarray(indices).shape[0] != int(permutations):
            raise ValueError(f"{np.asarray(indices).shape[0]} restricted draws for {permutations} permutations.")
        dist = restricted_null(
            latent, inputs.model_full, inputs.model_reduced, inputs.ls_means, contrast, indices,
            orientation_estimator=orientation_estimator, residual_correction=residual_correction,
        )
    null = _extract_null_distributions(*dist, shape_available=shape_available)
    p_values = {
        statistic: _empirical_p_value(null[statistic], value)
        for statistic, value in observed.items()
        if statistic in null and np.isfinite(value)
    }
    summary = _summarize_null_distributions(null)
    return {
        "statistics": observed,
        "p_values": p_values,
        "angle_null_q95": (summary.get("angle") or {}).get("q95"),
    }


def path_identity(
    paths: Mapping[str, Mapping[str, Any]], records: Mapping[str, Mapping[str, Any]], *,
    tolerance: float = PATH_IDENTITY_TOLERANCE,
) -> dict[str, Any]:
    """The path's statistics against the analysis's (to ``tolerance``) and its p-values (exactly), per arm."""

    arms: dict[str, Any] = {}
    worst = 0.0
    p_equal = True
    for arm, path in paths.items():
        record = records[arm]
        deviation = 0.0
        for statistic in STATISTICS:
            ours, theirs = _number(path["statistics"].get(statistic)), _number(record["statistics"].get(statistic))
            if ours is None and theirs is None:
                continue
            deviation = math.inf if ours is None or theirs is None else max(deviation, abs(ours - theirs))
        ours_p, theirs_p = path["p_values"], record["p_values"]
        differing = sorted(
            statistic for statistic in set(ours_p) | set(theirs_p)
            if _number(ours_p.get(statistic)) != _number(theirs_p.get(statistic))
        )
        passed = deviation <= tolerance and not differing
        arms[arm] = {"statistics_deviation": deviation, "p_values_differing": differing, "passed": passed}
        worst = max(worst, deviation)
        p_equal = p_equal and not differing
    applicable = bool(paths)
    return {
        "applicable": applicable,
        "passed": bool(applicable and worst <= tolerance and p_equal),
        "deviation": worst if applicable else None,
        "tolerance": tolerance,
        "p_values_equal": p_equal,
        "failing_arms": [arm for arm, entry in arms.items() if not entry["passed"]],
        "arms": arms,
    }


# ── Out-of-fold fitted stage indicators (design D6) ──────────────────────────


def stratified_folds(strata: Any, order: Sequence[str], seed: int, *, folds: int = FOLDS) -> np.ndarray:
    """One fold id per row: each stratum's rows shuffled and dealt round-robin, the start carried across strata."""

    labels = np.asarray(strata).astype(str)
    rng = np.random.default_rng(int(seed))
    assignment = np.full(labels.size, -1, dtype=int)
    offset = 0
    for stratum in order:
        rows = np.flatnonzero(labels == str(stratum))
        shuffled = rows[rng.permutation(rows.size)]
        assignment[shuffled] = (np.arange(rows.size) + offset) % int(folds)
        offset = (offset + rows.size) % int(folds)
    if (assignment < 0).any():
        raise ValueError(f"Rows outside the declared strata {list(order)}.")
    return assignment


def fitted_indicators(
    features: np.ndarray, labels: np.ndarray, train: np.ndarray, rows: np.ndarray, *, rank: int, scale: bool
) -> tuple[np.ndarray, np.ndarray]:
    """``transform(X_rows) · y_loadingsᵀ`` and the predicted indicators of a stage-PLS fitted on ``train``."""

    features = np.asarray(features, dtype=float)
    labels = np.asarray(labels).astype(str)
    model = fit_plsda_model(features[train], labels[train], n_components=int(rank), scale=bool(scale))
    loadings = np.asarray(model.y_loadings_, dtype=float)
    held = features[rows]
    return np.asarray(model.transform(held), dtype=float) @ loadings.T, np.asarray(model.predict(held), dtype=float)


def out_of_fold_fitted(
    features: np.ndarray, labels: Any, assignment: np.ndarray, *, rank: int, scale: bool
) -> dict[str, Any]:
    """The n × K out-of-fold fitted stage indicators, the held-out predictions and the training stage counts."""

    features = np.asarray(features, dtype=float)
    labels = np.asarray(labels).astype(str)
    categories = sorted(set(labels.tolist()))
    matrix = np.full((labels.size, len(categories)), np.nan)
    predictions = np.full((labels.size, len(categories)), np.nan)
    training_counts = []
    for fold in sorted(set(np.asarray(assignment).tolist())):
        rows = np.asarray(assignment) == fold
        train = ~rows
        present = sorted(set(labels[train].tolist()))
        if present != categories:
            raise ValueError(f"Fold {fold}'s training set lacks the stage(s) {sorted(set(categories) - set(present))}.")
        matrix[rows], predictions[rows] = fitted_indicators(features, labels, train, rows, rank=rank, scale=scale)
        training_counts.append([int(np.sum(labels[train] == category)) for category in categories])
    if not np.isfinite(matrix).all():
        raise ValueError("Some rows were never held out.")
    return {"matrix": matrix, "predictions": predictions, "training_stage_counts": training_counts,
            "categories": categories}


def explained_share(labels: Any, predictions: np.ndarray) -> float:
    """The share of the one-hot stage indicators' variance ``predictions`` explain (categories sorted)."""

    names = np.asarray(labels).astype(str)
    categories = sorted(set(names.tolist()))
    onehot = (names[:, None] == np.asarray(categories)[None, :]).astype(float)
    centred = onehot - onehot.mean(axis=0)
    return float(1.0 - np.sum((onehot - np.asarray(predictions, dtype=float)) ** 2) / np.sum(centred**2))


def stage_dispersion(latent: Any, model_reduced: Any, strata: Any, order: Sequence[str]) -> dict[str, Any]:
    """Each stage's mean squared residual norm under the reduced model, over the pooled one, and max / min."""

    values = np.asarray(latent, dtype=float)
    reduced = np.asarray(model_reduced, dtype=float)
    residuals = values - reduced @ np.asarray(estimate_betas(reduced, values), dtype=float)
    squared = np.sum(residuals**2, axis=1)
    labels = np.asarray(strata).astype(str)
    pooled = float(squared.mean())
    ratios = {f"stage_{index}": float(squared[labels == str(stage)].mean() / pooled) for index, stage in
              enumerate(order)}
    values_ = list(ratios.values())
    return {**ratios, "max_min_ratio": float(max(values_) / min(values_)) if min(values_) > 0 else None}


# ── One unit: fourteen arms, the checks and the readings ─────────────────────


def analyse_without_covariates(blocks: Mapping[str, pd.DataFrame], metadata: pd.DataFrame, export: Any,
                               params: Any) -> tuple[Any, Any]:
    """``plasmode_check.analyse`` with no covariates in the trajectory model."""

    from motco.analysis import OmicsBlock, TrajectoryDataset, analyze_trajectories

    dataset = TrajectoryDataset(
        blocks=[OmicsBlock(name, blocks[name], "identity") for name in pc.BLOCKS],
        metadata=metadata,
        sample_col=export.sample_col,
        group_col=export.group_col,
        stage_col=export.stage_col,
        stage_order=list(export.stage_order),
    )
    return dataset, analyze_trajectories(dataset, params)


def reference_latent(result: Any) -> pd.DataFrame:
    """The latent matrix an analysis measured in (its ``latent_scores``)."""

    return result.latent_scores


def _number(value: Any) -> float | None:
    return plsc._number(value)


def _path_record(arm: str, inferred: Mapping[str, Any], *, rank: int, dimension: int, seconds: float,
                 fit: str) -> dict[str, Any]:
    return {
        "status": "completed",
        **ARM_SPECS[arm],
        "rank": int(rank),
        "measured_dimension": int(dimension),
        "fit": fit,
        "pair_statistics": lsc._finite(inferred["statistics"]),
        "p_values": lsc._finite(inferred["p_values"]),
        "angle_null_q95": inferred["angle_null_q95"],
        "seconds": seconds,
    }


def _analysis_record(arm: str, result: Any, seconds: float, *, fit: str) -> dict[str, Any]:
    record = psc.analysis_arm_record(result, seconds)
    record.pop("span_capture", None)
    record.pop("population_latent", None)
    return {**record, **ARM_SPECS[arm], "fit": fit}


def _failed_arm(arm: str, exc: Exception, seconds: float) -> dict[str, Any]:
    return {**psc._failed_arm(exc, seconds), **ARM_SPECS[arm]}


def evaluate_unit(unit: lsc.Unit, settings: Mapping[str, Any], export: Any) -> dict[str, Any]:
    """Build the unit once, fit production once and evaluate the fourteen arms; one record of seeds and statistics."""

    start = perf_counter()
    pilot = settings["pilot"]
    parent = pilot["plasmode"]
    permutations = int(pilot["permutations"])
    seeds = pc.unit_seeds(int(parent["base_seed"]), str(parent["stream"]), unit.replicate)
    anchor = unit.cell == ANCHOR_CELL
    base: dict[str, Any] = {
        "design": unit.design,
        "cell": unit.cell,
        "trajectory_mode": "none" if anchor else plsc.PLANT_MODE,
        "effect_size": plsc.plasmode_effect(pilot, unit.cell),
        "replicate": unit.replicate,
        "stream": parent["stream"],
        "seeds": seeds,
        "rrpp_seed": seeds["rrpp"],
        "fold_seed": fold_seed(unit.replicate),
        "settings_hash": lsc.settings_hash(settings),
    }
    timing: dict[str, float] = {}
    try:
        mark = perf_counter()
        built = plsc.build_plasmode_unit(unit, pilot, export)
        base["label_redraws"] = built["label_redraws"]
        timing["generate"] = perf_counter() - mark
        production_params = psc.plasmode_production_params(pilot, built["seeds"])
        reference_params = (
            dataclasses.replace(production_params, forced_components=int(settings["anchor_rank"]))
            if anchor else production_params
        )
        mark = perf_counter()
        dataset, production = pc.analyse(built["blocks"], built["metadata"], export, reference_params)
        timing["production"] = perf_counter() - mark
        latent = production.record["latent"]
        rank = int(latent["selected_lv"])
        inputs = analysis._trajectory_inputs(dataset, production_params)
        order = [str(index) for index in range(len(inputs.stage_levels))]
        strata = inputs.working[analysis._STAGE].astype(str).to_numpy()
        artifacts = production.latent_artifacts
        model = artifacts.model
        features = artifacts.features.to_numpy(dtype=float)
        scale = bool(getattr(model, "scale", True))
    except Exception as exc:  # recorded, never aborts the shard
        return lsc._jsonable(psc._failed_unit(base, exc, start, timing))

    seed = int(seeds["rrpp"])
    n_stages = len(order)
    fit = "forced_at_fidelity_rank" if anchor else "forced_at_unit_rank"
    base.update({
        "rank": rank,
        "component_selection": latent.get("component_selection"),
        "cv_mean_auroc": latent.get("cv_mean_auroc"),
    })
    arms: dict[str, Any] = {
        REFERENCE_ARM: _analysis_record(
            REFERENCE_ARM, production, timing["production"],
            fit="forced_at_fidelity_rank" if anchor else "cross_validated",
        )
    }
    results: dict[str, Any] = {REFERENCE_ARM: production}
    for arm in PREDICTIVE_ARMS:
        mark = perf_counter()
        try:
            params = psc.plasmode_arm_params(arm, production_params, n_stages=n_stages, production_rank=rank)
            _, results[arm] = pc.analyse(built["blocks"], built["metadata"], export, params)
            arms[arm] = _analysis_record(arm, results[arm], perf_counter() - mark, fit=fit)
        except Exception as exc:  # recorded on the arm, the unit goes on
            arms[arm] = _failed_arm(arm, exc, perf_counter() - mark)

    # The checks that every predictive arm measures the production fit (the pilot's five).
    latents = {arm: reference_latent(result) for arm, result in results.items()}
    checks = psc.consistency_checks(
        model=model,
        scores=production.latent_scores.to_numpy(dtype=float),
        n_stages=n_stages,
        candidates={arm: latents[arm].to_numpy(dtype=float) for arm in PREDICTIVE_ARMS if arm in latents},
        arms=arms,
        with_capture=False,
        models={arm: results[arm].latent_artifacts.model for arm in PREDICTIVE_ARMS if arm in results},
    )

    # Path identity: the driver's path under unrestricted draws against each reference analysis.
    mark = perf_counter()
    paths: dict[str, Any] = {}
    for arm in REFERENCE_ARMS:
        if arm in latents:
            try:
                paths[arm] = infer(latents[arm], inputs, permutations=permutations, seed=seed)
            except Exception as exc:  # an inference path that raises fails its own identity
                paths[arm] = {"statistics": {}, "p_values": {}, "error": f"{type(exc).__name__}: {exc}"}
    checks["path_identity"] = path_identity(paths, {arm: results[arm].record for arm in paths})
    timing["path_identity"] = perf_counter() - mark

    # Within-stage draws: one restricted sequence, shared by every within-stage arm of the unit.
    indices = within_stage_indices(strata, order, permutations, seed)
    for source in REFERENCE_ARMS:
        arm = f"{source}{WITHIN_STAGE}"
        mark = perf_counter()
        try:
            if source not in latents:
                raise RuntimeError(f"the reference arm {source} did not complete")
            inferred = infer(latents[source], inputs, permutations=permutations, seed=seed, indices=indices)
            arms[arm] = _path_record(arm, inferred, rank=rank, dimension=latents[source].shape[1],
                                     seconds=perf_counter() - mark, fit=fit if source != REFERENCE_ARM else
                                     arms[REFERENCE_ARM]["fit"])
        except Exception as exc:  # recorded on the arm, the unit goes on
            arms[arm] = _failed_arm(arm, exc, perf_counter() - mark)

    # Out-of-fold fitted stage indicators, under both draws.
    mark = perf_counter()
    stage_fit: dict[str, Any] = {"in_sample_r2": psc.stage_r2(model, features, strata.tolist())}
    folds_record: dict[str, Any] = {"seed": base["fold_seed"], "count": FOLDS}
    out_of_fold: pd.DataFrame | None = None
    out_of_fold_error: Exception | None = None
    try:
        assignment = stratified_folds(strata, order, base["fold_seed"])
        built_oof = out_of_fold_fitted(features, strata, assignment, rank=rank, scale=scale)
        out_of_fold = pd.DataFrame(
            built_oof["matrix"], index=latents[REFERENCE_ARM].index,
            columns=[f"fitted_{index}" for index in range(built_oof["matrix"].shape[1])],
        )
        stage_fit["out_of_fold_share"] = explained_share(strata, built_oof["predictions"])
        folds_record["training_stage_counts"] = built_oof["training_stage_counts"]
    except Exception as exc:  # both out-of-fold arms fail on this unit
        out_of_fold_error = exc
        stage_fit["out_of_fold_share"] = None
    timing["out_of_fold_fit"] = perf_counter() - mark
    for arm, restricted in ((OUT_OF_FOLD_ARM, None), (f"{OUT_OF_FOLD_ARM}{WITHIN_STAGE}", indices)):
        mark = perf_counter()
        try:
            if out_of_fold is None:
                raise RuntimeError(f"the out-of-fold coordinates failed: {out_of_fold_error}")
            inferred = infer(out_of_fold, inputs, permutations=permutations, seed=seed, indices=restricted)
            arms[arm] = _path_record(arm, inferred, rank=rank, dimension=out_of_fold.shape[1],
                                     seconds=perf_counter() - mark, fit="out_of_fold_at_unit_rank")
        except Exception as exc:  # recorded on the arm, the unit goes on
            arms[arm] = _failed_arm(arm, exc, perf_counter() - mark)

    # The covariate-free arms: the same fit (covariate-blind), no covariates in the trajectory model.
    free_models: dict[str, Any] = {}
    for source in REFERENCE_ARMS:
        arm = f"{source}{NO_COVARIATES}"
        mark = perf_counter()
        try:
            params = (
                dataclasses.replace(production_params, forced_components=rank)
                if source == REFERENCE_ARM
                else psc.plasmode_arm_params(source, production_params, n_stages=n_stages, production_rank=rank)
            )
            _, result = analyse_without_covariates(built["blocks"], built["metadata"], export, params)
            arms[arm] = _analysis_record(arm, result, perf_counter() - mark, fit=fit)
            free_models[arm] = result.latent_artifacts.model
            if source in latents:
                difference = np.abs(result.latent_scores.to_numpy(dtype=float) - latents[source].to_numpy(dtype=float))
                arms[arm]["latent_deviation"] = float(difference.max()) if difference.size else 0.0
        except Exception as exc:  # recorded on the arm, the unit goes on
            arms[arm] = _failed_arm(arm, exc, perf_counter() - mark)
    free_deviation = 0.0
    for arm, fitted_model in free_models.items():
        free_deviation = max(free_deviation, psc._model_deviation(model, fitted_model))
        if arms[arm].get("rank") != rank:
            free_deviation = math.inf
    checks["same_fit_covariate_free"] = psc._check(
        free_deviation if free_models else None, CONSISTENCY_TOLERANCE, applicable=bool(free_models)
    )

    # Readings beside the arms: residual dispersion by stage per coordinate system.
    systems: dict[str, Any] = {}
    for system, matrix in ((SCORES, latents.get(REFERENCE_ARM)),
                           *((arm, latents.get(arm)) for arm in PREDICTIVE_ARMS), (OUT_OF_FOLD, out_of_fold)):
        if matrix is not None:
            systems[system] = {"dispersion": stage_dispersion(matrix, inputs.model_reduced, strata, order)}

    failing = [name for name in CONSISTENCY_CHECKS if not checks[name]["passed"]]
    timing["unit"] = perf_counter() - start
    record = {
        **base,
        "status": "inconsistent" if failing else "completed",
        "inconsistent_checks": failing,
        "consistency": checks,
        "stage_fit": stage_fit,
        "folds": folds_record,
        "coordinate_systems": systems,
        "arms": {arm: arms[arm] for arm in ARMS if arm in arms},
        "runtime_seconds": timing,
    }
    return lsc._jsonable(record)


# ── Pre-flight (design D3) ───────────────────────────────────────────────────


def preflight_unit(replicate: int, settings: Mapping[str, Any], export: Any) -> dict[str, Any]:
    """An anchor unit's cross-validated ``pls_cv`` against the forced call: rank, AUROC, statistics, p-values."""

    pilot = settings["pilot"]
    unit = lsc.Unit(DESIGN, ANCHOR_CELL, int(replicate))
    built = plsc.build_plasmode_unit(unit, pilot, export)
    production = psc.plasmode_production_params(pilot, built["seeds"])
    mark = perf_counter()
    _, cross_validated = pc.analyse(built["blocks"], built["metadata"], export, production)
    cv_seconds = perf_counter() - mark
    mark = perf_counter()
    _, forced = pc.analyse(built["blocks"], built["metadata"], export,
                           dataclasses.replace(production, forced_components=int(settings["anchor_rank"])))
    forced_seconds = perf_counter() - mark
    fidelity = plsc.fidelity_outcome(cross_validated, pilot["fidelity"])
    differing = {
        key: sorted(
            statistic
            for statistic in set(cross_validated.record[key]) | set(forced.record[key])
            if cross_validated.record[key].get(statistic) != forced.record[key].get(statistic)
        )
        for key in ("statistics", "p_values")
    }
    forced_rank = int(forced.record["latent"]["selected_lv"])
    passed = bool(fidelity["passed"]) and forced_rank == int(settings["anchor_rank"]) and not any(differing.values())
    return {
        "unit": unit_key(*unit.key),
        "rrpp_seed": built["seeds"]["rrpp"],
        "fidelity": fidelity,
        "forced_rank": forced_rank,
        "statistics_differing": differing["statistics"],
        "p_values_differing": differing["p_values"],
        "cross_validated": {"statistics": cross_validated.record["statistics"],
                            "p_values": cross_validated.record["p_values"]},
        "forced": {"statistics": forced.record["statistics"], "p_values": forced.record["p_values"]},
        "seconds": {"cross_validated": cv_seconds, "forced": forced_seconds},
        "passed": passed,
    }


def preflight(out_dir: Path, settings: Mapping[str, Any], export: Any) -> dict[str, Any]:
    """Write the settings and ``preflight.json``; the run waits on its ``passed``."""

    out_dir = Path(out_dir)
    if (out_dir / "settings.json").exists():
        verify(lsc.load_settings(out_dir))
    verify(settings, export)
    lsc.ensure_settings(out_dir, settings)
    units = []
    for replicate in settings["preflight_replicates"]:
        try:
            units.append(preflight_unit(int(replicate), settings, export))
        except Exception as exc:  # a pre-flight that raises fails
            units.append({"unit": unit_key(DESIGN, ANCHOR_CELL, int(replicate)), "passed": False,
                          "error": f"{type(exc).__name__}: {exc}"})
    result = {
        "settings_hash": lsc.settings_hash(settings),
        "rule": "the cross-validated pls_cv selects the fidelity rank and AUROC, and its statistics and p-values "
                "equal the forced call's exactly",
        "expected": {"rank": int(settings["anchor_rank"]), **settings["pilot"]["fidelity"]},
        "units": units,
        "passed": bool(units) and all(unit["passed"] for unit in units),
    }
    lsc._write_json(out_dir / PREFLIGHT_NAME, result)
    return result


def load_preflight(out_dir: Path, settings: Mapping[str, Any]) -> dict[str, Any] | None:
    """The pre-flight record written under these settings, or ``None``."""

    path = Path(out_dir) / PREFLIGHT_NAME
    if not path.exists():
        return None
    record = json.loads(path.read_text(encoding="utf-8"))
    return record if record.get("settings_hash") == lsc.settings_hash(settings) else None


# ── Shards, resume and the run ───────────────────────────────────────────────


def read_records(paths: Iterable[Path]) -> list[dict[str, Any]]:
    return psc.read_records(paths)


_WORKER_EXPORT: Any = None


def _init_worker(directory: str, files: Mapping[str, str]) -> None:
    global _WORKER_EXPORT
    export = pc.load_export(Path(directory))
    if export.files != dict(files):
        raise SystemExit(f"{directory} changed after the run started (file hashes differ).")
    _WORKER_EXPORT = export


def _evaluate_pending(unit: lsc.Unit, settings: Mapping[str, Any]) -> dict[str, Any]:
    return evaluate_unit(unit, settings, _WORKER_EXPORT)


def run_shard(
    out_dir: Path,
    settings: Mapping[str, Any],
    index: int,
    count: int,
    *,
    export: Any,
    limit: int | None = None,
    workers: int = 1,
) -> int:
    """Evaluate this shard's missing units into its file; returns how many were evaluated."""

    out_dir = Path(out_dir)
    if (out_dir / "settings.json").exists():
        verify(lsc.load_settings(out_dir))
    verify(settings, export)
    lsc.ensure_settings(out_dir, settings)
    checked = load_preflight(out_dir, settings)
    if checked is None:
        raise SystemExit(f"{out_dir / PREFLIGHT_NAME} is missing or was written under other settings; run "
                         "`preflight` first.")
    if not checked.get("passed"):
        raise SystemExit(f"The pre-flight failed ({out_dir / PREFLIGHT_NAME}); the run stops.")
    path = lsc.shard_path(out_dir, index, count)
    path.parent.mkdir(parents=True, exist_ok=True)
    done = {record_key(record): record for record in (read_records([path]) if path.exists() else [])}
    # Rewritten with complete units only, so a torn line is re-run, not duplicated.
    with path.open("w", encoding="utf-8") as handle:
        handle.write("".join(lsc._line(done[key]) for key in sorted(done)))
    pending = [
        unit
        for position, unit in enumerate(enumerate_units(settings))
        if position % count == index and unit_key(*unit.key) not in done
    ]
    if limit is not None:
        pending = pending[:limit]

    def write(position: int, unit: lsc.Unit, record: Mapping[str, Any]) -> None:
        with path.open("a", encoding="utf-8") as handle:
            handle.write(lsc._line(record))
            handle.flush()
        print(f"[{index}/{count}] {position}/{len(pending)} {unit.key} {record['status']}", file=sys.stderr, flush=True)

    if workers <= 1:
        for position, unit in enumerate(pending, start=1):
            write(position, unit, evaluate_unit(unit, settings, export))
    else:
        with ProcessPoolExecutor(
            max_workers=workers, initializer=_init_worker, initargs=(str(export.directory), dict(export.files))
        ) as pool:
            futures = {pool.submit(_evaluate_pending, unit, settings): unit for unit in pending}
            for position, future in enumerate(as_completed(futures), start=1):
                write(position, futures[future], future.result())
    return len(pending)


# ── Merge: rate, paired and cell tables, accounting, reproduction ────────────


RATE_COLUMNS = (
    "design", "cell", "trajectory_mode", "effect_size", "arm", "statistic", "rate", "se", "count", "rejections",
    "failures", "inconsistent",
)
PAIRED_COLUMNS = (
    "design", "cell", "trajectory_mode", "effect_size", "statistic", "arm", "reference_arm", "relation", "n_paired",
    "rate_reference", "rate_arm", "difference", "b", "c", "p_fewer", "p_more", "p_two_sided",
)
STAGE_COLUMNS: tuple[str, ...] = tuple(f"stage_{index}" for index in range(len(pc.STAGE_ORDER)))
CELL_COLUMNS = (
    "design", "cell", "trajectory_mode", "effect_size", "arm", "coordinates", "draws", "covariates", "n",
    "rank_median", "rank_q10", "rank_q90",
    *(f"dispersion_{stage}_median" for stage in STAGE_COLUMNS), "dispersion_ratio_median",
    "stage_r2_median", "out_of_fold_share_median",
    *(f"observed_{statistic}_median" for statistic in STATISTICS),
    "angle_null_q95_median",
)


def pairs() -> list[tuple[str, str, str]]:
    """(arm, reference arm, relation): every non-reference arm against its counterpart, every arm against pls_cv."""

    out = [(arm, COUNTERPARTS[arm], "counterpart") for arm in ARMS if arm in COUNTERPARTS]
    out += [
        (arm, REFERENCE_ARM, "production")
        for arm in ARMS
        if arm != REFERENCE_ARM and COUNTERPARTS.get(arm) != REFERENCE_ARM
    ]
    return out


def _cell_fields(record: Mapping[str, Any]) -> dict[str, Any]:
    return psc._cell_fields(record)


def _p(record: Mapping[str, Any], arm: str, statistic: str) -> Any:
    return psc._p(record, arm, statistic)


def rate_table(records: Sequence[Mapping[str, Any]], settings: Mapping[str, Any]) -> list[dict[str, Any]]:
    """Rejections (p ≤ α) per cell × arm × statistic; failed and inconsistent units excluded."""

    return psc.rate_table(records, {"alpha": settings["alpha"], "rejection": settings["rejection"],
                                    "arms": settings["arms"]})


def paired_table(records: Sequence[Mapping[str, Any]], settings: Mapping[str, Any]) -> list[dict[str, Any]]:
    """Every pair of :func:`pairs`, per cell and statistic, over units available in both arms."""

    alpha, rule = float(settings["alpha"]), str(settings["rejection"])
    rows = []
    for _, group in sorted(psc._by_cell(records).items()):
        fields = _cell_fields(group[0])
        ordered = sorted(group, key=lambda record: int(record["replicate"]))
        for statistic in STATISTICS:
            for arm, reference, relation in pairs():
                decisions = [
                    (plsc._rejects(_p(record, reference, statistic), alpha, rule),
                     plsc._rejects(_p(record, arm, statistic), alpha, rule))
                    for record in ordered
                ]
                both = [(ref, cand) for ref, cand in decisions if ref is not None and cand is not None]
                n = len(both)
                b = sum(1 for ref, cand in both if ref and not cand)
                c = sum(1 for ref, cand in both if cand and not ref)
                exact = plsc.paired_p_values(b, c)
                rows.append({
                    **fields, "statistic": statistic, "arm": arm, "reference_arm": reference, "relation": relation,
                    "n_paired": n,
                    "rate_reference": sum(ref for ref, _ in both) / n if n else None,
                    "rate_arm": sum(cand for _, cand in both) / n if n else None,
                    "difference": (c - b) / n if n else None,
                    "b": b, "c": c,
                    "p_fewer": exact["p_loss"], "p_more": exact["p_gain"], "p_two_sided": exact["p_two_sided"],
                })
    return rows


def cell_table(records: Sequence[Mapping[str, Any]], settings: Mapping[str, Any]) -> list[dict[str, Any]]:
    """Per cell and arm: rank, the arm's coordinate system's dispersion by stage, the stage fit, null width."""

    quantile = plsc._quantile
    rows = []
    for _, group in sorted(psc._by_cell(records).items()):
        fields = _cell_fields(group[0])
        units = [record for record in group if record.get("status") == "completed"]
        for arm in settings["arms"]:
            spec = ARM_SPECS[arm]
            done = [(record, entry) for record in units if (entry := psc._arm(record, arm)) is not None]
            dispersions = [
                ((record.get("coordinate_systems") or {}).get(spec["coordinates"]) or {}).get("dispersion") or {}
                for record, _ in done
            ]
            row: dict[str, Any] = {**fields, "arm": arm, **{key: spec[key] for key in ("coordinates", "draws",
                                                                                       "covariates")},
                                   "n": len(done)}
            ranks = [entry.get("rank") for _, entry in done]
            for suffix, q in (("median", 0.5), ("q10", 0.1), ("q90", 0.9)):
                row[f"rank_{suffix}"] = quantile(ranks, q)
            for stage in STAGE_COLUMNS:
                row[f"dispersion_{stage}_median"] = quantile((d.get(stage) for d in dispersions), 0.5)
            row["dispersion_ratio_median"] = quantile((d.get("max_min_ratio") for d in dispersions), 0.5)
            row["stage_r2_median"] = quantile(((r.get("stage_fit") or {}).get("in_sample_r2") for r, _ in done), 0.5)
            row["out_of_fold_share_median"] = quantile(
                ((r.get("stage_fit") or {}).get("out_of_fold_share") for r, _ in done), 0.5
            )
            for statistic in STATISTICS:
                row[f"observed_{statistic}_median"] = quantile(
                    ((entry.get("pair_statistics") or {}).get(statistic) for _, entry in done), 0.5
                )
            row["angle_null_q95_median"] = quantile((entry.get("angle_null_q95") for _, entry in done), 0.5)
            rows.append(row)
    return rows


def accounting(
    records: Sequence[Mapping[str, Any]], settings: Mapping[str, Any], *, expected: int, duplicates: int
) -> dict[str, Any]:
    """Units, failures and inconsistent units (check, arm, deviation) by cell, arm failures, timing."""

    summary = psc.accounting(records, {"profile": "plasmode", "arms": settings["arms"], **settings},
                             expected=expected, duplicates=duplicates)
    summary.pop("profile", None)
    summary["settings_hash"] = lsc.settings_hash(settings)
    details = []
    production: list[float] = []
    for record in records:
        timing = record.get("runtime_seconds") or {}
        if timing.get("production") is not None:
            production.append(float(timing["production"]))
        if record.get("status") != "inconsistent":
            continue
        entry: dict[str, Any] = {"unit": record_key(record), "checks": {}}
        for name in record.get("inconsistent_checks", []):
            check = record["consistency"][name]
            entry["checks"][name] = {"deviation": check.get("deviation")}
            if name == "path_identity":
                entry["checks"][name]["arms"] = {
                    arm: check["arms"][arm] for arm in check.get("failing_arms", [])
                }
        details.append(entry)
    summary["inconsistent_unit_checks"] = details
    summary.pop("fit_seconds_median", None)
    summary["production_seconds_median"] = plsc._quantile(production, 0.5)
    return summary


def _committed_counts(rows: Sequence[Mapping[str, Any]]) -> dict[tuple[str, str, str], dict[str, int]]:
    """The pilot's committed reference counts per (cell, arm, statistic)."""

    return {
        (str(row["cell"]), str(row["arm"]), str(row["statistic"])): {
            "count": int(float(row["count"])), "rejections": int(float(row["rejections"]))
        }
        for row in rows
        if str(row["design"]) == DESIGN and str(row["arm"]) in REFERENCE_ARMS
    }


def reproduction(
    rates: Sequence[Mapping[str, Any]], settings: Mapping[str, Any], committed_rows: Sequence[Mapping[str, Any]]
) -> dict[str, Any]:
    """The four reference arms' available and rejected counts per cell and statistic against the pilot's, exact."""

    committed = _committed_counts(committed_rows)
    ours = _committed_counts(rates)
    cells = [str(name) for name, _ in settings["pilot"]["cells"]]
    keys = [(cell, arm, statistic) for cell in cells for arm in REFERENCE_ARMS for statistic in STATISTICS]
    rows, differing = [], []
    for key in keys:
        theirs, mine = committed.get(key), ours.get(key)
        entry = {"cell": key[0], "arm": key[1], "statistic": key[2], "committed": theirs, "reproduced": mine,
                 "matches": theirs is not None and theirs == mine}
        rows.append(entry)
        if not entry["matches"]:
            differing.append(entry)
    return {
        "check": "pilot_reference_counts",
        "source": dict(settings["parents"]["pilot_plasmode_rates"]),
        "rule": "each reference arm's available and rejected counts (p <= alpha) per cell and statistic equal the "
                "pilot's committed rate_table.csv",
        "rows_compared": len(rows),
        "rows_matching": len(rows) - len(differing),
        "passed": bool(rows) and not differing,
        "differing": differing,
        "rows": rows,
    }


def merge(out_dir: Path) -> dict[str, Any]:
    """Merge the shards and write the tables, the accounting and the reproduction."""

    out_dir = Path(out_dir)
    settings = lsc.load_settings(out_dir)
    verify(settings)
    raw = read_records(sorted((out_dir / "records").glob("shard_*.jsonl")))
    stale = sorted({record_key(r) for r in raw if r.get("settings_hash") != lsc.settings_hash(settings)})
    if stale:
        raise SystemExit(f"{len(stale)} unit(s) were written under other settings, e.g. {stale[:3]}; refusing.")
    seen: dict[str, int] = {}
    by_key: dict[str, dict[str, Any]] = {}
    for record in raw:
        key = record_key(record)
        seen[key] = seen.get(key, 0) + 1
        by_key[key] = record
    expected = [unit_key(*unit.key) for unit in enumerate_units(settings)]
    records = [by_key[key] for key in expected if key in by_key]
    with (out_dir / "merged.jsonl").open("w", encoding="utf-8") as handle:
        handle.write("".join(lsc._line(record) for record in records))
    rates = rate_table(records, settings)
    lsc._write_csv(out_dir / "rate_table.csv", rates, RATE_COLUMNS)
    lsc._write_csv(out_dir / "paired_table.csv", paired_table(records, settings), PAIRED_COLUMNS)
    lsc._write_csv(out_dir / "cell_table.csv", cell_table(records, settings), CELL_COLUMNS)
    summary = accounting(records, settings, expected=len(expected), duplicates=sum(n - 1 for n in seen.values()))
    committed = lsc._read_csv(plsc._resolve(settings["parents"]["pilot_plasmode_rates"]["path"]))
    reproduced = reproduction(lsc._read_csv(out_dir / "rate_table.csv"), settings, committed)
    lsc._write_json(out_dir / "accounting.json", summary)
    lsc._write_json(out_dir / "reproduction.json", reproduced)
    return {
        **{key: summary[key] for key in ("complete_units", "missing_units", "failed_units", "inconsistent_units",
                                         "duplicate_records", "recorded_core_hours")},
        "reproduction_passed": reproduced["passed"],
        "reproduction_differing": len(reproduced["differing"]),
    }


# ── Decision (pure functions, design D8 and D9) ──────────────────────────────


def _index_rates(rows: Iterable[Mapping[str, Any]]) -> dict[tuple[str, str, str], Mapping[str, Any]]:
    return {(str(row["cell"]), str(row["arm"]), str(row["statistic"])): row for row in rows}


def _index_paired(rows: Iterable[Mapping[str, Any]]) -> dict[tuple[str, str, str, str], Mapping[str, Any]]:
    return {
        (str(row["cell"]), str(row["statistic"]), str(row["arm"]), str(row["reference_arm"])): row for row in rows
    }


def observation(rates: Mapping[tuple[str, str, str], Any], cell: str, arm: str, statistic: str) -> dict[str, Any]:
    """One reading for ``judge_rates``: available and rejected counts, or ``status: unavailable``."""

    entry = {"cell": cell, "arm": arm, "statistic": statistic}
    row = rates.get((cell, arm, statistic))
    if row is None or int(float(row["count"])) <= 0:
        return {**entry, "status": "unavailable"}
    return {**entry, "n_available": int(float(row["count"])), "n_rejected": int(float(row["rejections"]))}


def paired_reading(
    paired: Mapping[tuple[str, str, str, str], Any], cell: str, statistic: str, arm: str, reference: str
) -> dict[str, Any]:
    """One exact paired comparison, recomputed from the discordant counts (the table's p-values are rounded)."""

    entry = {"cell": cell, "statistic": statistic, "arm": arm, "reference_arm": reference}
    row = paired.get((cell, statistic, arm, reference))
    if row is None or int(float(row["n_paired"])) <= 0:
        return {**entry, "status": "unavailable"}
    b, c = int(float(row["b"])), int(float(row["c"]))
    exact = plsc.paired_p_values(b, c)
    return {
        **entry,
        "n_paired": int(float(row["n_paired"])),
        "rate_reference": _number(row.get("rate_reference")),
        "rate_arm": _number(row.get("rate_arm")),
        "b": b,
        "c": c,
        "p_fewer": exact["p_loss"],
        "p_more": exact["p_gain"],
        "p_two_sided": exact["p_two_sided"],
    }


def grade_arm(
    arm: str, rates: Mapping[tuple[str, str, str], Any], paired: Mapping[tuple[str, str, str, str], Any]
) -> dict[str, Any]:
    """``restored``, ``reduced`` or ``unchanged`` from the arm's three anchor readings and its paired tests."""

    readings = [observation(rates, ANCHOR_CELL, arm, statistic) for statistic in STATISTICS]
    unavailable = [row["statistic"] for row in readings if row.get("status") == "unavailable"]
    judged = grade.judge_rates(readings)
    counterpart = COUNTERPARTS[arm]
    tests = {
        statistic: paired_reading(paired, ANCHOR_CELL, statistic, arm, counterpart)
        for statistic in REDUCED_STATISTICS
    }
    fewer = [
        statistic for statistic, test in tests.items()
        if test.get("status") != "unavailable" and float(test["p_fewer"]) <= REDUCED_ALPHA
    ]
    if unavailable:
        verdict, reason = "unchanged", f"unavailable anchor reading(s): {unavailable}"
    elif judged["holds"]:
        verdict, reason = "restored", "the three anchor readings hold the count-plus-ceiling rule (m = 3)"
    elif fewer:
        named = ", ".join(f"{s} (b = {tests[s]['b']}, c = {tests[s]['c']}, p = {tests[s]['p_fewer']:.3g})"
                          for s in fewer)
        verdict, reason = "reduced", f"fewer rejections than {counterpart} on {named}"
    else:
        verdict, reason = "unchanged", "the size rule fails and no paired reduction reaches 0.025"
    return {
        "grade": verdict,
        "reason": reason,
        "counterpart": counterpart,
        "unavailable": unavailable,
        "reduced_statistics": fewer if verdict == "reduced" else [],
        "size_rule": judged,
        "paired": tests,
    }


def weakest(grades: Iterable[str]) -> str:
    """The weakest grade in the order ``unchanged`` < ``reduced`` < ``restored``."""

    values = list(grades)
    return min(values, key=GRADES.index) if values else "unchanged"


def verdict_of(w: str, o: str) -> str:
    """The verdict table of design D8 (the gate aside)."""

    if w == "restored" and o == "restored":
        return "both_restore"
    if w == "restored":
        return "null_restores"
    if o == "restored":
        return "out_of_fold_restores"
    if "reduced" in (w, o):
        return "partly_explained"
    return "unexplained"


def gate(
    reproduced: Mapping[str, Any], summary: Mapping[str, Any], checked: Mapping[str, Any] | None
) -> dict[str, Any]:
    """``reproduction_failed`` on a failed reproduction, any inconsistent unit, or a missing or failed pre-flight."""

    failures = []
    if not reproduced.get("passed"):
        failures.append(f"reproduction failed ({len(reproduced.get('differing', []))} differing row(s))")
    inconsistent = int(summary.get("inconsistent_units", 0) or 0)
    if inconsistent:
        failures.append(f"{inconsistent} inconsistent unit(s)")
    if checked is None:
        failures.append("the pre-flight record is missing")
    elif not checked.get("passed"):
        failures.append("the pre-flight did not pass")
    return {"passed": not failures, "failures": failures}


def decision(
    tables: Mapping[str, Sequence[Mapping[str, Any]]],
    *,
    reproduced: Mapping[str, Any],
    summary: Mapping[str, Any],
    checked: Mapping[str, Any] | None,
) -> dict[str, Any]:
    """The predeclared verdict (rule version 1) from the rate and paired tables; the descriptive readings beside it."""

    rates, paired = _index_rates(tables["rates"]), _index_paired(tables["paired"])
    arms = {arm: grade_arm(arm, rates, paired) for arm in DECISIVE_ARMS}
    checks = {
        name: {"grade": weakest(arms[arm]["grade"] for arm in members), "arms": list(members)}
        for name, members in CHECKS.items()
    }
    ungated = verdict_of(checks["W"]["grade"], checks["O"]["grade"])
    gated = gate(reproduced, summary, checked)
    verdict = ungated if gated["passed"] else "reproduction_failed"
    reduced_checks = [name for name, check in checks.items() if check["grade"] == "reduced"]
    return {
        "rule_version": RULE_VERSION,
        "verdict": verdict,
        "consequence": CONSEQUENCES[verdict],
        "consequences": dict(CONSEQUENCES),
        "gate": gated,
        "checks": checks,
        "reduced_checks": reduced_checks if verdict == "partly_explained" else [],
        "arms": arms,
        "rules": {
            "restored": "the arm's anchor delta, angle and shape rejections hold judge_rates (m = 3)",
            "reduced": f"not restored, and the one-sided exact paired p-value for fewer rejections than the "
                       f"counterpart is at most {REDUCED_ALPHA} on delta or shape",
            "unchanged": "otherwise; an unavailable anchor reading is unchanged",
            "check_W": "the weakest grade of its three arms (unchanged < reduced < restored)",
            "check_O": "the grade of fitted_out_of_fold",
        },
        "alpha": ALPHA,
        "rejection": REJECTION,
        "counterparts": {arm: COUNTERPARTS[arm] for arm in DECISIVE_ARMS},
        "descriptive": descriptive_tables(tables, arms),
        "notes": [
            "Only the four decisive arms are graded; every table under 'descriptive' is a reading beside the "
            "verdict and never changes it.",
            "The plant cells are a magnitude edit along the in-sample stage means at R = 100; they say whether a "
            "repair keeps the real-noise delta gain, not whether it holds size.",
        ],
    }


def _rate(rates: Mapping[tuple[str, str, str], Any], cell: str, arm: str, statistic: str) -> float | None:
    row = rates.get((cell, arm, statistic))
    return None if row is None else _number(row.get("rate"))


def descriptive_tables(
    tables: Mapping[str, Sequence[Mapping[str, Any]]], arms: Mapping[str, Mapping[str, Any]]
) -> dict[str, Any]:
    """Plant cells, the descriptive arms, the covariate-free arms, the dispersion and the null widths (D9)."""

    rates, paired = _index_rates(tables["rates"]), _index_paired(tables["paired"])
    plant = []
    for cell in PLANT_CELLS:
        for arm in ARMS:
            for statistic in STATISTICS:
                row = observation(rates, cell, arm, statistic)
                plant.append({
                    **row,
                    "rate": _rate(rates, cell, arm, statistic),
                    "against_pls_cv": None if arm == REFERENCE_ARM
                    else paired_reading(paired, cell, statistic, arm, REFERENCE_ARM),
                    "against_counterpart": None if arm not in COUNTERPARTS
                    else paired_reading(paired, cell, statistic, arm, COUNTERPARTS[arm]),
                })
    restored_plant = []
    for arm, graded in arms.items():
        if graded["grade"] != "restored":
            continue
        counterpart = COUNTERPARTS[arm]
        cells = {
            cell: {
                "rate_arm": _rate(rates, cell, arm, "delta"),
                "rate_pls_cv": _rate(rates, cell, REFERENCE_ARM, "delta"),
                "rate_counterpart": _rate(rates, cell, counterpart, "delta"),
            }
            for cell in PLANT_CELLS
        }
        at = cells.get(PRODUCTION_PLANT_CELL, {})
        arm_rate, production_rate = at.get("rate_arm"), at.get("rate_pls_cv")
        above = None if arm_rate is None or production_rate is None else arm_rate > production_rate
        restored_plant.append({
            "arm": arm, "counterpart": counterpart, "delta_rates": cells,
            "statement_cell": PRODUCTION_PLANT_CELL, "delta_above_production": above,
        })
    anchor_descriptive = []
    for arm in DESCRIPTIVE_ARMS[:2]:
        readings = [observation(rates, ANCHOR_CELL, arm, statistic) for statistic in STATISTICS]
        anchor_descriptive.append({
            "arm": arm,
            "counterpart": COUNTERPARTS[arm],
            "size_rule": grade.judge_rates(readings),
            "paired": {statistic: paired_reading(paired, ANCHOR_CELL, statistic, arm, COUNTERPARTS[arm])
                       for statistic in STATISTICS},
        })
    cells = [str(name) for name in (ANCHOR_CELL, *PLANT_CELLS)]
    covariate_free = [
        {
            "cell": cell, "arm": arm, "counterpart": COUNTERPARTS[arm], "statistic": statistic,
            "rate_arm": _rate(rates, cell, arm, statistic),
            "rate_counterpart": _rate(rates, cell, COUNTERPARTS[arm], statistic),
            "paired": paired_reading(paired, cell, statistic, arm, COUNTERPARTS[arm]),
        }
        for cell in cells for arm in COVARIATE_FREE_ARMS for statistic in STATISTICS
    ]
    seen: set[tuple[str, str]] = set()
    dispersion = []
    null_widths = []
    for row in tables.get("cells", []):
        cell, arm = str(row["cell"]), str(row["arm"])
        width = _number(row.get("angle_null_q95_median"))
        null_widths.append({"cell": cell, "arm": arm, "angle_null_q95_median": width})
        system = str(row.get("coordinates"))
        if (cell, system) in seen or str(row.get("draws")) != "unrestricted" or str(row.get("covariates")) != "six":
            continue
        seen.add((cell, system))
        dispersion.append({
            "cell": cell, "coordinates": system,
            **{f"dispersion_{stage}_median": _number(row.get(f"dispersion_{stage}_median")) for stage in STAGE_COLUMNS},
            "dispersion_ratio_median": _number(row.get("dispersion_ratio_median")),
            "stage_r2_median": _number(row.get("stage_r2_median")),
            "out_of_fold_share_median": _number(row.get("out_of_fold_share_median")),
        })
    return {
        "plant_cells": plant,
        "restored_arms_plant_delta": restored_plant,
        "anchor_descriptive_arms": anchor_descriptive,
        "covariate_free": covariate_free,
        "dispersion": dispersion,
        "null_widths": null_widths,
    }


_TABLE_FILES: dict[str, str] = {"rates": "rate_table.csv", "paired": "paired_table.csv", "cells": "cell_table.csv"}


def _read_json(path: Path) -> dict[str, Any]:
    if not path.exists():
        raise SystemExit(f"{path} does not exist; run `merge` first.")
    return json.loads(path.read_text(encoding="utf-8"))


def decide(out_dir: Path) -> dict[str, Any]:
    """Read the merged tables, the reproduction and the pre-flight; write the decision."""

    out_dir = Path(out_dir)
    settings = lsc.load_settings(out_dir)
    verify(settings)
    summary = _read_json(out_dir / "accounting.json")
    if summary.get("settings_hash") != lsc.settings_hash(settings):
        raise SystemExit(f"{out_dir}/accounting.json was merged under other settings; run `merge` again.")
    tables = {}
    for name, filename in _TABLE_FILES.items():
        path = out_dir / filename
        if not path.exists():
            raise SystemExit(f"{path} does not exist; run `merge` first.")
        tables[name] = lsc._read_csv(path)
    reproduced = _read_json(out_dir / "reproduction.json")
    checked = load_preflight(out_dir, settings)
    result = decision(tables, reproduced=reproduced, summary=summary, checked=checked)
    result["complete"] = {
        key: summary[key]
        for key in ("expected_units", "complete_units", "missing_units", "failed_units", "inconsistent_units")
    }
    result["reproduction"] = {key: reproduced.get(key) for key in ("passed", "rows_compared", "rows_matching")}
    result["preflight"] = None if checked is None else {
        "passed": checked.get("passed"),
        "units": [
            {key: unit.get(key) for key in ("unit", "passed", "forced_rank", "statistics_differing",
                                            "p_values_differing")}
            | {"fidelity": {key: (unit.get("fidelity") or {}).get(key) for key in ("rank", "cv_mean_auroc_rounded",
                                                                                   "passed")}}
            for unit in checked.get("units", [])
        ],
    }
    inputs = {name: plsc.sha256_file(out_dir / name)
              for name in ("settings.json", *_TABLE_FILES.values(), "accounting.json", "reproduction.json")}
    if (out_dir / PREFLIGHT_NAME).exists():
        inputs[PREFLIGHT_NAME] = plsc.sha256_file(out_dir / PREFLIGHT_NAME)
    result["inputs"] = inputs
    lsc._write_json(out_dir / DECISION_NAME, result)
    return result


# ── CLI ──────────────────────────────────────────────────────────────────────


def _unit_arguments(parser: argparse.ArgumentParser) -> None:
    parser.add_argument("--export-dir", type=Path, required=True, help="The SEA-AD cohort_30 export")
    parser.add_argument("--out-dir", type=Path, required=True)
    parser.add_argument("--reps", type=int, default=ANCHOR_REPLICATES,
                        help=f"Anchor replicates (default {ANCHOR_REPLICATES})")
    parser.add_argument("--plant-reps", type=int, default=PLANT_REPLICATES,
                        help=f"Replicates of each plant cell (default {PLANT_REPLICATES})")
    parser.add_argument("--perms", type=int, default=PERMUTATIONS, help=f"Permutations (default {PERMUTATIONS})")


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = parser.add_subparsers(dest="command", required=True)
    pre = sub.add_parser("preflight", help="Write the settings and check the forced anchor fit against the CV")
    _unit_arguments(pre)
    run = sub.add_parser("run", help="Evaluate (a shard of) the units")
    _unit_arguments(run)
    run.add_argument("--shard", type=lsc.parse_shard, default=(0, 1), help="i/N, round-robin over the units")
    run.add_argument("--limit", type=int, default=None, help="Evaluate at most this many pending units")
    run.add_argument("--workers", type=int, default=1, help="Units evaluated in parallel (not part of the settings)")
    merge_parser = sub.add_parser("merge", help="Merge the shards")
    merge_parser.add_argument("--out-dir", type=Path, required=True)
    decide_parser = sub.add_parser("decide", help=f"Write {DECISION_NAME}")
    decide_parser.add_argument("--out-dir", type=Path, required=True)
    return parser


def settings_from_args(args: argparse.Namespace) -> tuple[dict[str, Any], Any]:
    try:
        export = pc.load_export(args.export_dir)
    except pc.ExportContractError as exc:
        raise pc.contract_exit(exc) from exc
    settings = build_settings(
        export, replicates=args.reps, plant_replicates=args.plant_reps, permutations=args.perms
    )
    return settings, export


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    if args.command == "preflight":
        settings, export = settings_from_args(args)
        result = preflight(args.out_dir, settings, export)
        print(json.dumps({"passed": result["passed"], "units": [
            {key: unit.get(key) for key in ("unit", "passed", "forced_rank", "statistics_differing",
                                            "p_values_differing", "error")}
            for unit in result["units"]
        ]}, indent=1))
        return 0 if result["passed"] else 1
    if args.command == "run":
        index, count = args.shard
        settings, export = settings_from_args(args)
        evaluated = run_shard(args.out_dir, settings, index, count, export=export, limit=args.limit,
                              workers=args.workers)
        print(f"evaluated {evaluated} unit(s)", file=sys.stderr)
        return 0
    if args.command == "merge":
        print(json.dumps(merge(args.out_dir), indent=1))
        return 0
    result = decide(args.out_dir)
    print(json.dumps({"verdict": result["verdict"], "checks": {k: v["grade"] for k, v in result["checks"].items()},
                      "arms": {arm: entry["grade"] for arm, entry in result["arms"].items()},
                      "gate": result["gate"], "complete": result["complete"]}, indent=1))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
