#!/usr/bin/env python3
"""Latent retrofit check: an edit made in a replicate's own fitted PLS basis, pulled back and re-measured.

Descriptive companion of the stage-span pilot (``openspec/changes/archive/2026-11-04-stage-span-constructions/``,
design D12). For every replicate of a pilot profile's zero-effect anchor at one design column it

1. rebuilds the anchor dataset from the profile's matched seeds and fits the
   production PLS latent space once (the evaluation's own preprocessor and rank);
2. edits group B's **latent** stage means with the stage-span edits
   (:func:`~motco.simulations.semisynthetic.edit_stage_configuration`: size,
   orientation, shape and translation, in the latent coordinates);
3. pulls each edit back to the features twice — the **loadings** preimage
   ``Δ Pᵀ`` (``Pᵀ R = I``) and the **minimum-norm** preimage ``Δ R⁺`` — in the
   model's input coordinates, and adds it to group B's samples at each stage
   (methylation in M-values);
4. re-measures the edited dataset through the **frozen** basis (the null fit's
   preprocessor and model: the recovered latent stage means must equal the
   edited ones) and through a full **refit** (the harness from scratch: rank,
   statistics, RRPP p-values).

Beside the geometry of the edit itself (B against edited B), the intended
geometry (A against edited B, which carries the null sample's own A − B
difference) and the frozen and refit geometry it reports the footprint of each
pulled-back edit in pooled-standardized feature units — the effective number of
features (participation ratio), the share on the truth's stage-active features,
the methylation share and the constant-over-stages share — and its cosine with
the stage-span construction's population plant for the same mode and effect on
the same replicate's B configuration. It issues no decision: it documents that
an edit in the basis is recovered, and how close a basis-defined plant is to the
population-defined one the construction uses (the circular alternative design D1
rejects).

Outputs under ``--out-dir``: ``latent_retrofit_check.csv`` (one row per replicate,
mode and preimage; :data:`CSV_COLUMNS`) and ``summary.json`` (settings and, per
mode and preimage, medians, the largest frozen recovery error and the refit
rejection rates).

Pin BLAS to one thread before a ``--workers N > 1`` run
(``OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1``).

Example:

    python scripts/latent_retrofit_check.py \\
        --config examples/trajectory_power_study/stage_span_pilot_surgery_geometry.json \\
        --column '{"generator.n_samples": 200}' --replicates 10 --permutations 199 \\
        --workers 8 --out-dir results/stage-span-pilot-<date>/retrofit
"""

from __future__ import annotations

import argparse
import json
import math
from concurrent.futures import ProcessPoolExecutor
from dataclasses import replace
from pathlib import Path
from typing import Any, Mapping, Sequence

import numpy as np
import pandas as pd

from motco.simulations.diagnostics import geometry_from_means
from motco.simulations.evaluation import (
    SimulationEvaluationParams,
    evaluate_semisynthetic_trajectory,
    integrate_semisynthetic_dataset,
)
from motco.simulations.generator import logit, rev_logit
from motco.simulations.grid import SimulationCell, derive_replicate_seed
from motco.simulations.preprocessing import concatenate_blocks, fit_omics_preprocessor, selected_layers
from motco.simulations.semisynthetic import (
    SemiSyntheticTrajectoryDataset,
    SemiSyntheticTrajectoryParams,
    edit_stage_configuration,
    generate_semisynthetic_trajectory,
)
from motco.simulations.study import enumerate_study, load_study_config

EDIT_MODES: tuple[str, ...] = ("magnitude", "orientation", "shape", "translation")
PREIMAGES: tuple[str, ...] = ("loadings", "minimum_norm")
#: Default edit sizes: a 25 % size change, a 30° rotation, a half-step bend and a half-step offset.
DEFAULT_EFFECTS: dict[str, float] = {"magnitude": 0.25, "orientation": 1.0 / 3.0, "shape": 0.5, "translation": 0.5}
STATISTICS: tuple[str, ...] = ("delta", "angle", "shape")
ALPHA = 0.05
CSV_NAME = "latent_retrofit_check.csv"
SUMMARY_NAME = "summary.json"
CSV_COLUMNS: tuple[str, ...] = (
    "replicate_index",
    "generator_seed",
    "trajectory_mode",
    "effect",
    "preimage",
    "rank",
    "n_features",
    *(f"edit_{s}" for s in STATISTICS),
    *(f"intended_{s}" for s in STATISTICS),
    *(f"frozen_{s}" for s in STATISTICS),
    "frozen_recovery_error",
    "refit_selected_lv",
    *(f"refit_{s}" for s in STATISTICS),
    *(f"refit_p_{s}" for s in STATISTICS),
    "footprint_norm",
    "footprint_effective_features",
    "footprint_stage_active_share",
    "footprint_methylation_share",
    "footprint_constant_share",
    "cosine_population_plant",
)


class RetrofitCheckError(ValueError):
    """Raised when the check's inputs are unusable."""


# --------------------------------------------------------------------------- #
# The column's anchor and its replicates
# --------------------------------------------------------------------------- #


def anchor_cell(config_path: str | Path, column: Mapping[str, Any] | None = None) -> SimulationCell:
    """The profile's zero-effect anchor at the column whose design point holds every ``column`` value.

    ``None`` (or ``{}``) selects the baseline column (the ``power_primary`` anchor).
    """

    cells = [
        cell
        for cell in enumerate_study(load_study_config(config_path)).cells
        if cell.metadata.get("zero_effect_anchor") and cell.phase in ("power_primary", "power_design")
    ]
    wanted = dict(column or {})
    if not wanted:
        matches = [cell for cell in cells if cell.phase == "power_primary"]
    else:
        matches = [
            cell
            for cell in cells
            if all(
                json.dumps((cell.metadata.get("design_point") or {}).get(axis)) == json.dumps(value)
                for axis, value in wanted.items()
            )
        ]
    if len(matches) != 1:
        raise RetrofitCheckError(f"--column {json.dumps(wanted)} selects {len(matches)} anchor cells; name one column.")
    return matches[0]


def replicate_inputs(
    cell: SimulationCell, index: int, permutations: int
) -> tuple[SemiSyntheticTrajectoryParams, SimulationEvaluationParams]:
    """The replicate's generator parameters and the evaluation it runs under (as the study runs it)."""

    seed = derive_replicate_seed(cell, index)
    evaluation = cell.evaluation_params
    if evaluation.integration_method != "pls":
        raise RetrofitCheckError(
            f"The retrofit needs a PLS evaluation; the profile uses {evaluation.integration_method!r}."
        )
    evaluation = replace(
        evaluation,
        permutations=permutations,
        seed=evaluation.seed if evaluation.seed is not None else seed,
        span_capture=None,
    )
    return replace(cell.generator_params, seed=seed), evaluation


# --------------------------------------------------------------------------- #
# Latent edit, preimages and the edited dataset
# --------------------------------------------------------------------------- #


def stage_means(frame: pd.DataFrame, metadata: pd.DataFrame, group: str, stages: Sequence[str]) -> np.ndarray:
    """Group ``group``'s stage means of ``frame`` (rows aligned to ``metadata``), stages in ``stages`` order."""

    groups = metadata["group"].astype(str).to_numpy()
    labels = metadata["stage"].astype(str).to_numpy()
    values = np.asarray(frame, dtype=float)
    return np.vstack([values[(groups == group) & (labels == stage)].mean(axis=0) for stage in stages])


def preimage(model: Any, delta: np.ndarray, kind: str) -> np.ndarray:
    """Feature shifts in the model's input coordinates whose scores shift by ``delta`` (stages × rank)."""

    rotations = np.asarray(model.x_rotations_, dtype=float)
    if kind == "loadings":
        # ``x_rotations_ = W (PᵀW)⁻¹``, so ``Pᵀ R = I`` and ``Δ Pᵀ`` scores to ``Δ``.
        return delta @ np.asarray(model.x_loadings_, dtype=float).T
    if kind == "minimum_norm":
        return delta @ np.linalg.pinv(rotations)
    raise RetrofitCheckError(f"Unknown preimage {kind!r}; expected one of {list(PREIMAGES)}.")


def edited_dataset(
    dataset: SemiSyntheticTrajectoryDataset,
    shifts: np.ndarray,
    *,
    layers: Sequence[str],
    scales: Mapping[str, np.ndarray],
    stages: Sequence[str],
    group: str = "B",
) -> SemiSyntheticTrajectoryDataset:
    """Add per-stage shifts in pooled-standardized units to ``group``'s samples, in analysed units.

    ``shifts`` is stages × joint features in ``layers`` order; methylation is
    shifted in M-values. Population trajectories are dropped: the edit is
    defined on the sample, not on the population.
    """

    metadata = dataset.metadata
    groups = metadata["group"].astype(str).to_numpy()
    labels = metadata["stage"].astype(str).to_numpy()
    blocks: dict[str, pd.DataFrame] = {}
    offset = 0
    for layer in layers:
        frame = getattr(dataset, layer)
        width = frame.shape[1]
        units = shifts[:, offset : offset + width] * np.asarray(scales[layer], dtype=float)
        offset += width
        values = frame.to_numpy(dtype=float).copy()
        analysed = logit(values) if layer == "methylation" else values
        for index, stage in enumerate(stages):
            rows = (groups == group) & (labels == stage)
            analysed[rows] = analysed[rows] + units[index]
        blocks[layer] = pd.DataFrame(
            rev_logit(analysed) if layer == "methylation" else analysed, index=frame.index, columns=frame.columns
        )
    return replace(dataset, population_trajectories=None, **blocks)


# --------------------------------------------------------------------------- #
# Footprints
# --------------------------------------------------------------------------- #


def stage_active_mask(dataset: SemiSyntheticTrajectoryDataset, layers: Sequence[str], group: str = "A") -> np.ndarray:
    """Joint features differential at any stage of ``group``'s truth indicators."""

    indicators = dataset.truth["indicators"][group]
    return np.concatenate([np.asarray(indicators[layer]).sum(axis=1) > 0 for layer in layers])


def footprint(shifts: np.ndarray, *, active: np.ndarray, methylation_width: int) -> dict[str, float]:
    """Where a per-stage plant (stages × features, pooled-standardized) puts its energy."""

    energy = np.sum(shifts**2, axis=0)
    total = float(energy.sum())
    if total <= 0:
        return {
            "footprint_norm": 0.0,
            "footprint_effective_features": math.nan,
            "footprint_stage_active_share": math.nan,
            "footprint_methylation_share": math.nan,
            "footprint_constant_share": math.nan,
        }
    constant = shifts.shape[0] * float(np.sum(shifts.mean(axis=0) ** 2))
    return {
        "footprint_norm": math.sqrt(total),
        "footprint_effective_features": total**2 / float(np.sum(energy**2)),
        "footprint_stage_active_share": float(energy[active].sum()) / total,
        "footprint_methylation_share": float(energy[:methylation_width].sum()) / total,
        "footprint_constant_share": constant / total,
    }


def population_plant(
    generator: SemiSyntheticTrajectoryParams,
    null: SemiSyntheticTrajectoryDataset,
    *,
    mode: str,
    effect: float,
    layers: Sequence[str],
    scales: Mapping[str, np.ndarray],
    group: str = "B",
) -> np.ndarray:
    """The stage-span construction's population plant for ``mode``/``effect``, pooled-standardized.

    Same seed, so the same B configuration as ``null``; stages × joint features.
    """

    stage_span = generate_semisynthetic_trajectory(
        replace(
            generator,
            trajectory_mode=mode,  # type: ignore[arg-type]
            group_effect_size=effect,
            construction="stage_span",
            span_layers=tuple(layers),  # type: ignore[arg-type]
        )
    )
    assert stage_span.population_trajectories is not None and null.population_trajectories is not None
    parts = []
    for layer in layers:
        planted = stage_span.population_trajectories.layers[layer].loc[group].to_numpy(dtype=float)  # type: ignore[index]
        base = null.population_trajectories.layers[layer].loc[group].to_numpy(dtype=float)  # type: ignore[index]
        parts.append((planted - base) / np.asarray(scales[layer], dtype=float))
    return np.hstack(parts)


def cosine(left: np.ndarray, right: np.ndarray) -> float:
    a, b = np.ravel(left), np.ravel(right)
    denominator = float(np.linalg.norm(a) * np.linalg.norm(b))
    return float(a @ b) / denominator if denominator > 0 else math.nan


# --------------------------------------------------------------------------- #
# One replicate
# --------------------------------------------------------------------------- #


def _geometry(means: np.ndarray, stages: Sequence[str], estimator: str) -> dict[str, float | None]:
    measured = geometry_from_means(means, ["A", "B"], list(stages), orientation_estimator=estimator)
    return {"delta": measured.delta, "angle": measured.angle, "shape": measured.shape}


def evaluate_replicate(task: tuple[str, dict[str, Any], int, int, dict[str, float]]) -> list[dict[str, Any]]:
    """Every (mode, preimage) row of one replicate."""

    config_path, column, index, permutations, effects = task
    cell = anchor_cell(config_path, column)
    generator, evaluation = replicate_inputs(cell, index, permutations)
    estimator = evaluation.resolved_orientation_estimator
    layers = selected_layers(evaluation.integration_params)

    null = generate_semisynthetic_trajectory(generator)
    if set(null.metadata["group"].astype(str)) != {"A", "B"}:
        raise RetrofitCheckError("The retrofit expects the generator's default group labels A and B.")
    preprocessor = fit_omics_preprocessor(null, layers)
    observed = preprocessor.transform_dataset(null)
    latent = integrate_semisynthetic_dataset(null, evaluation, preprocessor=preprocessor, observed_blocks=observed)
    assert latent.artifacts is not None
    model = latent.artifacts.model
    x_std = np.asarray(model._x_std, dtype=float)
    scales = {layer: preprocessor.scalers[layer].scale for layer in layers}
    stages = sorted(null.metadata["stage"].astype(str).unique().tolist())
    latent_a = stage_means(latent.matrix, null.metadata, "A", stages)
    latent_b = stage_means(latent.matrix, null.metadata, "B", stages)
    active = stage_active_mask(null, layers)
    methylation_width = null.methylation.shape[1] if "methylation" in layers else 0
    rank = int(latent.matrix.shape[1])

    rows: list[dict[str, Any]] = []
    for mode in EDIT_MODES:
        effect = float(effects[mode])
        edited_b, _ = edit_stage_configuration(latent_b, mode, effect)
        delta = edited_b - latent_b
        edit = _geometry(np.vstack([latent_b, edited_b]), stages, estimator)
        intended = _geometry(np.vstack([latent_a, edited_b]), stages, estimator)
        plant = population_plant(generator, null, mode=mode, effect=effect, layers=layers, scales=scales)
        for kind in PREIMAGES:
            shifts = preimage(model, delta, kind) * x_std  # input → pooled-standardized coordinates
            dataset = edited_dataset(null, shifts, layers=layers, scales=scales, stages=stages)
            frozen_scores = np.asarray(
                model.transform(concatenate_blocks(preprocessor.transform_dataset(dataset)).to_numpy(dtype=float))
            )
            frozen_b = stage_means(frozen_scores, dataset.metadata, "B", stages)
            frozen_a = stage_means(frozen_scores, dataset.metadata, "A", stages)
            frozen = _geometry(np.vstack([frozen_a, frozen_b]), stages, estimator)
            refit = evaluate_semisynthetic_trajectory(dataset, evaluation)
            row: dict[str, Any] = {
                "replicate_index": index,
                "generator_seed": generator.seed,
                "trajectory_mode": mode,
                "effect": effect,
                "preimage": kind,
                "rank": rank,
                "n_features": int(x_std.size),
                **{f"edit_{s}": edit[s] for s in STATISTICS},
                **{f"intended_{s}": intended[s] for s in STATISTICS},
                **{f"frozen_{s}": frozen[s] for s in STATISTICS},
                "frozen_recovery_error": float(np.max(np.abs(frozen_b - edited_b))),
                "refit_selected_lv": refit.latent_matrix_metadata.get("selected_lv"),
                **{f"refit_{s}": refit.pair_statistics.get(s) for s in STATISTICS},
                **{f"refit_p_{s}": refit.p_values.get(s) for s in STATISTICS},
                **footprint(shifts, active=active, methylation_width=methylation_width),
                "cosine_population_plant": cosine(shifts, plant),
            }
            rows.append(row)
    return rows


# --------------------------------------------------------------------------- #
# Summary and outputs
# --------------------------------------------------------------------------- #


def _median(values: Sequence[Any]) -> float | None:
    finite = [float(v) for v in values if v is not None and math.isfinite(float(v))]
    return float(np.median(finite)) if finite else None


def summarize(rows: Sequence[Mapping[str, Any]], settings: Mapping[str, Any]) -> dict[str, Any]:
    """Settings and, per mode and preimage, medians, the largest frozen error and refit rejection rates."""

    frame = pd.DataFrame(list(rows), columns=list(CSV_COLUMNS))
    cells: list[dict[str, Any]] = []
    for (mode, kind), group in frame.groupby(["trajectory_mode", "preimage"], sort=False):
        entry: dict[str, Any] = {
            "trajectory_mode": mode,
            "preimage": kind,
            "effect": float(group["effect"].iloc[0]),
            "n_replicates": int(len(group)),
            "max_frozen_recovery_error": float(group["frozen_recovery_error"].max()),
            "median_rank": _median(group["rank"].tolist()),
            "median_refit_selected_lv": _median(group["refit_selected_lv"].tolist()),
        }
        for column in (
            *(f"edit_{s}" for s in STATISTICS),
            *(f"intended_{s}" for s in STATISTICS),
            *(f"frozen_{s}" for s in STATISTICS),
            *(f"refit_{s}" for s in STATISTICS),
            "footprint_effective_features",
            "footprint_stage_active_share",
            "footprint_methylation_share",
            "footprint_constant_share",
            "cosine_population_plant",
        ):
            entry[f"median_{column}"] = _median(group[column].tolist())
        for statistic in STATISTICS:
            p = [float(v) for v in group[f"refit_p_{statistic}"] if v is not None and math.isfinite(float(v))]
            entry[f"refit_rejection_rate_{statistic}"] = (sum(v < ALPHA for v in p) / len(p)) if p else None
            entry[f"n_refit_p_{statistic}"] = len(p)
        cells.append(entry)
    return {
        "settings": dict(settings),
        "alpha": ALPHA,
        "n_rows": int(len(frame)),
        "max_frozen_recovery_error": float(frame["frozen_recovery_error"].max()) if len(frame) else None,
        "cells": cells,
        "units": "latent: the replicate's fitted PLS scores; footprints and cosines: pooled-standardized features",
        "decision": None,
    }


def write_outputs(rows: Sequence[Mapping[str, Any]], settings: Mapping[str, Any], out_dir: Path) -> dict[str, Path]:
    out_dir.mkdir(parents=True, exist_ok=True)
    paths = {"csv": out_dir / CSV_NAME, "summary": out_dir / SUMMARY_NAME}
    frame = pd.DataFrame(list(rows), columns=list(CSV_COLUMNS))
    frame.sort_values(["replicate_index", "trajectory_mode", "preimage"]).to_csv(paths["csv"], index=False)
    summary = summarize(rows, settings)
    paths["summary"].write_text(json.dumps(summary, indent=1, sort_keys=True) + "\n", encoding="utf-8")
    return paths


def run_check(
    config: str | Path,
    *,
    column: Mapping[str, Any] | None,
    replicates: Sequence[int],
    permutations: int,
    effects: Mapping[str, float],
    workers: int = 1,
) -> list[dict[str, Any]]:
    anchor_cell(config, column)  # fail on a bad column before any work
    tasks = [(str(config), dict(column or {}), int(index), permutations, dict(effects)) for index in replicates]
    if workers > 1:
        with ProcessPoolExecutor(max_workers=workers) as pool:
            batches = list(pool.map(evaluate_replicate, tasks))
    else:
        batches = [evaluate_replicate(task) for task in tasks]
    return [row for batch in batches for row in batch]


def parse_effects(clauses: Sequence[str] | None) -> dict[str, float]:
    effects = dict(DEFAULT_EFFECTS)
    for clause in clauses or []:
        mode, sep, raw = clause.partition("=")
        if not sep or mode not in EDIT_MODES:
            raise RetrofitCheckError(f"--effect must read MODE=VALUE with MODE in {list(EDIT_MODES)}; got {clause!r}.")
        effects[mode] = float(raw)
    return effects


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--config", type=Path, required=True, help="A stage-span pilot profile.")
    parser.add_argument(
        "--column", default=None, help="JSON design-point subset selecting the column (default: the baseline)."
    )
    parser.add_argument("--replicates", type=int, default=10, help="Replicate indices 0 … N − 1 (default 10).")
    parser.add_argument("--permutations", type=int, default=199, help="RRPP permutations of each refit (default 199).")
    parser.add_argument(
        "--effect", action="append", default=None, help="MODE=VALUE edit size (repeatable; defaults in the docstring)."
    )
    parser.add_argument("--workers", type=int, default=1, help="Worker processes, one replicate each (default 1).")
    parser.add_argument("--out-dir", type=Path, required=True)
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    column = json.loads(args.column) if args.column else None
    effects = parse_effects(args.effect)
    rows = run_check(
        args.config,
        column=column,
        replicates=range(args.replicates),
        permutations=args.permutations,
        effects=effects,
        workers=args.workers,
    )
    settings = {
        "config": str(args.config),
        "column": column or {},
        "replicates": args.replicates,
        "permutations": args.permutations,
        "effects": effects,
        "preimages": list(PREIMAGES),
    }
    paths = write_outputs(rows, settings, args.out_dir)
    summary = json.loads(paths["summary"].read_text(encoding="utf-8"))
    print(f"{summary['n_rows']} rows; largest frozen recovery error {summary['max_frozen_recovery_error']:.3g}")
    for cell in summary["cells"]:
        rates = ", ".join(f"{s} {cell[f'refit_rejection_rate_{s}']}" for s in STATISTICS)
        print(
            f"  {cell['trajectory_mode']:<11} {cell['preimage']:<12} refit rejection [{rates}]; "
            f"cosine with the population plant {cell['median_cosine_population_plant']}"
        )
    print(f"Wrote {paths['csv']} and {paths['summary']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
