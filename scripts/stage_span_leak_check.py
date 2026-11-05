#!/usr/bin/env python3
"""Stage-span leak decomposition (``stage-span-paper-grade``, design D7 and D8): descriptive, no verdict.

Runs a profile's stage-span cells unit by unit and reports, per unit, the
signed B − A path-length difference (and both path lengths) at six checkpoints,
so each route by which a stage-span plant leaks into ``delta`` can be read off
one table:

- ``pop_null``    population stage means / the null pooled SD (the
  construction's units under the null scale; 0 for the bend there)
- ``pop_planted`` population stage means / the exact planted pooled SD
  (``planted_pooled_scale``, the sample scaler's population analogue; 0 for the
  bend, the rotation and the translation under ``span_scale: "planted"``)
- ``pop_std``     population means through the replicate's sample scaler
- ``pop_span``    those rows projected on an orthonormal basis of the fitted span
- ``pop_lat``     those rows in the fitted PLS score coordinates
  (``model.transform``, the ``population_latent`` checkpoint)
- ``obs_lat``     the sample LS means in the latent space (what the test measures)

plus the selected rank, the capture share and, with permutations, the three
p-values. Each unit rebuilds the cell's dataset from the profile's matched
seeds (seeds depend only on the family and the replicate index, so an
overridden effect or scale reads the same baseline draw) and fits the
production latent space once.

Subcommands:

- ``run``   — one shard of the units: ``--config``, ``--shard``, ``--n-shards``,
  ``--span-scale null|planted`` (default: the profile's), ``--columns``
  (default: every column), ``--modes`` (default:
  the profile's), ``--effects`` (default: the profile's non-zero effects),
  ``--replicates``, ``--permutations``, ``--out-dir``. Every column's anchor is
  included. Writes ``shard_<i>.csv`` (resumable: units already complete are
  skipped) and ``settings.json`` (a resume under other settings is refused).
- ``merge`` — ``leak_units.csv``, ``leak_cells.csv`` (per column, mode and
  effect: median of every difference, the differences as a share of A's path
  length, median capture and rank, rejection rates with p < α) and
  ``summary.json`` (settings, units, failures, the config's sha256).
- ``gate``  — reads a planted arm and a null arm and writes
  ``preflight_gate.json`` (design D8): G1 the planted arm's units all
  generated and their planted-scale difference is at most 1e-8 of A's path;
  G2 the planted orientation and shape off-target rates pass the decision's
  bound rule; G3 at every column the planted arm's median score-coordinate
  difference for ``shape`` at e = 0.5 is smaller in magnitude than the null
  arm's. ``gate_holds`` | ``gate_fails``, naming every failing condition and cell.

Example:

    python scripts/stage_span_leak_check.py run \\
        --config examples/trajectory_power_study/stage_span_paper_grade_span_geometry.json \\
        --span-scale planted --modes orientation shape --effects 0.1 0.25 0.5 \\
        --replicates 100 --permutations 199 --shard 0 --n-shards 100 --out-dir results/leak/planted
    python scripts/stage_span_leak_check.py merge --out-dir results/leak/planted
    python scripts/stage_span_leak_check.py gate --planted results/leak/planted --null results/leak/null \\
        --out results/leak/preflight_gate.json
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import sys
import time
from collections.abc import Mapping, Sequence
from dataclasses import replace
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

from motco.simulations.evaluation import evaluate_semisynthetic_trajectory, integrate_semisynthetic_dataset
from motco.simulations.grid import SimulationCell, derive_replicate_seed
from motco.simulations.preprocessing import concatenate_blocks, fit_omics_preprocessor, selected_layers
from motco.simulations.reference import load_reference
from motco.simulations.semisynthetic import (
    generate_semisynthetic_trajectory,
    null_pooled_scale,
    planted_pooled_scale,
)
from motco.simulations.study import enumerate_study, load_study_config

sys.path.insert(0, str(Path(__file__).resolve().parent))
try:
    from stage_span_paper_grade_decision import judge_rates
finally:
    sys.path.pop(0)

ALPHA = 0.05
STATISTICS: tuple[str, ...] = ("delta", "angle", "shape")
TARGETS: dict[str, str] = {"magnitude": "delta", "orientation": "angle", "shape": "shape"}
CHECKPOINTS: tuple[str, ...] = ("pop_null", "pop_planted", "pop_std", "pop_span", "pop_lat", "obs_lat")
EXACT_TOLERANCE = 1e-8
GATE_SHAPE_EFFECT = 0.5
SPAN_SCALES = ("null", "planted")
_DESIGN_PHASES = ("power_primary", "power_design")
SEAAD_TABLE = ((11, 10, 28), (9, 10, 12))

#: Unit key: (column, mode, effect, replicate).
UnitKey = tuple[str, str, float, int]


class LeakCheckError(ValueError):
    """Raised when the check's inputs are inconsistent."""


# --------------------------------------------------------------------------- #
# Units
# --------------------------------------------------------------------------- #


def column_of(cell: SimulationCell) -> str:
    """The canonical column name of a cell: ``n200``, ``n80``, ``seaad``, ``…_three_blocks``, ``…_four_stages``."""

    generator, evaluation = cell.generator_params, cell.evaluation_params
    if generator.group_stage_sizes is not None:
        table = tuple(tuple(int(count) for count in group) for group in generator.group_stage_sizes)
        name = "seaad" if table == SEAAD_TABLE else f"table{table}"
    else:
        name = f"n{generator.n_samples}"
    if (evaluation.integration_params or {}).get("layers") is None:
        name += "_three_blocks"
    if generator.n_stages == 4:
        name += "_four_stages"
    return name


def column_cells(config_path: Path) -> dict[str, dict[str, SimulationCell]]:
    """Per column, one primary-family cell per mode (``none`` = the anchor)."""

    columns: dict[str, dict[str, SimulationCell]] = {}
    for cell in enumerate_study(load_study_config(config_path)).cells:
        if cell.phase not in _DESIGN_PHASES:
            continue
        effect = float(cell.generator_params.group_effect_size)
        mode = "none" if effect == 0 else cell.generator_params.trajectory_mode
        columns.setdefault(column_of(cell), {}).setdefault(mode, cell)
    return columns


def unit_keys(
    columns: Mapping[str, Mapping[str, SimulationCell]],
    modes: Sequence[str],
    effects: Sequence[float],
    replicates: int,
) -> list[UnitKey]:
    keys: list[UnitKey] = []
    for column in sorted(columns):
        cells = [("none", 0.0)] + [(mode, float(effect)) for mode in modes for effect in effects]
        for mode, effect in cells:
            if mode not in columns[column]:
                raise LeakCheckError(f"Column {column!r} has no {mode!r} cell in the profile.")
            keys.extend((column, mode, effect, index) for index in range(replicates))
    return keys


def unit_params(
    columns: Mapping[str, Mapping[str, SimulationCell]],
    key: UnitKey,
    *,
    span_scale: str | None,
    permutations: int,
) -> tuple[Any, Any]:
    column, mode, effect, index = key
    cell = columns[column][mode]
    generator = replace(cell.generator_params, seed=derive_replicate_seed(cell, index))
    if mode != "none":
        generator = replace(generator, group_effect_size=effect)
        if span_scale is not None:
            generator = replace(generator, span_scale=None if span_scale == "null" else span_scale)
    evaluation = replace(cell.evaluation_params, permutations=permutations, span_capture=True)
    return generator, evaluation


def _stage_order(stage: Any) -> tuple[int, Any]:
    text = str(stage)
    return (0, int(text)) if text.isdigit() else (1, text)


def _path(rows: np.ndarray) -> float:
    return float(np.linalg.norm(np.diff(rows, axis=0), axis=1).sum())


def run_unit(generator: Any, evaluation: Any, reference: Any) -> dict[str, Any]:
    """The six checkpoints, rank, capture and p-values of one unit."""

    dataset = generate_semisynthetic_trajectory(generator, reference=reference)
    layers = selected_layers(evaluation.integration_params)
    preprocessor = fit_omics_preprocessor(dataset, layers)
    latent = integrate_semisynthetic_dataset(
        dataset, evaluation, preprocessor=preprocessor, observed_blocks=preprocessor.transform_dataset(dataset)
    )
    result = evaluate_semisynthetic_trajectory(dataset, evaluation, latent=latent)
    latent_map = latent.linear_map
    if latent_map is None:
        raise LeakCheckError("The latent space carries no linear map; the leak check needs a PLS measurement.")

    label_a, label_b = generator.group_labels
    population = dataset.population_trajectories
    frame = concatenate_blocks({layer: population.layers[layer] for layer in layers})
    stages = sorted({stage for _, stage in frame.index}, key=_stage_order)
    means = np.vstack(
        [np.asarray(frame.loc[(group, stage)], dtype=float) for group in (label_a, label_b) for stage in stages]
    )
    k = len(stages)
    counts = dataset.metadata.assign(stage=dataset.metadata["stage"].astype(str)).groupby(["group", "stage"]).size()
    sizes_a = np.array([counts.loc[(label_a, str(stage))] for stage in stages], dtype=float)
    sizes_b = np.array([counts.loc[(label_b, str(stage))] for stage in stages], dtype=float)
    widths = [population.layers[layer].shape[1] for layer in layers]
    bounds = np.cumsum([0, *widths])
    means_a = {layer: means[:k, bounds[i] : bounds[i + 1]] for i, layer in enumerate(layers)}
    means_b = {layer: means[k:, bounds[i] : bounds[i + 1]] for i, layer in enumerate(layers)}
    s_null = null_pooled_scale(means_a, sizes_a + sizes_b, layers=layers, reference=reference)
    s_planted = planted_pooled_scale(means_a, means_b, sizes_a, sizes_b, layers=layers, reference=reference)

    row: dict[str, Any] = {"rank": int(latent.matrix.shape[1]), "n_features": int(latent_map.basis.shape[0])}

    def put(name: str, rows: np.ndarray) -> None:
        a, b = _path(rows[:k]), _path(rows[k:])
        row[f"{name}_A"], row[f"{name}_B"], row[f"{name}_d"] = a, b, b - a

    put("pop_null", means / s_null)
    put("pop_planted", means / s_planted)
    standardized = concatenate_blocks(preprocessor.transform_population(population))
    x = np.vstack(
        [np.asarray(standardized.loc[(group, stage)], dtype=float) for group in (label_a, label_b) for stage in stages]
    )
    put("pop_std", x)
    put("pop_span", latent_map.inputs(x) @ latent_map.basis)
    put("pop_lat", np.asarray(latent_map.transform(x), dtype=float))
    observed = result.realized_geometry["checkpoints"]["pls_latent"]["joint"]["path_lengths"]
    row["obs_lat_A"], row["obs_lat_B"] = observed[label_a], observed[label_b]
    row["obs_lat_d"] = observed[label_b] - observed[label_a]
    row["capture"] = (result.realized_geometry.get("span_capture") or {}).get("group_difference")
    for statistic in STATISTICS:
        row[f"p_{statistic}"] = (result.p_values or {}).get(statistic)
    stage_span = dataset.truth.get("stage_span") or {}
    row["scale_iterations"] = stage_span.get("scale_iterations")
    row["scale_residual"] = stage_span.get("scale_residual")
    return row


# --------------------------------------------------------------------------- #
# run
# --------------------------------------------------------------------------- #


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _settings(args: argparse.Namespace, modes: Sequence[str], effects: Sequence[float]) -> dict[str, Any]:
    return {
        "config": str(args.config),
        "config_sha256": _sha256(args.config),
        "span_scale": args.span_scale,
        "modes": list(modes),
        "effects": [float(effect) for effect in effects],
        "replicates": int(args.replicates),
        "permutations": int(args.permutations),
        "columns": sorted(args.columns) if args.columns else None,
        "n_shards": int(args.n_shards),
        "alpha": ALPHA,
    }


def _bind_settings(out_dir: Path, settings: Mapping[str, Any]) -> None:
    path = out_dir / "settings.json"
    if path.exists():
        existing = json.loads(path.read_text())
        differing = sorted(key for key in set(existing) | set(settings) if existing.get(key) != settings.get(key))
        if differing:
            raise LeakCheckError(f"{path} holds other settings ({differing}); use a new --out-dir.")
        return
    tmp = path.with_suffix(f".{time.time_ns()}.tmp")
    tmp.write_text(json.dumps(dict(settings), indent=1, sort_keys=True) + "\n")
    tmp.replace(path)


def command_run(args: argparse.Namespace) -> int:
    config = load_study_config(args.config)
    modes = list(args.modes or [mode for mode in config.trajectory_modes])
    effects = [float(effect) for effect in (args.effects or [e for e in config.effect_sizes if float(e) != 0.0])]
    args.out_dir.mkdir(parents=True, exist_ok=True)
    _bind_settings(args.out_dir, _settings(args, modes, effects))
    columns = column_cells(args.config)
    if args.columns:
        unknown = sorted(set(args.columns) - set(columns))
        if unknown:
            raise LeakCheckError(f"--columns names unknown column(s) {unknown}; the profile has {sorted(columns)}.")
        columns = {name: cells for name, cells in columns.items() if name in args.columns}
    keys = unit_keys(columns, modes, effects, args.replicates)
    mine = keys[args.shard :: args.n_shards]
    out = args.out_dir / f"shard_{args.shard:03d}.csv"
    rows = pd.read_csv(out).to_dict("records") if out.exists() else []
    done = {
        (row["column"], row["mode"], float(row["effect"]), int(row["replicate"]))
        for row in rows
        if not isinstance(row.get("error"), str)
    }
    rows = [row for row in rows if (row["column"], row["mode"], float(row["effect"]), int(row["replicate"])) in done]
    reference = load_reference()
    span_scale = None if args.span_scale is None else args.span_scale
    for key in mine:
        if key in done:
            continue
        generator, evaluation = unit_params(columns, key, span_scale=span_scale, permutations=args.permutations)
        column, mode, effect, index = key
        row: dict[str, Any] = {
            "column": column,
            "mode": mode,
            "effect": effect,
            "replicate": index,
            "seed": generator.seed,
            "span_scale": generator.span_scale or "null_scale",  # "null" would read back as NaN
        }
        started = time.perf_counter()
        try:
            row.update(run_unit(generator, evaluation, reference))
        except Exception as error:  # recorded per unit, as the study's error policy does
            row["error"] = f"{type(error).__name__}: {error}"
        row["seconds"] = time.perf_counter() - started
        rows.append(row)
        pd.DataFrame(rows).to_csv(out, index=False)
    print(f"shard {args.shard}: {len(mine)} of {len(keys)} units", flush=True)
    return 0


# --------------------------------------------------------------------------- #
# merge
# --------------------------------------------------------------------------- #


def _median(values: pd.Series) -> float | None:
    finite = pd.to_numeric(values, errors="coerce").dropna()
    return float(finite.median()) if len(finite) else None


def summarize_cells(units: pd.DataFrame, *, alpha: float = ALPHA) -> pd.DataFrame:
    """Per (column, mode, effect): medians, shares of A's path and rejection rates (p < α)."""

    rows = []
    failed_mask = units["error"].notna() if "error" in units else pd.Series(False, index=units.index)
    for (column, mode, effect), group in units.groupby(["column", "mode", "effect"], sort=True):
        ok = group[~failed_mask.loc[group.index]]
        row: dict[str, Any] = {
            "column": column,
            "mode": mode,
            "effect": float(effect),
            "n_units": int(len(group)),
            "n_failed": int(len(group) - len(ok)),
        }
        for name in CHECKPOINTS:
            row[f"median_{name}_d"] = _median(ok[f"{name}_d"]) if len(ok) else None
            row[f"median_{name}_share"] = _median(ok[f"{name}_d"] / ok[f"{name}_A"]) if len(ok) else None
        row["median_capture"] = _median(ok["capture"]) if len(ok) and "capture" in ok else None
        row["median_rank"] = _median(ok["rank"]) if len(ok) else None
        for statistic in STATISTICS:
            column_name = f"p_{statistic}"
            values = pd.to_numeric(ok[column_name], errors="coerce").dropna() if column_name in ok else pd.Series()
            row[f"n_available_{statistic}"] = int(len(values))
            row[f"n_rejected_{statistic}"] = int((values < alpha).sum())
            row[f"rate_{statistic}"] = float((values < alpha).mean()) if len(values) else None
        rows.append(row)
    return pd.DataFrame(rows)


def command_merge(args: argparse.Namespace) -> int:
    out_dir: Path = args.out_dir
    settings_path = out_dir / "settings.json"
    if not settings_path.exists():
        raise LeakCheckError(f"Missing {settings_path}; merge reads a run directory.")
    settings = json.loads(settings_path.read_text())
    shards = sorted(out_dir.glob("shard_*.csv"))
    if not shards:
        raise LeakCheckError(f"No shard_*.csv under {out_dir}.")
    units = pd.concat([pd.read_csv(path) for path in shards], ignore_index=True)
    if "error" not in units:
        units["error"] = np.nan
    keys = units[["column", "mode", "effect", "replicate"]].apply(tuple, axis=1)
    duplicates = int(keys.duplicated().sum())
    units = units.sort_values(["column", "mode", "effect", "replicate"]).reset_index(drop=True)
    cells = summarize_cells(units, alpha=float(settings.get("alpha", ALPHA)))
    units.to_csv(out_dir / "leak_units.csv", index=False)
    cells.to_csv(out_dir / "leak_cells.csv", index=False)
    columns = column_cells(Path(settings["config"]))
    if settings.get("columns"):
        columns = {name: cells for name, cells in columns.items() if name in settings["columns"]}
    expected = len(unit_keys(columns, settings["modes"], settings["effects"], settings["replicates"]))
    failures = units[units["error"].notna()]
    summary = {
        "settings": settings,
        "config_sha256_now": _sha256(Path(settings["config"])),
        "n_shards_found": len(shards),
        "n_units": int(len(units)),
        "n_expected": expected,
        "n_duplicates": duplicates,
        "n_failed": int(len(failures)),
        "failures": [
            {key: row[key] for key in ("column", "mode", "effect", "replicate", "error")}
            for row in failures.to_dict("records")
        ],
        "median_seconds": _median(units["seconds"]) if "seconds" in units else None,
        "core_hours": float(pd.to_numeric(units["seconds"], errors="coerce").sum() / 3600)
        if "seconds" in units
        else None,
    }
    (out_dir / "summary.json").write_text(json.dumps(summary, indent=1, sort_keys=True, default=str) + "\n")
    print(
        f"merged {len(units)} units ({expected} expected, {duplicates} duplicates, {len(failures)} failed) "
        f"-> {out_dir / 'leak_cells.csv'}"
    )
    return 0


# --------------------------------------------------------------------------- #
# gate
# --------------------------------------------------------------------------- #


def gate(planted_units: pd.DataFrame, planted_cells: pd.DataFrame, null_cells: pd.DataFrame) -> dict[str, Any]:
    """The pre-flight conditions G1–G3 (design D8) on the two arms' merged tables."""

    failures: list[dict[str, Any]] = []
    errors = planted_units["error"].notna() if "error" in planted_units else pd.Series(False, index=planted_units.index)
    g1_rows = []
    for row in planted_units.to_dict("records"):
        key = {k: row[k] for k in ("column", "mode", "effect", "replicate")}
        if isinstance(row.get("error"), str):
            g1_rows.append({**key, "status": "failed", "error": row["error"]})
            continue
        gap = abs(float(row["pop_planted_d"])) / float(row["pop_planted_A"])
        if row["mode"] == "magnitude":
            continue
        if gap > EXACT_TOLERANCE:
            g1_rows.append({**key, "status": "inexact", "relative_gap": gap})
    g1 = {
        "holds": not g1_rows,
        "n_units": int(len(planted_units)),
        "n_failed": int(errors.sum()),
        "tolerance": EXACT_TOLERANCE,
        "violations": g1_rows,
    }
    failures.extend({"condition": "G1", **row} for row in g1_rows)

    observations = []
    for row in planted_cells.to_dict("records"):
        if row["mode"] not in ("orientation", "shape"):
            continue
        for statistic in STATISTICS:
            if statistic == TARGETS[row["mode"]]:
                continue
            cell = {"column": row["column"], "trajectory_mode": row["mode"], "effect_size": row["effect"]}
            n_available = int(row[f"n_available_{statistic}"])
            if n_available <= 0:
                observations.append({**cell, "statistic": statistic, "status": "unavailable"})
            else:
                observations.append(
                    {
                        **cell,
                        "statistic": statistic,
                        "n_available": n_available,
                        "n_rejected": int(row[f"n_rejected_{statistic}"]),
                    }
                )
    g2 = judge_rates(observations)
    failures.extend({"condition": "G2", **row} for row in g2["observations"] if row["status"] != "holds")

    g3_rows = []
    for column in sorted(set(planted_cells["column"])):

        def pick(cells: pd.DataFrame) -> float | None:
            match = cells[
                (cells["column"] == column)
                & (cells["mode"] == "shape")
                & np.isclose(cells["effect"].astype(float), GATE_SHAPE_EFFECT)
            ]
            value = match["median_pop_lat_d"].iloc[0] if len(match) else None
            return None if value is None or (isinstance(value, float) and math.isnan(value)) else float(value)

        planted, null = pick(planted_cells), pick(null_cells)
        if planted is None or null is None:
            status = "unavailable"
        else:
            status = "shrinks" if abs(planted) < abs(null) else "does_not_shrink"
        g3_rows.append(
            {"column": column, "planted_median_pop_lat_d": planted, "null_median_pop_lat_d": null, "status": status}
        )
    g3 = {"holds": all(row["status"] == "shrinks" for row in g3_rows), "effect": GATE_SHAPE_EFFECT, "columns": g3_rows}
    failures.extend({"condition": "G3", **row} for row in g3_rows if row["status"] != "shrinks")

    conditions = {"G1": g1, "G2": g2, "G3": g3}
    holds = all(condition["holds"] for condition in conditions.values())
    return {
        "verdict": "gate_holds" if holds else "gate_fails",
        "failing_conditions": [name for name, condition in conditions.items() if not condition["holds"]],
        "conditions": conditions,
        "failures": failures,
    }


def command_gate(args: argparse.Namespace) -> int:
    for arm in (args.planted, args.null):
        for name in ("leak_units.csv", "leak_cells.csv"):
            if not (arm / name).exists():
                raise LeakCheckError(f"Missing {arm / name}; merge each arm first.")
    result = gate(
        pd.read_csv(args.planted / "leak_units.csv"),
        pd.read_csv(args.planted / "leak_cells.csv"),
        pd.read_csv(args.null / "leak_cells.csv"),
    )
    result["inputs"] = {
        str(path): _sha256(path)
        for arm in (args.planted, args.null)
        for path in (arm / "leak_units.csv", arm / "leak_cells.csv", arm / "settings.json")
        if path.exists()
    }
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(result, indent=1, sort_keys=True, default=str) + "\n")
    print(f"pre-flight gate: {result['verdict']} (failing: {result['failing_conditions']})")
    return 0


# --------------------------------------------------------------------------- #
# CLI
# --------------------------------------------------------------------------- #


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = parser.add_subparsers(dest="command", required=True)
    run = sub.add_parser("run", help="Run one shard of the units.")
    run.add_argument("--config", type=Path, required=True)
    run.add_argument("--shard", type=int, required=True)
    run.add_argument("--n-shards", type=int, required=True)
    run.add_argument("--span-scale", choices=SPAN_SCALES, default=None, help="Override the profile's span scale.")
    run.add_argument("--columns", nargs="+", default=None, help="Restrict to these columns (default: all).")
    run.add_argument("--modes", nargs="+", default=None)
    run.add_argument("--effects", nargs="+", type=float, default=None)
    run.add_argument("--replicates", type=int, required=True)
    run.add_argument("--permutations", type=int, required=True)
    run.add_argument("--out-dir", type=Path, required=True)
    merge = sub.add_parser("merge", help="Merge the shards into leak_units.csv, leak_cells.csv and summary.json.")
    merge.add_argument("--out-dir", type=Path, required=True)
    gate_parser = sub.add_parser("gate", help="Read the planted and null arms into preflight_gate.json.")
    gate_parser.add_argument("--planted", type=Path, required=True)
    gate_parser.add_argument("--null", type=Path, required=True)
    gate_parser.add_argument("--out", type=Path, required=True)
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    try:
        if args.command == "run":
            if not 0 <= args.shard < args.n_shards:
                raise LeakCheckError(f"--shard must be in [0, {args.n_shards - 1}]; got {args.shard}.")
            return command_run(args)
        if args.command == "merge":
            return command_merge(args)
        return command_gate(args)
    except LeakCheckError as exc:
        raise SystemExit(str(exc)) from exc


if __name__ == "__main__":
    raise SystemExit(main())
