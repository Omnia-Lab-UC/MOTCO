#!/usr/bin/env python3
"""Covariate-adjustment Type I check at the Phase 6 SEA-AD design (no R).

The gate behind ``openspec/changes/covariate-adjusted-trajectories/`` (design
D10, spec ``covariate-adjustment-type-i``). The case study may report
covariate-adjusted ``motco analyze`` results only after this check returns
``covariates_validated``.

Each **replicate** generates a null semi-synthetic dataset at the Phase 6
small-n design (group × stage sizes ``[[11, 10, 28], [9, 10, 12]]``, three
stages, methylation + expression, joint magnitude, ``group_effect_size = 0``),
then:

- draws a numeric covariate ``z = a_g·g + a_s·s + a_gs·g·s + ε`` (``g`` the
  group indicator, ``s`` the stage index, ``ε ~ N(0, 1)``) and a three-level
  categorical covariate from a multinomial logit with level logits
  ``k·c·η`` (``k = 0, 1, 2``, ``c`` the recorded logit scale, which keeps
  the factor mixed within cells) on the same linear predictor ``η``. The
  ``a_gs`` interaction makes the covariate means differ between the groups'
  *trajectories*, not only their locations;
- injects the covariate effect into every measured feature,
  ``X_f += s_z·sd_f·w_f·z + s_c·sd_f·Σ_k u_{k,f}·1[cat = k]``, where ``sd_f``
  is the feature's pooled SD and ``w``/``u_k`` are unit directions over a
  random fraction of the joint features drawn once from the base seed (so the
  confounding structure is fixed across replicates). Methylation is converted
  to M-values first and analysed with the ``identity`` transform, so the
  injection is additive on the analysed scale;
- analyses the dataset with :func:`motco.analysis.analyze_trajectories` twice,
  with both covariates (``adjusted``) and without (``unadjusted``), with the
  same CV settings and RRPP seed. The two latent spaces must be identical; the
  replicate fails loudly otherwise.

**Calibration (2026-10-16).** The design's first strengths (``a_s = 0.5``,
``a_gs = 0.75``, ``s_z = 0.5``, ``s_c = 0.3``) did not confound: the unadjusted
arm rejected at about α in the R = 20 smoke run. Raising ``s_z`` alone (to 30)
does not help, because the injection carries ``z``'s within-cell noise along
``w`` in the same proportion as its cell-mean pattern, so RRPP's null widens
with the signal. The confounding signal-to-noise is ``a_gs`` relative to the
covariate noise, so the defaults are ``a_s = 1``, ``a_gs = 2``, ``s_z = 3``,
``s_c = 2`` and a categorical logit scale of 0.5 (which keeps the factor mixed
within cells, away from the full-rank guard). Every value is recorded.

The **decision** (pure, :func:`covariate_type_i_decision`): per statistic, the
adjusted arm's rejection rate (``p ≤ alpha``) over its R available replicates
must be at or below ``alpha + k·√(alpha(1 − alpha)/R)``; every statistic
holding gives ``covariates_validated``, otherwise ``diagnose_covariate_excess``.
The unadjusted arm is reported only, as evidence the construction confounds.

Outputs (``--out-dir``):

- ``settings.json``: settings and their hash; a resume into a directory
  written with different settings is refused;
- ``records/shard_<i>_of_<N>.jsonl``: one line per replicate (both arms);
  resume skips recorded replicates;
- after ``merge`` (automatic for an unsharded run): ``merged.jsonl``,
  ``covariate_type_i_pvalues.csv``, ``covariate_type_i_rates.csv`` and
  ``covariate_type_i_decision.json``.

Examples:

    # smoke profile
    python scripts/covariate_type_i_check.py run --reps 20 --perms 99 --out-dir /tmp/cov-smoke

    # one SLURM array task, then the merge
    python scripts/covariate_type_i_check.py run --shard 3/50 --out-dir results/cov
    python scripts/covariate_type_i_check.py merge --out-dir results/cov
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import sys
from pathlib import Path
from time import perf_counter
from typing import Any, Iterable, Mapping, Sequence

import numpy as np
import pandas as pd

SCRIPT_VERSION = 1
STATISTICS: tuple[str, ...] = ("delta", "angle", "shape")
ARMS: tuple[str, ...] = ("adjusted", "unadjusted")
CATEGORY_LEVELS: tuple[str, ...] = ("c0", "c1", "c2")
LAYERS: tuple[str, ...] = ("methylation", "expression")

#: Generator settings of ``examples/trajectory_power_study/phase6_small_n_study.json`` at e = 0.
PHASE6_GENERATOR: dict[str, Any] = {
    "trajectory_mode": "magnitude",
    "group_effect_size": 0.0,
    "n_stages": 3,
    "group_stage_sizes": [[11, 10, 28], [9, 10, 12]],
    "p_dmp": 0.1,
    "baseline_continuity": 0.0,
    "magnitude_kind": "joint",
}

#: Stage-PLS cross-validation knobs of the Phase 6 study configs (= the analysis defaults).
PLS_CV_PARAMS: dict[str, Any] = {
    "cv1_splits": 3,
    "cv2_splits": 4,
    "n_repeats": 5,
    "max_components": 20,
    "random_state": 1203,
}


# ── Settings ─────────────────────────────────────────────────────────────────


def build_settings(args: argparse.Namespace) -> dict[str, Any]:
    return {
        "script_version": SCRIPT_VERSION,
        "replicates": int(args.reps),
        "permutations": int(args.perms),
        "base_seed": int(args.base_seed),
        "alpha": float(args.alpha),
        "k": float(args.k),
        "generator": dict(PHASE6_GENERATOR),
        "layers": list(LAYERS),
        "methylation_units": "mvalue",
        "cv": dict(PLS_CV_PARAMS),
        "orientation_estimator": "net_displacement",
        "covariates": {
            "numeric": {"a_g": float(args.a_g), "a_s": float(args.a_s), "a_gs": float(args.a_gs), "noise_sd": 1.0},
            "categorical": {
                "levels": list(CATEGORY_LEVELS),
                "logits": "k * logit_scale * eta",
                "logit_scale": float(args.categorical_logit_scale),
            },
        },
        "injection": {
            "s_z": float(args.s_z),
            "s_c": float(args.s_c),
            "feature_fraction": float(args.feature_fraction),
            "scale": "pooled feature SD on the analysed scale",
        },
    }


def settings_hash(settings: Mapping[str, Any]) -> str:
    encoded = json.dumps(settings, sort_keys=True, separators=(",", ":")).encode("utf-8")
    return hashlib.sha256(encoded).hexdigest()


def ensure_settings(out_dir: Path, settings: Mapping[str, Any]) -> None:
    """Write the settings header, or refuse a directory written with different settings."""

    path = out_dir / "settings.json"
    if path.exists():
        stored = json.loads(path.read_text(encoding="utf-8"))
        if stored.get("hash") != settings_hash(settings):
            old, new = stored.get("settings", {}), settings
            differing = sorted(key for key in set(old) | set(new) if old.get(key) != new.get(key))
            raise SystemExit(
                f"{path} was written with different settings; differing: {', '.join(differing)}. "
                "Use a fresh --out-dir."
            )
        return
    out_dir.mkdir(parents=True, exist_ok=True)
    payload = {"hash": settings_hash(settings), "settings": settings}
    path.write_text(json.dumps(payload, indent=1, sort_keys=True) + "\n", encoding="utf-8")


def load_settings(out_dir: Path) -> dict[str, Any]:
    path = out_dir / "settings.json"
    if not path.exists():
        raise SystemExit(f"{path} does not exist; nothing to merge.")
    stored = json.loads(path.read_text(encoding="utf-8"))
    if stored.get("hash") != settings_hash(stored["settings"]):
        raise SystemExit(f"{path} is corrupt: its hash does not match its settings.")
    return stored["settings"]


def _seed(*parts: Any) -> int:
    encoded = json.dumps(["covariate-type-i", *parts], separators=(",", ":")).encode("utf-8")
    return int(hashlib.sha256(encoded).hexdigest()[:8], 16) & 0x7FFFFFFF


def replicate_seeds(base_seed: int, replicate: int) -> dict[str, int]:
    return {
        "generator": _seed(base_seed, "generator", replicate),
        "covariates": _seed(base_seed, "covariates", replicate),
        "rrpp": _seed(base_seed, "rrpp", replicate),
    }


# ── Construction ─────────────────────────────────────────────────────────────


def draw_covariates(
    group: np.ndarray, stage: np.ndarray, settings: Mapping[str, Any], rng: np.random.Generator
) -> tuple[np.ndarray, np.ndarray]:
    """Numeric ``z`` and categorical labels whose distributions depend on group × stage."""

    numeric = settings["covariates"]["numeric"]
    eta = numeric["a_g"] * group + numeric["a_s"] * stage + numeric["a_gs"] * group * stage
    z = eta + numeric["noise_sd"] * rng.normal(size=group.size)
    scale = settings["covariates"]["categorical"]["logit_scale"]
    logits = np.outer(scale * eta, np.arange(len(CATEGORY_LEVELS), dtype=float))
    probabilities = np.exp(logits - logits.max(axis=1, keepdims=True))
    probabilities /= probabilities.sum(axis=1, keepdims=True)
    draws = np.array([rng.choice(len(CATEGORY_LEVELS), p=row) for row in probabilities])
    return z, np.asarray(CATEGORY_LEVELS, dtype=object)[draws]


def feature_directions(n_features: int, settings: Mapping[str, Any]) -> tuple[np.ndarray, np.ndarray]:
    """Unit directions ``w`` (numeric) and ``u_k`` (one per level), each on a random feature subset.

    Drawn from the base seed only, so every replicate shares the same directions.
    """

    rng = np.random.default_rng(_seed(settings["base_seed"], "directions", n_features))
    size = max(1, int(round(settings["injection"]["feature_fraction"] * n_features)))

    def direction() -> np.ndarray:
        vector = np.zeros(n_features)
        support = rng.choice(n_features, size=size, replace=False)
        vector[support] = rng.normal(size=size)
        return vector / np.linalg.norm(vector)

    w = direction()
    u = np.vstack([direction() for _ in CATEGORY_LEVELS])
    return w, u


def build_replicate(replicate: int, settings: Mapping[str, Any], reference: Any = None) -> dict[str, Any]:
    """Generate the null dataset, draw the covariates and inject their effect."""

    from motco.simulations.preprocessing import apply_unit_transform
    from motco.simulations.reference import load_reference
    from motco.simulations.semisynthetic import SemiSyntheticTrajectoryParams, generate_semisynthetic_trajectory

    seeds = replicate_seeds(int(settings["base_seed"]), replicate)
    generator = dict(settings["generator"])
    generator["group_stage_sizes"] = tuple(tuple(row) for row in generator["group_stage_sizes"])
    dataset = generate_semisynthetic_trajectory(
        SemiSyntheticTrajectoryParams(seed=seeds["generator"], **generator),
        reference=reference if reference is not None else load_reference(),
    )
    metadata = dataset.metadata.reset_index(drop=True)
    group_labels = sorted(set(metadata["group"].astype(str)))
    group = (metadata["group"].astype(str) == group_labels[1]).to_numpy(dtype=float)
    stage = metadata["stage"].to_numpy(dtype=float)
    z, category = draw_covariates(group, stage, settings, np.random.default_rng(seeds["covariates"]))

    raw = {
        "methylation": apply_unit_transform("methylation", dataset.methylation.to_numpy(dtype=float), "beta_to_mvalue"),
        "expression": dataset.expression.to_numpy(dtype=float),
    }
    frames = {"methylation": dataset.methylation, "expression": dataset.expression}
    n_features = sum(raw[name].shape[1] for name in settings["layers"])
    w, u = feature_directions(n_features, settings)
    indicators = np.vstack([(category == level).astype(float) for level in CATEGORY_LEVELS])  # levels × n
    s_z, s_c = settings["injection"]["s_z"], settings["injection"]["s_c"]
    blocks: dict[str, pd.DataFrame] = {}
    offset = 0
    for name in settings["layers"]:
        values = raw[name]
        width = values.shape[1]
        sd = values.std(axis=0, ddof=1)
        effect = s_z * np.outer(z, w[offset : offset + width]) + s_c * (indicators.T @ u[:, offset : offset + width])
        blocks[name] = pd.DataFrame(values + effect * sd, index=frames[name].index, columns=frames[name].columns)
        offset += width
    metadata = metadata.assign(z=z, category=category)
    return {"seeds": seeds, "metadata": metadata, "blocks": blocks, "group_labels": group_labels}


def covariate_cell_summary(metadata: pd.DataFrame, group_labels: Sequence[str]) -> dict[str, Any]:
    """Per-cell covariate means/proportions and the numeric group × stage interaction contrast."""

    stages = sorted(set(metadata["stage"].tolist()))
    means: dict[str, dict[str, float]] = {}
    proportions: dict[str, dict[str, dict[str, float]]] = {}
    for label in group_labels:
        rows = metadata["group"].astype(str) == label
        means[label] = {}
        proportions[label] = {}
        for stage in stages:
            cell = metadata[rows & (metadata["stage"] == stage)]
            means[label][str(stage)] = float(cell["z"].mean())
            proportions[label][str(stage)] = {
                level: float((cell["category"] == level).mean()) for level in CATEGORY_LEVELS
            }
    first, last = str(stages[0]), str(stages[-1])
    a, b = group_labels
    interaction = (means[b][last] - means[b][first]) - (means[a][last] - means[a][first])
    return {"numeric_means": means, "categorical_proportions": proportions, "numeric_interaction_contrast": interaction}


# ── Replicate evaluation ─────────────────────────────────────────────────────


def _jsonable(value: Any) -> Any:
    if isinstance(value, Mapping):
        return {str(k): _jsonable(v) for k, v in value.items()}
    if isinstance(value, (list, tuple)):
        return [_jsonable(v) for v in value]
    if isinstance(value, np.generic):
        return _jsonable(value.item())
    if isinstance(value, float) and not math.isfinite(value):
        return None
    return value


def evaluate_replicate(
    replicate: int, settings: Mapping[str, Any], *, n_jobs: int = 1, reference: Any = None
) -> dict[str, Any]:
    from motco.analysis import OmicsBlock, TrajectoryAnalysisParams, TrajectoryDataset, analyze_trajectories

    start = perf_counter()
    built = build_replicate(replicate, settings, reference=reference)
    metadata = built["metadata"]
    blocks = [OmicsBlock(name, built["blocks"][name], "identity") for name in settings["layers"]]
    stage_order = [str(stage) for stage in sorted(set(metadata["stage"].tolist()))]
    params = TrajectoryAnalysisParams(
        permutations=int(settings["permutations"]),
        seed=built["seeds"]["rrpp"],
        n_jobs=n_jobs,
        orientation_estimator=settings["orientation_estimator"],
        **settings["cv"],
    )
    arms: dict[str, Any] = {}
    results = {}
    for arm in ARMS:
        covariate_kwargs: dict[str, Any] = (
            {"covariate_cols": ["z", "category"], "categorical_covariates": ["category"]} if arm == "adjusted" else {}
        )
        dataset = TrajectoryDataset(
            blocks=blocks,
            metadata=metadata,
            sample_col="sample_id",
            group_col="group",
            stage_col="stage",
            stage_order=stage_order,
            **covariate_kwargs,
        )
        result = analyze_trajectories(dataset, params)
        results[arm] = result
        record = result.record
        arms[arm] = {
            "statistics": record["statistics"],
            "p_values": record["p_values"],
            "selected_lv": record["latent"]["selected_lv"],
            "has_covariate_block": "covariates" in record,
        }
        if "covariates" in record:
            arms[arm]["covariate_design_columns"] = record["covariates"]["design_columns"]
    adjusted, unadjusted = results["adjusted"], results["unadjusted"]
    if not (
        np.array_equal(adjusted.latent_scores.to_numpy(), unadjusted.latent_scores.to_numpy())
        and arms["adjusted"]["selected_lv"] == arms["unadjusted"]["selected_lv"]
    ):
        raise RuntimeError(f"replicate {replicate}: the adjusted and unadjusted latent spaces differ.")
    return _jsonable(
        {
            "replicate": replicate,
            "settings_hash": settings_hash(settings),
            "seeds": built["seeds"],
            "group_labels": built["group_labels"],
            "covariate_cells": covariate_cell_summary(metadata, built["group_labels"]),
            "arms": arms,
            "runtime_seconds": perf_counter() - start,
        }
    )


# ── Shards, resume, merge ────────────────────────────────────────────────────


def parse_shard(text: str) -> tuple[int, int]:
    try:
        index_text, count_text = text.split("/")
        index, count = int(index_text), int(count_text)
    except ValueError as exc:
        raise argparse.ArgumentTypeError(f"--shard must be i/N, got {text!r}.") from exc
    if count < 1 or not 0 <= index < count:
        raise argparse.ArgumentTypeError(f"--shard index must be in [0, N), got {text!r}.")
    return index, count


def shard_path(out_dir: Path, index: int, count: int) -> Path:
    return out_dir / "records" / f"shard_{index:03d}_of_{count:03d}.jsonl"


def read_records(paths: Iterable[Path]) -> dict[int, dict[str, Any]]:
    """Records keyed by replicate; torn lines (an interrupted write) are skipped and re-run."""

    records: dict[int, dict[str, Any]] = {}
    for path in paths:
        with path.open(encoding="utf-8") as handle:
            for line in handle:
                try:
                    record = json.loads(line)
                except json.JSONDecodeError:
                    continue
                if set(record.get("arms", {})) == set(ARMS):
                    records[int(record["replicate"])] = record
    return records


def run_shard(
    out_dir: Path, settings: Mapping[str, Any], index: int, count: int, *, n_jobs: int = 1, limit: int | None = None
) -> int:
    """Evaluate this shard's missing replicates; returns how many were evaluated.

    A replicate recorded in any shard file under ``records/`` counts as done, so
    a re-run with a different shard count (e.g. moving a run to a cluster array)
    skips work already recorded elsewhere; ``merge`` keeps one record per replicate.
    """

    path = shard_path(out_dir, index, count)
    path.parent.mkdir(parents=True, exist_ok=True)
    own = read_records([path]) if path.exists() else {}
    # Rewrite this shard with complete records only, so a torn line is re-run instead of duplicated.
    with path.open("w", encoding="utf-8") as handle:
        for replicate in sorted(own):
            handle.write(json.dumps(own[replicate], sort_keys=True, allow_nan=False) + "\n")
    others = [other for other in sorted(path.parent.glob("shard_*.jsonl")) if other != path]
    done = {**read_records(others), **own}
    pending = [r for r in range(int(settings["replicates"])) if r % count == index and r not in done]
    if limit is not None:
        pending = pending[:limit]
    from motco.simulations.reference import load_reference

    reference = load_reference()
    for position, replicate in enumerate(pending, start=1):
        record = evaluate_replicate(replicate, settings, n_jobs=n_jobs, reference=reference)
        with path.open("a", encoding="utf-8") as handle:
            handle.write(json.dumps(record, sort_keys=True, allow_nan=False) + "\n")
            handle.flush()
        print(f"[shard {index}/{count}] {position}/{len(pending)} replicate {replicate}", file=sys.stderr, flush=True)
    return len(pending)


def p_value_rows(records: Mapping[int, Mapping[str, Any]]) -> list[dict[str, Any]]:
    rows = []
    for replicate in sorted(records):
        for arm in ARMS:
            p_values = records[replicate]["arms"][arm]["p_values"]
            for statistic in STATISTICS:
                rows.append(
                    {"replicate": replicate, "arm": arm, "statistic": statistic, "p_value": p_values.get(statistic)}
                )
    return rows


# ── Decision (pure) ──────────────────────────────────────────────────────────


RATE_COLUMNS = ("arm", "statistic", "replicates", "rejections", "rate", "bound", "verdict")


def anchor_bound(alpha: float, replicates: int, k: float) -> float:
    return alpha + k * math.sqrt(alpha * (1.0 - alpha) / replicates)


def rejection_rates(rows: Iterable[Mapping[str, Any]], *, alpha: float, k: float) -> list[dict[str, Any]]:
    """Per arm × statistic: available replicates R, rejections (``p ≤ alpha``), rate and bound.

    Only the adjusted arm carries a verdict; the unadjusted arm is descriptive.
    """

    tallies: dict[tuple[str, str], list[int]] = {}
    for row in rows:
        p_value = row.get("p_value")
        if p_value is None or not math.isfinite(float(p_value)):
            continue
        tally = tallies.setdefault((str(row["arm"]), str(row["statistic"])), [0, 0])
        tally[0] += 1
        tally[1] += int(float(p_value) <= alpha)
    out = []
    for arm in ARMS:
        for statistic in STATISTICS:
            replicates, rejections = tallies.get((arm, statistic), [0, 0])
            rate = rejections / replicates if replicates else None
            bound = anchor_bound(alpha, replicates, k) if replicates else None
            if arm != "adjusted":
                verdict = "descriptive"
            elif rate is None or bound is None:
                verdict = "unavailable"
            else:
                verdict = "holds" if rate <= bound else "exceeds"
            out.append(
                {
                    "arm": arm,
                    "statistic": statistic,
                    "replicates": replicates,
                    "rejections": rejections,
                    "rate": rate,
                    "bound": bound,
                    "verdict": verdict,
                }
            )
    return out


def covariate_type_i_decision(rows: Iterable[Mapping[str, Any]], *, alpha: float, k: float) -> dict[str, Any]:
    """``covariates_validated`` when the adjusted arm holds the bound for every statistic."""

    rates = rejection_rates(rows, alpha=alpha, k=k)
    adjusted = {row["statistic"]: row for row in rates if row["arm"] == "adjusted"}
    unadjusted = {row["statistic"]: row for row in rates if row["arm"] == "unadjusted"}
    verdicts = {statistic: adjusted[statistic]["verdict"] for statistic in STATISTICS}
    validated = all(verdict == "holds" for verdict in verdicts.values())
    return {
        "rule": "adjusted-arm rate (p <= alpha) <= alpha + k * sqrt(alpha * (1 - alpha) / R) for every statistic",
        "alpha": alpha,
        "k": k,
        "replicates": {statistic: adjusted[statistic]["replicates"] for statistic in STATISTICS},
        "bound": {statistic: adjusted[statistic]["bound"] for statistic in STATISTICS},
        "adjusted": {statistic: adjusted[statistic]["rate"] for statistic in STATISTICS},
        "unadjusted": {statistic: unadjusted[statistic]["rate"] for statistic in STATISTICS},
        "verdicts": verdicts,
        "unadjusted_role": "descriptive (evidence the construction confounds); not part of the decision",
        "decision": "covariates_validated" if validated else "diagnose_covariate_excess",
    }


def _write_csv(path: Path, rows: Sequence[Mapping[str, Any]], columns: Sequence[str]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(columns), lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow({column: "" if row.get(column) is None else row[column] for column in columns})


def write_decision(out_dir: Path, rows: Sequence[Mapping[str, Any]], settings: Mapping[str, Any]) -> dict[str, Any]:
    alpha, k = float(settings["alpha"]), float(settings["k"])
    decision = covariate_type_i_decision(rows, alpha=alpha, k=k)
    decision["settings_hash"] = settings_hash(settings)
    decision["construction"] = {
        "generator": settings["generator"],
        "layers": settings["layers"],
        "covariates": settings["covariates"],
        "injection": settings["injection"],
        "permutations": settings["permutations"],
        "base_seed": settings["base_seed"],
    }
    (out_dir / "covariate_type_i_decision.json").write_text(
        json.dumps(_jsonable(decision), indent=2, allow_nan=False) + "\n", encoding="utf-8"
    )
    _write_csv(out_dir / "covariate_type_i_rates.csv", rejection_rates(rows, alpha=alpha, k=k), RATE_COLUMNS)
    return decision


def merge(out_dir: Path, settings: Mapping[str, Any]) -> dict[str, Any]:
    """Merge every shard, then write the p-value table, the rates and the decision."""

    records = read_records(sorted((out_dir / "records").glob("shard_*.jsonl")))
    expected = int(settings["replicates"])
    with (out_dir / "merged.jsonl").open("w", encoding="utf-8") as handle:
        for replicate in sorted(records):
            handle.write(json.dumps(records[replicate], sort_keys=True, allow_nan=False) + "\n")
    rows = p_value_rows(records)
    _write_csv(out_dir / "covariate_type_i_pvalues.csv", rows, ("replicate", "arm", "statistic", "p_value"))
    decision = write_decision(out_dir, rows, settings)
    contrasts = [records[r]["covariate_cells"]["numeric_interaction_contrast"] for r in records]
    return {
        "replicates": len(records),
        "expected": expected,
        "missing": sorted(set(range(expected)) - set(records))[:20],
        "mean_numeric_interaction_contrast": float(np.mean(contrasts)) if contrasts else None,
        "decision": decision["decision"],
        "adjusted": decision["adjusted"],
        "unadjusted": decision["unadjusted"],
        "bound": decision["bound"],
    }


# ── CLI ──────────────────────────────────────────────────────────────────────


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Covariate-adjustment Type I check at the Phase 6 design (no R).")
    sub = parser.add_subparsers(dest="command", required=True)
    run = sub.add_parser("run", help="Evaluate (a shard of) the replicates; an unsharded run also merges")
    run.add_argument("--out-dir", type=Path, required=True, help="Output directory")
    run.add_argument("--reps", type=int, default=500, help="Replicates (default: 500)")
    run.add_argument("--perms", type=int, default=999, help="RRPP permutations (default: 999)")
    run.add_argument("--base-seed", type=int, default=2026, help="Parent seed (default: 2026)")
    run.add_argument("--alpha", type=float, default=0.05, help="Rejection threshold, p <= alpha (default: 0.05)")
    run.add_argument("--k", type=float, default=2.0, help="Bound multiplier k (default: 2)")
    run.add_argument("--a-g", type=float, default=0.5, help="Covariate group effect a_g (default: 0.5)")
    run.add_argument("--a-s", type=float, default=1.0, help="Covariate stage effect a_s (default: 1.0)")
    run.add_argument("--a-gs", type=float, default=2.0, help="Covariate group x stage effect a_gs (default: 2.0)")
    run.add_argument(
        "--categorical-logit-scale", type=float, default=0.5, help="Categorical logit scale c (default: 0.5)"
    )
    run.add_argument("--s-z", type=float, default=3.0, help="Numeric covariate injection strength (default: 3.0)")
    run.add_argument("--s-c", type=float, default=2.0, help="Categorical covariate injection strength (default: 2.0)")
    run.add_argument(
        "--feature-fraction", type=float, default=0.2, help="Fraction of features per direction (default: 0.2)"
    )
    run.add_argument("--shard", type=parse_shard, default=(0, 1), help="Run shard i of N (i/N; default 0/1)")
    run.add_argument("--n-jobs", type=int, default=1, help="Jobs for CV and RRPP (not part of the settings)")
    run.add_argument("--limit", type=int, default=None, help="Evaluate at most this many pending replicates")
    merge_parser = sub.add_parser("merge", help="Merge every shard and write the decision (reads settings.json)")
    merge_parser.add_argument("--out-dir", type=Path, required=True, help="Output directory")
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    if args.command == "merge":
        settings = load_settings(args.out_dir)
    else:
        settings = build_settings(args)
        ensure_settings(args.out_dir, settings)
        index, count = args.shard
        evaluated = run_shard(args.out_dir, settings, index, count, n_jobs=args.n_jobs, limit=args.limit)
        print(f"evaluated {evaluated} replicate(s) into {shard_path(args.out_dir, index, count)}", file=sys.stderr)
        if count > 1 or args.limit is not None:
            return 0
    print(json.dumps(merge(args.out_dir, settings), indent=1))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
