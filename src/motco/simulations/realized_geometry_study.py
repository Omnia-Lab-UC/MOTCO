"""Matched-seed Phase 2 characterization of realized generator geometry."""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any, Iterable, Mapping

import pandas as pd

from motco.simulations.diagnostics import flatten_geometry_diagnostics
from motco.simulations.evaluation import SimulationEvaluationParams, evaluate_semisynthetic_trajectory
from motco.simulations.reference import IntersimReference, load_reference
from motco.simulations.semisynthetic import SemiSyntheticTrajectoryParams, generate_semisynthetic_trajectory

PHASE2_CONSTRUCTIONS: tuple[tuple[str, str], ...] = (
    ("none", "relocate"),
    ("translation", "relocate"),
    ("magnitude", "relocate"),
    ("orientation", "relocate"),
    ("shape", "relocate"),
    ("shape", "magnitude"),
)


@dataclass(frozen=True)
class Phase2CharacterizationConfig:
    """Reproducible parameter grid for the realized-geometry study."""

    effect_sizes: tuple[float, ...] = (0.0, 0.25, 0.5, 0.75, 1.0)
    seeds: tuple[int, ...] = (0, 1, 2)
    n_samples: int = 300
    n_stages: int = 4
    p_dmp: float = 0.2
    integration_params: Mapping[str, Any] = field(default_factory=dict)


def run_phase2_characterization(
    config: Phase2CharacterizationConfig,
    *,
    reference: IntersimReference | None = None,
) -> list[dict[str, Any]]:
    """Run every construction/effect cell on matched generator seeds."""

    ref = reference if reference is not None else load_reference()
    rows: list[dict[str, Any]] = []
    evaluation = SimulationEvaluationParams(
        integration_method="pls",
        integration_params=dict(config.integration_params),
        permutations=0,
        orientation_estimator="pc1",  # historical results measured PC1 angle
        min_components=1,  # historical results ran the unfloored rank rule
    )
    for mode, shape_kind in PHASE2_CONSTRUCTIONS:
        for effect_size in config.effect_sizes:
            for seed in config.seeds:
                params = SemiSyntheticTrajectoryParams(
                    seed=seed,
                    trajectory_mode=mode,  # type: ignore[arg-type]
                    n_samples=config.n_samples,
                    n_stages=config.n_stages,
                    group_effect_size=effect_size,
                    p_dmp=config.p_dmp,
                    shape_kind=shape_kind,  # type: ignore[arg-type]
                    # Committed runs predate the joint default; keep their construction.
                    magnitude_kind="all",
                )
                result = evaluate_semisynthetic_trajectory(
                    generate_semisynthetic_trajectory(params, reference=ref),
                    evaluation,
                )
                row: dict[str, Any] = {
                    "seed": seed,
                    "mode": mode,
                    "shape_kind": shape_kind if mode == "shape" else None,
                    "effect_size": effect_size,
                    "n_samples": config.n_samples,
                    "n_stages": config.n_stages,
                    "p_dmp": config.p_dmp,
                    "selected_lv": result.latent_matrix_metadata.get("integration_params", {}).get(
                        "selected_lv"
                    ),
                }
                row.update(flatten_geometry_diagnostics(result.realized_geometry))
                rows.append(row)
    return rows


def summarize_phase2_characterization(rows: Iterable[Mapping[str, Any]]) -> pd.DataFrame:
    """Return long-form mean/SD summaries for every diagnostic scalar."""

    records = list(rows)
    long_rows: list[dict[str, Any]] = []
    identifiers = {"seed", "mode", "shape_kind", "effect_size", "n_samples", "n_stages", "p_dmp"}
    for row in records:
        for key, value in row.items():
            if key in identifiers or key in {"selected_lv", "trajectory_mode", "group_effect_size"}:
                continue
            parts = key.split(".")
            if len(parts) < 3 or parts[0] == "diagnostic_schema_version":
                continue
            checkpoint, scope = parts[0], parts[1]
            if parts[2] == "path_length" and len(parts) == 4:
                statistic = f"path_length_{parts[3]}"
            elif len(parts) == 3:
                statistic = parts[2]
            else:
                continue
            long_rows.append(
                {
                    "mode": row["mode"],
                    "shape_kind": row.get("shape_kind"),
                    "effect_size": float(row["effect_size"]),
                    "checkpoint": checkpoint,
                    "scope": scope,
                    "statistic": statistic,
                    "value": value,
                }
            )
    frame = pd.DataFrame(long_rows)
    if frame.empty:
        return pd.DataFrame(
            columns=[
                "mode",
                "shape_kind",
                "effect_size",
                "checkpoint",
                "scope",
                "statistic",
                "mean",
                "std",
                "n_available",
            ]
        )
    frame["value"] = pd.to_numeric(frame["value"], errors="coerce")
    grouped = frame.groupby(
        ["mode", "shape_kind", "effect_size", "checkpoint", "scope", "statistic"],
        dropna=False,
        sort=True,
    )["value"]
    return grouped.agg(mean="mean", std="std", n_available="count").reset_index()


def monotonicity_summary(summary: pd.DataFrame) -> pd.DataFrame:
    """Summarize requested-effect monotonicity using rank correlation."""

    rows = []
    keys = ["mode", "shape_kind", "checkpoint", "scope", "statistic"]
    for key, group in summary.groupby(keys, dropna=False, sort=True):
        valid = group[["effect_size", "mean"]].dropna().sort_values("effect_size")
        correlation = (
            float(valid["effect_size"].corr(valid["mean"], method="spearman"))
            if len(valid) >= 2 and valid["mean"].nunique() >= 2
            else None
        )
        rows.append(dict(zip(keys, key, strict=True), spearman=correlation, n_effects=len(valid)))
    return pd.DataFrame(rows)
