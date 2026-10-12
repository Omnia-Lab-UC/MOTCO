"""Realized trajectory-geometry decomposition for semi-synthetic studies."""

from __future__ import annotations

from dataclasses import asdict, dataclass
from typing import Any, Mapping

import numpy as np
import pandas as pd

from motco.simulations.preprocessing import (
    FittedOmicsPreprocessor,
    concatenate_blocks,
)
from motco.simulations.semisynthetic import OmicsLayer, SemiSyntheticTrajectoryDataset
from motco.stats.trajectory import estimate_difference, get_observed_vectors

DIAGNOSTIC_SCHEMA_VERSION = 1
_PATH_TOLERANCE = 1e-15


@dataclass(frozen=True)
class GeometryDiagnostic:
    """JSON-safe scalar geometry for one checkpoint and scope."""

    path_lengths: dict[str, float]
    delta: float
    angle: float | None
    shape: float | None
    availability: dict[str, bool]

    def to_dict(self) -> dict[str, Any]:
        return asdict(self)


@dataclass(frozen=True)
class RealizedGeometryDiagnostics:
    """Nested checkpoint decomposition with construction attribution."""

    schema_version: int
    requested: dict[str, Any]
    checkpoints: dict[str, dict[str, GeometryDiagnostic]]

    def to_dict(self) -> dict[str, Any]:
        return {
            "schema_version": self.schema_version,
            "requested": self.requested,
            "checkpoints": {
                checkpoint: {scope: diagnostic.to_dict() for scope, diagnostic in scopes.items()}
                for checkpoint, scopes in self.checkpoints.items()
            },
        }


def calculate_realized_geometry(
    dataset: SemiSyntheticTrajectoryDataset,
    *,
    preprocessor: FittedOmicsPreprocessor,
    observed_blocks: Mapping[OmicsLayer, pd.DataFrame],
    latent_matrix: pd.DataFrame | None = None,
    group_col: str = "group",
    stage_col: str = "stage",
    latent_checkpoint: str = "pls_latent",
    orientation_estimator: str = "pc1",
) -> RealizedGeometryDiagnostics:
    """Calculate all applicable geometry checkpoints for one replicate.

    Only the layers the preprocessor was fitted on are measured: per-block
    scopes cover those layers, and every joint scope concatenates them alone.
    ``latent_matrix``, when given, is measured under the checkpoint named
    ``latent_checkpoint``. Every checkpoint ``angle`` uses
    ``orientation_estimator``, so the checkpoints measure the tested estimand.
    """

    def measure(means: pd.DataFrame | np.ndarray) -> GeometryDiagnostic:
        return geometry_from_means(means, group_levels, stage_levels, orientation_estimator=orientation_estimator)

    layers = preprocessor.layers

    group_levels = sorted(pd.unique(dataset.metadata[group_col].astype(str)).tolist())
    stage_levels = sorted(pd.unique(dataset.metadata[stage_col].astype(str)).tolist())
    checkpoints: dict[str, dict[str, GeometryDiagnostic]] = {}

    population = dataset.population_trajectories
    if population is not None:
        checkpoints["population_native"] = {
            layer: measure(population.layers[layer])
            for layer in layers
        }
        standardized_population = preprocessor.transform_population(population)
        checkpoints["population_standardized"] = {
            layer: measure(standardized_population[layer])
            for layer in layers
        }
        checkpoints["population_standardized"]["joint"] = measure(concatenate_blocks(standardized_population))

    observed_scopes: dict[str, GeometryDiagnostic] = {}
    for layer in layers:
        means = get_observed_vectors(
            dataset.metadata,
            observed_blocks[layer],
            group_col=group_col,
            level_col=stage_col,
        )
        observed_scopes[layer] = measure(means)
    observed_joint = concatenate_blocks({layer: observed_blocks[layer] for layer in layers})
    joint_means = get_observed_vectors(
        dataset.metadata,
        observed_joint,
        group_col=group_col,
        level_col=stage_col,
    )
    observed_scopes["joint"] = measure(joint_means)
    checkpoints["observed_standardized"] = observed_scopes

    if latent_matrix is not None:
        latent_means = get_observed_vectors(
            dataset.metadata,
            latent_matrix,
            group_col=group_col,
            level_col=stage_col,
        )
        checkpoints[latent_checkpoint] = {"joint": measure(latent_means)}

    requested = {
        "trajectory_mode": dataset.truth.get("trajectory_mode"),
        "group_effect_size": dataset.truth.get("group_effect_size"),
        "transform": dataset.truth.get("transform", {}),
    }
    return RealizedGeometryDiagnostics(
        schema_version=DIAGNOSTIC_SCHEMA_VERSION,
        requested=requested,
        checkpoints=checkpoints,
    )


def geometry_from_means(
    means: pd.DataFrame | np.ndarray,
    group_levels: list[str],
    stage_levels: list[str],
    *,
    orientation_estimator: str = "pc1",
) -> GeometryDiagnostic:
    """Measure production trajectory estimands from ordered group-stage means.

    ``angle`` is unavailable when a group's path has zero length, or when the
    selected orientation estimator leaves it undefined (a zero net
    displacement under ``"net_displacement"``).
    """

    if isinstance(means, pd.DataFrame) and isinstance(means.index, pd.MultiIndex):
        ordered = np.vstack(
            [np.asarray(means.loc[(group, stage)], dtype=float) for group in group_levels for stage in stage_levels]
        )
    else:
        ordered = np.asarray(means, dtype=float)
    n_groups = len(group_levels)
    n_stages = len(stage_levels)
    expected_rows = n_groups * n_stages
    if ordered.shape[0] != expected_rows:
        raise ValueError(f"Expected {expected_rows} group-stage means, got {ordered.shape[0]}.")

    contrast = [list(range(group * n_stages, (group + 1) * n_stages)) for group in range(n_groups)]
    identity = np.eye(expected_rows)
    deltas, angles, shapes = estimate_difference(
        ordered, identity, identity, contrast, orientation_estimator=orientation_estimator
    )
    path_lengths = {
        group: _path_length(ordered[index * n_stages : (index + 1) * n_stages])
        for index, group in enumerate(group_levels)
    }
    angle_available = all(length > _PATH_TOLERANCE for length in path_lengths.values()) and bool(
        np.isfinite(angles[0, 1])
    )
    shape_available = n_stages >= 3
    return GeometryDiagnostic(
        path_lengths=path_lengths,
        delta=float(deltas[0, 1]),
        angle=float(angles[0, 1]) if angle_available else None,
        shape=float(shapes[0, 1]) if shape_available else None,
        availability={"delta": True, "angle": angle_available, "shape": shape_available},
    )


def flatten_geometry_diagnostics(diagnostics: Mapping[str, Any]) -> dict[str, float | int | str | None]:
    """Flatten scalar diagnostics for tabular Phase 2 summaries."""

    flattened: dict[str, float | int | str | None] = {
        "diagnostic_schema_version": diagnostics.get("schema_version"),
    }
    requested = diagnostics.get("requested", {})
    flattened["trajectory_mode"] = requested.get("trajectory_mode")
    flattened["group_effect_size"] = requested.get("group_effect_size")
    for checkpoint, scopes in diagnostics.get("checkpoints", {}).items():
        for scope, values in scopes.items():
            prefix = f"{checkpoint}.{scope}"
            for group, length in values.get("path_lengths", {}).items():
                flattened[f"{prefix}.path_length.{group}"] = length
            for statistic in ("delta", "angle", "shape"):
                flattened[f"{prefix}.{statistic}"] = values.get(statistic)
    return flattened


def _path_length(trajectory: np.ndarray) -> float:
    if trajectory.shape[0] < 2:
        return 0.0
    return float(np.linalg.norm(np.diff(trajectory, axis=0), axis=1).sum())
