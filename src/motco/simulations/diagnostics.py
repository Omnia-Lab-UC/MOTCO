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

#: Checkpoint of the population means seen through the replicate's fitted PLS basis.
POPULATION_LATENT_CHECKPOINT = "population_latent"

#: A population group difference whose stage-varying part is at most this
#: fraction of the configuration's scale is constant across stages (the anchor,
#: a translation): its span-capture share is undefined and reported as null.
_CONSTANT_DIFFERENCE_TOLERANCE = 1e-9


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
class SpanCapture:
    """Shares of the population trajectories inside the fitted PLS span.

    Measured in the fitted model's input coordinates with an orthonormal basis
    of its rotations. ``group_difference`` is the share of
    ``Σ_s ‖d_s‖²`` in the span, ``d_s`` the population B − A difference at
    stage ``s`` minus its mean over stages (``None`` when the difference is the
    same at every stage); ``pooled_stage`` the share of the centred average
    stage means. ``random_share`` = ``rank / n_features`` is the expected share
    of a random direction.
    """

    group_difference: float | None
    pooled_stage: float | None
    rank: int
    n_features: int
    random_share: float

    def to_dict(self) -> dict[str, Any]:
        return asdict(self)


@dataclass(frozen=True)
class RealizedGeometryDiagnostics:
    """Nested checkpoint decomposition with construction attribution.

    ``span_capture`` is present only when it was requested for a PLS
    evaluation; without it the record is unchanged.
    """

    schema_version: int
    requested: dict[str, Any]
    checkpoints: dict[str, dict[str, GeometryDiagnostic]]
    span_capture: SpanCapture | None = None

    def to_dict(self) -> dict[str, Any]:
        payload: dict[str, Any] = {
            "schema_version": self.schema_version,
            "requested": self.requested,
            "checkpoints": {
                checkpoint: {scope: diagnostic.to_dict() for scope, diagnostic in scopes.items()}
                for checkpoint, scopes in self.checkpoints.items()
            },
        }
        if self.span_capture is not None:
            payload["span_capture"] = self.span_capture.to_dict()
        return payload


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
    latent_model: Any | None = None,
) -> RealizedGeometryDiagnostics:
    """Calculate all applicable geometry checkpoints for one replicate.

    Only the layers the preprocessor was fitted on are measured: per-block
    scopes cover those layers, and every joint scope concatenates them alone.
    ``latent_matrix``, when given, is measured under the checkpoint named
    ``latent_checkpoint``. Every checkpoint ``angle`` uses
    ``orientation_estimator``, so the checkpoints measure the tested estimand.

    ``latent_model`` is the replicate's fitted PLS model, passed only when span
    capture is requested. With population trajectories it adds the
    ``population_latent`` checkpoint (the population means through the fitted
    preprocessor and ``model.transform``) and the :class:`SpanCapture` shares.
    Nothing is refitted.
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

    span_capture: SpanCapture | None = None
    if latent_model is not None and population is not None:
        population_joint = concatenate_blocks(preprocessor.transform_population(population))
        ordered = np.vstack(
            [
                np.asarray(population_joint.loc[(group, stage)], dtype=float)
                for group in group_levels
                for stage in stage_levels
            ]
        )
        scores = np.asarray(latent_model.transform(ordered), dtype=float)
        checkpoints[POPULATION_LATENT_CHECKPOINT] = {"joint": measure(scores)}
        span_capture = measure_span_capture(latent_model, ordered, n_groups=len(group_levels))

    requested = {
        "trajectory_mode": dataset.truth.get("trajectory_mode"),
        "group_effect_size": dataset.truth.get("group_effect_size"),
        "transform": dataset.truth.get("transform", {}),
    }
    return RealizedGeometryDiagnostics(
        schema_version=DIAGNOSTIC_SCHEMA_VERSION,
        requested=requested,
        checkpoints=checkpoints,
        span_capture=span_capture,
    )


def measure_span_capture(model: Any, means: np.ndarray, *, n_groups: int = 2) -> SpanCapture:
    """Shares of population group-stage means inside a fitted PLS model's span.

    ``means`` holds the group-stage means in the model's feature order, group
    major and stage minor (the first two groups are compared, B − A). They are
    mapped to the model's input coordinates, ``(x − x_mean) / x_std``, and
    measured against an orthonormal basis of ``model.x_rotations_`` (QR). For a
    single direction ``d`` the share equals the square of the plasmode
    ``latent_capture.py`` norm share ``‖Qᵀd‖ / ‖d‖``.
    """

    rotations = np.asarray(model.x_rotations_, dtype=float)
    basis, _ = np.linalg.qr(rotations)
    inputs = (np.asarray(means, dtype=float) - np.asarray(model._x_mean)) / np.asarray(model._x_std)
    n_features = int(rotations.shape[0])
    rank = int(rotations.shape[1])
    n_stages = inputs.shape[0] // n_groups
    group_a = inputs[:n_stages]
    group_b = inputs[n_stages : 2 * n_stages]
    difference = group_b - group_a
    varying = difference - difference.mean(axis=0)
    pooled = (group_a + group_b) / 2.0
    pooled = pooled - pooled.mean(axis=0)
    scale = max(float(np.linalg.norm(difference)), float(np.linalg.norm(pooled)), np.finfo(float).tiny)
    constant = float(np.linalg.norm(varying)) <= _CONSTANT_DIFFERENCE_TOLERANCE * scale
    return SpanCapture(
        group_difference=None if constant else span_share(varying, basis),
        pooled_stage=span_share(pooled, basis) if float(np.linalg.norm(pooled)) > 0 else None,
        rank=rank,
        n_features=n_features,
        random_share=rank / n_features,
    )


def span_share(vectors: np.ndarray, basis: np.ndarray) -> float:
    """Share of ``Σ ‖v‖²`` over the rows of ``vectors`` inside ``span(basis)``.

    ``basis`` has orthonormal columns. The result is clipped to ``[0, 1]``
    against rounding.
    """

    rows = np.atleast_2d(np.asarray(vectors, dtype=float))
    total = float(np.sum(rows**2))
    inside = float(np.sum((rows @ basis) ** 2))
    return float(min(max(inside / total, 0.0), 1.0))


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
