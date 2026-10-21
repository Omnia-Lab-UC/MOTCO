"""Evaluation harness for semi-synthetic trajectory datasets.

Latent-space architecture
--------------------------
The trajectory pipeline is **features -> latent (molecular) space -> trajectory
measurement**. The integration step constructs the low-dimensional *molecular
latent space*, and trajectory geometry (``delta``/``angle``/``shape``) is
estimated **within** that space — it is the measurement substrate, not a
visualization artifact. The integration methods are:

- ``pls`` (the production measurement space): the transform of the omic
  features into the subspace that maximizes covariance with the **stage**
  label, with latent dimensionality selected by double nested cross-validation
  (a stable, non-overfitted molecular space).
- ``snf``: graph-spectral embedding of the fused per-omic similarity networks.
  Supported, but its use with the Euclidean trajectory statistics is deferred
  pending graph-native statistics — ``docs/roadmap.md`` is the authoritative
  status.

``concat`` is a **baseline/diagnostic** path only: it standardizes and
column-binds the omic blocks, leaving the data in rescaled *feature* space, not
a constructed latent space. The viz down-projection (``viz.plot_trajectory_from_*``)
is display-only and distinct from this measurement space.

Two further **diagnostic** methods exist for measurement-design comparison
(the orientation measurement probe) and are not production latent spaces:

- ``pca``: label-free principal-component scores of the ``concat`` matrix, at
  a required ``integration_params["n_components"]``.
- ``pls_within_group``: a stage-supervised PLS fitted separately within each
  group (``components_per_group``, default ``n_stages - 1``), whose stacked
  weight vectors are orthonormalized into one basis spanning both groups'
  stage subspaces. The basis is fitted once per evaluation, never inside RRPP.
"""

from __future__ import annotations

import hashlib
from dataclasses import dataclass, field
from time import perf_counter
from typing import Any, Literal, Mapping, Sequence

import numpy as np
import pandas as pd

from motco.simulations.diagnostics import calculate_realized_geometry
from motco.simulations.generator import logit
from motco.simulations.preprocessing import (
    OMIC_LAYERS,
    BlockScaler,
    FittedOmicsPreprocessor,
    concatenate_blocks,
    fit_omics_preprocessor,
    selected_layers,
)
from motco.simulations.semisynthetic import OmicsLayer, SemiSyntheticTrajectoryDataset
from motco.stats.design import _sort_levels, build_ls_means, get_model_matrix
from motco.stats.permutation import RRPP
from motco.stats.pls import (
    NuisanceDesign,
    _modal_int_with_parsimony,
    fit_nuisance_coefficients,
    fit_plsda_model,
    nuisance_partial_r2,
    plsda_doubleCV,
    residualize,
)
from motco.stats.snf import SNF, get_affinity_matrix, get_spectral
from motco.stats.trajectory import estimate_difference

IntegrationMethod = Literal["concat", "snf", "pls", "pca", "pls_within_group"]

#: Integration methods the harness accepts.
INTEGRATION_METHODS: frozenset[str] = frozenset({"concat", "snf", "pls", "pca", "pls_within_group"})

#: Probe-only latent spaces: labelled ``integration_role="diagnostic"``, never
#: production, and excluded from attribution and the forced-rank override.
DIAGNOSTIC_INTEGRATION_METHODS: frozenset[str] = frozenset({"pca", "pls_within_group"})

#: Methods whose latent matrix is measured as its own realized-geometry
#: checkpoint, keyed to that checkpoint's name.
_LATENT_CHECKPOINTS: dict[str, str] = {
    "pls": "pls_latent",
    "pca": "pca_latent",
    "pls_within_group": "pls_within_group_latent",
}

#: Orientation estimators an evaluation may select. ``None`` (the default)
#: means PC1, and an explicit ``"pc1"`` is normalized to ``None``.
EVALUATION_ORIENTATION_ESTIMATORS: tuple[str, ...] = ("pc1", "net_displacement")


class SimulationEvaluationError(ValueError):
    """Raised when simulation evaluation inputs or parameters are invalid."""


@dataclass(frozen=True)
class AttributionDiagnosticSettings:
    """Optional orientation-attribution diagnostic settings for one evaluation.

    Attribution is **disabled by default**: the zero-argument instance carries
    ``enabled=False`` and no bootstrap work, so evaluations that do not ask for
    diagnostics behave exactly as before. When enabled the settings are frozen
    into the cell signature, so a diagnostic change cannot silently resume into
    an existing shard.
    """

    enabled: bool = False
    bootstrap_replicates: int = 0
    bootstrap_seed: int | None = 0
    top_k: int = 20
    zero_tolerance: float = 1e-12


@dataclass(frozen=True)
class SimulationEvaluationParams:
    """Parameters for evaluating one semi-synthetic trajectory dataset."""

    integration_method: IntegrationMethod = "concat"
    integration_params: Mapping[str, Any] = field(default_factory=dict)
    permutations: int = 0
    n_jobs: int | None = 1
    seed: int | None = None
    progress: bool = False
    include_null_distributions: bool = False
    group_col: str = "group"
    stage_col: str = "stage"
    attribution: AttributionDiagnosticSettings = field(default_factory=AttributionDiagnosticSettings)
    # Orientation estimator for the observed ``angle``, every RRPP draw and every
    # realized-geometry checkpoint. The default is net displacement; it is hashed
    # into cell ids and signatures and recorded in ``runtime_metadata``. ``None``
    # is the hash-omitted spelling of PC1 (``grid._SIGNATURE_OMIT_WHEN_NONE``) and
    # an explicit ``"pc1"`` normalizes to it, so PC1 cell ids, signatures and
    # records equal those produced before the field existed.
    orientation_estimator: str | None = "net_displacement"

    def __post_init__(self) -> None:
        estimator = self.orientation_estimator
        if estimator is None:
            return
        if estimator not in EVALUATION_ORIENTATION_ESTIMATORS:
            raise SimulationEvaluationError(
                f"Unknown orientation_estimator {estimator!r}; allowed: "
                f"{list(EVALUATION_ORIENTATION_ESTIMATORS)}."
            )
        if estimator == "pc1":
            object.__setattr__(self, "orientation_estimator", None)

    @property
    def resolved_orientation_estimator(self) -> str:
        """The estimator name the statistics are computed with."""

        return self.orientation_estimator or "pc1"


#: Version of the JSON-safe integration metadata contract persisted per replicate.
INTEGRATION_METADATA_VERSION = 1


#: Quantiles of the per-feature partial R² recorded for a latent adjustment.
_R2_QUANTILES: tuple[tuple[str, float], ...] = (
    ("q05", 0.05),
    ("q25", 0.25),
    ("q50", 0.50),
    ("q75", 0.75),
    ("q95", 0.95),
)


@dataclass(frozen=True)
class LatentAdjustmentSummary:
    """How the features were residualized before a covariate-aware PLS basis was fitted.

    ``design_columns`` are the residualization design columns beyond the
    intercept (conditioning, then covariates); ``residual_df`` is the
    full-sample fit's ``n − 1 − k − q``; ``cv_min_training_df`` the smallest
    training-set residual df seen during cross-validation (``None`` when the
    rank was forced and no CV ran); ``feature_r2`` summarizes the per-feature
    partial R² of the covariates given the intercept and conditioning design.
    """

    conditioning: str | None
    design_columns: tuple[str, ...]
    n_samples: int
    residual_df: int
    cv_min_training_df: int | None
    feature_scale: str
    feature_r2: dict[str, float]

    @classmethod
    def from_partial_r2(
        cls,
        r2: np.ndarray,
        *,
        nuisance: NuisanceDesign,
        residual_df: int,
        cv_min_training_df: int | None,
        feature_scale: str = "pooled_standardized",
    ) -> "LatentAdjustmentSummary":
        values = np.asarray(r2, dtype=float)
        summary = {name: float(np.quantile(values, q)) for name, q in _R2_QUANTILES}
        summary["mean"] = float(values.mean())
        summary["fraction_above_0.5"] = float(np.mean(values > 0.5))
        return cls(
            conditioning=nuisance.kind,
            design_columns=nuisance.column_names,
            n_samples=nuisance.n_samples,
            residual_df=residual_df,
            cv_min_training_df=cv_min_training_df,
            feature_scale=feature_scale,
            feature_r2=summary,
        )


@dataclass(frozen=True)
class PLSIntegrationArtifacts:
    """Ephemeral fitted artifacts behind one PLS latent space.

    These objects are the analysis boundary attribution conditions on: the exact
    standardized joint feature matrix the latent scores are the projection of
    (the matrix the model was fitted to, unless the basis was fitted on
    covariate-residualized features), that fitted estimator, and the per-feature
    scales needed to express effects in original units. ``latent_adjustment``
    describes the residualization when there was one. They are deliberately
    **not** part of the persisted evaluation result — a fitted sklearn estimator
    is neither JSON-safe nor version-stable — and are discarded once the compact
    diagnostics are built.
    """

    model: Any
    features: pd.DataFrame
    original_scales: np.ndarray
    methylation_units: str
    latent_adjustment: LatentAdjustmentSummary | None = None


@dataclass(frozen=True)
class LatentIntegrationResult:
    """Integrated latent/outcome matrix and metadata.

    ``artifacts`` carries the ephemeral fitted objects for the PLS path and is
    never serialized; only ``metadata`` reaches the persisted record.
    """

    matrix: pd.DataFrame
    metadata: dict[str, Any]
    artifacts: PLSIntegrationArtifacts | None = None


@dataclass(frozen=True)
class SimulationTrajectoryDesign:
    """Trajectory design objects derived from sample metadata."""

    model_full: np.ndarray
    model_reduced: np.ndarray
    ls_means: np.ndarray
    contrast: list[list[int]]
    group_levels: list[str]
    stage_levels: list[str]


@dataclass(frozen=True)
class SimulationEvaluationResult:
    """Result from evaluating one semi-synthetic trajectory dataset."""

    observed_deltas: np.ndarray
    observed_angles: np.ndarray
    observed_shapes: np.ndarray
    pair_statistics: dict[str, float]
    p_values: dict[str, float]
    latent_matrix_metadata: dict[str, Any]
    truth_metadata: dict[str, Any]
    runtime_metadata: dict[str, Any]
    evaluation_params: SimulationEvaluationParams
    group_levels: list[str]
    stage_levels: list[str]
    contrast: list[list[int]]
    null_distributions: dict[str, list[float]] | None = None
    realized_geometry: dict[str, Any] = field(default_factory=dict)
    attribution_diagnostics: dict[str, Any] = field(default_factory=dict)
    # Compact per-statistic description of the permutation null. Produced on
    # every RRPP run — unlike ``null_distributions``, which stays opt-in — so a
    # record can locate its own observed statistic against its own critical
    # value without retaining the full draw vectors.
    null_summary: dict[str, dict[str, float]] = field(default_factory=dict)
    # Eigenspectrum of the centered stage-mean configuration in the *evaluated
    # latent space*, pooled and per contrast group, plus — whenever permutations
    # ran — a summary of the pooled relative eigengap over the permutation draws.
    # A recorded covariate: the geometry audit
    # (``docs/reports/geometry-audit-2026-09-01.md``) found the relative eigengap
    # predicts how wide a replicate's own ``angle`` permutation null will be, and
    # therefore whether its orientation is resolvable at all. It qualifies
    # reporting; it enters no statistic and no decision rule.
    config_spectrum: dict[str, Any] = field(default_factory=dict)


def evaluate_semisynthetic_trajectory(
    dataset: SemiSyntheticTrajectoryDataset,
    params: SimulationEvaluationParams | None = None,
    *,
    latent: LatentIntegrationResult | None = None,
) -> SimulationEvaluationResult:
    """Evaluate one semi-synthetic trajectory dataset through MOTCO routines.

    ``latent`` optionally supplies an integration already built for this
    dataset with the same ``integration_method`` and ``integration_params``
    (from :func:`integrate_semisynthetic_dataset`), so evaluations that differ
    only in measurement settings — the orientation estimator, say — share one
    space fit. Omitted (the default), the integration is built here.
    """

    params = params or SimulationEvaluationParams()
    _validate_evaluation_params(params)
    start = perf_counter()

    preprocessor = fit_omics_preprocessor(dataset, _measured_layers(params.integration_params))
    observed_blocks = preprocessor.transform_dataset(dataset)
    if latent is None:
        latent = integrate_semisynthetic_dataset(
            dataset,
            params,
            preprocessor=preprocessor,
            observed_blocks=observed_blocks,
        )
    elif latent.metadata.get("integration_method") != params.integration_method:
        raise SimulationEvaluationError(
            f"Supplied latent integration was built with {latent.metadata.get('integration_method')!r}, "
            f"not integration_method={params.integration_method!r}."
        )
    design = build_simulation_trajectory_design(
        dataset.metadata,
        group_col=params.group_col,
        stage_col=params.stage_col,
    )

    observed_deltas, observed_angles, observed_shapes, observed_spectra = estimate_difference(
        Y=latent.matrix,
        model_matrix=design.model_full,
        LS_means=design.ls_means,
        contrast=design.contrast,
        return_spectra=True,
        orientation_estimator=params.resolved_orientation_estimator,
    )
    config_spectrum = _config_spectrum_block(observed_spectra, design.group_levels)
    shape_available = len(design.stage_levels) >= 3
    pair_statistics = _extract_pair_statistics(
        observed_deltas,
        observed_angles,
        observed_shapes,
        shape_available=shape_available,
    )

    p_values: dict[str, float] = {}
    null_distributions: dict[str, list[float]] | None = None
    null_summary: dict[str, dict[str, float]] = {}
    if params.permutations > 0:
        dist_delta, dist_angle, dist_shape, dist_eigengap = RRPP(
            latent.matrix,
            design.model_full,
            design.model_reduced,
            design.ls_means,
            design.contrast,
            permutations=params.permutations,
            n_jobs=params.n_jobs,
            progress=params.progress,
            seed=params.seed,
            return_eigengaps=True,
            orientation_estimator=params.resolved_orientation_estimator,
        )
        null_values = _extract_null_distributions(
            dist_delta,
            dist_angle,
            dist_shape,
            shape_available=shape_available,
        )
        p_values = {
            statistic: _empirical_p_value(null_values[statistic], observed)
            for statistic, observed in pair_statistics.items()
            if statistic in null_values and np.isfinite(observed)
        }
        null_summary = _summarize_null_distributions(null_values)
        # Only the summary survives: the per-permutation eigengaps are discarded
        # here, exactly as the per-permutation statistics are.
        config_spectrum["permutation_pooled_eigengap"] = _summarize_permutation_eigengaps(
            dist_eigengap
        )
        if params.include_null_distributions:
            null_distributions = null_values

    runtime_seconds = perf_counter() - start
    runtime_metadata: dict[str, Any] = {
        "runtime_seconds": runtime_seconds,
        "permutations": params.permutations,
        "n_jobs": params.n_jobs,
        "seed": params.seed,
        "progress": params.progress,
        "shape_available": shape_available,
    }
    if params.orientation_estimator is not None:
        runtime_metadata["orientation_estimator"] = params.orientation_estimator
    latent_checkpoint = _LATENT_CHECKPOINTS.get(params.integration_method)
    realized_geometry = calculate_realized_geometry(
        dataset,
        preprocessor=preprocessor,
        observed_blocks=observed_blocks,
        latent_matrix=latent.matrix if latent_checkpoint is not None else None,
        group_col=params.group_col,
        stage_col=params.stage_col,
        latent_checkpoint=latent_checkpoint or "pls_latent",
        orientation_estimator=params.resolved_orientation_estimator,
    ).to_dict()
    # The fitted estimator and standardized joint matrix live only inside this
    # call: attribution consumes them here, and only the compact JSON-safe
    # record survives into the result.
    attribution_diagnostics = _attribution_diagnostics(dataset, params, latent)

    return SimulationEvaluationResult(
        observed_deltas=observed_deltas,
        observed_angles=observed_angles,
        observed_shapes=observed_shapes,
        pair_statistics=pair_statistics,
        p_values=p_values,
        latent_matrix_metadata=latent.metadata,
        # Drop the raw indicator arrays — they are not JSON-serializable for the
        # study's JSONL persistence; the per-stage/group counts and transform
        # notes (which are) are retained for downstream characterization.
        truth_metadata={k: v for k, v in dataset.truth.items() if k != "indicators"},
        runtime_metadata=runtime_metadata,
        evaluation_params=params,
        group_levels=design.group_levels,
        stage_levels=design.stage_levels,
        contrast=design.contrast,
        null_distributions=null_distributions,
        realized_geometry=realized_geometry,
        attribution_diagnostics=attribution_diagnostics,
        null_summary=null_summary,
        config_spectrum=config_spectrum,
    )


def _attribution_diagnostics(
    dataset: SemiSyntheticTrajectoryDataset,
    params: SimulationEvaluationParams,
    latent: LatentIntegrationResult,
) -> dict[str, Any]:
    """Compute compact orientation-attribution diagnostics when they are enabled."""

    from motco.simulations.attribution_diagnostics import (
        AttributionDiagnosticError,
        compute_attribution_diagnostics,
        unavailable_record,
    )

    settings = params.attribution
    if not settings.enabled:
        return unavailable_record("attribution diagnostics were not requested for this cell")
    artifacts = latent.artifacts
    if artifacts is None:
        raise SimulationEvaluationError(
            "Attribution diagnostics were requested but the integration produced no fitted "
            f"PLS artifacts (integration_method={params.integration_method!r})."
        )
    try:
        return compute_attribution_diagnostics(
            dataset,
            model=artifacts.model,
            features=artifacts.features,
            original_scales=artifacts.original_scales,
            settings=settings,
            group_col=params.group_col,
            stage_col=params.stage_col,
            selected_components=latent.metadata.get("selected_lv"),
            feature_order_signature=latent.metadata.get("feature_order_signature"),
            methylation_units=artifacts.methylation_units,
            layers=_measured_layers(params.integration_params),
        )
    except AttributionDiagnosticError as exc:
        # A diagnostic is auxiliary evidence: a replicate whose attribution
        # cannot be computed still contributes a valid trajectory measurement,
        # so the failure is recorded rather than allowed to abort the run.
        record = unavailable_record(str(exc), status="failed")
        record["error_type"] = type(exc).__name__
        return record


def integrate_semisynthetic_dataset(
    dataset: SemiSyntheticTrajectoryDataset,
    params: SimulationEvaluationParams,
    *,
    preprocessor: FittedOmicsPreprocessor | None = None,
    observed_blocks: Mapping[OmicsLayer, pd.DataFrame] | None = None,
) -> LatentIntegrationResult:
    """Construct the molecular latent space from aligned omics layers.

    ``snf`` and ``pls`` build genuine latent spaces (the measurement substrate
    for trajectory geometry); ``concat`` is a standardized-feature-concatenation
    baseline, not a constructed latent space. See the module docstring for the
    latent-space architecture.
    """

    _validate_dataset(dataset, params.group_col, params.stage_col)
    layers = _measured_layers(params.integration_params)
    if preprocessor is None:
        preprocessor = fit_omics_preprocessor(dataset, layers)
    else:
        try:
            preprocessor.require_layers(layers)
        except ValueError as exc:
            raise SimulationEvaluationError(str(exc)) from exc
    if observed_blocks is None:
        observed_blocks = preprocessor.transform_dataset(dataset)
    method = params.integration_method
    if method == "concat":
        return _concat_integration(dataset, params.integration_params, observed_blocks)
    if method == "snf":
        return _snf_integration(dataset, params.integration_params)
    if method == "pls":
        return _pls_integration(
            dataset,
            params.integration_params,
            stage_col=params.stage_col,
            observed_blocks=observed_blocks,
            preprocessor=preprocessor,
        )
    if method == "pca":
        return _pca_integration(dataset, params.integration_params, observed_blocks)
    if method == "pls_within_group":
        return _pls_within_group_integration(
            dataset,
            params.integration_params,
            group_col=params.group_col,
            stage_col=params.stage_col,
            observed_blocks=observed_blocks,
        )
    raise SimulationEvaluationError(f"Unsupported integration_method: {method!r}.")


def build_simulation_trajectory_design(
    metadata: pd.DataFrame,
    *,
    group_col: str = "group",
    stage_col: str = "stage",
) -> SimulationTrajectoryDesign:
    """Build full/reduced model matrices, LS means, and two-group contrast."""

    _validate_metadata(metadata, group_col, stage_col)
    design_frame = metadata[[group_col, stage_col]].reset_index(drop=True).copy()
    design_frame[group_col] = design_frame[group_col].astype(str)
    design_frame[stage_col] = design_frame[stage_col].astype(str)
    group_levels = _sort_levels(pd.unique(design_frame[group_col]).tolist())
    stage_levels = _sort_levels(pd.unique(design_frame[stage_col]).tolist())
    if len(group_levels) != 2:
        raise SimulationEvaluationError(f"Expected exactly two groups, found {len(group_levels)}: {group_levels}.")
    if len(stage_levels) < 2:
        raise SimulationEvaluationError(
            f"Expected at least two stages, found {len(stage_levels)}: {stage_levels}."
        )
    _validate_group_stage_combinations(design_frame, group_col, stage_col, group_levels, stage_levels)

    model_full = get_model_matrix(design_frame, group_col=group_col, level_col=stage_col, full=True)
    model_reduced = get_model_matrix(design_frame, group_col=group_col, level_col=stage_col, full=False)
    ls_means = build_ls_means(group_levels, stage_levels, full=True)
    n_stages = len(stage_levels)
    contrast = [[group_index * n_stages + stage_index for stage_index in range(n_stages)] for group_index in range(2)]
    return SimulationTrajectoryDesign(
        model_full=model_full,
        model_reduced=model_reduced,
        ls_means=ls_means,
        contrast=contrast,
        group_levels=group_levels,
        stage_levels=stage_levels,
    )


def _validate_evaluation_params(params: SimulationEvaluationParams) -> None:
    if params.integration_method not in INTEGRATION_METHODS:
        raise SimulationEvaluationError(f"Unsupported integration_method: {params.integration_method!r}.")
    if (
        params.integration_method in DIAGNOSTIC_INTEGRATION_METHODS
        and params.integration_params.get("forced_components") is not None
    ):
        raise SimulationEvaluationError(
            "integration_params['forced_components'] is specific to integration_method='pls'; "
            f"got it with the diagnostic method {params.integration_method!r}."
        )
    if params.permutations < 0:
        raise SimulationEvaluationError("permutations must be greater than or equal to 0.")
    if params.n_jobs == 0:
        raise SimulationEvaluationError("n_jobs must be None, -1, or a non-zero integer.")
    _measured_layers(params.integration_params)
    _validate_attribution_settings(params.attribution, params.integration_method)


def _measured_layers(integration_params: Mapping[str, Any]) -> tuple[OmicsLayer, ...]:
    """The omic blocks this evaluation measures (``integration_params["layers"]``)."""

    try:
        return selected_layers(integration_params)
    except ValueError as exc:
        raise SimulationEvaluationError(str(exc)) from exc


def _layers_metadata(integration_params: Mapping[str, Any]) -> dict[str, Any]:
    """The ``layers`` metadata entry — present only when a selection was declared.

    Absent selection adds nothing, so the metadata of an evaluation that never
    set the key is unchanged.
    """

    if integration_params.get("layers") is None:
        return {}
    return {"layers": list(_measured_layers(integration_params))}


def _validate_attribution_settings(
    settings: AttributionDiagnosticSettings,
    integration_method: IntegrationMethod,
) -> None:
    if not settings.enabled:
        return
    if integration_method != "pls":
        raise SimulationEvaluationError(
            "Orientation attribution diagnostics require integration_method='pls'; "
            f"got {integration_method!r}. Attribution conditions on the fitted production PLS "
            "estimator, which no other integration method supplies."
        )
    if settings.bootstrap_replicates < 0:
        raise SimulationEvaluationError("attribution.bootstrap_replicates must be greater than or equal to 0.")
    if settings.top_k < 1:
        raise SimulationEvaluationError("attribution.top_k must be at least 1.")
    if not (settings.zero_tolerance > 0):
        raise SimulationEvaluationError("attribution.zero_tolerance must be positive.")


def _validate_dataset(dataset: SemiSyntheticTrajectoryDataset, group_col: str, stage_col: str) -> None:
    _validate_metadata(dataset.metadata, group_col, stage_col)
    if "sample_id" not in dataset.metadata.columns:
        raise SimulationEvaluationError("metadata must contain required column 'sample_id'.")
    expected = dataset.metadata["sample_id"].astype(str).tolist()
    if len(expected) == 0:
        raise SimulationEvaluationError("metadata must contain at least one sample.")
    for layer in OMIC_LAYERS:
        matrix = getattr(dataset, layer)
        if matrix.shape[0] != len(expected):
            raise SimulationEvaluationError(f"{layer} rows must match metadata rows.")
        if matrix.index.astype(str).tolist() != expected:
            raise SimulationEvaluationError(f"{layer} rows are not aligned to metadata sample_id order.")
        if matrix.shape[1] == 0:
            raise SimulationEvaluationError(f"{layer} must contain at least one feature.")
        values = matrix.to_numpy(dtype=float)
        if not np.all(np.isfinite(values)):
            raise SimulationEvaluationError(f"{layer} contains NaN or Inf values.")


def _validate_metadata(metadata: pd.DataFrame, group_col: str, stage_col: str) -> None:
    required = [group_col, stage_col]
    missing = [column for column in required if column not in metadata.columns]
    if missing:
        raise SimulationEvaluationError(f"metadata is missing required column(s): {missing}.")


def _validate_group_stage_combinations(
    metadata: pd.DataFrame,
    group_col: str,
    stage_col: str,
    group_levels: list[str],
    stage_levels: list[str],
) -> None:
    present = set(zip(metadata[group_col].astype(str), metadata[stage_col].astype(str)))
    missing = [(group, stage) for group in group_levels for stage in stage_levels if (group, stage) not in present]
    if missing:
        raise SimulationEvaluationError(f"metadata is missing group/stage combination(s): {missing}.")


def _concat_integration(
    dataset: SemiSyntheticTrajectoryDataset,
    integration_params: Mapping[str, Any],
    observed_blocks: Mapping[OmicsLayer, pd.DataFrame],
) -> LatentIntegrationResult:
    standardize = bool(integration_params.get("standardize", True))
    layers = _measured_layers(integration_params)
    layer_feature_counts = {layer: int(getattr(dataset, layer).shape[1]) for layer in layers}
    if standardize:
        latent = concatenate_blocks(observed_blocks)
    else:
        frames = []
        for layer in layers:
            matrix = getattr(dataset, layer).astype(float)
            values = matrix.to_numpy(dtype=float)
            if layer == "methylation":
                values = logit(values)
            columns = [f"{layer}__{column}" for column in matrix.columns.astype(str)]
            frames.append(pd.DataFrame(values, index=matrix.index.astype(str), columns=columns))
        latent = pd.concat(frames, axis=1)
    return LatentIntegrationResult(
        matrix=latent,
        metadata={
            "integration_method": "concat",
            "integration_role": "baseline",
            "integration_params": {"standardize": standardize, **_layers_metadata(integration_params)},
            "shape": tuple(latent.shape),
            "n_samples": int(latent.shape[0]),
            "n_features": int(latent.shape[1]),
            "layer_feature_counts": layer_feature_counts,
        },
    )


def _snf_integration(
    dataset: SemiSyntheticTrajectoryDataset,
    integration_params: Mapping[str, Any],
) -> LatentIntegrationResult:
    n_samples = int(dataset.metadata.shape[0])
    K = _bounded_int_param(integration_params, "K", default=min(20, n_samples - 1), minimum=1, maximum=n_samples - 1)
    k = _bounded_int_param(integration_params, "k", default=min(20, n_samples - 1), minimum=1, maximum=n_samples - 1)
    t = _bounded_int_param(integration_params, "t", default=20, minimum=1)
    n_components = _bounded_int_param(
        integration_params,
        "spectral_components",
        default=min(10, n_samples - 1),
        minimum=1,
        maximum=n_samples - 1,
    )
    eps = float(integration_params.get("eps", 0.5))
    if eps <= 0:
        raise SimulationEvaluationError("SNF integration parameter 'eps' must be positive.")

    layers = [
        logit(getattr(dataset, layer).to_numpy(dtype=float)) if layer == "methylation"
        else getattr(dataset, layer).to_numpy(dtype=float)
        for layer in _measured_layers(integration_params)
    ]
    affinities = get_affinity_matrix(layers, K=K, eps=eps)
    fused = SNF(affinities, k=k, t=t)
    embedding = get_spectral(fused, n_components=n_components)
    latent = pd.DataFrame(
        embedding,
        index=dataset.metadata["sample_id"].astype(str).tolist(),
        columns=[f"snf_{i}" for i in range(embedding.shape[1])],
    )
    return LatentIntegrationResult(
        matrix=latent,
        metadata={
            "integration_method": "snf",
            "integration_role": "latent_space",
            "integration_params": {
                "K": K,
                "eps": eps,
                "k": k,
                "t": t,
                "spectral_components": n_components,
                **_layers_metadata(integration_params),
            },
            "shape": tuple(latent.shape),
            "n_samples": int(latent.shape[0]),
            "n_features": int(latent.shape[1]),
            "fused_shape": tuple(fused.shape),
        },
    )


def _pls_integration(
    dataset: SemiSyntheticTrajectoryDataset,
    integration_params: Mapping[str, Any],
    *,
    stage_col: str,
    observed_blocks: Mapping[OmicsLayer, pd.DataFrame],
    preprocessor: FittedOmicsPreprocessor,
) -> LatentIntegrationResult:
    """Build the PLS latent space for a semi-synthetic dataset.

    Adapter over :func:`fit_pls_latent_space`: reads the dataset's blocks,
    stage labels and sample ids, and hands them to the shared core.
    """
    if stage_col not in dataset.metadata.columns:
        raise SimulationEvaluationError(
            f"PLS integration requires stage column {stage_col!r} in dataset metadata."
        )

    X = concatenate_blocks(observed_blocks).reset_index(drop=True)
    return fit_pls_latent_space(
        X,
        dataset.metadata[stage_col].reset_index(drop=True).astype(str),
        integration_params,
        stage_col=stage_col,
        sample_ids=dataset.metadata["sample_id"].astype(str).tolist(),
        layer_feature_counts={
            layer: int(getattr(dataset, layer).shape[1]) for layer in _measured_layers(integration_params)
        },
        layers_metadata=_layers_metadata(integration_params),
        original_scales=_joint_original_scales(preprocessor.scalers, preprocessor.layers, X.columns),
        methylation_units=preprocessor.methylation_units,
    )


def fit_pls_latent_space(
    X: pd.DataFrame,
    y: pd.Series,
    integration_params: Mapping[str, Any],
    *,
    stage_col: str,
    sample_ids: Sequence[str],
    layer_feature_counts: Mapping[str, int],
    layers_metadata: Mapping[str, Any],
    original_scales: np.ndarray,
    methylation_units: str,
    nuisance: NuisanceDesign | None = None,
) -> LatentIntegrationResult:
    """Build the PLS molecular latent space (stage-conditioned, double-CV-sized).

    ``X`` is the pooled standardized joint feature matrix, row-aligned to ``y``
    and ``sample_ids``. Fit PLS-DA conditioned on the stage label ``y`` and
    return the X-score matrix as the latent space. The number of latent
    variables is selected by ``plsda_doubleCV`` (modal LV across repeats,
    parsimony tie-break) so the molecular space is stable and non-overfitted —
    the space in which trajectory geometry is then measured.

    ``integration_params["forced_components"]`` overrides that selection for
    rank diagnostics only: cross-validation is skipped entirely and the pooled
    model is fitted at the requested rank. It is rejected rather than clamped
    when infeasible — a diagnostic that silently moved the rank it was asked to
    measure would defeat its purpose — and the metadata records
    ``component_selection="forced"`` so no consumer can read a forced-rank run
    as a production evaluation. Absent (the default), behavior is unchanged and
    the metadata records ``component_selection="cv"``.

    ``nuisance`` (row-aligned to ``X``; the real-data latent adjustment)
    fits the basis on covariate-residualized features: the covariate
    coefficients are estimated on all samples conditional on the nuisance
    design's conditioning columns, the cross-validation residualizes every
    training set with its own coefficients (``plsda_doubleCV(...,
    nuisance=)``) on folds stratified on the conditioning cells (so every
    training set contains every cell; under stage conditioning the folds are
    those of the unadjusted fit), and every PLS fit uses ``scale=False`` so the
    residuals keep the pooled scale. The latent matrix is the raw ``X``
    projected onto that basis (``model.transform(X)``), the metadata gains
    ``integration_params["latent_adjustment"]`` (the conditioning) and the
    artifacts carry a :class:`LatentAdjustmentSummary`. Without ``nuisance``
    nothing changes.

    This is the production measurement shared by the simulation harness and
    the real-data entry point (:mod:`motco.analysis`); it never reads a dataset.
    """

    n_samples, n_features = X.shape
    stage_counts = y.value_counts()
    if len(stage_counts) < 2 or int(stage_counts.min()) < 2:
        raise SimulationEvaluationError(
            "PLS integration requires at least two stages with >= 2 samples each; "
            f"got stage counts {stage_counts.to_dict()}."
        )
    min_stage_count = int(stage_counts.min())
    max_feasible_components = min(n_features, n_samples)

    X_fit = X
    scale = True
    cv_nuisance: dict[str, Any] = {}
    full_residual_df = 0
    if nuisance is not None:
        if nuisance.n_samples != n_samples:
            raise SimulationEvaluationError(
                f"PLS integration: the nuisance design has {nuisance.n_samples} rows but X has {n_samples}."
            )
        try:
            coefficients, full_residual_df = fit_nuisance_coefficients(X, nuisance)
        except ValueError as exc:
            raise SimulationEvaluationError(f"PLS integration latent adjustment failed: {exc}") from exc
        X_fit = residualize(X, nuisance.covariates, coefficients)
        scale = False
        # Stratify the folds on the conditioning cells so every training set can estimate them.
        cv_nuisance = {"scale": False, "nuisance": nuisance, "strata": nuisance.cells()}

    forced_components = _forced_components(integration_params, maximum=max_feasible_components)
    cv_params: dict[str, Any]
    if forced_components is not None:
        # Diagnostic path: no cross-validation runs, so none of its knobs are
        # resolved or recorded.
        selected_lv = forced_components
        mean_auroc = None
        cv_repeats_completed = 0
        cv_params = {"forced_components": forced_components}
    else:
        # Cross-validation knobs (configurable; clamped to a feasible range so the
        # stratified folds stay valid on small samples). Effective values are
        # recorded in the result metadata.
        max_components = _clamp_int(
            integration_params.get("max_components", 20), minimum=2, maximum=max_feasible_components
        )
        cv2_splits = _clamp_int(integration_params.get("cv2_splits", 4), minimum=2, maximum=min_stage_count)
        # Samples per stage left after removing one outer test fold (conservative).
        rest_min_per_stage = min_stage_count - -(-min_stage_count // cv2_splits)
        if rest_min_per_stage < 2:
            raise SimulationEvaluationError(
                "PLS integration cross-validation is infeasible: too few samples per stage "
                f"(min {min_stage_count}) for cv2_splits={cv2_splits}."
            )
        cv1_splits = _clamp_int(integration_params.get("cv1_splits", 3), minimum=2, maximum=rest_min_per_stage)
        n_repeats = _clamp_int(integration_params.get("n_repeats", 5), minimum=1, maximum=1000)
        random_state = int(integration_params.get("random_state", 1203))
        n_jobs = int(integration_params.get("n_jobs", 1))

        try:
            cv_result = plsda_doubleCV(
                X,
                y,
                cv1_splits=cv1_splits,
                cv2_splits=cv2_splits,
                n_repeats=n_repeats,
                max_components=max_components,
                random_state=random_state,
                n_jobs=n_jobs,
                progress=False,
                **cv_nuisance,
            )
        except ValueError as exc:
            raise SimulationEvaluationError(f"PLS integration cross-validation failed: {exc}") from exc

        selected_lv = _modal_int_with_parsimony([int(v) for v in cv_result["table"]["LV"].tolist()])
        mean_auroc = float(cv_result["table"]["AUROC"].mean())
        cv_repeats_completed = int(cv_result["table"].shape[0])
        cv_params = {
            "cv1_splits": cv1_splits,
            "cv2_splits": cv2_splits,
            "n_repeats": n_repeats,
            "max_components": max_components,
            "random_state": random_state,
        }

    # One final pooled fit. Its training scores are the trajectory measurement
    # matrix *and* the same fitted estimator is handed to attribution, so both
    # necessarily share one coordinate system; no second fit or second component
    # selection ever happens. With a nuisance design the basis is fitted on the
    # residualized features and the raw features are projected onto it.
    model = fit_plsda_model(X_fit, y, n_components=selected_lv, scale=scale)
    adjustment: LatentAdjustmentSummary | None = None
    adjustment_params: dict[str, Any] = {}
    if nuisance is not None:
        scores = np.asarray(model.transform(X.to_numpy(dtype=float)))
        adjustment = LatentAdjustmentSummary.from_partial_r2(
            nuisance_partial_r2(X, nuisance),
            nuisance=nuisance,
            residual_df=full_residual_df,
            cv_min_training_df=(
                None if forced_components is not None else int(cv_result["nuisance_min_training_df"])
            ),
        )
        adjustment_params = {"latent_adjustment": nuisance.kind}
    else:
        scores = np.asarray(model.x_scores_)
    latent = pd.DataFrame(
        scores,
        index=list(sample_ids),
        columns=[f"pls_{i}" for i in range(scores.shape[1])],
    )
    return LatentIntegrationResult(
        matrix=latent,
        metadata={
            "integration_method": "pls",
            "integration_role": "latent_space",
            "integration_metadata_version": INTEGRATION_METADATA_VERSION,
            "integration_params": {
                "stage_col": stage_col,
                "selected_lv": selected_lv,
                **cv_params,
                **layers_metadata,
                **adjustment_params,
            },
            "component_selection": "forced" if forced_components is not None else "cv",
            "cv_mean_auroc": mean_auroc,
            "cv_repeats_completed": cv_repeats_completed,
            "selected_lv": selected_lv,
            "shape": tuple(latent.shape),
            "n_samples": int(latent.shape[0]),
            "n_features": int(latent.shape[1]),
            "layer_feature_counts": dict(layer_feature_counts),
            "feature_order_signature": _feature_order_signature(X.columns),
        },
        artifacts=PLSIntegrationArtifacts(
            model=model,
            features=X,
            original_scales=original_scales,
            methylation_units=methylation_units,
            latent_adjustment=adjustment,
        ),
    )


def _pca_integration(
    dataset: SemiSyntheticTrajectoryDataset,
    integration_params: Mapping[str, Any],
    observed_blocks: Mapping[OmicsLayer, pd.DataFrame],
) -> LatentIntegrationResult:
    """Diagnostic label-free PCA of exactly the matrix ``concat`` produces.

    The blocks are already pooled-standardized, so no further scaling is
    applied: the scores are a plain SVD of the column-centered matrix. The rank
    is a required ``integration_params["n_components"]`` — there is no rank
    rule, so the space type is not confounded with one.
    """

    layers = _measured_layers(integration_params)
    layer_feature_counts = {layer: int(getattr(dataset, layer).shape[1]) for layer in layers}
    X = concatenate_blocks(observed_blocks)
    n_samples, n_features = X.shape
    n_components = _required_rank_param(
        integration_params, "n_components", method="pca", maximum=min(n_samples, n_features)
    )
    values = X.to_numpy(dtype=float)
    centered = values - values.mean(axis=0, keepdims=True)
    U, singular_values, _ = np.linalg.svd(centered, full_matrices=False)
    scores = U[:, :n_components] * singular_values[:n_components]
    latent = pd.DataFrame(
        scores,
        index=dataset.metadata["sample_id"].astype(str).tolist(),
        columns=[f"pca_{i}" for i in range(n_components)],
    )
    return LatentIntegrationResult(
        matrix=latent,
        metadata={
            "integration_method": "pca",
            "integration_role": "diagnostic",
            "integration_metadata_version": INTEGRATION_METADATA_VERSION,
            "integration_params": {"n_components": n_components, **_layers_metadata(integration_params)},
            "retained_dimensionality": n_components,
            "shape": tuple(latent.shape),
            "n_samples": int(latent.shape[0]),
            "n_features": int(latent.shape[1]),
            "layer_feature_counts": layer_feature_counts,
        },
    )


#: Relative tolerance below which a pivoted-QR diagonal counts as rank-deficient.
_BASIS_RANK_TOLERANCE = 1e-10


def _pls_within_group_integration(
    dataset: SemiSyntheticTrajectoryDataset,
    integration_params: Mapping[str, Any],
    *,
    group_col: str,
    stage_col: str,
    observed_blocks: Mapping[OmicsLayer, pd.DataFrame],
) -> LatentIntegrationResult:
    """Diagnostic per-group stage-supervised PLS, merged into one basis.

    Fitted once from the observed data and labels; RRPP never refits it.
    """

    layers = _measured_layers(integration_params)
    layer_feature_counts = {layer: int(getattr(dataset, layer).shape[1]) for layer in layers}
    X = concatenate_blocks(observed_blocks)
    groups = dataset.metadata[group_col].astype(str).to_numpy()
    stages = dataset.metadata[stage_col].astype(str).to_numpy()
    n_stages = len(pd.unique(stages))
    raw = integration_params.get("components_per_group")
    components_per_group = _required_rank_param(
        {"components_per_group": n_stages - 1 if raw is None else raw},
        "components_per_group",
        method="pls_within_group",
        maximum=int(X.shape[1]),
    )
    basis = within_group_pls_basis(
        X.to_numpy(dtype=float), groups, stages, components_per_group=components_per_group
    )
    values = X.to_numpy(dtype=float)
    scores = (values - values.mean(axis=0, keepdims=True)) @ basis
    rank = int(basis.shape[1])
    latent = pd.DataFrame(
        scores,
        index=dataset.metadata["sample_id"].astype(str).tolist(),
        columns=[f"plswg_{i}" for i in range(rank)],
    )
    return LatentIntegrationResult(
        matrix=latent,
        metadata={
            "integration_method": "pls_within_group",
            "integration_role": "diagnostic",
            "integration_metadata_version": INTEGRATION_METADATA_VERSION,
            "integration_params": {
                "stage_col": stage_col,
                "components_per_group": components_per_group,
                **_layers_metadata(integration_params),
            },
            "components_per_group": components_per_group,
            "basis_rank": rank,
            "retained_dimensionality": rank,
            "shape": tuple(latent.shape),
            "n_samples": int(latent.shape[0]),
            "n_features": int(latent.shape[1]),
            "layer_feature_counts": layer_feature_counts,
        },
    )


def within_group_pls_basis(
    X: np.ndarray,
    groups: Sequence[str] | np.ndarray,
    stages: Sequence[str] | np.ndarray,
    *,
    components_per_group: int,
) -> np.ndarray:
    """Orthonormal basis spanning every group's stage-supervised PLS subspace.

    For each group (sorted), a PLS-DA on stage is fitted with ``scale=False`` to
    that group's rows of ``X`` — per-group rescaling would give the groups
    different metrics. The groups' ``x_rotations_`` are stacked column-wise and
    orthonormalized by pivoted QR; columns whose diagonal falls below a relative
    tolerance are dropped, so the returned ``p x r`` basis has rank ``r <=
    n_groups * components_per_group``.
    """

    from scipy.linalg import qr

    values = np.asarray(X, dtype=float)
    group_labels = np.asarray(groups).astype(str)
    stage_labels = np.asarray(stages).astype(str)
    rotations = []
    for group in _sort_levels(pd.unique(group_labels).tolist()):
        mask = group_labels == group
        counts = pd.Series(stage_labels[mask]).value_counts()
        all_stages = _sort_levels(pd.unique(stage_labels).tolist())
        short = [stage for stage in all_stages if int(counts.get(stage, 0)) < 2]
        if short:
            raise SimulationEvaluationError(
                f"pls_within_group integration requires at least two samples per stage in every group; "
                f"group {group!r} has fewer than two in stage(s) {short}."
            )
        if components_per_group > int(mask.sum()):
            raise SimulationEvaluationError(
                f"pls_within_group integration: components_per_group={components_per_group} exceeds "
                f"group {group!r}'s {int(mask.sum())} samples."
            )
        model = fit_plsda_model(values[mask], stage_labels[mask], n_components=components_per_group, scale=False)
        rotations.append(np.asarray(model.x_rotations_, dtype=float))
    stacked = np.hstack(rotations)
    Q, R, _ = qr(stacked, mode="economic", pivoting=True)
    diagonal = np.abs(np.diag(R))
    keep = diagonal > _BASIS_RANK_TOLERANCE * (diagonal[0] if diagonal.size else 0.0)
    return np.asarray(Q[:, keep])


def _required_rank_param(params: Mapping[str, Any], key: str, *, method: str, maximum: int) -> int:
    """Resolve a required integer rank in ``[1, maximum]``, rejecting rather than clamping."""

    raw = params.get(key)
    if raw is None:
        raise SimulationEvaluationError(f"{method} integration requires integration_params[{key!r}].")
    try:
        value = int(raw)
        exact = value == raw and not isinstance(raw, bool)
    except (TypeError, ValueError) as exc:
        raise SimulationEvaluationError(f"{method} integration: {key} must be an integer, got {raw!r}.") from exc
    if not exact:
        raise SimulationEvaluationError(f"{method} integration: {key} must be an integer, got {raw!r}.")
    if value < 1 or value > maximum:
        raise SimulationEvaluationError(
            f"{method} integration: {key}={value} is outside the feasible range [1, {maximum}] for this sample."
        )
    return value


def _feature_order_signature(columns: pd.Index) -> str:
    """Stable digest of the canonical joint feature order.

    Persisted so top-k truncation and feature alignment stay auditable without
    storing the full feature list in every record.
    """

    encoded = "\n".join(str(column) for column in columns).encode("utf-8")
    return hashlib.sha256(encoded).hexdigest()[:16]


def _joint_original_scales(
    scalers: Mapping[Any, BlockScaler],
    order: Sequence[Any],
    columns: pd.Index,
) -> np.ndarray:
    """Per-feature original-unit scales in joint feature order.

    ``order`` is the block order the joint matrix was concatenated in.
    Methylation scales are M-value scales, matching the fitted preprocessing
    contract; the unit basis is labeled wherever these are used.
    """

    scales = np.concatenate([scalers[name].scale for name in order])
    if scales.shape[0] != len(columns):
        raise SimulationEvaluationError(
            f"Fitted scales ({scales.shape[0]}) do not align to the joint feature matrix ({len(columns)})."
        )
    return scales


def _forced_components(params: Mapping[str, Any], *, maximum: int) -> int | None:
    """Resolve the diagnostic ``forced_components`` override, or ``None`` when absent.

    Unlike the cross-validation knobs, an out-of-range value is rejected rather
    than clamped: the override exists to measure the latent space at a named
    rank, so silently measuring a different one would be worse than failing.
    """

    raw = params.get("forced_components")
    if raw is None:
        return None
    try:
        value = int(raw)
        exact = value == raw
    except (TypeError, ValueError) as exc:
        raise SimulationEvaluationError(
            f"PLS integration: forced_components must be an integer, got {raw!r}."
        ) from exc
    if not exact:
        raise SimulationEvaluationError(
            f"PLS integration: forced_components must be an integer, got {raw!r}."
        )
    if value < 2 or value > maximum:
        raise SimulationEvaluationError(
            f"PLS integration: forced_components={value} is outside the feasible range "
            f"[2, {maximum}] for this sample."
        )
    return value


def _clamp_int(value: Any, *, minimum: int, maximum: int) -> int:
    if maximum < minimum:
        raise SimulationEvaluationError(
            f"PLS integration: infeasible cross-validation bound (minimum {minimum} > maximum {maximum})."
        )
    return max(minimum, min(int(value), maximum))


def _bounded_int_param(
    params: Mapping[str, Any],
    key: str,
    *,
    default: int,
    minimum: int,
    maximum: int | None = None,
) -> int:
    value = int(params.get(key, default))
    if value < minimum:
        raise SimulationEvaluationError(f"SNF integration parameter {key!r} must be at least {minimum}.")
    if maximum is not None and value > maximum:
        raise SimulationEvaluationError(f"SNF integration parameter {key!r} must be at most {maximum}.")
    return value


def _extract_pair_statistics(
    deltas: np.ndarray,
    angles: np.ndarray,
    shapes: np.ndarray,
    *,
    shape_available: bool,
) -> dict[str, float]:
    stats = {
        "delta": float(deltas[0, 1]),
        "angle": float(angles[0, 1]),
        "shape": float(shapes[0, 1]) if shape_available else float("nan"),
    }
    return stats


def _extract_null_distributions(
    dist_delta: list[np.ndarray],
    dist_angle: list[np.ndarray],
    dist_shape: list[np.ndarray],
    *,
    shape_available: bool,
) -> dict[str, list[float]]:
    out = {
        "delta": [float(matrix[0, 1]) for matrix in dist_delta],
        "angle": [float(matrix[0, 1]) for matrix in dist_angle],
    }
    if shape_available:
        out["shape"] = [float(matrix[0, 1]) for matrix in dist_shape]
    return out


def _empirical_p_value(null_values: list[float], observed: float) -> float:
    values = np.asarray(null_values, dtype=float)
    return float((1.0 + np.sum(values >= observed)) / (1.0 + values.size))


#: Null quantiles retained per statistic. ``q95`` is the alpha-level upper
#: critical value the test actually compares against; the RRPP null for an angle
#: is bounded and skewed, so it is not recoverable from the mean and sd.
_NULL_SUMMARY_QUANTILES: tuple[tuple[str, float], ...] = (
    ("q50", 0.50),
    ("q90", 0.90),
    ("q95", 0.95),
    ("q99", 0.99),
)


def _summarize_null_distributions(
    null_values: Mapping[str, list[float]],
) -> dict[str, dict[str, float]]:
    """Describe each permutation null compactly enough to persist per replicate.

    Non-finite draws are dropped and ``count`` reports how many survived, so an
    exclusion is visible rather than silently folded into the moments. Every
    value is a plain JSON-safe float: a statistic with too few retained draws to
    define a moment omits that key instead of emitting ``NaN``, which ``json``
    would write as a non-conforming literal.
    """

    summary: dict[str, dict[str, float]] = {}
    for statistic, values in null_values.items():
        draws = np.asarray(values, dtype=float)
        retained = draws[np.isfinite(draws)]
        entry: dict[str, float] = {"count": float(retained.size)}
        if retained.size >= 1:
            entry["mean"] = float(retained.mean())
            for name, q in _NULL_SUMMARY_QUANTILES:
                entry[name] = float(np.quantile(retained, q))
        if retained.size >= 2:
            entry["sd"] = float(retained.std(ddof=1))
        summary[statistic] = entry
    return summary


def _config_spectrum_block(
    spectra: Mapping[str, Any],
    group_levels: Sequence[str],
) -> dict[str, Any]:
    """Name the per-group spectra by group level and keep the block JSON-safe.

    ``configuration_spectra`` returns groups positionally, in contrast order;
    the contrast is built group-major from ``group_levels``, so the two line up
    and the persisted record can be read without knowing the contrast layout.
    """

    groups = list(spectra.get("groups") or [])
    if len(groups) != len(group_levels):
        raise SimulationEvaluationError(
            f"Spectrum block has {len(groups)} group configurations but the design carries "
            f"{len(group_levels)} groups."
        )
    return {
        "version": int(spectra["version"]),
        "pooled": spectra.get("pooled"),
        "groups": {str(level): entry for level, entry in zip(group_levels, groups)},
    }


#: Quantiles retained for the pooled eigengap over the permutation draws. The
#: observed eigengap is located against this distribution, so the tails matter.
_EIGENGAP_SUMMARY_QUANTILES: tuple[tuple[str, float], ...] = (
    ("q05", 0.05),
    ("q50", 0.50),
    ("q95", 0.95),
)


def _summarize_permutation_eigengaps(values: Sequence[float | None]) -> dict[str, float]:
    """Summarize the pooled relative eigengap across the permutation draws.

    Draws whose configuration had no defined eigengap are dropped and ``count``
    reports how many survived, matching the null-summary convention: an
    exclusion stays visible rather than being folded into the moments. Only this
    summary is retained — the per-permutation values never reach the result.
    """

    draws = np.asarray([np.nan if value is None else float(value) for value in values], dtype=float)
    retained = draws[np.isfinite(draws)] if draws.size else draws
    entry: dict[str, float] = {"count": float(retained.size)}
    if retained.size >= 1:
        entry["mean"] = float(retained.mean())
        for name, q in _EIGENGAP_SUMMARY_QUANTILES:
            entry[name] = float(np.quantile(retained, q))
    if retained.size >= 2:
        entry["sd"] = float(retained.std(ddof=1))
    return entry
