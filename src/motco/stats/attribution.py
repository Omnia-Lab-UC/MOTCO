"""Feature attribution for two-group trajectories in a shared PLS space.

The functions in this module deliberately condition on a fitted PLS estimator.
They explain an already detected orientation difference; they do not fit a
model, calculate a significance test, or replace observed features with their
low-rank reconstruction.

Two decompositions, and what each one explains
----------------------------------------------
Every analysis reports both, and they answer different questions:

- The **principal-orientation** block decomposes the *tested* ``angle``
  estimand — principal-axis divergence, the contrast between the groups' signed
  leading principal axes. It applies the very functional the statistic is built
  from (:func:`motco.stats.trajectory.principal_orientation`) to the group
  stage-mean configurations, so the quantity explained here and the quantity
  tested upstream cannot drift apart.
- The **per-adjacent-transition** blocks describe *where along the trajectory*
  the two groups' step directions diverge. This is a per-step description, not
  the tested estimand: the two coincide only for a two-stage design, where the
  configuration is rank one and PC1 is the transition direction.

Both decompositions are reported in the same observed / PLS-captured / residual
structure, and both appear in the tabular outputs under the ``transition_id``
column — the transitions under ``from->to`` identifiers, the principal block
under a reserved identifier validated against them.

Degeneracy qualifier
--------------------
A principal axis is defined for any configuration with variance, but it is only
*resolvable* when one axis dominates. Each contributing configuration's relative
eigengap is therefore reported, and the principal contrast is flagged
``degenerate`` when any of them falls below ``eigengap_threshold``. The flagged
contrast is still returned — degeneracy is a property of the data, so callers
stratify on the flag rather than lose rows. That flag is distinct from
*unavailability*: a trajectory whose net displacement vanishes is closed, and
no direction exists to report at all.

Neither the flag nor the attribution as a whole is an inference. Attribution
runs after the calling workflow has made its orientation decision.

Covariate-adjusted means
------------------------
By default each group × stage mean is the arithmetic mean of its rows. With
``covariates`` (row-aligned numeric design columns ``C``), each mean is instead
the least-squares mean of the cell-means model ``X ~ cells + C`` evaluated at a
reference covariate row ``c_ref`` (by default the column means)::

    m_gs = mean_gs(X) - gamma' (mean_gs(C) - c_ref)

with ``gamma`` the covariate coefficients estimated within cells
(:func:`motco.stats.pls.fit_nuisance_coefficients` on the cell indicators). The
reference row shifts every cell by the same vector ``gamma' c_ref``, so the
transitions, path lengths and both contrast blocks do not depend on it; only
the means tables do. The bootstrap draws the same within-cell resamples as
without covariates and refits ``gamma`` in each; a resample whose covariates
are collinear within cells is an invalid replicate, and those are counted.
Covariates enter additively and linearly, and adjustment removes nothing that
was not declared.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Any, Mapping, Sequence, cast

import numpy as np
import pandas as pd

from motco.stats.pls import NuisanceDesign, fit_nuisance_coefficients
from motco.stats.trajectory import configuration_spectrum, principal_orientation

#: Boundary statement recorded with every covariate-adjusted result.
ADJUSTMENT_BOUNDARY = (
    "Adjusted effects are group × stage least-squares means of an additive linear covariate model evaluated at "
    "the reference covariate row; bootstrap replicates refit the covariate coefficients. Adjustment removes only "
    "the declared covariates' linear contribution and does not remove unmeasured confounding."
)


class AttributionError(ValueError):
    """Raised when attribution inputs cannot satisfy the analysis contract."""


@dataclass(frozen=True)
class AttributionConfig:
    """Effective settings and coordinate-system metadata for an analysis.

    The three covariate fields are ``None`` unless the means are covariate
    adjusted; then they hold the design columns in order, the reference row
    (column → value) and the number of rank-deficient bootstrap replicates.
    """

    group_column: str
    stage_column: str
    feature_names: tuple[str, ...]
    groups: tuple[str, str]
    stages: tuple[str, ...]
    mean_source: str
    input_units: str
    original_units_available: bool
    zero_tolerance: float
    bootstrap_replicates: int
    bootstrap_seed: int | None
    top_k: int
    model_components: int
    label_namespace: str | None
    label_source: str | None
    eigengap_threshold: float
    principal_component_id: str
    covariate_columns: tuple[str, ...] | None = None
    covariate_reference: dict[str, float] | None = None
    rank_deficient_replicates: int | None = None


@dataclass(frozen=True)
class ComponentAttribution:
    """Raw and directional values for one transition/component pair."""

    component: str
    group_transitions: tuple[np.ndarray, np.ndarray]
    path_lengths: tuple[float, float]
    unit_directions: tuple[np.ndarray | None, np.ndarray | None]
    directional_contrast: np.ndarray | None
    contrast_available: bool


@dataclass(frozen=True)
class TransitionAttribution:
    """Observed, PLS-captured, and residual attribution for one transition."""

    transition_id: str
    from_stage: str
    to_stage: str
    observed: ComponentAttribution
    pls_captured: ComponentAttribution
    residual: ComponentAttribution


@dataclass(frozen=True)
class PrincipalOrientationAttribution:
    """Attribution of the tested principal-axis divergence.

    Unlike :class:`TransitionAttribution`, which describes one adjacent stage
    step, this block describes the whole stage-mean configuration: per group,
    the leading principal axis of the centered configuration signed by net
    displacement — the orientation functional the tested ``angle`` statistic is
    built from (:func:`motco.stats.trajectory.principal_orientation`).

    Each component reuses :class:`ComponentAttribution` with the configuration
    reading of its fields: ``group_transitions`` are the groups' net
    displacements (last stage minus first), ``path_lengths`` their norms,
    ``unit_directions`` the signed principal axes, and ``directional_contrast``
    the second ordered group's axis minus the first's.

    ``degenerate`` flags a near-isotropic contributing configuration — the axis
    is defined but unstable — and is distinct from an unavailable orientation,
    where a vanishing net displacement leaves no direction defined at all.
    """

    component_id: str
    from_stage: str
    to_stage: str
    observed: ComponentAttribution
    pls_captured: ComponentAttribution
    residual: ComponentAttribution
    observed_eigengaps: tuple[float | None, float | None]
    reconstructed_eigengaps: tuple[float | None, float | None]
    eigengap_threshold: float
    degenerate: bool
    degeneracy_sources: tuple[str, ...]
    orientation_available: tuple[bool, bool]


@dataclass(frozen=True)
class FeatureEffect:
    """One signed feature contribution in a transition/component."""

    transition_id: str
    component: str
    feature: str
    effect_standardized: float
    effect_original: float | None


@dataclass(frozen=True)
class AggregateEffect:
    """One caller-label aggregate in a transition/component."""

    transition_id: str
    component: str
    label: str
    signed_effect: float
    absolute_effect: float
    feature_count: int
    label_namespace: str


@dataclass(frozen=True)
class BootstrapSummary:
    """Per-feature stability summary for one transition/component."""

    transition_id: str
    component: str
    requested_replicates: int
    valid_replicates: int
    top_k: int
    seed: int | None
    sign_stability: np.ndarray
    rank_stability: np.ndarray
    top_k_selection_frequency: np.ndarray
    nonzero_sign_replicates: np.ndarray


@dataclass(frozen=True)
class InterpretationMetadata:
    """Interpretation boundary recorded with every result.

    ``adjustment_boundary`` is set only for covariate-adjusted means.
    """

    statement: str
    causal_boundary: str
    significance_boundary: str
    reconstruction_boundary: str
    bootstrap_boundary: str
    label_source: str | None
    decomposition_boundary: str
    adjustment_boundary: str | None = None


@dataclass(frozen=True)
class OrientationAttributionResult:
    """Structured output of :func:`analyze_orientation_attribution`."""

    config: AttributionConfig
    means: pd.DataFrame
    reconstructed_means: pd.DataFrame
    transitions: tuple[TransitionAttribution, ...]
    principal_orientation: PrincipalOrientationAttribution
    feature_effects: pd.DataFrame
    transition_summaries: pd.DataFrame
    transition_vectors: pd.DataFrame
    principal_summaries: pd.DataFrame
    principal_vectors: pd.DataFrame
    aggregate_effects: pd.DataFrame
    bootstrap_summaries: pd.DataFrame
    interpretation: InterpretationMetadata
    feature_records: tuple[FeatureEffect, ...]
    aggregate_records: tuple[AggregateEffect, ...]
    bootstrap_records: tuple[BootstrapSummary, ...]

    @property
    def groups(self) -> tuple[str, str]:
        """Ordered group labels."""

        return self.config.groups

    @property
    def stages(self) -> tuple[str, ...]:
        """Ordered stage labels."""

        return self.config.stages

    @property
    def feature_names(self) -> tuple[str, ...]:
        """Canonical feature order used by the model and result tables."""

        return self.config.feature_names

    def frames(self) -> dict[str, pd.DataFrame]:
        """Return findings-ready tabular views without prose parsing."""

        return attribution_frames(self)


AttributionResult = OrientationAttributionResult


def analyze_orientation_attribution(
    X: pd.DataFrame,
    metadata: pd.DataFrame,
    model: Any,
    *,
    group_col: str = "group",
    stage_col: str = "stage",
    feature_cols: Sequence[str] | None = None,
    groups: Sequence[str] | None = None,
    stages: Sequence[str] | None = None,
    mean_table: pd.DataFrame | None = None,
    precomputed_means: pd.DataFrame | None = None,
    covariates: pd.DataFrame | np.ndarray | None = None,
    covariate_reference: Sequence[float] | pd.Series | Mapping[str, float] | None = None,
    original_scales: Sequence[float] | pd.Series | Mapping[str, float] | None = None,
    feature_scales: Sequence[float] | pd.Series | Mapping[str, float] | None = None,
    feature_groups: Mapping[str, str] | pd.Series | pd.DataFrame | None = None,
    group_labels: Mapping[str, str] | pd.Series | pd.DataFrame | None = None,
    label_namespace: str = "module",
    bootstrap_replicates: int = 0,
    bootstrap_seed: int | None = 0,
    top_k: int = 10,
    zero_tolerance: float = 1e-12,
    eigengap_threshold: float = 0.05,
    principal_component_id: str = "principal",
    sample_id_col: str | None = None,
) -> OrientationAttributionResult:
    """Attribute two-group trajectory orientation differences.

    Reports both decompositions described in the module docstring: the
    principal-orientation contrast, which decomposes the tested ``angle``
    estimand and carries a relative-eigengap degeneracy qualifier, and the
    per-adjacent-transition contrasts, which describe where along the trajectory
    the group directions diverge. The two coincide only for a two-stage design.

    Parameters
    ----------
    X:
        Aligned feature matrix in the exact coordinate system consumed by
        ``model``. Observed values are never overwritten by reconstructions.
    metadata:
        Row-aligned metadata containing ``group_col`` and ``stage_col``.
    model:
        One already fitted pooled PLS estimator exposing ``transform`` and
        ``inverse_transform``. It is reused for all groups and replicates.
    mean_table, precomputed_means:
        Optional complete group-by-stage table. It may have a two-level
        MultiIndex or explicit group/stage columns. When supplied, bootstrap
        resampling is unavailable because row-level sampling semantics are not
        defined for precomputed means. Exclusive with ``covariates``.
    covariates:
        Optional numeric covariate design columns, row-aligned with ``X`` (a
        DataFrame must share its index; an array's columns are named
        ``covariate_0`` …). Every cell mean is then the least-squares mean of
        ``X ~ cells + covariates`` at ``covariate_reference`` (see the module
        docstring), and each bootstrap replicate refits the covariate
        coefficients on its resampled rows. Categorical covariates are the
        caller's to code as indicator columns. Covariates collinear within the
        group × stage cells are an error.
    covariate_reference:
        The covariate row the adjusted means are evaluated at: a mapping or
        Series over the covariate columns, or a sequence in their order. The
        default is the column means. It moves the means tables only.
    original_scales, feature_scales:
        Optional positive per-feature scale in canonical feature order. The
        aliases are accepted for readability at call sites; supplying both is
        an error.
    feature_groups, group_labels:
        Optional one-to-one caller-supplied feature-to-label mapping. The
        aliases are accepted for readability at call sites.
    eigengap_threshold:
        Relative-eigengap level at or above which a stage-mean configuration is
        treated as having a well-separated principal axis. A contributing
        configuration below it flags the principal-orientation contrast as
        degenerate; the contrast is still returned, so callers stratify rather
        than lose rows. The default (0.05) separates the audit's resolvable
        replicates from its near-isotropic independent-baseline draws.
    principal_component_id:
        Reserved identifier carrying the principal-orientation component in the
        tabular outputs' ``transition_id`` column. It is validated against the
        stage-derived transition ids so the two can never collide.

    Returns
    -------
    OrientationAttributionResult
        Frozen typed records and DataFrame views for means, transitions, the
        principal-orientation block, feature effects, aggregates, and bootstrap
        stability.
    """
    if mean_table is not None and precomputed_means is not None:
        raise AttributionError("Provide only one of mean_table and precomputed_means.")
    if covariates is not None and (mean_table is not None or precomputed_means is not None):
        raise AttributionError(
            "Provide only one mean source: a precomputed mean table, or covariates for adjusted means."
        )
    if original_scales is not None and feature_scales is not None:
        raise AttributionError("Provide only one of original_scales and feature_scales.")
    if feature_groups is not None and group_labels is not None:
        raise AttributionError("Provide only one of feature_groups and group_labels.")
    if not isinstance(label_namespace, str) or not label_namespace.strip():
        raise AttributionError("label_namespace must be a non-empty string.")
    if precomputed_means is not None:
        mean_table = precomputed_means
    if feature_scales is not None:
        original_scales = feature_scales
    if group_labels is not None:
        feature_groups = group_labels

    feature_frame, metadata_frame, feature_names = _validate_inputs(
        X,
        metadata,
        model,
        group_col,
        stage_col,
        feature_cols,
        sample_id_col,
    )
    ordered_groups = cast(
        tuple[str, str], _ordered_levels(metadata_frame[group_col], groups, 2, "groups")
    )
    ordered_stages = _ordered_levels(metadata_frame[stage_col], stages, 2, "stages")
    _validate_cells(metadata_frame, group_col, stage_col, ordered_groups, ordered_stages)
    scale_values = _validate_scales(original_scales, feature_names)
    label_values = _validate_feature_groups(feature_groups, feature_names)
    adjustment = _validate_covariates(covariates, covariate_reference, feature_frame.index)

    if adjustment is not None:
        cells = _cell_codes(metadata_frame, group_col, stage_col, ordered_groups, ordered_stages)
        try:
            means_array = _adjusted_means(
                feature_frame.to_numpy(dtype=float),
                adjustment.matrix,
                cells,
                adjustment,
                ordered_groups,
                ordered_stages,
            )
        except ValueError as exc:
            raise AttributionError(
                f"Covariates {list(adjustment.names)} are collinear with the group × stage cells or with one "
                f"another, so covariate-adjusted cell means are not identified: {exc}"
            ) from exc
        mean_source = "covariate_adjusted"
    elif mean_table is None:
        means_array = _arithmetic_means(
            feature_frame,
            metadata_frame,
            group_col,
            stage_col,
            ordered_groups,
            ordered_stages,
        )
        mean_source = "arithmetic"
    else:
        means_array = _precomputed_mean_array(
            mean_table,
            feature_names,
            ordered_groups,
            ordered_stages,
            group_col,
            stage_col,
        )
        mean_source = "precomputed"
        if bootstrap_replicates:
            raise AttributionError(
                "bootstrap_replicates must be zero when mean_table contains precomputed means; "
                "row-level bootstrap semantics are unavailable."
            )

    _validate_bootstrap_settings(bootstrap_replicates, top_k, zero_tolerance)
    _validate_principal_settings(eigengap_threshold, principal_component_id, ordered_stages)
    reconstructed_array = _reconstruct_means(means_array, model, feature_names)
    transitions = _build_transitions(
        means_array,
        reconstructed_array,
        ordered_groups,
        ordered_stages,
        zero_tolerance,
    )
    principal = _build_principal_orientation(
        means_array,
        reconstructed_array,
        ordered_groups,
        ordered_stages,
        zero_tolerance,
        eigengap_threshold,
        principal_component_id,
    )
    model_components = _model_components(model)
    bootstrap_records, rank_deficient = _bootstrap(
        feature_frame,
        metadata_frame,
        group_col,
        stage_col,
        ordered_groups,
        ordered_stages,
        model,
        feature_names,
        scale_values,
        bootstrap_replicates,
        bootstrap_seed,
        top_k,
        zero_tolerance,
        principal_component_id,
        adjustment,
    )
    config = AttributionConfig(
        group_column=group_col,
        stage_column=stage_col,
        feature_names=feature_names,
        groups=ordered_groups,
        stages=ordered_stages,
        mean_source=mean_source,
        input_units="standardized_input",
        original_units_available=scale_values is not None,
        zero_tolerance=float(zero_tolerance),
        bootstrap_replicates=int(bootstrap_replicates),
        bootstrap_seed=bootstrap_seed,
        top_k=int(top_k),
        model_components=model_components,
        label_namespace=label_namespace if label_values is not None else None,
        label_source="caller-supplied" if label_values is not None else None,
        eigengap_threshold=float(eigengap_threshold),
        principal_component_id=principal_component_id,
        covariate_columns=adjustment.names if adjustment is not None else None,
        covariate_reference=(
            dict(zip(adjustment.names, map(float, adjustment.reference))) if adjustment is not None else None
        ),
        rank_deficient_replicates=rank_deficient if adjustment is not None else None,
    )
    bootstrap_lookup = {
        (record.transition_id, record.component): record for record in bootstrap_records
    }
    feature_effects, feature_records = _feature_effect_table(
        transitions,
        principal,
        feature_names,
        scale_values,
        bootstrap_lookup,
        label_values,
        label_namespace,
        bootstrap_replicates,
    )
    aggregate_effects, aggregate_records = _aggregate_effect_table(
        feature_effects,
        feature_names,
        label_values,
        label_namespace,
        bootstrap_replicates,
    )
    transition_summaries, transition_vectors = _transition_tables(
        transitions,
        ordered_groups,
        feature_names,
    )
    principal_summaries, principal_vectors = _principal_tables(
        principal,
        ordered_groups,
        feature_names,
    )
    bootstrap_table = _bootstrap_table(bootstrap_records, feature_names)
    means_table = _means_table(means_array, ordered_groups, ordered_stages, feature_names)
    reconstructed_table = _means_table(
        reconstructed_array, ordered_groups, ordered_stages, feature_names
    )
    interpretation = InterpretationMetadata(
        statement="Attribution describes associations implied by one fitted shared PLS representation.",
        causal_boundary="Feature, module, and pathway effects are not causal claims.",
        significance_boundary=(
            "No p-values or significance decisions are produced; external orientation inference remains the caller's "
            "responsibility."
        ),
        reconstruction_boundary=(
            "PLS-captured effects are a lossy inverse-transform view and are not a unique inverse of measured features."
        ),
        bootstrap_boundary=(
            "Bootstrap stability is conditional on the fitted preprocessing contract and frozen PLS model."
        ),
        label_source="caller-supplied" if label_values is not None else None,
        decomposition_boundary=(
            f"The '{principal_component_id}' component decomposes the tested angle estimand — principal-axis "
            "divergence between the groups' signed leading principal axes — while per-adjacent-transition "
            "contrasts describe where along the trajectory the group directions diverge; the two coincide only "
            "for two-stage designs. A degenerate principal contrast is reported with its flag and MUST NOT be "
            "read as a well-defined orientation."
        ),
        adjustment_boundary=ADJUSTMENT_BOUNDARY if adjustment is not None else None,
    )
    return OrientationAttributionResult(
        config=config,
        means=means_table,
        reconstructed_means=reconstructed_table,
        transitions=tuple(transitions),
        principal_orientation=principal,
        feature_effects=feature_effects,
        transition_summaries=transition_summaries,
        transition_vectors=transition_vectors,
        principal_summaries=principal_summaries,
        principal_vectors=principal_vectors,
        aggregate_effects=aggregate_effects,
        bootstrap_summaries=bootstrap_table,
        interpretation=interpretation,
        feature_records=tuple(feature_records),
        aggregate_records=tuple(aggregate_records),
        bootstrap_records=tuple(bootstrap_records),
    )


run_orientation_attribution = analyze_orientation_attribution


def attribution_frames(result: OrientationAttributionResult) -> dict[str, pd.DataFrame]:
    """Return all findings-ready tables for ``result``.

    The covariate rows of ``configuration`` and the ``adjustment_boundary`` row
    of ``interpretation`` appear only for covariate-adjusted means.
    """

    adjustment_settings: list[tuple[str, Any]] = []
    if result.config.covariate_columns is not None:
        adjustment_settings = [
            ("covariate_columns", list(result.config.covariate_columns)),
            ("covariate_reference", dict(result.config.covariate_reference or {})),
            ("bootstrap_rank_deficient_replicates", result.config.rank_deficient_replicates),
        ]
    adjustment_fields: list[tuple[str, Any]] = []
    if result.interpretation.adjustment_boundary is not None:
        adjustment_fields = [("adjustment_boundary", result.interpretation.adjustment_boundary)]
    configuration = pd.DataFrame(
        [
            ("group_column", result.config.group_column),
            ("stage_column", result.config.stage_column),
            ("feature_names", list(result.config.feature_names)),
            ("groups", list(result.config.groups)),
            ("stages", list(result.config.stages)),
            ("mean_source", result.config.mean_source),
            ("input_units", result.config.input_units),
            ("original_units_available", result.config.original_units_available),
            ("zero_tolerance", result.config.zero_tolerance),
            ("bootstrap_replicates", result.config.bootstrap_replicates),
            ("bootstrap_seed", result.config.bootstrap_seed),
            ("top_k", result.config.top_k),
            ("model_components", result.config.model_components),
            ("label_namespace", result.config.label_namespace),
            ("label_source", result.config.label_source),
            ("eigengap_threshold", result.config.eigengap_threshold),
            ("principal_component_id", result.config.principal_component_id),
            *adjustment_settings,
        ],
        columns=["setting", "value"],
    )
    interpretation = pd.DataFrame(
        [
            ("statement", result.interpretation.statement),
            ("causal_boundary", result.interpretation.causal_boundary),
            ("significance_boundary", result.interpretation.significance_boundary),
            ("reconstruction_boundary", result.interpretation.reconstruction_boundary),
            ("bootstrap_boundary", result.interpretation.bootstrap_boundary),
            ("label_source", result.interpretation.label_source),
            ("decomposition_boundary", result.interpretation.decomposition_boundary),
            *adjustment_fields,
        ],
        columns=["field", "value"],
    )
    return {
        "means": result.means.copy(),
        "reconstructed_means": result.reconstructed_means.copy(),
        "feature_effects": result.feature_effects.copy(),
        "transition_summaries": result.transition_summaries.copy(),
        "transition_vectors": result.transition_vectors.copy(),
        "principal_orientation": result.principal_summaries.copy(),
        "principal_orientation_vectors": result.principal_vectors.copy(),
        "aggregate_effects": result.aggregate_effects.copy(),
        "bootstrap_summaries": result.bootstrap_summaries.copy(),
        "configuration": configuration,
        "interpretation": interpretation,
    }


def write_attribution_outputs(
    result: OrientationAttributionResult,
    output_dir: str | Path,
) -> dict[str, Path]:
    """Write attribution tables as CSV files and return their paths."""

    directory = Path(output_dir)
    directory.mkdir(parents=True, exist_ok=True)
    paths: dict[str, Path] = {}
    for name, frame in attribution_frames(result).items():
        path = directory / f"{name}.csv"
        frame.to_csv(path, index=False)
        paths[name] = path
    return paths


def _validate_inputs(
    X: pd.DataFrame,
    metadata: pd.DataFrame,
    model: Any,
    group_col: str,
    stage_col: str,
    feature_cols: Sequence[str] | None,
    sample_id_col: str | None,
) -> tuple[pd.DataFrame, pd.DataFrame, tuple[str, ...]]:
    if not isinstance(X, pd.DataFrame):
        raise AttributionError("X must be a pandas DataFrame with canonical feature columns.")
    if not isinstance(metadata, pd.DataFrame):
        raise AttributionError("metadata must be a pandas DataFrame aligned to X rows.")
    if not X.index.equals(metadata.index):
        raise AttributionError("X and metadata rows or sample identifiers are not aligned.")
    if X.index.has_duplicates or metadata.index.has_duplicates:
        raise AttributionError("X and metadata sample identifiers must be unique.")
    if group_col not in metadata.columns or stage_col not in metadata.columns:
        raise AttributionError(
            f"metadata must contain group column '{group_col}' and stage column '{stage_col}'."
        )
    if sample_id_col is not None:
        if sample_id_col not in metadata.columns:
            raise AttributionError(f"sample_id_col='{sample_id_col}' is not present in metadata.")
        sample_ids = metadata[sample_id_col]
        if sample_ids.isna().any() or sample_ids.duplicated().any():
            raise AttributionError("sample identifiers must be non-missing and unique.")
        if not pd.Index(sample_ids.tolist()).equals(X.index):
            raise AttributionError("sample identifiers do not match the feature matrix index.")

    names = tuple(map(str, feature_cols if feature_cols is not None else X.columns.tolist()))
    if len(names) == 0:
        raise AttributionError("At least one feature column is required.")
    if len(set(names)) != len(names):
        raise AttributionError("feature names must be unique.")
    actual_names = tuple(map(str, X.columns.tolist()))
    if actual_names != names:
        missing = [name for name in names if name not in actual_names]
        if missing:
            raise AttributionError(f"feature_cols contains unknown feature(s): {missing}.")
        positions = [actual_names.index(name) for name in names]
        X = X.iloc[:, positions]
    try:
        feature_frame = X.copy()
        feature_frame.columns = names
        feature_frame = feature_frame.astype(float)
    except (TypeError, ValueError) as exc:
        raise AttributionError("X feature columns must be numeric.") from exc
    if not np.isfinite(feature_frame.to_numpy()).all():
        raise AttributionError("X contains NaN or Inf values.")

    n_features = len(names)
    model_features = getattr(model, "n_features_in_", None)
    if model_features is not None and int(model_features) != n_features:
        raise AttributionError(
            f"fitted model expects {int(model_features)} features but X has {n_features}."
        )
    model_names = getattr(model, "feature_names_in_", None)
    if model_names is not None and tuple(map(str, model_names)) != names:
        raise AttributionError("X feature order does not match the fitted model feature order.")
    if not callable(getattr(model, "transform", None)) or not callable(
        getattr(model, "inverse_transform", None)
    ):
        raise AttributionError("fitted model must expose transform and inverse_transform methods.")
    return feature_frame, metadata.copy(), names


def _ordered_levels(
    values: pd.Series,
    explicit: Sequence[str] | None,
    minimum: int,
    name: str,
) -> tuple[str, ...]:
    if values.isna().any():
        raise AttributionError(f"{name} contain missing values.")
    observed = {str(value) for value in pd.unique(values)}
    if explicit is None:
        ordered = tuple(sorted(observed))
    else:
        ordered = tuple(map(str, explicit))
        if len(set(ordered)) != len(ordered):
            raise AttributionError(f"explicit {name} must not contain duplicates.")
        if set(ordered) != observed:
            raise AttributionError(
                f"explicit {name} {ordered} do not match observed {name} {tuple(sorted(observed))}."
            )
    if len(ordered) < minimum:
        raise AttributionError(f"expected at least {minimum} {name}; found {len(ordered)}: {ordered}.")
    if name == "groups" and len(ordered) != 2:
        raise AttributionError(f"expected exactly 2 groups; found {len(ordered)}: {ordered}.")
    return ordered


def _validate_cells(
    metadata: pd.DataFrame,
    group_col: str,
    stage_col: str,
    groups: tuple[str, str],
    stages: tuple[str, ...],
) -> None:
    group_values = metadata[group_col].astype(str)
    stage_values = metadata[stage_col].astype(str)
    missing = [
        (group, stage)
        for group in groups
        for stage in stages
        if not ((group_values == group) & (stage_values == stage)).any()
    ]
    if missing:
        raise AttributionError(f"missing complete group-by-stage cells: {missing}.")


def _validate_scales(
    scales: Sequence[float] | pd.Series | Mapping[str, float] | None,
    feature_names: tuple[str, ...],
) -> np.ndarray | None:
    if scales is None:
        return None
    if isinstance(scales, pd.Series):
        if not scales.index.is_unique or tuple(map(str, scales.index)) != feature_names:
            raise AttributionError("feature scales must have unique index in the input feature order.")
        values = scales.to_numpy(dtype=float)
    elif isinstance(scales, Mapping):
        keys = tuple(map(str, scales.keys()))
        if keys != feature_names:
            raise AttributionError("feature scale mapping keys must match the input feature order.")
        values = np.asarray([scales[name] for name in feature_names], dtype=float)
    else:
        values = np.asarray(scales, dtype=float)
        if values.ndim != 1 or len(values) != len(feature_names):
            raise AttributionError("feature scales must contain exactly one value per input feature.")
    if not np.isfinite(values).all() or np.any(values <= 0):
        raise AttributionError("feature scales must be finite and strictly positive.")
    return values


def _validate_feature_groups(
    labels: Mapping[str, str] | pd.Series | pd.DataFrame | None,
    feature_names: tuple[str, ...],
) -> dict[str, str] | None:
    if labels is None:
        return None
    pairs: list[tuple[Any, Any]]
    if isinstance(labels, pd.Series):
        pairs = list(zip(labels.index.tolist(), labels.tolist()))
    elif isinstance(labels, Mapping):
        pairs = list(labels.items())
    elif isinstance(labels, pd.DataFrame):
        if labels.shape[1] != 2:
            raise AttributionError("feature-group DataFrame must contain exactly feature and label columns.")
        pairs = list(zip(labels.iloc[:, 0].tolist(), labels.iloc[:, 1].tolist()))
    else:
        raise AttributionError("feature groups must be a mapping, Series, or two-column DataFrame.")
    result: dict[str, str] = {}
    for feature, label in pairs:
        feature_name = str(feature)
        if feature_name not in feature_names:
            raise AttributionError(f"feature-group mapping contains unknown feature '{feature_name}'.")
        if pd.isna(label) or str(label).strip() == "":
            raise AttributionError(f"feature-group mapping has an empty label for '{feature_name}'.")
        label_name = str(label)
        if feature_name in result and result[feature_name] != label_name:
            raise AttributionError(f"feature '{feature_name}' has conflicting group labels.")
        result[feature_name] = label_name
    missing = [feature for feature in feature_names if feature not in result]
    if missing:
        raise AttributionError(f"feature-group mapping omits input feature(s): {missing}.")
    return result


def _arithmetic_means(
    X: pd.DataFrame,
    metadata: pd.DataFrame,
    group_col: str,
    stage_col: str,
    groups: tuple[str, str],
    stages: tuple[str, ...],
) -> np.ndarray:
    group_values = metadata[group_col].astype(str)
    stage_values = metadata[stage_col].astype(str)
    values = X.to_numpy(dtype=float)
    result = np.empty((len(groups), len(stages), X.shape[1]), dtype=float)
    for gi, group in enumerate(groups):
        for si, stage in enumerate(stages):
            mask = ((group_values == group) & (stage_values == stage)).to_numpy()
            result[gi, si] = values[mask].mean(axis=0)
    return result


@dataclass(frozen=True)
class _CovariateAdjustment:
    """Validated covariate design columns and the reference row the adjusted means use."""

    matrix: np.ndarray
    names: tuple[str, ...]
    reference: np.ndarray


def _validate_covariates(
    covariates: pd.DataFrame | np.ndarray | None,
    reference: Sequence[float] | pd.Series | Mapping[str, float] | None,
    index: pd.Index,
) -> _CovariateAdjustment | None:
    if covariates is None:
        if reference is not None:
            raise AttributionError("covariate_reference requires covariates.")
        return None
    if isinstance(covariates, pd.DataFrame):
        if covariates.shape[0] != len(index):
            raise AttributionError(f"covariates have {covariates.shape[0]} rows but X has {len(index)}.")
        if not covariates.index.equals(index):
            raise AttributionError("covariate rows are not aligned to X rows; their indexes must be equal.")
        names = tuple(map(str, covariates.columns))
        raw: Any = covariates
    else:
        raw = np.asarray(covariates)
        if raw.ndim != 2:
            raise AttributionError("covariates must be a two-dimensional table of design columns.")
        if raw.shape[0] != len(index):
            raise AttributionError(f"covariates have {raw.shape[0]} rows but X has {len(index)}.")
        names = tuple(f"covariate_{column}" for column in range(raw.shape[1]))
    try:
        matrix = np.asarray(raw, dtype=float)
    except (TypeError, ValueError) as exc:
        raise AttributionError(
            "covariates must be numeric design columns; code categorical covariates as indicator columns."
        ) from exc
    if matrix.shape[1] == 0:
        raise AttributionError("covariates must contain at least one design column.")
    if len(set(names)) != len(names):
        raise AttributionError(f"covariate column names must be unique; got {list(names)}.")
    if not np.isfinite(matrix).all():
        raise AttributionError("covariates contain missing or non-finite values.")
    return _CovariateAdjustment(matrix, names, _validate_covariate_reference(reference, names, matrix))


def _validate_covariate_reference(
    reference: Sequence[float] | pd.Series | Mapping[str, float] | None,
    names: tuple[str, ...],
    matrix: np.ndarray,
) -> np.ndarray:
    if reference is None:
        return matrix.mean(axis=0)
    try:
        if isinstance(reference, (pd.Series, Mapping)):
            keys = [str(key) for key in reference.keys()]
            if len(keys) != len(names) or set(keys) != set(names):
                raise AttributionError(
                    f"covariate_reference names {keys} but the covariate columns are {list(names)}."
                )
            lookup = {str(key): value for key, value in reference.items()}
            values = np.asarray([lookup[name] for name in names], dtype=float)
        else:
            values = np.asarray(reference, dtype=float)
            if values.ndim != 1 or values.shape[0] != len(names):
                raise AttributionError(
                    f"covariate_reference must hold one value per covariate column ({len(names)})."
                )
    except (TypeError, ValueError) as exc:
        if isinstance(exc, AttributionError):
            raise
        raise AttributionError("covariate_reference values must be numeric.") from exc
    if not np.isfinite(values).all():
        raise AttributionError("covariate_reference contains a non-finite value.")
    return values


def _cell_codes(
    metadata: pd.DataFrame,
    group_col: str,
    stage_col: str,
    groups: tuple[str, str],
    stages: tuple[str, ...],
) -> np.ndarray:
    """Each row's group × stage cell, coded group-major, stage-minor (the order of every mean array)."""

    group_values = metadata[group_col].astype(str).to_numpy()
    stage_values = metadata[stage_col].astype(str).to_numpy()
    codes = np.empty(len(metadata), dtype=int)
    for gi, group in enumerate(groups):
        for si, stage in enumerate(stages):
            codes[(group_values == group) & (stage_values == stage)] = gi * len(stages) + si
    return codes


def _adjusted_means(
    values: np.ndarray,
    covariates: np.ndarray,
    cells: np.ndarray,
    adjustment: _CovariateAdjustment,
    groups: tuple[str, str],
    stages: tuple[str, ...],
) -> np.ndarray:
    """Least-squares cell means of ``values ~ cells + covariates`` at the reference covariate row.

    The covariate coefficients are estimated within cells (drop-first cell
    indicators as the conditioning design), so each mean is the cell's
    arithmetic mean moved along them from the cell's covariate mean to the
    reference row. Raises ``ValueError`` when the covariates are collinear
    within cells.
    """

    labels = [f"{group}|{stage}" for group in groups for stage in stages]
    conditioning = (cells[:, None] == np.arange(1, len(labels))[None, :]).astype(float)
    design = NuisanceDesign(conditioning, covariates, tuple(labels[1:]), adjustment.names)
    gamma, _ = fit_nuisance_coefficients(values, design)
    result = np.empty((len(groups), len(stages), values.shape[1]), dtype=float)
    for code in range(len(labels)):
        rows = cells == code
        shift = (covariates[rows].mean(axis=0) - adjustment.reference) @ gamma
        result[code // len(stages), code % len(stages)] = values[rows].mean(axis=0) - shift
    return result


def _precomputed_mean_array(
    table: pd.DataFrame,
    feature_names: tuple[str, ...],
    groups: tuple[str, str],
    stages: tuple[str, ...],
    group_col: str,
    stage_col: str,
) -> np.ndarray:
    if not isinstance(table, pd.DataFrame):
        raise AttributionError("mean_table must be a pandas DataFrame.")
    frame = table.copy()
    if group_col in frame.columns and stage_col in frame.columns:
        group_values = frame[group_col].astype(str)
        stage_values = frame[stage_col].astype(str)
        table_features = [column for column in frame.columns if column not in {group_col, stage_col}]
        if tuple(map(str, table_features)) != feature_names:
            raise AttributionError("mean_table feature columns must match the input feature order exactly.")
        values = frame.loc[:, table_features].to_numpy(dtype=float)
    elif isinstance(frame.index, pd.MultiIndex) and frame.index.nlevels == 2:
        group_values = pd.Series(frame.index.get_level_values(0).astype(str), index=frame.index)
        stage_values = pd.Series(frame.index.get_level_values(1).astype(str), index=frame.index)
        if tuple(map(str, frame.columns)) != feature_names:
            raise AttributionError("mean_table feature columns must match the input feature order exactly.")
        values = frame.loc[:, list(feature_names)].to_numpy(dtype=float)
    else:
        raise AttributionError(
            f"mean_table must have '{group_col}'/'{stage_col}' columns or a two-level MultiIndex."
        )
    if not np.isfinite(values).all():
        raise AttributionError("mean_table contains NaN or Inf values.")
    result = np.empty((len(groups), len(stages), len(feature_names)), dtype=float)
    for gi, group in enumerate(groups):
        for si, stage in enumerate(stages):
            mask = ((group_values == group) & (stage_values == stage)).to_numpy()
            if int(mask.sum()) != 1:
                raise AttributionError(
                    f"mean_table must contain exactly one row for group-stage cell ({group}, {stage})."
                )
            result[gi, si] = values[mask][0]
    if set(zip(group_values, stage_values)) != {(g, s) for g in groups for s in stages}:
        raise AttributionError("mean_table contains unknown or incomplete group-by-stage cells.")
    return result


def _reconstruct_means(means: np.ndarray, model: Any, feature_names: tuple[str, ...]) -> np.ndarray:
    flat = means.reshape((-1, len(feature_names)))
    try:
        scores = np.asarray(model.transform(flat), dtype=float)
        reconstructed = np.asarray(model.inverse_transform(scores), dtype=float)
    except Exception as exc:  # sklearn estimators expose varied exception types
        raise AttributionError(
            "fitted model could not transform and inverse-transform the aligned feature means."
        ) from exc
    if reconstructed.shape != flat.shape or not np.isfinite(reconstructed).all():
        raise AttributionError("fitted model returned invalid inverse-transformed mean dimensions or values.")
    return reconstructed.reshape(means.shape)


def _model_components(model: Any) -> int:
    value = getattr(model, "n_components_", getattr(model, "n_components", None))
    if value is None:
        return 0
    try:
        return int(value)
    except (TypeError, ValueError):
        return 0


def _component(
    name: str,
    transitions: tuple[np.ndarray, np.ndarray],
    tolerance: float,
    contrast: np.ndarray | None = None,
    derive_contrast: bool = True,
) -> ComponentAttribution:
    lengths: tuple[float, float] = (
        float(np.linalg.norm(transitions[0])),
        float(np.linalg.norm(transitions[1])),
    )
    directions: tuple[np.ndarray | None, np.ndarray | None] = (
        transitions[0] / lengths[0] if lengths[0] > tolerance else None,
        transitions[1] / lengths[1] if lengths[1] > tolerance else None,
    )
    directions_available = directions[0] is not None and directions[1] is not None
    if contrast is None and derive_contrast and directions_available:
        assert directions[0] is not None and directions[1] is not None
        contrast = directions[1] - directions[0]
    available = contrast is not None
    return ComponentAttribution(
        component=name,
        group_transitions=transitions,
        path_lengths=lengths,
        unit_directions=directions,
        directional_contrast=contrast,
        contrast_available=available,
    )


def _build_transitions(
    observed: np.ndarray,
    reconstructed: np.ndarray,
    groups: tuple[str, str],
    stages: tuple[str, ...],
    tolerance: float,
) -> list[TransitionAttribution]:
    del groups  # The ordered group contract is encoded by tuple position.
    records: list[TransitionAttribution] = []
    for stage_index in range(len(stages) - 1):
        from_stage, to_stage = stages[stage_index : stage_index + 2]
        transition_id = f"{from_stage}->{to_stage}"
        observed_d = tuple(observed[group, stage_index + 1] - observed[group, stage_index] for group in range(2))
        captured_d = tuple(
            reconstructed[group, stage_index + 1] - reconstructed[group, stage_index] for group in range(2)
        )
        residual_d = tuple(observed_d[index] - captured_d[index] for index in range(2))
        observed_component = _component("observed", observed_d, tolerance)
        captured_component = _component("pls_captured", captured_d, tolerance)
        residual_contrast = None
        if observed_component.directional_contrast is not None and captured_component.directional_contrast is not None:
            residual_contrast = observed_component.directional_contrast - captured_component.directional_contrast
        residual_component = _component(
            "residual", residual_d, tolerance, residual_contrast, derive_contrast=False
        )
        records.append(
            TransitionAttribution(
                transition_id=transition_id,
                from_stage=from_stage,
                to_stage=to_stage,
                observed=observed_component,
                pls_captured=captured_component,
                residual=residual_component,
            )
        )
    return records


def _principal_component(
    name: str,
    configurations: tuple[np.ndarray, np.ndarray],
    tolerance: float,
    contrast: np.ndarray | None = None,
    derive_contrast: bool = True,
) -> ComponentAttribution:
    """One component of the principal-orientation block.

    ``configurations`` are the two groups' ``k x p`` stage-mean configurations
    for this component. Each group's net displacement (last stage minus first)
    both anchors the principal axis's sign and decides availability: at or below
    ``tolerance`` the trajectory is closed and no direction is defined.
    """

    displacements = tuple(configuration[-1] - configuration[0] for configuration in configurations)
    lengths: tuple[float, float] = (
        float(np.linalg.norm(displacements[0])),
        float(np.linalg.norm(displacements[1])),
    )
    axes: tuple[np.ndarray | None, np.ndarray | None] = (
        principal_orientation(configurations[0]) if lengths[0] > tolerance else None,
        principal_orientation(configurations[1]) if lengths[1] > tolerance else None,
    )
    axes_available = axes[0] is not None and axes[1] is not None
    if contrast is None and derive_contrast and axes_available:
        assert axes[0] is not None and axes[1] is not None
        contrast = axes[1] - axes[0]
    return ComponentAttribution(
        component=name,
        group_transitions=cast("tuple[np.ndarray, np.ndarray]", displacements),
        path_lengths=lengths,
        unit_directions=axes,
        directional_contrast=contrast,
        contrast_available=contrast is not None,
    )


def _relative_eigengap(configuration: np.ndarray) -> float | None:
    value = configuration_spectrum(configuration)["relative_eigengap"]
    return None if value is None else float(value)


def _build_principal_orientation(
    observed: np.ndarray,
    reconstructed: np.ndarray,
    groups: tuple[str, str],
    stages: tuple[str, ...],
    tolerance: float,
    eigengap_threshold: float,
    component_id: str,
) -> PrincipalOrientationAttribution:
    """Principal-axis divergence of the two groups' stage-mean configurations.

    The observed and PLS-reconstructed configurations each contribute a signed
    principal axis per group; the residual component applies the same functional
    to the residual configuration, with its contrast defined — exactly as for
    transitions — as observed minus captured.
    """

    observed_configs = (observed[0], observed[1])
    captured_configs = (reconstructed[0], reconstructed[1])
    residual_configs = (observed[0] - reconstructed[0], observed[1] - reconstructed[1])

    observed_component = _principal_component("observed", observed_configs, tolerance)
    captured_component = _principal_component("pls_captured", captured_configs, tolerance)
    residual_contrast = None
    if (
        observed_component.directional_contrast is not None
        and captured_component.directional_contrast is not None
    ):
        residual_contrast = (
            observed_component.directional_contrast - captured_component.directional_contrast
        )
    residual_component = _principal_component(
        "residual", residual_configs, tolerance, residual_contrast, derive_contrast=False
    )

    observed_eigengaps = (_relative_eigengap(observed_configs[0]), _relative_eigengap(observed_configs[1]))
    reconstructed_eigengaps = (
        _relative_eigengap(captured_configs[0]),
        _relative_eigengap(captured_configs[1]),
    )
    sources: list[str] = []
    for label, gaps in (("observed", observed_eigengaps), ("reconstructed", reconstructed_eigengaps)):
        for group_index, gap in enumerate(gaps):
            if gap is None or gap < eigengap_threshold:
                sources.append(f"{label}:{groups[group_index]}")
    return PrincipalOrientationAttribution(
        component_id=component_id,
        from_stage=stages[0],
        to_stage=stages[-1],
        observed=observed_component,
        pls_captured=captured_component,
        residual=residual_component,
        observed_eigengaps=observed_eigengaps,
        reconstructed_eigengaps=reconstructed_eigengaps,
        eigengap_threshold=float(eigengap_threshold),
        degenerate=bool(sources),
        degeneracy_sources=tuple(sources),
        orientation_available=(
            observed_component.unit_directions[0] is not None,
            observed_component.unit_directions[1] is not None,
        ),
    )


def _validate_principal_settings(
    eigengap_threshold: float,
    component_id: str,
    stages: tuple[str, ...],
) -> None:
    if isinstance(eigengap_threshold, bool) or not np.isfinite(eigengap_threshold):
        raise AttributionError("eigengap_threshold must be a finite number in [0, 1].")
    if not 0.0 <= float(eigengap_threshold) <= 1.0:
        raise AttributionError("eigengap_threshold must be a finite number in [0, 1].")
    if not isinstance(component_id, str) or not component_id.strip():
        raise AttributionError("principal_component_id must be a non-empty string.")
    collisions = [
        f"{stages[index]}->{stages[index + 1]}"
        for index in range(len(stages) - 1)
        if f"{stages[index]}->{stages[index + 1]}" == component_id
    ]
    if collisions:
        raise AttributionError(
            f"principal_component_id='{component_id}' collides with stage-derived transition id(s) "
            f"{collisions}; choose a different reserved identifier."
        )


def _bootstrap(
    X: pd.DataFrame,
    metadata: pd.DataFrame,
    group_col: str,
    stage_col: str,
    groups: tuple[str, str],
    stages: tuple[str, ...],
    model: Any,
    feature_names: tuple[str, ...],
    scales: np.ndarray | None,
    replicates: int,
    seed: int | None,
    top_k: int,
    tolerance: float,
    principal_component_id: str,
    adjustment: _CovariateAdjustment | None = None,
) -> tuple[tuple[BootstrapSummary, ...], int]:
    """Stratified bootstrap stability for every transition and the principal block.

    Each replicate recomputes the group-by-stage means, rebuilds the transitions
    and the principal-orientation block, and contributes its directional
    contrasts. The principal axis is signed by *that replicate's own* net
    displacement — the shared functional's convention — so resampling noise
    cannot flip the axis and no post-hoc alignment to the point estimate is
    needed. A replicate whose net displacement falls at or below ``tolerance``
    contributes no contrast and counts as invalid, exactly like a replicate with
    a zero transition.

    With an ``adjustment`` the draws are the same, and each replicate's means
    are the adjusted means of its resampled rows: the covariate coefficients
    are refitted on them and evaluated at the fixed reference row. A resample
    whose covariates are collinear within cells contributes no contrast to any
    block; the second return value counts those replicates.
    """

    del scales  # Bootstrap stability is summarized in the fixed input units.
    block_ids = [f"{stages[index]}->{stages[index + 1]}" for index in range(len(stages) - 1)]
    block_ids.append(principal_component_id)
    if replicates == 0:
        return tuple(
            BootstrapSummary(
                transition_id=block_id,
                component=component,
                requested_replicates=0,
                valid_replicates=0,
                top_k=top_k,
                seed=seed,
                sign_stability=np.full(len(feature_names), np.nan),
                rank_stability=np.full(len(feature_names), np.nan),
                top_k_selection_frequency=np.full(len(feature_names), np.nan),
                nonzero_sign_replicates=np.zeros(len(feature_names), dtype=int),
            )
            for block_id in block_ids
            for component in ("observed", "pls_captured", "residual")
        ), 0
    group_values = metadata[group_col].astype(str).to_numpy()
    stage_values = metadata[stage_col].astype(str).to_numpy()
    strata = [
        np.flatnonzero((group_values == group) & (stage_values == stage))
        for group in groups
        for stage in stages
    ]
    rng = np.random.default_rng(seed)
    per_component: dict[tuple[str, str], list[np.ndarray]] = {}
    rank_deficient = 0
    if adjustment is not None:
        values = X.to_numpy(dtype=float)
        # Drawn rows are concatenated in stratum order, so their cell codes are fixed.
        drawn_cells = np.repeat(np.arange(len(strata)), [len(stratum) for stratum in strata])
    for _ in range(replicates):
        sampled = np.empty((len(groups), len(stages), X.shape[1]), dtype=float)
        drawn: list[np.ndarray] = []
        for cell_index, source_indices in enumerate(strata):
            selected = rng.choice(source_indices, size=len(source_indices), replace=True)
            if adjustment is None:
                sampled[cell_index // len(stages), cell_index % len(stages)] = (
                    X.iloc[selected].mean(axis=0).to_numpy()
                )
            else:
                drawn.append(selected)
        if adjustment is not None:
            rows = np.concatenate(drawn)
            try:
                sampled = _adjusted_means(
                    values[rows], adjustment.matrix[rows], drawn_cells, adjustment, groups, stages
                )
            except ValueError:
                rank_deficient += 1
                for block_id in block_ids:
                    for component_name in ("observed", "pls_captured", "residual"):
                        per_component.setdefault((block_id, component_name), []).append(
                            np.full(len(feature_names), np.nan)
                        )
                continue
        reconstructed = _reconstruct_means(sampled, model, feature_names)
        transitions = _build_transitions(sampled, reconstructed, groups, stages, tolerance)
        # Per-replicate degeneracy flags are not summarized — stability is the
        # summary — so the threshold argument is inert here.
        principal = _build_principal_orientation(
            sampled,
            reconstructed,
            groups,
            stages,
            tolerance,
            0.0,
            principal_component_id,
        )
        blocks: list[tuple[str, Any]] = [(record.transition_id, record) for record in transitions]
        blocks.append((principal.component_id, principal))
        for block_id, block in blocks:
            for component_name in ("observed", "pls_captured", "residual"):
                component = getattr(block, component_name)
                key = (block_id, component_name)
                per_component.setdefault(key, []).append(
                    component.directional_contrast
                    if component.directional_contrast is not None
                    else np.full(len(feature_names), np.nan)
                )
    summaries: list[BootstrapSummary] = []
    for transition_id in block_ids:
        for component_name in ("observed", "pls_captured", "residual"):
            values = np.asarray(per_component[(transition_id, component_name)], dtype=float)
            valid_mask = np.isfinite(values).all(axis=1)
            valid_values = values[valid_mask]
            valid_count = len(valid_values)
            sign_stability = np.full(len(feature_names), np.nan)
            rank_stability = np.full(len(feature_names), np.nan)
            top_frequency = np.full(len(feature_names), np.nan)
            nonzero_counts = np.zeros(len(feature_names), dtype=int)
            selected_counts = np.zeros(len(feature_names), dtype=int)
            for vector in valid_values:
                nonzero = np.abs(vector) > tolerance
                nonzero_counts += nonzero
                if nonzero.any():
                    selected = np.argsort(-np.where(nonzero, np.abs(vector), -np.inf), kind="stable")[:top_k]
                    selected = selected[nonzero[selected]]
                    selected_counts[selected] += 1
            if valid_count:
                top_frequency = selected_counts / valid_count
                rank_stability = top_frequency.copy()
                for feature_index in range(len(feature_names)):
                    nonzero_vectors = valid_values[:, feature_index]
                    nonzero_vectors = nonzero_vectors[np.abs(nonzero_vectors) > tolerance]
                    nonzero_counts[feature_index] = len(nonzero_vectors)
                    if len(nonzero_vectors):
                        modal_sign = 1.0 if np.sum(nonzero_vectors > 0) >= np.sum(nonzero_vectors < 0) else -1.0
                        sign_stability[feature_index] = float(
                            np.sum(np.sign(nonzero_vectors) == modal_sign) / len(nonzero_vectors)
                        )
            summaries.append(
                BootstrapSummary(
                    transition_id=transition_id,
                    component=component_name,
                    requested_replicates=replicates,
                    valid_replicates=valid_count,
                    top_k=top_k,
                    seed=seed,
                    sign_stability=sign_stability,
                    rank_stability=rank_stability,
                    top_k_selection_frequency=top_frequency,
                    nonzero_sign_replicates=nonzero_counts,
                )
            )
    return tuple(summaries), rank_deficient


def _validate_bootstrap_settings(replicates: int, top_k: int, tolerance: float) -> None:
    if isinstance(replicates, bool) or replicates < 0:
        raise AttributionError("bootstrap_replicates must be a non-negative integer.")
    if int(replicates) != replicates:
        raise AttributionError("bootstrap_replicates must be a non-negative integer.")
    if isinstance(top_k, bool) or top_k < 1:
        raise AttributionError("top_k must be a positive integer.")
    if int(top_k) != top_k:
        raise AttributionError("top_k must be a positive integer.")
    if not np.isfinite(tolerance) or tolerance < 0:
        raise AttributionError("zero_tolerance must be finite and non-negative.")


def _component_effect(component: ComponentAttribution, p: int) -> np.ndarray:
    if component.directional_contrast is None:
        return np.full(p, np.nan)
    return component.directional_contrast


def _feature_effect_table(
    transitions: Sequence[TransitionAttribution],
    principal: PrincipalOrientationAttribution,
    feature_names: tuple[str, ...],
    scales: np.ndarray | None,
    bootstrap: Mapping[tuple[str, str], BootstrapSummary],
    labels: Mapping[str, str] | None,
    label_namespace: str,
    requested_replicates: int,
) -> tuple[pd.DataFrame, list[FeatureEffect]]:
    """Per-feature signed effects for every transition and the principal block.

    The principal-orientation rows carry the reserved component identifier in
    ``transition_id`` and the first/last stage in ``from_stage``/``to_stage``,
    so the frame schema is unchanged and a consumer that groups by
    ``transition_id`` sees one extra group rather than a new column layout.
    """

    rows: list[dict[str, Any]] = []
    records: list[FeatureEffect] = []
    blocks: list[tuple[str, str, str, Any]] = [
        (transition.transition_id, transition.from_stage, transition.to_stage, transition)
        for transition in transitions
    ]
    blocks.append((principal.component_id, principal.from_stage, principal.to_stage, principal))
    for transition_id, from_stage, to_stage, block in blocks:
        for component_name in ("observed", "pls_captured", "residual"):
            component = getattr(block, component_name)
            effects = _component_effect(component, len(feature_names))
            summary = bootstrap[(transition_id, component_name)]
            for index, feature in enumerate(feature_names):
                standardized = float(effects[index])
                original = (
                    None
                    if scales is None or not np.isfinite(standardized)
                    else float(standardized * scales[index])
                )
                records.append(
                    FeatureEffect(transition_id, component_name, feature, standardized, original)
                )
                rows.append(
                    {
                        "transition_id": transition_id,
                        "from_stage": from_stage,
                        "to_stage": to_stage,
                        "component": component_name,
                        "feature": feature,
                        "effect_standardized": standardized,
                        "effect_original": original,
                        "effect_unit": "standardized_input",
                        "feature_group": labels[feature] if labels is not None else None,
                        "label_namespace": label_namespace if labels is not None else None,
                        "label_source": "caller-supplied" if labels is not None else None,
                        "sign_stability": _array_value(summary.sign_stability, index),
                        "rank_stability": _array_value(summary.rank_stability, index),
                        "top_k_selection_frequency": _array_value(summary.top_k_selection_frequency, index),
                        "valid_bootstrap_replicates": summary.valid_replicates,
                        "nonzero_sign_replicates": int(summary.nonzero_sign_replicates[index]),
                        "requested_bootstrap_replicates": requested_replicates,
                    }
                )
    return pd.DataFrame(rows), records


def _aggregate_effect_table(
    feature_effects: pd.DataFrame,
    feature_names: tuple[str, ...],
    labels: Mapping[str, str] | None,
    label_namespace: str,
    requested_replicates: int,
) -> tuple[pd.DataFrame, list[AggregateEffect]]:
    del feature_names
    if labels is None:
        return pd.DataFrame(), []
    rows: list[dict[str, Any]] = []
    records: list[AggregateEffect] = []
    for (transition_id, component, label), group in feature_effects.groupby(
        ["transition_id", "component", "feature_group"], sort=False, dropna=False
    ):
        standardized = group["effect_standardized"].to_numpy(dtype=float)
        original = group["effect_original"].to_numpy(dtype=float)
        signed = float(np.nansum(standardized)) if np.isfinite(standardized).any() else np.nan
        absolute = float(np.nansum(np.abs(standardized))) if np.isfinite(standardized).any() else np.nan
        original_signed = float(np.nansum(original)) if np.isfinite(original).any() else np.nan
        original_absolute = float(np.nansum(np.abs(original))) if np.isfinite(original).any() else np.nan
        def mean_available(column: str) -> float:
            values = group[column].to_numpy(dtype=float)
            return float(np.nanmean(values)) if np.isfinite(values).any() else np.nan

        row = {
            "transition_id": transition_id,
            "component": component,
            "label": label,
            "label_namespace": label_namespace,
            "label_source": "caller-supplied",
            "signed_effect_standardized": signed,
            "absolute_effect_standardized": absolute,
            "signed_effect_original": original_signed,
            "absolute_effect_original": original_absolute,
            "feature_count": int(len(group)),
            "sign_stability": mean_available("sign_stability"),
            "rank_stability": mean_available("rank_stability"),
            "top_k_selection_frequency": mean_available("top_k_selection_frequency"),
            "valid_bootstrap_replicates": int(group["valid_bootstrap_replicates"].max()),
            "requested_bootstrap_replicates": requested_replicates,
        }
        rows.append(row)
        records.append(
            AggregateEffect(transition_id, component, str(label), signed, absolute, len(group), label_namespace)
        )
    return pd.DataFrame(rows), records


def _transition_tables(
    transitions: Sequence[TransitionAttribution],
    groups: tuple[str, str],
    feature_names: tuple[str, ...],
) -> tuple[pd.DataFrame, pd.DataFrame]:
    summary_rows: list[dict[str, Any]] = []
    vector_rows: list[dict[str, Any]] = []
    for transition in transitions:
        for component_name in ("observed", "pls_captured", "residual"):
            component = getattr(transition, component_name)
            summary_rows.append(
                {
                    "transition_id": transition.transition_id,
                    "from_stage": transition.from_stage,
                    "to_stage": transition.to_stage,
                    "component": component_name,
                    "group_1": groups[0],
                    "group_2": groups[1],
                    "path_length_group_1": component.path_lengths[0],
                    "path_length_group_2": component.path_lengths[1],
                    "contrast_available": component.contrast_available,
                    "contrast_norm": (
                        float(np.linalg.norm(component.directional_contrast))
                        if component.directional_contrast is not None
                        else np.nan
                    ),
                }
            )
            for group_index, group in enumerate(groups):
                for feature_index, feature in enumerate(feature_names):
                    vector_rows.append(
                        {
                            "transition_id": transition.transition_id,
                            "from_stage": transition.from_stage,
                            "to_stage": transition.to_stage,
                            "component": component_name,
                            "group": group,
                            "feature": feature,
                            "raw_transition": float(component.group_transitions[group_index][feature_index]),
                        }
                    )
    return pd.DataFrame(summary_rows), pd.DataFrame(vector_rows)


def _principal_tables(
    principal: PrincipalOrientationAttribution,
    groups: tuple[str, str],
    feature_names: tuple[str, ...],
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Summary and per-feature-axis views of the principal-orientation block.

    The pair mirrors ``transition_summaries`` / ``transition_vectors``: one
    summary row per component carrying the degeneracy flag, the per-group
    relative eigengaps, and the contrast norm; and one long row per component,
    group, and feature carrying the signed principal axis together with the net
    displacement it was signed by.
    """

    summary_rows: list[dict[str, Any]] = []
    vector_rows: list[dict[str, Any]] = []
    for component_name in ("observed", "pls_captured", "residual"):
        component = getattr(principal, component_name)
        summary_rows.append(
            {
                "transition_id": principal.component_id,
                "from_stage": principal.from_stage,
                "to_stage": principal.to_stage,
                "component": component_name,
                "group_1": groups[0],
                "group_2": groups[1],
                "net_displacement_group_1": component.path_lengths[0],
                "net_displacement_group_2": component.path_lengths[1],
                "orientation_available_group_1": bool(component.unit_directions[0] is not None),
                "orientation_available_group_2": bool(component.unit_directions[1] is not None),
                "contrast_available": component.contrast_available,
                "contrast_norm": (
                    float(np.linalg.norm(component.directional_contrast))
                    if component.directional_contrast is not None
                    else np.nan
                ),
                "relative_eigengap_observed_group_1": _optional_float(principal.observed_eigengaps[0]),
                "relative_eigengap_observed_group_2": _optional_float(principal.observed_eigengaps[1]),
                "relative_eigengap_reconstructed_group_1": _optional_float(
                    principal.reconstructed_eigengaps[0]
                ),
                "relative_eigengap_reconstructed_group_2": _optional_float(
                    principal.reconstructed_eigengaps[1]
                ),
                "eigengap_threshold": principal.eigengap_threshold,
                "degenerate": principal.degenerate,
                "degeneracy_sources": ";".join(principal.degeneracy_sources),
            }
        )
        for group_index, group in enumerate(groups):
            axis = component.unit_directions[group_index]
            for feature_index, feature in enumerate(feature_names):
                vector_rows.append(
                    {
                        "transition_id": principal.component_id,
                        "from_stage": principal.from_stage,
                        "to_stage": principal.to_stage,
                        "component": component_name,
                        "group": group,
                        "feature": feature,
                        "principal_axis": (
                            float(axis[feature_index]) if axis is not None else np.nan
                        ),
                        "net_displacement": float(
                            component.group_transitions[group_index][feature_index]
                        ),
                    }
                )
    return pd.DataFrame(summary_rows), pd.DataFrame(vector_rows)


def _optional_float(value: float | None) -> float:
    return np.nan if value is None else float(value)


def _bootstrap_table(records: Sequence[BootstrapSummary], feature_names: tuple[str, ...]) -> pd.DataFrame:
    rows: list[dict[str, Any]] = []
    for record in records:
        for index, feature in enumerate(feature_names):
            rows.append(
                {
                    "transition_id": record.transition_id,
                    "component": record.component,
                    "feature": feature,
                    "requested_replicates": record.requested_replicates,
                    "valid_replicates": record.valid_replicates,
                    "top_k": record.top_k,
                    "seed": record.seed,
                    "sign_stability": _array_value(record.sign_stability, index),
                    "rank_stability": _array_value(record.rank_stability, index),
                    "top_k_selection_frequency": _array_value(record.top_k_selection_frequency, index),
                    "nonzero_sign_replicates": int(record.nonzero_sign_replicates[index]),
                }
            )
    return pd.DataFrame(rows)


def _means_table(
    values: np.ndarray,
    groups: tuple[str, str],
    stages: tuple[str, ...],
    feature_names: tuple[str, ...],
) -> pd.DataFrame:
    rows: list[dict[str, Any]] = []
    for group_index, group in enumerate(groups):
        for stage_index, stage in enumerate(stages):
            row: dict[str, Any] = {"group": group, "stage": stage}
            row.update(
                {
                    feature: float(values[group_index, stage_index, index])
                    for index, feature in enumerate(feature_names)
                }
            )
            rows.append(row)
    return pd.DataFrame(rows, columns=["group", "stage", *feature_names])


def _array_value(values: np.ndarray, index: int) -> float:
    value = values[index]
    return float(value) if np.isfinite(value) else np.nan


__all__ = [
    "ADJUSTMENT_BOUNDARY",
    "AggregateEffect",
    "AttributionConfig",
    "AttributionError",
    "AttributionResult",
    "BootstrapSummary",
    "ComponentAttribution",
    "FeatureEffect",
    "InterpretationMetadata",
    "OrientationAttributionResult",
    "PrincipalOrientationAttribution",
    "TransitionAttribution",
    "analyze_orientation_attribution",
    "attribution_frames",
    "run_orientation_attribution",
    "write_attribution_outputs",
]
