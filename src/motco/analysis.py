"""Real-data trajectory analysis: the validated MOTCO measurement on your own data.

:func:`analyze_trajectories` runs the production pipeline the simulation harness
characterizes — pooled per-feature standardization, a stage-supervised PLS-DA
latent space whose rank is chosen by double cross-validation, ``delta`` /
``angle`` / ``shape`` between two groups in that space, and RRPP p-values — on
caller-supplied omic blocks and sample metadata. It calls the same cores as
:func:`motco.simulations.evaluate_semisynthetic_trajectory` (block scaling in
:mod:`motco.simulations.preprocessing`, :func:`~motco.simulations.evaluation.fit_pls_latent_space`,
the trajectory design and the record helpers), so the measurement exists once.

Two inputs are load-bearing and therefore explicit:

- **Stage order.** Net displacement runs from the first to the last declared
  stage, and ``shape`` and the spectrum follow that order. It is never inferred
  by sorting labels.
- **Units.** Each block declares its unit transform (``identity`` or
  ``beta_to_mvalue``). MOTCO does not normalize counts, filter features or
  impute values; blocks are analysed in the units supplied.

Declared covariates enter the trajectory model. By default the latent space is
built without them; ``TrajectoryAnalysisParams.covariate_latent_adjustment``
(``"stage"`` or ``"group_stage"``) instead fits the basis on features
residualized on the covariates conditional on that design, and measures the
raw features in it.

:func:`write_analysis_outputs` writes the result as plain JSON/CSV (record,
latent scores, PLS loadings, latent stage means and, when requested, the
attribution tables) — the fitted estimator is never pickled.
"""

from __future__ import annotations

import json
import platform
from dataclasses import asdict, dataclass, field
from pathlib import Path
from time import perf_counter
from typing import Any, Literal, Mapping, Sequence

import numpy as np
import pandas as pd

from motco.simulations.attribution_diagnostics import (
    AttributionDiagnosticError,
    TruthDrivers,
    run_attribution_diagnostics,
    unavailable_record,
)
from motco.simulations.evaluation import (
    EVALUATION_ORIENTATION_ESTIMATORS,
    AttributionDiagnosticSettings,
    LatentAdjustmentSummary,
    SimulationEvaluationError,
    _config_spectrum_block,
    _empirical_p_value,
    _extract_null_distributions,
    _extract_pair_statistics,
    _joint_original_scales,
    _summarize_null_distributions,
    _summarize_permutation_eigengaps,
    build_simulation_trajectory_design,
    fit_pls_latent_space,
)
from motco.simulations.preprocessing import (
    UNIT_TRANSFORMS,
    concatenate_blocks,
    fit_block_scalers,
    transform_blocks,
)
from motco.stats.attribution import OrientationAttributionResult, write_attribution_outputs
from motco.stats.design import (
    CONDITIONING_KINDS,
    CovariateDesign,
    build_conditioning_matrix,
    build_covariate_matrix,
)
from motco.stats.permutation import RRPP
from motco.stats.pls import NuisanceDesign
from motco.stats.trajectory import estimate_betas, estimate_difference

UnitTransform = Literal["identity", "beta_to_mvalue"]
MissingCovariatePolicy = Literal["error", "complete_case"]

#: Allowed missing-covariate policies.
MISSING_COVARIATE_POLICIES: tuple[str, ...] = ("error", "complete_case")

#: Version of the JSON-safe analysis record contract.
ANALYSIS_SCHEMA_VERSION = 1

#: Units each declared transform yields, as recorded per block.
_TRANSFORM_UNITS: dict[str, str] = {"identity": "as_supplied", "beta_to_mvalue": "mvalue"}

_NO_TRUTH_REASON = "real data carry no generator truth"

# Internal design column names on the recoded working frame.
_GROUP = "group"
_STAGE = "stage"


class TrajectoryAnalysisError(ValueError):
    """Raised when real-data analysis inputs or parameters are invalid."""


@dataclass(frozen=True)
class OmicsBlock:
    """One named omic block: a samples × features table and its unit transform."""

    name: str
    values: pd.DataFrame
    transform: UnitTransform = "identity"


@dataclass(frozen=True)
class TrajectoryDataset:
    """Aligned omic blocks plus sample metadata, validated on construction.

    Every block's rows must equal ``metadata[sample_col]`` in order. The group
    column must hold exactly two labels; ``stage_order`` must list every
    observed stage label exactly once (labels are compared as strings), and
    every group × stage cell must be observed.

    ``covariate_cols`` names metadata columns that enter the trajectory model
    additively, in declared order. A covariate is categorical when listed in
    ``categorical_covariates`` or when its dtype is not numeric (booleans are
    categorical); otherwise numeric. ``missing_covariates`` decides what a
    missing covariate value does: ``"error"`` rejects it, ``"complete_case"``
    drops the sample from the metadata and every block before anything else is
    validated or fitted (``dropped_sample_ids`` lists them). Covariates are
    never imputed.
    """

    blocks: Sequence[OmicsBlock]
    metadata: pd.DataFrame
    sample_col: str
    group_col: str
    stage_col: str
    stage_order: Sequence[str]
    group_order: Sequence[str] | None = None
    covariate_cols: Sequence[str] = ()
    categorical_covariates: Sequence[str] = ()
    missing_covariates: MissingCovariatePolicy = "error"
    _dropped_sample_ids: tuple[str, ...] = field(default=(), init=False, repr=False, compare=False)
    _covariate_design: CovariateDesign | None = field(default=None, init=False, repr=False, compare=False)

    def __post_init__(self) -> None:
        object.__setattr__(self, "blocks", tuple(self.blocks))
        object.__setattr__(self, "stage_order", tuple(str(stage) for stage in self.stage_order))
        if self.group_order is not None:
            object.__setattr__(self, "group_order", tuple(str(group) for group in self.group_order))
        object.__setattr__(self, "covariate_cols", tuple(self.covariate_cols))
        object.__setattr__(self, "categorical_covariates", tuple(self.categorical_covariates))
        self._validate_metadata()
        self._validate_covariate_declarations()
        cells_before = self._observed_cells()
        self._apply_missing_covariate_policy()
        self._validate_blocks(self.sample_ids)
        self._build_covariate_design()
        self._validate_design(cells_before)

    @property
    def sample_ids(self) -> list[str]:
        return self.metadata[self.sample_col].astype(str).tolist()

    @property
    def groups(self) -> tuple[str, ...]:
        """Group labels in analysis order (declared, else sorted as strings)."""

        if self.group_order is not None:
            return tuple(self.group_order)
        return tuple(sorted(set(self.metadata[self.group_col].astype(str))))

    @property
    def stages(self) -> tuple[str, ...]:
        return tuple(self.stage_order)

    @property
    def dropped_sample_ids(self) -> tuple[str, ...]:
        """Samples removed by complete-case deletion (empty under ``"error"``)."""

        return self._dropped_sample_ids

    @property
    def covariate_design(self) -> CovariateDesign | None:
        """The covariate design columns on the retained samples, or ``None`` without covariates."""

        return self._covariate_design

    def _observed_cells(self) -> set[tuple[str, str]]:
        return set(zip(self.metadata[self.group_col].astype(str), self.metadata[self.stage_col].astype(str)))

    def _validate_covariate_declarations(self) -> None:
        if self.missing_covariates not in MISSING_COVARIATE_POLICIES:
            raise TrajectoryAnalysisError(
                f"Unknown missing_covariates policy {self.missing_covariates!r}; allowed: "
                f"{list(MISSING_COVARIATE_POLICIES)}."
            )
        covariates = list(self.covariate_cols)
        repeated = sorted({name for name in covariates if covariates.count(name) > 1})
        if repeated:
            raise TrajectoryAnalysisError(f"Covariates must be declared once; repeated: {repeated}.")
        reserved = [name for name in covariates if name in {self.sample_col, self.group_col, self.stage_col}]
        if reserved:
            raise TrajectoryAnalysisError(
                f"Covariates must differ from the sample, group and stage columns; got {reserved}."
            )
        absent = [name for name in covariates if name not in self.metadata.columns]
        if absent:
            raise TrajectoryAnalysisError(f"metadata is missing covariate column(s): {absent}.")
        undeclared = [name for name in self.categorical_covariates if name not in covariates]
        if undeclared:
            raise TrajectoryAnalysisError(
                f"categorical_covariates must also be declared in covariate_cols; not declared: {undeclared}."
            )

    def _apply_missing_covariate_policy(self) -> None:
        if not self.covariate_cols:
            return
        missing = self.metadata[list(self.covariate_cols)].isna()
        if not missing.to_numpy().any():
            return
        sample_ids = self.metadata[self.sample_col].astype(str)
        if self.missing_covariates == "error":
            column = next(name for name in self.covariate_cols if missing[name].any())
            affected = sample_ids[missing[column].to_numpy()].tolist()
            raise TrajectoryAnalysisError(
                f"Covariate {column!r} has missing values for {len(affected)} sample(s): {affected[:10]}. "
                "Pass missing_covariates='complete_case' to drop them; covariates are never imputed."
            )
        drop = missing.any(axis=1).to_numpy()
        keep = ~drop
        object.__setattr__(self, "_dropped_sample_ids", tuple(sample_ids[drop].tolist()))
        object.__setattr__(self, "metadata", self.metadata.loc[keep].reset_index(drop=True))
        retained = []
        for block in self.blocks:
            values = block.values
            if isinstance(values, pd.DataFrame) and values.shape[0] == keep.size:
                values = values.loc[keep]
            retained.append(OmicsBlock(block.name, values, block.transform))
        object.__setattr__(self, "blocks", tuple(retained))

    def _build_covariate_design(self) -> None:
        if not self.covariate_cols:
            return
        try:
            design = build_covariate_matrix(self.metadata, self.covariate_cols, self.categorical_covariates)
        except ValueError as exc:
            raise TrajectoryAnalysisError(str(exc)) from exc
        object.__setattr__(self, "_covariate_design", design)

    def _validate_metadata(self) -> None:
        if not isinstance(self.metadata, pd.DataFrame):
            raise TrajectoryAnalysisError("metadata must be a pandas DataFrame.")
        columns = {"sample_col": self.sample_col, "group_col": self.group_col, "stage_col": self.stage_col}
        missing = [f"{role}={name!r}" for role, name in columns.items() if name not in self.metadata.columns]
        if missing:
            raise TrajectoryAnalysisError(f"metadata is missing column(s): {', '.join(missing)}.")
        if len({self.sample_col, self.group_col, self.stage_col}) != 3:
            raise TrajectoryAnalysisError("sample_col, group_col and stage_col must name three different columns.")
        if self.metadata.shape[0] == 0:
            raise TrajectoryAnalysisError("metadata must contain at least one sample.")
        for name in columns.values():
            if self.metadata[name].isna().any():
                raise TrajectoryAnalysisError(f"metadata column {name!r} contains missing values.")
        sample_ids = self.sample_ids
        duplicated = sorted({sid for sid in sample_ids if sample_ids.count(sid) > 1})
        if duplicated:
            raise TrajectoryAnalysisError(f"metadata sample ids are not unique: {duplicated[:10]}.")

    def _validate_blocks(self, sample_ids: list[str]) -> None:
        if not self.blocks:
            raise TrajectoryAnalysisError("At least one omic block is required.")
        names = [block.name for block in self.blocks]
        repeated = sorted({name for name in names if names.count(name) > 1})
        if repeated:
            raise TrajectoryAnalysisError(f"Block names must be unique; repeated: {repeated}.")
        for block in self.blocks:
            if not isinstance(block.name, str) or not block.name:
                raise TrajectoryAnalysisError(f"Block names must be non-empty strings, got {block.name!r}.")
            if block.transform not in UNIT_TRANSFORMS:
                raise TrajectoryAnalysisError(
                    f"Block {block.name!r} declares unknown transform {block.transform!r}; "
                    f"allowed: {list(UNIT_TRANSFORMS)}."
                )
            values = block.values
            if not isinstance(values, pd.DataFrame):
                raise TrajectoryAnalysisError(f"Block {block.name!r} must be a pandas DataFrame.")
            if values.shape[0] != len(sample_ids):
                raise TrajectoryAnalysisError(
                    f"Block {block.name!r} has {values.shape[0]} rows but metadata has {len(sample_ids)} samples."
                )
            if values.index.astype(str).tolist() != sample_ids:
                raise TrajectoryAnalysisError(
                    f"Block {block.name!r} rows are not aligned to metadata column {self.sample_col!r} "
                    "(same sample ids in the same order are required)."
                )
            if values.shape[1] == 0:
                raise TrajectoryAnalysisError(f"Block {block.name!r} must contain at least one feature.")
            features = values.columns.astype(str).tolist()
            if len(set(features)) != len(features):
                raise TrajectoryAnalysisError(f"Block {block.name!r} has duplicate feature names.")
            try:
                matrix = values.to_numpy(dtype=float)
            except (TypeError, ValueError) as exc:
                raise TrajectoryAnalysisError(f"Block {block.name!r} holds non-numeric values.") from exc
            if not np.all(np.isfinite(matrix)):
                raise TrajectoryAnalysisError(
                    f"Block {block.name!r} contains missing or infinite values; complete matrices are "
                    "required and no imputation is performed."
                )
            if block.transform == "beta_to_mvalue" and not np.all((matrix > 0.0) & (matrix < 1.0)):
                raise TrajectoryAnalysisError(
                    f"Block {block.name!r} declares beta_to_mvalue but holds values outside the open "
                    "interval (0, 1)."
                )

    def _validate_design(self, cells_before: set[tuple[str, str]]) -> None:
        emptied = sorted(cells_before - self._observed_cells())
        if emptied:
            raise TrajectoryAnalysisError(
                f"Complete-case deletion of {len(self._dropped_sample_ids)} sample(s) with missing covariates "
                f"emptied group/stage combination(s): {emptied}."
            )
        observed_stages = sorted(set(self.metadata[self.stage_col].astype(str)))
        declared = list(self.stage_order)
        if len(set(declared)) != len(declared) or set(declared) != set(observed_stages):
            raise TrajectoryAnalysisError(
                f"stage_order must list every observed stage label exactly once; declared {declared}, "
                f"observed {observed_stages}."
            )
        if len(declared) < 2:
            raise TrajectoryAnalysisError(f"At least two stages are required; observed {observed_stages}.")
        observed_groups = sorted(set(self.metadata[self.group_col].astype(str)))
        if len(observed_groups) != 2:
            raise TrajectoryAnalysisError(
                f"Exactly two groups are required in {self.group_col!r}; found {len(observed_groups)}: "
                f"{observed_groups}."
            )
        if self.group_order is not None:
            group_order = list(self.group_order)
            if len(set(group_order)) != len(group_order) or set(group_order) != set(observed_groups):
                raise TrajectoryAnalysisError(
                    f"group_order must list both observed groups exactly once; declared {group_order}, "
                    f"observed {observed_groups}."
                )
        present = self._observed_cells()
        missing = [(group, stage) for group in self.groups for stage in declared if (group, stage) not in present]
        if missing:
            raise TrajectoryAnalysisError(f"metadata has no sample for group/stage combination(s): {missing}.")


@dataclass(frozen=True)
class TrajectoryAnalysisParams:
    """Measurement settings; the defaults are the simulation harness's.

    ``forced_components`` fixes the latent rank (dimensionality sensitivity) and
    skips cross-validation; ``None`` selects it by double CV. ``seed`` seeds
    RRPP exactly as the harness does.

    ``covariate_latent_adjustment`` (``None``, ``"stage"`` or
    ``"group_stage"``) opts into covariate-aware latent construction: the
    basis is fitted on features residualized on the declared covariates,
    conditional on the stage or the group × stage indicators. It requires
    declared covariates and is omitted from :meth:`to_dict` while ``None``.
    """

    permutations: int = 0
    seed: int | None = None
    n_jobs: int = 1
    orientation_estimator: str = "net_displacement"
    forced_components: int | None = None
    cv1_splits: int = 3
    cv2_splits: int = 4
    n_repeats: int = 5
    max_components: int = 20
    random_state: int = 1203
    attribution: AttributionDiagnosticSettings = field(default_factory=AttributionDiagnosticSettings)
    progress: bool = False
    covariate_latent_adjustment: str | None = None

    def __post_init__(self) -> None:
        if self.covariate_latent_adjustment is not None and self.covariate_latent_adjustment not in CONDITIONING_KINDS:
            raise TrajectoryAnalysisError(
                f"Unknown covariate_latent_adjustment {self.covariate_latent_adjustment!r}; allowed: "
                f"{list(CONDITIONING_KINDS)} or None."
            )
        if self.orientation_estimator not in EVALUATION_ORIENTATION_ESTIMATORS:
            raise TrajectoryAnalysisError(
                f"Unknown orientation_estimator {self.orientation_estimator!r}; allowed: "
                f"{list(EVALUATION_ORIENTATION_ESTIMATORS)}."
            )
        if self.permutations < 0:
            raise TrajectoryAnalysisError("permutations must be greater than or equal to 0.")
        if self.n_jobs == 0:
            raise TrajectoryAnalysisError("n_jobs must be -1 or a non-zero integer.")
        settings = self.attribution
        if settings.enabled:
            if settings.bootstrap_replicates < 0:
                raise TrajectoryAnalysisError("attribution.bootstrap_replicates must be greater than or equal to 0.")
            if settings.top_k < 1:
                raise TrajectoryAnalysisError("attribution.top_k must be at least 1.")

    def to_dict(self) -> dict[str, Any]:
        data = asdict(self)
        if self.covariate_latent_adjustment is None:
            # Absent means the 0.8.0 analysis: its resolved parameters stay byte-identical.
            data.pop("covariate_latent_adjustment")
        return data


@dataclass(frozen=True)
class TrajectoryAnalysisResult:
    """The analysis record plus the fitted representation.

    ``record`` is JSON-safe. ``latent_scores`` are the per-sample PLS scores,
    ``pls_loadings`` the joint-feature rotations (``block``/``feature`` columns
    first), ``stage_means`` the fitted latent LS-mean vector per group × stage
    in declared stage order, and ``attribution`` the full attribution result
    when it was requested.
    """

    record: dict[str, Any]
    latent_scores: pd.DataFrame
    pls_loadings: pd.DataFrame
    stage_means: pd.DataFrame
    attribution: OrientationAttributionResult | None = None

    def to_json(self, **kwargs: Any) -> str:
        return json.dumps(self.record, **kwargs)


def analyze_trajectories(
    dataset: TrajectoryDataset,
    params: TrajectoryAnalysisParams | None = None,
) -> TrajectoryAnalysisResult:
    """Measure two groups' trajectories across ordered stages on real data.

    Standardizes every block with pooled per-feature scalers (after its unit
    transform), concatenates the blocks in declared order as
    ``f"{block}__{feature}"``, builds the stage-supervised PLS-DA latent space
    (double-CV rank, or ``params.forced_components``), estimates ``delta``,
    ``angle`` (net displacement by default) and ``shape`` between the groups
    in that space and, when ``params.permutations > 0``, their RRPP p-values.
    Records the selected rank, null summary, stage-mean spectrum, cell sizes,
    block units and software versions beside the statistics. Optionally runs
    frozen-model orientation attribution on the fitted estimator.

    Declared covariates (``dataset.covariate_cols``) are appended additively to
    the full (group × stage) and reduced (group + stage) trajectory models, so
    RRPP permutes residuals of the reduced model with covariates. Every LS-mean
    is evaluated at the pooled covariate mean; ``delta``, ``angle`` and
    ``shape`` do not depend on that reference. The latent space is fitted
    without the covariates unless ``params.covariate_latent_adjustment`` is
    set: then the pooled-standardized features are residualized on the
    covariates conditional on the stage (``"stage"``) or group × stage
    (``"group_stage"``) indicators before the basis is fitted (per training set
    inside the cross-validation, at the pooled scale), the raw features are
    projected onto that basis, and the trajectory model, LS-means and RRPP are
    unchanged. Attribution is not available with covariates.
    """

    params = params or TrajectoryAnalysisParams()
    covariates = dataset.covariate_design
    if params.covariate_latent_adjustment is not None and covariates is None:
        raise TrajectoryAnalysisError(
            f"Latent adjustment ({params.covariate_latent_adjustment!r}) requires declared covariates; "
            "declare covariate_cols on the dataset or drop the setting."
        )
    if covariates is not None and params.attribution.enabled:
        raise TrajectoryAnalysisError(
            "Attribution is not available with covariates: it works from unadjusted cell means. "
            "Disable attribution or drop the covariates."
        )
    start = perf_counter()
    blocks = {block.name: block.values for block in dataset.blocks}
    order = [block.name for block in dataset.blocks]
    transforms = {block.name: block.transform for block in dataset.blocks}

    scalers = fit_block_scalers(blocks, transforms)
    standardized = transform_blocks(blocks, scalers, transforms)
    X = concatenate_blocks(standardized, order).reset_index(drop=True)

    groups = list(dataset.groups)
    stages = list(dataset.stages)
    group_codes = {label: str(index) for index, label in enumerate(groups)}
    stage_codes = {label: str(index) for index, label in enumerate(stages)}
    working = pd.DataFrame(
        {
            _GROUP: dataset.metadata[dataset.group_col].astype(str).map(group_codes).to_numpy(),
            _STAGE: dataset.metadata[dataset.stage_col].astype(str).map(stage_codes).to_numpy(),
        }
    )

    integration_params: dict[str, Any] = {
        "cv1_splits": params.cv1_splits,
        "cv2_splits": params.cv2_splits,
        "n_repeats": params.n_repeats,
        "max_components": params.max_components,
        "random_state": params.random_state,
        "n_jobs": params.n_jobs,
    }
    if params.forced_components is not None:
        integration_params["forced_components"] = params.forced_components
    nuisance = _latent_nuisance(working, covariates, params.covariate_latent_adjustment)
    sample_ids = dataset.sample_ids
    try:
        latent = fit_pls_latent_space(
            X,
            working[_STAGE],
            integration_params,
            stage_col=dataset.stage_col,
            sample_ids=sample_ids,
            layer_feature_counts={name: int(frame.shape[1]) for name, frame in blocks.items()},
            layers_metadata={},
            original_scales=_joint_original_scales(scalers, order, X.columns),
            methylation_units=_methylation_units(transforms),
            nuisance=nuisance,
        )
        design = build_simulation_trajectory_design(working, group_col=_GROUP, stage_col=_STAGE)
    except SimulationEvaluationError as exc:
        raise TrajectoryAnalysisError(str(exc)) from exc
    group_levels = [groups[int(code)] for code in design.group_levels]
    stage_levels = [stages[int(code)] for code in design.stage_levels]
    if stage_levels != stages or group_levels != groups:
        raise TrajectoryAnalysisError(
            f"Internal level mapping mismatch: design levels {group_levels}/{stage_levels}, "
            f"declared {groups}/{stages}."
        )
    model_full, model_reduced, ls_means = _covariate_adjusted_design(
        design.model_full, design.model_reduced, design.ls_means, covariates
    )

    observed_deltas, observed_angles, observed_shapes, observed_spectra = estimate_difference(
        Y=latent.matrix,
        model_matrix=model_full,
        LS_means=ls_means,
        contrast=design.contrast,
        return_spectra=True,
        orientation_estimator=params.orientation_estimator,
    )
    config_spectrum = _config_spectrum_block(observed_spectra, group_levels)
    shape_available = len(stage_levels) >= 3
    pair_statistics = _extract_pair_statistics(
        observed_deltas, observed_angles, observed_shapes, shape_available=shape_available
    )

    p_values: dict[str, float] = {}
    null_summary: dict[str, dict[str, float]] = {}
    if params.permutations > 0:
        dist_delta, dist_angle, dist_shape, dist_eigengap = RRPP(
            latent.matrix,
            model_full,
            model_reduced,
            ls_means,
            design.contrast,
            permutations=params.permutations,
            n_jobs=params.n_jobs,
            progress=params.progress,
            seed=params.seed,
            return_eigengaps=True,
            orientation_estimator=params.orientation_estimator,
        )
        null_values = _extract_null_distributions(
            dist_delta, dist_angle, dist_shape, shape_available=shape_available
        )
        p_values = {
            statistic: _empirical_p_value(null_values[statistic], observed)
            for statistic, observed in pair_statistics.items()
            if statistic in null_values and np.isfinite(observed)
        }
        null_summary = _summarize_null_distributions(null_values)
        config_spectrum["permutation_pooled_eigengap"] = _summarize_permutation_eigengaps(dist_eigengap)

    latent_scores = latent.matrix.copy()
    latent_scores.index.name = dataset.sample_col
    artifacts = latent.artifacts
    assert artifacts is not None  # the PLS core always returns its fitted artifacts
    pls_loadings = _loadings_frame(artifacts.model, X.columns, list(latent.matrix.columns))
    betas = np.asarray(estimate_betas(model_full, latent.matrix.to_numpy(dtype=float)), dtype=float)
    stage_means = _stage_means_frame(
        ls_means,
        betas,
        group_levels,
        stage_levels,
        list(latent.matrix.columns),
    )

    attribution_record, attribution_result = _run_attribution(dataset, params, latent, X, transforms)
    cell_sizes = {
        group: {
            stage: int(
                (
                    (dataset.metadata[dataset.group_col].astype(str) == group)
                    & (dataset.metadata[dataset.stage_col].astype(str) == stage)
                ).sum()
            )
            for stage in stage_levels
        }
        for group in group_levels
    }
    record: dict[str, Any] = {
        "schema_version": ANALYSIS_SCHEMA_VERSION,
        "statistics": pair_statistics,
        "p_values": p_values,
        "null_summary": null_summary,
        "config_spectrum": config_spectrum,
        "latent": latent.metadata,
        "design": {
            "group_col": dataset.group_col,
            "stage_col": dataset.stage_col,
            "sample_col": dataset.sample_col,
            "group_levels": group_levels,
            "stage_levels": stage_levels,
            "cell_sizes": cell_sizes,
            "contrast": design.contrast,
            "n_samples": len(sample_ids),
        },
        **(
            {"covariates": _covariate_record(dataset, covariates, betas, list(latent.matrix.columns))}
            if covariates is not None
            else {}
        ),
        **(
            {"latent_adjustment": _latent_adjustment_record(dataset, artifacts.latent_adjustment)}
            if artifacts.latent_adjustment is not None
            else {}
        ),
        "blocks": [
            {
                "name": block.name,
                "n_features": int(block.values.shape[1]),
                "transform": block.transform,
                "units": _TRANSFORM_UNITS[block.transform],
            }
            for block in dataset.blocks
        ],
        "measurement": {
            "orientation_estimator": params.orientation_estimator,
            "selected_lv": latent.metadata["selected_lv"],
            "component_selection": latent.metadata["component_selection"],
            "permutations": params.permutations,
            "seed": params.seed,
            "n_jobs": params.n_jobs,
            "shape_available": shape_available,
            "runtime_seconds": perf_counter() - start,
        },
        "params": params.to_dict(),
        "software": _software_versions(),
        "attribution": attribution_record,
    }
    return TrajectoryAnalysisResult(
        record=_json_safe(record),
        latent_scores=latent_scores,
        pls_loadings=pls_loadings,
        stage_means=stage_means,
        attribution=attribution_result,
    )


def write_analysis_outputs(
    result: TrajectoryAnalysisResult,
    out_dir: str | Path,
    *,
    overwrite: bool = False,
) -> dict[str, Path]:
    """Write the analysis artifact to ``out_dir`` and return the written paths.

    Writes ``result.json``, ``latent_scores.csv``, ``pls_loadings.csv``,
    ``stage_means.csv`` and, when attribution ran, its tables under
    ``attribution/``. Refuses to touch a directory that already holds
    ``result.json`` unless ``overwrite`` is set.
    """

    directory = Path(out_dir)
    result_path = directory / "result.json"
    if result_path.exists() and not overwrite:
        raise FileExistsError(f"{result_path} already exists; pass overwrite=True to replace the analysis.")
    directory.mkdir(parents=True, exist_ok=True)
    paths: dict[str, Path] = {}
    result_path.write_text(json.dumps(result.record, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    paths["result"] = result_path
    paths["latent_scores"] = directory / "latent_scores.csv"
    result.latent_scores.to_csv(paths["latent_scores"])
    paths["pls_loadings"] = directory / "pls_loadings.csv"
    result.pls_loadings.to_csv(paths["pls_loadings"], index=False)
    paths["stage_means"] = directory / "stage_means.csv"
    result.stage_means.to_csv(paths["stage_means"], index=False)
    if result.attribution is not None:
        for name, path in write_attribution_outputs(result.attribution, directory / "attribution").items():
            paths[f"attribution/{name}"] = path
    return paths


def _run_attribution(
    dataset: TrajectoryDataset,
    params: TrajectoryAnalysisParams,
    latent: Any,
    X: pd.DataFrame,
    transforms: Mapping[str, str],
) -> tuple[dict[str, Any], OrientationAttributionResult | None]:
    settings = params.attribution
    if not settings.enabled:
        return unavailable_record("attribution was not requested"), None
    artifacts = latent.artifacts
    metadata = pd.DataFrame(
        {
            dataset.group_col: dataset.metadata[dataset.group_col].astype(str).to_numpy(),
            dataset.stage_col: dataset.metadata[dataset.stage_col].astype(str).to_numpy(),
        }
    )
    truth = TruthDrivers(
        feature_names=tuple(str(column) for column in X.columns),
        global_mask=np.zeros(X.shape[1], dtype=bool),
        per_transition={},
        available=False,
        definition=f"unavailable: {_NO_TRUTH_REASON}",
    )
    try:
        return run_attribution_diagnostics(
            model=artifacts.model,
            features=artifacts.features,
            original_scales=artifacts.original_scales,
            settings=settings,
            group_col=dataset.group_col,
            stage_col=dataset.stage_col,
            selected_components=latent.metadata.get("selected_lv"),
            feature_order_signature=latent.metadata.get("feature_order_signature"),
            methylation_units=artifacts.methylation_units,
            metadata=metadata,
            groups=list(dataset.groups),
            stages=list(dataset.stages),
            truth=truth,
            original_basis={name: _TRANSFORM_UNITS[kind] for name, kind in transforms.items()},
        )
    except AttributionDiagnosticError as exc:
        raise TrajectoryAnalysisError(str(exc)) from exc


def _covariate_adjusted_design(
    model_full: np.ndarray,
    model_reduced: np.ndarray,
    ls_means: np.ndarray,
    covariates: CovariateDesign | None,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Append the covariate columns to both models and the reference row to every LS-mean row.

    The covariate columns go last, so the group, stage and interaction
    positions (and therefore the contrast) are unchanged. Without covariates
    the matrices are returned as they are.
    """

    if covariates is None:
        return model_full, model_reduced, ls_means
    full = np.hstack([np.asarray(model_full, dtype=float), covariates.matrix])
    rank = int(np.linalg.matrix_rank(full))
    if rank < full.shape[1]:
        raise TrajectoryAnalysisError(
            f"The trajectory model with covariates {list(covariates.columns)} is not of full rank "
            f"(rank {rank} of {full.shape[1]} columns); a covariate is collinear with the group × stage "
            "design or with another covariate."
        )
    reduced = np.hstack([np.asarray(model_reduced, dtype=float), covariates.matrix])
    reference = np.tile(covariates.reference_row, (np.asarray(ls_means).shape[0], 1))
    return full, reduced, np.hstack([np.asarray(ls_means, dtype=float), reference])


def _latent_nuisance(
    working: pd.DataFrame,
    covariates: CovariateDesign | None,
    kind: str | None,
) -> NuisanceDesign | None:
    """The residualization design for latent adjustment, or ``None`` when it is off.

    The conditioning indicators are built on the recoded working frame, so their
    names use the ordinal codes (``stage[1]``) that map back through the
    declared orders. The full-sample design must have full column rank.
    """

    if kind is None or covariates is None:
        return None
    conditioning = build_conditioning_matrix(working, _GROUP, _STAGE, kind)
    nuisance = NuisanceDesign(
        conditioning.matrix, covariates.matrix, conditioning.column_names, covariates.column_names, kind
    )
    full = nuisance.matrix()
    rank = int(np.linalg.matrix_rank(full))
    if rank < full.shape[1]:
        raise TrajectoryAnalysisError(
            f"The latent-adjustment residualization design with covariates {list(covariates.columns)} under "
            f"{kind!r} conditioning is not of full rank (rank {rank} of {full.shape[1]} columns: intercept, "
            f"{', '.join(nuisance.column_names)}); a covariate is collinear with the conditioning design or "
            "with another covariate."
        )
    return nuisance


def _latent_adjustment_record(dataset: TrajectoryDataset, summary: LatentAdjustmentSummary) -> dict[str, Any]:
    return {
        "conditioning": summary.conditioning,
        "covariates": list(dataset.covariate_cols),
        "design_columns": list(summary.design_columns),
        "n_samples": summary.n_samples,
        "residual_df": summary.residual_df,
        "cv_min_training_df": summary.cv_min_training_df,
        "feature_scale": summary.feature_scale,
        "feature_r2": dict(summary.feature_r2),
    }


def _covariate_record(
    dataset: TrajectoryDataset,
    covariates: CovariateDesign,
    betas: np.ndarray,
    components: list[str],
) -> dict[str, Any]:
    q = len(covariates.column_names)
    coefficients = betas[betas.shape[0] - q :, :]
    return {
        "columns": list(covariates.columns),
        "kinds": dict(covariates.kinds),
        "levels": {name: list(levels) for name, levels in covariates.levels.items()},
        "reference_levels": dict(covariates.reference_levels),
        "design_columns": list(covariates.column_names),
        "reference_row": dict(zip(covariates.column_names, covariates.reference_row)),
        "lsmean_reference": "pooled_mean",
        "missing_policy": dataset.missing_covariates,
        "n_dropped": len(dataset.dropped_sample_ids),
        "dropped_sample_ids": list(dataset.dropped_sample_ids),
        "coefficients": {
            name: dict(zip(components, coefficients[index])) for index, name in enumerate(covariates.column_names)
        },
    }


def _methylation_units(transforms: Mapping[str, str]) -> str:
    return "mvalue" if "beta_to_mvalue" in transforms.values() else "not_applicable"


def _loadings_frame(model: Any, columns: pd.Index, components: list[str]) -> pd.DataFrame:
    rotations = np.asarray(model.x_rotations_, dtype=float)
    split = [str(column).split("__", 1) for column in columns]
    frame = pd.DataFrame(rotations, columns=components)
    frame.insert(0, "feature", [parts[1] for parts in split])
    frame.insert(0, "block", [parts[0] for parts in split])
    return frame


def _stage_means_frame(
    ls_means: np.ndarray,
    betas: np.ndarray,
    group_levels: list[str],
    stage_levels: list[str],
    components: list[str],
) -> pd.DataFrame:
    means = np.asarray(ls_means, dtype=float) @ betas
    frame = pd.DataFrame(means, columns=components)
    frame.insert(0, "stage", [stage for _ in group_levels for stage in stage_levels])
    frame.insert(0, "group", [group for group in group_levels for _ in stage_levels])
    return frame


def _software_versions() -> dict[str, str]:
    import scipy
    import sklearn

    from motco import __version__

    return {
        "motco": __version__,
        "numpy": np.__version__,
        "scipy": scipy.__version__,
        "sklearn": sklearn.__version__,
        "pandas": pd.__version__,
        "python": platform.python_version(),
    }


def _json_safe(value: Any) -> Any:
    """Convert numpy scalars/arrays and tuples; non-finite floats become ``None``."""

    if isinstance(value, Mapping):
        return {str(key): _json_safe(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [_json_safe(item) for item in value]
    if isinstance(value, np.ndarray):
        return _json_safe(value.tolist())
    if isinstance(value, (bool, np.bool_)):
        return bool(value)
    if isinstance(value, (int, np.integer)):
        return int(value)
    if isinstance(value, (float, np.floating)):
        numeric = float(value)
        return numeric if np.isfinite(numeric) else None
    return value


__all__ = [
    "ANALYSIS_SCHEMA_VERSION",
    "MISSING_COVARIATE_POLICIES",
    "OmicsBlock",
    "TrajectoryAnalysisError",
    "TrajectoryAnalysisParams",
    "TrajectoryAnalysisResult",
    "TrajectoryDataset",
    "analyze_trajectories",
    "write_analysis_outputs",
]
