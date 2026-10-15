# Design

## Context

See proposal.md for motivation. These facts shape the approach:

- The validated pipeline lives in `simulations/evaluation.py` and `simulations/preprocessing.py` and is
  bound to `SemiSyntheticTrajectoryDataset`: `OMIC_LAYERS = ("methylation", "expression", "proteomics")`
  is the canonical block order, `OmicsLayer` is a `Literal` of those names, `fit_omics_preprocessor` and
  `FittedOmicsPreprocessor.transform_dataset` read blocks with `getattr(dataset, layer)`, and
  `_integration_values` applies `logit` to the block *named* `methylation`. `concatenate_blocks` takes
  canonical order and prefixes `f"{layer}__{column}"`. `_pls_integration` holds the production
  measurement: CV-knob clamping, the feasibility check, `plsda_doubleCV`, `_modal_int_with_parsimony`,
  one `fit_plsda_model`, the metadata record, and `PLSIntegrationArtifacts` for attribution.
- `build_simulation_trajectory_design` orders groups and stages with `design._sort_levels` (numeric when
  every label is an integer, else lexicographic); `get_model_matrix` and `build_ls_means` do the same and
  take no explicit level order. The harness's integer stage labels sort correctly; real stage labels do not.
- `estimate_difference`, `RRPP`, `_extract_pair_statistics`, `_empirical_p_value`,
  `_summarize_null_distributions`, `_config_spectrum_block` and `_summarize_permutation_eigengaps` are
  already dataset-free. `calculate_realized_geometry` and `derive_truth_driver_features` read generator
  truth (`population_trajectories`, `truth["indicators"]`) and degrade gracefully when it is absent.
- `analyze_orientation_attribution` takes the standardized joint matrix, row-aligned metadata, the fitted
  estimator, and explicit `groups=`/`stages=` orders validated against the observed labels.
  `compute_attribution_diagnostics` wraps it into the compact record but takes the synthetic dataset for
  truth derivation.
- Identity is pinned by `tests/test_orientation_probe_byte_identity.py`, `tests/test_phase6_byte_identity.py`
  (recorded statistics), `tests/test_committed_config_identity.py`, `tests/test_magnitude_default_pins.py`
  and `tests/test_orientation_default_pins.py` (signatures, cell ids, seeds), and the R-reference regression
  tests (primitives). Signatures hash parameters, not code, so an internal refactor that preserves outputs
  is safe and these tests prove it.
- `cmd_de` calls `estimate_difference`/`RRPP` with no estimator argument; the primitives default to PC1 and
  `trajectory-orientation-invariance` requires that default to stay.
- The Phase 6 design has n = 80 with a minimum stage count of 20 (10 + 10 at Intermediate); the harness's
  CV defaults (`cv2_splits=4`, `cv1_splits=3`, `n_repeats=5`, `max_components=20`, `random_state=1203`)
  were what the small-n operating study ran, after the same clamping the core applies.

## Goals / Non-Goals

**Goals:**

- One implementation of the production measurement, called by both the harness and the real-data path,
  with a test that the two agree on identical data and the existing identity tests proving the harness did
  not move.
- An input contract that cannot silently measure the wrong thing: stage order is explicit, block units are
  declared, alignment and finiteness are validated, and the record states what was measured.
- Everything the Phase 6 interpretability statement needs beside a result — selected rank, estimator, null
  summary, spectrum — is in the record without the caller asking.

**Non-Goals:**

- A general-purpose preprocessing library (count normalization, filtering, batch correction, imputation).
- Covariate-adjusted designs. Out of scope per the proposal; see R3.
- Exposing `concat`, `pca`, `snf` or `pls_within_group` on real data. They are simulation diagnostics.
- Serializing the sklearn estimator. The artifact stores scores, loadings and the fitted contract instead.
- Moving `OMIC_LAYERS`/`OmicsLayer` out of the simulation modules. The synthetic path keeps its fixed names.

## Decisions

### D1. New top-level module `motco.analysis`, not a submodule of `simulations`

Real-data analysis is a product surface, not a simulation. `src/motco/analysis.py` holds the public API:

- `OmicsBlock(name, values: pd.DataFrame, transform: Literal["identity", "beta_to_mvalue"] = "identity")`.
- `TrajectoryDataset(blocks: Sequence[OmicsBlock], metadata, *, sample_col, group_col, stage_col,
  stage_order, group_order=None)` — frozen dataclass with validation in `__post_init__`
  (alignment, finiteness, duplicate block names, beta range for `beta_to_mvalue`, stage-order
  completeness, exactly two groups, complete group × stage table).
- `TrajectoryAnalysisParams(permutations=0, seed=None, n_jobs=1, orientation_estimator="net_displacement",
  forced_components=None, cv1_splits=3, cv2_splits=4, n_repeats=5, max_components=20,
  random_state=1203, attribution: AttributionDiagnosticSettings = disabled, progress=False)`.
  Reusing `AttributionDiagnosticSettings` keeps the knobs identical to the harness's.
- `analyze_trajectories(dataset, params) -> TrajectoryAnalysisResult`.
- `write_analysis_outputs(result, out_dir, *, overwrite=False) -> dict[str, Path]`.
- `TrajectoryAnalysisResult` with `record` (JSON-safe dict), `latent_scores`, `pls_loadings`,
  `stage_means`, `attribution` (the `OrientationAttributionResult` or `None`) and `to_json()`.

Alternative considered: extend `SimulationEvaluationParams`/`evaluate_semisynthetic_trajectory` with a
dataset protocol. Rejected: the harness API is hashed into cell signatures and pinned by five identity
tests; growing it for a real-data concern risks the thing the tests protect and would still leave the
synthetic-only fields (`truth_metadata`, `realized_geometry`) in a real-data result.

### D2. Dataset-free cores extracted from the harness; the harness becomes a thin caller

Two behaviour-preserving extractions, each gated by the identity tests:

- `preprocessing.fit_block_scalers(blocks: Mapping[str, pd.DataFrame], transforms: Mapping[str, str])
  -> dict[str, BlockScaler]` and `transform_blocks(blocks, scalers, transforms)`; `concatenate_blocks`
  gains an `order: Sequence[str] | None` parameter (default: canonical `OMIC_LAYERS` filtering, exactly as
  today). `fit_omics_preprocessor`/`transform_dataset` call the cores with `transforms = {"methylation":
  "beta_to_mvalue"}` and the dataset's blocks, so the synthetic path is byte-identical.
- `evaluation.fit_pls_latent_space(X: pd.DataFrame, y: pd.Series, integration_params, *,
  sample_ids, layer_feature_counts, layers_metadata) -> LatentIntegrationResult` — the body of
  `_pls_integration` from the stage-count check onward, with the dataset reads hoisted into
  `_pls_integration`, which becomes a six-line adapter. The metadata dict and `PLSIntegrationArtifacts`
  are produced by the core, so the real-data result records the same `selected_lv`,
  `component_selection`, `cv_mean_auroc`, `cv_repeats_completed`, `feature_order_signature` and CV knobs.
  `original_scales` are computed by a small helper that takes a scaler mapping instead of the preprocessor.

Alternative considered: re-implement the glue in `analysis.py` and rely on the equivalence test alone.
Rejected: two copies of the clamping and rank rule is precisely the drift the proposal exists to prevent;
the equivalence test stays as a second guard, not the only one.

### D3. Explicit stage order is honoured by ordinal recoding, not by changing `design.py`

`get_model_matrix`/`build_ls_means` sort levels. Rather than add an order parameter to the primitives (a
public-contract change in `module-contracts` for a single caller), `analyze_trajectories` builds a working
metadata frame whose group and stage columns are recoded to zero-padded ordinal strings
(`"0"`, `"1"`, … in the declared order), calls `build_simulation_trajectory_design` on it, and maps the
returned `group_levels`/`stage_levels` back to the original labels for every report. The LS-mean row order
is then group-major, stage-minor in *declared* order, which is what net displacement, `shape` and the
spectrum need. The equivalence test exercises this on the synthetic dataset whose labels are already
ordinal, and a dedicated test uses non-sorting labels (`Low` < `Intermediate` < `High`) and checks that
`angle` equals the primitives' answer on explicitly ordered LS means.

Attribution receives the *original* metadata with `groups=`/`stages=` set to the declared orders, since it
supports explicit levels natively, so its tables carry the caller's labels directly.

Alternative considered: require the caller to pre-code stages as integers. Rejected: it moves the
load-bearing step out of the validated code and into every case-study script.

### D4. Record layout mirrors the harness's `SimulationEvaluationResult`, minus the synthetic blocks

`record` keys: `schema_version` (1), `statistics` (pair statistics), `p_values`, `null_summary`,
`config_spectrum` (incl. `permutation_pooled_eigengap` when permutations ran), `latent` (the integration
metadata dict from the core: `selected_lv`, `component_selection`, `cv_*`, `integration_params`,
`feature_order_signature`, `layer_feature_counts`), `design` (`group_levels`, `stage_levels`,
`cell_sizes` as nested `{group: {stage: n}}`, `contrast`), `blocks` (`[{name, n_features, transform,
units}]`), `measurement` (`orientation_estimator`, `permutations`, `seed`, `n_jobs`, `shape_available`,
`runtime_seconds`), `params` (the resolved `TrajectoryAnalysisParams` as a dict), `software`
(`motco`, `numpy`, `scipy`, `sklearn`, `pandas`, `python`), and `attribution` (the compact record, or
`unavailable_record("not requested")`). No `truth_metadata`, no `realized_geometry`: a real-data result
must not carry fields that only make sense with generator truth.

`compute_attribution_diagnostics` is generalized to take `truth: TruthDrivers | None` (default: derive
from a synthetic dataset as today when a dataset is passed; `analysis.py` passes an explicit unavailable
`TruthDrivers` with the reason "real data carry no generator truth"). The synthetic call site is unchanged.

### D5. The written artifact is plain CSV/JSON; the estimator is not pickled

`write_analysis_outputs` writes `result.json`, `latent_scores.csv` (index = sample ids, columns `pls_0..`),
`pls_loadings.csv` (columns `block`, `feature`, `pls_0..` from `model.x_rotations_`, rows in joint order),
`stage_means.csv` (columns `group`, `stage`, `pls_0..`, rows group-major in declared stage order; the
input for `viz.plot_trajectories`), and under `attribution/` the frames from
`motco.stats.attribution.write_attribution_outputs` when attribution ran. Refuses to overwrite an existing
`result.json` unless `overwrite=True`. Loadings plus the recorded scalers' contract (feature order
signature, units) are enough to audit and to re-project; a pickled `PLSRegression` is version-fragile and
the roadmap item asks for a *versioned* artifact.

### D6. `motco analyze` parses `NAME=PATH` and `NAME=KIND` pairs; `motco de` gains one flag

`--block` values are split on the first `=`; a value without `=` or a repeated name is a usage error. CSVs
are read with the first column as the index (sample id). `--transform` names must match declared blocks.
`--stage-order` is required and comma-split; `--group-order` is optional. CV knobs are exposed as
`--cv1-splits`, `--cv2-splits`, `--n-repeats`, `--max-components`, `--cv-random-state` with the harness
defaults. `ValueError` (including `TrajectoryAnalysisError`, a `ValueError` subclass) becomes
`SystemExit("Error: …")` as in the other subcommands. On success the command prints a short summary and
the output paths.

`cmd_de` gains `--orientation-estimator {pc1,net_displacement}` (default `pc1`) passed to both
`estimate_difference` and `RRPP`, and the JSON gains `"orientation_estimator"`. Adding a key to the JSON
is the one visible change to default `de` output; the spec scenario "default output is unchanged" is read
as unchanged *statistics*, and the test compares `deltas`/`angles`/`shapes`. (If strict byte identity of
the JSON is preferred, the key is written only when the flag is given — decided at implementation per the
existing CLI test's expectations; either reading satisfies the roadmap's purpose.)

### D7. Version bump to 0.7.0

New public module and subcommand → minor bump in `src/motco/__init__.py`. Tagging is the maintainer's.

## Risks / Trade-offs

- [Refactoring `_pls_integration` could perturb recorded statistics] → The byte-identity fixtures
  (`orientation_probe`, `phase6`) and the committed-config identity tests run in the fast suite; the
  extraction is done first, alone, and gated before any new code lands (tasks 1.x).
- [Ordinal recoding hides a mismatch between recoded and original labels] → The design builder's returned
  levels are mapped back through the same table used to recode, and a test asserts the reported
  `stage_levels` equal `stage_order` and the LS-mean row order matches `cell_sizes`.
- [Small-n CV infeasibility on a real cohort] → Same clamping and the same explicit
  `SimulationEvaluationError` text as the harness, surfaced as `TrajectoryAnalysisError` with the stage
  counts, so the caller sees why rather than a silent rank.
- [A caller passes raw counts with `identity`] → The record states each block's declared units, the docs
  state that MOTCO does not normalize, and the CLI help says so; the operating study's blocks were
  continuous, so this is a documented precondition, not something MOTCO can check.
- [R3: no covariates] → Documented limitation in `docs/api/analysis.md` and the roadmap: the operating
  characteristics were measured on a group × stage design, and adjusting for covariates would need its own
  validation before any claim.
- [`shape` with net displacement reads only `angle`; `shape` and `delta` unchanged] → Record states the
  estimator; the equivalence test checks `delta`/`shape` match the harness regardless of estimator.

## Migration Plan

Additive. No committed `results/` or `docs/reports/` change. Existing `motco de` scripts run unchanged.
The case-study repository adopts `motco analyze` (or `analyze_trajectories`) pinned to the release that
contains this change. Rollback is removing the module and subcommand; the refactor cores can stay.
