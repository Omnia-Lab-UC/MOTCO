# Proposal

## Why

The joint feature matrix is built by z-scoring every feature and concatenating the blocks unweighted, so a block's weight in every fitted space grows with its feature count: in the simulations methylation (367 features) outweighs expression (131) by 2.8 to 1 in total variance, and the SEA-AD export works around the same thing by forcing 5,000 + 5,000 features. Roadmap item 5 names block weighting as the one preprocessing step that belongs in MOTCO, because it is a property of the joint measurement rather than of one omic; the latent-space comparison (2026-11-04, §12) listed it as a limitation of every space it compared. It is also the first arm of the production latent-space selection (`production-latent-space-selection`), which needs it to exist before any regularized candidate can be attributed to its method rather than to the weights it runs on.

## What Changes

- **An opt-in block weighting in the shared preprocessing.** After pooled per-feature standardization, every feature of block *b* is multiplied by one label-blind scalar *w_b*, fitted once on the pooled sample like the scalers and applied to the observed rows and to the population means alike. Two kinds:
  - `equal_variance`: *w_b* = 1/√*p_b*, so every block has unit total variance (the multi-block-PLS block scaling; after per-feature standardization it is the same thing as "unit total block variance", the roadmap's other candidate);
  - `mfa`: *w_b* = √n/σ₁(*Z_b*) = 1/√λ₁, the inverse standard deviation of the standardized block's leading principal component, so every block's dominant direction has unit variance (multiple-factor-analysis scaling).

  Absent (the default), weights are 1 and every output is byte-identical.
- **PLS fits in the weighted metric.** When weighting is set, every PLS fit of the production latent space (cross-validation, outer refits, the final model) runs with `scale=False`, as the latent-adjustment path already does, because `scale=True` would re-standardize each column and undo the weights exactly. Without weights the `scale=True` path is untouched.
- **Harness setting.** `integration_params["block_weighting"]` (`None` | `"equal_variance"` | `"mfa"`; an unknown value rejected before any preprocessing) applies to every integration method that reads the shared standardized blocks (`snf` and the unstandardized `concat` do not, and refuse it) and is recorded in the integration metadata only when set. Because it is an integration parameter, a study configuration can name it as the nested axis `evaluation.integration_params.block_weighting` with no engine change, and absent keys are already left out of cell ids and signatures.
- **Real-data setting.** `TrajectoryAnalysisParams.block_weighting` (same values; omitted from `to_dict()` while `None`) and `motco analyze --block-weighting {equal_variance,mfa}`. The record gains a `block_weighting` block after `blocks` (kind, per-block weight and feature count) and `latent.integration_params.block_weighting`, only when set. Attribution's original-scale mapping and `latent_artifacts.original_scales` include the weight, so reconstructed features are in their original units.
- **Documentation.** `docs/api/simulations.md`, `docs/api/analysis.md`, `docs/api/cli.md`, roadmap item 5, `AGENTS.md`. Version 0.21.0.

**Non-goals.**
- **No Type I or power run.** Roadmap item 5 requires a size and power check before the case study may use weighting; that check is the `weighted` arm of `production-latent-space-selection`, measured against production on matched data. Until that reading, the setting is available and documented as unvalidated.
- **No default change.** Flipping the default follows only from the selection's verdict and is its own change, with pins for every committed configuration.
- **Per-feature or data-driven weights beyond the two kinds**, and weights fitted with any label.
- **The SEA-AD export's feature split** stays `sead_sexdiff`'s decision.

## Capabilities

### New Capabilities
- `block-weighting`: the label-blind per-block weighting of the pooled-standardized blocks: its kinds and their definitions, that it is fitted once on the pooled sample and applied to observed rows and population means alike, its invisibility while absent, and the PLS metric it implies.

### Modified Capabilities
- `simulation-evaluation-harness`: the preprocessing accepts the weighting through `integration_params`, every integration method runs on the weighted blocks, the PLS space fits in the weighted metric, the metadata records the kind only when set, and an unknown kind is rejected before any preprocessing.
- `real-data-trajectory-analysis`: the analysis offers the same setting, records it only when set, maps attribution and the exposed original scales through the weight, and `motco analyze` exposes it.

## Impact

- **Code:**
  - `src/motco/simulations/preprocessing.py`: `BlockScaler.weight` (default 1.0), `fit_block_scalers(..., weighting=)`, `transform_blocks` multiplying by the weight, `FittedOmicsPreprocessor.block_weighting`, `fit_omics_preprocessor(..., weighting=)`, the weighting kinds constant;
  - `src/motco/simulations/evaluation.py`: reading the key, passing it to the preprocessor, `scale=False` in `fit_pls_latent_space` when set, metadata, `_joint_original_scales` through the weight;
  - `src/motco/analysis.py`: the parameter, `to_dict`, the integration-params key, the record block, original scales;
  - `src/motco/cli.py`: `--block-weighting`;
  - `src/motco/__init__.py`: 0.21.0.
- **Tests:** weight values for both kinds on a known block; byte identity while absent (the committed cell ids, signatures, the core/harness identity fixtures, `tests/test_analysis_identity.py`); that `equal_variance` gives every block unit total variance and `mfa` unit leading variance on the transformed blocks; that the population means are weighted like the observed rows; that a weighted PLS fit differs from a `scale=True` fit on the same rows (the weights survive); that the rank and CV table of a weighted fit with weights of 1 equal production's (the `scale=False` metric selects the same rank); record and `to_dict` contents; CLI round trip; attribution's reconstructed features in original units under weighting.
- **Compatibility:** every committed configuration, fixture, record and report is unchanged; the setting is absent everywhere.
- **Compute:** none beyond the test suite.
