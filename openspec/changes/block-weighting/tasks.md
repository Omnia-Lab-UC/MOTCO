# Tasks

## 1. Weighting in the shared preprocessing

- [ ] 1.1 Add `BLOCK_WEIGHTINGS = ("equal_variance", "mfa")`, `BlockScaler.weight` (default 1.0) and `fit_block_scalers(..., weighting=None)` computing `1/√p_b` or `1/σ₁(Z_b)` on the standardized block (design D1); reject an unknown kind naming it and the allowed kinds before transforming anything. `transform_blocks` multiplies by the weight. Verify in `tests/test_block_weighting.py` that `equal_variance` gives each transformed block unit total variance and `mfa` unit leading singular value on a two-block fixture, that unknown kinds are rejected, and that `transform_blocks` without a weighting is byte-identical to a pinned pre-change output.
- [ ] 1.2 Add `FittedOmicsPreprocessor.block_weighting` and `fit_omics_preprocessor(dataset, layers, weighting=None)`; `transform_dataset` and `transform_population` go through the weighted scalers. Verify that the population means are weighted like the observed rows (each block equals the unweighted transform times the weight) and that `require_layers` is unchanged.

## 2. Harness

- [ ] 2.1 Read `integration_params["block_weighting"]` in `integrate_semisynthetic_dataset` (and wherever a preprocessor is fitted for an evaluation), validate it before any preprocessing, pass it to `fit_omics_preprocessor`, and refuse a prebuilt preprocessor whose weighting differs from the request, naming both. Verify one test per rejection and that `concat`, `pca` and `pls` evaluations with `equal_variance` measure weighted blocks (the standardized `concat` matrix has unit total variance per block).
- [ ] 2.2 In `fit_pls_latent_space`, switch every PLS fit to `scale=False` when the key is set and record `cv_params["block_weighting"]` (design D2); route `_joint_original_scales` through `scale / weight`. Verify that a weighted fit's `_x_std` is unity and its first weight vector differs from a `scale=True` fit's on the same rows; that a `scale=False` fit with unit weights selects the same rank and CV table as production on a pinned dataset; and that the committed-config identity tests (`tests/test_committed_config_identity.py`, the pin tests, the core/harness identity fixtures) pass unchanged.
- [ ] 2.3 Verify under `span_capture` that `population_latent` and the span-capture shares are computed on the weighted inputs (the population checkpoint of a weighted `pls` evaluation equals the unweighted one only when the weights are 1), and that an evaluation without the key is byte-identical to the pre-change result on a pinned dataset.

## 3. Real-data analysis and CLI

- [ ] 3.1 Add `TrajectoryAnalysisParams.block_weighting` (validated in `__post_init__`, popped from `to_dict()` while `None`), pass the key into the analysis's `integration_params`, fit the scalers with the kind, and route attribution's original scales and `latent_artifacts.original_scales` through the weight (design D5). Verify in `tests/test_analysis_block_weighting.py` that attribution's reconstructed features of a weighted analysis are in original units (equal to the unweighted reconstruction of the same latent means), that an unknown kind is rejected before fitting, and that `tests/test_analysis_identity.py` passes unchanged.
- [ ] 3.2 Write the `block_weighting` record block after `block_units` (kind; per block weight and feature count) and `latent.integration_params.block_weighting`, only when set. Verify the record contents and their absence without the setting.
- [ ] 3.3 Add `motco analyze --block-weighting {equal_variance,mfa}` and print the kind beside the selected rank when given. Verify a CLI round trip on the test dataset writes `block_weighting.kind`, and that the help text names both kinds.

## 4. Documentation and version

- [ ] 4.1 Document the setting in `docs/api/simulations.md` (preprocessing and `integration_params`), `docs/api/analysis.md` (parameter, record block, `delta` in the weighted metric, "unvalidated until the production latent-space selection reads it") and `docs/api/cli.md`; update roadmap item 5 (the setting exists; the check is `production-latent-space-selection`) and the `preprocessing`/`evaluation.py`/`analysis.py`/`cli.py` bullets of `AGENTS.md`. Verify with `uv run mkdocs build` (no new warnings).
- [ ] 4.2 Bump `src/motco/__init__.py` to 0.21.0 and run the pre-commit gate: `uv run ruff check src/ tests/ && uv run mypy src/motco/ && MOTCO_TEST_PERMS=99 uv run pytest tests/ -m "not slow" --tb=short`.

## Workflow follow-up

- Sync the delta specs and archive the change once merged.
- The size and power check of weighting is the `weighted` arm of `production-latent-space-selection`; a default flip follows only from its verdict, as its own change.
