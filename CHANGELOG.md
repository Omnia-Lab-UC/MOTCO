# Changelog

All notable changes to the `motco` package are recorded here. The format follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/). Before 1.0, every feature change bumps the minor version.

Some versions were set on `main` but never tagged or published to PyPI. They are marked *not published*, with the release that first shipped them. Their numbers still matter, because result records carry the `motco` version that produced them.

Every opt-in option listed below leaves all outputs byte-identical while it is unset. Validation studies, their drivers under `scripts/` and their findings are not listed here: see `docs/reports/` and the [roadmap](docs/roadmap.md).

## [Unreleased]

## 0.24.0 - 2026-11-12

### Added

- `SimulationEvaluationParams.latent_coordinates` takes three predictive-subspace values: `predictive`, `predictive_orthonormal` and `predictive_fitted`. Each measures the final PLS model on the row space of its Y loadings, which has `min(rank, n_stages − 1)` dimensions. The values are the scores projected on that subspace, the input coordinates on an orthonormal basis of it, and the fitted stage indicators. The fit, the cross-validation and the rank are unchanged. The values are rejected with integration methods other than `pls`, with attribution and with fewer than three stages.
- `pls_coordinates_from_pls` builds any non-score coordinate system from an existing score-coordinate `pls` integration without refitting. `pls_orthonormal_from_pls` is now its orthonormal case.
- `LinearLatentMap.from_pls_predictive`, `predictive_directions` and `PREDICTIVE_COORDINATES` in `motco.simulations.diagnostics`.
- `TrajectoryAnalysisParams.latent_coordinates` takes the same three values in `motco.analysis`. The record gains `measurement.latent_coordinates`, which gives the value and the measured dimension. The field is rejected with attribution, with the split interaction test and with fewer than three stages. It has no CLI flag.

## [0.23.0] - 2026-11-11

### Added

- `interaction_feature_scan(..., blocks=)` runs the Westfall–Young step-down within each block and combines the tested blocks by Bonferroni (`INTERACTION_SCAN_COMBINATIONS`: `joint` | `block_bonferroni`). A single block label reproduces the joint scan exactly.
- `TrajectoryAnalysisParams.interaction_scan_by_block` and `motco analyze --interaction-scan-by-block` pass the declared blocks as the scan's labels. They require the interaction scan.

### Changed

- `pandas-stubs` and `scipy-stubs` are no longer runtime dependencies; they moved to the `test` extra.

### Removed

- `scripts/compare_lsmeans.py`, `scripts/lsmeans_obs.R` and `scripts/time_estimate_difference.py`, which no longer ran after the `motco.stats.sd` split in 0.3.1.

## 0.22.0 - 2026-11-11

*Not published; first released in 0.23.0.*

### Added

- `motco.stats.sparse_pls.SparsePLSDA`, a sparse stage-PLS-DA that soft-thresholds the X weights within each feature group. At the full budget it calls `PLSRegression.fit` unchanged. It is checked against mixOmics 6.34.0 `splsda`.
- `SparsityDesign` (scope `joint` | `block`, default keep-fraction ladder 1.0 / 0.5 / 0.2 / 0.1 / 0.05) and `plsda_doubleCV(..., sparsity=)`, which cross-validates rank and keep fraction together and adds `keep_fraction` to the table.
- `TrajectoryAnalysisParams.sparsity` (`SparsitySettings`), `motco analyze --sparsity` and `--sparsity-ladder`, and the harness option `integration_params["sparsity"]`.

## 0.21.0 - 2026-11-11

*Not published; first released in 0.23.0.*

### Added

- Label-blind block weighting of the standardized blocks (`BLOCK_WEIGHTINGS`: `equal_variance` = 1/√p_b, `mfa` = √n/σ₁(Z_b)). It is available as `TrajectoryAnalysisParams.block_weighting`, `motco analyze --block-weighting` and `integration_params["block_weighting"]`. The latent space is fitted in the weighted metric; attribution stays in original units.

## [0.20.0] - 2026-11-10

### Added

- `motco.stats.interaction.interaction_feature_scan` tests each feature's group × stage interaction F under Freedman–Lane permutation. It gives a max-T global p-value and Westfall–Young step-down adjusted p-values.
- `TrajectoryAnalysisParams.interaction_scan` and `motco analyze --interaction-scan` add an `interaction_scan` record block and `interaction_scan.csv`.

## 0.19.0 - 2026-11-09

*Not published; first released in 0.20.0.*

### Added

- `latent_coordinates="orthonormal"` in the harness measures on an orthonormal basis of the fitted PLS span, with the same fit, rank and permutation draws. It is a diagnostic: production stays on the PLS scores.

## [0.18.0] - 2026-11-05

### Added

- `motco.stats.interaction.interaction_split_test` tests the group × stage interaction in total, inside the fitted latent span and in its complement.
- `TrajectoryAnalysisParams.span_split_test` and `motco analyze --span-split-test` add a `span_split` record block.

## 0.17.0 - 2026-11-05

*Not published; first released in 0.18.0.*

### Added

- `span_scale="planted"` for the stage-span construction makes the planted edit exact under the planted pooled standard deviation.

## 0.16.0 - 2026-11-05

*Not published; first released in 0.18.0.*

### Added

- `span_alignment` and `span_support` tilt the stage-span orientation plant toward a direction the null stage program does not use.
- Span capture through `LinearLatentMap` for `pca`, `pls_residual_pca` and standardized `concat` as well as PLS.
- The diagnostic integration method `pls_residual_pca`: the production fit plus stage-residual components orthogonal to its span.

## 0.15.0 - 2026-11-04

*Not published; first released in 0.18.0.*

### Added

- The stage-span construction (`construction="stage_span"`, `span_layers`) plants magnitude, orientation, shape or translation edits directly on group B's null stage means.
- `span_capture` in the harness records how much of the population group difference and of the stage means falls inside the fitted span.

## [0.14.0] - 2026-11-02

### Added

- Covariate-adjusted orientation attribution: `analyze_orientation_attribution(..., covariates=, covariate_reference=)` attributes from adjusted cell means and refits the covariate coefficients in each bootstrap resample.

### Changed

- `motco analyze --attribution` now works with covariates and adds an `adjustment` block; before, the combination was rejected.

## [0.13.0] - 2026-10-31

### Added

- `RRPP(..., residual_correction=)` (`df` | `leverage`) rescales the reduced-model residuals before permuting. It is available as `TrajectoryAnalysisParams.residual_correction` and `motco analyze --residual-correction`. Only p-values and the null summary change.

## 0.12.0 - 2026-10-29

*Not published; first released in 0.13.0.*

### Added

- Study `type_i` acceptance entries can be scoped to a `trajectory_mode` and `statistic`.

## 0.11.0 - 2026-10-29

*Not published; first released in 0.13.0.*

### Added

- Study reports gain a stage-signal table (`design_point_stage_signal.csv`) for design-grid runs.

## [0.10.0] - 2026-10-26

### Added

- `plsda_doubleCV(..., min_components=)`, `motco plsr --min-components` (default 1) and `motco analyze --min-components`.

### Changed

- The PLS latent rank is floored at 2 by default in `motco analyze` and the simulation harness. `min_components=1` restores the unfloored rule. Results whose selected rank was 1 change.

## [0.9.0] - 2026-10-20

### Added

- Covariate-aware latent construction: `TrajectoryAnalysisParams.covariate_latent_adjustment` (`stage` | `group_stage`) and `motco analyze --latent-adjustment`. These residualize the features on the covariates inside every cross-validation training set.
- `plsda_doubleCV(..., scale=, nuisance=, strata=)`, `NuisanceDesign`, `fit_nuisance_coefficients` and `residualize`.

## [0.8.0] - 2026-10-16

### Added

- Covariates in `motco.analysis` and `motco analyze` (`--covariate`, `--categorical-covariate`, `--missing-covariates`). They enter the trajectory models; the latent space stays covariate-blind. With `complete_case`, samples with missing covariates are dropped and recorded, never imputed.

## [0.7.0] - 2026-10-15

### Added

- `motco.analysis` (`analyze_trajectories`, `write_analysis_outputs`) and `motco analyze`: the validated pipeline on real data. It runs pooled standardization, a stage-supervised PLS-DA latent space with a double-CV rank, `delta` / `angle` / `shape` and RRPP p-values, with optional orientation attribution.
- `orientation_estimator` (`pc1` | `net_displacement`) on `estimate_difference` and `RRPP`, and `motco de --orientation-estimator` (default `pc1`).
- A numpy reimplementation of InterSIM's generative model, so simulation no longer needs R at runtime. It is validated against InterSIM; see `src/motco/simulations/FIDELITY.md`.
- The `pls` integration method, the production latent space, and the diagnostic methods `pca` and `pls_within_group`.
- Orientation driver attribution (`motco.stats.attribution`).
- The declarative power-study engine `motco.simulations.study`: configs, sharding, merging, reports, design grids and acceptance rules.
- Generator options `magnitude_kind`, `surgery_censoring`, `baseline_continuity` and explicit group × stage cell sizes.

### Changed

- `shape` is now the strict pairwise ordinary Procrustes distance, with reflections aligned away at every dimension. Its values differ from 0.6.0.
- PC1's sign is anchored on the net displacement (last stage minus first) instead of the centred first stage. This can change `angle` where the old anchor was near zero.
- The simulation harness defaults to `magnitude_kind="joint"` and `orientation_estimator="net_displacement"`. The primitives `estimate_difference` and `RRPP` keep PC1 as their default.
- Methylation is converted to M-values at the integration step.
- Shape and orientation estimation use a thin-QR path when there are more dimensions than stages, which is much faster in high dimensions.

### Removed

- The R-backed `generate_semisynthetic_trajectory_from_intersim`. `motco.simulations.intersim` remains as a legacy bridge that nothing calls at runtime.

## [0.6.0] - 2026-05-08

### Added

- InterSIM cluster mean-shift controls in the simulation bridge, and a retuned toy dataset.

### Fixed

- `plsda_doubleCV` aggregates the nested CV by fold means.

## [0.5.0] - 2026-05-07

### Added

- `motco simulate`, which then needed R InterSIM.
- `motco plsr --metadata` and `--out-scores` for PLS-DA latent-space output.
- `motco.viz`: `plot_trajectories` and `plot_trajectory_from_data`.
- An end-to-end example notebook.

## [0.4.0] - 2026-05-06

### Added

- The simulation framework: an InterSIM bridge, a semi-synthetic trajectory generator, an evaluation harness and grid orchestration.

### Fixed

- `get_observed_vectors` no longer returns NaN.

## [0.3.1] - 2026-05-05

### Added

- Input validation with clear errors in `estimate_difference`, `RRPP`, `get_model_matrix`, `center_matrix`, `plsda_doubleCV`, `get_affinity_matrix` and `SNF`.
- `motco --version` and `motco de --out-observed`.
- `calculate_vips` and `motco plsr --out-vips`.
- `get_spectral(n_components=)` and `motco snf --spectral-components`.
- Progress bars for `RRPP` and `plsda_doubleCV`, and logging of numerical fallbacks.
- All public functions are exported from `motco.stats`.

### Changed

- `motco.stats.sd` was split into `design`, `trajectory` and `permutation`. Code importing `motco.stats.sd` must import from those modules or from `motco.stats`.

### Fixed

- Multiprocessing fragility in `RRPP` and the PLS cross-validation.
- A NumPy deprecation warning in `calculate_vips`.

## [0.3.0] - 2025-12-18

### Changed

- The trajectory statistics were refactored for speed and numerical stability.

## [0.2.1] - 2025-12-09

### Added

- Parallel `RRPP` permutations.

## [0.2.0] - 2025-12-09

### Changed

- `get_model_matrix` and the LS-mean construction (`build_ls_means`) take the group and level columns explicitly. They were checked against reference datasets.

## [0.1.1] - 2025-12-05

First published release: PLS-DA with double cross-validation, similarity network fusion, and trajectory size, orientation and shape differences with RRPP. The CLI has `plsr`, `snf` and `de`.

## 0.1.0 - 2025-12-05

*Not published; superseded the same day by 0.1.1.*

[Unreleased]: https://github.com/Omnia-Lab-UC/MOTCO/compare/v0.23.0...HEAD
[0.23.0]: https://github.com/Omnia-Lab-UC/MOTCO/compare/v0.20.0...v0.23.0
[0.20.0]: https://github.com/Omnia-Lab-UC/MOTCO/compare/v0.18.0...v0.20.0
[0.18.0]: https://github.com/Omnia-Lab-UC/MOTCO/compare/v0.14.0...v0.18.0
[0.14.0]: https://github.com/Omnia-Lab-UC/MOTCO/compare/v0.13.0...v0.14.0
[0.13.0]: https://github.com/Omnia-Lab-UC/MOTCO/compare/v0.10.0...v0.13.0
[0.10.0]: https://github.com/Omnia-Lab-UC/MOTCO/compare/v0.9.0...v0.10.0
[0.9.0]: https://github.com/Omnia-Lab-UC/MOTCO/compare/v0.8.0...v0.9.0
[0.8.0]: https://github.com/Omnia-Lab-UC/MOTCO/compare/v0.7.0...v0.8.0
[0.7.0]: https://github.com/Omnia-Lab-UC/MOTCO/compare/v0.6.0...v0.7.0
[0.6.0]: https://github.com/Omnia-Lab-UC/MOTCO/compare/v0.5.0...v0.6.0
[0.5.0]: https://github.com/Omnia-Lab-UC/MOTCO/compare/v0.4.0...v0.5.0
[0.4.0]: https://github.com/Omnia-Lab-UC/MOTCO/compare/v0.3.1...v0.4.0
[0.3.1]: https://github.com/Omnia-Lab-UC/MOTCO/compare/v0.3.0...v0.3.1
[0.3.0]: https://github.com/Omnia-Lab-UC/MOTCO/compare/v0.2.1...v0.3.0
[0.2.1]: https://github.com/Omnia-Lab-UC/MOTCO/compare/v0.2.0...v0.2.1
[0.2.0]: https://github.com/Omnia-Lab-UC/MOTCO/compare/v0.1.1...v0.2.0
[0.1.1]: https://github.com/Omnia-Lab-UC/MOTCO/releases/tag/v0.1.1
