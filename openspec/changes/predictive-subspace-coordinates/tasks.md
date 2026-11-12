# Tasks

## 1. Predictive-subspace coordinates in the harness

- [x] 1.1 Extend `LATENT_COORDINATES` and `SimulationEvaluationParams.latent_coordinates` with `"predictive"`, `"predictive_orthonormal"` and `"predictive_fitted"` (design D2, D4). `"scores"` still normalizes to `None`, and an unknown value is still rejected in `__post_init__` naming it and the allowed values. Verify in `tests/test_predictive_coordinates.py` the three accepted values and the rejection, and that `tests/test_orthonormal_coordinates.py`, `tests/test_committed_config_identity.py`, `tests/test_orientation_default_pins.py`, `tests/test_rank_floor_default_pins.py` and `tests/test_post_floor_config_identity.py` pass unchanged.
- [x] 1.2 Add the predictive step to `fit_pls_latent_space(coordinates=)` after the final fit (design D1, D2), with a `LinearLatentMap` constructor in `diagnostics.py` whose basis is `qr(x_rotations_·V)` and whose `transform` is in the measured coordinates:
  - `V` from the SVD of `y_loadings_`, numerical rank at 1e-10 of the largest singular value, an error naming the rank when it is below min(fit rank, stages − 1);
  - `predictive` is `scores·V`, `predictive_orthonormal` is `inputs(X)·basis`, `predictive_fitted` is `scores·y_loadings_ᵀ`, on the training rows of the latent-adjusted path too;
  - `integration_params.latent_coordinates` recorded; the scores and orthonormal paths byte for byte as before.

  Verify in the same test file, on a three-stage dataset at forced ranks 2 and 5:
  - the column counts and ranks (2, 2 and 3 of rank 2);
  - the fitted stage indicators recomputed from the subspace equal the model's to 1e-8, and `predictive_fitted` equals the model's predictions centred and rescaled;
  - the three matrices are fixed linear images of one another to 1e-8;
  - the rank and CV table equal the score measurement's;
  - `tests/test_analysis_identity.py` and the core/harness identity fixtures pass unchanged.
- [x] 1.3 Reject a predictive value before any preprocessing or fit (design D4) with a method other than `pls`, with attribution, and on a dataset with fewer than three stages, each error naming the offending value. Verify one test per rejection, and that the existing orthonormal rejections still fire.
- [x] 1.4 Verify the measurement-level contracts in the same test file:
  - at a fit rank of 2 on three stages, `predictive` gives `delta`/`angle`/`shape` within 1e-10 of the scores and equal RRPP p-values;
  - rotating the `predictive_orthonormal` matrix by a random orthogonal matrix leaves the statistics within 1e-10 and the p-values equal;
  - under `span_capture` the three values report identical shares, each at most the score measurement's, and `population_latent` is in the measured coordinates;
  - the `pls_latent` checkpoint and the configuration spectrum come from the predictive matrix;
  - a forced rank, block weighting, sparsity and the latent adjustment each run under `predictive`;
  - an evaluation without the setting is byte-identical to the pre-change result on a pinned dataset.
- [x] 1.5 Generalize the prebuilt-fit path (spec: prebuilt fit): a helper that builds the integration of any non-score coordinate value from a score-coordinate `pls` integration, with `pls_orthonormal_from_pls` kept as its orthonormal case. It refuses a non-`pls` integration, one without artifacts, or one not in score coordinates; `evaluate_semisynthetic_trajectory(..., latent=)` refuses a prebuilt integration whose coordinates differ from its params, naming both. Verify that:
  - each predictive value built from a prebuilt fit equals a from-scratch fit (matrix and metadata) and fits no PLS model (patched `fit_plsda_model` call count);
  - the three values from one prebuilt fit make one fit for four measurements;
  - each refusal fires.
- [x] 1.6 Document the three values in `docs/api/simulations.md` (definition, measured dimension, rejections, invisibility while absent, the prebuilt path) and in the `evaluation.py` bullets of `AGENTS.md`. Verify that `uv run mkdocs build` adds no warnings beyond those on `main`.

## 2. The setting in the real-data analysis

- [x] 2.1 Add `TrajectoryAnalysisParams.latent_coordinates` (default `None`; the three predictive values; anything else, including `"orthonormal"` and `"scores"`, rejected in `__post_init__` naming the value and the allowed values), omitted from `to_dict()` while unset, passed to `fit_pls_latent_space(coordinates=)` (design D4). Reject it before fitting with attribution, with `span_split_test` and on fewer than three stages, naming both settings or the stage count. Verify in `tests/test_analysis_predictive_coordinates.py` each rejection, the omission from `to_dict()`, and that `tests/test_analysis_identity.py`, `tests/test_attribution_identity.py` and the other analysis identity tests pass unchanged.
- [x] 2.2 Record `measurement.latent_coordinates` (the value and the measured dimension) and `latent.integration_params.latent_coordinates` when set; write `latent_scores.csv` and `stage_means.csv` in the measured coordinates and `pls_loadings.csv` from the fitted model. Verify in the same test file:
  - the record block and its absence without the setting;
  - on a three-stage analysis at a forced rank of 5, two score columns under `predictive`, the same selected rank and AUROC as without the setting, and a five-component loadings table;
  - a `test_matches_harness_*` case under `predictive`: identical latent matrices, statistics, p-values, rank and null summary to the harness on the same data;
  - covariates, `residual_correction` and a forced rank each run under the setting.
- [x] 2.3 Document the field in `docs/api/analysis.md` and add its row to the opt-in table of `AGENTS.md`, stating that it has no CLI flag until a verdict. Bump `src/motco/__init__.py` to 0.24.0 and add the `CHANGELOG.md` entry (package changes only). Verify with the pre-commit gate: `uv run ruff check src/ tests/ && uv run mypy src/motco/ && MOTCO_TEST_PERMS=99 uv run pytest tests/ -m "not slow" --tb=short`, including `tests/test_release_notes.py`.

## 3. Pilot driver: units and records

- [x] 3.1 Create `scripts/predictive_subspace_check.py` with `run --profile simulation|weak_size|plasmode`, `merge --profile` and `decide` (design D9): `settings.json` bound by hash; round-robin shards; resume of complete unit keys; refusal of a resume or merge under other settings, naming the keys; the parents' sha256 recorded, with a refusal naming any changed file. Verify in `tests/test_predictive_subspace_check.py` the unit counts (7,800; 1,000; 800), resume, the settings-mismatch refusal and a refused parent whose hash differs.
- [x] 3.2 Implement the `simulation` builder (design D6) over `latent_space_comparison`'s `design_generators`, `generator_seed`, `rrpp_seed`, `unit_dataset` and `cell_spec`, imported unchanged, at base seed 2036 with `span_scale="planted"`: the seven designs, the eleven cells and the surgery cell on its own design entry. Verify that:
  - the designs' generator parameters differ from the comparison's only in the stage effect, the planted scale and, for the surgery design, the construction;
  - the weak δ equals the committed calibration's;
  - no dataset seed equals the comparison's for any shared design name, cell and replicate;
  - cells of one design and replicate share their null dataset.
- [x] 3.3 Implement the `weak_size` builder over `span_split_check` / `latent_space_comparison`'s `seaad_weak` anchor and translation units at replicates 0–499 and 999 permutations, and the `plasmode` builder over `production_latent_space_check` and `plasmode_check` (export contract, `paper` stream, `span_plant`, the analysis-parameter builder, the fidelity values), all imported unchanged: the anchor at replicates 0–499 and the three plant cells at 0–99. Verify on sampled keys that the generator parameters and both seeds equal the span-split check's, and with a synthetic export that the plasmode units' labels, plant and analysis parameters equal the production selection's for the same key and that no record field holds a donor id or a label vector.
- [x] 3.4 Implement the unit (design D5): one production fit, the three predictive integrations from it without refitting, the two forced-rank arms, every arm evaluated serially at the unit's RRPP seed; the five consistency checks with their deviations; a failure recorded without aborting the shard. On the plasmode, one `analyze_trajectories` call per arm with the production selection's parameters plus the arm's setting, the predictive arms at `forced_components` = the unit's `pls_cv` rank (amended 2026-11-12; `notes.md`). Verify on a small design that:
  - the `pls_cv` record equals a plain production evaluation of the same unit;
  - each predictive record equals an evaluation fitted from scratch with that value;
  - one cross-validated fit and two forced fits are made (patched call counts);
  - the record carries every field the spec names;
  - a patched basis is flagged as inconsistent with the check's name and deviation, and the unit is left out of the rates.
- [x] 3.5 Add `scripts/predictive_subspace_array.sbatch` (a profile variable and a dry-run variable that prints the command; absolute paths, BLAS pinned to one thread, a throttled array, `PYTHONPATH` on the worktree's `src`). Verify with `bash -n` and a dry run of its command lines against a temporary directory.

## 4. Merge and reproduction

- [x] 4.1 Implement `merge`: the rate table, the paired table (discordant counts over units available in both arms; exact one-sided loss, one-sided gain and two-sided p-values, 1 when there is no discordant unit), the cell table (rank quantiles, AUROC, measured dimension, capture shares, median observed and `population_latent` statistics, `angle` null 95th percentile) and `accounting.json`. Verify with tests on synthetic records for the columns, failed and inconsistent units counted by cell and excluded from rates, and duplicate handling.
- [x] 4.2 Implement the reproduction checks, written to `reproduction.json`, each listing every differing row with both values:
  - `weak_size` `pls_cv` against `results/span-split-2026-11-05/simulation/production_table.csv`, exact on count and rejections;
  - `plasmode` `pls_cv` anchor against `results/plasmode-2026-11-02/plasmode_summary.csv` (rejections of `delta`, `angle` and `shape` over 500), exact;
  - `plasmode` `pls_cv` plant rows against `results/production-latent-space-2026-11-10/plasmode/rate_table.csv`, exact on count and rejections;
  - the fidelity values on every anchor unit.

  Verify with tests for both outcomes of each check, and that the committed files' columns resolve every compared row.

## 5. Decision

- [x] 5.1 Implement `decide`, which writes `predictive_subspace_decision.json` (rule version 1, design D7), per candidate arm:
  - the gate: `reproduction_failed` on any failed reproduction or fidelity check or any inconsistent unit;
  - S1: three families through `stage_span_paper_grade_decision.judge_rates`, imported unchanged (m = 42, 6, 3);
  - S2: the one-sided exact paired loss p-value per target cell under the count-plus-ceiling rule at m = 57;
  - S3: the δ = 2.0 family through `judge_rates` (m = 36) and the band family by the one-sided exact paired excess p-value under the count-plus-ceiling rule (m = 72), with the count of off-target readings above the bound over both;
  - G: a gain p-value at or below 0.05/18 in at least one δ = 2.0 target cell, with the band's, the weak design's and the plasmode's gains recorded.

  Write the per-arm verdicts in precedence order, the overall verdict with the named arm and the reason it was preferred, every criterion recorded whatever the verdict, an unavailable reading failing its criterion, and the consequence statement for each verdict. Verify with tests on synthetic tables for each per-arm verdict with the failing readings named, precedence, the gate still recording S1–G, both S3 families failing separately, a gain only outside δ = 2.0 failing G, the preference rule and its tie-break, and the allowed counts at m = 57 and 72.
- [x] 5.2 Implement the descriptive tables (design D8): the control arm, the rank-plus-one arm, observed against population statistics per arm, the surgery cell, the weak design's target cells, production's band off-target rates, and rank and AUROC per design. Write the driver's module docstring (profiles, arms, seeds, parents, outputs, decision) and add the driver's row to the validation-drivers table of `AGENTS.md`. Verify with tests that each table has one row per cell and arm it names, and with `uv run ruff check scripts/predictive_subspace_check.py` and the pre-commit gate.

## 6. Runs and report

- [x] 6.1 Run a local smoke of the three profiles at 2 replicates per cell with BLAS pinned (every design of `simulation`, including the planted scale in the band; a few `weak_size` and `plasmode` units at 999 permutations). Record per-unit timing by profile and arm and any failed planted fixed point by design, and confirm that `merge`, the reproduction checks and `decide` run, mismatching only because of the smoke's size. Verify that the timing, the core-hour projection, the proposed cap and any amendment are written into `notes.md` and put to the user for approval before any array is submitted.
- [x] 6.2 Run the three profiles on `ing` through the sbatch array from a clean worktree at a pushed revision, then `merge` each. Verify 0 inconsistent units, every reproduction and the fidelity check passed, and failures only where the smoke found them. Record the job ids, `sacct` output and package versions.
- [x] 6.3 Run `decide`, then commit aggregates only under `results/predictive-subspace-<date>/`: settings, rate, paired and cell tables, accounting, reproduction, the decision, `PROVENANCE.txt` and `sacct.txt`. Write `docs/reports/predictive-subspace-<date>.md` with the verdict, S1–G per arm, every reproduction, the descriptive comparisons and what follows from the verdict. Verify that every number the report quotes can be read from the committed aggregates.
- [x] 6.4 Add the measured-dimension section to `docs/measurement-visibility.md` (spec `latent-space-metric-compatibility`): the fit rank against the measured dimension, the predictive subspace and its three coordinate systems, the control arm, what is given up, the verdict and its readings, and that production stays on the fitted span's scores until a separate change. Update roadmap open item 1 with the verdict and the follow-up it permits, and summarize the run in `AGENTS.md`. Verify with `uv run mkdocs build` and the pre-commit gate.

## Workflow follow-up

- Sync the delta specs and archive the change once the report is reviewed.
- Under `propose_paper_grade`, a separate change runs the named arm at 500 × 999 on the stage-span paper grade's columns and the band. Only its verdict opens the default flip, the CLI flag, attribution and the split test on the measured subspace, and the SEA-AD re-issue.
- Under `keep_production`, the three values stay opt-in diagnostics and the measurement-visibility page records the cost of the measured dimension as a limitation.
