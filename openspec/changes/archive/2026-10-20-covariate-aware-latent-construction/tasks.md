# Tasks

## 1. Byte-identity baseline (before any code change)

- [x] 1.1 Add a capture helper in `tests/test_analysis_identity.py` (beside the 0.7.0 one) that builds a seeded two-block, three-stage fixture with one numeric and one categorical covariate, runs the 0.8.0 `analyze_trajectories` (99 permutations, fixed seed, attribution off), and writes `tests/data/analysis_no_latent_adjustment_identity.json` (record without `measurement.runtime_seconds` and `software`, plus SHA-256 digests of the three CSVs). Verify by committing the fixture from the unmodified tree and confirming the comparison test passes on it alongside the existing no-covariate test.

## 2. Nuisance residualization in the PLS core

- [x] 2.1 Add `NuisanceDesign`, `fit_nuisance_coefficients` (lstsq on `[1, D, C]` after a full-rank check, returning the covariate block and the residual df) and `residualize` to `src/motco/stats/pls.py`. Verify with unit tests: the stage-collinear scenario (within-stage coefficient recovered, stage effect kept), the covariate-removal scenario (within-design correlation ≈ 0), and a rank-deficient design raising an error naming the columns.
- [x] 2.2 Add keyword-only `scale=True` and `nuisance=None` to `plsda_doubleCV` and `scale` to `_plsda_auroc`; residualize the outer rest/test sets, the inner training/validation sets and the per-repeat full refits with coefficients from the training rows only; collect the minimum training-set residual df; raise naming the repeat and fold on a rank-deficient training design. Verify with `tests/test_pls.py`: defaults reproduce a stored table and model scores on a fixed input and `random_state`; permuting held-out labels leaves every training fit and held-out residualized matrix unchanged (hook the residualizer to record them); a fold whose training set lacks a factor level raises naming the fold.
- [x] 2.3 Document the two arguments in the `plsda_doubleCV` docstring and `docs/api/pls.md`. Verify with `uv run mkdocs build --strict`.
- [x] 2.4 Add keyword-only `strata=None` to `plsda_doubleCV` (fold stratification labels, default `y`) and `NuisanceDesign.cells()`; the residualized path in `fit_pls_latent_space` stratifies both CV levels on the conditioning cells (design D1, added during implementation: on group-sorted rows the stage-only inner folds empty a group × stage cell in every Phase 6 replicate). Verify with `tests/test_pls.py`: same-partition strata reproduce the default table; `group_stage` on group-sorted rows raises without strata and completes with them; and a core test that `group_stage` CV runs on group-sorted rows.

## 3. Conditioning design and the latent core

- [x] 3.1 Add `ConditioningDesign` and `build_conditioning_matrix(frame, group_col, stage_col, kind)` to `src/motco/stats/design.py` (`stage` → drop-first stage indicators; `group_stage` → group, stage and interaction indicators; `_sort_levels` ordering; `name[level]` column names). Verify with unit tests on a three-stage, two-group frame for both kinds, including column names and order.
- [x] 3.2 Add `nuisance=` to `fit_pls_latent_space`: full-sample coefficients and R² per feature, residualized CV with `scale=False`, final `fit_plsda_model(X_res, ..., scale=False)`, `latent.matrix = model.transform(X)`, the `latent_adjustment` entry in `integration_params`, and a `LatentAdjustmentSummary` on the artifacts (residual df, CV min training df, R² summary, design columns, feature scale). Honour `forced_components` with `nuisance`. Verify with unit tests on the core: without `nuisance` the output is unchanged (compare with a stored result); with `nuisance` the scores equal `transform(X)` and differ from `x_scores_` by a term in the column space of `C`; the summary values match a direct computation.
- [x] 3.3 Add a core-level test that projecting `X` and projecting `X_res` onto the same basis give identical `estimate_difference` statistics and identical `RRPP` null draws when `C` is in both models (the equivalence the design relies on). Verify the test passes at both orientation estimators.

## 4. Analysis entry point and record

- [x] 4.1 Add `covariate_latent_adjustment: str | None = None` to `TrajectoryAnalysisParams`, validate it, omit it from `to_dict()` while `None`, and raise in `analyze_trajectories` before any fitting when it is set without covariates. Verify with tests: invalid value rejected; setting without covariates raises before `fit_pls_latent_space` is called (monkeypatch); `to_dict()` has no key while `None` and has it when set.
- [x] 4.2 Build the conditioning design on the recoded `working` frame, assemble the `NuisanceDesign` from it and `dataset.covariate_design`, pass it to `fit_pls_latent_space`, and keep `_covariate_adjusted_design`, `estimate_difference` and `RRPP` unchanged. Verify with tests: the two analyses (with and without the setting) share `model_full`, `model_reduced`, `ls_means` and `contrast` and differ in latent scores; a `stage`-conditioned run with a covariate equal to the stage index raises naming the covariate; the same seed gives identical p-values across two runs.
- [x] 4.3 Emit the `latent_adjustment` record block (D6) after `covariates` and only when the setting is on; keep the `covariates` block as in 0.8.0. Verify with a record-layout test (conditioning, design columns, `residual_df` = n − 1 − k − q, `cv_min_training_df`, `feature_scale`, R² summary keys) and with the byte-identity tests from 1.1 and the existing 0.7.0 fixture on the no-setting path.
- [x] 4.4 Confirm `test_matches_harness_*`, `tests/test_analysis_identity.py`, `tests/test_committed_config_identity.py` and the orientation-default pins pass unchanged. Verify with the fast pytest command.

## 5. CLI

- [x] 5.1 Add `--latent-adjustment {stage,group_stage}` to `motco analyze`, wire it to the params field, reject it without a covariate flag, and print the conditioning in the summary. Verify with `tests/test_analysis_cli.py`: a `--latent-adjustment stage` run writes the block and prints the conditioning; the flag without covariates exits non-zero naming the requirement. Update `docs/api/cli.md`.

## 6. Type I check at the Phase 6 design

- [x] 6.1 Write `scripts/latent_adjustment_type_i_check.py` (D8): import the D10 construction, three arms per replicate with shared covariates, CV `random_state` and RRPP seed, per-record selected rank, first-component Spearman with the numeric covariate and the latent-adjustment block; `--shard`, resume and `merge`. Verify with a fast test that runs 2 replicates × 9 permutations twice and gets identical records, and that each residualized arm's record carries the block while the baseline's does not.
- [x] 6.2 Add the decision step writing `latent_adjustment_type_i_decision.json` and `latent_adjustment_type_i_rates.csv` (per-arm rates, bound, R, verdicts, admissible conditionings, `adopt_stage` | `adopt_group_stage` | `diagnose_latent_excess`, the baseline arm's correlation median and distribution). Verify with a unit test on synthetic tables covering the three outcomes and checking the baseline arm does not affect the decision.
- [x] 6.3 Add `scripts/latent_adjustment_type_i_array.sbatch` modeled on `covariate_type_i_array.sbatch`. Run the smoke profile (R = 20, 99 permutations) locally and verify the calibration criterion (baseline arm median |ρ| ≥ 0.8); if it fails, raise `s_z`, record it in the settings, and re-run the smoke profile. Report the calibration and the cost estimate to the user before any cluster submission.
- [x] 6.4 After the user's approval, run the paper-grade check (R = 500, 999 permutations) on ing, merge the shards, and write `docs/reports/latent-adjustment-type-i-<date>.md` with the rates, bound, calibration evidence, decision and the mediator limitation. Verify that the decision file exists and the report cites it.

## 7. Power comparison against the covariate-blind basis

- [x] 7.1 Write `scripts/latent_adjustment_power_comparison.py` (D9): the `nuisance` and `clean` columns, matched datasets and RRPP seeds across arms, effect levels per mode, anchors and translation controls, selected rank per record; `--shard`, resume and `merge`. Verify with a fast test that the arms of one replicate share the generated blocks and covariates and differ only in the setting, and that the rates table carries the rank distribution and the rank-1 fraction.
- [x] 7.2 Add the decision step writing `latent_adjustment_power_decision.json` (`candidate_qualifies` | `candidate_fails` | `unavailable` per arm, gain/no-loss/control reads and thresholds, with the baseline arm measured only). Verify with a unit test on synthetic tables covering the gain-and-no-loss, loss and control-excess cases.
- [x] 7.3 Add `scripts/latent_adjustment_power_array.sbatch`. Run the smoke profile locally, settle the grid within the D9 cost envelope, and report the grid and cost to the user before any cluster submission.
- [x] 7.4 After the user's approval, run the paper-grade comparison on ing, merge, and write `docs/reports/latent-adjustment-power-<date>.md` with the rates, rank distributions, verdicts and the mediator limitation. Verify that the decision file exists and the report cites it.

## 8. Docs, version and gate

- [x] 8.1 Update `docs/api/analysis.md`: a latent-adjustment section (the setting, the two conditionings and their FWL reading, leak-free CV, the kept scale, the raw-feature projection and its equivalence, the guards, the record block, the equality of residualization and model covariates, the mediator limitation), replacing the "covariate-blind basis" limitation paragraph. Update the `analysis.py` paragraph in `AGENTS.md`. Update `docs/roadmap.md`: Phase 6 case-study note and a Next changes item recording the decision outcomes and the separate rank-1 follow-up. Verify with `uv run mkdocs build --strict`.
- [x] 8.2 Bump `src/motco/__init__.py` to 0.9.0. Verify that the package version reports 0.9.0.
- [x] 8.3 Run the pre-commit gate (`ruff check`, `mypy src/motco/`, fast pytest) and verify all three pass.

## Workflow follow-up

- Archive the change after the Type I report is merged; the power report may land in a follow-up commit against the same specs.
- Open the separate rank-1 `angle` change with the rank distributions from 7.4 as its evidence.
