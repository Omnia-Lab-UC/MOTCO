# Tasks

## 1. Driver: settings, units and the reference arms

- [ ] 1.1 Create `scripts/predictive_size_excess_check.py` with `run`, `merge` and `decide` (design D10): `settings.json` bound by hash; round-robin shards; resume of complete unit keys; refusal of a resume or merge under other settings, naming the keys; the sha256 of the pilot's plasmode `settings.json` and `rate_table.csv` recorded beside the parents `predictive_subspace_check.plasmode_settings` already binds, with a refusal naming any changed file. Verify in `tests/test_predictive_size_excess_check.py` the resume, the settings refusal and the parent refusal, on the plasmode tests' synthetic stand-in export.
- [ ] 1.2 Build the units from `predictive_subspace_check.plasmode_settings` / `plasmode_units` and `plasmode_check.unit_seeds`, imported unchanged (design D2). Verify that the cells, replicates, label seeds, plant seeds and RRPP seeds equal the pilot's plasmode profile for every unit, that the count is 800 at the committed replicates, and that a changed export file is refused before any unit.
- [ ] 1.3 Implement the reference arms of the unit (design D3): on the anchor, `pls_cv` and the three predictive arms through `analyze_trajectories` at `forced_components` equal to the settings' fidelity rank; on a plant cell, `pls_cv` with the cross-validation and the predictive arms forced at its rank; the pilot's `same_fit` check on every unit; a unit that raises recorded as failed without aborting the shard. Verify that:
  - an anchor unit makes no cross-validation and a plant unit exactly one (patched call counts);
  - on the stand-in export the four reference records equal the pilot's `evaluate_plasmode_unit` records for the same unit (statistics and p-values);
  - a record holds seeds and statistics only.

## 2. The inference path and the two checks

- [ ] 2.1 Implement the driver's inference path (design D4): the covariate-adjusted models, LS-mean rows and contrast from `motco.analysis._trajectory_inputs`; observed statistics from `estimate_difference` under net displacement; the unrestricted null from serial `RRPP` with `residual_correction="df"`; p-values by `_empirical_p_value`. Add the per-unit `path_identity` check on the four reference arms (statistics to 1e-10, p-values exact) and record a failing unit as inconsistent with the arm and deviation. Verify identity on stand-in units of both cell kinds, and that a perturbed latent matrix is flagged and the unit leaves every rate.
- [ ] 2.2 Implement the within-stage draw loop (design D5): per draw and per stage in declared order, a permutation of that stage's rows; the reduced fit, the `df` scale and the per-draw statistics as in `RRPP`; one restricted sequence per unit shared by its within-stage arms. Verify that:
  - every drawn index stays in its row's stage, over many draws;
  - with one stratum the null distributions equal `RRPP`'s at the same seed exactly, with and without the `df` correction;
  - the permuted residuals sum to zero within each stage;
  - two arms of one unit receive the same index sequence.
- [ ] 2.3 Implement the out-of-fold coordinates (design D6): four stage-stratified shuffled folds from a seed derived by sha256 from the driver's name, base seed 2037, the replicate and `folds`; `fit_plsda_model` on the training rows of the analysis's PLS input features at the unit's rank and scaling; `transform(X_held) · y_loadingsᵀ` for the held-out rows. Verify that:
  - every donor is held out once and is in no training set that produced its row, and every training set holds every stage;
  - with every donor in the training set the result equals the `predictive_fitted` matrix to 1e-8;
  - the same replicate gives the same folds and matrix, and two replicates give different folds;
  - the record holds the fold seed and no assignment.
- [ ] 2.4 Assemble the fourteen arms of the unit (design D7): the three predictive within-stage arms and `pls_cv__within_stage` on the reference arms' latent matrices; `fitted_out_of_fold` and `fitted_out_of_fold__within_stage`; the four covariate-free arms as `analyze_trajectories` calls on a dataset without covariates, with `df` kept. Record per arm the observed statistics, p-values and `angle` null 95th percentile, and per coordinate system the residual dispersion by stage relative to the pooled one, the in-sample stage R² and the out-of-fold explained share (design D9). An arm that raises is recorded as failed on its unit. Verify the record's fields, that a patched failing arm leaves the other thirteen, and that the covariate-free arms' latent matrices equal the reference arms'.
- [ ] 2.5 Add `scripts/predictive_size_excess_array.sbatch` (a dry-run variable that prints the command; absolute paths, BLAS pinned to one thread, a throttled array, `PYTHONPATH` on the worktree's `src`). Verify with `bash -n` and a dry run of its command lines against a temporary directory.

## 3. Merge and reproduction

- [ ] 3.1 Implement `merge`: the rate table; the paired table with each arm's counterpart as the spec defines it and against `pls_cv` (discordant counts over units available in both arms; exact one-sided fewer, one-sided more and two-sided p-values, 1 when there is no discordant unit, through `production_latent_space_check`'s exact tests imported unchanged); the cell table (rank, the dispersion medians by stage and their ratio, both stage-fit readings, the median `angle` null 95th percentile); `accounting.json`. Verify with tests on synthetic records for each table, the counterpart mapping of all ten non-reference arms, and the no-discordant case.
- [ ] 3.2 Implement the reproduction check, written to `reproduction.json`: the four reference arms' available and rejected counts per cell and statistic against `results/predictive-subspace-2026-11-12/plasmode/rate_table.csv`, exact, listing every differing row with both values. Verify with tests for both outcomes and that the committed file's columns resolve all 48 compared rows.

## 4. Decision

- [ ] 4.1 Implement `decide`, which writes `predictive_size_excess_decision.json` (rule version 1, design D8):
  - the gate `reproduction_failed` on a failed reproduction, any inconsistent unit or a missing pre-flight record, with the grades still recorded;
  - per decisive arm, `restored` through `stage_span_paper_grade_decision.judge_rates` (imported unchanged) on its three anchor readings, else `reduced` when the one-sided exact paired p-value for fewer rejections than its counterpart is at most 0.025 on `delta` or `shape`, else `unchanged`; an unavailable reading graded `unchanged` and named;
  - the within-stage check as the weakest of its three arms, the out-of-fold check as its arm;
  - the verdict and its consequence statement.

  Verify with tests on synthetic tables for each of the six verdicts, the weakest-arm rule, both routes into `partly_explained`, the gate with grades recorded, and that the consequence statements are those of design D8.
- [ ] 4.2 Implement the descriptive tables (design D9) inside the decision: plant-cell rates and paired tests per arm, with the statement for each `restored` arm against production's `delta` at e = 0.5; `pls_cv__within_stage` and `fitted_out_of_fold__within_stage` on the anchor; the covariate-free arms; the dispersion table; null widths. Verify that changing any descriptive reading leaves the verdict unchanged. Write the driver's module docstring (units, arms, seeds, parents, checks, outputs, decision) and add its row to the validation-drivers table of `AGENTS.md`. Verify with the pre-commit gate: `uv run ruff check src/ tests/ && uv run mypy src/motco/ && MOTCO_TEST_PERMS=99 uv run pytest tests/ -m "not slow" --tb=short`.

## 5. Runs and report

- [ ] 5.1 Run the pre-flight and a local smoke on the real export with BLAS pinned (design D3): the cross-validated `pls_cv` on anchor replicates 0 and 1 against the forced call; then two anchor units and two units per plant cell at 999 permutations through `run`, `merge` and `decide`. Record in `notes.md` the pre-flight outcome, per-arm and per-unit timing, and the projected core-hours for the laptop and for `ing`. Stop and give the projection to the user; the full run waits on their approval and their choice of where it runs.
- [ ] 5.2 Run the 800 units where approved, from a clean worktree at a pushed revision, then `merge`. Verify 0 failed and 0 inconsistent units and a passed reproduction. Record the package versions and, on `ing`, the job ids and `sacct` output.
- [ ] 5.3 Run `decide`, then commit aggregates only under `results/predictive-size-excess-<date>/`: `settings.json`, the rate, paired and cell tables, `accounting.json`, `reproduction.json`, the decision and `PROVENANCE.txt`. Write `docs/reports/predictive-size-excess-<date>.md` with the verdict, each decisive arm's grade and counts beside the pilot's, the reproduction, every descriptive reading, what follows from the verdict and the limitations. Verify that every number the report quotes is in a committed aggregate and that no committed file holds a donor identifier, a label vector or a fold assignment.
- [ ] 5.4 Add the diagnosis to the measured-dimension section of `docs/measurement-visibility.md` (spec `latent-space-metric-compatibility`): what each check tests, the verdict, the decisive arms' anchor rates beside the pilot's, and that the predictive values stay opt-in diagnostics. Update roadmap open items 1 and 2 and the "Real-noise evidence rests on one export" entry with the verdict and its consequence, and add the Done-list entry. Verify that `uv run mkdocs build` adds no warnings beyond those on `main` and that the roadmap's statements match the decision file.

## Workflow follow-up

- Sync the delta specs and archive the change once the report is reviewed.
- Under `null_restores` or `both_restore`, a separate change proposes within-stage draws as an opt-in option of `RRPP`, the harness and `motco analyze`, byte-identical when unset, and reads their size and power on the pilot's simulated designs, the weak design and this anchor.
- Under `out_of_fold_restores`, `partly_explained` or `unexplained`, no package change follows; the roadmap records the lead for open item 2's cohort readings.
