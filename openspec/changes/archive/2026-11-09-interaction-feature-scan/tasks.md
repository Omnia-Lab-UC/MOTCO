# Tasks

## 1. The per-feature scan (`motco.stats.interaction`)

- [x] 1.1 Factor the split test's input validation and column-space setup into a shared helper in `src/motco/stats/interaction.py` (design D4), leaving its messages and results unchanged. Verify that `tests/test_interaction_split.py` passes unchanged.
- [x] 1.2 Implement `interaction_feature_scan` and `InteractionScanResult` with observed quantities only (design D1, D3).
  - Per feature: SS_int, RSS, F, a degenerate flag.
  - The maximum F and its column; df_int, df_res; the degenerate count and columns.
  - Validation before any computation: rows, finite values, nested models, positive df, a non-negative permutation count, and at least one non-degenerate feature.

  Verify in `tests/test_interaction_scan.py` that:
  - every F matches direct hat-matrix computation (1e-10);
  - each column's F is unchanged under a·y + b;
  - the per-feature sums add to the split test's `total` sums;
  - a constant column is flagged degenerate and leaves the other results equal to a scan without it;
  - an all-degenerate outcome and every validation error raise with their messages.
- [x] 1.3 Add the permutation path (design D2, D3).
  - RRPP's serial draws at the seed, with corrections through `residual_correction_scale`, applied before any draw.
  - Per draw, one W'ᵀ·E product, the per-feature F*, the null maximum and its column, and the step-down suffix-maximum counts.
  - The global p-value, the step-down adjusted p-values with their running maximum, and the B × p null F only with `return_null`.

  Verify in the same file that:
  - at every draw of a small case the null F equal a brute-force refit of Ŷ_R + E[π], with π from `np.random.default_rng(seed)`;
  - the permuted outcomes match RRPP's serial path at the same seed;
  - the adjusted p-values equal a direct step-down implementation on the returned null, are monotone in F, and the top one equals the global p-value;
  - every null maximum's column attains that draw's maximum;
  - the F and every p-value are identical with and without `df`;
  - `leverage` scales rows by 1/√(1 − h_ii) and leaves the observed F unchanged;
  - an unknown correction is rejected before any draw;
  - zero permutations give no p-values;
  - a 10,000-feature, 80-sample case with 12 model columns runs 999 draws without forming an n × p permuted matrix, and its runtime is recorded in the test output.
- [x] 1.4 Export `interaction_feature_scan` and `InteractionScanResult` from `motco.stats`, and document the function:
  - its docstring and the `motco.stats` module list;
  - a section in `docs/api/sd.md`;
  - the `stats/interaction.py` bullet in `AGENTS.md`.

  Verify with a test that the direct and flat imports return the same object, and that `uv run mkdocs build` adds no warnings beyond those already on `main`.

## 2. The opt-in in `motco analyze`

- [x] 2.1 Add `TrajectoryAnalysisParams.interaction_scan` (default `False`, omitted from `to_dict()` while false). Verify with tests for the default and for the setting's absence and presence in the resolved parameters.
- [x] 2.2 Run the scan in `analyze_trajectories` on `X` with the covariate-adjusted models, the seed, the permutations and the residual correction, the correction applied only when permutations run (design D5).
  - Write the `interaction_scan` block after `measurement` and any `span_split` block.
  - Build the per-feature table with the per-stage group differences from `estimate_betas(model_full, X)` and the LS-mean rows, and expose it as `TrajectoryAnalysisResult.interaction_scan`.

  Verify in `tests/test_analysis_interaction_scan.py` that:
  - latent scores, rank, statistics, p-values, null summary, stage means, the `span_split` block and attribution outputs are identical with and without the setting;
  - the block and table equal the primitive applied to the same inputs (six covariates, 99 permutations, `df`);
  - the differences equal the full model's LS-mean differences per stage at the reference row;
  - the scan is identical with and without latent adjustment `stage`;
  - 0 permutations record observed values only;
  - `draws_match_rrpp` is false with two jobs;
  - a constant feature is flagged and ignored;
  - `tests/test_analysis_identity.py` passes unchanged.
- [x] 2.3 Write `interaction_scan.csv` from `write_analysis_outputs` when the scan ran, and nothing new otherwise. Verify with tests on both cases: the file's columns, row count and joint feature order, and the absence of the file without the scan.
- [x] 2.4 Add `--interaction-scan` to `motco analyze` and print one summary line: the maximum F and its feature, the global p-value and the hits at 0.05. Verify with `tests/test_analysis_cli.py` cases for the flag, the `result.json` block, the CSV, the printed line, and a run with `--latent-adjustment stage`.
- [x] 2.5 Document the option, the record block and the CSV:
  - `docs/api/analysis.md`: the option, a record table row and a pointer to the measurement-visibility reading;
  - the `analysis.py` and `cli.py` bullets of `AGENTS.md`.

  Verify that `uv run mkdocs build` adds no new warnings.

## 3. Validation driver: simulation profile

- [x] 3.1 Create `scripts/interaction_scan_check.py` with `run --profile simulation|plasmode`, `merge --profile …` and `decide`, outputs under `--out-dir/<profile>/`.
  - Settings bound to a hash, round-robin shards, resume of complete units, and refusal of a resume or merge under other settings, naming the keys.
  - The simulation settings are rebuilt through `span_split_check.simulation_settings` and must equal `results/span-split-2026-11-05/simulation/settings.json`, whose sha256 is recorded. The parent tables are bound by sha256 too.

  Verify in `tests/test_interaction_scan_check.py`:
  - the unit counts and replicate ranges equal the parent's;
  - resume works;
  - a settings mismatch is refused, naming the keys;
  - a parent file whose hash differs is refused by name.
- [x] 3.2 Implement the simulation unit (design D6). It composes the parent's steps imported from `span_split_check` (dataset, `pls_cv` fit and evaluation, production record, split inputs and split record), then runs the scan on the split inputs at the unit's RRPP seed.
  - It records the scan summary, the per-block hits and null-maximum shares, and the hits on population-null features (design D8).
  - A failure is recorded without aborting the shard.

  Verify on a small design that:
  - the production and split fields equal `span_split_check.evaluate_simulation_unit`'s record for the same unit;
  - the scan fields equal the primitive applied to the unit's split inputs;
  - the population-null count matches a direct computation from `population_trajectories`.
- [x] 3.3 Implement the simulation merge. It writes:
  - the scan table (global rate, SE and count per design and cell, median hits);
  - the block-balance table;
  - the false-hit table;
  - the accounting;
  - the reproduction check against the parent's `split_table.csv` pseudo-F rows and `production_table.csv`, exact in count and rejections under p < α, listing differing rows.

  Verify with tests on synthetic records for the table columns, failed units counted by cell, and both reproduction outcomes.

## 4. Validation driver: plasmode profile

- [x] 4.1 Add the plasmode settings, rebuilt through `span_split_check.plasmode_settings`.
  - They must equal `results/span-split-2026-11-05/plasmode/settings.json`.
  - The export's file hashes must equal the plasmode's committed settings, and the fidelity gate must pass.
  - Cells: the anchor, e = 2.0 and 3.0 at 500 × 999 in the `paper` stream.

  Verify on the plasmode tests' synthetic export the unit enumeration and the refusal of an export whose hashes differ, naming the files.
- [x] 4.2 Implement the plasmode unit: `plasmode_check.build_unit`, then one `pc.analyse` with the parent's parameters plus `interaction_scan=True`, then the consistency and anchor latent checks.
  - The planted counts are computed from the planted names against the scan table: hits, hits inside and outside K, top-20 by F inside K, in total and per block.
  - The record holds seeds, statistics and counts only.

  Verify on the synthetic export that:
  - the record holds no donor identifier, label vector, feature list or donor-level value;
  - the scan fields equal the analysis's `interaction_scan` block;
  - the trajectory p-values and `span_split` block equal `span_split_check.evaluate_plasmode_unit`'s for the same unit.
- [x] 4.3 Implement the plasmode merge. It writes:
  - the scan table per cell;
  - the false-hit rate outside K and the mean recall, hit precision and top-20 precision;
  - the block-balance table;
  - the accounting;
  - the reproduction check against the parent's `split_table.csv` pseudo-F rows and `trajectory_table.csv`, exact per cell and statistic under p ≤ α.

  Verify with tests on synthetic records for both reproduction outcomes and for the recovery metrics on hand-built counts.
- [x] 4.4 Add `scripts/interaction_scan_check_array.sbatch` for both profiles: absolute paths, BLAS pinned to one thread, a throttled array, one job per unit's analysis. Verify with `bash -n` and a dry run of its command lines against a temporary directory.

## 5. Decision

- [x] 5.1 Implement `decide`, which writes `interaction_scan_decision.json` (rule version 1, design D7).
  - The reproduction gate (`reproduction_failed`).
  - S1 (m = 6), S2 (m = 1) and S3 (m = 2) through `stage_span_paper_grade_decision.judge_rates`, imported.
  - G at e = 3.0 against the largest of the six existing rates on the same units, by more than 2 pooled SE.
  - Verdict precedence, the consequence statement (including that Phase 6 is not reopened), and every descriptive table (design D8), each row naming its source run and grade.

  Verify with tests on synthetic tables that exercise:
  - each verdict, with the failing rates and cells named;
  - `size_excess` over a failing S3;
  - `fwer_excess` over a failing G;
  - the absence of criteria under `reproduction_failed`.
- [x] 5.2 Write the driver's module docstring (profiles, units, seeds, outputs, decision) and add the driver to the `AGENTS.md` architecture beside the other validation scripts. Verify that `uv run ruff check scripts/interaction_scan_check.py` passes.
- [x] 5.3 Bump `src/motco/__init__.py` to 0.20.0. Verify with the pre-commit gate: `uv run ruff check src/ tests/ && uv run mypy src/motco/ && MOTCO_TEST_PERMS=99 uv run pytest tests/ -m "not slow" --tb=short`.

## 6. Runs and report

- [x] 6.1 Run a local smoke of both profiles with BLAS pinned: 2 replicates × 19 permutations per cell in each, the plasmode profile on the local export. Record per-unit timing by profile and the scan's share of it, and confirm the reproduction check runs, mismatching only because of the smoke's size. Verify that the core-hour projection is written into `notes.md` before any cluster submission.
- [x] 6.2 Run the simulation profile on `ing` from a clean worktree through the sbatch array, then `merge`. Verify 0 failed units and a passed reproduction in the merge output, and record the job ids and `sacct` output.
- [x] 6.3 Run the plasmode profile on `ing` from the same worktree, then `merge`. Verify that the fidelity gate passed, 0 failed units and a passed reproduction, and record the job ids and `sacct` output.
- [x] 6.4 Run `decide`, then commit aggregates only under `results/interaction-scan-<date>/`: settings, scan, block-balance and false-hit tables, accounting, reproduction checks, the decision, `PROVENANCE.txt` and `sacct.txt`. Write `docs/reports/interaction-scan-<date>.md` with:
  - the verdict and the four criteria;
  - plasmode power and recovery beside attribution's committed precision;
  - block balance;
  - simulation localization against the split parts and production statistics;
  - simulation false hits.

  Verify that every number the report quotes can be read from the committed aggregates.
- [x] 6.5 Add the scan section to `docs/measurement-visibility.md` (spec `latent-space-metric-compatibility`): what a hit and a non-rejection mean with the measured power, why the maximum, the contrast with attribution and the split test, and what the scan does not provide. Record the verdict under roadmap item 16's thread 3 and summarize the run in `AGENTS.md`. Verify with `uv run mkdocs build` and the pre-commit gate.

## Workflow follow-up

- Sync the delta specs and archive the change once the report is reviewed.
- Under `scan_validated`, `sead_sexdiff` applies the scan to both SEA-AD cell types in its own change, pinning MOTCO 0.20.0 and fixing its multiplicity policy across the scan, the split test and the trajectory statistics before any real-label analysis.
- Depending on the measured power at e = 2.0 and the block balance, a moderated statistic or per-block maxima may be proposed as their own changes.
