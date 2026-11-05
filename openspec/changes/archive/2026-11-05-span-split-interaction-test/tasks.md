# Tasks

## 1. The split interaction test (`motco.stats.interaction`)

- [x] 1.1 Implement `interaction_split_test` and its result types in `src/motco/stats/interaction.py` (design D1, D2), with observed quantities only.
  - Per part (`total`, `span`, `complement`): the dimension, SS_int, RSS, df_int, df_res and the pseudo-F.
  - The observed in-span share, absent when the total SS_int is 0.
  - Input validation before any computation: row and column agreement, finite values, orthonormal basis (1e-8, naming the largest deviation), non-empty complement, nested models, positive df_int and df_res.

  Verify in `tests/test_interaction_split.py` that:
  - the definitions match direct hat-matrix computation (1e-10);
  - span + complement = total for both sums of squares;
  - the result is unchanged under a rotation Q·O of the basis;
  - a noise-free interaction inside the span gives share 1 and a zero complement SS;
  - every validation error fires with its message.
- [x] 1.2 Add the permutation path (design D4): Gram matrices formed from E·Q and E − (E·Q)·Qᵀ, RRPP's serial draw sequence at the seed, residual corrections through `residual_correction_scale`, the pseudo-F and SS p-values (1 + #{null ≥ obs}) / (1 + B), and optional null distributions. Verify in the same file that:
  - at every draw the null F and SS equal a brute-force refit of Ŷ_R + E[π], with π from `np.random.default_rng(seed)`, for a small case;
  - the permuted outcomes match RRPP's serial path at the same seed;
  - null span + complement = total at every draw;
  - F p-values are identical with and without `df`, and the SS nulls scale by κ²;
  - `leverage` scales rows by 1/√(1 − h_ii) and leaves the observed values unchanged;
  - an unknown correction is rejected before any draw;
  - zero permutations give no p-values;
  - a 10,000-feature, 80-sample case runs 999 draws without forming any n × p permuted matrix.
- [x] 1.3 Export `interaction_split_test` and its result type from `motco.stats` and document the module (docstring, `motco.stats` module list, a section in `docs/api/sd.md`, the `stats` architecture bullets in `AGENTS.md`). Verify:
  - a test that the direct and flat imports return the same object;
  - `uv run mkdocs build` adds no warnings beyond the eight already on `main`.

## 2. The opt-in in `motco analyze`

- [x] 2.1 Add `TrajectoryAnalysisParams.span_split_test` (default `False`, omitted from `to_dict()` while false) and reject it with `covariate_latent_adjustment` before fitting, naming both settings (design D6). Verify with tests for the default, the resolved-parameter omission and presence, and the rejection.
- [x] 2.2 Run the test in `analyze_trajectories` after the fit and write the `span_split` block after `measurement`: parts, `span_share`, `residual_correction`, `permutations`, `seed`, `draws_match_rrpp`.
  - Inputs: `latent_artifacts`' features in the model's input coordinates, `qr(x_rotations_)`, and the covariate-adjusted full and reduced models, with the analysis's seed, permutations and residual correction.

  Verify in `tests/test_analysis_span_split.py` that:
  - latent scores, rank, statistics, p-values, null summary and stage means are identical with and without the setting;
  - the block equals the primitive applied to the same inputs (six covariates, 99 permutations, `df`);
  - 0 permutations record observed values only;
  - `draws_match_rrpp` is false with two jobs;
  - `tests/test_analysis_identity.py` passes unchanged.
- [x] 2.3 Add `--span-split-test` to `motco analyze`, print one summary line (three pseudo-F p-values and the in-span share), and exit cleanly with `--latent-adjustment`. Verify with `tests/test_analysis_cli.py` cases for the flag, the `result.json` block, the printed line and the clean error.
- [x] 2.4 Document the option and the record block in `docs/api/analysis.md` (option, record table row, reading pointer to the measurement-visibility page) and in the `analysis.py` and `cli.py` bullets of `AGENTS.md`. Verify that `uv run mkdocs build` adds no new warnings.

## 3. Validation driver: simulation profile

- [x] 3.1 Create `scripts/span_split_check.py` with `run --profile simulation|plasmode`, `merge` and `decide`. It holds settings bound to a hash, round-robin shards, resume of complete units, and refusal of a resume or merge under other settings, naming the keys.
  - The simulation settings import `latent_space_comparison`'s designs, cells and seed functions unchanged.
  - The weak δ is read from `results/latent-space-comparison-2026-11-04/calibration/calibration.json`, with its sha256 recorded.
  - Units: the anchor and translation at replicates 0–499, the other twelve cells at 0–99, 999 permutations (design D7).

  Verify in `tests/test_span_split_check.py` the unit counts and replicate ranges, resume, the settings-mismatch refusal and a refused calibration file whose hash differs.
- [x] 3.2 Implement the simulation unit: the comparison's dataset, one `pls_cv` integration and its evaluation (999 permutations, null distributions kept, so the p-values from the first 199 draws are recorded), and the split test on the integration's linear map at the unit's RRPP seed. Record a failure without aborting the shard. Verify on a small design that:
  - the split record equals the primitive applied to the unit's fit;
  - the `pls_cv` record's statistics and capture equal `latent_space_comparison.evaluate_unit`'s `pls_cv` record for the same unit;
  - the first-199 p-values equal a 199-permutation evaluation at the same seed;
  - the record carries every field named in the spec.
- [x] 3.3 Implement the simulation merge. It writes:
  - the split table (rate, SE and count per design, cell, part and statistic: pseudo-F, uncorrected SS and `df` SS), plus the median in-span share;
  - the production rates and the median population capture;
  - the accounting;
  - the reproduction check against the committed `rate_table.csv`, exact on count and rejections for every shared design, cell and statistic, listing differing rows.

  Verify with tests on synthetic records for the table columns, failed units counted by cell, and both reproduction outcomes.

## 4. Validation driver: plasmode profile

- [x] 4.1 Add the plasmode settings (design D8). They import `plasmode_check`'s export loader, settings, unit builder and analysis parameters unchanged. The export's file hashes must equal `results/plasmode-2026-11-02/settings.json`'s, and the fidelity gate must pass. Cells: the anchor, e = 2.0 and e = 3.0 at 500 × 999 in the `paper` stream. Verify, on the plasmode tests' synthetic export, the unit enumeration and the refusal of an export whose hashes differ, naming the files.
- [x] 4.2 Implement the plasmode unit: `plasmode_check.build_unit`, then one `analyze_trajectories` with the plasmode's case-study parameters, attribution off and the split test on, then the anchor latent check. The record holds seeds and statistics only. Verify on the synthetic export that:
  - the record holds no donor identifier, label vector or donor-level value;
  - the split fields equal the analysis's `span_split` block;
  - the trajectory p-values equal the plasmode's own analysis of the same unit with attribution on.
- [x] 4.3 Implement the plasmode merge: the split table per cell and part, the trajectory rates, the accounting and the reproduction check against the committed `plasmode_summary.csv` (exact per cell and statistic). Verify with tests on synthetic records for both reproduction outcomes.
- [x] 4.4 Add `scripts/span_split_check_array.sbatch` for both profiles: absolute paths, BLAS pinned to one thread, a throttled array, one job per unit's analysis. Verify with `bash -n` and a dry run of its command lines against a temporary directory.

## 5. Decision

- [x] 5.1 Implement `decide`, which writes `span_split_decision.json` (rule version 1, design D10).
  - The reproduction gate (`reproduction_failed`).
  - C1 (m = 18) and C2 (m = 3) through the stage-span paper grade's count-plus-ceiling implementation, imported.
  - C3 at both δ = 2.0 designs against `pls_cv`'s `angle`, by more than 2 pooled SE.
  - Verdict precedence, the consequence statement, and every descriptive table, each row naming its source run and grade.

  Verify with tests on synthetic tables that exercise:
  - each verdict, with the failing rates and designs named;
  - the precedence of `size_excess` over a failing C3;
  - the absence of criteria under `reproduction_failed`.
- [x] 5.2 Write the driver's module docstring (profiles, units, seeds, outputs, decision) and add the driver to the `AGENTS.md` architecture beside the other validation scripts. Verify that `uv run ruff check scripts/span_split_check.py` passes.
- [x] 5.3 Bump `src/motco/__init__.py` to 0.18.0. Verify with the pre-commit gate: `uv run ruff check src/ tests/ && uv run mypy src/motco/ && MOTCO_TEST_PERMS=99 uv run pytest tests/ -m "not slow" --tb=short`.

## 6. Runs and report

- [x] 6.1 Run a local smoke of both profiles: the simulation profile at 2 replicates × 19 permutations per cell, and the plasmode profile at 2 replicates × 19 permutations on the local export, with BLAS pinned. Record per-unit timing by profile, and confirm the reproduction check runs, mismatching only because of the smoke's size. Verify that the core-hour projection is written into `notes.md`.
- [x] 6.2 Run the simulation profile on `ing` from a clean worktree through the sbatch array, then `merge`. Verify 0 failed units and a passed reproduction in the merge output, and record the job ids and `sacct` output.
- [x] 6.3 Run the plasmode profile on `ing` from the same worktree, then `merge`. Verify that the fidelity gate passed, 0 failed units and a passed reproduction, and record the job ids and `sacct` output.
- [x] 6.4 Run `decide`, then commit aggregates only under `results/span-split-<date>/`: settings, split and production tables, accounting, reproduction checks, the decision, `PROVENANCE.txt` and `sacct.txt`. Write `docs/reports/span-split-<date>.md` with:
  - the verdict and the three criteria;
  - localization across φ and modes, and the share against capture;
  - the span part against the production statistics, and the complement against the comparison's arms;
  - plasmode power;
  - the SS diagnostic.

  Verify that every number the report quotes can be read from the committed aggregates.
- [x] 6.5 Add the split-test section to `docs/measurement-visibility.md` (spec `latent-space-metric-compatibility`): the parts, the fitted-span caveat with the weak-design capture, the reading of rejections and non-rejections with the plasmode power, the ratio's rationale with the diagnostic, and what the test does not provide. Add roadmap item 16 with the verdict and its follow-ups (threads 2 and 3), and summarize the run in `AGENTS.md`. Verify with `uv run mkdocs build` and the pre-commit gate.

## Workflow follow-up

- Sync the delta specs and archive the change once the report is reviewed.
- Under `split_validated`, `sead_sexdiff` applies the test to both SEA-AD cell types in its own change, pinning MOTCO 0.18.0 and fixing its multiplicity policy (design Open Questions) before any real-label analysis.
- The result sets the bar for thread 2 (a group-aware basis refitted per permutation) and informs thread 3 (a sparse interaction scan), each its own proposal.
