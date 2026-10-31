# Tasks

## 1. Baseline on the unmodified tree

- [x] 1.1 Confirm the starting state before adding code, and record the package version (0.12.0) the change starts from. Verify that these pass on a clean `main`:
  - the pre-commit gate (`uv run ruff check src/ tests/ && uv run mypy src/motco/ && MOTCO_TEST_PERMS=99 uv run pytest tests/ -m "not slow" --tb=short`);
  - `uv run pytest tests/test_analysis_identity.py tests/test_permutation.py tests/test_case_study_type_i_check.py tests/test_case_study_diagnosis.py -q`.

## 2. The residual correction in `RRPP`

- [x] 2.1 Add `RESIDUAL_CORRECTIONS = ("df", "leverage")` and the public `residual_correction_scale(model_reduced, correction)` to `src/motco/stats/permutation.py` (design D2):
  - the thin SVD and numpy's rank tolerance;
  - κ = √(n / (n − r)) for `df`, and 1 / √(1 − h_ii) from `U[:, :r]` for `leverage`;
  - `ValueError` for an unknown value, for n − r ≤ 0 (stating n and r), and for 1 − h_ii ≤ 1e-10 (naming the row).

  Verify with a new `tests/test_residual_correction.py`:
  - κ = √(80/70) at an 80-row, rank-10 model;
  - a duplicated column uses the rank;
  - leverage multipliers equal 1/√(1 − diag(X(XᵀX)⁻¹Xᵀ)) on an unbalanced design with a numeric covariate, and the h_ii sum to r;
  - a balanced one-factor design gives multipliers equal to κ;
  - each error case raises with the stated content.
- [x] 2.2 Add `residual_correction: str | None = None` to `RRPP` (design D1). Validate it before any draw. Rescale `y_res_np` once by the multipliers, before the serial and parallel paths, and leave the RNG calls and `_RRPPWorker` unchanged. Extend the `RRPP` docstring. Verify in `tests/test_residual_correction.py`:
  - RRPP with the setting absent equals RRPP with it explicitly `None`, and `tests/test_permutation.py` passes unchanged;
  - `df` equals RRPP without a correction on `ŷ_red + κ·e_red`, the diagnosis's construction, cross-checked against `scripts/case_study_diagnosis.rescaled_residuals`;
  - `leverage` equals a manual loop over `default_rng(seed).permutation(n)` of `ŷ_red + (s ⊙ e_red)[perm]` through `estimate_difference`, which also shows the indices are unchanged;
  - centring the leverage-rescaled residuals leaves all three null distributions unchanged to 1e-10;
  - two parallel runs with a correction, a seed and `n_jobs=2` are identical;
  - unknown, saturated and leverage-one inputs raise before any permutation.
- [x] 2.3 Export `residual_correction_scale` in the module contract (spec `module-contracts`), and check that `docs/api/sd.md`'s `RRPP` entry renders the new keyword. Verify:
  - a test imports `residual_correction_scale` from `motco.stats.permutation` and calls it;
  - `uv run mkdocs build --strict` reports no warnings beyond those already present on `main`.

## 3. The analysis, the CLI and the release version

- [x] 3.1 Add `residual_correction: str | None = None` to `TrajectoryAnalysisParams` (design D3):
  - validate it in `__post_init__` against `RESIDUAL_CORRECTIONS` with a `TrajectoryAnalysisError` naming the allowed values;
  - drop it from `to_dict()` while `None`;
  - pass it to `RRPP` in `analyze_trajectories`;
  - when it is set and `permutations > 0`, add `measurement.residual_correction` from `residual_correction_scale(model_reduced, …)`: kind, n and reduced rank, plus κ for `df`, or the leverage and multiplier minimum and maximum for `leverage`.

  Verify in `tests/test_residual_correction.py`:
  - no correction vs `df` on a covariate-adjusted fixture gives identical latent scores, rank, statistics and stage means;
  - the analysis's `leverage` p-values equal `RRPP` with `leverage` on its own latent scores, `_covariate_adjusted_design` models, LS-means, contrast and seed, through the analysis's p-value rule;
  - `"uniform"` raises before `fit_pls_latent_space` is called;
  - a no-covariate analysis with `df` records reduced rank 4;
  - the `df` record states the kind, n, rank and κ;
  - `permutations=0` with `df` keeps the parameter and writes no measurement entry;
  - `tests/test_analysis_identity.py` and `tests/test_analysis.py` pass unchanged.
- [x] 3.2 Add `--residual-correction {df,leverage}` (default absent) to `motco analyze` in `src/motco/cli.py`, pass it to the parameters, and print `RRPP null: residual correction …` with κ or the leverage range when it is applied. Verify with new cases in `tests/test_analysis_cli.py`:
  - a covariate-adjusted run with `--residual-correction df` exits zero, `result.json` holds the correction and κ, and stdout names both;
  - `--residual-correction hc3` exits non-zero with the allowed values in the message;
  - a run without the flag prints no correction line.
- [x] 3.3 Document the setting:
  - a "Residual correction" section in `docs/api/analysis.md` covering the mechanism, both corrections, the record, the CLI example, and the validation status "pending the corrected-permutation check";
  - the regenerated `motco analyze --help` block in `docs/api/cli.md`;
  - the `RRPP` and `analysis.py` paragraphs of `AGENTS.md` (opt-in, byte-identical when absent).

  Verify that `uv run mkdocs build --strict` passes as in 2.3, and that the help block equals `motco analyze --help`.
- [x] 3.4 Bump `src/motco/__init__.py` to 0.13.0 before any paper-grade unit runs (design D10). Verify that `motco --version`, or `python -c "import motco; print(motco.__version__)"`, prints 0.13.0, and that the pre-commit gate passes.

## 4. The validation script

- [x] 4.1 Write the settings and construction of `scripts/corrected_permutation_check.py` (design D5, D7):
  - import `case_study_type_i_check as cs` and `covariate_type_i_check as d10`;
  - `--profile {type_i,power}`, `--table {astrocyte,microglia}` (the two `cohort_30` tables from the design, with their `result.json` sha256 recorded), `--delta`, base seed 2031, `s_z = 10`, five noise covariates;
  - the cells per profile, `cs.check_headroom`, the hashed `settings.json`, and refusal of a resume under other settings.

  Verify with a new `tests/test_corrected_permutation_check.py` (small CV, two replicates, nine permutations):
  - the four Type I runs and the power run have the stated cells and tables;
  - seeds depend only on base seed and replicate, are shared across a replicate's cells, and differ from base seed 2026's at the same replicate;
  - a resume with another table, δ, base seed or permutation count is refused;
  - an over-headroom cell refuses to start.
- [x] 4.2 Write the unit evaluation (design D6). It runs three `analyze_trajectories` calls (`None`, `"df"`, `"leverage"`) with identical parameters and the unit's RRPP seed, and asserts exact equality of the latent scores, selected rank and statistics, raising on any difference. Each arm records the statistics, p-values, null summary, rank, component selection, `cv_mean_auroc`, first-component Spearman ρ with `z`, `measurement.residual_correction` and a latent-score sha256. Add sharded, resumable `run` with the parent's shard helpers. Verify with tests:
  - the arms of a fixture unit share digest, rank and statistics, and differ only in p-values, null summaries and the correction block;
  - a perturbed arm (monkeypatched seed) raises and writes no line;
  - resume skips recorded units;
  - two shards cover every unit exactly once.
- [x] 4.3 Write `merge` (design D11). It produces `corrected_permutation_pvalues.csv`, `corrected_permutation_rates.csv` (R, rejections, rate and rank distribution per arm × cell × statistic) and, for the power profile, `corrected_permutation_contrasts.csv`. The contrasts cover each correction vs none and `leverage` vs `df`, per cell and statistic: both rates, discordant counts and one-sided replicate-clustered sign-flip p-values from `case_study_diagnosis.sign_flip_p_value`. Verify with tests:
  - rates and discordant counts on synthetic records match hand counts;
  - the clustered p-value pools a replicate's cells;
  - the merge is idempotent.
- [x] 4.4 Write the pure `corrected_permutation_decision` and `decide` (design D8–D10). They cover:
  - the pooled tier per table × statistic (α + 2·√(α(1 − α)/R_pooled)) and the per-reading tier at k = 3 over the Type I and power controls;
  - the descriptive uncorrected outcome and `uncorrected_excess_reproduced`;
  - the tie-break by clustered sign-flip over each target pair's four effects;
  - the verdict with its consequence, the calibration criteria and every settings hash;
  - `indicative` below 999 permutations;
  - `size_readings.csv`.

  Verify with tests on synthetic rows for:
  - the bounds 0.0792 and 0.0597 at R = 500 and 2,000;
  - `adopt_df` when only `df` holds;
  - `adopt_leverage` when only `leverage` holds;
  - `adopt_leverage` on a significant gain with no loss;
  - `adopt_df` on no significant difference;
  - `adopt_df` when `leverage` gains on one pair and loses on another;
  - `no_correction_holds`, listing a pooled microglia `delta` of 0.062;
  - the uncorrected flag;
  - an indicative smoke decision.
- [x] 4.5 Write `scripts/corrected_permutation_array.sbatch`, modelled on `scripts/case_study_diagnosis_array.sbatch`: profile, table, δ and shard count as arguments; partition `512x1024`; one CPU and BLAS threads set to 1; n11 excluded. Verify that `bash -n` passes and that a dry render of the `sbatch` command lines for the five runs names the intended output directories.

- [x] 4.6 (Added 2026-10-30 at the user's decision after the smoke.) Report the pooled tier with a replicate-clustered standard error beside D8, descriptively (`clustered_pooled_rows`, `clustered_sensitivity`, tier `pooled_clustered` in `size_readings.csv`). Verify with tests: duplicated readings give a design effect of 4R/(R − 1) and leave the verdict unchanged; disjoint readings give a design effect at most 1.

## 5. Smoke, timing and approval

- [x] 5.1 Run the smoke profiles locally:
  - every run at R = 20 × 99 permutations into `results/corrected-permutation-<date>/smoke-*/`, then `merge` and `decide`;
  - one Type I run at R = 4 × 999 into a scratch directory.

  Verify: every unit completes with zero consistency failures, the decision is marked indicative, and the 999-permutation units run end to end.
- [x] 5.2 Time the units at 999 permutations per profile and table. Append a dated "Smoke and shard plan" section to `design.md` with the smoke rates, the timing, the shards per run and the planned core-hours. Present it to the user and obtain approval for the five arrays. Verify that the section and the recorded approval exist in `design.md`.

## 6. Paper-grade run on `ing`

- [x] 6.1 Push the branch, create a clean worktree on `ing`, confirm that its tree hash equals the local tree, and run `tests/test_residual_correction.py` and `tests/test_corrected_permutation_check.py` there. Verify: identical tree hashes and passing tests on `ing`.
- [x] 6.2 Submit the four Type I arrays (50 shards each) and the power array (shard count from 5.2) at 500 replicates × 999 permutations, into `results/corrected-permutation-<date>/<run>/`. Monitor them to completion, and resubmit any stalled shard. Verify: every task `COMPLETED 0:0`, 1,000 units per Type I run and 7,000 for power, and zero consistency failures.
- [x] 6.3 Run `merge` per run and `decide` on `ing`, and copy back the rate, contrast and size tables and the decision. Commit each CSV under about 1 MB; otherwise name it in provenance. Write `PROVENANCE.txt` per run: command, jobs, tree hash, settings hash, versions (MOTCO 0.13.0), file sha256s and core-hours. Verify that the decision file's hashes equal the five settings files' hashes, and that its verdict is not indicative.

## 7. Report, roadmap and docs

- [x] 7.1 Write `docs/reports/corrected-permutation-<date>.md`. It covers:
  - the construction (fresh family, both tables, the power design);
  - the rates per run × control × arm × statistic against both tiers, and the pooled rates;
  - `uncorrected_excess_reproduced`;
  - the power cost of each correction against none, with `angle` and `shape` stated;
  - the tie-break contrasts;
  - the verdict with its consequence for `sead_sexdiff` and MOTCO;
  - limitations: the `cohort_50` tables, the additive-linear covariates, the structural-analogue blocks, and calibration on the astrocyte table only.

  Verify that every number is read from the committed decision, rate, contrast and size files, and that the report names them.
- [x] 7.2 Update `docs/roadmap.md`: item 8 is marked done with its verdict and report link, item 1's condition on the case study changes under the verdict, and a follow-up is recorded if the verdict names one. Replace the validation status in `docs/api/analysis.md` with the verdict, and add the check (script, verdict, date) to the `analysis.py` paragraph of `AGENTS.md`. Append a "Paper-grade run" record to `design.md`. Verify with `uv run mkdocs build --strict` as in 2.3, and by reading the three documents against the decision file.

## 8. Gate

- [x] 8.1 Run the pre-commit gate (`uv run ruff check src/ tests/ && uv run mypy src/motco/ && MOTCO_TEST_PERMS=99 uv run pytest tests/ -m "not slow" --tb=short`) plus `uv run ruff check scripts/corrected_permutation_check.py`. Confirm that all pass, that `tests/test_analysis_identity.py` and the committed fixtures are unchanged (`git diff --stat main -- tests/data` is empty), and that the version is 0.13.0.

## Workflow follow-up

- Open the pull request, merge, tag `v0.13.0`, sync the `rrpp-residual-correction`, `corrected-permutation-validation`, `real-data-trajectory-analysis` and `module-contracts` specs into `openspec/specs/`, and archive the change.
- Under `adopt_df` or `adopt_leverage`: propose the primary re-issue in `sead_sexdiff` (both cell types, `--residual-correction <adopted>`, pin 0.13.0) in its own change.
- Under `no_correction_holds`: propose a diagnosis of the failing readings before any re-issue.
- Whether the adopted correction should become the analysis default is a separate, later change.
