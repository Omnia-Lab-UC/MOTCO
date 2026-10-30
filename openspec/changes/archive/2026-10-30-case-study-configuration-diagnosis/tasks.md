# Tasks

## 1. Baseline on the unmodified tree

- [x] 1.1 Confirm the starting state before adding code, and record the package version (0.12.0) the diagnosis runs under. Verify:
  - `MOTCO_TEST_PERMS=99 uv run pytest tests/test_case_study_type_i_check.py tests/test_covariate_type_i_check.py -q` passes;
  - `case_study_type_i_check.settings_hash` of each committed `results/case-study-type-i-2026-10-29/delta*/settings.json` equals its recorded hash, and `build_settings` rebuilds it from the recorded arguments.

## 2. Parent binding, construction and the oracle

- [x] 2.1 Write the parent binding in `scripts/case_study_diagnosis.py` (design D1). Import `case_study_type_i_check as cs`. `--parent <run dir>` loads `settings.json`, checks the hash and refuses on mismatch, and rebuilds through `cs.build_settings` from the recorded arguments, refusing a different hash. The diagnosis `settings.json` records the parent path and hash and refuses a resume under another parent. `--reps` takes a prefix of the parent's replicates. `--perms` below the parent's marks the run as smoke. Verify with a new `tests/test_case_study_diagnosis.py` (small CV, two replicates, nine permutations, a temporary parent written by the parent check's own `run`); it asserts:
  - a tampered parent settings file is refused;
  - a resume under another parent is refused;
  - units regenerate the parent's datasets, covariates and RRPP seeds.
- [x] 2.2 Write the oracle (design D3). Build `X0` with D10's builder at `s_z = 0` and assert `X − X0` is rank one in `z` per block. Reconstruct the fitted covariate-blind projection from `fit_block_scalers`/`transform_blocks`, the PLS centring and `ddof=1` scaling, and `pls_loadings`, and assert it reproduces `latent_scores` to 1e-9 (a hard error otherwise). Form `Y0` by projecting `X0`. Verify with tests:
  - the reconstruction matches `latent_scores` on a fixture unit;
  - `Y − Y0` equals `outer(z, v)` for a fixed vector `v`;
  - a deliberately wrong scaling makes the unit raise.

  If the reconstruction cannot be made to match, implement D3's fallback (a direct `fit_pls_latent_space` call) instead, record that in `design.md`, and test it the same way.

## 3. Trajectory-model arms and the correction

- [x] 3.1 Write the script's inference path and the factorial (design D2, D4). One `analyze_trajectories` call per unit with the six covariates (`z_noise5`) supplies the record, `latent_scores` and `pls_loadings`. The design is built as the analysis builds it, covariate subsets come from `TrajectoryDataset(...).covariate_design`, and `_covariate_adjusted_design`, `estimate_difference` and serial `RRPP` run with the unit's RRPP seed. Run `z_noise{0,2,5}` on `Y` and `oracle_noise{0,2,5}` on `Y0`, with noise columns as prefixes of the parent's draw. Assert on every unit that the script's `z_noise5` equals the record's statistics and p-values. Record `p_red` per arm. Verify with tests:
  - the identity assertion passes on fixture units and raises when the design is perturbed;
  - `z_noise2`'s noise columns equal the first two of `z_noise5`'s;
  - all six arms share the latent scores and the permutation indices;
  - `p_red` is 4, 6, 9, 5, 7 and 10.
- [x] 3.2 Write the corrected arms (design D5): `RRPP` on `ŷ_red + κ·e_red` with κ = √(n / (n − p_red)), the arm's models and the unit's seed, observed statistics from the uncorrected `Y`, and p-values by the analysis's rule. Record `kappa`. Verify with tests:
  - observed statistics are equal between each corrected arm and its uncorrected arm;
  - at κ = 1 the corrected p-values equal the uncorrected ones exactly (same draws);
  - `kappa` is larger for `z_noise5` than for `oracle_noise0`.
- [x] 3.3 Extend the reference-row exclusion to the p-values (spec, reference-row requirement): add a test beside `tests/test_analysis.py::test_statistics_do_not_depend_on_the_reference_row` showing that `RRPP` p-values at one seed are identical with the pooled-mean and a shifted reference row. Verify the test passes.

## 4. Latent arms, records, merge and the decision

- [x] 4.1 Write the latent arms (design D6): `latent_stage_cv` (floor 2), `latent_stage_floor1` (`min_components=1`) and `latent_stage_rank2` (`forced_components=2`), all with the six covariates, the `stage` adjustment and the unit's seed. A `TrajectoryAnalysisError` makes that arm unavailable while the others proceed. Verify with tests:
  - `latent_stage_rank2` records a forced selection at rank 2 and no CV AUROC;
  - `latent_stage_floor1` can record rank 1;
  - a monkeypatched failure in one arm leaves the other arms recorded.
- [x] 4.2 Write records, sharding and resume (design D9). One JSONL line per unit carries the arms, the reconstruction error and the rank-one residual. `run` takes every N-th unit, skips units recorded in any shard file and drops torn lines. Verify with tests:
  - a re-run after a torn line re-evaluates only the missing unit;
  - a different shard count skips recorded units.
- [x] 4.3 Write `merge` (design D9, spec rates and reproduction requirements). It writes `diagnosis_pvalues.csv`, `diagnosis_rates.csv` (R, rejections, rate, the parent's one-sided bound and the verdict per run × control × arm × statistic), `diagnosis_rank_table.csv` (anchor buckets 2, 3 and ≥ 4 for `z_noise5` and the latent arms) and `diagnosis_reproduction.json`. The reproduction file compares `z_noise5` and `latent_stage_cv` with the parent's `case_study_type_i_pvalues.csv` as k/1000 integers and is `not_compared` for smoke runs. Verify with tests on synthetic records:
  - 500 available replicates give the bound 0.0695 and the 34 / 35 verdict boundary;
  - every anchor unit falls in exactly one rank bucket;
  - one altered parent p-value is reported as a mismatch.
- [x] 4.4 Write the contrasts and the pure decision (design D7, D8). The replicate-clustered sign-flip test computes exact one-sided p-values by convolution over the per-replicate `d_r`, for C1–C4 on `delta` (and descriptively on `angle`/`shape`), with discordant counts and per-reading tables. `decide` reads both merged runs and writes `diagnosis_contrasts.csv` and `case_study_diagnosis_decision.json`, giving each verdict, its consequence and the reproduction gate. Verify with synthetic-row tests:
  - the sign-flip p-value equals a brute-force enumeration on a small case;
  - one replicate's four readings contribute one `d_r`;
  - each of `residual_df_shrinkage`, `collinear_covariate`, `unresolved`, `rank_selection`, `residualized_basis` and `reproduction_failed` is reached;
  - smoke runs yield verdicts marked indicative.
- [x] 4.5 Add `scripts/case_study_diagnosis_array.sbatch`, modelled on `scripts/case_study_type_i_array.sbatch` (out dir, parent, shard count, single-threaded BLAS, resumable). Write the script's module docstring (purpose, parent binding, arms, correction, contrasts, verdicts, outputs, examples). Verify:
  - `bash -n` passes;
  - one shard runs locally through the template's command line (`--shard 0/2 --limit 1`);
  - `uv run ruff check scripts/case_study_diagnosis.py tests/test_case_study_diagnosis.py` and `MOTCO_TEST_PERMS=99 uv run pytest tests/test_case_study_diagnosis.py -q` pass.

## 5. Smoke, timing and approval

- [x] 5.1 Run the smoke profiles locally:
  - R = 20 at 99 permutations against both parents, into `results/case-study-diagnosis-<date>/smoke-delta*/`, then `merge` and `decide`;
  - R = 4 at 999 permutations against both parents, into a scratch directory, to exercise reproduction.

  Verify: 40 / 40 units per parent with 0 reconstruction or identity failures, `not_compared` and indicative verdicts at 99 permutations, and zero mismatches over the 8 × 2 units at 999 permutations.
- [x] 5.2 Time units at 999 permutations (`--limit`) and write the shard plan as a dated "Smoke and shard plan" section at the end of `design.md`: smoke rates, timing, shards per parent and planned core-hours. Present it to the user and obtain approval for the two arrays. Verify: the section and the recorded approval exist in `design.md`.

## 6. Paper-grade run on `ing`

- [x] 6.1 Push the branch, create the clean worktree on `ing`, confirm the tree hash matches the local tree, and run `tests/test_case_study_diagnosis.py` there. Verify: identical tree hashes and passing tests on `ing`.
- [x] 6.2 Submit one 50-shard array per parent (500 replicates, 999 permutations) into `results/case-study-diagnosis-<date>/delta*/`, monitoring to completion and excluding or resubmitting shards stalled on n11. Verify: every task `COMPLETED 0:0`, 1,000 units per parent, and arm failures counted.
- [x] 6.3 Run `merge` per parent and `decide` on `ing`, copy the tables and the decision back (committing each CSV under about 1 MB, else naming it in provenance), and write `PROVENANCE.txt` per run: command, jobs, tree hash, parent and diagnosis hashes, versions, file sha256s and core-hours. Verify: `diagnosis_reproduction.json` shows zero mismatches over 3,000 + 3,000 p-values per parent, and the decision file's hashes equal the settings files'.

## 7. Report, roadmap and docs

- [x] 7.1 Write `docs/reports/case-study-diagnosis-<date>.md`. It covers:
  - the construction (the parent binding) and the reproduction counts;
  - rates per run × control × arm × statistic against the bound;
  - the contrasts with discordant counts and sign-flip p-values;
  - the rank tables;
  - both verdicts with their consequences for `sead_sexdiff` and MOTCO;
  - the parent's four limitations, plus the rank-tail / per-training-set confound of `latent_stage_rank2`.

  Verify: every number is read from the committed decision, rate, contrast and rank files, and the report names them.
- [x] 7.2 Update `docs/roadmap.md`: item 1's condition on the case study changes per the verdicts, item 7(c) gains the diagnosis with its report link, and a follow-up item is added when the verdict names one (the corrected permutation, or the R = 1,000 replication). Append a "Paper-grade run" record to `design.md`. Verify: `uv run mkdocs build --strict` reports no warnings beyond those already present on `main`.
- [x] 7.3 Add the diagnosis to `docs/api/analysis.md` beside the case-study check's validation paragraph, and to the `analysis.py` paragraph of `AGENTS.md` (script, verdicts, date). Verify by reading both against the decision file and with `uv run mkdocs build --strict` as in 7.2.

## 8. Gate

- [x] 8.1 Run the pre-commit gate (`uv run ruff check src/ tests/ && uv run mypy src/motco/ && MOTCO_TEST_PERMS=99 uv run pytest tests/ -m "not slow" --tb=short`) plus `uv run ruff check scripts/case_study_diagnosis.py`, and confirm all pass. Confirm `src/motco/__init__.py` is unchanged at 0.12.0 (design D10).

## Workflow follow-up

- Open the pull request, merge, sync the `case-study-configuration-diagnosis` spec into `openspec/specs/` and archive the change.
- Per the verdicts: propose the corrected-permutation change in MOTCO (`residual_df_shrinkage`), the R = 1,000 replication (`unresolved`) or the confounding-strength diagnostic (`collinear_covariate`). Propose the `stage`-at-rank-2 re-issue in `sead_sexdiff` only under `rank_selection` and with the user's agreement.
