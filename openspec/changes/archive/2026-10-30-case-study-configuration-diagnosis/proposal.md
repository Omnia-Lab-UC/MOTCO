# case-study-configuration-diagnosis

## Why

The case-study configuration Type I check (2026-10-29, [report](../../../docs/reports/case-study-type-i-2026-10-29.md)) returned `diagnose_configuration_excess` under its predeclared rule. Two findings caused it:
- The covariate-adjusted model, which is the case study's primary, reads `delta` at 0.060–0.072 in every control reading (0.066 pooled over 2,000 units, about 3 Monte Carlo SE above α). It exceeds the bound on translation `delta` at δ = 0.1.
- The `stage` latent-adjusted fallback exceeds at the anchor on `angle` (0.084–0.090) and `shape` (0.078–0.092).

The report withholds the `sead_sexdiff` primary pending a diagnosis and names candidates in its §4. This change runs that diagnosis.

The case study's conclusions do not hang on it. Its primary rejects nothing anywhere, and a slightly liberal `delta` test does not weaken a non-rejection. What depends on the diagnosis is the paper's size claim at the configuration it actually runs, and the reading of any future rejection.

The §4 candidates for the `adjusted` excess, in order:
1. The degrees-of-freedom cost of six covariate columns at n = 80.
2. The LS-mean reference row and the RRPP reduced model under covariates.
3. The Monte Carlo bound at R = 500.

Candidate 2's reference row is already excluded. `tests/test_analysis.py::test_statistics_do_not_depend_on_the_reference_row` shows the three statistics do not depend on it, and RRPP evaluates the same LS-mean rows in every draw.

That leaves a specific mechanism worth testing, which is how the df cost would act. RRPP permutes the reduced model's residuals (the Freedman–Lane scheme). Their mean square is RSS/n ≈ σ²(n − p_red)/n, not σ². The reduced model here has p_red = 10 columns at n = 80, so the permuted noise is about 6 % narrower in SD than the true noise. A statistic that scales with the noise, like `delta`, then gets a null that is too narrow, and the test is liberal. A scale-free statistic such as `angle` does not. That is the pattern observed: `delta` 0.066 pooled, `shape` 0.062, `angle` 0.046. A rough normal-tail estimate gives 0.067 at six columns and 0.056 with none.

The `latent_stage` candidates are:
- the CV-selected rank tail of the residualized basis (rank ≥ 4 in 91 / 500 anchor units, maximum 11), including the per-training-set nuisance coefficients at 31 residual df;
- the residualized basis interacting with the rank floor, which the 2026-10-20 check never measured because it ran at floor 1.

## What Changes

- **A diagnosis script on the parent check's exact data.** The new `scripts/case_study_diagnosis.py` imports `case_study_type_i_check` for its construction. It reads its construction from the committed run directories (δ ∈ {0.1, 0.05}, `s_z = 10`, base seed 2026, 500 replicates × 999 permutations, both controls). It refuses a parent whose settings the current code does not rebuild to the recorded hash. Every unit regenerates the same dataset, covariates, CV folds and RRPP seed. The six-column `adjusted` arm must reproduce the committed p-values on every unit, and the merge reports this as its reproduction check.
- **Trajectory-model arms on the fixed covariate-blind basis.** These need only RRPP, about 0.2 s per 999 permutations. The script fits the basis once per unit and runs a 2 × 3 factorial on it, every arm sharing the RRPP seed:
  - **Confounded covariate `z` estimated, or removed exactly.** "Removed" means the oracle: the known injection is subtracted in the fixed basis, which is the analysis as if `z`'s effect were known. Comparing the two factor levels isolates the cost of estimating one covariate collinear with the group × stage design.
  - **Independent noise columns: 0, 2 or 5.** Comparing these isolates the pure column cost. The noise covariates are nested: the same draws, the first k columns.

  This gives 1, 3 and 6 design columns with `z` (the last is the case study's model) and 0, 2 and 5 under the oracle.
- **A residual-df correction arm for each factorial cell.** It is the same Freedman–Lane test with the permuted residuals rescaled by √(n / (n − p_red)), the same permutation indices and the observed statistic unchanged. It tests the shrinkage mechanism directly: the mechanism predicts the uncorrected `delta` rate rises with p_red while the corrected one stays flat near α.
- **Latent-adjustment arms, each with its own basis fit.**
  - `latent_stage` at the floor of 2 under CV re-runs the committed arm, so its p-values must reproduce.
  - `latent_stage` at floor 1 is the 2026-10-20 rule at this configuration.
  - `latent_stage` at forced rank 2 has no CV and therefore neither a rank tail nor per-training-set coefficients.
  - A tabulation of every arm's anchor rates by selected rank (2, 3, ≥ 4).
- **Paired reading and a predeclared decision.** Rates are reported per run × control × arm × statistic against the 0.0695 bound, as in the parent check. Arm contrasts are paired on the same units, using exact one-sided sign-flip tests on discordant rejections. Each replicate's four readings (two runs × two controls share its seeds) count as one cluster. A pure function returns one verdict per finding:
  - The `adjusted` `delta`: `residual_df_shrinkage`, `collinear_covariate` or `unresolved`.
  - The `latent_stage` `angle`/`shape`: `rank_selection` or `residualized_basis`, with the floor-1 reading descriptive.
  Each verdict carries its consequence. `residual_df_shrinkage` names a follow-up MOTCO change: an opt-in df-corrected permutation, validated by its own Type I and power check before the case study may use it. `unresolved` names the deferred R = 1,000 replication as the next step.
- **Run, report and docs.** A local smoke run, a timing measurement and an approval gate, then two 50-shard arrays on `ing` from a clean worktree (about 35–40 core-hours). After that come a dated report, `results/case-study-diagnosis-<date>/`, and updates to `docs/roadmap.md` (item 1's condition and item 7(c)'s follow-up), `docs/api/analysis.md` and `AGENTS.md`.
- **Out of scope.**
  - No change under `src/motco/`. The correction is computed in the script from the public `RRPP` primitive, by running it on `ŷ_red + κ·e_red`; adopting it in the package is the follow-up change.
  - The R = 1,000 replication, deferred to `unresolved`.
  - The microglia cell sizes and the `group_stage` conditioning.
  - Any change in `sead_sexdiff`.

## Capabilities

### New Capabilities

- `case-study-configuration-diagnosis`: diagnoses the size excess found at the case study's configuration. It covers:
  - reproducing the parent check's data and `adjusted` p-values;
  - the trajectory-model factorial on the fixed basis (covariate `z` estimated or oracle-removed × 0/2/5 noise columns) with and without the residual-df correction;
  - the latent-adjustment rank arms (CV at floor 2, floor 1, forced rank 2);
  - paired contrasts and the rank tabulation;
  - the predeclared verdicts and the consequence each carries for the case study and for MOTCO.

### Modified Capabilities

None. `case-study-configuration-validation` keeps its requirements and its decision. The diagnosis reads that decision and adds its own.

## Impact

- **Code:** new `scripts/case_study_diagnosis.py` and `scripts/case_study_diagnosis_array.sbatch` (modelled on `scripts/case_study_type_i_array.sbatch`). The script imports `case_study_type_i_check` (construction, settings, sharding) and calls `analyze_trajectories` and the public `RRPP`/`estimate_difference` primitives. It also calls `motco.analysis._covariate_adjusted_design`, which `tests/test_analysis.py` already imports, and checks on every unit that its inference path reproduces `analyze_trajectories`. No change under `src/motco/`, so there is no release and the case study's pin is untouched.
- **Tests:** new `tests/test_case_study_diagnosis.py` covers:
  - the refusal of an unreproducible parent;
  - the nested noise columns;
  - the oracle removal against a direct construction without the injection;
  - the script's inference path against `analyze_trajectories`;
  - the correction's permutation indices and the observed statistic it leaves unchanged;
  - reference-row invariance of RRPP p-values, extending the existing statistics test;
  - the clustered sign-flip contrasts;
  - each verdict of the pure decision function;
  - resume and merge.
- **Results and docs:** `results/case-study-diagnosis-<date>/` (settings, PROVENANCE, rates, contrasts, rank tables, decision; `records/` stay on `ing`), `docs/reports/case-study-diagnosis-<date>.md`, `docs/roadmap.md`, `docs/api/analysis.md` and `AGENTS.md`.
- **Downstream:** `sead_sexdiff` reads the verdicts. Nothing there changes in this change. A df-corrected re-issue would follow the MOTCO change that adopts the correction.
- **Compute:** about 1,000 units per strength, each with three CV basis fits (`adjusted`, `latent_stage` at floors 2 and 1), one forced-rank fit and about 13 RRPP-only arms. That is roughly 55–70 s per unit on an `ing` EPYC core, about 17 core-hours per strength and 35–40 core-hours for both. The plan is re-timed in the smoke run before submission.
