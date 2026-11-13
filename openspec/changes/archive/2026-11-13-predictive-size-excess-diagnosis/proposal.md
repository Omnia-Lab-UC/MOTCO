# Proposal

## Why

The predictive-subspace pilot (roadmap item 18, 2026-11-12) read `keep_production` on one finding: on the SEA-AD plasmode anchor all three predictive coordinate systems reject `delta` in 65–68 of 500 replicates and `shape` in 38–39, where production rejects 21 and 13 on the same units. The same arms hold size at every simulated design, so the generator lacks whatever drives the excess, and the report left the cause unidentified (§4).

That excess is now what stands between the roadmap and two of its open items:
- **It vetoed gains that are real on simulated noise.** `predictive_orthonormal` lost no protected power, gained in 10 of 18 δ = 2.0 target cells and passed both specificity families. The roadmap reopens it only on "an understood and fixed excess" or on size held on another cohort.
- **Open item 4 inherits it.** Shrinkage LDA estimates the same stage-discriminant subspace under another regularizer, so its real-noise size has to be read too.

Roadmap open item 2 names a guess and two checks on the SEA-AD anchor that need no new export. The stage fit is mostly in-sample at this design (stage R² 0.840 against a held-out AUROC of 0.636), and a two-component fit is conservative on the same anchor (7 / 10 / 5 rejections), so the excess belongs to the two predictive directions of the rank-5 fit. The guess is that those directions hold the in-sample fit, and that its residuals are not exchangeable across stages, which RRPP's draws assume. Nothing has tested it.

## What Changes

- **A diagnosis driver, and no package change.** `scripts/predictive_size_excess_check.py` (`run` | `merge` | `decide`; sbatch array `scripts/predictive_size_excess_array.sbatch`) re-reads the pilot's 800 plasmode units (the anchor at the `paper` stream's replicates 0–499, the real-noise span plant at e ∈ {0.25, 0.5, 1.0} at replicates 0–99, 999 permutations, the case-study configuration) with the pilot's labels and RRPP seeds, so every arm is paired with the committed readings. `motco` stays at 0.24.0.
- **Check W: RRPP's draws restricted within stage.** The reduced model's residuals are permuted among the donors of each stage, never across stages, with the same fitted values, `df` correction and statistics. It is read on the three in-sample predictive coordinate systems. The driver carries the draw loop; with one stratum it must equal `RRPP` exactly.
- **Check O: predictive coordinates fitted out of fold.** Each donor's fitted stage indicators come from a stage-PLS model fitted without that donor: four stage-stratified folds at the unit's production rank, read under the unrestricted draws. The fitted indicators are the one predictive coordinate system whose axes mean the same thing in every fold.
- **Four reference arms that must reproduce the pilot.** `pls_cv`, `predictive`, `predictive_orthonormal` and `predictive_fitted` under the unrestricted draws must reproduce the pilot's committed counts in all four cells exactly. On each unit the driver's own inference path must also equal `analyze_trajectories`' statistics and p-values for those four arms.
- **A predeclared decision on the anchor.** `predictive_size_excess_decision.json` (rule version 1) grades each decisive arm:
  - `restored`: its anchor `delta` / `angle` / `shape` hold the count-plus-ceiling rule (`judge_rates`, imported unchanged; m = 3 at R = 500), the rule the pilot's S1 applied to this family;
  - `reduced`: size does not hold, but a one-sided exact paired test against the arm's unrestricted in-sample counterpart shows fewer rejections on `delta` or `shape` at 0.05/2;
  - `unchanged`: neither.

  Check W takes the weakest grade of its three arms, check O the grade of its one arm. Verdict: `both_restore` | `null_restores` | `out_of_fold_restores` | `partly_explained` | `unexplained`; `reproduction_failed` gates everything. Each verdict's consequence for the roadmap is written before any unit runs.
- **Descriptive readings beside it, which never change the verdict.**
  - **The plant cells:** whether an arm that restores size keeps the real-noise `delta` gain (0.51 → 0.79–0.82 at e = 0.5).
  - **Production under the within-stage draws,** and the out-of-fold coordinates under them.
  - **Covariate-free arms:** the four in-sample coordinate systems with no covariates in the trajectory model, for the report's other named candidate (the six-covariate design with its `df` correction).
  - **A dispersion table:** per coordinate system, the reduced model's residual dispersion by stage, and the in-sample against the out-of-fold stage fit.
- **Documentation and records.** A dated report `docs/reports/predictive-size-excess-<date>.md`; aggregates under `results/predictive-size-excess-<date>/`; the verdict in roadmap open items 1 and 2 and the Done list; the measured-dimension section of `docs/measurement-visibility.md`; one row in the validation-drivers table of `AGENTS.md`.

**Non-goals.**
- **A within-stage option in `RRPP`, the harness or `motco analyze`.** A follow-up change proposes it only on `null_restores` or `both_restore`, with its own size and power readings in simulation. This change reads one anchor and validates no test.
- **Reopening `predictive_orthonormal`.** No verdict here changes the pilot's `keep_production`.
- **Out-of-fold coordinates as a measurement.** The arm is a diagnostic: its coordinates mix four fits, and it has no population counterpart.
- **Other cohorts.** The guess predicts an excess that shrinks as donors per feature grow and is absent on ADNI. That reading needs ADNI's export contract and stays with open item 2.
- **Why simulation does not reproduce the excess.** A generator field for it is considered only once the cause is named.
- **A label-permutation test.** Permuting group labels within stage is exact on the plasmode, but it tests the group main effect too and would reject a pure translation.

## Capabilities

### New Capabilities
- `predictive-size-excess-diagnosis`: the diagnosis of the predictive-subspace coordinates' size excess on the SEA-AD plasmode anchor. It covers the units and arms, the within-stage draws and the out-of-fold coordinates, the reproduction of the pilot and the per-unit consistency checks, the tables, the predeclared decision with its descriptive readings, and the dated report.

### Modified Capabilities
- `latent-space-metric-compatibility`: the measurement-visibility documentation reports the diagnosis of the predictive subspace's real-noise size excess and what its verdict means for the opt-in predictive values.

## Impact

- **Code:** none under `src/motco/`. No version bump and no `CHANGELOG.md` entry.
- **Scripts:** `scripts/predictive_size_excess_check.py` and `scripts/predictive_size_excess_array.sbatch`. They import unchanged `predictive_subspace_check` (the plasmode settings, units, arm parameters and records), `production_latent_space_check` (unit construction, the span plant, the exact paired tests), `plasmode_check` (the export contract, label permutation, seeds, the analysis call), `latent_space_comparison` (settings hashing, units) and `stage_span_paper_grade_decision` (`judge_rates`). They bind by sha256 the pilot's plasmode `settings.json` and `rate_table.csv` and, through them, the plasmode's and the production selection's committed files, and refuse a changed file by name.
- **Private helpers:** the driver's inference path reads `motco.analysis._trajectory_inputs` and the p-value helpers of `motco.simulations.evaluation`, as `case_study_diagnosis.py` does. The per-unit identity check is what keeps that path honest.
- **Tests:** `tests/test_predictive_size_excess_check.py`: the within-stage loop (strata respected, one stratum equal to `RRPP`, the `df` scale), the out-of-fold coordinates (every donor held out once, no donor in its own training set, the all-donor limit equal to `predictive_fitted`), path identity, units and seeds equal to the pilot's, parents, resume, merge, every reproduction outcome and every decision path.
- **Compute:** about 13 core-hours before the smoke. The anchor needs no cross-validation (about 1.5 core-hours for 500 units); the 300 plant units run one production CV each (about 11). That is roughly two hours on the laptop at seven workers with BLAS pinned, or one small array on `ing`. A local smoke replaces the estimate, and launch waits on the user's approval.
- **Data:** the SEA-AD astrocyte `cohort_30` export is read in place through the plasmode's contract. Only within-stage-permuted labels are analysed. No donor identifier, label vector or fold assignment is written.
- **Compatibility:** every committed configuration, fixture, record and report is unchanged.
