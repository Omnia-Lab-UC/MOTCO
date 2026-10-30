# residual-df-corrected-permutation

## Why

The case-study configuration diagnosis (2026-10-30, [report](../../../docs/reports/case-study-diagnosis-2026-10-30.md)) returned `residual_df_shrinkage` for the covariate-adjusted `delta`. RRPP permutes the reduced model's residuals (Freedman–Lane), and their mean square is RSS/n rather than the unbiased RSS/(n − p_red). At the case study's model (n = 80, six covariate columns, p_red = 10) the permuted noise is about 6 % narrower in SD than the real noise. The null of a scale-dependent statistic is then too narrow: `delta` reads 0.066 pooled against α = 0.05, and the rate rises with the reduced model's column count (C2 p = 0.015). The diagnosis computed a remedy in its script, Freedman–Lane on `ŷ_red + κ·e_red` with κ = √(n / (n − p_red)), and that corrected test held size in every reading.

The `sead_sexdiff` primary is withheld until MOTCO offers that correction and validates it (roadmap item 8). This change does both. The case study's conclusions do not depend on it, since its primary rejects nothing. What depends on it is the paper's size claim at the configuration the case study actually runs, and how any future rejection is read.

## What Changes

- **An opt-in residual correction on `RRPP`.** New keyword `residual_correction`, with values `None` (default), `"df"` or `"leverage"`:
  - `"df"` multiplies every reduced-model residual by κ = √(n / (n − r)), where r is the reduced model's rank.
  - `"leverage"` divides each residual by √(1 − h_ii), where h is the reduced model's hat matrix.
  - Either way the permutation indices, the fitted values and the full-model estimation are unchanged.
  - A public helper `residual_correction_scale(model_reduced, correction)` returns the per-row multipliers.
  - With `None`, every output is byte-identical to 0.12.0.
- **The same setting on the real-data analysis.**
  - `TrajectoryAnalysisParams.residual_correction` and `motco analyze --residual-correction {df,leverage}`.
  - The observed statistics are unchanged, and only the permutation null moves.
  - The record gains `measurement.residual_correction`, with the kind, n, the reduced rank and κ or the leverage range. The resolved parameters carry the setting only when it is set.
  - Without it, every record and table equals 0.12.0's output, checked against the committed identity fixtures.
- **Validation at the case study's configuration: `scripts/corrected_permutation_check.py`.** Every unit is analysed three times through `analyze_trajectories`, with no correction, `df` and `leverage`. The three analyses share the data, the latent basis and the RRPP seed. There are two profiles:
  - **`type_i`**, on a fresh seed family (base seed 2031) rather than the replicates the diagnosis chose the remedy on:
    - the case-study construction: six numeric covariate columns, floor 2, net displacement, weak stage signal δ ∈ {0.1, 0.05}, `s_z = 10`;
    - two cohort cell tables: astrocyte `cohort_30` `[[11, 10, 28], [9, 10, 12]]` (n = 80) and microglia `cohort_30` `[[11, 10, 24], [9, 11, 11]]` (n = 76);
    - the anchor and translation controls at 500 × 999, i.e. 4 runs and 4,000 units.
  - **`power`**: the Phase 6 astrocyte table at the generator's full stage effect (δ = 2.0) with the same six covariates:
    - cells: anchor, translation, and magnitude, orientation and shape at the latent-adjustment power grid;
    - 500 × 999, i.e. 7,000 units.
- **A predeclared decision.** A pure function produces `corrected_permutation_decision.json` with one of three verdicts:
  - `adopt_df` or `adopt_leverage`: the named correction holds size, read against a pooled bound per cohort table × statistic and a per-reading bound.
  - `no_correction_holds`: neither correction holds.
  - Tie-break: when both hold, `leverage` is adopted only if it shows a significant paired power gain over `df` and no significant loss.
  - The decision also records whether the fresh family reproduces the uncorrected `delta` excess, and the `angle`/`shape` cost of each correction against no correction.
- **Run, report, release.** Local smoke, a timing measurement and an approval gate, then SLURM arrays on `ing` from a clean worktree (about 140 core-hours, re-timed at smoke). Then a dated report, `results/corrected-permutation-<date>/`, and updates to `docs/roadmap.md` (item 8, item 1's condition), `docs/api/analysis.md`, `AGENTS.md` and the version (0.13.0).
- **Out of scope.**
  - Changing any default. The correction stays opt-in whatever the verdict; a default flip would be its own change.
  - The simulation harness (`SimulationEvaluationParams`), `motco de` and `motco plsr`.
  - The `cohort_50` sensitivity tables (n = 63 microglia), whose larger p_red / n is stated as a limitation.
  - The `stage` latent adjustment's `residualized_basis` finding, which the diagnosis left without a MOTCO change.
  - Any change in `sead_sexdiff`; its re-issue is its own change after the release.

## Capabilities

### New Capabilities

- `rrpp-residual-correction`: the residual-corrected permutation in the `RRPP` primitive. It covers the two corrections and their multipliers, the unchanged permutation draws and fitted values, serial/parallel identity, the errors for an undefined correction, and byte identity when absent.
- `corrected-permutation-validation`: the Type I and power check of the corrected permutation at the case study's configuration. It covers the fresh-seed construction on both cohort tables, the three-arm units and their consistency checks, records, merge, the predeclared size and choice rules, and the verdicts with their consequences.

### Modified Capabilities

- `real-data-trajectory-analysis`: the analysis parameters and `motco analyze` accept the residual correction. The record states it only when set, and outputs without it stay identical.
- `module-contracts`: `motco.stats.permutation` also exports `residual_correction_scale`.

## Impact

- **Code:**
  - `src/motco/stats/permutation.py`: `RRPP` keyword, `residual_correction_scale`, `RESIDUAL_CORRECTIONS`.
  - `src/motco/analysis.py`: params, validation, record.
  - `src/motco/cli.py`: `--residual-correction`, and the summary line.
  - `src/motco/__init__.py`: version 0.13.0.
  - New `scripts/corrected_permutation_check.py`, which imports `case_study_type_i_check` and `covariate_type_i_check` for construction, settings and sharding, and `case_study_diagnosis` for the clustered sign-flip test. New `scripts/corrected_permutation_array.sbatch`.
- **Tests:**
  - new `tests/test_residual_correction.py` (primitive and analysis);
  - CLI cases in `tests/test_analysis_cli.py`;
  - new `tests/test_corrected_permutation_check.py`;
  - the existing identity tests (`tests/test_analysis_identity.py`, `tests/test_permutation.py`) must pass unchanged.
- **API:** additive only, with no breaking change. The `RRPP` return value is unchanged.
- **Results and docs:**
  - `results/corrected-permutation-<date>/`: settings, PROVENANCE, rates, contrasts, decision; `records/` stay on `ing`.
  - `docs/reports/corrected-permutation-<date>.md`.
  - `docs/roadmap.md`, `docs/api/analysis.md` and `AGENTS.md`.
- **Downstream:** `sead_sexdiff` reads the verdict. Under `adopt_*` it re-issues both primaries with the adopted correction and pins 0.13.0, in its own change.
- **Compute:** about 4,000 Type I units at about 17 s per analysis and 7,000 power units at about 14 s per analysis on an `ing` EPYC core, with three analyses per unit. That totals roughly 140 core-hours, re-timed at smoke before the approval gate.
