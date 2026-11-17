# Proposal

## Why

Every real-noise reading MOTCO has comes from one export: SEA-AD astrocytes, 80 donors on 10,000 features (roadmap, "Real-noise evidence rests on one export"). That one anchor carries more than it can:

- **It vetoed the predictive subspace.** `predictive_orthonormal` gains in simulation and loses no protected power, but on the SEA-AD anchor it rejects `delta` in 67 of 500 replicates and `shape` in 38, where production rejects 21 and 13 (roadmap item 18).
- **The diagnosis left a prediction that nothing has read.** Item 19 traced the excess to a stage fit that is mostly in-sample (stage R² 0.84 against a held-out AUROC of 0.636), read under RRPP's unrestricted draws. If that is the cause, the excess shrinks as donors per feature grow and is absent on ADNI, at 0.2 features per donor against SEA-AD's 125.
- **Production has no second reading either.** Its size on real noise is known at SEA-AD's configuration only. ADNI's design was measured in simulation (item 9), never on ADNI's blocks.

Roadmap open item 1 reopens `predictive_orthonormal` on size held on an ADNI or ROSMAP anchor together with a gain at that design, "under a rule that scopes admissibility by design and is declared before the reading". Open item 2 names ADNI as the cohort to read first. Its export exists (`adni_sexdiff`, 1,501 participants on 298 metabolites), and its size lets participants be subsampled into many datasets, each with its own fit, where SEA-AD gives one.

A label-free probe of the stage fit on ADNI subsamples (2026-11-16, `notes.md`; it reads no rejection rate) shaped the design:

- **The cohort's fit reproduces.** Rank 8 at a held-out AUROC of 0.599514, the case study's values, under MOTCO 0.25.0.
- **The CV rank hides the mechanism at small n.** At 80 participants the CV selects rank 2 in 10 of 16 subsamples. There the predictive subspace is the scores' own two dimensions, so a ladder read only at the CV rank would hardly test item 19's guess.
- **A fixed rank exposes it.** At rank 19 the median in-sample stage R² is 0.93 / 0.62 / 0.31 / 0.16 at n = 80 / 200 / 500 / 1,501, while the held-out AUROC stays near 0.55 to 0.60. That brackets SEA-AD's 0.84.

## What Changes

- **A driver, and no package change.** `scripts/adni_anchor_check.py` (`preflight` | `run` | `merge` | `decide`) reads the `adni_sexdiff` export in place through an input contract. Every arm is one `analyze_trajectories` call on options that already exist. `motco` stays at 0.25.0.
- **Four design points on ADNI's blocks, the global null only.** Participants are subsampled within sex × stage cells at n = 80, 200 and 500, in the cohort's proportions, and the fourth point is the full cohort of 1,501. Each has 500 replicates at 999 permutations: 2,000 units. A subsampled unit has its own participants and its own fit. The cohort's 500 units share one fit and differ in the labels, as SEA-AD's anchor does.
- **Only permuted labels are analysed.** Each unit permutes the sex labels within stage, the plasmode's null-making step, with blocks and covariates left on their rows. No analysis runs under the real labels.
- **Operating arms, which decide.** Production's scores (`pls_cv`) and `predictive_orthonormal` on the same fit at the unit's CV rank, each under RRPP's unrestricted draws and under the within-stage draws (`__within_stage`). They run at the case-study configuration with the `df` residual correction.
- **Fixed-rank arms, which read the mechanism.** The same four combinations at rank 8 (the cohort's CV rank) and at rank 19 (the highest rank the CV can select). They put the over-fitted stage fit the CV avoids at small n in front of the test, with ADNI's noise.
- **One control arm.** `pls_cv__uncorrected` drops the `df` correction. At the cohort it is the configuration the case study was issued under (MOTCO 0.10.0).
- **A pre-flight that pins the cohort's fit.** The cross-validated fit on the full cohort must select rank 8 at AUROC 0.599514 and equal the forced-rank call exactly, before any unit runs.
- **A size rule per arm, and a size scope.** An arm's *size scope* is the largest set of design points, reaching up to the full cohort with none skipped, on which its `delta` / `angle` / `shape` readings hold the count-plus-ceiling rule (`judge_rates`, imported unchanged) as one family.
- **One predeclared verdict** in `adni_anchor_decision.json` (rule version 1), read on the operating arms under the unrestricted draws:
  - `gate_failed`: the pre-flight failed or a unit is inconsistent. Nothing is read.
  - `production_excess`: production does not hold size on the whole ladder.
  - `excess_persists`: `predictive_orthonormal` fails at the full cohort.
  - `excess_scoped`: it holds at the full cohort and fails lower on the ladder.
  - `excess_absent`: it holds at all four points.

  Each verdict's consequence for the roadmap is written before any unit runs. The within-stage arms' scopes are recorded beside it.
- **A second predeclared reading, of the fixed ranks**: `reference_exceeds` | `absent` | `reproduced` | `reproduced_unrestored`. It says whether ADNI's noise gives the predictive subspace SEA-AD's excess once the fit is as over-fitted, and whether the within-stage draws remove it there too. It never changes the verdict.
- **The scoping rule the roadmap asks for.** A real-noise size reading speaks for the export and the design points it was read at. A hold on one cohort does not lift a failure on another, and size alone admits nothing: reopening an option at a design also needs a gain and specificity there, each under its own predeclared decision. This change reads the size part for ADNI.
- **Descriptive readings beside it:** the rank distribution and rates by rank class per point, the in-sample against the out-of-fold stage fit, the residual dispersion by stage, exact paired comparisons, and SEA-AD's committed anchor counts in one cross-cohort table.
- **Documentation and records.** A dated report `docs/reports/adni-real-noise-anchor-<date>.md`; aggregates under `results/adni-real-noise-anchor-<date>/`; the verdict in roadmap open items 1 and 2, in "Real-noise evidence rests on one export" and in the Done list; a section of `docs/measurement-visibility.md`; one row in the validation-drivers table of `AGENTS.md`.

**Non-goals.**

- **Power and specificity on ADNI's noise.** Plants in a held-out reference span, the out-of-span plant and their off-target readings are the next slice of open item 2.
- **The gain at ADNI's design in simulation,** the other half of open item 1's reopening condition.
- **Reopening `predictive_orthonormal` or changing a default.** No verdict here does either.
- **ROSMAP.** It has no export yet.
- **The other two predictive coordinate systems.** On SEA-AD the three agreed to within three rejections per statistic, and only `predictive_orthonormal` has a path to reopening.
- **The covariate-free lead of item 19.** It is SEA-AD's. On ADNI the p180 run offset would make a covariate-free arm a different question.
- **A feature-subsampling ladder on SEA-AD's export,** the reading that would vary features per donor on SEA-AD's own noise. The `absent` reading of the fixed ranks names it as the lead.
- **Any change in `adni_sexdiff`.** The case study's record is not re-issued.

## Capabilities

### New Capabilities

- `cohort-real-noise-validation`: size readings of the trajectory tests on a cohort's exported blocks under within-stage-permuted labels. It covers the export contract, the design points drawn by subsampling participants within cells, the operating, fixed-rank and control arms, the pre-flight on the cohort's fit, the per-unit checks, the tables, the size scope, the predeclared verdict with its second reading, the rule that scopes a size reading by design, and the dated report. This change adds the ADNI anchor; plants and further cohorts extend it.

### Modified Capabilities

- `latent-space-metric-compatibility`: the measurement-visibility documentation reports the real-noise size readings on a second cohort, what they say about the predictive subspace's excess and about production, and the scope a size reading has.

## Impact

- **Code:** none under `src/motco/`. No version bump and no `CHANGELOG.md` entry.
- **Scripts:** `scripts/adni_anchor_check.py`. It imports unchanged `plasmode_check` (the contract's error classes, file hashing, the within-stage label permutation), `predictive_subspace_check` (arm parameters and records, the fit checks, stage R²), `predictive_size_excess_check` (out-of-fold fitted indicators, explained share, dispersion by stage), `production_latent_space_check` (the exact paired tests), `latent_space_comparison` (settings hashing, units) and `stage_span_paper_grade_decision` (`judge_rates`). It binds by sha256 the within-stage run's committed plasmode `rate_table.csv`, the source of SEA-AD's counts in the cross-cohort table, and refuses a changed file by name.
- **No sbatch array.** `adni_sexdiff` forbids copying its export between machines, and the run is about 7 core-hours. It runs on the workstation, where the export is.
- **Tests:** `tests/test_adni_anchor_check.py` on a synthetic stand-in export in ADNI's layout: every contract violation, the design-point tables, the subsample and label draws, seeds, the arms and their fit counts, the checks, records free of participant data, resume and refusals, the three tables, every scope, every verdict and every fixed-rank reading.
- **Compute:** about 7 core-hours before the smoke, about an hour at seven workers with BLAS pinned. A local smoke replaces the estimate, and the full run waits on the user's approval.
- **Data:** ADNI is controlled-access. The export is read in place and never copied. Records hold seeds and statistics only: no participant identifier, label vector, subsample membership or fold assignment. Only aggregates are committed.
- **Compatibility:** every committed configuration, fixture, record and report is unchanged.
