# Spec Delta

## Purpose

Establishes the size of the real-data trajectory analysis at the SEA-AD case study's configuration — the rank floor of 2, six numeric covariate design columns, a weak stage signal and a nuisance covariate the covariate-blind basis follows — so the case study's primary analysis, or its latent-adjusted fallback, rests on a measurement at the regime it is run in.

## ADDED Requirements

### Requirement: Null datasets at the case study's configuration

The check SHALL generate null semi-synthetic datasets at the Phase 6 small-n design (three stages, sizes `[[11, 10, 28], [9, 10, 12]]`, methylation and expression, joint magnitude) with the per-omic stage effect weakened to a recorded strength. Each replicate SHALL draw one numeric covariate confounded with group, stage and their interaction, inject its effect into the features, and declare five independent, never-injected numeric covariates: six numeric design columns, no factor.

#### Scenario: Six numeric design columns

- **WHEN** a replicate is analysed with the covariates declared
- **THEN** the record's covariate block lists exactly six design columns, one per declared covariate, and no factor level columns

#### Scenario: Replicates are reproducible

- **WHEN** the check is run twice with the same parent seed and unit indices
- **THEN** the generated datasets, covariates and analysis records are identical

#### Scenario: Controls share the baseline draw

- **WHEN** the anchor and the translation control are generated for the same replicate index
- **THEN** they share the generator seed, the covariate values and the RRPP seed, and differ only in the translation surgery

### Requirement: Calibration places the construction in the case study's regime

Before any paper-grade run, a recorded smoke run SHALL meet three criteria at the anchor: the `adjusted` arm's mean held-out CV stage AUROC lies inside a band declared before the run (default 0.60–0.75); the `unadjusted` arm's median first-component |Spearman ρ| with the confounded covariate is at least 0.8; and the `unadjusted` arm rejects at least one statistic above α. Every iteration SHALL be appended to a calibration record with its strengths.

#### Scenario: Criteria met admit the paper-grade run

- **WHEN** a smoke run's adjusted-arm mean AUROC is inside the band, its unadjusted median |ρ| is at least 0.8 and its unadjusted arm rejects above α on some statistic
- **THEN** the calibration record marks the iteration as met and the paper-grade run may use its strengths

#### Scenario: Unmet criterion blocks the paper-grade run

- **WHEN** any of the three criteria fails on a smoke run
- **THEN** the stage-effect or injection strength is changed, the iteration is recorded with its values, and the smoke run is repeated before any paper-grade run

#### Scenario: Realized placement is reported

- **WHEN** the paper-grade run completes
- **THEN** the decision file reports the adjusted arm's mean AUROC with its quantiles, the unadjusted arm's median first-component |ρ|, and the strengths that produced them

### Requirement: Three arms per unit on matched data

Each unit SHALL be analysed three times at the rank floor of 2 with the net-displacement estimator, the same CV settings, permutation count and RRPP seed: six covariates in the trajectory model only (`adjusted`), no covariates (`unadjusted`), and six covariates with the `stage` latent adjustment (`latent_stage`). Each arm SHALL record p-values, rank, CV AUROC, first-component correlation and null summary; a failed arm SHALL be unavailable. Runs SHALL be resumable and sharded.

#### Scenario: Model-only arms share the basis

- **WHEN** a unit is analysed in the `adjusted` and `unadjusted` arms
- **THEN** both record the same latent scores and selected rank, and only the `adjusted` arm records a covariate block

#### Scenario: Latent-adjusted arm fits its own basis

- **WHEN** a unit is analysed in the `latent_stage` arm
- **THEN** it records a latent-adjustment block with `stage` conditioning and six covariate columns, and its selected rank is at least 2

#### Scenario: Failed arm is unavailable, not a rejection

- **WHEN** the analysis raises an error in one arm of a unit
- **THEN** that arm's record carries the error and no p-values, and the unit's other arms are recorded normally

#### Scenario: Interrupted run resumes

- **WHEN** a shard is re-run after an interruption
- **THEN** it skips units already recorded in any shard file and appends only the missing ones

### Requirement: Judged arms must hold the one-sided anchor bound at both controls

For the zero-effect anchor and the `translation` control at e = 1.0, the check SHALL compute each arm's rejection rate at α = 0.05 for `delta`, `angle` and `shape` over the R available replicates and compare the `adjusted` and `latent_stage` rates with `α + k·√(α(1 − α)/R)` at k = 2. An arm SHALL be admissible when every statistic's rate at both controls is at or below its bound. The `unadjusted` arm's rates SHALL be reported and SHALL NOT enter the decision.

#### Scenario: Bound is computed from available replicates

- **WHEN** some units are unavailable in an arm
- **THEN** that arm's R is the count of available replicates for that control and statistic, and its bound is computed from it

#### Scenario: Unadjusted inflation is reported but not judged

- **WHEN** the `unadjusted` arm rejects above the bound on any statistic
- **THEN** its rates appear in the decision file and the decision depends on the judged arms alone

### Requirement: Decision names the admissible configuration

The decision file SHALL state the rates, bounds, R and verdicts (`holds` | `exceeds` | `unavailable`) per arm, control and statistic, the realized calibration values, the construction and the decision: `primary_validated` when the `adjusted` arm is admissible; otherwise `adopt_latent_stage` when the `latent_stage` arm is; otherwise `diagnose_configuration_excess` naming every exceeding arm, control and statistic. Over several stage-signal strengths the overall decision SHALL be the conjunction.

#### Scenario: Primary holds

- **WHEN** the `adjusted` arm's rates over 500 replicates are at or below 0.0695 for every statistic at both controls
- **THEN** the decision is `primary_validated`

#### Scenario: Primary exceeds and the latent arm holds

- **WHEN** the `adjusted` arm's `shape` rate at the anchor exceeds its bound and the `latent_stage` arm holds every statistic at both controls
- **THEN** the decision is `adopt_latent_stage` and the `adjusted` arm's anchor `shape` verdict is `exceeds`

#### Scenario: Both judged arms exceed

- **WHEN** the `adjusted` and `latent_stage` arms each exceed the bound on at least one statistic
- **THEN** the decision is `diagnose_configuration_excess` and the file names each exceeding arm, control and statistic

#### Scenario: Conjunction over strengths

- **WHEN** the paper-grade run comprises two stage-signal strengths and the `adjusted` arm is admissible in one but not the other
- **THEN** the overall decision is not `primary_validated`, and each run's own decision is reported beside it

### Requirement: Paper-grade run is reported with its consequence and limits

The paper-grade run SHALL use 500 replicates and 999 permutations and SHALL be written as a dated report under `docs/reports/` citing its decision file. The report SHALL state the consequence for the case study under its decision and SHALL carry four limitations: the astrocyte primary cohort's cell sizes only, additive-linear covariate effects, one confounded covariate beside five independent ones, and structural-analogue blocks.

#### Scenario: Report names the consequence

- **WHEN** the report is written
- **THEN** it names the decision file, the bound, R and the per-arm rates it reads; states that the primary stands as issued (`primary_validated`), must move to the `stage` latent adjustment (`adopt_latent_stage`) or is withheld pending diagnosis (`diagnose_configuration_excess`); and lists the four limitations
