# Spec Delta

## Purpose

Establishes the operating characteristics of the covariate-residualized latent space against the covariate-blind basis at the Phase 6 SEA-AD design: size under a calibrated covariate-aligned nuisance axis, and power with and without such an axis. It is the gate the case study must pass before it reports latent-adjusted results.

## ADDED Requirements

### Requirement: Null datasets carry a covariate-aligned nuisance axis

The Type I check SHALL generate null semi-synthetic datasets at the Phase 6 small-n design with confounded numeric and categorical covariates and a feature injection, as the covariate-adjustment Type I check does, and SHALL record the construction, its strengths and the seeds.

#### Scenario: Replicates are reproducible

- **WHEN** the check is run twice with the same parent seed and replicate indices
- **THEN** the generated datasets, covariates and analysis records are identical

### Requirement: Nuisance strength is calibrated on the covariate-blind basis

The injection strength SHALL be calibrated on a recorded smoke run so that the covariate-blind basis's first component tracks the numeric covariate: the median replicate's absolute Spearman correlation between that component's scores and the covariate SHALL be at least 0.8 before any paper-grade run.

#### Scenario: Baseline basis is spent on the nuisance axis

- **WHEN** the paper-grade run completes
- **THEN** the decision file reports, for the covariate-blind arm, the per-replicate correlation of the first component with the numeric covariate and its median at or above 0.8

#### Scenario: Calibration below the threshold blocks the paper-grade run

- **WHEN** the smoke run's median correlation is below 0.8
- **THEN** the strengths are raised and recorded, and the smoke run is repeated before any paper-grade run

### Requirement: Three arms per replicate, each with its own basis

Each replicate SHALL be analysed through the real-data analysis three times with the same covariates declared, the same cross-validation settings, permutation count and RRPP seed: with no latent adjustment, with `stage` conditioning and with `group_stage` conditioning. Each arm SHALL fit its own basis and record its selected rank. The check SHALL support resumable, sharded execution and merge shard records into one table of per-arm, per-replicate p-values and selected ranks.

#### Scenario: Arms differ in basis and rank

- **WHEN** a replicate is analysed in the three arms
- **THEN** each arm records its own latent scores and selected rank, and only the residualized arms record a latent-adjustment block

#### Scenario: Interrupted run resumes

- **WHEN** a shard is re-run after an interruption
- **THEN** it skips replicates already recorded and appends only the missing ones

### Requirement: Residualized arms must hold the one-sided anchor bound

For each arm and statistic, the check SHALL compute the rejection rate at α = 0.05 over the R available replicates and compare the residualized arms' rates with the bound `α + k·√(α(1 − α)/R)` at k = 2. A conditioning SHALL be admissible when every statistic's rate is at or below the bound. The covariate-blind arm SHALL be reported and SHALL NOT enter the decision.

#### Scenario: Baseline arm does not enter the decision

- **WHEN** the covariate-blind arm rejects above the bound
- **THEN** its rates are reported and the decision depends on the residualized arms alone

### Requirement: Type I decision names the admissible conditionings

The decision file SHALL state the per-arm rates, the bound, R, the per-statistic verdicts, the admissible conditionings and the decision: `adopt_stage` when `stage` is admissible, `adopt_group_stage` when only `group_stage` is, otherwise `diagnose_latent_excess`.

#### Scenario: Both conditionings hold

- **WHEN** both residualized arms' rates over 500 replicates are at or below 0.0695 for `delta`, `angle` and `shape`
- **THEN** the decision is `adopt_stage` and both conditionings are listed as admissible

#### Scenario: Only group × stage holds

- **WHEN** the `stage` arm's `shape` rate exceeds the bound and the `group_stage` arm holds every statistic
- **THEN** the decision is `adopt_group_stage` and the `stage` arm's `shape` verdict is `exceeds`

### Requirement: Power comparison against the covariate-blind basis

The power comparison SHALL run magnitude, orientation and shape effects at the Phase 6 effect levels, with the same covariates and arms as the Type I check, in two columns: one with the calibrated nuisance injection and one with the covariates present but the injection at zero. It SHALL use matched datasets and permutation draws across arms, record the CV-selected rank per arm and replicate, and write per-column, per-arm rejection rates and rank distributions.

#### Scenario: Matched datasets across arms

- **WHEN** a replicate is analysed in every arm of one column
- **THEN** the arms share the generated blocks, covariates and RRPP seed, and differ only in the latent-adjustment setting

#### Scenario: Selected rank is recorded

- **WHEN** the comparison completes
- **THEN** its rate table carries, per column, arm and effect level, the distribution of the selected rank and the fraction of replicates at rank 1

### Requirement: Power acceptance rule is written before any run

Per residualized arm against the covariate-blind arm, the rule SHALL require: in the nuisance column, a gain at the top effect of the target mode larger than the pooled-SE gain multiplier; in the zero-injection column, no loss larger than the pooled-SE loss multiplier at any effect; and every control statistic within the anchor bound. The rule and its multipliers SHALL be recorded before the run.

#### Scenario: Gain with the nuisance axis and no loss without it

- **WHEN** the `stage` arm's orientation `angle` rate at the top effect exceeds the baseline's by more than the gain threshold in the nuisance column, no rate falls below the baseline's beyond the loss threshold in the zero-injection column, and every control holds
- **THEN** the `stage` arm's verdict is `candidate_qualifies`

### Requirement: Power decision is advisory per arm

The power decision file SHALL state `candidate_qualifies`, `candidate_fails` or `unavailable` per residualized arm, with the rates and thresholds it read; the covariate-blind arm SHALL be measured, not judged.

#### Scenario: Loss without the nuisance axis fails the candidate

- **WHEN** the `group_stage` arm's magnitude `delta` rate in the zero-injection column falls below the baseline's by more than the loss threshold
- **THEN** the `group_stage` arm's verdict is `candidate_fails`

### Requirement: Paper-grade runs are reported with their limitation

Each paper-grade run SHALL use 500 replicates and 999 permutations and be written as a dated report under `docs/reports/`, citing its decision file. The report SHALL state that residualizing a variable that mediates the stage effect removes real signal, and that no simulation distinguishes a confounder from a mediator, so the choice of covariates is the case study's to justify.

#### Scenario: Report cites the decision file

- **WHEN** a paper-grade run's report is written
- **THEN** it names the decision file, the bound, R and the per-arm rates it reads, and carries the mediator limitation
