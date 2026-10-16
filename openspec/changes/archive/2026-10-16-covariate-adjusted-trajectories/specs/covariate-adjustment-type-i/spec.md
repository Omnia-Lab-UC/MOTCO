# Spec Delta

## Purpose

Establishes that covariate-adjusted real-data analysis holds its nominal Type I error at the Phase 6 SEA-AD design when covariates are confounded with group, stage and their interaction. It is the gate the case study must pass before it reports covariate-adjusted results.

## ADDED Requirements

### Requirement: Confounded null datasets at the Phase 6 design

The check SHALL generate null semi-synthetic datasets (group effect size 0) at the Phase 6 small-n design: three stages, group × stage sizes `[[11, 10, 28], [9, 10, 12]]`, the methylation and expression blocks, and joint magnitude. Each replicate SHALL draw one numeric and one three-level categorical covariate per sample. Their distributions SHALL depend on group, stage and the group × stage interaction, so that the covariate means differ between the groups' trajectories. Each replicate SHALL add a covariate effect to the standardized-unit features of the measured blocks, along fixed random feature directions with a declared strength. Each replicate SHALL be reproducible from a declared parent seed, and the covariate construction, its strengths and the seeds SHALL be recorded.

#### Scenario: Covariates are confounded by construction

- **WHEN** the check generates its replicates
- **THEN** the recorded per-cell covariate means differ across group × stage cells, with a non-zero group × stage interaction contrast for the numeric covariate

#### Scenario: Replicates are reproducible

- **WHEN** the check is run twice with the same parent seed and replicate indices
- **THEN** the generated datasets, covariates and analysis records are identical

### Requirement: Paired adjusted and unadjusted arms share one latent fit

Each replicate SHALL be analysed through the real-data analysis twice: once with both covariates declared (the adjusted arm) and once without (the unadjusted arm). Both arms SHALL use the default net-displacement estimator, cross-validated rank selection with the harness settings, the same RRPP permutation count and the same RRPP seed. The two arms SHALL share the latent space, so they differ only in the trajectory model. The check SHALL support resumable, sharded execution and SHALL merge shard records into one table of per-replicate p-values for `delta`, `angle` and `shape` in each arm.

#### Scenario: Arms differ only in the model

- **WHEN** a replicate is analysed in both arms
- **THEN** both arms record the same latent scores and selected rank, and only the adjusted arm records a covariate block

#### Scenario: Interrupted run resumes

- **WHEN** a shard is re-run after an interruption
- **THEN** it skips replicates already recorded and appends only the missing ones

### Requirement: Adjusted arm must hold the one-sided anchor bound

For each statistic, the check SHALL compute the rejection rate at α = 0.05 in each arm over the R available replicates. It SHALL compare the adjusted arm's rate with the bound `α + k·√(α(1 − α)/R)` at k = 2. The adjusted arm SHALL hold when every statistic's rate is at or below the bound. The check SHALL write `covariate_type_i_decision.json` stating the per-statistic rates in both arms, the bound, R, the per-statistic verdicts (`holds` | `exceeds`), and the overall decision: `covariates_validated` when the adjusted arm holds for every statistic, otherwise `diagnose_covariate_excess`. The unadjusted arm SHALL be reported descriptively, as evidence that the construction confounds, and SHALL NOT enter the decision. The paper-grade run SHALL use 500 replicates and 999 permutations, and its findings SHALL be written as a dated report under `docs/reports/`.

#### Scenario: Adjusted arm within bound validates covariates

- **WHEN** the adjusted arm's rejection rates over 500 replicates are at or below 0.0695 for `delta`, `angle` and `shape`
- **THEN** the decision is `covariates_validated`

#### Scenario: Excess in any statistic blocks validation

- **WHEN** the adjusted arm's `angle` rejection rate exceeds the bound
- **THEN** the decision is `diagnose_covariate_excess` and the `angle` verdict is `exceeds`

#### Scenario: Unadjusted inflation is reported but not judged

- **WHEN** the unadjusted arm rejects at a rate well above α
- **THEN** the decision file reports that rate and the decision depends on the adjusted arm alone
