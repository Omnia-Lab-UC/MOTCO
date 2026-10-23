# rank-floor-validation Specification

## Purpose
Establishes the size and power of the rank-2 floor against the unfloored rank rule on a small-n construction calibrated to select rank 1, so the default floor is adopted on evidence from the regime where it changes the measurement.

## Requirements

### Requirement: Null datasets are calibrated to select rank 1 without the floor

The check SHALL generate null semi-synthetic datasets at the Phase 6 small-n design with the stage signal weakened so that, on a recorded smoke run, the unfloored rule selects rank 1 in at least half of the replicates before any paper-grade run. The construction, its strengths, the calibration result and the seeds SHALL be recorded.

#### Scenario: Calibration is reported

- **WHEN** the paper-grade run completes
- **THEN** the decision file reports the fraction of unfloored null replicates that selected rank 1 and the construction settings that produced it

#### Scenario: Calibration below the threshold blocks the paper-grade run

- **WHEN** the smoke run's rank-1 fraction in the unfloored arm is below one half
- **THEN** the stage signal is weakened further, the settings are recorded, and the smoke run is repeated before any paper-grade run

#### Scenario: Replicates are reproducible

- **WHEN** the check is run twice with the same parent seed and replicate indices
- **THEN** the generated datasets and the analysis records are identical

### Requirement: Two arms per replicate on matched data

Each replicate SHALL be analysed through the real-data analysis twice, with the same cross-validation settings, permutation count and RRPP seed: at a floor of 1 and at a floor of 2. Each arm SHALL record its p-values, selected rank and null summary. The check SHALL support resumable, sharded execution and merge shard records into one table of per-arm, per-replicate p-values and selected ranks.

#### Scenario: Arms share the data and differ in the floor

- **WHEN** a replicate is analysed in both arms
- **THEN** both arms read the same blocks and metadata, the floor-2 arm's selected rank is at least 2, and the floor-1 arm's selected rank is whatever the unfloored rule chose

#### Scenario: Rank-1 null is reported as degenerate

- **WHEN** the floor-1 arm selects rank 1
- **THEN** its `angle` null summary is recorded, and the merged table marks the replicate's `angle` as measured in one dimension

#### Scenario: Interrupted run resumes

- **WHEN** a shard is re-run after an interruption
- **THEN** it skips replicates already recorded and appends only the missing ones

### Requirement: Floor-2 arm must hold the one-sided anchor bound

For the zero-effect anchor and the `translation` negative control, the check SHALL compute the floor-2 arm's rejection rate at α = 0.05 for `delta`, `angle` and `shape` over the R available replicates and compare each with `α + k·√(α(1 − α)/R)` at k = 2. The floor-2 arm SHALL be admissible when every rate is at or below the bound. The floor-1 arm's rates SHALL be reported and SHALL NOT enter the decision.

#### Scenario: Floor-1 arm does not enter the decision

- **WHEN** the floor-1 arm rejects above the bound on any statistic
- **THEN** its rates are reported and the decision depends on the floor-2 arm alone

#### Scenario: Bound is computed from available replicates

- **WHEN** some replicates are unavailable in the floor-2 arm
- **THEN** R is the count of available replicates and the bound is computed from it

### Requirement: Power is measured in both arms

The check SHALL run magnitude, orientation and shape effects at the Phase 6 power profile's effect levels in both arms, on matched datasets and permutation draws, and SHALL write per-cell rejection rates per arm, their difference with its pooled Monte Carlo standard error, and the selected-rank distribution per arm and cell. The power rows are descriptive and SHALL NOT enter the decision.

#### Scenario: Rates and rank distributions are tabulated

- **WHEN** the paper-grade run completes
- **THEN** the rates table holds, for every mode, effect and statistic, both arms' rates, their difference, the pooled standard error and each arm's rank distribution

### Requirement: Decision names adopt_floor or diagnose_floor_excess

The decision file SHALL state the per-arm rates at the anchor and the translation control, the bound, R, the per-statistic verdicts, the calibration result, and the decision: `adopt_floor` when the floor-2 arm is admissible, otherwise `diagnose_floor_excess`.

#### Scenario: Floor holds

- **WHEN** the floor-2 arm's anchor and translation rates over 500 replicates are at or below 0.0695 for `delta`, `angle` and `shape`
- **THEN** the decision is `adopt_floor`

#### Scenario: Floor exceeds

- **WHEN** the floor-2 arm's anchor `angle` rate exceeds the bound
- **THEN** the decision is `diagnose_floor_excess` and the `angle` verdict is `exceeds`

### Requirement: Findings are reported with the limitation

Each paper-grade run SHALL produce a dated report under `docs/reports/` stating the construction and calibration, the rates and bound, the rank distributions, the power differences, the decision, and the limitation that the floor guarantees a two-dimensional space, not a second axis with signal.

#### Scenario: Report cites the decision

- **WHEN** the report is written
- **THEN** it cites the decision file and states the decision verbatim
