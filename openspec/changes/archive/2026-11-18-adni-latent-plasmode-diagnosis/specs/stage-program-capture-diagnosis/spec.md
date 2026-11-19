# Spec Delta

## Purpose

Diagnoses a stage-program capture estimator's calibration failure from the run's own records, without new units. It says whether the bias lies in the estimator's parts or in the ratio, whether the spread is the cohort's or the units', and whether the cohort's own interval for the capture is bounded. The reading is stated beside the parent's verdict and does not change it.

## ADDED Requirements

### Requirement: The diagnosis reads the run's records under the parent's bindings

The diagnosis SHALL read the latent-plasmode run's merged records for both profiles where the run wrote them, and SHALL bind the parent's committed settings, calibration table, capture table and decision by sha256. It SHALL refuse a changed parent file, a record written under other settings, and a missing or duplicated unit, before any reading.

#### Scenario: A changed parent file is refused
- **WHEN** a committed parent file's hash differs from the one in the diagnosis's settings
- **THEN** `run` and `decide` exit naming the file

#### Scenario: A record under other settings is refused
- **WHEN** a record's settings hash differs from its profile's committed settings
- **THEN** `run` exits before any reading

#### Scenario: A missing unit is refused
- **WHEN** a profile's records lack an expected unit or hold one twice
- **THEN** `run` exits before any reading

### Requirement: The records reproduce the committed tables

Before any reading, the diagnosis SHALL rebuild the parent's calibration table and the cohort's capture table from the records with the parent's own builders and writer, and SHALL compare them with the committed files byte for byte.

#### Scenario: Reproduction holds
- **WHEN** both rebuilt tables equal the committed files
- **THEN** the reading proceeds and the records' sha256 are recorded

#### Scenario: Reproduction fails
- **WHEN** either rebuilt table differs
- **THEN** the verdict is `reproduction_failed`

### Requirement: The diagnosis records the calibration's parts, aggregate and spread

For every calibration cell, space and stage reading, the diagnosis SHALL record the bias of the numerator and of the denominator against their truths, the ratio of sums over all units against the truth's, and the within- and between-cohort covariances of the numerator and denominator.

#### Scenario: Unbiased parts
- **WHEN** the numerator and denominator are unbiased for their truths
- **THEN** their mean differences lie within the stated multiple of their cohort-level standard errors

#### Scenario: The aggregate
- **WHEN** the ratio of sums over all of a cell's units is read
- **THEN** it is recorded with the truth's ratio and a cluster-bootstrap standard error and interval over cohorts

#### Scenario: A spread fixed by the cohort
- **WHEN** most of the variance of a cohort's mean denominator is between cohorts
- **THEN** the between share is at least the stated limit

#### Scenario: Aggregates only
- **WHEN** the outputs are written
- **THEN** they hold no participant identifier, no per-unit value and no per-cohort value

### Requirement: The diagnosis records the cohort's own interval for the capture

For every cohort point, space and stage reading, the diagnosis SHALL record Fieller's 95 % interval for the capture from the cohort's mean numerator and denominator, under four stand-ins for the cohort's between-cohort covariance.

#### Scenario: An unbounded interval
- **WHEN** the cohort's mean denominator is not more than 1.96 of its standard deviations from zero
- **THEN** its interval is recorded as `exclusive` or `whole_line`, not `bounded`

#### Scenario: Every stand-in is recorded
- **WHEN** the reading is written
- **THEN** each point carries the interval's kind and `t` under all four stand-ins, and whether they agree

### Requirement: The reading is stated with its order and does not change the parent's verdict

The diagnosis SHALL issue one of `reproduction_failed`, `estimator_biased`, `resolvable`, `units_limited` or `information_limited`, in that order of precedence, with what each means for a rule version 2. The decision file SHALL state that the rule was written after an exploratory pass, and SHALL carry the parent's verdict unchanged.

#### Scenario: Information-limited
- **WHEN** the parts and the aggregate are unbiased, the spread is mostly between cohorts, and the cohort's interval is unbounded at every point under the primary stand-in
- **THEN** the verdict is `information_limited`, and no rule version 2 is proposed on that cohort

#### Scenario: Resolvable
- **WHEN** the parts and the aggregate are unbiased and the cohort's interval is bounded at some point under the primary stand-in
- **THEN** the verdict is `resolvable`, and a rule version 2 is left to its own change
