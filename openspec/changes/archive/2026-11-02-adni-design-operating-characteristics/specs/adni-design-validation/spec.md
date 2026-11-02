# Spec Delta

## Purpose

Establishes the size and the power of the real-data trajectory analysis at the ADNI case study's design — n = 1,501 in ADNI's sex × stage cells, two blocks, six covariate design columns of ADNI's kinds including a phase-like binary nuisance the covariate-blind basis follows, a weak stage signal and the rank floor — so the case study's reading of its non-rejections rests on a measurement at the regime it is run in instead of on the absence of one.

## ADDED Requirements

### Requirement: Datasets at the ADNI design

The check SHALL generate semi-synthetic datasets at the cell sizes `[[237, 310, 125], [223, 448, 158]]`, three stages, methylation and expression blocks, joint magnitude and a recorded weakened stage effect. Each replicate SHALL draw a numeric covariate confounded with group × stage and inject it, draw a three-level categorical covariate on the same linear predictor without injecting it, draw the phase-like covariate, and draw two independent never-injected numeric covariates.

#### Scenario: Six design columns of ADNI's kinds

- **WHEN** a unit is analysed with the covariates declared
- **THEN** the record's covariate block lists exactly six design columns: the confounded numeric, two drop-first columns of the categorical, one column of the phase-like covariate and the two noise numerics

#### Scenario: Units are reproducible

- **WHEN** the check is run twice with the same base seed and unit indices
- **THEN** the generated datasets, covariates and analysis records are identical

#### Scenario: Every cell of a replicate shares its seeds

- **WHEN** the anchor, the translation control and a power cell are generated for the same replicate index
- **THEN** they share the generator seed, the covariate values, the phase draw and the RRPP seed, and differ only in the surgery

#### Scenario: Seeds are fresh

- **WHEN** units are generated for the paper-grade run
- **THEN** no generator, covariate or RRPP seed equals one used by a committed check at the same replicate index

### Requirement: The phase-like nuisance is drawn per cell and injected into one block

The phase-like covariate SHALL be binary, drawn per sex × stage cell at the ADNI2GO proportions the ADNI cohort recorded (F 0.603 / 0.623 / 0.424, M 0.561 / 0.507 / 0.481), and injected as a dense per-feature offset, drawn once from the base seed and scaled by each feature's pooled SD and a recorded strength, into the expression block only. The methylation block SHALL be untouched by it.

#### Scenario: Cell proportions follow ADNI

- **WHEN** the phase-like covariate is drawn for many replicates at the ADNI table
- **THEN** its mean per sex × stage cell is within Monte Carlo error of the recorded proportion for that cell

#### Scenario: Only the expression block carries the offset

- **WHEN** one replicate is built with the phase strength set to zero and again at the calibrated strength, from the same seeds
- **THEN** the methylation block is identical between the two builds and every expression feature differs by the covariate times its recorded offset

#### Scenario: The realized offset is recorded

- **WHEN** a unit is built
- **THEN** its record carries the median absolute standardized phase difference over the expression features and the expression block's first-principal-component correlation with the phase-like covariate

### Requirement: Calibration places the construction in the ADNI regime

Before any paper-grade run, a recorded calibration on anchor units SHALL meet three criteria: the `adjusted` arm's mean held-out CV stage AUROC lies inside a band declared before the run (default 0.55–0.65); its median first-component |Spearman ρ| with the phase-like covariate is at least a declared threshold (default 0.4); and the `unadjusted` arm rejects at least one statistic above α. Every iteration SHALL be appended to a calibration record with its strengths.

#### Scenario: Criteria met admit the paper-grade run

- **WHEN** a calibration's adjusted-arm mean AUROC is inside the band, its median |ρ| with the phase-like covariate is at least the threshold and its unadjusted arm rejects above α on some statistic
- **THEN** the calibration record marks the iteration as met and the paper-grade run may use its strengths

#### Scenario: Unmet criterion blocks the paper-grade run

- **WHEN** any criterion fails on a calibration
- **THEN** the stage-effect, numeric-injection or phase strength is changed, the iteration is recorded with its values, and the calibration is repeated before any paper-grade run

#### Scenario: Placement descriptives accompany the criteria

- **WHEN** a calibration iteration is recorded
- **THEN** it also carries the selected-rank distribution, the median |ρ| with the confounded numeric covariate, the realized phase offset, the pooled relative eigengap and the `angle` null q95, none of which is judged

### Requirement: One primary analysis per unit on matched seeds

Each paper-grade unit SHALL be analysed once through the real-data analysis with the six covariates in the trajectory model, the rank floor of 2, net displacement, the harness CV settings, 999 permutations and the replicate's RRPP seed, and no latent adjustment or residual correction. The record SHALL carry the statistics, p-values, null summary, selected rank, CV AUROC, first-component correlations and eigengap; a failed unit SHALL be recorded with its error and be unavailable.

#### Scenario: Unit record

- **WHEN** a unit completes
- **THEN** its record line holds the unit key, the seeds, the statistics, p-values, null summary, selected rank, component selection, mean held-out CV AUROC, the first component's |ρ| with the phase-like and the confounded covariates, and the pooled relative eigengap

#### Scenario: Failed unit is unavailable, not a rejection

- **WHEN** the analysis raises an error for a unit
- **THEN** the record carries the error and no p-values, and the unit is excluded from that cell's available replicates

#### Scenario: Interrupted run resumes

- **WHEN** a shard is re-run after an interruption
- **THEN** it skips units already recorded in any shard file, re-runs a torn line, and appends only the missing ones

#### Scenario: Resume under other settings is refused

- **WHEN** an output directory is resumed with a different table, strength, base seed, cell set or permutation count
- **THEN** the check exits with an error and appends nothing

#### Scenario: Headroom is checked before running

- **WHEN** a pool-limited cell would exceed the generator's surgery headroom at the ADNI table
- **THEN** the check refuses to start and names the cell

### Requirement: Fourteen cells at one calibrated strength

The paper-grade run SHALL comprise the zero-effect anchor, the `translation` control at e = 1.0, magnitude at e ∈ {0.1, 0.25, 0.5, 1.0}, and orientation and shape at e ∈ {0.25, 0.5, 0.75, 1.0}, at 500 replicates and 999 permutations by default, all at the calibrated stage-signal strength.

#### Scenario: Cell enumeration

- **WHEN** the run's units are enumerated
- **THEN** there are 14 cells × 500 replicates, replicate-major, and each shard takes every N-th unit

### Requirement: Size is judged at both controls against the one-sided anchor bound

For the anchor and the translation control, the check SHALL compute the rejection rate at α = 0.05 of `delta`, `angle` and `shape` over the R available replicates and compare each with α + 2·√(α(1 − α)/R). Size SHALL hold when every one of the six readings is at or below its bound; a reading with no available replicates SHALL be `unavailable`.

#### Scenario: Bound at the default replicates

- **WHEN** a reading has 500 available replicates
- **THEN** its bound is 0.05 + 2·√(0.0475 / 500) ≈ 0.0695

#### Scenario: One exceeding reading fails size

- **WHEN** the translation control's `delta` rate is 0.072 over 500 replicates and every other reading holds
- **THEN** size does not hold and the decision names the control and statistic

#### Scenario: Bound is computed from available replicates

- **WHEN** some units of a control are unavailable
- **THEN** R is the count of available replicates for that control and statistic, and its bound is computed from it

### Requirement: Power is reported per statistic with an informative effect

For each mode's target statistic (magnitude/`delta`, orientation/`angle`, shape/`shape`), the check SHALL report the rejection rate at every effect with its binomial standard error and the selected-rank distribution, and SHALL name `informative_from`: the smallest effect at which the rate is at least 0.80, or null when none is. Off-target rates and rank-stratified rates SHALL be reported descriptively and SHALL NOT enter the verdict.

#### Scenario: Informative effect found

- **WHEN** orientation/`angle` rejects at 0.42, 0.71, 0.86 and 0.95 at e = 0.25, 0.5, 0.75 and 1.0
- **THEN** `angle`'s `informative_from` is 0.75

#### Scenario: No effect reaches the floor

- **WHEN** shape/`shape` rejects at most 0.31 at any effect
- **THEN** `shape`'s `informative_from` is null

#### Scenario: Rates are stratified by selected rank

- **WHEN** the power cells have been merged
- **THEN** a table lists, per cell and statistic, the rate within each selected-rank bucket (2–3, 4–7, 8 and above) beside the pooled rate

### Requirement: The decision issues one verdict with a measured statement

A pure function over the merged rows SHALL write a decision file with the verdict `power_statement` (size holds and a target statistic has an informative effect), `non_rejections_uninformative` (size holds and none has) or `diagnose_design_excess` (a size reading exceeds), every size reading with its bound, every power rate, the calibration descriptives, the settings hash and, per statistic, a measured statement for the case study. A run below 999 permutations SHALL mark the verdict indicative.

#### Scenario: Power statement

- **WHEN** every size reading holds and `delta`'s `informative_from` is 0.5 while `angle`'s and `shape`'s are null
- **THEN** the verdict is `power_statement`, `delta`'s statement says a non-rejection argues against size differences of e ≥ 0.5 at this design, and `angle`'s and `shape`'s say their non-rejections stay uninformative up to the construction maximum

#### Scenario: Non-rejections uninformative

- **WHEN** every size reading holds and no target statistic reaches 0.80 at any effect
- **THEN** the verdict is `non_rejections_uninformative` and the statements confirm rule 3 with the measured rates

#### Scenario: Size excess

- **WHEN** any size reading exceeds its bound
- **THEN** the verdict is `diagnose_design_excess`, the file names each exceeding control and statistic, and no power statement is issued

#### Scenario: Smoke verdict is indicative

- **WHEN** the decision is computed from a run at 99 permutations
- **THEN** the decision file marks the verdict indicative

### Requirement: Paper-grade run is reported with its consequence and limits

The paper-grade run SHALL be written as a dated report under `docs/reports/` citing its decision file. The report SHALL state the consequence for the case study under its verdict and SHALL carry its limitations: omic analogue blocks for metabolomics (498 features against 298), the intersection cohort's table only, one calibrated strength, additive-linear covariates with one confounded numeric and one confounded factor, and the primary analysis alone.

#### Scenario: Report names the consequence

- **WHEN** the report is written
- **THEN** it names the decision file, the bound, R and the six size readings, lists each statistic's power curve beside its selected-rank distribution, states the replacement or confirmation of the case study's rule 3 (or the diagnosis to run), and lists the limitations

#### Scenario: Limitations are specific

- **WHEN** the report lists its limitations
- **THEN** it states the block split (367 methylation / 131 expression against 133 p180 / 165 NMR) and names the arms it did not run: latent-adjusted, residual-corrected and forced-rank
