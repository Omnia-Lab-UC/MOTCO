# cohort-stage-program-capture Specification

## Purpose
Estimates the share of a cohort's real stage program that lies inside a fitted measurement space, from independent held-out halves of the cohort's exported omic blocks. The estimator is calibrated on synthetic cohorts with the truth known, and a predeclared decision says whether the capture is measured.

## Requirements

### Requirement: Units are the cohort anchor's samples

The check SHALL read the ADNI export in place through the cohort anchor's input contract and SHALL draw each unit's sample exactly as the cohort real-noise anchor does, at 80, 200, 500 and 750 participants. Nothing SHALL be written to the export's directory.

#### Scenario: A unit is an earlier run's unit
- **WHEN** a unit is built at 80, 200 or 500 participants at a replicate
- **THEN** its sample equals the cohort real-noise anchor's unit at that point and replicate

#### Scenario: A contract violation is named
- **WHEN** the export violates the anchor's input contract
- **THEN** the check exits naming the violation and the file, and no unit is evaluated

#### Scenario: The export is not copied
- **WHEN** a run, a merge or a decision completes
- **THEN** no file of the export has been copied or modified

### Requirement: Each unit has three disjoint parts

A unit SHALL consist of its sample, on which every space is fitted, and of two halves of the participants outside it, split within stage. The check SHALL draw 20 such splits per unit from seeds of its own family. No participant SHALL be in two parts of one split.

#### Scenario: The parts partition the cohort
- **WHEN** a unit's sample and one of its splits are built
- **THEN** the sample and the two halves are pairwise disjoint and together hold every participant of the cohort

#### Scenario: A split keeps the stages
- **WHEN** a split is drawn
- **THEN** each stage's participants outside the sample are divided between the halves with sizes differing by at most one

### Requirement: The stage program is group-blind and covariate-adjusted

Each half's stage configuration SHALL be the stage coefficients of its features on stage indicators and the case study's covariates, centred over stages, in the sample's standardized coordinates. The group label SHALL NOT enter the configuration or any fit, and no analysis SHALL run under the real labels.

#### Scenario: The labels do not change a reading
- **WHEN** a unit is evaluated under two different within-stage permutations of the group labels
- **THEN** every recorded reading is equal

#### Scenario: A rank-deficient half fails its unit
- **WHEN** a half's stage and covariate design lacks full column rank
- **THEN** the unit is recorded as failed with the reason and enters no reading

### Requirement: The capture is a ratio of summed cross-products

The capture of a space at a design point SHALL be the sum over units of the halves' cross-product inside the space divided by the sum over units of their cross-product in total. Per-unit ratios SHALL be reported for description only, and the estimate SHALL NOT be clipped.

#### Scenario: Unbiased with the truth known
- **WHEN** the estimator is applied to independent Gaussian halves around a known configuration and a fixed space
- **THEN** its mean over repetitions equals the configuration's true share inside the space within sampling error, and the share of a single noisy configuration does not

#### Scenario: A noise-free configuration is recovered exactly
- **WHEN** both halves equal one configuration
- **THEN** the capture equals that configuration's share inside the space

#### Scenario: An estimate outside the unit interval is kept
- **WHEN** a pooled estimate is below 0 or above 1
- **THEN** it is reported unchanged

### Requirement: The capture is read by stage contrast

The check SHALL report the capture in total and for two stage contrasts whose energies sum to the total: the net displacement (last stage minus first) and the interior of the centred stage space. It SHALL report the program's energy and the eigenvalues of the halves' stage-space cross-product.

#### Scenario: The contrasts partition the energy
- **WHEN** the total, net-displacement and interior energies are computed at a point
- **THEN** the two contrast energies sum to the total

### Requirement: Estimation noise is read beside the program

For every space the check SHALL report the share of the halves' difference inside the space, and the share of the whole held-out configuration and of the sample's own configuration.

#### Scenario: The difference holds no program
- **WHEN** both halves equal one configuration plus independent noise
- **THEN** the difference's share does not depend on the configuration

### Requirement: Production's span is the designated space

The check SHALL fit production's space at the case study's configuration, cross-validated once per unit, and SHALL read its fitted span as the designated space. It SHALL also record the predictive subspace of the same fit, the fit forced at ranks 2, 8 and 19, and a principal-component space of the sample at the unit's cross-validated rank. No permutation SHALL be drawn.

#### Scenario: One cross-validation per unit
- **WHEN** a unit is evaluated
- **THEN** exactly one cross-validated fit is made and every other fit is at a forced rank

#### Scenario: The predictive subspace shares the fit
- **WHEN** a unit's fitted span and predictive subspace are read
- **THEN** both come from the same fitted model, and the predictive subspace has one dimension fewer than the number of stages

### Requirement: The calibration plants a known program

The calibration profile SHALL evaluate synthetic cohorts that keep every participant's stage, group and covariates and draw features from a Gaussian model with the cohort's covariate coefficients and residual covariance and the cohort's stage configuration at strengths 0, 0.67 and 1. Each unit SHALL be evaluated as a cohort unit is, and its true capture SHALL be recorded.

#### Scenario: The truth is exact
- **WHEN** a synthetic cohort is drawn at strength 1 with its noise scaled to zero
- **THEN** every space's estimated capture equals its true capture

#### Scenario: Strength 0 plants nothing
- **WHEN** a synthetic cohort is drawn at strength 0
- **THEN** its planted stage configuration is zero and no true capture is defined

#### Scenario: Synthetic data is not kept
- **WHEN** a calibration run completes
- **THEN** no synthetic feature value is in a record or a written file

### Requirement: The calibration reports bias, spread and the null level

The calibration SHALL report, per strength, point, space and reading, the mean estimate, the mean truth, their difference and the standard deviation of the pooled estimate across cohorts, and at strength 0 the standard deviation of each energy.

#### Scenario: The spread is across cohorts
- **WHEN** the calibration table is written
- **THEN** each standard deviation is taken over one pooled estimate per synthetic cohort

### Requirement: Units are checked and earlier runs are reproduced

A unit that fails a consistency check SHALL be marked inconsistent and SHALL enter no reading. A run SHALL NOT start without a passed pre-flight. The merge SHALL compare the cross-validated rank counts at 80, 200 and 500 participants with the cohort anchor's committed counts, and the median naive captures at 750 with the gain-and-specificity run's committed values.

#### Scenario: An inconsistent unit leaves the readings
- **WHEN** a unit's parts overlap, a basis is not orthonormal or a cross-product disagrees with its matrix
- **THEN** the unit is marked inconsistent with the failing check and is excluded from every estimate

#### Scenario: A differing reproduction is recorded
- **WHEN** a rank count or a naive capture differs from its committed value
- **THEN** the reproduction is recorded as failed with the point, both values and the quantity

#### Scenario: A run without a pre-flight is refused
- **WHEN** a run is started for a profile with no passed pre-flight under the same settings
- **THEN** it exits without evaluating a unit

### Requirement: Settings and parent files are bound

Each profile's settings SHALL be bound by hash and SHALL carry the sha256 of the export's files and of the committed parent tables. A resume, merge or decision under other settings, or with a changed export or parent file, SHALL be refused.

#### Scenario: Other settings are refused
- **WHEN** a run is resumed with a different replicate count
- **THEN** it exits naming the differing keys

#### Scenario: A changed parent is refused by name
- **WHEN** a committed parent table differs from its bound hash
- **THEN** every subcommand exits naming the file

### Requirement: The decision is declared before any unit runs

The decision SHALL carry a rule version and SHALL return one of `gate_failed`, `no_program`, `unresolved`, `partly_measured` or `measured`. It SHALL give each designated capture an interval of twice the calibration's standard deviation at strength 0.67 and SHALL mark a decision below the committed replicate counts as indicative.

#### Scenario: A biased estimator fails the gate
- **WHEN** the calibration's absolute bias of the designated capture exceeds 0.05 at a point at strength 0.67 or 1
- **THEN** the verdict is `gate_failed` and no capture is stated

#### Scenario: A failed reproduction fails the gate
- **WHEN** a reproduction is recorded as failed
- **THEN** the verdict is `gate_failed`

#### Scenario: No program
- **WHEN** the gate holds and the total energy at 80 participants is at most twice its standard deviation at strength 0
- **THEN** the verdict is `no_program`

#### Scenario: Resolution decides the verdict
- **WHEN** the gate holds, the program is present, and the designated capture's half-width is at most 0.15 at four, at one to three, or at none of the points
- **THEN** the verdict is `measured`, `partly_measured` or `unresolved`

#### Scenario: Contrast grades
- **WHEN** a contrast's energy exceeds twice its standard deviation at strength 0 at a point
- **THEN** the contrast is graded present there, and captured, missed or undetermined as its capture's interval lies above 0.5, below 0.5 or across it

#### Scenario: A verdict states what it licenses
- **WHEN** a decision is written
- **THEN** it names the captures that may be stated and whether a follow-up with plants is licensed for each contrast

### Requirement: Records and committed files hold no participant-level value

Records SHALL hold seeds, ranks, scalars and stage-by-stage matrices only. No record or committed file SHALL hold a participant identifier, a label vector, the membership of a sample or of a half, a stage configuration or a synthetic feature.

#### Scenario: A record is searched
- **WHEN** a unit's record is written
- **THEN** it contains no identifier of the export and no field with one value per participant or per feature
