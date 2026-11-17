# Spec Delta

## ADDED Requirements

### Requirement: Real-noise size on a second cohort is documented with its scope

The measurement-visibility documentation SHALL report the size readings on ADNI's blocks under within-stage-permuted labels: the design points, the verdict on the production scores and on `predictive_orthonormal` under each null, the fixed-rank reading against the in-sample stage fit, and the scope a real-noise size reading has.

#### Scenario: Reader asks whether the predictive subspace's excess is specific to SEA-AD
- **WHEN** the reader consults the measured-dimension section
- **THEN** it states the ADNI anchor's verdict, `predictive_orthonormal`'s size scope beside production's, and the SEA-AD anchor's counts beside each ADNI design point's

#### Scenario: Reader asks what produces the excess
- **WHEN** the reader asks whether an over-fitted stage fit is enough to produce it
- **THEN** the documentation states the fixed-rank grade, the readings at ranks 8 and 19 against the in-sample stage R² at each design point, and whether the within-stage draws remove the excess there

#### Scenario: Reader asks whether production's size has a second real-noise reading
- **WHEN** the reader asks on which real data the production test's size was read
- **THEN** the documentation names both exports, states production's size scope on ADNI's design points, and states that the full-cohort rates are conditional on one fit

#### Scenario: Reader asks what a size reading on one cohort licenses
- **WHEN** the reader asks whether a hold on ADNI admits the predictive values elsewhere
- **THEN** the documentation states that a size reading speaks for its export and design points, that it does not lift the SEA-AD reading, and that reopening an option at a design also needs a gain and specificity there
