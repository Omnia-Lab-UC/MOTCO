# Spec Delta

## ADDED Requirements

### Requirement: The comparison with the default per-feature analysis is documented

The measurement-visibility documentation SHALL describe how the production trajectory test compares with a group × stage interaction test on each feature with no integration: what each side tests, which differences each detects, what only the trajectory test can claim, and the limits of the comparison. It SHALL cite the validation's verdict, its characterization grade and its measured rates.

#### Scenario: Reader asks whether integration is needed to detect a difference
- **WHEN** the reader asks whether the trajectory test detects more than a per-feature test
- **THEN** the documentation states the validation's verdict and, by plant family and design, which side leads, with the rates
- **AND** it states that production entered as one test at 0.05 / 3 and the baseline as one test at 0.05

#### Scenario: Reader asks what the per-feature test cannot say
- **WHEN** the reader asks what the trajectory test adds beyond detection
- **THEN** the documentation states that a per-feature test has no geometry and cannot say whether a difference is one of size, orientation or shape, and it states the validation's characterization grade for production

#### Scenario: Reader asks which baseline was used
- **WHEN** the reader asks how the per-feature tests were combined
- **THEN** the documentation names the maximum, the false-discovery reading and the pooled features, states which held size, and states that a moderated statistic was not an arm

#### Scenario: Reader asks how far the result carries
- **WHEN** the reader asks whether the result holds for a real cohort
- **THEN** the documentation states that it was read in simulation at 658 features with Gaussian noise, that real blocks hold far more features, and that no real-noise reading exists on methylation, expression and proteomics
