# Spec Delta

## ADDED Requirements

### Requirement: The two permutation nulls and their reading are documented

The measurement-visibility documentation SHALL state what each permutation null of the trajectory test assumes and holds fixed, and SHALL report the within-stage pilot: the verdict, the size readings of both nulls where stage dispersion is equal and where it is unequal, the paired power and specificity readings on the production scores, and the re-read of `predictive_orthonormal` under the same draws.

#### Scenario: Reader asks what the two nulls assume

- **WHEN** the reader consults the section on the permutation null
- **THEN** it states that the unrestricted draws assume the reduced model's residuals are exchangeable across all samples, that the within-stage draws assume it within each stage only, and that the within-stage draws keep every stage's pooled mean at its fitted value

#### Scenario: Reader asks when the choice matters

- **WHEN** the reader asks in which data the two nulls differ
- **THEN** the documentation reports both nulls' rates on the dispersion cells, with the direction of the unrestricted null's error when the largest stage is the narrowest and when it is the widest

#### Scenario: Reader asks which null to use

- **WHEN** the reader asks whether a case study may report within-stage p-values
- **THEN** the documentation states the pilot's verdict and its consequence, that the unrestricted draws remain the default, and that the split interaction test and the scan are drawn unrestricted

#### Scenario: Reader asks whether the predictive coordinates are reopened

- **WHEN** the reader asks what the pilot means for the predictive values
- **THEN** the documentation states the candidate's reading against production under the same draws, and whether a paper grade is proposed
