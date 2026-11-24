# Spec Delta

## ADDED Requirements

### Requirement: What a weak stage effect costs the trajectory test is documented

The measurement-visibility documentation SHALL describe how the production trajectory test behaves as the stage effect weakens and the sample grows: the CV rank, the power on in-span plants, whether a rejection names the kind of a difference, and how it compares with the per-feature test there. It SHALL cite the stage-effect scaling validation's verdict, its grade by design and its measured rates.

#### Scenario: Reader asks whether more samples make up for weak staging
- **WHEN** the reader asks whether a weakly staged cohort is helped by a larger sample
- **THEN** the documentation states the validation's verdict and, by plant family, whether a design with half the stage effect and four times the samples read as its match did, with the rates
- **AND** it states what the per-feature control arm read in the same pairs

#### Scenario: Reader asks where a rejection names the kind of difference
- **WHEN** the reader asks at which designs the rejecting statistic can be read as the kind of difference
- **THEN** the documentation gives the grade by design with each design's stage effect, sample size and CV rank

#### Scenario: Reader places a cohort
- **WHEN** the reader has a cohort's held-out stage AUROC, selected rank and sample size
- **THEN** the documentation gives each simulated design's held-out AUROC, rank and δ√n, and states what was detectable at the designs calibrated to the case studies' AUROC

#### Scenario: Reader asks how far the result carries
- **WHEN** the reader asks whether the result holds for a real cohort
- **THEN** the documentation states that it was read in simulation at 658 features with Gaussian noise, equal groups and three stages, that the calibrated designs match the case studies' AUROC and not their cell tables, covariates or rank, and that no real-noise reading exists on methylation, expression and proteomics
