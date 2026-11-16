# Spec Delta

## ADDED Requirements

### Requirement: Real-data analysis offers within-stage permutation draws

The analysis parameters SHALL accept a permutation strata setting (`permutation_strata`: absent or `stage`). When it is `stage`, every permutation draw of the trajectory test SHALL exchange the reduced model's residuals only among samples of the same stage, as the simulation harness defines it, with the covariate-adjusted models and the residual correction of the analysis. The latent fit, the selected rank, the observed statistics, the latent scores and the stage means SHALL be those of the analysis without the setting. The setting SHALL combine with covariates, the latent adjustment, the residual correction, attribution, block weighting, sparsity and the measurement coordinates. An unknown value SHALL be rejected before any fitting. The setting SHALL be omitted from the resolved parameters while absent.

#### Scenario: Only the null moves

- **WHEN** one dataset is analysed with the same parameters and seed, once without the setting and once with `stage`
- **THEN** the selected rank, the held-out AUROC, the observed statistics, the latent scores, the loadings and the stage means are identical
- **AND** the p-values and the null summary of the second analysis are those of draws restricted within stage

#### Scenario: Equivalence with the harness

- **WHEN** the same data are analysed through the real-data entry point and through the simulation harness under `stage`, with the same parameters and seed
- **THEN** the latent matrices, statistics, p-values, rank and null summary are identical

#### Scenario: The case-study configuration equals the diagnosis's construction

- **WHEN** a dataset with covariates is analysed under `stage` and the `df` residual correction, and the predictive size-excess diagnosis's within-stage loop runs on that analysis's latent matrix, covariate-adjusted models and seed
- **THEN** the analysis's p-values equal the loop's exactly

#### Scenario: Absent setting

- **WHEN** no permutation strata are requested
- **THEN** every output equals the pre-change output, and the resolved parameters carry no strata key

#### Scenario: Unknown value

- **WHEN** the analysis parameters name the strata `group`
- **THEN** the analysis raises a validation error naming the value and the allowed values before any fitting

### Requirement: Real-data analysis records the permutation strata

When the permutation strata setting is set and permutations are drawn, the record SHALL carry `measurement.permutation_strata`, naming `stage` and each stage's sample count in the declared stage order. The key SHALL be absent while the setting is absent.

#### Scenario: Record contents

- **WHEN** an analysis of 80 samples in stages of 20, 20 and 40 runs under `stage` with permutations
- **THEN** `measurement.permutation_strata` names `stage` and the counts 20, 20 and 40 in the declared order

#### Scenario: No permutations

- **WHEN** the setting is `stage` and no permutation is requested
- **THEN** the analysis completes, the resolved parameters carry the setting, and the record has no `measurement.permutation_strata`

### Requirement: Within-stage draws are rejected with the split interaction test and the interaction scan

The analysis SHALL reject the permutation strata setting requested together with the split interaction test or with the per-feature interaction scan, before any fitting, naming both settings. Those tests draw their own unrestricted nulls, which is how they were validated.

#### Scenario: Split interaction test

- **WHEN** `stage` is set and the split interaction test is requested
- **THEN** the analysis raises a validation error naming both settings before any fitting

#### Scenario: Interaction scan

- **WHEN** `stage` is set and the interaction scan is requested
- **THEN** the analysis raises a validation error naming both settings before any fitting

### Requirement: motco analyze exposes within-stage draws

`motco analyze` SHALL accept `--permutation-strata stage`, absent by default, and pass it to the analysis. The command's summary SHALL name the strata when they are set. A rejected combination SHALL exit with the analysis's message.

#### Scenario: Flag given

- **WHEN** `motco analyze` runs with permutations and `--permutation-strata stage`
- **THEN** `result.json` carries `measurement.permutation_strata`, and its p-values equal those of the analysis called with the same setting and seed

#### Scenario: Flag absent

- **WHEN** `motco analyze` runs without the flag
- **THEN** every written file equals the pre-change output

#### Scenario: Flag with the interaction scan

- **WHEN** `motco analyze` runs with `--permutation-strata stage` and `--interaction-scan`
- **THEN** it exits with a message naming both settings and writes no output
