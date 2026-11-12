# Spec Delta

## ADDED Requirements

### Requirement: Real-data analysis offers predictive-subspace coordinates

The analysis parameters SHALL accept a measurement coordinate setting (`latent_coordinates`: absent, `predictive`, `predictive_orthonormal` or `predictive_fitted`). When set, the latent space is fitted exactly as without it, and the trajectory statistics, the permutation draws, the latent scores and the latent stage means are computed in that coordinate system on the fitted model's predictive subspace, as the simulation harness defines it. An unknown value, and a value set on fewer than three stages, SHALL be rejected before any fitting. The setting SHALL be omitted from the resolved parameters while absent.

#### Scenario: Predictive analysis
- **WHEN** a three-stage analysis whose cross-validation selects five components requests `predictive`
- **THEN** the selected rank and held-out AUROC equal those of the analysis without the setting
- **AND** the latent scores and stage means have two columns, and the PLS loadings table is that of the five-component fit

#### Scenario: Equivalence with the harness
- **WHEN** the same data are analysed through the real-data entry point and through the simulation harness under the same predictive value, parameters and seed
- **THEN** the latent matrices, statistics, p-values, rank and null summary are identical

#### Scenario: Absent setting
- **WHEN** no coordinate setting is requested
- **THEN** every output equals the pre-change output, and the resolved parameters carry no coordinate key

#### Scenario: Unknown value
- **WHEN** the analysis parameters name an unknown coordinate system, or `orthonormal`
- **THEN** the analysis raises a validation error naming the value and the allowed values before any fitting

#### Scenario: Two stages
- **WHEN** a predictive value is set on a dataset with two stages
- **THEN** the analysis raises a validation error naming the value and the stage count before any fitting

### Requirement: Real-data analysis records the measurement coordinates

When the coordinate setting is set, the record SHALL carry `measurement.latent_coordinates` with the value and the measured dimension, and `latent.integration_params.latent_coordinates`. Neither SHALL appear while absent.

#### Scenario: Record contents
- **WHEN** a three-stage analysis under `predictive` is recorded
- **THEN** `measurement.latent_coordinates` names `predictive` and a measured dimension of two, beside the selected rank

### Requirement: Predictive-subspace coordinates are rejected with attribution and the split interaction test

The analysis SHALL reject a predictive coordinate setting requested together with attribution or with the split interaction test, before any fitting, naming both settings. Both are defined on the fitted span.

#### Scenario: Attribution
- **WHEN** a predictive value is set and attribution is enabled
- **THEN** the analysis raises a validation error naming both settings before any fitting

#### Scenario: Split interaction test
- **WHEN** a predictive value is set and the split interaction test is requested
- **THEN** the analysis raises a validation error naming both settings before any fitting
