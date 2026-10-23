# Spec Delta

## MODIFIED Requirements

### Requirement: Harness supports initial integration methods
The harness SHALL construct the molecular latent space — the measurement substrate in which trajectory geometry is estimated — via a selectable integration method, operating on M-value-converted methylation and raw expression/proteomics. The production latent-space methods are **SNF** (graph-spectral embedding) and **PLS** (the transform of the omic features into the subspace that maximizes covariance with the stage label). `concat` is retained as a **baseline/diagnostic** path (standardized feature concatenation), not a constructed latent space. The viz down-projection (`plot_trajectory_from_*`) is display-only and distinct from this measurement space.

An optional integration parameter SHALL select which omic blocks are measured: a non-empty, duplicate-free subset of methylation, expression, and proteomics. Only the selected blocks SHALL enter pooled preprocessing, the integration method, orientation attribution, and every realized-geometry checkpoint (population, standardized population, observed standardized, and latent). Every joint scope SHALL be computed over the selected blocks only, and per-block scopes SHALL be reported for the selected blocks only. Blocks that are not selected SHALL be generated but not measured. When the parameter is absent, all three blocks SHALL be measured, and the evaluation result SHALL be identical to the harness's output before the parameter existed.

#### Scenario: Concatenated baseline integration
- **WHEN** the caller selects `concat` integration
- **THEN** the harness converts methylation to M-values, standardises all selected layers, and concatenates them into the outcome matrix
- **AND** the result metadata identifies `concat` as a baseline rather than a constructed latent space

#### Scenario: SNF integration
- **WHEN** the caller selects `snf` integration
- **THEN** the harness converts methylation to M-values and creates the latent space from SNF fusion of the selected layers and spectral embedding

#### Scenario: PLS integration
- **WHEN** the caller selects `pls` integration
- **THEN** the harness converts methylation to M-values, standardises the selected layers, fits PLS-DA conditioned on the stage label, and returns the PLS X-score matrix as the latent space
- **AND** the number of latent variables is selected by the double nested cross-validation (modal LV across repeats, parsimony tie-break) from the candidate ranks at or above the rank floor (2 by default) to secure a stable, non-overfitted molecular space in which orientation and shape are defined
- **AND** the result metadata records the selected number of latent variables and the cross-validation parameters

#### Scenario: PLS integration is infeasible
- **WHEN** the caller selects `pls` integration but the sample provides too few observations per stage for the cross-validation
- **THEN** the harness raises a clear validation error

#### Scenario: Unsupported integration method
- **WHEN** the caller selects an unsupported integration method
- **THEN** the harness raises a clear validation error

#### Scenario: A block subset is measured
- **WHEN** the caller selects the layers methylation and expression
- **THEN** the outcome matrix, the latent space, every realized-geometry checkpoint's joint scope, and every attribution feature record contain methylation and expression features only
- **AND** realized geometry reports no proteomics per-block scope
- **AND** the result metadata records the selected layers in canonical order

#### Scenario: Absent layer selection is byte-identical
- **WHEN** the evaluation parameters omit the layer selection
- **THEN** the evaluation result equals the result produced before the selection existed, for the same dataset and parameters

#### Scenario: Invalid layer selection is rejected
- **WHEN** the layer selection is empty, names an unknown layer, or repeats a layer
- **THEN** the harness raises a clear validation error naming the allowed layers

### Requirement: PLS integration supports a diagnostic forced-component override

The harness SHALL accept an optional integration parameter that forces the number of PLS latent variables, bypassing cross-validated component selection, for rank diagnostics only. When the override is absent the behavior MUST be unchanged. When the override is set, the harness MUST skip the double nested cross-validation, fit the single pooled PLS estimator at the forced component count, and record in the integration metadata both the forced count and an explicit marker that component selection was overridden, so no downstream consumer can mistake a forced-rank run for a production evaluation. The override MUST be validated against the feasible component range, which is at least 2 and at most the feasible maximum regardless of the cross-validation rank floor, and rejected with a clear error otherwise.

#### Scenario: Default behavior is unchanged

- **WHEN** a PLS evaluation runs without the forced-component parameter
- **THEN** the number of latent variables is selected by the double nested cross-validation exactly as before
- **AND** the integration metadata carries no override marker

#### Scenario: Forced component count is honored and recorded

- **WHEN** a PLS evaluation runs with the forced-component parameter set to a feasible value
- **THEN** the latent space has exactly that many components and no cross-validation is performed
- **AND** the integration metadata records the forced count and marks component selection as overridden

#### Scenario: Infeasible forced count is rejected

- **WHEN** the forced-component parameter exceeds the feasible range for the sample
- **THEN** the harness raises a clear validation error

#### Scenario: Forced rank is independent of the floor

- **WHEN** a PLS evaluation runs with the forced-component parameter set to 2 and a cross-validation rank floor of 3
- **THEN** the latent space has exactly 2 components, no cross-validation is performed, and the floor is not applied

## ADDED Requirements

### Requirement: PLS rank selection respects a floor of 2 by default

The evaluation parameters SHALL accept an integer rank floor for PLS integration, default 2. Cross-validation SHALL consider only candidates at or above the floor and below the maximum, so every outer-fold selection, the modal rank and the pooled fit are at or above it. A floor below 1, or a non-integer, SHALL be rejected with a clear error. The floor is part of every evaluation's identity and applies to PLS only.

#### Scenario: Default floor lifts a rank-1 selection

- **WHEN** a PLS evaluation runs without naming a floor on an input where the unfloored rule selected rank 1
- **THEN** the selected rank is at least 2 and the latent space has at least 2 components

#### Scenario: Named floor is honoured

- **WHEN** a PLS evaluation names a floor of 3
- **THEN** every cross-validation selection and the final pooled fit have at least 3 components

#### Scenario: Invalid floor is rejected

- **WHEN** the evaluation parameters name a floor of 0 or a non-integer floor
- **THEN** the harness raises a clear validation error

### Requirement: Maximum rank is clamped above the floor

The cross-validation maximum SHALL be clamped to the range from the floor plus one to the feasible maximum (the smaller of the feature and sample counts). When the feasible maximum is below the floor plus one, the harness SHALL raise a clear error naming the floor and the counts rather than fit at a lower rank.

#### Scenario: Maximum at the floor is raised

- **WHEN** a PLS evaluation names a maximum of 2 at the default floor on a sample with at least 3 features and samples
- **THEN** the cross-validation considers rank 2 only, the recorded maximum is 3, and the selected rank is 2

#### Scenario: Sample too small for the floor

- **WHEN** a PLS evaluation runs at the default floor on a sample with 2 features
- **THEN** the harness raises a clear validation error naming the floor and the feature count

### Requirement: Floor of 1 is invisible; a floor above 1 is recorded

When the floor is 1, it SHALL be omitted from the persisted record, the parameter signature and the cell identity, so that floor-1 evaluation results, signatures and cell ids equal those produced before the floor existed. When the floor is above 1, by default or by name, the record's runtime metadata and the latent integration parameters SHALL state it.

#### Scenario: Floor of 1 reproduces the pre-floor evaluation

- **WHEN** the evaluation parameters name a floor of 1
- **THEN** the evaluation result, its persisted record, its parameter signature and its cell id equal those produced before the floor existed

#### Scenario: Default floor is recorded

- **WHEN** a PLS evaluation runs at the default floor
- **THEN** the runtime metadata states the floor and the latent integration parameters state it beside the selected rank
