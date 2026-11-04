# Spec Delta

## ADDED Requirements

### Requirement: Harness supports the diagnostic stage-plus-residual integration

The harness SHALL accept a diagnostic integration method `pls_residual_pca`. It SHALL return the production stage-PLS scores, at the CV rank with the rank floor, followed by r principal-component scores of the stage-residualized features, taken orthogonal to the PLS span. r is an integration parameter, default `n_stages − 1`. The basis SHALL use the features and stage labels, never the group labels, and SHALL be fitted once, not inside RRPP.

#### Scenario: The leading coordinates are the production scores
- **WHEN** a dataset is evaluated with `pls_residual_pca`
- **THEN** its first columns equal the latent matrix of a `pls` evaluation of the same dataset and parameters, and its last r columns are the residual component scores

#### Scenario: The basis does not use the group labels
- **WHEN** the group labels are permuted within each stage and the dataset is integrated again
- **THEN** the `pls_residual_pca` latent matrix is unchanged

#### Scenario: Residual components are orthogonal to the PLS span
- **WHEN** a dataset is integrated with `pls_residual_pca`
- **THEN** each residual component direction is orthogonal to the span of the PLS rotations in the model's input coordinates, to numerical tolerance

#### Scenario: Infeasible residual rank is rejected
- **WHEN** r is below 1 or above the number of features less the PLS rank
- **THEN** the harness raises a clear validation error naming r and the feasible range

### Requirement: The stage-plus-residual integration is labelled diagnostic

`pls_residual_pca` SHALL be labelled diagnostic in its metadata, which SHALL record the PLS rank, r and the retained dimensionality. Realized geometry SHALL include its latent checkpoint. Attribution and the forced-component override SHALL reject it.

#### Scenario: Metadata and checkpoint
- **WHEN** an evaluation uses `pls_residual_pca`
- **THEN** the metadata identifies the integration as diagnostic and records the PLS rank, r and the retained dimensionality, and realized geometry carries a latent checkpoint named for the method

#### Scenario: Attribution is rejected
- **WHEN** attribution or a forced component count is requested with `pls_residual_pca`
- **THEN** the harness raises a clear validation error

### Requirement: Harness accepts a prebuilt stage-PLS fit for the residual integration

The harness SHALL let a caller build the `pls_residual_pca` integration from an existing `pls` integration of the same dataset and parameters, without refitting the PLS model. The result SHALL equal the integration fitted from scratch.

#### Scenario: Reusing the PLS fit changes nothing
- **WHEN** `pls_residual_pca` is built from a prebuilt `pls` integration and, separately, fitted from scratch on the same dataset and parameters
- **THEN** the two latent matrices and metadata are identical, and the reuse path fits no PLS model

## MODIFIED Requirements

### Requirement: Harness measures span capture on request

The evaluation parameters SHALL accept `span_capture`, absent by default, with `false` normalized to absent. When it is on and the integration is a linear method (`pls`, `pca`, `pls_residual_pca` or `concat`), the harness SHALL compute the `population_latent` checkpoint and the `span_capture` shares from the same fitted basis and preprocessor the trajectory test uses, without refitting. While it is absent, cell identifiers, signatures, records and runtime metadata SHALL be unchanged.

#### Scenario: Span capture reuses the evaluation's fit
- **WHEN** an evaluation runs with `span_capture` on and PLS integration
- **THEN** the PLS model is fitted once, and the realized geometry carries `population_latent` and `span_capture` computed from that model

#### Scenario: Span capture covers the linear diagnostic spaces
- **WHEN** an evaluation runs with `span_capture` on and `pca`, `pls_residual_pca` or `concat` integration
- **THEN** the realized geometry carries `population_latent` and `span_capture` computed from that evaluation's own basis, and for `concat` every share is 1 and the rank equals the number of features

#### Scenario: PLS span capture records are unchanged
- **WHEN** an evaluation runs with `span_capture` on and PLS integration
- **THEN** its `population_latent` checkpoint and `span_capture` block are byte-identical to those produced before span capture covered other methods

#### Scenario: Span capture is inert for other integration methods
- **WHEN** an evaluation runs with `span_capture` on and `snf` or `pls_within_group` integration
- **THEN** the evaluation completes and records no `population_latent` checkpoint and no `span_capture` block

#### Scenario: Absent span capture changes nothing
- **WHEN** an evaluation runs without `span_capture`, or with it set to `false`
- **THEN** its result and cell identity are byte-identical to those produced before this change

#### Scenario: Statistics do not depend on span capture
- **WHEN** the same dataset is evaluated with and without `span_capture`
- **THEN** the latent matrix, the selected rank, the statistics and the p-values are identical
