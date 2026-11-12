# Spec Delta

## ADDED Requirements

### Requirement: Harness measures PLS in the predictive subspace on request

Under the `latent_coordinates` values `predictive`, `predictive_orthonormal` and `predictive_fitted`, the PLS latent matrix SHALL be a coordinate system on the predictive subspace of the final model: the part of its score space spanned by the rows of its response loadings, which has the rank of those loadings (the number of stages minus one, or the fit rank when that is smaller).
- **`predictive`:** the scores projected on an orthonormal basis of that row space.
- **`predictive_orthonormal`:** the model's input coordinates projected on an orthonormal basis of the same subspace in the input space.
- **`predictive_fitted`:** the fitted stage indicators in the model's response units.

Every trajectory statistic, RRPP draw, realized-geometry checkpoint and configuration spectrum SHALL be computed in those coordinates. The model, the cross-validation and the selected or forced rank SHALL be those of the score measurement.

#### Scenario: Measured dimension
- **WHEN** a three-stage dataset is evaluated under `predictive` or `predictive_orthonormal` at a fit rank of five
- **THEN** the latent matrix has two columns
- **AND** under `predictive_fitted` it has three columns and rank two

#### Scenario: The fitted stage labels are kept
- **WHEN** a dataset is evaluated under any predictive value
- **THEN** the fitted stage indicators recomputed from the predictive subspace equal the final model's fitted stage indicators to 1e-8

#### Scenario: Equal to the score measurement at a low fit rank
- **WHEN** the fit rank is at most the number of stages minus one and the dataset is measured under `predictive` and in scores at one RRPP seed
- **THEN** `delta`, `angle` and `shape` agree to 1e-10 and their p-values are equal

#### Scenario: One subspace, three coordinate systems
- **WHEN** the same dataset is evaluated under the three predictive values
- **THEN** the `predictive_orthonormal` and `predictive_fitted` matrices each equal the `predictive` matrix times one fixed matrix of rank equal to the measured dimension, to 1e-8
- **AND** their `span_capture` shares are identical and at most those of the score measurement

#### Scenario: Fit and rank are unchanged
- **WHEN** the same dataset is evaluated in scores and under a predictive value
- **THEN** the selected rank, the cross-validation table and the held-out AUROC are identical

#### Scenario: Statistics do not depend on the orthonormal basis
- **WHEN** the `predictive_orthonormal` latent matrix is multiplied by any orthogonal matrix and measured with the same design and RRPP seed
- **THEN** `delta`, `angle` and `shape` agree with the unrotated measurement to 1e-10, and their p-values are equal

#### Scenario: Population checkpoint in the same coordinates
- **WHEN** span capture is on under a predictive value
- **THEN** the `population_latent` checkpoint is measured in that value's coordinates

#### Scenario: A forced rank
- **WHEN** a predictive value is set together with a forced rank
- **THEN** the predictive subspace is that of the model fitted at the forced rank

#### Scenario: The setting is recorded
- **WHEN** an evaluation runs under a predictive value
- **THEN** its integration parameters and runtime metadata state the value

#### Scenario: A degenerate fit
- **WHEN** the response loadings of the final model have a rank below both the fit rank and the number of stages minus one
- **THEN** the harness raises an evaluation error naming the rank

### Requirement: Harness builds a predictive-subspace integration from a prebuilt PLS fit

The harness SHALL let a caller build the integration of any predictive value from an existing score-coordinate `pls` integration of the same dataset and parameters, without refitting. The result SHALL equal the integration fitted from scratch with the setting. An evaluation SHALL refuse a prebuilt integration whose coordinates differ from its own setting.

#### Scenario: Reusing the fit changes nothing
- **WHEN** a predictive integration is built from a prebuilt `pls` integration and, separately, fitted from scratch on the same dataset and parameters
- **THEN** the two latent matrices and metadata are identical, and the reuse path fits no PLS model

#### Scenario: Three values from one fit
- **WHEN** the three predictive integrations are built from one prebuilt `pls` integration
- **THEN** exactly one PLS model was fitted for the four measurements

#### Scenario: Mismatched prebuilt integration
- **WHEN** an evaluation set to one coordinate system receives a prebuilt integration in another
- **THEN** the harness raises a validation error naming both

#### Scenario: Not a score-coordinate PLS fit
- **WHEN** the prebuilt integration is not a `pls` integration with fitted artifacts, or is not in score coordinates
- **THEN** the harness raises a validation error naming what it received

## MODIFIED Requirements

### Requirement: The coordinate setting is normalized and validated before fitting

`latent_coordinates` SHALL accept absent (the score coordinates), `scores` (normalized to absent), `orthonormal`, `predictive`, `predictive_orthonormal` and `predictive_fitted`. The harness SHALL reject the setting before any fit when:
- the value is unknown;
- the integration method is not `pls`;
- attribution is requested together with any value other than the score coordinates;
- a predictive value is set on a dataset with fewer than three stages.

Each rejection SHALL name the offending value.

#### Scenario: Scores normalize to absent
- **WHEN** the evaluation parameters name `scores`
- **THEN** the parameters equal those with the setting absent

#### Scenario: Unknown value
- **WHEN** the evaluation parameters name an unknown coordinate system
- **THEN** the harness raises a validation error naming the value and the allowed values

#### Scenario: Other integration methods
- **WHEN** `orthonormal` or a predictive value is set with `pca`, `concat`, `snf`, `pls_within_group` or `pls_residual_pca` integration
- **THEN** the harness raises a validation error naming the method before any preprocessing or fit

#### Scenario: Attribution
- **WHEN** `orthonormal` or a predictive value is set and attribution is requested
- **THEN** the harness raises a validation error before any fit, stating that attribution is defined only in score coordinates

#### Scenario: Two stages
- **WHEN** a predictive value is set on a two-stage dataset
- **THEN** the harness raises a validation error before any fit, naming the value and the stage count

### Requirement: PLS fitting remains single-pass within an evaluation

The harness SHALL fit one final pooled PLS estimator after component selection. It SHALL measure trajectories on that estimator's training scores or, under another coordinate setting, on a fixed linear image of that estimator's training inputs: their projection on an orthonormal basis of its rotations' span, or a coordinate system on its predictive subspace. It SHALL use that same estimator for attribution projection and reconstruction.

#### Scenario: Scores and attribution share one final fit

- **WHEN** a PLS evaluation with attribution succeeds
- **THEN** the score matrix equals the training scores of the estimator supplied to attribution
- **AND** the selected component count and score values remain equivalent to a PLS evaluation with attribution disabled under the same inputs and seeds

#### Scenario: Orthonormal coordinates reuse the final fit

- **WHEN** a PLS evaluation runs under orthonormal span coordinates
- **THEN** exactly one final PLS estimator is fitted, and the latent matrix is that estimator's training inputs projected on the orthonormal basis of its rotations' span

#### Scenario: Predictive coordinates reuse the final fit

- **WHEN** a PLS evaluation runs under a predictive value
- **THEN** exactly one final PLS estimator is fitted, and the latent matrix is computed from that estimator's training scores, response loadings and rotations alone
