# Spec Delta

## ADDED Requirements

### Requirement: Harness measures PLS in orthonormal span coordinates on request

The evaluation parameters SHALL accept `latent_coordinates`. Under `orthonormal`, the PLS latent matrix SHALL be the final model's input coordinates projected on an orthonormal basis of its rotations' span, and every trajectory statistic, RRPP draw, realized-geometry checkpoint and configuration spectrum SHALL be computed in those coordinates. The model, the cross-validation and the selected or forced rank SHALL be those of the score measurement.

#### Scenario: Same span, orthonormal coordinates
- **WHEN** a dataset is evaluated with `latent_coordinates` set to `orthonormal`
- **THEN** the latent matrix spans the same subspace as the score matrix of the same fit, and equals the input coordinates times the PLS weights up to column signs, to 1e-8

#### Scenario: Fit and rank are unchanged
- **WHEN** the same dataset is evaluated in score and in orthonormal coordinates
- **THEN** the selected rank, the cross-validation table and the held-out AUROC are identical

#### Scenario: Statistics do not depend on the orthonormal basis
- **WHEN** the orthonormal latent matrix is multiplied by any orthogonal matrix and measured with the same design and RRPP seed
- **THEN** `delta`, `angle` and `shape` agree with the unrotated measurement to 1e-10, and their p-values are equal

#### Scenario: Population checkpoint in the same coordinates
- **WHEN** span capture is on under orthonormal coordinates
- **THEN** the `population_latent` checkpoint is measured in the orthonormal coordinates
- **AND** the `span_capture` shares are identical to those of the score measurement

#### Scenario: The setting is recorded
- **WHEN** an evaluation runs under orthonormal coordinates
- **THEN** its integration parameters and runtime metadata state `latent_coordinates: "orthonormal"`

### Requirement: The coordinate setting is normalized and validated before fitting

`latent_coordinates` SHALL accept absent (the score coordinates), `scores` (normalized to absent) and `orthonormal`. The harness SHALL reject the setting before any fit when:
- the value is unknown;
- the integration method is not `pls`;
- attribution is requested together with `orthonormal`.

Each rejection SHALL name the offending value.

#### Scenario: Scores normalize to absent
- **WHEN** the evaluation parameters name `scores`
- **THEN** the parameters equal those with the setting absent

#### Scenario: Unknown value
- **WHEN** the evaluation parameters name an unknown coordinate system
- **THEN** the harness raises a validation error naming the value and the allowed values

#### Scenario: Other integration methods
- **WHEN** `orthonormal` is set with `pca`, `concat`, `snf`, `pls_within_group` or `pls_residual_pca` integration
- **THEN** the harness raises a validation error naming the method before any preprocessing or fit

#### Scenario: Attribution
- **WHEN** `orthonormal` is set and attribution is requested
- **THEN** the harness raises a validation error before any fit, stating that attribution is defined only in score coordinates

### Requirement: An absent coordinate setting is invisible

While `latent_coordinates` is absent, it SHALL be omitted from the persisted record, the runtime metadata, the parameter signature and the cell identity. Every evaluation result, signature and cell id SHALL be identical to those produced before the setting existed.

#### Scenario: Committed identities are unchanged
- **WHEN** every committed study configuration is enumerated after this change
- **THEN** its cell ids, signatures and matched seeds equal the committed identity fixtures

#### Scenario: Results are unchanged
- **WHEN** a dataset is evaluated without the setting, or with `scores`
- **THEN** the result is byte-identical to the result produced before this change

### Requirement: Harness builds the orthonormal integration from a prebuilt PLS fit

The harness SHALL let a caller build the orthonormal-coordinate integration from an existing score-coordinate `pls` integration of the same dataset and parameters, without refitting. The result SHALL equal the integration fitted from scratch with the setting. An evaluation SHALL refuse a prebuilt integration whose coordinates differ from its own setting.

#### Scenario: Reusing the fit changes nothing
- **WHEN** the orthonormal integration is built from a prebuilt `pls` integration and, separately, fitted from scratch on the same dataset and parameters
- **THEN** the two latent matrices and metadata are identical, and the reuse path fits no PLS model

#### Scenario: Mismatched prebuilt integration
- **WHEN** an evaluation set to one coordinate system receives a prebuilt integration in the other
- **THEN** the harness raises a validation error naming both

#### Scenario: Not a PLS fit
- **WHEN** the prebuilt integration is not a `pls` integration with fitted artifacts, or already carries orthonormal coordinates
- **THEN** the harness raises a validation error naming what it received

## MODIFIED Requirements

### Requirement: PLS fitting remains single-pass within an evaluation

The harness SHALL fit one final pooled PLS estimator after component selection. It SHALL measure trajectories on that estimator's training scores or, under orthonormal span coordinates, on its training inputs projected on an orthonormal basis of its rotations' span. It SHALL use that same estimator for attribution projection and reconstruction.

#### Scenario: Scores and attribution share one final fit

- **WHEN** a PLS evaluation with attribution succeeds
- **THEN** the score matrix equals the training scores of the estimator supplied to attribution
- **AND** the selected component count and score values remain equivalent to a PLS evaluation with attribution disabled under the same inputs and seeds

#### Scenario: Orthonormal coordinates reuse the final fit

- **WHEN** a PLS evaluation runs under orthonormal span coordinates
- **THEN** exactly one final PLS estimator is fitted, and the latent matrix is that estimator's training inputs projected on the orthonormal basis of its rotations' span
