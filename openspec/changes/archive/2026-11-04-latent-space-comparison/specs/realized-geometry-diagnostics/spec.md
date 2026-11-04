# Spec Delta

## MODIFIED Requirements

### Requirement: Population latent checkpoint on request

When span capture is requested for an evaluation with a linear integration method (`pls`, `pca`, `pls_residual_pca` or `concat`) of a dataset with population trajectories, the diagnostics SHALL add a `population_latent` checkpoint. It is the population cell means standardized by the replicate's fitted preprocessor, concatenated over the measured blocks, mapped by the evaluation's own fitted linear map to its latent coordinates, and measured with the evaluation's orientation estimator.

#### Scenario: Checkpoint measures the truth through the fitted basis
- **WHEN** span capture is on, the integration is linear and the dataset carries population trajectories
- **THEN** realized geometry contains `population_latent` with each group's path length and the pairwise `delta`, `angle` and `shape`, under the same availability rules as every other checkpoint

#### Scenario: PLS checkpoint is unchanged
- **WHEN** span capture is on and the integration is PLS
- **THEN** the `population_latent` checkpoint equals the population means projected through the fitted PLS model, byte-identical to the checkpoint produced before other methods were covered

#### Scenario: Checkpoint is absent otherwise
- **WHEN** span capture is off, the integration is `snf` or `pls_within_group`, or the dataset has no population trajectories
- **THEN** realized geometry contains no `population_latent` checkpoint

### Requirement: Span capture shares on request

With span capture on, realized geometry SHALL carry a `span_capture` block, computed in the fitted map's input coordinates with an orthonormal basis of the directions the map keeps: the PLS rotations, the principal axes, the PLS rotations with the residual axes, or every feature for `concat`. `group_difference` is the share of Σ_s‖d_s‖² inside the span, with d_s the population B − A difference at stage s minus its stage mean. `pooled_stage` is the share of the centred average stage means. The block also records `rank`, `n_features` and `random_share` = rank / n_features.

#### Scenario: Shares are reported for a planted difference
- **WHEN** span capture is on for a linear evaluation of a dataset whose population B − A difference varies across stages
- **THEN** `span_capture.group_difference` and `span_capture.pooled_stage` are numbers in [0, 1], beside `rank`, `n_features` and `random_share`

#### Scenario: Shares agree with the definition for any basis
- **WHEN** span capture is computed for a fixed set of population means against a known orthonormal basis
- **THEN** each share equals the squared norm of the projection divided by the squared norm of the vectors, whichever integration method supplied the basis

#### Scenario: Concat keeps everything
- **WHEN** span capture is on for a `concat` evaluation
- **THEN** `group_difference` (when defined), `pooled_stage` and `random_share` are 1 and `rank` equals `n_features`

#### Scenario: Share is unavailable when the difference is constant
- **WHEN** the population B − A difference is the same at every stage, as at the anchor or under a translation
- **THEN** `span_capture.group_difference` is null and the other fields are still reported

#### Scenario: Records are unchanged without span capture
- **WHEN** span capture is off
- **THEN** the realized geometry record is byte-identical to the record produced before this change
