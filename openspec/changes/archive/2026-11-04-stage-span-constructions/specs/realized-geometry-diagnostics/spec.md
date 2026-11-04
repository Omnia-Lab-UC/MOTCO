# Spec Delta

## ADDED Requirements

### Requirement: Population latent checkpoint on request

When span capture is requested for a PLS evaluation of a dataset with population trajectories, the diagnostics SHALL add a `population_latent` checkpoint. It is the population cell means standardized by the replicate's fitted preprocessor, concatenated over the measured blocks, projected through the replicate's fitted PLS model, and measured with the evaluation's orientation estimator.

#### Scenario: Checkpoint measures the truth through the fitted basis
- **WHEN** span capture is on, the integration is PLS and the dataset carries population trajectories
- **THEN** realized geometry contains `population_latent` with each group's path length and the pairwise `delta`, `angle` and `shape`, under the same availability rules as every other checkpoint

#### Scenario: Checkpoint is absent otherwise
- **WHEN** span capture is off, the integration is not PLS, or the dataset has no population trajectories
- **THEN** realized geometry contains no `population_latent` checkpoint

### Requirement: Span capture shares on request

With span capture on, realized geometry SHALL carry a `span_capture` block, computed in the fitted model's input coordinates with the orthonormal span of its rotations. `group_difference` is the share of Σ_s‖d_s‖² inside the span, with d_s the population B − A difference at stage s minus its stage mean. `pooled_stage` is the share of the centred average stage means. The block also records `rank`, `n_features` and `random_share` = rank / n_features.

#### Scenario: Shares are reported for a planted difference
- **WHEN** span capture is on for a PLS evaluation of a dataset whose population B − A difference varies across stages
- **THEN** `span_capture.group_difference` and `span_capture.pooled_stage` are numbers in [0, 1], beside `rank`, `n_features` and `random_share`

#### Scenario: Share is unavailable when the difference is constant
- **WHEN** the population B − A difference is the same at every stage, as at the anchor or under a translation
- **THEN** `span_capture.group_difference` is null and the other fields are still reported

#### Scenario: Records are unchanged without span capture
- **WHEN** span capture is off
- **THEN** the realized geometry record is byte-identical to the record produced before this change
