# Spec Delta

## ADDED Requirements

### Requirement: Stage-span scale is an opt-in generator setting

The generator SHALL accept a `span_scale` setting, `null` (the default) or `planted`, valid only under `construction: "stage_span"`. `null` SHALL normalize to absent, and while it is absent every dataset, truth record, cell identifier, signature and matched seed SHALL equal those produced before the setting existed.

#### Scenario: Default scale is byte-identical
- **WHEN** a stage-span dataset is generated without `span_scale`, or with `span_scale` set to `null`
- **THEN** its samples, truth, cell identifier and signature are byte-identical to those produced before this change at the same parameters and seed

#### Scenario: Scale without the stage-span construction is rejected
- **WHEN** `span_scale` is `planted` while `construction` is absent or `surgery`
- **THEN** the generator raises a validation error naming the setting

#### Scenario: Unknown scale is rejected
- **WHEN** `span_scale` is any value other than `null` or `planted`
- **THEN** the generator raises a validation error naming the offending value

#### Scenario: Planted scale with an alignment is rejected
- **WHEN** `span_scale` is `planted` and `span_alignment` is set
- **THEN** the generator raises a validation error stating that the two settings are not supported together

### Requirement: Planted-scale stage-span edits are exact in the planted pooled scale

Under `span_scale: "planted"`, the stage-span edit SHALL be made in coordinates where each span feature is divided by the planted population's pooled SD. That SD combines the within-cell variance with the between-cell variance of the planted cell means of both groups at the realized cell sizes. Because it depends on the plant, the plant SHALL be solved as a fixed point.

#### Scenario: The edit holds exactly in the planted coordinates
- **WHEN** a stage-span dataset is generated under `planted` for `magnitude`, `orientation`, `shape` or `translation` with effect e
- **THEN** with A's and B's population stage means each divided by the planted pooled SD, B's configuration is the mode's edit of A's configuration
- **AND** the intended `delta`, `angle` and `shape` computed in those coordinates equal the mode's definition to within 1e-8 relative to A's path length

#### Scenario: The fixed point converges or fails loudly
- **WHEN** the fixed-point iteration does not reach a relative change in the planted SD of at most 1e-12 within its iteration limit
- **THEN** the generator raises a validation error naming the mode, the effect and the last residual, rather than returning an inexact plant

#### Scenario: The planted plant consumes no randomness
- **WHEN** a planted-scale stage-span dataset and a `none` dataset are generated at the same seed and design
- **THEN** group A's samples are identical
- **AND** group B's samples in analysed units differ from the null dataset's by exactly the plant of their stage, to floating-point precision

#### Scenario: Zero effect yields the null dataset
- **WHEN** `trajectory_mode` is `none` or the effect is 0 under `planted`
- **THEN** the dataset equals the null dataset at the same seed

### Requirement: Planted-scale truth names its scale and convergence

Truth for a planted-scale stage-span dataset SHALL record `span_scale: "planted"`, state that its stage-span units are planted-standardized, and record the fixed point's iteration count and final residual. The intended geometry, path length, mean step and plant norm SHALL be reported in the planted-standardized coordinates.

#### Scenario: Truth records the planted scale
- **WHEN** a planted-scale stage-span dataset is generated
- **THEN** its truth names `span_scale: "planted"`, the planted-standardized units, the iteration count and the residual
- **AND** its intended geometry is the geometry of the planted-standardized population configurations of A and B

#### Scenario: Truth without the setting is unchanged
- **WHEN** a stage-span dataset is generated without `span_scale`
- **THEN** its truth carries no `span_scale` entry and reports its geometry in the null-standardized units, as before this change
