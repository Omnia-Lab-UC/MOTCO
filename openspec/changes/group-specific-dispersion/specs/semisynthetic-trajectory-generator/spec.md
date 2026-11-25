# Spec Delta

## ADDED Requirements

### Requirement: Group dispersion scale is an opt-in generator setting

The generator SHALL accept `group_dispersion_scale`, a positive finite number s that multiplies the noise standard deviation of the second group in `group_labels`, in every omic block and at every stage. The first group SHALL keep the reference noise. The setting SHALL be accepted under every trajectory mode and construction. A value of 1 SHALL normalize to absent, and any other value that is not a positive finite number SHALL be rejected with an error naming the setting.

#### Scenario: The second group's noise is scaled and the first group's is not
- **WHEN** a dataset is generated with `group_dispersion_scale` 2 and another without the setting, at the same seed and parameters
- **THEN** the first group's samples are identical in the two datasets
- **AND** every sample of the second group, in analysed units (methylation in M-values), equals its population cell mean plus 2 times the deviation of the same sample of the other dataset from that mean, to floating-point precision

#### Scenario: Population means do not move
- **WHEN** a dataset is generated with a dispersion scale
- **THEN** its population trajectories and its differential indicators equal those of the dataset generated without the setting at the same seed

#### Scenario: The scale combines with an effect
- **WHEN** a dataset is generated with a non-zero group effect and a dispersion scale
- **THEN** the second group's population means are those of the same dataset without the scale, and its noise is scaled about them

#### Scenario: Invalid scale is rejected
- **WHEN** `group_dispersion_scale` is zero, negative, not finite or not a number
- **THEN** the generator raises a validation error naming the setting

### Requirement: The dispersion scale is recorded and invisible while absent

Truth SHALL record the dispersion scale and the group it applies to when the setting is present. While it is absent, or 1, every dataset, truth record, cell identifier, signature and matched seed SHALL equal those produced before the setting existed. The scale is part of the design and not of the group transform, so a zero-effect cell SHALL keep it.

#### Scenario: Absent scale is byte-identical
- **WHEN** a dataset is generated without `group_dispersion_scale`, or with it set to 1
- **THEN** its samples, truth, cell identifier and signature are byte-identical to those produced before this change at the same parameters and seed

#### Scenario: Truth names the scale
- **WHEN** a dataset is generated with `group_dispersion_scale` 1.5 and group labels A and B
- **THEN** its truth records the scale 1.5, the covariance multiplier 2.25 and the group B it applies to

#### Scenario: A set scale changes the cell identity
- **WHEN** two cells differ only in `group_dispersion_scale`, absent in one and 1.5 in the other
- **THEN** their cell identifiers and signatures differ

#### Scenario: A zero-effect cell keeps the scale
- **WHEN** a null cell is derived from generator parameters that set the scale
- **THEN** the null cell's parameters carry the same scale, and its dataset has the second group's noise scaled

### Requirement: Stage-span scales use each group's own noise variance

Under `stage_span` with a dispersion scale s, the null and the planted population pooled standard deviations SHALL take each feature's within-cell variance as the reference variance in the first group's cells and s² times it in the second group's, weighted by the realized cell sizes. Without a dispersion scale both pooled deviations SHALL equal those computed before the setting existed.

#### Scenario: The within-cell term is weighted by group
- **WHEN** a stage-span dataset is generated with dispersion scale s, n_A samples in the first group and n_B in the second
- **THEN** the within-cell term of every feature's pooled variance is (n_A + s²·n_B) / (n_A + n_B) times its reference variance

#### Scenario: A null plant is the null dataset
- **WHEN** a stage-span dataset is generated with mode `none`, or with an effect of 0, and a dispersion scale
- **THEN** it equals the dataset generated without a construction at the same seed, parameters and scale

#### Scenario: Stage-span datasets without the scale are unchanged
- **WHEN** a stage-span dataset is generated without `group_dispersion_scale` under either span scale
- **THEN** its plant, samples and truth are byte-identical to those produced before this change
