# Spec Delta

## ADDED Requirements

### Requirement: Stage dispersion scale is an opt-in generator setting

The generator SHALL accept `stage_dispersion_scale`, one positive finite multiplier of the noise standard deviation per stage, applied to both groups in every block. The population trajectories and indicators SHALL NOT move, and no random number SHALL be drawn for it. It SHALL compose with every mode and construction, and a zero-effect cell SHALL keep it.

#### Scenario: The dataset is the equal-dispersion one, stretched by stage
- **WHEN** a dataset is generated with the setting at a seed
- **THEN** it equals the dataset without the setting at that seed with each stage's samples stretched about their exact cell means by that stage's multiplier
- **AND** the population trajectories equal those without the setting

#### Scenario: It composes with the group dispersion scale
- **WHEN** both `stage_dispersion_scale` and `group_dispersion_scale` are set
- **THEN** each cell's noise standard deviation is the reference one times the product of its stage's and its group's multipliers

#### Scenario: Invalid values are rejected
- **WHEN** the setting's length differs from the number of stages, or a value is not a positive finite number
- **THEN** generation raises naming the setting and the value, before any draw

### Requirement: The stage dispersion scale is recorded and invisible while absent

When set, the truth record SHALL gain a `stage_dispersion` block with the multiplier and the covariance multiplier per stage. A setting of all ones SHALL normalize to the absent setting. While absent, every dataset, truth record, cell id, signature, matched seed and dumped configuration SHALL be unchanged.

#### Scenario: Truth names the scales
- **WHEN** a dataset is generated with the setting
- **THEN** its truth holds each stage's multiplier and its square

#### Scenario: All ones is the absent setting
- **WHEN** the setting is given as all ones
- **THEN** the parameters equal those without it, and the cell id and signature are those without it

#### Scenario: Committed identities are unchanged
- **WHEN** the committed study configurations are enumerated
- **THEN** every cell id, signature and matched seed equals its committed value

### Requirement: Stage-span scales use each cell's own noise variance under a stage dispersion scale

Under the stage-span construction with a stage dispersion scale, the null and planted pooled standard deviations SHALL take the within-cell variance of each feature as the reference variance times the mean, over cells weighted by their realized sizes, of each cell's squared multiplier. The planted edit SHALL then be exact in the scale the generated data have.

#### Scenario: The pooled scale matches the population
- **WHEN** a stage-span dataset is generated with stage multipliers and unequal stage sizes
- **THEN** the pooled scale used by the plant equals the population pooled standard deviation of the generated cells

#### Scenario: Equal multipliers reduce to a uniform scale
- **WHEN** every stage carries the same multiplier s
- **THEN** the within-cell term is s² times the reference variance
