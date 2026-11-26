# Spec Delta

## ADDED Requirements

### Requirement: Sampling accepts an opt-in noise scale

The generator SHALL accept an optional positive noise scale s for a sampling call. Every sample's deviation from its exact cell mean SHALL be multiplied by s in the units the noise is drawn in (methylation M-values, before the inverse logit), so the noise covariance is s² times the reference covariance and the cell means are unchanged. The scale SHALL draw no random numbers. While it is absent, sampling SHALL be byte-identical to sampling before the option existed.

#### Scenario: Scaled samples are the unscaled samples stretched about their cell means
- **WHEN** one sampling call is made twice at the same seed, once without a noise scale and once with scale s
- **THEN** in analysed units every sample of the second call equals its cell mean plus s times the deviation of the same sample of the first call from that mean, to floating-point precision
- **AND** the random generator is at the same stream position after both calls

#### Scenario: Methylation is scaled before the inverse logit
- **WHEN** methylation is sampled with a noise scale
- **THEN** the deviation is scaled in M-value space and the inverse logit is applied last, so every value stays in `(0, 1)`

#### Scenario: Absent scale is byte-identical
- **WHEN** a sampling call omits the noise scale
- **THEN** the returned matrices equal those produced before the option existed, at the same seed and parameters

#### Scenario: Invalid scale is rejected
- **WHEN** the noise scale is zero, negative or not finite
- **THEN** the generator raises a validation error naming the option
