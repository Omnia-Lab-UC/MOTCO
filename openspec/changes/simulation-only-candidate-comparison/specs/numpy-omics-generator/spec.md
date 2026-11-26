# Spec Delta

## ADDED Requirements

### Requirement: Sampling accepts one noise scale per cell

A sampling call SHALL accept, in place of a single noise scale, one positive finite scale per cell. Each cell's samples SHALL be stretched about that cell's exact mean by its scale, in the units the noise is drawn in, after the draw. No random number SHALL be drawn for it, and a scale of 1 for every cell SHALL give the samples drawn without a scale.

#### Scenario: Per-cell scales stretch each cell about its mean
- **WHEN** a call is made with one scale per cell at a fixed generator state
- **THEN** each cell's samples equal the unscaled call's samples moved to mean + scale × (sample − mean) for that cell, methylation in M-values before the inverse logit
- **AND** the generator ends at the same stream position as the unscaled call

#### Scenario: Equal scales match the single scale
- **WHEN** every cell carries the same scale s
- **THEN** the samples equal those of the call with the single noise scale s

#### Scenario: Invalid scales are rejected
- **WHEN** the scales' count differs from the number of cells, or one is not a positive finite number
- **THEN** the call raises naming the scales, before any draw
