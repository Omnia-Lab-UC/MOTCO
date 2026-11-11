# Spec Delta

## ADDED Requirements

### Requirement: The scan by block and its reading are documented

The measurement-visibility documentation SHALL describe the by-block scan: when to use it, what a by-block hit claims, and how detection moves between blocks. It SHALL cite the by-block validation's verdict and measured rates.

#### Scenario: Reader asks when to scan by block
- **WHEN** the reader asks whether to scan jointly or by block
- **THEN** the documentation states that the joint scan's familywise bar is set by the block that most often gives the largest null F, cites the measured null-maximum shares, and recommends the by-block scan when the blocks differ in size or noise

#### Scenario: Reader asks what a by-block hit means
- **WHEN** the reader consults the documentation about a feature with by-block adjusted p ≤ 0.05
- **THEN** it states that the familywise error over all scanned features is still controlled at 0.05, now split equally across the blocks, so each block has the same chance of a false hit whatever its size

#### Scenario: Reader asks what the by-block scan costs
- **WHEN** the reader asks what the by-block scan gives up
- **THEN** the documentation states that detection moves from the block with more features or heavier-tailed noise to the others, cites the validation's per-block recall in both arms, and states that only two blocks were measured
