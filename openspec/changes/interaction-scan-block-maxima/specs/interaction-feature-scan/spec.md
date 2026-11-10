# Spec Delta

## MODIFIED Requirements

### Requirement: The global p-value is the permutation p-value of the maximum F

With permutations and without block labels, the scan SHALL report the maximum F over the non-degenerate features, the feature attaining it, and the global p-value (1 + #{draws whose maximum null F ≥ the observed maximum}) / (1 + permutations).

#### Scenario: Global p-value
- **WHEN** the scan runs with B permutations
- **THEN** the global p-value equals (1 + the number of draws whose largest null F over the non-degenerate features is at least the observed maximum F) / (1 + B)

### Requirement: Per-feature p-values are step-down max-T adjusted

With permutations and without block labels, each non-degenerate feature SHALL receive a Westfall–Young step-down max-T adjusted p-value. Features are ordered by observed F, descending, with ties broken by column index. At rank j, the raw value is (1 + #{draws whose maximum null F over the features ranked j or lower ≥ the j-th observed F}) / (1 + B). The adjusted p-value is the running maximum of these raw values down the order.

#### Scenario: Matches a brute-force step-down
- **WHEN** the scan runs on a small case with null distributions returned
- **THEN** every adjusted p-value equals a direct implementation of the step-down rule on the returned null F

#### Scenario: Adjusted p-values are monotone in F
- **WHEN** feature a has a larger observed F than feature b
- **THEN** a's adjusted p-value is at most b's

#### Scenario: The top feature's adjusted p-value is the global p-value
- **WHEN** the scan runs with permutations
- **THEN** the adjusted p-value of the feature with the maximum F equals the global p-value

## ADDED Requirements

### Requirement: Blocks share the familywise error equally on request

Given one block label per column, the scan SHALL apply the step-down rule within each block on the same draws, and SHALL give each feature the adjusted p-value K times its within-block value, capped at 1, with K the number of blocks holding a non-degenerate feature. The global p-value SHALL be K times the smallest within-block global p-value, capped at 1.

#### Scenario: Matches a brute-force per-block step-down
- **WHEN** the scan runs with block labels on a small case of two blocks with null distributions returned
- **THEN** every adjusted p-value equals min(1, 2 × the step-down rule applied to that feature's block alone on the returned null F)

#### Scenario: The global p-value is the smallest adjusted p-value
- **WHEN** the scan runs with block labels and permutations
- **THEN** the global p-value equals the smallest adjusted p-value, and equals min(1, K × the smallest block global p-value)

#### Scenario: Adjusted p-values are monotone in F within a block
- **WHEN** features a and b are in the same block and a has the larger observed F
- **THEN** a's adjusted p-value is at most b's

### Requirement: One block reproduces the joint scan

The block-wise combination with a single block SHALL give the joint scan's global p-value, adjusted p-values and hits, and the scan without block labels SHALL be unchanged.

#### Scenario: One label
- **WHEN** the same scan runs once without labels and once with the same label on every column, at the same seed
- **THEN** the global p-values, every adjusted p-value and the hits are identical

#### Scenario: No labels
- **WHEN** the scan runs without labels
- **THEN** every field and its serialization equal those of the scan before block labels existed, and no block field is present

### Requirement: The block-wise result reports every block

With block labels, the result SHALL give the combination, K, and per block in order of first appearance: the maximum F and its column, and with permutations the block's own global p-value and every draw's block maximum. A block whose features are all degenerate SHALL have no maximum or p-value and SHALL not count in K. The joint null maxima SHALL still be returned.

#### Scenario: Per-block maxima
- **WHEN** the scan runs with two block labels and B permutations
- **THEN** each block reports its maximum F, the column attaining it, its global p-value and B null maxima, each the largest null F over that block's non-degenerate columns at that draw

#### Scenario: A wholly degenerate block
- **WHEN** every column of one of three blocks is constant
- **THEN** that block has no maximum or p-value, K is 2, and the other blocks' adjusted p-values are 2 times their within-block values, capped at 1

### Requirement: Malformed block labels are rejected before any draw

The scan SHALL raise an error naming both lengths when the block labels do not give exactly one label per outcome column.

#### Scenario: Wrong length
- **WHEN** 9 labels are given for 10 outcome columns
- **THEN** the scan raises an error naming 9 and 10, before any draw
