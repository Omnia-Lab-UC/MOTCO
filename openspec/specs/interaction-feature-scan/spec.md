# interaction-feature-scan Specification

## Purpose
Tests the group × stage interaction one feature at a time and combines the features by their maximum, so a difference confined to a few features can be detected and the features carrying it named under familywise error control.

## Requirements

### Requirement: The scan computes an interaction F for every feature

Given an outcome matrix and a full and a nested reduced model matrix, the scan SHALL report for each outcome column its interaction sum of squares, its full-model residual sum of squares and its F = (SS_int / df_int) / (RSS / df_res), with df_int = rank(F) − rank(R) and df_res = n − rank(F).

#### Scenario: Definitions match direct computation
- **WHEN** the scan runs on an outcome Y, full model F and reduced model R
- **THEN** for every column y the interaction sum of squares equals ‖(H_F − H_R)·y‖², the residual sum of squares equals ‖(I − H_F)·y‖², and the F equals their ratio scaled by the degrees of freedom, to a relative tolerance of 1e-10

#### Scenario: Each feature's F ignores its location and scale
- **WHEN** one outcome column is replaced by a·y + b with a ≠ 0
- **THEN** that column's F is unchanged to a relative tolerance of 1e-10, and every other column's F is unchanged

#### Scenario: The per-feature sums add to the split test's total part
- **WHEN** the scan and the split interaction test run on the same outcome and models
- **THEN** the scan's interaction and residual sums of squares, summed over the features, equal the split test's `total` interaction and residual sums of squares to a relative tolerance of 1e-10

### Requirement: Degenerate features are excluded and reported

A feature whose full-model residual sum of squares is zero to a relative tolerance of 1e-12 of its centred sum of squares (a constant feature included) SHALL be marked degenerate. It SHALL have no F and no p-value, and it SHALL be excluded from the observed and null maxima. The result SHALL report the number of degenerate features and which they are.

#### Scenario: Constant feature
- **WHEN** one outcome column is constant
- **THEN** that column is reported as degenerate with no F or p-value, and the maximum and every other column's results equal those of the scan without that column

#### Scenario: Every feature degenerate
- **WHEN** every outcome column is degenerate
- **THEN** the scan raises an error stating that no feature has a defined F, before any draw

### Requirement: Permutation follows RRPP's residual randomization

With a positive permutation count, the scan SHALL permute the reduced model's residuals by Freedman–Lane, drawing the same index sequence as RRPP's serial path at the same seed, and SHALL evaluate every feature's F at every draw from the same permuted residuals. With zero permutations, it SHALL report the observed F only.

#### Scenario: Draws match RRPP
- **WHEN** the scan and RRPP run on the same outcome and models with the same seed and permutation count
- **THEN** every permuted outcome of the scan equals the reduced fitted values plus the reduced residuals in the row order RRPP's serial path uses at that draw

#### Scenario: Null F match a brute-force refit
- **WHEN** null distributions are returned for a small case
- **THEN** each draw's per-feature null F equals the F computed by refitting both models to Ŷ_R + E[π], with π that draw's permutation, to a relative tolerance of 1e-10

#### Scenario: No permutations
- **WHEN** the permutation count is 0
- **THEN** the result carries the per-feature F, the maximum and its feature, and no p-values

### Requirement: The global p-value is the permutation p-value of the maximum F

With permutations, the scan SHALL report the maximum F over the non-degenerate features, the feature attaining it, and the global p-value (1 + #{draws whose maximum null F ≥ the observed maximum}) / (1 + permutations).

#### Scenario: Global p-value
- **WHEN** the scan runs with B permutations
- **THEN** the global p-value equals (1 + the number of draws whose largest null F over the non-degenerate features is at least the observed maximum F) / (1 + B)

### Requirement: Per-feature p-values are step-down max-T adjusted

With permutations, each non-degenerate feature SHALL receive a Westfall–Young step-down max-T adjusted p-value. Features are ordered by observed F, descending, with ties broken by column index. At rank j, the raw value is (1 + #{draws whose maximum null F over the features ranked j or lower ≥ the j-th observed F}) / (1 + B). The adjusted p-value is the running maximum of these raw values down the order.

#### Scenario: Matches a brute-force step-down
- **WHEN** the scan runs on a small case with null distributions returned
- **THEN** every adjusted p-value equals a direct implementation of the step-down rule on the returned null F

#### Scenario: Adjusted p-values are monotone in F
- **WHEN** feature a has a larger observed F than feature b
- **THEN** a's adjusted p-value is at most b's

#### Scenario: The top feature's adjusted p-value is the global p-value
- **WHEN** the scan runs with permutations
- **THEN** the adjusted p-value of the feature with the maximum F equals the global p-value

### Requirement: The null maxima are returned

With permutations, the scan SHALL return the maximum null F of every draw and the column index attaining it, so a caller can check how the null maxima fall across blocks of features. The full per-feature null F SHALL be returned only on request.

#### Scenario: Null maxima returned
- **WHEN** the scan runs with B permutations
- **THEN** the result carries B null maxima and B column indices, each index naming a non-degenerate column whose null F at that draw equals that draw's maximum
- **AND** it carries the B × p null F only when null distributions were requested

### Requirement: Residual corrections rescale the permuted residuals

The scan SHALL accept RRPP's residual corrections (`df`, `leverage`, or none) and apply them to the reduced residuals before any draw. The observed F SHALL be unchanged by the correction. An unknown or undefined correction SHALL be rejected before any draw.

#### Scenario: df leaves every p-value unchanged
- **WHEN** the same scan runs with no correction and with `df` at the same seed
- **THEN** every null F, the global p-value and every adjusted p-value are identical

#### Scenario: Leverage correction
- **WHEN** the scan runs with `leverage`
- **THEN** each reduced residual row is divided by √(1 − h_ii) of the reduced model before permutation, and the observed F equal the uncorrected run's

#### Scenario: Unknown correction
- **WHEN** the correction `uniform` is given
- **THEN** the scan raises an error naming the value and the allowed corrections before any draw

### Requirement: Invalid inputs are rejected before any draw

The scan SHALL raise an error naming the problem when the outcome and the models disagree in rows, when any of them is non-finite, when the reduced model does not lie in the full model's column space, when df_int or df_res is not positive, or when the permutation count is negative.

#### Scenario: Reduced model not nested
- **WHEN** a reduced-model column lies outside the full model's column space
- **THEN** the scan raises an error stating that the models are not nested

#### Scenario: Non-finite outcome
- **WHEN** the outcome contains a NaN
- **THEN** the scan raises an error naming the outcome before any computation
