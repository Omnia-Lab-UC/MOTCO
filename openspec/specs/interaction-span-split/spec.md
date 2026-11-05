# interaction-span-split Specification

## Purpose
Tests the group × stage interaction of a multivariate outcome separately inside a given linear span and in its orthogonal complement, so a difference the stage-supervised measurement cannot see can still be tested, and located relative to that span.

## Requirements

### Requirement: The test measures interaction and residual energy in three parts

Given an outcome matrix, a full and a nested reduced model matrix and a basis with orthonormal columns, the test SHALL compute three parts: `total` (all coordinates), `span` (the projection onto the basis) and `complement` (the projection onto its orthogonal complement). For each part it SHALL report the dimension, the interaction sum of squares, the full model's residual sum of squares and the pseudo-F.

#### Scenario: Definitions match direct computation
- **WHEN** the test runs on an outcome Y, full model F, reduced model R and basis Q
- **THEN** each part's interaction sum of squares equals ‖(H_F − H_R)·Y_part‖², its residual sum of squares equals ‖(I − H_F)·Y_part‖², and its pseudo-F equals (interaction SS / df_int) / (residual SS / df_res) to a relative tolerance of 1e-10
- **AND** H are the models' hat matrices, Y_part is Y, Y·Q·Qᵀ or Y·(I − Q·Qᵀ), df_int = rank(F) − rank(R) and df_res = n − rank(F)

#### Scenario: Dimensions are reported
- **WHEN** the outcome has p columns and the basis r columns
- **THEN** the parts report dimensions p, r and p − r

#### Scenario: The basis matters only through its span
- **WHEN** the basis is replaced by Q·O for an orthogonal r × r matrix O
- **THEN** every reported quantity is unchanged to a relative tolerance of 1e-10

### Requirement: The parts decompose the total exactly

The `span` and `complement` interaction sums of squares SHALL add to the `total` one, and likewise the residual sums of squares. The test SHALL report the observed in-span share, `span` interaction SS over `total` interaction SS. The share SHALL be absent when the total interaction SS is zero.

#### Scenario: Exact decomposition
- **WHEN** the test runs on any valid input
- **THEN** span + complement equals total for both sums of squares to a relative tolerance of 1e-10, and the in-span share lies in [0, 1]

#### Scenario: Interaction energy inside the span
- **WHEN** the outcome's group × stage interaction lies entirely in the span of the basis and carries no noise
- **THEN** the in-span share is 1 and the complement's interaction SS is 0 to numerical tolerance

### Requirement: Permutation follows RRPP's residual randomization

With a positive permutation count, the test SHALL permute the reduced model's residuals by Freedman–Lane, drawing the same index sequence as RRPP's serial path at the same seed. Every draw SHALL be evaluated for all three parts. Each part SHALL report the pseudo-F p-value and the interaction-SS p-value, (1 + #{null ≥ observed}) / (1 + permutations). With zero permutations, only the observed quantities SHALL be reported.

#### Scenario: Draws match RRPP
- **WHEN** the test and RRPP run on the same outcome and models with the same seed and permutation count
- **THEN** every permuted outcome of the test equals the reduced fitted values plus the reduced residuals in the row order RRPP's serial path uses at that draw

#### Scenario: Parts share each draw
- **WHEN** null distributions are returned
- **THEN** the three parts' null values at draw k come from the same permuted outcome, and span + complement equals total for the null interaction SS at every draw

#### Scenario: No permutations
- **WHEN** the permutation count is 0
- **THEN** the result carries the observed quantities and no p-values

### Requirement: Residual corrections rescale the permuted residuals

The test SHALL accept RRPP's residual corrections (`df`, `leverage`, or none) and apply them to the reduced residuals before any draw. The observed quantities SHALL be unchanged by the correction. An unknown or undefined correction SHALL be rejected before any draw.

#### Scenario: df leaves the pseudo-F null unchanged
- **WHEN** the same test runs with no correction and with `df` at the same seed
- **THEN** the pseudo-F p-values are identical, and every null interaction SS under `df` equals the uncorrected one times κ² = n / (n − rank(R))

#### Scenario: Leverage correction
- **WHEN** the test runs with `leverage`
- **THEN** each reduced residual row is divided by √(1 − h_ii) of the reduced model before permutation, and the observed quantities equal the uncorrected run's

#### Scenario: Unknown correction
- **WHEN** the correction `uniform` is given
- **THEN** the test raises an error naming the value and the allowed corrections before any draw

### Requirement: Invalid inputs are rejected before any draw

The test SHALL raise an error naming the problem when the outcome, the models and the basis disagree in rows or columns, when any of them is non-finite, when the basis's columns are not orthonormal (tolerance 1e-8), when the complement is empty, when the reduced model does not lie in the full model's column space, or when df_int or df_res is not positive.

#### Scenario: Non-orthonormal basis
- **WHEN** the basis has columns whose Gram matrix differs from the identity by more than 1e-8
- **THEN** the test raises an error naming the largest deviation

#### Scenario: Empty complement
- **WHEN** the basis spans all p coordinates
- **THEN** the test raises an error stating that the complement is empty

#### Scenario: Reduced model not nested
- **WHEN** a reduced-model column lies outside the full model's column space
- **THEN** the test raises an error stating that the models are not nested
