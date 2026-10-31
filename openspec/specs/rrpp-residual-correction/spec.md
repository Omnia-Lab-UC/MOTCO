# rrpp-residual-correction Specification

## Purpose
Defines the residual-corrected permutation in the RRPP primitive: an opt-in rescaling of the reduced model's residuals before they are permuted, which restores the residual variance the Freedman–Lane null loses to the reduced model's degrees of freedom.

## Requirements

### Requirement: RRPP accepts an opt-in residual correction

RRPP SHALL accept a residual correction that is absent by default or names one of `df` and `leverage`. When it is absent, every returned null distribution SHALL be identical to RRPP without the setting. Any other value SHALL be rejected before any permutation is drawn.

#### Scenario: Absent correction is byte-identical

- **WHEN** RRPP runs on fixed inputs, permutation count and seed, once without the setting and once with it explicitly absent
- **THEN** both runs return identical `delta`, `angle` and `shape` null distributions, and they equal those of the previous release

#### Scenario: Unknown correction is rejected

- **WHEN** RRPP is called with the correction `hc3`
- **THEN** it raises an error naming the value and the allowed corrections, and draws no permutation

### Requirement: The df correction scales every residual uniformly

Under `df`, each reduced-model residual SHALL be multiplied by κ = √(n / (n − r)), where n is the number of samples and r the numerical rank of the reduced model matrix.

#### Scenario: Multiplier at the case study's model

- **WHEN** the reduced model has 80 rows and rank 10
- **THEN** every multiplier equals √(80 / 70)

#### Scenario: df correction equals the diagnosis's construction

- **WHEN** RRPP runs with `df` on a response matrix, and separately without a correction on the reduced fitted values plus κ times the reduced residuals of that matrix, with the same models and seed
- **THEN** the two runs return identical null distributions up to floating-point tolerance

#### Scenario: Rank-deficient reduced model uses its rank

- **WHEN** the reduced model carries a duplicated column, so its rank is one less than its column count
- **THEN** κ is computed from the rank, not from the column count

### Requirement: The leverage correction scales each residual by its own leverage

Under `leverage`, each reduced-model residual SHALL be divided by √(1 − h_ii), where h_ii is that sample's diagonal element of the reduced model's hat matrix.

#### Scenario: Multipliers match the hat matrix

- **WHEN** the leverage multipliers are computed for a reduced model with unequal cell sizes and a numeric covariate
- **THEN** each equals 1 / √(1 − h_ii), with h computed as X (XᵀX)⁻¹ Xᵀ for that model, and the h_ii sum to the model's rank

#### Scenario: Balanced design reduces to the df correction

- **WHEN** the reduced model is a balanced intercept-plus-factor design in which every sample has the same leverage
- **THEN** the leverage multipliers all equal the df multiplier κ

### Requirement: The correction moves only the permuted residuals

With either correction, the permutation indices SHALL be drawn exactly as without a correction, and the reduced fitted values and the full-model estimation SHALL be unchanged. The serial and parallel paths SHALL return identical null distributions for the same seed and worker count as without the correction.

#### Scenario: Same permutation indices

- **WHEN** RRPP runs without a correction, with `df` and with `leverage` on the same inputs and seed
- **THEN** each permutation draw uses the same sample order in all three runs

#### Scenario: Statistics do not depend on the residuals' mean

- **WHEN** the leverage-rescaled residuals, which need not sum to zero, are centred before permutation
- **THEN** the `delta`, `angle` and `shape` null distributions are unchanged up to floating-point tolerance

#### Scenario: Parallel path matches the corrected serial contract

- **WHEN** RRPP runs with a correction in parallel with a fixed seed and worker count, twice
- **THEN** both runs return identical null distributions, and their values are drawn from the corrected residuals

### Requirement: An undefined correction is an error

RRPP SHALL raise an error before any permutation when a requested correction is undefined. This covers `df` when n − r is not positive, and `leverage` when any sample's leverage is 1 to numerical tolerance. The error SHALL name the counts or the sample index.

#### Scenario: Saturated reduced model

- **WHEN** `df` is requested and the reduced model's rank equals the number of samples
- **THEN** RRPP raises an error stating n and the rank

#### Scenario: Sample with leverage one

- **WHEN** `leverage` is requested and one sample is the only member of a reduced-model category, so its leverage is 1
- **THEN** RRPP raises an error naming that sample's row index

### Requirement: The multipliers are available without running RRPP

The permutation module SHALL provide the per-sample multipliers of a correction for a given reduced model, so a caller can record them. The multipliers SHALL be those RRPP applies, and SHALL raise the same errors for an undefined correction.

#### Scenario: Helper agrees with RRPP

- **WHEN** a caller computes the `leverage` multipliers for a reduced model and, for each permutation order drawn from RRPP's seed, evaluates the statistics on the reduced fitted values plus the permuted multiplied reduced residuals
- **THEN** the values equal RRPP's `leverage` null distributions on the original response
