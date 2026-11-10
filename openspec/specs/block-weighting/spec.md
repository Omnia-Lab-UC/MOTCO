# block-weighting Specification

## Purpose
Defines the opt-in, label-blind weighting of the pooled-standardized omic blocks that every MOTCO latent space may be fitted on, so that a block's weight in the joint measurement is a declared choice rather than its feature count.

## Requirements

### Requirement: Block weighting multiplies each standardized block by one label-blind scalar

When block weighting is requested, every feature of a block SHALL be multiplied, after pooled per-feature standardization, by one scalar fitted on the pooled sample without any group or stage label. Two kinds SHALL be offered:
- `equal_variance`: the inverse square root of the block's feature count, so every block has unit total variance;
- `mfa`: the inverse standard deviation of the standardized block's leading principal component, so every block's dominant direction has unit variance.

An unknown kind SHALL be rejected, naming the value and the allowed kinds, before any block is transformed.

#### Scenario: Equal variance
- **WHEN** blocks of 367 and 131 standardized features are weighted with `equal_variance`
- **THEN** the sum of the feature variances of each transformed block is 1 (within floating-point tolerance, counting constant features at zero variance against the block's count)

#### Scenario: MFA
- **WHEN** blocks are weighted with `mfa`
- **THEN** the leading principal component of each transformed block has unit variance

#### Scenario: Unknown kind
- **WHEN** a weighting other than `equal_variance` or `mfa` is requested
- **THEN** the request is rejected naming the value and the allowed kinds, and no block is transformed

### Requirement: The weights apply alike to observed rows and population means

The fitted weight of a block SHALL be applied to every transformation through the fitted preprocessing: observed rows, held-out rows and analytic population means.

#### Scenario: Population means
- **WHEN** a weighted preprocessor transforms a dataset's population stage means
- **THEN** each block of the result equals the unweighted transformation times the block's weight

### Requirement: Weighting is invisible while absent

Without a requested weighting, every weight SHALL be 1 and every output — transformed blocks, latent spaces, statistics, records, cell identifiers and signatures — SHALL be byte-identical to the output before weighting existed. No metadata entry SHALL mention weighting.

#### Scenario: Absent weighting
- **WHEN** a dataset is preprocessed and evaluated without a weighting
- **THEN** every output equals the pre-change output byte for byte, and no record or metadata carries a weighting key

### Requirement: Weighted PLS spaces are fitted in the weighted metric

When weighting is set, every PLS fit of the production latent space — cross-validation, per-repeat refits and the final model — SHALL be fitted without per-feature rescaling, so the weights determine the metric the space is fitted in.

#### Scenario: Weights survive the fit
- **WHEN** a PLS latent space is fitted on weighted blocks
- **THEN** its input scaling is unity and the fitted weights of a block are those of the weighted features, not of features rescaled to unit variance

#### Scenario: Unit weights reproduce production at a fixed rank
- **WHEN** a PLS latent space is fitted at a fixed rank in the weighted metric with every weight equal to 1
- **THEN** its weight vectors equal those of the production fit on the same rows and its scores equal production's divided by the uniform factor √((n − 1)/n)

#### Scenario: Unit weights select production's rank on a pinned dataset
- **WHEN** the cross-validation runs in the weighted metric with every weight equal to 1 on the pinned dataset
- **THEN** the selected rank and the cross-validation table equal those of the production fit on the same rows (the cross-validation differs from production only by its per-fold rescaling, which need not vanish on other datasets)
