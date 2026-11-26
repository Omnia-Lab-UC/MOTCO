# numpy-omics-generator Specification

## Purpose
Defines MOTCO's numpy-native multi-omic generator: it reproduces InterSIM's mean-shift generative model from a cached reference, runs reproducibly by seed with no R at runtime, exposes per-stage differential indicators, and has its realism validated against InterSIM.
## Requirements
### Requirement: numpy generator reproduces InterSIM's mean-shift model
MOTCO SHALL provide a numpy-native multi-omic generator that produces methylation, gene expression, and proteomics matrices using InterSIM's generative model: each stage mean is `μ = base + δ · v` for a per-stage differential indicator `v ∈ {0,1}ᵖ`, with samples drawn from a multivariate normal with the reference covariance, and methylation passed through the inverse-logit (`rev.logit`) transform after the additive shift.

#### Scenario: Generates the three aligned omic matrices
- **WHEN** the generator is invoked with valid parameters
- **THEN** it returns methylation, expression, and proteomics matrices with the same number of rows, aligned to a shared sample order

#### Scenario: Methylation shift is applied in logit space
- **WHEN** a differential shift is applied to methylation features
- **THEN** the shift is added to the pre-logit (M-value) mean and the inverse-logit transform is applied last, so methylation values remain in `(0, 1)`

#### Scenario: Cross-omic coupling is preserved
- **WHEN** differential expression and proteomics features are derived by default from the differential methylation features
- **THEN** the generator maps differential CpGs to genes and genes to proteins using the cached cross-omic maps

### Requirement: Generation uses cached InterSIM reference data with provenance
The generator SHALL source its reference means, covariances, and cross-omic maps from a cached artifact committed to the repository, recorded with provenance, and reproducible by a one-time export.

#### Scenario: Reference cache is loaded without R
- **WHEN** the generator runs
- **THEN** it loads the reference means, covariances, and maps from the cached artifact without invoking R

#### Scenario: Provenance is recorded
- **WHEN** the reference cache is produced
- **THEN** it records the InterSIM version, export date, and the export procedure used to create it

### Requirement: Generation is reproducible by seed and requires no R at runtime
The generator SHALL be deterministic given its seed and parameters, and SHALL NOT require an R runtime.

#### Scenario: Same seed produces identical matrices
- **WHEN** the generator is invoked twice with the same seed and parameters
- **THEN** the returned matrices are identical

#### Scenario: Runs without Rscript on PATH
- **WHEN** the generator is invoked in an environment with no `Rscript`
- **THEN** generation succeeds using only the cached reference data

### Requirement: Generator realism is validated against InterSIM
The generator's output distributions SHALL be validated to match InterSIM's across a swept, replicate-based protocol — not a single degenerate point — as a guard against reimplementation drift and as paper-supportable evidence of fidelity. The validation SHALL run InterSIM and the numpy generator each multiple times per parameter cell over a grid of `delta` and `p.DMP`, and SHALL compare, per cell: per-omic marginal moments/quantiles, cluster separation (η²), differential-feature rates (the DMP→DEG→DEP coupling), covariance structure, and cross-omic coupling. Agreement SHALL be judged against InterSIM's own sampling distribution (the numpy statistic falling within InterSIM's documented central interval), so the criterion accounts for InterSIM's RNG variability. The InterSIM side SHALL be captured as committed fixtures so the validation runs without R, and SHALL be reproducible from a committed R script with recorded provenance.

#### Scenario: Fidelity holds across the parameter sweep
- **WHEN** the validation is run over the `delta` × `p.DMP` grid
- **THEN** for each cell and each compared statistic, the numpy generator's value falls within InterSIM's documented central interval for that statistic

#### Scenario: Effect injection and cross-omic coupling are validated at non-zero effect
- **WHEN** the validation is run at `delta > 0`
- **THEN** cluster separation (η²) and differential-feature rates (DMP→DEG→DEP) for the numpy generator agree with InterSIM's distribution, exercising the effect injection and the cross-omic coupling that a `delta = 0` check cannot

#### Scenario: Validation runs without R from committed fixtures
- **WHEN** the validation runs in CI with no `Rscript` available
- **THEN** it compares the numpy generator against the committed InterSIM summary fixtures and passes without invoking R

#### Scenario: Fixtures are reproducible with recorded provenance
- **WHEN** the InterSIM fixtures are regenerated
- **THEN** the committed R script reproduces them, and the fixtures record the InterSIM version, generation date, seeds, and the parameter grid

#### Scenario: A reproducible supplementary artifact is produced
- **WHEN** the supplementary artifact generator is run against the committed fixtures
- **THEN** it produces a paper-ready table and figure summarizing numpy-vs-InterSIM fidelity across the grid

#### Scenario: A qualitative visual supplement is produced
- **WHEN** the visual-supplement generator is run against regenerated InterSIM raw data (produced with InterSIM via the dev flake) with the numpy side generated live
- **THEN** it renders side-by-side InterSIM-vs-numpy figures — per-omic marginal densities, per-modality clustermap heatmaps (with sample/feature dendrograms and a cluster colour bar), per-modality PCA, per-feature mean/variance agreement scatter, and a cross-omic coupling correlation block
- **AND** the InterSIM raw matrices are not committed to the repository; the rendering code is exercised R-free in CI via a synthetic stand-in fixture

### Requirement: Generator exposes per-stage differential indicators
The generator SHALL expose, as output, the differential-feature indicator used for each stage (and group) so downstream consumers can use the ground truth.

#### Scenario: Differential indicators are returned per stage
- **WHEN** generation succeeds
- **THEN** the per-stage differential-feature indicators for each omic layer are available in the returned truth structure

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
