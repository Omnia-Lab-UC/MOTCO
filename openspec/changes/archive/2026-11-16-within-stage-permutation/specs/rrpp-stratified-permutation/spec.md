# Spec Delta

## Purpose

Defines the stratified permutation in the RRPP primitive: an opt-in restriction of the residual permutation to exchangeability strata, so that the reduced model's residuals are exchanged only among samples of the same stratum.

## ADDED Requirements

### Requirement: RRPP accepts opt-in permutation strata

RRPP SHALL accept permutation strata that are absent by default or give one label per sample. When they are absent, every returned null distribution SHALL be identical to RRPP without the setting. Strata whose length differs from the number of samples, or that hold a missing label, SHALL be rejected before any permutation is drawn.

#### Scenario: Absent strata are byte-identical

- **WHEN** RRPP runs on fixed inputs, permutation count and seed, once without the setting and once with it explicitly absent
- **THEN** both runs return identical `delta`, `angle` and `shape` null distributions, and they equal those of the previous release

#### Scenario: Wrong length is rejected

- **WHEN** RRPP is called on 80 samples with 79 strata labels
- **THEN** it raises an error naming both counts, and draws no permutation

#### Scenario: Missing label is rejected

- **WHEN** one strata label is missing
- **THEN** RRPP raises an error naming the row, and draws no permutation

### Requirement: Each draw exchanges residuals only within a stratum

With strata, every permutation draw SHALL replace each sample's reduced-model residual by the residual of a sample of the same stratum, and SHALL use every residual exactly once.

#### Scenario: Draws stay inside the strata

- **WHEN** RRPP draws many permutations under strata of unequal sizes
- **THEN** in every draw each sample receives the residual of a sample with its own strata label, and the draw is a permutation of the samples

#### Scenario: A one-sample stratum is held fixed

- **WHEN** one sample is the only member of its stratum
- **THEN** it keeps its own residual in every draw

#### Scenario: Stratum means stay at the reduced fit

- **WHEN** the reduced model contains the strata indicators, and no residual correction is requested
- **THEN** in every draw the permuted response has the same mean within each stratum as the reduced model's fitted values

### Requirement: The stratified draw sequence is defined

With strata, each draw SHALL visit the strata in ascending order of their labels and draw one permutation of each stratum's samples, all from the one generator the seed defines. The sequence SHALL depend on the labels and not on the order in which the strata first appear among the samples. With every sample in one stratum the draws SHALL be those of RRPP without strata.

#### Scenario: One stratum equals the unrestricted draws

- **WHEN** RRPP runs with every sample in one stratum, and separately without strata, on the same inputs and seed
- **THEN** the two runs return identical null distributions

#### Scenario: Strata equal the diagnosis's construction

- **WHEN** RRPP runs with the stage of each sample as its stratum, coded by the stage's position in the declared order, and the predictive size-excess diagnosis's within-stage loop runs on the same response, models, residual correction and seed
- **THEN** the two return identical null distributions

#### Scenario: Label order, not row order

- **WHEN** the samples are reordered together with their strata labels, so another stratum appears first
- **THEN** each stratum is still visited in ascending label order

### Requirement: Strata move only the permutation draws

With strata, the reduced fitted values, the residual correction's multipliers and the full-model estimation SHALL be those of RRPP without strata. The strata SHALL combine with either residual correction, with either orientation estimator and with the per-draw eigengap return, none of which changes the draws.

#### Scenario: Residual correction under strata

- **WHEN** RRPP runs with strata and `df` on a response matrix, and separately with the same strata and no correction on the reduced fitted values plus κ times the reduced residuals of that matrix, with the same models and seed
- **THEN** the two runs return identical null distributions up to floating-point tolerance

#### Scenario: Eigengap return does not change the draws

- **WHEN** RRPP runs with strata, once with and once without the per-draw eigengap return
- **THEN** the `delta`, `angle` and `shape` null distributions are identical

### Requirement: The parallel path honours the strata

With strata and more than one worker, every draw of every worker SHALL exchange residuals only within a stratum, and two runs with the same seed and worker count SHALL return identical null distributions.

#### Scenario: Parallel stratified run is reproducible

- **WHEN** RRPP runs with strata in parallel with a fixed seed and worker count, twice
- **THEN** both runs return identical null distributions

#### Scenario: Parallel draws stay inside the strata

- **WHEN** RRPP runs with strata in parallel on a response whose residuals are constant within each stratum
- **THEN** every draw returns the statistics of the unpermuted response
