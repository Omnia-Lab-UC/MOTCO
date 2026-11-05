# Spec Delta

## ADDED Requirements

### Requirement: Real-data analysis offers the split interaction test

The analysis parameters SHALL accept an opt-in split interaction test, off by default. When on, the analysis SHALL run it on its pooled-standardized features in the fitted PLS model's input coordinates, against an orthonormal basis of the model's rotations. It SHALL use the full and reduced models with any covariates, the analysis's seed, permutation count and residual correction, and leave every other output unchanged.

#### Scenario: Other outputs are unchanged
- **WHEN** the same analysis runs with and without the split test on the same inputs, permutation count and seed
- **THEN** the latent scores, selected rank, statistics, p-values, null summary and stage means are identical

#### Scenario: Analysis uses the primitive
- **WHEN** an analysis with six covariates runs with the split test, 99 permutations and `df`
- **THEN** its split results equal the split interaction test applied to the analysis's standardized features in the model's input coordinates, the orthonormalized rotations, the full and reduced models with covariates, the seed and `df`

#### Scenario: Draws follow the serial sequence at the seed
- **WHEN** the analysis runs with more than one job
- **THEN** the split test still draws the serial index sequence at the analysis's seed, and the record states that its draws match RRPP's only when one job is used

### Requirement: Real-data analysis records the split interaction test

When the split test is on, the record SHALL carry a `span_split` block after `measurement`. Per part (`total`, `span`, `complement`) it SHALL give the dimension, the interaction and residual sums of squares, the degrees of freedom, the pseudo-F and, with permutations, the pseudo-F and interaction-SS p-values. It SHALL also give the observed in-span share and the residual correction applied. The resolved parameters SHALL carry the setting only when it is on.

#### Scenario: Block is recorded
- **WHEN** an analysis of 10,000 features at rank 5 runs with the split test and 999 permutations
- **THEN** the `span_split` block gives dimensions 10,000, 5 and 9,995, the observed quantities and both p-values per part, and the in-span share, and the resolved parameters carry `span_split_test: true`

#### Scenario: No setting means no change
- **WHEN** the committed identity analyses run without the split test
- **THEN** their records, minus the runtime, and their tables equal the committed fixtures, and neither the resolved parameters nor the record carries a split-test key

#### Scenario: Without permutations only the observed quantities are recorded
- **WHEN** an analysis runs with the split test and 0 permutations
- **THEN** the `span_split` block gives the observed quantities and the in-span share and no p-values

### Requirement: The split interaction test is rejected with latent adjustment

An analysis that requests both the split test and covariate latent adjustment SHALL raise an error before fitting the latent space, stating that the split test is defined on the covariate-blind basis.

#### Scenario: Both requested
- **WHEN** the parameters set the split test and latent adjustment `stage` on a dataset with covariates
- **THEN** the analysis raises an error naming both settings before any fitting

### Requirement: motco analyze exposes the split interaction test

`motco analyze` SHALL accept `--span-split-test`, off by default, which sets the split test. When it is on, the printed summary SHALL add one line with the three parts' pseudo-F p-values and the observed in-span share. Combining it with `--latent-adjustment` SHALL exit with a clean error message.

#### Scenario: Flag runs the test
- **WHEN** `motco analyze` runs with `--span-split-test` and `--permutations 99`
- **THEN** `result.json` carries the `span_split` block and the printed summary shows the three p-values and the in-span share

#### Scenario: Flag combined with latent adjustment
- **WHEN** `motco analyze` runs with `--span-split-test` and `--latent-adjustment stage`
- **THEN** it exits with an error message naming both options, without a traceback
