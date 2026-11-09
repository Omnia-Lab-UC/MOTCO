# Spec Delta

## ADDED Requirements

### Requirement: Real-data analysis offers the per-feature interaction scan

The analysis parameters SHALL accept an opt-in per-feature interaction scan, off by default. When on, the analysis SHALL run it on its standardized features with the full and reduced models including any covariates, and with the analysis's seed, permutation count and residual correction. Every other output SHALL be unchanged. The scan SHALL be available with latent adjustment and with attribution.

#### Scenario: Other outputs are unchanged
- **WHEN** the same analysis runs with and without the scan on the same inputs, permutation count and seed
- **THEN** the latent scores, selected rank, statistics, p-values, null summary, stage means, any split-test block and any attribution output are identical

#### Scenario: Analysis uses the primitive
- **WHEN** an analysis with six covariates runs with the scan, 99 permutations and `df`
- **THEN** its scan results equal the per-feature interaction scan applied to the analysis's standardized features, the full and reduced models with covariates, the seed and `df`

#### Scenario: Latent adjustment does not change the scan
- **WHEN** the same analysis runs with the scan, once without latent adjustment and once with latent adjustment `stage`
- **THEN** the two scan blocks and per-feature tables are identical

#### Scenario: Draws follow the serial sequence at the seed
- **WHEN** the analysis runs with more than one job
- **THEN** the scan still draws the serial index sequence at the analysis's seed, and the record states that its draws match RRPP's only when one job is used

### Requirement: Real-data analysis records the per-feature interaction scan

When the scan is on, the record SHALL carry an `interaction_scan` block after `measurement` and after any `span_split` block, giving the feature and degenerate counts, df, the maximum F with its feature and block, and the residual correction. The resolved parameters SHALL carry the setting only when it is on.

#### Scenario: Observed quantities are recorded
- **WHEN** an analysis of two blocks of 5,000 features runs with the scan
- **THEN** the `interaction_scan` block gives 10,000 features, the degenerate count, df and the maximum F with its feature and block
- **AND** the resolved parameters carry `interaction_scan: true`

#### Scenario: No setting means no change
- **WHEN** the committed identity analyses run without the scan
- **THEN** their records, minus the runtime, and their tables equal the committed fixtures, and neither the resolved parameters nor the record carries a scan key

### Requirement: The scan record carries the permutation results

With permutations, the `interaction_scan` block SHALL add the global p-value, the features with adjusted p ≤ 0.05 (block, feature, F, adjusted p) ordered by F, per block the hit count and the share of null maxima, and whether its draws match RRPP's. Without permutations it SHALL carry none of these.

#### Scenario: Permutation results are recorded
- **WHEN** an analysis runs with the scan and 999 permutations
- **THEN** the block gives the global p-value, the hit list ordered by F, and per block the hit count and the share of null maxima, the shares summing to 1

#### Scenario: Without permutations only the observed quantities are recorded
- **WHEN** an analysis runs with the scan and 0 permutations
- **THEN** the `interaction_scan` block gives the counts, df and the maximum F with its feature, and no p-values, hit list or null-maximum shares

### Requirement: Real-data analysis tabulates the scan for every feature

When the scan is on, the analysis result SHALL carry a per-feature table in the joint feature order. Each row SHALL give the block, the original feature name, the interaction and residual sums of squares, the F, the adjusted p-value (with permutations) and a degenerate flag.

#### Scenario: One row per joint feature
- **WHEN** an analysis of two blocks of 5,000 features runs with the scan and permutations
- **THEN** the table has 10,000 rows in the joint feature order, and each non-degenerate row's F and adjusted p-value equal the primitive's for that column

#### Scenario: Degenerate feature
- **WHEN** one feature is constant across samples
- **THEN** its row is flagged degenerate with no F or adjusted p-value, and the scan's maximum and hit list ignore it

### Requirement: The scan table gives each feature's group difference at each stage

Each row of the scan table SHALL give the group difference at each stage, in declared stage order: the second declared group minus the first, from the full model's LS means at the covariate reference row, in standardized units.

#### Scenario: Group differences are the model's LS-mean differences
- **WHEN** an analysis with covariates runs with the scan
- **THEN** each feature's group difference at each stage equals the difference of the full model's fitted LS means for that stage at the record's covariate reference row, computed on that feature's standardized values

### Requirement: motco analyze exposes the per-feature interaction scan

`motco analyze` SHALL accept `--interaction-scan`, off by default, which sets the scan. When it is on, the printed summary SHALL add one line with the maximum F and its feature, the global p-value and the number of features with adjusted p ≤ 0.05, and the output directory SHALL hold `interaction_scan.csv`.

#### Scenario: Flag runs the scan
- **WHEN** `motco analyze` runs with `--interaction-scan` and `--permutations 99`
- **THEN** `result.json` carries the `interaction_scan` block, the output directory holds `interaction_scan.csv`, and the printed summary shows the maximum F, its feature, the global p-value and the hit count

#### Scenario: Flag combined with latent adjustment
- **WHEN** `motco analyze` runs with `--interaction-scan` and `--latent-adjustment stage` on a dataset with covariates
- **THEN** the analysis runs and the scan line is printed

## MODIFIED Requirements

### Requirement: Real-data analysis writes a reproducible analysis artifact

The analysis SHALL be writable to a directory as: `result.json` (the JSON-safe record and resolved parameters), `latent_scores.csv` (samples × components, indexed by sample identifier), `pls_loadings.csv` (joint features × components, with block name and original feature name columns), `stage_means.csv` (group × stage latent LS-mean vectors, in declared stage order), when attribution ran, its tables, and, when the interaction scan ran, `interaction_scan.csv` (its per-feature table). The writer SHALL create the directory if absent and SHALL refuse to overwrite an existing `result.json` unless overwriting is requested.

#### Scenario: Outputs are complete and loadable

- **WHEN** a result with permutations is written to an empty directory
- **THEN** the directory contains `result.json`, `latent_scores.csv`, `pls_loadings.csv` and `stage_means.csv`, and `result.json` parses as JSON holding the p-values and the selected rank

#### Scenario: Existing result is protected

- **WHEN** a result is written to a directory that already holds `result.json` without requesting overwrite
- **THEN** the writer raises an error and leaves the existing files unchanged

#### Scenario: Scan table is written only when the scan ran

- **WHEN** a result with the scan is written, and a result without it is written to another directory
- **THEN** the first directory also contains `interaction_scan.csv` with one row per joint feature, and the second contains no such file
