# Spec Delta

## ADDED Requirements

### Requirement: Real-data analysis offers the scan by block

The analysis parameters SHALL accept an opt-in by-block scan, off by default, that runs the per-feature interaction scan with the declared omic blocks as its block labels. It SHALL require the scan and be rejected without it before any fitting, naming both settings. Every output other than the scan's p-values and hits SHALL be unchanged.

#### Scenario: Scan by block
- **WHEN** an analysis of two declared blocks runs with the scan by block and permutations
- **THEN** its scan results equal the per-feature interaction scan with the two blocks as labels on the analysis's standardized features, models, seed and residual correction

#### Scenario: By block without the scan
- **WHEN** the parameters request the scan by block without the scan
- **THEN** they are rejected before any fitting, naming both settings

#### Scenario: Other outputs are unchanged
- **WHEN** the same analysis runs with the scan, once jointly and once by block, on the same inputs and seed
- **THEN** the latent scores, selected rank, statistics, p-values, stage means, any split-test block, and the scan's per-feature F and group differences are identical

#### Scenario: Absent setting
- **WHEN** the committed identity analyses run with or without the joint scan and without the by-block setting
- **THEN** their records, minus the runtime, and their tables equal the pre-change outputs, and the resolved parameters carry no by-block key

### Requirement: Real-data analysis records the scan by block

With the by-block scan, the `interaction_scan` record block SHALL give the combination and the number of blocks tested, and per block the maximum F with its feature and, with permutations, the block's own global p-value. Its global p-value, hit list and the table's adjusted p-values SHALL be the by-block ones.

#### Scenario: Record contents
- **WHEN** an analysis of two blocks runs with the scan by block and 999 permutations
- **THEN** the block names the combination, gives 2 blocks tested, and each block entry gives its maximum F, that feature and its global p-value, beside its hit count and share of null maxima

#### Scenario: Table holds the by-block adjusted p-values
- **WHEN** the same analysis writes `interaction_scan.csv`
- **THEN** its columns are unchanged and each adjusted p-value equals the primitive's by-block value for that feature

### Requirement: motco analyze exposes the scan by block

`motco analyze` SHALL accept `--interaction-scan-by-block`, off by default, exit with a clean error when it is given without `--interaction-scan`, and name the combination and the number of blocks in the printed scan line when it is on.

#### Scenario: Flag runs the scan by block
- **WHEN** `motco analyze` runs with `--interaction-scan --interaction-scan-by-block --permutations 99`
- **THEN** `result.json`'s `interaction_scan` block names the by-block combination and the printed scan line names it and the number of blocks

#### Scenario: Flag without the scan
- **WHEN** `motco analyze` runs with `--interaction-scan-by-block` and without `--interaction-scan`
- **THEN** it exits non-zero with a message naming both flags, before any fitting
