# Spec Delta

## ADDED Requirements

### Requirement: Real-data analysis offers the residual-corrected permutation

The analysis parameters SHALL accept a residual correction, absent by default or one of `df` and `leverage`, and pass it to RRPP with the reduced model including any covariates. The observed statistics, latent space, models and LS-means SHALL be unchanged by it. An unknown value SHALL be rejected before any fitting.

#### Scenario: Correction moves only the p-values

- **WHEN** the same covariate-adjusted analysis runs with no correction and with `df`, on the same inputs, permutation count and seed
- **THEN** the latent scores, selected rank, statistics and stage means are identical, and only the p-values and null summary may differ

#### Scenario: Analysis uses the primitive's corrected null

- **WHEN** an analysis with six covariates runs with `leverage`
- **THEN** its p-values equal those computed from RRPP with `leverage` on the analysis's latent scores, full and reduced models with covariates, LS-means, contrast and seed

#### Scenario: Unknown correction is rejected before fitting

- **WHEN** the parameters name the correction `uniform`
- **THEN** the analysis raises an error naming the value and the allowed corrections before fitting the latent space

#### Scenario: Correction without covariates is allowed

- **WHEN** an analysis without covariates runs with `df`
- **THEN** it completes, and κ is computed from the reduced model of intercept, group and stage

### Requirement: Real-data analysis records the residual correction

When a correction is set and permutations run, the measurement block SHALL state the correction's kind, the number of samples and the reduced model's rank. For `df` it SHALL also state κ, and for `leverage` the minimum and maximum leverage and multiplier. The resolved parameters SHALL carry the correction only when it is set. Without a correction, every output SHALL equal 0.12.0's except the runtime.

#### Scenario: df correction is recorded

- **WHEN** an analysis with 80 samples and a reduced model of rank 10 runs with `df` and 99 permutations
- **THEN** the record's measurement block states kind `df`, 80 samples, rank 10 and κ = √(80 / 70), and the resolved parameters carry `residual_correction: "df"`

#### Scenario: No correction means no change

- **WHEN** the committed identity analyses run without a correction
- **THEN** their records, minus the runtime, and their tables equal the committed fixtures, and neither the resolved parameters nor the measurement block carries a residual-correction key

#### Scenario: Correction without permutations is not reported as applied

- **WHEN** an analysis runs with `df` and 0 permutations
- **THEN** the resolved parameters carry the setting and the measurement block carries no residual-correction entry

## MODIFIED Requirements

### Requirement: motco analyze subcommand

`motco analyze` SHALL be a CLI subcommand over the real-data analysis. It SHALL accept: a repeatable `--block NAME=PATH` (CSV, first column the sample identifier); a repeatable `--transform NAME=KIND` (default `identity` for any block not named); `--metadata PATH`, `--sample-col`, `--group-col` and `--stage-col`; a required `--stage-order` (comma-separated labels) and an optional `--group-order`; a repeatable `--covariate COL`; a repeatable `--categorical-covariate COL`, which also declares the column a covariate; `--missing-covariates` (`error` default, `complete_case` selectable); `--latent-adjustment` (`stage` or `group_stage`; absent by default, and rejected without a covariate); `--permutations` (default 0), `--seed`, `--n-jobs` and `--forced-components`; `--residual-correction` (`df` or `leverage`; absent by default); `--min-components` (default 2; `1` restores the unfloored rule); `--orientation-estimator` (`net_displacement` default, `pc1` selectable); `--attribution` with `--bootstrap-replicates`, `--bootstrap-seed` and `--top-k`; `--out-dir` (required); and `--overwrite`. Covariates SHALL be declared in the order the flags are given. Cross-validation settings SHALL be exposed with the harness defaults. Input and validation errors SHALL exit non-zero with a descriptive message. A successful run SHALL write the analysis artifact and print the pair statistics, p-values, selected rank, rank floor and estimator. When covariates are declared, it SHALL also print the covariate names and the number of dropped samples. When latent adjustment is requested, it SHALL also print the conditioning. When a residual correction is applied, it SHALL also print the correction and its κ or leverage range.

#### Scenario: Two-block analysis from CSV files

- **WHEN** `motco analyze --block rna=rna.csv --block atac=atac.csv --metadata meta.csv --sample-col donor --group-col Sex --stage-col ADNC --stage-order "Low,Intermediate,High" --permutations 99 --seed 1 --out-dir out/` runs on aligned files
- **THEN** the command exits zero, `out/result.json` exists with p-values for `delta`, `angle` and `shape`, and the printed summary names the selected rank, the floor of 2 and `net_displacement`

#### Scenario: Unknown block in a transform is rejected

- **WHEN** `--transform methyl=beta_to_mvalue` names a block that was not declared with `--block`
- **THEN** the command exits non-zero and the message names the undeclared block

#### Scenario: Missing stage order is rejected

- **WHEN** `--stage-order` is omitted
- **THEN** the command exits non-zero with a message stating the stage order is required

#### Scenario: Covariate-adjusted analysis from the CLI

- **WHEN** the two-block command above is run with `--covariate age_at_death --covariate pmi --categorical-covariate education_level --missing-covariates complete_case`
- **THEN** the command exits zero, `out/result.json` holds the covariate block naming the three covariates in that order with `education_level` categorical, and the printed summary names the covariates and the dropped-sample count

#### Scenario: Missing covariate column is rejected

- **WHEN** `--covariate rin` names a column absent from the metadata
- **THEN** the command exits non-zero and the message names `rin`

#### Scenario: Latent-adjusted analysis from the CLI

- **WHEN** the covariate-adjusted command above is run with `--latent-adjustment stage`
- **THEN** the command exits zero, `out/result.json` holds a latent-adjustment block with conditioning `stage`, and the printed summary names the conditioning

#### Scenario: Latent adjustment without covariates is rejected

- **WHEN** `--latent-adjustment stage` is given and no `--covariate` or `--categorical-covariate`
- **THEN** the command exits non-zero with a message stating latent adjustment requires declared covariates

#### Scenario: Unfloored rule from the CLI

- **WHEN** the two-block command above is run with `--min-components 1`
- **THEN** the command exits zero and `out/result.json` equals the artifact the 0.9.0 command wrote for the same inputs

#### Scenario: Invalid floor is rejected

- **WHEN** `--min-components 0` is given
- **THEN** the command exits non-zero with a message stating the floor must be at least 1

#### Scenario: Residual-corrected analysis from the CLI

- **WHEN** the covariate-adjusted command above is run with `--residual-correction df`
- **THEN** the command exits zero, `out/result.json` states the `df` correction with its κ in the measurement block, and the printed summary names the correction and κ

#### Scenario: Unknown residual correction is rejected

- **WHEN** `--residual-correction hc3` is given
- **THEN** the command exits non-zero and the message names the allowed corrections
