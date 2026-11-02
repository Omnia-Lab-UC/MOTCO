# Spec Delta

## REMOVED Requirements

### Requirement: Real-data analysis can run frozen-model orientation attribution

**Reason**: The requirement rejected attribution whenever covariates were declared, through its scenario "Attribution with covariates is rejected". Covariate-adjusted attribution replaces that rejection, and a MODIFIED block cannot retire a scenario. Its other content is carried, unchanged, into "Real-data analysis runs frozen-model orientation attribution". The one rejection that remains is "Attribution is rejected with latent adjustment", and the covariate behaviour is "Real-data attribution adjusts for declared covariates".

**Migration**: Nothing to change for analyses without covariates. An analysis that declares covariates and requests attribution, which used to raise an error, now returns covariate-adjusted attribution. An analysis that requests attribution together with latent adjustment still raises an error before fitting.

## MODIFIED Requirements

### Requirement: motco analyze subcommand

`motco analyze` SHALL be a CLI subcommand over the real-data analysis. It SHALL accept: a repeatable `--block NAME=PATH` (CSV, first column the sample identifier); a repeatable `--transform NAME=KIND` (default `identity` for any block not named); `--metadata PATH`, `--sample-col`, `--group-col` and `--stage-col`; a required `--stage-order` (comma-separated labels) and an optional `--group-order`; a repeatable `--covariate COL`; a repeatable `--categorical-covariate COL`, which also declares the column a covariate; `--missing-covariates` (`error` default, `complete_case` selectable); `--latent-adjustment` (`stage` or `group_stage`; absent by default, and rejected without a covariate); `--permutations` (default 0), `--seed`, `--n-jobs` and `--forced-components`; `--residual-correction` (`df` or `leverage`; absent by default); `--min-components` (default 2; `1` restores the unfloored rule); `--orientation-estimator` (`net_displacement` default, `pc1` selectable); `--attribution` with `--bootstrap-replicates`, `--bootstrap-seed` and `--top-k`; `--out-dir` (required); and `--overwrite`. Covariates SHALL be declared in the order the flags are given. Cross-validation settings SHALL be exposed with the harness defaults. Input and validation errors SHALL exit non-zero with a descriptive message. A successful run SHALL write the analysis artifact and print the pair statistics, p-values, selected rank, rank floor and estimator. When covariates are declared, it SHALL also print the covariate names and the number of dropped samples. When latent adjustment is requested, it SHALL also print the conditioning. When a residual correction is applied, it SHALL also print the correction and its κ or leverage range. When attribution ran with covariates, it SHALL also print that the attribution used covariate-adjusted means and the number of rank-deficient bootstrap replicates.

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

#### Scenario: Covariate-adjusted attribution from the CLI

- **WHEN** the covariate-adjusted command above is run with `--attribution --bootstrap-replicates 20 --bootstrap-seed 3`
- **THEN** the command exits zero, `out/attribution/configuration.csv` states mean source `covariate_adjusted`, `out/result.json` holds the attribution's adjustment block, and the printed summary names covariate-adjusted attribution and the rank-deficient replicate count

#### Scenario: Attribution with latent adjustment is rejected from the CLI

- **WHEN** the covariate-adjusted command above is run with `--latent-adjustment stage --attribution`
- **THEN** the command exits non-zero with a message stating attribution is not available with latent adjustment

## ADDED Requirements

### Requirement: Real-data analysis runs frozen-model orientation attribution

When requested, the analysis SHALL run orientation attribution on the fitted pooled PLS estimator and the standardized joint feature matrix it was fitted to. It SHALL use the caller's stage and group order, a caller-chosen bootstrap count and seed, and a top-k size, and SHALL never refit the latent space. Attribution SHALL be disabled by default. The compact record SHALL mark truth recovery unavailable, and the full tables SHALL be available for writing.

#### Scenario: Attribution reuses the fitted estimator

- **WHEN** attribution is requested
- **THEN** the attribution's model component count equals the analysis's selected rank, and the estimator is not refitted

#### Scenario: Truth recovery is unavailable on real data

- **WHEN** attribution is requested
- **THEN** the compact record marks truth recovery unavailable with a reason, and no driver set is invented

### Requirement: Attribution is rejected with latent adjustment

The analysis SHALL NOT combine attribution with latent adjustment. That basis is fitted on covariate-residualized features, and its `angle` and `shape` are not admissible at the case study's configuration. A request for both SHALL raise an error before any fitting.

#### Scenario: Attribution with latent adjustment is rejected

- **WHEN** attribution is requested together with a latent adjustment
- **THEN** the analysis raises an error stating attribution is not available with latent adjustment, and no latent space is fitted

### Requirement: Real-data attribution adjusts for declared covariates

With declared covariates and attribution requested, the analysis SHALL attribute from covariate-adjusted cell means. These SHALL be evaluated at the record's covariate reference row, from the trajectory models' covariate design columns, with the bootstrap that refits the adjustment. The compact attribution record SHALL then carry an adjustment block. Without covariates it SHALL carry none, and every attribution output SHALL be unchanged.

#### Scenario: Attribution means project to the stage means

- **WHEN** a covariate-adjusted analysis runs with attribution
- **THEN** the attribution's observed group × stage means, transformed by the fitted estimator, equal the analysis's latent stage means to numerical tolerance

#### Scenario: Adjustment block is recorded

- **WHEN** a covariate-adjusted analysis with attribution is written out
- **THEN** the record's attribution block holds an adjustment block with mean source `covariate_adjusted`, the covariate design columns in model order, a reference row equal to the covariate block's, and the number of rank-deficient bootstrap replicates

#### Scenario: Attribution leaves the measurement unchanged

- **WHEN** the same covariate-adjusted analysis runs with and without attribution, on the same inputs, permutation count and seed
- **THEN** the latent scores, selected rank, statistics, p-values, null summaries and stage means are identical

#### Scenario: No covariates means no change in attribution

- **WHEN** an analysis with attribution and without covariates runs on a fixed dataset and settings
- **THEN** its record's attribution block, minus runtimes, and its attribution tables equal those produced before covariate-adjusted attribution existed, and the attribution block has no adjustment key

### Requirement: Real-data analysis exposes its frozen measurement in memory

The result SHALL expose, in memory, the fitted estimator, the standardized joint feature matrix whose projection the latent scores are, and the per-feature original scales. The written artifact SHALL NOT contain them, no estimator SHALL be serialized, and the record SHALL be unchanged by their presence.

#### Scenario: Exposed measurement reproduces the latent scores

- **WHEN** an analysis completes
- **THEN** transforming the exposed feature matrix with the exposed estimator returns the latent scores to within 1e-10 of the largest absolute score (the default path's scores are the fit's training scores, which agree with the transform only to rounding)

#### Scenario: Exposed measurement is not written

- **WHEN** a result is written to an empty directory
- **THEN** the directory holds the same files as before the measurement was exposed, and none of them serializes the estimator
