# Spec Delta

## MODIFIED Requirements

### Requirement: Covariates enter the full and reduced trajectory models

When covariates are declared, the analysis SHALL append the covariate design columns, unchanged, to both the full model (intercept, group, stage, group × stage) and the reduced model (intercept, group, stage). Covariates SHALL enter additively, with no interaction with group or stage. The full model with covariates MUST have full column rank; otherwise the analysis SHALL raise an error naming the covariates, for example when a covariate is collinear with the group × stage design. Every group × stage LS-mean SHALL be evaluated at one reference covariate row: the mean of each covariate design column over the analysed samples, which for a categorical covariate is its observed level proportions. RRPP SHALL permute the residuals of the reduced model with covariates, so every permuted dataset keeps the fitted covariate effects. Unless latent adjustment is requested, the latent space SHALL be built exactly as without covariates, from the standardized features and the stage labels alone. With latent adjustment requested, the latent space SHALL be built as the covariate-residualized latent space capability specifies, and the trajectory models above SHALL be unchanged.

#### Scenario: Statistics do not depend on the reference value

- **WHEN** the same covariate-adjusted analysis is evaluated with its LS-means at the pooled covariate mean and at any other fixed covariate row
- **THEN** `delta`, `angle` and `shape` are equal up to floating-point tolerance, and only the latent stage means differ, by the same vector for every cell

#### Scenario: Adjustment removes a known covariate effect

- **WHEN** a latent matrix is built as a pure group × stage cell-mean structure plus a covariate times a fixed vector, with the covariate distribution differing across group × stage cells
- **THEN** the adjusted statistics equal those of the cell-mean structure alone, and the unadjusted statistics differ from them

#### Scenario: Collinear covariate is rejected

- **WHEN** a numeric covariate equals 1 for every sample of one group and 0 for the other
- **THEN** the analysis raises an error stating the model with covariates is not of full rank and naming the covariates

#### Scenario: Latent space is unchanged by covariates

- **WHEN** the same data are analysed with and without covariates, the same rank selection settings and no latent adjustment
- **THEN** the latent scores, PLS loadings and selected rank are identical

#### Scenario: Latent adjustment changes the basis only

- **WHEN** the same data and covariates are analysed with and without latent adjustment
- **THEN** the two analyses use the same full model, reduced model, LS-mean reference row and contrast, and differ in the latent scores, loadings and possibly the selected rank

#### Scenario: Adjusted permutation p-values are reproducible

- **WHEN** a covariate-adjusted analysis runs twice with the same inputs, permutation count and seed
- **THEN** both runs return identical p-values and null summaries

### Requirement: Real-data analysis records the covariate design

When covariates are declared, the record SHALL carry a covariate block stating: each covariate's name and kind; each categorical covariate's ordered levels and reference level; the design column names in model order; the reference row used for the LS-means; the missing-covariate policy; the number and identifiers of samples it dropped; and the fitted covariate coefficients per latent component. The design block SHALL state the number of analysed samples. The resolved analysis parameters SHALL carry the latent-adjustment setting only when it is set. When no covariates are declared, the record SHALL contain no covariate block, the schema version SHALL be unchanged, and the record and every written output SHALL be identical to an analysis without covariate support, except for the measured runtime. When latent adjustment is not requested, the record and every written output SHALL be identical to an analysis without latent-adjustment support, except for the measured runtime.

#### Scenario: Covariate design is in result.json

- **WHEN** an analysis with covariates `age` and `education_level` is written out
- **THEN** `result.json` holds the covariate block with both names and kinds, the design column names, the reference row and the missing policy

#### Scenario: No covariates means no change

- **WHEN** an analysis without covariates runs on a fixed dataset and settings
- **THEN** its record, minus the runtime, and its CSV outputs equal those produced before covariate support existed, and the record has no covariate key

#### Scenario: No latent adjustment means no change

- **WHEN** an analysis with covariates and without latent adjustment runs on a fixed dataset and settings
- **THEN** its record, minus the runtime, and its CSV outputs equal those produced before latent-adjustment support existed, and the resolved parameters carry no latent-adjustment key

### Requirement: motco analyze subcommand

`motco analyze` SHALL be a CLI subcommand over the real-data analysis. It SHALL accept: a repeatable `--block NAME=PATH` (CSV, first column the sample identifier); a repeatable `--transform NAME=KIND` (default `identity` for any block not named); `--metadata PATH`, `--sample-col`, `--group-col` and `--stage-col`; a required `--stage-order` (comma-separated labels) and an optional `--group-order`; a repeatable `--covariate COL`; a repeatable `--categorical-covariate COL`, which also declares the column a covariate; `--missing-covariates` (`error` default, `complete_case` selectable); `--latent-adjustment` (`stage` or `group_stage`; absent by default, and rejected without a covariate); `--permutations` (default 0), `--seed`, `--n-jobs` and `--forced-components`; `--orientation-estimator` (`net_displacement` default, `pc1` selectable); `--attribution` with `--bootstrap-replicates`, `--bootstrap-seed` and `--top-k`; `--out-dir` (required); and `--overwrite`. Covariates SHALL be declared in the order the flags are given. Cross-validation settings SHALL be exposed with the harness defaults. Input and validation errors SHALL exit non-zero with a descriptive message. A successful run SHALL write the analysis artifact and print the pair statistics, p-values, selected rank and estimator. When covariates are declared, it SHALL also print the covariate names and the number of dropped samples. When latent adjustment is requested, it SHALL also print the conditioning.

#### Scenario: Two-block analysis from CSV files

- **WHEN** `motco analyze --block rna=rna.csv --block atac=atac.csv --metadata meta.csv --sample-col donor --group-col Sex --stage-col ADNC --stage-order "Low,Intermediate,High" --permutations 99 --seed 1 --out-dir out/` runs on aligned files
- **THEN** the command exits zero, `out/result.json` exists with p-values for `delta`, `angle` and `shape`, and the printed summary names the selected rank and `net_displacement`

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
