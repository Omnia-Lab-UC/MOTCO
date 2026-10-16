# Spec Delta

## ADDED Requirements

### Requirement: Real-data analysis accepts declared covariates

The analysis SHALL accept an optional ordered list of covariate columns from the sample metadata. A covariate whose column has a numeric (non-boolean) dtype SHALL be numeric unless the caller declares it categorical. Every other covariate SHALL be categorical. A numeric covariate SHALL contribute one design column holding its values. A categorical covariate SHALL contribute one indicator column per observed level except the first, with levels ordered by sorted string representation and the first level as reference. A covariate column MUST be distinct from the sample, group and stage columns and MUST appear once. A column the caller declares categorical MUST also be declared a covariate. The analysis SHALL reject the following, each with an error naming the covariate: a constant numeric covariate, a categorical covariate with one observed level, and a non-numeric value in a numeric covariate. All of these checks SHALL apply to the samples that are analysed after the missing-covariate policy has run. When no covariates are declared, the analysis SHALL behave exactly as without this requirement.

#### Scenario: Numeric and categorical covariates build the design columns

- **WHEN** the covariates `age` (numeric) and `education_level` (levels `College`, `Graduate`, `HighSchool`) are declared
- **THEN** the covariate design has three columns: `age`, `education_level[Graduate]` and `education_level[HighSchool]`, with `College` as the reference level

#### Scenario: Numeric-coded factor is declared categorical

- **WHEN** a covariate column holds the integers 1, 2 and 3 and the caller declares it categorical
- **THEN** it contributes two indicator columns and is recorded as categorical

#### Scenario: Constant covariate is rejected

- **WHEN** a declared numeric covariate takes the same value in every analysed sample
- **THEN** the analysis raises an error naming the covariate as constant

#### Scenario: Covariate naming the group column is rejected

- **WHEN** the group column is also declared a covariate
- **THEN** the analysis raises an error stating covariates must differ from the sample, group and stage columns

### Requirement: Real-data analysis applies an explicit missing-covariate policy

The caller SHALL choose a missing-covariate policy, either `error` (the default) or `complete_case`. Under `error`, a missing value in any declared covariate SHALL raise an error naming the column and listing the affected sample identifiers. Under `complete_case`, the analysis SHALL drop every sample with a missing value in any declared covariate from every block and from the metadata before standardization, the latent fit and the design. Every later step SHALL run on the retained samples only. After deletion, every group × stage cell MUST still be observed; otherwise the analysis SHALL raise an error naming the emptied cells. Covariate values SHALL never be imputed. The finite-values requirement on blocks is unchanged. The policy SHALL have no effect when no covariates are declared.

#### Scenario: Missing covariate is an error by default

- **WHEN** one sample has no `pmi` value and the policy is not given
- **THEN** the analysis raises an error naming `pmi` and that sample's identifier

#### Scenario: Complete-case deletion drops the sample everywhere

- **WHEN** two of 80 samples miss `rin` and the policy is `complete_case`
- **THEN** the latent scores have 78 rows, the recorded cell sizes sum to 78, and the record lists the two dropped sample identifiers

#### Scenario: Deletion that empties a cell is rejected

- **WHEN** complete-case deletion removes every sample of one group at one stage
- **THEN** the analysis raises an error naming that group × stage combination

### Requirement: Covariates enter the full and reduced trajectory models

When covariates are declared, the analysis SHALL append the covariate design columns, unchanged, to both the full model (intercept, group, stage, group × stage) and the reduced model (intercept, group, stage). Covariates SHALL enter additively, with no interaction with group or stage. The full model with covariates MUST have full column rank; otherwise the analysis SHALL raise an error naming the covariates, for example when a covariate is collinear with the group × stage design. Every group × stage LS-mean SHALL be evaluated at one reference covariate row: the mean of each covariate design column over the analysed samples, which for a categorical covariate is its observed level proportions. RRPP SHALL permute the residuals of the reduced model with covariates, so every permuted dataset keeps the fitted covariate effects. The latent space SHALL be built exactly as without covariates, from the standardized features and the stage labels alone.

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

- **WHEN** the same data are analysed with and without covariates and the same rank selection settings
- **THEN** the latent scores, PLS loadings and selected rank are identical

#### Scenario: Adjusted permutation p-values are reproducible

- **WHEN** a covariate-adjusted analysis runs twice with the same inputs, permutation count and seed
- **THEN** both runs return identical p-values and null summaries

### Requirement: Real-data analysis records the covariate design

When covariates are declared, the record SHALL carry a covariate block stating: each covariate's name and kind; each categorical covariate's ordered levels and reference level; the design column names in model order; the reference row used for the LS-means; the missing-covariate policy; the number and identifiers of samples it dropped; and the fitted covariate coefficients per latent component. The design block SHALL state the number of analysed samples. When no covariates are declared, the record SHALL contain no covariate block, the schema version SHALL be unchanged, and the record and every written output SHALL be identical to an analysis without covariate support, except for the measured runtime.

#### Scenario: Covariate design is in result.json

- **WHEN** an analysis with covariates `age` and `education_level` is written out
- **THEN** `result.json` holds the covariate block with both names and kinds, the design column names, the reference row and the missing policy

#### Scenario: No covariates means no change

- **WHEN** an analysis without covariates runs on a fixed dataset and settings
- **THEN** its record, minus the runtime, and its CSV outputs equal those produced before covariate support existed, and the record has no covariate key

## MODIFIED Requirements

### Requirement: Real-data analysis can run frozen-model orientation attribution

When requested, the analysis SHALL run orientation attribution on the fitted pooled PLS estimator and the standardized joint feature matrix it was fitted to, with the caller's stage order and group labels, a caller-chosen bootstrap replicate count and seed, and a top-k size; it SHALL never refit the latent space to do so. The compact attribution record SHALL report the truth-recovery block as unavailable, because real data carry no generator truth, and the full attribution tables SHALL be available for writing. Attribution SHALL be disabled by default. Attribution SHALL NOT be combined with declared covariates: it works from unadjusted cell means, so a request for both SHALL raise an error before any fitting.

#### Scenario: Attribution reuses the fitted estimator

- **WHEN** attribution is requested
- **THEN** the attribution's model component count equals the analysis's selected rank, and the estimator is not refitted

#### Scenario: Truth recovery is unavailable on real data

- **WHEN** attribution is requested
- **THEN** the compact record marks truth recovery unavailable with a reason, and no driver set is invented

#### Scenario: Attribution with covariates is rejected

- **WHEN** attribution is requested and at least one covariate is declared
- **THEN** the analysis raises an error stating attribution is not available with covariates, and no latent space is fitted

### Requirement: motco analyze subcommand

`motco analyze` SHALL be a CLI subcommand over the real-data analysis. It SHALL accept: a repeatable `--block NAME=PATH` (CSV, first column the sample identifier); a repeatable `--transform NAME=KIND` (default `identity` for any block not named); `--metadata PATH`, `--sample-col`, `--group-col` and `--stage-col`; a required `--stage-order` (comma-separated labels) and an optional `--group-order`; a repeatable `--covariate COL`; a repeatable `--categorical-covariate COL`, which also declares the column a covariate; `--missing-covariates` (`error` default, `complete_case` selectable); `--permutations` (default 0), `--seed`, `--n-jobs` and `--forced-components`; `--orientation-estimator` (`net_displacement` default, `pc1` selectable); `--attribution` with `--bootstrap-replicates`, `--bootstrap-seed` and `--top-k`; `--out-dir` (required); and `--overwrite`. Covariates SHALL be declared in the order the flags are given. Cross-validation settings SHALL be exposed with the harness defaults. Input and validation errors SHALL exit non-zero with a descriptive message. A successful run SHALL write the analysis artifact and print the pair statistics, p-values, selected rank and estimator. When covariates are declared, it SHALL also print the covariate names and the number of dropped samples.

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
