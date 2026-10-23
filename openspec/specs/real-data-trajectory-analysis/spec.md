# real-data-trajectory-analysis Specification

## Purpose

Defines the real-data entry point of MOTCO: how caller-supplied multi-omic blocks and sample metadata are analysed with the same validated pipeline the simulation harness measures — pooled standardization, stage-supervised PLS latent space, net-displacement trajectory statistics with RRPP — what the analysis records, how it is written down, and the `motco analyze` command that exposes it.

## Requirements

### Requirement: Real-data analysis accepts named omic blocks with declared unit transforms

The analysis SHALL accept one or more omic blocks, each identified by a caller-chosen name and holding a samples × features table with string feature names, together with a per-block unit transform that is either `identity` or `beta_to_mvalue` (the logit of methylation beta values). Blocks SHALL be processed in the order the caller declares them, and that order SHALL define the joint feature order, with every joint feature name prefixed by its block name so names cannot collide across blocks. Every block MUST have the same number of rows as the metadata, rows aligned to the metadata's sample identifiers in the same order, at least one feature, and only finite values; the analysis SHALL reject a violation with an error naming the block and the problem. The `beta_to_mvalue` transform SHALL reject values outside the open interval (0, 1). The analysis SHALL NOT normalize counts, filter features or impute values: the blocks are taken in the units the caller supplies.

#### Scenario: Two blocks are concatenated in declared order

- **WHEN** blocks named `rna` (200 features) and `atac` (300 features) are declared in that order
- **THEN** the joint feature matrix has 500 columns, the first 200 named `rna__<feature>` and the next 300 named `atac__<feature>`

#### Scenario: Methylation beta values are analysed on the M-value scale

- **WHEN** a block declares the `beta_to_mvalue` transform
- **THEN** every value is mapped through the logit before pooled standardization, and the record states that block's units as M-values

#### Scenario: Misaligned block is rejected

- **WHEN** a block's row order differs from the metadata's sample identifiers
- **THEN** the analysis raises an error naming that block and stating the rows are not aligned

#### Scenario: Non-finite value is rejected

- **WHEN** a block contains a missing or infinite value
- **THEN** the analysis raises an error naming the block, and no imputation is performed

### Requirement: Real-data analysis takes an explicit stage order and exactly two groups

The analysis SHALL require the metadata to carry a sample-identifier column, a group column and a stage column, all named by the caller. The caller MUST supply the stage order explicitly as a sequence of the distinct stage labels; the analysis SHALL NOT infer stage order by sorting labels. The supplied order MUST list every observed stage label exactly once and nothing else, and SHALL define the first and last stages used by the net-displacement orientation estimator and the order of stages reported everywhere. The group column MUST have exactly two distinct labels; the caller MAY supply their order, and otherwise they are ordered by sorted string representation. Every group × stage cell MUST be observed. Reported levels SHALL be the caller's original labels.

#### Scenario: Stage order drives the orientation measurement

- **WHEN** the stages are labelled `Low`, `Intermediate`, `High` and the caller declares the order `["Low", "Intermediate", "High"]`
- **THEN** net displacement is measured from the `Low` stage mean to the `High` stage mean, and the record lists the stages in that order

#### Scenario: Incomplete stage order is rejected

- **WHEN** the declared stage order omits an observed label or names an unobserved one
- **THEN** the analysis raises an error listing the declared and the observed labels

#### Scenario: Three groups are rejected

- **WHEN** the group column has three distinct labels
- **THEN** the analysis raises an error stating exactly two groups are required and naming the labels found

#### Scenario: Missing group-by-stage cell is rejected

- **WHEN** some group has no sample at some stage
- **THEN** the analysis raises an error naming the missing combination

### Requirement: Real-data analysis measures with the validated production pipeline

The analysis SHALL standardize every feature with pooled per-feature location and scale fitted on all samples (a feature with zero variance is left centered and unscaled), concatenate the blocks in declared order, and construct the latent space by stage-supervised PLS-DA whose rank is selected by double cross-validation with the harness's rank rule — the modal selected rank across outer repeats with the smaller rank breaking ties, from the candidate ranks at or above the rank floor, 2 by default — then fitted once on all samples at that rank. The caller MAY instead force a rank; the analysis SHALL then skip cross-validation and record the component selection as forced. The analysis SHALL build the group × stage design with interactions, estimate `delta`, `angle` and `shape` between the two groups in the latent space, and, when a positive permutation count is given, estimate their p-values by residual randomization under the reduced (no-interaction) model with a caller-supplied seed. The orientation estimator SHALL default to net displacement; the PC1 estimator SHALL be selectable by name. `shape` SHALL be reported unavailable with fewer than three stages.

#### Scenario: Cross-validated rank is recorded with its mode

- **WHEN** the analysis runs without a forced rank
- **THEN** the record states the selected rank, the component selection as cross-validated, the cross-validation settings and the mean cross-validated AUROC

#### Scenario: Cross-validated rank is at least 2 by default

- **WHEN** the analysis runs without a forced rank and without naming a floor
- **THEN** the selected rank is at least 2, the latent scores have at least 2 columns, and `angle` is measured in a space of at least two dimensions

#### Scenario: Forced rank skips cross-validation

- **WHEN** the caller forces a rank of 2
- **THEN** no cross-validation is run, the latent space has exactly 2 components, and the record states the component selection as forced

#### Scenario: Net displacement is the default estimator

- **WHEN** the analysis runs without naming an estimator
- **THEN** `angle` is the angle between the groups' latent net-displacement vectors from the first to the last declared stage, and the record states the estimator as net displacement

#### Scenario: Permutation p-values are reproducible

- **WHEN** the analysis runs twice with the same inputs, permutation count and seed
- **THEN** both runs return identical p-values and null summaries

### Requirement: Real-data analysis is equivalent to the simulation harness on identical data

On the same blocks, metadata, stage order and settings, including the rank floor, the real-data analysis SHALL produce the same latent matrix, observed `delta`/`angle`/`shape`, p-values, selected rank and null summary as the simulation evaluation harness's `pls` integration. The validated measurement SHALL exist in one implementation shared by both entry points, and the harness's outputs for every committed fixture and configuration SHALL be unchanged by the change, every committed configuration and fixture naming the floor of 1 under which it was produced.

#### Scenario: Harness and real-data analysis agree on a semi-synthetic dataset

- **WHEN** a semi-synthetic dataset is evaluated by the harness with `pls` integration, and its blocks and metadata are passed to the real-data analysis with the same stage order, cross-validation settings, rank floor, permutation count and seed
- **THEN** the latent matrices are equal, and the observed statistics, p-values, selected rank and null summaries are identical

#### Scenario: Harness and real-data analysis agree at their defaults

- **WHEN** the harness and the real-data analysis both run without naming a floor
- **THEN** both select the rank from candidates at or above 2 and their latent matrices, statistics, p-values and selected rank are identical

#### Scenario: Committed harness identity is preserved

- **WHEN** the committed byte-identity fixtures and committed configuration identity tests run after the change, at the floor of 1 they name
- **THEN** every recorded statistic, signature, cell id and matched seed is unchanged

### Requirement: Real-data analysis records the covariates the interpretability statement requires

The result SHALL carry, in a JSON-safe record: the pair statistics; the p-values; the per-statistic null summary (retained draw count, mean, sd and the 50/90/95/99 % quantiles); the stage-mean configuration spectrum in the latent space, pooled and per group, and the permutation eigengap summary when permutations ran; the selected rank and component selection; the orientation estimator; the group levels and the ordered stage levels; the sample counts per group × stage cell; each block's name, feature count and declared units; the joint feature order signature; the resolved analysis parameters; and the versions of MOTCO and its numerical dependencies. The result SHALL also expose the latent scores per sample, the PLS loadings per joint feature, and the fitted latent stage-mean vectors per group × stage cell.

#### Scenario: Selected rank accompanies the angle result

- **WHEN** an analysis completes
- **THEN** its record states the selected rank and the component-selection mode beside the `angle` statistic and p-value

#### Scenario: Software versions are recorded

- **WHEN** an analysis completes
- **THEN** its record states the MOTCO, numpy, scipy, scikit-learn and pandas versions it ran under

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

### Requirement: Real-data analysis writes a reproducible analysis artifact

The analysis SHALL be writable to a directory as: `result.json` (the JSON-safe record and resolved parameters), `latent_scores.csv` (samples × components, indexed by sample identifier), `pls_loadings.csv` (joint features × components, with block name and original feature name columns), `stage_means.csv` (group × stage latent LS-mean vectors, in declared stage order), and, when attribution ran, its tables. The writer SHALL create the directory if absent and SHALL refuse to overwrite an existing `result.json` unless overwriting is requested.

#### Scenario: Outputs are complete and loadable

- **WHEN** a result with permutations is written to an empty directory
- **THEN** the directory contains `result.json`, `latent_scores.csv`, `pls_loadings.csv` and `stage_means.csv`, and `result.json` parses as JSON holding the p-values and the selected rank

#### Scenario: Existing result is protected

- **WHEN** a result is written to a directory that already holds `result.json` without requesting overwrite
- **THEN** the writer raises an error and leaves the existing files unchanged

### Requirement: motco analyze subcommand

`motco analyze` SHALL be a CLI subcommand over the real-data analysis. It SHALL accept: a repeatable `--block NAME=PATH` (CSV, first column the sample identifier); a repeatable `--transform NAME=KIND` (default `identity` for any block not named); `--metadata PATH`, `--sample-col`, `--group-col` and `--stage-col`; a required `--stage-order` (comma-separated labels) and an optional `--group-order`; a repeatable `--covariate COL`; a repeatable `--categorical-covariate COL`, which also declares the column a covariate; `--missing-covariates` (`error` default, `complete_case` selectable); `--latent-adjustment` (`stage` or `group_stage`; absent by default, and rejected without a covariate); `--permutations` (default 0), `--seed`, `--n-jobs` and `--forced-components`; `--min-components` (default 2; `1` restores the unfloored rule); `--orientation-estimator` (`net_displacement` default, `pc1` selectable); `--attribution` with `--bootstrap-replicates`, `--bootstrap-seed` and `--top-k`; `--out-dir` (required); and `--overwrite`. Covariates SHALL be declared in the order the flags are given. Cross-validation settings SHALL be exposed with the harness defaults. Input and validation errors SHALL exit non-zero with a descriptive message. A successful run SHALL write the analysis artifact and print the pair statistics, p-values, selected rank, rank floor and estimator. When covariates are declared, it SHALL also print the covariate names and the number of dropped samples. When latent adjustment is requested, it SHALL also print the conditioning.

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

### Requirement: Real-data analysis records the rank floor

The analysis parameters SHALL accept an integer rank floor, default 2, with 1 restoring the unfloored rule. A floor of 1 SHALL be omitted from the resolved parameters and from the latent record, so a floor-1 analysis reproduces every 0.9.0 record and table byte for byte. A floor above 1 SHALL appear in the resolved parameters and in the latent block beside the selected rank. A floor below 1 SHALL be rejected before any fitting.

#### Scenario: Floor of 1 is byte-identical to 0.9.0

- **WHEN** the committed identity analyses run with a floor of 1
- **THEN** their records, including the resolved parameters, and their tables equal the committed 0.9.0 fixtures

#### Scenario: Default floor is recorded

- **WHEN** the analysis runs without naming a floor
- **THEN** the resolved parameters and the latent block state a floor of 2

#### Scenario: Floor below 1 is rejected

- **WHEN** the parameters name a floor of 0
- **THEN** the analysis raises an error before fitting the latent space
