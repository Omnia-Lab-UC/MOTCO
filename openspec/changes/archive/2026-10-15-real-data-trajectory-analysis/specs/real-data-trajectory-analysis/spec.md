# Spec Delta

## Purpose

Defines the real-data entry point of MOTCO: how caller-supplied multi-omic blocks and sample metadata are analysed with the same validated pipeline the simulation harness measures — pooled standardization, stage-supervised PLS latent space, net-displacement trajectory statistics with RRPP — what the analysis records, how it is written down, and the `motco analyze` command that exposes it.

## ADDED Requirements

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

The analysis SHALL standardize every feature with pooled per-feature location and scale fitted on all samples (a feature with zero variance is left centered and unscaled), concatenate the blocks in declared order, and construct the latent space by stage-supervised PLS-DA whose rank is selected by double cross-validation with the harness's rank rule — the modal selected rank across outer repeats with the smaller rank breaking ties — then fitted once on all samples at that rank. The caller MAY instead force a rank; the analysis SHALL then skip cross-validation and record the component selection as forced. The analysis SHALL build the group × stage design with interactions, estimate `delta`, `angle` and `shape` between the two groups in the latent space, and, when a positive permutation count is given, estimate their p-values by residual randomization under the reduced (no-interaction) model with a caller-supplied seed. The orientation estimator SHALL default to net displacement; the PC1 estimator SHALL be selectable by name. `shape` SHALL be reported unavailable with fewer than three stages.

#### Scenario: Cross-validated rank is recorded with its mode

- **WHEN** the analysis runs without a forced rank
- **THEN** the record states the selected rank, the component selection as cross-validated, the cross-validation settings and the mean cross-validated AUROC

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

On the same blocks, metadata, stage order and settings, the real-data analysis SHALL produce the same latent matrix, observed `delta`/`angle`/`shape`, p-values, selected rank and null summary as the simulation evaluation harness's `pls` integration. The validated measurement SHALL exist in one implementation shared by both entry points, and the harness's outputs for every committed fixture and configuration SHALL be unchanged by the change.

#### Scenario: Harness and real-data analysis agree on a semi-synthetic dataset

- **WHEN** a semi-synthetic dataset is evaluated by the harness with `pls` integration, and its blocks and metadata are passed to the real-data analysis with the same stage order, cross-validation settings, permutation count and seed
- **THEN** the latent matrices are equal, and the observed statistics, p-values, selected rank and null summaries are identical

#### Scenario: Committed harness identity is preserved

- **WHEN** the committed byte-identity fixtures and committed configuration identity tests run after the change
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

When requested, the analysis SHALL run orientation attribution on the fitted pooled PLS estimator and the standardized joint feature matrix it was fitted to, with the caller's stage order and group labels, a caller-chosen bootstrap replicate count and seed, and a top-k size; it SHALL never refit the latent space to do so. The compact attribution record SHALL report the truth-recovery block as unavailable, because real data carry no generator truth, and the full attribution tables SHALL be available for writing. Attribution SHALL be disabled by default.

#### Scenario: Attribution reuses the fitted estimator

- **WHEN** attribution is requested
- **THEN** the attribution's model component count equals the analysis's selected rank, and the estimator is not refitted

#### Scenario: Truth recovery is unavailable on real data

- **WHEN** attribution is requested
- **THEN** the compact record marks truth recovery unavailable with a reason, and no driver set is invented

### Requirement: Real-data analysis writes a reproducible analysis artifact

The analysis SHALL be writable to a directory as: `result.json` (the JSON-safe record and resolved parameters), `latent_scores.csv` (samples × components, indexed by sample identifier), `pls_loadings.csv` (joint features × components, with block name and original feature name columns), `stage_means.csv` (group × stage latent LS-mean vectors, in declared stage order), and, when attribution ran, its tables. The writer SHALL create the directory if absent and SHALL refuse to overwrite an existing `result.json` unless overwriting is requested.

#### Scenario: Outputs are complete and loadable

- **WHEN** a result with permutations is written to an empty directory
- **THEN** the directory contains `result.json`, `latent_scores.csv`, `pls_loadings.csv` and `stage_means.csv`, and `result.json` parses as JSON holding the p-values and the selected rank

#### Scenario: Existing result is protected

- **WHEN** a result is written to a directory that already holds `result.json` without requesting overwrite
- **THEN** the writer raises an error and leaves the existing files unchanged

### Requirement: motco analyze subcommand

`motco analyze` SHALL be a CLI subcommand over the real-data analysis. It SHALL accept a repeatable `--block NAME=PATH` (CSV, first column the sample identifier), a repeatable `--transform NAME=KIND` (default `identity` for any block not named), `--metadata PATH`, `--sample-col`, `--group-col`, `--stage-col`, a required `--stage-order` (comma-separated labels), an optional `--group-order`, `--permutations` (default 0), `--seed`, `--n-jobs`, `--forced-components`, `--orientation-estimator` (`net_displacement` default, `pc1` selectable), `--attribution` with `--bootstrap-replicates`, `--bootstrap-seed` and `--top-k`, `--out-dir` (required) and `--overwrite`. Cross-validation settings SHALL be exposed with the harness defaults. Input and validation errors SHALL exit non-zero with a descriptive message; a successful run SHALL write the analysis artifact and print the pair statistics, p-values, selected rank and estimator.

#### Scenario: Two-block analysis from CSV files

- **WHEN** `motco analyze --block rna=rna.csv --block atac=atac.csv --metadata meta.csv --sample-col donor --group-col Sex --stage-col ADNC --stage-order "Low,Intermediate,High" --permutations 99 --seed 1 --out-dir out/` runs on aligned files
- **THEN** the command exits zero, `out/result.json` exists with p-values for `delta`, `angle` and `shape`, and the printed summary names the selected rank and `net_displacement`

#### Scenario: Unknown block in a transform is rejected

- **WHEN** `--transform methyl=beta_to_mvalue` names a block that was not declared with `--block`
- **THEN** the command exits non-zero and the message names the undeclared block

#### Scenario: Missing stage order is rejected

- **WHEN** `--stage-order` is omitted
- **THEN** the command exits non-zero with a message stating the stage order is required
