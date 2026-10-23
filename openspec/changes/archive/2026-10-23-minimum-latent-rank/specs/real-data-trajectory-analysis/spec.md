# Spec Delta

## MODIFIED Requirements

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

## ADDED Requirements

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
