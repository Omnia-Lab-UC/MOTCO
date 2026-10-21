# covariate-residualized-latent-space Specification

## Purpose
Defines the opt-in covariate-aware construction of the PLS-DA latent space in the real-data analysis: the features are residualized on the declared covariates, conditional on a chosen design, before the basis is fitted, so the basis is not spent on an axis aligned with the covariates.

## Requirements

### Requirement: Latent adjustment is opt-in and requires declared covariates

The analysis SHALL accept a latent-adjustment setting with the values `stage`, `group_stage` or absent (the default). The setting SHALL be rejected with an error when no covariates are declared. When the setting is absent, the latent space SHALL be built from the standardized features and the stage labels alone, exactly as without this capability.

#### Scenario: Setting without covariates is rejected

- **WHEN** the caller requests latent adjustment and declares no covariates
- **THEN** the analysis raises an error before any fitting, stating that latent adjustment requires declared covariates

#### Scenario: Absent setting leaves the basis covariate-blind

- **WHEN** covariates are declared and the setting is absent
- **THEN** the latent scores, loadings and selected rank equal those of an analysis without latent adjustment support

### Requirement: Features are residualized on the covariates conditional on the chosen design

With the setting on, the analysis SHALL regress every pooled-standardized joint feature on an intercept, the conditioning design and the declared covariate design columns, and subtract the fitted covariate part. Under `stage` the conditioning design SHALL be the stage indicators; under `group_stage` the group, stage and group × stage indicators. The residualization covariates SHALL be exactly the declared covariates.

#### Scenario: Stage-collinear covariate signal is kept

- **WHEN** a feature equals a stage effect plus a covariate effect, and the covariate is correlated with stage
- **THEN** the residualized feature retains the stage effect and its covariate coefficient equals the within-stage effect, not the marginal one

#### Scenario: Covariate effect is removed

- **WHEN** a feature equals a fixed vector times a covariate plus noise
- **THEN** the residualized feature has a within-design correlation with the covariate of zero up to floating-point tolerance

### Requirement: Residualized features keep the pooled scale

The residualized features SHALL keep the pooled per-feature scale fitted before residualization and SHALL NOT be re-standardized before or inside the basis fit.

#### Scenario: Scale is kept

- **WHEN** a feature's variance is mostly explained by the covariates
- **THEN** its residualized variance is smaller than one and the basis fit does not rescale it to unit variance

### Requirement: Residualization is estimated on training samples only inside cross-validation

During rank selection, every training set at both levels of the double cross-validation SHALL estimate its own residualization coefficients, residualize its own features with them, and residualize the held-out features with the same coefficients and the held-out samples' covariates alone. Held-out labels SHALL NOT be read. The final basis SHALL use coefficients estimated on all analysed samples.

#### Scenario: Held-out labels are not read

- **WHEN** the held-out samples' labels are permuted before rank selection while their features and covariates are kept
- **THEN** every training-set fit and every held-out residualized feature matrix is unchanged

#### Scenario: Fold coefficients differ from the full-sample coefficients

- **WHEN** a training set excludes samples whose covariates differ from the rest
- **THEN** that fold residualizes with its own coefficients, not the full-sample ones

### Requirement: Cross-validation folds are stratified on the conditioning cells

With the setting on, both levels of the double cross-validation SHALL stratify their folds on the conditioning cells (the stages under `stage`, the group × stage cells under `group_stage`), so every training set contains every conditioning cell.

#### Scenario: Group-sorted rows keep every cell in training

- **WHEN** `group_stage` conditioning is requested on rows sorted by group within stage
- **THEN** every cross-validation training set contains every group × stage cell, and rank selection completes

#### Scenario: Stage-conditioned folds equal the unadjusted folds

- **WHEN** `stage` conditioning is requested
- **THEN** the cross-validation folds are those of the analysis without latent adjustment

### Requirement: Every residualization fit must have full column rank

The residualization design (intercept, conditioning and covariate columns) MUST have full column rank in the full-sample fit and in every cross-validation training set. Otherwise the analysis SHALL raise an error naming the covariates, and for a training set, the repeat and fold.

#### Scenario: Collinear covariate fails the full-sample fit

- **WHEN** a numeric covariate equals the stage index under `stage` conditioning
- **THEN** the analysis raises an error stating the residualization design is not of full rank and naming the covariate

#### Scenario: Factor level absent from a training set fails that fold

- **WHEN** every sample with one categorical covariate level lands in the held-out part of one fold
- **THEN** the analysis raises an error naming that repeat and fold

### Requirement: Trajectory statistics are measured on the raw features in the residualized basis

With the setting on, the latent scores SHALL be the pooled-standardized features, centered at the residualized features' mean, projected onto the basis fitted on the residualized features. The trajectory model, LS-means and permutation scheme SHALL be those of the covariate-adjusted analysis with the same covariates. Projecting the residualized features instead SHALL give the same statistics and p-values.

#### Scenario: Raw and residualized projections agree in inference

- **WHEN** the same analysis is run once with the raw features projected and once with the residualized features projected, with the same seed
- **THEN** `delta`, `angle`, `shape`, their p-values and null summaries are equal up to floating-point tolerance

#### Scenario: Covariate coefficients remain in the record

- **WHEN** a latent-adjusted analysis completes
- **THEN** the covariate block's coefficients are the fitted covariate effects on the projected raw features and are not identically zero

### Requirement: Latent adjustment is recorded

With the setting on, the record SHALL carry a latent-adjustment block stating the conditioning, the residualization design column names, the residual degrees of freedom of the full-sample fit, the smallest training-set residual degrees of freedom seen during cross-validation, the feature scale, and a per-feature R² summary (the 5/25/50/75/95 % quantiles, the mean and the fraction above 0.5). The latent block SHALL name the conditioning. With the setting absent, neither key SHALL appear.

#### Scenario: Block is present with the setting on

- **WHEN** a `stage`-conditioned analysis with three covariates and three stages on 80 samples completes
- **THEN** the latent-adjustment block states conditioning `stage`, five design columns beyond the intercept, residual degrees of freedom 74, and the R² summary

#### Scenario: Block is absent with the setting off

- **WHEN** an analysis runs without the setting
- **THEN** the record has no latent-adjustment key and its latent block has no conditioning entry
