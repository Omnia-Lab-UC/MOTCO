# Spec Delta

## ADDED Requirements

### Requirement: Attribution can use covariate-adjusted cell means

When the caller supplies row-aligned numeric covariate design columns, the system SHALL compute each group × stage feature mean as the least-squares mean of a model with one indicator per group × stage cell plus the covariates, evaluated at a reference covariate row, by default the covariate column means. It SHALL use those means for the observed, PLS-captured and residual components. Covariates enter additively.

#### Scenario: Known covariate effect is removed

- **WHEN** the features are a pure group × stage cell-mean structure plus a covariate times a fixed feature vector, with the covariate's mean differing across cells and no within-cell variation correlated with it
- **THEN** the adjusted transitions, path lengths, transition contrasts and principal-orientation contrast equal those of the cell-mean structure alone, and the contrasts from unadjusted means differ from them

#### Scenario: Adjusted means are the model's least-squares means

- **WHEN** attribution runs with covariates on any valid input
- **THEN** each reported mean equals, to numerical tolerance, the fitted cell coefficient plus the covariate coefficients times the reference row, from an independent least-squares fit of the features on the cell indicators and the covariates

#### Scenario: Contrasts do not depend on the reference row

- **WHEN** the same adjusted attribution is evaluated at the covariate column means and at any other fixed covariate row
- **THEN** the transitions, path lengths, transition contrasts and principal-orientation contrasts are equal up to floating-point tolerance, and the observed means differ by the same vector in every cell

#### Scenario: Covariates collinear within cells are rejected

- **WHEN** a covariate takes one value within every group × stage cell, so it is collinear with the cell indicators
- **THEN** the system raises a descriptive attribution error naming the covariate columns before calculating driver values

#### Scenario: Invalid covariate input is rejected

- **WHEN** the covariate table has a different number of rows from the features, contains a missing or non-finite value, has duplicated column names, or the reference row names different columns or holds a non-finite value
- **THEN** the system raises a descriptive validation error

#### Scenario: Precomputed means and covariates are exclusive

- **WHEN** both a precomputed mean table and covariates are supplied
- **THEN** the system raises an error stating only one mean source may be given

### Requirement: Covariate-adjusted bootstrap refits the adjustment in each resample

When bootstrap resampling is enabled with covariates, each replicate SHALL draw the same stratified within-cell resample of rows that attribution without covariates draws for the same seed. It SHALL refit the covariate coefficients on the resampled rows and evaluate the means at the fixed reference row, keeping the preprocessing and PLS model frozen. A resample whose covariates are collinear within cells SHALL be an invalid replicate, and invalid replicates SHALL be counted.

#### Scenario: Adjusted and unadjusted bootstraps resample the same rows

- **WHEN** attribution runs on the same data, bootstrap count and seed, once with covariates and once without
- **THEN** both resample the same row indices in every replicate

#### Scenario: Adjusted bootstrap is reproducible

- **WHEN** adjusted attribution is run twice with the same data, covariates, reference row, bootstrap count and seed
- **THEN** all bootstrap summaries are identical

#### Scenario: The adjustment is refitted per replicate

- **WHEN** a bootstrap replicate resamples the rows
- **THEN** its means use covariate coefficients estimated from that replicate's rows, not the full-sample coefficients, at the full-sample reference row

#### Scenario: Rank-deficient resample is invalid

- **WHEN** a covariate's indicator is non-zero in only one sample of one cell, so some resamples omit it
- **THEN** those replicates contribute no contrast, every component's valid replicate count excludes them, and the result records how many replicates were rank-deficient

### Requirement: Covariate-adjusted attribution is recorded and interpretation-bounded

With covariates, the configuration SHALL state the mean source as covariate-adjusted, the covariate columns in order, the reference row and the number of rank-deficient bootstrap replicates. The interpretation SHALL state that adjusted effects assume an additive linear covariate model, are evaluated at the reference row, and do not remove unmeasured confounding. Without covariates, every output SHALL equal the output produced before covariate support existed.

#### Scenario: Adjusted configuration is recorded

- **WHEN** attribution runs with the covariates `age` and `pmi` and 50 bootstrap replicates
- **THEN** the configuration table states mean source `covariate_adjusted`, the columns `age` and `pmi`, their reference values and the count of rank-deficient replicates, and the interpretation table holds the adjustment boundary statement

#### Scenario: No covariates means no change

- **WHEN** attribution runs without covariates on a fixed dataset, model and settings
- **THEN** every table, including the configuration and interpretation tables, equals the one produced before covariate support existed, and no adjustment row appears
