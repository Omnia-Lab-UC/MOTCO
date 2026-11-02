# orientation-driver-attribution Specification

## Purpose

This capability explains a detected group trajectory orientation difference in a shared PLS molecular space by identifying signed feature contrasts, the portion retained by PLS, omitted residual structure, and stability of the resulting drivers.

## Requirements

### Requirement: Attribution uses one shared PLS coordinate system

The system SHALL accept one aligned feature matrix, group and stage metadata, and one fitted pooled PLS model, and SHALL use that same model for both groups and every transition. It MUST NOT fit separate group-specific PLS models during attribution.

#### Scenario: Valid shared-model attribution

- **WHEN** the feature rows, metadata rows, feature order, and fitted model input dimensions agree and metadata contains exactly two groups and at least two stages
- **THEN** the system returns an attribution result using the shared model

#### Scenario: Invalid alignment is rejected

- **WHEN** feature rows, metadata rows, sample identifiers, or model input dimensions do not agree
- **THEN** the system raises a descriptive validation error before calculating driver values

#### Scenario: Unsupported design is rejected

- **WHEN** metadata contains fewer than two groups, more than two groups, fewer than two stages, or a missing group-by-stage cell
- **THEN** the system raises a descriptive validation error naming the invalid design

### Requirement: Attribution separates direction from transition magnitude

The system SHALL calculate, for each group and every adjacent ordered stage transition, the feature transition vector, its Euclidean norm, its unit direction when the norm is nonzero, and the signed directional contrast between the two groups. The group contrast SHALL use a documented deterministic group ordering and SHALL preserve feature signs.

#### Scenario: Two-stage directional contrast

- **WHEN** two groups have nonzero transitions between two stages
- **THEN** the result contains one directional contrast equal to the second ordered group's unit transition minus the first ordered group's unit transition, together with both group path lengths

#### Scenario: Multi-stage adjacent transitions

- **WHEN** the design contains three or more ordered stages
- **THEN** the result contains a separate directional contrast and path length pair for every adjacent stage transition

#### Scenario: Zero transition is unavailable

- **WHEN** either group's transition norm is zero or below the configured numerical tolerance
- **THEN** direction-dependent values for that transition are marked unavailable and the result retains the raw transition and path length information

### Requirement: Attribution decomposes observed direction into PLS-captured and residual components

The system SHALL compute observed group-stage feature means, project those means through the shared PLS model, reconstruct the retained feature-space component, and report observed, PLS-captured, and observed-minus-captured directional contrasts for every applicable transition. The decomposition SHALL use a fixed model and SHALL not replace observed features with low-rank reconstructions.

#### Scenario: PLS-captured contrast is reported

- **WHEN** a valid nonzero transition is available for both groups
- **THEN** the result contains observed and reconstructed PLS directional contrasts, their feature-level signed values, and a residual contrast defined as observed minus reconstructed

#### Scenario: Model remains frozen

- **WHEN** attribution calculates observed values, PLS reconstructions, or bootstrap replicates
- **THEN** the fitted PLS model and its coordinate system remain unchanged

#### Scenario: Reconstruction is unavailable for an invalid model

- **WHEN** the fitted model cannot transform or inverse-transform the aligned feature matrix
- **THEN** the system raises a descriptive attribution error rather than returning partial driver values

### Requirement: Results distinguish standardized and original-unit effects

The system SHALL report feature-level attribution in the input coordinate units and SHALL optionally report a second set of observed, PLS-captured, and residual effects in original feature units when valid per-feature scales are supplied. Unit labels and feature order SHALL be retained in the result.

#### Scenario: Standardized effects are returned

- **WHEN** attribution is run on the pooled standardized feature matrix
- **THEN** every feature-level effect identifies that it is in standardized input units

#### Scenario: Original-unit effects are requested

- **WHEN** positive finite per-feature scales matching the input feature order are supplied
- **THEN** the result reports corresponding original-unit effects without changing the standardized PLS geometry

#### Scenario: Invalid unit conversion is rejected

- **WHEN** supplied feature scales are missing, non-finite, non-positive, duplicated, or differently ordered from the input features
- **THEN** the system raises a descriptive validation error

### Requirement: Bootstrap estimates driver stability within design strata

The system SHALL optionally resample subjects with replacement independently within each group-by-stage stratum, while keeping the fitted preprocessing contract and PLS model fixed, and SHALL summarize feature sign stability, rank stability, and top-k selection frequency for each transition and attribution component.

#### Scenario: Stratified bootstrap is reproducible

- **WHEN** attribution is run twice with the same data, bootstrap count, and seed
- **THEN** the bootstrap summaries are identical

#### Scenario: Bootstrap preserves design strata

- **WHEN** a bootstrap replicate is generated
- **THEN** every group-by-stage stratum contributes the same number of resampled observations as the source stratum and no observation is sampled across strata

#### Scenario: Stability summaries are returned

- **WHEN** bootstrap resampling is enabled with a positive replicate count
- **THEN** the result includes per-feature sign stability, rank stability, top-k selection frequency, bootstrap configuration, and the number of valid replicates

#### Scenario: Bootstrap is disabled explicitly

- **WHEN** the bootstrap count is zero
- **THEN** no resampling is performed and stability fields are marked unavailable rather than treated as zero stability

### Requirement: Caller-supplied feature groups can be aggregated

The system SHALL optionally accept validated caller-supplied feature-to-group labels for correlated-feature modules or biological pathways and SHALL return group-level summaries of signed and absolute attribution, feature counts, and bootstrap stability when available. The system MUST NOT imply that the labels came from an internal pathway database or that aggregation establishes causality.

#### Scenario: Module aggregation

- **WHEN** multiple features share a supplied module label
- **THEN** the result contains one summary row per module and transition for observed, PLS-captured, and residual attribution

#### Scenario: Pathway labels are preserved

- **WHEN** supplied labels identify pathways
- **THEN** the result retains the caller-provided pathway names and records the label source as caller-supplied

#### Scenario: Incomplete mapping is rejected

- **WHEN** a supplied feature-group mapping omits an input feature, contains an unknown feature, or assigns conflicting labels
- **THEN** the system raises a descriptive validation error

### Requirement: Attribution output is machine-readable and interpretation-bounded

The system SHALL return structured results that retain group and stage ordering, transition identifiers, feature names, model/component metadata, effective bootstrap parameters, and unit metadata. The output SHALL state that attribution describes associations implied by the fitted shared PLS representation and SHALL not claim a unique inverse, biological causality, or inferential significance. The interpretation metadata SHALL state which quantity each decomposition explains: the principal-orientation contrast decomposes the tested `angle` estimand (principal-axis divergence), while per-adjacent-transition contrasts describe where along the trajectory the group directions diverge and coincide with the tested estimand only for two-stage designs.

#### Scenario: Result can be consumed without parsing prose

- **WHEN** a valid attribution analysis completes
- **THEN** feature-level, transition-level, principal-orientation, aggregate, and configuration data are available as structured records or tables

#### Scenario: Significance is not fabricated

- **WHEN** attribution is run without an external hypothesis-test result
- **THEN** the result contains no synthetic p-value or significance decision and leaves selection of significant global results to the calling workflow

#### Scenario: Decomposition targets are stated

- **WHEN** a valid attribution analysis completes
- **THEN** the interpretation metadata distinguishes the principal-orientation contrast as the decomposition of the tested `angle` estimand from the per-transition contrasts as a per-step description

### Requirement: Attribution decomposes the tested principal-orientation contrast

The system SHALL compute, for each group, the principal orientation of its centered stage-mean configuration in feature space — the leading principal axis, signed so it points along the trajectory's net displacement from first stage to last, using the same convention as the tested `angle` estimator — for both the observed group-stage means and their PLS-reconstructed counterparts. The system SHALL report the signed per-feature principal-orientation contrast between the two ordered groups as observed, PLS-captured, and residual (observed minus captured) components, beside the existing per-transition contrasts. The group ordering and sign conventions SHALL be identical to those used for transition contrasts.

#### Scenario: Principal contrast is estimator-consistent

- **WHEN** the principal orientation is computed from a group's observed stage means
- **THEN** it equals, within numerical tolerance, the orientation the trajectory `angle` estimator returns for the same stage-mean configuration

#### Scenario: Two-stage principal contrast equals the transition contrast

- **WHEN** the design contains exactly two stages and both groups have nonzero transitions
- **THEN** the observed, PLS-captured, and residual principal-orientation contrasts equal the corresponding single-transition contrasts exactly

#### Scenario: Multi-stage principal contrast is reported beside transitions

- **WHEN** the design contains three or more ordered stages
- **THEN** the result contains one principal-orientation contrast per attribution component in addition to the per-adjacent-transition contrasts, in the same structured and tabular forms

### Requirement: Principal-orientation degeneracy is qualified, not silently reported

The system SHALL report, for each group and for both the observed and PLS-reconstructed configurations, the relative eigengap of the centered stage-mean configuration, and SHALL mark the principal-orientation contrast as degenerate when any contributing configuration's relative eigengap falls below a documented threshold. A degenerate principal contrast SHALL remain available in the output but MUST carry an explicit degeneracy flag; the system MUST NOT present a near-isotropic configuration's principal axis as a well-defined orientation.

#### Scenario: Well-separated configuration passes the gate

- **WHEN** every contributing stage-mean configuration has a relative eigengap at or above the documented threshold
- **THEN** the principal-orientation contrast is reported without a degeneracy flag and the per-group eigengaps are included in the result

#### Scenario: Near-isotropic configuration is flagged

- **WHEN** a contributing stage-mean configuration has a relative eigengap below the documented threshold
- **THEN** the principal-orientation contrast carries a degeneracy flag identifying which group and configuration triggered it

#### Scenario: Closed trajectory is unavailable

- **WHEN** a group's net displacement from first to last stage is zero or below the configured numerical tolerance
- **THEN** the principal orientation for that group is marked unavailable, mirroring the zero-transition handling

### Requirement: Bootstrap stability covers the principal-orientation contrast

When bootstrap resampling is enabled, the system SHALL summarize feature sign stability, rank stability, and top-k selection frequency for the principal-orientation contrast with the same stratified resampling contract as transition contrasts, and SHALL anchor each replicate's principal-axis sign by that replicate's net displacement so resampling noise cannot flip the axis arbitrarily.

#### Scenario: Principal-orientation bootstrap is reproducible

- **WHEN** attribution is run twice with the same data, bootstrap count, and seed
- **THEN** the principal-orientation stability summaries are identical

#### Scenario: Replicate sign is anchored

- **WHEN** a bootstrap replicate's stage-mean configuration is resampled
- **THEN** its principal axis is signed by that replicate's own net displacement before contrasting, using the same convention as the point estimate

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
