# simulation-evaluation-harness Specification

## Purpose
Define the per-replicate evaluation layer that runs one semi-synthetic trajectory dataset through MOTCO integration, trajectory statistics, optional RRPP p-value estimation, and result metadata capture.

## Requirements

### Requirement: Harness evaluates one semi-synthetic trajectory dataset
MOTCO SHALL provide a simulation evaluation harness that accepts one `SemiSyntheticTrajectoryDataset` and returns a structured evaluation result.

#### Scenario: Successful single-dataset evaluation
- **WHEN** a caller provides a valid semi-synthetic trajectory dataset and evaluation parameters
- **THEN** the harness returns observed trajectory statistics and evaluation metadata

#### Scenario: Result includes generator truth
- **WHEN** evaluation succeeds
- **THEN** the result includes the generator truth metadata from the input dataset

### Requirement: Harness supports initial integration methods
The harness SHALL construct the molecular latent space — the measurement substrate in which trajectory geometry is estimated — via a selectable integration method, operating on M-value-converted methylation and raw expression/proteomics. The production latent-space methods are **SNF** (graph-spectral embedding) and **PLS** (the transform of the omic features into the subspace that maximizes covariance with the stage label). `concat` is retained as a **baseline/diagnostic** path (standardized feature concatenation), not a constructed latent space. The viz down-projection (`plot_trajectory_from_*`) is display-only and distinct from this measurement space.

An optional integration parameter SHALL select which omic blocks are measured: a non-empty, duplicate-free subset of methylation, expression, and proteomics. Only the selected blocks SHALL enter pooled preprocessing, the integration method, orientation attribution, and every realized-geometry checkpoint (population, standardized population, observed standardized, and latent). Every joint scope SHALL be computed over the selected blocks only, and per-block scopes SHALL be reported for the selected blocks only. Blocks that are not selected SHALL be generated but not measured. When the parameter is absent, all three blocks SHALL be measured, and the evaluation result SHALL be identical to the harness's output before the parameter existed.

#### Scenario: Concatenated baseline integration
- **WHEN** the caller selects `concat` integration
- **THEN** the harness converts methylation to M-values, standardises all selected layers, and concatenates them into the outcome matrix
- **AND** the result metadata identifies `concat` as a baseline rather than a constructed latent space

#### Scenario: SNF integration
- **WHEN** the caller selects `snf` integration
- **THEN** the harness converts methylation to M-values and creates the latent space from SNF fusion of the selected layers and spectral embedding

#### Scenario: PLS integration
- **WHEN** the caller selects `pls` integration
- **THEN** the harness converts methylation to M-values, standardises the selected layers, fits PLS-DA conditioned on the stage label, and returns the PLS X-score matrix as the latent space
- **AND** the number of latent variables is selected by the double nested cross-validation (modal LV across repeats, parsimony tie-break) from the candidate ranks at or above the rank floor (2 by default) to secure a stable, non-overfitted molecular space in which orientation and shape are defined
- **AND** the result metadata records the selected number of latent variables and the cross-validation parameters

#### Scenario: PLS integration is infeasible
- **WHEN** the caller selects `pls` integration but the sample provides too few observations per stage for the cross-validation
- **THEN** the harness raises a clear validation error

#### Scenario: Unsupported integration method
- **WHEN** the caller selects an unsupported integration method
- **THEN** the harness raises a clear validation error

#### Scenario: A block subset is measured
- **WHEN** the caller selects the layers methylation and expression
- **THEN** the outcome matrix, the latent space, every realized-geometry checkpoint's joint scope, and every attribution feature record contain methylation and expression features only
- **AND** realized geometry reports no proteomics per-block scope
- **AND** the result metadata records the selected layers in canonical order

#### Scenario: Absent layer selection is byte-identical
- **WHEN** the evaluation parameters omit the layer selection
- **THEN** the evaluation result equals the result produced before the selection existed, for the same dataset and parameters

#### Scenario: Invalid layer selection is rejected
- **WHEN** the layer selection is empty, names an unknown layer, or repeats a layer
- **THEN** the harness raises a clear validation error naming the allowed layers

### Requirement: Harness builds MOTCO trajectory design objects
The harness SHALL construct model matrices, LS means, and trajectory contrasts from generated sample metadata. Level ordering SHALL be numeric-aware: when every distinct label of a factor (group or stage) parses as an integer, the levels are ordered by numeric value; otherwise the levels keep the deterministic lexicographic order of their string representations. The same ordering rule MUST govern the model matrix columns, the LS-mean rows, and the contrast indices so they stay mutually consistent.

#### Scenario: Design objects are derived from metadata
- **WHEN** the dataset metadata contains valid `group` and `stage` columns
- **THEN** the harness builds full and reduced model matrices, LS means, and a two-group trajectory contrast

#### Scenario: Missing metadata columns are rejected
- **WHEN** required metadata columns are missing
- **THEN** the harness raises a clear validation error

#### Scenario: Ten or more integer-labeled stages order numerically
- **WHEN** the metadata's stage labels are the strings "0" through "11"
- **THEN** the design's stage-level order is 0, 1, 2, ..., 11 (numeric), not the lexicographic order that places "10" and "11" before "2", and LS-mean rows and contrast indices follow the same order

#### Scenario: Non-numeric labels keep lexicographic order
- **WHEN** any factor level label does not parse as an integer
- **THEN** that factor's levels are ordered lexicographically by string representation, identical to the previous behavior

### Requirement: Harness estimates observed trajectory differences
The harness SHALL estimate observed `deltas`, `angles`, and `shapes` using MOTCO trajectory routines.

#### Scenario: Observed statistics are returned
- **WHEN** evaluation succeeds
- **THEN** the result includes observed `deltas`, `angles`, and `shapes`

#### Scenario: Pairwise group statistic is exposed
- **WHEN** evaluation succeeds for two groups
- **THEN** the result includes scalar statistics for the generated group comparison

### Requirement: Harness optionally runs RRPP p-value estimation
The harness SHALL optionally run RRPP and compute empirical p-values for observed trajectory statistics.

#### Scenario: RRPP disabled
- **WHEN** the caller sets permutation count to 0
- **THEN** the harness returns observed statistics without RRPP p-values

#### Scenario: RRPP enabled
- **WHEN** the caller sets permutation count greater than 0
- **THEN** the harness runs RRPP and returns empirical p-values for available statistics

#### Scenario: P-values use plus-one correction
- **WHEN** RRPP p-values are computed
- **THEN** the harness uses `(1 + count(null >= observed)) / (1 + n_permutations)`

### Requirement: Harness records evaluation metadata
The harness SHALL record parameters and runtime metadata needed for later grid aggregation.

#### Scenario: Evaluation metadata is returned
- **WHEN** evaluation succeeds
- **THEN** the result includes integration method, integration parameters, permutation count, runtime seconds, and evaluation parameters

### Requirement: Harness can compute orientation attribution from the evaluated PLS representation

The harness SHALL optionally compute orientation-attribution diagnostics for a PLS evaluation using the exact pooled standardized feature matrix and fitted PLS estimator that produced the trajectory score matrix. It MUST keep preprocessing, feature order, selected component count, and the PLS estimator fixed while calculating observed, PLS-captured, residual, and bootstrap values.

#### Scenario: Attribution is enabled for PLS

- **WHEN** a caller enables attribution diagnostics for a valid PLS evaluation
- **THEN** the harness returns attribution diagnostics calculated from the same standardized features and fitted PLS estimator used for trajectory measurement
- **AND** no second PLS fit or component selection is performed for attribution

#### Scenario: Attribution is disabled

- **WHEN** attribution diagnostics are not enabled
- **THEN** evaluation returns the existing statistics, p-values, geometry, and metadata without attribution work

#### Scenario: Attribution is requested for a non-PLS method

- **WHEN** a caller enables PLS orientation attribution with `concat` or `snf` integration
- **THEN** the harness rejects the incompatible configuration with a descriptive validation error

### Requirement: Attribution diagnostics are bounded and machine-readable

The harness SHALL return a JSON-safe diagnostic summary containing the effective attribution configuration, ordered transitions, observed/PLS-captured/residual transition metrics, configured top-k signed feature records in pooled standardized units and in original units wherever the fitted preprocessor supplies a positive per-feature scale, bootstrap sign and selection stability, and available recovery metrics against generator truth. It MUST NOT place a fitted estimator, full standardized matrix, full bootstrap matrix, or unrestricted feature table in the persisted diagnostic payload.

#### Scenario: Compact diagnostic is produced

- **WHEN** attribution completes for an eligible replicate
- **THEN** the result contains compact transition and top-k feature records with feature identifiers, component labels, signs, ranks, and stability values
- **AND** each feature record carries its standardized effect and, when a per-feature scale is available, its original-unit effect, with the unit basis labeled so M-value methylation is not read as beta values
- **AND** the payload records the attribution seed, bootstrap count, top-k value, selected PLS component count, feature order signature, and diagnostic schema version

#### Scenario: Generator truth is available

- **WHEN** the generated dataset identifies features whose group-stage differential pattern changed
- **THEN** the diagnostic reports top-k precision, recall, and selection counts against that truth using a documented truth definition that includes propagated omic effects

#### Scenario: A transition is unavailable

- **WHEN** either group has a zero or near-zero transition for an attribution component
- **THEN** the diagnostic marks the affected metrics unavailable and preserves transition identity and path lengths without treating unavailable values as zero

### Requirement: PLS fitting remains single-pass within an evaluation

The harness SHALL fit one final pooled PLS estimator after component selection. It SHALL measure trajectories on that estimator's training scores or, under orthonormal span coordinates, on its training inputs projected on an orthonormal basis of its rotations' span. It SHALL use that same estimator for attribution projection and reconstruction.

#### Scenario: Scores and attribution share one final fit

- **WHEN** a PLS evaluation with attribution succeeds
- **THEN** the score matrix equals the training scores of the estimator supplied to attribution
- **AND** the selected component count and score values remain equivalent to a PLS evaluation with attribution disabled under the same inputs and seeds

#### Scenario: Orthonormal coordinates reuse the final fit

- **WHEN** a PLS evaluation runs under orthonormal span coordinates
- **THEN** exactly one final PLS estimator is fitted, and the latent matrix is that estimator's training inputs projected on the orthonormal basis of its rotations' span

### Requirement: Harness summarizes each permutation null

Whenever the harness runs RRPP, it SHALL return a compact per-statistic summary of the permutation null distribution alongside the observed statistics and p-values. The summary MUST locate the observed statistic against its own null: for each statistic with a null distribution it carries the number of retained permutation draws, the null mean, the null standard deviation, and null quantiles including the median and the alpha-level upper critical value. The summary MUST be JSON-safe scalars only, and non-finite draws MUST be excluded from the summary with the retained count reported.

The summary is separate from the existing opt-in retention of full null distributions: it is produced whenever permutations are run, whereas the full null vectors remain available only when the caller asks for them.

#### Scenario: Summary accompanies every RRPP evaluation

- **WHEN** a caller evaluates a dataset with permutation count greater than 0
- **THEN** the result carries a null summary for each statistic that has a null distribution
- **AND** each summary reports the retained draw count, mean, standard deviation, and quantiles for that statistic

#### Scenario: Summary is independent of full-distribution retention

- **WHEN** a caller runs RRPP without requesting retention of the full null distributions
- **THEN** the result still carries the compact null summary
- **AND** the full null vectors are absent

#### Scenario: No permutations means no summary

- **WHEN** a caller sets permutation count to 0
- **THEN** the result carries no null summary and the existing observed statistics are unchanged

#### Scenario: Summary does not alter the test

- **WHEN** the null summary is produced for an evaluation
- **THEN** the observed statistics, p-values, permutation draws, and every other field of the result are identical to those the same inputs produced before the summary existed

### Requirement: Harness records the latent stage-mean configuration spectrum

On every evaluation, the harness SHALL compute and return a JSON-safe spectrum description of the centered stage-mean configuration in the evaluated latent space — for the pooled configuration and for each contrast group — derived from the same fitted LS-mean vectors the trajectory statistics use. Each description MUST carry the normalized eigenvalue spectrum and the relative eigengap (λ₁−λ₂)/Σλ as a named scalar. A configuration whose eigengap is undefined (zero total variance) MUST record that explicitly rather than a non-finite value, and a two-stage configuration's trivially saturated eigengap MUST be recorded as computed.

#### Scenario: Spectrum accompanies every evaluation

- **WHEN** a caller evaluates a dataset through the harness
- **THEN** the result carries a configuration-spectrum block with one entry for the pooled configuration and one per contrast group
- **AND** each entry reports the normalized spectrum and the relative eigengap as JSON-safe scalars

#### Scenario: Spectrum matches the measured configurations

- **WHEN** the spectrum block is computed for an evaluation
- **THEN** its values equal an independent singular value decomposition of the centered LS-mean configurations implied by the evaluation's fitted betas, contrast, and stage ordering

### Requirement: Harness summarizes the permutation eigengap distribution

Whenever the harness runs RRPP, it SHALL additionally summarize the pooled-configuration relative eigengap across the permutation draws — retained draw count, mean, standard deviation, and quantiles — so a replicate can locate its observed eigengap against its own permutation distribution. Full per-permutation spectra MUST NOT be retained in the result. The spectrum computation MUST consume no randomness: with identical inputs, seed, and worker count, the permutation draws, null distributions, observed statistics, and p-values MUST be identical to those produced before this capability existed.

#### Scenario: Permutation eigengap summary accompanies RRPP runs

- **WHEN** a caller evaluates a dataset with permutation count greater than 0
- **THEN** the result carries a summary of the pooled eigengap over the permutation draws with retained count, mean, standard deviation, and quantiles
- **AND** no per-permutation spectrum vectors are present in the result

#### Scenario: No permutations means no permutation summary

- **WHEN** a caller sets permutation count to 0
- **THEN** the result carries the observed configuration-spectrum block but no permutation eigengap summary

#### Scenario: Spectrum recording does not alter the test

- **WHEN** the same evaluation inputs are run with and without spectrum recording available
- **THEN** the observed statistics, p-values, permutation draws, and every pre-existing result field are identical

### Requirement: PLS integration supports a diagnostic forced-component override

The harness SHALL accept an optional integration parameter that forces the number of PLS latent variables, bypassing cross-validated component selection, for rank diagnostics only. When the override is absent the behavior MUST be unchanged. When the override is set, the harness MUST skip the double nested cross-validation, fit the single pooled PLS estimator at the forced component count, and record in the integration metadata both the forced count and an explicit marker that component selection was overridden, so no downstream consumer can mistake a forced-rank run for a production evaluation. The override MUST be validated against the feasible component range, which is at least 2 and at most the feasible maximum regardless of the cross-validation rank floor, and rejected with a clear error otherwise.

#### Scenario: Default behavior is unchanged

- **WHEN** a PLS evaluation runs without the forced-component parameter
- **THEN** the number of latent variables is selected by the double nested cross-validation exactly as before
- **AND** the integration metadata carries no override marker

#### Scenario: Forced component count is honored and recorded

- **WHEN** a PLS evaluation runs with the forced-component parameter set to a feasible value
- **THEN** the latent space has exactly that many components and no cross-validation is performed
- **AND** the integration metadata records the forced count and marks component selection as overridden

#### Scenario: Infeasible forced count is rejected

- **WHEN** the forced-component parameter exceeds the feasible range for the sample
- **THEN** the harness raises a clear validation error

#### Scenario: Forced rank is independent of the floor

- **WHEN** a PLS evaluation runs with the forced-component parameter set to 2 and a cross-validation rank floor of 3
- **THEN** the latent space has exactly 2 components, no cross-validation is performed, and the floor is not applied

### Requirement: Harness supports diagnostic PCA and within-group PLS integration

The harness SHALL accept two additional **diagnostic** integration methods for measurement-design comparison. Both operate on the same preprocessing and the same layer selection as `pls` (M-value methylation, pooled standardization, selected blocks only). Neither is a production latent space, and the result metadata SHALL identify each as diagnostic.

- `pca` SHALL return the scores of the leading principal components of the standardized, concatenated selected blocks. It uses no labels. The number of components SHALL be a required, explicit integration parameter.
- `pls_within_group` SHALL fit a stage-supervised PLS separately within each group, at a declared number of components per group (default `n_stages − 1`). It SHALL orthonormalize the groups' stacked weight vectors into one basis and return the pooled standardized samples projected onto that basis. The basis SHALL be fitted once per evaluation, from the observed data and labels, and SHALL NOT be refitted inside RRPP.

The result metadata SHALL record the method, the retained dimensionality, and, for `pls_within_group`, the per-group component count and the rank of the orthonormalized basis. Orientation attribution SHALL remain restricted to `pls` and SHALL reject either diagnostic method with a clear validation error. The forced-component override SHALL remain specific to `pls`.

#### Scenario: PCA integration is label-free

- **WHEN** the caller selects `pca` with a component count
- **THEN** the latent matrix equals the leading principal-component scores of the standardized selected blocks, and permuting the group or stage labels leaves it unchanged

#### Scenario: PCA without a component count is rejected

- **WHEN** the caller selects `pca` without a component count, or with one outside the feasible range
- **THEN** the harness raises a clear validation error

#### Scenario: Within-group PLS spans both groups' stage subspaces

- **WHEN** the caller selects `pls_within_group` on a dataset whose two groups progress along orthogonal feature directions
- **THEN** the returned basis has rank up to twice the per-group component count, and both groups' stage-mean displacements lie in its span within numerical tolerance

#### Scenario: Within-group PLS is infeasible

- **WHEN** a group has fewer than two samples in any stage
- **THEN** the harness raises a clear validation error naming the group and stage

#### Scenario: Diagnostic methods are labelled

- **WHEN** an evaluation uses `pca` or `pls_within_group`
- **THEN** the result metadata identifies the integration as diagnostic and records its retained dimensionality

#### Scenario: Attribution rejects diagnostic methods

- **WHEN** attribution is enabled with `pca` or `pls_within_group`
- **THEN** the harness raises a clear validation error

### Requirement: Harness records the latent checkpoint for the diagnostic latent spaces

Realized-geometry diagnostics SHALL include the latent checkpoint (the trajectory geometry measured on the evaluated latent matrix) for `pls`, `pca`, and `pls_within_group`. For `pls`, the checkpoint's scope name and values SHALL be unchanged. `snf` and `concat` evaluations SHALL be unchanged. `concat`'s geometry is already reported by the observed standardized checkpoint, and adding a checkpoint to `snf` would change existing records.

#### Scenario: Diagnostic latent spaces report the latent checkpoint

- **WHEN** an evaluation uses `pca` or `pls_within_group`
- **THEN** the realized geometry includes a latent checkpoint computed on that evaluation's latent matrix, named for its method

#### Scenario: PLS latent checkpoint is unchanged

- **WHEN** an evaluation uses `pls`
- **THEN** the latent checkpoint's scope name and values equal those produced before this requirement existed

### Requirement: Harness measures orientation with net displacement by default

The evaluation parameters SHALL accept an orientation estimator selection, either the net-displacement estimator (the default) or the PC1 estimator. The selection SHALL be used for the observed `angle`, for every RRPP permutation, and for every realized-geometry checkpoint `angle`. When the PC1 estimator is selected, the selection SHALL be omitted from the persisted record, the parameter signature and the cell identity, so that PC1 evaluation results, signatures and cell ids are identical to those produced before the selection existed. When the selection holds its default, or names the net-displacement estimator explicitly, the persisted record SHALL state the estimator and the parameter signature and cell identity SHALL include it. An unknown estimator SHALL be rejected with a clear validation error.

#### Scenario: Net displacement is the default

- **WHEN** the evaluation parameters omit the estimator
- **THEN** the observed `angle`, the RRPP null and the realized-geometry `angle` values use the net-displacement estimator, and the record states `net_displacement` as the selected estimator

#### Scenario: PC1 estimator is invisible

- **WHEN** the evaluation parameters name the PC1 estimator
- **THEN** the evaluation result, its persisted record, its parameter signature and its cell id equal those produced before the selection existed

#### Scenario: Net-displacement estimator is recorded

- **WHEN** the evaluation parameters select the net-displacement estimator, by default or by name
- **THEN** the observed `angle`, the RRPP null and the realized-geometry `angle` values use it, and the record states the selected estimator

#### Scenario: Unknown estimator is rejected

- **WHEN** the evaluation parameters name an unknown estimator
- **THEN** the harness raises a clear validation error listing the allowed estimators

### Requirement: PLS rank selection respects a floor of 2 by default

The evaluation parameters SHALL accept an integer rank floor for PLS integration, default 2. Cross-validation SHALL consider only candidates at or above the floor and below the maximum, so every outer-fold selection, the modal rank and the pooled fit are at or above it. A floor below 1, or a non-integer, SHALL be rejected with a clear error. The floor is part of every evaluation's identity and applies to PLS only.

#### Scenario: Default floor lifts a rank-1 selection

- **WHEN** a PLS evaluation runs without naming a floor on an input where the unfloored rule selected rank 1
- **THEN** the selected rank is at least 2 and the latent space has at least 2 components

#### Scenario: Named floor is honoured

- **WHEN** a PLS evaluation names a floor of 3
- **THEN** every cross-validation selection and the final pooled fit have at least 3 components

#### Scenario: Invalid floor is rejected

- **WHEN** the evaluation parameters name a floor of 0 or a non-integer floor
- **THEN** the harness raises a clear validation error

### Requirement: Maximum rank is clamped above the floor

The cross-validation maximum SHALL be clamped to the range from the floor plus one to the feasible maximum (the smaller of the feature and sample counts). When the feasible maximum is below the floor plus one, the harness SHALL raise a clear error naming the floor and the counts rather than fit at a lower rank.

#### Scenario: Maximum at the floor is raised

- **WHEN** a PLS evaluation names a maximum of 2 at the default floor on a sample with at least 3 features and samples
- **THEN** the cross-validation considers rank 2 only, the recorded maximum is 3, and the selected rank is 2

#### Scenario: Sample too small for the floor

- **WHEN** a PLS evaluation runs at the default floor on a sample with 2 features
- **THEN** the harness raises a clear validation error naming the floor and the feature count

### Requirement: Floor of 1 is invisible; a floor above 1 is recorded

When the floor is 1, it SHALL be omitted from the persisted record, the parameter signature and the cell identity, so that floor-1 evaluation results, signatures and cell ids equal those produced before the floor existed. When the floor is above 1, by default or by name, the record's runtime metadata and the latent integration parameters SHALL state it.

#### Scenario: Floor of 1 reproduces the pre-floor evaluation

- **WHEN** the evaluation parameters name a floor of 1
- **THEN** the evaluation result, its persisted record, its parameter signature and its cell id equal those produced before the floor existed

#### Scenario: Default floor is recorded

- **WHEN** a PLS evaluation runs at the default floor
- **THEN** the runtime metadata states the floor and the latent integration parameters state it beside the selected rank

### Requirement: Harness measures span capture on request

The evaluation parameters SHALL accept `span_capture`, absent by default, with `false` normalized to absent. When it is on and the integration is a linear method (`pls`, `pca`, `pls_residual_pca` or `concat`), the harness SHALL compute the `population_latent` checkpoint and the `span_capture` shares from the same fitted basis and preprocessor the trajectory test uses, without refitting. While it is absent, cell identifiers, signatures, records and runtime metadata SHALL be unchanged.

#### Scenario: Span capture reuses the evaluation's fit
- **WHEN** an evaluation runs with `span_capture` on and PLS integration
- **THEN** the PLS model is fitted once, and the realized geometry carries `population_latent` and `span_capture` computed from that model

#### Scenario: Span capture covers the linear diagnostic spaces
- **WHEN** an evaluation runs with `span_capture` on and `pca`, `pls_residual_pca` or `concat` integration
- **THEN** the realized geometry carries `population_latent` and `span_capture` computed from that evaluation's own basis, and for `concat` every share is 1 and the rank equals the number of features

#### Scenario: PLS span capture records are unchanged
- **WHEN** an evaluation runs with `span_capture` on and PLS integration
- **THEN** its `population_latent` checkpoint and `span_capture` block are byte-identical to those produced before span capture covered other methods

#### Scenario: Span capture is inert for other integration methods
- **WHEN** an evaluation runs with `span_capture` on and `snf` or `pls_within_group` integration
- **THEN** the evaluation completes and records no `population_latent` checkpoint and no `span_capture` block

#### Scenario: Absent span capture changes nothing
- **WHEN** an evaluation runs without `span_capture`, or with it set to `false`
- **THEN** its result and cell identity are byte-identical to those produced before this change

#### Scenario: Statistics do not depend on span capture
- **WHEN** the same dataset is evaluated with and without `span_capture`
- **THEN** the latent matrix, the selected rank, the statistics and the p-values are identical

### Requirement: Harness refuses a stage-span dataset measured on other blocks

The harness SHALL refuse, before any fitting, to evaluate a stage-span dataset whose recorded span layers differ from the blocks the evaluation measures.

#### Scenario: Span layers and measured blocks differ
- **WHEN** a dataset generated with `construction: "stage_span"` and span layers methylation + expression is evaluated on all three blocks
- **THEN** the harness raises an evaluation error naming both block sets, and no PLS model is fitted

#### Scenario: Matching blocks are accepted
- **WHEN** the evaluation measures exactly the dataset's span layers
- **THEN** the evaluation proceeds

### Requirement: Harness supports the diagnostic stage-plus-residual integration

The harness SHALL accept a diagnostic integration method `pls_residual_pca`. It SHALL return the production stage-PLS scores, at the CV rank with the rank floor, followed by r principal-component scores of the stage-residualized features, taken orthogonal to the PLS span. r is an integration parameter, default `n_stages − 1`. The basis SHALL use the features and stage labels, never the group labels, and SHALL be fitted once, not inside RRPP.

#### Scenario: The leading coordinates are the production scores
- **WHEN** a dataset is evaluated with `pls_residual_pca`
- **THEN** its first columns equal the latent matrix of a `pls` evaluation of the same dataset and parameters, and its last r columns are the residual component scores

#### Scenario: The basis does not use the group labels
- **WHEN** the group labels are permuted within each stage and the dataset is integrated again
- **THEN** the `pls_residual_pca` latent matrix is unchanged

#### Scenario: Residual components are orthogonal to the PLS span
- **WHEN** a dataset is integrated with `pls_residual_pca`
- **THEN** each residual component direction is orthogonal to the span of the PLS rotations in the model's input coordinates, to numerical tolerance

#### Scenario: Infeasible residual rank is rejected
- **WHEN** r is below 1 or above the number of features less the PLS rank
- **THEN** the harness raises a clear validation error naming r and the feasible range

### Requirement: The stage-plus-residual integration is labelled diagnostic

`pls_residual_pca` SHALL be labelled diagnostic in its metadata, which SHALL record the PLS rank, r and the retained dimensionality. Realized geometry SHALL include its latent checkpoint. Attribution and the forced-component override SHALL reject it.

#### Scenario: Metadata and checkpoint
- **WHEN** an evaluation uses `pls_residual_pca`
- **THEN** the metadata identifies the integration as diagnostic and records the PLS rank, r and the retained dimensionality, and realized geometry carries a latent checkpoint named for the method

#### Scenario: Attribution is rejected
- **WHEN** attribution or a forced component count is requested with `pls_residual_pca`
- **THEN** the harness raises a clear validation error

### Requirement: Harness accepts a prebuilt stage-PLS fit for the residual integration

The harness SHALL let a caller build the `pls_residual_pca` integration from an existing `pls` integration of the same dataset and parameters, without refitting the PLS model. The result SHALL equal the integration fitted from scratch.

#### Scenario: Reusing the PLS fit changes nothing
- **WHEN** `pls_residual_pca` is built from a prebuilt `pls` integration and, separately, fitted from scratch on the same dataset and parameters
- **THEN** the two latent matrices and metadata are identical, and the reuse path fits no PLS model

### Requirement: Harness measures PLS in orthonormal span coordinates on request

The evaluation parameters SHALL accept `latent_coordinates`. Under `orthonormal`, the PLS latent matrix SHALL be the final model's input coordinates projected on an orthonormal basis of its rotations' span, and every trajectory statistic, RRPP draw, realized-geometry checkpoint and configuration spectrum SHALL be computed in those coordinates. The model, the cross-validation and the selected or forced rank SHALL be those of the score measurement.

#### Scenario: Same span, orthonormal coordinates
- **WHEN** a dataset is evaluated with `latent_coordinates` set to `orthonormal`
- **THEN** the latent matrix spans the same subspace as the score matrix of the same fit, and equals the input coordinates times the PLS weights up to column signs, to 1e-8

#### Scenario: Fit and rank are unchanged
- **WHEN** the same dataset is evaluated in score and in orthonormal coordinates
- **THEN** the selected rank, the cross-validation table and the held-out AUROC are identical

#### Scenario: Statistics do not depend on the orthonormal basis
- **WHEN** the orthonormal latent matrix is multiplied by any orthogonal matrix and measured with the same design and RRPP seed
- **THEN** `delta`, `angle` and `shape` agree with the unrotated measurement to 1e-10, and their p-values are equal

#### Scenario: Population checkpoint in the same coordinates
- **WHEN** span capture is on under orthonormal coordinates
- **THEN** the `population_latent` checkpoint is measured in the orthonormal coordinates
- **AND** the `span_capture` shares are identical to those of the score measurement

#### Scenario: The setting is recorded
- **WHEN** an evaluation runs under orthonormal coordinates
- **THEN** its integration parameters and runtime metadata state `latent_coordinates: "orthonormal"`

### Requirement: The coordinate setting is normalized and validated before fitting

`latent_coordinates` SHALL accept absent (the score coordinates), `scores` (normalized to absent) and `orthonormal`. The harness SHALL reject the setting before any fit when:
- the value is unknown;
- the integration method is not `pls`;
- attribution is requested together with `orthonormal`.

Each rejection SHALL name the offending value.

#### Scenario: Scores normalize to absent
- **WHEN** the evaluation parameters name `scores`
- **THEN** the parameters equal those with the setting absent

#### Scenario: Unknown value
- **WHEN** the evaluation parameters name an unknown coordinate system
- **THEN** the harness raises a validation error naming the value and the allowed values

#### Scenario: Other integration methods
- **WHEN** `orthonormal` is set with `pca`, `concat`, `snf`, `pls_within_group` or `pls_residual_pca` integration
- **THEN** the harness raises a validation error naming the method before any preprocessing or fit

#### Scenario: Attribution
- **WHEN** `orthonormal` is set and attribution is requested
- **THEN** the harness raises a validation error before any fit, stating that attribution is defined only in score coordinates

### Requirement: An absent coordinate setting is invisible

While `latent_coordinates` is absent, it SHALL be omitted from the persisted record, the runtime metadata, the parameter signature and the cell identity. Every evaluation result, signature and cell id SHALL be identical to those produced before the setting existed.

#### Scenario: Committed identities are unchanged
- **WHEN** every committed study configuration is enumerated after this change
- **THEN** its cell ids, signatures and matched seeds equal the committed identity fixtures

#### Scenario: Results are unchanged
- **WHEN** a dataset is evaluated without the setting, or with `scores`
- **THEN** the result is byte-identical to the result produced before this change

### Requirement: Harness builds the orthonormal integration from a prebuilt PLS fit

The harness SHALL let a caller build the orthonormal-coordinate integration from an existing score-coordinate `pls` integration of the same dataset and parameters, without refitting. The result SHALL equal the integration fitted from scratch with the setting. An evaluation SHALL refuse a prebuilt integration whose coordinates differ from its own setting.

#### Scenario: Reusing the fit changes nothing
- **WHEN** the orthonormal integration is built from a prebuilt `pls` integration and, separately, fitted from scratch on the same dataset and parameters
- **THEN** the two latent matrices and metadata are identical, and the reuse path fits no PLS model

#### Scenario: Mismatched prebuilt integration
- **WHEN** an evaluation set to one coordinate system receives a prebuilt integration in the other
- **THEN** the harness raises a validation error naming both

#### Scenario: Not a PLS fit
- **WHEN** the prebuilt integration is not a `pls` integration with fitted artifacts, or already carries orthonormal coordinates
- **THEN** the harness raises a validation error naming what it received
