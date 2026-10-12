## ADDED Requirements

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

### Requirement: Harness accepts an orientation estimator selection

The evaluation parameters SHALL accept an orientation estimator selection, either the PC1 estimator (default) or the net-displacement estimator. The selection SHALL be used for the observed `angle`, for every RRPP permutation, and for every realized-geometry checkpoint `angle`. When the selection holds its default, it SHALL be omitted from the persisted record, the parameter signature and the cell identity. Evaluation results, signatures and cell ids SHALL then be identical to those produced before the selection existed. An unknown estimator SHALL be rejected with a clear validation error.

#### Scenario: Default estimator is invisible

- **WHEN** the evaluation parameters omit the estimator or name the PC1 estimator
- **THEN** the evaluation result, its parameter signature and its cell id equal those produced before the selection existed

#### Scenario: Net-displacement estimator is recorded

- **WHEN** the evaluation parameters select the net-displacement estimator
- **THEN** the observed `angle`, the RRPP null and the realized-geometry `angle` values use it, and the record states the selected estimator

#### Scenario: Unknown estimator is rejected

- **WHEN** the evaluation parameters name an unknown estimator
- **THEN** the harness raises a clear validation error listing the allowed estimators

### Requirement: Harness records the latent checkpoint for the diagnostic latent spaces

Realized-geometry diagnostics SHALL include the latent checkpoint (the trajectory geometry measured on the evaluated latent matrix) for `pls`, `pca`, and `pls_within_group`. For `pls`, the checkpoint's scope name and values SHALL be unchanged. `snf` and `concat` evaluations SHALL be unchanged. `concat`'s geometry is already reported by the observed standardized checkpoint, and adding a checkpoint to `snf` would change existing records.

#### Scenario: Diagnostic latent spaces report the latent checkpoint

- **WHEN** an evaluation uses `pca` or `pls_within_group`
- **THEN** the realized geometry includes a latent checkpoint computed on that evaluation's latent matrix, named for its method

#### Scenario: PLS latent checkpoint is unchanged

- **WHEN** an evaluation uses `pls`
- **THEN** the latent checkpoint's scope name and values equal those produced before this requirement existed
