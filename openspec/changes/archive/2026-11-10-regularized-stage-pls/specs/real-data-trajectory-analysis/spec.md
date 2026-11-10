# Spec Delta

## ADDED Requirements

### Requirement: Real-data analysis offers sparse stage-PLS

The analysis parameters SHALL accept a sparsity setting (scope `joint` | `block`; a ladder of keep fractions in (0, 1], distinct, defaulting to 1.0, 0.5, 0.2, 0.1, 0.05). When set, the latent space is the sparse stage-PLS-DA cross-validated over the rank and the ladder jointly, with block membership from the declared blocks, and attribution, the split interaction test and the interaction scan run unchanged on it. A malformed setting SHALL be rejected before any fitting, and the setting SHALL be omitted from the resolved parameters while absent.

#### Scenario: Sparse analysis
- **WHEN** an analysis declares two blocks and requests scope `block`
- **THEN** every component of the fitted space loads on the selected fraction of each block, and the written loadings show the dropped features at zero

#### Scenario: Absent sparsity
- **WHEN** no sparsity is requested
- **THEN** every output equals the pre-change output, and the resolved parameters carry no sparsity key

### Requirement: Real-data analysis records the sparsity setting

When set, the record SHALL carry `latent.integration_params.sparsity` and a `sparsity` block after `latent` with the scope, the ladder, the selected keep fraction and, per component, the support size per block. Neither SHALL appear while absent.

#### Scenario: Record contents
- **WHEN** a sparse analysis is recorded
- **THEN** the `sparsity` block names the scope, lists the ladder, the selected fraction and each component's support per block

### Requirement: motco analyze exposes sparse stage-PLS

`motco analyze` SHALL accept `--sparsity {joint,block}` and an optional `--sparsity-ladder` (comma-separated fractions), print the selected fraction beside the selected rank when given, and reject `--sparsity-ladder` without `--sparsity`.

#### Scenario: CLI round trip
- **WHEN** `motco analyze --sparsity block --sparsity-ladder 1.0,0.2` runs
- **THEN** the written record's `sparsity.scope` is `block` and its ladder is 1.0, 0.2
