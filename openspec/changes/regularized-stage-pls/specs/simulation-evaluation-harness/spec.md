# Spec Delta

## ADDED Requirements

### Requirement: The stage-PLS space takes a sparsity budget on request

The harness SHALL accept a sparsity setting in the integration parameters (`sparsity`: absent, or a scope `joint` | `block` with a ladder of keep fractions). When set, the production PLS space is cross-validated over the rank and the ladder jointly, fitted at the selected pair with block membership taken from the measured blocks in concatenation order, and its metadata records the scope, the ladder, the selected fraction and the per-component, per-block support sizes. The setting SHALL compose with the rank floor, latent adjustment, latent coordinates and block weighting. A forced rank SHALL be accepted only with a one-element ladder, else rejected naming both. The setting SHALL be rejected before any preprocessing with an integration method other than `pls` or `pls_residual_pca`, and for a malformed value, naming the offence. Absent, nothing changes and no metadata entry is written.

#### Scenario: Sparse evaluation
- **WHEN** an evaluation requests scope `block` with the default ladder on methylation and expression
- **THEN** the latent space is fitted at the selected rank and fraction, and the metadata records the scope, the ladder, the selected fraction and each component's support size per block

#### Scenario: Forced rank needs a fixed budget
- **WHEN** an evaluation requests a forced rank with a ladder of more than one fraction
- **THEN** it is rejected before any fit, naming the forced rank and the ladder

#### Scenario: Other methods reject the setting
- **WHEN** `pca`, `concat`, `snf` or `pls_within_group` integration is requested with a sparsity setting
- **THEN** the evaluation is rejected before any preprocessing, naming the method

#### Scenario: Absent sparsity is byte-identical
- **WHEN** an evaluation runs without the key
- **THEN** its record, cell identifier and signature equal those from before the key existed
