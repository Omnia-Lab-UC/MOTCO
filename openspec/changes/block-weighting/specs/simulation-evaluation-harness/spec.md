# Spec Delta

## ADDED Requirements

### Requirement: Harness weights blocks on request

The harness SHALL accept a block weighting kind in the integration parameters (`block_weighting`: absent, `equal_variance` or `mfa`). When set, the shared preprocessing weights every measured block by the kind's scalar before any integration method runs, the PLS space is fitted in the weighted metric, and the integration metadata records the kind. An unknown kind SHALL be rejected before any preprocessing, naming the value and the allowed kinds. Absent, nothing changes and no metadata entry is written.

#### Scenario: Weighted evaluation
- **WHEN** an evaluation requests `equal_variance` weighting on methylation and expression
- **THEN** both transformed blocks have unit total variance, the latent space is fitted on them, and the integration metadata records `block_weighting: equal_variance`

#### Scenario: Every method sees the weighted blocks
- **WHEN** `concat`, `pca` or `pls` integration runs with a weighting
- **THEN** each measures the weighted blocks, and the span-capture inputs and population checkpoints use the same weights

#### Scenario: Unknown weighting
- **WHEN** an evaluation requests an unknown weighting kind
- **THEN** it is rejected before any preprocessing, naming the value and the allowed kinds

#### Scenario: Absent weighting is byte-identical
- **WHEN** an evaluation runs without the key
- **THEN** its record, cell identifier and signature equal those from before the key existed
