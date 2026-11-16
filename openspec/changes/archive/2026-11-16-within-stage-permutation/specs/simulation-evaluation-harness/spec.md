# Spec Delta

## ADDED Requirements

### Requirement: Harness draws the permutation null within stage on request

The evaluation parameters SHALL accept a permutation strata setting (`permutation_strata`: absent or `stage`). When it is `stage`, every permutation draw SHALL exchange the reduced model's residuals only among samples of the same stage, with the strata the stage's position in the design's stage order, as the RRPP primitive defines stratified draws. The integration, the selected or forced rank, the measurement coordinates and the observed statistics SHALL be those of the evaluation without the setting. The setting SHALL apply with every integration method and every measurement coordinate system. Any other value SHALL be rejected before any fit, naming the value and the allowed values.

#### Scenario: Only the null moves

- **WHEN** one dataset is evaluated with the same parameters and seed, once without the setting and once with `stage`
- **THEN** the latent matrix, the latent metadata and the observed `delta`, `angle` and `shape` are identical
- **AND** the p-values and the null summary of the second evaluation are those of draws restricted within stage

#### Scenario: Agreement with the primitive

- **WHEN** an evaluation runs with `stage`, and RRPP is called on that evaluation's latent matrix, models and seed with each sample's stage position as its stratum
- **THEN** the evaluation's null summary and p-values are those of that call

#### Scenario: Predictive coordinates under within-stage draws

- **WHEN** an evaluation requests `predictive_orthonormal` coordinates and `stage`
- **THEN** it measures the two predictive dimensions of the unchanged fit, and every draw is restricted within stage

#### Scenario: A prebuilt integration

- **WHEN** two evaluations of one dataset share a prebuilt integration and differ only in the setting
- **THEN** no model is refitted, and their observed statistics are identical

#### Scenario: Unknown value

- **WHEN** the evaluation parameters name the strata `group`
- **THEN** a validation error names the value and the allowed values before any fit

### Requirement: An absent strata setting is invisible

While the permutation strata setting is absent, cell identifiers, parameter signatures, replicate records and reports SHALL be those produced before the setting existed. When it is set, it SHALL enter the cell identifier and the parameter signature, and the runtime metadata SHALL name it.

#### Scenario: Committed configurations are unchanged

- **WHEN** every committed study configuration is enumerated after the change
- **THEN** its cell identifiers, signatures and matched seeds equal the pinned ones

#### Scenario: A set value is hashed and recorded

- **WHEN** two cells differ only in the setting, absent in one and `stage` in the other
- **THEN** their identifiers and signatures differ, and only the second's runtime metadata carries `permutation_strata`
