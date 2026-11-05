# Spec Delta

## MODIFIED Requirements

### Requirement: Study configurations accept the stage-span construction and span capture

The study configuration SHALL accept `generator.construction`, `generator.span_layers`, `generator.span_scale` and `evaluation.span_capture` in its baseline blocks. `span_layers` may also be a design axis. The loader SHALL reject a stage-span configuration whose span layers differ from the measured blocks at any column, where absent means all three blocks, and SHALL reject `generator.construction` and `generator.span_scale` as design axes. Configurations that use none of the keys SHALL enumerate exactly as before.

#### Scenario: Stage-span configuration loads and enumerates
- **WHEN** a configuration sets `generator.construction` to `stage_span`, `generator.span_layers` to methylation + expression, `evaluation.integration_params.layers` to the same blocks and `evaluation.span_capture` to `true`
- **THEN** it loads, enumerates, and reloading the same file yields identical cell identifiers, signatures and matched seeds

#### Scenario: Planted-scale configuration loads and enumerates
- **WHEN** a stage-span configuration also sets `generator.span_scale` to `planted`
- **THEN** it loads and enumerates, its stage-span cells carry `span_scale: "planted"`, and their cell identifiers and signatures differ from the same configuration without the key
- **AND** a configuration that sets `generator.span_scale` to `null` enumerates identically to one that omits the key

#### Scenario: Mismatched span layers are rejected at load
- **WHEN** a stage-span configuration's span layers differ from the measured blocks at the baseline or at any design column
- **THEN** loading fails with an error naming the column and both block sets

#### Scenario: Construction is not a design axis
- **WHEN** a design grid declares `generator.construction`
- **THEN** loading fails with an error stating that constructions are compared as sibling profiles in one seed family

#### Scenario: Stage-span scale is not a design axis
- **WHEN** a design grid or an OFAT axis declares `generator.span_scale`
- **THEN** loading fails with an error stating that stage-span scales are compared as sibling profiles in one seed family

#### Scenario: Committed configurations are unchanged
- **WHEN** any configuration committed before this change is enumerated
- **THEN** its cell identifiers, parameter signatures and matched seeds are identical to those enumerated before this change

### Requirement: Stage-span profiles share the surgery anchor

The zero-effect anchor of a stage-span profile SHALL carry no construction, no span layers and no span scale. It is therefore the same cell, with the same identifier and the same dataset, as the anchor of a surgery profile with otherwise equal parameters in the same seed family.

#### Scenario: Sibling anchors have one identity
- **WHEN** a stage-span profile and its surgery sibling are enumerated in one seed family
- **THEN** at every design column their anchor cells have identical cell identifiers, signatures and matched seeds

#### Scenario: Planted-scale anchors are the surgery anchors
- **WHEN** a planted-scale stage-span profile and a surgery profile with otherwise equal parameters are enumerated in one seed family
- **THEN** at every design column they share their anchor cell, which carries no `span_scale`
