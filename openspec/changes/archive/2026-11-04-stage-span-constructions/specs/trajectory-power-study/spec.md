# Spec Delta

## ADDED Requirements

### Requirement: Study configurations accept the stage-span construction and span capture

The study configuration SHALL accept `generator.construction`, `generator.span_layers` and `evaluation.span_capture` in its baseline blocks. `span_layers` may also be a design axis. The loader SHALL reject a stage-span configuration whose span layers differ from the measured blocks at any column, where absent means all three blocks, and SHALL reject `generator.construction` as a design axis. Configurations that use none of the keys SHALL enumerate exactly as before.

#### Scenario: Stage-span configuration loads and enumerates
- **WHEN** a configuration sets `generator.construction` to `stage_span`, `generator.span_layers` to methylation + expression, `evaluation.integration_params.layers` to the same blocks and `evaluation.span_capture` to `true`
- **THEN** it loads, enumerates, and reloading the same file yields identical cell identifiers, signatures and matched seeds

#### Scenario: Mismatched span layers are rejected at load
- **WHEN** a stage-span configuration's span layers differ from the measured blocks at the baseline or at any design column
- **THEN** loading fails with an error naming the column and both block sets

#### Scenario: Construction is not a design axis
- **WHEN** a design grid declares `generator.construction`
- **THEN** loading fails with an error stating that constructions are compared as sibling profiles in one seed family

#### Scenario: Committed configurations are unchanged
- **WHEN** any configuration committed before this change is enumerated
- **THEN** its cell identifiers, parameter signatures and matched seeds are identical to those enumerated before this change

### Requirement: Stage-span profiles share the surgery anchor

The zero-effect anchor of a stage-span profile SHALL carry no construction and no span layers. It is therefore the same cell, with the same identifier and the same dataset, as the anchor of a surgery profile with otherwise equal parameters in the same seed family.

#### Scenario: Sibling anchors have one identity
- **WHEN** a stage-span profile and its surgery sibling are enumerated in one seed family
- **THEN** at every design column their anchor cells have identical cell identifiers, signatures and matched seeds

### Requirement: Stage-span cells are exempt from surgery headroom

Enumeration SHALL apply the surgery headroom check only to surgery cells. Stage-span cells draw from no feature pool, so no headroom applies, and realized-surgery summaries SHALL report no surgery for them.

#### Scenario: Large stage-span effects enumerate
- **WHEN** a stage-span profile requests `orientation`, `shape` or `translation` at effects a surgery's pool could not supply
- **THEN** enumeration completes without a headroom failure, subject only to the stage-span construction's own validation

### Requirement: Reporting summarizes span capture

When records carry span capture, the report SHALL write `span_capture.csv`. It has one row per cell: design column, mode, effect, construction, the number of available replicates, the median and 10th and 90th percentiles of the `group_difference` and `pooled_stage` shares, the median rank and the median `population_latent` target statistic. Without span capture the file SHALL NOT be written and every other report output SHALL be unchanged.

#### Scenario: Capture table is written for span-capture records
- **WHEN** a study whose evaluation sets `span_capture` is reported
- **THEN** `report/span_capture.csv` exists with one row per cell and the columns above, and a share that is null in every replicate reports zero available replicates and empty summaries

#### Scenario: Reports are unchanged without span capture
- **WHEN** a study without `span_capture` is reported
- **THEN** no `span_capture.csv` is written and every other report output is byte-identical to the output produced before this change
