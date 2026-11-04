# Spec Delta

## ADDED Requirements

### Requirement: Harness measures span capture on request

The evaluation parameters SHALL accept `span_capture`, absent by default, with `false` normalized to absent. When it is on and the integration is PLS, the harness SHALL compute the `population_latent` checkpoint and the `span_capture` shares from the same fitted PLS model and preprocessor the trajectory test uses, without refitting. While it is absent, cell identifiers, signatures, records and runtime metadata SHALL be unchanged.

#### Scenario: Span capture reuses the evaluation's fit
- **WHEN** an evaluation runs with `span_capture` on and PLS integration
- **THEN** the PLS model is fitted once, and the realized geometry carries `population_latent` and `span_capture` computed from that model

#### Scenario: Span capture is inert for other integration methods
- **WHEN** an evaluation runs with `span_capture` on and an integration method other than PLS
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
