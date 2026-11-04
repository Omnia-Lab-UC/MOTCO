# stage-span-construction-validation Specification

## Purpose
Measures, on matched null datasets, how the stage-span construction behaves next to the feature-set surgeries: how much of each plant the fitted PLS basis sees, whether the stage-span modes stay specific to their target statistic, and whether the size control holds. An advisory reading decides whether a paper-grade follow-up is proposed.

## Requirements

### Requirement: Stage-span pilot profiles

The study SHALL provide four version-controlled pilot configurations: {surgery, stage-span} × {magnitude, geometry}. Each runs 100 replicates at 199 permutations in the `design-guidance` seed family at base seed 1028, with span capture on. Each names `net_displacement`, a rank floor of 2 and `joint` magnitude. Each profile has three design columns at three stages measured on methylation + expression: n = 200 (the baseline) and n = 80 balanced, and the SEA-AD table `[[11, 10, 28], [9, 10, 12]]`.

#### Scenario: Profiles enumerate their columns and cells
- **WHEN** each pilot profile is loaded and enumerated
- **THEN** it yields three design columns, each with one zero-effect anchor and one cell per (trajectory mode, nonzero effect)
- **AND** the magnitude profiles' only mode is `magnitude` at effects 0.02, 0.05 and 0.1
- **AND** the geometry profiles' modes are `orientation`, `shape` and `translation` at effects 0.25, 0.5 and 1.0

#### Scenario: Stage-span profiles declare their construction and span
- **WHEN** a stage-span profile is loaded
- **THEN** its generator declares `construction: "stage_span"` and `span_layers` equal to the measured blocks at every column
- **AND** its surgery sibling is identical except that it leaves `construction` and `span_layers` unset

#### Scenario: Profiles are post-floor pinned configurations
- **WHEN** the committed post-floor fixture is checked
- **THEN** the four pilot profiles are listed in it with their cell identifiers, signatures and matched seeds

### Requirement: Construction siblings share their anchors

Within one design column, the zero-effect anchor of a stage-span profile SHALL be the same cell, and evaluate the same dataset, as the anchor of its surgery sibling. Anchors evaluated by several sibling profiles SHALL never be pooled.

#### Scenario: Anchors reproduce across siblings
- **WHEN** the four pilot profiles have run
- **THEN** at every design column and replicate index, the anchor's statistics and p-values are identical in all four profiles

### Requirement: Surgery siblings reproduce the design-guidance pilot

At the n = 200 and n = 80 columns, every surgery-sibling cell SHALL evaluate the same dataset as the matching cell of the 2026-10-28 design-guidance pilot. Turning on span capture SHALL change none of the statistics or p-values.

#### Scenario: Shared columns reproduce the earlier pilot
- **WHEN** the surgery siblings' records at n = 200 and n = 80 are paired with the design-guidance pilot's records on (column, mode, effect, replicate index)
- **THEN** every pair has identical statistics and p-values, and the comparison reports the number of pairs and of mismatches

### Requirement: Latent retrofit check

The system SHALL provide a descriptive check. It edits group B's latent stage means in a replicate's own fitted PLS basis (size, orientation, shape, translation), pulls each edit back to the features through the loadings and the minimum-norm preimage, and re-measures it twice: through the frozen basis and through a full re-analysis. It issues no decision.

#### Scenario: Frozen basis reproduces the edit
- **WHEN** an edited feature matrix is projected through the replicate's frozen basis
- **THEN** the recovered latent stage means equal the edited ones to within 1e-8, and the check reports the intended and recovered geometry

#### Scenario: Re-analysis is reported beside the intended geometry
- **WHEN** an edited dataset is re-analysed from scratch
- **THEN** the check reports the selected rank, the three statistics and their p-values beside the intended geometry, for each edit and each preimage

#### Scenario: Footprint is compared with the population plant
- **WHEN** an edit is pulled back to the features
- **THEN** the check reports the footprint's effective number of features and its share on the truth's stage-active features
- **AND** it reports the footprint's cosine with the stage-span construction's population plant for the same mode and the same B configuration

### Requirement: Advisory pilot reading

A reading driver SHALL read the four pilot reports and write `stage_span_pilot_reading.json`, verdict `propose_paper_grade` or `revise_construction`. The bound is α + 2·√(α(1 − α)/R), with R the available replicates. The verdict is `propose_paper_grade` only when all three criteria below hold. Otherwise it is `revise_construction`, naming every failing criterion, cell and statistic.

#### Scenario: Criterion R1, the stage-span translation control holds size
- **WHEN** the reading evaluates the stage-span `translation` cells
- **THEN** R1 holds only if, at every column, the `delta`, `angle` and `shape` rejection rates are at or below the bound

#### Scenario: Criterion R2, the stage-span plant is captured
- **WHEN** the reading evaluates every stage-span `magnitude`, `orientation` and `shape` cell
- **THEN** R2 holds only if the median group-difference capture share is at least 0.90 in every cell

#### Scenario: Criterion R3, the stage-span modes are specific
- **WHEN** the reading evaluates every stage-span `magnitude`, `orientation` and `shape` cell
- **THEN** R3 holds only if the rejection rate of each statistic other than the mode's target is at or below the bound in every cell
- **AND** the targets are `magnitude`→`delta`, `orientation`→`angle` and `shape`→`shape`

#### Scenario: Descriptive comparisons accompany the verdict
- **WHEN** the reading is written
- **THEN** it reports, per column, mode and effect, both constructions' target power and median capture share
- **AND** it reports the stage-span minus surgery `magnitude` `delta` power difference with its pooled standard error
- **AND** it reports the anchor size of every column as read, without judging it again

### Requirement: Dated stage-span pilot report

The run SHALL be recorded in a dated report under `docs/reports/`, with aggregates only (reports, the reading, the retrofit check's tables, reproduction summaries, provenance) under `results/stage-span-pilot-<date>/`. The report SHALL state the verdict, the surgery capture shares per mode, and the effect grid the paper-grade follow-up would use.

#### Scenario: Report states the reading and its evidence
- **WHEN** the dated report is written
- **THEN** it quotes the reading's verdict and criteria, the capture shares of every mode under both constructions, the anchor and surgery reproduction results, and the retrofit check's recovery
- **AND** it names the paper-grade follow-up and the latent plasmode as separate changes
