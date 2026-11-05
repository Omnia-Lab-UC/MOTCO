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

### Requirement: Paper-grade stage-span profiles

The study SHALL provide three post-floor paper-grade profiles at 500 replicates × 999 permutations in the `design-guidance` seed family (base seed 1028), with span capture on, naming `net_displacement`, a rank floor of 2 and `joint` magnitude. Two are planted-scale stage-span profiles at five design columns; the third is a surgery profile at the two columns no earlier paper-grade run covers.

#### Scenario: Stage-span geometry profile
- **WHEN** the stage-span geometry profile is loaded and enumerated
- **THEN** its generator declares `construction: "stage_span"` and `span_scale: "planted"`, with span layers equal to the measured blocks at every column
- **AND** its modes are `orientation`, `shape` and `translation` at effects 0.1, 0.25 and 0.5, with one zero-effect anchor per column

#### Scenario: Stage-span magnitude profile
- **WHEN** the stage-span magnitude profile is loaded and enumerated
- **THEN** it is identical to the geometry profile except that its only mode is `magnitude`, at effects 0.02, 0.05 and 0.1

#### Scenario: The five stage-span columns
- **WHEN** either stage-span profile is enumerated
- **THEN** its design columns are: n = 200 (the baseline) and n = 80 at three stages on methylation + expression; the SEA-AD table `[[11, 10, 28], [9, 10, 12]]` at three stages on methylation + expression; n = 200 at three stages on all three blocks; and n = 200 at four stages on methylation + expression

#### Scenario: Surgery geometry profile
- **WHEN** the surgery geometry profile is loaded and enumerated
- **THEN** it leaves `construction`, `span_layers` and `span_scale` unset, runs `orientation`, `shape` and `translation` at effects 0.25, 0.5 and 1.0, and has two design columns: the SEA-AD table and n = 200 at four stages, both on methylation + expression

#### Scenario: Profiles are post-floor pinned configurations
- **WHEN** the committed post-floor fixture is checked
- **THEN** the three profiles are listed in it with their cell identifiers, signatures and matched seeds

### Requirement: Paper-grade siblings share anchors and reproduce earlier runs

Within a design column, the paper-grade profiles SHALL evaluate one anchor dataset (same generator parameters and seed), even where a column is one profile's baseline and another's design point. Paper-grade cells SHALL draw the same generator seeds as the earlier runs' matching cells, so the run can be checked against the design-guidance slice, the stage-span pilot and the design-guidance pilot.

#### Scenario: Sibling anchors are identical
- **WHEN** the paper-grade profiles have run
- **THEN** at every shared column and replicate index the anchor's statistics and p-values are identical across the profiles that run that column

#### Scenario: Anchors reproduce the design-guidance slice
- **WHEN** the anchors at n = 200 and n = 80 on two blocks and at n = 200 on three blocks are paired with the 2026-10-29 slice's anchors on (column, replicate index)
- **THEN** every pair has identical seeds, ranks, statistics and p-values, and the comparison reports the number of pairs and of mismatches

#### Scenario: Surgery cells reproduce the pilots' datasets
- **WHEN** the surgery profile's cells are paired with the 2026-11-04 stage-span pilot's SEA-AD surgery cells and with the 2026-10-28 design-guidance pilot's four-stage cells at replicates 0–99
- **THEN** every pair has identical generator seeds, ranks and statistics, and p-values are not compared because the permutation counts differ

#### Scenario: Planted-scale magnitude reproduces the pilot's stage-span magnitude
- **WHEN** the planted-scale magnitude cells at n = 200, n = 80 and the SEA-AD table are paired with the 2026-11-04 stage-span pilot's magnitude cells at replicates 0–99
- **THEN** every pair has identical generator seeds, ranks and statistics, because scaling about the centroid gives the same plant under either scale

### Requirement: Stage-span leak decomposition check

The system SHALL provide a descriptive check that runs a profile's stage-span cells and, per unit, reports the signed B − A path-length difference at six checkpoints. It SHALL accept overrides of the span scale, modes, effects, replicates and permutations, run in shards, and issue no verdict.

#### Scenario: The six checkpoints are reported
- **WHEN** the check runs one unit
- **THEN** it reports B's and A's path lengths and their difference for the population means divided by the null pooled SD; divided by the exact planted pooled SD; through the replicate's sample scaler; projected on an orthonormal basis of the fitted span; in the fitted PLS score coordinates; and for the sample latent stage means
- **AND** it reports the selected rank, the capture share and, when permutations are run, the three p-values

#### Scenario: Merged output summarizes each cell
- **WHEN** the shards are merged
- **THEN** a per-cell table reports the median of every difference, the differences as a share of A's path length, the median capture and rank, and the rejection rates at α = 0.05 with p < α as the rejection rule

#### Scenario: The null-scale arm reproduces the pilot
- **WHEN** the check runs the paper-grade geometry profile's three pilot columns with `span_scale` overridden to `null` at effects 0.25, 0.5 and 1.0, 100 replicates and 199 permutations
- **THEN** its rejection rates equal the 2026-11-04 stage-span pilot's rates at those cells

### Requirement: Pre-flight gate before the paper-grade array

Before the paper-grade array is submitted, the leak check SHALL run at 100 replicates × 199 permutations on the stage-span geometry profile's five columns: under `planted` at the paper-grade effects and under `null` at 0.25, 0.5 and 1.0, for the anchor, `orientation` and `shape`. The array SHALL be submitted only if the gate holds; otherwise the work stops and is returned to the user.

#### Scenario: The scale route is removed exactly
- **WHEN** the planted arm is read
- **THEN** the gate requires that every unit generated without error and that, in every unit, the planted-scale path-length difference is at most 1e-8 of A's path length

#### Scenario: Off-target rates stay at chance level
- **WHEN** the planted arm is read
- **THEN** the gate requires that the `orientation` and `shape` off-target rejection rates over every column and effect pass the decision's bound rule for the number of tests (at R = 100 and 60 rates: at most 4 above 0.0936, and none above the Bonferroni ceiling)

#### Scenario: The planted plant shrinks the score-coordinate residual
- **WHEN** both arms are read at e = 0.5 for `shape`
- **THEN** the gate requires that, at every column, the planted arm's median score-coordinate B − A path-length difference is smaller in magnitude than the null arm's

#### Scenario: A failing gate stops the run
- **WHEN** either gate condition fails
- **THEN** no paper-grade array is submitted, and the pre-flight tables and the failing cells are reported to the user

### Requirement: Predeclared paper-grade decision

A decision driver SHALL read the three paper-grade reports and write `stage_span_paper_grade_decision.json` (rule version 1). It judges size (P1), capture (P2) and specificity (P3) and returns `size_excess`, `capture_shortfall`, `specificity_excess` or `sensitivity_established`, the first failing criterion in that order, naming every failing criterion, cell and statistic.

#### Scenario: Criterion P1, size holds at every column
- **WHEN** the decision reads the anchor of every column and the stage-span `translation` cells
- **THEN** P1 is judged on their `delta`, `angle` and `shape` rejection rates under the bound rule

#### Scenario: Criterion P2, the planted difference is captured
- **WHEN** the decision reads every stage-span `magnitude`, `orientation` and `shape` cell
- **THEN** P2 holds only if the median group-difference capture share is at least 0.90 in every cell

#### Scenario: Criterion P3, the stage-span modes are specific
- **WHEN** the decision reads every stage-span `magnitude`, `orientation` and `shape` cell
- **THEN** P3 is judged on each statistic other than the mode's target (`magnitude`→`delta`, `orientation`→`angle`, `shape`→`shape`) under the bound rule

#### Scenario: The bound rule accounts for the number of tests
- **WHEN** a criterion judges m rejection rates
- **THEN** it holds only if the number of rates above α + 2·√(α(1 − α)/R) is at most the 95th percentile of a binomial count over m tests, each exceeding with the exact probability that a Binomial(R, α) rate exceeds that bound
- **AND** no rate exceeds the Bonferroni ceiling: the smallest rate that a Binomial(R, α) rate exceeds with probability at most 0.05/m
- **AND** the decision records m, the exceedance probability, the allowed count, the observed count and the ceiling for each criterion

#### Scenario: A missing reading fails its criterion
- **WHEN** a required row, rate or capture share is missing
- **THEN** it is recorded as unavailable and counts as a failure of its criterion

### Requirement: Paper-grade decision reports sensitivity descriptively

Beside its verdict, the decision SHALL report target power, the smallest informative effect, capture and realized geometry for every stage-span cell, and the visibility contrast with the surgeries. None of these enters the verdict.

#### Scenario: Target power and the informative effect
- **WHEN** the decision is written
- **THEN** it reports each stage-span cell's target rejection rate with its standard error
- **AND** for every column and mode, `informative_from`: the smallest effect whose target rate is at least 0.80, or null

#### Scenario: Capture and realized geometry
- **WHEN** the decision is written
- **THEN** it reports each stage-span cell's capture quantiles (q10, median, q90), median rank, and the median `population_latent` target statistic beside the intended one

#### Scenario: Visibility contrast with the surgeries
- **WHEN** the decision is written
- **THEN** it reports the surgery cells' target power and capture at the SEA-AD and four-stage columns from this run
- **AND** at the other three columns, the slice's surgery target power and the stage-span pilot's surgery capture where the pilot measured it, labelled with their source and grade, or unavailable

### Requirement: Dated stage-span paper-grade report

The run SHALL be recorded in a dated report under `docs/reports/`, with aggregates only under `results/stage-span-paper-grade-<date>/`: the three reports, the decision, the pre-flight leak tables, the reproduction summaries, provenance and accounting.

#### Scenario: Report states the decision and its evidence
- **WHEN** the dated report is written
- **THEN** it quotes the verdict and each criterion's counts, the target power and `informative_from` per column, the capture shares, and the visibility contrast
- **AND** it reports the leak decomposition under both scales, including the score-coordinate residual at each column, and the reproduction checks
