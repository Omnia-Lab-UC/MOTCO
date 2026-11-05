# Spec Delta

## ADDED Requirements

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
