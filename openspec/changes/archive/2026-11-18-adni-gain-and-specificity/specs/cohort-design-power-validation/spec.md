# Spec Delta

## Purpose

Reads the power and specificity of the trajectory tests at a cohort's design: on plants added to the cohort's exported omic blocks under within-stage-permuted labels, and on stage-span plants in a simulation of that design. A predeclared decision says whether a candidate measurement is reopened at the design and states what the production measurement itself detects there.

## ADDED Requirements

### Requirement: The plants profile builds its units as the cohort anchor does

The plants profile SHALL read the ADNI export in place through the cohort anchor's input contract and SHALL draw each unit's participants, within-stage label permutation and seeds exactly as the cohort real-noise anchor does. No analysis SHALL run under the real labels, and nothing SHALL be written to the export's directory.

#### Scenario: A cohort anchor unit is the anchor run's unit
- **WHEN** the anchor cell of the full-cohort point is built at a replicate
- **THEN** its labels and its permutation seed equal those of the cohort real-noise anchor's unit at that replicate

#### Scenario: A contract violation is named
- **WHEN** the export violates the anchor's input contract
- **THEN** the check exits naming the violation and the file, and no unit is evaluated

#### Scenario: The export is not copied
- **WHEN** a run, a merge or a decision completes
- **THEN** no file of the export has been copied or modified

### Requirement: Two design points, one of them with a held-out reference

The plants profile SHALL read two design points: 750 participants drawn without replacement within each sex × stage cell at the cohort's proportions rounded by largest remainder, and the full cohort of 1,501. At the 750-participant point the participants outside the unit's sample SHALL form its held-out reference.

#### Scenario: Cell tables
- **WHEN** the 750-participant point is built on the cohort's table (Female 237 / 310 / 125, Male 223 / 448 / 158 at CN / MCI / AD)
- **THEN** the sample's table is Female 118 / 155 / 63 and Male 111 / 224 / 79, and the held-out reference holds the other 751 participants

#### Scenario: The reference is disjoint from the sample
- **WHEN** a unit of the 750-participant point is built
- **THEN** no participant is in both its sample and its held-out reference

#### Scenario: The cohort point
- **WHEN** a unit of the full-cohort point is built
- **THEN** it holds every participant and has no held-out reference

### Requirement: A plant edits a covariate-adjusted stage configuration

A plant SHALL be the difference between an edited and an unedited stage configuration, added to every row of the second group at each stage under the unit's permuted labels. The configuration SHALL be the stage means adjusted for the case study's covariates, in the analysis sample's pooled-standardized coordinates, and the shift SHALL be mapped back to each block's units.

#### Scenario: Only the planted group moves
- **WHEN** a plant is applied
- **THEN** every row of the first group is unchanged, every row of the second group at one stage moves by the same vector, and covariates stay on their rows

#### Scenario: A stretch
- **WHEN** a stretch of effect e is planted
- **THEN** the second group's rows at each stage move by e times that stage's deviation from the configuration's centroid

#### Scenario: A rotation
- **WHEN** a rotation of effect e is planted
- **THEN** the shift is that of a rigid rotation of the configuration about its centroid by 90°·e in the plane of its net displacement and its second direction

#### Scenario: A bend
- **WHEN** a bend of effect e is planted
- **THEN** the shift is that of moving the middle stage off the chord between its neighbours by e mean steps and rescaling the configuration to its path length

#### Scenario: The adjusted configuration is the one the trajectory model measures
- **WHEN** the stage means of the planted group's rows are adjusted for the covariates before and after a plant
- **THEN** they differ by the plant's shift exactly

#### Scenario: Zero effect
- **WHEN** a cell's effect is zero
- **THEN** the blocks are returned untouched

### Requirement: A plant is defined on the sample or on held-out participants

Each plant SHALL be defined on one of two references: `sample`, the adjusted configuration of the unit's analysis sample, or `heldout`, the adjusted configuration of the participants outside it, standardized by the analysis sample's scaler. Held-out plants SHALL exist at the 750-participant point only.

#### Scenario: A held-out configuration uses no analysed participant
- **WHEN** a held-out plant is built
- **THEN** its stage means and covariate adjustment are computed from participants outside the unit's sample only

#### Scenario: The cohort point has no held-out cell
- **WHEN** the cells of the full-cohort point are enumerated
- **THEN** every plant cell is on the sample reference

### Requirement: Each design point holds two controls and a fixed set of plant cells

Each point SHALL hold an anchor and a translation of effect 1.0 along the sample's net displacement (500 replicates each), and plants at 200 replicates: on `sample` a stretch at 0.25, 0.5 and 1.0, a rotation at 0.25, 0.5 and 1.0 and a bend at 0.5 and 1.0; on `heldout` a stretch at 0.5, 1.0 and 2.0, a rotation at 0.5 and 1.0 and a bend at 1.0. Every unit SHALL run 999 permutations.

#### Scenario: Unit counts
- **WHEN** the plants profile is enumerated at the committed replicates
- **THEN** it holds 4,200 units at the 750-participant point and 3,000 at the full cohort, the out-of-span cells included

#### Scenario: Cells of one replicate are matched
- **WHEN** two cells of one design point are built at the same replicate
- **THEN** they share the sample, the labels and the permutation seed, and differ only in the plant

#### Scenario: The translation plants no interaction
- **WHEN** the translation control is planted
- **THEN** every row of the second group moves by the same vector at every stage

### Requirement: Out-of-span cells are recorded and never judged

Each design point SHALL also hold, at 200 replicates, a rotation at 0.25 and at 0.5 on the sample reference whose plane is spanned by the configuration's net displacement and a seeded unit direction orthogonal to the configuration's stage span. Both arms' rates on these cells SHALL be reported for the three statistics, and no criterion SHALL read them.

#### Scenario: The direction is outside the stage span
- **WHEN** the out-of-span direction of a unit is drawn
- **THEN** it has unit norm and is orthogonal to every centred stage mean of the configuration to 1e-10

#### Scenario: The direction is reproducible
- **WHEN** one unit is built twice
- **THEN** its out-of-span direction is the same, and it is the same at both effects of the replicate

#### Scenario: The cells do not decide
- **WHEN** an arm rejects or fails to reject on an out-of-span cell
- **THEN** every criterion and the verdict are unchanged

### Requirement: The simulation profile plants stage-span edits at the ADNI design

The simulation profile SHALL generate each dataset by the committed ADNI-design construction at its calibrated strengths and ADNI's cell table, with the stage-span construction under the planted scale over the measured blocks, from a seed family no committed check uses. Its cells SHALL be the anchor and a translation at 1.0 (500 replicates each) and, at 200 replicates, a stretch at 0.5, 1.0 and 2.0, a rotation at 0.25, 0.5 and 1.0 and a bend at 1.0 and 2.0, all at 999 permutations.

#### Scenario: Unit count
- **WHEN** the simulation profile is enumerated at the committed replicates
- **THEN** it holds 2,600 units

#### Scenario: Cells of one replicate share everything but the plant
- **WHEN** two cells are generated at the same replicate
- **THEN** they share the baseline draw, the covariates, the phase draw and the permutation seed

#### Scenario: The intended geometry is recorded
- **WHEN** a rotation of effect e is generated
- **THEN** its truth records a population `angle` of 90°·e with `delta` and `shape` of zero in the planted coordinates

#### Scenario: The anchor is the construction's null dataset
- **WHEN** the anchor is generated at a seed
- **THEN** its blocks equal those of the ADNI-design construction's anchor at that seed

### Requirement: Each unit measures one production fit in two coordinate systems

On each unit the check SHALL fit the production stage-PLS model once at its cross-validated rank and evaluate on that fit `pls_cv` and `predictive_orthonormal` under RRPP's unrestricted draws and both under within-stage draws, at the case study's configuration with the `df` residual correction and one permutation seed. A full-cohort anchor unit SHALL fit at the pre-flight's rank without a cross-validation.

#### Scenario: One cross-validation per unit
- **WHEN** a unit other than a full-cohort anchor runs, with or without a plant, on ADNI's blocks or in simulation
- **THEN** exactly one cross-validation is made on the unit's data, and the other arms are fitted at the rank it selects

#### Scenario: A cohort anchor unit
- **WHEN** an anchor unit of the full-cohort point runs
- **THEN** no cross-validation is made and every arm reports the pre-flight's rank

#### Scenario: A failed arm
- **WHEN** one arm raises
- **THEN** the arm is recorded as failed with the error and the unit's other arms are evaluated

#### Scenario: A failed unit
- **WHEN** building the unit or its production fit raises
- **THEN** the unit is recorded as failed with the error and the shard continues

### Requirement: Each unit checks its arms and its plant

Each unit SHALL check the rank of every forced arm, that every arm fitted the production model, the dimension the predictive arms measure, that each within-stage arm observes its unrestricted counterpart's statistics, the realized cell table and, on ADNI's blocks, that the realized plant is the declared one. A unit failing a check SHALL be recorded as inconsistent and excluded from every rate.

#### Scenario: The plant is the declared one
- **WHEN** a planted unit of the plants profile is checked
- **THEN** the second group's stage means moved by the declared shift and the first group's did not, to 1e-8

#### Scenario: Same fit
- **WHEN** an arm is fitted at the unit's rank
- **THEN** its rotations and response loadings equal the production model's to 1e-8

#### Scenario: Measured dimension
- **WHEN** a predictive arm is evaluated on a three-stage fit of rank 2 or more
- **THEN** it measures two dimensions

#### Scenario: Inconsistent unit
- **WHEN** a check fails
- **THEN** the unit is recorded as inconsistent with the check's name, the arm and the deviation, and it is excluded from every rate

### Requirement: Pre-flights and a reproduction gate tie both profiles to committed runs

Before any plants unit runs, the cohort's cross-validated fit SHALL select rank 8 at a mean held-out stage AUROC of 0.599514 and equal the forced call exactly. Before any simulated unit runs, two anchor units of the committed ADNI-design run SHALL reproduce its committed p-values exactly. At merge the full-cohort anchor's counts SHALL equal the cohort real-noise anchor's committed counts for the four arms.

#### Scenario: The pre-flights pass
- **WHEN** the cohort's fit and the two simulated anchor units reproduce their committed values
- **THEN** each pre-flight is recorded as passed with the compared values

#### Scenario: A pre-flight fails
- **WHEN** the rank, the AUROC or a p-value differs from its committed value
- **THEN** the pre-flight is recorded as failed, naming both values, and `run` exits for that profile

#### Scenario: The anchor does not reproduce
- **WHEN** a full-cohort anchor count differs from the committed one
- **THEN** the merge records the reproduction as failed, naming the arm, the statistic and both counts

### Requirement: Records hold no participant data

A unit's record SHALL hold seeds, counts, ranks, statistics, p-values, summary readings of the fit and of the plant, and timings only. No record, table or committed file SHALL hold a participant identifier, a label vector, a subsample's membership or a plant's shift vector.

#### Scenario: A recorded unit
- **WHEN** a unit of the plants profile is recorded
- **THEN** no identifier of the export appears in the record, and no field of the record holds one entry per participant or per feature

### Requirement: Runs are resumable, sharded and bound to their settings

Each profile SHALL write its settings with their hash, the hashes of the export's files where it reads the export, and the sha256 of every committed parent file it reads. A resume, merge or decision under other settings SHALL exit naming the differing keys, a changed parent file or export file SHALL be refused by name, shards SHALL partition the units round-robin, and a resume SHALL skip complete units.

#### Scenario: Resume skips completed units
- **WHEN** a shard is rerun
- **THEN** units with a complete record are not evaluated again

#### Scenario: Settings mismatch is refused
- **WHEN** a run resumes in a directory whose settings differ
- **THEN** it exits naming the differing keys

#### Scenario: A changed parent is refused
- **WHEN** a committed parent table's hash differs from the recorded one
- **THEN** `run`, `merge` and `decide` exit naming the file

### Requirement: Merge writes the rate, paired and cell tables

The merge of a profile SHALL write a rate table (per design point, cell, arm and statistic: available and rejected counts, rate, Monte Carlo SE, failures), a paired table (`predictive_orthonormal` against `pls_cv` under each null, and each within-stage arm against its counterpart: both rates, the discordant counts and exact p-values), a cell table and an accounting of units, failures, inconsistent units and redraws.

#### Scenario: Rejection rule
- **WHEN** a rate is computed
- **THEN** a unit rejects a statistic when its p-value is at most 0.05, and failed or inconsistent units are not counted

#### Scenario: No discordant units
- **WHEN** an arm and its reference reject on the same units
- **THEN** every paired p-value is 1

#### Scenario: Cell table
- **WHEN** the cell table of the plants profile is written
- **THEN** it holds, per design point and cell, the units by cross-validated rank, the median held-out AUROC, the plant's median norm per stage, the largest relative change of a pooled SD, and the shares of the plant and of its reference configuration inside the fitted span and inside the predictive subspace
- **AND** per design point with a held-out reference, the median reliability of the held-out configuration between two halves of the reference

### Requirement: The decision judges size on the controls this change adds

The decision SHALL judge `pls_cv` and `predictive_orthonormal` for size by the count-plus-ceiling rule of the stage-span paper grade on two families: on ADNI's blocks, the anchor at 750 participants and the translation at both points (nine readings); in simulation, the anchor and the translation (six readings). The full-cohort anchor is the anchor run's reading and SHALL enter through the reproduction only.

#### Scenario: Size holds
- **WHEN** an arm's nine readings on ADNI's blocks have at most one above 34 of 500 and none above 38, and its six simulated readings at most one above 34 and none above 37
- **THEN** the arm holds size

#### Scenario: An excess
- **WHEN** a family has two readings above the bound or one above its ceiling
- **THEN** the arm does not hold size and the readings are named

#### Scenario: Unavailable reading
- **WHEN** a reading of a family has no available unit
- **THEN** the arm does not hold on that family and the reading is named

### Requirement: The decision judges protected power on every plant's target statistic

The decision SHALL judge protected power on one family: the target statistic of every simulated, held-out and in-sample plant cell (`delta` for a stretch, `angle` for a rotation, `shape` for a bend). Each reading SHALL be the one-sided exact paired test of a loss of `predictive_orthonormal` against `pls_cv` on the discordant units, and the family SHALL be judged by the count-plus-ceiling rule on those p-values.

#### Scenario: Power is protected
- **WHEN** at most two of the thirty readings have a loss p-value below 0.025 and none is below 0.05 / 30
- **THEN** the criterion holds

#### Scenario: A loss
- **WHEN** three readings are below 0.025, or one is below 0.05 / 30
- **THEN** the criterion fails and the cells are named with both rates and the discordant counts

### Requirement: The decision judges specificity where the off-target truth is exact

The decision SHALL judge specificity on two families of off-target readings, each the one-sided exact paired test of an excess of `predictive_orthonormal` over `pls_cv`, under the count-plus-ceiling rule: the simulated plants (sixteen readings) and the in-sample plants on ADNI's blocks (thirty-two). Off-target rates of the held-out plants SHALL be reported and not judged.

#### Scenario: Specificity holds
- **WHEN** each family has at most two excess p-values below 0.025 and none below 0.05 divided by its size
- **THEN** the criterion holds

#### Scenario: An excess in one family
- **WHEN** either family fails the rule
- **THEN** the criterion fails, naming the family, the cells and the statistics

#### Scenario: A held-out plant's off-target rate
- **WHEN** an arm rejects a held-out stretch by `angle` or `shape`
- **THEN** the rate is reported beside the reference's reliability and no criterion reads it

### Requirement: The decision judges gain in simulation at the cohort's table

The decision SHALL judge gain on the target statistics of the eight simulated plant cells: the criterion holds when at least one exact paired gain p-value of `predictive_orthonormal` over `pls_cv` is at most 0.05 / 8. Gains on the held-out and on the in-sample plants SHALL be judged the same way, recorded beside it, and SHALL NOT qualify the candidate.

#### Scenario: A simulated gain
- **WHEN** one simulated target cell has a gain p-value at most 0.05 / 8
- **THEN** the criterion holds and the gaining cells are named

#### Scenario: A real-noise gain without a simulated one
- **WHEN** no simulated cell gains and a held-out cell does
- **THEN** the criterion fails, and the held-out gain is recorded

### Requirement: The decision writes one predeclared verdict on the candidate

The decision SHALL write `adni_gain_specificity_decision.json` with a rule version and one verdict, the first that applies: `gate_failed`; `production_excess`; `size_excess`; `power_loss`; `specificity_excess`; `no_gain`; `reopened`. Every criterion SHALL be evaluated whichever fails first, and each verdict SHALL carry the consequence statement fixed for it.

#### Scenario: Gate
- **WHEN** a pre-flight is missing or failed, the anchor reproduction failed or a unit is inconsistent
- **THEN** the verdict is `gate_failed`, and the criteria are still recorded

#### Scenario: Production exceeds
- **WHEN** `pls_cv` does not hold size on a family
- **THEN** the verdict is `production_excess`, naming the readings, whatever the candidate reads

#### Scenario: The first failing criterion names the verdict
- **WHEN** `pls_cv` holds size and `predictive_orthonormal` fails specificity and gain
- **THEN** the verdict is `specificity_excess` and both failing criteria are listed

#### Scenario: Reopened
- **WHEN** production holds size and the candidate holds size, protected power and specificity and gains in simulation
- **THEN** the verdict is `reopened`

#### Scenario: Below grade
- **WHEN** a control has fewer than 500 replicates, a plant fewer than 200, or a unit fewer than 999 permutations
- **THEN** the decision is marked `indicative`

### Requirement: The decision states what the verdict licenses

The decision SHALL carry, with its rule version, what its verdict licenses: it speaks for ADNI's export, the two design points and the simulated table; `reopened` permits a later change to propose the candidate's paper-grade run with a rank rule and changes no default; the SEA-AD reading stands for SEA-AD's design; and no held-out plant was read at the full cohort.

#### Scenario: Reopened changes no default
- **WHEN** the verdict is `reopened`
- **THEN** the consequence states that `predictive_orthonormal` remains an opt-in diagnostic until a later paper-grade run

#### Scenario: A closed verdict
- **WHEN** the verdict is `size_excess`, `power_loss`, `specificity_excess` or `no_gain`
- **THEN** the consequence states that open item 1's second reopening condition is not met for ADNI and names the criterion and its readings

### Requirement: The decision states what production detects on the cohort's noise

The decision SHALL state, for `pls_cv` alone on ADNI's blocks, the smallest effect at which each target statistic reaches a rejection rate of 0.80 on the held-out and on the in-sample plants at each design point, or that none does, and whether its off-target rates on the in-sample plants hold the count-plus-ceiling rule as one family. The statement SHALL NOT change the verdict.

#### Scenario: An informative effect
- **WHEN** production's `delta` reaches 0.80 on a held-out stretch from some effect upward
- **THEN** the statement names the design point, the reference and that effect with its rate

#### Scenario: No informative effect
- **WHEN** a target statistic stays below 0.80 at every effect of a reference
- **THEN** the statement gives its largest rate and the effect it was read at

#### Scenario: Production's own off-target rates
- **WHEN** production's thirty-two off-target readings on the in-sample plants have more than two above 16 of 200 or one above 20
- **THEN** the statement names the cells and statistics that exceed

### Requirement: Descriptive readings never change the verdict

Beside the verdict the decision SHALL report: both arms' rates in every cell; the within-stage arms; the out-of-span cells; each target rate on the sample reference beside the held-out one at the same mode and effect; rates among units whose rank exceeds the number of stages minus one; and the cell table's capture and reliability readings.

#### Scenario: Descriptive readings do not decide
- **WHEN** a within-stage arm, a rank class, an out-of-span cell or a capture share is high or low
- **THEN** the verdict and production's statement are those the rules give

### Requirement: The run is recorded in a dated report

The run SHALL be recorded in a dated findings report and committed aggregates: each profile's settings, pre-flight record, three tables and accounting, the decision and the provenance. The roadmap SHALL carry the verdict, production's statement and what the verdict licenses.

#### Scenario: Report contents
- **WHEN** the report is written
- **THEN** it states the verdict with each criterion's readings, production's statement, the bracket between the two references, the out-of-span cells, every descriptive reading, what follows and the limitations
- **AND** every number it quotes is in a committed aggregate

#### Scenario: Roadmap
- **WHEN** the run is recorded
- **THEN** roadmap open items 1 and 2, the entry on real-noise power and specificity and the Done list state the verdict and its consequence

#### Scenario: Nothing participant-level is committed
- **WHEN** the aggregates and the report are committed
- **THEN** no committed file holds a participant identifier, a label vector, a subsample's membership or a plant's shift vector
