# cohort-real-noise-validation Specification

## Purpose
Reads the size of the trajectory tests on a cohort's exported omic blocks, with the group labels permuted within stage so that the null is true on real noise. Participants are subsampled within cells into design points, and a predeclared decision says at which design points each measurement holds size and what that reading does and does not license.

## Requirements

### Requirement: The cohort export is read in place through an input contract

The check SHALL read the ADNI export in place: a manifest recording each participant-level file's hash, the metadata (participant, sex, stage), the covariates, and the `p180` and `nmr` blocks. It SHALL verify the files against the manifest, the two groups, the stage order CN / MCI / AD, the cell table, aligned participants, finite block values and complete covariates at their declared levels. A violation SHALL raise a named error before any analysis.

#### Scenario: A violation is named
- **WHEN** the export lacks a declared file, a block lists other participants than the metadata, a stage label is outside the declared order, a cell is empty, a block value is not finite or a covariate value is missing
- **THEN** the check exits naming the violation and the file, and no unit is evaluated

#### Scenario: A file that differs from its manifest is refused
- **WHEN** a participant-level file's hash differs from the one the export's manifest records
- **THEN** the check exits naming the file

#### Scenario: Another export is refused on resume
- **WHEN** a run resumes with an export whose file hashes differ from those in its settings
- **THEN** it exits naming the differing files before any unit is evaluated

#### Scenario: The export is not copied
- **WHEN** a run, a merge or a decision completes
- **THEN** no file of the export has been copied or modified, and nothing has been written to the export's directory

### Requirement: Design points subsample participants within group × stage cells

The check SHALL read four design points: 80, 200 and 500 participants drawn without replacement within each sex × stage cell, at the cohort's proportions rounded by largest remainder, and the full cohort. Each point SHALL have 500 replicates at 999 permutations. A subsample in which a level of a categorical covariate is absent SHALL be refused and redrawn from the replicate's next seed, and the redraw SHALL be counted.

#### Scenario: Cell tables
- **WHEN** the design points are built on the cohort's table (Female 237 / 310 / 125, Male 223 / 448 / 158 at CN / MCI / AD)
- **THEN** the subsampled tables are Female 13 / 16 / 7 and Male 12 / 24 / 8 at n = 80, Female 31 / 41 / 17 and Male 30 / 60 / 21 at n = 200, and Female 79 / 103 / 42 and Male 74 / 149 / 53 at n = 500

#### Scenario: A subsample holds each participant once
- **WHEN** a subsample is drawn
- **THEN** every participant appears at most once, in the cell of its real sex and stage

#### Scenario: Subsamples are reproducible and differ between replicates
- **WHEN** one replicate is drawn twice and two replicates are drawn once each
- **THEN** the first two subsamples are identical and the last two differ

#### Scenario: An absent covariate level
- **WHEN** a drawn subsample holds no participant at some level of a categorical covariate
- **THEN** it is refused, the next seed's subsample is drawn, and the unit records one redraw

#### Scenario: The cohort point
- **WHEN** a unit of the full-cohort point is built
- **THEN** it holds every participant and draws no subsample

### Requirement: Only within-stage-permuted labels are analysed

Each unit SHALL permute the group labels among the participants of each stage of its sample, which keeps the cell table, and SHALL leave each participant's blocks and covariates on its own row. A draw equal to the sample's real labels SHALL be refused and redrawn, and the redraw counted. No analysis SHALL run under the real labels.

#### Scenario: The cell table is kept
- **WHEN** the labels of a unit are drawn
- **THEN** every group × stage cell has the size of the design point's table

#### Scenario: A draw equal to the real labels
- **WHEN** a within-stage permutation reproduces the sample's real labels
- **THEN** it is refused and the next seed's draw is used

### Requirement: Seeds derive from the design point and the replicate

Every seed SHALL derive from the driver's name, base seed 2039, the design point, the replicate and a purpose: the subsample, the labels, the permutation draws and the folds. Every arm of a unit SHALL share the unit's sample, labels and permutation seed.

#### Scenario: Arms of one unit are paired
- **WHEN** two arms of one unit are evaluated under the same null
- **THEN** they use the same participants, the same labels and the same sequence of permutation draws

#### Scenario: Design points do not share seeds
- **WHEN** the same replicate index is evaluated at two design points
- **THEN** its subsample, label and permutation seeds differ between them

### Requirement: Operating arms measure the production fit in two coordinate systems under both nulls

On each unit the check SHALL fit the production stage-PLS model at its cross-validated rank and evaluate on that fit `pls_cv`, `pls_cv__within_stage`, `predictive_orthonormal` and `predictive_orthonormal__within_stage`. Every arm SHALL use the case-study configuration: the five covariates in the trajectory model, rank floor 2, net displacement and the `df` residual correction. A cohort unit SHALL fit at the pre-flight's rank without a cross-validation.

#### Scenario: A subsampled unit
- **WHEN** a unit of a subsampled design point runs
- **THEN** exactly one cross-validation is made, and the other operating arms are fitted at the rank it selects

#### Scenario: A cohort unit
- **WHEN** a unit of the full-cohort point runs
- **THEN** no cross-validation is made and every operating arm reports the pre-flight's rank

#### Scenario: A failed arm
- **WHEN** one arm raises
- **THEN** the arm is recorded as failed with the error and the unit's other arms are evaluated

#### Scenario: A failed unit
- **WHEN** building the unit or its production fit raises
- **THEN** the unit is recorded as failed with the error and the shard continues

### Requirement: Fixed-rank arms read the same sample at ranks 8 and 19

Each unit SHALL also evaluate the four combinations of coordinate system (scores, `predictive_orthonormal`) and null (unrestricted, within stage) on stage-PLS fits forced at rank 8 and at rank 19, on the unit's sample, labels and permutation seed and at the operating arms' configuration. At the full-cohort point the rank-8 combinations are the operating arms and are not evaluated twice.

#### Scenario: A subsampled unit
- **WHEN** a unit of a subsampled design point completes
- **THEN** it holds eight fixed-rank arms, four at each rank, and none of them made a cross-validation

#### Scenario: A cohort unit
- **WHEN** a unit of the full-cohort point completes
- **THEN** it holds four fixed-rank arms at rank 19

### Requirement: A control arm drops the residual correction

Each unit SHALL evaluate `pls_cv__uncorrected`: the production scores of the unit's fit under the unrestricted null with no residual correction and every other setting of `pls_cv`.

#### Scenario: Only the null moves
- **WHEN** the control arm of a unit is compared with its `pls_cv` arm
- **THEN** their observed `delta`, `angle` and `shape` are equal exactly

### Requirement: The pre-flight pins the cohort's fit to the case study's

Before any unit runs, the cross-validated production analysis of the full cohort under a permuted label draw SHALL select rank 8 at a mean held-out stage AUROC of 0.599514 to six decimals, the case study's primary values, and its statistics and p-values SHALL equal those of the analysis forced at that rank exactly. `run` SHALL refuse to start without a passed pre-flight under the same settings.

#### Scenario: The pre-flight passes
- **WHEN** the cross-validated fit selects rank 8 at AUROC 0.599514 on both pre-flight replicates and equals the forced call
- **THEN** the pre-flight is recorded as passed with the rank, the AUROC and both comparisons

#### Scenario: The fit differs
- **WHEN** the rank or the AUROC differs from the case study's
- **THEN** the pre-flight is recorded as failed, naming both values, and `run` exits

### Requirement: Each unit checks its arms against one another

Each unit SHALL check its arms against one another: the rank of every forced arm, the fitted model of every arm at one rank, the dimension each predictive arm measures, the observed statistics under the two nulls and the realized cell table. A unit failing a check SHALL be recorded as inconsistent and excluded from every rate.

#### Scenario: Forced rank
- **WHEN** an arm is forced at a rank
- **THEN** its record reports that rank

#### Scenario: Same fit
- **WHEN** an operating arm is fitted at the unit's rank
- **THEN** its rotations and response loadings equal the production model's to 1e-8

#### Scenario: Measured dimension
- **WHEN** a predictive arm is evaluated on a three-stage fit of rank 2 or more
- **THEN** it measures two dimensions

#### Scenario: Cell table
- **WHEN** a unit is evaluated
- **THEN** its realized group × stage table equals the design point's

#### Scenario: The two nulls read one measurement
- **WHEN** a within-stage arm is compared with the same coordinate system and rank under unrestricted draws
- **THEN** their observed statistics are equal exactly

#### Scenario: Inconsistent unit
- **WHEN** a check fails
- **THEN** the unit is recorded as inconsistent with the check's name, the arm and the deviation, and it is excluded from every rate

### Requirement: Records hold no participant data

A unit's record SHALL hold seeds, counts, ranks, statistics, p-values, summary readings of the fit and timings only. No record, table or committed file SHALL hold a participant identifier, a label vector, a subsample's membership or a fold assignment.

#### Scenario: A recorded unit
- **WHEN** a unit is recorded
- **THEN** no identifier of the export appears in the record, and no field of the record holds one entry per participant

### Requirement: Runs are resumable, sharded and bound to their settings

A run SHALL write its settings with their hash, the hashes of the export's files and the sha256 of every committed parent file it reads. A resume, merge or decision under other settings SHALL exit naming the differing keys, and a changed parent file SHALL be refused by name. Shards SHALL partition the units round-robin, and a resume SHALL skip complete units.

#### Scenario: Resume skips completed units
- **WHEN** a shard is rerun
- **THEN** units with a complete record are not evaluated again

#### Scenario: Settings mismatch is refused
- **WHEN** a run resumes in a directory whose settings differ
- **THEN** it exits naming the differing keys

#### Scenario: A changed parent is refused
- **WHEN** the committed SEA-AD anchor rate table's hash differs from the recorded one
- **THEN** `run`, `merge` and `decide` exit naming the file

### Requirement: Merge writes the rate, paired and point tables

The merge SHALL write a rate table (per design point, arm and statistic: available and rejected counts, rate, Monte Carlo SE, failures), a paired table (each within-stage and control arm against its counterpart, and each predictive arm against the scores at the same rank and null: both rates, the discordant counts and exact p-values), a point table and an accounting of units, failures, inconsistent units and redraws.

#### Scenario: Rejection rule
- **WHEN** a rate is computed
- **THEN** a unit rejects a statistic when its p-value is at most 0.05, and failed or inconsistent units are not counted

#### Scenario: No discordant units
- **WHEN** an arm and its counterpart reject on the same units
- **THEN** every paired p-value is 1

#### Scenario: Point table
- **WHEN** the point table is written
- **THEN** it holds, per design point, the counts of units by cross-validated rank and the share at the number of stages minus one, the median held-out AUROC, and per fit rank the median in-sample stage R² beside the median share of the stage indicators' variance explained out of fold
- **AND** per coordinate system and fit rank, the median residual dispersion of each stage relative to the pooled one

### Requirement: An arm's size scope is the largest upper set of design points on which it holds size

The decision SHALL judge an arm on a set of design points by the count-plus-ceiling rule of the stage-span paper grade over that set's `delta`, `angle` and `shape` readings as one family. The arm's size scope SHALL be the first of these sets on which it holds: all four points; 200 participants and above; 500 and above; the full cohort alone. An arm whose full-cohort readings fail SHALL have no scope.

#### Scenario: Size held on the whole ladder
- **WHEN** an arm's twelve readings hold the rule
- **THEN** its scope starts at 80 participants

#### Scenario: An excess at the smallest point only
- **WHEN** the twelve readings fail the rule and the nine readings at 200 participants and above hold it
- **THEN** the arm's scope starts at 200 participants

#### Scenario: A point is not skipped
- **WHEN** an arm exceeds at 200 participants and holds at 80, at 500 and at the full cohort
- **THEN** its scope does not include 80 participants

#### Scenario: The cohort fails
- **WHEN** the arm's three full-cohort readings fail the rule
- **THEN** the arm has no scope

#### Scenario: Unavailable reading
- **WHEN** a reading of a set has no available unit
- **THEN** the arm does not hold on that set and the reading is named

### Requirement: The decision writes one predeclared verdict on the operating arms

The decision SHALL write `adni_anchor_decision.json` with a rule version and one verdict, the first that applies: `gate_failed` when the pre-flight did not pass or a unit is inconsistent; `production_excess` when `pls_cv` does not hold on all four points; `excess_persists` when `predictive_orthonormal` has no scope; `excess_scoped` when its scope starts above 80 participants; `excess_absent` otherwise. Each verdict SHALL carry the consequence statement fixed for it.

#### Scenario: Gate
- **WHEN** the pre-flight record is missing or failed, or a unit is inconsistent
- **THEN** the verdict is `gate_failed`, and the scopes are still recorded

#### Scenario: Production exceeds
- **WHEN** `pls_cv`'s scope does not start at 80 participants
- **THEN** the verdict is `production_excess`, naming the readings that exceed, whatever the other arms read

#### Scenario: The excess persists
- **WHEN** `pls_cv` holds on all four points and `predictive_orthonormal` has no scope
- **THEN** the verdict is `excess_persists`

#### Scenario: The excess is scoped
- **WHEN** `pls_cv` holds on all four points and `predictive_orthonormal`'s scope starts at 200 participants, at 500 or at the full cohort
- **THEN** the verdict is `excess_scoped`, naming where the scope starts

#### Scenario: The excess is absent
- **WHEN** both arms hold on all four points
- **THEN** the verdict is `excess_absent`

#### Scenario: The within-stage arms are recorded beside the verdict
- **WHEN** a within-stage operating arm does not hold on all four points
- **THEN** the decision names the arm and its exceeding readings, and the verdict is the one the unrestricted arms give

### Requirement: The decision grades the fixed-rank reading

The decision SHALL grade the fixed-rank arms on one family of 21 readings per arm: three statistics at ranks 8 and 19 on the three subsampled points and at rank 19 on the full cohort. The grade SHALL be the first that applies of `reference_exceeds`, `absent`, `reproduced` and `reproduced_unrestored`, with the statement fixed for it, and SHALL NOT change the verdict.

#### Scenario: The excess is reproduced and removed
- **WHEN** the scores hold, `predictive_orthonormal` fails under unrestricted draws and holds under within-stage draws
- **THEN** the grade is `reproduced`, naming the exceeding readings with their fit rank, design point and in-sample stage R²

#### Scenario: The excess is reproduced and stays
- **WHEN** the scores hold and `predictive_orthonormal` fails under both nulls
- **THEN** the grade is `reproduced_unrestored`

#### Scenario: No excess at any fixed rank
- **WHEN** the scores and `predictive_orthonormal` both hold under unrestricted draws
- **THEN** the grade is `absent`

#### Scenario: The reference exceeds
- **WHEN** the scores fail under unrestricted draws
- **THEN** the grade is `reference_exceeds`, and the predictive arms' readings are reported without a grade

### Requirement: The decision states what a size reading licenses

The decision SHALL carry, with its rule version, the rule that scopes a real-noise size reading: it speaks for the export and the design points it was read at; an option's scope on a cohort is its size scope; a hold on one cohort does not lift a failure on another, and does not extend to a cohort without its own reading; and size alone does not reopen an option, which also needs a gain and specificity at that design under their own predeclared decisions.

#### Scenario: A hold on ADNI and the SEA-AD veto
- **WHEN** the verdict is `excess_absent` or `excess_scoped`
- **THEN** the consequence states that the SEA-AD anchor's reading stands for SEA-AD's design and that `predictive_orthonormal` remains an opt-in diagnostic

#### Scenario: Size alone
- **WHEN** an arm has a scope
- **THEN** the decision states that the scope is the size part of the reopening condition only, and names the readings still missing

#### Scenario: Descriptive covariates do not extend a scope
- **WHEN** the decision reports features per participant, the selected ranks or the in-sample stage fit beside a scope
- **THEN** it states that they describe the design points and license no reading at another design

### Requirement: Descriptive readings never change the verdict

Beside the verdict and the grade the decision SHALL report: every arm's rates per design point with their status under the rule; the operating rates among units whose rank exceeds the number of stages minus one; the control arm against `pls_cv`; the paired comparisons; the point table; each arm's median `angle` null 95th percentile; and SEA-AD's committed anchor counts beside ADNI's in one table.

#### Scenario: Descriptive readings do not decide
- **WHEN** the control arm, a rank class or a paired comparison exceeds or holds
- **THEN** the verdict and the grade are those the rule gives

#### Scenario: Cross-cohort table
- **WHEN** the decision is written
- **THEN** it holds, for the scores and `predictive_orthonormal` under both nulls, the SEA-AD anchor's committed counts beside each ADNI design point's, with features per participant, the fit rank and the in-sample stage R²

### Requirement: The run is recorded in a dated report

The run SHALL be recorded in a dated findings report and committed aggregates: the settings, the pre-flight record, the three tables, the accounting, the decision and the provenance. The roadmap SHALL carry the verdict, the grade and the scoping rule.

#### Scenario: Report contents
- **WHEN** the report is written
- **THEN** it states the verdict, each operating arm's scope and counts, the fixed-rank grade with its readings against the in-sample stage fit, every descriptive reading, what follows and the limitations
- **AND** every number it quotes is in a committed aggregate

#### Scenario: Roadmap
- **WHEN** the run is recorded
- **THEN** roadmap open items 1 and 2, the entry on real-noise evidence and the Done list state the verdict and its consequence

#### Scenario: Nothing participant-level is committed
- **WHEN** the aggregates and the report are committed
- **THEN** no committed file holds a participant identifier, a label vector or a subsample's membership
