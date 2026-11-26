# group-dispersion-validation Specification

## Purpose
Reads the size of the production trajectory test, and of the per-feature interaction tests on the same units, when the two groups differ in noise dispersion and in size but not in trajectory. A predeclared decision states whether size needs equal group dispersion in balanced and in unbalanced designs, and a map gives the scale at which each statistic first fails.

## Requirements

### Requirement: Every unit is a null unit evaluated as the univariate baseline evaluates one

Every unit SHALL carry no group difference in trajectory: its cell is the zero-effect anchor or the translation control. It SHALL be evaluated as the univariate baseline evaluates one: one dataset, one production stage-PLS fit with net displacement and rank floor 2, and every arm read at the unit's RRPP seed. A failure SHALL be recorded for the unit without aborting the shard.

#### Scenario: One cross-validation per unit
- **WHEN** a unit runs
- **THEN** exactly one cross-validated PLS fit is made, and the record of every arm the univariate baseline defines equals that driver's record for the same generator parameters, cell and seeds

#### Scenario: Failed unit
- **WHEN** generation, fitting or an arm fails for a unit
- **THEN** the unit is recorded as failed with the error, and the shard continues

### Requirement: The core profile crosses the dispersion scale with group balance

The `core` profile SHALL run designs on methylation, expression and proteomics at n = 200, three stages and stage effect δ = 2.0. They cross the second group's dispersion scale s ∈ {1.25, 1.5, 2} with a first-group share of 0.5, 0.65, 0.8, 0.35 and 0.2, and add the equal-dispersion designs at shares 0.65 and 0.8. Each design SHALL have the anchor and the translation at e = 1.0 under the planted scale, at 500 replicates and 999 permutations.

#### Scenario: Designs and units
- **WHEN** the profile enumerates its new designs
- **THEN** there are 17 of them and 17,000 units

#### Scenario: Only dispersion and balance vary
- **WHEN** a new design's generator parameters are compared with those of the univariate baseline's design at n = 200 and δ = 2.0
- **THEN** they differ only in `group_ratio` and `group_dispersion_scale`

#### Scenario: The share sets which group is wider
- **WHEN** a design has a dispersion scale above 1 and a first-group share of 0.65 or 0.8
- **THEN** the wider group is the smaller one, and at a share of 0.35 or 0.2 it is the larger one

#### Scenario: New datasets
- **WHEN** a new design is named
- **THEN** its name differs from every design of the univariate baseline and of the stage-effect scaling run, so its dataset seeds are derived apart from theirs
- **AND** the two cells of one design and replicate share their null dataset

### Requirement: The columns profile reads a reduced cross at the other committed designs

The `columns` profile SHALL run four committed designs: the univariate baseline's n = 80 and n = 500 at δ = 2.0 and n = 200 at δ = 0.5, and the stage-effect scaling run's `band_n320` (δ = 0.05 at n = 320). At each it SHALL read the anchor at 500 replicates and 999 permutations for s ∈ {1.5, 2} at first-group shares 0.5, 0.8 and 0.2, and for equal dispersion at share 0.8. Each column's balanced equal-dispersion reading SHALL be the committed one.

#### Scenario: Designs and units
- **WHEN** the profile enumerates its designs
- **THEN** there are 28 of them, seven per column, and 14,000 units

#### Scenario: Each column differs from its committed design in dispersion and balance only
- **WHEN** a column design's generator parameters are compared with its committed design's
- **THEN** they differ only in `group_ratio` and `group_dispersion_scale`

#### Scenario: The band column's stage effect is the committed one
- **WHEN** the `band_n320` column is built
- **THEN** its stage effect is read from the stage-effect scaling run's committed settings, and a changed settings file is refused by name

### Requirement: The committed balanced design is rerun and must reproduce

The `core` profile SHALL rerun the anchor and the translation of the univariate baseline's design at n = 200 and δ = 2.0 with that run's name, seeds and grade. For both cells and every arm the available and rejected counts SHALL equal the committed rate table exactly. The committed settings, rate tables and cell tables of both parent runs SHALL be bound by sha256 and refused by name if changed.

#### Scenario: Reproduction holds
- **WHEN** every row of the rerun cells matches the committed counts
- **THEN** the merge records the reproduction as passed with the number of rows compared

#### Scenario: Reproduction fails
- **WHEN** any row differs
- **THEN** the merge records the reproduction as failed and lists every differing row with both counts

#### Scenario: A changed parent is refused
- **WHEN** a bound parent file's sha256 differs from the settings'
- **THEN** `run`, `merge` and `decide` exit naming the file

### Requirement: Within-group draws are recorded on the same fit

For every unit the driver SHALL also read `delta`, `angle` and `shape` under RRPP draws that permute the reduced model's residuals within each group, on the unit's production fit, observed statistics and RRPP seed. The arm SHALL be recorded and SHALL NOT bear on the verdict.

#### Scenario: Same fit and statistics
- **WHEN** the within-group arm is read for a unit
- **THEN** its observed statistics equal the production arm's exactly
- **AND** a unit where they differ is marked inconsistent for that arm and leaves its rates

#### Scenario: One stratum reproduces production
- **WHEN** the arm's inference path is given a single stratum holding every sample
- **THEN** its p-values equal the production arm's exactly

#### Scenario: Recorded only
- **WHEN** the decision is computed with the within-group arm's rates changed
- **THEN** the verdict is unchanged

### Requirement: Size is judged in four families

Size SHALL be judged on `delta`, `angle` and `shape`, each at p < 0.05, under the stage-span paper grade's count-plus-ceiling rule, in four families over both profiles: `unbalanced_equal`, the equal-dispersion designs with unequal groups (m = 24); `balanced_unequal`, the designs with a scale and equal groups (m = 42); `wider_smaller`, the designs whose wider group is the smaller (m = 60); and `wider_larger`, those whose wider group is the larger (m = 60).

#### Scenario: Family membership
- **WHEN** the families are assembled
- **THEN** every new design's control cells enter exactly one family with three readings each, and the rerun design enters none

#### Scenario: A family holds
- **WHEN** a family's readings stay within the rule's allowed exceedances and none is above its ceiling
- **THEN** the family holds

#### Scenario: An unavailable reading fails its family
- **WHEN** a reading of a family has no available unit
- **THEN** the family does not hold and the reading is listed as unavailable

### Requirement: The verdict follows a fixed order

The decision SHALL be `reproduction_failed` when the rerun design differs from its committed counts; otherwise `unbalance_excess` when `unbalanced_equal` fails; otherwise `dispersion_excess` when `balanced_unequal` fails; otherwise `unbalanced_dispersion_excess` when `wider_smaller` or `wider_larger` fails; otherwise `size_holds`. Every family SHALL be recorded whatever the verdict, with the rule version.

#### Scenario: Balanced designs lose size
- **WHEN** the reproduction passes, `unbalanced_equal` holds and `balanced_unequal` fails
- **THEN** the verdict is `dispersion_excess`

#### Scenario: Only unbalanced designs lose size
- **WHEN** the reproduction passes, `unbalanced_equal` and `balanced_unequal` hold and `wider_smaller` fails
- **THEN** the verdict is `unbalanced_dispersion_excess`

#### Scenario: Unbalance alone loses size
- **WHEN** the reproduction passes and `unbalanced_equal` fails
- **THEN** the verdict is `unbalance_excess`, and the other three families are still recorded

#### Scenario: Size holds
- **WHEN** the reproduction passes and the four families hold
- **THEN** the verdict is `size_holds`

#### Scenario: The gate comes first
- **WHEN** the reproduction fails
- **THEN** the verdict is `reproduction_failed`, and the families are still recorded

### Requirement: A map by statistic is reported beside the verdict

The decision SHALL report, for each statistic, column and first-group share, the smallest scale whose reading is above its family's ceiling, and SHALL flag every reading below 0.05 minus two standard errors as conservative. It SHALL report the selected rank and the held-out AUROC by design. None of these SHALL change the verdict.

#### Scenario: First failing scale
- **WHEN** `delta` at n = 200, δ = 2.0 and share 0.5 is within the ceiling at s = 1.25 and above it at s = 1.5 and 2
- **THEN** the map gives 1.5 for that statistic, column and share

#### Scenario: No failing scale
- **WHEN** no reading of a statistic at a column and share is above the ceiling
- **THEN** the map states that none failed up to the largest scale read there

#### Scenario: Conservative readings are flagged
- **WHEN** a reading's rate is below 0.05 minus two standard errors
- **THEN** it is flagged conservative, and its family's judgement is unchanged

### Requirement: Recorded arms are labelled and not judged

`production` (the omnibus at 0.05 / 3), `scan`, `fdr`, `total`, `span`, `complement`, `pca_matched` and the within-group statistics SHALL each be given a size label per family under the same rule, with every arm's rate per design and cell. They SHALL NOT bear on the verdict.

#### Scenario: Labels per family
- **WHEN** the decision is computed
- **THEN** each recorded arm has one label per family, `holds` or `exceeds`, with its exceedances and its readings above the ceiling

#### Scenario: Recorded arms are outside the decision
- **WHEN** the decision is computed with a recorded arm's rates changed
- **THEN** the verdict and the map are unchanged

### Requirement: Runs are resumable and bound to their settings

The driver SHALL provide `run`, `merge` and `decide`. Settings SHALL be written with a hash, and a resume or a merge under different settings SHALL exit naming the differing keys. Units SHALL be sharded round-robin, complete units SHALL be skipped on resume, and records SHALL hold seeds, statistics and counts only.

#### Scenario: Resume skips complete units
- **WHEN** a shard is rerun after an interruption
- **THEN** units whose records are complete are not evaluated again

#### Scenario: Other settings are refused
- **WHEN** `run` or `merge` is invoked on a directory whose stored settings differ
- **THEN** it exits naming the differing keys

### Requirement: Findings are reported with provenance and limits

The run SHALL produce a dated report and committed aggregates: settings, rate and cell tables, accounting, the reproduction, the decision and provenance. Every number the report quotes SHALL be readable from a committed aggregate. The report SHALL state that the readings hold for simulated Gaussian noise, one scale applied to every feature of one group, three stages and no covariates, and that no real cohort was read.

#### Scenario: Numbers trace to aggregates
- **WHEN** the report quotes a rate, a count or a label
- **THEN** the value is in a committed CSV or JSON file under the run's results directory

#### Scenario: Limits are stated
- **WHEN** the report states the verdict
- **THEN** it states the limits beside it and what the verdict changes in the documentation and the roadmap
