# unbalance-diagnosis-validation Specification

## Purpose
Reads the size of the production trajectory test at equal group dispersion across the smaller group's cell size and group balance, to say whether the group-dispersion validation's one failing equal-dispersion reading replicates and what it follows. A predeclared decision gives the verdict and states how the group-dispersion validation's families are read under it.

## Requirements

### Requirement: Every unit is an equal-dispersion null unit evaluated as the group-dispersion validation evaluates one

Every unit SHALL be the zero-effect anchor of a design with equal group dispersion on methylation, expression and proteomics, three stages and stage effect δ = 2.0. It SHALL be evaluated as the group-dispersion validation evaluates one, with 999 permutations. A failure SHALL be recorded for the unit without aborting the shard.

#### Scenario: Same unit as the parent
- **WHEN** a unit runs
- **THEN** the record of every arm the group-dispersion validation defines equals that driver's record for the same generator parameters, cell and seeds

#### Scenario: Failed unit
- **WHEN** generation, fitting or an arm fails for a unit
- **THEN** the unit is recorded as failed with the error, and the shard continues

### Requirement: The designs cross the smaller group's cell size with balance

The run SHALL read seven new designs at 1,000 replicates each: n = 80 at first-group shares 0.8, 0.74, 0.65 and 0.5, n = 200 at shares 0.925 and 0.895, and n = 42 at share 0.5. Their smaller group SHALL hold 5, 7, 9 and 13 samples per stage at n = 80, 5 and 7 at n = 200, and 7 at n = 42.

#### Scenario: Designs and units
- **WHEN** the run enumerates its units
- **THEN** there are seven new designs and the rerun design, and 7,500 units

#### Scenario: Realized cell sizes
- **WHEN** a new design's dataset is generated
- **THEN** the smaller group's samples per stage equal the design's declared cell size

#### Scenario: New datasets
- **WHEN** a new design is named
- **THEN** its name differs from every design of the group-dispersion validation, the univariate baseline and the stage-effect scaling run, so its dataset seeds are derived apart from theirs

#### Scenario: Only sample size and balance vary
- **WHEN** a new design's generator parameters are compared with those of the committed design at n = 80 and share 0.8
- **THEN** they differ only in the sample size and the first-group share

### Requirement: The committed failing design is rerun and must reproduce

The run SHALL rerun the group-dispersion validation's equal-dispersion design at n = 80 and share 0.8 with that run's name, seeds and grade. For every arm the available and rejected counts SHALL equal the committed tables exactly. The committed settings, tables and decision of the group-dispersion validation SHALL be bound by sha256 and refused by name if changed.

#### Scenario: Reproduction holds
- **WHEN** every row of the rerun design matches the committed counts
- **THEN** the merge records the reproduction as passed with the number of rows compared

#### Scenario: Reproduction fails
- **WHEN** any row differs
- **THEN** the merge records the reproduction as failed and lists every differing row with both counts

#### Scenario: A changed parent is refused
- **WHEN** a bound parent file's sha256 differs from the settings'
- **THEN** `run`, `merge` and `decide` exit naming the file

### Requirement: Residual corrections are recorded on the same fit

For every unit the driver SHALL also read `delta`, `angle` and `shape` under the `df` and the `leverage` residual corrections, with unrestricted draws, on the unit's production fit, observed statistics and RRPP seed. These arms SHALL be recorded and SHALL NOT bear on the verdict.

#### Scenario: Same fit and statistics
- **WHEN** a corrected arm is read for a unit
- **THEN** its observed statistics equal the production arm's exactly
- **AND** a unit where they differ is marked inconsistent for that arm and leaves its rates

#### Scenario: No correction reproduces production
- **WHEN** the arm's inference path is run with no correction
- **THEN** its p-values equal the production arm's exactly

### Requirement: Size is judged in three families

Size SHALL be judged on `delta`, `angle` and `shape`, each at p < 0.05, under the stage-span paper grade's count-plus-ceiling rule, in three families: `balanced_small`, the balanced design at 7 samples per stage (m = 3); `unbalanced_small`, the four unbalanced designs whose smaller group holds 5 or 7 samples per stage (m = 12); and `reference`, the n = 80 designs at 9 and 13 samples per stage (m = 6). The rerun design SHALL enter no family.

#### Scenario: Family membership
- **WHEN** the families are assembled
- **THEN** every new design enters exactly one family with three readings

#### Scenario: A design short of units fails its family
- **WHEN** a design has fewer than 900 available units for a statistic
- **THEN** that reading is listed as unavailable and its family does not hold

### Requirement: The verdict follows a fixed order

The decision SHALL be `reproduction_failed` when the rerun design differs from its committed counts; otherwise `small_cells` when `balanced_small` fails; otherwise `small_unbalanced_cells` when `unbalanced_small` fails; otherwise `moderate_excess` when `reference` fails; otherwise `not_replicated`. Every family SHALL be recorded whatever the verdict, with the rule version.

#### Scenario: Small cells alone lose size
- **WHEN** the reproduction passes and `balanced_small` fails
- **THEN** the verdict is `small_cells`

#### Scenario: A small group in an unbalanced design loses size
- **WHEN** the reproduction passes, `balanced_small` holds and `unbalanced_small` fails
- **THEN** the verdict is `small_unbalanced_cells`

#### Scenario: The excess is not confined to small cells
- **WHEN** the reproduction passes, the two small-cell families hold and `reference` fails
- **THEN** the verdict is `moderate_excess`

#### Scenario: Nothing fails
- **WHEN** the reproduction passes and the three families hold
- **THEN** the verdict is `not_replicated`, and the decision reports the pooled `delta` count of the committed design and its fresh replicate

### Requirement: The group-dispersion families are re-read by a rule fixed in advance

The decision SHALL re-read the group-dispersion validation's four families from its committed decision, with no new unit, by the same count-plus-ceiling rule and that validation's verdict order. Which readings enter the re-read SHALL depend only on this run's verdict, as the scenarios state. The group-dispersion validation's own files SHALL NOT be changed.

#### Scenario: Small-cell designs are set aside
- **WHEN** the verdict is `small_cells` or `small_unbalanced_cells`
- **THEN** the re-read removes the five group-dispersion designs at n = 80 with shares 0.8 and 0.2, judges the four families on the remaining committed counts, and records the re-read verdict with the designs set aside

#### Scenario: The committed reading is replaced
- **WHEN** the verdict is `not_replicated`
- **THEN** the re-read judges `unbalanced_equal` with the fresh replicate's three readings in place of the committed design's, judges the other families in full, and records the re-read verdict

#### Scenario: No re-read
- **WHEN** the verdict is `moderate_excess` or `reproduction_failed`
- **THEN** the decision records that nothing is re-read and that the group-dispersion verdict stands

#### Scenario: The parent decision is bound
- **WHEN** the committed group-dispersion decision's sha256 differs from the settings'
- **THEN** `decide` exits naming the file

### Requirement: Recorded arms are labelled and not judged

The omnibus, `scan`, `fdr`, `total`, `span`, `complement`, `pca_matched`, the within-group statistics and the corrected statistics SHALL each be given a size label per family under the same rule, with every arm's rate per design. They SHALL NOT bear on the verdict or on the re-read.

#### Scenario: Recorded arms are outside the decision
- **WHEN** the decision is computed with a recorded arm's rates changed
- **THEN** the verdict and the re-read verdict are unchanged

### Requirement: Runs are resumable and bound to their settings

The driver SHALL provide `run`, `merge` and `decide`. Settings SHALL be written with a hash, and a resume or a merge under different settings SHALL exit naming the differing keys. Units SHALL be sharded round-robin, complete units SHALL be skipped on resume, and records SHALL hold seeds, statistics and counts only.

#### Scenario: Resume skips complete units
- **WHEN** a shard is rerun after an interruption
- **THEN** units whose records are complete are not evaluated again

#### Scenario: Other settings are refused
- **WHEN** `run` or `merge` is invoked on a directory whose stored settings differ
- **THEN** it exits naming the differing keys

### Requirement: Findings are reported with provenance and limits

The run SHALL produce a dated report and committed aggregates: settings, rate and cell tables, accounting, the reproduction, the decision and provenance. Every number the report quotes SHALL be readable from a committed aggregate, and the report SHALL state its limits beside the verdict.

#### Scenario: Numbers trace to aggregates
- **WHEN** the report quotes a rate, a count or a label
- **THEN** the value is in a committed CSV or JSON file under the run's results directory

#### Scenario: Limits are stated
- **WHEN** the report states the verdict or the re-read verdict
- **THEN** it states that the readings hold for simulated Gaussian noise at equal dispersion, δ = 2.0, three stages and no covariates, that no balanced design below 7 samples per stage was read, and that the re-read rule was written after the group-dispersion families' rates were known
- **AND** it states what the verdict changes in the documentation and the roadmap
