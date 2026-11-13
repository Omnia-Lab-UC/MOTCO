# Spec Delta

## Purpose

Diagnoses why the predictive-subspace coordinates exceed size on the SEA-AD plasmode anchor, by re-reading the predictive-subspace pilot's plasmode units under within-stage permutation draws and under out-of-fold coordinates, with a predeclared decision that says which repair, if any, restores size.

## ADDED Requirements

### Requirement: The diagnosis re-reads the pilot's plasmode units

The diagnosis SHALL evaluate the plasmode units of the predictive-subspace pilot: the anchor at the `paper` stream's replicates 0–499 and the real-noise span plant at e ∈ {0.25, 0.5, 1.0} at replicates 0–99, at 999 permutations and the case-study configuration, 800 units. Each unit SHALL use the pilot's label draw, plant and RRPP seed. The export SHALL be read in place through the plasmode's input contract.

#### Scenario: Units equal the pilot's
- **WHEN** the units are enumerated
- **THEN** their cells, replicates and seeds equal those of the pilot's plasmode profile

#### Scenario: A changed export file is refused
- **WHEN** an export file's hash differs from the one in the plasmode's committed settings
- **THEN** the run exits naming the file before any unit is evaluated

#### Scenario: Records hold no donor data
- **WHEN** a unit is recorded
- **THEN** the record holds seeds and statistics only, with no donor identifier, no label vector and no fold assignment

#### Scenario: Failed unit
- **WHEN** building or fitting a unit raises
- **THEN** the unit is recorded as failed with the error and the shard continues

### Requirement: Every in-sample arm of a unit measures the unit's production fit

Every in-sample arm of a unit SHALL be measured on that unit's production stage-PLS model, at one rank. An anchor unit SHALL fit at the rank of the plasmode's committed fidelity values without cross-validation. A plant unit SHALL run the production cross-validation once and fit every other arm at the rank it selects.

#### Scenario: Anchor unit
- **WHEN** an anchor unit runs
- **THEN** no cross-validation is made and every arm reports the committed fidelity rank

#### Scenario: Same fit
- **WHEN** a predictive or covariate-free arm is fitted
- **THEN** its fitted model equals the unit's production model: the same rank, and rotations and response loadings equal to 1e-8

#### Scenario: Plant unit
- **WHEN** a plant unit runs
- **THEN** exactly one cross-validation is made and every other arm is fitted at the rank it selected

#### Scenario: Pre-flight
- **WHEN** the cross-validated production analysis is run on anchor replicates 0 and 1 before the full run
- **THEN** it selects the committed fidelity rank and AUROC, and its statistics and p-values equal the forced-rank analysis's exactly
- **AND** a mismatch stops the run

### Requirement: Within-stage draws permute residuals only among donors of one stage

The within-stage null SHALL permute the reduced model's residuals among the rows of each stage and never across stages. It SHALL keep the reduced fit, the residual correction, the statistics, the unit's RRPP seed and the p-value rule of the unrestricted null.

#### Scenario: Strata are respected
- **WHEN** a within-stage draw is made
- **THEN** every row receives the residual of a row of its own stage

#### Scenario: One stratum is the unrestricted null
- **WHEN** every row belongs to one stratum
- **THEN** the null distributions equal those of the package's permutation procedure at the same seed exactly

#### Scenario: The residual correction is kept
- **WHEN** the within-stage null runs with the `df` correction
- **THEN** the permuted residuals carry the same multiplier the unrestricted null applies

### Requirement: Out-of-fold coordinates are fitted without the donors they measure

The out-of-fold coordinates SHALL be the fitted stage indicators of stage-PLS models fitted on stage-stratified folds: each donor's row comes from the model fitted on the folds that exclude it, at the unit's production rank and scaling, in that model's response units. Four folds SHALL be drawn from a seed derived per replicate, and the seed SHALL be recorded.

#### Scenario: Every donor is held out once
- **WHEN** the out-of-fold matrix is built
- **THEN** every donor is in exactly one held-out fold and in no training set that produced its row
- **AND** every training set contains every stage

#### Scenario: All-donor limit
- **WHEN** the same computation is made with every donor in the training set
- **THEN** it equals the in-sample `predictive_fitted` coordinates to 1e-8

#### Scenario: Reproducible folds
- **WHEN** the same replicate is evaluated twice
- **THEN** the fold seed and the out-of-fold matrix are the same

### Requirement: Each unit evaluates reference, check and descriptive arms

Each unit SHALL evaluate four reference arms (`pls_cv`, `predictive`, `predictive_orthonormal`, `predictive_fitted` under unrestricted draws), four decisive arms (the three predictive coordinate systems under within-stage draws, and the out-of-fold coordinates under unrestricted draws) and six descriptive arms. An arm that raises SHALL be recorded as failed on its unit without failing the unit.

#### Scenario: Decisive arms
- **WHEN** a unit completes
- **THEN** it holds `predictive__within_stage`, `predictive_orthonormal__within_stage`, `predictive_fitted__within_stage` and `fitted_out_of_fold`, each with its observed `delta`, `angle` and `shape`, their p-values and the `angle` null 95th percentile

#### Scenario: Descriptive arms
- **WHEN** a unit completes
- **THEN** it also holds `pls_cv__within_stage`, `fitted_out_of_fold__within_stage` and the four reference coordinate systems analysed with no covariates in the trajectory model

#### Scenario: Within-stage arms share their draws
- **WHEN** the within-stage arms of one unit are evaluated
- **THEN** they use the same sequence of restricted permutations

#### Scenario: A failed arm
- **WHEN** one arm raises
- **THEN** the arm is recorded as failed with the error and the unit's other arms are evaluated

### Requirement: Each unit checks its own inference path against the analysis

Each unit SHALL recompute, through the path the check arms use, the statistics and unrestricted p-values of the four reference arms, and SHALL compare them with the real-data analysis's. It SHALL also check that every predictive reference arm fitted the production model. A unit failing a check SHALL be recorded as inconsistent, naming the check, the arm and the deviation, and SHALL be excluded from every rate.

#### Scenario: Path identity holds
- **WHEN** the path is run on a reference arm's latent matrix under unrestricted draws
- **THEN** its statistics equal the analysis's to 1e-10 and its p-values equal them exactly

#### Scenario: Inconsistent unit
- **WHEN** a check fails
- **THEN** the unit is recorded as inconsistent with the check's name and deviation, and it is excluded from every rate

### Requirement: The reference arms must reproduce the pilot

The merge SHALL compare the four reference arms' available and rejected counts, per cell and statistic, with the pilot's committed plasmode rate table, and SHALL write the comparison with every differing row and both values.

#### Scenario: Reproduction holds
- **WHEN** every reference count equals the committed one in all four cells
- **THEN** the reproduction is recorded as passed

#### Scenario: Reproduction fails
- **WHEN** a count differs
- **THEN** the reproduction is recorded as failed, listing the cell, arm, statistic and both counts

### Requirement: Runs are resumable, sharded and bound to their settings and parents

A run SHALL write its settings with their hash and the sha256 of every committed parent file it reads. A resume or merge under other settings SHALL exit naming the differing keys, and a changed parent file SHALL be refused by name. Shards SHALL partition the units round-robin, and a resume SHALL skip complete units.

#### Scenario: Resume skips completed units
- **WHEN** a shard is rerun
- **THEN** units with a complete record are not evaluated again

#### Scenario: Settings mismatch is refused
- **WHEN** a run resumes in a directory whose settings differ
- **THEN** it exits naming the differing keys

#### Scenario: A changed parent is refused
- **WHEN** a committed parent file's hash differs from the recorded one
- **THEN** `run`, `merge` and `decide` exit naming the file

### Requirement: Merge writes the rate, paired and cell tables

The merge SHALL write, per cell, a rate table (per arm and statistic: available and rejected counts, rate, Monte Carlo SE, failures), a paired table (each non-reference arm against its counterpart and against `pls_cv`: both rates, the discordant counts and exact one-sided and two-sided p-values) and a cell table, with an accounting of units, failures and inconsistent units.

#### Scenario: Counterparts
- **WHEN** the paired table is written
- **THEN** a within-stage in-sample arm or a covariate-free arm is paired with the same coordinate system under unrestricted draws and six covariates, and each out-of-fold arm with the `predictive_fitted` arm under the same draws

#### Scenario: No discordant units
- **WHEN** an arm and its counterpart reject on the same units
- **THEN** every paired p-value is 1

#### Scenario: Dispersion by stage
- **WHEN** the cell table is written
- **THEN** it holds, per coordinate system, the median over units of each stage's residual dispersion relative to the pooled one and of the largest-to-smallest ratio
- **AND** the in-sample stage R² beside the share of the stage indicators' variance the out-of-fold predictions explain

### Requirement: The decision grades each decisive arm on the anchor

The decision SHALL grade each decisive arm from its anchor readings. `restored`: its `delta`, `angle` and `shape` rejections hold the count-plus-ceiling rule of the stage-span paper grade over those three readings. `reduced`: not `restored`, and the one-sided exact paired p-value for fewer rejections than its counterpart is at most 0.025 on `delta` or on `shape`. `unchanged`: otherwise.

#### Scenario: Restored
- **WHEN** an arm's three anchor readings hold the count-plus-ceiling rule
- **THEN** the arm is graded `restored`, with its counts, bound and ceiling recorded

#### Scenario: Reduced
- **WHEN** an arm's readings do not hold the rule and its paired p-value on `delta` is at most 0.025
- **THEN** the arm is graded `reduced`, naming the statistic and the discordant counts

#### Scenario: Unchanged
- **WHEN** an arm's readings do not hold the rule and neither paired p-value is at most 0.025
- **THEN** the arm is graded `unchanged`

#### Scenario: Unavailable reading
- **WHEN** a decisive arm has no available anchor reading for a statistic
- **THEN** the arm is graded `unchanged` and the reading is named

### Requirement: The decision writes one predeclared verdict

The decision SHALL write `predictive_size_excess_decision.json` with a rule version. The within-stage check SHALL take the weakest grade of its three arms and the out-of-fold check the grade of its arm. The verdict SHALL be `both_restore`, `null_restores`, `out_of_fold_restores`, `partly_explained` or `unexplained`, with the consequence statement fixed for it, and `reproduction_failed` SHALL gate all of them.

#### Scenario: Both checks restore
- **WHEN** both checks are `restored`
- **THEN** the verdict is `both_restore`

#### Scenario: Only the draws restore
- **WHEN** the within-stage check is `restored` and the out-of-fold check is not
- **THEN** the verdict is `null_restores`

#### Scenario: Only the out-of-fold coordinates restore
- **WHEN** the out-of-fold check is `restored` and the within-stage check is not
- **THEN** the verdict is `out_of_fold_restores`

#### Scenario: Partly explained
- **WHEN** neither check is `restored` and at least one is `reduced`
- **THEN** the verdict is `partly_explained`, naming the reduced check

#### Scenario: Unexplained
- **WHEN** both checks are `unchanged`
- **THEN** the verdict is `unexplained`

#### Scenario: One metric is not enough
- **WHEN** two within-stage arms are `restored` and the third is `unchanged`
- **THEN** the within-stage check is `unchanged`

#### Scenario: Gate
- **WHEN** the reproduction failed, a unit is inconsistent or the pre-flight was not passed
- **THEN** the verdict is `reproduction_failed`, and the grades are still recorded

### Requirement: The decision reports descriptive readings that never change the verdict

Beside the verdict the decision SHALL report: every arm's plant-cell rates with its paired gain against `pls_cv` and paired loss against its counterpart; `pls_cv__within_stage` and `fitted_out_of_fold__within_stage` on the anchor; the covariate-free arms against their counterparts; the dispersion by stage per coordinate system; and each arm's median `angle` null 95th percentile.

#### Scenario: Plant cells under a restoring arm
- **WHEN** an arm is graded `restored`
- **THEN** the decision states its plant-cell `delta` rates beside production's and the reference arm's

#### Scenario: Descriptive readings do not decide
- **WHEN** a descriptive arm exceeds or holds the size rule
- **THEN** the verdict is the one the four decisive arms give

### Requirement: The diagnosis is recorded in a dated report

The run SHALL be recorded in a dated findings report and committed aggregates: the settings, the three tables, the accounting, the reproduction, the decision and the provenance. The roadmap SHALL carry the verdict.

#### Scenario: Report contents
- **WHEN** the report is written
- **THEN** it states the verdict, each decisive arm's grade and counts, the reproduction, every descriptive reading and what follows from the verdict
- **AND** every number it quotes is in a committed aggregate

#### Scenario: Roadmap
- **WHEN** the run is recorded
- **THEN** roadmap open items 1 and 2 and the Done list state the verdict and its consequence
