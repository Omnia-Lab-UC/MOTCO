# interaction-scan-block-validation Specification

## Purpose
Validates the per-feature interaction scan's by-block combination against the joint scan, on the interaction-scan validation's simulated and real-noise units and on a real-noise export with blocks of unequal size, with a predeclared decision.

## Requirements

### Requirement: Every unit runs both combinations on the same inputs and draws

Every unit SHALL run the per-feature interaction scan twice on the same features, models, seed and residual correction: `joint`, without block labels, and `by_block`, with each feature's block as its label. No unit SHALL fit a latent space. A scan rejects globally, and a feature is a hit, at adjusted p ≤ 0.05.

#### Scenario: Paired arms
- **WHEN** a unit is evaluated
- **THEN** both arms' per-feature F are identical, and the arms differ only in their p-values and hits

#### Scenario: A failing unit is recorded
- **WHEN** a unit raises an error
- **THEN** it is recorded as failed with the error type and message, and the shard continues

### Requirement: The simulation profile reruns the interaction-scan validation's simulation units

The `simulation` profile SHALL rebuild the interaction-scan validation's simulation settings and require them to equal the committed ones, naming any differing key. Each unit SHALL scan the harness's pooled-standardized features on the measured blocks with the trajectory models and the unit's RRPP seed, without a residual correction.

#### Scenario: Settings match the parent
- **WHEN** the profile's settings are built
- **THEN** they equal `results/interaction-scan-2026-11-09/simulation/settings.json` in designs, cells, replicates, permutations and seeds, or the driver exits naming the differing keys

#### Scenario: Records hold per-block counts
- **WHEN** a simulation unit completes
- **THEN** its record gives, per arm, the global p-value, the hits per block and the hits on population-null features per block, and the joint arm's per-block null-maximum shares

### Requirement: The plasmode profile reruns the interaction-scan validation's plasmode units

The `plasmode` profile SHALL rebuild the interaction-scan validation's plasmode settings and require them to equal the committed ones, and SHALL require the export's file hashes to equal the plasmode's committed settings, naming any differing file. Each unit SHALL scan the analysis's own scan inputs at the parent's parameters, with `df`, and SHALL record per arm the planted hits and the hits outside the planted set per block.

#### Scenario: Export hashes are checked
- **WHEN** an export file differs from the plasmode's committed hash
- **THEN** the run exits naming the file before any unit

#### Scenario: The scan inputs are the analysis's
- **WHEN** a plasmode unit builds its scan inputs
- **THEN** a joint scan on them equals the scan `analyze_trajectories` runs on the same unit at the parent's parameters

### Requirement: The unequal profile reduces one block of the real export

The `plasmode_unequal` profile SHALL reduce the export's RNA block to 500 of its 5,000 features, chosen by a seeded permutation fixed for the profile and recorded by the sha256 of the feature list, and SHALL run the plasmode's anchor, e = 2.0 and e = 3.0 at 500 replicates and 999 permutations, with 10 planted features per block drawn from the reduced export.

#### Scenario: The subset is fixed and recorded
- **WHEN** the profile's settings are built twice
- **THEN** both name the same 500 RNA features by the same sha256, and every unit uses that subset

#### Scenario: Planted features come from the reduced block
- **WHEN** a unit of the e = 2.0 cell is built
- **THEN** its 10 planted RNA features are among the 500 retained ones

### Requirement: Records hold counts and seeds only

No record SHALL hold a feature name, donor identifier, label vector or donor-level value. Each unit SHALL also record the 20 largest F's planted count (plasmode profiles), the correlation of the block maxima across draws, and the global p-value a minP combination between blocks would give.

#### Scenario: Record contents
- **WHEN** a plasmode unit completes
- **THEN** its record holds seeds, statistics and counts, and no feature or donor name appears in it

### Requirement: Runs are resumable, sharded and bound to their settings

Each profile SHALL write a settings file with a hash and one record per unit to its shard, skip complete units on a resume, and refuse a resume or merge whose settings differ, naming the differing keys. The sha256 of every parent file SHALL be recorded, and `run`, `merge` and `decide` SHALL refuse a changed parent file by name.

#### Scenario: Resume skips completed units
- **WHEN** a shard is restarted after some units completed
- **THEN** only the remaining units are evaluated, and the merged output has one record per unit

#### Scenario: A changed parent file is refused
- **WHEN** a parent table's sha256 differs from the recorded one
- **THEN** the driver exits naming the file

### Requirement: The joint arm must reproduce the interaction-scan validation

Merging the `simulation` and `plasmode` profiles SHALL compare the `joint` arm with the parent's committed tables and write the result. The comparison SHALL be exact, in count, rejections and hits for the scan table, in mean hits and null-maximum share for the block-balance table, and in the false-hit or recovery table's counts, rejections and means.

#### Scenario: Exact reproduction
- **WHEN** both profiles are merged after a complete run
- **THEN** every compared row equals the parent's, and the reproduction passes

#### Scenario: A mismatch is named
- **WHEN** one cell's rejection count differs from the parent's
- **THEN** the reproduction fails, naming the profile, the cell and both counts

### Requirement: Merge writes the validation tables

Merging SHALL write the run's accounting (units, failures, duplicates, recorded core-hours) and per-cell, per-arm tables of global rejections, per-block hits, hits outside the planted set, paired discordance and the minP comparison.

#### Scenario: Tables cover every cell and arm
- **WHEN** a full run is merged
- **THEN** the tables have a row for every profile, design, cell and arm, and failed units are counted by cell rather than dropped

#### Scenario: Table contents
- **WHEN** a full run is merged
- **THEN** the tables give per profile, design, cell and arm the global rejection rate, Monte Carlo SE, count and median hit count
- **AND** per block and arm the rate of units with a hit in the block and the mean hits, with recall and hit precision in the plasmode profiles
- **AND** per planted cell and arm the rate of units with a hit outside the planted set
- **AND** per cell the paired discordant counts of the two arms' global rejections, and the minP and Bonferroni global rates with the median correlation of block maxima

### Requirement: Reproduction gates the decision

`decide` SHALL read the three profiles' merged tables and write `interaction_scan_blocks_decision.json` (rule version 1). When a reproduction failed, the verdict SHALL be `reproduction_failed`, naming the failing rows, and no criterion SHALL be read.

#### Scenario: Failed reproduction
- **WHEN** the simulation profile's reproduction failed
- **THEN** the verdict is `reproduction_failed` and the criteria are absent

### Requirement: Size and familywise control are judged by the count-plus-ceiling rule

The by-block arm SHALL be judged by the stage-span paper grade's count-plus-ceiling rule on: S1, the global rate at the anchor and translation of every simulation design (m = 6); S2, the global rate at both exports' anchors (m = 2); S3, the rate of units with a hit outside the planted set at e = 2.0 and 3.0 on both exports (m = 4).

#### Scenario: Size holds
- **WHEN** every S1 and S2 rate is within the rule
- **THEN** S1 and S2 hold, and each rate is reported beside its bound

#### Scenario: Familywise excess
- **WHEN** the unequal export's rate of units with a hit outside the planted set at e = 3.0 is above the Bonferroni ceiling
- **THEN** S3 fails, naming that cell

### Requirement: Power is protected and the gain is judged on the unequal export

P SHALL hold when, at e = 2.0 and 3.0 on both exports, no exact one-sided paired loss p-value of the by-block arm's global rejections against the joint arm's is below 0.05/4. G SHALL hold when, on the unequal export, the exact one-sided sign test of the per-unit difference in planted hits in the reduced block (by block minus joint) gives p ≤ 0.025 at e = 2.0 or at e = 3.0.

#### Scenario: Power loss
- **WHEN** at e = 2.0 on the balanced export 30 units reject only jointly and 5 only by block
- **THEN** the loss p-value is below 0.05/4 and P fails, naming the cell and both counts

#### Scenario: Allocation gain
- **WHEN** on the unequal export at e = 2.0 more units gain planted reduced-block hits by block than lose them, with a sign-test p-value of at most 0.025
- **THEN** G holds

### Requirement: The decision issues one verdict with its consequence

After a passed reproduction, the verdict SHALL be `size_excess` when S1 or S2 fails, else `fwer_excess` when S3 fails, else `power_loss` when P fails, else `no_allocation_gain` when G fails, else `by_block_validated`, naming every failing rate. The consequence SHALL state that the joint scan stays the default under every verdict.

#### Scenario: Validated
- **WHEN** S1, S2, S3, P and G hold
- **THEN** the verdict is `by_block_validated`, and the consequence documents the option as validated for blocks of unequal size or noise, with the joint scan the default

#### Scenario: Size excess takes precedence
- **WHEN** S2 fails and P fails
- **THEN** the verdict is `size_excess`

### Requirement: The decision writes descriptive tables that never change the verdict

After a passed reproduction, the decision SHALL include descriptive tables, each row naming its source run and grade: block balance per arm, recall and hit precision per block per arm beside the parent's committed recovery, the minP and Bonferroni global rates with the correlation of block maxima, and simulation localization and population-null false hits per arm. No descriptive table SHALL enter a criterion.

#### Scenario: Descriptive tables do not change the verdict
- **WHEN** a descriptive rate changes while S1, S2, S3, P and G keep their outcomes
- **THEN** the verdict is unchanged

### Requirement: The validation is reported in a dated report

A dated report under `docs/reports/` SHALL record the configuration and provenance, the reproduction and accounting, the five criteria and the verdict, and the descriptive tables. Aggregates SHALL be committed under `results/interaction-scan-blocks-<date>/`, and per-unit records SHALL stay off the repository.

#### Scenario: Report numbers match the aggregates
- **WHEN** the report is written
- **THEN** every rate, count and verdict it quotes can be read from the committed aggregates and decision
