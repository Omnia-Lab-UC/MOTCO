# Spec Delta

## Purpose

Validates the per-feature interaction scan before a case study may use it: its size on semi-synthetic designs and on real noise, its familywise control when a few features truly differ, and its power for the sparse difference the other tests miss, on units that reproduce the split-test validation.

## ADDED Requirements

### Requirement: The simulation profile reruns the split-test validation's simulation units

The `simulation` profile SHALL use the split-test validation's simulation designs, cells, replicate counts, seeds and permutation count, read from its committed settings, whose sha256 the profile's settings record. The anchor and translation SHALL run replicates 0–499 and every other cell 0–99, all at 999 permutations.

#### Scenario: Units equal the parent's
- **WHEN** the simulation profile enumerates its units
- **THEN** its designs, cells, replicates, generator seeds and RRPP seeds equal the split-test validation's simulation units

#### Scenario: Changed parent settings are refused
- **WHEN** the committed parent settings file's sha256 differs from the one the profile's settings record
- **THEN** the profile exits naming the file and evaluates no unit

### Requirement: Each simulation unit runs the scan beside the split test at one seed

Every simulation unit SHALL evaluate the parent unit (the production fit and evaluation and the split test) unchanged. It SHALL run the scan on the features and the full and reduced models that the split test received, at the unit's RRPP seed and 999 permutations. A failure SHALL be recorded for the unit without aborting the shard.

#### Scenario: Scan equals the primitive on the unit's inputs
- **WHEN** a unit's scan record is compared with the scan applied to the same unit's split-test features, models and seed
- **THEN** the maximum F, its feature, the global p-value and the hit count are identical

#### Scenario: Record contents
- **WHEN** a unit completes
- **THEN** its record holds the parent unit's record fields
- **AND** it holds the scan's feature and degenerate counts, the maximum F and its block, the global p-value, the number of features with adjusted p ≤ 0.05 in total and per block, and the per-block share of null maxima
- **AND** where the dataset carries population stage means, it holds the number of hits on features whose population group difference is constant across stages

### Requirement: The plasmode profile reruns the split-test validation's plasmode units

The `plasmode` profile SHALL read the SEA-AD export through the plasmode's input contract, refuse an export whose file hashes differ from the plasmode's committed settings, and pass the plasmode's fidelity gate. It SHALL run the anchor, e = 2.0 and e = 3.0 at 500 replicates × 999 permutations with the plasmode's null-making step, plant and seeds (base seed 2034, `paper` stream).

#### Scenario: Another export is refused
- **WHEN** the export's file hashes differ from those in the plasmode's committed settings
- **THEN** the profile exits naming the differing files and runs no analysis

#### Scenario: Each unit is one case-study analysis with the scan
- **WHEN** a plasmode unit runs
- **THEN** it is one real-data analysis with the six covariates, rank floor 2, net displacement, the `df` residual correction, the harness CV defaults, the split test on, the scan on and attribution off
- **AND** each anchor unit passes the plasmode's latent check

#### Scenario: Planted-set counts are recorded
- **WHEN** a plasmode unit at a planted cell completes
- **THEN** its record holds the number of hits, the hits inside and outside the planted set, and the number of planted features among the 20 features with the largest F, in total and per block

#### Scenario: No donor-level data is recorded
- **WHEN** a plasmode record is written
- **THEN** it holds seeds, statistics and counts, and no donor identifier, label vector, feature list or donor-level value

### Requirement: Each profile must reproduce the split-test validation

For each profile the merge SHALL compare, per cell (and design), the split test's pseudo-F rejection counts per part and the trajectory statistics' rejection counts with the split-test validation's committed tables, each under that run's rejection rule. They SHALL match exactly.

#### Scenario: Reproduction holds
- **WHEN** every compared row matches the committed count and rejections
- **THEN** the merge records the reproduction as passed with the number of rows compared

#### Scenario: Reproduction fails
- **WHEN** any row differs
- **THEN** the merge records the reproduction as failed and lists every differing row

### Requirement: Runs are resumable, sharded and bound to their settings

Each profile SHALL write a settings file with a hash and one record per unit to its shard, skip complete units on a resume, and refuse a resume or merge whose settings differ, naming the differing keys.

#### Scenario: Resume skips completed units
- **WHEN** a shard is restarted after some units completed
- **THEN** only the remaining units are evaluated, and the merged output has one record per unit

#### Scenario: Settings mismatch is refused
- **WHEN** a resume or merge is attempted with settings whose hash differs from the stored one
- **THEN** the driver exits listing the differing keys

### Requirement: Merge writes the validation tables

Merging SHALL write the following, with the run's accounting (units, failures, duplicates, recorded core-hours):
- per profile, design and cell: the scan's global rejection rate, Monte Carlo SE and count, and the median hit count;
- for the plasmode's planted cells: the rate of units with a hit outside the planted set, and the mean recall, hit precision and top-20 precision;
- per block: the mean share of null maxima beside the block's share of features.

#### Scenario: Tables cover every cell
- **WHEN** a full run is merged
- **THEN** the scan table has a row for every profile, design and cell, and failed units are counted by cell rather than dropped

### Requirement: Reproduction gates the decision

`decide` SHALL read both profiles' merged tables and write `interaction_scan_decision.json`. When either profile's reproduction failed, the verdict SHALL be `reproduction_failed`, naming the failing rows, and no criterion SHALL be read.

#### Scenario: Failed reproduction
- **WHEN** the plasmode profile's reproduction failed
- **THEN** the verdict is `reproduction_failed` and the criteria are absent

### Requirement: Size is judged by the count-plus-ceiling rule

The scan's global rejection is adjusted p ≤ 0.05. Size criteria SHALL use the stage-span paper grade's rule:
- **count:** at most the 95th percentile of Binomial(m, p₂) rates above α + 2·√(α(1 − α)/R), with p₂ the exact Binomial(R, α) probability of exceeding that bound at the median R;
- **ceiling:** no rejected count above the Bonferroni ceiling at 0.05/m;
- **availability:** an unavailable rate fails.

#### Scenario: Simulation size (S1)
- **WHEN** S1 is evaluated
- **THEN** it judges the scan's global rates at the anchor and translation of all three designs together (m = 6)

#### Scenario: Real-noise size (S2)
- **WHEN** S2 is evaluated
- **THEN** it judges the scan's global rate at the plasmode anchor (m = 1)

### Requirement: Strong familywise control is judged under the plasmode's partial null

S3 SHALL judge, at the plasmode's e = 2.0 and e = 3.0 cells together (m = 2), the rate of units with at least one feature outside the planted set at adjusted p ≤ 0.05, under the same count-plus-ceiling rule.

#### Scenario: False hits within the bound
- **WHEN** the rates of units with a hit outside the planted set are within the rule at both cells
- **THEN** S3 holds, and the decision records both rates with their bounds

#### Scenario: False hits above the bound
- **WHEN** the rate at e = 3.0 exceeds the Bonferroni ceiling
- **THEN** S3 fails, naming the cell and its rate

### Requirement: Sparse gain is judged against every existing test

G SHALL hold when, at the plasmode's e = 3.0 cell, the scan's global rejection rate exceeds the largest rejection rate among `delta`, `angle`, `shape` and the split test's three parts on the same units by more than 2 pooled SE.

#### Scenario: Gain over the best existing test
- **WHEN** the scan's rate exceeds the largest existing rate by more than 2 pooled SE
- **THEN** G holds, and the decision records the comparator, both rates and the difference with its SE

#### Scenario: No gain
- **WHEN** the scan's rate exceeds the largest existing rate by less than 2 pooled SE
- **THEN** G fails, naming the comparator and the difference

### Requirement: The decision issues one verdict with its consequence

After a passed reproduction, the verdict SHALL be `size_excess` when S1 or S2 fails, else `fwer_excess` when S3 fails, else `no_sparse_gain` when G fails, else `scan_validated`, naming every failing rate. The decision SHALL state the consequence: only under `scan_validated` may a case study report the scan, and no verdict reopens Phase 6 or changes its attribution leg.

#### Scenario: Validated
- **WHEN** S1, S2, S3 and G hold
- **THEN** the verdict is `scan_validated`, and the consequence permits case-study use with the reading in the measurement-visibility documentation

#### Scenario: Size excess takes precedence
- **WHEN** S2 fails and S3 fails
- **THEN** the verdict is `size_excess`, naming the plasmode anchor's rate

#### Scenario: Familywise excess before gain
- **WHEN** S1 and S2 hold, S3 fails and G holds
- **THEN** the verdict is `fwer_excess`

### Requirement: The decision writes descriptive tables that never change the verdict

After a passed reproduction, the decision SHALL include descriptive tables, each row naming its source run and grade. No descriptive table SHALL enter a criterion.

#### Scenario: Descriptive tables are present
- **WHEN** the decision is written after a passed reproduction
- **THEN** it holds these tables:
  - **plasmode power and recovery:** the scan's global rate, recall, hit precision and top-20 precision at e = 2.0 and e = 3.0, beside the six existing tests' rates and the plasmode's committed planted-set precision of attribution;
  - **block balance:** per profile and cell, the share of null maxima in each block beside its share of features;
  - **simulation localization:** the scan's rate for every design and cell beside the split test's three parts and the production statistics on the same units;
  - **simulation false hits:** per cell, the rate of units with a hit on a feature whose population group difference is constant across stages.

#### Scenario: Descriptive tables do not change the verdict
- **WHEN** a descriptive rate changes while S1, S2, S3 and G keep their outcomes
- **THEN** the verdict is unchanged

### Requirement: The validation is reported in a dated report

A dated report under `docs/reports/` SHALL record:
- the configuration and provenance, the reproduction checks and the accounting;
- the four criteria and the verdict;
- plasmode power and recovery, block balance and simulation localization.

Aggregates SHALL be committed under `results/interaction-scan-<date>/`, and per-unit records SHALL stay off the repository.

#### Scenario: Report numbers match the aggregates
- **WHEN** the report is written
- **THEN** every rate, count and verdict it quotes can be read from the committed aggregates and decision
