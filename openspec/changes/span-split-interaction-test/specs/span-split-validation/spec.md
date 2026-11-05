# Spec Delta

## Purpose

Validates the split interaction test before a case study may use it: its size on semi-synthetic designs and on real noise, its power for a difference outside the fitted stage span, and where it places differences relative to that span, on units that reproduce two committed runs.

## ADDED Requirements

### Requirement: The simulation profile reuses the latent-space comparison's units

The `simulation` profile SHALL use the latent-space comparison's three designs, fourteen cells, base seed 2035 and its generator and RRPP seed derivations. The weak design's δ SHALL come from the comparison's committed calibration, whose sha256 the settings record. The anchor and translation cells SHALL run 500 replicates, every other cell 100, all at 999 permutations.

#### Scenario: Units and grade
- **WHEN** the simulation profile enumerates its units
- **THEN** each design has the anchor and translation at replicates 0–499 and the other twelve cells at replicates 0–99, all at 999 permutations
- **AND** the weak design's δ equals the comparison's adopted δ

#### Scenario: Shared units are the comparison's datasets
- **WHEN** a unit at replicate 0–99 is generated
- **THEN** its dataset and RRPP seed equal the latent-space comparison's for the same design, cell and replicate

### Requirement: Each simulation unit evaluates the production arm and the split test at one seed

Every simulation unit SHALL fit the production stage-PLS space once and evaluate it with net displacement, rank floor 2 and span capture. It SHALL run the split test on that fit's standardized features, in the fit's input coordinates and against its orthonormal span basis, with the full and reduced trajectory models, at the unit's RRPP seed and 999 permutations. A failure SHALL be recorded for the unit without aborting the shard.

#### Scenario: Split test equals the primitive on the unit's fit
- **WHEN** a unit's split record is compared with the split test applied to the same dataset's production fit, models and seed
- **THEN** every part's sums of squares, pseudo-F and p-values are identical

#### Scenario: Record contents
- **WHEN** a unit completes
- **THEN** its record names the design, cell, replicate and both seeds
- **AND** it holds the production arm's rank, held-out AUROC, statistics, p-values at 999 permutations, p-values from its first 199 draws and span capture
- **AND** per split part it holds the dimension, sums of squares, pseudo-F, the pseudo-F p-value and the interaction-SS p-value both uncorrected and under `df`, together with the observed in-span share

### Requirement: The simulation profile must reproduce the comparison

On the units the profile shares with the latent-space comparison, the production arm's rejection counts SHALL be computed from the first 199 draws of each unit's permutation sequence. For every design, cell and statistic they SHALL equal the comparison's committed `rate_table.csv` exactly.

#### Scenario: Reproduction holds
- **WHEN** every shared design, cell and statistic matches the committed count and rejections
- **THEN** the merge records the reproduction as passed with the number of rows compared

#### Scenario: Reproduction fails
- **WHEN** any row differs
- **THEN** the merge records the reproduction as failed and lists every differing row

### Requirement: The plasmode profile reuses the plasmode demonstration's units

The `plasmode` profile SHALL read the SEA-AD export through the plasmode's input contract and refuse an export whose file hashes differ from the plasmode's committed settings. It SHALL pass the plasmode's fidelity gate, then run the anchor, e = 2.0 and e = 3.0 at 500 replicates × 999 permutations with the plasmode's null-making step, plant and seeds (base seed 2034, `paper` stream).

#### Scenario: Another export is refused
- **WHEN** the export's file hashes differ from those in the plasmode's committed settings
- **THEN** the profile exits naming the differing files and runs no analysis

#### Scenario: Each unit is one case-study analysis
- **WHEN** a plasmode unit runs
- **THEN** it is one real-data analysis with the six covariates, rank floor 2, net displacement, the `df` residual correction, the harness CV defaults, the split test on and attribution off
- **AND** each anchor unit passes the plasmode's latent check

#### Scenario: No donor-level data is recorded
- **WHEN** a plasmode record is written
- **THEN** it holds seeds and statistics, and no donor identifier, label vector or donor-level value

### Requirement: The plasmode profile must reproduce the plasmode demonstration

For each plasmode cell the trajectory statistics' rejection rates (`delta`, `angle`, `shape`) SHALL equal the plasmode's committed `plasmode_summary.csv` exactly.

#### Scenario: Plasmode reproduction
- **WHEN** the plasmode profile is merged
- **THEN** the merge records, per cell and statistic, the committed and reproduced rates and whether every one matches

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
- per profile, cell (and design) and split part: the rejection rate, Monte Carlo SE and count for the pseudo-F and each interaction-SS diagnostic, and the median observed in-span share;
- per cell: the production arm's rates;
- for the simulation profile, the median population capture.

#### Scenario: Tables cover every cell and part
- **WHEN** a full run is merged
- **THEN** the split table has a row for every profile, design, cell, part and statistic, and failed units are counted by cell rather than dropped

### Requirement: Reproduction gates the decision

`decide` SHALL read both profiles' merged tables and write `span_split_decision.json`. When either profile's reproduction failed, the verdict SHALL be `reproduction_failed`, naming the failing rows, and no criterion SHALL be read.

#### Scenario: Failed reproduction
- **WHEN** the simulation profile's reproduction failed
- **THEN** the verdict is `reproduction_failed` and the criteria are absent

### Requirement: Size is judged by the count-plus-ceiling rule

Size criteria SHALL use the pseudo-F p-values at α = 0.05 under the count-plus-ceiling rule:
- **count:** at most the 95th percentile of Binomial(m, p₂) rates above α + 2·√(α(1 − α)/R), with p₂ the exact Binomial(R, α) probability of exceeding that bound at the median R;
- **ceiling:** no rejected count above the Bonferroni ceiling at 0.05/m;
- **availability:** an unavailable rate fails.

#### Scenario: Simulation size (C1)
- **WHEN** C1 is evaluated
- **THEN** it judges the `total`, `span` and `complement` rates at the anchor and translation of all three designs together (m = 18)

#### Scenario: Real-noise size (C2)
- **WHEN** C2 is evaluated
- **THEN** it judges the three parts' rates at the plasmode anchor together (m = 3)

### Requirement: Out-of-span gain is judged against the production arm

C3 SHALL hold when, at both δ = 2.0 designs, the `complement` pseudo-F rejection rate for orientation e = 0.5 at φ = 90 with dense support exceeds the production arm's `angle` rate on the same units by more than 2 pooled SE.

#### Scenario: Gain at both designs
- **WHEN** the complement rate exceeds the production `angle` rate by more than 2 pooled SE at `n200` and at `seaad`
- **THEN** C3 holds, and the decision records both differences with their SE

#### Scenario: Gain at one design only
- **WHEN** the gain holds at `n200` but not at `seaad`
- **THEN** C3 fails, naming `seaad`

### Requirement: The decision issues one verdict with its consequence

After a passed reproduction, the verdict SHALL be `size_excess` when C1 or C2 fails (naming every failing rate), else `no_out_of_span_gain` when C3 fails, else `split_validated`. The decision SHALL state the consequence: only under `split_validated` may a case study report the split test, with the documented reading of its parts.

#### Scenario: Validated
- **WHEN** C1, C2 and C3 hold
- **THEN** the verdict is `split_validated`, and the consequence permits case-study use with the reading in the measurement-visibility documentation

#### Scenario: Size excess takes precedence
- **WHEN** C2 fails and C3 holds
- **THEN** the verdict is `size_excess`, naming the failing parts at the plasmode anchor

### Requirement: The decision writes descriptive tables that never change the verdict

After a passed reproduction, the decision SHALL include descriptive tables of localization, power and the interaction-SS diagnostic, each row naming its source run and grade. No descriptive table SHALL enter a criterion.

#### Scenario: Descriptive tables are present
- **WHEN** the decision is written after a passed reproduction
- **THEN** it holds these tables:
  - **localization:** `span` and `complement` rates for every cell and design;
  - **share against capture:** the median observed in-span share beside the production arm's median population capture;
  - **span part against production:** the `span` part's rates beside the production arm's three statistics;
  - **against the comparison:** the `complement` rates beside the comparison's committed `concat` and `pls_cv` `angle` rates;
  - **plasmode power:** the three parts' rates at e = 2.0 and e = 3.0 beside the plasmode's trajectory rates;
  - **diagnostic:** the interaction-SS sizes, uncorrected and under `df`, by part and design.

#### Scenario: Descriptive tables do not change the verdict
- **WHEN** a descriptive rate changes while C1, C2 and C3 keep their outcomes
- **THEN** the verdict is unchanged

### Requirement: The validation is reported in a dated report

A dated report under `docs/reports/` SHALL record:
- the configuration and provenance, the reproduction checks and the accounting;
- the three criteria and the verdict;
- localization across φ and modes, plasmode power and the diagnostic.

Aggregates SHALL be committed under `results/span-split-<date>/`, and per-unit records SHALL stay off the repository.

#### Scenario: Report numbers match the aggregates
- **WHEN** the report is written
- **THEN** every rate, share and verdict it quotes can be read from the committed aggregates and decision
