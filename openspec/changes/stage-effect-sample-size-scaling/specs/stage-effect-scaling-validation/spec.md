# Spec Delta

## Purpose

Reads how a weaker stage effect trades against sample size for the production trajectory test, on simulated designs matched in δ√n to committed ones and on designs calibrated to the case studies' held-out stage AUROC. A predeclared decision states whether a weak stage effect costs more than sample size buys back, and a grade per design states where a rejection names the kind of a difference.

## ADDED Requirements

### Requirement: Each unit is the univariate baseline's unit

Every unit SHALL be evaluated as the univariate baseline evaluates one: one dataset, one production stage-PLS fit with net displacement and rank floor 2, and every arm read at the unit's RRPP seed. A failure SHALL be recorded for the unit without aborting the shard.

#### Scenario: One cross-validation per unit
- **WHEN** a unit runs
- **THEN** exactly one cross-validated PLS fit is made, and every arm's record equals the univariate baseline's record for the same dataset, parameters and seed

#### Scenario: Failed unit
- **WHEN** generation, fitting or an arm fails for a unit
- **THEN** the unit is recorded as failed with the error, and the shard continues

### Requirement: Three arms are judged and the others recorded

The judged arms SHALL be `production`, which rejects when the smallest of the `delta`, `angle` and `shape` p-values is at or below 0.05 / 3, and the per-feature arms `fdr` and `total`. `scan`, the three single statistics, `span`, `complement`, `production_minp` and `pca_matched` SHALL be recorded and SHALL NOT bear on the verdict.

#### Scenario: Each arm keeps its recorded rule
- **WHEN** a unit's rejections are recorded
- **THEN** each arm uses the rejection rule the univariate baseline recorded for it

#### Scenario: Recorded arms are outside the decision
- **WHEN** the decision is computed with a recorded arm's rates changed
- **THEN** the verdict is unchanged

### Requirement: The saturated profile's new designs are matched in δ√n to committed designs

The `saturated` profile SHALL add four designs on methylation, expression and proteomics, three stages and equal groups, each differing from a committed design of the univariate baseline only in stage effect δ and sample size n, with the same δ√n: δ = 1.0 at n = 320 and δ = 0.5 at n = 1,280 against δ = 2.0 at n = 80; δ = 1.0 at n = 800 against δ = 2.0 at n = 200; and δ = 0.25 at n = 800 against δ = 0.5 at n = 200.

#### Scenario: Matched designs
- **WHEN** a new design's generator parameters are compared with its match's
- **THEN** they differ only in the per-omic stage effects and the sample size, and δ√n is equal for the two

#### Scenario: Cells and grade
- **WHEN** the profile enumerates a new design's units
- **THEN** it has the univariate baseline's thirteen stage-span cells with their scales, the anchor and translation at 500 replicates and the eleven other cells at 200, at 999 permutations

#### Scenario: New datasets
- **WHEN** a new design's dataset seeds are derived
- **THEN** none equals a dataset seed of the univariate baseline, and cells of one design and replicate share their null dataset

### Requirement: The committed weak-stage design is rerun and must reproduce

The `saturated` profile SHALL rerun the univariate baseline's design at δ = 0.5 and n = 200 with that run's seeds, cells and grade. For every cell and arm the available and rejected counts SHALL equal the committed rate table exactly. The committed settings, rate table and cell table SHALL be bound by sha256 and refused by name if changed.

#### Scenario: Reproduction holds
- **WHEN** every row of the rerun design matches the committed counts
- **THEN** the merge records the reproduction as passed with the number of rows compared

#### Scenario: Reproduction fails
- **WHEN** any row differs
- **THEN** the merge records the reproduction as failed and lists every differing row with both counts

#### Scenario: A changed parent is refused
- **WHEN** a bound parent file's sha256 differs from the one in the stored settings
- **THEN** the driver exits naming the file and runs no unit

### Requirement: The band's stage effect is calibrated to the held-out AUROC

A calibration step SHALL fit the production space, with no permutation, on 40 anchor datasets per ladder value of δ at n = 320, and SHALL choose the value whose median held-out stage AUROC is nearest the latent-space comparison's target, 0.636. The chosen median SHALL lie in that comparison's band, 0.60–0.67; otherwise the `band` profile SHALL NOT be built.

#### Scenario: A value inside the band
- **WHEN** the ladder value nearest the target has a median AUROC inside the band
- **THEN** the calibration records it as δ_c with every ladder value's AUROC and rank quantiles

#### Scenario: No value inside the band
- **WHEN** the nearest ladder value's median AUROC is outside the band
- **THEN** the driver refuses to build the `band` settings and names the ladder's readings

#### Scenario: The run is bound to its calibration
- **WHEN** the `band` profile is run or merged
- **THEN** its settings carry the calibration's hash, and a changed calibration is refused

### Requirement: The band profile reads three designs below the AUROC ceiling

The `band` profile SHALL run three designs on the same blocks, stages and groups: δ_c at n = 320; δ_c / 2 at n = 1,280, its matched design; and δ_c at n = 1,280. Each SHALL have eleven cells under the planted scale: the anchor, a translation at e = 1.0, stretches at e ∈ {0.5, 1.0, 2.0}, and in-span rotations and bends at e ∈ {0.25, 0.5, 1.0}.

#### Scenario: Units and grade
- **WHEN** the profile enumerates its units
- **THEN** each design has its two controls at 500 replicates and its nine other cells at 200, 8,400 units in all, at 999 permutations

#### Scenario: The matched design is not calibrated
- **WHEN** the design at δ_c / 2 and n = 1,280 is merged
- **THEN** its held-out AUROC is reported beside the calibrated design's and is not required to lie in the band

### Requirement: Runs are resumable, sharded and bound to their settings

Each profile SHALL write a settings file with a hash and one record per unit to its shard. It SHALL skip complete units on a resume and refuse a resume or merge whose settings differ, naming the differing keys. A record SHALL hold seeds, statistics and counts only, and no per-feature array.

#### Scenario: Resume skips completed units
- **WHEN** a shard is restarted after some units completed
- **THEN** only the remaining units are evaluated, and the merged output has one record per unit

#### Scenario: Settings mismatch is refused
- **WHEN** a resume or merge is attempted with settings whose hash differs from the stored one
- **THEN** the driver exits listing the differing keys

### Requirement: Merge writes the rate, paired and cell tables

Merging a profile SHALL write the univariate baseline's tables for its designs: a rate table (per design, cell and arm: available and rejected counts, the rate and its Monte Carlo SE), a paired table (`production` against each per-feature arm), a cell table (rank quantiles, held-out AUROC, capture shares, hits and discoveries) and accounting (units, failures by cell, duplicates, recorded core-hours).

#### Scenario: Tables per design
- **WHEN** a profile is merged
- **THEN** every design, cell and arm has one rate row, and every design and cell has one cell row

#### Scenario: Failed units are excluded
- **WHEN** a cell has failed units
- **THEN** they are counted in the accounting and left out of every rate

### Requirement: The decision judges size per arm on the new designs

The decision SHALL apply the stage-span paper grade's count-plus-ceiling bound rule, imported unchanged, to the anchor and translation rates of `production`, `fdr` and `total` in two families: the four new `saturated` designs (m = 8) and the three `band` designs (m = 6). A per-feature arm that fails a family SHALL be reported as inadmissible and SHALL leave the map.

#### Scenario: Production fails size
- **WHEN** `production` fails either family
- **THEN** the verdict is `production_size_excess`, naming the family and every rate above its bound or ceiling

#### Scenario: A per-feature arm fails size
- **WHEN** `fdr` or `total` fails a family
- **THEN** the arm is recorded as inadmissible with every rate above its bound or ceiling, and no map label uses it

### Requirement: The decision labels each matched pair-cell by exact unpaired tests

For each signal cell of each matched pair, the decision SHALL compare `production`'s rejected counts in the two designs by Fisher's exact test, one-sided in each direction. The label SHALL be `weak_costs` when the weaker-staged design's rate is lower by at least 0.10 with a p-value at or below 0.025 / 53, `weak_gains` when it is higher by the same margin and level, and `unresolved` otherwise.

#### Scenario: The pairs
- **WHEN** the pair-cells are enumerated
- **THEN** they are the eleven signal cells of each of the four `saturated` pairs and the nine of the `band` design at δ_c / 2 and n = 1,280 against the one at δ_c and n = 320, 53 in all

#### Scenario: A cost
- **WHEN** the weaker-staged design rejects in 82 of 200 units and its match in 164 of 200
- **THEN** the pair-cell is labelled `weak_costs`

#### Scenario: A significant difference below the margin
- **WHEN** the exact p-value passes its level and the two rates differ by less than 0.10
- **THEN** the pair-cell is labelled `unresolved`

#### Scenario: The reference counts are the committed ones
- **WHEN** a `saturated` pair-cell is compared
- **THEN** the match's counts are read from the committed rate table

#### Scenario: Unavailable readings
- **WHEN** a pair-cell has no available unit in either design
- **THEN** it is labelled `unresolved` and named as unavailable

### Requirement: The per-feature false-discovery arm controls the pairing

The decision SHALL apply the same test and labels to the `fdr` arm in every pair-cell. A `production` label SHALL count toward the verdict unless `fdr` carries the same label in that pair-cell. Every label of both arms SHALL be recorded, counted or not.

#### Scenario: The control shares the cost
- **WHEN** `production` and `fdr` are both labelled `weak_costs` in a pair-cell
- **THEN** the label is recorded and does not count toward the verdict

#### Scenario: The control moves the other way
- **WHEN** `production` is labelled `weak_costs` and `fdr` is `unresolved` or `weak_gains` in a pair-cell
- **THEN** the label counts

### Requirement: The decision writes one predeclared verdict

The decision SHALL write `stage_effect_scaling_decision.json` (rule version 1) with one verdict, the first that applies in this order: `reproduction_failed`, `production_size_excess`, then `weak_stage_costs`, `weak_stage_gains`, `mixed` or `scaling_unrefuted` by which counted labels occur in at least one pair-cell. Every criterion SHALL be recorded whatever the verdict, with each verdict's consequence statement.

#### Scenario: Costs only
- **WHEN** at least one counted label is `weak_costs` and none is `weak_gains`
- **THEN** the verdict is `weak_stage_costs`, naming the pair-cells

#### Scenario: Both
- **WHEN** counted labels of both kinds occur
- **THEN** the verdict is `mixed`, naming the pair-cells of each kind

#### Scenario: Neither
- **WHEN** no counted label occurs
- **THEN** the verdict is `scaling_unrefuted`, and its consequence states that a difference below the margin is not excluded

#### Scenario: Gate
- **WHEN** the reproduction fails
- **THEN** the verdict is `reproduction_failed`, and size, the labels, the grades and the map are still recorded

### Requirement: The decision grades per design whether production names the kind of difference

For each design, the four committed ones included, the decision SHALL take the in-span cells under the planted scale where `production` rejects in at least 0.80 of units and judge the two off-target statistics of each, as one family per design, under the bound rule of the size criterion. The grade SHALL be `names_kind` when the family holds, `detects_only` when it fails and `no_detected_cell` when no cell qualifies. The grades SHALL be reported beside the verdict and SHALL NOT change it.

#### Scenario: Committed designs are graded the same way
- **WHEN** a committed design is graded
- **THEN** its family is built from the committed rate table by the same rule as a new design's

#### Scenario: A design that detects without naming
- **WHEN** a design's off-target family fails the bound rule
- **THEN** its grade is `detects_only`, naming every off-target rate above its bound or ceiling

#### Scenario: Grades are outside the verdict
- **WHEN** the decision is computed with the off-target rates changed
- **THEN** the verdict is unchanged

### Requirement: The decision reports the map against the per-feature arms at each new design

Over the 71 signal cells of the new designs, the decision SHALL label each cell `production_leads`, `baseline_leads` or `unresolved` by the univariate baseline's paired exact tests against the admissible per-feature arms, with its levels and its 0.10 margin, and SHALL report the count of cells under each label by plant family and design. The map SHALL NOT change the verdict.

#### Scenario: The map
- **WHEN** the decision is written
- **THEN** it has one row per plant family and new design with the number of cells under each label

#### Scenario: No admissible arm
- **WHEN** neither `fdr` nor `total` is admissible
- **THEN** every map cell is `unresolved` and the reason is recorded

### Requirement: The decision reports descriptive tables that never change the verdict

Beside the verdict the decision SHALL report: per design, the committed ones included, δ, n, δ√n, the held-out AUROC, the CV rank quantiles and the capture shares; every arm's rate per cell; the rates of `production`, each statistic and `fdr` on the stretch cells against e·δ·√n; for each pair-cell labelled `weak_costs`, the weaker-staged design's rates by selected rank; and the calibration's ladder.

#### Scenario: A design can be placed
- **WHEN** the decision is written
- **THEN** every design has one index row with its held-out AUROC, rank and δ√n

#### Scenario: Descriptive tables are outside the verdict
- **WHEN** the decision is computed with a recorded arm's rates changed
- **THEN** the verdict, the grades and the map are unchanged

### Requirement: The run is recorded in a dated report

The run SHALL be recorded in a dated report under `docs/reports/`, with aggregates only under `results/stage-effect-scaling-<date>/`: the merged tables, the calibration, reproduction and accounting, the decision, and provenance (revision, tree hash, job ids, core-hours, package versions). Per-unit records SHALL stay on the cluster.

#### Scenario: Report contents
- **WHEN** the report is written
- **THEN** it states the verdict, each arm's size reading, the reproduction, the labels by pair and plant family with the control's, the grade by design, the map, the descriptive tables and what follows from the verdict
- **AND** every number it quotes can be read from the committed aggregates

#### Scenario: The limits and the probe are stated
- **WHEN** the report states the verdict
- **THEN** it states in the same section that the reading holds at 658 simulated features with Gaussian noise, equal groups and three stages, and it states that the decision was written after a timing probe of null datasets
