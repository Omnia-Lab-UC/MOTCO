# Spec Delta

## Purpose

Validates measuring the trajectory statistics in orthonormal coordinates of the fitted PLS span against the production score coordinates, on matched units that reproduce two committed runs, under a predeclared decision on size, specificity and power that gates any change to the production measurement.

## ADDED Requirements

### Requirement: The paper-grade profile reuses the stage-span paper grade's units

The `paper_grade` profile SHALL rebuild the stage-span paper grade's units from its committed profiles' matched seeds, under the planted scale, at its five columns. Its cells are the anchor, translation, orientation and shape at e ∈ {0.1, 0.25, 0.5}, and magnitude at e ∈ {0.02, 0.05, 0.1}, each at 500 replicates and 999 permutations. The settings SHALL record the profiles' sha256.

#### Scenario: Units and grade
- **WHEN** the paper-grade profile enumerates its units
- **THEN** each of the five columns has 13 cells at replicates 0–499, 65 cells and 32,500 units in all, at 999 permutations

#### Scenario: Units are the paper grade's datasets
- **WHEN** a unit is generated
- **THEN** its generator parameters, dataset seed and RRPP seed equal the stage-span paper grade's for the same column, mode, effect and replicate

#### Scenario: A changed profile is refused
- **WHEN** a committed profile's sha256 differs from the one in the stored settings
- **THEN** the driver exits naming the profile and runs no unit

### Requirement: The weak profile reuses the latent-space comparison's weak design

The `weak` profile SHALL rebuild the latent-space comparison's `seaad_weak` anchor and translation (e = 1.0) units at replicates 0–499 and 999 permutations. It SHALL use the comparison's designs, base seed 2035, seed derivations and committed calibrated δ, whose sha256 the settings record.

#### Scenario: Weak units
- **WHEN** the weak profile enumerates its units
- **THEN** it has 1,000 units, and each unit's dataset and RRPP seed equal the span-split check's simulation unit for the same design, cell and replicate

### Requirement: Each unit measures one fit in both coordinate systems at one seed

Every unit SHALL generate its dataset once and fit the production stage-PLS space once, with net displacement, rank floor 2 and span capture. It SHALL evaluate that fit twice at the unit's RRPP seed: in score coordinates (`scores`) and in orthonormal span coordinates (`orthonormal`). A failure SHALL be recorded for the unit without aborting the shard.

#### Scenario: One fit, two measurements
- **WHEN** a unit runs
- **THEN** one PLS model is fitted, and both arms report the same selected rank, held-out AUROC and span-capture shares

#### Scenario: Record contents
- **WHEN** a unit completes
- **THEN** its record names the profile, column or design, mode, effect, replicate and both seeds
- **AND** per arm it holds the observed `delta`, `angle` and `shape`, their p-values, the pooled eigengap, the `angle` null 95th percentile, and the `population_latent` path lengths of A and B

#### Scenario: Failed unit
- **WHEN** generation or fitting fails for a unit
- **THEN** the unit is recorded as failed with the error, and the shard continues

### Requirement: Each unit checks that the two arms measure one span

Each unit SHALL run two checks, each to 1e-8:
- **Orthonormal basis:** the basis of the orthonormal arm is orthonormal.
- **Same span:** the score matrix equals the orthonormal arm's latent matrix times one square matrix, the basis transposed times the rotations.

A unit failing a check SHALL be recorded as inconsistent, naming the check and the deviation.

#### Scenario: Consistent unit
- **WHEN** a unit's two latent matrices are compared
- **THEN** the largest deviation of each check is recorded and is at most 1e-8

### Requirement: The paper-grade profile must reproduce the paper grade

For every column, cell and statistic of the paper-grade profile, the `scores` arm's available and rejected counts at p < 0.05 SHALL equal the `n_available` and `n_rejected` of the committed stage-span paper-grade `design_point_operating.csv` exactly. That includes the cells where four-stage orientation units failed.

#### Scenario: Reproduction holds
- **WHEN** every row matches the committed counts
- **THEN** the merge records the reproduction as passed with the number of rows compared

#### Scenario: Reproduction fails
- **WHEN** any row differs
- **THEN** the merge records the reproduction as failed and lists every differing row with both counts

### Requirement: The weak profile must reproduce the span-split check

For the weak design's anchor and translation, the `scores` arm's count and rejections for each statistic SHALL equal the `pls_cv` rows of the span-split check's committed `production_table.csv` exactly.

#### Scenario: Weak reproduction
- **WHEN** the weak profile is merged
- **THEN** the merge records, per cell and statistic, the committed and reproduced counts and whether every one matches

### Requirement: Runs are resumable, sharded and bound to their settings

Each profile SHALL write a settings file with a hash and one record per unit to its shard. It SHALL skip complete units on a resume and refuse a resume or merge whose settings differ, naming the differing keys.

#### Scenario: Resume skips completed units
- **WHEN** a shard is restarted after some units completed
- **THEN** only the remaining units are evaluated, and the merged output has one record per unit

#### Scenario: Settings mismatch is refused
- **WHEN** a resume or merge is attempted with settings whose hash differs from the stored one
- **THEN** the driver exits listing the differing keys

### Requirement: Merge writes the paired rate tables

Merging SHALL write one row per profile, column or design, cell, statistic and arm, with the available and rejected counts, the rate and its Monte Carlo SE. It SHALL also write a paired table and the run's accounting.
- **Paired table:** both rates, their difference, the discordant counts and the exact two-sided paired p-value.
- **Accounting:** units, failures, inconsistent units, duplicates and recorded core-hours.

#### Scenario: Paired table
- **WHEN** a profile is merged
- **THEN** each cell and statistic has one paired row whose discordant counts are taken over units available in both arms

### Requirement: The decision judges size of the orthonormal arm (O1)

The decision SHALL apply the stage-span paper grade's count-plus-ceiling bound rule to the `orthonormal` arm's control rates, as two families. Both SHALL hold.
- **Paper-grade family:** the anchor and the three translation cells at each of the five columns, for `delta`, `angle` and `shape` (m = 60).
- **Weak family:** the weak design's anchor and translation, for the same statistics (m = 6).

#### Scenario: Size holds
- **WHEN** both families pass the bound rule
- **THEN** O1 is recorded as holding, with every observation's rate, bound and status

#### Scenario: Size fails
- **WHEN** either family fails the bound rule
- **THEN** O1 fails, naming the family and every rate above its bound

### Requirement: The decision judges specificity of the orthonormal arm (O2)

The decision SHALL apply the same bound rule to the `orthonormal` arm's off-target rates, m = 90:
- orientation: `delta` and `shape`;
- shape: `delta` and `angle`;
- magnitude: `angle` and `shape`.

Each mode is read at its three effects in each of the five columns.

#### Scenario: Specificity fails
- **WHEN** the off-target family fails the bound rule
- **THEN** O2 fails, naming every rate above its bound or ceiling

### Requirement: The decision judges target power against the scores (O3)

For each of the 45 target cells (orientation `angle`, shape `shape` and magnitude `delta` at three effects in five columns), the decision SHALL compute a one-sided exact paired p-value for a loss of `orthonormal` against `scores`. O3 SHALL hold when no more cells have a loss p-value below 0.025 than the 95th percentile of Binomial(45, 0.025), and none is below 0.05/45.

#### Scenario: No discordant units
- **WHEN** a target cell has no unit where exactly one arm rejects
- **THEN** its loss p-value is 1

#### Scenario: Power loss
- **WHEN** more target cells than allowed have a loss p-value below 0.025, or any is below 0.05/45
- **THEN** O3 fails, naming each such cell with both rates and the discordant counts

### Requirement: The decision writes one predeclared verdict

The decision SHALL write `orthonormal_coordinates_decision.json` (rule version 1). If either reproduction fails or any unit is inconsistent, the verdict is `reproduction_failed`. Otherwise it is the first failing criterion's verdict in the order O1 `size_excess`, O2 `specificity_excess`, O3 `power_loss`, else `adopt_orthonormal`. Every criterion SHALL be evaluated and recorded whatever the verdict.

#### Scenario: Adoption
- **WHEN** both reproductions pass, every unit is consistent, and O1, O2 and O3 hold
- **THEN** the verdict is `adopt_orthonormal`

#### Scenario: Gate
- **WHEN** a reproduction fails or a unit is inconsistent
- **THEN** the verdict is `reproduction_failed`, and O1–O3 are still recorded

#### Scenario: Unavailable readings fail
- **WHEN** a reading a criterion needs has no available units
- **THEN** that criterion fails, naming the reading

### Requirement: The decision reports descriptive comparisons beside the verdict

Beside the verdict the decision SHALL report, without bearing on it:
- the paired off-target and target differences per cell, naming the 45° rotation's `shape` at n = 200 and at four stages;
- per cell and arm, the median population-level B − A path-length difference as a share of A's path;
- the median pooled eigengap and `angle` null 95th percentile per cell and arm.

#### Scenario: Watched cells
- **WHEN** the decision is written
- **THEN** it reports both arms' `shape` rate for the 45° rotation at n = 200 and at four stages, with the paired difference and its exact two-sided p-value

#### Scenario: Population residual
- **WHEN** the decision is written
- **THEN** for the shape cells at e = 0.5 it reports each arm's median population path-length difference per column

### Requirement: The run is recorded in a dated report

The run SHALL be recorded in a dated report under `docs/reports/`, with aggregates only under `results/orthonormal-coordinates-<date>/`:
- the merged tables;
- reproduction and accounting;
- the decision;
- provenance: revision, tree hash, job ids, core-hours.

Per-unit records SHALL stay on the cluster.

#### Scenario: Report contents
- **WHEN** the report is written
- **THEN** it states the verdict, every criterion's reading, both reproductions, the watched cells and the population residuals, and names what follows from the verdict
