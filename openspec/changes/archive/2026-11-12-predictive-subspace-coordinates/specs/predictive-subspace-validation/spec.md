# Spec Delta

## Purpose

Validates measuring the trajectory statistics in the predictive subspace of the production stage-PLS fit against the production score measurement, on matched units that are fresh or reproduce committed runs, under a predeclared decision on size, power, specificity and gain that gates any paper grade of a change to the production measurement.

## ADDED Requirements

### Requirement: Each unit measures one production fit under six arms at one seed

Every unit SHALL generate or assemble its dataset once and fit the production stage-PLS space once, with net displacement, rank floor 2 and, in simulation, span capture. It SHALL evaluate these arms at the unit's RRPP seed with identical permutation draws:
- `pls_cv`: the production score measurement;
- `predictive`, `predictive_orthonormal` and `predictive_fitted`: the three predictive-subspace coordinate systems of that same fit, built without refitting;
- `forced_k_minus_1`: the fit forced to the number of stages minus one components, in scores;
- `predictive_rank_plus_one`: the `predictive` coordinates of the fit forced to one component more than the production rank, in the simulation profile only.

On the plasmode profile each arm SHALL be one real-data analysis. Each predictive arm SHALL run at a forced rank equal to the unit's `pls_cv` rank, so that it fits the production model without repeating its cross-validation, and its record SHALL carry `pls_cv`'s held-out AUROC.

A failure SHALL be recorded for the unit without aborting the shard.

#### Scenario: One cross-validation per unit
- **WHEN** a simulation unit runs
- **THEN** exactly one cross-validated PLS fit and two forced-rank fits are made, and `pls_cv` and the three predictive arms report the same selected rank and held-out AUROC

#### Scenario: One cross-validation per plasmode unit
- **WHEN** a plasmode unit runs
- **THEN** exactly one cross-validation is made, and the three predictive arms' statistics and p-values equal those of a cross-validated analysis under the same setting

#### Scenario: Record contents
- **WHEN** a unit completes
- **THEN** its record names the profile, design, cell, mode, effect, replicate and both seeds
- **AND** per arm it holds the observed `delta`, `angle` and `shape`, their p-values, the measured dimension, the `angle` null 95th percentile and, in simulation, the span-capture shares and the `population_latent` statistics

#### Scenario: Records hold no donor data
- **WHEN** a plasmode unit is recorded
- **THEN** the record holds seeds and statistics only, with no donor identifier and no label vector

#### Scenario: Failed unit
- **WHEN** generation or fitting fails for a unit
- **THEN** the unit is recorded as failed with the error, and the shard continues

### Requirement: Each unit checks that the predictive arms measure the production fit

Each unit SHALL run these checks and record a failing unit as inconsistent, naming the check and its deviation:
- **Rank:** the response loadings of the production fit have rank equal to the smaller of the fit rank and the number of stages minus one.
- **Fitted labels:** the fitted stage indicators recomputed from the predictive subspace equal the production fit's to 1e-8.
- **Same fit:** every predictive arm's fitted model is the production fit: the same rank, and rotations and response loadings equal to 1e-8.
- **One subspace:** the three predictive arms report equal span-capture shares, in simulation.
- **Identity at a low rank:** when the production rank is at most the number of stages minus one, the `predictive` arm's p-values equal `pls_cv`'s.

#### Scenario: Consistent unit
- **WHEN** a unit's arms are compared
- **THEN** each check's deviation is recorded and within its tolerance

#### Scenario: Inconsistent unit
- **WHEN** a check fails
- **THEN** the unit is recorded as inconsistent with the check's name and deviation, and it is excluded from every rate

### Requirement: The simulation profile is a fresh seed family across stage effects

The `simulation` profile SHALL build its units from the latent-space comparison's design generators and seed derivations, imported unchanged, at base seed 2036, under the planted scale, at 100 replicates and 199 permutations. Its designs are:
- n = 200 and the SEA-AD table at δ = 2.0;
- n = 200 at δ = 1.0, 0.5 and 0.25, and the SEA-AD table at δ = 0.5;
- the SEA-AD table at the comparison's committed weak δ, whose calibration sha256 the settings record.

Each design has eleven stage-span cells: the anchor, translation at e = 1.0, magnitude at e ∈ {0.05, 0.1, 0.25}, orientation at e ∈ {0.1, 0.25, 0.5} with alignment 0, and shape at e ∈ {0.1, 0.25, 0.5}. One further cell applies the default feature-set orientation surgery at e = 1.0 on the SEA-AD table at δ = 2.0.

#### Scenario: Units and grade
- **WHEN** the simulation profile enumerates its units
- **THEN** it has seven designs of eleven cells and one surgery cell at replicates 0–99, 7,800 units in all, at 199 permutations

#### Scenario: A fresh family
- **WHEN** a unit's dataset seed is derived
- **THEN** it differs from the latent-space comparison's seed for the same design name, cell and replicate, and no unit shares a dataset seed with the comparison's or the production selection's runs

#### Scenario: Matched across cells
- **WHEN** two cells of one design and replicate are generated
- **THEN** they share their null dataset, and every arm of a unit shares its RRPP seed

### Requirement: The weak-size profile reuses the comparison's weak controls and must reproduce the span-split check

The `weak_size` profile SHALL rebuild the latent-space comparison's `seaad_weak` anchor and translation units at replicates 0–499 and 999 permutations, with the comparison's base seed, seed derivations and committed weak δ. For each cell and statistic, the `pls_cv` arm's count and rejections SHALL equal the `pls_cv` rows of the span-split check's committed `production_table.csv` exactly.

#### Scenario: Weak units
- **WHEN** the weak-size profile enumerates its units
- **THEN** it has 1,000 units, and each unit's dataset and RRPP seed equal the span-split check's simulation unit for the same design, cell and replicate

#### Scenario: Reproduction holds
- **WHEN** every row matches the committed counts
- **THEN** the merge records the reproduction as passed with the number of rows compared

#### Scenario: Reproduction fails
- **WHEN** any row differs
- **THEN** the merge records the reproduction as failed and lists every differing row with both counts

### Requirement: The plasmode profile reads the export through its contract and must reproduce two committed runs

The `plasmode` profile SHALL read the SEA-AD astrocyte export through the plasmode's input contract, refusing a file whose sha256 differs from the committed plasmode settings'. It SHALL permute the sex labels within stage with the plasmode's `paper` stream and analyse each unit at the case-study configuration at 999 permutations. Its cells are:
- the anchor at replicates 0–499;
- the production selection's real-noise span plant, imported unchanged, at e ∈ {0.25, 0.5, 1.0} and replicates 0–99.

The `pls_cv` arm SHALL reproduce, exactly:
- the plasmode's committed anchor rejection counts over 500 replicates for `delta`, `angle` and `shape`;
- the production selection's committed `pls_cv` counts and rejections for the three plant cells and every statistic;
- the label-free fidelity values (selected rank and held-out AUROC) on every anchor unit.

#### Scenario: Plasmode units
- **WHEN** the plasmode profile enumerates its units
- **THEN** it has 500 anchor units and 300 plant units, 800 in all

#### Scenario: A changed export file is refused
- **WHEN** an export file's sha256 differs from the committed one
- **THEN** the driver exits naming the file and runs no unit

#### Scenario: Reproduction and fidelity
- **WHEN** the plasmode profile is merged
- **THEN** the merge records, per check, the committed and reproduced values and whether every one matches

### Requirement: Runs are resumable, sharded and bound to their settings

Each profile SHALL write a settings file with a hash, the sha256 of every parent file it binds to, and one record per unit to its shard. It SHALL skip complete units on a resume, refuse a resume or merge whose settings differ, naming the differing keys, and refuse a changed parent file by name.

#### Scenario: Resume skips completed units
- **WHEN** a shard is restarted after some units completed
- **THEN** only the remaining units are evaluated, and the merged output has one record per unit

#### Scenario: Settings mismatch is refused
- **WHEN** a resume or merge is attempted with settings whose hash differs from the stored one
- **THEN** the driver exits listing the differing keys

#### Scenario: A changed parent is refused
- **WHEN** a committed parent file's sha256 differs from the one in the stored settings
- **THEN** the driver exits naming the file and runs no unit

### Requirement: Merge writes the rate, paired and cell tables

Merging a profile SHALL write:
- **a rate table:** one row per design, cell, statistic and arm with the available and rejected counts, the rate and its Monte Carlo SE;
- **a paired table:** per design, cell, statistic and candidate or control arm against `pls_cv`, both rates, the discordant counts over units available in both arms, and the exact one-sided loss, one-sided gain and two-sided paired p-values;
- **a cell table:** per design, cell and arm, the selected rank (median, 10th and 90th percentiles), the held-out AUROC, the measured dimension, the span-capture shares, the median observed and `population_latent` statistics, and the median `angle` null 95th percentile;
- **accounting:** units, failures by cell, inconsistent units, duplicates and recorded core-hours.

#### Scenario: Paired table
- **WHEN** a profile is merged
- **THEN** each design, cell, statistic and non-production arm has one paired row whose discordant counts are taken over units available in both arms

#### Scenario: No discordant units
- **WHEN** a paired row has no unit where exactly one arm rejects
- **THEN** its three paired p-values are 1

### Requirement: The decision judges size of each candidate arm (S1)

For each of the three predictive arms the decision SHALL apply the stage-span paper grade's count-plus-ceiling bound rule, imported unchanged, to the arm's control rates of `delta`, `angle` and `shape`, as three families. All three SHALL hold.
- **Simulation family:** the anchor and translation of the seven designs (m = 42).
- **Weak-size family:** the two `weak_size` cells (m = 6).
- **Plasmode family:** the plasmode anchor (m = 3).

#### Scenario: Size holds
- **WHEN** the three families pass the bound rule for an arm
- **THEN** S1 is recorded as holding for that arm, with every observation's rate, bound and status

#### Scenario: Size fails
- **WHEN** any family fails the bound rule for an arm
- **THEN** S1 fails for that arm, naming the family and every rate above its bound or ceiling

### Requirement: The decision judges protected power of each candidate arm (S2)

For each candidate arm the decision SHALL compute a one-sided exact paired loss p-value against `pls_cv` for each of 57 target cells: magnitude `delta`, orientation `angle` and shape `shape` at their three effects in the six simulation designs at δ ≥ 0.25, and `delta` at the plasmode's three plant cells. S2 SHALL hold when no more cells have a loss p-value below 0.025 than the 95th percentile of Binomial(57, 0.025), and none is below 0.05/57.

#### Scenario: Power loss
- **WHEN** more target cells than allowed have a loss p-value below 0.025, or any is below 0.05/57
- **THEN** S2 fails for that arm, naming each such cell with both rates and the discordant counts

### Requirement: The decision judges specificity of each candidate arm (S3)

For each candidate arm the decision SHALL judge the off-target rates of the stage-span target cells (for orientation, `delta` and `shape`; for shape, `delta` and `angle`; for magnitude, `angle` and `shape`) in two families. Both SHALL hold.
- **δ = 2.0 family (m = 36):** the two δ = 2.0 designs, under the bound rule of S1.
- **Band family (m = 72):** the four designs at δ ∈ {1.0, 0.5, 0.25}, by a one-sided exact paired excess p-value against `pls_cv` on the off-target rejections. It holds when no more readings have an excess p-value below 0.025 than the 95th percentile of Binomial(72, 0.025), and none is below 0.05/72.

The decision SHALL count, per arm, the off-target readings of both families that are above the S1 bound.

#### Scenario: Specificity fails at δ = 2.0
- **WHEN** the δ = 2.0 family fails the bound rule for an arm
- **THEN** S3 fails for that arm, naming every rate above its bound or ceiling

#### Scenario: Specificity fails in the band
- **WHEN** more band readings than allowed reject off target in excess of production, or any does below 0.05/72
- **THEN** S3 fails for that arm, naming each such reading with both rates and the discordant counts

#### Scenario: Production's own band rates are reported
- **WHEN** the decision is written
- **THEN** it reports `pls_cv`'s off-target rate for every band reading beside each candidate's

### Requirement: The decision judges gain of each candidate arm (G)

For each candidate arm the decision SHALL compute a one-sided exact paired gain p-value against `pls_cv` for the 18 target cells of the two δ = 2.0 designs. G SHALL hold when at least one is at or below 0.05/18. The same p-value SHALL be recorded, without bearing on G, for the 36 target cells of the band, the 9 of the weak design and the plasmode's three plant cells.

#### Scenario: Gain holds
- **WHEN** an arm's gain p-value is at or below 0.05/18 in at least one δ = 2.0 target cell
- **THEN** G holds for that arm, naming the cells that pass

#### Scenario: Gain only outside δ = 2.0
- **WHEN** an arm gains in the band or on the plasmode but in no δ = 2.0 target cell
- **THEN** G fails for that arm, and the recorded gains are reported

### Requirement: The decision writes one predeclared verdict

The decision SHALL write `predictive_subspace_decision.json` (rule version 1). If any reproduction or fidelity check fails or any unit is inconsistent, the verdict is `reproduction_failed`. Otherwise each candidate arm receives the first failing criterion's verdict in the order S1 `size_excess`, S2 `power_loss`, S3 `specificity_excess`, G `no_gain`, else `qualifies`. The overall verdict is `propose_paper_grade`, naming one qualifying arm, or `keep_production` when none qualifies. Among qualifying arms the named one SHALL be the arm with the fewest off-target readings above the bound; ties are broken in the order `predictive`, `predictive_orthonormal`, `predictive_fitted`. Every criterion SHALL be evaluated and recorded for every arm whatever the verdict.

#### Scenario: A paper grade is proposed
- **WHEN** every reproduction passes, every unit is consistent and at least one arm qualifies
- **THEN** the verdict is `propose_paper_grade` and names one arm with the reason it was preferred

#### Scenario: Production is kept
- **WHEN** every reproduction passes and no arm qualifies
- **THEN** the verdict is `keep_production`, with each arm's verdict and failing readings

#### Scenario: Gate
- **WHEN** a reproduction or fidelity check fails or a unit is inconsistent
- **THEN** the verdict is `reproduction_failed`, and S1–G are still recorded for every arm

#### Scenario: Unavailable readings fail
- **WHEN** a reading a criterion needs has no available units
- **THEN** that criterion fails for the arm, naming the reading

#### Scenario: Preference among qualifying arms
- **WHEN** two arms qualify with different counts of off-target readings above the bound
- **THEN** the arm with the smaller count is named

### Requirement: The decision reports descriptive comparisons beside the verdict

Beside the verdict the decision SHALL report, without bearing on it:
- the `forced_k_minus_1` arm's rates, paired differences against `pls_cv` and against `predictive`, and `pooled_stage` capture per cell;
- the `predictive_rank_plus_one` arm's target and off-target rates per cell beside `predictive`'s;
- per cell and arm, the median observed `angle` and path-length difference beside the `population_latent` values;
- the surgery cell's rates per arm and statistic, with the `group_difference` capture;
- the weak design's target-cell rates per arm;
- the selected rank and held-out AUROC per design.

#### Scenario: The control arm
- **WHEN** the decision is written
- **THEN** for every target cell it reports `forced_k_minus_1`'s rate and capture beside `pls_cv`'s and `predictive`'s

#### Scenario: Rank sensitivity
- **WHEN** the decision is written
- **THEN** for every orientation cell at δ ≥ 0.25 it reports the off-target `shape` rate of `predictive` and of `predictive_rank_plus_one`, with the paired two-sided p-value

#### Scenario: Estimation
- **WHEN** the decision is written
- **THEN** for every orientation cell it reports, per arm, the median observed angle and the median population angle through the same map

### Requirement: The run is recorded in a dated report

The run SHALL be recorded in a dated report under `docs/reports/`, with aggregates only under `results/predictive-subspace-<date>/`: the merged tables, reproduction and accounting, the decision, and provenance (revision, tree hash, job ids, core-hours, package versions). Per-unit records SHALL stay on the cluster.

#### Scenario: Report contents
- **WHEN** the report is written
- **THEN** it states the verdict, every criterion's reading per arm, every reproduction, the descriptive comparisons and what follows from the verdict
- **AND** every number it quotes can be read from the committed aggregates
