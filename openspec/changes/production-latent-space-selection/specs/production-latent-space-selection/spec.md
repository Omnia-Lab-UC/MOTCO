# Spec Delta

## Purpose

Validates the production stage-PLS latent space against its block-weighted and sparse variants on matched units, in simulation and on real noise, under a predeclared decision on size, protected power, specificity and gain that gates any paper grade and any change to the production default.

## ADDED Requirements

### Requirement: Each unit evaluates four stage-supervised arms on one dataset and one set of draws

Every unit SHALL evaluate one dataset under four arms, each a stage-PLS latent space fitted with net displacement, rank floor 2 and span capture, at the unit's RRPP seed with serial draws: production; block weighting `equal_variance`; sparse stage-PLS with the block scope and the default ladder; and both together. Each arm SHALL be its own cross-validated fit. A failing arm SHALL be recorded as failed for that unit without aborting the unit or the shard.

#### Scenario: Four arms, one dataset
- **WHEN** a unit runs
- **THEN** the dataset is generated or built once, every arm's record names the same generator and RRPP seeds, and the arms differ only in the weighting and sparsity settings

#### Scenario: Record contents
- **WHEN** a unit completes
- **THEN** per arm its record holds the observed `delta`, `angle` and `shape`, their p-values, the selected rank, the selected keep fraction where sparsity is set, the held-out stage AUROC, the span-capture shares and the population target where the dataset has population trajectories, and the per-block support sizes where sparsity is set

### Requirement: The simulation profile reuses the latent-space comparison's designs and in-span cells

The `simulation` profile SHALL build its units from the latent-space comparison's three designs, base seed, seed derivations and committed weak δ, at 100 replicates and 199 permutations. Its cells SHALL be the comparison's seven in-span cells by name — the anchor, translation at e = 1.0, magnitude at e = 0.05, orientation at e ∈ {0.25, 0.5} with φ = 0, shape at e ∈ {0.25, 0.5} — and three further cells at every design: magnitude at e = 0.25, orientation at e = 1.0 with φ = 0, and shape at e = 1.0. The settings SHALL record the sha256 of the comparison's settings, calibration and rate table, and a changed file SHALL be refused by name.

#### Scenario: Units
- **WHEN** the simulation profile enumerates its units
- **THEN** it has 10 cells at each of 3 designs at replicates 0–99, 3,000 units in all

#### Scenario: Shared cells are the comparison's units
- **WHEN** a unit of one of the seven shared cells is built
- **THEN** its generator parameters, dataset seed and RRPP seed equal the latent-space comparison's for the same design, cell and replicate

#### Scenario: Changed parent is refused
- **WHEN** a committed parent file's sha256 differs from the stored one
- **THEN** the driver exits naming the file and runs no unit

### Requirement: The plasmode profile analyses the SEA-AD export with a real-noise span plant

The `plasmode` profile SHALL read the SEA-AD export through the plasmode's input contract, with file hashes equal to the committed plasmode settings', permute the sex labels within stage as the plasmode does for replicates 0–99 of its `paper` stream, and analyse each unit through the real-data entry point at the case-study configuration (the six covariates, rank floor 2, net displacement, the `df` correction, attribution off) at 999 permutations. Its cells SHALL be the anchor and a span plant at e ∈ {0.25, 0.5, 1.0}: every donor labelled `Male` at stage *s* is shifted by e times the difference between the pooled stage-*s* mean and the centroid of the pooled stage means, in the analysis's standardized units mapped back to the block's units, so that the planted group's centred stage configuration gains e times the pooled centred one (the pooled one stretched by 1 + e about its centroid, up to the permuted labels' sampling noise). The real labels SHALL never be analysed and no record SHALL hold a donor identifier or label vector.

#### Scenario: Plant geometry
- **WHEN** the span plant at e is applied to a block
- **THEN** in standardized units the planted group's stage means minus their centroid equal its unplanted ones plus e times the pooled stage means minus their centroid (to 1e-8), and rows not labelled `Male` are unchanged

#### Scenario: Anchor is unplanted
- **WHEN** the anchor cell is built
- **THEN** the blocks equal the export's blocks exactly

#### Scenario: Units
- **WHEN** the plasmode profile enumerates its units
- **THEN** it has 4 cells at replicates 0–99, 400 units in all

### Requirement: The production arm must reproduce both parents

For every shared simulation cell, design and statistic, the production arm's count, rejections and failures SHALL equal the latent-space comparison's committed rate table exactly. For every plasmode anchor unit, the production arm's selected rank and held-out AUROC SHALL equal the plasmode's committed label-free fidelity values exactly.

#### Scenario: Simulation reproduction
- **WHEN** the simulation profile is merged
- **THEN** the reproduction lists the 63 compared rows and whether every one matches, naming each differing row with both counts

#### Scenario: Plasmode fidelity
- **WHEN** a plasmode anchor unit completes
- **THEN** its production arm's rank is 5 and its AUROC equals 0.635556, or the unit is recorded as failing the fidelity gate

### Requirement: Runs are resumable, sharded and bound to their settings

Each profile SHALL write a settings file with a hash and one record per unit and arm to its shard, skip complete units on resume, and refuse a resume or merge whose settings differ, naming the differing keys.

#### Scenario: Resume
- **WHEN** a shard restarts after some units completed
- **THEN** only the remaining units run and the merged output has one record per unit and arm

### Requirement: Merge writes the paired comparison tables

Merging SHALL write, per profile: a rate table (one row per design or cell, statistic and arm: count, rejections, failures, rate, Monte Carlo SE); a paired table (per cell, statistic and candidate arm against production: both rates, the difference, the discordant counts over units available in both arms, exact one-sided loss and gain p-values and the two-sided p-value); a cell table (per cell and arm the median and q10/q90 of the selected rank, keep fraction, held-out AUROC and `pooled_stage` capture, the median population target where available, and on the plasmode the median support per block); accounting (units, failures by arm and cell, fidelity failures, duplicates, core-hours, per-arm timing); and the reproduction.

#### Scenario: Paired table
- **WHEN** a profile is merged
- **THEN** every cell and statistic has one paired row per candidate arm whose discordant counts are taken over units available in both arms

### Requirement: The decision judges size of every candidate arm (S1)

For each candidate arm, the anchor and translation rates of `delta`, `angle` and `shape` at the three simulation designs (m = 18) and the plasmode anchor rates (m = 3) SHALL each be judged as one family under the stage-span paper grade's count-plus-ceiling rule on the bound α + 2·√(α(1 − α)/R), R the available replicates. S1 holds when both families hold. An unavailable reading fails its family.

#### Scenario: Count and ceiling
- **WHEN** a family of m rates is judged
- **THEN** it holds when at most the 95th percentile of Binomial(m, p₂) of its rates exceed the bound, p₂ the exact Binomial(R, α) exceedance probability at the median R, and none exceeds the Bonferroni ceiling at 0.05/m

#### Scenario: Size excess
- **WHEN** a family of an arm has more rates above the bound than allowed, or any rate above its ceiling
- **THEN** S1 fails for that arm, naming every such rate with its bound and ceiling

#### Scenario: One chance exceedance
- **WHEN** one simulation control rate of an arm is above the bound, below its ceiling, and the family allows one exceedance
- **THEN** that family holds, and the rate is still listed as exceeding

### Requirement: The decision judges protected power (S2) and specificity (S3)

For each candidate arm at the two δ = 2.0 designs:
- **S2:** for each of the 16 target cells (magnitude `delta` at e ∈ {0.05, 0.25}, orientation `angle` at e ∈ {0.25, 0.5, 1.0}, shape `shape` at e ∈ {0.25, 0.5, 1.0}, at both designs) the one-sided exact paired loss p-value against production SHALL be computed; S2 holds when no more cells than the 95th percentile of Binomial(16, 0.025) have p < 0.025 and none has p < 0.05/16;
- **S3:** the off-target rates of those cells at e ≤ 0.5 (m = 24) SHALL be judged as one family under the S1 rule. The e = 1.0 cells' off-target rates SHALL be reported and SHALL NOT enter S3: production's own construction leak there exceeds the bound (the stage-span pilot, 2026-11-04).

#### Scenario: No discordant units
- **WHEN** a target cell has no unit where exactly one arm rejects
- **THEN** its loss p-value is 1

#### Scenario: Power loss
- **WHEN** more cells than allowed have a loss p-value below 0.025, or any is below 0.05/16
- **THEN** S2 fails, naming each such cell with both rates and the discordant counts

### Requirement: The decision judges gain where production has no power (G)

For each candidate arm, the one-sided exact paired gain p-value against production SHALL be computed for the weak design's 8 target cells and the plasmode's 3 span-plant `delta` cells. G holds when at least one of the 11 cells has a gain p-value at most 0.05/11. The pooled gain over the 11 cells SHALL be reported beside it.

#### Scenario: Gain
- **WHEN** one gain cell's exact paired p-value is at most 0.05/11
- **THEN** G holds for that arm, naming the cell, both rates and the discordant counts

#### Scenario: No gain
- **WHEN** no gain cell reaches 0.05/11
- **THEN** G fails, and every gain cell's rates and p-value are still reported

### Requirement: The decision writes one predeclared verdict

The decision SHALL write `production_latent_space_decision.json` (rule version 1). If either reproduction fails, the verdict is `reproduction_failed`. Otherwise each candidate arm receives the first failing criterion's verdict in the order S1 `size_excess`, S2 `power_loss`, S3 `specificity_excess`, G `no_gain`, else `qualifies`; the overall verdict is `propose_paper_grade` naming the first qualifying arm in the order `weighted`, `sparse_block`, `weighted_sparse_block`, else `keep_production`. Every criterion SHALL be evaluated and recorded for every arm whatever the verdict, with the descriptive tables (paired differences for every cell, rank and fraction distributions, AUROC and capture per arm, the `delta` metric note) beside it.

#### Scenario: Proposal
- **WHEN** both reproductions pass and `weighted` and `weighted_sparse_block` qualify
- **THEN** the verdict is `propose_paper_grade` naming `weighted`, and both arms' criteria are recorded

#### Scenario: Gate
- **WHEN** a reproduction fails
- **THEN** the verdict is `reproduction_failed`, and S1–G are still recorded for every arm

### Requirement: The run is recorded in a dated report

The run SHALL be recorded in a dated report under `docs/reports/`, with aggregates only under `results/production-latent-space-<date>/`: the merged tables, reproduction and accounting, the decision, and provenance (revision, tree hash, job ids, core-hours, package version). Per-unit records SHALL stay on the cluster.

#### Scenario: Report contents
- **WHEN** the report is written
- **THEN** it states the verdict, every arm's criteria, both reproductions, the gain cells with their paired readings, and what follows from the verdict
