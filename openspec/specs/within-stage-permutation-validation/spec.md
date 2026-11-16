# within-stage-permutation-validation Specification

## Purpose
Validates the within-stage permutation null against RRPP's unrestricted null, on matched units that are fresh or reproduce committed runs, under a predeclared decision on size, protected power and specificity for the production measurement, and re-reads the `predictive_orthonormal` coordinates against production under the same draws.

## Requirements

### Requirement: Each fitted unit measures one production fit in four coordinate systems under both nulls

Every unit of the `simulation`, `weak_size` and `plasmode` profiles SHALL generate or assemble its dataset once and fit the production stage-PLS space once, with net displacement and rank floor 2. It SHALL evaluate eight arms at the unit's RRPP seed:
- `pls_cv`, `predictive`, `predictive_orthonormal` and `predictive_fitted`: the production scores and the three predictive-subspace coordinate systems of that same fit, under the unrestricted draws;
- the same four with the suffix `__within_stage`: the same matrices under draws restricted within stage, requested through the package's permutation strata setting.

In simulation each arm SHALL be one harness evaluation on the prebuilt fit. On the plasmode each arm SHALL be one real-data analysis at the case-study configuration (six covariates, the `df` residual correction), forced at the unit's rank. A failure SHALL be recorded for the unit without aborting the shard.

#### Scenario: One fit per unit
- **WHEN** a simulation unit runs
- **THEN** exactly one cross-validated PLS fit is made, and all eight arms report its selected rank and held-out AUROC

#### Scenario: Unrestricted arms equal the pilot's
- **WHEN** a unit of the predictive-subspace pilot is evaluated by this driver
- **THEN** its four unrestricted arms' statistics and p-values equal the pilot's records for that unit

#### Scenario: Record contents
- **WHEN** a unit completes
- **THEN** its record holds the profile, design, cell, replicate and seeds, the rank and AUROC, and per arm the observed statistics, p-values, `angle` null 95th percentile and timing
- **AND** it holds no latent matrix, no label vector and no donor identifier

### Requirement: Each unit checks that the two nulls read the same measurement

Each unit SHALL run these checks and record a failing unit as inconsistent, naming the check and its deviation:
- **The pilot's checks:** the rank, fitted-label, same-fit, one-subspace and low-rank identity checks of the predictive-subspace pilot, on the four unrestricted arms.
- **Same statistics:** each within-stage arm's observed `delta`, `angle` and `shape` equal its unrestricted counterpart's exactly.
- **Construction identity** (plasmode): each within-stage arm's p-values equal those of the predictive size-excess diagnosis's within-stage loop on that arm's latent matrix, models, residual correction and seed, exactly.

#### Scenario: Consistent unit
- **WHEN** a unit's arms are compared
- **THEN** each applicable check's deviation is recorded and within its tolerance

#### Scenario: Inconsistent unit
- **WHEN** a check fails
- **THEN** the unit is recorded as inconsistent with the check's name and deviation, and it is excluded from every rate

### Requirement: The simulation profile is a fresh seed family on the pilot's designs

The `simulation` profile SHALL build the predictive-subspace pilot's seven stage-span designs with their eleven cells and its surgery cell, through the pilot's builders imported unchanged, at base seed 2038, 100 replicates and 199 permutations, with rejection at p < α. No dataset seed SHALL equal one of the pilot (base seed 2036), of the latent-space comparison (2035) or of the probes recorded in the change's notes.

#### Scenario: Unit count
- **WHEN** the profile's units are enumerated at the committed replicates
- **THEN** there are 7,800, and every cell of a design and replicate shares its null dataset

#### Scenario: Fresh datasets
- **WHEN** the profile's dataset seeds are compared with the pilot's and the comparison's
- **THEN** none coincides

### Requirement: The weak-size profile reuses the pilot's units and must reproduce its table

The `weak_size` profile SHALL evaluate the predictive-subspace pilot's 1,000 weak-design control units, built by the pilot's functions imported unchanged, at 999 permutations. The available and rejected counts of its four unrestricted arms SHALL equal the pilot's committed weak-size rate table in every cell and statistic.

#### Scenario: Reproduction passes
- **WHEN** the merged counts of the four unrestricted arms equal the committed table's
- **THEN** the reproduction is recorded as passed with the number of compared rows

#### Scenario: Reproduction fails
- **WHEN** any count differs
- **THEN** the reproduction is recorded as failed, listing every differing row with both values

### Requirement: The plasmode profile reads the anchor on reproduced and on fresh label draws

The `plasmode` profile SHALL read the SEA-AD astrocyte export in place through the plasmode's input contract and analyse only labels permuted within stage. It SHALL evaluate the anchor at replicates 0–999 of the plasmode's `paper` label stream at 999 permutations, with rejection at p ≤ α, each unit fitted at the plasmode's committed fidelity rank without cross-validation.
- **Pre-flight:** before any unit, the cross-validated production fit on replicates 0 and 500 SHALL select the fidelity rank at the committed held-out AUROC and give the forced fit's statistics and p-values exactly.
- **Reproduced draws (replicates 0–499, cell `anchor`):** the four unrestricted arms' counts SHALL equal the predictive-subspace pilot's committed plasmode anchor rows, and the four within-stage arms' counts SHALL equal the predictive size-excess diagnosis's committed within-stage anchor rows.
- **Fresh draws (replicates 500–999, cell `anchor_fresh`):** label draws that no committed run has read.

No record SHALL hold a donor identifier, a label vector or a matrix.

#### Scenario: Both parents reproduced
- **WHEN** the merged counts of cell `anchor` equal both committed tables in all 24 rows
- **THEN** the reproduction is recorded as passed

#### Scenario: Fresh replicates are unread
- **WHEN** the committed settings of every run that read the `paper` stream are inspected
- **THEN** none evaluated a replicate at or above 500

#### Scenario: A changed export is refused
- **WHEN** an export file's sha256 differs from the committed plasmode settings'
- **THEN** the driver exits naming the file and runs no unit

#### Scenario: Failed pre-flight
- **WHEN** the cross-validated fit selects another rank or differs from the forced fit on either pre-flight replicate
- **THEN** the pre-flight is recorded as failed and the decision reads `reproduction_failed`

### Requirement: The dispersion profile reads both nulls where stage dispersion is unequal

The `dispersion` profile SHALL draw Gaussian latent matrices with no group difference on the SEA-AD cell table (F 11 / 10 / 28, M 9 / 10 / 12), with stage means at fixed centroids and a noise variance that depends on the stage only, and read `delta`, `angle` and `shape` under net displacement by the RRPP primitive, without covariates or residual correction, under the unrestricted draws (`matrix`) and under draws restricted within stage (`matrix__within_stage`) at one seed. It SHALL cross five variance patterns with 2 and 5 dimensions, at 1,000 replicates and 999 permutations, with rejection at p < α and seeds derived from base seed 2038:
- `equal`: 1, 1, 1;
- `predictive`: 1.120, 1.278, 0.802, the diagnosis's median relative dispersion by stage in `predictive_orthonormal`;
- `scores`: 0.805, 0.921, 1.137, its medians in the scores;
- `predictive_doubled` and `scores_doubled`: the same patterns with every multiplier squared.

The centroids SHALL give a population stage share of variance of 0.84 at two dimensions under equal dispersion, the anchor's in-sample stage R².

#### Scenario: Unit count
- **WHEN** the profile's units are enumerated at the committed replicates
- **THEN** there are 10,000, and both arms of a unit read one matrix at one seed

#### Scenario: No fit
- **WHEN** a dispersion unit runs
- **THEN** no latent space is fitted, and its record holds seeds, statistics and p-values only

### Requirement: Runs are resumable, sharded and bound to their settings

Each profile SHALL write a settings file with a hash, the sha256 of every parent file it binds to, and one record per unit to its shard. It SHALL skip complete units on a resume, refuse a resume or merge whose settings differ, naming the differing keys, and refuse a changed parent file by name. The parents SHALL include the predictive-subspace pilot's settings and rate tables and the predictive size-excess diagnosis's settings and rate table.

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
- **a paired table:** per design, cell and statistic, each within-stage arm against its unrestricted counterpart, and each predictive arm against the production arm under the same draws: both rates, the discordant counts over units available in both arms, and the exact one-sided loss, one-sided gain and two-sided paired p-values, 1 when no unit is discordant;
- **a cell table:** per design, cell and arm, the selected rank, the held-out AUROC, the measured dimension, the span-capture shares where measured, and the median `angle` null 95th percentile;
- **accounting:** units, failures by cell, inconsistent units, duplicates, recorded core-hours and per-arm timing;
- **reproduction:** every comparison with a committed parent, with every differing row.

#### Scenario: Paired table
- **WHEN** a fitted profile is merged
- **THEN** each design, cell and statistic has four rows pairing the nulls and three rows per null pairing the predictive arms with production

#### Scenario: No discordant units
- **WHEN** a paired row has no unit where exactly one arm rejects
- **THEN** its three paired p-values are 1

### Requirement: The decision judges the within-stage null on the production scores

The decision SHALL judge the within-stage null by three criteria, each evaluated and recorded whatever the others read.
- **Size.** The stage-span paper grade's count-plus-ceiling bound rule, imported unchanged, on the `delta`, `angle` and `shape` control rates in four families, all of which SHALL hold: `pls_cv__within_stage` on the simulation's anchors and translations (m = 42), on the two `weak_size` cells (m = 6) and on the plasmode's `anchor_fresh` cell (m = 3), and `matrix__within_stage` on the ten dispersion cells (m = 30).
- **Protected power.** A one-sided exact paired loss p-value of `pls_cv__within_stage` against `pls_cv` in each of the 54 target cells of the six simulation designs at δ ≥ 0.25 (magnitude `delta`, orientation `angle` and shape `shape` at their three effects). It SHALL hold when no more cells have a loss p-value below 0.025 than the 95th percentile of Binomial(54, 0.025), and none is below 0.05/54.
- **Specificity.** The off-target rates of `pls_cv__within_stage` in the stage-span target cells, in two families, both of which SHALL hold. At the two δ = 2.0 designs (m = 36), the bound rule of the size criterion. In the four band designs (m = 72), a reading counts as an excess when its rate is above the size bound and its one-sided exact paired excess p-value against `pls_cv` is below 0.025; the family holds when no more readings count than the 95th percentile of Binomial(72, 0.025), and no reading above the bound has an excess p-value below 0.05/72.

The one-sided exact paired gain p-value against `pls_cv` SHALL be recorded for every target cell without bearing on the criteria.

#### Scenario: Size fails where dispersion is unequal
- **WHEN** the dispersion family fails the bound rule for `matrix__within_stage`
- **THEN** the size criterion fails, naming the family and every rate above its bound or ceiling

#### Scenario: A restored rate is not an excess
- **WHEN** a band off-target reading is rejected more often under within-stage draws than under unrestricted ones, and its rate is at or below the size bound
- **THEN** it does not count against specificity, and both rates are reported

#### Scenario: Power loss
- **WHEN** more target cells than allowed have a loss p-value below 0.025, or any is below 0.05/54
- **THEN** protected power fails, naming each such cell with both rates and the discordant counts

#### Scenario: Unavailable readings fail
- **WHEN** a reading a criterion needs has no available units
- **THEN** that criterion fails, naming the reading

### Requirement: The decision re-reads predictive_orthonormal under the same draws

The decision SHALL judge `predictive_orthonormal__within_stage` against `pls_cv__within_stage` by the predictive-subspace pilot's four criteria with that reference:
- **S1 size:** the bound rule on its control rates in the simulation (m = 42), `weak_size` (m = 6) and `anchor_fresh` (m = 3) families;
- **S2 protected power:** the paired loss rule on the 54 target cells;
- **S3 specificity:** the bound rule on the δ = 2.0 off-target readings (m = 36), and the pilot's paired excess rule on the band's (m = 72);
- **G gain:** at least one one-sided exact paired gain p-value at or below 0.05/18 among the 18 target cells of the two δ = 2.0 designs.

The arm SHALL receive the first failing criterion's reading in the order `size_excess`, `power_loss`, `specificity_excess`, `no_gain`, else `qualifies`. The reading SHALL be recorded whatever the null's criteria read, and SHALL bear on the verdict only when all three of the null's criteria hold.

#### Scenario: A gain that is the null's
- **WHEN** `predictive_orthonormal__within_stage` rejects a δ = 2.0 target cell more often than `pls_cv` but not more often than `pls_cv__within_stage`
- **THEN** that cell does not count toward G

#### Scenario: Candidate recorded under a failed null
- **WHEN** the null fails its size criterion
- **THEN** the candidate's four criteria are still recorded, and the verdict does not depend on them

### Requirement: The decision writes one predeclared verdict

The decision SHALL write `within_stage_permutation_decision.json` (rule version 1) with one verdict and the consequence statement fixed with the rule:
- `reproduction_failed`: a reproduction or the pre-flight failed, or a unit is inconsistent;
- `size_excess`, `power_loss` or `specificity_excess`: the null's first failing criterion, in that order;
- `propose_paper_grade`: the null's three criteria hold and the candidate qualifies;
- `within_stage_validated`: the null's three criteria hold and the candidate does not qualify, naming its reading.

#### Scenario: Gate
- **WHEN** a reproduction fails or a unit is inconsistent
- **THEN** the verdict is `reproduction_failed`, and every criterion is still recorded

#### Scenario: The null is validated and production keeps its coordinates
- **WHEN** every reproduction passes, the null's three criteria hold and the candidate reads `no_gain`
- **THEN** the verdict is `within_stage_validated`, naming `no_gain`

#### Scenario: A paper grade is proposed
- **WHEN** every reproduction passes, the null's three criteria hold and the candidate qualifies
- **THEN** the verdict is `propose_paper_grade`, naming `predictive_orthonormal` under within-stage draws

#### Scenario: The null fails
- **WHEN** the null's size criterion holds and its protected power fails
- **THEN** the verdict is `power_loss`, with the failing cells

### Requirement: The decision reports descriptive comparisons beside the verdict

Beside the verdict the decision SHALL report, without bearing on it:
- the `predictive__within_stage` and `predictive_fitted__within_stage` arms under the candidate's four criteria;
- the unrestricted arms' rates on the dispersion cells and on `anchor_fresh`, beside the within-stage arms';
- each arm's rates over the pooled 1,000 anchor draws;
- the gain p-values of `pls_cv__within_stage` against `pls_cv` per target cell, the weak design's target cells included;
- the surgery cell's rates per arm and statistic;
- the median `angle` null 95th percentile per arm and cell;
- the diagnosis's committed plant-cell rates under both nulls, quoted from its rate table.

#### Scenario: The mechanism's direction
- **WHEN** the decision is written
- **THEN** for every dispersion cell it reports both arms' rates per statistic, with the pattern and the dimension

#### Scenario: Fresh draws beside reproduced ones
- **WHEN** the decision is written
- **THEN** for every arm it reports the `anchor`, `anchor_fresh` and pooled counts per statistic

### Requirement: The run is recorded in a dated report

The run SHALL be recorded in a dated report under `docs/reports/`, with aggregates only under `results/within-stage-permutation-<date>/`: the merged tables, reproduction and accounting, the pre-flight, the decision, and provenance (revision, tree hash, job ids, core-hours, package versions). Per-unit records SHALL stay off the repository.

#### Scenario: Report contents
- **WHEN** the report is written
- **THEN** it states the verdict, every criterion's reading for the null and for the candidate, every reproduction, the descriptive comparisons, the limitations and what follows from the verdict
- **AND** every number it quotes can be read from the committed aggregates
