# Spec Delta

## Purpose

Judges the production trajectory test against the default analysis with no integration, a group × stage interaction test on each feature, on matched simulated units over the three generator blocks. A predeclared decision states which of the two detects a planted difference, and a separate grade states whether production names the kind of difference.

## ADDED Requirements

### Requirement: Each unit reads every arm on one dataset, one production fit and one seed

Every unit SHALL generate its dataset once and fit the production stage-PLS space once, with net displacement and rank floor 2. It SHALL read every arm on that dataset at the unit's RRPP seed, so that all arms see the same permutation draws. A failure SHALL be recorded for the unit without aborting the shard.

#### Scenario: One cross-validation per unit
- **WHEN** a unit runs
- **THEN** exactly one cross-validated PLS fit is made, and the `pls_cv` record equals a plain production evaluation of the same dataset at the same seed

#### Scenario: The baseline does not depend on the fit
- **WHEN** the per-feature interaction F is computed from features standardized without any latent fit
- **THEN** it equals the F the unit records, and so does the scan's global p-value

#### Scenario: Failed unit
- **WHEN** generation, fitting or an arm fails for a unit
- **THEN** the unit is recorded as failed with the error, and the shard continues

### Requirement: Production enters the comparison as one test

The `production` arm SHALL reject when the smallest of the unit's `delta`, `angle` and `shape` p-values is at or below 0.05 / 3. The three statistics SHALL also be recorded singly at 0.05, and `production_minp` SHALL be recorded as the smallest of the three calibrated on the unit's joint permutation draws. Only `production` is judged.

#### Scenario: Omnibus rejection
- **WHEN** a unit's smallest trajectory p-value is 0.012
- **THEN** `production` rejects

#### Scenario: A single statistic at 0.05 does not reject the omnibus
- **WHEN** a unit's trajectory p-values are 0.030, 0.400 and 0.700
- **THEN** `production` does not reject, and the first statistic is recorded as rejecting singly

### Requirement: The baseline is the per-feature interaction F combined three ways

Each baseline arm SHALL be one test at 0.05 on the per-feature group × stage interaction F of the standardized joint features, under the trajectory test's full and reduced models:
- `scan` rejects on the permutation p-value of the largest F;
- `fdr` rejects when Benjamini–Hochberg at 0.05 over the parametric F p-values of the non-degenerate features has a discovery;
- `total` rejects on the split interaction test's `total` part.

#### Scenario: The scan arm is the package's scan
- **WHEN** a unit's `scan` arm is compared with the interaction scan run on the same features, models and seed
- **THEN** the largest F and the global p-value are identical

#### Scenario: The false-discovery arm on a known case
- **WHEN** the per-feature F values and degrees of freedom of a small case are given
- **THEN** the arm's adjusted p-values equal hand-computed Benjamini–Hochberg values of the upper-tail F probabilities, and it rejects exactly when the smallest is at or below 0.05

#### Scenario: Degenerate features
- **WHEN** a feature is degenerate for the scan
- **THEN** it is left out of the `fdr` arm's family and counted in the record

### Requirement: Further arms are recorded and decide nothing

Each unit SHALL record, without bearing on the verdict, the split test's `span` and `complement` parts and the trajectory statistics measured on `concat` and on PCA at the unit's production rank, each singly and as the omnibus of the `production` arm.

#### Scenario: Matched PCA uses the production rank
- **WHEN** a unit is evaluated
- **THEN** the `pca_matched` arm's dimension equals the `pls_cv` arm's selected rank on that unit

#### Scenario: Recorded arms are outside the decision
- **WHEN** the decision is computed with the recorded arms' rates changed
- **THEN** the verdict and the grade are unchanged

### Requirement: The unit record holds seeds and statistics only

A completed unit's record SHALL name its profile, design, cell, mode, effect, replicate and both seeds. Per arm it SHALL hold the p-values and the rejection under the arm's rule. It SHALL also hold the selected rank, the held-out AUROC, the span-capture shares, the number of scan hits and of `fdr` discoveries, and the count of each on population-null features where that set is defined. It SHALL hold no per-feature array.

#### Scenario: Record contents
- **WHEN** a unit completes
- **THEN** its record carries every field above and no array with one entry per feature

#### Scenario: Each arm keeps its validated rule
- **WHEN** a unit's rejections are recorded
- **THEN** `scan` uses p ≤ 0.05, the split parts and the single trajectory statistics use p < 0.05, `fdr` uses an adjusted p ≤ 0.05, and `production` uses p ≤ 0.05 / 3

### Requirement: The parent profile reruns the scan validation's units and must reproduce it

The `parent` profile SHALL rebuild the scan validation's simulation units for the `n200` and SEA-AD-table designs, with that run's settings, seeds and grade. For each design, cell and test, the counts and rejections of the scan, the three split parts and the three trajectory statistics SHALL equal the committed comparator table exactly. The profile's other readings are descriptive.

#### Scenario: Parent units
- **WHEN** the parent profile enumerates its units
- **THEN** it has two designs of fourteen cells, the anchor and translation at replicates 0–499 and the other cells at 0–99, 4,400 units at 999 permutations

#### Scenario: Reproduction holds
- **WHEN** every row matches the committed counts
- **THEN** the merge records the reproduction as passed with the number of rows compared

#### Scenario: Reproduction fails
- **WHEN** any row differs
- **THEN** the merge records the reproduction as failed and lists every differing row with both counts

### Requirement: The stage-span profile is a fresh seed family on three blocks

The `stage_span` profile SHALL build its units from the latent-space comparison's seed derivations at base seed 2043, at three stages, measuring methylation, expression and proteomics with the span layers to match. Its designs are balanced n = 80, 200 and 500 at the generator's stage effect δ = 2.0, and n = 200 at δ = 0.5.

#### Scenario: Cells per design
- **WHEN** the profile enumerates a design's cells
- **THEN** it has the anchor; translation at e = 1.0; magnitude at e ∈ {0.05, 0.1, 0.25}; orientation at e ∈ {0.1, 0.25, 0.5} with alignment 0; shape at e ∈ {0.1, 0.25, 0.5}; orientation at e = 0.5 with alignment 90; and that cell on a 20-feature support

#### Scenario: Scales
- **WHEN** a cell's plant is built
- **THEN** the cells at alignment 0 use the planted scale and the two cells at alignment 90 use the null scale

#### Scenario: Units and grade
- **WHEN** the profile enumerates its units
- **THEN** each design has its anchor and translation at 500 replicates and its eleven other cells at 200, 12,800 units in all, at 999 permutations

#### Scenario: A fresh family
- **WHEN** a unit's dataset seed is derived
- **THEN** it equals no dataset seed of the latent-space comparison, the split-test validation, the scan validation or the predictive-subspace pilot

#### Scenario: Matched across cells
- **WHEN** two cells of one design and replicate are generated
- **THEN** they share their null dataset

### Requirement: The surgery profile reruns the design-guidance slice's three-block columns and must reproduce it

The `surgery` profile SHALL rebuild, from the committed study configurations, the design-guidance paper-grade slice's three-block columns at n = 80, 200 and 500, with the slice's seeds, at 500 replicates and 999 permutations. For each column, mode, effect and statistic, the `pls_cv` arm's available and rejected counts SHALL equal the slice's committed operating tables exactly.

#### Scenario: Cells per column
- **WHEN** the profile enumerates a column's cells
- **THEN** it has the anchor; translation at e = 1.0; magnitude at e ∈ {0.02, 0.05, 0.1}; and orientation and shape at e ∈ {0.25, 0.5, 1.0}, 16,500 units over the three columns

#### Scenario: Units are the slice's
- **WHEN** a unit's generator parameters and seeds are built
- **THEN** they equal the study engine's for the same column, mode, effect and replicate

#### Scenario: Reproduction
- **WHEN** the profile is merged
- **THEN** the merge records whether every compared row matches and lists every differing row with both counts

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

Merging a profile SHALL write a rate table (per design, cell and arm: available and rejected counts, the rate and its Monte Carlo SE), a paired table (per design, cell and baseline arm against `production`: both rates, the discordant counts and the exact one-sided gain and loss p-values), a cell table (rank quantiles, held-out AUROC, capture shares, hits and discoveries) and accounting (units, failures by cell, duplicates, recorded core-hours).

#### Scenario: Paired table
- **WHEN** a profile is merged
- **THEN** each design, cell and baseline arm has one paired row whose discordant counts are taken over units available in both arms

#### Scenario: No discordant units
- **WHEN** a paired row has no unit where exactly one arm rejects
- **THEN** its paired p-values are 1

#### Scenario: Failed units are excluded
- **WHEN** a cell has failed units
- **THEN** they are counted in the accounting and left out of every rate

### Requirement: The decision judges size per arm and admits only arms that hold it

The decision SHALL apply the stage-span paper grade's count-plus-ceiling bound rule, imported unchanged, to the anchor and translation rates of `production`, `scan`, `fdr` and `total`, in two families: the four `stage_span` designs (m = 8) and the three `surgery` columns (m = 6). An arm SHALL be admissible when both families hold. A baseline arm that is not admissible SHALL leave the detection comparison and be reported.

#### Scenario: An admissible arm
- **WHEN** both families pass the bound rule for an arm
- **THEN** the arm is recorded as admissible, with every rate, bound and status

#### Scenario: A baseline arm fails size
- **WHEN** either family fails for a baseline arm
- **THEN** the arm is recorded as inadmissible with every rate above its bound or ceiling, and no detection label uses it

### Requirement: The decision labels each signal cell by paired exact tests

For each of the 71 signal cells of `stage_span` and `surgery`, and with A admissible baseline arms, the decision SHALL label the cell:
- `production_leads` when, against every admissible arm, the gain p-value is at or below 0.025 / 71 and production's rate is at least 0.10 higher;
- `baseline_leads` when, against at least one admissible arm, the loss p-value is at or below 0.025 / (71·A) and that arm's rate is at least 0.10 higher;
- `unresolved` otherwise.

#### Scenario: Production must beat every admissible arm
- **WHEN** production gains against two admissible arms and not the third
- **THEN** the cell is not labelled `production_leads`

#### Scenario: One arm is enough for the baseline
- **WHEN** one admissible arm beats production at the divided level and by the margin
- **THEN** the cell is labelled `baseline_leads`, naming the arm

#### Scenario: A significant difference below the margin
- **WHEN** a paired p-value passes its level and the two rates differ by less than 0.10
- **THEN** the cell is labelled `unresolved`

#### Scenario: Unavailable readings
- **WHEN** a signal cell has no available unit for an arm the label needs
- **THEN** the cell is labelled `unresolved` and named as unavailable

### Requirement: The decision writes one predeclared verdict

The decision SHALL write `univariate_baseline_decision.json` (rule version 1) with one verdict, the first that applies in this order: `reproduction_failed`, `production_size_excess`, `no_admissible_baseline`, then `production_leads`, `baseline_leads`, `complementary` or `unresolved` by which labels occur in at least one cell. Every criterion SHALL be recorded whatever the verdict, with each verdict's consequence statement.

#### Scenario: Production fails size
- **WHEN** every reproduction passes and `production` is inadmissible
- **THEN** the verdict is `production_size_excess`, naming the family and every rate above its bound or ceiling

#### Scenario: No baseline arm holds size
- **WHEN** every reproduction passes, `production` is admissible and no baseline arm is
- **THEN** the verdict is `no_admissible_baseline`

#### Scenario: Complementary
- **WHEN** at least one cell is labelled `production_leads` and at least one `baseline_leads`
- **THEN** the verdict is `complementary`, naming the cells of each label

#### Scenario: One side only
- **WHEN** cells are labelled `baseline_leads` and none `production_leads`
- **THEN** the verdict is `baseline_leads`

#### Scenario: Neither side
- **WHEN** every signal cell is `unresolved`
- **THEN** the verdict is `unresolved`

#### Scenario: Gate
- **WHEN** a reproduction fails
- **THEN** the verdict is `reproduction_failed`, and size, the labels and the grade are still recorded

### Requirement: The decision grades whether production names the kind of difference

Over the 36 `stage_span` cells at alignment 0 with a magnitude, orientation or shape plant, the decision SHALL take those where `production` rejects in at least 0.80 of units and judge the two off-target statistics of each, as one family, under the bound rule of the size criterion. The grade SHALL be `names_kind` when the family holds, `detects_only` when it fails and `no_detected_cell` when no cell qualifies. The grade SHALL be reported beside the verdict and SHALL NOT change it.

#### Scenario: Off-target statistics
- **WHEN** a detected cell's off-target statistics are taken
- **THEN** they are `angle` and `shape` for magnitude, `delta` and `shape` for orientation, and `delta` and `angle` for shape

#### Scenario: Production detects without naming
- **WHEN** the off-target family fails the bound rule
- **THEN** the grade is `detects_only`, naming every off-target rate above its bound or ceiling

#### Scenario: Surgery cells are not graded
- **WHEN** the grade is computed
- **THEN** no `surgery` cell and no cell at alignment 90 enters it, and their off-target rates are reported

### Requirement: The decision reports descriptive tables that never change the verdict

Beside the verdict the decision SHALL report: the count of cells under each label by plant family and design; every arm's rate per cell, the recorded arms included; per cell the median number of scan hits and of `fdr` discoveries and the realized false discovery proportion where it is defined; the selected rank, held-out AUROC and capture shares per design; and, for the `parent` profile, `fdr` and `production` beside the seven reproduced tests.

#### Scenario: The map
- **WHEN** the decision is written
- **THEN** it has one row per plant family and design with the number of cells labelled `production_leads`, `baseline_leads` and `unresolved`

#### Scenario: The earlier baselines in the same table
- **WHEN** the decision is written
- **THEN** every signal cell reports the omnibus rate of `concat` and of `pca_matched` beside `production`'s and the three baseline arms'

### Requirement: The run is recorded in a dated report

The run SHALL be recorded in a dated report under `docs/reports/`, with aggregates only under `results/univariate-baseline-<date>/`: the merged tables, reproduction and accounting, the decision, and provenance (revision, tree hash, job ids, core-hours, package versions). Per-unit records SHALL stay on the cluster.

#### Scenario: Report contents
- **WHEN** the report is written
- **THEN** it states the verdict, the grade, each arm's size reading, every reproduction, the map, the descriptive tables and what follows from the verdict
- **AND** every number it quotes can be read from the committed aggregates

#### Scenario: The limits are stated beside the claim
- **WHEN** the report states the verdict
- **THEN** it states in the same section that the reading holds at 658 simulated features with Gaussian noise and that no real-noise reading exists on the intended blocks
