# latent-space-comparison Specification

## Purpose

Compares candidate measurement spaces across the alignment of a group difference with the stage span. For each space it measures what the space sees, whether it holds size and what it costs on what the production space already sees, so that a space is chosen for the question it answers rather than for the construction it was tested on.

## Requirements

### Requirement: Comparison evaluates five arms on matched data

Every unit SHALL generate one dataset and evaluate five arms on it with one RRPP seed:
- `pls_cv`, the production measurement;
- `pca_matched`, PCA at that unit's `pls_cv` rank;
- `pca_fixed`, PCA at rank 2(K − 1);
- `pls_residual_pca`;
- `concat`.

Each arm SHALL use net displacement, span capture on and the declared permutations. The generator seed SHALL depend only on the design and replicate.

#### Scenario: Arms differ only in measurement
- **WHEN** a unit is evaluated
- **THEN** all five arms read the same dataset and the same RRPP seed, and every record names its arm, design, cell, replicate, generator seed and RRPP seed

#### Scenario: Production arm equals a direct evaluation
- **WHEN** the `pls_cv` record of a unit is compared with a production evaluation of the same dataset at the same RRPP seed
- **THEN** their latent rank, statistics and p-values are identical

#### Scenario: Matched PCA uses the production rank
- **WHEN** a unit is evaluated
- **THEN** the `pca_matched` arm's retained dimensionality equals the `pls_cv` arm's selected rank on that unit

#### Scenario: Cells at a replicate share their null dataset
- **WHEN** two cells of the same design and replicate are evaluated
- **THEN** their datasets are the same null dataset plus their own plant

### Requirement: Comparison covers three designs and fourteen cells

The comparison SHALL run three designs at three stages on methylation + expression with span layers to match: n = 200 balanced at the generator's stage effect δ = 2.0, the SEA-AD table `[[11, 10, 28], [9, 10, 12]]` at δ = 2.0, and that table at the calibrated δ. Each design SHALL run the fourteen cells of design D10 under the stage-span construction.

#### Scenario: Cells per design
- **WHEN** the comparison enumerates its units
- **THEN** each design has the anchor; translation e = 1.0; magnitude e = 0.05; orientation e ∈ {0.25, 0.5} at φ ∈ {0, 30, 60, 90}; shape e ∈ {0.25, 0.5} at φ = 0; and orientation e = 0.5 at φ = 90 with a 20-feature support

#### Scenario: Pilot grade
- **WHEN** the comparison runs at its declared grade
- **THEN** every cell has 100 replicates at 199 permutations, from base seed 2035

### Requirement: The weak design's stage signal is calibrated to the case study

The comparison SHALL choose the weak design's δ before the run, from a predeclared ascending ladder. Each step runs anchor-only `pls_cv` fits on the SEA-AD table without permutations. The first step whose median held-out stage AUROC lies in [0.60, 0.67] SHALL be adopted. If no step lands there, the step closest to 0.636 SHALL be adopted and marked `band_missed`.

#### Scenario: Calibration is recorded
- **WHEN** calibration runs
- **THEN** every step's δ, replicate count, median and quartiles of the held-out AUROC and median rank are written to the calibration output, with the adopted δ and how it was chosen

#### Scenario: Run requires a calibration
- **WHEN** the run is started without a calibration output or an explicit δ
- **THEN** the driver refuses to start and names the missing input

### Requirement: Runs are resumable, sharded and bound to their settings

The driver SHALL write a settings file with a hash and one record per (unit, arm) to its shard. On a resume it SHALL skip complete units. It SHALL refuse a resume or a merge whose settings differ, naming the differing keys. A failure in one arm SHALL be recorded for that arm without aborting the unit.

#### Scenario: Resume skips completed units
- **WHEN** a shard is restarted after some units completed
- **THEN** only the remaining units are evaluated, and the merged output has one record per unit and arm

#### Scenario: Settings mismatch is refused
- **WHEN** a resume or merge is attempted with settings whose hash differs from the stored one
- **THEN** the driver exits with an error listing the differing keys

### Requirement: Merge writes the comparison tables

Merging SHALL write, per arm, design, cell and statistic, the rejection rate with its Monte Carlo standard error, the median capture shares, the median retained rank and the median `population_latent` target statistic. It SHALL also write the run's accounting: units, failures, duplicates and recorded core-hours.

#### Scenario: Tables cover every arm and cell
- **WHEN** the merge completes on a full run
- **THEN** the rate table has a row for every arm, design, cell and statistic, and the capture table a row for every arm, design and cell

#### Scenario: Failures are counted, not dropped
- **WHEN** some arm records failed
- **THEN** the accounting reports them by arm and cell, and their rates are computed over the available records with the count shown

### Requirement: An advisory reading decides whether a paper grade is proposed

`decide` SHALL read the merged tables and write `latent_space_comparison_reading.json`. It SHALL apply the criteria to `pca_matched`, `pca_fixed` and `pls_residual_pca`, each against `pls_cv`; `concat` is a reference and is never eligible. The bound is α + 2·√(α(1 − α)/R), and comparisons use 2 pooled SE.

#### Scenario: Size criterion
- **WHEN** the reading evaluates an arm's size
- **THEN** it requires the anchor and translation rejection rates of `delta`, `angle` and `shape` at every design to be at or below the bound

#### Scenario: Out-of-span gain criterion
- **WHEN** the reading evaluates an arm's gain
- **THEN** it requires, at both δ = 2.0 designs, the arm's `angle` rate for orientation e = 0.5 at φ = 90 (dense) to exceed `pls_cv`'s by at least 2 pooled SE

#### Scenario: Protected criterion
- **WHEN** the reading evaluates an arm's cost
- **THEN** it requires, at both δ = 2.0 designs, the arm's rate to be at least `pls_cv`'s minus 2 pooled SE for orientation e = 0.5 at φ = 0 on `angle`, shape e = 0.5 on `shape` and magnitude e = 0.05 on `delta`

#### Scenario: Verdicts
- **WHEN** the criteria are evaluated
- **THEN** each eligible arm reads `qualifies` (all three hold), `trade_off` (size and gain hold, protected fails), `no_gain` (size holds, gain fails) or `size_excess` (size fails), naming every failing design, cell and statistic
- **AND** the overall verdict is `propose_paper_grade` when at least one arm qualifies, else `document_tradeoff`

#### Scenario: Descriptive tables
- **WHEN** the reading is written
- **THEN** it includes, for every arm and design, capture and target power against φ, `pls_cv`'s own control rates, the calibrated design's tables and the sparse-minus-dense difference in capture and `angle` power at φ = 90 with its pooled SE

### Requirement: The comparison is reported in a dated report

A dated report under `docs/reports/` SHALL record:
- the configuration and provenance;
- the calibration, accounting and identity checks;
- each arm's visibility profile against φ, its size and its protected power;
- the sparse-versus-dense contrast;
- the reading's verdict and the proposed follow-up.

Aggregates SHALL be committed under `results/latent-space-comparison-<date>/`, and per-unit records SHALL stay off the repository.

#### Scenario: Report numbers match the aggregates
- **WHEN** the report is written
- **THEN** every rate, share and verdict it quotes can be read from the committed aggregates and reading
