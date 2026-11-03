# Spec Delta

## Purpose

Demonstrates, on a real case-study export, that a planted group × stage orientation difference is detected by `motco analyze` at the case-study configuration and traced by its covariate-adjusted attribution back to the planted features. A predeclared rule either lifts Phase 6's withheld attribution leg or reports it failed.

## ADDED Requirements

### Requirement: The check reads the case-study export through a declared input contract

The check SHALL read an export directory holding `manifest.json`, an `rna` and an `atac` block (`identity` transform), `metadata.csv` (sample, group, stage columns, two groups, explicit stage order) and `covariates.csv`. It SHALL derive the primary's six covariate columns and record every input file's sha256 in the settings. A violation SHALL stop the check before any analysis.

#### Scenario: The primary's covariates are derived

- **WHEN** the export is loaded
- **THEN** the covariate columns, in model order, are `age_at_death`, `pmi`, `rin`, `years_education`, `log10_n_rna` and `log10_n_atac`, the last two the base-10 logarithms of `n_rna` and `n_atac`

#### Scenario: Contract violation stops the check

- **WHEN** a block's donors differ from the metadata's, a stage label is outside the declared order, there are not exactly two groups, or a covariate value is missing
- **THEN** the check exits with an error naming the violation and runs no analysis

#### Scenario: Resume on another export is refused

- **WHEN** an output directory is resumed with an export whose file hashes differ from the recorded ones
- **THEN** the check exits with an error and appends nothing

### Requirement: The check never analyses the real group labels

Every analysis SHALL use group labels permuted within stage. The check SHALL NOT compute, record or read any statistic, p-value or attribution under the export's real group labels.

#### Scenario: A draw equal to the real labels is refused

- **WHEN** a within-stage permutation reproduces the real label vector
- **THEN** the draw is discarded and redrawn from the next seed in the replicate's stream, and the redraw is counted in the record

#### Scenario: No record holds a label vector

- **WHEN** a record line is written
- **THEN** it holds the permutation's seed and no donor identifier, label vector or donor-level value

### Requirement: The latent space must reproduce the case-study primary before any unit

Before calibration or any unit, one analysis on permuted labels with no permutations and no attribution SHALL reproduce the primary's label-free latent space: CV-selected rank 5 and mean held-out stage AUROC 0.635556 to six decimals. Every anchor unit SHALL be held to the same check.

#### Scenario: Fidelity gate passes

- **WHEN** the gate's analysis selects rank 5 with mean held-out AUROC 0.635556
- **THEN** the gate's result is written to the output directory and calibration may start

#### Scenario: Fidelity gate fails

- **WHEN** the gate's analysis selects another rank or AUROC
- **THEN** the check exits with an error naming both values, and no calibration or unit runs

#### Scenario: Anchor unit drifts

- **WHEN** an anchor unit's latent rank or AUROC differs from the gate's
- **THEN** the unit raises an error naming the check and no record line is written for it

### Requirement: The null-making step permutes group labels within stage

Each replicate SHALL permute the group labels within each stage, keeping every group × stage cell size. Each donor's blocks and covariates SHALL stay on the donor's own row. Everything the replicate draws (permutation, planted set, signs, RRPP seed, bootstrap seed) SHALL derive from the base seed and the replicate index only.

#### Scenario: Cell sizes are kept

- **WHEN** a replicate's labels are permuted on the astrocyte `cohort_30` export
- **THEN** the group × stage table equals `[[11, 10, 28], [9, 10, 12]]`

#### Scenario: Cells of a replicate are matched

- **WHEN** the anchor and planted units of one replicate are built
- **THEN** they share the permuted labels, the planted set, the signs, the RRPP seed and the bootstrap seed, and differ only in the strength

#### Scenario: Calibration units are not paper-grade units

- **WHEN** calibration and paper-grade units are enumerated
- **THEN** their seeds come from separate streams, so no paper-grade unit repeats a calibration draw

### Requirement: The planted orientation is a declared ramp on a known feature set

The planted set K SHALL be 10 features drawn per replicate from each block, each with a random sign `s_f`. In the second declared group under the permuted labels, each feature f in K SHALL be shifted by `e · c(stage) · σ_f · s_f`. Here `σ_f` is the feature's pooled within-stage SD in export units, and `c` is the stage's position, from 0 at the first stage to 1 at the last.

#### Scenario: Shift formula

- **WHEN** a donor in the planted group at the middle of three stages is built at e = 2.0
- **THEN** each planted feature of that donor equals its export value plus `σ_f · s_f`, and every other feature and donor equals the export

#### Scenario: Zero strength is the export

- **WHEN** a unit is built at e = 0
- **THEN** its blocks equal the export's byte for byte

#### Scenario: The planted set is recorded

- **WHEN** a unit completes
- **THEN** its record holds K's feature names by block and their signs

### Requirement: The strength is calibrated on detection only

Before any paper-grade unit, the check SHALL climb the ladder e ∈ {0.5, 1.0, 1.5, 2.0, 3.0}. Each step SHALL run 40 calibration units at 99 permutations with attribution off. It SHALL adopt as `e*` the first step whose `angle` rejection rate is at least 0.90. Each iteration SHALL be recorded with its rate and outcome, and the settings SHALL name the adopted one.

#### Scenario: First qualifying step is adopted

- **WHEN** `angle` rejects in 0.80 of the calibration units at e = 1.5 and in 0.925 at e = 2.0
- **THEN** `e*` is 2.0, the descriptive cell is 1.5, and the ladder is not climbed further

#### Scenario: Cap reached

- **WHEN** no step reaches 0.90
- **THEN** `e*` is 3.0, the calibration file records `cap_reached`, and the paper-grade run proceeds at the cap

#### Scenario: Lowest step qualifies

- **WHEN** e = 0.5 qualifies
- **THEN** the descriptive cell is e = 0.25

### Requirement: Each unit is analysed once at the case-study configuration

Each unit SHALL be analysed once by the real-data analysis with the six covariates, attribution, rank floor 2, net displacement, the `df` residual correction and the harness CV defaults. That analysis is the primary arm. An unadjusted attribution on the analysis's exposed estimator and features, with the same bootstrap seed, B and top-k, SHALL be recorded as a descriptive arm.

#### Scenario: Paper-grade cells

- **WHEN** units are enumerated with default settings
- **THEN** there are three cells (the anchor, `e*` and the descriptive step) × 500 replicates, at 999 permutations, 200 bootstrap replicates and top-k 20

#### Scenario: Means project to the analysis

- **WHEN** a unit is evaluated
- **THEN** the primary arm's observed attribution means projected through the exposed estimator equal the analysis's stage means to a relative tolerance of 1e-8, or the unit raises an error and writes no record line

### Requirement: Each unit records detection, recovery and stability

On every arm, component and block of each unit, the check SHALL record:
- the planted-set precision |top-k ∩ K| / min(k, |K|);
- the mean bootstrap sign stability and top-k frequency over K;
- the share of recovered planted features whose observed sign is the planted sign.

Each record SHALL also hold the analysis's statistics, p-values, selected rank, AUROC, eigengap flag and valid bootstrap count.

#### Scenario: Merge writes the summary

- **WHEN** a run is merged
- **THEN** it writes a per-cell summary of every metric's mean and SE, over all units and over the units whose `angle` rejects, and repeating the merge yields identical files

#### Scenario: Anchor precision is the chance reference

- **WHEN** the anchor cell is summarized
- **THEN** its planted-set precision is reported as the false-recovery reference beside the planted cells

### Requirement: Criterion (i) — the unplanted construction holds size

At the anchor, the rejection rate (p ≤ 0.05) of each of `delta`, `angle` and `shape` SHALL be at most α + 2·√(α(1 − α)/R), with R the anchor's available units.

#### Scenario: Size holds

- **WHEN** the anchor rates are 0.046, 0.058 and 0.040 at R = 500
- **THEN** criterion (i) holds (bound 0.0695)

#### Scenario: Size fails on one statistic

- **WHEN** the anchor `shape` rate is 0.072 at R = 500
- **THEN** criterion (i) fails and names `shape`

### Requirement: Criterion (ii) — the planted difference is detected

At `e*`, the primary arm's `angle` rejection rate SHALL be at least 0.80.

#### Scenario: Detection fails at the cap

- **WHEN** `e*` is the cap and `angle` rejects in 0.62 of its units
- **THEN** criterion (ii) fails and names the cell and the rate

### Requirement: Criterion (iii) — the planted features are recovered with stable signs

On the units at `e*` whose `angle` rejects, on the primary arm's observed principal-orientation contrast, the lower bound mean − 2·SE SHALL clear each threshold:
- planted-set precision: 0.50;
- sign stability over K: 0.90;
- top-k frequency over K: 0.50.

#### Scenario: Recovery holds

- **WHEN** the means are 0.71, 0.95 and 0.62 with SEs 0.01, 0.005 and 0.01
- **THEN** criterion (iii) holds

#### Scenario: Stability fails

- **WHEN** the mean top-k frequency over K is 0.51 with an SE of 0.01
- **THEN** criterion (iii) fails and names top-k frequency (0.49 < 0.50)

### Requirement: The decision issues one verdict

A pure decision function SHALL write `plasmode_decision.json`, holding every criterion (value, bound, SE, R, verdict), the settings hash, the adopted calibration iteration and one verdict: `attribution_leg_lifted` when (i), (ii) and (iii) all hold, else `attribution_leg_failed` naming every failing criterion and cell.

#### Scenario: Lifted

- **WHEN** criteria (i), (ii) and (iii) hold at paper grade
- **THEN** the verdict is `attribution_leg_lifted`

#### Scenario: Failed on several criteria

- **WHEN** criterion (ii) and the sign-stability part of (iii) fail
- **THEN** the verdict is `attribution_leg_failed` and its failing list names both

#### Scenario: Smoke decision is indicative

- **WHEN** the decision is computed below 500 replicates, 999 permutations or 200 bootstrap replicates
- **THEN** the file is marked indicative and its verdict is not read as the check's result

#### Scenario: Descriptive readings do not move the verdict

- **WHEN** the decision is written
- **THEN** it reports the unadjusted arm, the PLS-captured and residual components, the per-transition blocks, the unconditional readings, `delta` and `shape` at `e*` and the descriptive cell, and none of them changes the verdict

### Requirement: Only aggregates are committed

The committed results SHALL be the settings, the calibration file, the fidelity gate's result, the per-cell summary, the merge and decision files and a PROVENANCE file. Per-unit records SHALL stay on the run machine, named in PROVENANCE with their sha256.

#### Scenario: Results directory content

- **WHEN** the run's results directory is prepared for commit
- **THEN** it holds no per-unit record and no file derived from donor-level values other than per-cell aggregates
