# adjusted-attribution-validation Specification

## Purpose
Validates covariate-adjusted orientation attribution at the case study's configuration. On semi-synthetic data with a confounded covariate and a planted orientation, the check tests whether the adjusted attribution removes the covariate's contribution from the drivers, and whether it recovers and stabilizes the planted drivers as well as an oracle that knows the injection.

## Requirements

### Requirement: Datasets at the case study's configuration with a planted orientation

The check SHALL build the case-study construction: the Phase 6 replicate on the astrocyte `cohort_30` table at the full stage effect (δ = 2.0), a covariate `z` confounded with group × stage and injected into a fixed feature subset, and five noise covariates. Seeds SHALL come from a base seed no committed check uses. The cells SHALL be the anchor and orientation at e = 0.5 and 1.0. Settings SHALL be hashed.

#### Scenario: Three cells at paper grade

- **WHEN** the check's units are enumerated with default settings
- **THEN** there are three cells (the anchor, orientation at e = 0.5 and orientation at e = 1.0) × 500 replicates, at 999 permutations, 200 bootstrap replicates and a top-k size of 20

#### Scenario: Seeds are fresh and matched

- **WHEN** units are generated
- **THEN** their generator, covariate, injection-direction, RRPP and bootstrap seeds derive from base seed 2033 and the replicate index only, so the three cells of a replicate share them, and no seed equals one used by a committed check at the same replicate

#### Scenario: Resume under other settings is refused

- **WHEN** an output directory is resumed with a different injection strength, base seed, permutation count, bootstrap count or top-k size
- **THEN** the check exits with an error and appends nothing

#### Scenario: Headroom is checked before running

- **WHEN** an orientation cell would exceed the generator's surgery headroom at the table and δ
- **THEN** the check refuses to start and names the cell

### Requirement: Calibration places the construction where attribution can be judged

Before any paper-grade unit, calibration units SHALL be run at the anchor and at orientation e = 1.0. The injection strength `s_z` SHALL be the only calibrated quantity. It SHALL be raised until the covariate measurably confounds unadjusted attribution. Every iteration SHALL be recorded, with the version of the rule it was judged under. If the oracle cannot recover the planted drivers, the check SHALL stop for the user's decision.

#### Scenario: Confounding criterion

- **WHEN** a calibration iteration completes
- **THEN** it records the mean paired difference in top-k contamination between the unadjusted and oracle arms at the anchor, and the mean paired difference in the contrast's projection on the injection direction in both cells, and the criterion is met when the contamination difference is at least 0.25 and the absolute projection difference is at least 0.2 in both cells

#### Scenario: Recoverability criterion stops the check

- **WHEN** the oracle arm's mean top-k truth precision at orientation e = 1.0 is below 0.5 in a calibration iteration that meets the contamination criterion
- **THEN** the check records the iteration, starts no paper-grade unit, and reports the failed criterion

#### Scenario: Iterations are kept

- **WHEN** calibration is run more than once
- **THEN** the calibration file lists every iteration with its `s_z`, its criteria and whether each was met, and the paper-grade settings name the iteration they adopt

### Requirement: One analysis and three attributions per unit

Each unit SHALL be analysed once through the real-data analysis with the six covariates, attribution, rank floor 2, net displacement and `df`. This is the `adjusted` arm. Its exposed estimator and features SHALL feed two more attributions with the same bootstrap settings: `oracle`, on the features minus the exact injection in analysed coordinates, and `unadjusted`, on the features as analysed. A failed consistency check SHALL stop the unit.

#### Scenario: Oracle removes the injection exactly

- **WHEN** the oracle features are built
- **THEN** they equal the analysed features minus `z` times the injection vector divided by each feature's pooled analysed scale, so the oracle's cell means are the uninjected cell means in the analysed coordinates

#### Scenario: Consistency checks

- **WHEN** a unit is evaluated
- **THEN** the exposed estimator applied to the exposed features reproduces the latent scores to within 1e-10 of the largest absolute score, and the adjusted arm's observed means projected through the estimator equal the analysis's stage means to a relative tolerance of 1e-8

#### Scenario: Failed check stops the unit

- **WHEN** either consistency check fails
- **THEN** the unit raises an error naming the check, and no record line is written for it

### Requirement: Attribution quality is measured against truth, the injection and the oracle

For every arm, component and block, each unit SHALL record:
- top-k truth precision, over min(k, |truth|);
- top-k contamination, the share of injected features that are not truth drivers;
- the contrast's projection on the unit injection direction;
- cosine and top-k overlap with the oracle;
- mean bootstrap sign stability and top-k frequency over the arm's own top-k;
- valid and rank-deficient bootstrap counts.

#### Scenario: Truth is the generator's driver set

- **WHEN** truth precision is computed for the principal-orientation block
- **THEN** the truth set is the union of the generator's per-transition driver features, and the record holds its size

#### Scenario: Analysis outcome is recorded beside the attribution

- **WHEN** a unit completes
- **THEN** its record line also holds the unit key, seeds, the adjusted analysis's statistics, p-values, selected rank and mean held-out CV stage AUROC

#### Scenario: Merge writes the metric tables

- **WHEN** a run is merged
- **THEN** it writes a per-unit metric table and a per-cell summary with, per arm, component and block, the mean and standard error of every metric and of every paired difference against the oracle, and repeating the merge yields identical files

### Requirement: Adjustment must remove the confounded covariate from the drivers

On the observed principal-orientation contrast, the adjusted arm SHALL hold two bounds at every cell, with paired standard errors over the cell's R available units:
- contamination: mean(adjusted − oracle) + 2·SE ≤ 0.05;
- residual confounding: |mean projection difference (adjusted − oracle)| + 2·SE ≤ 0.10 × |mean projection difference (unadjusted − oracle)|.

#### Scenario: Contamination bound

- **WHEN** the adjusted arm's mean top-k contamination exceeds the oracle's by 0.04 with a paired SE of 0.01 at the anchor
- **THEN** the contamination criterion fails at the anchor (0.06 > 0.05)

#### Scenario: Residual-confounding bound

- **WHEN** the unadjusted arm's mean projection exceeds the oracle's by 0.40 and the adjusted arm's by 0.02 with a paired SE of 0.005
- **THEN** the residual-confounding criterion holds (0.03 ≤ 0.04)

### Requirement: Adjustment must cost no recovery or stability against the oracle

On the observed principal-orientation contrast, the adjusted arm SHALL hold three non-inferiority bounds against the oracle at both orientation cells, with paired standard errors:
- truth precision: mean(adjusted − oracle) − 2·SE ≥ −0.05;
- bootstrap sign stability: mean(adjusted − oracle) − 2·SE ≥ −0.05;
- top-k frequency: mean(adjusted − oracle) − 2·SE ≥ −0.10.

#### Scenario: Recovery bound

- **WHEN** the adjusted arm's mean truth precision is 0.03 below the oracle's with a paired SE of 0.005 at e = 1.0
- **THEN** the recovery criterion holds at that cell (−0.04 ≥ −0.05)

#### Scenario: Stability bound fails

- **WHEN** the adjusted arm's mean top-k frequency is 0.10 below the oracle's with a paired SE of 0.01 at e = 0.5
- **THEN** the stability criterion fails at that cell and the failure names top-k frequency

### Requirement: The check must show the covariate confounds unadjusted attribution

At paper grade the decision SHALL first test two confounding gates:
- at the anchor, mean top-k contamination (unadjusted − oracle) − 2·SE ≥ 0.25;
- at every cell, |mean projection on the injection direction (unadjusted − oracle)| − 2·SE ≥ 0.2.

A failed gate SHALL make the verdict `construction_uninformative`, whatever the criteria read.

#### Scenario: Weak confounding is uninformative

- **WHEN** the unadjusted arm's contamination exceeds the oracle's by only 0.10 at the anchor
- **THEN** the verdict is `construction_uninformative` and names the contamination gate and the anchor

#### Scenario: Confounding that misses the contrast is uninformative

- **WHEN** the unadjusted arm's projection on the injection direction differs from the oracle's by −0.15 with a paired SE of 0.01 at orientation e = 1.0
- **THEN** the verdict is `construction_uninformative` and names the projection gate and the cell (0.13 < 0.2)

#### Scenario: The contamination gate does not read the orientation cells

- **WHEN** the unadjusted and oracle arms have equal top-k contamination at both orientation cells and every other gate and criterion holds
- **THEN** the verdict is `adjusted_attribution_validated`

### Requirement: The check must show the planted drivers are recoverable

At paper grade the decision SHALL also test that the oracle's mean top-k truth precision at orientation e = 1.0 is at least 0.5. A failed gate SHALL make the verdict `construction_uninformative`, whatever the criteria read.

#### Scenario: Unrecoverable drivers are uninformative

- **WHEN** the oracle's mean truth precision at orientation e = 1.0 is 0.4 and the adjusted arm fails the contamination criterion at the anchor
- **THEN** the verdict is `construction_uninformative` and names the recoverability gate, not the contamination criterion

### Requirement: The decision issues one verdict with its consequence

A pure decision function SHALL write `adjusted_attribution_decision.json`. It SHALL hold every gate and criterion (value, bound, SE, R, verdict), the settings hash, the adopted calibration iteration, and one verdict with its consequence. The verdict SHALL be indicative below 999 permutations, 200 bootstrap replicates or 500 replicates.

#### Scenario: Validated

- **WHEN** both gates and every criterion hold at paper grade
- **THEN** the verdict is `adjusted_attribution_validated`, and its consequence states that item 12 may start on the released version

#### Scenario: Diagnose

- **WHEN** both gates hold and the contamination criterion fails at orientation e = 0.5
- **THEN** the verdict is `diagnose_adjusted_attribution`, names the failing criterion and cell, and its consequence states that the attribution leg stays withheld and a diagnosis precedes item 12

#### Scenario: Smoke decision is indicative

- **WHEN** the decision is computed on 20 replicates at 99 permutations
- **THEN** the decision file is marked indicative and its verdict is not read as the check's result

#### Scenario: Unadjusted and other components are descriptive

- **WHEN** the decision is written
- **THEN** it reports the unadjusted arm's criteria, the PLS-captured and residual components, the per-transition blocks and the readings conditional on `angle` rejection, and none of them changes the verdict
