# corrected-permutation-validation Specification

## Purpose
Validates the residual-corrected permutation at the case study's configuration. It covers the size of both corrections on fresh data for both case-study cohort tables, and their power cost at the Phase 6 design. A predeclared rule then names the correction the case study adopts, or states that none holds.

## Requirements

### Requirement: Type I profile reproduces the case-study construction on a fresh seed family

The Type I profile SHALL build the case-study construction from a base seed that no committed check uses. That construction is the Phase 6 replicate, a confounded covariate `z` injected at `s_z = 10`, five noise covariates, floor 2, net displacement and the weak stage signal. It SHALL run each cohort table at δ = 0.1 and 0.05 with the anchor and translation controls, and every hashed settings file SHALL record the table, δ and base seed.

#### Scenario: Four Type I runs

- **WHEN** the Type I profile is configured for the astrocyte table `[[11, 10, 28], [9, 10, 12]]` and the microglia table `[[11, 10, 24], [9, 11, 11]]` at δ = 0.1 and 0.05
- **THEN** four runs are defined, each with the anchor and the translation control at e = 1.0, 500 replicates and 999 permutations by default

#### Scenario: Seeds are fresh and matched

- **WHEN** units are generated for one run
- **THEN** their generator, covariate and RRPP seeds derive from the profile's base seed and the replicate index only, so the two controls of a replicate share them, and no seed equals one used by the committed case-study check at the same replicate

#### Scenario: Resume under other settings is refused

- **WHEN** an output directory written for one run is resumed with a different table, δ, base seed or permutation count
- **THEN** the check exits with an error and appends nothing

#### Scenario: Headroom is checked before running

- **WHEN** a run's translation control would exceed the generator's surgery headroom at its table and δ
- **THEN** the check refuses to start and names the cell

### Requirement: Power profile measures the correction's cost at the Phase 6 design

The power profile SHALL use the astrocyte table at the generator's full stage effect (δ = 2.0), with the same six covariates and seed rule as the Type I profile. Its cells SHALL be the anchor, the translation control at e = 1.0, magnitude at 0.02, 0.05, 0.1 and 0.25, and orientation and shape at 0.25, 0.5, 0.75 and 1.0. It SHALL use 500 replicates and 999 permutations by default.

#### Scenario: Fourteen cells

- **WHEN** the power profile's units are enumerated
- **THEN** there are 14 cells × 500 replicates, and every cell of one replicate shares that replicate's seeds

### Requirement: Every unit is analysed under all three nulls through the real-data analysis

Each unit SHALL be analysed three times through the real-data analysis, under no correction, `df` and `leverage`. The three analyses SHALL use the same covariates, rank floor, estimator, CV settings, permutation count and RRPP seed. Their latent scores, selected rank and statistics SHALL be identical; otherwise the unit SHALL fail with an error naming the arms.

#### Scenario: Arms share the measurement

- **WHEN** a unit is evaluated
- **THEN** its three arms record identical latent-score digests, selected ranks and statistics, and differ only in p-values, null summaries and the recorded correction

#### Scenario: Arm disagreement stops the unit

- **WHEN** the three arms of a unit differ in their selected rank or statistics
- **THEN** the unit raises an error, and no record line is written for it

#### Scenario: Unit record

- **WHEN** a unit completes
- **THEN** its record line holds the unit key, seeds and, per arm, the statistics, p-values, null summary, selected rank, mean held-out CV stage AUROC and the recorded correction block

### Requirement: Size is judged by a pooled bound and a per-reading bound

A reading is one run × control × statistic. A correction SHALL hold size when two conditions are met:
- for each cohort table and statistic, its rate pooled over that table's four readings is at most α + 2·√(α(1 − α)/R_pooled);
- every reading, including the power profile's controls, is at most α + 3·√(α(1 − α)/R).

R is the number of available replicates, and α = 0.05.

#### Scenario: Bounds at the default replicates

- **WHEN** every reading has 500 available replicates
- **THEN** the per-reading bound is 0.05 + 3·√(0.0475 / 500) ≈ 0.0792, and the pooled bound per table and statistic is 0.05 + 2·√(0.0475 / 2000) ≈ 0.0597

#### Scenario: Pooled excess fails the correction

- **WHEN** a correction's `delta` pooled over the microglia table's four readings is 0.062 and every reading is below 0.0792
- **THEN** the correction does not hold size, and the decision names the table and statistic

#### Scenario: Uncorrected excess is reported

- **WHEN** the rule is applied to the uncorrected arm
- **THEN** its outcome is recorded as descriptive, with a flag stating whether the fresh family reproduces a pooled `delta` excess in either table, and it never enters the verdict

### Requirement: A replicate-clustered pooled sensitivity is reported

Because a cohort table's Type I readings share replicate seeds, the decision SHALL also report, per table, arm and statistic, the pooled rate against α + 2·√(α(1 − α)/R_eff), where R_eff is the pooled reading count divided by the design effect of the replicate-clustered rate. This sensitivity SHALL be descriptive and SHALL NOT change the verdict.

#### Scenario: Duplicated readings form one cluster

- **WHEN** every reading of a table rejects on the same 31 of 500 replicates, so its pooled rate is 0.062
- **THEN** D8's pooled tier fails, the clustered sensitivity reports a design effect of 4·500/499 and a bound of α + 2·√(α(1 − α)/499), under which the rate holds, and the verdict is unchanged

### Requirement: A predeclared rule chooses between two size-valid corrections

When both corrections hold size, the check SHALL pair them on the power profile's units. For each target pair (magnitude/`delta`, orientation/`angle`, shape/`shape`) it SHALL pool the mode's four effects and compute a replicate-clustered exact one-sided sign-flip test in each direction. `leverage` SHALL be chosen only with a gain at p ≤ 0.05 on at least one target pair and no loss at p ≤ 0.05 on any; otherwise `df`.

#### Scenario: Leverage gains on orientation

- **WHEN** both corrections hold size, `leverage` rejects orientation/`angle` more often than `df` with a clustered one-sided p of 0.01, and no target pair shows a significant loss
- **THEN** the choice is `leverage`

#### Scenario: No significant difference

- **WHEN** both corrections hold size and no target pair differs significantly in either direction
- **THEN** the choice is `df`

### Requirement: The decision issues one verdict with its consequence

A pure function over the merged rows SHALL write `corrected_permutation_decision.json` with:
- the verdict `adopt_df`, `adopt_leverage` or `no_correction_holds`;
- every reading, pooled rate and bound;
- the paired contrasts and the uncorrected flag;
- the settings hashes, and the consequence for `sead_sexdiff` and MOTCO.

A run at fewer permutations than its settings' default SHALL mark the verdict indicative.

#### Scenario: Only df holds

- **WHEN** `df` holds size and `leverage` does not
- **THEN** the verdict is `adopt_df`, with no paired choice computed, and the consequence states that `sead_sexdiff` re-issues its primaries with `df` under the release that contains it

#### Scenario: Neither holds

- **WHEN** neither correction holds size
- **THEN** the verdict is `no_correction_holds`, it lists every failing table, statistic and reading per correction, and the consequence states that the case study keeps reading `delta` against its realized size and the correction stays unvalidated

#### Scenario: Smoke verdict is indicative

- **WHEN** the decision is computed from runs at 99 permutations
- **THEN** the decision file marks the verdict indicative

### Requirement: The power cost of each correction is reported

For each correction against no correction, and for `leverage` against `df`, the merge SHALL write per power cell and statistic the rejection rates of both arms, the discordant counts and the clustered sign-flip p-values. These SHALL be descriptive and SHALL NOT enter the verdict except through the choice rule.

#### Scenario: Cost table

- **WHEN** the power profile has been merged
- **THEN** a table lists, for every cell, statistic and arm pair, both rates, the units only one arm rejects in each direction, and the one-sided p-values
