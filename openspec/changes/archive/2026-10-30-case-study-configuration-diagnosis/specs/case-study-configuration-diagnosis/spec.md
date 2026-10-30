# Spec Delta

## Purpose

Diagnoses the size excess that the case-study configuration Type I check found: the covariate-adjusted model's `delta` and the `stage` latent adjustment's anchor `angle` and `shape`. It reruns the check's exact data under arms that separate the candidate causes, and issues predeclared verdicts with their consequence for the case study and for MOTCO.

## ADDED Requirements

### Requirement: Diagnosis runs on the parent check's data

The diagnosis SHALL take its construction from a named committed run of the case-study Type I check, one per stage-signal strength. Every unit SHALL regenerate that run's data, covariates, CV settings and RRPP seed. It SHALL refuse a parent whose settings file fails its recorded hash, or whose settings the current construction code does not rebuild to the same hash.

#### Scenario: Unreproducible parent is refused

- **WHEN** the diagnosis is started against a parent run whose settings file fails its hash, or whose settings the current code rebuilds to a different hash
- **THEN** it exits with an error naming the parent and the mismatch, and writes no records

#### Scenario: Resume under another parent is refused

- **WHEN** an output directory written against one parent run is resumed against a different parent
- **THEN** the diagnosis exits with an error and appends nothing

#### Scenario: Units match the parent check

- **WHEN** a unit is generated for a given control and replicate index
- **THEN** its datasets, covariate values and RRPP seed are identical to those the parent check generated for that unit

### Requirement: Judged arms reproduce the committed p-values

For every unit, the six-column covariate-adjusted arm and the `stage` latent-adjusted arm at the floor of 2 SHALL be compared with the committed per-unit p-values of the parent run. The merge SHALL report the number of compared and mismatching p-values per arm and run, and any mismatch SHALL block every verdict.

#### Scenario: Reproduction holds

- **WHEN** every compared p-value equals its committed value
- **THEN** the decision file reports zero mismatches for both arms in both runs and issues its verdicts

#### Scenario: Reproduction fails

- **WHEN** any compared p-value differs from its committed value
- **THEN** the decision is `reproduction_failed`, the file lists the mismatching units, and no diagnostic verdict is issued

#### Scenario: Smoke runs at fewer permutations are not compared

- **WHEN** a run uses fewer permutations than its parent
- **THEN** its reproduction is reported as `not_compared` and its verdicts are marked indicative

### Requirement: Trajectory-model factorial on the fixed covariate-blind basis

On each unit's covariate-blind basis, the diagnosis SHALL run six trajectory-model arms sharing the RRPP seed: the confounded covariate estimated (`z_noise0`, `z_noise2`, `z_noise5`) or its injected effect removed exactly from the latent scores (`oracle_noise0`, `oracle_noise2`, `oracle_noise5`). Each is crossed with 0, 2 or 5 independent noise covariates, nested so that fewer columns are a prefix of more.

#### Scenario: Case-study model is one cell

- **WHEN** a unit is evaluated
- **THEN** `z_noise5` has the six design columns of the parent check's covariate-adjusted arm, and its statistics and p-values equal that arm's analysis of the same unit

#### Scenario: Noise covariates are nested

- **WHEN** `z_noise2` and `z_noise5` are evaluated on one unit
- **THEN** the two noise covariates of `z_noise2` are identical to the first two of `z_noise5`

#### Scenario: Oracle removes exactly the injected effect

- **WHEN** an oracle arm's latent scores are formed for a unit
- **THEN** they equal the projection, through the unit's fitted scalers and covariate-blind basis, of the same features without the covariate injection, and the unit fails with an error if the reconstruction of the fitted projection does not reproduce the analysis's latent scores

#### Scenario: Basis is shared

- **WHEN** the six factorial arms of a unit are evaluated
- **THEN** they use the same latent basis, the same selected rank and the same permutation draws

### Requirement: Residual-df corrected permutation arms

For each factorial arm the diagnosis SHALL also run a corrected arm. It permutes the reduced model's residuals scaled by √(n / (n − p_red)), where p_red is the reduced model's column count, with the same permutation indices and the uncorrected observed statistic.

#### Scenario: Observed statistic is unchanged

- **WHEN** a corrected arm and its uncorrected arm are evaluated on one unit
- **THEN** their observed `delta`, `angle` and `shape` are equal and only their null distributions differ

#### Scenario: Correction scales with the reduced model

- **WHEN** a corrected arm's scale factor is recorded
- **THEN** it equals √(n / (n − p_red)) for that arm's reduced model, and it is larger for `z_noise5` than for `oracle_noise0`

### Requirement: Latent-adjustment rank arms

Each unit SHALL also be analysed with the `stage` latent adjustment and the six covariates in three arms that share the data and RRPP seed: CV rank at the floor of 2 (`latent_stage_cv`, the parent check's arm), CV rank at the floor of 1 (`latent_stage_floor1`) and a forced rank of 2 without cross-validation (`latent_stage_rank2`). Each SHALL record its selected rank and component selection.

#### Scenario: Forced rank skips cross-validation

- **WHEN** a unit is analysed in `latent_stage_rank2`
- **THEN** its record states a forced component selection at rank 2 and no CV AUROC

#### Scenario: Floor-1 arm may select rank 1

- **WHEN** `latent_stage_floor1` selects rank 1 on a unit
- **THEN** the unit's rank is recorded as 1, and its `angle` p-value is recorded and counted in the arm's rates

#### Scenario: Failed arm is unavailable

- **WHEN** an analysis raises an error in one arm of a unit
- **THEN** that arm's record carries the error and no p-values, and the unit's other arms are recorded normally

### Requirement: Reference row is excluded by invariance

The diagnosis SHALL not run arms for the LS-mean reference row. A test SHALL show that the covariate-adjusted model's statistics and RRPP p-values are unchanged when the reference row is shifted.

#### Scenario: Shifted reference row

- **WHEN** a unit's covariate-adjusted statistics and RRPP p-values are computed with the pooled-mean reference row and with a shifted reference row, using the same seed
- **THEN** the statistics agree to numerical tolerance and the p-values are identical

### Requirement: Rates, paired contrasts and rank tables

The merge SHALL report every arm's rejection rate at α = 0.05, with R, bound and verdict, per run × control × statistic, using the parent check's one-sided bound. Each predeclared contrast SHALL be an exact one-sided paired sign-flip test on `delta` rejections, clustered by replicate over both runs and controls. Anchor rates by selected rank (2, 3, ≥ 4) SHALL be reported for `z_noise5` and the latent arms.

#### Scenario: Bound matches the parent check

- **WHEN** an arm has 500 available replicates in a run and control
- **THEN** its bound is 0.0695 and its verdict is `holds` at 34 or fewer rejections and `exceeds` above

#### Scenario: Paired contrast counts discordant units

- **WHEN** a contrast compares arm B with arm A
- **THEN** it reports the readings rejected by B only and by A only, and the exact one-sided p-value of the sign-flip test that B rejects more often than A, in which each replicate's readings across runs and controls flip together

#### Scenario: Correlated readings are one cluster

- **WHEN** the anchor and translation readings of one replicate index in both runs enter a contrast
- **THEN** they contribute a single signed discordance count to the test, not four independent pairs

#### Scenario: Rank table covers every selected rank

- **WHEN** the rank table is written for an arm
- **THEN** every available anchor unit falls in exactly one rank bucket, and each bucket reports its count, rejections and rate per statistic

### Requirement: Verdict on the covariate-adjusted delta

A pure function SHALL issue one verdict. `residual_df_shrinkage`: an uncorrected column contrast (`oracle_noise5` vs `oracle_noise0`, or `z_noise5` vs `z_noise0`) is significant at one-sided 0.05, corrected `z_noise5` holds every statistic at both controls in both runs, and the corrected `z_noise5` vs `z_noise0` contrast is not significant. Else `collinear_covariate` if `z_noise0` vs `oracle_noise0` is significant. Else `unresolved`.

#### Scenario: Shrinkage explains the excess

- **WHEN** `z_noise5` rejects `delta` significantly more often than `z_noise0`, the corrected `z_noise5` holds the bound everywhere, and the corrected contrast is not significant
- **THEN** the verdict is `residual_df_shrinkage`

#### Scenario: Estimating the confounded covariate costs size

- **WHEN** the shrinkage conditions fail and `z_noise0` rejects `delta` significantly more often than `oracle_noise0`
- **THEN** the verdict is `collinear_covariate`

#### Scenario: Nothing separates

- **WHEN** neither the shrinkage conditions nor the collinear contrast hold
- **THEN** the verdict is `unresolved`, and the decision file names the R = 1,000 replication as the next step

### Requirement: Verdict on the latent-adjusted angle and shape

The same function SHALL issue `rank_selection` when `latent_stage_rank2` holds the bound on every statistic at both controls in both runs, and `residualized_basis` otherwise. The floor-1 arm's rates and its rank-1 fraction SHALL be reported as descriptive and SHALL NOT enter the verdict.

#### Scenario: Fixed rank holds

- **WHEN** `latent_stage_rank2` is at or below the bound for every statistic at both controls in both runs while `latent_stage_cv` exceeds at the anchor
- **THEN** the verdict is `rank_selection`

#### Scenario: Fixed rank still exceeds

- **WHEN** `latent_stage_rank2` exceeds the bound on `angle` or `shape` at the anchor in either run
- **THEN** the verdict is `residualized_basis`, and the file names each exceeding control, statistic and run

### Requirement: Diagnosis is reported with its consequences

The paper-grade diagnosis SHALL be written as a dated report under `docs/reports/` that cites its decision file and states each verdict's consequence. `residual_df_shrinkage` calls for a MOTCO change making the corrected permutation available, validated before the case study uses it. `collinear_covariate` and `unresolved` name the next diagnostic. The latent verdict states whether the `stage` arm is admissible at a fixed rank.

#### Scenario: Report states the case-study condition

- **WHEN** the report is written
- **THEN** it reads both verdicts and the reproduction counts from the decision file, states what changes for the `sead_sexdiff` primary and its `stage` sensitivity arm under them, and carries the parent check's four limitations
