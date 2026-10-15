# Spec Delta

## ADDED Requirements

### Requirement: Column decision is predeclared in the configuration and evaluated against the baseline column

A study configuration MAY declare an advisory column decision rule that judges one **candidate** design column against the baseline column on the matched power grid. The rule SHALL name the candidate by design-axis coordinates that resolve, after the same value normalization the design grid applies, to exactly one non-baseline design point; a target mode × statistic pair with a criterion kind of `gain` or `no_loss`; a gain standard-error multiplier and a loss standard-error multiplier; a Type I bound (significance level and standard-error tolerance); a non-empty list of **mandatory control** mode × statistic pairs; and a list of **protected** mode × statistic pairs. Every named mode MUST be enumerated by the study, the target pair MUST NOT also be a control or protected pair, and the loader MUST reject a candidate that matches no design point, more than one, or the baseline.

When declared, the report SHALL measure, in the baseline column and the candidate column, each named mode's rejection rate at the largest enumerated effect of that mode and each column's zero-effect anchor rate, with Monte Carlo standard errors and available replicate counts, and SHALL evaluate four criteria for the candidate:

- **target** — under `gain`, the candidate's target rate exceeds the baseline's by more than the gain multiplier times the pooled standard error; under `no_loss`, the baseline's target rate exceeds the candidate's by at most the loss multiplier times the pooled standard error;
- **anchor** — the candidate's zero-effect rate is at or below α + k·√(α(1 − α)/R) on `delta`, `angle` and `shape`, with α and k the rule's Type I bound and R the available anchor replicates;
- **controls** — every mandatory control pair's rate in the candidate column, at that mode's largest effect, is at or below the same one-sided bound computed with that pair's available replicates;
- **protected** — no protected pair's rate in the candidate column falls below the baseline's by more than the loss multiplier times the pooled standard error.

The verdict SHALL be `candidate_qualifies` when every criterion passes, `candidate_fails` when any criterion fails, naming the failing criteria and pairs, `unavailable` when no criterion fails but one cannot be evaluated, `no_candidate` when no record carries the candidate's coordinates, and `no_design_grid` when no design-point record is present. The output MUST list every rate, standard error, replicate count, threshold and bound each criterion used, MUST report the same measures for every other design column as descriptive context without a verdict, and MUST be written as a JSON document and a CSV table. The rule MUST NOT feed the Phase 4 gate or the acceptance-target report, and a configuration without the rule MUST produce outputs byte-identical to before the rule existed.

#### Scenario: Candidate qualifies on a gain criterion

- **WHEN** the rule declares `gain` on orientation/`angle` with multiplier 2, the candidate's rate at the top orientation effect is 0.670 and the baseline's is 0.300 over 500 replicates each, the candidate's anchor is at or below 0.0695 on every statistic, every mandatory control is at or below its bound, and no protected pair loses more than 2 pooled SE
- **THEN** the target criterion passes (the gain 0.370 exceeds 2 × 0.0294) and the verdict is `candidate_qualifies`

#### Scenario: Mandatory control above the bound fails the candidate

- **WHEN** the rule names translation/`angle` as a mandatory control and the candidate's translation `angle` rate at the top translation effect is 0.080 over 500 available replicates at α = 0.05 and k = 2
- **THEN** the controls criterion fails naming translation/`angle` against the bound 0.0695, and the verdict is `candidate_fails` whatever the other criteria show

#### Scenario: Candidate qualifies on a no-loss criterion

- **WHEN** the rule declares `no_loss` on orientation/`angle` with loss multiplier 2, the baseline's rate at the top orientation effect is 0.860 and the candidate's is 0.840 over 500 replicates each, and the other three criteria pass
- **THEN** the target criterion passes (the loss 0.020 is within 2 pooled SE) and the verdict is `candidate_qualifies`

#### Scenario: Candidate coordinates must resolve to one non-baseline point

- **WHEN** a configuration declares a column decision whose candidate coordinates match no design point, or match the baseline point, or match more than one point
- **THEN** loading the configuration raises a clear validation error naming the candidate coordinates and the problem

#### Scenario: Rule absent leaves every output unchanged

- **WHEN** a configuration declares no column decision rule
- **THEN** no column decision output is written and every other report output is byte-identical to before the rule existed

### Requirement: Study provides fixed net-displacement power comparison profiles

The study SHALL provide two version-controlled net-displacement power comparison profiles that measure the four trajectory modes at paper-grade size across two measurement columns on matched datasets. Both SHALL use the numpy generator, `p_dmp = 0.1`, baseline continuity `0.0`, the joint magnitude construction, the default fail-loud surgery-censoring policy, the trajectory modes `magnitude`, `orientation`, `shape` and `translation` with effect sizes `0.0`, `0.25`, `0.5`, `0.75` and `1.0`, 500 replicates per cell, 999 RRPP permutations, matched seeds with a shared zero-effect anchor, attribution disabled, the Phase 4 gate disabled, no power target, and Type I and specificity targets (translation × `delta`/`angle`/`shape`, magnitude × `angle`/`shape`) at α = 0.05 with a two-SE tolerance in the **upper** bound form. The design grid SHALL be declared in the explicit form over the axes `evaluation.integration_method` and `evaluation.orientation_estimator`, with the points: the baseline (cross-validated stage PLS, PC1) and cross-validated stage PLS with net displacement (the candidate). Both SHALL declare the column decision rule with the candidate at the net-displacement stage-PLS point, target orientation/`angle`, mandatory controls translation/`angle`, magnitude/`angle` and shape/`angle`, protected magnitude/`delta` and shape/`shape`, Type I bound α = 0.05 with k = 2, and multipliers 2. The **Phase 6** profile SHALL use three stages, `group_stage_sizes = [[11, 10, 28], [9, 10, 12]]`, the methylation and expression blocks, the Phase 6 small-n study's base seed and matched-seed family, and the `gain` criterion. The **Phase 5** profile SHALL use `n_samples = 1200`, four stages, `group_ratio = 0.5`, all three blocks, the Phase 5 primary base seed and matched-seed family, and the `no_loss` criterion. Neither profile may set `surgery_censoring` to `clamp`.

#### Scenario: Phase 6 comparison profile is loaded

- **WHEN** the committed Phase 6 net-displacement power comparison profile is loaded and enumerated
- **THEN** it deterministically enumerates two Type I baseline cells, the primary anchor and sixteen primary power cells, and one further design column of one anchor and sixteen power cells, thirty-six cells in all, at 500 replicates and 999 permutations
- **AND** its baseline cells draw the same generator seeds as the Phase 6 small-n study's cells of the same mode and effect at every replicate index
- **AND** every pool-limited cell fits its expected surgery headroom

#### Scenario: Phase 5 comparison profile is loaded

- **WHEN** the committed Phase 5 net-displacement power comparison profile is loaded and enumerated
- **THEN** it deterministically enumerates the same thirty-six-cell structure at `n_samples = 1200` and four stages
- **AND** its baseline cells draw the same generator seeds as the Phase 5 primary family's cells of the same mode and effect at every replicate index, and the first replicates' orientation datasets equal those generated from the Phase 5 paper-grade profile at the same effect

#### Scenario: Committed profiles are identity-pinned

- **WHEN** either profile is enumerated
- **THEN** its cell identifiers, parameter signatures and matched seeds equal the committed identity fixture, and every previously committed configuration's fixture entry is unchanged

### Requirement: Comparison selection is computed from both profiles' column decisions

The study SHALL provide a comparison selection that reads the column decision outputs of the Phase 6 and Phase 5 net-displacement power comparison profiles and writes a machine-readable decision. The verdict SHALL be `adopt_net_displacement` when the candidate's verdict is `candidate_qualifies` at both designs and `keep_pc1` otherwise. The output MUST name, per design, the candidate's verdict and the failing criteria and pairs when it does not qualify, and MUST state that the selection is advisory, the rule version, and the two input files it read with their content hashes.

#### Scenario: Candidate qualifies at both designs

- **WHEN** the candidate's verdict is `candidate_qualifies` in both inputs
- **THEN** the verdict is `adopt_net_displacement`

#### Scenario: Candidate fails at one design

- **WHEN** the candidate's verdict is `candidate_fails` or `unavailable` in either input
- **THEN** the verdict is `keep_pc1` and the output names, per design, the failing criteria and pairs

### Requirement: Study produces a dated net-displacement power comparison findings report

The study SHALL produce a dated findings report for the net-displacement power comparison run, tied to the exact configurations, code revision, parameter signatures, software versions, record counts, failure counts and reproduction commands. The report SHALL present, per design: each mode's `angle` rejection rate at every effect under both estimators with Monte Carlo standard errors, with the probe's concatenation ceiling cited as descriptive context; the anchor and every mandatory control against the bound; the recorded eigengap and `angle` null-width covariates per column; the column decision and the comparison verdict. It SHALL state that `delta` and `shape` cannot differ between the two stage-PLS columns, because the estimator enters only the orientation statistic, and report whether their p-values are identical record for record. It SHALL state that the baseline column re-evaluates the datasets of the corresponding committed studies with the same evaluation parameters and RRPP seed, and report whether its records reproduce those studies' recorded p-values exactly wherever a parent record of the same mode and effect exists. When the verdict is `adopt_net_displacement`, the report SHALL state the revised case-study interpretability statement for `angle` that the measured operating characteristics support and SHALL name the production-default change it selects; otherwise it SHALL state that the Phase 6 interpretability statement and the production default are unchanged.

#### Scenario: Report is tied to its run

- **WHEN** the findings report is written
- **THEN** it names both configurations, the code revision, the software versions, the record and failure counts, the committed results directories, and the commands that reproduce the merge, the report and the comparison selection

#### Scenario: Report names the branch taken

- **WHEN** the comparison selection has been computed
- **THEN** the report quotes its verdict, states the interpretability consequence for the case study, and the roadmap entry for item 2 is updated to the branch taken
