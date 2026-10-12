# orientation-measurement-probe Specification

## Purpose

An exploratory, matched-seed probe that compares orientation measurement designs (latent space × orientation estimator) at the Phase 6 small-n and Phase 5 large-n designs. It decides which designs advance to a paper-grade comparison study.

## Requirements

### Requirement: Probe crosses measurement spaces with orientation estimators on matched data

The system SHALL provide a probe driver that evaluates every arm on the same generated datasets. An arm is a measurement space × an orientation estimator.

- **Spaces:** stage-supervised PLS with CV rank (the production reference); stage-supervised PLS at forced rank `2(n_stages − 1)`; `pls_within_group` at `n_stages − 1` components per group; `pca` at `2(n_stages − 1)` components; and `concat`.
- **Estimators:** PC1 and net displacement.

For each design, replicate and cell, the driver SHALL generate the dataset once and evaluate every arm on it. It SHALL use one RRPP seed across arms, so that all arms see the same permutation draws. `concat` SHALL be reported as a ceiling reference and SHALL never be eligible to advance.

#### Scenario: Arms share data and permutations

- **WHEN** the probe evaluates one replicate of one cell
- **THEN** every arm's record carries the same generator seed and the same RRPP seed

#### Scenario: Concat is reference-only

- **WHEN** the advancement decision is computed
- **THEN** `concat` arms are reported with their rates but are never listed as advancing

### Requirement: Probe covers two designs and five cells

The probe SHALL run two designs. Each design's generator parameters SHALL be recorded in the outputs.

- **Phase 6 SEA-AD design:** three stages, group × stage sizes `[[11, 10, 28], [9, 10, 12]]`, methylation and expression measured, joint magnitude construction, `p_dmp = 0.1`, ρ = 0.
- **Phase 5 design:** n = 1200, four stages, all three blocks measured, joint magnitude construction, `p_dmp = 0.1`, ρ = 0.

For each design it SHALL run five cells: the zero-effect anchor, orientation at e = 0.5 and e = 1.0, magnitude at e = 1.0, and shape at e = 1.0. The defaults SHALL be 100 replicates and 199 permutations per cell. Surgery censoring SHALL be `error`.

#### Scenario: Default probe size

- **WHEN** the probe runs with default settings
- **THEN** each design × cell × arm has 100 replicates evaluated with 199 permutations

#### Scenario: Phase 5 deviation is recorded

- **WHEN** the probe writes its outputs
- **THEN** they state that the Phase 5 design uses the joint magnitude construction, not the methylation-only construction pinned in the Phase 5 paper-grade config

### Requirement: Probe applies a predeclared advisory advancement rule

The probe SHALL compute a machine-readable advancement decision. The reference arm is stage-PLS with CV rank × the PC1 estimator. The pooled Monte Carlo SE of a difference SHALL be computed from the two arms' binomial rates. The anchor bound SHALL be α + 2·sqrt(α(1 − α)/R) at α = 0.05 with R replicates.

A non-reference, non-`concat` arm SHALL advance only if all of the following hold:

1. **Anchor.** At both designs, its anchor rejection rate for each of `delta`, `angle` and `shape` is at or below the bound.
2. **Gain.** At the Phase 6 design, its orientation `angle` rate at e = 1.0 exceeds the reference's by more than 2 pooled SE.
3. **No loss.** At the Phase 5 design, its orientation `angle` rate at e = 1.0 is not below the reference's by more than 2 pooled SE.
4. **Protected rates.** At both designs, neither magnitude/`delta` nor shape/`shape` at e = 1.0 is below the reference by more than 2 pooled SE.
5. **Specificity.** At both designs, neither magnitude/`angle` nor shape/`angle` at e = 1.0 exceeds the reference by more than 2 pooled SE.

The decision SHALL list, for every arm, each criterion's observed values, thresholds and pass/fail status. It SHALL be advisory: it informs the follow-up study's arm selection and changes no production default.

#### Scenario: An arm that inflates the anchor does not advance

- **WHEN** an arm's anchor `angle` rejection rate at either design exceeds the bound
- **THEN** the decision marks the anchor criterion failed and the arm as not advancing, whatever its power

#### Scenario: An arm that buys orientation with specificity does not advance

- **WHEN** an arm raises orientation `angle` power at the Phase 6 design but its magnitude/`angle` rate exceeds the reference by more than 2 pooled SE
- **THEN** the decision marks the specificity criterion failed and the arm as not advancing

#### Scenario: No arm advances

- **WHEN** no arm satisfies every criterion
- **THEN** the decision records an empty advancing set and the per-criterion reasons for every arm

### Requirement: Probe outputs are reproducible and self-describing

The probe SHALL write the following to an output directory:

- **Per-replicate records:** generator seed, RRPP seed, arm, cell, design, retained dimensionality, observed statistics, p-values, and the realized `angle` at the population, observed-standardized and latent checkpoints.
- **Per-arm rate table:** design, cell and arm, with rates and Monte Carlo SEs.
- **The advancement decision**, in machine-readable form.
- **A Markdown summary**, including the checkpoint table per arm.

Runs SHALL be resumable at replicate granularity: a resumed run SHALL skip completed (design, cell, replicate) units and SHALL refuse to resume into an output directory written with different probe settings. The committed results directory SHALL include a provenance file with the command, code revision, software versions, host and unit timings.

#### Scenario: Interrupted run resumes

- **WHEN** the probe is re-run into an output directory with completed units and the same settings
- **THEN** it evaluates only the missing units, and the final outputs equal those of an uninterrupted run

#### Scenario: Mismatched resume is refused

- **WHEN** the probe is re-run into an output directory written with different settings
- **THEN** it exits with an error naming the differing settings

### Requirement: Motivating evidence is reproducible from committed study records

The system SHALL provide a script that reads a trajectory power study `merged.jsonl` and reproduces the motivating evidence. It SHALL produce two tables:

- **Checkpoint table:** the median realized orientation `angle` at the population, standardized population, observed standardized and PLS latent checkpoints per orientation cell, split by selected rank above or at most `n_stages`, with that split's `angle` rejection rate (rejection at p < α, the study convention).
- **Flip-tail table:** the anchor `angle` quantiles, and the fraction of replicates above the anchor's 95th percentile for the signed and folded angles in the observed and latent spaces.

#### Scenario: Phase 6 evidence reproduces

- **WHEN** the script is run on the Phase 6 small-n orientation profile's `merged.jsonl`, baseline column
- **THEN** it reports population angles of arccos(1 − e) within 0.5° at e = 0.25 / 0.5 / 0.75 / 1.0, and the PLS latent medians and rank-split rejection rates recorded in the probe's findings report
