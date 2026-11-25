# Proposal

## Why

RRPP exchanges the reduced model's residuals between the two groups, so its null assumes that the groups share one noise dispersion. Every simulated reading draws both groups from one covariance, and nothing has tested the assumption (roadmap open item 8).

- **The literature predicts a failure where the case studies sit.** A permutation test of location can lose its size when the groups' dispersions differ and the cells are unbalanced, and keep it when they are balanced (Anderson and Walsh 2013, doi:10.1890/12-2010.1). Both case-study tables are unbalanced (SEA-AD F 11 / 10 / 28 against M 9 / 10 / 12; ADNI F 237 / 310 / 125 against M 223 / 448 / 158), and unbalance has been read at equal dispersion only: on those two tables, and at pilot grade for group ratios of 0.65 and 0.8.
- **Dispersion by stage already moves the null.** On Gaussian matrices the unrestricted draws read 0.024–0.102 (roadmap item 20). Dispersion by group has no reading.
- **The generator cannot produce it.** Both groups are sampled from the reference covariance. The design-guidance pilot deferred the covariance and correlation scaling rows of the originally planned grid to this change.

An exploratory probe on 2026-11-24 (`notes.md`; not a result) rescaled group B's noise about its exact cell means on null datasets of the committed generator. It read production's test at n = 200, three stages and two blocks, 100 replicates × 199 permutations per cell, at seeds outside every run. With s the multiplier of group B's noise standard deviation:

| stage effect δ | group A's share | s | `delta` | `angle` | `shape` | median CV rank |
|---|---|---|---|---|---|---|
| 2.0 | 0.5 | 1 | 0.03 | 0.01 | 0.02 | 2 |
| 2.0 | 0.5 | 2 | 0.17 | 0.02 | 0.01 | 3 |
| 2.0 | 0.8 | 2 | 0.44 | 0.16 | 0.10 | 3 |
| 2.0 | 0.2 | 2 | 0.01 | 0.00 | 0.00 | 3 |
| 0.5 | 0.5 | 1 | 0.05 | 0.04 | 0.03 | 4 |
| 0.5 | 0.5 | 2 | 0.65 | 0.16 | 0.07 | 9 |
| 0.5 | 0.8 | 2 | 0.71 | 0.38 | 0.10 | 5 |
| 0.5 | 0.2 | 2 | 0.10 | 0.01 | 0.00 | 12 |

- **Equal groups are not protected.** With group B at twice the noise standard deviation, `delta` rejects in 0.17 of balanced null datasets at δ = 2.0, where `angle` and `shape` hold. A group's estimated path is lengthened by its own noise, while a draw gives both groups the pooled noise.
- **It grows as the stage effect weakens.** At δ = 0.5, where the held-out stage AUROC is still above 0.99, the same scale makes `delta` reject in 0.65 of balanced null datasets, and the selected rank rises from 4 to 9.
- **The wider group being the smaller is worst.** With group B a fifth of the samples the three statistics read 0.44 / 0.16 / 0.10 at δ = 2.0 and 0.71 / 0.38 / 0.10 at δ = 0.5.
- **The wider group being the larger is conservative at δ = 2.0** (0.01 / 0.00 / 0.00), and not for `delta` at δ = 0.5 (0.10).
- **Within-group draws repair part of it.** A second probe (three blocks, a forced rank, 200 × 199) permuted the residuals within each group, which the RRPP primitive already allows. At δ = 2.0 with the wider group the smaller it brought 0.37 / 0.19 / 0.065 back to 0.03 / 0.07 / 0.00. At δ = 0.5 with equal groups it left `delta` where it was, 0.53 against 0.515. So part of the excess is the exchange of residuals between groups, and part is carried by the fixed scores, which no residual permutation refits.
- **Correlation scaling is a dispersion scale in the measured space.** Shrinking group B's correlations at fixed variances lowers the within-cell variance of its scores to 0.67 of group A's at half the correlations and to 0.28–0.35 with none (δ = 2.0), so the same setting covers it in total dispersion.

So a rejection, of `delta` above all, may be a difference in dispersion and not in trajectory. Neither case study rejects, so no committed claim is at stake. What is at stake is what a future rejection means and how unequal two groups may be before it matters. The setting also serves open item 3: its bias correction assumes one covariance for both groups, and its sketch names this generator field as the test of that limit.

## What Changes

- **An opt-in dispersion scale in the generator.** `SemiSyntheticTrajectoryParams.group_dispersion_scale` multiplies the noise standard deviation of the second group, in every block and at every stage. The means do not move.
  - It draws no random numbers: the dataset is the equal-dispersion dataset at the same seed with group B stretched about its exact cell means, methylation in M-values. `generate_omics` gains the opt-in `noise_scale` that does it.
  - `None` is the default and 1 normalizes to it. While absent, every dataset, truth record, cell id, signature, matched seed and dumped config is unchanged (`grid._SIGNATURE_OMIT_WHEN_NONE`). Zero-effect cells keep it.
  - Under `stage_span`, the null and planted pooled scales use each group's own noise variance, so the construction stays exact when one group is wider.
  - MOTCO 0.26.0.
- **A validation driver.** `scripts/group_dispersion_check.py` (`run --profile core|columns` | `merge --profile …` | `decide`; sbatch array beside it). Every unit is a null unit of the univariate baseline, imported unchanged: one dataset, one production fit, `delta`, `angle` and `shape`, then the split test and the scan on the same fit.
- **Profile `core`: the cross at n = 200, δ = 2.0.** Three blocks, three stages, base seed 2043 with new design names; the anchor and the planted translation at 500 × 999. 18,000 units.
  - **The gate.** The univariate baseline's `n200` controls are rerun and must reproduce their committed counts exactly, which also checks at run time that the generator is unchanged while the scale is absent.
  - **17 new designs.** s ∈ {1.25, 1.5, 2} (variance ratios 1.56, 2.25 and 4) at a first-group share of 0.5; of 0.65 and 0.8, the pilot's ratios, where the wider group is the smaller; and of 0.35 and 0.2, where it is the larger. Equal dispersion at shares 0.65 and 0.8, the paper-grade reading of unbalance alone.
- **Profile `columns`: a reduced cross at four other committed designs.** The univariate baseline's n = 80 and n = 500 at δ = 2.0 and n = 200 at δ = 0.5, and the stage-effect scaling run's `band_n320` (δ = 0.05 at n = 320, held-out AUROC 0.664, the case studies' regime). Per column: s ∈ {1.5, 2} at shares 0.5, 0.8 and 0.2, and equal dispersion at share 0.8; the anchor at 500 × 999. Each column's balanced equal-dispersion reading is the committed one. 14,000 units.
- **A recorded arm: within-group draws.** On every unit, `delta`, `angle` and `shape` under `RRPP(..., strata=<group>)` on the same fit, statistics and seed, through the primitive's existing option. It decides nothing; it says how much of each column's excess belongs to the null.
- **A predeclared decision.** `group_dispersion_decision.json` (rule version 1):
  - **Gate.** `reproduction_failed` if the rerun design differs from its committed counts.
  - **Size, in four families** of `delta`, `angle` and `shape` at p < 0.05 under the count-plus-ceiling rule (`judge_rates`, imported unchanged): `unbalanced_equal` (m = 24), `balanced_unequal` (m = 42), `wider_smaller` (m = 60) and `wider_larger` (m = 60).
  - **Verdict**, the first that applies: `unbalance_excess` | `dispersion_excess` (balanced designs fail) | `unbalanced_dispersion_excess` | `size_holds`.
  - **Beside it, changing nothing:** for each statistic, column and share, the smallest scale that fails; every conservative reading; a label per family for `production`, `scan`, `fdr`, `total`, `span`, `complement`, `pca_matched` and the within-group statistics; the rank and AUROC by design.
- **Documentation and records.** A dated report `docs/reports/group-dispersion-<date>.md`; aggregates under `results/group-dispersion-<date>/`; a section of `docs/measurement-visibility.md` on what the tests assume about group dispersion and how a cohort is placed; `docs/api/simulations.md`; roadmap open item 8 and the interpretability statement; `AGENTS.md`; `CHANGELOG.md`.

**Non-goals.**
- **Power.** Every unit is a null unit. What a conservative null costs, and power under unequal dispersion, are not read.
- **A remedy as an option.** Within-group draws stay a call to the primitive. Adding them to the harness or to `motco analyze` is a later change with its own power and specificity readings, as the within-stage null had.
- **Any default.** Nothing `motco analyze` or the harness runs by default changes.
- **A correlation setting.** In the scores it acts as a dispersion scale, which the cross covers in both directions; what it adds, a ratio that differs by component, stays open (design D3).
- **Dispersion by stage, by feature, or by group and stage together,** and `motco simulate` as a user of the setting.
- **A dispersion estimate or test in `motco analyze`.** The documentation says what a case study would compute.
- **Matching a cohort.** No design reproduces SEA-AD's or ADNI's cell table, covariates or blocks, and neither cohort's dispersion is read.
- **Real noise and feature count** (roadmap open items 2 and 15).

## Capabilities

### New Capabilities
- `group-dispersion-validation`: the judged null-only reading of the trajectory test under group-specific dispersion crossed with group balance. It covers the two profiles and their reproduction of the committed parent, the recorded within-group arm, the four size families and the verdict, the map by statistic, the recorded arms' labels and the dated report.

### Modified Capabilities
- `numpy-omics-generator`: a sampling call gains an opt-in noise scale that stretches the drawn noise about the exact cell means and draws nothing.
- `semisynthetic-trajectory-generator`: the generator gains the opt-in `group_dispersion_scale`, with its validation, its truth record, its invisibility while absent and the stage-span scales under it.
- `trajectory-power-study`: study configurations accept `generator.group_dispersion_scale`, hash-omitted while absent and kept by null cells.
- `latent-space-metric-compatibility`: the measurement-visibility documentation states what the tests assume about group dispersion, the verdict, the map and how a cohort is placed.

## Impact

- **Code:**
  - `src/motco/simulations/generator.py`: `noise_scale` in `generate_omics` and its sampling helper;
  - `src/motco/simulations/semisynthetic.py`: the field, its validation and normalization, the scale passed to group B's call, the truth record, and the within-cell multiplier of `null_pooled_scale` and `planted_pooled_scale`;
  - `src/motco/simulations/grid.py` and `src/motco/simulations/study/config.py`: `group_dispersion_scale` in the omit-when-`None` sets;
  - `src/motco/__init__.py` and `CHANGELOG.md`: 0.26.0.

  `motco.stats`, the evaluation harness, `motco.analysis` and the CLI are not touched.
- **Scripts:** `scripts/group_dispersion_check.py` and `scripts/group_dispersion_array.sbatch`. They import unchanged `univariate_baseline_check` (the unit and its arms, the rate and cell tables, the accounting), `latent_space_comparison` (seeds, settings binding) and `stage_span_paper_grade_decision` (`judge_rates`).
- **Committed parents bound by sha256:** `settings.json`, `rate_table.csv` and `cell_table.csv` under `results/univariate-baseline-2026-11-23/stage_span/` and under `results/stage-effect-scaling-2026-11-24/band/`.
- **Tests:** `tests/test_group_dispersion_scale.py` (the stretch about the cell means, the stream position, validation, truth, identity while absent and at 1, the stage-span scales, hash omission, null cells keeping the scale, study configs) and `tests/test_group_dispersion_check.py` (designs, unit counts and seeds per profile, the unit against the parent's, the within-group arm and its checks, resume and settings refusals, a refused parent, both outcomes of the reproduction and every decision path on synthetic tables). The committed identity fixtures must pass unchanged.
- **Dependencies:** none added.
- **Compute:** 32,000 units. The stage-effect scaling run measured 25 s per unit at n = 320 and 40 s at n = 800 on the workstation, and its arrays ran at about 1.4 times the unpacked projection. That puts the estimate before the smoke at 250–350 core-hours on `ing`, the within-group arm included. A local smoke replaces it, and launch waits on the user's approval with a cap.
- **Compatibility:** every committed configuration, fixture, record and report is unchanged. The setting is absent by default everywhere, and the parents' committed files are read, hashed and refused by name if changed.
