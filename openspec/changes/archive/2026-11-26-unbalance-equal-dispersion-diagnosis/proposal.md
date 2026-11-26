# Proposal

## Why

The group-dispersion run (roadmap item 26, [report](../../../docs/reports/group-dispersion-2026-11-24.md)) read `unbalance_excess` on one reading of 24. At equal dispersion, `delta` at n = 80 with a first-group share of 0.8 rejected 44 times of 500 against a ceiling of 40. That design leaves 5 samples per stage in the smaller group. By the predeclared order the verdict stopped there, and the three dispersion families, which fail by wide margins, are reported and not interpreted until unbalance at equal dispersion is diagnosed (roadmap open item 16).

- **The committed readings point at the smallest cell.** `delta` at equal dispersion reads 44 of 500 at 5 per stage (n = 80, share 0.8), 27 and 34 at 13 per stage (n = 200, share 0.8), 30 at 13 per stage with balanced groups (n = 80) and 19 at 33 per stage (n = 500, share 0.8).
- **They cannot separate small cells from unbalance.** The one failing design is both the most unbalanced and the smallest. No committed design has small cells with balanced groups, and none has a small group inside a large sample.
- **Chance is not excluded.** The reading is one of 24, on one seed stream.
- **Much waits on it.** With balanced groups `delta` exceeds from s = 1.25 at n = 200 (86 of 500), at 33 samples per stage. Those readings carry no interpretation today.

## What Changes

- **A validation driver.** `scripts/unbalance_diagnosis_check.py` (`run` | `merge` | `decide`, an sbatch array beside it). No package change. Every unit is the group-dispersion run's unit, imported unchanged, on an equal-dispersion null dataset: three blocks, three stages, δ = 2.0, the anchor cell, base seed 2043, 999 permutations.
- **The gate.** The committed `n80_r0.8` design is rerun at its name and seeds (500 units) and must reproduce the group-dispersion run's committed counts exactly.
- **Seven new designs at 1,000 replicates** (samples per stage, larger / smaller group):

  | design | n | share | cells | what it separates |
  |---|---|---|---|---|
  | `n80_r0.8_b` | 80 | 0.8 | 22 / 5 | a fresh replicate of the failing design |
  | `n200_r0.925` | 200 | 0.925 | 62 / 5 | a 5-per-stage group inside a large sample |
  | `n42_r0.5` | 42 | 0.5 | 7 / 7 | small cells with balanced groups |
  | `n80_r0.74` | 80 | 0.74 | 20 / 7 | 7 per stage, unbalanced, at n = 80 |
  | `n200_r0.895` | 200 | 0.895 | 60 / 7 | 7 per stage inside a large sample |
  | `n80_r0.65` | 80 | 0.65 | 18 / 9 | SEA-AD's smallest cell size |
  | `n80_r0.5_b` | 80 | 0.5 | 14 / 13 | a fresh balanced reference at n = 80 |

  7,500 units in all. A balanced design at 5 per stage (n = 30) cannot run under the production cross-validation, so 7 per stage is the smallest balanced cell read (design D2).
- **Recorded arms.** The parent unit's arms, the within-group draws, and the two existing residual corrections (`df`, `leverage`) on the same fit, statistics and seed. They decide nothing.
- **A predeclared decision.** `unbalance_diagnosis_decision.json` (rule version 1). `delta`, `angle` and `shape` at p < 0.05, judged by `judge_rates` in three families:
  - `balanced_small` (`n42_r0.5`, m = 3);
  - `unbalanced_small` (the four designs with a 5- or 7-per-stage group, m = 12);
  - `reference` (`n80_r0.65` and `n80_r0.5_b`, m = 6).

  The verdict is the first that applies: `reproduction_failed` | `small_cells` (`balanced_small` fails) | `small_unbalanced_cells` (`unbalanced_small` fails) | `moderate_excess` (`reference` fails) | `not_replicated`.
- **A predeclared re-read of the group-dispersion families,** from that run's committed decision and with no new unit:
  - under `small_cells` or `small_unbalanced_cells`, the designs with a 5-per-stage group (n = 80 at shares 0.8 and 0.2) are set aside, the four families are judged again by the same rule, and the group-dispersion order gives the re-read verdict;
  - under `not_replicated`, the three readings of `n80_r0.8` are replaced by the fresh replicate's, and the families are judged in full;
  - under `moderate_excess` or `reproduction_failed`, nothing is re-read and `unbalance_excess` stands.
- **Documentation and records.** A dated report `docs/reports/unbalance-diagnosis-<date>.md`; aggregates under `results/unbalance-diagnosis-<date>/`; the measurement-visibility section on group dispersion; roadmap open item 16, the interpretability statement and the "Not yet established" entries; `AGENTS.md`.

**Non-goals.**
- **Any new dispersion unit.** Every design is at equal dispersion.
- **A remedy.** The corrections and the within-group draws are recorded only.
- **Any default or package change.** `motco analyze` and the harness run what they ran.
- **Balanced cells below 7 per stage,** which need a change to the cross-validation's candidate ranks.
- **Other stage effects, four stages, covariates, or a cohort's own table.** SEA-AD's table was read at its own design (roadmap items 7 (c) and 8).
- **Power.** Every unit is a null unit.

## Capabilities

### New Capabilities
- `unbalance-diagnosis-validation`: the judged null-only reading of the trajectory test at equal dispersion across the smaller group's cell size and group balance. It covers the reproduction of the committed failing design, the three families and the verdict, the recorded arms, the predeclared re-read of the group-dispersion families and the dated report.

### Modified Capabilities
- `latent-space-metric-compatibility`: the measurement-visibility documentation states from which cell size the tests lose size at equal dispersion, and the status of the group-dispersion readings after the re-read.

## Impact

- **Scripts:** `scripts/unbalance_diagnosis_check.py` and `scripts/unbalance_diagnosis_array.sbatch`. They import unchanged `group_dispersion_check` (the unit, the within-group arm, the merge helpers, the family judgement), `univariate_baseline_check`, `latent_space_comparison` (seeds, settings binding) and `stage_span_paper_grade_decision` (`judge_rates`).
- **Committed parents bound by sha256:** `columns/settings.json`, `columns/rate_table.csv`, `columns/cell_table.csv`, `columns/within_group_table.csv` and `group_dispersion_decision.json` under `results/group-dispersion-2026-11-24/`.
- **Tests:** `tests/test_unbalance_diagnosis_check.py`.
- **Package code:** none. MOTCO stays at 0.26.0.
- **Branch:** the change stacks on `feat/group-specific-dispersion`, which holds the parent driver and its committed results and is not merged yet.
- **Compute:** 7,500 units at n ≤ 200. The group-dispersion run measured a median of 14.8 s per unit at n = 200 on `ing`, which puts the estimate near 40 core-hours. A local smoke replaces it, and launch waits on the user's approval with a cap.
- **Compatibility:** every committed configuration, fixture, record and report is unchanged.
