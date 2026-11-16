# Notes

## Exploratory probe (2026-11-16; not a result)

Run on the workstation from `main` at `bb6004c` (MOTCO 0.24.0) before the proposal was written, with BLAS pinned to one thread and 10 workers. The script was a scratch file outside the repository and is not committed. Its seeds are spent: no unit of this change's `dispersion` profile may reuse them (the profile derives its seeds from base seed 2038 through the driver's own derivation).

**Question.** The diagnosis (roadmap item 19) described unequal residual dispersion by stage, in opposite directions in the scores and in the predictive coordinates, and did not test whether that alone mis-sizes RRPP's unrestricted draws. Does it, and do within-stage draws hold?

**Construction.**
- The SEA-AD cell table: F 11 / 10 / 28, M 9 / 10 / 12 (80 samples, three stages of 20 / 20 / 40).
- `Y = c[stage] + ε`, `ε ~ N(0, v_stage · I_d)`, d ∈ {2, 5}, no group difference. Centroids (0, 0), (4, 0) and (2, 3.6) in the first two dimensions, a population stage share of variance of about 0.72 at d = 2.
- The design from `build_simulation_trajectory_design`; no covariates, no residual correction; net displacement.
- Both nulls on the same matrix at one seed and 199 draws: `RRPP` (unrestricted), and a loop equal to `predictive_size_excess_check.within_stage_indices` / `restricted_null` (within stage). Rejection at p < 0.05.
- 2,000 replicates per cell, seeded from `[20261116, pattern index, d, 400, replicate]`. The Monte Carlo SE at a rate of 0.05 is 0.005.

**Variance multipliers by stage** (stages of 20 / 20 / 40 samples):

| pattern | multipliers | source |
|---|---|---|
| `equal` | 1 / 1 / 1 | — |
| `largest_narrow` | 1.12 / 1.25 / 0.82 | rounded from the diagnosis's predictive medians (its report, section 9) |
| `largest_wide` | 0.80 / 0.92 / 1.14 | rounded from its score medians |
| `largest_narrow_x2` | 1.5 / 1.5 / 0.5 | chosen by hand |
| `largest_wide_x2` | 0.5 / 0.5 / 1.5 | chosen by hand |

**Rejection rates** (`delta` / `angle` / `shape`):

| pattern | d | unrestricted | within stage |
|---|---|---|---|
| `equal` | 2 | 0.051 / 0.045 / 0.046 | 0.049 / 0.041 / 0.043 |
| `largest_narrow` | 2 | 0.066 / 0.052 / 0.057 | 0.045 / 0.049 / 0.039 |
| `largest_wide` | 2 | 0.043 / 0.043 / 0.034 | 0.048 / 0.047 / 0.043 |
| `largest_narrow_x2` | 2 | 0.095 / 0.064 / 0.097 | 0.052 / 0.054 / 0.049 |
| `largest_wide_x2` | 2 | 0.014 / 0.034 / 0.019 | 0.038 / 0.045 / 0.050 |
| `equal` | 5 | 0.053 / 0.049 / 0.058 | 0.050 / 0.049 / 0.052 |
| `largest_narrow` | 5 | 0.070 / 0.050 / 0.074 | 0.052 / 0.051 / 0.050 |
| `largest_wide` | 5 | 0.040 / 0.037 / 0.042 | 0.049 / 0.046 / 0.048 |
| `largest_narrow_x2` | 5 | 0.092 / 0.076 / 0.096 | 0.040 / 0.051 / 0.041 |
| `largest_wide_x2` | 5 | 0.015 / 0.024 / 0.016 | 0.050 / 0.048 / 0.044 |

**What it suggests.**
- With equal dispersion both nulls read near 0.05.
- When the largest stage is the narrowest, the unrestricted draws exceed on `delta` and `shape`, and more as the inequality grows. This is the direction of the predictive coordinates' excess on the anchor.
- When the largest stage is the widest, the unrestricted draws are conservative. This is the scores' pattern, and it would account for production's anchor rates under the unrestricted draws (`delta` 0.042, `shape` 0.026) and for its plant-power gain under within-stage draws (0.51 → 0.63 at e = 0.5).
- The within-stage draws read 0.038–0.054 in every cell.
- At the measured inequality the unrestricted `delta` is 0.066–0.070, a fifth to a quarter of the anchor's excess over 0.05 (0.130–0.136). Unequal isotropic dispersion gives the direction and not the size.

**What it is not.** One cell table, isotropic Gaussian noise, no fit, no covariates, 199 draws, and hand-chosen doubled patterns. It motivated the `dispersion` profile (design D11), which redraws every cell on fresh seeds at 999 permutations under the predeclared rule.

## Pre-flight and smoke (2026-11-16, task 6.1)

Run on the workstation from the worktree's uncommitted state at `bb6004c` plus this change (MOTCO 0.25.0; Python 3.11.15, numpy 2.3.5, scipy 1.16.3, scikit-learn 1.8.0, pandas 2.3.3, the pilot's and the diagnosis's package versions), BLAS pinned to one thread. The export is `sead_sexdiff/results/motco_input/cohort_30`, read in place. Output went to `/tmp/wsp-smoke` and is not committed. A first pass had run on a Python 3.14 venv that `uv run` had created in the worktree; the venv was rebuilt on 3.11.15 and everything below is from the second pass, whose timings were within a few percent of the first.

**Pre-flight: passed.** On replicates 0 and 500 (the first fresh draw) the cross-validated `pls_cv` selected rank 5 at held-out AUROC 0.635556, the plasmode's fidelity values, and its statistics and p-values equalled the forced call's exactly. The CV took 67.9 s and 67.2 s.

**Smoke: every unit completed and consistent.**

| profile | units | workers | wall | failed | inconsistent |
|---|---|---|---|---|---|
| `plasmode` (replicates 0, 1, 500, 501) | 4 | 4 | 6 s | 0 | 0 |
| `weak_size` (replicates 0, 1 × anchor, translation) | 4 | 4 | 14 s | 0 | 0 |
| `dispersion` (20 per cell) | 200 | 7 | 15 s | 0 | 0 |
| `simulation` (`n200`, `seaad`; 2 per cell) | 44 | 7 | 92 s | 0 | 0 |

- **Checks.** On every fitted unit the pilot's five checks held, `same_statistics` held exactly, and on every plasmode unit `construction_identity` held: the package's within-stage p-values equal the diagnosis's loop for all four within-stage arms.
- **A known value.** Replicate 0's production `delta` p-value is 0.036 under the unrestricted draws, the diagnosis's pre-flight value.
- **Merge and decide** ran on all four profiles. The reproductions fail on every compared row (24 weak-size rows, 12 + 12 plasmode rows), because the counts are 2 to 4 against the committed 500, so `decide` returns `reproduction_failed`, as it must at this size. Its criteria carry no information at smoke size.

**Timing** (per unit, under the smoke's parallel load):

| quantity | simulation | weak_size | plasmode | dispersion |
|---|---|---|---|---|
| unit (median) | 12.6 s | 13.0 s | 4.3 s | 0.47 s |
| production CV | 10.4 s | 8.9 s | none (forced) | none |
| each arm | 0.25 s | 0.47–0.51 s | 0.39–0.49 s | 0.23–0.24 s |

**Core-hour projection.**

- **`ing`.** From the pilot's `ing` accounting (simulation unit median 20.6 s with six 0.4 s arms; weak-size 24.9 s with five 0.97 s arms) plus this driver's two extra arms: simulation 7,800 × 21.4 s ≈ 46 core-hours, weak size 1,000 × 27.8 s ≈ 7.7. The plasmode and dispersion profiles at the workstation-to-`ing` ratio of the pilot's forced arms (about 1.5): about 1.8 and 2.1. About 58 core-hours in all.
- **Workstation.** Simulation 7,800 × 12.6 s ≈ 27, weak size 3.6, plasmode 1.2, dispersion 1.3: about 33 core-hours, about 4.8 h of wall clock at 7 workers.
- **Proposed split (design D13).** `plasmode` and `dispersion` on the workstation (about 2.5 core-hours, under half an hour at 7 workers); `simulation` and `weak_size` on `ing` from a clean worktree at a pushed revision (about 54 core-hours; a 100-shard simulation array and a 20-shard weak-size array).
- **Cap.** 85 recorded core-hours on `ing` (about 1.5 times its projection), or 50 if everything runs on the workstation.

The full run waits on the user's approval of the cap and the host.
