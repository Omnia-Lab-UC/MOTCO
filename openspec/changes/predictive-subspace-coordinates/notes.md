# Notes

## Exploration probe (2026-11-12)

An exploratory probe run while roadmap open item 1 was being explored, before this proposal. **Nothing below is a result.** It was run with throwaway scripts on the workstation; its records were kept in a scratch directory and are not committed. The pilot re-measures every reading on fresh datasets.

**Setup.** MOTCO 0.23.0 (package code as at `ca41a36`), Python 3.11.15, numpy 2.3.5, scipy 1.16.3, scikit-learn 1.8.0, BLAS pinned to one thread.
- **Units.** `latent_space_comparison`'s `design_generators`, `generator_seed`, `rrpp_seed` and `unit_dataset`, imported unchanged, at base seed 4242 and the comparison's null scale: `n200`, `seaad` and `seaad_weak` (δ = 0.11); n = 200 at δ = 1.0 and 0.25; and the SEA-AD table under the default feature-set orientation surgery at e = 1.0. 199 permutations, rejection at p < 0.05. Base seed 4242 is spent.
- **Arms**, all on one dataset and one RRPP seed per unit:
  - production (`pls_cv`, floor 2, net displacement) and a forced rank of 2;
  - four coordinate systems on the predictive subspace of the production fit: the three this change proposes, and the oblique projection `T·Qᵀ(Q·D·Qᵀ)⁺Q·D` with `D = diag(TᵀT)`, whose dropped scores are uncorrelated with the kept ones;
  - four stricter rank rules, each read as a forced rank chosen on a replica of the production CV's inner AUROC curves: a one-standard-error rule and an absolute tolerance of 0.01 inside each outer fold, and the same two on the curve pooled over all folds.
- **Sanity, over the 1,080 units with a production CV:** the forced fit at the production rank gave production's p-values on every unit; the CV replica's argmax rule gave the production rank on every unit; the Y loadings had rank 2 on every unit.

**A. Target rates against production** (60 units per cell, 20 where noted; `predictive` with its discordant counts b/c against production, b = only production rejects).

| design | cell, statistic | production | forced rank 2 | tolerance rule | `predictive` |
|---|---|---|---|---|---|
| SEA-AD table, δ = 2.0 | magnitude e = 0.05, `delta` | 0.23 | 0.20 | 0.17 | 0.53 (0/18) |
| | 22.5° rotation, `angle` | 0.45 | 0.93 | 0.78 | 1.00 (0/33) |
| | bend e = 0.25, `shape` | 0.68 | 0.60 | 0.55 | 0.88 (0/12) |
| | 45° rotation at φ = 60, `angle` | 0.42 | 0.97 | 0.78 | 1.00 (0/35) |
| n = 200, δ = 2.0 | magnitude e = 0.05, `delta` | 0.63 | 0.65 | 0.65 | 0.80 (0/10) |
| | 22.5° rotation, `angle` | 0.82 | 1.00 | 1.00 | 1.00 (0/11) |
| n = 200, δ = 1.0 (20) | magnitude e = 0.25, `delta` | 0.85 | 0.95 | 0.85 | 1.00 |
| n = 200, δ = 0.25 (20) | magnitude e = 0.25, `delta` | 0.25 | 0.20 | 0.30 | 0.95 |
| | 45° rotation, `angle` | 0.55 | 0.15 | 0.55 | 1.00 |
| SEA-AD table, δ = 0.11 | 45° rotation, `angle` | 0.03 | 0.02 | 0.03 | 0.10 (1/5) |
| | magnitude e = 0.05, `delta` | 0.05 | 0.05 | 0.03 | 0.07 (0/1) |

Production's δ = 2.0 rates agree with the production selection's `pls_cv` on other datasets (0.23 / 0.46 / 0.70 on the SEA-AD table; 0.74 / 0.84 at n = 200).

**B. The rank rules.**
- **At δ = 2.0 the AUROC is at its ceiling.** On the SEA-AD anchor the mean inner AUROC is 0.9960 at rank 2 and 1.0000 at rank 3, and the best candidate's fold SE is exactly 0 in every outer fold. The inner one-standard-error rule therefore selected production's rank on every unit of the anchor and target cells at δ = 2.0 and at δ = 1.0 (480 of 480; 30 / 30 at ranks 2 / 3 on the SEA-AD anchor). The translation cell, whose AUROC is 0.93, is the exception (28 of 60 equal). On the SEA-AD anchor the tolerance rule moved the split to 49 / 11, and the pooled one-standard-error rule selected rank 3 in 58 of 60.
- **At δ = 0.25** production selected ranks 8–19 on the anchor and the inner rules 7–9, with no change in any rate (magnitude `delta` 0.30 against 0.25; rotation `angle` 0.50–0.55 against 0.55).
- **At the weak design** production selected rank 2 in 111 of 240 anchor units and ranks 3–11 in the rest; the inner one-standard-error rule selected rank 2 in 165. No rate moved.

**C. Fit rank with the measured dimension fixed** (40 units per cell, forced fit ranks, no CV). Entries are "measured in all A dimensions / in the predictive subspace"; at rank 2 the two are the same measurement.

SEA-AD table, δ = 2.0:

| fit rank | in-sample R² of the stage labels (anchor) | magnitude e = 0.05, `delta` | 22.5° rotation, `angle` | bend e = 0.25, `shape` | rotation's off-target `shape` | `pooled_stage` capture of the predictive subspace (rotation cell) |
|---|---|---|---|---|---|---|
| 2 | 0.946 | 0.25 | 0.90 | 0.60 | 0.05 | 0.905 |
| 3 | 0.993 | 0.12 / 0.85 | 0.20 / 1.00 | 0.55 / 1.00 | 0.10 / 0.03 | 0.935 |
| 4 | 0.996 | 0.12 / 0.85 | 0.33 / 1.00 | 0.47 / 1.00 | 0.07 / 0.85 | 0.878 |
| 6 | 0.999 | 0.12 / 0.95 | 0.70 / 1.00 | 0.38 / 1.00 | 0.05 / 0.53 | 0.729 |
| 10 | 1.000 | 0.12 / 0.45 | 0.90 / 0.72 | 0.28 / 0.70 | 0.07 / 0.05 | 0.510 |
| 16 | 1.000 | 0.12 / 0.00 | 0.93 / 0.00 | 0.28 / 0.05 | 0.10 / 0.00 | 0.433 |

n = 200, δ = 0.25:

| fit rank | in-sample R² of the stage labels (anchor) | magnitude e = 0.25, `delta` | 45° rotation, `angle` | rotation's off-target `shape` |
|---|---|---|---|---|
| 2 | 0.358 | 0.17 | 0.20 | 0.65 |
| 3 | 0.719 | 0.35 / 0.57 | 0.42 / 0.97 | 0.50 / 0.68 |
| 4 | 0.820 | 0.30 / 0.72 | 0.45 / 1.00 | 0.55 / 0.75 |
| 6 | 0.907 | 0.38 / 0.97 | 0.40 / 1.00 | 0.57 / 0.82 |
| 10 | 0.966 | 0.30 / 0.97 | 0.65 / 1.00 | 0.55 / 0.68 |
| 16 | 0.989 | 0.28 / 0.93 | 0.85 / 1.00 | 0.57 / 0.42 |

The anchors of this ladder read 0.00–0.07 for the predictive subspace at fit ranks 3–16 on both designs, and 0.00–0.10 at the weak design.

**D. Controls** (`delta` / `angle` / `shape`).

| controls | production | forced rank 2 | `predictive` | `predictive_orthonormal` | `predictive_fitted` |
|---|---|---|---|---|---|
| SEA-AD table at δ = 2.0, anchor and translation pooled (120 readings) | 0.033 / 0.050 / 0.017 | 0.075 / 0.017 / 0.092 | 0.017 / 0.042 / 0.008 | 0.017 / 0.042 / 0.008 | 0.017 / 0.042 / 0.025 |
| weak anchor (300 units) | 0.017 / 0.043 / 0.050 | 0.017 / 0.047 / 0.047 | 0.033 / 0.057 / 0.077 | 0.020 / 0.060 / 0.063 | 0.047 / 0.063 / 0.040 |

The anchors at n = 200 and δ = 1.0 and 0.25 have 20 units each and read 0.00–0.15 in every arm.

**E. Off-target rates under a pure rotation** (`delta` / `shape`).

| design, rotation (units) | production | `predictive` | oblique | `predictive_fitted` | `predictive_orthonormal` |
|---|---|---|---|---|---|
| SEA-AD table, δ = 2.0, 22.5° (60) | 0.02 / 0.07 | 0.00 / 0.02 | 0.00 / 0.07 | 0.00 / 0.03 | 0.00 / 0.00 |
| n = 200, δ = 2.0, 22.5° (60) | 0.02 / 0.08 | 0.00 / 0.07 | 0.00 / 0.12 | 0.00 / 0.07 | 0.00 / 0.00 |
| n = 200, δ = 1.0, 45° (20) | 0.05 / 0.50 | 0.05 / 0.50 | 0.10 / 0.80 | 0.00 / 0.00 | 0.00 / 0.00 |
| n = 200, δ = 0.25, 45° (20) | 0.20 / 0.60 | 0.30 / 0.80 | 0.50 / 1.00 | 0.00 / 0.00 | 0.00 / 0.25 |

Target `angle` power is 1.00 in all four predictive coordinate systems in every row. The probe's plants are at the null scale, so the band's off-target rates, production's included, cannot be attributed to the measurement.

**F. Observed against population** (median observed `angle`, then the population's through the same map, in degrees).

| design, rotation | production | `predictive` | `predictive_fitted` | `predictive_orthonormal` |
|---|---|---|---|---|
| SEA-AD table, δ = 2.0, 22.5° | 22.4 / 22.2 | 20.2 / 22.3 | 20.3 / 22.1 | 20.7 / 22.5 |
| n = 200, δ = 2.0, 22.5° | 22.4 / 22.6 | 21.9 / 22.5 | 22.2 / 22.8 | 21.9 / 22.5 |
| n = 200, δ = 1.0, 45° | 44.2 / 45.7 | 43.1 / 45.1 | 42.2 / 44.6 | 41.5 / 45.0 |
| n = 200, δ = 0.25, 45° | 32.9 / 41.6 | 16.3 / 40.6 | 19.0 / 46.0 | 16.6 / 43.8 |

**G. Capture and the surgery cell.**
- **`pooled_stage` capture** (median): production / forced rank 2 / predictive subspace read 0.93 / 0.90 / 0.93 on the SEA-AD table at δ = 2.0, 0.67–0.77 / 0.30–0.35 / 0.61–0.70 at n = 200 and δ = 0.25, and 0.08–0.10 / 0.07 / 0.08–0.10 at the weak design. The largest principal angle between the predictive subspace and the rank-2 span has a median of 12–18° at δ = 2.0 (units where the CV chose 3), 88° at δ = 0.25 and 84° at the weak design.
- **The Phase 6 feature-set orientation surgery at e = 1.0** (60 units; production rank 3 in 34, 4–7 in 26): `delta` / `angle` / `shape` read 0.32 / 0.58 / 0.38 for production (`group_difference` capture 0.081), 0.57 / 0.35 / 0.53 for the forced rank 2 (0.015) and 0.53 / 0.37 / 0.60 for `predictive` (0.004).

**What the probe does not show.** Rates carry an SE of up to 0.065 at 60 units and 0.11 at 20. Nothing was predeclared. The band's cells are at the null scale. The plasmode export was not read. The fit-rank ladder pools units the CV would have placed at different ranks.

## Smoke (2026-11-12, task 6.1)

A local smoke of the three profiles at the worktree's uncommitted state (MOTCO 0.24.0; Python 3.11.15, numpy 2.3.5,
scipy 1.16.3, scikit-learn 1.8.0), BLAS pinned to one thread, on the 16-thread workstation with 11 processes at once
(simulation on 6 workers, plasmode on 4, `weak_size` serial). Output in `/tmp/psc_smoke/run` (not committed).

- **Units.** `simulation` at 2 replicates per cell (156 units: every design, the planted scale in the band
  included), `weak_size` at 2 replicates per cell (4 units, 999 permutations), `plasmode` at 2 anchor and 2 plant
  replicates per cell (8 units, 999 permutations). Wall clock 863 s.
- **Outcome.** 168 of 168 units completed: none failed and none inconsistent. **No planted fixed point failed in any
  design**, so the band needs no amendment. Largest consistency deviations: `fitted_labels` 5.4e-14, every other check
  0. `identity_at_low_rank` applied on 34 of 156 simulation units and 2 of 4 `weak_size` units, and never on the
  plasmode (rank 5–7). The plasmode fidelity values (rank 5, AUROC 0.635556) held on all 4 anchor units.
- **Merge, reproduction and decide** ran. The gate reads `reproduction_failed` only because the smoke's counts are 2
  against the committed 500 (`weak_size`, plasmode anchor) and 100 (plasmode plant); every differing row is a count
  difference. Every arm read `no_gain`, which carries no information at 2 replicates.
- **CV rank by design** (production arm, min–max over the smoke's units, median): `n200` 2–3 (2); `seaad` 2–3 (3);
  `n200_delta1` 3 (3); `n200_delta0.5` 3–5 (3); `n200_delta0.25` 5–15 (9.5); `seaad_delta0.5` 3–8 (4); `seaad_weak`
  2–5 (3); `seaad_surgery` 3–5 (4). The CV rank inflates in the band, as the design expects.

**Timing** (per unit, under the smoke's load):

| profile | unit median (mean) | where the time goes |
|---|---|---|
| `simulation` | 18.1 s (19.0 s); `n200` designs 19–23 s, SEA-AD designs 13–15 s | the one production CV (12–20 s); every arm, the two forced fits included, 0.3–0.4 s |
| `weak_size` | 24.0 s (24.1 s) | the production CV (19.5 s); each arm about 1 s at 999 permutations |
| `plasmode` | 410 s (421 s) | four cross-validated `analyze_trajectories` calls (`pls_cv` 106 s; `predictive` 96 s, `predictive_orthonormal` 92 s, `predictive_fitted` 100 s); `forced_k_minus_1` 0.4 s |

**Core-hour projection.** Local core-hours: `simulation` 7,800 × 19.0 s = 41.2, `weak_size` 1,000 × 24.1 s = 6.7 and
`plasmode` 800 × 421 s = 93.6, 141.5 in all. That is above the proposal's 100–140 estimate, because each plasmode arm
repeats the label-free CV. The production selection's plasmode `pls_cv` arm took 134 s per unit on `ing` against
106 s here, so on `ing` the run would take about 1.26 × 141.5 ≈ 178 core-hours.

**Amendment (approved and applied on 2026-11-12, before any array unit).** On the plasmode the three predictive arms repeat a CV whose
result is known from `pls_cv`: the fit depends on the features and the stage labels only. An `analyze_trajectories`
call at `forced_components` = the unit's `pls_cv` rank fits the same final model with the same `fit_plsda_model` call,
so it gives the same latent matrix. On the smoke's first two plasmode units (anchor rank 5, `span_e0.25` rank 7) the
forced call reproduced each predictive arm's statistics, p-values, measured dimension and `angle` null quantile exactly,
in 0.4 s instead of about 100 s.
- **What it changes.** D5's plasmode clause would read "`pls_cv` with its CV, then each predictive arm at
  `forced_components` = the `pls_cv` rank". `same_fit` would compare the predictive arms' fitted model with
  `pls_cv`'s (rank and a hash of `x_rotations_`), not their AUROC, because a forced analysis records none. The arms'
  records would carry `pls_cv`'s AUROC.
- **What it saves.** The plasmode unit falls to about 110 s: 800 × 110 s ≈ 24 core-hours locally. The whole run falls
  to about 72 local core-hours, about 91 on `ing`.
- **What it does not change.** The simulation and `weak_size` profiles already fit once per unit, and no rule, unit,
  seed or reproduction changes.

**Cap.** 220 core-hours on `ing` as designed, or 120 core-hours with the amendment, about 1.25 × the projection
either way. The user approved the amendment, which sets the cap at 120 core-hours. The other lever the design names, fewer plasmode anchor replicates, would break the exact
reproduction of the plasmode's committed counts over 500 replicates, so it is not proposed.

**Applied.** In the driver, `plasmode_arm_params(..., production_rank=)` runs each predictive arm at the `pls_cv`
rank. The arm's record carries `fit: "forced_at_pls_cv_rank"` and `pls_cv`'s AUROC, and `same_fit` compares every
predictive arm's rotations and Y loadings with production's. The plasmode settings record `predictive_arm_fit`.
`tests/test_predictive_subspace_check.py` checks two things on the stand-in export: the forced arms equal a
cross-validated analysis under the same setting, and a plasmode unit runs one CV.
