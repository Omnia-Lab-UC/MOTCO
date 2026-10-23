# Rank-floor check — 2026-10-23

> Paper-grade size and power check of the PLS latent-space rank floor (`min_components=2`, MOTCO 0.10.0) against the
> unfloored 0.9.0 rule at the Phase 6 SEA-AD design (`openspec/changes/archive/2026-10-23-minimum-latent-rank/`, design D7, spec
> `rank-floor-validation`). It is the gate behind the default flip: `motco analyze` and the simulation harness select
> the PLS rank from candidates at or above 2 only once this check names the floor-2 arm admissible. Every number below
> is read from `results/rank-floor-2026-10-23/` (`delta0.05/` and `delta0.1/`: `rank_floor_decision.json`,
> `rank_floor_rates.csv`, `rank_floor_replicates.csv`, `rank_floor_pvalues.csv`, `merged.jsonl`, the last three kept on ing as
> `PROVENANCE.txt` records; `calibration/`;
> `smoke-delta0.05/`, `smoke-delta0.1/`). Rates are rejection rates at p ≤ 0.05 over the available replicates of 500,
> at 999 permutations.

**Run:** two SLURM arrays of 100 shards on `ing` (partition `512x1024`): job 925691 for δ = 0.05 and job 925791 for
δ = 0.1, 7,000 units each, 100 / 100 tasks completed, 0 arm failures, 53.0 + 45.8 recorded
core-hours (29 and 25 s per unit, two analyses at 999 permutations). The code is local `main` at `0a15d44` plus this
change's working tree, rsynced over the cluster checkout and verified by tree hash; file hashes are in each
`PROVENANCE.txt`.

**Decision:** `adopt_floor` at δ = 0.05 and `adopt_floor` at δ = 0.1, hence **`adopt_floor`** under the predeclared
conjunction. The floor-2 arm holds the bound 0.0695 for every statistic at both controls in both runs.

**In one paragraph.** At the Phase 6 SEA-AD design with the stage signal weakened until the unfloored rank rule
selects rank 1 in 58 % (δ = 0.05) and 42 % (δ = 0.1) of the null replicates (held-out CV AUROC
0.51 and 0.57), each of 7,000 units per strength was analysed twice on the same data and permutation draws,
with the 0.9.0 rule and with the rank floor of 2. In the floored arm the zero-effect anchor rejects `delta` 0.042 / 0.044,
`angle` 0.054 / 0.060 and `shape` 0.036 / 0.038 (δ = 0.05 / 0.1), and the translation
control 0.036 / 0.036, 0.046 / 0.042 and 0.026 / 0.030, all at or below
0.0695 at R = 500. The floor lifted 4,257 and 3,228 rank-1 selections, mostly to rank 2; in the unfloored arm every rank-1 unit measured
`angle` as exactly 0° or 180° with a two-point null that never rejected. Power is at the null in both arms at
these strengths, so the check establishes size, not a gain: the floor defines `angle` and `shape` where the 0.9.0
rule left them degenerate, and the record's spectrum and null summary say how resolvable they are.

## 1. Construction and calibration

- **Why a new construction.** The 0.9.0 latent-adjustment power comparison recorded the CV-selected rank per unit at
  the Phase 6 construction: no unit in any arm selected rank 1 ([roadmap](../roadmap.md), Next changes item 4), while
  both SEA-AD cohorts did (held-out AUROC 0.63 / 0.66). The floor changes the measurement only where the unfloored
  rule selects rank 1, so it is validated where that happens.
- **Data.** `generate_semisynthetic_trajectory` at the Phase 6 small-n design: three stages,
  `group_stage_sizes = [[11, 10, 28], [9, 10, 12]]` (n = 80), methylation (M-values) and expression measured, 498
  features, joint magnitude, `p_dmp = 0.1`, `baseline_continuity = 0`, with the **stage signal weakened** through the
  per-omic effect sizes `delta_methyl = delta_expr = delta_protein = δ`, at the Phase 6 value 2.0 for every committed
  study. Two strengths were run, δ = 0.05 and δ = 0.1 (below). `p_dmp` stays at 0.1 so every pool-limited surgery
  keeps its Phase 6 headroom (translation 98 sites available, orientation and shape 242; every cell fits).
- **Cells.** The zero-effect anchor; the `translation` negative control at e = 1.0; `magnitude` at e = 0.02, 0.05,
  0.1, 0.25 (the Phase 6 magnitude profile below its saturation); `orientation` and `shape` at e = 0.25, 0.5, 0.75,
  1.0 (the Phase 6 study profile). 14 cells × 500 replicates = 7,000 units per strength. The generator and RRPP seeds
  depend on the replicate index only (base seed 2026), so every cell and both arms of replicate `r` share the
  baseline draw and the permutation draws.
- **Arms.** Each unit is analysed with `analyze_trajectories` twice, with the harness CV defaults
  (`cv2_splits = 4`, `cv1_splits = 3`, `n_repeats = 5`, `max_components = 20`, `random_state = 1203`), 999
  permutations, the same RRPP seed, no covariates and the net-displacement estimator:
  - `floor_1`: `min_components = 1`, the 0.9.0 rule (measured, not judged);
  - `floor_2`: `min_components = 2`, the candidate default.
- **Rule** (predeclared in design D7, pure function `rank_floor_decision`). For the anchor and the translation
  control, the `floor_2` arm's rate over its R available replicates must be at most α + k·√(α(1 − α)/R) with
  α = 0.05 and k = 2 (0.0695 at R = 500) for `delta`, `angle` and `shape`. Every statistic holding at both controls
  gives `adopt_floor`, otherwise `diagnose_floor_excess`. With two strengths, `adopt_floor` requires both runs to hold
  (settled with the user before submission). The `floor_1` arm's rates are reported and never enter the decision;
  power rows are descriptive.

### Calibration

The calibration criterion is that the unfloored rule selects rank 1 in at least half of the null (anchor) replicates
on a recorded smoke run before any paper-grade run (`scripts/rank_floor_check.py calibrate`; every iteration is in
`results/rank-floor-2026-10-23/calibration/calibration.json`). Weakening δ alone first *raises* the selected rank:
at δ = 1.0 and 0.5 the held-out AUROC stays at 1.0 and the CV picks 3–6 components, at δ = 0.25 it is 0.92 with
4–11 components, because the stage signal thins across more directions before it fades. It reaches rank 1 only once
the signal is near chance.

| δ | R | rank-1 fraction | rank distribution | mean CV AUROC |
|---|---|---|---|---|
| 2.0 (Phase 6) | 20 | 0.00 | 2: 8, 3: 12 | 1.000 |
| 1.0 | 20 | 0.00 | 3: 20 | 1.000 |
| 0.5 | 20 | 0.00 | 3: 1, 4: 8, 5: 9, 6: 2 | 1.000 |
| 0.25 | 20 | 0.00 | 4–11 | 0.920 |
| 0.15 | 50 | 0.08 | 1: 4, 2: 12, 3: 8, 4–12: 26 | 0.702 |
| 0.1 | 20 / 50 | 0.55 / 0.44 | 1: 22, 2: 13, 3: 5, 4–10: 10 (R = 50) | 0.579 / 0.583 |
| 0.05 | 50 | 0.58 | 1: 29, 2: 14, 3: 4, 4: 3 | 0.525 |

Both strengths meet the criterion on their smoke replicates (δ = 0.05: 0.70; δ = 0.1: 0.55), and they bracket the
SEA-AD regime from the two sides: δ = 0.1 has the cohorts' held-out AUROC (0.58 against 0.63–0.66) with a rank-1
fraction near one half, δ = 0.05 is robustly in the rank-1 regime with almost no stage signal left (AUROC 0.52).

### Smoke profiles

R = 20 at 99 permutations, 14 cells × 2 arms = 560 analyses per strength (`smoke-delta0.05/`, `smoke-delta0.1/`). No
arm failed; the `floor_2` arm selected rank ≥ 2 in every unit while `floor_1` measured `angle` in one dimension in
183 / 280 (δ = 0.05) and 147 / 280 (δ = 0.1) units. Each smoke decision reads `diagnose_floor_excess` on one
statistic at 3 of 20 rejections (0.150 against the R = 20 bound 0.147): anchor `angle` at δ = 0.05, translation
`delta` at δ = 0.1. At R = 20 one rejection moves a rate by 0.05, so the smoke decisions are calibration and cost
evidence, not verdicts. Cost: 26–53 s per unit at 99 permutations locally, 58 s at 999 permutations.

## 2. Rates against the bound

Rate (Monte Carlo SE) with the floor-2 verdict against α + 2·√(α(1 − α)/500) = 0.0695; the floor-1 arm is descriptive.

| control | statistic | floor-2 δ = 0.05 | floor-2 δ = 0.1 | floor-1 δ = 0.05 | floor-1 δ = 0.1 |
|---|---|---|---|---|---|
| anchor | `delta` | 0.042 (0.0090) holds | 0.044 (0.0092) holds | 0.044 (0.0092) | 0.042 (0.0090) |
| anchor | `angle` | 0.054 (0.0101) holds | 0.060 (0.0106) holds | 0.030 (0.0076) | 0.038 (0.0086) |
| anchor | `shape` | 0.036 (0.0083) holds | 0.038 (0.0086) holds | 0.048 (0.0096) | 0.046 (0.0094) |
| translation | `delta` | 0.036 (0.0083) holds | 0.036 (0.0083) holds | 0.038 (0.0086) | 0.036 (0.0083) |
| translation | `angle` | 0.046 (0.0094) holds | 0.042 (0.0090) holds | 0.016 (0.0056) | 0.026 (0.0071) |
| translation | `shape` | 0.026 (0.0071) holds | 0.030 (0.0076) holds | 0.034 (0.0081) | 0.026 (0.0071) |

At p < 0.05 instead of p ≤ 0.05 the floor-2 anchor rates are `delta` 0.042 / 0.044, `angle` 0.052 / 0.058, `shape`
0.036 / 0.038 and the translation rates 0.036 / 0.036, 0.046 / 0.042, 0.024 / 0.030; the decision is unchanged.

**Stability.** Floor-2 rates over five consecutive blocks of 100 replicates:

- δ = 0.05, anchor: `delta` 0.06, 0.03, 0.05, 0.05, 0.02; `angle` 0.05, 0.04, 0.07, 0.05, 0.06; `shape` 0.04, 0.01, 0.07, 0.04, 0.02.
  Translation: `delta` 0.03, 0.03, 0.02, 0.05, 0.05; `angle` 0.02, 0.04, 0.06, 0.04, 0.07; `shape` 0.04, 0.01, 0.03, 0.02, 0.03.
- δ = 0.1, anchor: `delta` 0.07, 0.02, 0.05, 0.04, 0.04; `angle` 0.04, 0.05, 0.07, 0.05, 0.09; `shape` 0.06, 0.01, 0.07, 0.02, 0.03.
  Translation: `delta` 0.03, 0.03, 0.02, 0.05, 0.05; `angle` 0.02, 0.03, 0.06, 0.03, 0.07; `shape` 0.05, 0.01, 0.05, 0.01, 0.03.

No block exceeds 0.09, and no statistic is near the bound in both runs at once (`angle` at the anchor is the closest,
0.054 and 0.060). Every unit was available in both arms: no CV training fold failed, which the 0.9.0 checks saw only
with categorical covariates.

## 3. What the floor changes

- **Calibration realized.** The unfloored anchor arm selected rank 1 in 289 / 500 replicates at δ = 0.05
  (0.578; criterion met) and 210 / 500 at δ = 0.1 (0.420; the smoke replicates' 0.55 was the high side of the
  R = 50 estimate 0.44, so the paper-grade fraction sits just under one half, as reported before submission). The anchor
  rank distributions are 1: 289, 2: 120, 3: 59, 4: 18, 5: 10, ≥6: 4 (floor 1) against 2: 346, 3: 91, 4: 27, 5: 23, 6: 7, ≥7: 6 (floor 2) at δ = 0.05, and
  1: 210, 2: 131, 3: 66, 4: 27, 5: 25, ≥6: 41 against 2: 275, 3: 104, 4: 37, 5: 39, 6: 22, ≥7: 23 at δ = 0.1.
- **The lift.** Over all 7,000 units the unfloored rule selected rank 1 in 4,257 (δ = 0.05) and 3,228 (δ = 0.1); the floor
  moved those to rank 2 in 3,289 / 2,187, rank 3 in 558 / 499, rank 4 in 184 / 180 and higher in the rest. The two arms
  select the same rank in 2,635 / 3,633 units; in 108 / 139 units the unfloored modal rank was already ≥ 2 but
  moved, because outer folds that had chosen rank 1 chose another candidate under the floor. The floor-2 arm
  selected rank 2 in 4,795 / 3,909 units and rank ≥ 3 in the remainder; it never fell back to a lower rank.
- **`angle` at rank 1.** In every unfloored anchor unit that selected rank 1 the observed `angle` is exactly 0°
  (246 and 182 units) or 180° (43 and 28), its null's 95 % quantile is 180° (a two-point distribution) and
  its smallest p-value is 0.068 / 0.054: the rank-1 `angle` test never rejects. The unfloored arm's anchor `angle`
  rates (0.030 / 0.038) therefore come entirely from its rank ≥ 2 units, where they are 15 / 211 = 0.071 and
  19 / 290 = 0.066. Under the 0.9.0 rule a real-data `angle` p-value of 1.0 at rank 1 says nothing about
  orientation; this is the degeneracy the floor removes.
- **Null width.** The median over anchor replicates of the `angle` null's 95 % quantile is 84.0° / 61.2° in the
  unfloored arm (180° in its rank-1 units, 55.9° / 52.4° in its others) against 61.2° / 57.0° with the floor; the
  null's median SD is 34.0° / 21.0° against 17.8° / 16.3°. In a two-dimensional space with no stage signal the
  `angle` null is wide (two near-random directions), which is what a resolvable-orientation reader should see in
  `null_summary.angle`.
- **The signal.** Mean held-out CV AUROC at the anchor is 0.513 / 0.571 (floor 1) and 0.516 / 0.580 (floor 2):
  the floor does not change what the CV can find, only the smallest rank it may report.

## 4. Power in both arms

Rejection rates at p ≤ 0.05 over 500 replicates, floor 1 / floor 2, with the `angle` difference (floor 2 − floor 1)
and its pooled Monte Carlo SE; the last column is the share of units whose unfloored rule selected rank 1. Descriptive.

δ = 0.05:

| mode | e | `delta` 1 / 2 | `angle` 1 / 2 (diff, pooled SE) | `shape` 1 / 2 | rank-1 share, floor 1 |
|---|---|---|---|---|---|
| magnitude | 0.02 | 0.044 / 0.042 | 0.030 / 0.056 (+0.026, 0.013) | 0.050 / 0.038 | 0.58 |
| magnitude | 0.05 | 0.044 / 0.042 | 0.030 / 0.054 (+0.024, 0.013) | 0.046 / 0.036 | 0.57 |
| magnitude | 0.1 | 0.046 / 0.038 | 0.030 / 0.054 (+0.024, 0.013) | 0.044 / 0.036 | 0.58 |
| magnitude | 0.25 | 0.046 / 0.044 | 0.032 / 0.058 (+0.026, 0.013) | 0.044 / 0.040 | 0.58 |
| orientation | 0.25 | 0.030 / 0.022 | 0.020 / 0.032 (+0.012, 0.010) | 0.038 / 0.046 | 0.65 |
| orientation | 0.5 | 0.034 / 0.032 | 0.024 / 0.040 (+0.016, 0.011) | 0.024 / 0.024 | 0.60 |
| orientation | 0.75 | 0.044 / 0.038 | 0.010 / 0.042 (+0.032, 0.010) | 0.038 / 0.040 | 0.64 |
| orientation | 1 | 0.030 / 0.028 | 0.020 / 0.036 (+0.016, 0.010) | 0.032 / 0.040 | 0.64 |
| shape | 0.25 | 0.046 / 0.036 | 0.024 / 0.060 (+0.036, 0.013) | 0.044 / 0.054 | 0.62 |
| shape | 0.5 | 0.054 / 0.048 | 0.030 / 0.052 (+0.022, 0.013) | 0.048 / 0.046 | 0.61 |
| shape | 0.75 | 0.036 / 0.038 | 0.028 / 0.056 (+0.028, 0.013) | 0.048 / 0.036 | 0.62 |
| shape | 1 | 0.046 / 0.040 | 0.026 / 0.060 (+0.034, 0.013) | 0.044 / 0.032 | 0.64 |

δ = 0.1:

| mode | e | `delta` 1 / 2 | `angle` 1 / 2 (diff, pooled SE) | `shape` 1 / 2 | rank-1 share, floor 1 |
|---|---|---|---|---|---|
| magnitude | 0.02 | 0.046 / 0.048 | 0.038 / 0.060 (+0.022, 0.014) | 0.046 / 0.040 | 0.42 |
| magnitude | 0.05 | 0.046 / 0.046 | 0.038 / 0.060 (+0.022, 0.014) | 0.044 / 0.038 | 0.41 |
| magnitude | 0.1 | 0.044 / 0.048 | 0.034 / 0.060 (+0.026, 0.013) | 0.044 / 0.038 | 0.40 |
| magnitude | 0.25 | 0.046 / 0.050 | 0.042 / 0.062 (+0.020, 0.014) | 0.038 / 0.034 | 0.38 |
| orientation | 0.25 | 0.032 / 0.026 | 0.024 / 0.034 (+0.010, 0.011) | 0.052 / 0.048 | 0.49 |
| orientation | 0.5 | 0.030 / 0.030 | 0.026 / 0.048 (+0.022, 0.012) | 0.036 / 0.028 | 0.51 |
| orientation | 0.75 | 0.048 / 0.046 | 0.018 / 0.042 (+0.024, 0.011) | 0.038 / 0.040 | 0.54 |
| orientation | 1 | 0.030 / 0.032 | 0.028 / 0.042 (+0.014, 0.012) | 0.034 / 0.042 | 0.58 |
| shape | 0.25 | 0.044 / 0.040 | 0.036 / 0.056 (+0.020, 0.013) | 0.058 / 0.062 | 0.43 |
| shape | 0.5 | 0.058 / 0.048 | 0.038 / 0.052 (+0.014, 0.013) | 0.048 / 0.042 | 0.43 |
| shape | 0.75 | 0.036 / 0.034 | 0.042 / 0.062 (+0.020, 0.014) | 0.040 / 0.040 | 0.48 |
| shape | 1 | 0.054 / 0.040 | 0.030 / 0.060 (+0.030, 0.013) | 0.050 / 0.042 | 0.51 |

At these strengths the surgeries scale a stage effect that is itself near the noise (magnitude multiplies δ by
1 + e; orientation, shape and translation relocate differential sites whose shift is δ), so every rate in both arms
sits at the null, 0.010–0.062. The construction therefore measures size, not power: nothing here says the floor
makes `angle` or `shape` powerful. The `angle` differences of +0.016 to +0.036, several above two pooled SEs, are
the unfloored arm's under-rejection at rank 1 (its `angle` rates of 0.010–0.042 are below α, as at the anchor), not
detection; `delta` and `shape` differ by at most ±0.014 between the arms.

## 5. Decision, limitation and scope

**Decision**, read verbatim from `results/rank-floor-2026-10-23/delta0.05/rank_floor_decision.json` and
`delta0.1/rank_floor_decision.json`: `"decision": "adopt_floor"` in both, with `"exceeding": []` and every
`floor_2` verdict `"holds"` at the anchor and the translation control. Under the predeclared conjunction the
check adopts the floor: MOTCO 0.10.0 keeps `min_components=2` as the default of `fit_pls_latent_space`,
`SimulationEvaluationParams` and `TrajectoryAnalysisParams` / `motco analyze`.

**Limitation: the floor guarantees two dimensions, not a second axis with signal.** When the stage signal supports
one component, or none, the second component is the best of the remaining candidates and carries whatever the inner
CV preferred among them; `angle` and `shape` are then defined but their nulls are wide, and at these strengths
their power is at α. A rank-2 selection under the floor is a statement about the measurement space, not evidence of
a second biological axis. Read the record's `config_spectrum` relative eigengap and `null_summary.angle` together:
they say how resolvable orientation is in each analysis, and the selected rank stays the covariate reported beside
any `angle` result. The 0.9.0 rule remains reachable as `min_components=1` for comparison.

**Scope.**

- The construction weakens the per-omic effect sizes uniformly; `p_dmp`, the cell sizes, the block set and the
  joint magnitude are the Phase 6 study's. It reproduces the SEA-AD cohorts' held-out AUROC regime, not their
  mechanism (a dominant technical axis aligned with stage), which the covariate-aware latent construction addresses.
- Both arms used the harness CV defaults and no covariates; the floor applies inside the same cross-validation
  with covariates and latent adjustment, but their interaction was not measured here.
- Every committed study config, fixture and dated report pins `min_components: 1` and is unchanged by the default
  (`tests/test_rank_floor_default_pins.py`); the sead_sexdiff cohorts are to be rerun at the floor in their own
  proposal, pinning MOTCO 0.10.0.
- Replicates are matched across cells and arms by seed; the two strengths are independent runs with their own
  decision files, and the report's conjunction rule was fixed with the user before submission.
