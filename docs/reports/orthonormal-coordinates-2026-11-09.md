# Orthonormal span coordinates — findings report (2026-11-09)

**Change:** `orthonormal-span-coordinates` (MOTCO 0.19.0). **Reads:** `results/orthonormal-coordinates-2026-11-09/` (`orthonormal_coordinates_decision.json`; `paper_grade/` and `weak/` rate, paired and cell tables, `accounting.json`, `reproduction.json`, `settings.json` and the post-hoc `population_geometry.csv`; `PROVENANCE.txt`, `sacct.txt`). **Concept page:** [measurement visibility](../measurement-visibility.md#score-coordinates-the-span-is-seen-through-a-shear).

The production statistics are measured in PLS score coordinates, `T = X_in·R` with `R = W·(PᵀW)⁻¹`, a sheared image of the fitted span. The stage-span paper grade's leak decomposition ([report](stage-span-paper-grade-2026-11-04.md#6-the-leak-decomposition)) found that the shear moves a planted bend's path length by −0.48 % on the SEA-AD table and −0.59 % at four stages. On an orthonormal basis of the same span the bend kept its length to within 0.11 %.

MOTCO 0.19.0 adds the opt-in setting `SimulationEvaluationParams.latent_coordinates = "orthonormal"`. It measures on `X_in·Q`, `Q = qr(x_rotations_)`, with the same fit, rank and draws. `delta`, `angle` and `shape` do not change under any orthogonal map, so every orthonormal basis of the span gives the same statistics. This run compared the two coordinate systems on matched units under a decision fixed in advance (`design.md` D7).

**In one paragraph.** **The verdict is `power_loss`.** The orthonormal coordinates hold size and remove the score shear's off-target leak. But they lose target power at the small, unequal and four-stage designs, so production stays on the scores.

- **Size holds.** The 60 paper-grade controls read 0.000–0.066 and the 6 weak-design controls 0.010–0.050. None is above its bound.
- **Specificity holds, and the leak goes.** At the 45° rotation, off-target `shape` falls from 0.070 to 0.000 at n = 200 and from 0.086 to 0.000 at four stages. 9 of 90 off-target cells read lower with orthonormal coordinates, and none higher, at two-sided paired p < 0.05.
- **Power is lost.** 19 of 45 target cells lose, with one-sided paired p < 0.025; 3 losses are allowed, and 16 are below the Bonferroni ceiling. The largest losses:
  - the SEA-AD table's 45° rotation, `angle` 0.962 → 0.760;
  - the 22.5° rotation at four stages, 0.744 → 0.519;
  - the e = 0.25 bend, `shape` 0.724 → 0.518 at n = 80 and 0.678 → 0.482 on the SEA-AD table;
  - the e = 0.1 magnitude on the SEA-AD table, `delta` 0.542 → 0.416.

  No target cell gains. The three-block column loses nothing.
- **The scores do not win by enlarging the planted effect, except for magnitude.** The population-level target `angle` and `shape` are the same in both coordinate systems: 43.7°–45.2° for the 45° rotation, and 0.1255–0.1303 for the e = 0.25 bend. The scores' advantage there lies on the sampling side. Where the orthonormal coordinates lose, their `angle` null is wider: at the SEA-AD anchor its 95th percentile is 34.6° against the scores' 26.1°. For magnitude, the scores do stretch the planted size, by 1.6–5.0 %.

## 1. Configuration and provenance

- **The unit.** Each unit generates one dataset and fits stage-supervised PLS once (double-CV rank, floor 2). The two arms share that fit:
  - **`scores`:** the fit's production evaluation;
  - **`orthonormal`:** `pls_orthonormal_from_pls`, the same fit in orthonormal coordinates, with no refit.

  Both arms run at the unit's RRPP seed with 999 permutations, serially, under net displacement, with span capture on. Each unit checks to 1e-8 that `QᵀQ = I` and that `T = Z·(QᵀR)`. A unit failing either check is inconsistent, and an inconsistent unit fails the decision's gate.
- **`paper_grade` profile.** It reuses the stage-span paper grade's units under the planted scale, through `stage_span_leak_check.column_cells` / `unit_params`, imported unchanged:
  - **Columns (5):** n = 200, n = 80 and the SEA-AD table (all on two blocks), n = 200 on three blocks, and n = 200 at four stages.
  - **Cells (13):** the anchor; translation, orientation and shape at e 0.1 / 0.25 / 0.5; magnitude at e 0.02 / 0.05 / 0.1.
  - **Size:** 500 replicates per cell, 32,500 units.
- **`weak` profile.** The latent-space comparison's `seaad_weak` anchor and translation (e = 1.0), at δ = 0.11 from its committed calibration and 500 replicates: 1,000 units, the span-split check's units. It reads size only.
- **Rejection rule:** p < 0.05, as the parents count it.
- **Code.** Revision `d776366` on `feat/orthonormal-span-coordinates`, run from a clean detached worktree on `ing`. Python 3.11.16, numpy 2.3.5, scipy 1.16.3 and scikit-learn 1.8.0, with BLAS pinned to one thread.
- **Cluster.** Partition `512x1024`, n11 excluded:
  - array 936125 (`paper_grade`): 650 shards, at most 96 at once;
  - array 936126 (`weak`): 20 shards.
- **Pre-flight.**
  - The change's test modules and the identity suites passed on a compute node (95 tests).
  - A smoke of both profiles ran.
  - Timing at 999 permutations gave a projection of about 206 core-hours. The user approved a 300 core-hour cap before submission (`notes.md`).

Details are in `PROVENANCE.txt`.

## 2. Accounting and reproduction

- **Units.**
  - `paper_grade`: 32,500 of 32,500 records, of which 12 failed with `SemiSyntheticTrajectoryError`: 9 in the four-stage orientation cell at e = 0.5 and 3 at e = 0.25. These are the planted-scale fixed point the stage-span paper grade also failed on, in the same cells and counts. A failed unit fails before the fit, so it is missing from both arms, and the pairing uses the 497 and 491 units available there.
  - `weak`: 1,000 of 1,000.
  - Both profiles: 0 inconsistent units and 0 duplicates. All 670 array tasks ended `COMPLETED 0:0` (`sacct.txt`).
- **Consistency.** The largest deviations were 2.0e-15 for the orthonormal basis and 3.4e-13 for the same span, against the 1e-8 tolerance.
- **Cost.** 173.3 + 2.5 = 175.8 recorded core-hours. The median unit took 18.2 s (paper grade) and 8.9 s (weak); the PLS fit took 16.6 s, and each arm's RRPP pass 0.70 s. That compares with the 206 projected and the 300 cap.
- **Reproduction (the decision's gate): both passed exactly.**
  - **`paper_grade`:** the `scores` arm's available and rejected counts equal the stage-span paper grade's committed `design_point_operating.csv` in all 195 rows (5 columns × 13 cells × 3 statistics).
  - **`weak`:** the counts and rejections equal the span-split check's `production_table.csv` in all 6 rows.
  - **So** every `scores` rate below is the parents' rate. Every `orthonormal` rate pairs with it on the same datasets, fits and permutation draws.

## 3. The decision

`orthonormal_coordinates_decision.json` (rule version 1) returns **`power_loss`**. O1 and O2 apply the stage-span paper grade's count-plus-ceiling rule (`judge_rates`, imported unchanged) to the `orthonormal` arm.

| criterion | rule | result |
|---|---|---|
| gate | both reproductions pass and no unit is inconsistent | **passed** |
| O1 size | anchor and translation, `delta`/`angle`/`shape`: m = 60 at the five columns (allowed 4 above the bound 0.0695, ceiling 0.084), and m = 6 at the weak design (allowed 1, ceiling 0.074) | **holds**: 0 of 60 (0.000–0.066) and 0 of 6 (0.010–0.050) |
| O2 specificity | off-target `delta`/`shape` of orientation, `delta`/`angle` of shape and `angle`/`shape` of magnitude, m = 90 (allowed 6, ceiling 0.084) | **holds**: 1 of 90 above the bound, the SEA-AD bend at e = 0.5 on `angle` at 0.072; range 0.000–0.072 |
| O3 power | per target cell (m = 45), one-sided exact paired loss p = P(Binomial(b + c, ½) ≥ b). Holds if at most 3 cells fall below 0.025 (the 95th percentile of Binomial(45, 0.025)) and none below 0.05/45 = 0.0011 | **fails**: 19 below 0.025, 16 of them below the ceiling |

- **The same families in the `scores` arm** read 0.002–0.062 at the paper-grade controls and 0.026–0.052 at the weak design. They have 3 of 90 off-target rates above the bound: the 45° rotation's `shape` at 0.070 (n = 200) and 0.086 (four stages), and the SEA-AD bend's `angle` at 0.078. These are the paper grade's P3 readings.
- **The weak design's controls.** With orthonormal coordinates `delta` is conservative there: 0.012 (anchor) and 0.010 (translation), against the scores' 0.034 and 0.036. `angle` reads 0.050 and 0.048, against 0.038 and 0.026.

**Consequence** (the decision's text): "Production stays on scores: the orthonormal coordinates lose target power against the scores at the failing cells. The report documents the trade-off and the setting stays an opt-in diagnostic."

## 4. Where the power is lost

Each entry gives the `scores` → `orthonormal` rate over 500 paired units (497 and 491 in the four-stage orientation cells at e = 0.25 and 0.5). In brackets is b/c: b counts units only the scores reject, c units only the orthonormal coordinates reject. **‡** marks a loss below the ceiling 0.0011 and **†** a loss below 0.025 that is above the ceiling.

**Orientation, `angle`:**

| column | e = 0.1 | e = 0.25 | e = 0.5 |
|---|---|---|---|
| n = 200 | 0.738 → 0.736 (3/2) | 0.880 → 0.720 (83/3) **‡** | 1.000 → 1.000 (0/0) |
| n = 80 | 0.400 → 0.396 (4/2) | 0.412 → 0.366 (28/5) **‡** | 0.992 → 0.892 (50/0) **‡** |
| SEA-AD table | 0.382 → 0.382 (1/1) | 0.376 → 0.350 (19/6) **†** | 0.962 → 0.760 (103/2) **‡** |
| three blocks | 0.984 → 0.984 (0/0) | 1.000 → 0.994 (3/0) | 1.000 → 1.000 (0/0) |
| four stages | 0.478 → 0.464 (8/1) **†** | 0.744 → 0.519 (117/5) **‡** | 1.000 → 1.000 (0/0) |

**Shape, `shape`:**

| column | e = 0.1 | e = 0.25 | e = 0.5 |
|---|---|---|---|
| n = 200 | 0.704 → 0.672 (46/30) | 0.974 → 0.966 (16/12) | 1.000 → 0.998 (1/0) |
| n = 80 | 0.210 → 0.174 (24/6) **‡** | 0.724 → 0.518 (150/47) **‡** | 0.972 → 0.898 (49/12) **‡** |
| SEA-AD table | 0.208 → 0.170 (22/3) **‡** | 0.678 → 0.482 (132/34) **‡** | 0.968 → 0.846 (76/15) **‡** |
| three blocks | 0.914 → 0.918 (2/4) | 1.000 → 0.998 (1/0) | 1.000 → 1.000 (0/0) |
| four stages | 0.478 → 0.328 (82/7) **‡** | 0.962 → 0.944 (28/19) | 1.000 → 1.000 (0/0) |

**Magnitude, `delta`:**

| column | e = 0.02 | e = 0.05 | e = 0.1 |
|---|---|---|---|
| n = 200 | 0.228 → 0.210 (15/6) | 0.694 → 0.656 (42/23) **†** | 0.936 → 0.918 (25/16) |
| n = 80 | 0.036 → 0.040 (3/5) | 0.236 → 0.212 (23/11) | 0.624 → 0.506 (95/36) **‡** |
| SEA-AD table | 0.056 → 0.058 (4/5) | 0.188 → 0.170 (22/13) | 0.542 → 0.416 (95/32) **‡** |
| three blocks | 0.430 → 0.428 (2/1) | 0.878 → 0.872 (4/1) | 0.992 → 0.988 (2/0) |
| four stages | 0.172 → 0.134 (21/2) **‡** | 0.648 → 0.560 (98/54) **‡** | 0.940 → 0.938 (28/27) |

- **The losses sit where the stage program is hardest to resolve.** They fall at n = 80, on the SEA-AD table and at four stages: 6, 6 and 5 of their 9 target cells. At n = 200 on two blocks there are 2, and on three blocks none. They are spread across the modes: 7 orientation, 7 shape and 5 magnitude cells.
- **The losses are one-sided.** In the 45 target cells c exceeds b only three times, by at most 2 units (n = 80 and SEA-AD magnitude at e = 0.02, the three-block bend at e = 0.1).
- **They are not small.** At 500 paired units the differences run up to 0.225 (four stages, 22.5°), 0.206 (n = 80, bend at e = 0.25) and 0.202 (SEA-AD, 45°).

## 5. Specificity: the watched cells

The watched cells are the 45° rotation's off-target `shape`. Before the run, they were the paper grade's two largest P3 readings that the shear could explain.

| column | `scores` | `orthonormal` | b / c | two-sided paired p |
|---|---|---|---|---|
| n = 200 | 0.070 | 0.000 | 35 / 0 | 5.8e-11 |
| four stages (491 units) | 0.086 | 0.000 | 42 / 0 | 4.5e-13 |

The orthonormal coordinates remove the shear's leak, and not only at the watched cells. Nine off-target cells read lower at two-sided p < 0.05, with b/c in brackets:

- **The rotation's `shape`:**
  - n = 200: 0.046 → 0.004 at e = 0.25 (21/0) and 0.070 → 0.000 at e = 0.5 (35/0);
  - SEA-AD table: 0.022 → 0.002 at both effects (11/1);
  - four stages: 0.030 → 0.000 at e = 0.25 (15/0) and 0.086 → 0.000 at e = 0.5 (42/0).
- **The rotation's `delta` at four stages:** 0.020 → 0.000 at e = 0.5 (10/0).
- **The bend's `delta` at e = 0.5:** 0.034 → 0.012 at n = 200 (11/0) and 0.026 → 0.008 at four stages (11/2).

No off-target cell reads higher. The O2 exceedance, the SEA-AD bend's `angle` at e = 0.5, reads 0.078 in scores and 0.072 in orthonormal coordinates (two-sided p 0.45). It is not a leak the shear caused.

## 6. Population geometry: the residuals and the planted effect

**The bend's residual.** For the planted-scale bend at e = 0.5, B's population path length equals A's. The table gives each arm's median `population_latent` B − A path-length difference, as a share of A's path (`cell_table.csv`, `median_population_share`):

| column | `scores` | `orthonormal` |
|---|---|---|
| n = 200 | +0.06 % | +0.06 % |
| n = 80 | −0.04 % | +0.11 % |
| SEA-AD table | **−0.38 %** | +0.02 % |
| three blocks | +0.09 % | +0.06 % |
| four stages | **−0.60 %** | +0.07 % |

- **The leak check's pattern reproduces at R = 500.** The leak check, at R = 100 and 199 permutations, read −0.48 % and −0.59 % in scores against −0.01 % and +0.06 % on the orthonormal span. With orthonormal coordinates the residual stays within 0.11 % at every column.
- **The 45° rotation's residual on the SEA-AD table** is −1.60 % in scores against +0.52 % with orthonormal coordinates. The leak check's 100 replicates read −3.35 % and +0.66 %. Its median moved with the replicate set, but its sign and its excess over the orthonormal coordinates did not.

**The planted effect is not larger in scores, except for magnitude.** `population_geometry.csv` is a post-hoc descriptive aggregate (`PROVENANCE.txt`). It gives the median population-level statistics per arm, beside the spectrum of the anchor (`cell_table.csv`):

| column | 45° rotation, population `angle` | bend e = 0.25, population `shape` | magnitude e = 0.1, population `delta` | anchor `angle` null q95 | anchor pooled eigengap |
|---|---|---|---|---|---|
| n = 200 | 44.98° / 44.97° | 0.1272 / 0.1272 | 3.377 / 3.323 (+1.6 %) | 2.59° / 2.58° | 0.0881 / 0.0793 |
| n = 80 | 45.15° / 44.86° | 0.1270 / 0.1272 | 3.402 / 3.263 (+4.3 %) | 23.94° / 31.93° | 0.0955 / 0.0820 |
| SEA-AD table | 43.71° / 45.04° | 0.1255 / 0.1257 | 3.561 / 3.391 (+5.0 %) | 26.08° / 34.61° | 0.1117 / 0.0984 |
| three blocks | 45.04° / 44.97° | 0.1274 / 0.1274 | 4.172 / 4.169 (+0.1 %) | 1.55° / 1.55° | 0.0892 / 0.0878 |
| four stages | 44.81° / 45.06° | 0.1289 / 0.1303 | 5.372 / 5.180 (+3.7 %) | 17.09° / 24.46° | 0.0562 / 0.0470 |

Each pair reads `scores` / `orthonormal`. The percentage is the scores' excess.

- **Orientation and shape.** The planted target is the same in both coordinate systems. On the SEA-AD table the scores' population angle is even the smaller (43.7° against 45.0°), yet the scores reject the rotation in 0.962 against 0.760. So the scores' advantage is not a stretched target.
- **Magnitude.** The scores stretch the planted size by 1.6–5.0 % at the columns that lose magnitude power, and by 0.1 % on three blocks, which loses none. The shear's stretch accounts for at least part of the magnitude losses.
- **The null side.** Where the orthonormal coordinates lose, their `angle` null is wider: its 95th percentile rises by 7–9° at n = 80, on the SEA-AD table and at four stages, and the pooled eigengap shrinks. At n = 200 and on three blocks both are nearly unchanged. At the weak design the null widens most, from 53.8° to 100.0° at the anchor. There the median observed null `angle` also rises, from 21.6° to 46.2° (`weak/population_geometry.csv`), and the orthonormal `angle` still holds size (0.050).
- **What the aggregates do not explain.** At n = 200 the 22.5° rotation loses 0.880 → 0.720, yet the medians agree: the observed `angle` is 22.22° / 22.26° and the null q95 4.84° / 4.88°. The loss there lies in the tails of the per-unit null, which the committed medians do not resolve.

## 7. What follows

- **Production and `motco analyze` stay on score coordinates.** The default flip and the `motco analyze` adoption were owed only to `adopt_orthonormal`. Neither is opened.
- **`latent_coordinates = "orthonormal"` stays an opt-in diagnostic in the harness.** It is not exposed in `motco analyze`, in study configurations or with attribution. Without it every output equals 0.18.0.
- **The score shear is a measured trade-off, not a defect to remove.**
  - **The cost.** In scores, the 45° rotation leaks into `shape` at 0.070–0.086 and the e = 0.5 bend into `delta` at 0.026–0.034. With orthonormal coordinates these read 0.000 and 0.008–0.012. Elsewhere the off-target rates are alike in both arms; the bend's `delta` peaks at 0.036 / 0.038 on three blocks at e = 0.1.
  - **The gain.** In exchange, the scores have up to 0.23 more target power at small n, unequal cells and four stages, with equal size.
  - **Reading a rejection.** An off-target `shape` or `delta` rejection beside a strong `angle` or `shape` rejection can be the shear's. At paper grade such a rate stayed within 0.086.
- **Roadmap item 15 is closed with this verdict.** The change does not permit a production follow-up. Why the score metric narrows the null at weak or unbalanced designs is open. The aggregates show the association with the null width and the eigengap, but not the mechanism. Any future work on the measurement metric should take the scores' power, not the orthonormal span, as the benchmark.

## 8. Limitations

- **One construction family.** Every plant is the stage-span construction under the planted scale, exact in planted-standardized coordinates. A real difference may load differently on the shear's directions, so the size of the trade-off is specific to these plants.
- **Medians.** The cell and population tables summarize each cell by medians. The n = 200 rotation loss (§6) shows that they can miss where a loss comes from.
- **`population_geometry.csv` was written after the decision.** It does not bear on the verdict, and it reads the same merged records the decision read.
- **The weak design reads size only.** The latent-space comparison found no power there for any plant, so it was not measured.
- **Units that cannot be planted.** The 12 four-stage orientation units that cannot be planted at the planted scale are missing from both arms. The four-stage orientation rates are therefore over 497 and 491 units.
