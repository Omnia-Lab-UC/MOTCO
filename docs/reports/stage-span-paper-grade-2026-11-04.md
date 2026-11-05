# Stage-span paper grade — findings report (2026-11-04)

**Change:** `stage-span-paper-grade` (MOTCO 0.17.0). **Reads:** `results/stage-span-paper-grade-2026-11-04/`. It holds:

- `report/` for each of the three profiles;
- `stage_span_paper_grade_decision.json`;
- `leak/`, the pre-flight arms and `preflight_gate.json`;
- `reproduction/`, `PROVENANCE.txt` and `sacct.txt`.

**Concept page:** [measurement visibility](../measurement-visibility.md).

The [stage-span pilot](stage-span-pilot-2026-11-04.md) showed that `angle` and `shape` are sensitive once a group difference lies inside the stage-supervised span, but only at pilot grade (100 × 199). It also left the bend's `delta` leak unexplained. This run states sensitivity given visibility at paper grade (500 × 999) at five columns, with a construction whose leaks are measured.

**In one paragraph.** The predeclared decision is **`sensitivity_established`**.

- **Size.** The anchors and the stage-span translation control hold size at every column: 0 of 60 rates above the 0.0695 bound.
- **Capture.** Every stage-span plant keeps 0.92–0.97 of its energy in the fitted span. The surgeries keep 0.2–0.8 % at e = 0.5 at the same columns.
- **Power given visibility.** A 45° stage-span rotation rejects `angle` in 0.96–1.00 of replicates and a half-step bend rejects `shape` in 0.97–1.00, at every column including n = 80 and the SEA-AD table. The surgeries reach 0.08–0.67 at the same effect.
- **Specificity.** 3 of 90 off-target rates are above the per-cell bound. The bound rule allows 6, and none is above the Bonferroni ceiling.
- **The planted scale.** It solves each edit as a fixed point so that it is exact where the sample scaler looks. That removes the bend's `delta` leak, now 0.016–0.036.
- **The residual.** What remains is the PLS score metric, a sheared map of the span: −0.48 % of A's path on the SEA-AD table and −0.59 % at four stages, for the bend at e = 0.5. It is now documented ([measurement visibility](../measurement-visibility.md#score-coordinates-the-span-is-seen-through-a-shear), roadmap item 15).
- **Two qualifications.**
  - The rotation's off-target `shape` rises with the angle, reaching 0.086 at four stages at 45°: exactly at that cell's ceiling.
  - 12 of the 1,000 four-stage orientation units have no planted plant and were reported as failures.

## 1. Configuration and provenance

| | stage-span geometry | stage-span magnitude | surgery geometry |
|---|---|---|---|
| config | `stage_span_paper_grade_span_geometry.json` | `stage_span_paper_grade_span_magnitude.json` | `stage_span_paper_grade_surgery_geometry.json` |
| config sha256 (first 16) | `6697d927854b3899` | `4e2ef7d3905438f7` | `e1cae479c0e42a9c` |
| construction | `stage_span`, `span_scale: "planted"` | same | surgery (unset) |
| modes × effects | `orientation`, `shape`, `translation` × {0.1, 0.25, 0.5} | `magnitude` × {0.02, 0.05, 0.1} | `orientation`, `shape`, `translation` × {0.25, 0.5, 1.0} |
| columns | five (below) | five | SEA-AD table, n = 200 at four stages |
| cells / units | 52 / 26,000 | 22 / 11,000 | 22 / 11,000 |
| SLURM job (shards) | 933810 (260) | 933811 (110) | 933812 (110) |

**The five columns**, each measured on its span layers:

- n = 200 (the baseline) and n = 80, at three stages on methylation + expression;
- the SEA-AD table `[[11, 10, 28], [9, 10, 12]]` on methylation + expression;
- n = 200 on all three blocks;
- n = 200 at four stages on methylation + expression.

**Shared settings:**

- 500 replicates × 999 permutations; `base_seed` 1028 in the `design-guidance` family; `evaluation.seed` 3, which is also the RRPP seed.
- The design-guidance generator: `p_dmp` 0.1, ρ = 0, every per-omic δ 2.0, `joint` magnitude.
- The production measurement: pooled PLS with the double-CV rank floored at 2, `net_displacement`, `span_capture` on.
- **Code:** revision `e41c1ae` in a clean detached worktree on `ing`, partition `512x1024` (n1, n3, n4; n2 and n11 excluded, reasons in `PROVENANCE.txt`), one CPU and 4 GB per task, BLAS pinned to one thread.

**Pre-flight before submission:**

- the change's suites (258 passed);
- a 2 × 19 smoke of all three profiles (192/192 units);
- 32 random units timed at 999 permutations (mean 17.3 s, about 230 core-hours projected).

**Order of work (user decision).** The workstation had to shut down before the leak-decomposition script existed. The user chose to submit the paper grade that night under a pre-approved cap of 450 projected core-hours. The leak check and its gate, which the design placed before the paper grade, ran the next morning as a descriptive companion (§6). Its gate holds, so the reordering changed nothing that would have stopped the run. The decision script was committed (`cbd3e57`) before any rate was read.

## 2. Unit, failure and reproduction accounting

**Units.** All 480 array tasks finished `COMPLETED 0:0` (22:25–01:13). There are 48,000 records, 48,000 distinct, with no duplicates and 229.0 recorded core-hours (160.9 / 38.1 / 30.0). Every completed record carries `net_displacement`, `min_components` 2, `n_jobs` 1 and `span_capture` true.

**Failures.** 12 units failed, all in stage-span `orientation` at n = 200 / four stages: replicates 311, 409 and 481 at e = 0.25, and 123, 266, 311, 368, 409, 425, 439, 477 and 481 at e = 0.5.

- **Cause.** The planted fixed point did not converge: the planted pooled SD cycled in two steps after 200 updates.
- **Diagnosis.** In those datasets the configuration with the net displacement removed has two nearly tied singular values (ratio 0.95–0.97, against 0.86 in a replicate that converges). The rotation's second direction is chosen from those two, so it swings with the scale.
- **No cheap fix.** Damping the iteration by half-steps rescued 2 of the 12, so the other 10 most likely have no consistent planted plant.
- **Determinism.** A local scan of all 1,000 units of the two cells finds exactly these 12.

The user chose to report them as failures, as the design's risk entry anticipated. The two cells' rates use the completed replicates (R = 497 and 491), and every other cell is complete. The failures select against near-tied configurations, so those two cells describe configurations with a well-defined rotation plane.

**Reproduction.** Run on `ing` with `scripts/anchor_reproduction.py`, from `reproduction/`:

| comparison | pairs | match |
|---|---|---|
| stage-span magnitude vs geometry anchors (5 columns, `--by-column`) | 2,500 | 2,499 (seeds, ranks, p-values 2,500) |
| surgery vs stage-span geometry anchors (2 columns) | 1,000 | 1,000 |
| anchors vs the 2026-10-29 slice (n = 200, n = 80, three blocks; 999 permutations both) | 1,500 | 1,499 (p-values 1,500) |
| SEA-AD surgery cells vs the stage-span pilot, replicates 0–99 | 1,000 | 1,000 (p-values not comparable) |
| four-stage surgery cells vs the design-guidance pilot, replicates 0–99 | 1,000 | 1,000 (p-values not comparable) |
| planted magnitude vs the pilot's stage-span magnitude (n = 200, n = 80, SEA-AD), replicates 0–99 | 1,200 | 1,200 |

Both misses are one anchor dataset (replicate 222 at n = 200), whose observed `angle` is 0.00019°. It differs by 1.9e-9° between runs, at the limit of the arccos precision near zero, and its p-value is 1.0 everywhere. The planted magnitude cells reproduce the null-scale pilot because scaling about the centroid gives the same plant under either scale.

## 3. The decision

`stage_span_paper_grade_decision.json`, rule version 1. **Bound rule.** For each criterion over m rates:

- A rate *exceeds* when it is above α + 2·√(α(1 − α)/R), with R the row's available replicates.
- The criterion fails if the number of exceedances passes the 95th percentile of Binomial(m, p₂), where p₂ is the exact probability that a Binomial(R, α)/R rate exceeds that bound.
- It also fails if any rate is above the Bonferroni ceiling at 0.05/m, or if any row is unavailable.

At R = 500, p₂ = 0.0303.

| criterion | what | m | allowed | observed | ceiling | holds |
|---|---|---|---|---|---|---|
| P1 size | anchor (5 columns × 3 statistics) and stage-span translation (5 × 3 effects × 3) | 60 | 4 | 0 | 0.084 | yes |
| P2 capture | median `group_difference` ≥ 0.90 in every stage-span magnitude, orientation and shape cell | 45 cells | — | 0 below | — | yes |
| P3 specificity | off-target statistics of those cells | 90 | 6 | 3 | 0.084 | yes |

Verdict: **`sensitivity_established`**. The three exceedances are listed in §8.

## 4. Visibility: capture

Median `group_difference` share per cell, with its q10–q90 range across the column's cells:

| column | magnitude | orientation | shape | median rank |
|---|---|---|---|---|
| n = 200 | 0.970–0.971 | 0.970–0.971 | 0.967–0.970 | 2 (3 at 45°) |
| n = 80 | 0.936–0.937 | 0.927–0.935 | 0.927–0.934 | 2–3 |
| SEA-AD table | 0.934–0.935 | 0.924–0.933 | 0.931–0.936 | 3 |
| three blocks | 0.972 | 0.970–0.971 | 0.967–0.971 | 2 |
| four stages | 0.969 | 0.966–0.969 | 0.968 | 3–4 |

Every cell's q10 is at least 0.903. Through the fitted basis the plants read close to their intended geometry. The rotation reads its angle to within 1.3°: 45° reads 43.7° on the SEA-AD table and 44.8–45.2° elsewhere. The bend's `shape` reads within 0.004, and the scaling's `delta` within 2.8 % (the `population_latent` columns of the decision's `target_power_and_capture`).

## 5. Size

Anchor and stage-span translation rejection rates (500 replicates; bound 0.0695):

| column | anchor `delta` / `angle` / `shape` | translation, all three effects |
|---|---|---|
| n = 200 | 0.040 / 0.060 / 0.034 | 0.012–0.062 |
| n = 80 | 0.022 / 0.036 / 0.010 | 0.006–0.036 |
| SEA-AD table | 0.018 / 0.052 / 0.012 | 0.004–0.052 |
| three blocks | 0.038 / 0.060 / 0.062 | 0.036–0.062 |
| four stages | 0.028 / 0.052 / 0.012 | 0.002–0.062 |

Every rate is at or below 0.062. The report's own advisory anchor decision reads `SIZE_HOLDS` in all three profiles. The surgery profile's anchors are the same datasets (§2).

## 6. The leak decomposition

`scripts/stage_span_leak_check.py` ran two arms at 100 × 199 on the five columns (jobs 934447 and 934448, revision `b40fb91`; 3,500 units each, 0 failed, 32.0 core-hours). It covered the anchor, `orientation` and `shape`, under `planted` at e ∈ {0.1, 0.25, 0.5} and under `null` at {0.25, 0.5, 1.0}. The null arm reproduces the stage-span pilot's rejection counts at its three columns (63 of 63).

For each unit the check reports B's path length minus A's at six checkpoints. The table gives the median for the bend at e = 0.5, as a share of A's path:

| column | arm | null scale | planted scale | sample scaler | orthonormal span | PLS scores | sample latent means |
|---|---|---|---|---|---|---|---|
| n = 200 | null | +0.00 % | −0.83 % | −0.82 % | −0.77 % | −0.79 % | −1.12 % |
| | planted | +0.86 % | −0.00 % | −0.02 % | +0.05 % | +0.03 % | −0.33 % |
| n = 80 | null | +0.00 % | −0.82 % | −0.80 % | −0.72 % | −0.82 % | −0.71 % |
| | planted | +0.85 % | −0.00 % | −0.00 % | +0.11 % | −0.02 % | +0.05 % |
| SEA-AD table | null | +0.00 % | −0.65 % | −0.65 % | −0.64 % | −1.08 % | +0.07 % |
| | planted | +0.67 % | −0.00 % | −0.01 % | −0.01 % | −0.48 % | +0.62 % |
| three blocks | null | +0.00 % | −0.90 % | −0.89 % | −0.84 % | −0.78 % | −0.83 % |
| | planted | +0.94 % | −0.00 % | −0.01 % | +0.06 % | +0.10 % | +0.00 % |
| four stages | null | +0.00 % | −1.06 % | −1.06 % | −0.99 % | −1.62 % | −2.41 % |
| | planted | +1.15 % | −0.00 % | −0.02 % | +0.06 % | −0.59 % | −1.40 % |

**Reading the routes:**

- **The scale route.** Under the null scale the bend keeps its path length in the null coordinates it is defined in. The planted data's own pooled SD shortens it by 0.65–1.06 %, and the sample scaler sees the same. The planted plant is exact there in every unit (G1), and the scaler follows it.
- **The score-metric route.** On an orthonormal basis of the fitted span the planted bend keeps its length to within 0.11 %. In the PLS score coordinates it does not on the SEA-AD table (−0.48 %) or at four stages (−0.59 %).
- **The rotation shows the shear more strongly.** The planted 45° rotation on the SEA-AD table reads −3.35 % in scores against +0.66 % on the orthonormal basis.
- **Size is unaffected by either route** (§5).

**Rates.** Under the null scale the bend's `delta` rate at R = 100 is 0.08 at n = 200 and 0.11 on three blocks at e = 0.5. At e = 1.0 it is 0.17 at n = 200 and 0.33 at four stages. Under the planted scale it is 0.00–0.04 at every column and effect.

**Pre-flight gate.** `preflight_gate.json` reads **`gate_holds`**:

- **G1.** Every planted unit generated, with a planted-scale gap ≤ 1e-8 of A's path.
- **G2.** 0 of 60 planted off-target rates above the R = 100 bound (allowed 4, ceiling 0.130).
- **G3.** The planted arm's median score-coordinate difference is smaller than the null arm's at every column (−0.28 → +0.01, −0.29 → −0.01, −0.36 → −0.16, −0.35 → +0.05, −0.92 → −0.33, in A's path units).

The gate ran after the paper grade (§1), so it is a diagnostic here, not a condition the run waited on.

## 7. Sensitivity given visibility: target power

Target rejection rate (± SE about 0.02 or less) and the smallest effect with power ≥ 0.80:

| column | orientation, 9° / 22.5° / 45° | shape, e = 0.1 / 0.25 / 0.5 | magnitude, e = 0.02 / 0.05 / 0.1 | informative from (orientation / shape / magnitude) |
|---|---|---|---|---|
| n = 200 | 0.738 / 0.880 / 1.000 | 0.704 / 0.974 / 1.000 | 0.228 / 0.694 / 0.936 | 0.25 / 0.25 / 0.1 |
| n = 80 | 0.400 / 0.412 / 0.992 | 0.210 / 0.724 / 0.972 | 0.036 / 0.236 / 0.624 | 0.5 / 0.5 / — |
| SEA-AD table | 0.382 / 0.376 / 0.962 | 0.208 / 0.678 / 0.968 | 0.056 / 0.188 / 0.542 | 0.5 / 0.5 / — |
| three blocks | 0.984 / 1.000 / 1.000 | 0.914 / 1.000 / 1.000 | 0.430 / 0.878 / 0.992 | 0.1 / 0.1 / 0.05 |
| four stages | 0.478 / 0.744 / 1.000 | 0.478 / 0.962 / 1.000 | 0.172 / 0.648 / 0.940 | 0.5 / 0.25 / 0.1 |

**What the table shows:**

- **A visible 45° rotation and a half-step bend are detected at every column.** That includes n = 80 and the unequal SEA-AD table.
- **Orientation power is flat between 9° and 22.5° at n = 80 and on the SEA-AD table** (0.38–0.41), then jumps at 45°. The design-guidance runs recorded the same flat response of `angle` in e at small n.
- **Magnitude is the hardest at small n.** A 10 % scaling reaches 0.54–0.62 at n = 80 and on the table.
- **A third block raises every statistic's power** at n = 200.

## 8. Specificity

Off-target rejection rates of the stage-span cells (500 replicates; bound 0.0695, ceiling 0.084 at R = 500):

| mode | statistic | n = 200 | n = 80 | SEA-AD | three blocks | four stages |
|---|---|---|---|---|---|---|
| magnitude | `angle` | 0.050–0.058 | 0.032–0.038 | 0.050–0.052 | 0.036–0.058 | 0.040–0.050 |
| magnitude | `shape` | 0.036–0.038 | 0.008–0.010 | 0.012 | 0.042–0.062 | 0.012–0.014 |
| orientation | `delta` | 0.002–0.018 | 0.000–0.014 | 0.000–0.012 | 0.000–0.018 | 0.018–0.032 |
| orientation | `shape` | 0.030 / 0.046 / **0.070** | 0.004–0.006 | 0.012–0.022 | 0.002–0.022 | 0.004 / 0.030 / **0.086** |
| shape | `angle` | 0.050–0.062 | 0.034–0.042 | 0.054 / 0.066 / **0.078** | 0.048–0.054 | 0.044–0.050 |
| shape | `delta` | 0.024–0.036 | 0.022–0.024 | 0.016–0.020 | 0.022–0.036 | 0.020–0.030 |

The three exceedances (bold) are:

- the 45° rotation's `shape` at n = 200 (35/500 = 0.070);
- the 45° rotation's `shape` at four stages (42/491 = 0.086). At R = 491 that cell's ceiling is 42/491, so it sits exactly at the ceiling; one more rejection would have failed P3;
- the half-step bend's `angle` on the SEA-AD table (39/500 = 0.078).

**Patterns.**

- The rotation's off-target `shape` rises steadily with the angle at n = 200 and at four stages, so it is a systematic leak, not noise. The leak check reads path lengths and does not decompose it. The plant keeps `shape` exactly in the planted coordinates, so the leak enters through the measurement. The score-metric shear is the candidate the decomposition supports, but this run does not isolate it.
- The bend's `delta`, the pilot's failing reading, is 0.016–0.036 at every column and effect.
- At n = 200 the anchor `angle` (0.060) and every `angle` cell read 0.05–0.06. That is the column's shared noise, which the pilot saw at 0.10 at R = 100.

## 9. The visibility contrast

The surgeries' target power and capture at the same columns (`descriptive.visibility_contrast`; each source is labelled in the JSON):

| column | source | orientation `angle`, e = 0.25 / 0.5 / 1.0 | shape `shape`, e = 0.25 / 0.5 / 1.0 | capture, e = 0.5 (orientation / shape) |
|---|---|---|---|---|
| n = 200 | slice power (paper); pilot capture (pilot) | 0.158 / 0.242 / 0.184 | 0.236 / 0.418 / 0.562 | 0.002 / 0.004 |
| n = 80 | slice; pilot | 0.100 / 0.110 / 0.150 | 0.044 / 0.104 / 0.204 | 0.003 / 0.004 |
| SEA-AD table | this run (paper) | 0.110 / 0.152 / 0.690 | 0.062 / 0.084 / 0.206 | 0.008 / 0.006 |
| three blocks | slice; capture unavailable | 0.406 / 0.494 / 0.492 | 0.486 / 0.668 / 0.856 | — |
| four stages | this run (paper) | 0.094 / 0.118 / 0.136 | 0.132 / 0.292 / 0.552 | 0.003 / 0.008 |

At e = 0.5 the surgeries keep 0.2–0.8 % of their energy in the span and reach `angle` 0.11–0.49 and `shape` 0.08–0.67. The stage-span plants at the same effect reach 0.96–1.00. The exception is the surgery's full relocation on the SEA-AD table (capture 0.217, `angle` 0.690), the unequal-cell "kite" route the pilot found. The surgeries' low power measures visibility, not the tests.

## 10. Limitations

- **The construction is circular by design.** It plants where the measurement looks, so it measures sensitivity given visibility. How much of a real difference the span sees is the latent plasmode's question (roadmap item 13 (b)).
- **The paper grade ran before its pre-flight**, on the user's decision. The gate held when it ran (§6), and the decision rule was committed before any rate was read.
- **Excluded configurations.** 12 four-stage orientation units have no planted plant (near-tied in-plane directions). The two cells read over 497 and 491 replicates and describe configurations with a well-defined rotation plane.
- **The score-coordinate residual remains.** It is −0.48 % and −0.59 % of A's path for the bend at e = 0.5 on the SEA-AD table and at four stages, and −3.35 % for the 45° rotation on the table. It did not produce a `delta` leak at these effects. Removing it would mean measuring on orthonormal span coordinates, which changes the production measurement and is unvalidated (roadmap item 15).
- **The rotation's off-target `shape` is at the per-cell bound at 45°**, and sits at the ceiling at four stages. Larger stage-span rotations were a non-goal and were not run; they would likely fail P3.
- **The count rule assumes independent rates.** Rates within a column share matched-seed noise, so exceedances can cluster (n = 200 reads 0.05–0.06 on every `angle` cell). Every exceeding row is listed.
- **The surgeries' capture is measured only at pilot grade** at n = 200 and n = 80, and not at all on three blocks.
- **Simulator scope.** The structural-analogue limitation applies: InterSIM blocks, 498 measured features, ρ = 0.
