# Stage-span pilot — findings report (2026-11-04)

**Change:** `stage-span-constructions` (MOTCO 0.15.0). **Reads:** `results/stage-span-pilot-2026-11-04/` (`report/` per profile, `stage_span_pilot_reading.json`, `retrofit/`, `reproduction/`, `PROVENANCE.txt`, `sacct.txt`). **Concept page:** [measurement visibility](../measurement-visibility.md).

The production measurement sees a group difference only through its component inside the stage-supervised PLS span ([measurement visibility](../measurement-visibility.md)). The simulator's power curves therefore mixed two questions: how much of a planted difference the measurement can see, and how sensitive the test is to what it sees. This pilot separates them. It pairs the feature-set surgeries with a new **stage-span construction**, which plants group B's difference inside the span of B's own stage configuration, on the same baseline programs and noise. Every replicate records how much of the population difference the fitted basis sees (`span_capture`). The design is four sibling profiles: {surgery, stage-span} × {magnitude, geometry} at 100 replicates × 199 permutations. Each has three columns at three stages on methylation + expression: n = 200 (baseline), n = 80, and the SEA-AD table `[[11, 10, 28], [9, 10, 12]]`. All four run in the design-guidance seed family, so the surgery siblings reproduce the 2026-10-28 pilot cell for cell.

**In one paragraph.** **Visibility, not sensitivity, limits the surgeries.**

- **Capture.** The orientation and shape surgeries keep 0.2–1.1 % of their population difference's energy inside the fitted span. That is the share a random direction keeps (rank / 498 = 0.4–0.6 %). The one exception is the full relocation on the unequal SEA-AD table (21.5 %, the "kite" route). The stage-span plants keep 85–97 %. Both magnitude constructions keep 93–97 %.
- **Power given visibility.** At e = 0.5 the stage-span orientation (a 45° rotation) rejects `angle` in 0.95–1.00 of replicates and the stage-span bend rejects `shape` in 0.98–1.00, at every column including n = 80. The surgeries reach 0.06–0.22 and 0.11–0.40 at the same effect. Seen through the fitted basis, the surgery's full relocation is 4.4° at n = 200, against a population 90°; the stage-span 90° rotation is 90.4°.
- **Magnitude.** The two magnitude constructions give the same power to ±0.01 (pooled SE 0.02–0.07).
- **Reading.** The predeclared reading returns **`revise_construction`**. Its failures sit at e = 1.0: the 90° rotation's capture falls to 0.85 at n = 80, and the full-strength bend leaks into `delta` (0.17 at n = 200). One further failure is a 0.10 `angle` reading at n = 200. That reading is inherited from those 100 replicates' shared noise: the shared anchor and every non-`angle` cell at that column read 0.07–0.10 under both constructions.
- **Proposal.** The follow-up should cap the stage-span geometry effects at e ≤ 0.5. Below that every criterion except that inherited reading holds, and the paper grade's 500 replicates should settle it.

## 1. Configuration and provenance

| | surgery magnitude | surgery geometry | stage-span magnitude | stage-span geometry |
|---|---|---|---|---|
| config | `stage_span_pilot_surgery_magnitude.json` | `stage_span_pilot_surgery_geometry.json` | `stage_span_pilot_span_magnitude.json` | `stage_span_pilot_span_geometry.json` |
| config sha256 (first 16) | `1e62a7c52f0de1a4` | `732de3404f76a1ba` | `d4fffa05dd9edb67` | `18e3a1f88fa49796` |
| construction | surgery (unset) | surgery (unset) | `stage_span`, span methylation + expression | same |
| modes × effects | `magnitude` × {0.02, 0.05, 0.1} | `orientation`, `shape`, `translation` × {0.25, 0.5, 1.0} | as surgery magnitude | as surgery geometry |
| cells / units | 14 / 1,400 | 32 / 3,200 | 14 / 1,400 | 32 / 3,200 |
| SLURM job (shards) | 933350 (8) | 933351 (16) | 933352 (8) | 933353 (16) |

**Shared settings (all four profiles):**

- 100 replicates × 199 permutations; `base_seed` 1028 in the `design-guidance` family; `evaluation.seed` 3, which is also the RRPP seed.
- The generator of the design-guidance core profiles: `p_dmp` 0.1, ρ = 0, every per-omic δ 2.0, `joint` magnitude.
- The production measurement: pooled PLS with the stage-supervised double-CV rank floored at 2, `net_displacement`, with the Phase 6 CV knobs.
- `evaluation.span_capture` on.
- **Code:** revision `7bd6bb9` on `feat/stage-span-constructions`, in a clean detached worktree on `ing`.
- **Cluster:** partition `512x1024`, every task on n1 (n2 and n11 excluded), one CPU and 4 GB per task, BLAS pinned to one thread.
- **Pre-flight:**
  - the fast suite: 2,121 passed;
  - the change's nine test modules on `ing`: 135 passed;
  - a 2 × 19 smoke of all four profiles: 184/184 units;
  - eight units timed at 199 permutations: 13–18 s each.

Details are in `PROVENANCE.txt` and `openspec/changes/stage-span-constructions/notes.md`.

## 2. Unit, failure and reproduction accounting

- **Units.** 9,200/9,200 records over 92 cells: no duplicates, no failed records, and all 48 array tasks `COMPLETED 0:0` (`sacct.txt`). Every record carries `orientation_estimator` net displacement, `min_components` 2, `n_jobs` 1, `span_capture` true and the `population_latent` checkpoint. The run used 45.5 recorded core-hours (median 16.9 s per unit) against a projection of 66.
- **Sibling anchors** (`reproduction/siblings_*`, `--by-column`). The stage-span profiles' anchors are the surgery profiles' anchor cells, because every zero-effect cell clears the construction. All three comparisons against the surgery geometry anchor pair 300/300 (column, replicate) records with identical seeds, ranks, statistics and p-values.
- **The design-guidance pilot** (`reproduction/surgery_*_vs_design_guidance_pilot*`, `--by-cell`). This check projects onto `n_samples` and filters the pilot to three stages and two blocks. The surgery siblings reproduce the 2026-10-28 pilot at every shared cell: 800/800 magnitude pairs (8 cells) and 2,000/2,000 geometry pairs (20 cells), with identical seeds, ranks, statistics and p-values. Span capture changed no statistic. The pilot's other columns and effects are listed with zero pairs.

## 3. Visibility: how much of each plant the fitted span sees

Median `group_difference` share over 100 replicates: the share of Σ_s‖d_s‖² inside the fitted span, with d_s the population B − A difference at stage s minus its stage mean. These are energy shares. A random direction keeps rank / 498 = 0.004 at rank 2 and 0.006 at rank 3 (`report/span_capture.csv`).

| mode | e | surgery n = 200 | surgery n = 80 | surgery SEA-AD | stage-span n = 200 | stage-span n = 80 | stage-span SEA-AD |
|---|---|---|---|---|---|---|---|
| magnitude | 0.02 / 0.05 / 0.1 | 0.970 | 0.938–0.939 | 0.934–0.936 | 0.970 | 0.938–0.939 | 0.934–0.936 |
| orientation | 0.25 | 0.002 | 0.003 | 0.005 | 0.971 | 0.935 | 0.934 |
| orientation | 0.5 | 0.002 | 0.003 | 0.008 | 0.969 | 0.928 | 0.924 |
| orientation | 1.0 | 0.004 | 0.006 | **0.215** | 0.946 | **0.849** | **0.850** |
| shape | 0.25 | 0.004 | 0.004 | 0.004 | 0.969 | 0.934 | 0.935 |
| shape | 0.5 | 0.004 | 0.004 | 0.006 | 0.967 | 0.929 | 0.932 |
| shape | 1.0 | 0.010 | 0.008 | 0.011 | 0.964 | 0.917 | 0.926 |

- **The surgeries are invisible.** The orientation surgery relocates stage-changing CpGs, and the shape surgery relocates one stage's. Both put their difference on features the pooled stage program does not use, so it lands outside the span at the random-direction share. The exception is the full relocation on the SEA-AD table (0.215, median rank 4): with unequal cells the pooled means are weighted and the rhombus becomes a kite.
- **The stage-span plants lie in the population span by construction**, and the fitted span keeps 85–97 % of them. Capture falls with n (0.97 at n = 200, 0.93 at n = 80) and with the size of the rotation: a 90° rotation moves B's configuration furthest from the pooled stage means the basis is fitted to.
- **Magnitude** scales the shared program, so both constructions sit inside the span (0.93–0.97). They are the same plant up to the cascade's per-feature pattern.
- **Pooled stage means** keep 0.97 (n = 200) and 0.93–0.94 (n = 80, SEA-AD) of their energy in the span under every construction. Translation has no stage-varying difference, so its `group_difference` share is null in every replicate.

**Seen through the fitted basis** (`population_latent`, median target statistic). The stage-span rotation reads 22.0–22.5° / 43.7–45.3° / 86.2–90.4° against an intended 22.5° / 45° / 90°. The stage-span bend reads `shape` 0.125–0.126 / 0.221–0.224 / 0.359–0.365 against an intended 0.125–0.126 / 0.218–0.221 / 0.346–0.351. The surgery's relocation reads 1.1–1.9° at e = 0.25 and 2.0–3.8° at e = 0.5 at every column. At e = 1.0 it reads 4.4–4.8° at the balanced columns and 53.8° on the SEA-AD table. Its population angle is arccos(1 − e) in the methylation block: 41°, 60° and 90° ([orientation probe](orientation-measurement-probe-2026-10-12.md) §1).

## 4. Sensitivity given visibility: target power

Rejection rate at α = 0.05 over 100 replicates (SE ≤ 0.05; `report/design_point_operating.csv`).

| mode / statistic | e | surgery n = 200 | surgery n = 80 | surgery SEA-AD | stage-span n = 200 | stage-span n = 80 | stage-span SEA-AD |
|---|---|---|---|---|---|---|---|
| orientation / `angle` | 0.25 | 0.17 | 0.05 | 0.09 | 0.85 | 0.43 | 0.32 |
| | 0.5 | 0.22 | 0.06 | 0.11 | 1.00 | 0.99 | 0.95 |
| | 1.0 | 0.12 | 0.15 | 0.67 | 1.00 | 1.00 | 1.00 |
| shape / `shape` | 0.25 | 0.24 | 0.02 | 0.04 | 0.96 | 0.74 | 0.65 |
| | 0.5 | 0.40 | 0.14 | 0.11 | 1.00 | 0.98 | 0.98 |
| | 1.0 | 0.55 | 0.20 | 0.17 | 1.00 | 1.00 | 0.99 |
| magnitude / `delta` | 0.02 | 0.21 | 0.04 | 0.02 | 0.21 | 0.04 | 0.03 |
| | 0.05 | 0.72 | 0.23 | 0.19 | 0.72 | 0.23 | 0.18 |
| | 0.1 | 0.94 | 0.59 | 0.50 | 0.94 | 0.60 | 0.49 |

- **`angle` and `shape` are not weak tests.** Given a difference inside the span they detect a 45° rotation and a half-step bend at n = 80 in 0.98–0.99 of replicates.
- **Their low power under the surgeries is visibility.** The surgeries' orientation power is flat or non-monotone in e at the balanced columns (0.17 / 0.22 / 0.12 at n = 200), the pattern the design-guidance pilot recorded. It follows the small out-of-span leak, not the effect. The SEA-AD table's 0.67 at e = 1.0 is the kite route, at capture 0.215.
- **Magnitude power is the same** under the two constructions: the stage-span − surgery difference is −0.01 to +0.01 at every column and effect, against a pooled SE of 0.022–0.071 (`magnitude_delta_difference` in the reading).

## 5. Specificity: the stage-span modes move their target, the surgeries a mix

Off-target rates (same tables):

- **Stage-span orientation** leaves `delta` at 0.00 and `shape` at 0.00–0.07 at every column and effect.
- **Stage-span shape** leaves `angle` at 0.04–0.12 and `delta` at 0.00–0.08, except at e = 1.0 and n = 200, where `delta` reads 0.17. The full-strength bend (rescale factor ≈ 0.5) adds enough between-cell variance that the sample scaler departs from the null population scale, so the sample path length moves.
- **Stage-span magnitude** moves only `delta`: `angle` 0.03–0.09 and `shape` 0.01–0.05, the same as the surgery's.
- **The surgeries register as a mix**, as [measurement visibility](../measurement-visibility.md) explains. The full relocation at n = 200 rejects `delta` 0.22, `angle` 0.12 and `shape` 0.21. On the SEA-AD table it rejects 0.31 / 0.67 / 0.52. The full shape surgery is detected more often by `delta` than by `shape`: 0.63 against 0.55 at n = 200, and 0.36 against 0.17 on the SEA-AD table.

## 6. Size

**The shared anchor** (`anchor_decision.json`, identical in all four reports, bound 0.0936 at R = 100) reads `delta` / `angle` / `shape`:

| column | `delta` | `angle` | `shape` |
|---|---|---|---|
| n = 200 | 0.02 | 0.10 | 0.05 |
| n = 80 | 0.01 | 0.04 | 0.01 |
| SEA-AD | 0.00 | 0.04 | 0.01 |

`anchor_decision` reads `excess` at n = 200 on `angle` in every profile, because it is the same cell evaluated on the same datasets. These are the design-guidance pilot's baseline anchor datasets at replicates 0–99: the pilot recorded the same 0.10 as one of its two chance-level excesses, and the reproduction above confirms identical p-values. The independent Type I-phase draw at n = 200 (`type_i_table.csv`, the same cell in all four profiles) reads 0.02 / 0.02 / 0.01. The paper-grade slice held the bound at this column at R = 500 (anchor ≤ 0.062).

**At n = 200, every cell that does not target `angle` reads it at 0.07–0.10 under both constructions.** That covers magnitude (0.09 at every effect), translation (0.05–0.10) and stage-span shape (0.09). Those cells share the anchor's noise stream at each replicate index, so an `angle` reading near the bound at this column is a property of these 100 replicates, not of a construction.

**Translation** (the negative control):

| construction | cells | `delta` | `angle` | `shape` |
|---|---|---|---|---|
| stage-span | design-point cells | 0.00–0.02 | 0.03–0.10 | 0.00–0.04 |
| surgery | design-point cells | 0.01–0.05 | 0.02–0.07 | 0.01–0.04 |
| stage-span | Type I phase, e = 0.1 / 1.0 (fresh draws) | 0.01 / 0.01 | 0.03 / 0.08 | 0.01 / 0.00 |
| surgery | Type I phase, e = 0.1 / 1.0 (fresh draws) | 0.05 / 0.01 | 0.03 / 0.04 | 0.02 / 0.01 |

## 7. The reading

`stage_span_pilot_reading.json` (rule version 1, bound α + 2·√(α(1 − α)/R) = 0.0936 at R = 100) returns **`revise_construction`**:

| criterion | observations | failures |
|---|---|---|
| R1: stage-span translation at or below the bound, every column, effect and statistic | 27 | 1: n = 200, e = 0.25, `angle` 0.10 |
| R2: median `group_difference` ≥ 0.90 in every stage-span magnitude, orientation and shape cell | 27 | 2: orientation e = 1.0 at n = 80 (0.849) and SEA-AD (0.850) |
| R3: off-target rates at or below the bound in those cells | 54 | 3: shape e = 1.0 `delta` 0.17 (n = 200), `angle` 0.12 (n = 80) and 0.10 (SEA-AD) |

**How to read the failures:**

- **R2's two failures and R3's `delta` 0.17 are real properties of the largest effects.** A 90° rotation moves B's configuration away from the pooled stage means the basis is fitted to, and capture falls to 0.85 at n = 80. The full-strength bend's between-cell variance shifts the sample scale enough to move `delta` (§5).
- **R1's failure and R3's two `angle` readings at 0.10–0.12 are at chance level.** At R = 100, each rate exceeds the bound with probability 0.028 when its true value is α, so 2.3 of R1's and R3's 81 bound tests are expected to exceed by chance. R1's failure sits at the n = 200 column, whose `angle` reads 0.07–0.10 in every non-target cell (§6).
- **Restricted to e ≤ 0.5:**
  - R2 holds at every column (capture ≥ 0.924).
  - R3 holds at every column (the largest off-target rate is 0.09, stage-span shape `angle` at n = 200).
  - R1's only failure is the inherited n = 200 reading.

**Descriptive tables** in the reading: target power and capture for both constructions (§3–§4), the magnitude `delta` difference with its pooled SE (§4), and each profile's anchor size as read (§6).

## 8. Latent retrofit check

`retrofit/` (`scripts/latent_retrofit_check.py`, the n = 200 column, replicates 0–9, 199 permutations): each replicate's own fitted basis, B's latent stage means edited by a 25 % scaling, a 30° rotation, a half-step bend and a half-step offset, each pulled back by the loadings and the minimum-norm preimage.

- **Frozen basis:** every edit is recovered to 7.8 × 10⁻¹⁴ in the latent coordinates.
- **Refit from scratch:** each edit is detected in its own statistic in 10/10 replicates (`delta`, `angle`, `shape` for magnitude, orientation, shape) and in no other (0/10). The offset is detected in no statistic, except `shape` in 1/10. Median ranks are 2 at the null fit and 2–3 at the refit.
- **Footprints:** 105–170 effective features, 96–97 % of the energy on the truth's stage-active features, 63–67 % on methylation.
- **Agreement with the population plant:** each footprint has cosine 0.965–0.984 with the stage-span construction's population plant for the same mode on the same replicate. The loadings preimage is the lower of the two in every mode.

A basis-defined plant and the population-defined plant the construction uses nearly coincide. The circular alternative (design D1, planting in each replicate's own basis) would therefore have bought little, and it would have aligned the plant with that replicate's noise.

## 9. Proposed paper-grade effect grid, and the revision

The reading is advisory, and the decision belongs to the follow-up change. The evidence supports one revision and one grid:

- **Revision.** Cap the stage-span geometry effects at e ≤ 0.5: rotations of at most 45° and bends of at most half a mean step. This is where R2 and R3 hold at every column. The e = 1.0 failures are the construction leaving the regime it was built for. They are not evidence against the measurement, and the paper should not quote the stage-span 90° rotation as a clean orientation-only plant at n ≤ 80.
- **Grid.**
  - Stage-span `orientation` and `shape` at e ∈ {0.1, 0.25, 0.5}: e = 0.5 already saturates at n ≥ 80, and the informative range is below it.
  - Stage-span `translation` at {0.25, 0.5, 1.0}, as the control.
  - `magnitude` at {0.02, 0.05, 0.1} under one construction only: the two are equivalent here.
  - Surgery siblings at the slice's {0.25, 0.5, 1.0}, as the visibility contrast.
  - All at 500 × 999 in the same seed family, so that the n = 200 `angle` reading near the bound is measured at R = 500 (bound 0.0695), where the slice already held size at that column.
- **Columns.** The pilot's three, plus the three-block and four-stage columns the proposal names, to see how capture moves with the span's dimension.

The **latent plasmode** on the SEA-AD export (a stage-span plant on real noise, which needs a split-sample basis because real data has no population means) is a separate change. So is the paper-grade decision.

## 10. Limitations

- **Pilot grade.** At R = 100, rates near the bound are not resolved (SE 0.022 at α). The n = 200 column's datasets carry an `angle` excess that the slice did not reproduce at R = 500.
- **The stage-span construction is circular by design.** It plants where the measurement looks. R2 is a premise check, not a finding, and the construction answers sensitivity given visibility, not how visible a biological difference is. It also bypasses the methylation → expression → protein cascade; the surgeries remain the biologically grounded constructions.
- **One generator regime.** These results are at δ = 2.0, three stages, two blocks and ρ = 0. The design-guidance arm showed that a weaker stage signal inflates the CV rank, and a larger span may capture more of the surgeries.
- **The SEA-AD column is the case-study table under the design-guidance generator**, not the case-study configuration (no covariates, saturated stage signal).
