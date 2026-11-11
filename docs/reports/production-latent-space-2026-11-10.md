# Production latent-space selection — findings report (2026-11-10)

**Change:** `production-latent-space-selection` (run on MOTCO 0.22.0; no package code changed). **Reads:** `results/production-latent-space-2026-11-10/` (`production_latent_space_decision.json`; `simulation/` and `plasmode/` rate, paired and cell tables, `accounting.json`, `reproduction.json`, `settings.json`; `PROVENANCE.txt`, `sacct.txt`). **Concept page:** [measurement visibility](../measurement-visibility.md#which-stage-supervised-space-weighting-and-sparsity).

Two opt-in settings change the production stage-PLS space without changing its objective:

- **block weighting** (`block_weighting: equal_variance`, MOTCO 0.21.0) multiplies each standardized block by 1/√p_b, so a block's weight no longer grows with its feature count;
- **sparse stage-PLS** (`sparsity: {scope: block}`, MOTCO 0.22.0) lets each component load on a budget of features per block, with the rank and the keep fraction chosen together by the held-out AUROC.

The latent-space comparison ([report](latent-space-comparison-2026-11-04.md)) kept the stage-PLS span as the production measurement, but it left open how much noise that span carries where real data sit: at SEA-AD's stage signal it keeps 0.10 of the pooled stage program. This check asked whether either setting, or both, holds size, loses no in-span power, and gains power where production has none. It compared four arms on matched units under a decision fixed in advance (`design.md` D4).

**In one paragraph.** **The verdict is `keep_production`.** No candidate gains power in the regime the decision reads (the weak design and real SEA-AD noise), so neither setting becomes a production candidate.

- **`weighted` reads `no_gain`.** Size, protected power and specificity hold. Its best gain p-values are 0.105 (the plasmode plant at e = 0.25, `delta` 0.24 → 0.30) and 0.125 (the weak design's 45° rotation, `angle` 0.02 → 0.05), against the threshold 0.05/11 = 0.0045.
- **`sparse_block` and `weighted_sparse_block` read `size_excess`.** On the plasmode anchor, `delta` reads 0.11, above the family's ceiling of 0.10 at R = 100. Production reads 0.11 there too: these 100 replicates hold 11 of the 21 `delta` rejections the 2026-11-02 plasmode found over 500. Both arms also fail the gain criterion. `sparse_block`'s best gain is the plasmode plant at e = 0.5, `delta` 0.51 → 0.62, p = 0.0096.
- **Where the candidates do gain, the decision does not look.** At the δ = 2.0 designs every candidate rejects more often than production where production is below saturation. The largest gain is the SEA-AD table's 22.5° rotation: `angle` 0.46 → 0.90 (`weighted`), 1.00 (`sparse_block`) and 0.99 (`weighted_sparse_block`). The candidates pick rank 2 where production's CV picks rank 3. These cells are protected-power readings (S2), so the gain counts only as "no loss".

## 1. Configuration and provenance

- **The unit.** Each unit builds one dataset and evaluates four arms on it, serially, at the unit's RRPP seed. Every arm is its own stage-PLS fit (double-CV rank, floor 2, net displacement, span capture):
  - **`pls_cv`:** production;
  - **`weighted`:** `block_weighting: equal_variance` (every PLS fit with `scale=False`, so `delta` is in the weighted metric);
  - **`sparse_block`:** `sparsity: {scope: block}` with the default ladder 1.0 / 0.5 / 0.2 / 0.1 / 0.05;
  - **`weighted_sparse_block`:** both.

  The arms share the data and the permutation indices, so the decision compares them by exact paired tests on the discordant units.
- **`simulation` profile.** The latent-space comparison's three designs, seed functions and committed weak δ (base seed 2035), at 100 replicates × 199 permutations, rejection p < 0.05:
  - **Designs:** `n200` (n = 200, δ = 2.0), `seaad` (the SEA-AD table, δ = 2.0) and `seaad_weak` (the SEA-AD table at δ = 0.11, held-out AUROC about 0.6).
  - **Cells (10):** the comparison's seven in-span cells by name (anchor; translation; magnitude e = 0.05; orientation e = 0.25 and 0.5 at φ = 0; shape e = 0.25 and 0.5), whose datasets and draws are the comparison's, plus `magnitude_e0.25`, `orient_e1.0_phi0` and `shape_e1.0` for a power gradient. Orientation e is in units of 90°, so e = 0.25 is a 22.5° rotation.
  - **Units:** 3,000.
- **`plasmode` profile.** The SEA-AD astrocyte `cohort_30` export (RNA and ATAC, 5,000 features each, 80 donors), read through the plasmode's contract with its file hashes checked. Each replicate permutes the sex labels within stage as the plasmode's `paper` stream does (base seed 2034, replicates 0–99). Each arm runs `analyze_trajectories` at the case-study configuration: six covariates, floor 2, net displacement, `df`, attribution off, 999 permutations, rejection p ≤ 0.05.
  - **Cells (4):** the anchor and a real-noise span plant at e = 0.25, 0.5 and 1.0: the donors labelled `Male` at stage s gain e · (x̄_s − x̄), the pooled stage-s mean minus the mean of the stage means. B's centred stage configuration thus gains e times the pooled one, a stage-span magnitude edit on real noise.
  - **Units:** 400.
- **Code.** Revision `d0a5ddb` on `proposal/production-latent-space` (MOTCO 0.22.0), run from a clean detached worktree on `ing`. Python 3.11.16, numpy 2.3.5, scipy 1.16.3, scikit-learn 1.8.0 and pandas 2.3.3, with BLAS pinned to one thread.
- **Cluster.** Partition `512x1024`, n11 excluded: array 937293 (`simulation`, 300 shards of 10 units, at most 96 at once) and array 937294 (`plasmode`, 100 shards of 4 units, at most 64 at once).
- **Pre-flight.** A workstation smoke of both profiles (68 units, 0 failed records) projected about 490 core-hours. The user approved the full ladder at 100 replicates with a cap of 650 core-hours (`openspec/changes/production-latent-space-selection/notes.md`).
- **Amendments before any unit** (2026-11-10, `notes.md`): S3 reads the e ≤ 0.5 cells only, and S1 and S3 use the stage-span paper grade's count-plus-ceiling rule.

Details are in `PROVENANCE.txt`.

## 2. Accounting and reproduction

- **Units.** `simulation`: 3,000 of 3,000 units, 12,000 records. `plasmode`: 400 of 400, 1,600 records. Both: 0 failed records, 0 duplicates, 0 fidelity failures. All 400 array tasks ended `COMPLETED 0:0` (`sacct.txt`).
- **Cost.** 184.0 + 325.0 = **509.0 recorded core-hours**, against the projection of 490 and the cap of 650.
  - **`simulation`:** median unit 208.5 s; per arm, `pls_cv` 13.2 s, `weighted` 13.2 s, `sparse_block` 91.5 s, `weighted_sparse_block` 90.7 s.
  - **`plasmode`:** median unit 2,852 s; per arm, 133.8 s, 133.7 s, 1,261 s and 1,268 s.

  The sparse arms are 87 % (simulation) and 91 % (plasmode) of the cost: every (rank, fraction) candidate is refitted in every fold.
- **Reproduction (the decision's gate): both passed.**
  - **`simulation`:** the `pls_cv` arm's count, rejections and failures equal the comparison's committed `rate_table.csv` in all 63 rows (3 designs × 7 shared cells × 3 statistics).
  - **`plasmode`:** every one of the 100 anchor units' `pls_cv` arm selected rank 5 at held-out AUROC 0.635556, the plasmode's fidelity values.
  - **So** every `pls_cv` rate on the shared cells is the comparison's rate, and every candidate rate pairs with it on the same datasets and draws.

## 3. The decision

`production_latent_space_decision.json` (rule version 1) returns **`keep_production`**: no arm qualifies.

| criterion | rule | `weighted` | `sparse_block` | `weighted_sparse_block` |
|---|---|---|---|---|
| gate | both reproductions pass | **passed** | | |
| S1 size, simulation | anchor and translation `delta`/`angle`/`shape` at three designs, m = 18: at most 2 above the bound 0.0936, none above the ceiling 0.12 | holds: 0 above (0.00–0.05) | holds: 0 above (0.01–0.06) | holds: 0 above (0.01–0.06) |
| S1 size, plasmode | the anchor's three statistics, m = 3: at most 1 above 0.0936, none above the ceiling 0.10 | holds: 0.05 / 0.05 / 0.04 | **fails**: `delta` 0.11 above the ceiling | **fails**: `delta` 0.11 above the ceiling |
| S2 protected power | 16 target cells at δ = 2.0, one-sided exact paired loss p: at most 2 below 0.025, none below 0.05/16 | holds: 0 losses | holds: 0 losses | holds: 0 losses |
| S3 specificity | off-target rates of the δ = 2.0 target cells at e ≤ 0.5, m = 24, under the S1 rule | holds: 0 above (0.00–0.08) | holds: 0 above (0.00–0.07) | holds: 0 above (0.00–0.07) |
| G gain | 11 cells (the weak design's 8 targets, the plasmode's 3 plant `delta` cells), one exact paired gain p ≤ 0.05/11 = 0.0045 | **fails**: best p 0.105 | **fails**: best p 0.0096 | **fails**: best p 0.125 |
| **arm verdict** | `size_excess` → `power_loss` → `specificity_excess` → `no_gain` → `qualifies` | **`no_gain`** | **`size_excess`** | **`size_excess`** |

**Consequence** (the decision's text): "Production keeps the stage-PLS space unchanged; block weighting and sparse stage-PLS stay opt-in and are documented as measured."

The `pls_cv` arm's own controls read 0.00–0.06 at the 18 simulation readings and 0.11 / 0.06 / 0.04 (`delta` / `angle` / `shape`) at the plasmode anchor. Its `delta` there is above the family's ceiling too (§4).

## 4. Size, and the plasmode anchor's 0.11

Every arm holds size in simulation: the 72 control readings (18 per arm) run from 0.00 to 0.06, and none is above the bound.

On the plasmode anchor, `delta` reads:

| arm | rate | b / c against `pls_cv` |
|---|---|---|
| `pls_cv` | 0.11 | — |
| `weighted` | 0.05 | 6 / 0 |
| `sparse_block` | 0.11 | 6 / 6 |
| `weighted_sparse_block` | 0.11 | 6 / 6 |

(b counts units only `pls_cv` rejects, c units only the candidate rejects; `plasmode/paired_table.csv`.)

- **Production reads 0.11 on these replicates.** The `pls_cv` arm's anchor p-values equal the 2026-11-02 plasmode's per-replicate records at all 100 replicates, with equal seeds (checked after the decision; `PROVENANCE.txt`). That run's replicates 0–99 hold 11 of its 21 `delta` rejections over 500 (0.042 over 500, [plasmode report](plasmode-2026-11-02.md)). The 0.11 is the replicate subset, not an excess the sparse fit adds.
- **The sparse arms share 5 of their 11 rejections with production.** Their anchor space is one label-free fit, the same in every unit (rank 4 at keep fraction 0.05, 250 features per block; `plasmode/cell_table.csv`), so they reject on different replicates than production at the same rate.
- **The verdict does not depend on it.** Without S1 the sparse arms would read `no_gain`, since every arm fails G. The S1 failure ties production on these permutations; it does not measure an excess of the sparse fit. At R = 100 a rate of 0.11 is about 2.8 SE above α, and the paper grade the decision would have opened is where size is resolved.
- **`weighted` reads 0.05.** Its rejections are a subset of production's (b = 6, c = 0). On the plasmode the weighting is not a rebalancing (§7), so this is the effect of fitting without sklearn's per-fold re-standardization, which also moves the anchor's CV rank from 5 to 7.

## 5. Gain: the regime the decision reads

G reads the cells where production rejects at or near the null rate (the weak design) and the plasmode plant. Each entry is `pls_cv` → candidate over 100 paired units, with b/c in brackets and the one-sided gain p after it.

**Weak design (SEA-AD table, δ = 0.11):**

| cell | statistic | `weighted` | `sparse_block` | `weighted_sparse_block` |
|---|---|---|---|---|
| magnitude e = 0.05 | `delta` | 0.05 → 0.04 (2/1) 0.875 | 0.05 → 0.05 (0/0) 1 | 0.05 → 0.07 (2/4) 0.344 |
| magnitude e = 0.25 | `delta` | 0.06 → 0.05 (4/3) 0.773 | 0.06 → 0.07 (1/2) 0.5 | 0.06 → 0.06 (4/4) 0.637 |
| orientation e = 0.25 | `angle` | 0.01 → 0.02 (0/1) 0.5 | 0.01 → 0.02 (0/1) 0.5 | 0.01 → 0.03 (0/2) 0.25 |
| orientation e = 0.5 | `angle` | 0.02 → 0.05 (0/3) 0.125 | 0.02 → 0.02 (0/0) 1 | 0.02 → 0.05 (0/3) 0.125 |
| orientation e = 1.0 | `angle` | 0.03 → 0.02 (2/1) 0.875 | 0.03 → 0.04 (0/1) 0.5 | 0.03 → 0.04 (2/3) 0.5 |
| shape e = 0.25 | `shape` | 0.05 → 0.04 (3/2) 0.813 | 0.05 → 0.05 (1/1) 0.75 | 0.05 → 0.01 (4/0) 1 |
| shape e = 0.5 | `shape` | 0.03 → 0.06 (1/4) 0.188 | 0.03 → 0.05 (1/3) 0.313 | 0.03 → 0.06 (1/4) 0.188 |
| shape e = 1.0 | `shape` | 0.02 → 0.04 (1/3) 0.313 | 0.02 → 0.04 (0/2) 0.25 | 0.02 → 0.05 (1/4) 0.188 |

**Plasmode span plant (`delta`):**

| e | `weighted` | `sparse_block` | `weighted_sparse_block` |
|---|---|---|---|
| 0.25 | 0.24 → 0.30 (5/11) 0.105 | 0.24 → 0.28 (4/8) 0.194 | 0.24 → 0.21 (11/8) 0.820 |
| 0.5 | 0.51 → 0.52 (5/6) 0.5 | 0.51 → 0.62 (4/15) **0.0096** | 0.51 → 0.49 (11/9) 0.748 |
| 1.0 | 1.00 → 1.00 (0/0) 1 | 1.00 → 0.99 (1/0) 1 | 1.00 → 0.99 (1/0) 1 |

- **The weak design stays at the null rate in every arm.** Production reads 0.01–0.06 at its eight target cells, and no candidate moves a cell by more than 0.03. Pooled over the 11 cells the candidates gain 0.019 (`sparse_block`, b/c 12/33), 0.011 (`weighted`, 23/35) and 0.004 (`weighted_sparse_block`, 37/41).
- **The one sizeable gain is `sparse_block` on the plasmode plant at e = 0.5** (0.51 → 0.62, p = 0.0096). It misses the Bonferroni threshold by a factor of two, and the arm fails S1 regardless.
- **The capture does not rise.** At the weak design's anchor the fitted span keeps a median 0.100 of the pooled stage program under `pls_cv`, 0.098 under `weighted`, 0.085 under `sparse_block` and 0.082 under `weighted_sparse_block` (`pooled_stage`, `simulation/cell_table.csv`). Over the design's ten cells the medians run 0.054–0.123 for production and 0.043–0.128 for the candidates, and the held-out stage AUROC 0.548–0.630 and 0.538–0.631. Neither setting finds more of the stage program in the weak regime, which is the mechanism G was meant to reward.

## 6. Protected power, and where the candidates gain

S2 reads the 16 target cells at the δ = 2.0 designs for losses. There are none: in every candidate and cell c ≥ b. Where production is below saturation the candidates gain:

| design | cell | statistic | `pls_cv` | `weighted` (b/c) | `sparse_block` (b/c) | `weighted_sparse_block` (b/c) |
|---|---|---|---|---|---|---|
| n = 200 | magnitude e = 0.05 | `delta` | 0.74 | 0.93 (3/22) | 0.81 (0/7) | 0.93 (3/22) |
| n = 200 | orientation e = 0.25 | `angle` | 0.84 | 1.00 (0/16) | 1.00 (0/16) | 1.00 (0/16) |
| SEA-AD table | magnitude e = 0.05 | `delta` | 0.23 | 0.31 (7/15) | 0.39 (4/20) | 0.31 (7/15) |
| SEA-AD table | orientation e = 0.25 | `angle` | 0.46 | 0.90 (1/45) | 1.00 (0/54) | 0.99 (1/54) |
| SEA-AD table | orientation e = 0.5 | `angle` | 0.96 | 1.00 (0/4) | 1.00 (0/4) | 1.00 (0/4) |
| SEA-AD table | shape e = 0.25 | `shape` | 0.70 | 0.87 (10/27) | 0.92 (1/23) | 0.88 (10/28) |
| SEA-AD table | shape e = 0.5 | `shape` | 0.95 | 0.98 (2/5) | 0.99 (0/4) | 0.99 (1/5) |
| SEA-AD table | magnitude e = 0.25 | `delta` | 0.95 | 0.99 (1/5) | 0.97 (0/2) | 0.99 (1/5) |

The other eight cells are at 1.00 in every arm. The two-sided paired p-value is below 0.05 for 13 of these 24 comparisons, and down to 1e-16 for the SEA-AD rotation (`descriptive.paired_differences`).

- **The gains follow the selected rank.** At these cells production's CV picks rank 3 in a large share of units (median 3 on the SEA-AD rotation and bend, q90 3 at n = 200), and every candidate picks rank 2 (median 2, q90 2–2.1 for `weighted`, 2 for the sparse arms; `simulation/cell_table.csv`). The capture of the pooled stage program is 0.92–0.98 in every arm there, so the candidates do not see more of the program; they measure it in fewer dimensions. The design-guidance slice found the same association between `angle` and the CV rank ([report](design-guidance-slice-2026-10-29.md) §8).
- **`weighted` and `weighted_sparse_block` often coincide in simulation.** Under weighting the sparse CV keeps the full budget in most units (median keep fraction 1.0 at every δ = 2.0 target cell except the 90° rotation at n = 200 and the SEA-AD 45° rotation), so the two weighted arms often fit the same space and read the same b/c at five of the eight cells above.
- **None of this enters the verdict.** The decision protects these cells against loss and looks for gains only where production has none. A gain at a stronger stage signal was not a predeclared reason to change production, and at R = 100 the paired gains are not size-checked at a matched paper grade.

## 7. Specificity, and what the plasmode's `weighted` arm measures

**S3** (off-target rates of the δ = 2.0 target cells at e ≤ 0.5) holds for every arm with no reading above the bound 0.0936. `pls_cv`'s own off-target rates there read 0.00–0.07. One off-target reading falls in all three candidates: the n = 200 rotation at e = 0.5 leaks into `shape` in 0.06 of production's units and in none of the candidates' (6/0, two-sided p 0.031).

**The e = 1.0 cells' off-target rates** are descriptive (`descriptive.large_plant_off_target`). Production's construction leak there is the stage-span pilot's: the n = 200 full bend's `delta` reads 0.17 (`pls_cv`), 0.11 (`weighted`), 0.15 (`sparse_block`) and 0.09 (`weighted_sparse_block`), and the SEA-AD bend's `angle` 0.12, 0.10, 0.08 and 0.10. No candidate adds a leak there.

**At the weak design** no off-target rate is read by a criterion. The largest is `weighted`'s `shape` under the 90° rotation, 0.10 against production's 0.03 (1/8, two-sided p 0.039); `weighted_sparse_block` reads 0.08 there.

**On the plasmode the weighting is not a rebalancing.** The export's blocks have 5,000 features each, so `equal_variance` gives them equal weights. The `weighted` arm is therefore production fitted with `scale=False`, without sklearn's per-fold re-standardization. That alone changes the label-free anchor space (rank 5 → 7, held-out AUROC 0.635556 → 0.638444). The plasmode's `weighted` readings measure that change, not block weighting; the simulation's blocks (367 and 131 features) carry the weighting question.

## 8. Rank, keep fraction and support

From `simulation/cell_table.csv` and `plasmode/cell_table.csv` (medians, q10–q90 in brackets):

- **δ = 2.0 designs.** Production's median rank is 2–3, the candidates' 2. `sparse_block` keeps half the features (median 0.5) at the SEA-AD table's rotations and bends, at both translations, and at n = 200's 45° rotation and full bend; elsewhere it keeps the full budget (median 1.0, q10 0.5 at most cells). AUROC is 1.000 in every arm.
- **Weak design.** Production's median rank is 2–3 (q90 5–6); `sparse_block`'s 2–3 (q90 4–5) and `weighted_sparse_block`'s 2. Both sparse arms keep a median 1.0 of the features with q10 0.05–0.1, so a tenth of the units drop to the smallest budgets.
- **Plasmode anchor.** The label-free space is identical in every unit: production rank 5, `weighted` 7, both sparse arms rank 4 at keep fraction 0.05 (250 features per block).
- **Plasmode plant.** Production's rank rises to a median 7 at e = 0.25 and 0.5 (q90 10.1 and 9) and falls back to 5 at e = 1.0. `sparse_block` keeps 0.2 of the features at e = 0.25 (1,000 per block) and 0.5 at e = 0.5 and 1.0 (2,500 per block); `weighted_sparse_block` keeps 0.5 at every effect.

## 9. What follows

- **Production and `motco analyze` keep the stage-PLS space unchanged.** The paper grade and the default flip were owed only to `propose_paper_grade`. Neither is opened.
- **Both settings stay opt-in and are now measured.**
  - **`block_weighting: equal_variance`** holds size, loses no protected power and adds no specificity excess on these designs. It gains no power at a weak stage signal or on the plasmode plant. A case study may declare it as a sensitivity analysis with this reading. Under it `delta` is in the weighted metric.
  - **`sparsity: {scope: block}`** holds size in simulation and reads the plasmode anchor's `delta` at production's rate on the same replicates; it loses no protected power. It gains no Bonferroni-protected power where the decision looks, and it costs about seven to nine production fits. The plasmode anchor's size is unresolved at R = 100 for every arm, production included.
- **The gain at δ = 2.0 is an open thread, not a finding for production.** Every candidate gains up to 0.54 in `angle` where production's CV selects rank 3 and the candidates select rank 2. That points at the rank rule, not at weighting or sparsity. The latent-rank ladder (2026-09-08) kept the CV rank at the Phase 5 design; whether a stricter rank rule recovers these gains at held size, and what it costs at the weak design, is a question for its own change.
- **Thread 2 of roadmap item 16** (the group-stratified stage space) was sequenced after this selection because it is paired against production. Production is unchanged, so it pairs against `pls_cv`.

## 10. Limitations

- **R = 100.** Every rate has an SE of about 0.02 at α, and the bound is 0.0936. The check decides whether a paper grade is worth its cost; it does not resolve rates near the bound.
- **One kind of each setting.** `mfa` weighting and joint-scope sparsity were not run (proposal, non-goals).
- **The plasmode's weighted arms do not weight.** With equal block sizes `equal_variance` is uniform (§7).
- **Weighted `delta` is in another metric.** Only p-values and rejections are compared across arms (`descriptive.delta_metric_note`).
- **Convergence.** The smoke traced rare `ConvergenceWarning`s in the weighted sparse arm at the smallest budget (17 of 45,360 power-method calls in one unit; `notes.md`). They were not counted in the run and do not change a selection measurably in the smoke.
- **The δ = 2.0 gains are paired at 199 permutations on the comparison's datasets.** They were not part of a predeclared criterion and are reported as read.
