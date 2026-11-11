# Interaction scan by block — findings report (2026-11-10)

**Change:** `interaction-scan-block-maxima` (MOTCO 0.23.0). **Reads:** `results/interaction-scan-blocks-2026-11-10/` (`interaction_scan_blocks_decision.json`, the `simulation/`, `plasmode/` and `plasmode_unequal/` scan, block, false-hit, paired, minP and allocation tables, `accounting.json`, `reproduction.json`, `settings.json`, `PROVENANCE.txt`, `sacct.txt`). **Parent:** the [interaction-scan validation](interaction-scan-2026-11-09.md). **Concept page:** [measurement visibility](../measurement-visibility.md#testing-feature-by-feature-the-interaction-scan).

The per-feature interaction scan (MOTCO 0.20.0) takes one maximum F over every feature, so its familywise bar is set by whichever block tends to give the largest null F. On the SEA-AD export each block holds half the features, but ATAC gave 0.618 of the null maxima ([parent](interaction-scan-2026-11-09.md), §5). With blocks of unequal size the larger block would set the bar almost alone.

The by-block combination (`interaction_feature_scan(..., blocks=)`, `motco analyze --interaction-scan-by-block`) runs the scan's step-down max-T within each block, on the same draws, and splits the familywise error equally over the K blocks that have a non-degenerate feature (Bonferroni). A feature's adjusted p-value is K times its within-block adjusted p-value, capped at 1. The global p-value is K times the smallest block global p-value, capped at 1. A hit is a feature with adjusted p ≤ 0.05, so at K = 2 each block is tested at 0.025.

**In one paragraph.** **The by-block scan is validated: `by_block_validated`.** It holds size and familywise control, it keeps the joint scan's global power on both exports, and on an export with a small block it names more of the small block's planted features.

- **Size holds everywhere.** The six simulation control rates read 0.046–0.058, against a per-cell bound of 0.0695. On real noise the global rates read 0.040 on the SEA-AD export and 0.044 on the unequal export.
- **Familywise control holds under the plasmode's partial null.** A feature outside the planted set was a hit in 0.040 (SEA-AD export) and 0.044 (unequal export) of units, at both planted strengths.
- **No global power is lost.** On the SEA-AD export both arms reject in 0.692 of units at e = 2.0 and in every unit at e = 3.0. On the unequal export the by-block scan rejects in 0.892 at e = 2.0, against the joint scan's 0.760.
- **The small block gains.** On the unequal export (RNA cut to 500 of 5,000 features) the by-block arm finds more planted RNA features in 296 of 500 units at e = 2.0 and fewer in none. RNA recall rises from 0.069 to 0.158 at e = 2.0 and from 0.570 to 0.755 at e = 3.0.
- **The large block pays.** ATAC recall falls from 0.076 to 0.058 at e = 2.0 and from 0.575 to 0.500 at e = 3.0 on the same export. On the SEA-AD export, with equal blocks, detection moves from ATAC to RNA at the same total.
- **Bonferroni costs little at K = 2.** The block maxima correlate weakly (median Spearman 0.18–0.20 on the exports, at most 0.05 in the simulation), and a minP combination would have rejected at most 2 more units out of 500 in any cell.

## 1. Configuration and provenance

- **Arms.** Every unit runs two scans on the same inputs and draws: `joint` (the 0.20.0 scan, no labels) and `by_block` (labels from `block__feature`). Both reject at adjusted p ≤ 0.05.
- **Scan-only units.** No unit fits a latent space. The scan reads the unweighted pooled-standardized features and the covariate-adjusted full and reduced models, built by the analysis's own input helper.
- **Simulation profile.** The interaction-scan validation's 6,600 simulation units: the latent-space comparison's three designs (`n200`, `seaad`, and `seaad_weak` at δ = 0.11), methylation + expression (498 features, 0.263 of them expression), fourteen cells, base seed 2035. The anchor and translation ran 500 replicates, every other cell 100, all at 999 permutations. The settings were rebuilt through `interaction_scan_check.simulation_settings` and equal the parent's committed ones.
- **Plasmode profile.** The interaction-scan validation's 1,500 SEA-AD plasmode units: the astrocyte `cohort_30` export (5,000 RNA + 5,000 ATAC features) under the plasmode's contract, within-stage label permutation, the 20-feature plant (10 per block), base seed 2034, `paper` stream, the case study's six covariates and `df`. Cells: the anchor, e = 2.0 and e = 3.0 at 500 × 999.
- **Unequal profile (`plasmode_unequal`).** The same export and units with RNA reduced to 500 of its 5,000 features: the first 500 positions of a seeded permutation (seed 455583459), kept in export order and recorded by the sha256 of the feature list. The plant still puts 10 features in each block. This is the unequal-size case the change is for, and it has no parent.
- **Records** hold seeds, statistics and counts only.
- **Code:** `fe87330` on `proposal/interaction-scan-block-maxima`, in clean worktrees on the workstation and on `ing`.
- **Compute.** The simulation profile ran on the workstation (AMD Ryzen 7 5800X, 8 workers). The two plasmode profiles ran on `ing` as arrays 937733 and 937734 (10 shards each, all on n4, n11 excluded, BLAS pinned to one thread), queued behind the production latent-space arrays. Job 937821 merged both and ran `decide`.
- **Pre-flight:** a local smoke of all three profiles at 2 replicates per cell × 99 permutations, a serial timing at 999 permutations (projection 2.0 core-hours), and the change's test modules on `ing` (164 passed).

Details are in `PROVENANCE.txt` and `openspec/changes/interaction-scan-block-maxima/notes.md`.

## 2. Accounting and reproduction

- **Units:** 9,600 of 9,600 complete, 0 failed, 0 duplicates (6,600 simulation, 1,500 per plasmode profile); all 20 array tasks and the merge job completed.
- **Compute:** 2.89 recorded core-hours (simulation 1.57, plasmode 0.94, unequal 0.38), against the projection of 2.0. Unit medians were 0.81 s (simulation), 2.23 s (plasmode) and 0.91 s (unequal). The parent run, which fitted a latent space in every unit, took 69.6.
- **Reproduction.** The joint arm reproduced the interaction-scan validation exactly at the written precision, under p ≤ 0.05:
  - **Simulation:** 168 rows, the scan table's count, rejections and median hits, the block balance's mean hits and null-maximum share, and the population-null false-hit count and rejections, for every design and cell.
  - **Plasmode:** 11 rows, the same scan and block-balance columns, and the recovery table's false-hit count and rejections, recall, hit precision and top-20 precision.

  The by-block arm therefore pairs unit by unit with the joint scan of the parent run, and through it with the split test and the plasmode.

## 3. The decision

`interaction_scan_blocks_decision.json`, rule version 1. S1–S3 use the stage-span paper grade's count-plus-ceiling rule. At R = 500 the per-cell bound is 0.0695, and the exceedance probability per rate is 0.030.

| Criterion | What | m | By-block rates | Exceedances (allowed) | Above ceiling | Holds |
|---|---|---|---|---|---|---|
| S1 | simulation size: global by-block test at the anchor and translation of every design | 6 | 0.046–0.058 | 0 (1) | 0 (ceiling 0.074) | yes |
| S2 | real-noise size: global by-block test at both exports' anchors | 2 | 0.040, 0.044 | 0 (1) | 0 (ceiling 0.070) | yes |
| S3 | strong familywise control: units with a hit outside the planted set at e = 2.0 and 3.0 on both exports | 4 | 0.040, 0.040, 0.044, 0.044 | 0 (1) | 0 (ceiling 0.072) | yes |
| P | protected power: no exact one-sided paired loss p-value below 0.05/4 at e = 2.0 and 3.0 on both exports | 4 | smallest p_loss 0.593 | — | — | yes |
| G | allocation gain: more planted reduced-block (RNA) hits by block, exact sign test at p ≤ 0.025 at e = 2.0 or 3.0 | — | p = 7.9 × 10⁻⁹⁰ and 7.2 × 10⁻¹³⁰ | — | — | yes |

**S1.** The six control rates, with the joint arm's on the same units beside them:

| Design | Anchor (by block / joint) | Translation (by block / joint) |
|---|---|---|
| `n200` | 0.052 / 0.056 | 0.052 / 0.052 |
| `seaad` | 0.046 / 0.054 | 0.048 / 0.048 |
| `seaad_weak` | 0.050 / 0.056 | 0.058 / 0.058 |

**S2 and S3.**

| Export | Anchor (by block / joint) | Hit outside the planted set, e = 2.0 (by block / joint) | e = 3.0 (by block / joint) |
|---|---|---|---|
| SEA-AD (`plasmode`) | 0.040 / 0.042 | 0.040 / 0.042 | 0.040 / 0.042 |
| unequal (`plasmode_unequal`) | 0.044 / 0.056 | 0.044 / 0.056 | 0.044 / 0.056 |

In each arm the count of units with an outside hit at e = 2.0 and 3.0 equals the anchor's count of global rejections (20 and 22 by block, 21 and 28 joint). The parent run found the same pattern and checked on its per-unit records that these are the same replicates: the labels and RRPP seed are shared across a replicate's cells and the plant changes only 20 features. Here only the counts are compared.

**P.** The paired global rejections (`paired_table.csv`):

| Export | Cell | Both | Joint only | By block only | Neither | p_loss |
|---|---|---|---|---|---|---|
| SEA-AD | e = 2.0 | 337 | 9 | 9 | 145 | 0.593 |
| SEA-AD | e = 3.0 | 500 | 0 | 0 | 0 | 1.000 |
| unequal | e = 2.0 | 372 | 8 | 74 | 46 | 1.000 |
| unequal | e = 3.0 | 500 | 0 | 0 | 0 | 1.000 |

On the unequal export at e = 2.0 the by-block arm gains 74 units and loses 8 (descriptive gain p-value 8.2 × 10⁻¹⁵). The global rate rises from 0.760 to 0.892.

**G.** On the unequal export, per unit, the by-block arm's planted RNA hits minus the joint arm's (`plasmode_unequal/allocation_table.csv`):

| Cell | Gains | Losses | Ties | Mean difference | p_sign |
|---|---|---|---|---|---|
| e = 2.0 | 296 | 0 | 204 | +0.894 | 7.9 × 10⁻⁹⁰ |
| e = 3.0 | 429 | 0 | 71 | +1.854 | 7.2 × 10⁻¹³⁰ |

**Verdict: `by_block_validated`.** The decision's consequence: the by-block scan is validated for exports whose blocks differ in size or noise. It controls the familywise error over all scanned features at 0.05, split equally across the blocks, keeps the joint scan's global power and names more features of the small block. Only two blocks were measured. The joint scan stays the default, and the by-block combination stays opt-in.

## 4. Where the hits go: recall per block

Recall is the planted features found in a block divided by its 10. Hit precision is the planted share of the block's hits, over units with a hit in the block (`block_table.csv`).

| Export | Cell | Block | Recall (joint → by block) | Mean hits (joint → by block) | Hit precision (joint / by block) |
|---|---|---|---|---|---|
| SEA-AD | e = 2.0 | ATAC | 0.063 → 0.058 | 0.67 → 0.61 | 0.948 / 0.951 |
| SEA-AD | e = 2.0 | RNA | 0.056 → 0.060 | 0.58 → 0.63 | 0.979 / 0.978 |
| SEA-AD | e = 3.0 | ATAC | 0.523 → 0.500 | 5.27 → 5.04 | 0.994 / 0.994 |
| SEA-AD | e = 3.0 | RNA | 0.527 → 0.547 | 5.30 → 5.50 | 0.997 / 0.997 |
| unequal | e = 2.0 | ATAC | 0.076 → 0.058 | 0.83 → 0.61 | 0.932 / 0.951 |
| unequal | e = 2.0 | RNA | 0.069 → 0.158 | 0.69 → 1.60 | 0.995 / 0.989 |
| unequal | e = 3.0 | ATAC | 0.575 → 0.500 | 5.82 → 5.04 | 0.991 / 0.994 |
| unequal | e = 3.0 | RNA | 0.570 → 0.755 | 5.70 → 7.57 | 0.999 / 0.997 |

- **On the SEA-AD export** the blocks are equal in size, so the change only redistributes. The by-block arm finds slightly fewer ATAC and slightly more RNA features. The total is unchanged: the mean of the two blocks' recall is 0.059 in both arms at e = 2.0, and 0.525 against 0.524 at e = 3.0. The per-unit allocation table agrees: RNA gains in 22 units and loses in none at e = 2.0, 89 against 6 at e = 3.0, while ATAC loses in 25 and 99 (`plasmode/allocation_table.csv`). The joint arm's totals equal the parent's committed recall (0.0594 and 0.525).
- **On the unequal export** the small block gains much more than the large block loses. The mean of the two blocks' recall rises from 0.07 to 0.11 at e = 2.0 and from 0.57 to 0.63 at e = 3.0. Under the joint scan, a 500-feature RNA block faced a bar set by 5,500 features, 0.92 of whose null maxima came from ATAC. By block, RNA faces its own bar at 0.025.
- **Hit precision stays at 0.93 or above** in both arms, in every block and cell.
- **The by-block ATAC rows are the same on both exports** in every recall, hit, precision and false-hit column. That is as it should be: the within-block step-down reads only the block's own features and the shared draws, K is 2 on both exports, and ATAC and its plant are the same.

## 5. Block balance

The joint arm's share of null maxima per block, against its share of features, and the rate of units with a false hit in each block at the anchors, per arm (`block_table.csv`). Under the by-block scan each block's false-hit rate should be near α/K = 0.025.

| Profile | Design | Block | Share of features | Share of null maxima (joint) | False-hit units, joint | False-hit units, by block |
|---|---|---|---|---|---|---|
| plasmode | — | ATAC | 0.500 | 0.618 | 0.034 | 0.030 |
| plasmode | — | RNA | 0.500 | 0.382 | 0.014 | 0.016 |
| plasmode_unequal | — | ATAC | 0.909 | 0.920 | 0.052 | 0.030 |
| plasmode_unequal | — | RNA | 0.091 | 0.080 | 0.004 | 0.016 |
| simulation | `n200` | methylation | 0.737 | 0.501 | 0.038 | 0.036 |
| simulation | `n200` | expression | 0.263 | 0.499 | 0.018 | 0.016 |
| simulation | `seaad` | methylation | 0.737 | 0.506 | 0.030 | 0.022 |
| simulation | `seaad` | expression | 0.263 | 0.494 | 0.024 | 0.024 |
| simulation | `seaad_weak` | methylation | 0.737 | 0.506 | 0.036 | 0.028 |
| simulation | `seaad_weak` | expression | 0.263 | 0.494 | 0.022 | 0.024 |

- **The joint scan spends its error unevenly.** On the unequal export nearly all of its false hits fall in ATAC (0.052 against RNA's 0.004). On the SEA-AD export ATAC takes 0.034 against RNA's 0.014.
- **By block, the error is shared.** Every by-block rate at the plasmode anchors is 0.016 or 0.030, on both exports.
- **In the simulation** the joint shares were already near even, since correlated methylation features act as fewer tests (methylation's share of null maxima is 0.497–0.523 in every cell). There the two arms differ little.

## 6. What dependence between blocks would have bought: minP

The decision rule uses Bonferroni across blocks, which ignores dependence between the block maxima. The run measured a minP combination descriptively: each block maximum's tail probability under the block's own null maxima, combined by the smallest (`minp_table.csv`).

| Profile | Cell | Joint | Bonferroni (by block) | minP | Median Spearman of block maxima |
|---|---|---|---|---|---|
| plasmode | anchor | 0.042 | 0.040 | 0.044 | 0.196 |
| plasmode | e = 2.0 | 0.692 | 0.692 | 0.694 | 0.196 |
| plasmode | e = 3.0 | 1.000 | 1.000 | 1.000 | 0.196 |
| plasmode_unequal | anchor | 0.056 | 0.044 | 0.048 | 0.184 |
| plasmode_unequal | e = 2.0 | 0.760 | 0.892 | 0.896 | 0.182 |
| plasmode_unequal | e = 3.0 | 1.000 | 1.000 | 1.000 | 0.183 |

- **On the exports** minP rejects at most 2 more units of 500 than Bonferroni (446 → 448 on the unequal export at e = 2.0), and 2 more at each anchor.
- **In the simulation** minP equals Bonferroni in 40 of 42 cells and adds one unit of 100 in the other two (`seaad` orientation e = 0.25 at φ = 60° and 90°). The block maxima are nearly independent there: median Spearman −0.007 to 0.050.
- **Consequence.** At K = 2, with weakly correlated block maxima, minP would not change any reading. The design named a minP global test as a follow-up only if the cost were material, and it is not.

## 7. Simulation: localization and false hits

The two arms' global rates on the simulation units (`simulation/scan_table.csv`, `simulation/paired_table.csv`):

- **They agree closely.** No cell's rates differ by more than 0.03 (3 units of 100: `n200` shape e = 0.25, 0.85 joint against 0.82 by block; `seaad_weak` orientation e = 0.25 at φ = 30°, 0.08 against 0.05). No cell has more than 6 discordant units.
- **Signal cells keep their power.** Each cell where the joint scan had power reads within 0.03 by block: for example 0.76 at `n200` orientation e = 0.25, φ = 90° in both arms, and 0.96 → 0.98 at `seaad` orientation e = 0.5, φ = 90°.
- **At the weak design**, where nothing has power, the by-block arm reads 0.05–0.07 in the non-control cells, against the joint arm's 0.07–0.08.

The rate of units with a hit on a population-null feature, over the eight signal cells per design that have such features, at R = 100 (`simulation/false_hit_table.csv`):

| Design | Joint | By block |
|---|---|---|
| `n200` | 0.02–0.04 | 0.02–0.05 |
| `seaad` | 0.00–0.03 | 0.00–0.01 |
| `seaad_weak` | 0.00–0.08 | 0.00–0.07 |

At the anchor and translation every feature is population-null, so the rate is the global rate (S1). These rates are descriptive and stay within α + 2 SE at R = 100 (0.094) in every cell.

## 8. What follows

- **Reading a by-block hit.** It names a feature whose covariate-adjusted group difference changes across stages, with the familywise error over all scanned features controlled at 0.05, as for the joint scan. The difference is in who gets the error: each of the K blocks is tested at 0.05/K, whatever its size or noise.
- **When to use it.** When the blocks differ in size, or their noise is known to differ. On the unequal export the by-block scan named more planted features of the small block in 296 of 500 units at e = 2.0 and lost none, and it raised the global rate from 0.760 to 0.892. With equal blocks it mostly redistributes detection between them.
- **What it costs.** The larger or noisier block's features face a stricter bar: ATAC recall fell from 0.575 to 0.500 at e = 3.0 on the unequal export. With many small blocks each would get only 0.05/K. Only K = 2 was measured.
- **Default.** The joint scan stays the default; the by-block combination is opt-in (`interaction_scan_by_block`, `--interaction-scan-by-block`). A default flip would change `sead_sexdiff`'s pinned reading and needs its own change.
- **Case-study use.** `sead_sexdiff` may read the scan by block in its own change. It must pin MOTCO 0.23.0 and fix its multiplicity policy across the scan, the split test and the trajectory statistics before any real-label analysis. Its SEA-AD exports have equal blocks, where §4 shows the by-block scan trades ATAC detection for RNA detection at the same total.
- **Not needed:** a minP global test at K = 2 (§6).

## 9. Limitations

- **Two blocks only.** The Bonferroni split's cost grows with K, and no export with more blocks was run.
- **The unequal export is a construction:** a seeded random 500 of the RNA features, with the real noise. It tests a count imbalance with RNA's lighter null tail. A small ATAC block was not run.
- **The plant is stylized:** 20 features, 10 per block, one group, a stage ramp. The allocation gain is a statement about that plant.
- **Real-noise size and familywise control are conditional** on one data matrix and its subset, and on two effect sizes. The simulation gives unconditional size, at 498 features.
- **Within-block subset pivotality is approximate** under Freedman–Lane with covariates, as for the joint scan. S3 measures strong control on real noise in both exports.
- **Hits are detections, not drivers**, and the scan has no geometry, as in the parent report.
