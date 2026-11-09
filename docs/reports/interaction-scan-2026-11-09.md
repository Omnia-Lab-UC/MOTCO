# Per-feature interaction scan — findings report (2026-11-09)

**Change:** `interaction-feature-scan` (MOTCO 0.20.0). **Reads:** `results/interaction-scan-2026-11-09/` (`interaction_scan_decision.json`, `simulation/` and `plasmode/` scan, block-balance, comparator, false-hit and recovery tables, `accounting.json`, `reproduction.json`, `settings.json`, `plasmode-fidelity/fidelity.json`, `PROVENANCE.txt`, `sacct.txt`). **Concept page:** [measurement visibility](../measurement-visibility.md#testing-feature-by-feature-the-interaction-scan).

Every test MOTCO had pooled the group × stage difference over many dimensions:
- `delta`, `angle` and `shape` see only the part of a difference inside the fitted stage-PLS span;
- the split interaction test sums interaction energy over its parts.

Neither sees a difference confined to a few features. The SEA-AD plasmode planted one, a ramp to 3 within-stage SDs on 20 of 10,000 features in one group, and all six tests rejected it at the null rate ([plasmode](plasmode-2026-11-02.md), [split test](span-split-2026-11-05.md)). `motco.stats.interaction.interaction_feature_scan` tests each feature's interaction F on its own and combines the features by their maximum. Its global p-value is the permutation p-value of the largest F. Each feature gets a Westfall–Young step-down max-T adjusted p-value, under RRPP's Freedman–Lane draws and residual corrections. A hit is a feature with adjusted p ≤ 0.05.

**In one paragraph.** **The scan is validated: `scan_validated`.** It holds size and familywise control, and it sees the sparse difference every other test misses.

- **Size holds everywhere.** The six simulation control rates read 0.048–0.058, against a per-cell bound of 0.0695. On real noise at the case study's configuration (six covariates, 10,000 features) the scan reads 0.042.
- **Familywise control holds under the plasmode's partial null.** At e = 2.0 and 3.0, a feature outside the planted set was a hit in 0.042 of units at both strengths.
- **It sees the sparse plant.** At e = 3.0 the scan rejects in every unit (1.000), against 0.016–0.048 for the six existing tests on the same units. At e = 2.0 it rejects in 0.692.
- **It names the planted features.** At e = 3.0, 0.995 of the hits are planted features and half the plant is found (recall 0.525). The 20 largest F hold 0.887 of the plant, against attribution's committed 0.432.
- **The null maxima are not balanced across blocks.** On the SEA-AD export each block holds half the features, but ATAC gives 0.618 of the null maxima. In the simulation, methylation holds 0.737 of the features and gives 0.50. The verdict is unaffected, but the familywise threshold leans on one block.

## 1. Configuration and provenance

- **The parent.** Both profiles rerun the split-test validation's units (`results/span-split-2026-11-05/`). Each profile's settings are rebuilt through `scripts/span_split_check.py` at the parent's grade and must equal the committed `settings.json`. The sha256 of that file and of the parent's split and production (or trajectory) tables are bound in the settings and re-checked at every step.
- **Simulation profile.** The latent-space comparison's three designs: `n200`, `seaad`, and `seaad_weak` at the calibrated δ = 0.11. Each has methylation + expression (498 features), fourteen cells and base seed 2035.
  - **Replicates.** The anchor and translation ran 500 replicates, every other cell 100, all at 999 permutations: 6,600 units.
  - **Each unit** is the parent's unit unchanged: one `pls_cv` fit and evaluation, and the split test. The scan then runs on the split test's features and models at the unit's RRPP seed. Those features are the fit's standardized features in the PLS input coordinates, a per-feature affine map the scan ignores.
  - **Recorded per unit:** the scan's counts, its maximum and block, the global p-value, hits per block, the per-block share of null maxima, and hits on population-null features. A feature is population-null when its population B − A difference is constant across stages.
- **Plasmode profile.** The SEA-AD plasmode's units: the astrocyte `cohort_30` export under the plasmode's contract, within-stage label permutation, and the 20-feature plant (10 per block) with its seeds (base seed 2034, `paper` stream).
  - **Cells:** the anchor, e = 2.0 and e = 3.0 at 500 × 999, 1,500 units.
  - **Each unit** is one `analyze_trajectories` at the case study's configuration: six covariates, rank floor 2, net displacement, `df`, harness CV. The split test and the scan are on and attribution is off. The plasmode's consistency check runs on every unit, and its latent check on the anchor.
  - **Planted counts** are recorded as counts only, with no feature names.
- **Rejection.** The scan rejects at adjusted p ≤ 0.05 in both profiles. The six existing tests use their parent's rule (p < 0.05 in the simulation, p ≤ 0.05 in the plasmode), so every rate pairs with the parent's.
- **Code:** `359f57f` on `feat/interaction-feature-scan`, in a clean detached worktree on `ing`.
- **Cluster:** partition `512x1024`, arrays 936835 (simulation, 200 shards) and 936836 (plasmode, 150 shards), each throttled to 64 tasks. They ran on n1–n4 with n11 excluded and BLAS pinned to one thread. The fidelity gate passed on the cluster's export copy (rank 5, AUROC 0.635556).
- **Pre-flight:**
  - the fast suite (2,528 passed);
  - a local smoke of both profiles at 2 × 19;
  - a timing run at 999 permutations: the scan alone takes 0.72 s at 80 × 10,000, under 1 % of a plasmode unit.

Details are in `PROVENANCE.txt` and `openspec/changes/archive/2026-11-09-interaction-feature-scan/notes.md`.

## 2. Accounting and reproduction

- **Units:** 8,100 of 8,100 complete, 0 failed, 0 duplicates; all 350 array tasks completed. Every plasmode anchor unit passed the latent check.
- **Compute:** 69.6 recorded core-hours (simulation 22.6, plasmode 47.0), against the projection of 63 (budget 75), in 71 minutes of wall time. Unit medians were 11.5 s (simulation; the scan takes 0.067 s, 0.66 % of unit seconds) and 86 s (plasmode).
- **Reproduction.** Both profiles reproduced the split-test validation exactly, in count and rejections:
  - **Simulation:** 252 rows, the three split parts' pseudo-F and `delta`/`angle`/`shape` for every design and cell, under p < 0.05.
  - **Plasmode:** 18 rows, the same per cell, under p ≤ 0.05.

  The scan therefore pairs unit by unit with the split test and the trajectory statistics, and through them with the latent-space comparison and the plasmode.

## 3. The decision

`interaction_scan_decision.json`, rule version 1. S1–S3 use the stage-span paper grade's count-plus-ceiling rule. At R = 500 the per-cell bound is 0.0695, and the exceedance probability per rate is 0.030.

| Criterion | What | m | Rates | Exceedances (allowed) | Above ceiling | Holds |
|---|---|---|---|---|---|---|
| S1 | simulation size: global scan at the anchor and translation of every design | 6 | 0.048–0.058 | 0 (1) | 0 | yes |
| S2 | real-noise size: global scan at the plasmode anchor | 1 | 0.042 | 0 (0) | 0 | yes |
| S3 | strong familywise control: units with a hit outside the planted set at e = 2.0 and 3.0 | 2 | 0.042, 0.042 | 0 (1) | 0 | yes |
| G | sparse gain at e = 3.0: the scan against the best of six existing tests, by more than 2 pooled SE | — | 1.000 against 0.048 (`span`) | — | — | yes |

**S1.** The six control rates:

| Design | Anchor | Translation |
|---|---|---|
| `n200` | 0.056 | 0.052 |
| `seaad` | 0.054 | 0.048 |
| `seaad_weak` | 0.056 | 0.058 |

Translation is a constant group offset, which has no interaction, so it is a null for the scan.

**S3.** The 21 of 500 units with a false hit at each planted strength are the same replicates as the anchor's 21 global rejections (checked unit by unit on the per-unit records). Within a replicate the labels and the RRPP seed are shared across cells, and the plant changes only 20 features. So the null features' maxima are the anchor's, and step-down removes the planted features from the suffix only once they are rejected.

**G.** At e = 3.0 the existing tests read `delta` 0.026, `angle` 0.036, `shape` 0.016, and the split test's `total` 0.042, `span` 0.048 and `complement` 0.018. The comparator is the largest of the six, chosen after the fact, which is the conservative choice. The difference is 0.952 with a pooled SE of 0.0096.

**Verdict: `scan_validated`.** The consequence: a case study may report the per-feature interaction scan, read as the [measurement-visibility](../measurement-visibility.md#testing-feature-by-feature-the-interaction-scan) page states. No verdict reopens Phase 6 or changes its attribution leg. The plasmode's `attribution_leg_failed` stands, and the scan is a new test with its own claim.

## 4. Real noise: power and recovery on the plasmode

| Cell | Scan | `delta` | `angle` | `shape` | `total` | `span` | `complement` | Recall | Hit precision | Top-20 precision | Attribution top-20 (adjusted: observed / residual) |
|---|---|---|---|---|---|---|---|---|---|---|---|
| e = 2.0 | 0.692 | 0.030 | 0.034 | 0.020 | 0.042 | 0.046 | 0.022 | 0.059 | 0.966 (346 units) | 0.473 | 0.168 / 0.227 |
| e = 3.0 | 1.000 | 0.026 | 0.036 | 0.016 | 0.042 | 0.048 | 0.018 | 0.525 | 0.995 (500 units) | 0.887 | 0.344 / 0.432 |

- **Source and grade.** Everything except the last column comes from this run, 500 × 999 on the same units. The last column is the plasmode's committed attribution planted-set precision (2026-11-02, 500 × 999, B = 200), on the same metric.
- **Definitions.** Recall is the planted features found, divided by 20. Hit precision is the planted share of hits, over units with any hit. Top-20 precision is the planted share of the 20 largest F.
- **The proposal's back-of-envelope calculation held.** It assumed Gaussian noise and a Bonferroni bar, and was not a measurement. It gave each planted feature power 0.06 at e = 2.0 and 0.51 at e = 3.0, with at least one detected in 0.71 and 1.00 of units. The measured recall, 0.059 and 0.525, is the per-feature power. The global rates are 0.692 and 1.000.
- **e = 2.0 is detectable but thin.** The median unit has one hit. The scan says that the groups differ by stage on some feature, and names about one of the 20.
- **Ranking.** The scan's ranking recovers the plant better than attribution's at both strengths, without a latent space and without needing a prior global rejection.

## 5. Block balance

If every feature were equally likely to give a draw's maximum, a block's share of null maxima would equal its share of features. It does not:

| Profile | Cell | Block | Share of features | Share of null maxima |
|---|---|---|---|---|
| plasmode | anchor | ATAC | 0.500 | 0.618 |
| plasmode | anchor | RNA | 0.500 | 0.382 |
| simulation | anchor (`n200`, `seaad`, `seaad_weak`) | methylation | 0.737 | 0.501, 0.506, 0.506 |
| simulation | anchor (`n200`, `seaad`, `seaad_weak`) | expression | 0.263 | 0.499, 0.494, 0.494 |

- **In the simulation**, the methylation share stays at 0.497–0.523 in every cell. Methylation features are strongly correlated through the generator's coupling, so they act as fewer independent tests and win the maximum less often.
- **On the SEA-AD export**, ATAC takes more maxima than its share. Its per-feature null F has the heavier tail, or its features are less correlated. The shares are near-identical across the plasmode's cells (0.6178, 0.6178, 0.6177), because the plant changes 20 of 10,000 features.
- **Consequence.** Validity is unaffected: S2 and S3 hold on this export. But the common familywise threshold is set more by ATAC's null than by RNA's, so a difference of a given size is somewhat harder to detect on RNA features than per-block thresholds would make it. At e = 3.0 the hits are still balanced (5.27 ATAC, 5.30 RNA per unit, against 10 planted per block).
- **Follow-up.** The design named per-block maxima as the follow-up if the balance were materially off. On the export it is (0.62 / 0.38).

## 6. Simulation localization

The global scan beside the split test's parts and the production statistics, on the same units. The rates are this run's: the scan is judged at p ≤ 0.05 and the others at p < 0.05. Anchors and translations are at 500 replicates, other cells at 100.

| Design | Cell | Scan | `delta` | `angle` | `shape` | `total` | `span` | `complement` |
|---|---|---|---|---|---|---|---|---|
| `n200` | magnitude e = 0.05 | 0.09 | 0.74 | 0.02 | 0.00 | 0.02 | 0.56 | 0.01 |
| `n200` | orientation e = 0.25, φ = 0 | 1.00 | 0.00 | 0.87 | 0.05 | 0.80 | 0.93 | 0.11 |
| `n200` | orientation e = 0.25, φ = 90 | 0.76 | 0.05 | 0.05 | 0.01 | 0.19 | 0.02 | 0.29 |
| `n200` | orientation e = 0.5, φ = 90 | 1.00 | 0.20 | 0.08 | 0.51 | 1.00 | 0.52 | 1.00 |
| `n200` | orientation e = 0.5, φ = 90, 20 features | 1.00 | 0.68 | 0.07 | 0.96 | 0.90 | 0.79 | 0.78 |
| `n200` | shape e = 0.25 | 0.85 | 0.03 | 0.03 | 1.00 | 0.05 | 0.80 | 0.04 |
| `seaad` | magnitude e = 0.05 | 0.01 | 0.21 | 0.05 | 0.02 | 0.04 | 0.09 | 0.05 |
| `seaad` | orientation e = 0.25, φ = 0 | 1.00 | 0.00 | 0.47 | 0.01 | 0.12 | 0.45 | 0.07 |
| `seaad` | orientation e = 0.25, φ = 90 | 0.16 | 0.01 | 0.06 | 0.03 | 0.07 | 0.04 | 0.31 |
| `seaad` | orientation e = 0.5, φ = 90 | 0.96 | 0.03 | 0.10 | 0.20 | 0.36 | 0.19 | 0.68 |
| `seaad` | orientation e = 0.5, φ = 90, 20 features | 1.00 | 0.13 | 0.09 | 0.56 | 0.16 | 0.36 | 0.58 |
| `seaad` | shape e = 0.25 | 0.20 | 0.02 | 0.06 | 0.70 | 0.04 | 0.39 | 0.07 |
| `seaad` | shape e = 0.5 | 0.83 | 0.00 | 0.08 | 0.94 | 0.06 | 0.43 | 0.07 |

The full table, every design and cell, is `simulation/comparator_table.csv`.

- **The design's prediction was half wrong.** It said the scan would trail the complement on the dense φ = 90 rotation, which spreads over hundreds of features. That holds at e = 0.25 on the SEA-AD table (0.16 against 0.31). At e = 0.5 the scan ties the complement at n = 200 (1.00) and beats it on the SEA-AD table (0.96 against 0.68). At n = 200 it also beats the complement at e = 0.25 (0.76 against 0.29). With 498 features the search costs little.
- **The scan is weak where a difference spreads thinly over many features.** Magnitude reads 0.09 / 0.01 against `delta`'s 0.74 / 0.21, and the weak SEA-AD bend 0.20 against `shape`'s 0.70.
- **It is strong on in-span rotation:** 1.00 against `angle`'s 0.47 on the SEA-AD table at e = 0.25.
- **The tests are complementary:** each targets a different spread of the difference across features.
- **At the weak design nothing has power, the scan included.** The non-control cells read 0.07–0.08 against controls of 0.056–0.058: the stage-span plants scale with a stage signal that barely exists there.

## 7. Simulation false hits

The rate of units with a hit on a population-null feature (a feature whose population B − A difference is constant across stages), among the units at cells that have such features (R = 100):

- **`n200`:** 0.02–0.04 over the eight signal cells with population-null features.
- **`seaad`:** 0.00–0.03 over the same eight cells.
- **`seaad_weak`:** 0.00–0.08, where the plant is negligible and the global rate itself is 0.07–0.08.

At the anchor and translation every feature is population-null, so the rate is the global rate (§3, S1). At the φ = 30° and 60° rotations no feature is population-null. These are the generator-noise analogue of S3. They are descriptive and stay within α + 2 SE at R = 100 (0.094) in every cell.

## 8. What follows

- **Case-study use.** Under `scan_validated`, `sead_sexdiff` may apply the scan to both SEA-AD cell types in its own change. It must pin MOTCO 0.20.0 and fix its multiplicity policy across the scan, the split test and the trajectory statistics before any real-label analysis.
- **Reading.** A hit names a feature whose covariate-adjusted group difference changes across stages, with the familywise error rate over the scanned features controlled at 0.05. A scan without hits is uninformative unless power is shown. At the plasmode's design that power is 0.69 for a 20-feature ramp to 2 SDs and 1.00 at 3 SDs.
- **Possible follow-ups**, each its own change if proposed:
  - per-block maxima, since the export's null maxima split 0.62 / 0.38;
  - a moderated statistic, if a case study needs more than the measured per-feature power: 0.06 at 2 SDs, 0.53 at 3 SDs.
- **Phase 6.** Its attribution leg stays failed (item 12). The scan answers a different question, which features differ and not why the latent trajectory differs, and it is not a rescue of that rule.

## 9. Limitations

- **The plant is stylized:** 20 features, one group, a stage ramp. G is a statement about the documented gap, not about realistic power.
- **Real-noise size and familywise control are conditional** on one data matrix (the `cohort_30` export) and two effect sizes. The simulation gives unconditional size, but at 498 features, not 10,000.
- **The F is the plain per-feature F**, not a variance-moderated one. At n = 80 it is noisy (§4).
- **The null maxima are imbalanced across blocks** (§5). The threshold leans on the block with the heavier null tail.
- **Hits are detections, not drivers.** Correlated neighbours of a differing feature can be hits too, and the scan says nothing about the latent trajectory.
- **The scan has no geometry.** It says which features differ by stage, never how the trajectories differ.
