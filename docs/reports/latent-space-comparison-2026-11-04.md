# Latent-space comparison — findings report (2026-11-04)

**Change:** `latent-space-comparison` (MOTCO 0.16.0). **Reads:** `results/latent-space-comparison-2026-11-04/` (`rate_table.csv`, `capture_table.csv`, `accounting.json`, `latent_space_comparison_reading.json`, `settings.json`, `calibration/`, `calibration-original-ladder/`, `PROVENANCE.txt`, `sacct.txt`). **Concept page:** [measurement visibility](../measurement-visibility.md).

The stage-span pilot ([report](stage-span-pilot-2026-11-04.md)) established that the production measurement sees a group difference only through its component inside the stage-supervised PLS span. Comparing spaces on any one construction is circular: a plant inside the span favours stage-PLS, and one outside it favours a label-free space. How real group differences sit relative to the stage program is the unknown. This pilot parametrizes it as the **alignment** φ of a stage-span orientation with the stage span, at a fixed intended geometry, and reads five measurement spaces across it on identical data:

- `pls_cv`, the production measurement (stage PLS-DA, double-CV rank floored at 2);
- `pca_matched`, label-free PCA at that unit's `pls_cv` rank;
- `pca_fixed`, PCA at rank 2(K − 1) = 4, the orientation probe's arm;
- `pls_residual_pca`, the `pls_cv` scores followed by K − 1 = 2 principal components of the stage-residualized features, taken orthogonal to the PLS span;
- `concat`, the 498 standardized features, the reference ceiling.

**In one paragraph.** **No space buys out-of-span visibility without losing what stage-PLS sees.** The advisory reading returns **`document_tradeoff`**, with every eligible arm `no_gain`.

- **Size holds for every eligible arm at every design** (largest control rate 0.08 against the bound 0.0936).
- **The low-dimensional spaces are blind to the out-of-span rotation.** At φ = 90 every PLS- and PCA-based arm keeps 0.006–0.075 of the difference's energy and rejects `angle` at 0.05–0.18. `pls_residual_pca` comes closest to a gain: 0.17 and 0.18 against `pls_cv`'s 0.08 and 0.11, short of 2 pooled SE at both δ = 2.0 designs. Its residual components do not find the difference (capture 0.024 and 0.075).
- **Only `concat` sees it**, rejecting the 45° out-of-span rotation in 1.00 (n = 200) and 0.90 (SEA-AD table) of replicates.
- **Added dimensions cost in-span power.** `pls_residual_pca` shares `pls_cv`'s coordinates and adds two; magnitude `delta` falls from 0.74 to 0.37 at n = 200 and from 0.23 to 0.07 on the SEA-AD table. `concat` and `pca_fixed` lose as much, and `pca_matched` at rank 2 loses everything, keeping only 0.53–0.54 of the pooled stage means.
- **At SEA-AD's stage signal no space has power** for these plants. The fitted stage-PLS span keeps only 0.10 of the population stage program there. `concat`'s `angle` null is liberal at that design (anchor 0.21).
- **Support.** At a fixed norm, a 20-feature out-of-span direction is more visible than a dense one: `pls_cv` keeps 0.108–0.116 of it against 0.008–0.026 (differences +0.09 to +0.10, SE ≤ 0.003). It is no easier to detect: `pls_cv`'s `angle` difference is −0.02, within 1 SE.

## 1. Configuration and provenance

- **Designs** (three stages, methylation + expression, `p_dmp` 0.1, ρ = 0, joint magnitude, `construction = stage_span` with span layers = the measured blocks):
  - `n200`: n = 200, group ratio 0.5, every per-omic δ 2.0 (cells A 34/34/33, B 33/33/33);
  - `seaad`: the SEA-AD table `[[11, 10, 28], [9, 10, 12]]`, δ = 2.0;
  - `seaad_weak`: the same table at δ = 0.11, calibrated to SEA-AD's held-out stage AUROC (§2).
- **Cells** (14 per design): the anchor; translation e = 1.0; magnitude e = 0.05; orientation e ∈ {0.25, 0.5} × φ ∈ {0°, 30°, 60°, 90°} with a dense out-of-span direction; shape e ∈ {0.25, 0.5}; and orientation e = 0.5 at φ = 90° on a 20-feature support. Orientation e = 0.25 / 0.5 is a 22.5° / 45° rotation at every φ, with intended `delta` = `shape` = 0.
- **Arms:** five per unit on one dataset and one RRPP seed, all with net displacement, rank floor 2 and span capture. One stage-PLS fit per unit is shared by `pls_cv`, `pca_matched`'s rank and `pls_residual_pca`.
- **Grade:** 100 replicates × 199 permutations, α = 0.05, base seed 2035 (generator seed per design and replicate, so the 14 cells of a replicate share their null dataset; RRPP seed per unit; calibration on its own stream). 4,200 units, 21,000 arm records.
- **Code:** `d2c167b` on `feat/latent-space-comparison` (the run and the amended calibration), in a clean detached worktree on `ing`; merge and reading at `36f8fdc`, which only adds the population `delta`/`angle`/`shape` medians to the capture table.
- **Cluster:** partition `512x1024`, array 933442 (200 shards, throttle 80) on n1, n3 and n4; n2 and n11 excluded. BLAS pinned to one thread.
- **Pre-flight:** the fast suite (2,212 passed), a local smoke of every design at 2 × 19 (84/84 units), and 28 units timed at 199 permutations.

Details are in `PROVENANCE.txt` and `openspec/changes/latent-space-comparison/notes.md`.

## 2. Calibration of the weak design

The weak design's δ was chosen before the run from an ascending ladder (all three per-omic δ together, 40 anchor replicates of the SEA-AD table per step, `pls_cv` fits only), adopting the first step whose median held-out stage AUROC lies in [0.60, 0.67].

- **The predeclared ladder missed the band.** Its first two steps read 0.595 (δ = 0.10) and 0.710 (δ = 0.15), so the rule would have adopted δ = 0.10 as `band_missed` (job 933435, `calibration-original-ladder/`).
- **The ladder was amended before any run unit** (design D9) with steps 0.11–0.14, under the same rule (job 933439, `calibration/`). The calibration seeds do not depend on δ, and every original step reproduced exactly.

| δ | 0.10 | 0.11 | 0.12 | 0.13 | 0.14 | 0.15 | 0.20 | 0.25 | 0.30 | 0.40 | 0.50 |
|---|---|---|---|---|---|---|---|---|---|---|---|
| median AUROC | 0.595 | **0.617** | 0.633 | 0.660 | 0.687 | 0.710 | 0.838 | 0.922 | 0.973 | 0.999 | 1.000 |
| median rank | 2.5 | 3 | 3 | 4 | 4 | 4 | 6 | 7 | 7 | 7 | 4.5 |

**Adopted: δ = 0.11, `first_in_band`** (median 0.617, quartiles 0.549–0.645). The rule takes the first in-band step: δ = 0.12 (0.633) is nearer SEA-AD's 0.636. SEA-AD's own CV rank is 5; the adopted step's median rank is 3.

## 3. Accounting and identity checks

- **Units.** 4,200/4,200 complete, 0 failed records, 0 duplicates. All 200 array tasks `COMPLETED 0:0` (`sacct.txt`).
- **Cost.** 37.1 recorded core-hours (median 38.2 s per unit; `concat`'s RRPP in 498 dimensions is 16.0 s of it, every other arm 0.43–0.45 s). The projection was 35 core-hours.
- **Identity checks** (tests and the re-merge):
  - A `pls_cv` record equals a production evaluation of the same dataset at the same RRPP seed (rank, statistics, p-values and span capture; `tests/test_latent_space_comparison.py`).
  - Without the new settings and method, datasets, cell ids, signatures and PLS span-capture records are byte-identical to 0.15.0 (`tests/test_span_alignment.py`, `tests/test_span_capture.py`).
  - Re-merging at `36f8fdc` left `rate_table.csv`, `accounting.json` and the merged records byte-identical.

## 4. Visibility against alignment

Median `group_difference` energy share inside each arm's fitted span (`capture_table.csv`), orientation e = 0.5. Random-direction shares are 0.004 (rank 2), 0.008 (rank 4) and 1 (`concat`).

| arm | n200 φ = 0 | 30 | 60 | 90 | 90, k = 20 | SEA-AD φ = 0 | 30 | 60 | 90 | 90, k = 20 |
|---|---|---|---|---|---|---|---|---|---|---|
| `pls_cv` | 0.971 | 0.818 | 0.395 | 0.008 | 0.108 | 0.927 | 0.794 | 0.412 | 0.026 | 0.116 |
| `pca_matched` | 0.971 | 0.851 | 0.217 | 0.006 | 0.031 | 0.928 | 0.827 | 0.452 | 0.032 | 0.122 |
| `pca_fixed` | 0.974 | 0.857 | 0.450 | 0.024 | 0.131 | 0.936 | 0.833 | 0.470 | 0.051 | 0.143 |
| `pls_residual_pca` | 0.974 | 0.829 | 0.412 | 0.024 | 0.116 | 0.932 | 0.813 | 0.450 | 0.075 | 0.131 |
| `concat` | 1.000 | 1.000 | 1.000 | 1.000 | 1.000 | 1.000 | 1.000 | 1.000 | 1.000 | 1.000 |

- **Every low-dimensional space has the same visibility profile.** Capture falls from 0.93–0.97 at φ = 0 to the random-direction share at φ = 90 (0.006–0.075). The axis behaves as constructed: at φ = 90 the difference is orthogonal to the pooled stage span, and no basis built from the stage structure or the leading variance holds it.
- **`pls_residual_pca` adds least where it was meant to add most.** Its two residual components raise capture at φ = 90 from 0.008 to 0.024 (n = 200) and from 0.026 to 0.075 (SEA-AD).
- **Seen through each basis** (`population_latent`, median `angle`, intended 45°): `pls_cv` reads 44.6° / 38.1° / 21.6° / 0.3° at φ = 0 / 30 / 60 / 90 at n = 200 and 43.5° / 37.4° / 22.1° / 1.8° on the SEA-AD table. `pls_residual_pca` reads 5.5° and 9.8° at φ = 90. `concat` reads 45.1° / 44.5° / 42.7° / 41.7° and 45.0° / 44.5° / 43.1° / 42.4°.

## 5. Out-of-span power

`angle` rejection rate at α = 0.05 (`rate_table.csv`; SE ≤ 0.05 at R = 100):

| arm | n200 e = 0.25: φ = 0 / 30 / 60 / 90 | n200 e = 0.5: φ = 0 / 30 / 60 / 90 / k = 20 | SEA-AD e = 0.5: φ = 0 / 30 / 60 / 90 / k = 20 |
|---|---|---|---|
| `pls_cv` | 0.84 / 0.78 / 0.82 / 0.05 | 1.00 / 1.00 / 0.79 / 0.08 / 0.06 | 0.96 / 0.73 / 0.43 / 0.11 / 0.09 |
| `pca_matched` | 0.18 / 0.12 / 0.03 / 0.04 | 0.89 / 0.79 / 0.26 / 0.05 / 0.05 | 0.72 / 0.54 / 0.10 / 0.06 / 0.07 |
| `pca_fixed` | 0.36 / 0.24 / 0.07 / 0.03 | 1.00 / 1.00 / 0.67 / 0.09 / 0.07 | 0.75 / 0.53 / 0.15 / 0.09 / 0.08 |
| `pls_residual_pca` | 0.41 / 0.28 / 0.05 / 0.03 | 1.00 / 1.00 / 0.54 / 0.17 / 0.07 | 0.93 / 0.79 / 0.35 / 0.18 / 0.09 |
| `concat` | 0.44 / 0.41 / 0.36 / 0.36 | 1.00 / 1.00 / 1.00 / 1.00 / 0.98 | 0.81 / 0.85 / 0.89 / 0.90 / 0.49 |

- **The production space's `angle` power follows the visible part of the rotation.** It holds while any of it remains in the span (0.79 at φ = 60 at n = 200, with capture 0.40) and falls to the null at φ = 90.
- **Only `concat` keeps its power across the axis.** At e = 0.25 it is flat at 0.36–0.44 (n = 200), below `pls_cv`'s in-span 0.84. At e = 0.5 on the SEA-AD table it is 0.81–0.90. Its sparse-support cell drops to 0.49 there (§8).
- **Off-target registration.** At φ = 90, e = 0.5, `pls_cv` rejects `shape` in 0.47 and `delta` in 0.21 of replicates at n = 200, and `shape` in 0.18 on the SEA-AD table; no other arm exceeds 0.10 there. The stage-span rotation is rigid in null-standardized coordinates, but every space measures it through the replicate's sample scaler, which absorbs the plant's between-cell spread. Through that scaler the population means carry a size and shape difference in every space: at φ = 90, n = 200, population `delta` / `shape` read 0.559 / 0.034 through `pls_cv` and 0.533 / 0.027 through `concat`. This is the "unequal feature scales" route on the [visibility page](../measurement-visibility.md). The pilot does not explain why only `pls_cv` rejects on it, and only at φ = 90: its population `delta` at φ = 0 is similar (0.464) and is rejected in 0.00.

## 6. Size

Anchor and translation rejection rates (`delta` / `angle` / `shape`), against the bound α + 2·√(α(1 − α)/R) = 0.0936:

| arm | n200 anchor | n200 translation | SEA-AD anchor | SEA-AD translation | weak anchor | weak translation |
|---|---|---|---|---|---|---|
| `pls_cv` | 0.05 / 0.04 / 0.00 | 0.03 / 0.04 / 0.00 | 0.02 / 0.04 / 0.03 | 0.00 / 0.06 / 0.01 | 0.06 / 0.02 / 0.05 | 0.05 / 0.02 / 0.04 |
| `pca_matched` | 0.00 / 0.05 / 0.01 | 0.00 / 0.05 / 0.00 | 0.00 / 0.06 / 0.01 | 0.00 / 0.06 / 0.00 | 0.01 / 0.07 / 0.06 | 0.02 / 0.07 / 0.05 |
| `pca_fixed` | 0.00 / 0.05 / 0.00 | 0.00 / 0.05 / 0.00 | 0.00 / 0.06 / 0.01 | 0.00 / 0.06 / 0.00 | 0.01 / 0.08 / 0.04 | 0.02 / 0.07 / 0.05 |
| `pls_residual_pca` | 0.00 / 0.05 / 0.00 | 0.00 / 0.04 / 0.00 | 0.00 / 0.06 / 0.01 | 0.00 / 0.06 / 0.00 | 0.05 / 0.03 / 0.06 | 0.04 / 0.02 / 0.05 |
| `concat` | 0.00 / 0.05 / 0.00 | 0.00 / 0.05 / 0.00 | 0.00 / 0.06 / 0.01 | 0.00 / 0.06 / 0.00 | 0.01 / **0.21** / 0.04 | 0.02 / **0.22** / 0.05 |

- **Every eligible arm holds size at every design.** The fixed group-blind bases (design D5) do not repeat `pls_within_group`'s failure.
- **`concat`'s `angle` is liberal at the weak design** (0.21 and 0.22). Its `angle` reads 0.19–0.26 in every cell there, so its weak-design `angle` rates measure that excess, not power. It is reference-only and never eligible.
- **`pls_cv` holds every control** (`pls_cv_controls` in the reading). The n = 200 anchor `angle` that read 0.10 in the stage-span pilot reads 0.04 here, on an independent seed family.

## 7. In-span cost: information against power

Target rates on the cells `pls_cv` already sees (the reading's protected criterion is orientation e = 0.5 φ = 0, shape e = 0.5 and magnitude):

| arm | n200 orient. 0.25 / 0.5 `angle` | shape 0.25 / 0.5 `shape` | magnitude `delta` | SEA-AD orient. 0.25 / 0.5 `angle` | shape 0.25 / 0.5 `shape` | magnitude `delta` |
|---|---|---|---|---|---|---|
| `pls_cv` | 0.84 / 1.00 | 1.00 / 1.00 | 0.74 | 0.46 / 0.96 | 0.70 / 0.95 | 0.23 |
| `pca_matched` | 0.18 / 0.89 | 0.20 / 0.31 | 0.08 | 0.07 / 0.72 | 0.10 / 0.47 | 0.00 |
| `pca_fixed` | 0.36 / 1.00 | 0.98 / 1.00 | 0.35 | 0.09 / 0.75 | 0.22 / 0.85 | 0.01 |
| `pls_residual_pca` | 0.41 / 1.00 | 1.00 / 1.00 | 0.37 | 0.14 / 0.93 | 0.33 / 0.93 | 0.07 |
| `concat` | 0.44 / 1.00 | 0.98 / 1.00 | 0.37 | 0.09 / 0.81 | 0.12 / 0.73 | 0.01 |

- **`pls_residual_pca` is a clean test of the power cost of extra dimensions.** Its first coordinates are exactly `pls_cv`'s, and it adds two residual components that carry almost none of the in-span plant. At n = 200 the 22.5° rotation's `angle` falls from 0.84 to 0.41, and magnitude `delta` from 0.74 to 0.37. On the SEA-AD table they fall from 0.46 to 0.14 and from 0.23 to 0.07. The saturated cells (e = 0.5) are unaffected.
- **`concat` contains everything `pls_cv` sees and has less power for it**, at both designs. This is the [visibility page](../measurement-visibility.md)'s information-against-power distinction: a projection cannot add information, but by discarding noise directions it can gain power.
- **`pca_matched` fails for a different reason.** At rank 2 its leading principal axes miss the stage program: they keep 0.538 (n = 200) and 0.527 (SEA-AD) of the pooled stage means' energy, against `pls_cv`'s 0.970 and 0.932 at the same rank. Supervision, not dimension, is what puts the stage program in a two-dimensional space.

## 8. Sparse against dense support

The reading's `sparse_minus_dense` table: orientation e = 0.5 at φ = 90, a 20-feature out-of-span direction minus the dense one (every span-layer feature the null program does not use), at the same plant norm. Differences with pooled SE:

| design | arm | capture | `angle` |
|---|---|---|---|
| n200 | `pls_cv` | +0.102 (0.002) | −0.02 (0.036) |
| n200 | `pca_fixed` | +0.101 (0.004) | −0.02 (0.038) |
| n200 | `pls_residual_pca` | +0.065 (0.010) | −0.10 (0.045) |
| n200 | `concat` | 0 | −0.02 (0.014) |
| SEA-AD | `pls_cv` | +0.091 (0.003) | −0.02 (0.042) |
| SEA-AD | `pca_fixed` | +0.080 (0.004) | −0.01 (0.039) |
| SEA-AD | `pls_residual_pca` | +0.033 (0.010) | −0.09 (0.048) |
| SEA-AD | `concat` | 0 | −0.41 (0.058) |

- **Capture depends on the support; power barely does.** The sparse direction keeps 0.108–0.116 of its energy in the stage-PLS span, against 0.008–0.026 for the dense one (4–13 times more). That is still far from visible: `pls_cv`'s and `pca_fixed`'s `angle` power is unchanged (−0.01 to −0.02, within 1 SE), and `pls_residual_pca`'s is 0.09–0.10 lower (about 2 SE). The visibility page's claim that blindness is a direction, not a set of features, holds for power at these designs. Its wording for capture needs qualifying: a direction on few unused features keeps several times the random share.
- **The sample scaler bends sparse plants.** Through `concat`, whose basis keeps everything, the sparse rotation reads 31.8° (n = 200) and 33.8° (SEA-AD) against an intended 45° (dense: 41.7° and 42.4°). Large shifts on 20 features inflate those features' pooled SDs, so standardization shrinks the plant (design risk R6). This plausibly also explains the sparse cell's higher capture in the low-dimensional arms, through the same unequal-scales route; the pilot does not test that link.

## 9. The weak design

At δ = 0.11 (median held-out stage AUROC 0.617) no arm has power for any plant: every target rate is 0.01–0.10, except `concat`'s liberal `angle` (§6).

- **The fitted stage-PLS span keeps only 0.100 of the population stage program's energy** (pooled-stage capture at the anchor, median rank 3, random share 0.006); `pls_residual_pca` keeps 0.114. The PCA arms sit at the random share (0.006 and 0.009).
- **The in-span plants are therefore mostly invisible too.** `pls_cv` keeps 0.082–0.086 of the φ = 0 rotation. Through its basis the 45° rotation reads 25.3°, and through the PCA arms 33–40° at every φ: a noise basis gives a noise angle.
- **What this says about the case study.** Stage-span plants are relative to the configuration (design non-goal), so they are as weak as the program they reshape, and this design is descriptive. It adds one number to the [measurement visibility](../measurement-visibility.md) argument for a SEA-AD-like stage signal: the basis fitted there is about 90 % noise directions, measured against the population stage means.

## 10. The reading

`latent_space_comparison_reading.json` (rule version 1, advisory, bound 0.0936, comparisons at 2 pooled SE) returns **`document_tradeoff`**:

| arm | C1 size (every design) | C2 out-of-span gain (n200, SEA-AD) | C3 protected (n200, SEA-AD) | verdict |
|---|---|---|---|---|
| `pca_matched` | holds | fails: 0.05 vs 0.08; 0.06 vs 0.11 | fails at all six checks | `no_gain` |
| `pca_fixed` | holds | fails: 0.09 vs 0.08; 0.09 vs 0.11 | fails: magnitude at both, orientation and shape at SEA-AD | `no_gain` |
| `pls_residual_pca` | holds | fails: +0.09 (threshold 0.093); +0.07 (0.099) | fails: magnitude at both (0.37 vs 0.74; 0.07 vs 0.23) | `no_gain` |

`concat` is never eligible. On the reading's own terms it would pass C2 (1.00 and 0.90), but fail C1 at the weak design (`angle` 0.21 / 0.22). It would also fail C3 on magnitude at both designs (0.37 vs 0.74; 0.01 vs 0.23), and on orientation (0.81 vs 0.96) and shape (0.73 vs 0.95) on the SEA-AD table.

## 11. What follows

The reading is advisory and proposes no paper grade. The evidence supports three statements and one follow-up question:

- **The production space stays.** No candidate sees an out-of-span orientation without losing in-span power, and the one that sees it (`concat`) is liberal at a SEA-AD-like stage signal.
- **The trade-off is between dimension and noise, not supervision and information.** A space that adds even two near-empty directions halves `delta` power at n = 200 (§7). A label-free space at the production rank misses the stage program (§7). A within-stage residual space does not find a dense out-of-span difference at these designs (§4).
- **At SEA-AD's stage signal the measurement span is mostly noise** (§9), so the case study's non-rejections carry less information than the simulator's δ = 2.0 results imply. This is already the case study's binding reading; the pilot measures it in the space itself.
- **Open:** an out-of-span question would need a space whose added directions are chosen for the group difference. A group-aware basis must be refitted inside every permutation to hold size (the `pls_within_group` lesson), which multiplies the cost by the permutation count. Whether any such space is worth that cost is a separate change, as is a four-stage design where an out-of-span bend is defined (design open question).

## 12. Limitations

- **Pilot grade.** At R = 100, rates near the bound are not resolved (SE 0.022 at α), and a 2-pooled-SE gain of about 0.09 is the smallest the reading can credit. `pls_residual_pca`'s +0.09 at n = 200 sits at that threshold.
- **The calibration ladder was amended after one look at its result**, before any run unit, under the predeclared rule (§2). The adopted step is below SEA-AD's AUROC (0.617 against 0.636) and its rank below SEA-AD's (3 against 5).
- **The alignment axis is defined for orientation only.** At three stages an out-of-span bend changes the intended `shape` (design D2), so `shape` is measured at φ = 0 only.
- **Naive spaces.** PCA weights methylation (367 features) over expression (131); block weighting is roadmap item 5. There is no data-driven PCA rank rule.
- **The construction is relative.** Stage-span plants scale with the configuration, so the weak design is descriptive; absolute-scale plants on real noise are the plasmode's question.
- **The sample scaler distorts the planted geometry** (§5, §8). Population geometry through each basis is reported beside the intended values, but how much of the off-target registration it explains is not resolved.
