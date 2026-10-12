# Orientation measurement probe — 2026-10-12

> Exploratory, matched-seed probe (`openspec/changes/archive/2026-10-12-orientation-measurement-probe/`). It ranks measurement
> arms for a paper-grade follow-up. It certifies none of them and changes no production default. Every rate
> below is read from `results/orientation-measurement-probe-2026-10-12/` (`rate_table.csv`,
> `checkpoint_table.csv`, `advancement_decision.json`, `summary.md`). The motivating evidence in §1 is read from
> `evidence/phase6_checkpoint_evidence.md`, produced by `scripts/orientation_checkpoint_evidence.py`. Rates are
> rejection rates at p < 0.05 over 100 replicates (Monte Carlo SE about 0.022 at 0.05 and 0.05 at 0.5).

**Run:** `results/orientation-measurement-probe-2026-10-12/` · **Driver:** `scripts/orientation_measurement_probe.py`
at `49134c3` · **Decision:** no arm advances (`advancement_decision.json`, advisory).

**In one paragraph.** Two levers were crossed on the same datasets: the measurement space (five spaces) and the
orientation estimator (PC1 or net displacement). The estimator turned out to be the larger lever. In the
production stage-PLS space, net displacement more than doubles orientation/`angle` power at the Phase 6 SEA-AD
design (0.29 → 0.67). At Phase 5 it removes the shape construction's leak into `angle` (0.89 → 0.02), and it costs
no other power or specificity. That arm fails the predeclared rule on one criterion only: its Phase 6 anchor
`angle` rate is 0.100, against a bound of 0.0936. Every net-displacement arm shows a similar anchor excess at both
designs (0.08–0.10, against 0.01–0.05 under PC1). At 100 replicates this is suggestive, not established, and it
is the question the follow-up must settle first. The alternative spaces do not help on their own terms:
- `pls_within_group` without a per-permutation refit is invalid (anchor `angle` up to 0.88).
- `pca` and `pls_forced` at rank 2(K − 1) recover the rotation, but they lose Phase 6 shape power (0.19 → 0.00–0.01) and leak `angle` into the shape or magnitude cells.

## 1. Motivating evidence

From `evidence/phase6_checkpoint_evidence.md` (Phase 6 small-n orientation profile, baseline column,
`results/phase6-small-n-2026-10-09/merged.jsonl`, sha256 `6a91a47e…`):

| e | population (native methylation) | observed std. | PLS latent, rank ≤ 3 | PLS latent, rank > 3 | `angle` rate (≤ 3 / > 3) |
|---|---|---|---|---|---|
| 0.25 | 41.4° | 58.3° | 21.6° | — | 0.034 / — |
| 0.50 | 60.0° | 71.5° | 24.8° | — | 0.058 / — |
| 0.75 | 75.7° | 81.8° | 28.5° | — | 0.090 / — |
| 1.00 | 90.0° | 89.8° | 32.1° (n = 243) | 76.0° (n = 257) | 0.082 / 0.506 |

The surgery realizes arccos(1 − e) exactly in the population, and the observed standardized space keeps it. The
stage-supervised PLS space collapses it unless CV keeps more than three components. The signed anchor `angle` has
a heavy tail: q95 = 137.1° in the observed space and 135.6° in the latent space. On that oracle threshold the
signed observed angle separates e = 1 from the anchor in 0.000 of replicates, and its folded counterpart in 0.996
(0.496 in the latent space).

## 2. Design and provenance

- **Arms:** five spaces × two estimators, giving ten arms:
  - Spaces: `pls_cv` (production stage-PLS, CV rank; with PC1 this is the **reference arm**); `pls_forced` (stage-PLS at rank 2(K − 1)); `pls_within_group` (K − 1 components per group, fixed basis, no refit inside RRPP); `pca` (rank 2(K − 1)); and `concat` (ceiling reference, never eligible).
  - Estimators: `pc1` and `net_displacement`.
- **Designs:**
  - Phase 6 SEA-AD: three stages, `group_stage_sizes = [[11, 10, 28], [9, 10, 12]]`, methylation + expression measured.
  - Phase 5: n = 1200, four stages, all three blocks measured.
  - Both use `p_dmp = 0.1`, ρ = 0, joint magnitude and `surgery_censoring = error`.
  - **Deviation:** the Phase 5 design uses the joint magnitude construction, not the methylation-only `all` pinned in `phase5_power_study.json`. The reference arm's rates are recomputed here, and every comparison is against them.
- **Cells:** the anchor; orientation at e = 0.5 and 1.0; magnitude and shape at e = 1.0. Each cell has 100 replicates × 199 permutations.
- **Seeds:** the generator seed is matched across cells per (design, replicate). One RRPP seed per unit is shared by all ten arms, so arms differ only in measurement.
- **Run:** SLURM array 912606 on `ing`, 50 shards. 1,000/1,000 units and 10,000 arm records, with 0 failures. It took 28.2 recorded core-hours, 28 % over the 22 core-hour budget written before launch (`PROVENANCE.txt`).

## 3. Rejection rates

From `rate_table.csv` / `summary.md`. "orient." is orientation; "mag." is magnitude.

**Phase 6 (n = 80, K = 3).**

| arm | anchor `angle` | orient. 0.5 `angle` | orient. 1.0 `angle` | mag. `delta` | mag. `angle` | shape `shape` | shape `angle` |
|---|---|---|---|---|---|---|---|
| `pls_cv__pc1` (ref.) | 0.03 | 0.07 | 0.29 | 1.00 | 0.01 | 0.19 | 0.10 |
| `pls_cv__net_displacement` | **0.10** | 0.08 | **0.67** | 1.00 | 0.05 | 0.19 | 0.06 |
| `pls_forced__pc1` | 0.03 | 0.14 | 0.40 | 1.00 | 0.01 | 0.01 | 0.34 |
| `pls_forced__net_displacement` | 0.09 | 0.74 | 0.99 | 1.00 | 0.15 | 0.01 | 0.06 |
| `pls_within_group__pc1` | 0.06 | 0.23 | 0.14 | 1.00 | 0.01 | 0.03 | 0.64 |
| `pls_within_group__net_displacement` | 0.75 | 1.00 | 1.00 | 1.00 | 0.42 | 0.03 | 0.90 |
| `pca__pc1` | 0.01 | 0.03 | 0.69 | 1.00 | 0.01 | 0.00 | 0.39 |
| `pca__net_displacement` | 0.10 | 0.45 | 0.99 | 1.00 | 0.08 | 0.00 | 0.06 |
| `concat__pc1` | 0.01 | 0.06 | 0.14 | 1.00 | 0.01 | 0.00 | 0.46 |
| `concat__net_displacement` | 0.10 | 1.00 | 1.00 | 1.00 | 0.09 | 0.00 | 0.06 |

**Phase 5 (n = 1200, K = 4).**

| arm | anchor `angle` | orient. 0.5 `angle` | orient. 1.0 `angle` | mag. `delta` | mag. `angle` | shape `shape` | shape `angle` |
|---|---|---|---|---|---|---|---|
| `pls_cv__pc1` (ref.) | 0.04 | 0.86 | 0.89 | 1.00 | 0.01 | 1.00 | 0.89 |
| `pls_cv__net_displacement` | 0.08 | 0.93 | 0.89 | 1.00 | 0.00 | 1.00 | **0.02** |
| `pls_forced__pc1` | 0.05 | 0.85 | 0.84 | 1.00 | 0.03 | 0.65 | 0.95 |
| `pls_forced__net_displacement` | 0.09 | 0.99 | 1.00 | 1.00 | 0.32 | 0.65 | 0.03 |
| `pls_within_group__pc1` | 0.48 | 0.88 | 0.82 | 1.00 | 0.02 | 1.00 | 0.95 |
| `pls_within_group__net_displacement` | 0.88 | 1.00 | 1.00 | 1.00 | 0.94 | 1.00 | 0.86 |
| `pca__pc1` | 0.02 | 0.85 | 0.84 | 1.00 | 0.02 | 1.00 | 0.94 |
| `pca__net_displacement` | 0.09 | 0.96 | 1.00 | 1.00 | 0.20 | 1.00 | 0.04 |
| `concat__pc1` | 0.04 | 0.85 | 0.80 | 1.00 | 0.03 | 1.00 | 0.95 |
| `concat__net_displacement` | 0.10 | 1.00 | 1.00 | 1.00 | 0.21 | 1.00 | 0.05 |

Anchor `delta` and `shape` are at or below 0.04 in every arm at both designs.

## 4. Where the rotation is lost

From `checkpoint_table.csv`: median realized `angle` at orientation e = 1.0, and the median retained rank.

| arm | Phase 6 latent | Phase 6 rank | Phase 5 latent | Phase 5 rank |
|---|---|---|---|---|
| `pls_cv__pc1` | 63.9° | 4 | 72.0° | 3 |
| `pls_cv__net_displacement` | 56.7° | 4 | 6.1° | 3 |
| `pls_forced__pc1` | 71.1° | 4 | 89.5° | 6 |
| `pca__pc1` | 90.5° | 4 | 90.2° | 6 |
| `pls_within_group__pc1` | 90.3° | 4 | 90.2° | 6 |

The population angle is 90.0° and the observed standardized angle is 89.7–90.3° in every arm.

- **Label-free and two-sided spaces keep the full rotation.** `pca` and `pls_within_group` reproduce the observed angle. The CV stage-PLS space shrinks it.
- **Phase 5 under net displacement.** The `pls_cv` latent angle is only 6.1°, yet the test still rejects in 0.89 of replicates, because its anchor median is 0.3°. In a rank-3 = K − 1 space shared by both groups, the two net displacements are nearly parallel, and a small, consistent rotation is still resolvable.

## 5. Advancement decision

`advancement_decision.json`, rule version 1. The anchor bound is 0.05 + 2·√(0.05·0.95/100) = 0.0936. The other
criteria use 2 pooled SE against `pls_cv__pc1`.

| arm | failed criteria |
|---|---|
| `pls_cv__net_displacement` | anchor (Phase 6 `angle` 0.100 > 0.0936; Phase 5 0.080 passes) |
| `pls_forced__pc1` | gain, protected, specificity |
| `pls_forced__net_displacement` | protected (Phase 6 shape 0.01 vs 0.19; Phase 5 shape 0.65 vs 1.00), specificity (magnitude/`angle` 0.15 and 0.32 vs 0.01) |
| `pca__pc1` | protected (Phase 6 shape 0.00 vs 0.19), specificity (Phase 6 shape/`angle` 0.39 vs 0.10) |
| `pca__net_displacement` | anchor, protected, specificity |
| `pls_within_group__pc1` | anchor, gain, protected, specificity |
| `pls_within_group__net_displacement` | anchor, protected, specificity |
| `concat__*` | reference-only by rule |

**Advancing set: empty.**

## 6. Findings

1. **The estimator, not the space, is the main lever.**
   - In the production `pls_cv` space, net displacement raises Phase 6 orientation power from 0.29 to 0.67. That is +0.38, against a 2-SE margin of 0.13.
   - It leaves `delta` and `shape` power unchanged at both designs.
   - It removes the largest specificity failure the probe saw: shape → `angle` at Phase 5, from 0.89 to 0.02.
   - This matches the mechanism in §1. The signed PC1 angle's heavy null tail comes from PC1 swapping or flipping on near-isotropic configurations, and net displacement does not have that failure mode.
2. **Net displacement shows an anchor excess at both designs, in every space.**
   - Anchor `angle` rates are 0.08–0.10 for every non-`pls_within_group` space, against 0.01–0.05 under PC1.
   - The four spaces share datasets and permutation draws, so their rates are not independent. In effect, the evidence is about 9/100 per design, where 5/100 is expected.
   - That is suggestive at 100 replicates, but not established. The `pls_cv` arm misses the bound by less than one rejection.
   - The mechanism is not identified. One hypothesis, not tested here: the estimator uses only the first and last stages. RRPP permutes reduced-model residuals across all stages, so if the endpoint stages are noisier than the interior ones, the permutation null would be too narrow.
3. **A fixed `pls_within_group` basis is not a valid test space without a refit.** Anchor `angle` reaches 0.48–0.88. The basis is fitted to the observed groups' own stage directions, and so it carries their noise. The permuted data are then projected onto that fixed basis, so the null never sees the fitting step. This is the risk design D4 anticipated, and the anchor criterion caught it.
4. **Rank 2(K − 1) recovers orientation but costs shape.** At Phase 6, both `pca` and `pls_forced` drop shape/`shape` from 0.19 to 0.00–0.01. They also leak `angle` into the shape cell (PC1) or the magnitude cell (net displacement). Giving the space more dimensions is not the fix.
5. **`concat`, the ceiling reference, behaves like the label-free spaces.** It reaches power 1.00 under net displacement and has the same anchor excess. The p-padding question remains open.

## 7. Recommendation for the follow-up study

1. **Lead candidate: `pls_cv__net_displacement`, conditional on Type I.** First run an anchor-only Type I study for every net-displacement arm at both designs, at paper-grade size (500 replicates × 999 permutations). If the excess holds, diagnose it before running any power comparison. The endpoint-heteroscedasticity hypothesis in §6 is the first thing to test. If size holds, compare `pls_cv__net_displacement` against `pls_cv__pc1` at paper grade.
2. **Drop fixed-basis `pls_within_group`.** A refit-per-permutation variant would be valid by construction, but it multiplies cost by the permutation count. Pursue it only if net displacement fails.
3. **Do not carry `pca` or `pls_forced` as replacements.** Their Phase 6 shape loss is disqualifying for a three-statistic claim.
4. **Keep `concat` as the ceiling.** Decide p-padding in the follow-up.
5. **No production change.** The PC1 estimator and the CV stage-PLS space remain the defaults. The Phase 6 interpretability statement for the case study stands unchanged.
