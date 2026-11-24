# Stage-effect sample-size scaling — report (2026-11-24)

**Change:** `stage-effect-sample-size-scaling` (MOTCO 0.25.0, no package change; roadmap open item 14). **Reads:**
`results/stage-effect-scaling-2026-11-24/` (`stage_effect_scaling_decision.json`; per profile `settings.json`,
`rate_table.csv`, `paired_table.csv`, `cell_table.csv`, `rank_table.csv`, `accounting.json`;
`saturated/reproduction.json`; `calibration/`; `PROVENANCE.txt`, `sacct.txt`). **Parent:** the
[univariate baseline](univariate-baseline-2026-11-23.md)'s `stage_span` profile (its settings, rate table and cell
table, bound by sha256). **Concept page:**
[measurement visibility](../measurement-visibility.md#a-weak-stage-effect-against-sample-size).

At a weak stage effect the CV rank inflates, `delta` and `shape` lose power, and a rotation is read as `shape`. This
run asks whether more samples buy that back. In the generator a design's information about the stage program grows
as δ√n. So, if a weak stage effect were only a smaller effective sample, half the stage effect with four times the
samples would read the same, cell by cell. Each new design is paired with a committed design at the same δ√n:

- δ = 1.0 and δ = 0.5 at n = 320 and 1,280 against δ = 2.0 at n = 80;
- δ = 1.0 at n = 800 against δ = 2.0 at n = 200;
- δ = 0.25 at n = 800 against δ = 0.5 at n = 200.

A fifth pair sits below the AUROC ceiling: a stage effect calibrated to the case studies' held-out AUROC (δ_c = 0.05
at n = 320) against δ_c / 2 at n = 1,280. The decision was fixed before any unit ran. It labels each of 53 pair-cells
by Fisher's exact test on the production omnibus, with the per-feature `fdr` arm as the control of the pairing.

**In one paragraph.** **The verdict is `mixed`**: 3 counted costs and 12 counted gains. A weak stage effect is not
only a smaller sample: for some plants it costs more than δ√n says, and for others less.

- **The costs are rotations, all at δ = 0.5, n = 1,280 against δ = 2.0, n = 80.** The smallest in-span rotation
  (e = 0.1) reads 0.250 against 0.820. The two out-of-span rotations, which production sees only through incidental
  routes, read 0.060 against 0.255 and 0.375 against 0.835.
- **The gains are stretches and bends, and rotations below δ = 0.5.**
    - Stretches and small bends at δ = 1.0 gain 0.105–0.300 over their δ = 2.0 matches.
    - Small rotations at δ = 0.25 gain over δ = 0.5: 0.310 against 0.090, and 0.765 against 0.410.
    - At the calibrated weak signal, the n = 1,280 design gains on every stretch, for example 0.805 against 0.330 at
      e = 1.
    - The same δ = 0.5 design that loses the e = 0.1 rotation gains the e = 0.25 one (0.985 against 0.855).
- **The control moved once.** `fdr` carries no label in 52 of 53 pair-cells, and its one label (a gain on an
  out-of-span rotation at δ = 0.25) is not one production carries. No label was set aside.
- **Sample size does not restore the naming or the rank.**
    - The δ = 2.0 designs grade `names_kind` at n = 80, 200 and 500. Every weaker design grades `detects_only`, the
      matched ones included. At δ = 0.5, n = 1,280, a 45° rotation is read as `shape` in 0.720 of units.
    - The CV rank does not come down with more samples. At the AUROC ceiling it is set by δ: 3 at δ = 0.5 at both
      n = 200 and 1,280. Below the ceiling, n = 1,280 drives it to the cap of 19 at both δ_c / 2 and δ_c, against a
      median of 10 at δ_c, n = 320.
- **Size holds.** Production, `fdr` and `total` hold size in both families. The rerun of δ = 0.5, n = 200 reproduces
  all 234 committed rows exactly. 24,400 of 24,400 units completed.
- **What it holds for.** The reading holds at 658 simulated features with Gaussian noise, equal groups and three
  stages. The band designs match the case studies' held-out AUROC, not their cell tables, covariates or rank, and no
  real-noise reading exists on methylation, expression and proteomics (roadmap open item 2). The decision was
  written after a timing probe of two null datasets per design, at a seed outside the run, which saw no rejection
  rate. No verdict changes what `motco analyze` runs.

## 1. Configuration and provenance

- **Profile `saturated`: four matched designs and the gate.** Methylation, expression and proteomics (658 features),
  three stages, equal groups, base seed 2043, the parent's thirteen cells:
    - the anchor and a translation at 500 replicates;
    - stretches at e = 0.05, 0.1 and 0.25, and in-span rotations (φ = 0) and bends at e = 0.1, 0.25 and 0.5, under the
      planted scale, each at 200;
    - the out-of-span rotation at e = 0.5, φ = 90 and its 20-feature form, under the null scale, each at 200.

  The new designs are δ = 1.0 at n = 320, δ = 0.5 at n = 1,280, δ = 1.0 at n = 800 and δ = 0.25 at n = 800. Each
  differs from its committed match only in the per-omic stage effects and n. The parent's δ = 0.5, n = 200 design is
  rerun at its own name and seeds as the gate. 16,000 units.
- **Profile `band`: three designs below the AUROC ceiling.**
    - `calibrate` fitted production on 40 anchor datasets per δ ∈ {0.03, …, 0.1} at n = 320 (280 fits, §8.5) and
      adopted δ_c = 0.05, whose median held-out AUROC of 0.664 is nearest the target 0.636 and inside 0.60–0.67.
    - The designs are δ_c at n = 320, δ_c / 2 = 0.025 at n = 1,280 (its δ√n match) and δ_c at n = 1,280.
    - Eleven cells under the planted scale: the anchor and a translation at 500, stretches at e = 0.5, 1 and 2, and
      rotations and bends at e = 0.25, 0.5 and 1 at 200. The out-of-span cells are left out, because production is
      blind there at any stage effect.
    - 8,400 units.
- **The unit is the parent's, imported unchanged.** One production fit (`pls_cv`, rank floor 2, net displacement, 999
  permutations). The judged arms are `production` (the smallest of `delta`, `angle` and `shape` at 0.05 / 3), `fdr`
  (Benjamini–Hochberg on the per-feature interaction F) and `total` (the split test's sum). `scan`, the three single
  statistics, `span`, `complement`, `production_minp` and `pca_matched` are recorded and never judged.
- **Code and host.** Revision `3afc721`, pushed before any array unit ran, from a clean detached worktree on `ing`.
  Python 3.11.16, numpy 2.3.5, scipy 1.16.3, scikit-learn 1.8.0, pandas 2.3.3, BLAS pinned to one thread.
- **Where it ran.** Arrays 958096 (`saturated`, 200 tasks) and 958097 (`band`, 140), partition `512x1024` with n11
  excluded: 340 of 340 `COMPLETED 0:0`, 11:18–15:25.
    - **409.49 recorded core-hours** (206.09 / 203.40) against the approved cap of 600; the local smoke had projected
      480–540.
    - The median unit took 34.7 s in `saturated` and 113.8 s in `band`, most of it the cross-validated fit.
    - 65 rerun-design units stalled on the hosts for 10–47 minutes, across every stage of the unit. They hold
      38.8 of the core-hours, and their counts reproduce exactly.

## 2. Accounting and the reproduction

- **Units.** 16,000 and 8,400 complete: no failed unit, no duplicate, no failed planted fixed point (the band's e = 2
  stretch at δ = 0.025 included).
- **Reproduction: passed.** All 234 rows of the rerun design (available and rejected counts per cell and arm) equal
  the committed `stage_span` rate table, none unresolved. The index row of the rerun (rank 3, AUROC 1.000, capture
  0.905) equals the committed one.
- **The decision regenerates.** `decide` on the committed aggregates alone, on the workstation, writes a decision
  byte-identical to the cluster's (sha256 `c5cd1310…`).

## 3. The decision

The order is gate, size, then the matched comparison.

- **Size.** Judged per arm in two families: the four new `saturated` designs (m = 8) and the three `band` designs
  (m = 6). Each family uses the count-plus-ceiling rule (`judge_rates`), with a per-cell bound of 0.0695 at 500
  replicates.
- **Labels.** Each pair-cell compares `production`'s rejected counts in the two designs by Fisher's exact test, one-sided
  each way, at 0.025 / 53 = 4.7 × 10⁻⁴:
    - `weak_costs`: the weaker-staged design's rate is lower by at least 0.10;
    - `weak_gains`: it is higher by at least 0.10;
    - `unresolved`: otherwise.

  The `saturated` matches' counts are the committed ones. `fdr` is labelled the same way, and a `production` label
  counts unless `fdr` carries the same label there.
- **Verdict.** By existence over the counted labels. With 200 replicates per design the difference of two rates has a
  standard error up to 0.05, so `unresolved` is not agreement.
- **Same section, same limits:** the reading holds at 658 simulated features with Gaussian noise, equal groups and
  three stages. The rule was written after a timing probe of null datasets (design, Context), which saw the rank,
  AUROC, capture and timing of two units per design and no rejection rate.

| criterion | reading |
|---|---|
| gate | reproduction passed (234 rows); 0 failed units |
| size | `production`, `fdr` and `total` hold in both families |
| labels | 3 counted `weak_costs`, 12 counted `weak_gains`, 38 `unresolved`; 0 shared with `fdr`, 0 unavailable |
| verdict | **`mixed`** |
| grades | `names_kind` at the three δ = 2.0 designs; `detects_only` at the other eight |

## 4. Size per arm

Rejections of 500 on the null cells (anchor; translation):

| profile | design | δ√n | production | `fdr` | `total` |
|---|---|---|---|---|---|
| `saturated` | δ = 1.0, n = 320 | 17.9 | 24; 19 | 29; 29 | 27; 28 |
| `saturated` | δ = 0.5, n = 1,280 | 17.9 | 11; 11 | 11; 11 | 17; 18 |
| `saturated` | δ = 1.0, n = 800 | 28.3 | 26; 20 | 16; 16 | 32; **35** |
| `saturated` | δ = 0.25, n = 800 | 7.1 | 16; 21 | 17; 17 | 24; 25 |
| `band` | δ = 0.05, n = 320 | 0.89 | 12; 8 | 14; 14 | 24; 25 |
| `band` | δ = 0.025, n = 1,280 | 0.89 | 21; 19 | 19; 19 | 26; 26 |
| `band` | δ = 0.05, n = 1,280 | 1.79 | 20; 24 | 20; 20 | **37**; 29 |

- **Every arm holds.** `total` exceeds the bound once in each family (35 / 500 = 0.070; 37 / 500 = 0.074, at the
  band's ceiling and not above it). Each family allows one exceedance.
- **Production's largest null rate is 0.052.** `fdr`'s is 0.058.
- **The rerun design is outside the families.** It reads 13; 15 for production, 20; 20 for `fdr` and 30; 30 for
  `total`, as committed.

## 5. The matched comparison

Labels by pair and plant family (costs / gains / unresolved):

| pair (weaker ~ match) | δ√n | stretch | in-span rotation | bend | out-of-span rotation |
|---|---|---|---|---|---|
| δ = 1.0, n = 320 ~ δ = 2.0, n = 80 | 17.9 | 0 / 2 / 1 | 0 / 0 / 3 | 0 / 2 / 1 | 0 / 0 / 2 |
| δ = 0.5, n = 1,280 ~ δ = 2.0, n = 80 | 17.9 | 0 / 0 / 3 | 1 / 1 / 1 | 0 / 0 / 3 | 2 / 0 / 0 |
| δ = 1.0, n = 800 ~ δ = 2.0, n = 200 | 28.3 | 0 / 1 / 2 | 0 / 0 / 3 | 0 / 1 / 2 | 0 / 0 / 2 |
| δ = 0.25, n = 800 ~ δ = 0.5, n = 200 | 7.1 | 0 / 0 / 3 | 0 / 2 / 1 | 0 / 0 / 3 | 0 / 0 / 2 |
| δ = 0.025, n = 1,280 ~ δ = 0.05, n = 320 | 0.89 | 0 / 3 / 0 | 0 / 0 / 3 | 0 / 0 / 3 | — |
| **all pairs** | | 0 / 6 / 9 | 1 / 3 / 11 | 0 / 3 / 12 | 2 / 0 / 6 |

**The control.** `fdr` is `unresolved` in 52 of the 53 pair-cells. Its one label is a `weak_gains` on the
out-of-span rotation at δ = 0.25 against δ = 0.5 (0.970 against 0.840), where production is unresolved (0.030 against
0.035). That is the direction the design expected for the per-feature test, which has more residual degrees of
freedom at the larger n. In 39 of the 53 pair-cells `fdr` has little room to move. It reads 0.97 or more on both sides
in 25: every out-of-span cell but that one, every in-span rotation at δ ≥ 0.5, and the largest stretches and bends.
It reads 0.10 or less on both sides in 14.

### 5.1 Costs (3 counted)

| pair | cell | weaker | match | difference | p (costs) | `fdr` weaker; match |
|---|---|---|---|---|---|---|
| δ = 0.5, n = 1,280 ~ δ = 2.0, n = 80 | rotation e = 0.1, φ = 0 | 50 / 200 (0.250) | 164 / 200 (0.820) | −0.570 | 8.1 × 10⁻³² | 0.990; 0.970 |
| δ = 0.5, n = 1,280 ~ δ = 2.0, n = 80 | rotation out of span | 12 / 200 (0.060) | 51 / 200 (0.255) | −0.195 | 4.0 × 10⁻⁸ | 1.000; 1.000 |
| δ = 0.5, n = 1,280 ~ δ = 2.0, n = 80 | rotation out of span, 20 features | 75 / 200 (0.375) | 167 / 200 (0.835) | −0.460 | 7.9 × 10⁻²² | 1.000; 1.000 |

- **The rotation cost grows as δ falls along the diagonal.** At δ√n = 17.9, the e = 0.1 rotation reads 0.820 at
  δ = 2.0, 0.690 at δ = 1.0 (−0.130, p = 0.0018, above the level: `unresolved`) and 0.250 at δ = 0.5.
- **It is confined to the smallest rotation.** At e = 0.25 the same δ = 0.5 design gains (0.985 against 0.855), and at
  e = 0.5 both read 1.000 and 0.980. The δ = 2.0 design's rotation rates are flat across e (0.820, 0.855, 0.980). The
  δ = 0.5 design's are steep: from 0.250 to 0.985 between e = 0.1 and 0.25.
- **The out-of-span costs are on cells production cannot see.** Production reaches them only through incidental
  routes (measurement visibility). The routes are weaker at δ = 0.5 than at δ = 2.0, at the same δ√n. These cells
  are under the null scale, so the plant is fixed relative to the null pooled SD, and that SD includes the
  between-stage variance, which grows as δ². That is one of the places where the design named the scaling as
  approximate; the run does not separate it from the others.

### 5.2 Gains (12 counted)

| pair | cell | weaker | match | difference | `fdr` weaker; match |
|---|---|---|---|---|---|
| δ = 1.0, n = 320 ~ δ = 2.0, n = 80 | stretch e = 0.05 | 0.625 | 0.345 | +0.280 | 0.070; 0.060 |
| | stretch e = 0.1 | 0.960 | 0.750 | +0.210 | 0.330; 0.310 |
| | bend e = 0.1 | 0.665 | 0.365 | +0.300 | 0.100; 0.105 |
| | bend e = 0.25 | 0.985 | 0.880 | +0.105 | 0.940; 0.885 |
| δ = 0.5, n = 1,280 ~ δ = 2.0, n = 80 | rotation e = 0.25, φ = 0 | 0.985 | 0.855 | +0.130 | 1.000; 1.000 |
| δ = 1.0, n = 800 ~ δ = 2.0, n = 200 | stretch e = 0.05 | 0.975 | 0.815 | +0.160 | 0.205; 0.235 |
| | bend e = 0.1 | 0.990 | 0.860 | +0.130 | 0.390; 0.420 |
| δ = 0.25, n = 800 ~ δ = 0.5, n = 200 | rotation e = 0.1, φ = 0 | 0.310 | 0.090 | +0.220 | 0.095; 0.095 |
| | rotation e = 0.25, φ = 0 | 0.765 | 0.410 | +0.355 | 0.980; 0.995 |
| δ = 0.025, n = 1,280 ~ δ = 0.05, n = 320 | stretch e = 0.5 | 0.330 | 0.070 | +0.260 | 0.040; 0.035 |
| | stretch e = 1 | 0.805 | 0.330 | +0.475 | 0.060; 0.045 |
| | stretch e = 2 | 1.000 | 0.835 | +0.165 | 0.415; 0.350 |

Every gain p-value is at most 1.4 × 10⁻⁵, against the level of 4.7 × 10⁻⁴.

- **Stretches and bends gain at δ = 1.0 and do not cost at δ = 0.5.** On the δ = 0.5, n = 1,280 design they match
  their δ = 2.0, n = 80 references within 0.055 (stretches 0.380 / 0.750 / 1.000 against 0.345 / 0.750 / 0.990).
- **Every gain is a rate `fdr` does not share.** On those cells `fdr` moves by at most 0.065 between the two designs.
- **At the calibrated weak signal the matched design reads more on every stretch.** The two designs share δ√n and
  nearly their AUROC (0.649 against 0.664), but the n = 1,280 design selects rank 19 in every unit, against a median
  of 10 at n = 320 (§8.1). Its fitted span holds more of the stage program (capture 0.376 against 0.288).

### 5.3 Unresolved (38)

- **Saturated.** 16 of the 38 have production at 0.88 or more in both designs.
- **Near the level.** The band's six rotations and bends read 0.08 or less in both designs (§8.4).
- **Close.** The other 16 differ by less than 0.10, or do not reach the level:
    - the δ = 0.25 pair's stretches and bends within 0.025, except its bend at e = 0.5;
    - the δ = 0.5, n = 1,280 design's two small stretches and its e = 0.1 bend within 0.055 of δ = 2.0, n = 80;
    - the δ = 1.0, n = 320 design's two small rotations;
    - five out-of-span rotations, where production reads 0.03–0.87 on both sides.
- **Three differences meet only one criterion.**
    - The δ = 0.25 pair's 45° rotation differs by +0.090 at p = 1.3 × 10⁻⁵: significant, but below the margin.
    - Its bend at e = 0.5 differs by −0.145 at p = 0.0016: above the margin, but not significant at the level.
    - The δ = 1.0, n = 320 rotation at e = 0.1 differs by −0.130 at p = 0.0018 (§5.1).

Every rate is in the decision's `criteria.pair_cells`.

## 6. The grade by design

Per design, among the in-span cells that production detects in at least 0.80 of units, the two off-target statistics
of each form one family under `judge_rates` (bound 0.081 at 200 replicates). The committed designs are graded from the
committed table by the same rule.

| design | δ√n | detected cells (of 9) | off-target readings above the bound (above the ceiling) | grade |
|---|---|---|---|---|
| δ = 2.0, n = 80 | 17.9 | 6 | 0 | `names_kind` |
| δ = 2.0, n = 200 | 28.3 | 9 | 1 (0) | `names_kind` |
| δ = 2.0, n = 500 | 44.7 | 9 | 0 | `names_kind` |
| δ = 0.5, n = 200 | 7.1 | 2 | 2 (2) | `detects_only` |
| δ = 1.0, n = 320 | 17.9 | 5 | 4 (3) | `detects_only` |
| δ = 0.5, n = 1,280 | 17.9 | 5 | 6 (6) | `detects_only` |
| δ = 1.0, n = 800 | 28.3 | 9 | 2 (1) | `detects_only` |
| δ = 0.25, n = 800 | 7.1 | 1 | 2 (2) | `detects_only` |
| δ = 0.05, n = 320 | 0.89 | 1 | 1 (1) | `detects_only` |
| δ = 0.025, n = 1,280 | 0.89 | 2 | 2 (2) | `detects_only` |
| δ = 0.05, n = 1,280 | 1.79 | 4 | 6 (6) | `detects_only` |

The readings above the ceiling:

- **Rotations read as `shape`, and as `delta`, as δ falls.**

  | design | rotation | `shape` | `delta` |
  |---|---|---|---|
  | δ = 1.0, n = 320 | e = 0.5 | 0.175 | |
  | δ = 1.0, n = 800 | e = 0.5 | 0.140 | |
  | δ = 0.5, n = 1,280 | e = 0.25 | 0.540 | 0.175 |
  | δ = 0.5, n = 1,280 | e = 0.5 | 0.720 | 0.340 |
  | δ = 0.25, n = 800 | e = 0.5 | 0.930 | 0.530 |
  | δ = 0.05, n = 1,280 | e = 1 | 0.535 | 0.110 |

  The committed δ = 0.5, n = 200 design reads 0.665 and 0.155 on its e = 0.5 rotation.
- **Bends read as `delta` at δ = 0.5.** 0.140 at e = 0.25 and 0.305 at e = 0.5, at n = 1,280.
- **Stretches and bends read as `angle`.**
    - At δ = 1.0, n = 320: 0.100–0.115 on the bends.
    - At the band's large stretches: 0.215 at e = 2 and n = 320, and 0.580 at e = 2 and δ_c, n = 1,280, where `shape`
      reads 0.185.

**Matched samples do not restore the naming.** At δ√n = 17.9, δ = 2.0 grades `names_kind` and both weaker designs
`detects_only`, and at 28.3 so does δ = 1.0, n = 800 against δ = 2.0, n = 200. The leak grows as δ falls at fixed δ√n
(on the 45° rotation, `shape` reads 0.175 at δ = 1.0 and 0.720 at δ = 0.5). It also grows at fixed δ as n rises:
0.665 at δ = 0.5, n = 200, against 0.720 at n = 1,280. The grade never changes the verdict.

## 7. The map

The parent's labels, by paired exact tests of production against the admissible per-feature arms (`fdr`, `total`) at
its levels for m = 71 and its 0.10 margin: production leads / the per-feature test leads / unresolved.

| design | stretch | in-span rotation | bend | out-of-span rotation |
|---|---|---|---|---|
| δ = 1.0, n = 320 | 2 / 0 / 1 | 0 / 2 / 1 | 1 / 0 / 2 | 0 / 2 / 0 |
| δ = 0.5, n = 1,280 | 2 / 0 / 1 | 0 / 1 / 2 | 1 / 0 / 2 | 0 / 2 / 0 |
| δ = 1.0, n = 800 | 1 / 0 / 2 | 0 / 0 / 3 | 1 / 0 / 2 | 0 / 1 / 1 |
| δ = 0.25, n = 800 | 2 / 0 / 1 | 1 / 1 / 1 | 2 / 0 / 1 | 0 / 2 / 0 |
| δ = 0.05, n = 320 | 2 / 0 / 1 | 0 / 0 / 3 | 0 / 0 / 3 | — |
| δ = 0.025, n = 1,280 | 3 / 0 / 0 | 0 / 0 / 3 | 0 / 0 / 3 | — |
| δ = 0.05, n = 1,280 | 2 / 0 / 1 | 1 / 0 / 2 | 2 / 0 / 1 | — |
| **all** | 14 / 0 / 7 | 2 / 4 / 15 | 7 / 0 / 14 | 0 / 7 / 1 |

23 cells where production leads, 11 where the per-feature test leads, 37 unresolved. **The parent's map holds where the stage fit is
weak.**

- **Production leads on stretches and bends at every design, the band included.** For example, at δ_c, n = 1,280 a
  stretch at e = 0.5 reads 0.875 against `fdr`'s 0.080, and a bend at e = 1 reads 0.505 against 0.075.
- **The per-feature test leads on every out-of-span rotation but one** (n = 800 at δ = 1.0, 20 features: both sides
  at 1.000).
- **In-span rotations go to the per-feature test where the fit is weakest.**
    - At δ = 1.0, n = 320, rotations at e = 0.1 and 0.25 read 0.690 and 0.765 against `fdr`'s 0.985 and 1.000.
    - At δ = 0.5, n = 1,280, the e = 0.1 rotation reads 0.250 against 0.990.
    - At δ = 0.25, n = 800, the e = 0.25 rotation reads 0.765 against 0.980.
- **Production leads on two rotations:** at δ = 0.25, e = 0.1 (0.310 against 0.095), and at δ_c, n = 1,280, e = 0.5
  (0.600 against 0.200).
- **The band has no baseline lead.** Its in-span rotations and bends at n = 320 and at the matched n = 1,280 read
  0.080 or less for production and 0.130 or less for `fdr`.

## 8. Descriptive readings

### 8.1 The index: where a cohort sits

| design | source | n | δ | δ√n | held-out AUROC | CV rank median (q10–q90) | pooled-stage capture |
|---|---|---|---|---|---|---|---|
| δ = 2.0, n = 80 | committed | 80 | 2.0 | 17.9 | 1.000 | 2 (2–2) | 0.934 |
| δ = 2.0, n = 200 | committed | 200 | 2.0 | 28.3 | 1.000 | 2 (2–2) | 0.973 |
| δ = 2.0, n = 500 | committed | 500 | 2.0 | 44.7 | 1.000 | 2 (2–2) | 0.989 |
| δ = 0.5, n = 200 | committed; rerun equal | 200 | 0.5 | 7.1 | 1.000 | 3 (3–3) | 0.905 |
| δ = 1.0, n = 320 | this run | 320 | 1.0 | 17.9 | 1.000 | 2 (2–3) | 0.971 |
| δ = 0.5, n = 1,280 | this run | 1,280 | 0.5 | 17.9 | 1.000 | 3 (3–3) | 0.984 |
| δ = 1.0, n = 800 | this run | 800 | 1.0 | 28.3 | 1.000 | 2 (2–2) | 0.987 |
| δ = 0.25, n = 800 | this run | 800 | 0.25 | 7.1 | 1.000 | 5 (4–6) | 0.934 |
| δ = 0.05, n = 320 | this run | 320 | 0.05 | 0.89 | 0.664 | 10 (5–19) | 0.288 |
| δ = 0.025, n = 1,280 | this run | 1,280 | 0.025 | 0.89 | 0.649 | 19 (19–19) | 0.376 |
| δ = 0.05, n = 1,280 | this run | 1,280 | 0.05 | 1.79 | 0.883 | 19 (19–19) | 0.647 |

The anchor cell's medians. A cohort can be placed by its held-out AUROC and CV rank:

- **At the ceiling,** the AUROC cannot tell the designs apart, and the rank indexes δ.
- **Below it,** the AUROC follows δ√n: 0.649 against 0.664 at δ√n = 0.89, and 0.883 at 1.79.
- **The rank does not follow δ√n.** At the calibrated signal, four times the samples drives it to the cap of 19, where
  ADNI's design simulation also sat (rank 19 at an AUROC of 0.606, roadmap item 9).
- **The capture rises with n at matched δ√n in every pair.**

### 8.2 Stretches against e·δ·√n

Production's rate on every stretch cell, by the index the generator's scaling predicts it to follow:

| e·δ·√n | rates (design) |
|---|---|
| 0.35 | 0.075 (δ = 0.5, n = 200); 0.100 (δ = 0.25, n = 800) |
| 0.45 | 0.070 (δ_c, n = 320); 0.330 (δ_c / 2, n = 1,280) |
| 0.71 | 0.285 (δ = 0.5, n = 200); 0.290 (δ = 0.25, n = 800) |
| 0.89 | 0.345 (δ = 2.0, n = 80); 0.625 (δ = 1.0, n = 320); 0.380 (δ = 0.5, n = 1,280); 0.330 (δ_c, n = 320); 0.805 (δ_c / 2, n = 1,280); 0.875 (δ_c, n = 1,280) |
| 1.41 | 0.815 (δ = 2.0, n = 200); 0.975 (δ = 1.0, n = 800) |
| 1.77–1.79 | 0.810 (δ = 0.5, n = 200); 0.795 (δ = 0.25, n = 800); 0.750 (δ = 2.0, n = 80); 0.960 (δ = 1.0, n = 320); 0.750 (δ = 0.5, n = 1,280); 0.835 (δ_c, n = 320); 1.000 (δ_c / 2 and δ_c, n = 1,280) |
| ≥ 2.2 | 0.985–1.000 |

- **Where the index lines up.** It orders the rates within each design. It lines up the δ = 0.5 / 2.0 pair and the
  δ = 0.25 / 0.5 pair within 0.035.
- **Where it does not.** The δ = 1.0 designs read above both their δ = 2.0 and δ = 0.5 neighbours, and at the
  calibrated signal the n = 1,280 designs read far above n = 320.
- **What changed from the roadmap.** Item 24's four designs had suggested that stretches line up on e·δ·√n. Across
  these eleven designs they do so only between some pairs.
- **The other statistics.** `delta` alone reads 0.025–0.160 above production on every stretch where production is
  below 0.95, and `fdr` reads at most 0.420 on every stretch at an index of 1.79 or less.

### 8.3 Rank strata of the costs

The weaker design's rates by its selected rank, on each pair-cell labelled `weak_costs`:

| cell | rank 2: production (units) | rank 3: production (units) | `fdr` at rank 3 |
|---|---|---|---|
| rotation e = 0.1, φ = 0 | 1.000 (10) | 0.211 (190) | 0.989 |
| rotation out of span | 0.200 (5) | 0.056 (195) | 1.000 |
| rotation out of span, 20 features | 1.000 (5) | 0.359 (195) | 1.000 |

The δ = 0.5, n = 1,280 design selects rank 3 in 95–97.5 % of these units. Its δ = 2.0 match selects rank 2 in at least
90 % of the units of each of these cells (the committed cell table's rank q90 is 2). Production's rejections concentrate in the few rank-2 units in two of the three cells. Those strata are 5 and
10 units chosen by the data, so they describe and do not attribute.

### 8.4 What is detectable at the case studies' AUROC

Production's rates on the band's signal cells (`fdr` in brackets):

| cell | δ_c, n = 320 | δ_c / 2, n = 1,280 | δ_c, n = 1,280 |
|---|---|---|---|
| stretch e = 0.5 | 0.070 (0.035) | 0.330 (0.040) | 0.875 (0.080) |
| stretch e = 1 | 0.330 (0.045) | 0.805 (0.060) | 1.000 (0.420) |
| stretch e = 2 | 0.835 (0.350) | 1.000 (0.415) | 1.000 (1.000) |
| rotation e = 0.25 / 0.5 / 1 | 0.025 / 0.035 / 0.060 (0.040 / 0.040 / 0.125) | 0.050 / 0.055 / 0.080 (0.035 / 0.060 / 0.130) | 0.150 / 0.600 / 0.920 (0.070 / 0.200 / 0.920) |
| bend e = 0.25 / 0.5 / 1 | 0.015 / 0.010 / 0.020 (0.040) | 0.045 / 0.050 / 0.060 (0.040) | 0.125 / 0.295 / 0.505 (0.055 / 0.060 / 0.075) |

- **At δ_c, n = 320 (AUROC 0.664), only the largest stretch reaches 0.80** (e = 2, which triples group B's stage
  configuration). On every rotation and bend up to e = 1, production reads at most 0.060 and `fdr` at most 0.125.
- **At the matched δ_c / 2, n = 1,280,** the stretches gain (§5.2) and the rotations and bends stay near the level.
- **At δ_c with four times the samples (AUROC 0.883),** a rotation at e = 1 reaches 0.920 and a bend 0.505. There
  the rotation is read as `shape` in 0.535 of units (§6).

### 8.5 The calibration

| δ | median held-out AUROC (q25–q75) | CV rank q10 / median / q90 |
|---|---|---|
| 0.03 | 0.549 (0.531–0.570) | 2 / 3 / 13 |
| 0.04 | 0.600 (0.587–0.620) | 2.9 / 8.5 / 19 |
| **0.05** | **0.664 (0.636–0.681)** | 5 / 9 / 19 |
| 0.06 | 0.728 (0.702–0.745) | 8 / 15.5 / 19 |
| 0.07 | 0.789 (0.764–0.808) | 8 / 11.5 / 19 |
| 0.08 | 0.847 (0.821–0.859) | 9.9 / 12.5 / 19 |
| 0.1 | 0.930 (0.910–0.937) | 9 / 13 / 19 |

40 anchor fits per value at n = 320, no permutation, 0 failures. δ_c = 0.05 is 0.028 from the target 0.636. δ = 0.04
(0.600) sits on the band's lower edge, 0.036 from the target.

### 8.6 The recorded arms

- **`production_minp`.** It is never below the omnibus and at most 0.050 above it, in all 98 cells.
- **`pca_matched`.** PCA at the production rank never exceeds production by 0.10 on a signal cell. In the band it runs
  above its level on the nulls: 0.064–0.072 at n = 320 and 0.094–0.106 at n = 1,280. Its band rates are therefore not
  powers.
- **`span` and `complement`.**
    - `complement` reads near its level on the `saturated` stretches and bends (0.030–0.085), 0.155–0.490 on the 45°
      rotations, and 0.99–1.00 on every out-of-span rotation.
    - In the band, where the fitted span holds 0.29–0.65 of the stage program, it reads 0.245–0.985 on the rotation
      at e = 1.

## 9. What follows

- **The documentation states the map of costs and gains by plant family** (design, Migration Plan, `mixed`). The
  [measurement-visibility](../measurement-visibility.md#a-weak-stage-effect-against-sample-size) section says it:
    - **Stretches and bends.** For a weakly staged design with a matched sample, power follows δ√n or does better.
    - **Small in-span rotations.** Below δ = 1.0 they cost: at δ = 0.5 against δ = 2.0, sample size does not buy back
      a rotation at e = 0.1.
    - **The incidental route to out-of-span rotations** weakens with δ at fixed δ√n.
    - **Power cannot be read off one number.** A cohort's power is indexed by its stage effect and its sample size
      separately, not by δ√n alone.
- **The naming is not bought back.** Every design below δ = 2.0 grades `detects_only`, at matched δ√n and at four
  times the samples. A weakly staged cohort's rejection names a group × stage difference; the statistic that rejects
  does not name its kind, whatever its sample size.
- **The rank does not come down with sample size.** At the AUROC ceiling δ sets it (3 at δ = 0.5 at both n = 200
  and 1,280; 5 at δ = 0.25). Below the ceiling, more samples raise it to the cap of 19.
- **The rank rule stays on hold.** The design reopens it under `weak_stage_costs` with the costs in the inflated-rank
  strata. The verdict is `mixed`. The strata of the one in-span cost are recorded (§8.3) and attribute nothing.
- **What does not change.** `motco analyze` runs what it ran; no default moves.
- **Follow-ups.**
    - The same reading on real noise belongs to the first multi-omic ROSMAP or ADNI export (roadmap open item 2).
    - Unequal group dispersion and feature count belong to roadmap open items 8 and 15.
    - The δ = 2.0, n = 500 match (δ = 1.0 at n = 2,000) was not run.

## 10. Limitations

- **Simulated features with Gaussian noise.** Every claim holds at 658 simulated features, equal groups and three
  stages, with one noise covariance for both groups. No real-noise reading exists (roadmap open item 2).
- **The decision followed a timing probe.** The probe ran two null datasets per design at base seed 90001, outside the
  run. It saw their rank, AUROC, capture and timing, and no rejection rate (design, Context).
- **The pairing is exact only in the generator's means.** The design named three places where the scaling is
  approximate:
    - the pooled SD's between-stage share, largest at δ = 2.0;
    - the residual degrees of freedom (74 at n = 80, 1,274 at n = 1,280);
    - the nonlinearity of rotations and bends in e.

  The run reads their sum. The `fdr` control moved once, but it is saturated or near its level in most pair-cells
  (§5), so it rules out little of the pairing's own drift.
- **Existence, at 200 replicates.** A rate difference below 0.10 is not excluded anywhere, and the 38 unresolved
  pair-cells are not agreements.
- **The band matches an AUROC, not a cohort.** It matches neither case study's cell table, covariates or rank. The
  calibration's δ_c sits near the band's upper edge.
- **The rank strata are conditional** on a rank the data choose, and the rank-2 strata hold 5–10 units.
- **One match was not run:** δ = 1.0 at n = 2,000, against δ = 2.0 at n = 500 (cost, and a saturated reference).
