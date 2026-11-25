# Group-specific dispersion — report (2026-11-24)

**Change:** `group-specific-dispersion` (MOTCO 0.26.0: the opt-in generator setting `group_dispersion_scale`;
roadmap open item 8). **Reads:** `results/group-dispersion-2026-11-24/` (`group_dispersion_decision.json`; per
profile `settings.json`, `rate_table.csv`, `cell_table.csv`, `within_group_table.csv`, `accounting.json`;
`core/reproduction.json`; `PROVENANCE.txt`, `sacct.txt`). **Parents:** the
[univariate baseline](univariate-baseline-2026-11-23.md)'s `stage_span` profile and the
[stage-effect scaling](stage-effect-scaling-2026-11-24.md) run's `band` profile (their settings, rate tables and cell
tables, bound by sha256). **Concept page:**
[measurement visibility](../measurement-visibility.md#unequal-group-dispersion-and-unbalance).

RRPP permutes the reduced model's residuals over all samples, so each draw gives both groups the pooled noise. The
null therefore assumes that the two groups share one noise dispersion. This run asks whether the trajectory test holds
its size when they do not. Every unit is a null unit: the groups share their trajectory and differ only in the spread
of their noise and in their size.

- **The setting.** s multiplies the second group's noise standard deviation in every feature, block and stage, about
  unchanged means: s ∈ {1.25, 1.5, 2}, covariance multipliers 1.56, 2.25 and 4.
- **The shares.** The first group's share is 0.5 (balanced), 0.65 or 0.8 (the wider group is the smaller), or 0.35
  or 0.2 (the wider group is the larger). Equal dispersion at shares 0.65 and 0.8 reads unbalance alone.
- **The decision**, fixed before any unit ran, judges `delta`, `angle` and `shape` in four families and takes the
  first that fails in a fixed order: unbalance at equal dispersion, then balanced designs with a scale, then
  unbalanced designs with a scale.

**In one paragraph.**

- **The verdict is `unbalance_excess`, on one reading.** At equal dispersion with unequal groups, 23 of the 24
  readings stay within their ceiling. The exception is `delta` at n = 80 with share 0.8, where the smaller group has
  5 samples per stage: 44 rejections of 500, against a ceiling of 40 (the committed balanced design reads 30).
  Production's omnibus holds in this family, and so do the per-feature arms.
- **What the verdict sets.** By the predeclared order, unbalance at equal dispersion is diagnosed first, as its own
  change. The three dispersion families are reported and not interpreted.
- **The dispersion families all fail, and far from the bound** (reported, §5).
    - **Balanced:** 20 of 42 readings are above the ceiling. `delta` is above it from s = 1.25 at n = 200 (86 and 93
      of 500) and from s = 1.5 at every other column (124–444). `angle` and `shape` stay within it at δ ≥ 0.5, apart
      from one reading each at s = 2, and fail in the weak-stage band.
    - **Wider group the smaller:** 48 of 60 are above. All three statistics fail by s = 1.5 at every column; at
      n = 200, share 0.8 and s = 2, `delta` / `angle` / `shape` read 412 / 94 / 118.
    - **Wider group the larger:** 19 of 60 are above: fifteen `delta` readings and four `angle` and `shape`
      readings in the band. At δ ≥ 0.5 `angle` and `shape` turn conservative (25 readings flagged).
- **The within-group draws (recorded)** bring `angle` and `shape` back where the wider group is the smaller. They
  leave balanced `delta` where it was (332 → 319 of 500 at n = 200, s = 2), and raise `delta` where the wider group
  is the larger (109 → 223).
- **The per-feature tests (recorded)** hold where the groups are balanced or the wider group is the larger. They fail
  where the wider group is the smaller: `fdr` reads 134 of 500 at n = 200, share 0.8 and s = 1.25, and 498 at s = 2.
- **The gate passed.** The rerun of `n200` reproduces all 36 committed rows exactly; 32,000 of 32,000 units
  completed; 186.8 recorded core-hours against the approved cap of 400.
- **What it holds for.** Simulated Gaussian noise at 658 features, one scale on every feature of one group, three
  stages and no covariates. A difference in correlation structure was not simulated, and no real cohort's dispersion
  was read. No verdict changes what `motco analyze` runs.

## 1. Configuration and provenance

- **Profile `core`: the full cross at n = 200, δ = 2.0.** Methylation, expression and proteomics (658 features),
  three stages, base seed 2043. The anchor and the translation at e = 1.0 under the planted scale, 500 replicates
  each.
    - **Seventeen new designs:** s ∈ {1.25, 1.5, 2} at shares 0.5, 0.65, 0.8, 0.35 and 0.2, and equal dispersion at
      shares 0.65 and 0.8. Each differs from the committed `n200` only in `group_ratio` and `group_dispersion_scale`.
    - **The gate:** the univariate baseline's `n200` controls, rerun at their committed name and seeds.
    - 18,000 units.
- **Profile `columns`: a reduced cross at four other committed designs.** The univariate baseline's `n80` and `n500`
  at δ = 2.0 and `n200_delta0.5`, and the stage-effect scaling run's `band_n320` (δ = 0.05 at n = 320, its stage
  effect read from the committed `band/settings.json`). Per column: s ∈ {1.5, 2} at shares 0.5, 0.8 and 0.2, and
  equal dispersion at share 0.8; the anchor at 500 replicates. Each column's balanced equal-dispersion reading is its
  parent's committed one. 14,000 units.
- **Cell sizes.** Stages hold about n / 3 samples and the first group `round(share × stage)` of each:

  | n | share 0.5 | 0.65 / 0.35 | 0.8 / 0.2 |
  |---|---|---|---|
  | 80 | 14 / 13 | 18 / 9 | 22 / 5 |
  | 200 | 34 / 33 | 44 / 23 | 54 / 13 |
  | 320 | 54 / 53 | 70 / 37 | 86 / 21 |
  | 500 | 84 / 83 | 109 / 58 | 134 / 33 |

  Samples per stage in the larger / smaller group (the first stage; the others differ by at most one). The
  balanced designs give the second group one sample fewer per stage at n = 200.
- **The unit is the univariate baseline's, imported unchanged.** One dataset and one production fit (`pls_cv`, rank
  floor 2, net displacement, 999 permutations). On the same fit and RRPP seed it reads `delta`, `angle` and `shape`,
  the omnibus, the split test (`total`, `span`, `complement`), the scan (`scan`, `fdr`) and PCA at the production
  rank (`pca_matched`).
- **The within-group arm.** After each unit the driver regenerates the dataset and refits production at the recorded
  rank with no cross-validation. It then reads the three statistics under `RRPP(..., strata=<group>)` at the unit's
  RRPP seed. The refit's observed statistics equalled the parent record's on every unit (`same_statistics`).
- **Code and host.** Revision `5aab30e`, pushed before any array unit ran, from a clean worktree on `ing`.
  Python 3.11.16, numpy 2.3.5, scipy 1.16.3, scikit-learn 1.8.0, pandas 2.3.3, BLAS pinned to one thread.
- **Where it ran.** Arrays 958479 (`core`, 200 tasks) and 958544 (`columns`, 140), partition `512x1024` with n11
  excluded: 340 of 340 `COMPLETED 0:0`, 20:05–22:17.
    - **186.83 recorded core-hours** (80.86 / 105.97) against the approved cap of 400; the smoke had projected 233
      before inflation.
    - The median unit took 14.8 s in `core` and 25.6 s in `columns`, most of it the cross-validated fit. The
      within-group arm added a median of 0.6–0.9 s.

## 2. Accounting and the reproduction

- **Units.** 18,000 and 14,000 complete: no failed unit, no missing unit, no duplicate, and no failed planted fixed
  point at any scale.
- **The within-group arm.** 32,000 readings, none failed, none inconsistent.
- **Reproduction: passed.** All 36 rows of the rerun design (available and rejected counts per cell and arm) equal
  the committed `stage_span` rate table, none unresolved. This also checks, at run time, that the generator is
  unchanged while the scale is absent.
- **The decision regenerates.** `decide` on the committed aggregates alone, on the workstation, writes a decision
  byte-identical to the cluster's (sha256 `017f8112…`).

## 3. The decision

- **The readings.** `delta`, `angle` and `shape`, each at p < 0.05, on every control cell of every new design. The
  rerun design enters no family.
- **The rule.** `judge_rates` of the stage-span paper grade, imported unchanged. A family holds when its readings
  above the per-cell bound (0.0695, so 35 or more of 500) stay within the count a run with exact size stays within,
  and none is above the family's Bonferroni ceiling.
- **The verdict** is the first that applies: `reproduction_failed`, `unbalance_excess` (`unbalanced_equal` fails),
  `dispersion_excess` (`balanced_unequal` fails), `unbalanced_dispersion_excess` (`wider_smaller` or `wider_larger`
  fails), `size_holds`. Every family is recorded whatever the verdict.

| family | designs | m | allowed above the bound | ceiling (of 500) | above the bound | above the ceiling | holds |
|---|---|---|---|---|---|---|---|
| `unbalanced_equal` | equal dispersion at shares 0.65 and 0.8 | 24 | 2 | 40 | 2 | **1** | no |
| `balanced_unequal` | a scale at share 0.5 | 42 | 3 | 41 | 22 | 20 | no |
| `wider_smaller` | a scale at shares 0.65 and 0.8 | 60 | 4 | 42 | 59 | 48 | no |
| `wider_larger` | a scale at shares 0.35 and 0.2 | 60 | 4 | 42 | 21 | 19 | no |

**Verdict: `unbalance_excess`.** The reproduction passed and `unbalanced_equal` failed, so the order stops there. The
other three families are recorded and, under this verdict, not interpreted (design, Migration Plan).

## 4. Unbalance at equal dispersion

Rejections of 500 at equal dispersion. The balanced reference is the parent's committed design at the same n and δ.

| design | smaller group per stage | cell | `delta` | `angle` | `shape` | omnibus | balanced: `delta` / `angle` / `shape` (omnibus) |
|---|---|---|---|---|---|---|---|
| n = 200, share 0.65 | 23 | anchor | 25 | 20 | 23 | 23 | 34 / 22 / 23 (24) |
| | | translation | 24 | 32 | 18 | 19 | 27 / 22 / 20 (17) |
| n = 200, share 0.8 | 13 | anchor | 27 | 28 | 30 | 20 | 34 / 22 / 23 (24) |
| | | translation | 34 | 28 | 23 | 18 | 27 / 22 / 20 (17) |
| n = 80, share 0.8 | 5 | anchor | **44** | 30 | 31 | 32 | 30 / 27 / 23 (25) |
| n = 500, share 0.8 | 33 | anchor | 19 | 23 | 23 | 19 | 22 / 21 / 20 (15) |
| n = 200, δ = 0.5, share 0.8 | 13 | anchor | 21 | 28 | 21 | 15 | 27 / 25 / 17 (13) |
| band, n = 320, share 0.8 | 21 | anchor | 24 | *37* | 30 | 26 | 21 / 25 / 20 (12) |

Bold is above the ceiling (40); italic is above the bound only.

- **The family fails on one reading.** `delta` at n = 80, share 0.8 reads 0.088 against a ceiling of 0.080. It sits
  at the smallest cells the run reads: 5 samples per stage in the smaller group, 15 in all. The family's second
  reading above the bound, `angle` in the band at share 0.8 (0.074), is within its ceiling. The family allows two
  readings above the bound, so the count alone would have held.
- **Everything else at equal dispersion holds.** At n = 200 the four readings per statistic run 0.036–0.068, and at
  n = 500 and at δ = 0.5 they run 0.038–0.056.
- **The recorded arms hold in this family** (§7). At n = 80 and share 0.8 the omnibus reads 32, `fdr` 23, `total`
  28, `scan` 30, and the within-group `delta` 35.
- **What the run cannot say.** It reads one design at 5 per stage, at one share and one replicate stream. It does not
  say whether the reading is a small-cell property of `delta`, of unbalance, or of both. That is the diagnosis the
  verdict opens (§10).

## 5. The dispersion families (reported, not interpreted)

Under `unbalance_excess` these readings are recorded, as the design requires, and no conclusion about dispersion is
drawn from them until unbalance at equal dispersion is diagnosed.

### 5.1 The core cross at n = 200, δ = 2.0

Rejections of 500 (anchor; translation):

| share | s | wider group | `delta` | `angle` | `shape` | omnibus (recorded) |
|---|---|---|---|---|---|---|
| 0.5 | 1 (committed) | — | 34; 27 | 22; 22 | 23; 20 | 24; 17 |
| 0.5 | 1.25 | — | **86; 93** | 27; 20 | 26; 25 | 48; 56 |
| 0.5 | 1.5 | — | **219; 211** | 24; 34 | 30; *35* | 144; 159 |
| 0.5 | 2 | — | **332; 314** | 27; *38* | 27; **62** | 259; 260 |
| 0.65 | 1.25 | smaller | **89; 68** | *40*; *37* | *40*; *40* | 65; 53 |
| 0.65 | 1.5 | smaller | **226; 196** | **45; 44** | *37*; **55** | 165; 160 |
| 0.65 | 2 | smaller | **389; 335** | *42*; *39* | **51; 71** | 339; 286 |
| 0.8 | 1.25 | smaller | **91; 86** | *38*; *42* | *41*; *41* | 70; 63 |
| 0.8 | 1.5 | smaller | **249; 249** | **56; 55** | **70; 67** | 217; 208 |
| 0.8 | 2 | smaller | **412; 411** | **94; 88** | **118; 120** | 386; 394 |
| 0.35 | 1.25 | larger | **47**; *39* | 16; 23 | 21; 21 | 29; 30 |
| 0.35 | 1.5 | larger | **137; 110** | 20; 21 | 17; 22 | 73; 65 |
| 0.35 | 2 | larger | **229; 194** | 8; 16 | 15; 21 | 157; 132 |
| 0.2 | 1.25 | larger | 20; 28 | 16; 12 | 10; 12 | 14; 16 |
| 0.2 | 1.5 | larger | **66; 59** | 6; 4 | 5; 6 | 21; 16 |
| 0.2 | 2 | larger | **109; 119** | 2; 9 | 0; 3 | 64; 63 |

Bold is above the family's ceiling (41 in the balanced family, 42 in the others); italic is above the bound only.

### 5.2 The columns

Rejections of 500 on the anchor (`delta` / `angle` / `shape`; omnibus in brackets):

| column | share 0.5, s = 1.5 | share 0.5, s = 2 | share 0.8, s = 1.5 | share 0.8, s = 2 | share 0.2, s = 1.5 | share 0.2, s = 2 |
|---|---|---|---|---|---|---|
| n = 80, δ = 2.0 | 220 / 26 / 18 (159) | 334 / 27 / 15 (269) | 271 / 58 / 56 (200) | 422 / 68 / 87 (387) | 77 / 10 / 7 (36) | 84 / 5 / 3 (48) |
| n = 500, δ = 2.0 | 141 / 22 / 16 (79) | 295 / 29 / 28 (207) | 176 / 49 / 63 (136) | 351 / 78 / 125 (321) | 42 / 6 / 4 (13) | 76 / 3 / 1 (26) |
| n = 200, δ = 0.5 | 124 / 21 / 21 (74) | 287 / 46 / 15 (237) | 164 / 59 / 43 (142) | 315 / 97 / 55 (264) | 18 / 12 / 8 (6) | 105 / 4 / 1 (33) |
| band, n = 320, δ = 0.05 | 444 / 91 / 69 (404) | 498 / 229 / 120 (490) | 445 / 70 / 34 (400) | 498 / 263 / 72 (492) | 185 / 61 / 70 (128) | 350 / 110 / 92 (286) |

- **`delta` is above its ceiling in every column at s = 1.5 with balanced groups** (124–444) and wherever the wider
  group is the smaller (164–498).
- **The band fails all three statistics in every direction.** There the stage fit is weak (AUROC 0.55–0.66, §8).
  Where the wider group is the larger, its `angle` and `shape` read 61–110, while every other column reads them
  conservative.
- **Away from the band**, `angle` and `shape` stay within their ceiling with balanced groups, except `angle` at δ = 0.5
  and s = 2 (46). They fail wherever the wider group is the smaller (43–125), and turn conservative where it is the
  larger.

### 5.3 The map

The smallest scale whose reading is above the family's ceiling, by statistic, column and share. The core reads
s ∈ {1.25, 1.5, 2} and the columns s ∈ {1.5, 2}, so 1.5 at a column means at most 1.5. "None" means no reading up to
s = 2.

| column | share | wider group | `delta` | `angle` | `shape` |
|---|---|---|---|---|---|
| n = 200, δ = 2.0 | 0.5 | — | 1.25 | none | 2 (translation only) |
| | 0.65 | smaller | 1.25 | 1.5 | 1.5 (anchor 2) |
| | 0.8 | smaller | 1.25 | 1.5 | 1.5 |
| | 0.35 | larger | 1.25 (translation 1.5) | none | none |
| | 0.2 | larger | 1.5 | none | none |
| n = 80, δ = 2.0 | 0.5 | — | 1.5 | none | none |
| | 0.8 | smaller | 1.5 | 1.5 | 1.5 |
| | 0.2 | larger | 1.5 | none | none |
| n = 500, δ = 2.0 | 0.5 | — | 1.5 | none | none |
| | 0.8 | smaller | 1.5 | 1.5 | 1.5 |
| | 0.2 | larger | 2 | none | none |
| n = 200, δ = 0.5 | 0.5 | — | 1.5 | 2 | none |
| | 0.8 | smaller | 1.5 | 1.5 | 1.5 |
| | 0.2 | larger | 2 | none | none |
| band, n = 320, δ = 0.05 | 0.5 | — | 1.5 | 1.5 | 1.5 |
| | 0.8 | smaller | 1.5 | 1.5 | 2 |
| | 0.2 | larger | 1.5 | 1.5 | 1.5 |

At share 0.8 the equal-dispersion reading is within the ceiling for every statistic and column, except `delta` at
n = 80 (§4). Every map entry is in the decision's `map`, with its rates by cell.

### 5.4 Conservative readings

27 readings are below 0.05 − 2 SE (15 or fewer of 500). 25 of them are `angle` and `shape` where the wider group is
the larger, at every column but the band. For example, at n = 200, share 0.2 and s = 2, `angle` reads 2 and `shape`
0 of 500 on the anchor. The other two are `shape` at 15 with balanced groups and s = 2, at n = 80 and at δ = 0.5. A
conservative reading fails no family; what it costs in power is not read.

## 6. The within-group draws (recorded)

Rejections of 500 under the unrestricted draws and under draws within each group, on the same fit, observed
statistics and seed (anchor; translation in the core, anchor in the columns):

| design | `delta` unrestricted → within | `angle` unrestricted → within | `shape` unrestricted → within |
|---|---|---|---|
| n = 200, share 0.5, s = 2 | 332; 314 → 319; 299 | 27; 38 → 29; 35 | 27; 62 → 28; 57 |
| n = 200, share 0.8, s = 2 | 412; 411 → 288; 301 | 94; 88 → 22; 17 | 118; 120 → 10; 23 |
| n = 200, share 0.2, s = 2 | 109; 119 → 223; 247 | 2; 9 → 21; 33 | 0; 3 → 30; 34 |
| n = 80, share 0.8, equal | 44 → 35 | 30 → 31 | 31 → 32 |
| n = 500, share 0.5, s = 2 | 295 → 291 | 29 → 28 | 28 → 28 |
| n = 200, δ = 0.5, share 0.8, s = 2 | 315 → 136 | 97 → 19 | 55 → 6 |
| band, share 0.5, s = 2 | 498 → 497 | 229 → 233 | 120 → 113 |
| band, share 0.8, s = 2 | 498 → 462 | 263 → 36 | 72 → 2 |
| band, share 0.2, s = 2 | 350 → 451 | 110 → 219 | 92 → 162 |

Labels per family (readings above the bound, above the ceiling, of m):

| arm | `unbalanced_equal` (m = 8) | `balanced_unequal` (14) | `wider_smaller` (20) | `wider_larger` (20) |
|---|---|---|---|---|
| within-group `delta` | holds (1, 0) | exceeds (14, 14) | exceeds (20, 20) | exceeds (20, 20) |
| within-group `angle` | holds (0, 0) | exceeds (4, 3) | holds (1, 0) | exceeds (5, 3) |
| within-group `shape` | holds (0, 0) | exceeds (3, 3) | holds (0, 0) | exceeds (4, 4) |

- **Where the wider group is the smaller**, the within-group draws bring `angle` and `shape` within their ceiling at
  every column, the band included. They lower `delta` but leave it far above (412 → 288 at n = 200, s = 2).
- **With balanced groups** they change almost nothing. `delta` at s = 2 moves from 332 to 319 at n = 200, and from
  295 to 291 at n = 500.
- **Where the wider group is the larger**, they remove the conservatism of `angle` and `shape` at δ ≥ 0.5, and
  overshoot in places (`angle` 4 → 62 at δ = 0.5, share 0.2, s = 2). They raise `delta` (109 → 223 at n = 200,
  s = 2). In the band they raise all three.
- **What the arm separates.** The part of an excess that the within-group draws remove comes from exchanging
  residuals between groups of unequal spread. Balanced `delta` and every column's remaining `delta` are not of that
  kind. The design's untested guess (D7) is that the stage fit follows the wider group's noise more, so that group's
  in-sample path is longer in scores that no residual permutation refits. The run records this and does not test it.
- **Not an option.** The arm calls the primitive's `strata`; neither the harness nor `motco analyze` exposes group
  strata.

## 7. The per-feature and other recorded arms

Labels per family under the same rule, one reading per design and cell (readings above the bound, above the ceiling):

| arm | rule | `unbalanced_equal` (m = 8) | `balanced_unequal` (14) | `wider_smaller` (20) | `wider_larger` (20) |
|---|---|---|---|---|---|
| omnibus (`production`) | min p ≤ 0.05 / 3 | holds (0, 0) | exceeds (14, 14) | exceeds (20, 20) | exceeds (10, 9) |
| `fdr` | min BH-adjusted p ≤ 0.05 | holds (0, 0) | holds (1, 0) | exceeds (20, 20) | holds (0, 0) |
| `total` | p < 0.05 | holds (0, 0) | holds (1, 0) | exceeds (20, 19) | holds (0, 0) |
| `scan` | p ≤ 0.05 | holds (0, 0) | exceeds (2, 1) | exceeds (20, 20) | holds (0, 0) |
| `span` | p < 0.05 | holds (0, 0) | exceeds (13, 11) | exceeds (20, 20) | exceeds (4, 4) |
| `complement` | p < 0.05 | holds (0, 0) | holds (0, 0) | exceeds (20, 20) | holds (0, 0) |
| `pca_matched` | min p ≤ 0.05 / 3 | holds (1, 0) | exceeds (4, 4) | exceeds (18, 16) | holds (0, 0) |

- **The per-feature tests do not share the balanced reading.** `fdr` reads 14–36 of 500 and `total` 19–36 across the
  balanced designs. `scan`'s one reading above its ceiling is at n = 80 and s = 2 (56). The max-T scan already
  exceeded on the committed balanced n = 80 design (38).
- **They share the reading where the wider group is the smaller**, from the smallest scale read. `fdr` reads 57 and
  134 of 500 at s = 1.25 (shares 0.65 and 0.8 at n = 200), and 496–498 at s = 2 and share 0.8 at every column. At
  n = 200, share 0.8 and s = 2 the median null dataset holds 34 BH discoveries (cell table).
- **Where the wider group is the larger they are conservative.** `fdr` and `scan` read 0 of 500 at share 0.2 and
  s ≥ 1.5 at every column.
- **`span` follows production; `complement` follows the per-feature tests.** `span` (the part in the fitted span)
  exceeds in the balanced family (11 readings above the ceiling). `complement` holds there, and exceeds only where
  the wider group is the smaller, reading 498–500 at δ = 0.5 and 451–488 in the band.
- **PCA at the production rank** exceeds less than production: balanced 4 readings above the ceiling (48–82 of 500),
  none where the wider group is the larger.

## 8. The fit by design

The anchor cell's medians:

| column | balanced, equal (committed) | share 0.5, s = 1.5 / 2 | share 0.8: equal / s = 1.5 / s = 2 | share 0.2: s = 1.5 / 2 |
|---|---|---|---|---|
| n = 200, δ = 2.0: rank | 2 | 2 / 2 (q90 3) | 2 / 2 / 2 (q90 3) | 2 / 2 (q90 3) |
| n = 200, δ = 2.0: capture | 0.973 | 0.968 / 0.962 | 0.972 / 0.971 / 0.968 | 0.965 / 0.959 |
| n = 200, δ = 0.5: rank (q10–q90) | 3 | 3 (3–4) / 5 (4–6) | 3 (3–3) / 3 (3–3) / 3 (3–4) | 4 (3–4) / 6 (5–8) |
| n = 200, δ = 0.5: capture | 0.905 | 0.870 / 0.853 | 0.905 / 0.890 / 0.876 | 0.855 / 0.821 |
| band: rank (q10–q90) | 10 (5–19) | 8 (3–19) / 9 (2–19) | 10 (5–19) / 10 (4–19) / 14 (4–19) | 6 (2–19) / 4 (2–19) |
| band: held-out AUROC | 0.664 | 0.612 / 0.581 | 0.664 / 0.641 / 0.628 | 0.581 / 0.550 |
| band: capture | 0.288 | 0.208 / 0.169 | 0.286 / 0.258 / 0.249 | 0.146 / 0.074 |

- **At δ = 2.0 the fit barely moves.** At n = 80, 200 and 500 the median CV rank stays at 2 (3 at n = 80, share
  0.2, s = 2), and the held-out AUROC stays at 1.000. The pooled-stage capture moves at most 0.026 from the
  committed design (0.934 → 0.908 in that same cell).
- **At δ = 0.5 the rank rises with the scale,** most where the wider group is the larger (6 at s = 2), and the capture
  falls.
- **In the band the scale weakens the stage fit.** The AUROC falls from 0.664 to 0.550 and the capture from 0.288 to
  0.074 where the wider group is the larger.
- The committed values are the parents' cell tables; the run's are in `cell_table.csv` and in the decision's `fit`.

## 9. The run against the probes

The scales, the columns and the within-group arm were chosen after two exploratory probes, at seeds outside the run.
The probes saw rejection rates. The rule did not take any threshold from them: it is the repository's standing size
rule. Against the run, on `delta`. The probe values are from the change's design (Context) and `notes.md`, not from
this run's aggregates.

| δ | share | s | probe 1 (two blocks, CV rank, 100 × 199) | probe 2 (three blocks, forced rank, 200 × 199) | run (three blocks, CV rank, 500 × 999) |
|---|---|---|---|---|---|
| 2.0 | 0.5 | 2 | 0.17 (rank 3) | 0.075 (rank 3) | 0.664 (rank 2) |
| 2.0 | 0.8 | 2 | 0.44 (rank 3) | 0.370 (rank 3) | 0.824 (rank 2) |
| 2.0 | 0.2 | 2 | 0.01 (rank 3) | 0.005 (rank 3) | 0.218 (rank 2) |
| 0.5 | 0.5 | 2 | 0.65 (rank 9) | 0.530 (rank 4) | 0.574 (rank 5) |
| 0.5 | 0.8 | 2 | 0.71 (rank 5) | 0.575 (rank 4) | 0.630 (rank 3) |

- **At δ = 0.5 the run agrees with the probes** within about 0.1.
- **At δ = 2.0 it reads far more.** The within-group draws read 0.638 in the run against 0.015 in probe 2 at share
  0.5.
- **The probes read at rank 3 there and the run at rank 2.** They also differ in blocks, replicates, permutations and
  seeds, and the aggregates hold no rate by rank, so the run does not say which difference matters.
- **The scale itself was checked.** The generator's tests pin the stretch about the cell means to floating-point
  precision, and the gate pins the generator with the scale absent.

## 10. What follows

- **Unbalance at equal dispersion is diagnosed first, as its own change** (design, Migration Plan). The roadmap opens
  it (open item 16). It starts from one reading: `delta` at n = 80 with 5 samples per stage in the smaller group,
  0.088 against a ceiling of 0.080, beside a committed balanced 0.060 at the same n. It needs to say whether that
  reading is a property of small cells, of unbalance or of both, and whether it carries to a case study's table.
  SEA-AD's smallest cell holds 9 donors (F 11 / 10 / 28 against M 9 / 10 / 12). It read size at its own design under
  the `df` correction (roadmap items 7 (c) and 8).
- **The dispersion families are reported, not interpreted.** The documentation states their rates and map with that
  status ([measurement visibility](../measurement-visibility.md#unequal-group-dispersion-and-unbalance)). The
  interpretability statement gains no dispersion clause. It records that the null assumes one dispersion for both
  groups and that the reading waits on the diagnosis. Once the diagnosis reads, the families of this run can be
  re-read under the next verdict in the order, without new units.
- **The recorded arms say how far the per-feature tests share the reading.** They hold where the groups are balanced
  and where the wider group is the larger. They fail where the wider group is the smaller.
- **The roadmap's "Not yet established" entries on unbalance and on noise** now cite this run.
- **What does not change.** `motco analyze` runs what it ran; no default moves; `group_dispersion_scale` stays an
  opt-in generator setting.
- **Still open.** A dispersion that differs by feature, by stage or in correlation structure, and the same reading on
  a real cohort (roadmap open items 2 and 15).

## 11. Limitations

- **Simulated Gaussian noise, one scale.** Every reading holds at 658 simulated features with one scale applied to
  every feature of one group, three stages and no covariates. A difference in correlation structure, which acts in
  the scores as a dispersion scale but not a uniform one (design D3), was not simulated. No real cohort's dispersion
  was read, and neither case study rejects.
- **The verdict rests on one reading.** `unbalanced_equal` fails on a single reading above its ceiling, at the
  smallest cells read. The rule was fixed before any unit ran and is applied as written.
- **The rule followed the probes.** The probes saw rejection rates at a scale of 2. They chose the scales, the
  columns and the recorded arm, not the rule's thresholds (§9).
- **The columns read s ≥ 1.5.** Below that only n = 200 at δ = 2.0 is read, so a column's first failing scale of 1.5
  is an upper bound.
- **Size only.** No unit carries a trajectory difference, so neither the power a conservative reading costs nor power
  under unequal dispersion is read.
- **The within-group arm is a reading on a refit at the recorded rank.** Its observed statistics equal production's
  on every unit, but its draws assume equal dispersion across the stages of a group, the mirror of the within-stage
  null's assumption.
