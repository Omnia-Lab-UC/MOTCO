# Unbalance at equal dispersion — report (2026-11-25)

**Change:** `unbalance-equal-dispersion-diagnosis` (roadmap open item 16; no package change, MOTCO 0.26.0).
**Reads:** `results/unbalance-diagnosis-2026-11-25/` (`unbalance_diagnosis_decision.json`, `settings.json`,
`rate_table.csv`, `cell_table.csv`, `within_group_table.csv`, `corrected_table.csv`, `accounting.json`,
`reproduction.json`, `PROVENANCE.txt`, `sacct.txt`). **Parent:** the
[group-dispersion run](group-dispersion-2026-11-24.md)'s `columns` profile and decision (its settings, rate, cell and
within-group tables and `group_dispersion_decision.json`, bound by sha256). **Concept page:**
[measurement visibility](../measurement-visibility.md#unequal-group-dispersion-and-unbalance).

The group-dispersion run read `unbalance_excess` on one reading of 24. At equal dispersion, `delta` at n = 80 with a
first-group share of 0.8 rejected 44 times of 500 against a ceiling of 40. That design leaves 5 samples per stage in
the smaller group, so it is both the smallest and the most unbalanced design read, and the run could not say which
of the two mattered. By the predeclared order its three dispersion families were reported and not interpreted. This
run asks whether the reading replicates and whether it follows small cells, unbalance or both. Every unit is an
equal-dispersion null unit.

- **The designs.** The committed `n80_r0.8` design, rerun at its name and seeds as a gate, and seven new designs at
  1,000 replicates that cross the smaller group's cell size (5, 7, 9, 13 per stage) with balance and total sample
  size.
- **The decision**, fixed before any unit ran, judges `delta`, `angle` and `shape` in three families and takes the
  first that fails in a fixed order: balanced small cells, then a small group in an unbalanced design, then n = 80 at
  9 or 13 per stage. It also fixes how the group-dispersion families are re-read under each verdict.

**In one paragraph.**

- **The verdict is `small_cells`, on one reading.** `balanced_small` fails: `angle` at n = 42 with 7 samples per
  stage in both groups rejects 76 times of 1,000, against a ceiling of 65. `delta` (54) and `shape` (47) hold there.
- **Unbalance is not needed, and the committed reading did not recur.**
    - Every reading with a small group in an unbalanced design holds (`unbalanced_small`, 12 of 12), at 5 samples per
      stage as at 7.
    - So does every reading at n = 80 with 9 or 13 per stage (`reference`, 6 of 6).
    - A fresh replicate of the failing design reads `delta` at 45 of 1,000, and so does a 5-per-stage group inside
      n = 200 (44). Pooled with the rerun, `delta` at `n80_r0.8` reads 89 of 1,500 (0.059).
- **The re-read gives `dispersion_excess`.** The five group-dispersion designs with a 5-per-stage group are set aside.
  Unbalance at equal dispersion then holds (21 readings, none above the ceiling). The balanced dispersion family fails
  as it did (20 of 42 above the ceiling).
- **What that permits.** Within readings at 13 samples per stage or more in every cell, a difference in noise
  dispersion alone takes the trajectory statistics above their size, with balanced groups too. A rejection there is
  evidence of a difference in trajectory or in dispersion.
- **The recorded arms.** The omnibus holds in every family, and so do the `df` and `leverage` corrections and the
  per-feature tests. At n = 42 the `angle` excess is shared by PCA at the same rank (77), the within-group draws (70)
  and the split test's `span` (68). The `df` correction reads 62.
- **The gate passed.** The rerun reproduces all 22 committed rows exactly. 7,500 of 7,500 units completed, with no
  inconsistent refit, in 35.4 recorded core-hours against the approved cap of 80.
- **What it holds for.** Simulated Gaussian noise at equal dispersion, 658 features on three blocks, δ = 2.0, three
  stages and no covariates. No balanced design below 7 samples per stage was read. The re-read rule was written after
  the group-dispersion families' rates were known. No verdict changes what `motco analyze` runs.

## 1. Configuration and provenance

- **One run, no profile.** Methylation, expression and proteomics (658 features), three stages, δ = 2.0, the anchor
  cell, equal dispersion, base seed 2043, 999 permutations. The generator parameters are read from the parent's
  committed `columns/settings.json`; every new design differs from `n80_r0.8` only in `n_samples` and `group_ratio`.
- **The designs.** Stages hold about n / 3 samples and the first group `round(share × stage)` of each:

  | design | role | n | share | larger / smaller per stage | replicates | family |
  |---|---|---|---|---|---|---|
  | `n80_r0.8` | the gate (committed name and seeds) | 80 | 0.8 | 22 / 5 | 500 | none |
  | `n80_r0.8_b` | a fresh replicate of the failing design | 80 | 0.8 | 22 / 5 | 1,000 | `unbalanced_small` |
  | `n200_r0.925` | a 5-per-stage group inside a large sample | 200 | 0.925 | 62 / 5 | 1,000 | `unbalanced_small` |
  | `n42_r0.5` | small cells with balanced groups | 42 | 0.5 | 7 / 7 | 1,000 | `balanced_small` |
  | `n80_r0.74` | 7 per stage, unbalanced, at n = 80 | 80 | 0.74 | 20 / 7 | 1,000 | `unbalanced_small` |
  | `n200_r0.895` | 7 per stage inside a large sample | 200 | 0.895 | 60 / 7 | 1,000 | `unbalanced_small` |
  | `n80_r0.65` | SEA-AD's smallest cell size | 80 | 0.65 | 18 / 9 | 1,000 | `reference` |
  | `n80_r0.5_b` | a fresh balanced reference at n = 80 | 80 | 0.5 | 14 / 13 | 1,000 | `reference` |

  7,500 units. A balanced design at 5 per stage (n = 30) cannot run under the production cross-validation, so 7 per
  stage is the smallest balanced cell read.
- **The unit is the group-dispersion run's, imported unchanged.** One dataset and one production fit (`pls_cv`, rank
  floor 2, net displacement). On the same fit and RRPP seed it reads `delta`, `angle` and `shape`, the omnibus, the
  split test (`total`, `span`, `complement`), the scan (`scan`, `fdr`), PCA at the production rank (`pca_matched`),
  and the within-group draws (`RRPP(..., strata=<group>)`) on a refit at the recorded rank.
- **The corrected arms.** After the parent's unit the driver regenerates the dataset, refits production at the
  recorded rank with no cross-validation, and reads the three statistics under the unrestricted draws with
  `residual_correction="df"` and with `"leverage"`, at the unit's RRPP seed. The refit's observed statistics equalled
  the parent record's on every unit (`same_statistics`).
- **Code and host.** Revision `ade26ad`, pushed before any array unit ran, from a clean worktree on `ing`.
  Python 3.11.16, numpy 2.3.5, scipy 1.16.3, scikit-learn 1.8.0, pandas 2.3.3, BLAS pinned to one thread. The new and
  the group-dispersion test files passed on a compute node first (107 passed).
- **Where it ran.** Array 958918 (150 tasks, `%48`), partition `512x1024` with n11 excluded: 150 of 150
  `COMPLETED 0:0`, 14:17–15:11, tasks 11.2–27.0 min.
    - **35.43 recorded core-hours** against the approved cap of 80; the smoke had projected 48 on a loaded
      workstation (the change's `notes.md`).
    - The median unit took 12.8 s (q90 27.7 s); the median cross-validated fit took 10.1 s. The within-group arm
      added a median of 0.5 s and the corrected arms 0.8 s.

## 2. Accounting and the reproduction

- **Units.** 7,500 of 7,500 complete: no failed unit, no missing unit, no duplicate. n = 42 failed no unit.
- **The refits.** The within-group arm and the corrected arms each completed 7,500 refits, none failed and none
  inconsistent; neither correction failed on any unit.
- **Reproduction: passed.** All 22 compared rows of `n80_r0.8` (available and rejected counts per arm, and its
  within-group row) equal the committed `columns` tables, none unresolved. The `concat` rows are recorded by neither
  run.
- **The decision regenerates.** `decide` on the committed aggregates alone, on the workstation, writes a decision
  byte-identical to the cluster's (sha256 `96b7bede…`).

## 3. The decision

- **The readings.** `delta`, `angle` and `shape`, each at p < 0.05, on the anchor of every new design. The rerun
  design enters no family.
- **The rule.** `judge_rates` of the stage-span paper grade, imported unchanged. A family holds when its readings
  above the per-cell bound (0.0638 at 1,000 replicates, so 64 or more) stay within the count a run with exact size
  stays within, and none is above the family's Bonferroni ceiling. A reading with fewer than 900 available units would
  be unavailable and fail its family; none was.
- **The verdict** is the first that applies: `reproduction_failed`, `small_cells` (`balanced_small` fails),
  `small_unbalanced_cells` (`unbalanced_small` fails), `moderate_excess` (`reference` fails), `not_replicated`.
  Every family is recorded whatever the verdict.

| family | designs | m | allowed above the bound | ceiling (of 1,000) | above the bound | above the ceiling | holds |
|---|---|---|---|---|---|---|---|
| `balanced_small` | `n42_r0.5` (7 / 7) | 3 | 1 | 65 | 1 | **1** | no |
| `unbalanced_small` | `n80_r0.8_b`, `n200_r0.925`, `n80_r0.74`, `n200_r0.895` | 12 | 1 | 69 | 1 | 0 | yes |
| `reference` | `n80_r0.65`, `n80_r0.5_b` | 6 | 1 | 67 | 0 | 0 | yes |

**Verdict: `small_cells`.** The reproduction passed and `balanced_small` failed, so the order stops there. Small cells
lose size with balanced groups, so unbalance is not needed to produce an excess.

- **The failing reading is `angle`, not `delta`.** The run was opened by a `delta` reading; the reading that decides
  it is `angle` at n = 42, 0.076. The design judged all three statistics in advance for this reason (D5).
- **`balanced_small` is the easiest family to fail.** It has three readings, so its ceiling is the Bonferroni one for
  m = 3; one reading above it fails the family. The count is 76 against 65.

## 4. Counts by cell size and balance

Rejections at equal dispersion (of 500 for the rerun, of 1,000 for the new designs):

| design | larger / smaller per stage | `delta` | `angle` | `shape` | omnibus |
|---|---|---|---|---|---|
| `n80_r0.8` (rerun, of 500) | 22 / 5 | 44 | 30 | 31 | 32 |
| `n80_r0.8_b` | 22 / 5 | 45 | 53 | 51 | 37 |
| `n200_r0.925` | 62 / 5 | 44 | 46 | 51 | 40 |
| `n42_r0.5` | 7 / 7 | 54 | **76** | 47 | 53 |
| `n80_r0.74` | 20 / 7 | 58 | 46 | *65* | 43 |
| `n200_r0.895` | 60 / 7 | 48 | 47 | 51 | 33 |
| `n80_r0.65` | 18 / 9 | 52 | 62 | 59 | 47 |
| `n80_r0.5_b` | 14 / 13 | 36 | 57 | 47 | 37 |

Bold is above the family's ceiling; italic is above the bound only. The omnibus (min p ≤ 0.05 / 3) is recorded, not
judged.

- **A small group in an unbalanced design holds.** At 5 per stage `delta` / `angle` / `shape` read 45 / 53 / 51 at
  n = 80 and 44 / 46 / 51 at n = 200. At 7 per stage they read 58 / 46 / 65 at n = 80 and 48 / 47 / 51 at n = 200.
  The one reading above the bound, `shape` at `n80_r0.74` (65), is within its ceiling of 69.
- **n = 80 holds at 9 and 13 per stage.** The largest reading is `angle` at `n80_r0.65`, 62.
- **The excess is at n = 42.** There both groups hold 7 per stage, the sample is the smallest read, and only `angle`
  exceeds. The same smaller-group cell size inside n = 80 or n = 200 holds.
- **The omnibus holds everywhere**, from 33 to 53 of 1,000, n = 42 included.

## 5. The committed reading

`delta` at `n80_r0.8`, the reading that opened the run:

| source | replicates | rejected | rate |
|---|---|---|---|
| the group-dispersion run (committed), reproduced by the rerun | 500 | 44 | 0.088 |
| fresh replicate `n80_r0.8_b` | 1,000 | 45 | 0.045 |
| pooled | 1,500 | 89 | 0.059 (SE 0.006) |

- **The rerun reproduces the committed count exactly**, so the generator, the fit and the environment are unchanged.
- **The fresh replicate reads at the level.** The same design on new datasets gives 45 of 1,000, and a 5-per-stage
  group inside n = 200 gives 44.
- **The pooled count is reported beside the verdict, not judged.** `not_replicated` would have applied had every
  family held; the order stopped at `small_cells` first.

## 6. The recorded arms

Rejections of 1,000 (of 500 for the rerun), on the same fit, observed statistics and seed. `delta` / `angle` /
`shape`:

| design | production | within group | `df` | `leverage` | PCA at the production rank |
|---|---|---|---|---|---|
| `n80_r0.8` (of 500) | 44 / 30 / 31 | 35 / 31 / 32 | 38 / 26 / 31 | 38 / 26 / 31 | 28 / 34 / 28 |
| `n80_r0.8_b` | 45 / 53 / 51 | 47 / 47 / 50 | 37 / 45 / 45 | 37 / 45 / 44 | 33 / 57 / 51 |
| `n200_r0.925` | 44 / 46 / 51 | 42 / 42 / 43 | 40 / 44 / 48 | 41 / 44 / 47 | 59 / 48 / 53 |
| `n42_r0.5` | 54 / 76 / 47 | 44 / 70 / 49 | 32 / 62 / 34 | 32 / 62 / 34 | 18 / 77 / 45 |
| `n80_r0.74` | 58 / 46 / 65 | 52 / 46 / 55 | 49 / 40 / 54 | 49 / 40 / 52 | 40 / 56 / 55 |
| `n200_r0.895` | 48 / 47 / 51 | 42 / 42 / 50 | 45 / 43 / 47 | 44 / 44 / 49 | 50 / 48 / 54 |
| `n80_r0.65` | 52 / 62 / 59 | 53 / 68 / 50 | 51 / 54 / 51 | 51 / 54 / 49 | 30 / 68 / 42 |
| `n80_r0.5_b` | 36 / 57 / 47 | 38 / 54 / 47 | 33 / 54 / 38 | 33 / 54 / 38 | 12 / 54 / 39 |

The per-feature and split-test arms, rejections of 1,000:

| design | `fdr` | `scan` | `total` | `span` | `complement` |
|---|---|---|---|---|---|
| `n80_r0.8` (of 500) | 23 | 30 | 28 | 32 | 25 |
| `n80_r0.8_b` | 43 | 56 | 53 | 48 | 47 |
| `n200_r0.925` | 48 | 54 | 52 | 51 | 53 |
| `n42_r0.5` | 44 | 54 | 59 | 68 | 46 |
| `n80_r0.74` | 41 | 53 | 50 | 58 | 42 |
| `n200_r0.895` | 35 | 51 | 58 | 39 | 58 |
| `n80_r0.65` | 48 | 60 | 64 | 59 | 49 |
| `n80_r0.5_b` | 39 | 49 | 45 | 47 | 40 |

Labels per family under the same rule, one reading per new design (readings above the bound, above the ceiling):

| arm | `balanced_small` (m = 1) | `unbalanced_small` (4) | `reference` (2) |
|---|---|---|---|
| omnibus (`production`), `pca_matched`, `fdr`, `scan`, `complement` | holds (0, 0) | holds (0, 0) | holds (0, 0) |
| `total` | holds (0, 0) | holds (0, 0) | holds (1, 0) |
| `span` | exceeds (1, 1) | holds (0, 0) | holds (0, 0) |
| within-group `delta`, `shape` | holds (0, 0) | holds (0, 0) | holds (0, 0) |
| within-group `angle` | exceeds (1, 1) | holds (0, 0) | exceeds (1, 1) |
| `df` and `leverage`, each statistic | holds (0, 0) | holds (0, 0) | holds (0, 0) |

- **The `angle` excess at n = 42 is not production's alone.** PCA at the production rank reads 77, the within-group
  draws 70 and `span` 68, on the same datasets. The within-group draws leave it, so exchanging residuals between the
  groups does not cause it.
- **The `df` correction brings it below the bound** (62), and lowers `delta` and `shape` at n = 42 to 32 and 34.
  Without covariates the reduced model holds four columns, which the design (D8) put at a scale of about 1.05 at
  n = 42. `leverage` reads as `df` to within two rejections in every design.
- **The per-feature tests hold everywhere.** `fdr` reads 35–48 of 1,000 across the new designs and `scan` 49–60.
- **The within-group `angle` also exceeds in `reference`,** at `n80_r0.65` (68, where production reads 62). These
  labels are recorded and decide nothing.

## 7. The fit by design

| design | rank (q10–q90) | held-out AUROC | pooled-stage capture |
|---|---|---|---|
| `n80_r0.8`, `n80_r0.8_b` | 2 (2–2) | 1.000 | 0.933, 0.934 |
| `n200_r0.925`, `n200_r0.895` | 2 (2–2) | 1.000 | 0.972 |
| `n42_r0.5` | 2 (2–2) | 1.000 | 0.881 |
| `n80_r0.74`, `n80_r0.65`, `n80_r0.5_b` | 2 (2–2) | 1.000 | 0.934, 0.933, 0.934 |

The CV rank is 2 and the held-out AUROC 1.000 in every design. The fitted span keeps less of the pooled stage program
at n = 42 (0.881) than at n = 80 (0.933–0.934) or n = 200 (0.972). The values are in `cell_table.csv` and in the
decision's `fit`.

## 8. The re-read of the group-dispersion families

The rule (design D7) was fixed before this run's units and after the group-dispersion families' rates were known.
Under `small_cells` it sets aside the five group-dispersion designs with a 5-per-stage group: `n80_r0.8`,
`n80_r0.8_s1.5`, `n80_r0.8_s2`, `n80_r0.2_s1.5` and `n80_r0.2_s2`. It judges the four families again by
`judge_rates` on the remaining committed counts and takes the group-dispersion order. No other group-dispersion design
has a cell below 13 per stage.

| family | committed: m, above the bound, above the ceiling | re-read: m, above the bound, above the ceiling | re-read ceiling (of 500) | holds |
|---|---|---|---|---|
| `unbalanced_equal` | 24, 2, 1 | 21, 1, 0 | 40 | yes |
| `balanced_unequal` | 42, 22, 20 | 42, 22, 20 | 41 | no |
| `wider_smaller` | 60, 59, 48 | 54, 53, 44 | 41 | no |
| `wider_larger` | 60, 21, 19 | 54, 19, 18 | 41 | no |

**Re-read verdict: `dispersion_excess`.** With the 5-per-stage designs set aside, unbalance at equal dispersion holds:
its one reading above the bound, `angle` in the band at share 0.8 (37 of 500), is within the ceiling. The balanced
dispersion family is unchanged by the set-aside and fails first.

- **Its scope.** Readings at 13 samples per stage or more in every cell: n = 200, 320 and 500 at every share read,
  and n = 80 with balanced groups.
- **What it was foreseen to give.** The design noted that `balanced_unequal` fails at 33 or more samples per stage,
  so under `small_cells` or `small_unbalanced_cells` the re-read was bound to give `dispersion_excess`. What was not
  known in advance was this run's verdict.
- **What it permits** (the group-dispersion design's Migration Plan). Size needs equal group dispersion even with
  equal groups. A rejection is evidence of a difference in trajectory or in dispersion, for the statistics and scales
  the group-dispersion map gives ([its report](group-dispersion-2026-11-24.md) §5.3, outside the set-aside rows), in
  balanced designs too:
    - `delta` at every column and share read, from s = 1.25, 1.5 or 2;
    - `angle` and `shape` wherever the wider group is the smaller, and with balanced groups in the weak-stage band
      and at s = 2 in one other column each.
- **What it does not change.** The group-dispersion files are not written to; their verdict stays `unbalance_excess`
  as committed, and this decision carries the re-read.

## 9. What follows

- **The documentation states the cell sizes at which size held at equal dispersion**
  ([measurement visibility](../measurement-visibility.md#unequal-group-dispersion-and-unbalance)):
    - every design at n ≥ 80, down to a smaller group of 5 per stage, balanced or not;
    - the one failure, `angle` at n = 42 with 7 per stage in both groups;
    - that no balanced design below 7 per stage was read.
- **The dispersion readings are interpreted, within the re-read's scope.** The interpretability statement of each
  statistic gains the dispersion clause. A case study that rejects reports the two groups' residual dispersion in the
  scores beside the rejection, by the placement steps the measurement-visibility page gives.
- **The roadmap.**
    - Open item 16 moves to Done.
    - A diagnosis of the `delta` excess the within-group draws leave is opened, and after it a change for a null or a
      fit that holds size under unequal dispersion, starting from the recorded within-group arm.
    - Open item 3's bias correction needs a per-group trace.
    - The "Not yet established" entries on unbalance and on noise cite this run.
- **What does not change.** `motco analyze` runs what it ran; no default moves; `df` stays opt-in.
- **Still open.** Whether the n = 42 `angle` excess follows the cell size or the total sample size: at n = 42 both
  groups are small and so is the sample, and no design separates the two. A dispersion that differs by feature, by
  stage or in correlation structure, and the same reading on a real cohort (roadmap open items 2 and 15).

## 10. Limitations

- **Simulated Gaussian noise at one stage effect.** Every reading holds at 658 simulated features on three blocks,
  equal group dispersion, δ = 2.0 (held-out AUROC 1.000), three stages and no covariates. A weaker stage effect, four
  stages and covariates were not read. A cohort's own table is read at its own design: SEA-AD's held size at its own
  design under the `df` correction (roadmap items 7 (c) and 8). What carries over from this run is only the cell size
  at which a simulated design at δ = 2.0 loses size.
- **No balanced design below 7 per stage.** A balanced design at 5 per stage (n = 30) cannot run under the production
  cross-validation. Whether balanced groups would hold at 5 per stage is unread.
- **The verdict rests on one reading.** `balanced_small` fails on one of its three readings, by 76 against a ceiling
  of 65. The rule was fixed before any unit ran and is applied as written.
- **The re-read rule was written after the group-dispersion rates were known**, and its outcome under this verdict was
  foreseeable (§8).
- **The design set was chosen after the committed rates were seen.** The rule is the standing size rule, with no
  threshold taken from them, and `not_replicated` and `moderate_excess` were both outcomes it could return.
- **Size only.** Every unit is a null unit; power at these designs is not read.
- **The recorded arms are readings on refits at the recorded rank.** Their observed statistics equal production's on
  every unit, and their labels decide nothing.
