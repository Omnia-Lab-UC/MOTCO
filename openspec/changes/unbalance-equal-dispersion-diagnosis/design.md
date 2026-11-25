# Design

## Context

See proposal.md (Why) for the motivation. The facts below shape the approach.

- **The committed readings at equal dispersion** (`results/group-dispersion-2026-11-24/`, rejections of 500 on the anchor, δ = 2.0 unless noted):

  | design | cells (larger / smaller per stage) | `delta` | `angle` | `shape` |
  |---|---|---|---|---|
  | n = 80, share 0.8 | 22 / 5 | 44 | 30 | 31 |
  | n = 80, balanced (univariate baseline) | 14 / 13 | 30 | 27 | 23 |
  | n = 200, share 0.8 | 54 / 13 | 27 | 28 | 30 |
  | n = 200, share 0.65 | 44 / 23 | 25 | 20 | 23 |
  | n = 500, share 0.8 | 134 / 33 | 19 | 23 | 23 |

  On `n80_r0.8` the omnibus reads 32, `fdr` 23, `total` 28 and the within-group `delta` 35.
- **The sizing.** Each stage holds about n / 3 samples and the first group `round(share × stage)` of them. The cell sizes in the proposal's table were computed with the generator's sizing functions.
- **The production cross-validation bounds the sample size.** It tries ranks up to 19 inside inner training folds of about half the sample. A probe on 2026-11-25 (scratch, one fit per design, no permutation, no rate read) found that the fit raises at n = 30 and n = 36 ("`n_components` upper bound is 14" and 18), and runs at n = 42 (7 / 7), at n = 80 with share 0.74 (20 / 7) and at n = 200 with shares 0.895 (60 / 7) and 0.925 (62 / 5).
- **The seeds.** A unit's dataset seed is a hash of the base seed, the design name and the replicate. A committed name gives the committed datasets, and a new name gives new ones.
- **The parent's unit.** `group_dispersion_check.evaluate_unit` runs the univariate baseline's unit and then the within-group arm on a refit at the recorded rank. Its `strata_reading` does not pass a residual correction to `RRPP`; the primitive accepts one.
- **The rule.** `judge_rates` takes rows of available and rejected counts and sets the bound, the allowed count and the Bonferroni ceiling from the family's size and the replicates.
- **What was seen before this rule was written.** Every rate of the group-dispersion run, including all three dispersion families.

## Goals / Non-Goals

**Goals:**
- Say, under a rule fixed before any unit, whether the committed `delta` excess replicates, and whether it follows small cells, unbalance or both.
- Say in advance how the group-dispersion run's families are read under each outcome.
- Record whether the existing residual corrections or the within-group draws change the reading.

**Non-Goals:**
- Naming a mechanism. The recorded arms bracket it and decide nothing.
- A balanced reading below 7 per stage.
- Any reading at unequal dispersion, at another stage effect or on a cohort's table.

## Decisions

**D1. Every design is at equal dispersion, δ = 2.0, on the anchor cell.**
- The failing reading is the anchor of a δ = 2.0 design. The translation control and the other columns held, so they are not rerun.
- **Alternative: add δ = 0.5 and the band.** Both held at share 0.8 with 13 and 21 per stage, so they would need their own small-cell designs. That triples the run for a question the verdict does not ask.

**D2. The smallest balanced cell is 7 per stage, and the unbalanced designs are matched to it.**
- A balanced design at 5 per stage is n = 30, which the production fit cannot run. `n42_r0.5` (7 / 7) is the smallest that does.
- `n80_r0.74` (20 / 7) and `n200_r0.895` (60 / 7) hold the same smallest cell, so 7 per stage is read balanced, unbalanced at n = 80 and unbalanced inside a large sample.
- 5 per stage is read unbalanced only: `n80_r0.8_b` (22 / 5) and `n200_r0.925` (62 / 5).
- **Alternative: lower the candidate ranks to fit n = 30.** That changes the production fit, and a difference could then be the cross-validation's.
- **What it costs.** If 5 per stage fails and 7 per stage holds everywhere, the run cannot say whether balanced groups would have protected 5 per stage. The verdict is then `small_unbalanced_cells`, and the report states the limit.

**D3. One fresh replicate, a rerun gate and a fresh balanced reference.**
- **The gate:** `n80_r0.8` at its committed name and seeds, 500 units. Every arm's counts must equal the committed `columns/rate_table.csv` and `columns/within_group_table.csv`.
- **The replicate:** `n80_r0.8_b`, the same generator parameters under a new name, so new datasets.
- **The reference:** `n80_r0.5_b`. The committed balanced n = 80 design reads `delta` at 30 of 500. A fresh one says whether n = 80 itself sits above the level.
- `n80_r0.65` (18 / 9) holds SEA-AD's smallest cell size at n = 80.

**D4. 1,000 replicates per new design.**
- The committed reading is 0.088 and the balanced reference 0.060. At 500 replicates the standard error near these rates is 0.011–0.013, as large as their difference from the bound. At 1,000 it is 0.008–0.009.
- Units are cheap here: n ≤ 200.
- **Alternative: 500, the standing grade.** It would leave a replicate at 0.07 undecided.

**D5. Three families, judged by `judge_rates` on `delta`, `angle` and `shape`.**

| family | designs | m | question |
|---|---|---|---|
| `balanced_small` | `n42_r0.5` | 3 | do small cells lose size without unbalance? |
| `unbalanced_small` | `n80_r0.8_b`, `n200_r0.925`, `n80_r0.74`, `n200_r0.895` | 12 | does a small group in an unbalanced design lose size? |
| `reference` | `n80_r0.65`, `n80_r0.5_b` | 6 | does n = 80 lose size at 9 or 13 per stage? |

- All three statistics are judged, as in the parent, although the committed excess is `delta`'s. The report gives the statistic and design of every failing reading.
- **Why families and not designs.** Seven designs judged one by one would give a run with exact size about seven chances of some failure.

**D6. The verdict follows a fixed order.**
1. `reproduction_failed`: the gate.
2. `small_cells`: `balanced_small` fails. Small cells lose size with balanced groups, so unbalance is not needed.
3. `small_unbalanced_cells`: `unbalanced_small` fails. A group of 5 or 7 per stage in an unbalanced design loses size, and balanced groups at 7 per stage hold.
4. `moderate_excess`: `reference` fails. The excess is not confined to small cells.
5. `not_replicated`: every family holds. The committed reading is not reproduced on fresh datasets.

Every family is recorded whatever the verdict. Beside it the decision reports, for each new design, the three counts, and for `delta` the pooled count of `n80_r0.8` and `n80_r0.8_b` (1,500 units).

**D7. The re-read of the group-dispersion families is fixed here, before any unit.**
- **Under `small_cells` or `small_unbalanced_cells`.** The excess at equal dispersion is scoped to designs with a small group. The group-dispersion designs with a 5-per-stage group are set aside: `n80_r0.8`, `n80_r0.8_s1.5`, `n80_r0.8_s2`, `n80_r0.2_s1.5` and `n80_r0.2_s2`. That removes 3 readings from `unbalanced_equal` and 6 each from `wider_smaller` and `wider_larger`. The four families are judged again by `judge_rates` on the remaining committed counts, and the group-dispersion order gives the re-read verdict. No other group-dispersion design has a cell below 13 per stage.
- **Under `not_replicated`.** The three `n80_r0.8` readings of `unbalanced_equal` are replaced by `n80_r0.8_b`'s, the families are judged in full, and the group-dispersion order gives the re-read verdict.
- **Under `moderate_excess` or `reproduction_failed`.** Nothing is re-read. `unbalance_excess` stands, and the next step is a diagnosis at n = 80.
- **Where it is written.** `decide` reads the committed `group_dispersion_decision.json`, bound by sha256, and writes the re-read families, the set-aside designs and the re-read verdict into this run's decision. The group-dispersion run's own files are not changed.
- **What the re-read verdict permits.** The consequence the group-dispersion design attached to that verdict (its Migration Plan), with the scope stated: readings at 13 samples per stage or more in every cell.
- **Alternative: leave the re-read to a later change.** It would be written after this run's outcome was known as well.

**D8. The recorded arms.**
- **From the parent, unchanged:** the omnibus, `scan`, `fdr`, `total`, `span`, `complement`, `pca_matched` and the within-group statistics.
- **New: `df` and `leverage`.** After the parent's unit, on the within-group arm's refit at the recorded rank, the driver reads `delta`, `angle` and `shape` under the unrestricted draws with `residual_correction="df"` and with `"leverage"`, at the unit's RRPP seed. The observed statistics must equal the production record's (`same_statistics`), and a unit that differs leaves these arms only.
- **Why record them.** The case-study configuration's `delta` excess at n = 80 was diagnosed as residual shrinkage (roadmap item 7 (c)), and `df` was adopted there (item 8). Without covariates the reduced model has four columns, so `df` rescales by about 1.03 at n = 80 and 1.05 at n = 42, and it may change little.
- Each recorded arm gets a label per family under the same rule. None changes the verdict.

**D9. Execution follows the parent.**
- One script with `run`, `merge` and `decide`; `settings.json` bound by hash with the parent files' sha256; round-robin shards; resume by complete unit key; JSONL records with seeds, statistics and counts only.
- An sbatch array with a dry-run variable, BLAS pinned to one thread, run on `ing` from a clean worktree at a pushed revision.
- Committed: aggregates, the decision, provenance and the report under `results/unbalance-diagnosis-<date>/`. Records stay on the cluster.

## Risks / Trade-offs

- **[The re-read rule was written after the dispersion families were seen]** → Its outcome under three of the five verdicts is foreseeable: `balanced_unequal` fails at 33 or more samples per stage, so the re-read gives `dispersion_excess`. What is not known is this run's verdict. The report says so beside the re-read.
- **[The design set was chosen after seeing the committed rates]** → The rule is the standing size rule, with no threshold taken from them. The designs test the reading those rates suggest, and `not_replicated` and `moderate_excess` are both outcomes the rule can return.
- **[Balanced 5 per stage is unread]** → Stated in D2 and in the report's limits.
- **[`balanced_small` has three readings]** → Its ceiling is the Bonferroni one for m = 3 at 1,000 replicates, so it is the easiest family to fail by a single reading. The report gives the counts beside the label.
- **[n = 42 is near the fit's limit]** → The smoke runs every design, and failed units are counted per design. A design whose units fail at a rate that leaves fewer than 900 available is reported as unavailable, which fails its family by the rule.
- **[The gate fails]** → The reading stops. A mismatch is diagnosed, never waived.
- **[The verdict is read as a finding about SEA-AD]** → SEA-AD's table has six covariates, two blocks and a weak stage effect, and it was read at its own design. The report says what carries over: only the cell size at which a simulated design at δ = 2.0 loses size.

## Migration Plan

No package change and no default moves. What each verdict changes in the documentation:
- **`small_cells` or `small_unbalanced_cells`:** the documentation states the smallest cell size from which `delta`, `angle` and `shape` held and the designs that failed, and applies the re-read verdict's consequence with its scope. Roadmap open item 16 moves to Done.
- **`not_replicated`:** the documentation states that the committed reading did not reproduce, gives the pooled count, and applies the re-read verdict's consequence.
- **`moderate_excess`:** `unbalance_excess` stands; the roadmap keeps a diagnosis at n = 80 open and the dispersion families stay uninterpreted.
- **Rollback:** removing the driver, its tests, the aggregates and the report restores the previous state.

## Open Questions

- The core-hour cap for the array is set from the smoke timing, at the approval gate.
