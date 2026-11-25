# Design

## Context

See proposal.md (Why) for the motivation. The facts below shape the approach.

- **How the generator samples.** `generate_omics` draws each group × stage cell of each omic from a multivariate normal with the cell's exact mean and the reference covariance (`rng.multivariate_normal(mean, cov, method="svd")`), methylation in M-values before the inverse logit. Each group is one call, group A first. The semi-synthetic layer has no setting that touches the covariance.
- **What the null assumes.** RRPP permutes the reduced model's residuals over all samples. A draw therefore gives both groups the pooled dispersion. The opt-in `strata` of the primitive (0.25.0) restricts a draw to permute within the rows of each label; the harness and `motco analyze` expose only `stage` as the label and reject any other value.
- **What the parent committed.** The univariate baseline's `stage_span` profile ran four designs on methylation, expression and proteomics (658 features), three stages, equal groups, base seed 2043: `n80`, `n200` and `n500` at δ = 2.0 and `n200_delta0.5`. Each has the anchor and a translation at e = 1.0 under the planted scale at 500 × 999. Its `rate_table.csv` holds the available and rejected counts of every arm. Its unit builds the generator parameters from a per-design dictionary, so a design may carry any generator setting.
- **The seeds.** A unit's dataset seed is a hash of the base seed, the design name and the replicate, and its RRPP seed adds the cell. A new design name gives new datasets at the same base seed, and a committed name gives the committed ones.
- **The stage-span cells need a scale.** The parent's controls are stage-span cells. Their edits are defined in coordinates divided by the population pooled standard deviation, whose within-cell term is the reference variance. With one group wider that term is no longer the reference variance.
- **The probes** (2026-11-24, `notes.md`; scratch scripts at seeds outside every run; not results). s multiplies group B's noise standard deviation and the share is group A's.
  - **Dispersion through the production fit.** Null datasets with group B stretched about its exact cell means; two blocks, n = 200, three stages, the cross-validated rank, 100 × 199:

    | δ | share | s | `delta` | `angle` | `shape` | median rank |
    |---|---|---|---|---|---|---|
    | 2.0 | 0.5 | 1 | 0.03 | 0.01 | 0.02 | 2 |
    | 2.0 | 0.5 | 2 | 0.17 | 0.02 | 0.01 | 3 |
    | 2.0 | 0.8 | 2 | 0.44 | 0.16 | 0.10 | 3 |
    | 2.0 | 0.2 | 2 | 0.01 | 0.00 | 0.00 | 3 |
    | 0.5 | 0.5 | 1 | 0.05 | 0.04 | 0.03 | 4 |
    | 0.5 | 0.5 | 2 | 0.65 | 0.16 | 0.07 | 9 |
    | 0.5 | 0.8 | 2 | 0.71 | 0.38 | 0.10 | 5 |
    | 0.5 | 0.2 | 2 | 0.10 | 0.01 | 0.00 | 12 |

    The held-out stage AUROC stayed above 0.99 in every unit.
  - **Within-group draws.** Both nulls on one fit and seed; three blocks, a forced rank, 200 × 199 (`delta` / `angle` / `shape`):

    | δ | share | s | unrestricted | within group |
    |---|---|---|---|---|
    | 2.0 | 0.5 | 2 | 0.075 / 0.060 / 0.010 | 0.015 / 0.055 / 0.000 |
    | 2.0 | 0.8 | 2 | 0.370 / 0.190 / 0.065 | 0.030 / 0.070 / 0.000 |
    | 2.0 | 0.2 | 2 | 0.005 / 0.015 / 0.005 | 0.045 / 0.035 / 0.090 |
    | 0.5 | 0.5 | 2 | 0.530 / 0.050 / 0.030 | 0.515 / 0.050 / 0.030 |
    | 0.5 | 0.8 | 2 | 0.575 / 0.240 / 0.095 | 0.230 / 0.020 / 0.005 |

  - **Correlation.** Group B redrawn with its correlations multiplied by c at fixed variances, 30 fits per cell. The within-cell variance of its scores over group A's was 0.67 at c = 0.5 and 0.28–0.35 at c = 0 (δ = 2.0, rank 2), and 0.52–0.55 and 0.04–0.05 at δ = 0.5 (rank 4).
  - **Timing.** A two-block unit at 199 permutations took 13–14 s with four workers and 22–27 s with twelve on eight cores. The cross-validated fit is most of it.

## Goals / Non-Goals

**Goals:**
- Let the generator give one group a wider or narrower noise dispersion without moving any mean, any seed stream or any committed identifier.
- Say, under a rule fixed before any unit, whether the trajectory test holds its size when the groups' dispersions differ: with equal groups, with the wider group the smaller, and with it the larger.
- Give each statistic the scale at which it first fails, by group share, sample size and stage effect, so that a cohort can be placed.
- Record, on the same units, whether the per-feature tests and the within-group draws hold.

**Non-Goals:**
- Power. No unit carries a difference in trajectory, so what a conservative null costs is not read.
- A remedy as a package option. The within-group draws are read through the primitive and decide nothing (D7).
- A dispersion that varies by stage, by feature or in correlation structure (D2, D3).
- A dispersion estimate or test in `motco analyze`.

## Decisions

**D1. The dispersion scale stretches the drawn noise; it does not change the covariance handed to the sampler.**
- **The mechanism.** A sampling call takes an optional noise scale s. After each block is drawn, every sample's deviation from its exact cell mean is multiplied by s, in the units the noise is drawn in. For a Gaussian draw that is exactly a draw from s² times the covariance.
- **Why after the draw.**
  - It consumes no random numbers, so the stream is the one of the equal-dispersion dataset, and nothing drawn later can shift.
  - The dataset at scale s is the dataset without the scale with group B stretched about its means. Identity while absent is then structural: the absent path is the old code path.
  - Two designs that differ only in the scale are paired sample by sample at one seed, which a test can check to floating-point precision.
- **Where.** `generate_omics(..., noise_scale=)`, applied to methylation in M-values before the inverse logit, and to expression and proteomics in their native units. The cell means, and so the population trajectories and the truth indicators, are untouched.
- **Alternative: pass s²·Σ to the sampler.** Same distribution, but the factorization of a scaled matrix differs in floating point, so the pairing is only distributional and the absent path needs its own proof.
- **Alternative: a scale per cell.** It would also give dispersion by stage within a group. Nothing in this change needs it; the scalar is its special case and can stay the default spelling if it is added later.

**D2. One scalar on the second group, with the direction set by the group share.**
- **The setting.** `SemiSyntheticTrajectoryParams.group_dispersion_scale`: a positive finite s that multiplies the noise standard deviation of the second label in `group_labels`, in every block and at every stage. `None` is the default and 1 normalizes to it, the hash-omitted spelling, as φ = 0 does for `span_alignment`.
- **Why a standard-deviation multiplier.** It is what the mechanism applies. Truth records the covariance multiplier s² beside it, and the documentation gives both.
- **Why one group.** Which group is wider is exchanged by the group share: a first-group share of 0.8 makes the scaled group the smaller, 0.2 the larger. A second setting for the first group would name the same designs twice.
- **Under every mode and construction.** The scale acts on noise and the modes on means, so they compose. A zero-effect cell keeps the scale: it is part of the design, like `group_stage_sizes`, and not one of the fields the group transform ignores at zero effect.
- **Alternative: a scale drawn per feature** (a share of features wider by a random factor). It is closer to a real cohort and needs two more parameters. The scalar is the case where the scores' dispersion ratio is s² in every direction, which is what a first reading should hold fixed.

**D3. Correlation scaling gets no setting of its own.**
- **What the probe measured.** Shrinking group B's correlations toward zero at fixed variances lowers the within-cell variance of its scores: to 0.67 of group A's at half the correlations and 0.28–0.35 with none (δ = 2.0, rank 2), and to 0.52–0.55 and 0.04–0.05 at δ = 0.5 (rank 4). A likely reason is that the fitted directions weight correlated features alike, so less correlation is less variance along them.
- **What follows.** In the measured space correlation scaling is a dispersion scale below one on group B, which is a scale above one on the other group. The cross of D5 reads both directions at variance ratios of 1.56, 2.25 and 4. The half-correlation case (ratios of 1.5 to 1.9 once inverted) and the no-correlation case at δ = 2.0 (2.9 to 3.6) fall in or at the edge of that range.
- **What it leaves out.** With no correlation the ratio differs by component (0.01–0.31 at δ = 0.5), which a scalar cannot produce. The per-feature F statistics see only the variances, so a correlation difference at fixed variances does not reach them at all.
- **Disposition of the planned grid's row.** The design-guidance pilot deferred "correlation scaling" to this change. It is recorded here as covered in total dispersion and open in anisotropy. A `group_correlation_scale` would be its own change, with a second sampling path.

**D4. Under `stage_span`, the pooled scales use each group's own noise variance.**
- **Why it is needed.** The study's cells are the parent's stage-span controls (D5). Their coordinates divide each feature by the population pooled standard deviation of the generated data. With group B wider, the within-cell term of that variance is (n_A + s²·n_B) / (n_A + n_B) times the reference variance, at the realized group sizes.
- **The change.** `null_pooled_scale` and `planted_pooled_scale` take an optional multiplier of the within-cell term, absent by default, and the stage-span plant and the planted-scale fixed point pass it. Without a dispersion scale nothing is passed and both functions return what they returned.
- **What it keeps true.** An edit under the null scale is exact in the null population's pooled scale, and one under the planted scale in the planted population's, whatever the dispersion. A null plant is still the null dataset.
- **Alternative: reject the scale under `stage_span`.** The translation control would then come from the feature-set surgery, and the units could no longer be the parent's cells or reproduce its counts.
- **Alternative: leave the scales unaware.** A translation would still be a pure location offset, but its size would be stated in a scale the data do not have, and a later power study under unequal dispersion would inherit the error.

**D5. The designs are the committed three-block designs, each with dispersion and balance changed.**
- **`core`: the full cross at n = 200, δ = 2.0.** The roadmap's study, at the design most committed claims sit at.
  - **Scales:** s ∈ {1.25, 1.5, 2}, variance ratios of 1.56, 2.25 and 4. The probe read a clear excess at 2; the ladder goes down to where a reading is expected near the bound.
  - **First-group shares:** 0.5; 0.65 and 0.8, the pilot's ratios, where the wider group is the smaller; and their mirrors 0.35 and 0.2, where it is the larger. At n = 200 the smaller group has 23 or 13 samples per stage.
  - **Equal dispersion at shares 0.65 and 0.8.** Unbalance alone, which was read at pilot grade only. Without it an excess at an unbalanced design could not be assigned to the dispersion. The mirrors are the same designs with the labels exchanged and are not run.
  - **Cells:** the anchor and the translation at e = 1.0 under the planted scale, at 500 × 999, the parent's controls.
  - 17 new designs, 17,000 units.
- **`columns`: a reduced cross at four other committed designs.** The univariate baseline's n = 80 and n = 500 at δ = 2.0 and n = 200 at δ = 0.5, and the stage-effect scaling run's `band_n320`: δ = 0.05 at n = 320, calibrated to the case studies' held-out AUROC (0.664 against SEA-AD's 0.636; median rank 10).
  - **Why.** The probe read `delta` at 0.17 at δ = 2.0 and at 0.65 at δ = 0.5 in balanced designs, so the excess depends on the stage effect, and by the same argument on the cell size and the rank. Whether a cohort is exposed depends on these, so the map needs them. `band_n320` is the committed design nearest the regime a case study is in.
  - **Per column:** s ∈ {1.5, 2} at shares 0.5, 0.8 and 0.2, and equal dispersion at share 0.8. The anchor only, at 500 × 999. The translation is read in the core.
  - **References.** The balanced equal-dispersion reading of each column is its parent's committed one. None is rerun: the gate of D6 covers the code path.
  - 28 designs, 14,000 units.
- **Names.** `<committed design>_r<share>` at equal dispersion and `<committed design>_r<share>_s<scale>` otherwise, so each is new at base seed 2043.
- **Not run.** The SEA-AD cell table, which is a cohort's design with its covariates and two blocks and would need its own simulation, as ADNI's had (roadmap item 9); four stages; two blocks; scales below 1.25 or above 2.
- **Alternative: the pilot report's cross** of spread {1, 1.5, 2, 3} by ratio {0.2, 0.5, 0.8} at one design. It has no moderate ratio and a scale of 3 where 2 already fails, and it cannot say how the excess moves with sample size or stage effect.
- **Alternative: the full cross at every column.** 85 designs and 51,000 units against 45 and 32,000, for readings the reduced cross brackets.

**D6. The unit is the parent's, unchanged, and one committed design is rerun as the gate.**
- **The unit.** `univariate_baseline_check.evaluate_unit`, handed a settings block of the shape its `stage_span` profile writes, with this driver's designs. One dataset, one cross-validated production fit, the trajectory statistics with their null distributions, then the split test and the scan on the fit's standardized features at the unit's RRPP seed. `pca_matched` is the one recorded space, as in the parent.
- **The gate.** The anchor and the translation of `n200` are rerun at the committed name and seeds, 1,000 units, and every arm's counts must equal the committed rate table. It ties this driver, the changed generator and the cluster's environment to the table the columns' references are read from, and it is the end-to-end check that the generator is unchanged while the scale is absent.
- **Bound by sha256:** `settings.json`, `rate_table.csv` and `cell_table.csv` of the univariate baseline's `stage_span` profile and of the stage-effect scaling run's `band` profile. The band column's stage effect is read from the second.
- **Alternative: rerun every committed design.** About five times the gate's cost for the same code path.

**D7. The within-group draws are recorded on every unit.**
- **What they are.** RRPP with the group label as `strata`: each draw permutes the reduced model's residuals among the samples of one group. A draw then keeps each group's own dispersion, which is what the unrestricted draws pool away.
- **How they are read.** After the parent's unit returns, the driver regenerates the unit's dataset, fits production at the recorded rank with no cross-validation, and calls the primitive with the group labels at the unit's RRPP seed. The p-values follow the harness's rule.
- **Two checks.**
  - `same_statistics`: the refit's observed `delta`, `angle` and `shape` equal the parent record's exactly. A unit that differs is marked inconsistent for this arm and leaves its rates; the parent's arms are unaffected.
  - In the tests, the same path with one stratum reproduces the production p-values exactly.
- **What the arm separates.** The probe suggests two sources of excess. One is the exchange of residuals between groups, which these draws remove: at δ = 2.0 with the wider group the smaller they brought 0.37 / 0.19 / 0.065 back to 0.03 / 0.07 / 0.00. The other survives them: at δ = 0.5 with equal groups `delta` read 0.53 under the unrestricted draws and 0.515 under these. A guess, not tested here, is that the stage fit follows the wider group's noise more than the other's, so its in-sample path is longer in scores that no residual permutation refits. The arm says how much of each column's excess is of the first kind.
- **Why recorded and not judged.** A null is adopted on size, power and specificity together: the within-stage null held size and failed specificity (roadmap item 20). This run has no plant. The draws also assume equal dispersion across the stages of a group, the mirror of what the within-stage draws assume, and the probe shows they are not a full repair.
- **What it does not cover.** The split test and the scan draw their own nulls and take no strata. The recorded `pca_matched` and `total` arms read the same units with no stage fit, under the unrestricted draws only.
- **Alternative: add `group` to the harness's `permutation_strata` now.** It would put a second unvalidated null beside the one item 20 declined to adopt, before its power is known.

**D8. Size is judged in four families, and the verdict is the first that fails.**
- **The readings.** `delta`, `angle` and `shape`, each at p < 0.05, on every control cell of every new design. The single statistics are judged, not the omnibus, because a case study reads each with its own statement and the mechanism differs by statistic.
- **The rule.** `judge_rates` of the stage-span paper grade, imported unchanged: a family holds when the readings above 0.05 + 2 SE stay within the count a null run stays within, and none is above the family's Bonferroni ceiling. At 500 replicates the bound is 34 rejections and the ceiling 40 to 42.
- **The families**, over both profiles:

  | family | designs | m | allowed above the bound | ceiling |
  |---|---|---|---|---|
  | `unbalanced_equal` | equal dispersion at shares 0.65 and 0.8 (core) and 0.8 (columns) | 12 + 12 = 24 | 2 | 40 |
  | `balanced_unequal` | a scale at share 0.5 | 18 + 24 = 42 | 3 | 41 |
  | `wider_smaller` | a scale at shares 0.65 and 0.8 (core) and 0.8 (columns) | 36 + 24 = 60 | 4 | 42 |
  | `wider_larger` | a scale at shares 0.35 and 0.2 (core) and 0.2 (columns) | 36 + 24 = 60 | 4 | 42 |

- **The verdict**, in order:
  1. `reproduction_failed`: the gate.
  2. `unbalance_excess`: `unbalanced_equal` fails. Unequal groups lose size before any dispersion difference, and the other families cannot be assigned to dispersion.
  3. `dispersion_excess`: `balanced_unequal` fails. Size needs equal dispersion even with equal groups.
  4. `unbalanced_dispersion_excess`: `wider_smaller` or `wider_larger` fails. Size needs equal dispersion when the groups are unequal.
  5. `size_holds`.
- Every family is recorded whatever the verdict.
- **Why the three statistics share a family.** Twelve families of one statistic each would give a run with exact size about three times the chance of some failure that four families give. The map of D9 says which statistic failed where.
- **Why by family and not by cell.** The verdict states a condition a reader can check on a cohort (balanced or not, which group is wider). The cells are in the map.
- **A conservative reading does not fail a family.** It costs power, which this run does not read (D9).

**D9. What is reported beside the verdict.**
- **The map.** For each statistic, column and share, the smallest scale whose reading is above its family's ceiling, or that none is up to the largest scale read there.
- **Conservative readings.** Every reading below 0.05 − 2 SE (15 or fewer rejections of 500) is flagged.
- **Recorded arms.** `production` (the omnibus at 0.05 / 3), `scan`, `fdr`, `total`, `span`, `complement`, `pca_matched` and the three within-group statistics: a label per family under the same rule, and every rate per design and cell.
- **The fit.** The selected rank and the held-out AUROC by design. The probe saw the rank rise with the scale.
- None of these changes the verdict.

**D10. The setting is invisible while absent, in every identifier.**
- `group_dispersion_scale` joins `grid._SIGNATURE_OMIT_WHEN_NONE` and the study engine's dump omission, so cell ids, signatures, matched seeds, records and dumped configs are unchanged while it is `None`.
- It is not added to the fields the study engine clears or ignores at zero effect.
- Study configs accept `generator.group_dispersion_scale` through the existing generic handling of generator settings. `motco simulate` does not expose it.
- Identity is pinned twice: by the existing fixtures (`tests/test_committed_config_identity.py`, the generator's and the stage-span tests), which must pass unchanged, and by the gate of D6 at run time.

**D11. Execution follows the parent.**
- One script, `scripts/group_dispersion_check.py`, with `run --profile core|columns`, `merge --profile` and `decide`; `settings.json` bound by hash with the parents' sha256; round-robin shards; resume by complete unit key; JSONL records with seeds, statistics and counts only; an sbatch array with a profile variable and a dry-run variable, BLAS pinned to one thread, run on `ing` from a clean worktree at a pushed revision.
- The merge uses the parent's rate and cell tables and accounting, imported unchanged, and adds the within-group arm's rows and its consistency counts.
- Committed: aggregates, the decision, provenance and the report under `results/group-dispersion-<date>/`. Records stay on the cluster.

## Risks / Trade-offs

- **[The rule was written after a probe that saw rejection rates]** → The probe ran two blocks, 100 replicates and 199 permutations at seeds outside the run, and the rule is the repository's standing size rule with no threshold chosen from the probe. The scales and columns were chosen from it, and the report says so.
- **[The verdict is the one the probe predicts]** → The product is the map: the scale at which each statistic fails, by balance, sample size and stage effect, with the per-feature arms and the within-group draws on the same units.
- **[The excess has a source the run does not isolate]** → The within-group arm, `pca_matched` and `total` bracket it: what the within-group draws remove belongs to the null, and what remains in production but not in the arms without a stage fit points at the fit. Naming the mechanism is a diagnosis for a later change, as the case-study configuration's was (roadmap item 7 (c)).
- **[A scalar scale on every feature is not what a cohort has]** → The report states it beside every claim, and gives the reading as a function of the scores' dispersion ratio, which a cohort can estimate. A heterogeneous or anisotropic dispersion is open (D2, D3).
- **[The planted-scale fixed point behaves differently with one group wider]** → The smoke runs the translation at every core design, and failed units are counted per cell.
- **[The within-group refit differs from the cross-validated fit]** → `same_statistics` on every unit; inconsistent units leave that arm only.
- **[The gate fails]** → Then the generator is not unchanged while the scale is absent, or the environment drifted. The reading stops; a mismatch is diagnosed, never waived.
- **[Compute overrun]** → A local smoke measures per-unit time by design. The levers are, in order: the n = 500 column, the translation cells at s = 1.25, and the within-group arm. The projection goes to the user with a cap.
- **[The verdict is read as a finding about SEA-AD or ADNI]** → Neither cohort's dispersion was read, and neither case study rejects. The report says what a case study would have to report for the map to apply.

## Migration Plan

The package gains one opt-in generator setting (0.26.0); no default moves and no committed identifier changes.
- **On `dispersion_excess`:** the documentation states that a rejection is evidence of a difference in trajectory or in dispersion, for the statistics and scales the map gives, in balanced designs too. The interpretability statement of each failing statistic gains that clause. A case study that rejects reports the two groups' residual dispersion in the scores beside the rejection. The roadmap opens a diagnosis of the excess the within-group draws leave, and after it a change for a null or a fit that holds size, and notes that open item 3's bias correction needs a per-group trace.
- **On `unbalanced_dispersion_excess`:** the same, scoped to unequal groups, with the shares and scales of the map; balanced designs are stated to hold up to the largest scale read.
- **On `size_holds`:** the documentation states the range of scales, shares, sample sizes and stage effects over which size held.
- **On `unbalance_excess`:** unbalance at equal dispersion is diagnosed first, as its own change; the dispersion families are reported and not interpreted.
- **Under any verdict:** the recorded arms say whether the per-feature tests share the reading, and the roadmap's "Not yet established" entries on unbalance and on noise are updated.
- **Rollback:** the setting is opt-in and absent everywhere; removing it, the driver, the tests and the report restores the previous state.

## Open Questions

- The core-hour cap for the arrays is set from the smoke timing, at the approval gate.
