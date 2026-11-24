# Design

## Context

See proposal.md (Why) for the motivation. The facts below shape the approach.

- **What the parent committed.** The univariate baseline's `stage_span` profile ran four designs on methylation, expression and proteomics (658 features), three stages, equal groups, base seed 2043: `n80`, `n200` and `n500` at δ = 2.0 and `n200_delta0.5`. Each has thirteen cells, the two controls at 500 × 999 and the eleven signal cells at 200 × 999. Its `rate_table.csv` holds the available and rejected counts per design, cell and arm, and its `cell_table.csv` the rank, AUROC and capture. The per-unit records no longer exist.
- **How the generator scales.** Each omic's cell means are a base plus δ times the differential indicators, drawn over a noise covariance that does not depend on δ. Under the planted scale a stretch multiplies group B's centred stage configuration by 1 + e, so the planted difference is e·δ times a fixed pattern.
- **Where the scaling is only approximate.**
  - The pooled standard deviation that standardizes a feature includes its between-stage variance, which grows as δ². It is largest at δ = 2.0.
  - The residual degrees of freedom differ between matched designs (74 at n = 80, 1,274 at n = 1,280). The per-feature arms read far tails, where that matters.
  - Rotations and bends are not linear in e.
- **The seeds.** A unit's dataset seed is a hash of the base seed, the design name and the replicate, and its RRPP seed adds the cell. A new design name gives new datasets at the same base seed.
- **The cross-validation.** Three inner and four outer folds, five repeats, at most 20 components, so the largest rank it can select is 19.
- **The timing probe** (2026-11-24). Two anchor units per design through the parent's unit, at 999 permutations and base seed 90001, on the workstation with seven workers and BLAS pinned:

  | design | δ√n | CV rank | held-out AUROC | pooled-stage capture | seconds per unit |
  |---|---|---|---|---|---|
  | δ = 1.0, n = 320 | 17.9 | 2, 3 | 1.000 | 0.974, 0.978 | 24–26 |
  | δ = 1.0, n = 800 | 28.3 | 2, 2 | 1.000 | 0.985, 0.988 | 40 |
  | δ = 0.5, n = 1,280 | 17.9 | 3, 3 | 1.000 | 0.984, 0.985 | 75–77 |
  | δ = 0.25, n = 800 | 7.1 | 5, 5 | 1.000 | 0.936, 0.937 | 39–42 |
  | δ = 0.1, n = 320 | 1.8 | 11, 10 | 0.902, 0.929 | 0.563, 0.595 | 25–26 |
  | δ = 0.1, n = 1,280 | 3.6 | 19, 19 | 0.998 | 0.846, 0.848 | 76 |
  | δ = 0.05, n = 1,280 | 1.8 | 19, 19 | 0.905, 0.873 | 0.670, 0.592 | 71–75 |

  The committed anchors read rank 2 and capture 0.934 at δ = 2.0, n = 80; rank 2 and 0.973 at n = 200; and rank 3 and 0.905 at δ = 0.5, n = 200. The cross-validated fit is 85–90 % of every unit. The probe saw null datasets only, so it carries nothing about power or naming.

## Goals / Non-Goals

**Goals:**
- Say whether a weaker stage effect costs the trajectory test more than sample size buys back, under a rule fixed before any unit.
- Give "restored" a reference: a committed design with the same δ√n, read on the same cells.
- Separate what the stage fit loses from what any test loses, with an arm that uses no fit.
- Report every design by held-out AUROC, rank and δ√n, so that a cohort can be placed among them, and read designs at the case studies' AUROC.

**Non-Goals:**
- A power curve in n at every δ. The designs are chosen as matched pairs, not as a rectangular grid.
- Showing equivalence. At these replicate counts a difference below the margin cannot be excluded (D6).
- Choosing or testing a rank rule.

## Decisions

**D1. The designs are matched pairs in δ√n, each new design against a committed one.**
- **Why pairs.** A grid in δ and n gives curves and no statement of what "restored" means. A pair gives each weakly staged design a strongly staged one that, by the generator's scaling, holds the same information on the same cells.
- **The step.** Halving δ against four times n keeps δ√n, and with the same effect e it keeps e·δ·√n too. The new sample sizes are therefore 320, 800 and 1,280.
- **The four pairs.**
  - **δ = 0.5 at n = 1,280 against δ = 2.0 at n = 80.** The committed weak column's stage effect, at the sample size that matches the smallest committed strong design. This is the question as asked.
  - **δ = 1.0 at n = 320 against the same design.** The intermediate step on that diagonal, so a cost can be seen to grow as δ falls.
  - **δ = 1.0 at n = 800 against δ = 2.0 at n = 200.** A second diagonal, at the sample size most committed claims sit at.
  - **δ = 0.25 at n = 800 against δ = 0.5 at n = 200.** It matches the committed weak column from below. If that column's leak into `shape` belongs to δ√n = 7.1 the pair agrees; if it belongs to δ the weaker design is worse.
- **Not run.** The match of δ = 2.0 at n = 500 is δ = 1.0 at n = 2,000: its reference is saturated in most cells and the unit is the costliest.
- **Alternative: a rectangular grid** of δ ∈ {1.0, 0.5, 0.25} by n ∈ {320, 800, 1,280}. It costs about twice as much, and five of its nine designs have no committed match.
- **Alternative: matching by held-out AUROC.** It is at its ceiling for every δ ≥ 0.25, so it cannot match these designs. It is read below the ceiling instead (D4).

**D2. The committed designs are the references, and one is rerun as the gate.**
- **Why the committed counts are enough.** The two designs of a pair are different datasets, so the comparison is unpaired and needs counts only.
- **The gate.** `n200_delta0.5` is rerun at its committed seed and name and must reproduce every committed count of every arm. It is both a reference and the weak column the question starts from. That ties this driver and its environment to the table the pairs are read against.
- **Bound by sha256:** the parent's `settings.json`, `rate_table.csv` and `cell_table.csv`.
- **Alternative: rerun all four references.** It adds about 45 core-hours and no information, since the code path is the same.

**D3. The cells are the parent's in `saturated`, and an extended ladder in `band`.**
- **`saturated`.** The thirteen cells unchanged, so every pair-cell compares the same plant: the anchor, a translation, stretches at e ∈ {0.05, 0.1, 0.25}, in-span rotations and bends at e ∈ {0.1, 0.25, 0.5} under the planted scale, and the out-of-span rotation and its 20-feature form under the null scale.
- **`band`.** At the case studies' AUROC the parent's ladder would sit at the null. Stretches at e ∈ {0.5, 1.0, 2.0}, rotations and bends at e ∈ {0.25, 0.5, 1.0}, the anchor and a translation, all under the planted scale. e = 2.0 is the largest stretch ADNI's gain run read. The out-of-span cells are left out: production is blind there at any stage effect.
- **Grade.** Controls at 500 × 999 and signal cells at 200 × 999, the parent's.

**D4. The band's stage effect is calibrated to the held-out AUROC.**
- **The step.** `calibrate` fits production, with no permutation, on 40 anchor datasets per value of δ ∈ {0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.1} at n = 320, from the comparison's calibration seed stream. δ_c is the value whose median held-out AUROC is nearest 0.636, and that median must lie in 0.60–0.67. Both numbers are the latent-space comparison's, which calibrated its weak design to SEA-AD.
- **Outside the band.** The driver refuses to build the `band` settings. The ladder is amended in `notes.md` before any `band` unit.
- **The three designs.**
  - **δ_c at n = 320.** The calibrated design.
  - **δ_c / 2 at n = 1,280.** Its matched design. Its AUROC is not calibrated: whether it lands on the first design's is a reading of whether the AUROC follows δ√n.
  - **δ_c at n = 1,280.** The same stage effect with four times the samples, the question as a cohort would put it.
- **Why calibrate at n = 320.** A weak held-out AUROC is what a cohort of a few hundred reports. Calibrating at 1,280 would put the smaller design below any signal.
- **What the band does not match.** Neither cohort's cell table, covariates or rank. ADNI's own design simulation selects rank 19 at an AUROC of 0.606 (roadmap item 9).

**D5. The unit is the parent's, unchanged.**
- One dataset, one cross-validated production fit, the trajectory statistics with their null distributions, then the split test and the scan on the fit's standardized features at the unit's RRPP seed. `production` is the smallest of the three p-values at 0.05 / 3. The recorded arms are `span`, `complement`, `production_minp` and `pca_matched`; `concat` stays off.
- Each arm keeps the rejection rule the parent recorded, which the gate needs.

**D6. The matched comparison.**
- **The test.** For a pair-cell and an arm, Fisher's exact test on the two designs' rejected and not-rejected counts, one-sided in each direction.
- **Labels.** With 53 pair-cells (four pairs of eleven signal cells and the band's pair of nine):
  - `weak_costs`: the weaker-staged design's rate is lower by at least 0.10 and the p-value is at or below 0.025 / 53;
  - `weak_gains`: it is higher by at least 0.10 at the same level;
  - `unresolved`: otherwise.
- **The control arm.** `fdr` uses no stage fit, was admissible in the parent and draws no permutation. A cost or gain that `fdr` shares in the same pair-cell belongs to the pairing, not to the trajectory test, so `production`'s label counts toward the verdict unless `fdr` carries the same label there.
  - **Where the control is expected to move.** The design with more samples has more residual degrees of freedom, which helps a far-tail test a little. That would show as `weak_gains` on `fdr`.
  - **Why not `total` or `scan`.** `total` sits near its level on the stretches and bends, and `scan` was inadmissible at n = 80.
- **The verdict is by existence**, over the counted labels:
  - `weak_stage_costs`: at least one `weak_costs` and no `weak_gains`;
  - `weak_stage_gains`: the reverse;
  - `mixed`: both occur;
  - `scaling_unrefuted`: neither occurs.
- **Why "unrefuted".** With 200 replicates in each design the difference of two rates has a standard error up to 0.05. A true difference below the margin cannot be excluded, so the absence of a label is not agreement.
- **The out-of-span cells stay in the pairs.** The scaling applies to them as to any plant.
- **Order.** Gate, size, the matched comparison. Every criterion is recorded whatever the verdict.
- **Alternative: an equivalence test** at a margin of 0.10. At this familywise level it would need about 900 replicates per cell in both designs.
- **Alternative: a regression of the rate on δ and n.** It needs the rectangular grid of D1 and a model for the curve.

**D7. Naming is graded per design.**
- **The grade.** The parent's: among a design's in-span cells that `production` detects in at least 0.80 of units, the two off-target statistics of each are judged as one family under `judge_rates`. `names_kind` when the family holds, `detects_only` when it fails, `no_detected_cell` when none qualifies.
- **Per design.** The parent pooled its four designs into one family. Here each design is its own family, and the four committed designs are graded again the same way from the committed table, so a pair compares like with like.
- **How it is read.** A table of the grade by design, with every off-target rate above its bound. It never changes the verdict.

**D8. The map is read at each new design.**
- The parent's labels (`production_leads` | `baseline_leads` | `unresolved`) by its paired exact tests against the admissible per-feature arms, with its levels and margin, over this run's 71 signal cells (44 in `saturated`, 27 in `band`).
- It answers whether the trajectory test still leads on stretches and small bends where the stage fit is weak. It is reported by plant family and design and does not change the verdict.

**D9. What is descriptive.**
- **The index table:** per design δ, n, δ√n, the held-out AUROC, the CV rank quantiles and the capture shares, the committed designs included.
- **Every arm's rate per cell**, the three single statistics and the recorded arms included.
- **Rates against e·δ·√n:** `production`, each statistic and `fdr` on the stretch cells of every design.
- **Rank strata:** per pair-cell with a `weak_costs` label, the weaker-staged design's rates by selected rank.
- **The calibration:** each ladder value's AUROC and rank quantiles.

**D10. Execution follows the parent.**
- One script with `calibrate`, `run --profile`, `merge --profile` and `decide`; `settings.json` bound by hash with the parent's sha256 and, for `band`, the calibration's; round-robin shards; resume by complete unit key; JSONL records with seeds and statistics only; an sbatch array with a profile variable and a dry-run variable, BLAS pinned to one thread, run on `ing` from a clean worktree at a pushed revision.
- Committed: aggregates, the calibration, the decision, provenance and the report under `results/stage-effect-scaling-<date>/`. Records stay on the cluster.

## Risks / Trade-offs

- **[The scaling fails for reasons that have nothing to do with the stage fit]** → The control arm marks those pair-cells, and their labels do not count. The report states how often it moved.
- **[The decision was written after a probe of the designs]** → The probe saw two null datasets per design at a seed outside the run: rank, AUROC, capture and timing, and no rejection rate. The report says so.
- **[The band's cells sit at the null for every arm]** → That is a reading: what is detectable at the case studies' AUROC. The ladder reaches the largest stretch read on ADNI's noise.
- **[The calibration lands outside the band]** → The driver refuses the `band` settings, and the ladder is amended before any unit.
- **[The planted fixed point fails at e = 2.0 or at a weak stage effect]** → Failed units are recorded and counted per cell, and the smoke runs every design. A design with more than a few failures is amended before the arrays.
- **[Compute overrun]** → A local smoke measures per-unit time by design. The levers are, in order: the two `band` designs at n = 1,280 (117 core-hours at the probe's cost), the δ = 1.0, n = 800 pair, and the signal replicates. The projection goes to the user with a cap.
- **[The verdict is read as a statement about a real cohort]** → The report states beside every claim that it holds at 658 simulated features with Gaussian noise, equal groups and three stages.
- **[Reproduction fails from environment drift]** → The gate stops the reading. A mismatch is diagnosed, never waived.

## Migration Plan

Nothing migrates. No package code changes and no default moves.
- **On `weak_stage_costs`:** the documentation states that a weak stage effect is not bought back by sample size at the matched rate, with the plant families and the size of the loss. Guidance is indexed by the stage effect and the sample size separately. If the costs sit in the strata where the rank inflates, the rank rule under a fixed measured dimension is reopened (roadmap, On hold).
- **On `scaling_unrefuted`:** the documentation may state that no cost beyond δ√n was resolved at these designs, with the margin, and may place a design by n·(δ / 2)² as an approximate effective sample size.
- **On `weak_stage_gains`:** the documentation states where a weakly staged design with a matched sample reads better than its match.
- **On `mixed`:** the documentation states the map of costs and gains by plant family.
- **Under any verdict:** the grade by design says where a rejection names the kind of a difference, and the `band` readings say what is detectable at the case studies' AUROC.
- **Rollback:** remove the driver, its test and its report.

## Open Questions

- The core-hour cap for the arrays is set from the smoke timing, at the approval gate.
