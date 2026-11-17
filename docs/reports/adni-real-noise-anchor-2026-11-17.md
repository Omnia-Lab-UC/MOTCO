# The ADNI real-noise anchor — report (2026-11-17)

**Change:** `adni-real-noise-anchor` (MOTCO 0.25.0, no package change). **Reads:**
`results/adni-real-noise-anchor-2026-11-17/` (`adni_anchor_decision.json`, `rate_table.csv`, `paired_table.csv`,
`point_table.csv`, `accounting.json`, `preflight.json`, `settings.json`, `PROVENANCE.txt`). **Parents:** the
[within-stage permutation pilot](within-stage-permutation-2026-11-16.md) (its SEA-AD anchor tables, bound by sha256),
the [predictive size-excess diagnosis](predictive-size-excess-2026-11-13.md) and the
[predictive-subspace pilot](predictive-subspace-2026-11-12.md). **Concept page:**
[measurement visibility](../measurement-visibility.md#a-second-cohort-the-adni-anchor).

Every real-noise reading MOTCO had came from one export, SEA-AD astrocytes at 80 donors on 10,000 features. On its
anchor `predictive_orthonormal` rejected `delta` in 67 of 500 replicates where production rejected 21, and that veto
kept the predictive subspace an opt-in diagnostic (roadmap item 18). The diagnosis (item 19) traced the excess to a
stage fit that is mostly in-sample, read under RRPP's unrestricted draws, and predicted that it would shrink as
participants per feature grow and be absent on ADNI. This change reads that prediction, and production's size, on
ADNI's own blocks, under a decision fixed before any unit ran.

**In one paragraph.** **The verdict is `excess_absent`.**

- **Size holds for both operating arms everywhere.** At 80, 200 and 500 participants subsampled within cells and at
  the full cohort of 1,501, neither production (`pls_cv`) nor `predictive_orthonormal` has a reading above the bound
  (34 of 500) under RRPP's unrestricted draws. Their largest counts are 32 and 29 of 500.
- **The consequence fixed with the rule.** SEA-AD's excess does not transfer to ADNI. The size part of open item 1's
  second reopening condition is met for ADNI. `predictive_orthonormal` stays an opt-in diagnostic until the gain and
  specificity readings at ADNI's design, which are the next change.
- **The fixed-rank grade is `reference_exceeds`.** Forced to rank 8 or 19, the scores themselves exceed at
  n = 200 (`shape` 41 and 38 of 500; 41 is above the family's ceiling of 40). So the predictive arms carry no grade,
  and this run does not say whether the hold is owed to ADNI's noise or to the low ranks its CV selects.
- **Recorded beside the verdict.** Production under the within-stage draws has no size scope: at the full cohort it
  rejects `delta` 35 and `angle` 38 of 500, above that point's ceiling of 36. `predictive_orthonormal` under the
  within-stage draws holds everywhere.

## 1. Configuration and provenance

- **Export.** `adni_sexdiff/results/motco_input`, read in place through `ADNI_CONTRACT` and never copied: 1,501
  participants, the `p180` (133) and `nmr` (165) metabolite blocks, Female 237 / 310 / 125 and Male 223 / 448 / 158
  at CN / MCI / AD. Its five file hashes are in `settings.json` and were unchanged after the run.
- **Design points** (500 replicates × 999 permutations each, rejection at p ≤ 0.05):

  | point | participants | Female CN / MCI / AD | Male CN / MCI / AD | features per participant |
  |---|---|---|---|---|
  | `n80` | 80 | 13 / 16 / 7 | 12 / 24 / 8 | 3.7 |
  | `n200` | 200 | 31 / 41 / 17 | 30 / 60 / 21 | 1.5 |
  | `n500` | 500 | 79 / 103 / 42 | 74 / 149 / 53 | 0.6 |
  | `cohort` | 1,501 | 237 / 310 / 125 | 223 / 448 / 158 | 0.2 |

  A subsampled unit draws its participants within cells and has its own fit. The cohort's 500 units share one fit and
  differ only in the labels, as SEA-AD's anchor does.
- **The null.** Each unit permutes the sex labels within stage and leaves blocks and covariates on their rows. No
  analysis ran under the real labels.
- **Arms.** Each arm is one `analyze_trajectories` call at the case study's configuration (five covariates, rank
  floor 2, net displacement, the harness CV defaults) with the `df` residual correction, all at the unit's RRPP seed.
  - **Operating:** `pls_cv` (the cross-validated fit; at the cohort, forced at rank 8) and `predictive_orthonormal` on
    the same fit, each also under the within-stage draws (`__within_stage`).
  - **Fixed rank:** the same four combinations forced at ranks 8 and 19. Rank 8 is not repeated at the cohort.
  - **Control:** `pls_cv__uncorrected`, production without the `df` correction.

  A subsampled unit holds 13 arms and one cross-validation; a cohort unit holds 9 and none.
- **Code and host.** Revision `2099db4` on `proposal/adni-real-noise-anchor`, pushed before any unit ran, from a clean
  detached worktree. The workstation, 7 workers, Python 3.11.15, numpy 2.3.5, scipy 1.16.3, scikit-learn 1.8.0,
  pandas 2.3.3, BLAS pinned to one thread.
- **Compute.** 8.12 recorded core-hours against an approved cap of 10. A service restart on the workstation killed the
  first launch at 658 units; their records were complete, and the resumed run evaluated the other 1,342 under the same
  settings and seeds (`PROVENANCE.txt`).

## 2. Accounting and the pre-flight

- **Units.** 2,000 of 2,000 complete: no failed unit, no failed arm, no inconsistent unit and no duplicate. No
  subsample needed a redraw for an absent covariate level, and no label draw reproduced the real labels.
- **Checks.** On every unit each forced arm reported its rank, every arm fitted at one rank had the same model, each
  predictive arm measured two dimensions, each within-stage arm and the control observed exactly their counterpart's
  statistics, and the realized cell table was the point's.
- **Pre-flight: passed.** On cohort replicates 0 and 1 the cross-validated fit selected rank 8 at a held-out stage
  AUROC of 0.599514, the values of `adni_sexdiff`'s primary record (MOTCO 0.10.0), and the call forced at rank 8 gave
  identical statistics and p-values. The export, the row order and the environment reproduce the case study's fit.

## 3. The decision

Rule version 1 (`design.md` D8–D11).

- **Size scope.** An arm is judged on the ladder's upper sets (all four points; `n200` and above; `n500` and above;
  the cohort alone) by the count-plus-ceiling rule over each set's `delta`, `angle` and `shape` readings as one family.
  Its scope is the first set that holds.
- **Verdict.** Read on the operating arms under the unrestricted draws: `gate_failed`, then `production_excess`
  (production's scope does not start at `n80`), then `excess_persists` (the candidate has no scope), then
  `excess_scoped`, else `excess_absent`.

| arm | all four points (m = 12) | `n200` and above (9) | `n500` and above (6) | cohort (3) | scope |
|---|---|---|---|---|---|
| `pls_cv` | 0 above the bound | 0 | 0 | 0 | `n80` |
| `predictive_orthonormal` | 0 above the bound | 0 | 0 | 0 | `n80` |
| `pls_cv__within_stage` (recorded) | 2 above (allowed 1) | 2 | 2, one above the ceiling | 2, one above the ceiling | none |
| `predictive_orthonormal__within_stage` (recorded) | 0 above the bound | 0 | 0 | 0 | `n80` |

**The verdict is `excess_absent`.** The gate passed and the decision is at paper grade (not `indicative`). The
consequence fixed with the rule: "predictive_orthonormal holds size on ADNI's blocks at every design point under
RRPP's unrestricted draws: SEA-AD's excess does not transfer to ADNI. The size part of open item 1's second
condition is met for ADNI. It stays an opt-in diagnostic until the gain and specificity readings at ADNI's design,
the next change. The fixed-rank grade says whether the hold is owed to ADNI's noise or to the ranks its CV selects.
In every verdict the SEA-AD anchor's reading stands for SEA-AD's design, and keep_production stands."

## 4. Size of the operating arms

Rejections of 500 for `delta` / `angle` / `shape`:

| arm | `n80` | `n200` | `n500` | `cohort` |
|---|---|---|---|---|
| `pls_cv` | 14 / 12 / 23 | 26 / 20 / 26 | 17 / 20 / 20 | 25 / 32 / 28 |
| `predictive_orthonormal` | 10 / 16 / 22 | 25 / 10 / 27 | 17 / 21 / 29 | 29 / 20 / 18 |
| `pls_cv__within_stage` | 16 / 11 / 20 | 30 / 24 / 25 | 25 / 28 / 21 | **35 / 38** / 28 |
| `predictive_orthonormal__within_stage` | 13 / 15 / 20 | 32 / 18 / 28 | 18 / 33 / 30 | 27 / 28 / 21 |
| `pls_cv__uncorrected` (control) | 21 / 16 / 31 | 29 / 23 / 32 | 18 / 24 / 22 | 26 / 33 / 28 |

- **Production holds size on a second cohort.** Its rates run from 0.024 to 0.064, and the case study's design
  (the cohort point) reads 0.050 / 0.064 / 0.056.
- **`predictive_orthonormal` shows no excess.** Its largest count is 29 (0.058), against SEA-AD's 67 at its anchor.
  In paired tests against production on the same units, the smallest two-sided p-value is 0.041 (`n200` `angle`:
  15 units production alone rejects, 5 the candidate alone).
- **The cohort's rates are conditional on one fit**, as SEA-AD's anchor rates are. `n500` is the largest point whose
  fits vary.

## 5. The within-stage arms (recorded beside the verdict)

- **`pls_cv__within_stage` has no size scope.** At the cohort it rejects `delta` 35 and `angle` 38 of 500, both above
  the bound and `angle` above the ceiling of 36 for a one-point family. Paired against the unrestricted draws on the
  same units, the within-stage draws reject more at the larger points: `delta` 0 units against 10 at the cohort
  (two-sided p = 0.002) and 0 against 8 at `n500` (p = 0.008).
- **`predictive_orthonormal__within_stage` holds everywhere**, its largest count 33 (`n500` `angle`).
- **What the rule does with it.** Nothing for the verdict, which reads the unrestricted draws (design D9). The
  documentation of the opt-in null records the failing readings
  ([measurement visibility](../measurement-visibility.md#the-permutation-null-unrestricted-or-within-stage)). This
  change explains no failure.

## 6. The fixed-rank grade: `reference_exceeds`

One family of 21 readings per combination of coordinates and null: three statistics at ranks 8 and 19 on the three
subsampled points, and at rank 19 on the cohort. It holds with at most two readings above the bound and none above
the ceiling (40 of 500). Rejections of 500, with the median in-sample stage R² of each fit:

| arm | `n80` | `n200` | `n500` | `cohort` |
|---|---|---|---|---|
| in-sample stage R², rank 8 | 0.596 | 0.320 | 0.164 | (operating) |
| in-sample stage R², rank 19 | 0.935 | 0.603 | 0.315 | 0.157 |
| `pls_rank8` | 18 / 23 / 27 | 23 / 26 / **41** | 15 / 28 / 32 | — |
| `pls_rank19` | 16 / 28 / 29 | 24 / 30 / **38** | 15 / 32 / 27 | 25 / 33 / 29 |
| `pls_rank8__within_stage` | 20 / 14 / 25 | 29 / 27 / 32 | 21 / 26 / 26 | — |
| `pls_rank19__within_stage` | 18 / 16 / 25 | 30 / 30 / 35 | 19 / 28 / 25 | 32 / 37 / 28 |
| `predictive_orthonormal_rank8` | 22 / 24 / 26 | 32 / 21 / 29 | 26 / 16 / 25 | — |
| `predictive_orthonormal_rank19` | 33 / 26 / 22 | 35 / 27 / 18 | 28 / 20 / 24 | 36 / 22 / 20 |
| `predictive_orthonormal_rank8__within_stage` | 22 / 23 / 24 | 29 / 25 / 27 | 29 / 19 / 30 | — |
| `predictive_orthonormal_rank19__within_stage` | 30 / 24 / 20 | 30 / 25 / 18 | 28 / 23 / 28 | 28 / 24 / 21 |

- **The grade.** The scores under unrestricted draws fail their family: `pls_rank8` `shape` at `n200` reads 41,
  above the ceiling, and `pls_rank19` `shape` there reads 38. Both are on the same 500 subsamples. The statement fixed
  with the rule: "An over-fitted rank mis-sizes the scores themselves on ADNI's noise. The predictive arms are
  reported without a grade, and the roadmap records that production's size depends on a rank its CV does not select
  here."
- **The predictive arms, reported without a grade.** Their unrestricted family has two readings above the bound
  (`delta` 35 at `n200` and 36 at the cohort, both at rank 19) and none above the ceiling; their within-stage family
  has none above the bound.
- **What the grade cannot say.** It was to tell whether the operating hold is owed to ADNI's noise or to the low
  ranks its CV selects. Because the reference fails first, it says neither. Whether an over-fitted stage fit on ADNI's
  noise gives the predictive subspace SEA-AD's excess stays unread. The in-sample fits these ranks produce, from 0.157
  to 0.935, bracket SEA-AD's 0.840; that describes the ladder and licenses no reading (S5).

## 7. Descriptive readings

None changes the verdict or the grade.

**The stage fit by point** (medians over each point's 500 units; the cohort's AUROC is the pre-flight's):

| point | CV rank 2 | CV rank median | held-out AUROC | in-sample R², CV / 8 / 19 | out-of-fold share, CV / 8 / 19 |
|---|---|---|---|---|---|
| `n80` | 287 of 500 | 2 | 0.521 | 0.152 / 0.596 / 0.935 | −0.174 / −0.595 / −1.242 |
| `n200` | 228 | 3 | 0.539 | 0.088 / 0.320 / 0.603 | −0.083 / −0.293 / −0.966 |
| `n500` | 135 | 5 | 0.565 | 0.085 / 0.164 / 0.315 | −0.032 / −0.097 / −0.335 |
| `cohort` | 0 | 8 | 0.5995 | 0.086 / (CV) / 0.157 | 0.005 / (CV) / −0.046 |

- **ADNI's stage signal is weak at every point.** The out-of-fold share is negative at every subsampled point: a fit
  predicts held-out participants' stages worse than their mean. At the CV rank the in-sample fit is
  0.085–0.152, against SEA-AD's 0.840.
- **The CV rank grows with n.** Rank 2, where the predictive subspace is the scores' own two dimensions, is selected
  in 57 % of units at `n80`, 46 % at `n200` and 27 % at `n500`. Rank 19 is selected in 1, 12 and 11 units.

**Rates by rank class** (operating arms; `delta` / `angle` / `shape` rejections among units at CV rank 2 and above
it):

| point | `pls_cv`, rank 2 | `pls_cv`, above | `predictive_orthonormal`, rank 2 | `predictive_orthonormal`, above |
|---|---|---|---|---|
| `n80` | 8 / 4 / 15 of 287 | 6 / 8 / 8 of 213 | 5 / 5 / 14 of 287 | 5 / 11 / 8 of 213 |
| `n200` | 13 / 5 / 14 of 228 | 13 / 15 / 12 of 272 | 12 / 5 / 12 of 228 | 13 / 5 / 15 of 272 |
| `n500` | 7 / 2 / 4 of 135 | 10 / 18 / 16 of 365 | 5 / 1 / 6 of 135 | 12 / 20 / 23 of 365 |

Above rank 2, where the two coordinate systems differ, `predictive_orthonormal` reads at most 23 of 365 (0.063,
`n500` `shape`).

**The control arm.** Without the `df` correction production rejects more at the small points: every discordant unit is
one the uncorrected arm rejects alone (`n80`: 7 / 4 / 8 units; the cohort: 1 / 1 / 0). It stays under the bound
everywhere, its largest count 33 at the cohort, the configuration the case study was issued under. On ADNI's noise
`df` makes the null more conservative without being needed for size (open item 9).

**Residual dispersion by stage** (median, each stage over the pooled one, at the CV fit; CN / MCI / AD):

| point | scores | `predictive_orthonormal` |
|---|---|---|
| `n80` | 1.037 / 1.016 / 0.877 | 1.031 / 1.018 / 0.871 |
| `cohort` | 1.018 / 1.016 / 0.927 | 1.000 / 1.017 / 0.956 |

The narrowest stage is AD, the smallest; MCI, the largest, sits at the pooled spread in both coordinate systems. On
SEA-AD the largest stage was the narrowest in the predictive coordinates, the pattern under which item 20's Gaussian
profiles over-reject.

**The `angle` null's 95th percentile** (median): production 94.3° / 77.2° / 55.9° / 37.5° and
`predictive_orthonormal` 108.0° / 89.6° / 44.7° / 21.9° from `n80` to the cohort.

**The two cohorts side by side** (rejections of 500, `delta` / `angle` / `shape`; SEA-AD's from the within-stage
run's committed table):

| cohort, point | features per participant | fit rank | in-sample stage R² | `pls_cv` | `predictive_orthonormal` | `pls_cv__within_stage` | `predictive_orthonormal__within_stage` |
|---|---|---|---|---|---|---|---|
| SEA-AD `anchor` | 125 | 5 | 0.840 | 21 / 24 / 13 | 67 / 9 / 38 | 19 / 27 / 19 | 29 / 17 / 32 |
| SEA-AD `anchor_fresh` | 125 | 5 | 0.840 | 32 / 24 / 12 | 69 / 4 / 37 | 32 / 26 / 23 | 36 / 15 / 26 |
| ADNI `n80` | 3.7 | 2 | 0.152 | 14 / 12 / 23 | 10 / 16 / 22 | 16 / 11 / 20 | 13 / 15 / 20 |
| ADNI `n200` | 1.5 | 3 | 0.088 | 26 / 20 / 26 | 25 / 10 / 27 | 30 / 24 / 25 | 32 / 18 / 28 |
| ADNI `n500` | 0.6 | 5 | 0.085 | 17 / 20 / 20 | 17 / 21 / 29 | 25 / 28 / 21 | 18 / 33 / 30 |
| ADNI `cohort` | 0.2 | 8 | 0.086 | 25 / 32 / 28 | 29 / 20 / 18 | 35 / 38 / 28 | 27 / 28 / 21 |

Features per participant, the ranks and the in-sample fit describe each point; they license no reading at a design
that was not read (S5).

## 8. What follows

- **The size part of open item 1's second reopening condition is met for ADNI.** `predictive_orthonormal` stays an
  opt-in diagnostic. The next change reads its gain over production at ADNI's design in simulation and its specificity
  there, in simulation and on plants in a held-out reference span (open item 2's next slice), each under its own
  predeclared decision.
- **The veto stands for SEA-AD's design.** By the scoping rule (S1–S5, written into the decision), a size reading
  speaks for its export and design points, and cohorts do not vouch for one another. `keep_production` stands.
- **Production has a second real-noise size reading.** It holds on ADNI's blocks at every point, at the case study's
  configuration with `df` and, descriptively, without it.
- **Production's size at a forced rank.** At the ranks its CV does not select on ADNI, the scores exceed at `n200`.
  The roadmap records this (the grade's statement). No case study forces a rank, and the cross-validated production
  arm holds.
- **The within-stage null fails at ADNI's cohort.** Its failing readings join the opt-in null's documentation.
- **Item 19's mechanism is not tested on ADNI.** The grade could not say whether an over-fitted fit on ADNI's noise
  reproduces SEA-AD's excess. Reading that needs a reference that holds at the forced ranks.

## 9. Limitations

- **Subsamples of one cohort are not independent cohorts.** Every rate is conditional on ADNI's 1,501 participants, and
  at `n500` two subsamples share about a third of their participants. Independent real noise needs ROSMAP.
- **The cohort's 500 units share one fit.** Its rates, the within-stage arm's failure included, are conditional on
  that fit.
- **At small n the CV mostly selects rank 2**, where the candidate measures the scores' own subspace. The rates by rank
  class show the candidate holding above rank 2 as well, on fewer units.
- **The family rule has a false-alarm rate.** A correctly sized arm fails the ladder family about 7 % of the time
  (design D8),
  and two arms decide. The verdict does not rest on a reading near the ceiling: both operating arms have no reading
  above the bound.
- **The fixed ranks are not what a case study runs.** Their grade characterizes a mechanism and cannot change the
  verdict.
- **One omic type.** ADNI's blocks are metabolites at 0.2 to 3.7 features per participant; SEA-AD's were RNA and ATAC
  at 125. Nothing here says what SEA-AD's blocks add.
