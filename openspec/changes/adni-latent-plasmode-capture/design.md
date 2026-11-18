# Design

## Context

See `proposal.md` (Why) and `notes.md` (the probe). The constraints that shape the design:

- **Driver conventions** (`AGENTS.md`, Validation drivers): later drivers import earlier ones unchanged and bind their committed files by sha256; settings are bound by hash; a decision carries a rule version fixed before any unit runs; records hold seeds and statistics only; the ADNI export is read in place through `adni_anchor_check`'s contract.
- **What exists.** `adni_anchor_check` builds a unit (subsample within cells, permuted labels, seeds). `adni_gain_specificity_check` has the covariate-adjusted stage configuration, the sample's scalers, the within-stage half split and the map to a fitted model's input coordinates. `LinearLatentMap` gives the fitted span's and the predictive subspace's orthonormal bases. `analyze_trajectories` with `permutations=0` fits the space and skips RRPP.
- **ADNI.** 1,501 participants on 298 features, three stages (CN, MCI, AD), table Female 237 / 310 / 125 and Male 223 / 448 / 158.

## Goals / Non-Goals

**Goals**

- An estimate of the share of ADNI's real pooled stage program inside production's fitted span, at four sample sizes, with a cohort-level interval.
- The same for the reference spaces, by stage contrast, and for pure estimation noise.
- A check of the estimator with the truth known, before any real reading is trusted.

**Non-Goals**

- A test, a rejection rate or a power statement. No RRPP draw runs.
- A group difference. The program is pooled over groups, and no analysis uses the real labels.
- The full cohort's own fit. With every participant in the sample no one is left to form the halves, so 1,501 is read only as the trend across the four points.
- New spaces. Open item 4's candidates import the estimator in their own change.

## Decisions

### D1. The estimand

For a unit, `Q` is a space's orthonormal basis fitted on the sample `A`, in the fitted model's input coordinates. `M` is the population's covariate-adjusted stage configuration in those coordinates, centred over stages (K × p). The capture is `tr(M·Q·Q'·M') / tr(M·M')`: the share of the stage program's energy inside the space.

- **Why the pooled stage program.** It is the quantity the simulations report as `pooled_stage` (0.10 at the weak design), and a stage-span plant is an edit of it, so its capture bounds what any such plant can show.
- **Why covariate-adjusted.** The trajectory test compares covariate-adjusted cell means, and on ADNI the covariates carry half of the unadjusted program (item 22, D3). The configuration is `adni_gain_specificity_check.adjusted_configuration`, imported.
- **It is conditional on the fit.** Each unit has its own `Q`. The reading at a point is the capture averaged over samples of that size.

### D2. Three disjoint parts

```
 cohort (1,501)
 +----------------+---------------------------------------+
 | A: the sample  | remainder, split within stage         |
 | fit, basis Q   |   B1 -> m1          B2 -> m2          |
 +----------------+---------------------------------------+
```

- `A` is the anchor driver's subsample at the point (within sex × stage cells, largest-remainder table).
- The remainder is split into two halves within stage by `adni_gain_specificity_check.split_within_stage`, imported. Each half's configuration is fitted on that half alone, in `A`'s scaler and the model's input coordinates.
- **Why it works.** `m_i = M + e_i`, with `e_1`, `e_2` independent of each other and of `Q`. So `E[tr(m1·Q·Q'·m2')] = tr(M·Q·Q'·M')` and `E[tr(m1·m2')] = tr(M·M')` given the fit.
- **Twenty splits per unit**, averaged. A split costs two least-squares fits and no model fit, and averaging removes most of the split's own noise.
- **Alternative: the whole remainder as one reference** (item 22). Rejected: its share mixes the program's and the noise's (proposal).
- **Alternative: cross-fit the roles of A and the halves.** Rejected: the fit must come from a sample of the point's size, and the halves' sizes differ from it.

### D3. Points and units

- **Points:** 80, 200, 500 and 750 participants. The halves hold about 710, 650, 500 and 375.
- **Units:** 500 per point, 2,000 in all, drawn by the anchor driver's functions under its seed family (base seed 2039), unchanged. A unit at 80, 200 or 500 is the anchor run's unit, and one at 750 is the gain run's.
- **Why those units.** Their fits are already pinned by committed tables (D9), and the capture joins the size and power readings on the same samples.
- **Own seeds** (base seed 2041, used by no committed driver): `reference_split` per unit and split index, and the calibration's `cohort_noise` and subsample seeds.
- **The labels do not matter.** The unit builder permutes them as the anchor does, and nothing this driver reads depends on them: the fit and the configuration are group-blind. A test pins that.

### D4. The estimate is a ratio of sums

The reading at a point is `sum over units of numerator / sum over units of denominator`.

- **Why not a per-unit ratio.** The denominator is a cross-product of two noisy configurations and is near zero in some units (the probe saw per-unit captures from −12.9 to 22.8).
- **Per-unit quantiles are recorded** for description only.
- **The estimate is not clipped.** A value outside [0, 1] is reported as it is and counts against resolution (D8).

### D5. What is read

| reading | definition |
|---|---|
| total capture | D4 on the full cross-products |
| by contrast | D4 on the stage-space cross-Gram `G = (m1·m2' + m2·m1')/2` and its in-space counterpart, projected on a stage contrast |
| noise capture | `‖(m1 − m2)·Q‖² / ‖m1 − m2‖²`, summed over units as D4 |
| naive captures | the whole remainder's configuration, and the sample's own (item 22's two readings) |
| program size | the summed denominator: the program's energy in pooled SDs, and the cross-Gram's eigenvalues |
| reliability | the halves' correlation |

- **Contrasts.** `net_displacement` is the unit contrast last stage minus first, the direction `angle` is built on. `interior` is the rest of the centred stage space (at three stages, MCI's deviation from the chord). Their energies sum to the total.
- **Noise capture.** `m1 − m2` holds no program, so its share inside a space is the share of stage-mean estimation noise the space holds. It is always defined. It says whether a space follows the program or the noise, which is the condition under which the roadmap reopens label-free factor models for open item 4.
- **Records carry the K × K Gram matrices,** so any contrast is computed at merge. They are functions of stage means' inner products and hold no participant-level value.

### D6. Spaces

| space | basis | fit |
|---|---|---|
| `pls_cv` | the fitted span (`LinearLatentMap.from_pls`) | the case study's configuration, cross-validated on the unit |
| `predictive` | the predictive subspace of the same fit (`from_pls_predictive`) | the same model |
| `pls_rank2`, `pls_rank8`, `pls_rank19` | the fitted span | forced at the rank |
| `pca` | the leading right singular vectors of the sample's model inputs, at the unit's CV rank | none |

- **`pls_cv` is the designated space.** The others are recorded.
- **Why the forced ranks.** They are the anchor's (8, the cohort's rank; 19, its simulation's) and K − 1. They show how the capture moves with rank on a real program, which the on-hold rank rule and item 21's `reference_exceeds` grade both ask.
- **Why a PCA.** Item 14's label-free reference kept 0.53 of a simulated program. On real noise it separates what stage supervision adds from what variance alone gives.
- **One cross-validation per unit.** The forced fits are single fits.

### D7. The calibration

Synthetic cohorts with the truth known, built from the export in memory:

- **Generator.** On the full cohort, the least-squares fit `features ~ stage indicators + covariate columns` gives the stage configuration `Z`, the covariate coefficients and the residual covariance `S`. A synthetic cohort keeps every participant's stage, sex and covariates and draws features as `grand mean + s·centred(Z)[stage] + covariate part + N(0, S)`.
- **Truth.** The planted centred configuration `s·centred(Z)`, so each unit's true capture is exact.
- **Strengths.** `s = 1` (the cohort's configuration as estimated), `s = 0.67` (its energy scaled by the configuration's extrapolated reliability, 0.443 from item 22's split-half 0.166 at 375), and `s = 0` (no program).
- **Size.** 24 cohorts per strength, 25 units per cohort and point: 7,200 units. The pipeline per unit is the cohort profile's.
- **What it gives.**
  - *Bias:* the mean over cohorts of estimate minus truth, per point, space and reading.
  - *The standard error:* the SD of the pooled estimate across cohorts at `s = 0.67`. Subsamples of one cohort share participants, so the spread over ADNI's own units understates it.
  - *The null level:* the SD of each energy across cohorts at `s = 0`.
- **Why Gaussian noise.** Independent rows make the three parts independent, which the estimator assumes. Resampling ADNI's residual rows would put one participant's noise in several parts.
- **Alternative: the ADNI-design simulation.** Rejected as the calibration: its blocks are InterSIM's and its CV selects rank 19 against ADNI's 8, so its spread would not be ADNI's.

### D8. The decision (rule version 1)

Fixed before any unit of either profile runs.

- **Gate.** All of:
  - no failed and no inconsistent unit beyond 1 % of a point's units, in either profile;
  - both reproductions pass (D9);
  - at `s = 0.67` and `s = 1`, the absolute bias of `pls_cv`'s total capture is at most 0.05 at every point.
- **Program present.** The total energy at `n80` (the largest halves) exceeds twice its SD at `s = 0`. Otherwise `no_program`.
- **Interval.** The estimate ± twice the SD at `s = 0.67`, at the same point, space and reading.
- **Resolved.** A point is resolved when the half-width of `pls_cv`'s total capture is at most 0.15.
- **Verdict.**

| verdict | condition |
|---|---|
| `gate_failed` | the gate fails |
| `no_program` | the gate holds and the program is not present |
| `unresolved` | no point is resolved |
| `partly_measured` | one to three points are resolved |
| `measured` | all four are resolved |

- **Contrast grades** (per point, recorded with the verdict). A contrast is `present` when its energy exceeds twice its SD at `s = 0`, and `captured` / `missed` / `undetermined` as its capture's interval lies above 0.5, below 0.5 or across it. A contrast that is not present has no capture grade.
- **What a verdict licenses.**
  - `measured` or `partly_measured`: the resolved captures are stated on the measurement-visibility page as the real-noise capture at ADNI's design points, beside the simulated 0.10. A follow-up with plants along the cross-validated program is proposed for the contrasts graded `present`; a rotation or a bend needs both.
  - `unresolved`: the driver's estimator stands, and the number waits for ROSMAP or a larger cohort. No capture is stated.
  - `no_program`: ADNI's blocks carry no stage program the halves agree on, and stage-span plants on ADNI are edits of noise.
  - `gate_failed`: nothing is read.
- **Below grade.** A decision on fewer than 500 units at a cohort point, or fewer than 24 cohorts or 25 units in the calibration, is marked `indicative`.
- **Why these numbers.** A bias of 0.05 is a third of the resolution. A half-width of 0.15 separates the probe's 0.12 at 80 participants from its 0.43 at 750.

### D9. Checks, pre-flight and reproductions

- **Per-unit checks** (a failing unit is `inconsistent` and leaves every reading):
  - `cell_table`, the anchor driver's, imported;
  - `parts`: the sample and each split's halves are pairwise disjoint and cover the cohort;
  - `basis`: each space's basis is orthonormal to 1e-8 at the expected dimension;
  - `gram`: the trace of each Gram matrix equals its scalar cross-product to 1e-8;
  - `truth` (calibration): the planted configuration is recovered exactly from the noise-free features.
- **Pre-flight** (`run` refuses to start without a pass). Cohort: the anchor driver's pre-flight units reproduce the case study's rank 8 and AUROC 0.599514. Calibration: at `s = 1` with the noise scaled to zero, every space's cross-product capture equals the true capture to 1e-8.
- **Reproductions at merge**, against parents bound by sha256:
  - the CV-rank counts at 80, 200 and 500 equal the anchor run's `point_table.csv`;
  - the medians of the naive captures at 750 equal the gain run's `cell_table.csv` (0.53603 for the remainder and 0.925475 for the sample in the fitted span; 0.108803 and 0.548427 in the predictive subspace).

### D10. Driver shape and outputs

- `scripts/adni_latent_plasmode_check.py`: `preflight`, `run`, `merge` with `--profile cohort | calibration`; `decide`. `--export-dir` is required by both profiles and has no default. Each profile writes under `<out>/<profile>/` with its own `settings.json` bound by hash, the export's five file hashes and the parents' hashes in it.
- **Shards and resume** as in the gain driver: round-robin, complete units resumed, other settings refused by key.
- **Merge writes** `capture_table.csv` (point, space, reading: estimate, numerator, denominator, units, per-unit quantiles), `point_table.csv` (ranks, AUROC, reliability, energies, eigenvalues), `accounting.json`, `reproduction.json`, and for the calibration `calibration_table.csv` (strength, point, space, reading: mean estimate, mean truth, bias, SD). `merged.jsonl` stays off the repository.
- **Decide writes** `adni_latent_plasmode_decision.json`.
- **No sbatch array.** Both profiles read the export and fit on one machine in about 25 core-hours, as the anchor did.
- **Records** hold seeds, ranks, scalars and K × K matrices: no identifier, label vector, subsample membership, split membership, configuration or synthetic feature.

## Risks / Trade-offs

- **[The interval rests on Gaussian synthetic cohorts]** → ADNI's noise has heavier tails and its participants are one draw. The report states the interval as parametric, and the units' own spread is given beside it.
- **[The calibration's program is the cohort's estimate, 0.56 of it noise in direction]** → The estimator's validity does not depend on the direction. The true capture of the calibration is not ADNI's, and the report does not compare them.
- **[The halves shrink as the sample grows]** → 375 participants per half at 750, where the denominator is noisiest. Resolution is judged per point, and `partly_measured` exists for that case.
- **[A ratio of sums is biased when the denominator is noisy]** → The gate measures it at ADNI's signal and at a weaker one.
- **[Subsamples of one cohort are not independent cohorts]** → The standard error comes from the calibration, not from the units.
- **[The probe read estimates before the rule was written]** → It ran on scratch seeds at 32 units, is recorded in `notes.md` as not a result, and the thresholds are stated with their reasons (D8).
- **[A stage statistic of real participants reaches a committed file]** → Only pooled scalars and K × K matrices are recorded; a test searches records and committed files for identifiers and per-participant or per-feature fields.
- **[The strength 0.67 is extrapolated]** → It is a setting fixed before the run. The pre-flight records the cohort configuration's own split-half reliability for comparison.
