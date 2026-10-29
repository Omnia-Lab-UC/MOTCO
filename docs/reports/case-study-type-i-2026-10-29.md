# Case-study configuration Type I check — 2026-10-29

> Paper-grade size check of the real-data analysis (`motco analyze`, MOTCO 0.12.0) at the SEA-AD case study's
> configuration — the rank floor of 2, six numeric covariate design columns, a weak stage signal and a nuisance
> covariate the covariate-blind basis follows (`openspec/changes/case-study-configuration-type-i/`, spec
> `case-study-configuration-validation`). It is the measurement the case-study section cites for its primary
> analysis and for its `stage` latent-adjusted sensitivity arm. Every number below is read from
> `results/case-study-type-i-2026-10-29/` (`delta0.1/` and `delta0.05/`: `case_study_type_i_decision.json`,
> `case_study_type_i_rates.csv`, `case_study_type_i_replicates.csv`, `case_study_type_i_pvalues.csv`, `merge.json`,
> `settings.json`, `PROVENANCE.txt`; `case_study_type_i_overall_decision.json`; `calibration/calibration.json`;
> `smoke-delta0.1-sz10/`, `smoke-delta0.05-sz10/`; `merged.jsonl` and `records/` stay on `ing` as each
> `PROVENANCE.txt` records). Rates are rejection rates at p ≤ 0.05 over the R = 500 available replicates at 999
> permutations; the bound is the one-sided `α + 2·√(α(1 − α)/R)` = 0.0695, i.e. at most 34 rejections in 500.

**Run:** two 50-shard SLURM arrays on `ing` (partition `512x1024`, n11 excluded) from the clean worktree of
`feat/case-study-configuration-type-i` at `ffd6249` (tree hash `ececcb83…`, identical to the workstation): job 930122
for δ = 0.1 (50 / 50 tasks completed, 15.2 core-hours) and job 930172 for δ = 0.05, whose stalled and slow shards
were moved twice (jobs 930224 and 930266; 35.4 recorded core-hours, about 15 of them productive). 1,000 units per
strength (the zero-effect anchor and the translation control at e = 1.0 × 500 replicates), three analyses per unit,
6,000 analyses in all, 0 arm failures, every `latent_stage` arm available.

**Decision:** `diagnose_configuration_excess` at δ = 0.1 and `primary_validated` at δ = 0.05, hence
**`diagnose_configuration_excess`** under the predeclared conjunction (`case_study_type_i_overall_decision.json`).
The `adjusted` arm — the case study's primary — holds the bound for every statistic at both controls at δ = 0.05 and
for five of six at δ = 0.1, exceeding it on the translation control's `delta` by two rejections (36 / 500 = 0.072
against 0.0695). The `latent_stage` arm — the case study's sensitivity arm and the predeclared fallback — exceeds the
bound at the anchor on `angle` (0.090 / 0.084) and `shape` (0.092 / 0.078) at both strengths, so it is not
admissible either.

**In one paragraph.** At the Phase 6 SEA-AD design with the per-omic stage effect weakened to δ = 0.1 and 0.05 and
a numeric covariate `z` confounded with group, stage and their interaction injected at `s_z = 10` — calibrated so
that the covariate-adjusted analysis reads a held-out CV stage AUROC of 0.70 / 0.69 (the cohorts read 0.64–0.70) and
the covariate-blind basis spends its first component on `z` (median |Spearman| 0.92; the real basis reads 0.84) —
each of 1,000 units per strength was analysed three times on the same data, CV folds and permutation draws at the
rank floor of 2 with the net-displacement estimator: with `z` and five pure-noise covariates in the trajectory model
(`adjusted`), with no covariates (`unadjusted`) and with the six covariates under the `stage` latent adjustment
(`latent_stage`). The unadjusted arm confirms that the construction confounds (anchor `delta` 0.672 / 0.686, `shape`
0.252 / 0.222). The adjusted model removes almost all of it: anchor `delta` 0.068 / 0.064, `angle` 0.048 / 0.048,
`shape` 0.056 / 0.054, translation `delta` 0.072 / 0.060, `angle` 0.042 / 0.044, `shape` 0.068 / 0.068. Its `delta`
sits at or just above the bound in every reading (0.066 pooled over the four control readings, 34–36 rejections
where 34 are allowed) and its `shape` at the bound at translation; `angle` is comfortably under it but its RRPP null
is wide (median q95 109–113°). The `stage` latent adjustment removes the nuisance axis from the basis (|ρ| 0.08) at
the price of a heavier rank tail (up to 11), a lower AUROC (0.66 / 0.62) and anchor `angle` and `shape` rates of
0.078–0.092. Under the rule written before the run, the case study's primary is **withheld pending diagnosis** of a
small `delta` inflation at six design columns and n = 80, and its `stage` sensitivity arm's `angle` and `shape`
rejections are not interpretable at this configuration.

## 1. Construction and calibration

The check composes two validated constructions by import (`scripts/case_study_type_i_check.py`, design D1):

- **The D10 covariate construction** (`scripts/covariate_type_i_check.py`, [covariate Type I
  report](covariate-type-i-2026-10-16.md)): a null semi-synthetic dataset at the Phase 6 small-n design (three
  stages, `group_stage_sizes = [[11, 10, 28], [9, 10, 12]]`, n = 80, methylation + expression, joint magnitude),
  a numeric covariate `z = 0.5·g + 1·s + 2·g·s + ε` and its effect injected into every feature along a fixed
  direction over 20 % of the joint features, `X_f += s_z·sd_f·w_f·z`, on the M-value scale. The three-level
  categorical covariate D10 draws stays in the RNG stream but is neither declared nor injected (`s_c = 0`): the
  case study enters six numerics and no factor, and without a factor no CV training set can lose a level (design D2).
- **The weakened stage signal of the rank-floor check** ([report](rank-floor-2026-10-23.md)): one `--delta` applied
  to `delta_methyl`, `delta_expr` and `delta_protein`, with `p_dmp` at the Phase 6 value 0.1 so the translation
  surgery keeps its pool headroom.

Beside `z` the check declares five pure-noise covariates `x1`–`x5` (i.i.d. N(0, 1) from their own seed, never
injected), so the trajectory model carries six numeric design columns, as the case study's does. Cells are the
zero-effect anchor and the `translation` negative control at e = 1.0 (design D4). Seeds depend on the replicate
index only, so both cells and all three arms of a replicate share the baseline draw, the covariate values and the
RRPP seed; the `adjusted` and `unadjusted` latent spaces are asserted identical per unit (they share the
covariate-blind basis), so an `adjusted` excess is attributable to the trajectory model and an `unadjusted` excess
to confounding. All three arms run at `min_components=2`, the harness CV knobs (3 × 4 folds, 5 repeats, `max_components`
20, `random_state` 1203), 999 permutations and the net-displacement estimator.

### Calibration

Three criteria were fixed before any run (design D3): the `adjusted` arm's mean held-out CV AUROC inside
`[0.60, 0.75]` (the SEA-AD cohorts read 0.636–0.697); the `unadjusted` arm's median first-component |Spearman ρ|
with `z` at least 0.8 (the real basis reads 0.84); and the `unadjusted` arm rejecting above α on at least one
statistic at the anchor, so the gate is not vacuous. The knobs were δ and `s_z`. Six iterations
(`calibration/calibration.json`; R = 20, 99 permutations, base seed 2026):

| # | δ | `s_z` | adjusted mean AUROC | unadjusted median \|ρ\| | unadjusted anchor rates (`delta` / `angle` / `shape`) | met |
|---|---|---|---|---|---|---|
| 1 | 0.10 | 3 | 0.691 | 0.565 | 0.55 / 0.05 / 0.30 | no (\|ρ\|) |
| 2 | 0.15 | 3 | 0.744 | 0.571 | 0.55 / 0.05 / 0.35 | no (\|ρ\|) |
| 3 | 0.10 | 10 | 0.699 | 0.876 | 0.55 / 0.10 / 0.25 | **yes** |
| 4 | 0.10 | 20 | 0.699 | 0.918 | 0.65 / 0.10 / 0.30 | yes |
| 5 | 0.05 | 10 | 0.692 | 0.873 | 0.60 / 0.05 / 0.25 | **yes** |
| 6 | 0.05 | 20 | 0.693 | 0.917 | 0.65 / 0.10 / 0.25 | yes |

At D10's `s_z = 3` the basis does not follow `z` (median |ρ| 0.57, which is the regime of the 2026-10-16 check). At
`s_z ≥ 10` it does, and the realized stage AUROC is then set by the injected covariate rather than by δ — 0.69–0.70
at both δ = 0.05 and 0.10 — because `z` is confounded with stage. The criteria did not conflict. The settled
strengths are `s_z = 10` (the smallest injection meeting the threshold, so the nuisance is the case study's rather
than a harsher one) at δ = 0.10 and δ = 0.05, the rank-floor check's pair; the user approved two paper-grade runs
with the decision as their conjunction.

Realized at paper grade (anchor, R = 500): adjusted mean AUROC 0.702 (q05 0.641, q95 0.771) at δ = 0.1 and 0.688
(0.628, 0.756) at δ = 0.05; unadjusted median |ρ| 0.917 (q05 0.697, q95 0.985) and 0.918 (0.693, 0.985); unadjusted
anchor rates `delta` 0.672 / 0.686, `angle` 0.104 / 0.098, `shape` 0.252 / 0.222. All three criteria hold in both
runs (`calibration.met: true` in each decision file).

### Smoke profiles

`smoke-delta0.1-sz10/` and `smoke-delta0.05-sz10/` (R = 20, 99 permutations, both cells, three arms): 40 / 40 units
at each strength, 0 arm failures, the `latent_stage` arm available in every unit. Both smoke decisions were
`primary_validated`; the `latent_stage` arm's anchor `shape` rate of 0.35 (7 / 20) at both strengths was recorded in
the design as the warning sign to carry into the paper-grade run, and the paper-grade run confirmed an excess on
`shape` (0.092 / 0.078) and on `angle` (0.090 / 0.084) at R = 500. One unit took 57 s locally at 999 permutations
and 47–55 s on an `ing` EPYC core.

## 2. Rates against the bound

Rejection rates (rejections / 500) per control × arm × statistic, with the Monte Carlo SE `√(r(1 − r)/500)`, read
from `case_study_type_i_rates.csv` of each run. The judged arms carry a verdict against 0.0695; the `unadjusted`
arm is descriptive.

| control | arm | statistic | δ = 0.1 rate (rejections) | SE | verdict | δ = 0.05 rate (rejections) | SE | verdict |
|---|---|---|---|---|---|---|---|---|
| anchor | adjusted | `delta` | 0.068 (34) | 0.011 | holds | 0.064 (32) | 0.011 | holds |
| anchor | adjusted | `angle` | 0.048 (24) | 0.010 | holds | 0.048 (24) | 0.010 | holds |
| anchor | adjusted | `shape` | 0.056 (28) | 0.010 | holds | 0.054 (27) | 0.010 | holds |
| anchor | latent_stage | `delta` | 0.060 (30) | 0.011 | holds | 0.064 (32) | 0.011 | holds |
| anchor | latent_stage | `angle` | **0.090** (45) | 0.013 | **exceeds** | **0.084** (42) | 0.012 | **exceeds** |
| anchor | latent_stage | `shape` | **0.092** (46) | 0.013 | **exceeds** | **0.078** (39) | 0.012 | **exceeds** |
| anchor | unadjusted | `delta` | 0.672 (336) | 0.021 | descriptive | 0.686 (343) | 0.021 | descriptive |
| anchor | unadjusted | `angle` | 0.104 (52) | 0.014 | descriptive | 0.098 (49) | 0.013 | descriptive |
| anchor | unadjusted | `shape` | 0.252 (126) | 0.019 | descriptive | 0.222 (111) | 0.019 | descriptive |
| translation | adjusted | `delta` | **0.072** (36) | 0.012 | **exceeds** | 0.060 (30) | 0.011 | holds |
| translation | adjusted | `angle` | 0.042 (21) | 0.009 | holds | 0.044 (22) | 0.009 | holds |
| translation | adjusted | `shape` | 0.068 (34) | 0.011 | holds | 0.068 (34) | 0.011 | holds |
| translation | latent_stage | `delta` | 0.048 (24) | 0.010 | holds | 0.052 (26) | 0.010 | holds |
| translation | latent_stage | `angle` | 0.064 (32) | 0.011 | holds | 0.064 (32) | 0.011 | holds |
| translation | latent_stage | `shape` | 0.056 (28) | 0.010 | holds | 0.056 (28) | 0.010 | holds |
| translation | unadjusted | `delta` | 0.718 (359) | 0.020 | descriptive | 0.728 (364) | 0.020 | descriptive |
| translation | unadjusted | `angle` | 0.102 (51) | 0.014 | descriptive | 0.106 (53) | 0.014 | descriptive |
| translation | unadjusted | `shape` | 0.260 (130) | 0.020 | descriptive | 0.242 (121) | 0.019 | descriptive |

Verdicts are read verbatim from the decision files: at δ = 0.1, `"exceeding"` lists `latent_stage` anchor `angle`
and `shape` and `adjusted` translation `delta`, `"admissible": []`, `"decision": "diagnose_configuration_excess"`;
at δ = 0.05, `"exceeding"` lists `latent_stage` anchor `angle` and `shape`, `"admissible": ["adjusted"]`,
`"decision": "primary_validated"`. No verdict is `unavailable` (0 failures). Pooled descriptively over the four
control readings (two controls × two strengths, 2,000 units), the adjusted arm rejects `delta` at 0.066 (132 / 2,000),
`shape` at 0.062 (123 / 2,000) and `angle` at 0.046 (91 / 2,000); the `latent_stage` arm rejects at the anchor
`angle` 0.087 (87 / 1,000) and `shape` 0.085 (85 / 1,000).

## 3. What the arms show

**The construction confounds, and the trajectory model removes almost all of it.** The unadjusted arm rejects
`delta` in two thirds of the null units and `shape` in a quarter; the adjusted arm, on the identical latent scores
(asserted per unit), brings `delta` to 0.060–0.072, `shape` to 0.054–0.068 and `angle` to 0.042–0.048. The
covariate Type I check of 2026-10-16 measured the same model at three design columns (`z` plus the two indicator
columns of a three-level factor), a saturated stage signal (AUROC 1.000), the unfloored rank rule and a nuisance the basis did not follow
(|ρ| 0.64): adjusted `delta` 0.040, `angle` 0.062, `shape` 0.034. Here, at six columns, the case study's AUROC and
a basis that follows the nuisance, `delta` reads 0.064–0.072 in every control reading, at or above the 34-rejection
limit, and `shape` reads 0.068 at translation in both runs; `angle` is lower than in 2026-10-16. The excess that
trips the rule — two rejections on one control at one strength — is one Monte Carlo SE wide, but the pooled `delta`
of 0.066 over 2,000 units, every reading above α, is what the diagnosis should take as its starting point: a small
`delta` inflation of the covariate-adjusted model at six design columns and n = 80, which the three-column check did
not show.

**The adjusted arm's `angle` null is wide.** The median q95 of the RRPP `angle` null is 109° (δ = 0.1) and 113°
(δ = 0.05) in the adjusted arm against 36° in the unadjusted arm on the same basis and 57–59° in the `latent_stage`
arm (`case_study_type_i_replicates.csv`, `angle_null_q95`). With six covariates in the reduced model the permuted
residuals carry little stage structure, so the permuted net displacements point almost anywhere; `angle` is
size-controlled at this configuration but, as the design-guidance δ arm found at δ ≤ 0.1, it has no power to spend.

**The `stage` latent adjustment is inadmissible at this configuration.** It does what it is built to do — the
first component's |ρ| with `z` falls from 0.92 to 0.08 — and its translation control holds (0.048–0.064), but at the
anchor it rejects `angle` 0.090 / 0.084 and `shape` 0.092 / 0.078, about three Monte Carlo SE above α in each of the
four readings. Beside it, the basis changes: the held-out AUROC falls to 0.659 / 0.622 (from 0.702 / 0.688), the CV
rank's tail lengthens (rank ≥ 4 in 91 / 500 and 61 / 500 anchor units against 6 / 500 and 1 / 500 in the adjusted
arm; maximum 11), and the smallest training-set residual df inside the CV is 31 (`cv_min_training_df`; full-sample
`residual_df` 71: n − 1 − 2 stage indicators − 6 covariates). The 2026-10-20 latent-adjustment check measured the
same conditioning at the unfloored rule, eight covariates including a factor and a saturated stage signal (δ = 2.0):
`stage` held `delta` 0.031, `angle` 0.067, `shape` 0.058 against 0.0706, while its covariate-blind baseline arm
rejected `shape` at 0.078. The excess here therefore belongs to the weak-signal, floor-2, six-numeric regime, not to
the conditioning as such; the per-training-set residualization at 31 residual df and the rank tail are the two
things to look at first.

**Rank and signal.** The adjusted arm selects rank 2 in 366 / 500 (δ = 0.1) and 407 / 500 (δ = 0.05) anchor units
and rank 3 in most of the rest; no unit selected rank 1 (the floor binds). The rank-floor check at the same
strengths without covariates measured the floor-2 arm's anchor at `delta` 0.044 / 0.042, `angle` 0.060 / 0.054 and
`shape` 0.038 / 0.036 (δ = 0.1 / 0.05), so the covariate-adjusted model, not the floor, carries the `delta` and
`shape` readings here.

## 4. Decision, consequence and what to diagnose

**Decision**, read verbatim from `results/case-study-type-i-2026-10-29/case_study_type_i_overall_decision.json`:
`"decision": "diagnose_configuration_excess"`, `"admissible": []`, runs `delta0.1 → diagnose_configuration_excess`
and `delta0.05 → primary_validated`. The rule, written before the run (proposal, design D5): `primary_validated` when
the `adjusted` arm holds every statistic at both controls in every run; otherwise `adopt_latent_stage` when the
`latent_stage` arm does; otherwise `diagnose_configuration_excess`. The `adjusted` arm fails the conjunction on one
reading (translation `delta` at δ = 0.1, 36 / 500) and the `latent_stage` arm fails on four (anchor `angle` and
`shape` at both strengths).

**Consequence for the case study.** Under `diagnose_configuration_excess` the primary analysis of `sead_sexdiff`
(covariates in the trajectory model, floor 2, net displacement, 999 permutations) is **withheld pending diagnosis**:
it does not rest on a validated size at its configuration, and it cannot move to the `stage` latent adjustment,
which is not admissible here either. Concretely for the readers of its results, until the diagnosis is done: the
primary's `delta` p-values should be read against a realized size of about 0.066 rather than 0.05 (a p-value near
0.05 is not a rejection at the nominal level; a p-value well below it is unaffected in kind); its `angle` and `shape`
p-values are size-controlled in this check (0.042–0.068), with `angle`'s null so wide that a non-rejection is
uninformative; and the `stage` sensitivity arm's `angle` and `shape` rejections are not interpretable at this
configuration. Nothing in `sead_sexdiff` changes under this report; the diagnosis is its own change there and here.

**What to diagnose** (named by the rule; none of it is done here):

1. The `adjusted` arm's `delta` (0.060–0.072, pooled 0.066). Candidates, in the order they should be excluded: the
   degrees-of-freedom cost of six columns at n = 80 (the same seeds and construction at one, three and six noise
   columns, which isolates the column count from the confounding); the LS-mean reference row and the RRPP reduced
   model under covariates at these cell sizes (the design's `reference_row` is the pooled mean; the reduced model
   keeps the covariates); and the bound itself at R = 500 (a replication at R = 1,000 on the same seed family would
   halve the SE, at about 30 core-hours per strength).
2. The `latent_stage` arm's anchor `angle` and `shape` (0.078–0.092). Candidates: the per-training-set nuisance
   coefficients at `cv_min_training_df` 31 (the residualized CV rank's tail reaches 11, and `angle`/`shape` at rank
   ≥ 4 should be tabulated against rank 2–3); and the interaction of the residualized basis with the rank floor,
   which the 2026-10-20 check never measured (it ran at floor 1).

## 5. Limitations and scope

Four limitations the spec requires the report to carry, then the scope:

- **Astrocyte primary cohort sizes only.** The check uses the astrocyte `cohort_30` cell sizes
  `[[11, 10, 28], [9, 10, 12]]` (n = 80) that every Phase 6 check uses; the microglia cohort (n = 76) and the other
  cohorts' tables were not run, and size at their cell sizes is inferred, not measured.
- **Additive-linear covariate effects.** `z` enters the features additively along one fixed direction and the
  trajectory model linearly; non-linear or interaction effects of real covariates on the features are outside the
  construction.
- **One confounded covariate beside five independent ones.** The confounding rests on `z` alone; the five noise
  covariates are i.i.d. normals, uncorrelated with one another and with the design, unlike age, education, PMI,
  RIN and library sizes. The configuration matched is the column count, the floor and the nuisance regime, not the
  covariates' joint distribution (design D2 records the alternatives considered).
- **Structural-analogue blocks.** The blocks are InterSIM-calibrated methylation and expression analogues (498
  measured features), not the case study's TMM log2-CPM blocks of 5,000 features each; the per-omic δ is a
  generator strength, and the regime is matched on the realized AUROC and |ρ|, not on the feature space.

Scope:

- No power is claimed at this configuration (the design-guidance pilot's δ arm measured none at δ ≤ 0.1), and the
  `group_stage` conditioning, a floor-1 arm and the microglia sizes were not run.
- No change under `src/motco/`: every record field the check reads exists in 0.12.0, so the release and the case
  study's pin are untouched (design D9).
- The two strengths are independent runs with their own decision files; the conjunction rule and the number of
  runs were fixed with the user after calibration and before submission (design D3 record). The δ = 0.05 run's
  shards were moved twice between nodes (`PROVENANCE.txt`); shards are resumable and a unit recorded in any shard
  file is never re-run, so the move changes which node computed a unit and nothing else.
- The replicate-level tables (`case_study_type_i_replicates.csv`, `case_study_type_i_pvalues.csv`) are committed
  for both runs; `merged.jsonl` and `records/` stay on `ing` under the worktree named in each `PROVENANCE.txt`, and
  `case_study_type_i_check.py merge` regenerates every table from `records/`.
