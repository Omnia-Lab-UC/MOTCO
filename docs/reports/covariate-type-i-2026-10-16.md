# Covariate-adjustment Type I check — 2026-10-16

> Paper-grade size check of covariate-adjusted `motco analyze` at the Phase 6 SEA-AD design
> (`openspec/changes/covariate-adjusted-trajectories/`, spec `covariate-adjustment-type-i`). It is the gate the case
> study must pass before it reports covariate-adjusted trajectory results. Every number below is read from
> `results/covariate-type-i-2026-10-16/` (`covariate_type_i_decision.json`, `covariate_type_i_rates.csv`,
> `covariate_type_i_pvalues.csv`, `merged.jsonl`). Rates are rejection rates at p ≤ 0.05 over 500 replicates × 999
> permutations (Monte Carlo SE about 0.0097 at 0.05).

**Run:** 180 replicates on the workstation, then the remaining 320 on `ing` as SLURM array 921167 (60 shards, all
`COMPLETED`). The code is `c3ee13a` plus this change's working tree, and the script and `analysis.py` hashes in
`PROVENANCE.txt` were verified identical on both machines. **Decision:** `covariates_validated`.

**In one paragraph.** On null datasets at the SEA-AD MTG design, a numeric and a categorical covariate were made to
depend on group, stage and their interaction, and their effect was injected into the features. Without adjustment
the analysis rejects too often: `shape` at 0.364, `delta` at 0.098 and `angle` at 0.096. With both covariates
declared, all three statistics hold the predeclared bound 0.05 + 2·√(0.05·0.95/500) = 0.0695: `delta` 0.040,
`angle` 0.062 and `shape` 0.034. The two arms share every latent space (asserted per replicate), so the difference
comes entirely from the trajectory model. The case study may report covariate-adjusted results under the conditions
in §5.

## 1. Construction

- **Data.** `generate_semisynthetic_trajectory` at the `phase6_small_n` settings: three stages,
  `group_stage_sizes = [[11, 10, 28], [9, 10, 12]]`, `p_dmp = 0.1`, ρ = 0, joint magnitude, `group_effect_size = 0`
  (null). Only the methylation (367 CpGs, converted to M-values) and expression (131 genes) blocks are measured;
  both are analysed with the `identity` transform. Base seed 2026. Generator, covariate and RRPP seeds are derived
  per replicate.
- **Covariates.** With `g` the group indicator and `s` the stage index (0–2):
  - numeric `z = a_g·g + a_s·s + a_gs·g·s + ε`, with ε ~ N(0, 1);
  - categorical: three levels from a multinomial logit with level logits `k·c·η` (k = 0, 1, 2) on the same
    linear predictor η.
  - Values: `a_g = 0.5`, `a_s = 1`, `a_gs = 2` and `c = 0.5`.
  - Across replicates, the numeric covariate's group × stage interaction contrast
    `(z̄_B,2 − z̄_B,0) − (z̄_A,2 − z̄_A,0)` has mean 3.98 (SD 0.58). The covariates therefore bend the groups'
    trajectories differently, not only shift them.
- **Injection.** `X_f += s_z·sd_f·w_f·z + s_c·sd_f·Σ_k u_{k,f}·1[cat = k]`.
  - `sd_f` is the feature's pooled SD on the analysed scale.
  - `w` and the `u_k` are unit directions, each on a random 20 % of the 498 joint features, drawn once from the base
    seed.
  - Strengths: `s_z = 3` and `s_c = 2`.
- **Arms.** Each replicate is analysed twice with `analyze_trajectories`, using its defaults: net displacement, CV
  rank with the harness knobs, and the same RRPP seed.
  - `adjusted`: `covariate_cols = ["z", "category"]`, `category` categorical. The design columns are `z`,
    `category[c1]` and `category[c2]`.
  - `unadjusted`: no covariates.
  - The latent scores and selected rank were identical in both arms for all 500 replicates. CV selected rank 2 in
    264 replicates and rank 3 in 236.
- **Rule.** Per statistic, the adjusted rate over the R available replicates must be ≤ α + k·√(α(1 − α)/R), with
  α = 0.05 and k = 2. `covariates_validated` requires all three statistics to hold. The unadjusted arm is reported
  only and does not enter the decision.

### Calibration

The design's first strengths (`a_s = 0.5`, `a_gs = 0.75`, `s_z = 0.5`, `s_c = 0.3`) did not confound: in the
R = 20 smoke run, the unadjusted arm rejected at about α, which would have made the gate vacuous.

Raising `s_z` alone, even to 30, did not help. The injection carries `z`'s within-cell noise along `w` in the same
proportion as its cell-mean pattern, so the RRPP null widens together with the signal. The confounding
signal-to-noise is `a_gs` relative to the covariate's noise SD.

The calibrated values above were chosen on probe replicates and confirmed by the smoke profile (R = 20, 99
permutations): unadjusted `shape` 0.45. The categorical logit scale keeps the factor mixed within cells, so
complete-cell alignment cannot trip the full-rank guard. All values are recorded in the decision file.

## 2. Rates against the bound

From `covariate_type_i_rates.csv`. Rate (MC SE); R = 500 for every statistic, no failed replicate; bound 0.0695.

| statistic | adjusted | verdict | unadjusted (descriptive) |
|---|---|---|---|
| `delta` | 0.040 (0.0088) | holds | 0.098 (0.0133) |
| `angle` | 0.062 (0.0108) | holds | 0.096 (0.0132) |
| `shape` | 0.034 (0.0081) | holds | 0.364 (0.0215) |

**Stability.** Split into five consecutive blocks of 100 replicates, the adjusted rates are:

- `delta`: 0.07, 0.04, 0.05, 0.00, 0.04;
- `angle`: 0.05, 0.05, 0.05, 0.11, 0.05;
- `shape`: 0.06, 0.02, 0.03, 0.02, 0.04.

The `angle` rate of 0.062 comes from one block at 0.11; the other four blocks sit at 0.05. At p < 0.05 instead of
p ≤ 0.05, the convention of the earlier Type I reports, the rates are 0.040, 0.060 and 0.034. The decision is the
same.

**Against the no-covariate anchor.** The parent Phase 6 study's `pls_cv × net_displacement` anchor (no covariates,
independent seeds; [net-displacement Type I report](net-displacement-type-i-2026-10-13.md)) had `delta` 0.026,
`angle` 0.050 and `shape` 0.020. The adjusted arm here is 0.01–0.015 higher on each statistic but within the bound.
That is plausible for a model with three more columns at n = 80, where RRPP is exact only asymptotically. It is
recorded here and not judged by the rule.

## 3. What the unadjusted arm shows

The construction confounds, mostly through `shape`: an interaction between the covariate and stage bends one
group's path relative to the other's. With 36 % false `shape` rejections, an unadjusted SEA-AD `shape` result under
comparable confounding could not be told apart from a donor-composition difference. `delta` and `angle` inflate
about two-fold. The adjustment removes all three inflations, which is the property the gate tests.

## 4. Scope and limits

- The check covers **additive, linear** covariate effects that are correctly specified in the model. A covariate
  whose true effect is non-linear, or interacts with group or stage, is outside it.
- The latent space is covariate-blind. In this construction the covariate effect is correlated with stage, so PLS
  may spend capacity on it; that case is measured here rather than assumed away. A covariate effect much larger
  than the stage signal is not.
- This is one design: n = 80, two blocks, 498 features, three covariate design columns. More covariate columns at
  the same n use more residual degrees of freedom; this check does not bound that cost.
- There is no power claim. Adjustment costs degrees of freedom and, in the presence of confounding, also removes
  confounded signal; the operating studies' power numbers were measured without covariates.

## 5. Consequence for the case study

The decision is `covariates_validated`. The SEA-AD case study may report covariate-adjusted `motco analyze` results
(MOTCO ≥ 0.8.0) under these conditions:

- keep the covariate set parsimonious (each factor level is a column);
- report the record's `covariates` block (design columns, missing policy and dropped donors) beside the results;
- carry the additive-linear and covariate-blind-basis limitations above;
- keep the earlier binding statements unchanged: the CV-selected rank is reported beside `angle`, and
  non-rejections are uninformative.
