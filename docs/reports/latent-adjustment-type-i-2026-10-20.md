# Latent-adjustment Type I check — 2026-10-20

> Paper-grade size check of covariate-aware latent construction (`covariate_latent_adjustment`, `motco analyze
> --latent-adjustment`) at the Phase 6 SEA-AD design (`openspec/changes/covariate-aware-latent-construction/`, design
> D8, spec `latent-adjustment-validation`). It is the gate the case study must pass before it reports
> latent-adjusted trajectory results. Every number below is read from `results/latent-adjustment-type-i-2026-10-20/`
> (`latent_adjustment_type_i_decision.json`, `latent_adjustment_type_i_rates.csv`,
> `latent_adjustment_type_i_pvalues.csv`, `latent_adjustment_type_i_replicates.csv`, `merged.jsonl`). Rates are
> rejection rates at p ≤ 0.05 over the available replicates of 500, at 999 permutations.

**Run:** 500 replicates on `ing` (SLURM arrays 923443 and 923539, 50 shards, 5.1 recorded core-hours), from
`979657c` plus this change's working tree; file hashes are in `PROVENANCE.txt`. **Decision:** `adopt_stage`. The
admissible conditionings are `stage` only.

**In one paragraph.** Null datasets at the SEA-AD MTG design carried two covariates confounded with group × stage,
whose effect was injected strongly enough that the covariate-blind basis spends its first latent component on the
numeric covariate: median |Spearman| 0.83 between that component and the covariate. Each replicate was analysed with
eight covariates declared, in three arms that differ only in how the basis is built. With `stage` conditioning, the
three statistics hold the predeclared bound: `delta` 0.031, `angle` 0.067 and `shape` 0.058 against 0.0706
(R = 448). With `group_stage` conditioning, `angle` exceeds it: 0.082 against 0.0703 (R = 462). Residualization
moved the first component off the covariate (median |ρ| 0.34). The case study may use `--latent-adjustment stage`
under the conditions in §5, and not `group_stage`.

## 1. Construction

- **Data and covariates.** The covariate-adjustment check's construction
  ([report](covariate-type-i-2026-10-16.md) §1), imported from `scripts/covariate_type_i_check.py` rather than
  copied. The data are `generate_semisynthetic_trajectory` at the `phase6_small_n` settings: three stages,
  `group_stage_sizes = [[11, 10, 28], [9, 10, 12]]`, methylation (M-values) and expression blocks, 498 features,
  `group_effect_size = 0`. A numeric `z` and a three-level factor depend on group, stage and their interaction
  (`a_g = 0.5`, `a_s = 1`, `a_gs = 2`, logit scale 0.5).
  - The seeds are the covariate check's (base seed 2026). Replicate `r` therefore draws the same generated blocks,
    `z` and factor as that check's replicate `r`.
- **Noise covariates.** Six i.i.d. N(0, 1) numeric covariates `x1`–`x6`, declared but never injected, are drawn
  from their own seed. Every arm therefore declares 8 covariates (9 design columns). This puts the residualization
  at the case study's 8–10 covariates and its inner CV folds at the case study's degrees of freedom (design, Risks).
- **Injection.** `X_f += s_z·sd_f·w_f·z + s_c·sd_f·Σ_k u_{k,f}·1[cat = k]` on a random 20 % of features per
  direction. The strengths are `s_z = 40` (calibrated, below) and `s_c = 2`.
- **Arms.** Each replicate runs `analyze_trajectories` three times with the same covariates, CV settings
  (harness defaults) and RRPP seed, and the net-displacement estimator:
  - `baseline`: covariates in the trajectory model only, with the covariate-blind basis (the 0.8.0 analysis);
  - `stage`: `covariate_latent_adjustment="stage"`;
  - `group_stage`: `covariate_latent_adjustment="group_stage"`.
- **Rule.** For each residualized arm and statistic, the rate over the R available replicates must be at most
  α + k·√(α(1 − α)/R), with α = 0.05 and k = 2. An arm is admissible when all three statistics hold. The decision
  is `adopt_stage` when `stage` is admissible, `adopt_group_stage` when only `group_stage` is, and otherwise
  `diagnose_latent_excess`. The baseline arm is reported only.

### Calibration

The gate needs the covariate-blind basis to be spent on the nuisance axis. The criterion is a median |Spearman| of
at least 0.8 between the first component's scores and `z`, checked on the smoke profile before any paper-grade run.
At the covariate check's `s_z = 3`, the median over the 20 smoke replicates was 0.64. It rises with `s_z` and then
saturates: 0.74 at 6, 0.78 at 10, 0.80 at 20, 0.81 at 30 and 0.815 at 40. The covariate-blind fit re-standardizes
every feature, so once the injected features are dominated by `z` a stronger injection adds nothing; only the
injected fifth of the features can carry `z`. `s_z = 40` was chosen for the most margin. The smoke profile (R = 20,
99 permutations) confirmed 0.815.

In the paper-grade run the baseline median |ρ| is **0.830** (5–95 % range 0.751–0.891; minimum 0.708), so the
criterion is met. The covariate-blind basis selected rank 2 in 265 replicates, rank 3 in 233 and rank 4 in 2.

## 2. Rates against the bound

From `latent_adjustment_type_i_rates.csv`. Rate (Monte Carlo SE).

| statistic | `stage` (R = 448, bound 0.0706) | `group_stage` (R = 462, bound 0.0703) | `baseline` (R = 500, descriptive) |
|---|---|---|---|
| `delta` | 0.031 (0.0082) holds | 0.048 (0.0099) holds | 0.048 (0.0096) |
| `angle` | 0.067 (0.0118) holds | **0.082 (0.0128) exceeds** | 0.062 (0.0108) |
| `shape` | 0.058 (0.0110) holds | 0.069 (0.0118) holds | 0.078 (0.0120) |

At p < 0.05 instead of p ≤ 0.05, every rate is the same except `group_stage` `angle`, which is 0.078 and still
exceeds. The decision is unchanged.

**Stability.** In five consecutive blocks of 100 replicates (rates over each block's available replicates):

- `stage`: `delta` 0.04, 0.01, 0.03, 0.03, 0.03; `angle` 0.03, 0.06, 0.07, 0.13, 0.04; `shape` 0.05, 0.07, 0.06,
  0.04, 0.07.
- `group_stage`: `delta` 0.04, 0.02, 0.09, 0.04, 0.04; `angle` 0.05, 0.10, 0.05, 0.14, 0.06; `shape` 0.07, 0.07,
  0.08, 0.07, 0.06.

In both arms the `angle` rate is driven by block 4, replicates 300–399. The covariate check's adjusted arm, on the
same generated datasets, also had its `angle` excess in block 4 (0.11). The `stage` arm's 0.067 sits 0.004 below its
bound; `group_stage`'s 0.082 sits 0.012 above, and its block 2 is also high.

**Unavailable replicates.** An arm is unavailable for a replicate when a CV training design loses full rank.

- That happened in 52 `stage` replicates and 38 `group_stage` replicates (22 in both), never in the baseline arm.
- Every one is an inner-fold training set short exactly one column, an absent level of the three-level factor.
  They are spread across the five CV repeats (repeat 1 most often).
- Excluding them does not flatter the residualized arms. On the replicates `stage` lost, the always-available
  baseline arm rejects `delta` 0.019, `angle` 0.019 and `shape` 0.096 (n = 52). On the replicates `stage` kept, it
  rejects `delta` 0.051, `angle` 0.067 and `shape` 0.076 (n = 448), close to its overall 0.048 / 0.062 / 0.078.
- The `group_stage` comparison is similar: 0.026 / 0.026 / 0.132 on n = 38, against 0.050 / 0.065 / 0.074.

## 3. What residualization changes

- **The basis.** The first component's median |ρ| with `z` falls from 0.830 (baseline) to 0.336 (`stage`) and
  0.355 (`group_stage`). The residualized bases no longer spend their first component on the nuisance axis; what
  correlation remains is the stage structure `z` shares by construction.
- **The rank.** CV selects rank 3 in 417 of 448 `stage` replicates and 457 of 462 `group_stage` replicates, against
  233 of 500 for the baseline. With the nuisance axis gone, the CV keeps more of the stage structure.
- **Degrees of freedom.** The full-sample residual df is 68 under `stage` (80 − 1 − 2 − 9) and 65 under
  `group_stage`. The smallest CV training-set df is 28 and 25 in every replicate.
- **What is removed.** The per-feature partial R² of the covariates given the conditioning design has a
  replicate-median of 0.13 (median feature), 0.27 (mean) and 0.99 (95th percentile). About 18 % of features have
  R² above 0.5, roughly the fifth of features that carry `z`.
- **The baseline arm.** Here the 0.8.0 covariate adjustment alone, on the covariate-blind basis, lets `shape` reach
  0.078. That is above the 0.0695 bound it would face at R = 500, and its rule does not judge it. The covariate
  check's adjusted arm held `shape` at 0.034 with `s_z = 3` and two covariates. With a nuisance strong enough to
  take the first component, the model adjustment alone is no longer enough for `shape`, and latent adjustment is.
  This is descriptive, not a verdict.

## 4. Scope and limits

- **Mediators.** Residualizing a variable that mediates the stage effect removes real signal. No simulation can
  tell a confounder from a mediator, because the data look the same; this check injects a pure confounder. Choosing
  and justifying the covariate set is the case study's job, in its own proposal. A technical variable that is itself
  downstream of disease stage, such as a cell-quality metric that worsens with pathology, is a mediator candidate.
- **Additive, linear effects.** The covariate effects are additive and linear, and correctly specified in the
  residualization and the model. Non-linear effects and covariate × design interactions are outside the check.
- **One design.** n = 80, two blocks, 498 features, 8 covariates with 9 design columns, and one injected direction
  per covariate.
- **`angle` is near the bound.** `stage` holds `angle` by 0.004, with one 100-replicate block at 0.13. Treat a
  marginal latent-adjusted `angle` rejection with the care the earlier reports ask of `angle` near the translation
  control.
- **Rare factor levels.** About 10 % of replicates lose an arm to a factor level absent from a CV training set, and
  which fold fails depends on `random_state`. Merge rare levels before analysis; the record names the failing repeat
  and fold.
- **No power claim.** This report covers size only. Power, with and without a nuisance axis, is the subject of the
  separate [power comparison](latent-adjustment-power-2026-10-20.md).

## 5. Consequence for the case study

The decision is `adopt_stage`. The SEA-AD case study may report latent-adjusted `motco analyze` results with
`--latent-adjustment stage` (MOTCO ≥ 0.9.0) under these conditions:

- its own proposal names and justifies the covariate set, including why no residualized variable mediates the
  stage effect;
- `group_stage` conditioning is not used: it exceeded the bound on `angle`;
- the record's `latent_adjustment` block (conditioning, design columns, residual df, smallest CV training df, the
  feature R² summary) and its `covariates` block are reported beside the results;
- categorical covariates have no level rare enough to drop out of a CV training set;
- the earlier binding statements stand: the CV-selected rank is reported beside `angle`, and non-rejections are
  uninformative;
- the [power comparison](latent-adjustment-power-2026-10-20.md) §6 applies: latent adjustment is a sensitivity
  analysis beside the covariate-adjusted analysis, not a replacement for it.
