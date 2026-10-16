# Design

## Context

`analyze_trajectories` (0.7.0) builds its trajectory design with the harness's `build_simulation_trajectory_design`. That function returns `model_full` (intercept, group, stage, group × stage), `model_reduced` (no interaction), the group-major LS-mean rows and the two-group contrast. It then passes these to `estimate_difference` and `RRPP`, which already accept arbitrary model matrices and LS-mean rows: betas come from the normal equations, LS-means are `ls_means @ betas`, and RRPP permutes `Y − X_r β_r`. Covariates therefore need no change to the stats primitives. They need extra columns on the three matrices, plus the input handling and record around them. The harness must stay untouched. `tests/test_analysis.py::test_matches_harness_*` and the committed byte-identity fixtures pin it. See proposal.md for the motivation and the specs for the contract.

## Goals / Non-Goals

**Goals:**

- Covariate adjustment confined to `motco.analysis`, plus one stateless design helper in `motco.stats.design`.
- With no covariates, output identical to 0.7.0, pinned by a fixture captured before any code change.
- A Type I gate at the Phase 6 design that the case study can cite.

**Non-Goals:**

- Covariates in the simulation harness, the study engine or `motco de`.
- Covariate × group or covariate × stage interactions, splines and other non-linear covariate terms.
- Covariate-aware latent construction (residualizing the features before PLS) and covariate-aware attribution.
- Caller-chosen LS-mean reference values. The statistics do not depend on them (see D3), so the pooled mean is the only reference.

## Decisions

### D1. The latent space stays covariate-blind; adjustment happens in the latent regression

PLS-DA is fitted exactly as in 0.7.0: on the standardized joint features, against stage. Covariates enter only the trajectory regression on the latent scores `T = X W`. OLS commutes with a fixed linear map: regressing `T` on `[D, C]` gives `B_T = B_X W`, where `B_X` comes from regressing `X` on `[D, C]`. Adjusting in the latent space is therefore exactly the feature-space adjustment, projected. Only the basis `W` is chosen without the covariates.

*Alternative considered:* residualizing `X` on `C`, conditional on group × stage, before PLS. This changes the latent space, needs its own power and Type I validation, and would break the guarantee that the latent space is shared with the harness. The user deferred it, and it is recorded as a limitation. A case-study sensitivity analysis can compare it against this approach. Keeping the latent space unchanged also lets the Type I check's two arms share one latent fit.

### D2. Design construction: append covariate columns to both models

A new helper `build_covariate_matrix(frame, columns, categorical) -> CovariateDesign` in `stats/design.py` returns:

- the `n × q` design matrix, numeric columns as-is and categorical columns as drop-first indicators with `_sort_levels` ordering (the same ordering helper `get_model_matrix` uses);
- the column names (`name` for numeric, `name[level]` for factor indicators);
- the per-covariate kind, levels and reference level;
- the reference row, the column means.

In `analysis.py`, after `build_simulation_trajectory_design` (unchanged), the design gets `model_full = [M_f | C]`, `model_reduced = [M_r | C]` and `ls_means = [L | 1 ⊗ c̄]`. The covariate columns are appended last. Group, stage and interaction positions are unchanged, so the contrast and the group-major row order hold as they are.

*Rank check:* `np.linalg.matrix_rank(model_full)` must equal its column count. Otherwise the analysis raises `TrajectoryAnalysisError`, naming the covariates. The check is explicit because `estimate_betas` falls back silently to `lstsq`, which would return a minimum-norm fit and LS-means that depend on the arbitrary null-space choice.

*Alternatives:* centering numeric covariates before they enter. Fits and statistics are invariant to it, so raw values are kept, which makes the recorded coefficients directly interpretable. A formula interface such as patsy was rejected as a new dependency for one additive term.

### D3. LS-means at the pooled covariate mean; statistics are invariant to that choice

The stage mean is `μ_gs = L_gs β + c̄ γ`, so the covariate term adds the same vector `c̄ γ` to every cell. `delta` is a path length over differences, `angle` uses differences (net displacement) or centered PC1, and `shape` is a Procrustes distance after centering. All three are translation-invariant, so the reference row changes only `stage_means.csv`. A test pins this by comparing two reference rows. The pooled mean of the indicator columns weights factor levels by their observed proportions, which is emmeans' "proportional" weighting. Equal weighting is not offered, because it would change nothing that is tested.

### D4. Covariate declarations live on `TrajectoryDataset`, not on the params

`covariate_cols: Sequence[str] = ()`, `categorical_covariates: Sequence[str] = ()` and `missing_covariates: Literal["error", "complete_case"] = "error"` are dataset fields. They describe the data, and validation already lives in `TrajectoryDataset.__post_init__`. The decisive reason is byte identity: `TrajectoryAnalysisParams.to_dict()` is written whole into `record["params"]`, so a new params field would add keys to every 0.7.0 record. The covariate block is emitted only when `covariate_cols` is non-empty.

### D5. Complete-case deletion happens at dataset construction

With `complete_case`, `__post_init__` computes the rows that have a missing value in any declared covariate. It replaces `metadata` and every block's `values` with the retained rows (frozen dataclass, via `object.__setattr__`, as the existing normalizations do) and stores the dropped ids in a private field. It does this *before* the existing metadata, block and design validators run, so cell emptiness is checked on the retained samples by the existing missing-cell check, extended to name the cells deletion emptied. Everything downstream (scalers, PLS CV, design, RRPP) then sees only retained samples, with no further special-casing. Under `error`, missing values are reported with the column and up to ten sample ids, matching the style of the existing errors.

*Alternative:* subsetting inside `analyze_trajectories`. This was rejected because the dataset's `sample_ids`, `groups` and cell sizes would then disagree with what was analysed.

### D6. Covariate kind inference

The kind is categorical if the caller declares it, or if the column's dtype is not numeric or is boolean (`pd.api.types.is_numeric_dtype` and not `is_bool_dtype`). Otherwise it is numeric. A numeric column with non-finite values (other than the NaN handled by the policy) is rejected. On the CLI, `pd.read_csv` makes string columns object dtype, so they become categorical without a flag. `--categorical-covariate` covers numeric-coded factors. Categorical values are compared as strings, like group and stage.

### D7. Record shape

```json
"covariates": {
  "columns": ["age_at_death", "education_level"],
  "kinds": {"age_at_death": "numeric", "education_level": "categorical"},
  "levels": {"education_level": ["College", "Graduate", "HighSchool"]},
  "reference_levels": {"education_level": "College"},
  "design_columns": ["age_at_death", "education_level[Graduate]", "education_level[HighSchool]"],
  "reference_row": {"age_at_death": 84.2, "education_level[Graduate]": 0.31, "education_level[HighSchool]": 0.22},
  "lsmean_reference": "pooled_mean",
  "missing_policy": "complete_case",
  "n_dropped": 2,
  "dropped_sample_ids": ["H20.33.004", "H21.33.017"],
  "coefficients": {"age_at_death": {"LV1": 0.012, "LV2": -0.004}, "...": {}}
}
```

This block is inserted after `design`. `design.n_samples` already reports the analysed count. `ANALYSIS_SCHEMA_VERSION` stays 1, because the key is optional and its absence is the 0.7.0 contract. The coefficients come from the observed full-model fit.

### D8. Byte-identity fixture captured before the change

The first task, before any code edit, runs the 0.7.0 `analyze_trajectories` on a small seeded fixture (two blocks, three stages, 99 permutations, attribution off) and commits `tests/data/analysis_no_covariate_identity.json`. The file holds the record without `runtime_seconds` and the `software` block, which varies by environment, plus SHA-256 digests of the three CSVs. A test compares the post-change output with it. This is the same pattern as `tests/data/orientation_default_identity.json`.

### D9. Attribution guard

`analyze_trajectories` raises before any fitting when `params.attribution.enabled` and `dataset.covariate_cols` are both set. Failing early avoids spending a CV fit on a request that cannot be honoured.

### D10. Type I check: construction and execution

`scripts/covariate_type_i_check.py`, a self-contained driver in the style of `net_displacement_follow_up.py`:

- **Data.** `generate_semisynthetic_trajectory` with the `phase6_small_n` generator settings, `group_effect_size = 0`, and the layers `methylation` + `expression` only. The methylation block is passed as M-values (logit applied by the script) with the `identity` transform, so the effect can be injected on one additive scale. The analysis treats the units as supplied either way.
- **Covariates.**
  - Numeric `z_i = a_g·g_i + a_s·s_i + a_gs·g_i·s_i + ε_i`, where `g ∈ {0, 1}`, `s` is the stage index and `ε ~ N(0, 1)`.
  - Categorical: three levels drawn from a multinomial logit on the same linear predictor.
  - The defaults `a_g = 0.5`, `a_s = 0.5` and `a_gs = 0.75` are recorded. The interaction term is what creates a spurious *trajectory* difference. Main effects alone would only translate or uniformly stretch the trajectories.
- **Injection.** For each measured block, `X_f += s_z·sd_f·w_f·z + s_c·sd_f·Σ_k u_{f,k}·1[cat = k]`. Here `sd_f` is the feature's pooled SD, `w` and `u` are unit random directions over a random 20 % of features drawn once from the parent seed, and the default strengths are `s_z = 0.5` and `s_c = 0.3` (recorded). The directions are fixed across replicates, so the confounding is the same structure in every replicate.
- **Arms.** Two `analyze_trajectories` calls per replicate, with and without covariates. Each call fits its own PLS on identical `X` with a fixed CV `random_state`, so the latent spaces are identical (asserted per replicate). The cost is two CV fits per replicate. That was accepted over adding a private latent-injection hook to the public module: at 498 features and n = 80 the CV is cheap next to 999 permutations.
- **Execution.** `--shard i/n`, with JSONL output per shard, resume by replicate index and a `merge` subcommand. Defaults are R = 500 and 999 permutations; a smoke profile uses R = 20 and 99. The cluster entry is an sbatch array like `orientation_probe_array.sbatch`.
- **Decision.** Per arm and statistic: the rate is the fraction of replicates with `p ≤ 0.05`, compared with the bound `0.05 + 2·√(0.05·0.95/R)` (0.0695 at R = 500). The decision covers the adjusted arm only (spec). The script writes `covariate_type_i_decision.json` and `covariate_type_i_rates.csv`.
- **Calibration of the construction.** Before the paper-grade run, a smoke run must show the unadjusted arm rejecting at well above α for at least one statistic. Otherwise the strengths are raised and recorded. A construction that does not confound would make the gate vacuous.

## Risks / Trade-offs

- [Covariate-blind basis: if a covariate effect is large and aligned with stage, PLS may spend a component on it] → The adjustment is still exact within the basis (D1). `selected_lv` is recorded, and the limitation is documented. The Type I check runs with covariate effects that are correlated with stage, so this case is measured rather than assumed.
- [Small cells with many covariate columns (Phase 6 cells have 9–28 samples) can make the full model near-singular or the null distribution unstable] → Exact rank deficiency is an error (D2). The record carries `design_columns`, so the case study can report the degrees of freedom used. The docs advise a parsimonious covariate set.
- [Complete-case deletion can unbalance cells and bias the analysed sample if covariates are missing not at random] → The dropped ids are recorded, and the policy is opt-in. Imputation is out of scope by decision.
- [RRPP with covariates is exact only asymptotically under the reduced model] → This is the standard RRPP property. The Type I check measures it at the actual small-n design.
- [The two-call arms in the Type I check rely on PLS CV determinism] → An identity assertion per replicate fails loudly if that ever breaks.

## Migration Plan

This is additive. The default dataset has no covariates and stays identical to 0.7.0. Release as 0.8.0 once the specs, tests and docs land. The sead_sexdiff case study pins 0.8.0 and uses covariates only after the Type I decision is `covariates_validated`. Rollback: pin 0.7.0, which leaves the case study unadjusted with the no-covariate limitation stated.

## Open Questions

- ~~The exact covariate strengths for the paper-grade Type I run depend on the smoke calibration in D10.~~ Resolved 2026-10-16. The initial strengths did not confound (the unadjusted arm rejected at about α at R = 20), and raising `s_z` alone, even to 30, did not help: the injection carries `z`'s within-cell noise along `w` in the same proportion as its cell-mean pattern, so the RRPP null widens with the signal. The confounding signal-to-noise is `a_gs` relative to the covariate noise. The calibrated construction is `a_g = 0.5`, `a_s = 1`, `a_gs = 2`, `s_z = 3`, `s_c = 2`, plus a categorical logit scale `c = 0.5` (level logits `k·c·η`) that keeps the factor mixed within cells, away from the full-rank guard. Smoke (R = 20, 99 permutations): unadjusted `shape` 0.45, `delta` 0.15. All values are recorded in `covariate_type_i_decision.json`.
