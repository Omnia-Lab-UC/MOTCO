# covariate-adjusted-trajectories

## Why

The Phase 6 case study in [Omnia-Lab-UC/sead_sexdiff](https://github.com/Omnia-Lab-UC/sead_sexdiff) now exports donor covariates beside the RNA and ATAC blocks: age at death with a `90+` flag, PMI, RIN, years and level of education, and astrocyte nuclei counts per omic. Several of them are plausibly associated with sex and with ADNC stage. Age at death in particular rises with stage. `motco analyze` (0.7.0) has no covariates, so an SEA-AD trajectory difference cannot be separated from a donor-composition difference. Covariates enter the model: they change the LS-mean fit, the RRPP reduced model and therefore what is permuted. Adjusting upstream (residualizing the blocks in sead_sexdiff) would silently fix a different null than the one MOTCO tests. The adjustment therefore belongs in MOTCO. This is roadmap Next changes item 2, and the case study waits on it.

## What Changes

- **Covariates in the real-data analysis.** `TrajectoryDataset` gains `covariate_cols` (metadata columns, in declared order), `categorical_covariates` (columns among them to treat as factors whatever their dtype) and `missing_covariates` (`"error"`, the default, or `"complete_case"`). Numeric covariates enter as one column each. Categorical covariates enter as drop-first dummies with sorted levels.
- **Covariates in both models.** The full model is `group × stage + covariates` and the reduced model is `group + stage + covariates`. Covariates enter additively. They have no interaction with group or stage. RRPP permutes the reduced model's residuals as before, so the null fixes the covariate effects.
- **LS-means at a fixed, recorded covariate value.** Every group × stage LS-mean is evaluated at the pooled mean of each covariate design column over the analysed samples. For a factor, those means are the observed level proportions. Because covariates are additive, the reference value translates every stage mean by the same vector, so `delta`, `angle` and `shape` do not depend on it. Only the reported `stage_means.csv` does, and the record says which value was used.
- **Explicit missing-covariate policy.** `"error"` rejects any missing covariate value and names the column and the samples. `"complete_case"` drops those samples from every block and from the metadata before scaling and the latent fit, and records the dropped ids. Covariates are never imputed. The finite-values contract for blocks is unchanged.
- **Design guards.** A constant covariate, a factor with one observed level, a non-numeric value in a numeric covariate, a covariate column that is also the sample, group or stage column, and a full model that is not of full column rank are rejected with an error naming the covariate(s). The last case covers a covariate that is collinear with group × stage. Complete-case deletion that empties a group × stage cell is rejected too.
- **Recorded covariate design.** When covariates are declared, `result.json` gains a `covariates` block: names, kinds, factor levels and reference level, design column names, reference values, missing policy, the samples dropped, and the covariate coefficients in the latent space. Without covariates the record, the artifact and every other output are identical to 0.7.0, except for `runtime_seconds`. The key is absent and the schema version stays 1.
- **Attribution is rejected with covariates.** Frozen-model attribution works from unadjusted cell means, so attribution together with covariates is an error in this change. Covariate-aware attribution is a follow-up.
- **The latent space stays covariate-blind.** PLS-DA is still fitted on the standardized features against stage alone. This is exact rather than approximate: OLS in a fixed linear projection equals OLS on the features followed by projection, so the adjustment of the stage means is the same either way. Only the basis is chosen without the covariates. That is documented as a limitation, together with how a sensitivity analysis can address it.
- **CLI.** `motco analyze` gains a repeatable `--covariate COL`, a repeatable `--categorical-covariate COL` (which implies `--covariate`), and `--missing-covariates {error,complete_case}`. The printed summary names the covariates and the dropped-sample count.
- **Type I validation at the Phase 6 design.** A new script, `scripts/covariate_type_i_check.py`, runs null semi-synthetic datasets at the SEA-AD MTG design (the `phase6_small_n` group × stage sizes, methylation + expression, joint magnitude, e = 0). It adds a numeric and a categorical covariate whose distributions depend on group, stage and their interaction, and whose effect is injected into the features. Each replicate is analysed with and without adjustment on a shared latent fit. The adjusted arm must hold the one-sided anchor bound `α + k·√(α(1 − α)/R)` for every statistic. The unadjusted arm is reported to show that the confounding was real. The run writes `covariate_type_i_decision.json` and a dated report. The case study may use covariates only after the adjusted arm holds.

## Capabilities

### New Capabilities

- `covariate-adjustment-type-i`: the Type I check for covariate-adjusted analysis at the Phase 6 design. It covers the confounded null construction, the paired adjusted and unadjusted arms, the acceptance bound and the decision artifact.

### Modified Capabilities

- `real-data-trajectory-analysis`: it accepts declared covariates with a missing-covariate policy, puts them in the full and reduced models, evaluates the LS-means at a recorded reference, records the covariate design, rejects attribution with covariates, and exposes the new `motco analyze` flags. The no-covariate output stays identical to 0.7.0.

## Impact

- **Code:** `src/motco/analysis.py` (dataset fields, covariate design, complete-case subsetting, record block, attribution guard); `src/motco/stats/design.py` (a covariate design-matrix helper beside `get_model_matrix`); `src/motco/cli.py` (`motco analyze` flags and summary); new `scripts/covariate_type_i_check.py`.
- **Unchanged:** the simulation harness (`evaluation.py`, `build_simulation_trajectory_design`), the stats primitives' signatures and defaults, `motco de`, every committed study config, and the byte-identity fixtures. `RRPP` and `estimate_difference` already accept arbitrary full and reduced model matrices and LS-mean rows, so they are reused as they are.
- **Tests:** `tests/test_analysis.py` and `tests/test_analysis_cli.py` gain covariate cases. The existing harness-equivalence tests must pass unchanged.
- **Docs:** `docs/api/analysis.md`, the CLI reference, the `docs/roadmap.md` item 2 status and the case-study item, and the `analysis.py` paragraph in `AGENTS.md` ("No covariates" becomes the covariate contract).
- **Release:** a minor version (0.8.0). sead_sexdiff pins it.
