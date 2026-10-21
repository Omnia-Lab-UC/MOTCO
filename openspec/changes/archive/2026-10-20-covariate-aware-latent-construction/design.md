# Design

## Context

See proposal.md (Why) for the motivation. The pieces this change touches:

- `analyze_trajectories` (0.8.0) standardizes the blocks with pooled per-feature scalers, concatenates them, calls `fit_pls_latent_space` (the core shared with the harness) with the stage labels, then builds the trajectory design, appends the declared covariates to both models through `_covariate_adjusted_design`, and runs `estimate_difference` and `RRPP` on `latent.matrix`.
- `fit_pls_latent_space` calls `plsda_doubleCV` for the rank and `fit_plsda_model` for the final fit, and returns the training scores `model.x_scores_` as the latent matrix. Every PLS fit, in the folds and at the end, is `PLSRegression(scale=True)`.
- `plsda_doubleCV` is a fixed loop: `RepeatedStratifiedKFold` outer, `StratifiedKFold` inner, `_plsda_auroc` per candidate rank on `(X_train, X_val)` sub-frames, then per-repeat refits on the full input.
- `CovariateDesign` (`build_covariate_matrix`) already provides the covariate matrix `C`, its column names and its reference row; complete-case deletion already happens at dataset construction, so `X`, `C` and the labels are row-aligned.
- Two invariants constrain the implementation: `tests/test_analysis.py::test_matches_harness_*` (the analysis equals the harness on identical data) and `tests/test_analysis_identity.py` (the no-covariate output equals 0.7.0). This change adds a third: the no-setting output equals 0.8.0.

Two facts, verified numerically before this proposal, shape the design:

- `PLSRegression(scale=True)` re-standardizes its input internally and is invariant to any per-feature rescaling of it. Passing residualized features into the unchanged machinery is therefore the "residualize then re-standardize" option, whether or not one intends it.
- With the covariates in both RRPP models, projecting the raw or the residualized features onto the same basis gives identical observed statistics and identical permutation draws, because the subtracted term `C γ W` lies in the column space of `C`.

## Goals / Non-Goals

**Goals:**

- The residualizer lives in the shared PLS core, so the harness's measurement function carries it, and the analysis path is one argument away from the harness path.
- Leak-free CV at both levels without duplicating the CV loop.
- Byte identity of the no-setting path at the level of `record["params"]`, the CSVs and the harness.
- Validation scripts that reuse the D10 construction and the study engine's acceptance conventions, and that can be sharded on SLURM.

**Non-Goals:**

- Covariates in `SimulationEvaluationParams`, the study config schema or the generator. The validation scripts inject covariates at script level, as D10 does.
- A degrees-of-freedom floor. Full rank is the only hard guard; the Phase 6 check at 8–10 covariates measures the rest.
- The rank-1 `angle` degeneracy, caller-chosen conditioning designs beyond the two named, and covariate-aware attribution.

## Decisions

### D1. The residualizer is a per-fold transformer inside `plsda_doubleCV`, driven by a `NuisanceDesign`

A small frozen dataclass `NuisanceDesign(conditioning: np.ndarray, covariates: np.ndarray, conditioning_names, covariate_names)` in `stats/pls.py` holds the row-aligned conditioning columns (without intercept) and covariate columns. `plsda_doubleCV(..., *, scale=True, nuisance=None)` gains two keyword-only arguments. With `nuisance` given, the loop residualizes at three places: the outer rest set (coefficients from the rest rows, applied to rest and test), each inner training set (coefficients from the training rows, applied to training and validation), and the per-repeat full refits (full-sample coefficients). `_plsda_auroc` gains `scale`. A helper `fit_nuisance_coefficients(X, design) -> γ` does one `lstsq` on `[1, D, C]` after a `matrix_rank` check and returns the `C` block; `residualize(X, C, γ)` subtracts.

A third keyword-only argument, `strata=None`, names the labels both CV levels stratify on instead of `y`. The residualized path passes the conditioning cells (`NuisanceDesign.cells()`, the distinct conditioning rows: the stages under `stage`, the group × stage cells under `group_stage`), so every training set contains every cell the conditioning design must estimate. *Found during implementation:* the inner `StratifiedKFold` is unshuffled and stratifies on stage, so on rows sorted by group within stage (the generator's order, and plausibly a real metadata table's) one inner validation fold takes a whole group × stage cell and the `group_stage` training design loses an interaction column — 100 of 100 Phase 6 replicates at the harness CV settings. Stratifying on the conditioning cells removes that structural failure (to the 8 of 100 rare-level failures both conditionings share; see Risks). Folds depend only on the partition the labels define (sklearn encodes classes by order of appearance), so under `stage` the folds are exactly those of the unadjusted fit; a test pins this.

*Alternative: a sibling function `plsda_doubleCV_residualized`.* Rejected: it would duplicate the loop that the `plsda-nested-cv` spec pins, and the two would drift. The spec is amended instead to allow keyword-only arguments whose defaults reproduce the old behaviour, which `tests/test_pls.py` checks on a fixed input.

*Alternative: residualize once on the full sample and pass `X_res` to the unchanged CV.* Rejected: the held-out samples' labels would have entered their own residualization through γ, making the CV AUROC and the rank selection optimistic (question 3 of the design note).

### D2. Conditioning designs come from a helper beside `build_covariate_matrix`

`build_conditioning_matrix(frame, group_col, stage_col, kind) -> ConditioningDesign` in `stats/design.py` returns the drop-first stage indicators (`stage`) or the drop-first group, stage and interaction indicators (`group_stage`), with names, using the same `_sort_levels` ordering as `get_model_matrix`. In `analysis.py` it is called on the recoded `working` frame, so the ordinal recoding and the level mapping are reused. The nuisance design for the CV is `NuisanceDesign(conditioning.matrix, covariates.matrix, ...)`.

*Alternative: slice `design.model_full` / `model_reduced` from `build_simulation_trajectory_design`.* Rejected: those are built after the latent fit today, and their intercept column and column order are harness details the residualizer should not depend on.

### D3. `fit_pls_latent_space` gains `nuisance=` and projects the raw features when it is given

With `nuisance` given, the function: fits `γ_full` on all samples (full-rank check, error naming the covariates); residualizes `X` to `X_res`; runs `plsda_doubleCV(X, ..., scale=False, nuisance=nuisance, strata=nuisance.cells())` for the rank; fits `fit_plsda_model(X_res, y, n_components=selected_lv, scale=False)`; and returns `latent.matrix = model.transform(X)` (centered at `X_res`'s mean, no scaling because `scale=False`, times `x_rotations_`). Without `nuisance` the function is untouched, including `scale=True` and `x_scores_`.

The forced-rank path (`forced_components`) is honoured with `nuisance`: no CV runs, the full-sample residualization still happens, and the record states `component_selection="forced"`.

The metadata gains `integration_params["latent_adjustment"] = conditioning` only when `nuisance` is given, and the function returns a `LatentAdjustmentSummary` (residual df, CV minimum training df, per-feature R² summary, design column names, feature scale) on the result's artifacts for `analysis.py` to record. The CV minimum training df is collected by the residualizer in `plsda_doubleCV` and returned in a new optional key of its result dict only when `nuisance` is given, so the table schema and positional contract are unchanged.

*Alternative: project the residualized features (`model.x_scores_`).* Rejected on readability, not inference: the two give identical statistics and p-values (Context), but `X W'` keeps the covariate coefficients in the `covariates` block meaningful and `stage_means.csv` on the same footing as 0.8.0. A test pins the equivalence so the choice stays free.

### D4. Keep the pooled scale: `scale=False` on the residualized path only

The residualized path fits every PLS with `scale=False` on the pooled-standardized residuals. Re-standardization would give a feature whose variance was mostly covariate the same weight as a clean one and would inflate exactly the features where γ is noisiest (the case study's R² reaches 0.53 at the 95th percentile on RNA; the Type I construction goes higher). With the pooled scale such features are down-weighted, and when the covariate effects are weak the basis converges to the 0.8.0 basis, which is the no-loss case the power comparison checks. The no-setting path keeps `scale=True` so it stays byte-identical; the folds therefore re-standardize on the training fold there and keep the pooled scale here, a difference the record names as `feature_scale: "pooled_standardized"`.

### D5. The setting is a params field omitted from `to_dict()` while `None`

`TrajectoryAnalysisParams.covariate_latent_adjustment: str | None = None`, validated in `__post_init__` against `{"stage", "group_stage"}`. `to_dict()` drops the key while it is `None`, following `grid._SIGNATURE_OMIT_WHEN_NONE`, so `record["params"]` of every 0.8.0 record is unchanged. The setting describes the analysis, not the data, so it does not belong on the dataset (D4 of the 0.8.0 design), but D4's byte-identity reason still binds, hence the omission. `analyze_trajectories` raises before any fitting when the setting is on and `dataset.covariate_cols` is empty.

### D6. Record block

```json
"latent_adjustment": {
  "conditioning": "stage",
  "covariates": ["age_at_death", "pmi", "education_level"],
  "design_columns": ["stage[1]", "stage[2]", "age_at_death", "pmi", "education_level[Graduate]", "education_level[HighSchool]"],
  "n_samples": 80,
  "residual_df": 73,
  "cv_min_training_df": 27,
  "feature_scale": "pooled_standardized",
  "feature_r2": {"q05": 0.04, "q25": 0.10, "q50": 0.19, "q75": 0.30, "q95": 0.53, "mean": 0.21, "fraction_above_0.5": 0.06}
}
```

Inserted after `covariates`, only when the setting is on. `latent.integration_params.latent_adjustment` carries the conditioning. `ANALYSIS_SCHEMA_VERSION` stays 1: both keys are optional and their absence is the 0.8.0 contract. The `covariates` block is written as in 0.8.0, from the full-model fit on the projected raw features.

### D7. Byte-identity fixture before any code change

The first task captures `tests/data/analysis_no_latent_adjustment_identity.json` from the 0.8.0 tree: a seeded two-block, three-stage fixture with two covariates (one numeric, one categorical), 99 permutations, attribution off. The file holds the record without `measurement.runtime_seconds` and `software`, plus SHA-256 digests of the three CSVs. The existing no-covariate fixture stays and keeps passing, so the identity chain is 0.7.0 → 0.8.0 → this change.

### D8. Type I check: three arms on the D10 construction, with a calibration gate

`scripts/latent_adjustment_type_i_check.py` imports the D10 construction from `scripts/covariate_type_i_check.py` (data, covariates, injection, seeds, sharding, merge) rather than copying it, and replaces `evaluate_replicate` with three `analyze_trajectories` calls: `None`, `"stage"`, `"group_stage"`, all with both covariates declared, the same CV `random_state` and RRPP seed. Each record keeps the selected rank, the first component's Spearman correlation with the numeric covariate, and the latent-adjustment block. The calibration criterion (median |ρ| ≥ 0.8 in the baseline arm) is evaluated on the smoke profile (R = 20, 99 permutations) and the strengths recorded; D10's calibrated strengths (`s_z = 3`, `s_c = 2`) are the starting point and may already satisfy it. The decision step writes `latent_adjustment_type_i_decision.json` and `latent_adjustment_type_i_rates.csv`. The sbatch array follows `covariate_type_i_array.sbatch`.

Cost: three CV fits and three RRPP runs per replicate, 1.5× the D10 run (which took about 60 shard-hours on ing).

### D9. Power comparison: two columns, matched seeds, `column_decision` rule

`scripts/latent_adjustment_power_comparison.py` reuses the same construction with `group_effect_size` at the Phase 6 effect levels for `magnitude`, `orientation` and `shape`, in two columns: `nuisance` (the calibrated injection) and `clean` (`s_z = s_c = 0`, covariates still declared). The same replicate index draws the same dataset, covariates and RRPP seed in every arm. The decision reuses the study engine's pooled-SE multipliers and verdict vocabulary: target pairs orientation/`angle`, magnitude/`delta`, shape/`shape` at the top effect in the `nuisance` column (`gain`), every statistic at every effect in the `clean` column (`no_loss`), and the zero-effect anchor and translation held to the anchor bound. It writes `latent_adjustment_power_decision.json` and a rates table with the rank distribution per column, arm and effect. Cost: about 3 arms × 2 columns × 3 modes × 4 effects × 500 replicates, plus anchors, which is of the order of the 2026-10-14 net-displacement power comparison (384 core-hours); the exact grid is settled on the smoke run and reported before submission. Nothing is submitted without the user's approval.

### D10. Harness equivalence stays a no-setting property

`test_matches_harness_*` compares the analysis with the harness on identical data without covariates; it is unchanged. The harness has no covariates, so the setting has no harness counterpart, and the spec's equivalence requirement is read as applying to the shared settings. The residualizer in `fit_pls_latent_space` is exercised by the analysis tests and by a direct unit test on the core.

## Risks / Trade-offs

- [A factor level absent from a CV training set makes the fold rank-deficient, and which fold depends on `random_state`] → The error names the repeat and fold. The docs advise merging rare levels upstream. A seed-dependent failure is preferred to a silent minimum-norm fit, whose held-out shift would be arbitrary. In the D10 construction the categorical covariate's lowest level is rare enough to do this in about 8 % of Phase 6 replicates (both conditionings); the Type I and power scripts record such an arm as unavailable for that replicate and report the count, so R is the number of available replicates.
- [A conditioning cell absent from a CV training set] → Prevented by stratifying both CV levels on the conditioning cells (D1).
- [`group_stage` conditioning lets the group label touch the basis through a pooled within-cell slope] → Measured, not assumed: it is an arm of the Type I check, and the decision adopts `stage` whenever `stage` holds.
- [`stage` conditioning leaves the group-correlated part of a covariate's effect in the features] → The basis is stage-supervised, so this affects which variance PLS sees, not its response; the power comparison's `nuisance` column measures whether enough of the axis is removed.
- [Residualizing a mediator of stage removes real signal] → Documented limitation in the API docs, the report and the roadmap. No simulation distinguishes a confounder from a mediator; the case study must justify its set.
- [Many covariates at n = 80 leave about 27 training-set df in the inner folds] → Recorded as `cv_min_training_df`; the Type I check runs at 8–10 covariates so the case study's design is covered.
- [`scale=False` in the folds changes the fold metric relative to the no-setting path] → By design (D4), recorded as `feature_scale`; the no-setting path is untouched and pinned by the fixture.
- [The case study's technical variables must also enter the RRPP model] → A consequence of the equality requirement that makes the projection choice immaterial; stated in the docs. The 0.8.0 full-rank guard covers the larger model.

## Migration Plan

Additive. The default setting is `None`, pinned to 0.8.0 by the new fixture and to 0.7.0 by the existing one. Release as 0.9.0 once the specs, tests, docs and the Type I decision land; the power comparison may follow in its own report without blocking the release. sead_sexdiff pins 0.9.0 and uses the setting only after its own proposal names a covariate set and the Type I decision is `adopt_stage` or `adopt_group_stage`. Rollback: pin 0.8.0.

## Open Questions

- Whether D10's calibrated injection already makes the covariate-blind basis spend its first component on the nuisance (median |ρ| ≥ 0.8). If not, the smoke run raises `s_z` and records it; the specs and tasks do not change.
- The exact power grid (effect levels per mode, anchors) within the cost envelope of D9, settled on the smoke run and reported before cluster submission.
