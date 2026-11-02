# Design

## Context

See proposal.md (Why) for the motivation. The requirements are in the specs: `orientation-driver-attribution`, `real-data-trajectory-analysis` and `adjusted-attribution-validation`. The change builds on these pieces of 0.13.0.

**Attribution (`motco.stats.attribution.analyze_orientation_attribution`).**
- The means come from one of two sources:
  - `_arithmetic_means`, the per-cell feature means of the input rows (`mean_source = "arithmetic"`);
  - a precomputed `mean_table` (`"precomputed"`), which forces `bootstrap_replicates = 0`.
- `_reconstruct_means` passes the means through `model.transform` and `model.inverse_transform`. Both are affine for `PLSRegression`.
- `_build_transitions` and `_build_principal_orientation` build the observed, PLS-captured and residual components from the means.
- `_bootstrap` draws `rng.choice(stratum, size=len(stratum), replace=True)` once per group × stage stratum, in group-major, stage-minor order, from `default_rng(seed)`. It then recomputes the arithmetic means and the contrasts. A replicate with no defined contrast contributes a NaN row and is not counted as valid.
- `attribution_frames` writes the `configuration` and `interpretation` tables from `AttributionConfig` and `InterpretationMetadata`.

**The compact record (`motco.simulations.attribution_diagnostics.run_attribution_diagnostics`).** It calls the primitive and reduces the result to a record at `ATTRIBUTION_SCHEMA_VERSION = 2`. The simulation harness and `motco analyze` both use it.

**The analysis (`motco.analysis.analyze_trajectories`).**
- It raises before fitting when covariates are declared and attribution is enabled (`analysis.py:519–523`).
- `_run_attribution` passes `latent.artifacts.model` and `latent.artifacts.features` (the pooled-standardized joint matrix) to the record builder, with an unavailable `TruthDrivers`.
- The covariate design (`stats.design.CovariateDesign`) carries `matrix` (n × q), `column_names` and `reference_row`, the column means.
- `_covariate_adjusted_design` appends the covariates to `[1, g, s, g×s]` and the reference row to every LS-mean row.
- `TrajectoryAnalysisResult` holds the record, the latent scores, the loadings, the stage means and the attribution result. It does not hold the estimator.

**The within-cell coefficient fit (`motco.stats.pls.fit_nuisance_coefficients`).**
- It regresses `X` on `[1, conditioning, covariates]` and returns the covariate block γ̂ (q × p) and the residual df.
- It raises a `ValueError` naming the columns on rank deficiency.
- With the group × stage cell indicators as conditioning, γ̂ is the covariate coefficient of the analysis's full model, by the Frisch–Waugh–Lovell theorem.

**The case-study construction (`scripts/case_study_type_i_check.py`).**
- `build_unit(mode, effect, replicate, settings)` wraps `covariate_type_i_check.build_replicate`. That function generates the dataset and injects `s_z·sd_f·w_f·z` into the analysed scale (methylation as M-values, `identity` transform), with `s_c = 0`. `build_unit` then adds five noise covariates.
- `build_replicate` returns the seeds, metadata, blocks and group labels. It does not return the generator dataset, so the truth drivers are currently out of reach.
- The corrected-permutation power profile ran this construction at δ = 2.0 on the astrocyte table at about 14 s per analysis.

**Seeds.** The committed base seeds are 100, 400, 500, 600, 800, 900, 1028, 2026, 2031 and 2032. 2033 appears only as a hash-change argument in `tests/test_adni_design_check.py`.

**Byte identity.** No committed fixture pins attribution outputs byte for byte. `tests/test_analysis_identity.py` runs without attribution, and the harness identity fixtures pin statistics, not attribution tables.

## Goals / Non-Goals

**Goals:**

- One adjusted-means definition shared by the primitive, the analysis and the bootstrap, equal to the latent LS-means the analysis already reports, pushed back into feature space.
- A bootstrap whose uncertainty includes the covariate-coefficient estimation, because at n = 80 with six covariates that is the extra uncertainty adjustment brings.
- Attribution outputs without covariates that are byte-identical, pinned by a fixture captured before any edit.
- A validation that runs the shipped analysis path and compares it with an oracle on identical data, bootstrap draws and frozen model.

**Non-Goals:**

- Attribution on a latent-adjusted basis (D6).
- Covariate × group or covariate × stage interactions. The trajectory models are additive, and attribution follows them.
- Changing the harness's attribution, its record schema version or any study config.
- Absolute driver-recovery thresholds on real data. They are item 12's lift criterion (iii), fixed in its proposal.

## Decisions

### D1. Adjusted means in the primitive, in the within-cell form

`analyze_orientation_attribution` gains two keyword-only arguments:
- `covariates`: a DataFrame or array of numeric design columns, row-aligned with `X`;
- `covariate_reference`: a mapping, Series or sequence over those columns. The default is the column means.

With `covariates`, every cell mean is

```text
m̂_gs = x̄_gs − γ̂ᵀ (c̄_gs − c_ref)
```

with γ̂ from `fit_nuisance_coefficients(X, NuisanceDesign(conditioning = drop-first group × stage cell indicators, covariates = C))`. It is exactly the least-squares mean of `X ~ cells + C` at `c_ref` (FWL), so it equals the analysis's full-model LS-mean, whose column space is the same. Three consequences:
- **Reference-row invariance.** `c_ref` enters every cell through the same γ̂ᵀc_ref, so transitions, path lengths and both contrast blocks are invariant. Only the means tables move.
- **Projection identity.** `transform(·)` is affine, and every LS-mean row has intercept 1. So `model.transform(m̂)` equals the analysis's latent `stage_means` up to floating point.
- **The validity check is the rank check.** `[cells, C]` has full rank exactly when the within-cell-centred `C` does. `fit_nuisance_coefficients` reports a failure, and the primitive re-raises it as an `AttributionError` naming the columns.

`mean_table` and `covariates` together are an error. The categorical coding stays the caller's job: the analysis passes `CovariateDesign.matrix` with its `column_names`.

*Alternative: compute the adjusted means in the analysis and pass `mean_table`.* Rejected. The bootstrap needs rows, which is exactly the guard at `attribution.py:368–372`.

*Alternative: residualize the features (X − Cγ̂) and take arithmetic means.* It gives the same contrasts, but the means sit at c_ref = 0, so the projection identity with `stage_means` would fail by a constant. It is equivalent once γ̂ᵀc_ref is added back, which is the form above.

### D2. The bootstrap draws the unadjusted resamples and refits γ̂ in each

The adjusted path reuses `_bootstrap`'s draw loop unchanged: the same `default_rng(seed)`, the same strata order, the same `rng.choice` calls. So for one seed it resamples the same row indices as attribution without covariates. Per replicate, it:
- forms the resampled `X*` and `C*`;
- fits γ̂* with `fit_nuisance_coefficients` on them;
- evaluates `m̂*` at the full-sample `c_ref`;
- reconstructs through the frozen model, and builds the contrasts with the replicate's own net-displacement sign anchor.

A `ValueError` from the fit, meaning within-cell collinearity in the resample, makes the replicate invalid. It contributes a NaN row to every block and component, exactly as a zero transition does, and increments one `rank_deficient_replicates` counter. Preprocessing and the PLS model stay frozen.

Why the same draws: the validation's three arms, and any later adjusted-versus-unadjusted comparison, are then paired at the bootstrap level too. Without covariates the draw loop is untouched (D3).

*Alternative: freeze γ̂ at its full-sample value.* Rejected. It would leave out the γ̂ estimation error, which at n = 80 with q = 6 (68 residual df in the full model) is the uncertainty adjustment adds. The exit review asks for a refit.

*Alternative: stratify on cell × categorical level, or drop an unidentified column per replicate.* Rejected. The first changes the resampling scheme and does not cover numeric degeneracy. The second changes the estimand replicate by replicate.

### D3. Byte identity without covariates, pinned by a new fixture

With `covariates=None`, the code path is the current one: `_arithmetic_means`, the current `_bootstrap` branch, and no new configuration or interpretation rows. The new fields are optional and default to `None`:
- `AttributionConfig.covariate_columns`, `covariate_reference` and `rank_deficient_replicates`;
- `InterpretationMetadata.adjustment_boundary`.

`attribution_frames` emits their rows only when they are set. `ATTRIBUTION_SCHEMA_VERSION` stays 2, and the compact record's `adjustment` block appears only when covariates are given. The harness never passes covariates, so its records cannot change.

A new fixture, `tests/data/attribution_no_covariate_identity.json`, is captured before any code edit by running a new `tests/test_attribution_identity.py` as a script, as `tests/test_analysis_identity.py` was captured. It holds two things:
- the primitive's frames on a fixed synthetic dataset and frozen model, with 25 bootstrap replicates;
- the attribution block of an `analyze_trajectories` record, minus runtimes, and its attribution CSVs, on the `test_analysis_identity` dataset.

The comparison uses `tests.test_phase6_byte_identity._assert_close`: structure, strings and integers exact, floats at a tight relative tolerance.

### D4. What the adjusted result records

- `AttributionConfig.mean_source = "covariate_adjusted"`, plus the covariate columns, the reference row and the rank-deficient count. They become configuration rows `covariate_columns`, `covariate_reference` and `bootstrap_rank_deficient_replicates`.
- `InterpretationMetadata.adjustment_boundary`, with this text: "Adjusted effects are group × stage least-squares means of an additive linear covariate model evaluated at the reference covariate row; bootstrap replicates refit the covariate coefficients. Adjustment removes only the declared covariates' linear contribution and does not remove unmeasured confounding."
- The compact record gains `adjustment` with these fields: `mean_source`, `design_columns`, `reference_row` (a column → value mapping), `rank_deficient_replicates` and `bootstrap_replicates`.

### D5. The analysis passes its covariate design and exposes its measurement

- **Guards.** The covariates guard in `analyze_trajectories` is replaced by a latent-adjustment guard (D6), also raised before fitting.
- **The attribution call.** `_run_attribution` receives the `CovariateDesign`. It passes `covariates = DataFrame(matrix, columns=column_names)` and `covariate_reference = reference_row` through `run_attribution_diagnostics`, which forwards them to the primitive. The rows are the analysed samples in metadata order, the same rows the features have, after any complete-case deletion.
- **`TrajectoryAnalysisResult.latent_artifacts`.** A new field holding a `PLSIntegrationArtifacts | None`, default `None`, set by `analyze_trajectories`. It is excluded from `record`, `to_json` and `write_analysis_outputs`. Its docstring states that it is in memory only and never serialized.

The validation needs the exact frozen estimator and features for its oracle and unadjusted arms, and item 12 needs them for its planted-signal arms.

*Alternative: the script refits the basis with `fit_pls_latent_space` at the selected rank.* Rejected. It duplicates the analysis's measurement and invites silent drift. The exposed object is the measurement itself, and the consistency check (`model.transform(features)` against the latent scores, to rounding) holds by construction.

- **CLI.** When the record's attribution carries `adjustment`, the CLI prints `Attribution: covariate-adjusted means at the reference row; <k> of <B> bootstrap replicates rank-deficient`.

### D6. Attribution with latent adjustment stays rejected

The latent-adjusted basis is fitted on covariate-residualized features with `scale=False`, and `latent.matrix = model.transform(X)` on the raw features. Its `angle` and `shape` are `residualized_basis`, which makes them inadmissible at n = 80. So no claim that attribution would explain is read there. Running attribution on that basis would also mix the residualized basis with adjusted raw-feature means in the captured component, a quantity no requirement defines and no check validates.

The combination is rejected before fitting, with the message "Attribution is not available with latent adjustment …". It was already unreachable, because latent adjustment requires covariates and covariates rejected attribution, so no behaviour that worked before is removed.

### D7. Validation construction: the case study's configuration at the Phase 6 stage effect

`scripts/adjusted_attribution_check.py` imports `case_study_type_i_check as cs` and `covariate_type_i_check as d10`, as `scripts/corrected_permutation_check.py` does. The construction:
- **Design:** the astrocyte `cohort_30` table `[[11, 10, 28], [9, 10, 12]]`; methylation + expression (498 joint features); joint magnitude.
- **Stage effect:** δ = 2.0, the Phase 6 stage effect.
- **Covariates:** `z` confounded with group × stage (D10 strengths) and injected on a 20 % feature subset, `s_c = 0`, plus five noise covariates (six design columns).
- **Analysis settings:** rank floor 2, net displacement, `df`, 999 permutations.
- **Seeds:** base seed 2033.
- **Cells:** the anchor, and orientation at e = 0.5 and 1.0, 500 replicates each. Bootstrap B = 200 with seed `_seed(base, "bootstrap", replicate)`, and top-k 20.

`d10.build_replicate` additionally returns the generated dataset under a new key. It is never written to a record, so the existing checks' records and tests are unchanged. `derive_truth_driver_features` reads the dataset, and the principal block is scored against the union of the per-transition drivers, as `attribution_diagnostics` does.

Why δ = 2.0 rather than the weak case-study signal: attribution explains a planted orientation. At δ = 0.05–0.1 the stage transitions are near zero, and every arm's contrast is noise. How the method behaves at SEA-AD's real stage signal is item 12's question, answered on the real blocks.

Why these cells:
- e = 1.0 is the construction maximum, where net-displacement `angle` power is 0.648 at this design.
- e = 0.5 is a weak planted orientation (power 0.084–0.178). It tests that adjustment does not degrade a weak contrast.
- The anchor tests that adjustment removes covariate-driven pseudo-drivers when no orientation is planted.

Why only the astrocyte table: it is item 12's cohort.

### D8. Three arms per unit on one measurement

1. **`adjusted`.** `analyze_trajectories` runs with the six covariates and `attribution = AttributionDiagnosticSettings(enabled=True, bootstrap_replicates=200, bootstrap_seed=…, top_k=20)`, `residual_correction="df"` and the unit's RRPP seed. This is the shipped path end to end.
2. **`oracle`.** The primitive runs on `X̃_or = X̃ − z ⊗ (v / σ)` with the same model, metadata, bootstrap seed, B and top-k, and no covariates. Here:
   - `v_f = s_z · sd_f · w_f` is the raw-unit injection, with `sd_f` the pre-injection SD that `build_replicate` used;
   - `σ_f` is the analysis's pooled scale, `latent_artifacts.original_scales`.

   The oracle's cell means are therefore the uninjected cell means in the analysed coordinates, up to a constant shift.
3. **`unadjusted`.** The primitive runs on `X̃` as analysed, with the same settings. This is the existing path.

All three share D2's row resamples.

Each unit runs two consistency checks before it records anything, and raises if either fails:
- `model.transform(latent_artifacts.features)` equals the latent scores to within 1e-10 of the largest absolute score. Not bitwise: without latent adjustment the scores are the fit's `x_scores_`, which agree with `transform` only to rounding (3.6e-15 on the identity dataset, found at implementation);
- `model.transform(adjusted observed means)` equals `stage_means` to a relative 1e-8.

*Alternative: the oracle through `analyze_trajectories` on uninjected blocks.* Rejected. The pooled scaler and the basis would both change, so the arms' feature contrasts would not be comparable feature by feature.

### D9. Metrics are read on the observed principal-orientation contrast

The decision reads the **observed** component of the **principal** block:
- The principal block decomposes the tested `angle`.
- The Phase 4 pilot found that the PLS-captured component retains about 6 % of the observed orientation contrast, so driver reports use the observed component.

The PLS-captured and residual components and the per-transition blocks are recorded and reported descriptively. The metric definitions are in the spec. Three definitions need comment:
- **Contamination** counts injected features that are *not* truth drivers, because ranking an injected genuine driver is correct.
- **The projection** is on `ŵ = (v/σ)/‖v/σ‖`, the injection direction in analysed coordinates, so it measures confounding left in the whole contrast vector rather than only in the top-k.
- **Truth precision** divides by min(k, |truth|), so a small driver set cannot cap it.

Per-feature mean differences against the oracle are not used as a bias test. The contrast q is a difference of unit vectors, and a noisier estimate of a transition is attenuated toward zero after normalization. A per-feature t-test at R = 500 would then flag noise as bias. The projection on ŵ and the contamination share target the confounder itself.

### D10. Predeclared margins and verdicts

All bounds are one-sided at 2 paired standard errors per cell, over the cell's R available units, the convention of the anchor bounds elsewhere. Units within a cell are independent replicates. The margins are fixed here, before any run, and changing one needs the user's approval recorded in this file:

| quantity | bound | why |
|---|---|---|
| contamination, adjusted − oracle | ≤ 0.05 | one feature in twenty |
| residual confounding along ŵ | ≤ 0.10 × the unadjusted arm's | adjustment removes at least 90 % of the confounding the unadjusted contrast carries |
| truth precision, adjusted − oracle | ≥ −0.05 | one feature in twenty |
| bootstrap sign stability, adjusted − oracle | ≥ −0.05 | as precision |
| top-k frequency, adjusted − oracle | ≥ −0.10 | two features in twenty. The adjusted arm estimates six extra coefficients on 68 residual df, which widens per-feature bootstrap spread more than the sign. |
| gate: contamination, unadjusted − oracle | ≥ 0.25 at every cell | the confounder must visibly reach the unadjusted top-k, or nothing is being removed |
| gate: oracle truth precision at e = 1.0 | ≥ 0.5 | the planted drivers must be recoverable at all, or non-inferiority to the oracle is vacuous |

The verdicts and their consequences are in the spec. In short:
- `adjusted_attribution_validated` lets item 12 start on 0.14.0;
- `diagnose_adjusted_attribution` keeps the leg withheld and calls for a diagnosis proposal;
- `construction_uninformative` keeps the leg withheld and returns the construction to the user.

A verdict below 999 permutations, B = 200 or R = 500 is indicative.

### D11. Calibration, timing and where it runs

Calibration runs R = 20 × 99 permutations × B = 50 at the anchor and orientation e = 1.0, starting from `s_z = 10`, the corrected-permutation power profile's value.
- `s_z` is raised in steps of 5 until the contamination criterion is met in both cells, and every iteration goes to `calibration.json`.
- `s_z` enters every arm's pooled scale identically, so the arms stay paired at any value. An injected feature's clean signal shrinks in standardized units as `s_z` grows, in all three arms alike.
- If the oracle recoverability criterion fails at a strength that meets contamination, no further `s_z` is tried. The check stops, the user decides, and the decision is recorded here.

**Timing estimate.** About 14 s per adjusted analysis at δ = 2.0. Three attributions at B = 200 and p = 498 add an estimated 1–3 s. Over 1,500 units that is about 6–7 core-hours.

**Where it runs.** The rule is fixed now and applied to the smoke timing:
- if the measured total is at most 12 core-hours, the paper-grade run executes locally with process workers from a clean worktree at the committed tree, with the tree hash in `PROVENANCE.txt`;
- otherwise it runs as arrays on `ing` with a new `scripts/adjusted_attribution_array.sbatch`, modelled on `scripts/corrected_permutation_array.sbatch`.

Either way the plan goes to the user for approval first.

### D12. Version and naming

The version is bumped to 0.14.0 before any paper-grade unit, so the records name the release that carries the feature. The names are:
- script `adjusted_attribution_check.py`, with subcommands `calibrate | run | merge | decide`;
- output directory `results/adjusted-attribution-<date>/`;
- decision file `adjusted_attribution_decision.json`;
- report `docs/reports/adjusted-attribution-<date>.md`.

## Risks / Trade-offs

- **[The oracle cannot recover the planted drivers at n = 80, so the recoverability gate fails.]** → Calibration surfaces this before any paper-grade unit and stops for the user. A larger table would no longer be the case-study configuration, so it is the user's call, not a silent adjustment.
- **[At e = 1.0 a strong orientation crowds the injected features out of the unadjusted top-k, so the contamination gate fails at that cell only.]** → Calibration checks both cells. Raising `s_z` restores the gap. The projection criterion does not depend on the top-k.
- **[The adjusted contrast is noisier than the oracle's, so its normalized contrast is attenuated.]** → The criteria are non-inferiority margins, not equality. No per-feature t-test is used (D9).
- **[Real covariates include factors whose rare levels make resamples rank-deficient.]** → Such replicates are invalid, counted, recorded and printed, and the valid counts say how many contributed. The validation uses numeric covariates only, so this is unit-tested, not measured, and the report states it as a limitation.
- **[The additive linear covariate model is misspecified on real data.]** → The interpretation boundary states the model. Interactions stay out of scope, as in the trajectory models.
- **[Exposing the estimator invites someone to pickle it into the artifact.]** → The field is excluded from the record and the writer, and a test pins the written file set.
- **[Attribution outputs were never pinned byte for byte, so a regression in the no-covariate path would pass unnoticed.]** → The fixture is captured first, on the unmodified tree (task 1.2).
- **[The InterSIM blocks (498 features) are not SEA-AD's 10,000.]** → The structural-analogue limitation is stated. Item 12 runs on the real blocks. Bootstrap cost is O(B · n · p), about 10 MFLOP per replicate at p = 10,000, so it does not block item 12.

## Migration Plan

The API change is additive. `analyze_orientation_attribution` keeps its defaults, and an analysis that raised before (covariates with attribution) now returns a result. `tests/test_analysis.py::test_attribution_with_covariates_is_rejected_before_fitting` becomes a latent-adjustment rejection test, and a covariate-adjusted attribution test is added beside it. Documentation that says attribution is unavailable with covariates is updated in the same change: `docs/api/analysis.md` and `AGENTS.md`. Rollback is a revert. No committed result, config or fixture changes.

## Open Questions

- The calibrated `s_z`. It is fixed by D11's procedure and recorded, and does not change the specs.
- Whether item 12 reuses this check's metric code. That is decided in item 12's proposal.
