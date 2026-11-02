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
| gate: contamination, unadjusted − oracle | ≥ 0.25 at the anchor (amended; declared at every cell) | the confounder must visibly reach the unadjusted top-k, or nothing is being removed |
| gate: \|projection on ŵ, unadjusted − oracle\| − 2·SE | ≥ 0.2 at every cell (added by the amendment) | the confounder must reach the unadjusted contrast in every cell the criteria read |
| gate: oracle truth precision at e = 1.0 | ≥ 0.5 | the planted drivers must be recoverable at all, or non-inferiority to the oracle is vacuous |

**Amendment of the gates, 2026-11-02 (approved by the user, recorded before any smoke or paper-grade unit).** Calibration under the declared rule ended `stop_exhausted` (D11, iterations 0–8): at orientation e = 1.0 the unadjusted − oracle contamination was 0.0 at every `s_z` from 10 to 50. A diagnostic (`results/adjusted-attribution-2026-11-02/calibration/probe/`, six replicates per row, indicative) showed why. At e = 1.0 the planted drivers are 302 of the 498 features, and 61 of the 100 injected features are among them, so only about 39 injected features can count as contaminants; their best rank in the unadjusted contrast was 30th–90th at δ = 2.0. Weakening the stage effect to δ = 0.25 did not change that (0.0–0.07). Meanwhile the unadjusted − oracle projection on ŵ was about 0.5 in every cell, against about 0 after adjustment, so the confounding reaches the contrast; it just does not reach the top-20. Under the declared gate the verdict was known in advance to be `construction_uninformative`. The amendment, chosen by the user over stopping with validation pending or running as declared:
- the contamination gate reads the anchor only, where no orientation competes;
- a projection gate, |mean projection difference (unadjusted − oracle)| − 2·SE ≥ 0.2, reads every cell, so the residual-confounding criterion is never judged where there is no confounding to remove;
- every criterion and its margin is unchanged, including the contamination criterion at every cell; δ = 2.0 and the `s_z` start of 10 are unchanged;
- the script's rule is version 2, recorded in `settings.json` and in every later calibration iteration.

The verdicts and their consequences are in the spec. In short:
- `adjusted_attribution_validated` lets item 12 start on 0.14.0;
- `diagnose_adjusted_attribution` keeps the leg withheld and calls for a diagnosis proposal;
- `construction_uninformative` keeps the leg withheld and returns the construction to the user.

A verdict below 999 permutations, B = 200 or R = 500 is indicative.

### D11. Calibration, timing and where it runs

Calibration runs R = 20 × 99 permutations × B = 50 at the anchor and orientation e = 1.0, starting from `s_z = 10`, the corrected-permutation power profile's value.
- `s_z` is raised in steps of 5 until the confounding criterion is met, and every iteration goes to `calibration.json`. Under rule version 2 that criterion is the anchor's contamination gap of at least 0.25 and an absolute projection gap of at least 0.2 in both cells (means, no SE); under the declared rule it was the contamination gap of at least 0.25 in both cells.
- `s_z` enters every arm's pooled scale identically, so the arms stay paired at any value. An injected feature's clean signal shrinks in standardized units as `s_z` grows, in all three arms alike.
- If the oracle recoverability criterion fails at a strength that meets the confounding criterion, no further `s_z` is tried. The check stops, the user decides, and the decision is recorded here.

**Calibration under the declared rule, 2026-11-02.** Iterations 0–8 of `results/adjusted-attribution-2026-11-02/calibration/calibration.json`, at `s_z` 10, 15, …, 50, read:

| `s_z` | contamination gap, anchor | contamination gap, e = 1.0 | oracle truth precision, e = 1.0 |
|---|---|---|---|
| 10 | 0.485 | 0.000 | 0.975 |
| 15 | 0.4225 | 0.000 | 0.975 |
| 20 | 0.3525 | 0.000 | 0.970 |
| 25 | 0.280 | 0.000 | 0.970 |
| 30 | 0.2675 | 0.000 | 0.9725 |
| 35 | 0.2375 | 0.000 | 0.975 |
| 40 | 0.205 | 0.000 | 0.9775 |
| 45 | 0.1825 | 0.000 | 0.985 |
| 50 | 0.165 | 0.000 | 0.985 |

The outcome was `stop_exhausted`, and the risk below was wrong: raising `s_z` lowered the anchor's gap and never opened one at e = 1.0, because in standardized units an injected feature converges to `z`'s own cell pattern as `s_z` grows. Iteration 0 ran with BLAS unpinned (12 workers, load 112 on 16 cores, 7.7 min); it was stopped and resumed at `s_z = 15` with one BLAS thread per worker (75 s per iteration). The user amended the gates (D10). Calibration then reruns once under rule version 2 from `s_z = 10`.

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
- **[At e = 1.0 a strong orientation crowds the injected features out of the unadjusted top-k, so the contamination gate fails at that cell only.]** → It happened, and raising `s_z` did not restore the gap (D11, calibration under the declared rule). The gates were amended on 2026-11-02 (D10): contamination is gated at the anchor, and a projection gate, which does not depend on the top-k, reads every cell.
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

## Smoke and run plan (2026-11-02)

**Calibration.** `results/adjusted-attribution-2026-11-02/calibration/calibration.json` holds ten iterations. Iterations 0–8 (declared rule, `s_z` 10–50) ended `stop_exhausted` (D11). After the amendment (D10), iteration 9 (rule version 2, `s_z = 10`, R = 20 × 99 × B = 50 at the anchor and e = 1.0) met every criterion and was adopted: anchor contamination gap 0.485, projection gap 0.563 at the anchor and 0.547 at e = 1.0, oracle truth precision at e = 1.0 0.975. Its readings equal iteration 0's, so pinning BLAS to one thread changed only the speed. Selected ranks: 2 (14 units), 3 (19), 4 (3), 5 (3), 6 (1).

**Smoke.** `results/adjusted-attribution-2026-11-02/smoke/`: R = 20 × 99 permutations × B = 50 over the three cells at `s_z = 10`, 60 units, 14 workers, 1 min 43 s. Zero consistency failures and zero rank-deficient bootstrap replicates. The decision is indicative, as it must be at this size, and reads `diagnose_adjusted_attribution` on one criterion:

| reading (observed principal block) | anchor | e = 0.5 | e = 1.0 | bound |
|---|---|---|---|---|
| gate: contamination gap − 2·SE | 0.374 | — | — | ≥ 0.25 |
| gate: \|projection gap\| − 2·SE | 0.499 | 0.503 | 0.472 | ≥ 0.2 |
| gate: oracle truth precision | — | — | 0.975 | ≥ 0.5 |
| contamination, adjusted − oracle + 2·SE | 0.045 | 0.002 | 0.008 | ≤ 0.05 |
| residual confounding, \|adjusted − oracle\| + 2·SE | 0.036 | 0.018 | 0.023 | ≤ 0.10 × unadjusted (0.056 / 0.057 / 0.055) |
| truth precision, adjusted − oracle − 2·SE | — | **−0.110** | −0.010 | ≥ −0.05 |
| sign stability, adjusted − oracle − 2·SE | — | −0.013 | −0.024 | ≥ −0.05 |
| top-k frequency, adjusted − oracle − 2·SE | — | −0.073 | −0.055 | ≥ −0.10 |

The failing reading has mean +0.010 (the adjusted arm at or above the oracle) and a paired SE of 0.060 at R = 20; at R = 500 the same spread gives an SE of about 0.012 and a reading of about −0.014. It is a small-sample reading of an indicative decision, not a finding. The paper-grade run decides.

**Timing.** Two units at 999 permutations × B = 200, run serially with one BLAS thread on this machine (AMD Ryzen AI 7 PRO 350, 8 cores and 16 threads), took 7.4 s and 8.4 s (`timing/timing.json`). The 1,500 paper-grade units are therefore about 3.5 core-hours serial-equivalent. In the smoke, 14 workers on 8 physical cores ran at about 0.58 units per second at 99 permutations, so the paper-grade wall time is estimated at 40–75 minutes.

**Where it runs.** 3.5 core-hours is under D11's 12-core-hour rule, which would have kept the run local. The user offered the server on 2026-11-02 ("you can also use the server"), so it runs on `ing` as the earlier paper-grade checks did, which keeps the workstation free and runs from a clean cluster checkout. The procedure, following `results/adni-design-2026-11-02/run/PROVENANCE.txt`:
1. Commit the calibration, smoke, `scripts/adjusted_attribution_array.sbatch` (modelled on `scripts/corrected_permutation_array.sbatch`; `bash -n` passes) and this plan; push the branch.
2. On `ing`, create a clean worktree `/home1/tgonza/MOTCO-adjusted-attribution` of the branch at that commit, check that its tree hash equals the workstation's, symlink `.venv` from `/home1/tgonza/MOTCO`, `mkdir logs`, and run `tests/test_adjusted_attribution.py` and `tests/test_adjusted_attribution_check.py` there with `PYTHONPATH` on the worktree's `src`.
3. Pre-create `settings.json` (`run --shard 0/50 --limit 0`), run one timing task (`--limit 1`) to catch a stalled node, then submit `--array=0-49` (30 units per task) on `512x1024`, excluding `n11` and any node short of memory, with `--export=ALL,AA_OUT=<worktree>/results/adjusted-attribution-2026-11-02/run,N_SHARDS=50,AA_CALIBRATION=results/adjusted-attribution-2026-11-02/calibration/calibration.json,PYTHONPATH=<worktree>/src`. Watch shard mtimes; resubmit any stalled task.
4. After every task completes `0:0`, run `merge` and `decide` on `ing`, copy back `settings.json`, the decision, the summary and metric CSVs, `merge.json`, `decide.json` and `sacct.txt`, and write `PROVENANCE.txt` (commands, jobs, tree hash, settings hash, versions, file sha256s, core-hours).

**Approval.** Approved by the user on 2026-11-02, on `ing` ("you can also use the server").
