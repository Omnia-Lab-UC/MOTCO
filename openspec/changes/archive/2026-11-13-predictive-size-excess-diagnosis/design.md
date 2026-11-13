# Design

## Context

See proposal.md (Why) for the motivation. The facts below shape the approach.

- **What the pilot read on the anchor** (`results/predictive-subspace-2026-11-12/plasmode/rate_table.csv`; rejections over 500 units at p ≤ 0.05):

  | arm | `delta` | `angle` | `shape` |
  |---|---|---|---|
  | `pls_cv` | 21 | 24 | 13 |
  | `predictive` | 68 | 11 | 38 |
  | `predictive_orthonormal` | 67 | 9 | 38 |
  | `predictive_fitted` | 65 | 10 | 39 |
  | `forced_k_minus_1` | 7 | 10 | 5 |

- **The size rule on this family.** Under `judge_rates` at m = 3 and R = 500 a reading exceeds above 34 rejections, the ceiling is 36, and one exceedance is allowed. Every predictive arm is above the ceiling on `delta` and `shape`.
- **The anchor's units share one latent matrix.** The fit is label-free and the anchor does not move the features, so every anchor unit fits the same model (rank 5, held-out AUROC 0.635556, in-sample stage R² 0.840). The 500 units differ only in the group labels.
- **How the plasmode draws its null.** The group labels are permuted within each stage; blocks and covariates stay on their donors' rows. Across the anchor's units the cell sizes are fixed (F 11 / 10 / 28, M 9 / 10 / 12), the pooled stage means do not move, and each stage keeps its own residual dispersion.
- **How RRPP draws.** `RRPP` permutes the reduced model's residuals over all n rows (`rng.permutation(n)`), after the optional `df` scale. It has no stratified draw. A draw therefore mixes residuals across stages, which is valid when they are exchangeable across stages.
- **Costs in the pilot.** `pls_cv` took 117 s per plasmode unit, almost all of it the cross-validation. Every forced-rank arm took 0.6 s.
- **Precedent.** The case-study diagnosis (`scripts/case_study_diagnosis.py`, 2026-10-30) read its candidate repair, the κ rescaling, inside the driver through its own inference path, with a per-unit identity check against the analysis record. The package gained `residual_correction` only in the change that validated it (roadmap item 8).

## Goals / Non-Goals

**Goals:**
- Decide, under a rule fixed before any unit runs, whether the anchor's excess is repaired by within-stage draws, by out-of-fold coordinates, by both or by neither.
- Keep every reading paired with the pilot's committed counts: the same units, labels and RRPP seeds.
- Record beside the verdict what the next change needs: the plant cells under each repair, the per-stage dispersions, and the covariate-free arms.

**Non-Goals:**
- Validating the within-stage draws as a test. No simulated design and no translation control is read here.
- Resolving rates on the plant cells, which stay at R = 100.
- Changing any committed driver. The parents are imported and bound, never edited.

## Decisions

**D1. The diagnosis lives in the driver; the package is not touched.**
- **Choice.** The within-stage draw loop and the out-of-fold coordinates are functions of `scripts/predictive_size_excess_check.py`. `motco` stays at 0.24.0.
- **Why.** Neither repair has a reading yet. A package option is a contract with an identity fixture, API documentation and a version; it is owed to a repair that has passed a diagnosis, as `residual_correction` was. The out-of-fold arm has no package path under either choice, so the driver needs its own inference path anyway.
- **Alternative: an opt-in `RRPP(strata=)`, byte-identical when unset.** It gives one draw code path and is what a follow-up would add. It is rejected here because it ships a statistical option before anything says it helps, and under `out_of_fold_restores` or `unexplained` it would stay as unvalidated surface on the core primitive.
- **What guards the driver's loop.** With every row in one stratum the loop must reproduce `RRPP`'s null distributions exactly at the same seed (a test), and on every unit the driver's path under unrestricted draws must reproduce the analysis (D4).

**D2. The units are the pilot's plasmode units.**
- **Units.** `predictive_subspace_check.plasmode_settings` and `plasmode_units`, imported unchanged: the anchor at the `paper` stream's replicates 0–499 and the span plant at e ∈ {0.25, 0.5, 1.0} at replicates 0–99, 999 permutations, 800 units. Labels, plant and RRPP seeds come from `plasmode_check.unit_seeds`.
- **Why the same units.** The excess was read on them. On new label draws the reference arms would have to re-establish it at R = 500 before any check could be graded, and the paired `reduced` grade would lose its counterpart.
- **Why the plant cells run.** A repair that restores size by removing all sensitivity is a different finding from one that keeps the gain. The plant cells are the only real-noise power reading that exists, and they need no new export. They are descriptive (D9).

**D3. One fit per anchor unit without a cross-validation; one cross-validation per plant unit.**
- **Anchor.** Every arm runs at `forced_components` equal to the committed fidelity rank (5, from the plasmode's settings). The final fit is the same `fit_plsda_model` call the CV ends with, as the pilot's amendment of 2026-11-12 established for its predictive arms.
- **Plant.** The plant moves the features, so the rank is the unit's: `pls_cv` runs the CV once, and every other arm runs forced at that rank, as in the pilot.
- **Pre-flight.** Before the full run, the cross-validated `pls_cv` on anchor replicates 0 and 1 must select rank 5 at AUROC 0.635556 and give the forced call's statistics and p-values exactly. A mismatch stops the run; the outcome goes in `notes.md` and `PROVENANCE.txt`.
- **Alternative: cross-validate every anchor unit.** About 16 core-hours to recompute one label-free integer that three committed runs pin. The exact reproduction of twelve committed anchor counts is the stronger check and costs nothing.

**D4. One inference path for every arm the analysis cannot run, checked against the analysis on every unit.**
- **The path.** Given a latent matrix in the dataset's row order: the covariate-adjusted full and reduced models, LS-mean rows and contrast from `motco.analysis._trajectory_inputs`; the observed statistics from `estimate_difference` under net displacement; the null from `RRPP` (unrestricted, serial, `residual_correction="df"`) or from the within-stage loop (D5); p-values by the analysis's rule (`_empirical_p_value`); rejection at p ≤ α, the pilot's plasmode rule.
- **Which arms use it.** The within-stage arms (on the latent matrix of the corresponding analysis call) and both out-of-fold arms. The reference and covariate-free arms are plain `analyze_trajectories` calls.
- **Path identity.** For `pls_cv` and the three in-sample predictive arms, the path under unrestricted draws must return the analysis's statistics to 1e-10 and its p-values exactly. A unit failing it is recorded as `inconsistent`, naming the arm and deviation, and leaves every rate.
- **Why private helpers are acceptable.** `case_study_diagnosis.py` does the same, and the identity check turns any drift in them into a gate instead of a silent difference.

**D5. Within-stage draws (check W).**
- **Definition.** Per draw, for each stage in the declared order, the residual rows of that stage are permuted among themselves: `idx[rows] = rows[rng.permutation(rows.size)]`, with `rng = default_rng(the unit's RRPP seed)`. Everything else is `RRPP`'s: the reduced fit, the `df` scale applied once before any draw, `ŷ_red + e*`, `estimate_difference` on each draw.
- **What it changes.** Two things, and only these:
  - each stage's residuals stay in that stage, so a stage with larger dispersion keeps it in every draw;
  - the reduced model contains the stage indicators, so its residuals sum to zero within each stage, and a within-stage draw keeps the pooled stage means at their fitted values. Unrestricted draws let them move.

  Both are true of the plasmode's own label draws. The first is the one the guess names; the dispersion table (D9) says how unequal the stages are in each coordinate system.
- **Why it is still a test of the interaction.** The group main effect stays in the reduced fit, so it is not a label permutation: a pure translation is not rejected by construction. It assumes exchangeability within stage only, which is weaker than `RRPP`'s assumption.
- **Arms.** The three in-sample predictive coordinate systems are decisive. `pls_cv__within_stage` is descriptive: it shows whether production's own rates move.
- **Draws are not shared with the reference.** A restricted draw cannot reuse the unrestricted index sequence. Each within-stage arm of a unit uses the same restricted sequence, so the three predictive arms are paired with one another.

**D6. Out-of-fold fitted stage indicators (check O).**
- **Definition.** The donors are split into four stage-stratified, shuffled folds. For each fold, `fit_plsda_model` is fitted on the other three folds' rows of the analysis's PLS input features at the unit's production rank and scaling, and the held-out rows get `transform(X_held) · y_loadingsᵀ`: the fitted stage indicators in that model's response units, which is how `predictive_fitted` is computed in sample. The result is an n × K matrix measured through D4's path.
- **Why the fitted indicators.** `predictive` and `predictive_orthonormal` are coordinates on a basis that each fold chooses freely within its own subspace, so stacking them needs an alignment step with its own in-sample fit. The response units need none. In sample the three coordinate systems agree to within three rejections per statistic, so one arm reads the subspace.
- **Fold-invariant units.** The response centring and scale depend only on the training stage proportions. Stratification holds them equal across folds up to the rounding of stage counts, and on this export (20 / 20 / 40 donors) they are exactly equal.
- **Four folds.** The production CV's outer split count. Each training set has 60 donors with the full sample's stage proportions, and its smallest stage (15) is well above the rank.
- **A fold seed per replicate.** Derived by sha256 from the driver's name, base seed 2037, the replicate and the purpose `folds`. On the anchor one split would make the reading conditional on one partition of 80 donors; 500 splits average over it. The seed is recorded; the assignment is not.
- **What is out of fold and what is not.** Only the stage-supervised fit. The block scalers stay pooled: they are label-free and are the analysis's own input. The rank is the full-sample CV's.
- **Limit check.** With every donor in the training set the builder must return the `predictive_fitted` matrix to 1e-8 (a test).
- **Arms.** `fitted_out_of_fold` under unrestricted draws is decisive. `fitted_out_of_fold__within_stage` is descriptive and completes the two-by-two.

**D7. Fourteen arms per unit.**

| arm | coordinates | draws | covariates | role |
|---|---|---|---|---|
| `pls_cv` | scores | unrestricted | six | reference, reproduction |
| `predictive`, `predictive_orthonormal`, `predictive_fitted` | in-sample predictive | unrestricted | six | reference, reproduction |
| `predictive__within_stage`, `predictive_orthonormal__within_stage`, `predictive_fitted__within_stage` | in-sample predictive | within stage | six | check W, decisive |
| `fitted_out_of_fold` | out-of-fold fitted indicators | unrestricted | six | check O, decisive |
| `pls_cv__within_stage` | scores | within stage | six | descriptive |
| `fitted_out_of_fold__within_stage` | out-of-fold fitted indicators | within stage | six | descriptive |
| `pls_cv__no_covariates`, `predictive__no_covariates`, `predictive_orthonormal__no_covariates`, `predictive_fitted__no_covariates` | as named | unrestricted | none | descriptive |

- **The covariate-free arms** drop the six covariates from the trajectory model and keep the `df` correction, whose κ then follows the four-column reduced model. The latent fit is covariate-blind, so the coordinates are the reference arms'.
- **`forced_k_minus_1` is not rerun.** Its reading is committed and no check bears on it.

**D8. The decision.**
- **Gate.** `reproduction_failed` when any reference arm's count in any cell differs from the pilot's committed `rate_table.csv`, when any unit is inconsistent, or when the pre-flight was not passed. A failed unit already breaks the exact counts. The grades below are still recorded.
- **Grades, per decisive arm, on the anchor.**
  - `restored`: the arm's `delta`, `angle` and `shape` rejections hold `judge_rates` at m = 3. This is the rule and the family the pilot's S1 failed, so a `restored` arm is one the pilot would not have stopped for size on this anchor.
  - `reduced`: not `restored`, and the one-sided exact paired p-value for fewer rejections than the arm's counterpart is at most 0.05/2 on `delta` or on `shape`, the two statistics that exceeded. The counterpart is the same coordinate system under unrestricted draws; for `fitted_out_of_fold` it is `predictive_fitted`. With b the units only the counterpart rejects and c the units only the arm rejects, p = P(Bin(b + c, ½) ≥ b), and 1 when b + c = 0.
  - `unchanged`: neither.
- **Check outcomes.** W takes the weakest grade of its three arms, in the order `unchanged`, `reduced`, `restored`: the excess is the subspace's, so a repair that works in one metric only has not repaired it. O takes its arm's grade.
- **Verdict.**

  | W | O | verdict |
  |---|---|---|
  | `restored` | `restored` | `both_restore` |
  | `restored` | other | `null_restores` |
  | other | `restored` | `out_of_fold_restores` |
  | neither `restored`, at least one `reduced` | | `partly_explained` |
  | `unchanged` | `unchanged` | `unexplained` |

- **Consequences, fixed with the rule.**
  - `null_restores`: the unrestricted draws are what fails in the predictive subspace on this export. A follow-up change may propose within-stage draws as an opt-in option and must read their size and power on the pilot's simulated designs, the weak design and this anchor before `predictive_orthonormal` is re-read.
  - `both_restore`: as above, and the out-of-fold reading confirms that the in-sample fit is where the excess lives.
  - `out_of_fold_restores`: the excess lives in the in-sample fit and a restricted null does not repair it. A repair needs a cross-fitted measurement, which belongs to open item 2's held-out half and open item 5's cross-fitted arm. The predictive values stay diagnostics at SEA-AD's design.
  - `partly_explained`: the guess is incomplete. The report names which check reduced and by how much; the remainder waits for open item 2's cohorts.
  - `unexplained`: the guess is rejected and the roadmap's prediction for ADNI is withdrawn. The covariate-free arms and the dispersion table are the leads.
- **Why the plant cells do not decide.** The question is size. The plant is a magnitude edit along the in-sample stage means, which favours in-sample coordinates, and R = 100.
- **Unavailable readings.** A decisive arm with no available anchor reading for a statistic is graded `unchanged` and named.

**D9. Descriptive readings.**
- **Plant cells.** Every arm's rates at e ∈ {0.25, 0.5, 1.0}; each arm's exact paired gain against `pls_cv` and paired loss against its counterpart. The report states, for each arm graded `restored`, whether `delta` at e = 0.5 stays above production's 0.51.
- **Production under within-stage draws**, and the out-of-fold coordinates under them.
- **Covariate-free arms:** rates and paired tests against their six-covariate counterparts.
- **Dispersion.** Per unit and coordinate system (scores, the three predictive systems, out of fold): the mean squared norm of the reduced model's residuals in each stage, divided by the pooled value, and the largest-to-smallest ratio. Stages are named by their position in the declared order. Beside it, the in-sample stage R² and the share of the stage indicators' variance the out-of-fold predictions explain. The merge reports medians per cell.
- **Null width.** The median `angle` null 95th percentile per arm, as in the pilot.

**D10. Execution follows the recent validation drivers.**
- `run` / `merge` / `decide`; `settings.json` bound by hash, carrying the parents' sha256 (the pilot's plasmode `settings.json` and `rate_table.csv`, and through `plasmode_settings` the plasmode's and the production selection's files); a resume or merge under other settings exits naming the keys; a changed parent is refused by name.
- Round-robin shards, resume by complete unit key, `--workers` with BLAS pinned to one thread, JSONL records holding seeds and statistics only.
- The sbatch array has a dry-run variable that prints the command. The run fits on the laptop, so the array is the fallback for the plant cells.
- Committed: `settings.json`, the rate, paired and cell tables, `accounting.json`, `reproduction.json`, the decision and `PROVENANCE.txt` under `results/predictive-size-excess-<date>/`. Records stay off the repository.

## Risks / Trade-offs

- **[The forced-rank anchor arms are not the cross-validated ones]** → The pre-flight compares them on two units, the gate requires all twelve committed anchor counts, and the pilot's `same_fit` check runs on every unit.
- **[The driver's inference path drifts from the analysis]** → Path identity on every unit for four arms; an inconsistent unit gates the verdict.
- **[Within-stage draws hold size here because they mirror how the anchor's labels were drawn]** → That is the reading: the check asks whether the extra assumption of the unrestricted draws is what fails. Whether within-stage draws hold size where residuals are exchangeable across stages, and what they cost in power, is the follow-up's simulation reading, and the verdict's consequence says so.
- **[A repair restores size by losing sensitivity]** → The plant cells are reported for every arm. The out-of-fold arm is expected to lose most of it: the held-out AUROC is 0.636 and the plant follows the in-sample means, so its plant reading is a lower bound.
- **[Both checks restore and the cause is not separated]** → `both_restore` states that either suffices. The dispersion table and the descriptive two-by-two arm say which mechanism is larger; the follow-up starts from the draws because they keep the coordinates and the pilot's gains.
- **[The reading is conditional on one fit of 80 donors]** → As the pilot's was. The out-of-fold arm averages over fold splits, not over donors. Another cohort is open item 2.
- **[The rank is still chosen on all donors]** → It is one integer selected by held-out AUROC, and the anchor's is pinned at 5.
- **[Reproduction fails from environment drift]** → The gate stops the reading. A mismatch is diagnosed, never waived.
- **[Fourteen arms invite reading the descriptive ones as tests]** → Only four arms are graded, the rule names them, and the report separates the two tables.

## Migration Plan

Nothing migrates: no package code, configuration, fixture or committed record changes.
- **On `null_restores` or `both_restore`:** a separate change proposes the within-stage draw as an opt-in option of `RRPP`, the harness and `motco analyze`, byte-identical when unset, and validates it. Only its verdict could reopen `predictive_orthonormal`, under the admissibility rule the roadmap requires.
- **On `out_of_fold_restores`:** no package change follows. Open item 2's held-out half carries the cross-fitted measurement.
- **On `partly_explained` or `unexplained`:** the roadmap records the remaining candidates and open item 2's cohort readings go first.
- **Rollback:** remove the driver, its test and its array script. Nothing imports them.

## Open Questions

- Whether the full run goes on the laptop or on `ing`, and its core-hour cap, is set from the smoke timing at the approval gate.
