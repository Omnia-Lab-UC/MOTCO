# Design

## Context

See proposal.md (Why) for the motivation and specs/case-study-configuration-diagnosis/spec.md for the requirements. The pieces this design builds on, all in 0.12.0:

- **`scripts/case_study_type_i_check.py` (the parent check).**
  - `build_settings`, `build_unit` (the D10 replicate plus five noise covariates from their own seed, drawn as one `(n, 5)` matrix), `settings_hash`, the shard and resume helpers, and the bound `anchor_bound`.
  - The committed runs are `results/case-study-type-i-2026-10-29/delta0.1/` and `delta0.05/`. Each holds `settings.json` (hash `899ce29a…` for δ = 0.1), `case_study_type_i_pvalues.csv` (`mode, effect, replicate, arm, statistic, p_value`, 9,000 rows; 999 permutations, so every p-value is k/1000) and `case_study_type_i_replicates.csv` (rank and AUROC per unit × arm).
- **`scripts/covariate_type_i_check.py` (D10).** `build_replicate` converts methylation to M-values, then adds `s_z·outer(z, w)·sd` per block, where `sd` is the pre-injection `ddof=1` SD and `w` is drawn from the base seed. At `s_z = 0` the returned blocks are the pre-injection values, with the same seeds and the same `z`.
- **`motco.analysis.analyze_trajectories`.**
  - It standardizes with `preprocessing.fit_block_scalers` (pooled mean, `ddof=0` SD) and fits the covariate-blind basis with `scale=True`. In sklearn that means centring by `x_mean_` and dividing by the `ddof=1` SD before `x_rotations_`.
  - It returns `latent_scores` and `pls_loadings` (the `x_rotations_`), but not the estimator or the scalers.
  - Covariates enter through `_covariate_adjusted_design`, which appends them to both models and the pooled-mean reference row to every LS-mean row.
  - Inference is `estimate_difference` plus `RRPP`, with p = (1 + #{null ≥ observed}) / (1 + B).
  - `forced_components` works with a nuisance design (`tests/test_latent_adjustment.py::test_forced_components_are_honoured_with_nuisance`).
- **`motco.stats.permutation.RRPP`.** It refits the reduced model on its `Y` and draws `ŷ_red + e_red[perm]`, with `perm` from `default_rng(seed).permutation(n)`. The identity is not among the draws. One 999-permutation run on an 80 × 2 latent with ten design columns takes 0.19 s on the workstation (measured 2026-10-30).
- **Reference-row invariance** is already tested for the statistics (`tests/test_analysis.py::test_statistics_do_not_depend_on_the_reference_row`). RRPP evaluates the same LS-mean rows in every draw, so the p-values inherit it.
- **The parent run's realized values.**
  - The `adjusted` arm selects rank 2 in 366 / 500 (δ = 0.1) and 407 / 500 (δ = 0.05) anchor units.
  - `latent_stage` selects rank ≥ 4 in 91 and 61, with a maximum of 11. Its `cv_min_training_df` is 31, and `residual_df` is 71 for the full sample.
  - One parent unit (three analyses) took 47–55 s on an `ing` EPYC core.

## Goals / Non-Goals

**Goals:**

- Separate the four candidate causes of the `adjusted` excess on the same units, with the basis held fixed: column count, collinear estimation, residual shrinkage, and Monte Carlo chance. The diagnosis reads paired contrasts, not differences between independent runs.
- Separate the `latent_stage` excess's rank selection from its residualized basis.
- Make every derived arm trustworthy by tying it to the production path: the parent's p-values reproduced, the script's inference path identical to `analyze_trajectories`, and the oracle reconstruction verified on every unit.

**Non-Goals:**

- Adopting a remedy. The corrected permutation is computed in the script, and adopting it is a later MOTCO change with its own Type I and power check.
- Separating the per-training-set nuisance coefficients from the rank tail inside `latent_stage`. That would need a core switch; see Risks.
- New data. The diagnosis adds no replicates, so it cannot say whether the parent's two-rejection excess was chance. The R = 1,000 replication answers that and stays deferred.

## Decisions

### D1. A new script bound to the parent run, not flags on the parent check

`scripts/case_study_diagnosis.py` imports `case_study_type_i_check as cs` (and through it D10). Each run is bound to one parent directory through `--parent`:
- It reads `settings.json` and checks the recorded hash with `cs.settings_hash`.
- It rebuilds the settings with `cs.build_settings` from the recorded arguments and refuses a different hash. That catches construction code that has drifted since the parent run.
- It writes its own `settings.json`, which carries the parent's path and hash, so a resume under another parent is refused (spec, first requirement).
- `--reps` may take a prefix of the parent's replicates. `--perms` below the parent's marks the run as a smoke run: reproduction is `not_compared` and the verdicts are indicative.

*Alternative: arms and flags added to `case_study_type_i_check.py`.* Rejected. Its settings hash identifies the committed runs, its tests are keyed on three arms, and its decision vocabulary is a validation's, not a diagnosis's. It is the same reasoning as the parent check's D1.

*Alternative: re-specifying the construction by flags (`--delta`, `--s-z`).* Rejected. Binding to the parent's directory makes a mismatched construction impossible to request, and leaves the reproduction check to catch code drift.

### D2. One CV basis fit per basis; trajectory-model arms reuse it through the script's inference path

Per unit the script calls `analyze_trajectories` once with the six covariates, giving `z_noise5`. That call yields the record (the reproduction target), `latent_scores` and `pls_loadings`. The other factorial arms are inference-only on those scores:
- the analysis's design, built from the recoded working frame as `analyze_trajectories` builds it;
- `TrajectoryDataset(...).covariate_design` for each covariate subset;
- `motco.analysis._covariate_adjusted_design`;
- `estimate_difference` and serial `RRPP` with the unit's RRPP seed.

On every unit, the script's own `z_noise5` must equal the analysis record's statistics and p-values exactly, or the unit fails with an error. That assertion makes the private helper and the recode safe to use, and it pins the inference path every other factorial arm uses.

*Alternative: `analyze_trajectories` for every arm.* Rejected. It refits the CV basis six times (about 100 s more per unit), cannot express the oracle or the correction, and the basis-identity assertion would have to be repeated per arm.

*Alternative: a public inference hook in `motco.analysis`.* Rejected for now. It is a package change and a release for a diagnosis, and the follow-up change that adopts a correction is the place to add a public inference surface.

### D3. The oracle subtracts the injection inside the fixed basis

The oracle's question is: with the basis exactly as the analysis fitted it, what is the size if `z`'s effect is known rather than estimated? The script works it out in three steps:
1. It builds the pre-injection blocks `X0` with D10's builder at `s_z = 0` (same seeds), and asserts that `X − X0` is `outer(z, c_b)` per block for a fixed vector `c_b`.
2. It reconstructs the fitted projection: `fit_block_scalers` and `transform_blocks` on `X`, which are the analysis's own calls, then the PLS centring and `ddof=1` scaling of the standardized matrix, then `pls_loadings` as the rotations. It asserts the reconstruction reproduces `latent_scores` to 1e-9.
3. It sets `Y0` to the projection of `X0` through that reconstructed map.

`Y − Y0` is then exactly `z` times a fixed latent vector, so the oracle arms run the factorial on `Y0` without `z`. With `z` in the model, `Y` and `Y0` give the same inference, because OLS removes any multiple of a model column. That is why the `z_noise*` arms run on `Y`.

*Alternative: analysing `X0` from scratch.* Rejected. It fits a different basis, one that does not follow `z`, and so answers the D10 / rank-floor question instead of this one.

*Alternative: subtracting `z·ĉ` with `ĉ` estimated.* Rejected. That is the adjusted model, not an oracle.

*Fallback if the reconstruction does not match* (zero-variance features, a convention change): call `fit_pls_latent_space` directly with the analysis's inputs and use the estimator's `transform`, asserting its scores equal `latent_scores`. That costs one more CV fit per unit and is used only if smoke shows the reconstruction failing.

### D4. Factorial: `z` estimated or oracle × 0 / 2 / 5 noise columns

| arm | covariates in both models | p_red |
|---|---|---|
| `oracle_noise0` | none (on `Y0`) | 4 |
| `oracle_noise2` | x1, x2 (on `Y0`) | 6 |
| `oracle_noise5` | x1–x5 (on `Y0`) | 9 |
| `z_noise0` | z | 5 |
| `z_noise2` | z, x1, x2 | 7 |
| `z_noise5` | z, x1–x5 (the case study's model) | 10 |

p_red counts the intercept, one group, two stage and the covariate columns. The noise columns are prefixes of the parent's `(n, 5)` draw, so the arms are nested. The rows of the table answer three questions:
- **Oracle rows:** the pure column cost, with no collinearity.
- **`z` rows:** the same column cost with one column collinear with the group × stage cells.
- **`z_noise0` vs `oracle_noise0`:** the cost of estimating `z`.

Report §4 asked for "one, three and six" columns. The `z` rows are exactly that, and the oracle rows add the matching counts without `z`.

### D5. The correction: Freedman–Lane on rescaled residuals, through the public `RRPP`

For an arm with reduced model `M_r` (p_red columns), the script computes `ŷ_red` and `e_red` on the arm's `Y` and calls `RRPP` on `Y' = ŷ_red + κ·e_red`, with κ = √(n / (n − p_red)), the arm's full model, LS-means and the unit's RRPP seed. `κ·e_red` is orthogonal to `M_r`'s column space, so `RRPP` recovers `ŷ_red` and draws `ŷ_red + κ·e_red[perm]` with the same `perm` as the uncorrected arm. The observed statistic comes from the uncorrected `Y`, and the p-value uses the analysis's rule.

Why uniform κ: the mean square of the permuted residuals is RSS/n, and κ restores the unbiased σ² = RSS/(n − p_red). That is the minimal test of the shrinkage mechanism. κ is 1.026 at p_red = 4 and 1.069 at p_red = 10.

The mechanism's prediction, which is what the verdict reads:
- The uncorrected `delta` rejection count rises with p_red in both families.
- The corrected count is flat across p_red.
- `angle` is little affected in either.

The rough normal-tail figure in the proposal (0.056 → 0.067) overshoots the rank-floor check's no-covariate anchor (`delta` 0.042–0.044 measured), so the baseline may be conservative for other reasons. The verdict therefore reads contrasts and the corrected arm's bound, not absolute predictions.

*Alternative: leverage-adjusted residuals `e_i / √(1 − h_ii)`.* Deferred to the follow-up change. They are non-uniform and change which units dominate the null. If the uniform correction explains the excess, the follow-up compares both.

*Alternative: a studentized `delta`.* Rejected here. It changes the statistic rather than diagnosing the current test.

### D6. Latent arms: CV at floor 2, CV at floor 1, forced rank 2

Three `analyze_trajectories` calls per unit with the six covariates and `covariate_latent_adjustment="stage"`: `min_components=2` (`latent_stage_cv`, the reproduction target), `min_components=1` (`latent_stage_floor1`) and `forced_components=2` (`latent_stage_rank2`, no CV). Forced rank 2 matches the covariate-blind basis's modal rank and removes the CV entirely, taking the rank tail and the per-training-set coefficients with it. Floor 1's rank-1 units have a degenerate `angle` (0° or 180°, a two-point null), as in the rank-floor check. They are counted, and the floor-1 arm is descriptive (spec).

*Alternative: forced rank 3, or a ladder of forced ranks.* Rejected for cost. Every forced rank needs its own basis fit, though without CV, and rank 2 is the configuration that the primary actually measures.

### D7. Contrasts: a replicate-clustered sign-flip test

The two controls of one replicate share the baseline draw, `z`, the noise covariates and the RRPP seed. The two runs share the same seeds at different δ. So the 2,000 paired readings form 500 replicate clusters. For contrast B vs A, each replicate r gets `d_r` = (readings where B rejects and A does not) − (readings where A rejects and B does not), with d_r ∈ {−4, …, 4}. The one-sided p-value is P(Σ ε_r d_r ≥ Σ d_r) over independent equiprobable signs ε_r, computed exactly by convolution over the integer sums. The table also gives the discordant counts and, descriptively, the per-reading (run × control) contrasts.

Predeclared contrasts, all on `delta`:
- C1: `oracle_noise5` vs `oracle_noise0`
- C2: `z_noise5` vs `z_noise0`
- C3: `z_noise0` vs `oracle_noise0`
- C4: `z_noise5` corrected vs `z_noise0` corrected

The same contrasts on `shape` and `angle` are written as descriptive.

*Alternative: McNemar over the 2,000 readings.* Rejected. It treats correlated readings as independent and overstates significance.

### D8. Verdicts and their order

The verdicts are those of the spec: `residual_df_shrinkage` → `collinear_covariate` → `unresolved` for the `adjusted` `delta`, and `rank_selection` | `residualized_basis` for the `latent_stage` arm, behind the `reproduction_failed` gate. The verdict function is pure over the merged p-value rows, the rank rows and the reproduction counts, and is tested on synthetic rows for every outcome.

The order means a correction that restores size wins even when C3 is also significant: in that case the correction covers `z`'s column too. The decision file records every contrast either way, so the report can name joint mechanisms.

Each verdict's consequence is recorded with it:
- `residual_df_shrinkage`: a MOTCO change offering the corrected permutation (opt-in, byte-identical when absent, its own Type I and power check). After that, a `sead_sexdiff` re-issue.
- `collinear_covariate`: the case study keeps reading `delta` against the realized size, and the next diagnostic is the confounding strength (`a_gs`) at one column.
- `unresolved`: the R = 1,000 replication.
- `rank_selection`: the `stage` arm is admissible at a fixed rank of 2 at this configuration, which `sead_sexdiff` may adopt in its own change.
- `residualized_basis`: the `stage` arm stays inadmissible here.

### D9. Records, layout and subcommands

`--out-dir` holds one subdirectory per parent run (`delta0.1/`, `delta0.05/`), each with `settings.json`, `records/shard_<i>_of_<N>.jsonl` and, after `merge`, its per-run tables.

Each unit is one JSONL line with `key`, `seeds`, the oracle and reconstruction checks (max abs reconstruction error, the rank-one injection residual) and `arms`. Per arm it records `statistics`, `p_values`, `null_q95`, `selected_lv`, `component_selection`, `cv_mean_auroc`, `p_red`, `kappa` and `error`.

The subcommands:
- `run` (sharded, resumable, as in the parent);
- `merge` per run, writing `diagnosis_pvalues.csv`, `diagnosis_rates.csv`, `diagnosis_rank_table.csv` and `diagnosis_reproduction.json`;
- `decide` across runs, writing `diagnosis_contrasts.csv` and `case_study_diagnosis_decision.json`.

### D10. Run, provenance and no release

1. Smoke locally:
   - R = 20 at 99 permutations on both parents, exercising every arm and the indicative verdicts;
   - R = 4 at 999 permutations, eight units per parent, which exercises reproduction against the committed p-values.
2. Time the units.
3. Approval gate: the user approves the arrays.
4. One 50-shard array per parent on `ing`, from a clean worktree of the pushed branch (`scripts/case_study_diagnosis_array.sbatch`, modelled on the parent's template, excluding n11), then `merge`, `decide` and `PROVENANCE.txt` (command, jobs, tree hash, parent hashes, file sha256s, core-hours).

`records/` stays on `ing`. The CSVs are committed if each is under about 1 MB; that is about 50,000 p-value rows per run, so the p-value table may have to stay on `ing` and be named in `PROVENANCE.txt`.

No change under `src/motco/`, so the version stays 0.12.0.

## Risks / Trade-offs

- [The reconstruction of the fitted projection misses, e.g. through a constant feature whose scale the scalers set to 1 while sklearn's convention differs] → The unit fails loudly at smoke, and D3's fallback (a direct `fit_pls_latent_space` call with the estimator's `transform`) is taken, costing one more CV fit per unit.
- [The script's inference path drifts from `analyze_trajectories`] → The per-unit `z_noise5` identity assertion makes any drift a hard error, not a silent bias.
- [The correction over-corrects and the corrected arms sit well below α] → Reported as a power cost. The verdict only requires the corrected case-study model to hold and to show no column trend, and power belongs to the follow-up change.
- [The contrasts are underpowered for a 0.01 rate difference] → Pairing on identical observed statistics (corrected vs uncorrected) and on one basis keeps discordances few and directional under the mechanism. If they are not significant, the verdict is `unresolved` and names the R = 1,000 replication, as the spec requires.
- [`latent_stage_rank2` confounds the rank tail with the per-training-set coefficients] → Stated in the report. If `rank_selection` comes back, separating them needs a core switch (full-sample coefficients inside the CV), which is out of scope here.
- [The verdict order hides a joint mechanism] → Every contrast is in the decision file, and the report reads them, not just the label.
- [Cost: three CV fits, one forced fit and about 13 RRPP calls per unit, planned at 55–70 s on `ing`, 35–40 core-hours] → Re-timed at smoke before the approval gate. Shards are resumable, and n11 is excluded.

## Migration Plan

No package or API change, so there is nothing to migrate or roll back. The run sequence is smoke → timing → approval → arrays on `ing` → merge and decide → `PROVENANCE.txt` → report → roadmap and docs. The results go to `results/case-study-diagnosis-<date>/` with one subdirectory per parent run.

## Open Questions

- Whether `diagnosis_pvalues.csv` is committed or kept on `ing` is decided by its size at merge (D10). Neither choice changes the decision file.

## Smoke and shard plan (2026-10-30)

Run locally on the workstation (24 cores, single-threaded BLAS) from the working tree of branch `feat/case-study-configuration-diagnosis` before its first commit, with `scripts/case_study_diagnosis.py` as committed with this change.

**Smoke at R = 20 × 99 permutations** (`results/case-study-diagnosis-2026-10-30/smoke-delta0.1/`, `smoke-delta0.05/`, decision in `smoke-decision/`; 10 shards per parent):
- 40 / 40 units per parent, 0 arm failures, 0 identity or reconstruction failures. The maximum reconstruction error is 7.1e-14 (tolerance 1e-9), and the maximum rank-one residual of `X − X0` and `Y − Y0` is 4.9e-16.
- Reproduction `not_compared` in both runs, and the verdicts are marked indicative: `unresolved` (C1–C4 one-sided p = 0.50, 0.88, 0.88, 0.94) and `residualized_basis` (`latent_stage_rank2` anchor `shape` 6 / 20 in both runs).
- At R = 20 the bound is 0.147, so 3 / 20 already exceeds. The smoke's judged arms match the parent's own 999-permutation p-values on the same replicates: `z_noise5` anchor `delta` 2 / 20 (parent `adjusted` 2 / 20) and `latent_stage_cv` anchor `shape` 7 / 20 (parent 7 / 20, against 46 / 500 over the full run). Replicates 0–19 run hot, so the smoke rates are not a size reading.
- The oracle arms read anchor `delta` 5–7 / 20 against 2–6 / 20 for the `z` arms. Two probes on 10–20 anchor units found no sign that `Y0` carries `z` through the basis fitted on `X`. Its within-cell multiple correlation with `z` (0.05–0.38) and the share of its group × stage interaction pattern along `z`'s (mean 0.22) are at or below those of `X0` projected through its own basis (0.07–0.34; 0.31). This is read as Monte Carlo noise at R = 20. The paper-grade report reads the oracle rows' absolute rates beside the contrasts, not as an α baseline.
- `latent_stage_floor1` selects rank 1 in 4 / 20 (δ = 0.1) and 6 / 20 (δ = 0.05) anchor units.

**Reproduction at R = 4 × 999 permutations** (scratch, `build/case-study-diagnosis-2026-10-30/`): 8 / 8 units per parent and status `holds`, with 0 mismatches over 24 `z_noise5` and 24 `latent_stage_cv` p-values per parent. One unit run through the sbatch template's command line (`--shard 0/2 --limit 1`) also reproduced its six committed p-values.

**Timing.** A diagnosis unit at 999 permutations takes 20.5 s on one workstation core against 18.5 s for the parent's unit, a ratio of 1.11. The parent measured 54.7 s per unit on `ing` (15.2 core-hours per 1,000 units at δ = 0.1), so the plan is about 61 s per unit on an `ing` EPYC core.

**Shard plan.**
- One 50-shard array per parent (20 units per shard, about 20 min per shard on an uncontended core) on partition `512x1024` with `--exclude=n11`.
- Code: a clean worktree of the pushed branch, with `PYTHONPATH` on its `src` and `.venv` symlinked, as for the parent. `settings.json` is pre-created with `run --shard 0/50 --limit 0`.
- Outputs go to `results/case-study-diagnosis-2026-10-30/delta0.1/` and `delta0.05/`.
- About 17 core-hours per parent and 34 for both. The 8 h task limit covers a node running 5× slow, as n3 and n4 did during the parent's δ = 0.05 run. Stalled shards are cancelled and resubmitted, excluding the slow node; the shards are resumable.
- Committed afterwards: `settings.json`, `merge.json`, `diagnosis_rates.csv`, `diagnosis_rank_table.csv`, `diagnosis_reproduction.json`, `PROVENANCE.txt` per run, plus `case_study_diagnosis_decision.json` and `diagnosis_contrasts.csv`. `diagnosis_pvalues.csv` (about 45,000 rows, about 2.2 MB per run), `records/` and `merged.jsonl` stay on `ing`, named with their sha256 in `PROVENANCE.txt`.

**Approval:** the user approved both arrays as planned on 2026-10-30.
