# Design

## Context

See proposal.md (Why) for the motivation. The requirements are in the specs: `rrpp-residual-correction`, `real-data-trajectory-analysis`, `module-contracts` and `corrected-permutation-validation`. This change builds on these pieces of 0.12.0:

- **`motco.stats.permutation.RRPP`.**
  - It fits the reduced model on `Y` with `estimate_betas` (Cholesky, then a direct solve, then the lstsq fallback), and forms `y_hat` and `y_res`.
  - It draws `y_hat + y_res[perm]`, where `perm` comes from `default_rng(seed).permutation(n)` (serial path) or from per-worker `SeedSequence` children (parallel path, `_RRPPWorker`).
  - Both paths receive `y_hat_np` and `y_res_np` as NumPy arrays computed once before any draw.
  - It returns null distributions only; the observed statistic and the p-value are the caller's.
- **`motco.analysis.analyze_trajectories`.** It calls `RRPP` once on `latent.matrix` with `_covariate_adjusted_design`'s full and reduced models, the LS-means, the contrast and `params.seed`.
  - The record's `measurement` block always exists and holds the RRPP settings (`permutations`, `seed`, `n_jobs`).
  - `TrajectoryAnalysisParams.to_dict()` drops opt-in fields while they are `None`; `covariate_latent_adjustment` and `min_components` already work this way.
  - `tests/test_analysis_identity.py` pins the record (minus `runtime_seconds` and `software`) and the CSVs against two committed fixtures.
- **The diagnosis's correction (its D5).**
  - `scripts/case_study_diagnosis.py` ran `RRPP` on `ŷ_red + κ·e_red`, using `rescaled_residuals` and `residual_df_kappa`.
  - `κ·e_red` is orthogonal to the reduced model, so `RRPP`'s own refit recovers `ŷ_red` and draws `ŷ_red + κ·e_red[perm]`.
  - That route cannot express a non-uniform rescaling: `e_i / √(1 − h_ii)` is not orthogonal to the reduced model, so the refit would absorb part of it.
- **The parent construction.**
  - `scripts/case_study_type_i_check.py` provides `build_settings`, `build_unit`, `generator_params`, `check_headroom`, `calibration_criteria`, the shard/resume helpers and `anchor_bound`.
  - Seeds come from `covariate_type_i_check.replicate_seeds(base_seed, replicate)` and `_seed(base_seed, "noise_covariates", replicate)`. So a new base seed gives a new family, including the injection directions `w`.
  - `build_settings` takes the cell table from `PHASE6_GENERATOR` and the stage effect from `--delta`.
- **Committed base seeds.** 100, 400, 500, 600, 800, 900, 1028 and 2026. The new family uses **2031**.
- **The cohort tables** (read 2026-10-30 from the `cohort_30` `result.json` files on `ing`, MOTCO 0.10.0, with six covariate columns each):
  - astrocyte `[[11, 10, 28], [9, 10, 12]]`, n = 80, sha256 `34a744c2…`;
  - microglia `[[11, 10, 24], [9, 11, 11]]`, n = 76, sha256 `0dba07a7…`.

  With intercept, group, two stage columns and six covariates, p_red = 10. So κ = √(80/70) = 1.069 and √(76/66) = 1.073.
- **Timing.**
  - Parent check: one unit of three analyses took 47–55 s on an `ing` EPYC core, i.e. about 17 s per covariate-adjusted analysis at the weak signal.
  - Latent-adjustment power comparison at δ = 2.0: about 13.5 s per analysis.
  - 999 RRPP permutations on an 80 × 2 latent take about 0.2 s.

## Goals / Non-Goals

**Goals:**

- Put the correction inside the primitive, so every caller of `RRPP` (the analysis, scripts, the case study) gets the same null from one code path.
- Ensure the analysis without the setting changes nothing, byte for byte. The committed fixtures stay as they are.
- Validate the shipped analysis path end to end, not a reimplementation of it.
- Predeclare the size rule, its multiplicity, and the tie-break between the two corrections, before any paper-grade unit runs.

**Non-Goals:**

- Correcting the observed statistic, or studentizing `delta`. Only the null moves.
- Exchangeability theory beyond the two corrections; for example, wild bootstrap or sign-flip nulls.
- Wiring the correction into the simulation harness, study configs or `motco de`. They have no covariates, and their committed signatures must not move.

## Decisions

### D1. The correction rescales the reduced residuals once, inside `RRPP`

`RRPP` gains the keyword `residual_correction: str | None = None`. After `y_res_np` is formed, and before either path, it computes `scale = residual_correction_scale(model_reduced, residual_correction)` and sets `y_res_np = y_res_np * scale[:, None]`. Nothing else changes:
- the serial loop and `_RRPPWorker` receive the rescaled array;
- the RNG calls, and therefore the permutation indices, are untouched;
- the full-model estimation in `estimate_difference` is untouched.

With `None` the code path is the current one, so byte identity needs no special casing.

*Alternative: callers transform `Y` before calling `RRPP`, as the diagnosis did.* Rejected. It only works for an orthogonal (uniform) rescaling, it duplicates the reduced fit, and every caller would have to repeat it.

*Alternative: a separate `RRPP_corrected` function.* Rejected. The parallel and serial machinery would be duplicated for a one-line difference.

### D2. Multipliers from the reduced model's thin SVD

`residual_correction_scale(model_reduced, correction)` is public and is used by both `RRPP` and the analysis record. It computes the thin SVD `X_r = U S Vᵀ` and r = #{s > s_max · max(n, p) · ε}, which is numpy's `matrix_rank` tolerance.

- **`df`:** κ = √(n / (n − r)) for every row. It is an error if n − r ≤ 0, and the message states n and r.
- **`leverage`:** h = row sums of `U[:, :r]²`, and the multiplier is 1 / √(1 − h_ii). It is an error if any 1 − h_ii ≤ 1e-10, and the message names the first such row index.
- Unknown values raise a `ValueError` naming `RESIDUAL_CORRECTIONS = ("df", "leverage")`. The check runs before any draw.

Why the rank and not the column count: `estimate_betas` already tolerates rank deficiency through `lstsq`, and the residual df is n − r. Why SVD: it gives the rank and an orthonormal basis of the column space in one step, and it costs nothing at n ≤ a few hundred.

The leverage-rescaled residuals are **not centred**. Adding a constant row vector to every row of `Y` shifts every LS-mean by the same vector, because each LS-mean row has intercept entry 1. So `delta`, `angle` and `shape` are invariant to it, and the permuted residuals' column mean is the same constant in every draw. A test pins this (spec scenario "Statistics do not depend on the residuals' mean").

*Alternative: the leverage-adjusted residuals `e_i / √(1 − h_ii)` centred and rescaled to the df variance (the bootstrap "modified residuals").* Rejected as a third variant. The extra rescaling is close to the `df` correction already offered, and the comparison the roadmap asks for is uniform against per-sample.

### D3. The analysis passes the setting through and records what it applied

- **Parameters.** `TrajectoryAnalysisParams.residual_correction: str | None = None` is validated in `__post_init__` against `RESIDUAL_CORRECTIONS`, so the error comes before any fit, and it is dropped from `to_dict()` while `None`.
- **RRPP call.** `analyze_trajectories` passes it to `RRPP`.
- **Record.** When it is set and `permutations > 0`, `measurement.residual_correction` is added. It is placed beside the other RRPP settings:
  - for `df`: `{"kind": "df", "n": 80, "reduced_rank": 10, "kappa": 1.069…}`;
  - for `leverage`: `{"kind": "leverage", "n": …, "reduced_rank": …, "leverage_min": …, "leverage_max": …, "scale_min": …, "scale_max": …}`.
- **Without permutations.** With `permutations = 0` the parameter is recorded, but no `measurement` entry is written: nothing was applied.
- **Covariates are not required.** The no-covariate reduced model also has r = 4, so κ = 1.026 at n = 80, and the option is meaningful there too.

The CLI adds `--residual-correction {df,leverage}` (default absent) and passes it to the parameters. When the correction is applied, it prints one extra line: `RRPP null: residual correction df (kappa = 1.069)`, or the leverage range for `leverage`.

*Alternative: a top-level `residual_correction` record block.* Rejected. The setting belongs to the permutation null, and `measurement` is where the null's settings live.

*Alternative: also correct the harness (`SimulationEvaluationParams`).* Deferred (Non-Goals). Adding a hashed field would need the hash-omission machinery and fixture pins, for no current consumer.

### D4. Naming

The parameter is `residual_correction`, with values `"df"` and `"leverage"`:
- `df` matches the diagnosis's verdict (`residual_df_shrinkage`) and roadmap item 8;
- `leverage` names the per-sample quantity.

*Alternative: `hc2` / `hc3`-style names.* Rejected. Those name heteroskedasticity-consistent covariance estimators, which is a different object, and the names would invite that confusion.

### D5. One check script with two profiles

`scripts/corrected_permutation_check.py` has these subcommands:
- `run --profile {type_i,power} [--table {astrocyte,microglia}] [--delta …] --shard i/N --out-dir <run dir>`;
- `merge --out-dir <run dir>`;
- `decide --runs <run dir>… --out-dir <dir>`.

It imports `case_study_type_i_check as cs` and `covariate_type_i_check as d10`. Settings come from `cs.build_settings` with these overrides:
- `base_seed = 2031`;
- the named cohort table in `generator.group_stage_sizes`;
- the profile's δ: 0.1 or 0.05 for `type_i`, 2.0 for `power`;
- `s_z = 10`;
- the profile's cells.

The settings also record `script_version`, `check = "corrected_permutation"`, the profile, the table name and its sha256 source, the arms, and the rule constants. A resume under other settings is refused with `cs.ensure_settings`.

One script, because the decision reads both profiles and their unit records are identical in shape. *Alternative: separate Type I and power scripts, as for latent adjustment.* Rejected. The unit code, the records and the merge would be duplicated, and `decide` would import across them anyway.

Cells:
- `type_i`: `anchor` (e = 0) and `translation` (e = 1.0);
- `power`: those two plus magnitude {0.02, 0.05, 0.1, 0.25} and orientation and shape {0.25, 0.5, 0.75, 1.0}, the latent-adjustment power grid.

`cs.check_headroom` runs on every cell before the first unit.

### D6. Three production analyses per unit, with a hard consistency check

Each unit calls `analyze_trajectories` three times. The calls use identical `TrajectoryAnalysisParams`: six covariates (`z`, `x1`–`x5`), `min_components = 2`, net displacement, the harness CV knobs, `permutations = 999`, the unit's RRPP seed and `n_jobs = 1`. They differ only in `residual_correction ∈ {None, "df", "leverage"}`.

The unit then asserts that the three `latent_scores` arrays are equal exactly, along with the selected rank and the three statistics. Any difference raises and fails the shard loudly; none is expected, because the CV fit is deterministic in process.

Each arm records:
- the statistics, p-values and null summary;
- the selected rank, the component selection and `cv_mean_auroc`;
- the Spearman correlation of the first latent component with `z`;
- the record's `measurement.residual_correction`;
- a sha256 digest of the latent scores.

*Alternative: one analysis plus RRPP-only corrected arms, rebuilding the design through `_covariate_adjusted_design` (the diagnosis's inference path).* Rejected. It is about 3× cheaper (one CV fit per unit instead of three), but it validates a reimplementation rather than the shipped path. At about 140 core-hours the production path is affordable on `ing`. Tests separately pin the equivalence: the analysis's corrected p-values equal `RRPP` with the correction on the analysis's own inputs.

### D7. Fresh seed family, construction kept, calibration descriptive

Base seed 2031 is not used by any committed check. Every seed, including the injection directions, the covariates, the noise covariates, the generator and RRPP, derives from it and the replicate index. The two controls of a replicate (and, in `power`, all 14 cells) share them.

The diagnosis chose the remedy on the parent's replicates (base seed 2026), so validating on them would be circular. A fresh family is a confirmatory replication of both the excess and the fix.

The construction's δ and `s_z` stay those calibrated on the astrocyte table (the parent's `calibration.json`); they are not re-tuned. Per run, the decision reports these criteria descriptively from the uncorrected arm:
- the mean held-out CV AUROC against the parent's 0.60–0.75 band;
- the median |ρ| of the first component with `z` against 0.8.

The parent's third criterion needs an unadjusted arm, which this check does not run.

### D8. Size rule: a pooled tier and a per-reading tier

A reading is a (run, control, statistic) rate over R available replicates, with α = 0.05. A correction holds size when both tiers hold:
1. **Pooled, per cohort table × statistic.** The rate over the table's four readings (2 strengths × 2 controls, R_pooled = 2,000) is ≤ α + 2·√(α(1 − α)/R_pooled), i.e. 0.0597 or at most 119 / 2,000 rejections.
2. **Per reading.** Every reading of the four Type I runs, and the power profile's two controls (30 readings), is ≤ α + 3·√(α(1 − α)/R), i.e. 0.0792 or at most 39 / 500.

Why two tiers and not the precedent's per-reading k = 2 conjunction:
- At exact size α, one reading exceeds 0.0695 with probability 0.030 (exact binomial). With 24 readings, a per-reading rule at k = 2 fails a valid correction up to 52 % of the time (the bound if the readings were independent). That is too often.
- The pooled tier keeps the sensitivity that matters. A true rate of 0.066, the uncorrected `delta` the diagnosis measured, exceeds 0.0597 at R = 2,000 with probability 0.87.
- The per-reading tier still catches a localized excess.
- At exact α the two tiers together fail a valid test at most 21 % of the time (6 pooled tests at 0.025 each, plus 30 readings at 0.0027 each, under independence). The diagnosis measured the `df` correction as conservative (`delta` 0.048–0.0545, `angle` 0.042, `shape` 0.0465). At a true size of 0.048 the risk is 8.7 %, and at 0.045 it is 1.8 %.

Pooling per table, not over both tables, keeps a microglia-specific excess from being diluted by the astrocyte readings.

The uncorrected arm is run through the same rule descriptively. `uncorrected_excess_reproduced` is true when its pooled `delta` exceeds the pooled bound in either table.

### D9. Tie-break: `df` unless `leverage` earns it on paired power

When both corrections hold size, `leverage` is compared with `df` on the power profile's units. For each target pair (magnitude/`delta`, orientation/`angle`, shape/`shape`):
- each replicate's d_r is the number of the mode's four effect cells that `leverage` rejects and `df` does not, minus the reverse;
- the one-sided p-value comes from the replicate-clustered exact sign-flip test, in each direction, using `case_study_diagnosis.sign_flip_p_value`.

`adopt_leverage` requires a gain at p ≤ 0.05 on at least one target pair and no loss at p ≤ 0.05 on any. Otherwise the verdict is `adopt_df`.

Why `df` is the default:
- it is a single scalar per analysis, easy to state and audit;
- it is the correction the diagnosis already measured;
- C4's residual net of 11 (p = 0.156) is the only hint that leverage structure matters, and that hint is weak.

Clustering by replicate is needed because a mode's four effect cells share the replicate's seeds.

*Alternative: choose by the larger total rejection rate.* Rejected. It would reward any difference, however small, and it ignores the pairing.

### D10. Verdicts, consequences and the release

The verdict is computed by a pure function `corrected_permutation_decision(size_rows, power_rows, settings)`.

| verdict | condition | consequence |
|---|---|---|
| `adopt_df` | `df` holds, and `leverage` does not or loses the tie-break | `sead_sexdiff` re-issues both cell types' primaries with `--residual-correction df`, pinned to 0.13.0, in its own change |
| `adopt_leverage` | both hold and `leverage` wins the tie-break, or only `leverage` holds | as above, with `leverage` |
| `no_correction_holds` | neither holds | the correction ships but is documented as not validated at this configuration; the case study keeps reading `delta` against its realized size; the roadmap records the failing readings and names a diagnosis as the next step |

The decision file also carries:
- every reading and pooled rate with its bound;
- the paired contrasts and the uncorrected flag;
- the descriptive calibration criteria;
- every run's settings hash.

A run at fewer permutations than 999 marks the verdict indicative.

**Version.** `src/motco/__init__.py` goes to 0.13.0 with the package change, before the `ing` run. That way the records' `software.motco` names the release that ships the code. The verdict does not gate the release of an opt-in setting, which is the 0.9.0 precedent. It gates what the docs claim and whether the case study re-issues.

### D11. Records, layout and outputs

`results/corrected-permutation-<date>/` has one directory per run:
- `type-i-astrocyte-delta0.1/`, `type-i-astrocyte-delta0.05/`, `type-i-microglia-delta0.1/` and `type-i-microglia-delta0.05/`;
- `power/`.

Each holds `settings.json`, `records/shard_<i>_of_<N>.jsonl` (one line per unit, all three arms) and, after `merge`:
- `corrected_permutation_pvalues.csv`, with columns `run, mode, effect, replicate, arm, statistic, p_value`;
- `corrected_permutation_rates.csv`, per arm × cell × statistic: R, rejections, rate and the rank distribution;
- `corrected_permutation_contrasts.csv`, power only: paired cost and tie-break contrasts per cell and statistic.

`decide` writes these at the top level:
- `corrected_permutation_decision.json`;
- `size_readings.csv`: every judged reading and pooled rate with its bound and verdict.

`records/` and any CSV over about 1 MB stay on `ing`, named with their sha256 in each run's `PROVENANCE.txt`.

### D12. Run plan

1. **Smoke locally.** Every run at R = 20 × 99 permutations, then `merge` and `decide` (indicative). In addition, R = 4 × 999 on one Type I run, to exercise the full permutation count.
2. **Time** units at 999 permutations per profile, and append a dated "Smoke and shard plan" section to this design.
3. **Approval gate.** The user approves the arrays.
4. **On `ing`,** from a clean worktree of the pushed branch, with matching tree hashes and tests passing there:
   - four 50-shard Type I arrays (1,000 units each);
   - one power array of about 140 shards (7,000 units);
   - all on partition `512x1024`, n11 excluded;
   - `scripts/corrected_permutation_array.sbatch`, modelled on `scripts/case_study_diagnosis_array.sbatch`.
5. **Finish:** `merge` per run, `decide`, `PROVENANCE.txt` per run, the report, and the docs.

## Risks / Trade-offs

- [The fresh family does not reproduce the uncorrected `delta` excess (`uncorrected_excess_reproduced: false`)] → The correction's size verdict still stands. The report states that the excess was not replicated and reads the diagnosis's C2 as the evidence for the mechanism. No rule is re-tuned.
- [A valid correction fails on a marginal pooled reading (at most 21 % at exact α, 8.7 % at 0.048, D8)] → The verdict is issued as predeclared. The report reads the failing table and statistic against the diagnosis's rates, and any follow-up is its own change.
- [`leverage` is undefined for some unit (leverage 1)] → It cannot happen with numeric covariates and every cell occupied, but if it does the analysis raises. The unit fails loudly and is counted, and the arm is unavailable for that unit. That would itself be a finding against `leverage` for real data.
- [The three CV fits of a unit diverge (nondeterminism)] → `n_jobs = 1` inside units. The exact-equality assertion turns divergence into a shard failure, not a silent bias.
- [The correction is conservative for `angle` and `shape`] → This is expected (the diagnosis measured `angle` 0.042 and `shape` 0.0465). The power profile quantifies the cost, and the report states it beside the verdict.
- [The microglia table at the astrocyte calibration reads outside the AUROC band] → The criteria are descriptive (D7). The size claim is about the cohort table, not the band.
- [The `cohort_50` tables (n = 63 microglia, κ = 1.090) are not run] → This is stated as a limitation. The mechanism scales with p_red / n, so the correction matters more there, not less.
- [Compute: 3 CV fits per unit, about 140 core-hours] → Re-timed at smoke before the approval gate. Shards are resumable.

## Migration Plan

The API change is additive and the setting is opt-in, so nothing migrates. Callers that never pass the setting get 0.12.0's outputs, checked against the committed identity fixtures. Rollback is to not pass the setting, or to pin 0.12.0.

The case study adopts the correction in its own change, after the release and under the verdict's consequence.

## Open Questions

- The number of power shards is set by the smoke timing (D12). It does not affect any rule.
- Whether `corrected_permutation_pvalues.csv` is committed depends on its size at merge (9,000 rows per Type I run and 63,000 for the power run): it is committed if it is under about 1 MB. Otherwise it is named in `PROVENANCE.txt`. The decision is unaffected.

## Smoke and shard plan (2026-10-30)

Run locally on the workstation (24 cores, single-threaded BLAS) from branch `feat/residual-df-corrected-permutation` at `a2f1a1f` (tree `4d8a22fb…`), MOTCO 0.13.0.

**Smoke at R = 20 × 99 permutations** (`results/corrected-permutation-2026-10-30/smoke-type-i-{astrocyte,microglia}-delta{0.1,0.05}/`, `smoke-power/`, decision in `smoke-decision/`; 2 shards per Type I run, 12 for power):
- 40 / 40 units per Type I run and 280 / 280 power units, 0 arm failures and 0 consistency failures (every unit's three arms shared their latent scores, rank and statistics).
- The decision is marked indicative (99 permutations in every run) and reads `no_correction_holds`. The only failing entries are the pooled `shape` readings, 9–10 / 80 against the R = 80 bound of 0.099, in every arm including the uncorrected one. At R = 20 per reading this is Monte Carlo noise, not a size reading; the replicates are shared with the paper-grade run, so replicates 0–19 will be read again there.
- Pooled `delta` (astrocyte / microglia): none 2 / 7, `df` 0 / 6, `leverage` 0 / 6 of 80, so both corrections move `delta` down, as the diagnosis measured. `uncorrected_excess_reproduced` is false at R = 20.
- Calibration (uncorrected arm, anchor): mean held-out CV AUROC 0.695–0.716 in every Type I run (band 0.60–0.75) and median first-component |ρ| with `z` 0.93–0.97 (threshold 0.8); CV rank 2 or 3. The power run reads AUROC 1.000 and |ρ| 0.81.
- Power at R = 20 (none / `df` / `leverage`, top effect): magnitude `delta` 19 / 19 / 19, orientation `angle` 16 / 16 / 16, shape `shape` 7 / 6 / 6. The tie-break contrasts are empty (one discordant unit on magnitude `delta`, p_gain 0.5).

**999 permutations** (scratch, `build/corrected-permutation-2026-10-30/`): `type-i-astrocyte-delta0.1` at R = 4 (8 / 8 units), `type-i-microglia-delta0.1` at R = 2 (4 / 4) and `power` at R = 1 (14 / 14), each end to end with merge and 0 failures.

**Timing.** At 999 permutations a unit (three analyses) takes 19.0 s (astrocyte), 19.5 s (microglia) and 18.2 s (power) on one uncontended workstation core, against 19.3 s for the parent check's unit (`case_study_type_i_check.py`, δ = 0.1, `s_z = 10`). The parent measured 54.7 s per unit on `ing`, so the plan is about 54 s per Type I unit and 52 s per power unit on an `ing` EPYC core (the diagnosis ran at 31 s on less loaded nodes).

**Shard plan.**
- Four 50-shard Type I arrays (20 units per shard, about 18 min) and one 140-shard power array (50 units per shard, about 45 min), 500 × 999, on partition `512x1024` with `--exclude=n11` (in the template), one CPU and BLAS threads at 1. The 8 h task limit covers a node running 10× slow.
- Code: the clean worktree `/home1/tgonza/MOTCO-corrected-permutation` of the pushed branch at `a2f1a1f`, whose tree hash equals the workstation's; `PYTHONPATH` on its `src`, `.venv` symlinked from `/home1/tgonza/MOTCO`. `tests/test_residual_correction.py` and `tests/test_corrected_permutation_check.py` passed there (60 tests; MOTCO 0.13.0 imported from the worktree). Each `settings.json` is pre-created with `run … --shard 0/N --limit 0` before its array is submitted (`/home1/tgonza/corrected_permutation_submit.sh`).
- Outputs: `results/corrected-permutation-2026-10-30/type-i-{astrocyte,microglia}-delta{0.1,0.05}/` and `power/` in that worktree.
- About 60 core-hours for the Type I runs and 100 for power, 160 in all (design estimate 140).

**Approval:** the user approved committing and pushing the branch and submitting the five arrays as planned on 2026-10-30 (10:52 local), before the timing above was complete; submission waited for the smoke and the 999-permutation units.

**Open item found at smoke (for the user, before `decide`).** A table's four Type I readings share replicate seeds, and the two strengths of one replicate are strongly dependent: the smoke's per-replicate p-values correlate at Spearman 0.77–0.96 between δ = 0.1 and 0.05 (anchor), while anchor and translation are nearly independent (−0.05 to 0.27). The pooled tier (D8) treats the four readings as R_pooled = 2,000 independent replicates, so its 0.0597 bound is tighter than its stated 2.5 % one-sided risk per test: at an effective R of about 1,100–1,200 a correction at exact size fails a pooled test with probability about 6 %, not 2.5 %. D8's 21 % overall figure is therefore optimistic. The data the arrays collect do not depend on the rule; whether to amend the pooled tier (for example, a replicate-clustered standard error) before the results are read is the user's decision, recorded here before any paper-grade result exists.
