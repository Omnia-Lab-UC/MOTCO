# Design

## Context

See proposal.md (Why) for the motivation. The pieces this change composes, all in 0.12.0:

- `scripts/covariate_type_i_check.py` (D10 of `covariate-adjusted-trajectories`) builds a null Phase 6 replicate from `settings["generator"]` (so per-omic `delta_*` can be added to that mapping), draws a numeric covariate `z = a_g·g + a_s·s + a_gs·g·s + ε` and a three-level categorical on the same linear predictor, converts methylation to M-values, injects `s_z·sd_f·w_f·z + s_c·sd_f·Σ u_{k,f}·1[cat = k]` along directions drawn once from the base seed, and returns metadata with `z` and `category`. Its seeds depend on the base seed and replicate index only. It exports the settings hash, `ensure_settings`/`load_settings`, `parse_shard`, `shard_path`, `anchor_bound` and the CSV writer; `scripts/latent_adjustment_type_i_check.py` already imports it (`import covariate_type_i_check as d10`) and adds `extra_covariates` pure-noise numerics from their own seed.
- `scripts/rank_floor_check.py` has the unit grid (cells × replicates, replicate-major, sharded by unit), the weakened stage signal (`PHASE6_STRENGTHS` with `delta_*` and `p_dmp` as settings), the `calibrate` subcommand appending to `calibration.json`, per-arm records with `cv_mean_auroc`, `null_summary` and `pooled_spectrum`, the `_failed_arm` convention, and a decision over the anchor and translation controls. Its calibration found, at the Phase 6 design without covariates, mean held-out CV AUROC 0.58 at δ = 0.1 and 0.70 at δ = 0.15 (rank-1 fractions 0.44–0.55 and 0.08).
- `motco.analysis.analyze_trajectories` records `latent.selected_lv`, `latent.cv_mean_auroc`, `null_summary`, `config_spectrum`, a `covariates` block (`design_columns`) when covariates are declared and a `latent_adjustment` block when `covariate_latent_adjustment` is set; failures raise `TrajectoryAnalysisError`. `TrajectoryAnalysisParams.min_components` defaults to 2.
- The case study (`sead_sexdiff`, 2026-10-26): `ANALYSIS_COVARIATES` = six numerics, `missing_covariates="error"`, floor 2, net displacement, 999 permutations, held-out CV AUROC 0.6356–0.6965 across the four cohorts, first latent component |Spearman| 0.84 with a library-quality covariate in the unadjusted basis, latent-adjusted `stage` arms as sensitivity analyses.
- Evidence on the arms: the 0.9.0 Type I check's descriptive baseline arm (model adjustment only, nuisance taking the first component, eight covariates, floor 1, δ = 2.0) rejected `shape` at 0.078; its `stage` arm held (0.031 / 0.067 / 0.058) and `group_stage` exceeded on `angle` (0.082). Its `stage` arm lost 52 replicates to a categorical level missing from a CV training set.
- Cluster: `ing`, partition `512x1024`, 18–30 s per `analyze_trajectories` call at 999 permutations on an EPYC core; runs go from a clean worktree of the pushed branch with `PROVENANCE.txt` carrying the tree hash and file hashes; node n11 has stalled imports under NFS load.

## Goals / Non-Goals

**Goals:**

- Measure size where the case study is actually run — weak stage signal, six design columns, floor 2, a nuisance the basis follows — with a rule and a fallback written before the run, so the report can say which of the case study's two analysis types is admissible.
- Reuse the two validated constructions by import so the new check inherits their reproducibility, sharding and settings-hash discipline, and so its calibration can be read against theirs.
- Place the construction on the axes a real analysis reports (held-out CV AUROC, first-component correlation), not on generator strengths alone.

**Non-Goals:**

- Power at this configuration, the microglia cell sizes, the `group_stage` conditioning (inadmissible at 0.9.0), a floor-1 arm (the rank-floor check covers the floor's effect), or any change to `src/motco/`.
- Deciding anything in `sead_sexdiff`: the report states the consequence; acting on it is that repository's change.

## Decisions

### D1. A new script composing D10 and the rank-floor pattern, not flags on either

`scripts/case_study_type_i_check.py` imports `covariate_type_i_check as d10` for the construction and helpers (as the latent-adjustment check does) and follows `rank_floor_check.py` for the unit grid, calibration record, per-arm record shape and control-based decision.

*Alternative: `--min-components`, `--noise-covariates` and `--delta` flags on `covariate_type_i_check.py`.* Rejected: its settings hash identifies the committed 2026-10-16 run, its `read_records` and tests are keyed on its two arms, and its decision names `covariates_validated`; a third arm, a second cell and a three-way decision would make it a different check wearing the old name.

*Alternative: covariates added to `rank_floor_check.py`.* Rejected: its 14-cell grid and two floor arms are the floor's validation, not this one's, and it has no covariate machinery.

### D2. Six numeric design columns: `z` plus five independent noise covariates, no categorical

The check declares `covariate_cols = ["z", "x1", …, "x5"]` with no `categorical_covariates`, draws the noise covariates i.i.d. N(0, 1) from a dedicated seed (the latent-adjustment check's `extra_covariates` pattern) and sets D10's categorical injection strength `s_c` to 0; the `category` column is still drawn (so D10's RNG stream and `z` are unchanged at equal strengths) but not declared or injected.

Why: the case study enters six numerics and no factor, so six numeric columns is its configuration; without a factor no CV training set can lose a level, so the `latent_stage` arm keeps R = 500 (the 0.9.0 check lost about 10 % of its `stage` replicates to that); and the confounding rests on `z` alone, which the calibration's third criterion verifies rather than assumes.

*Alternative: D10's pair (`z` and the categorical, three columns) plus three noise covariates.* Rejected: a factor the case study does not have, with the dropout failure mode, in exchange for confounding already established at δ = 2.0 that must be re-verified at the weak signal anyway.

*Alternative: six confounded covariates along six directions.* Rejected: it would test a configuration harsher than the case study's by construction; the five noise covariates measure the degrees-of-freedom cost the D10 report declined to bound, which is the open question.

### D3. Calibration on the regime, not on the generator strengths

Three criteria on a smoke run (R = 20–50, 99 permutations), each recorded per iteration in `calibration.json`: the `adjusted` arm's mean held-out CV AUROC inside the declared band `[0.60, 0.75]` (the cohorts read 0.636–0.697; the band is fixed here, before any run); the `unadjusted` arm's median first-component |ρ| with `z` at least 0.8 (the real basis reads 0.84; the 0.9.0 criterion); and the `unadjusted` arm rejecting above α on at least one statistic at the anchor (the D10 non-vacuity condition). The knobs are the common per-omic δ (one `--delta` applied to `delta_methyl`, `delta_expr` and `delta_protein`; `p_dmp` stays 0.1 so the translation surgery keeps its Phase 6 headroom) and D10's `s_z`; `a_g`, `a_s`, `a_gs` stay at D10's values.

Starting point: δ = 0.1 and 0.15 (0.58 and 0.70 without injection; `z` is confounded with stage, so the injection adds stage-predictive signal and the band is read on the adjusted arm's realized AUROC), `s_z = 3` (D10; the 0.9.0 check reached |ρ| 0.64 at 3 and 0.83 at 40 against a saturated stage signal, so at a weak signal the threshold should be met well below 40). One or two strengths run at paper grade is settled after calibration with the user; with two, the decision is the conjunction, as the rank-floor check's was.

*Alternative: reuse the rank-floor check's rank-1-fraction criterion.* Rejected: the question here is not whether the floor binds but whether the case study's analysis holds size at its stage signal; the realized rank distribution is recorded, not required.

### D4. Cells: the anchor and the translation control at e = 1.0

Both controls, as in the rank-floor check and the design-guidance slice's column-scoped rule; the net-displacement Type I study found translation's `angle` near the bound (0.068 at Phase 6), so it is the control most likely to move under a covariate-adjusted model. Seeds depend on the replicate index only, so both cells and all three arms of replicate r share the baseline draw, the covariates and the permutation draws.

*Alternative: anchor only (D10).* Rejected: cheap to add, and a report that judged only the anchor would be weaker than the two checks it sits between.

### D5. Arms and the decision order

`adjusted` (judged first; the case study's primary), `latent_stage` (judged second; the case study's sensitivity arm and the only conditioning admissible at 0.9.0), `unadjusted` (descriptive). All at `min_components=2`, net displacement, the harness CV knobs and D10's `random_state`. The `adjusted` and `unadjusted` latent spaces are asserted identical per unit (D10's check), so an `adjusted` excess is attributable to the trajectory model and an `unadjusted` excess to confounding. The decision is `primary_validated` → `adopt_latent_stage` → `diagnose_configuration_excess`, over both controls and all three statistics, with `holds` | `exceeds` | `unavailable` per (arm, control, statistic).

Why the fallback is predeclared: the 0.9.0 baseline arm's `shape` 0.078 is the warning that the primary may fail exactly on `shape` while `stage` held; writing the fallback into the rule now means the outcome, either way, is a decision rather than a diagnosis.

*Alternative: a `group_stage` arm.* Rejected: inadmissible at 0.9.0; not a case-study analysis.

*Alternative: a floor-1 `adjusted` arm.* Rejected: the rank-floor check already measures the floor's effect on these controls; the cost is a fourth CV fit per unit; the basis-identity assertion gives the attribution a floor-1 arm would.

### D6. Records and tables

One JSONL line per unit with `key` = `[mode, effect, replicate]`, `seeds`, the D10 covariate-cell summary, and `arms`: per arm `statistics`, `p_values`, `selected_lv`, `component_selection`, `cv_mean_auroc`, `first_component_spearman_z`, `null_summary`, `pooled_spectrum`, `covariate_design_columns` (adjusted, latent_stage), `latent_adjustment` (latent_stage: `residual_df`, `cv_min_training_df`, `feature_r2` summary), `error`. Merge writes `merged.jsonl`, `case_study_type_i_pvalues.csv` (unit × arm × statistic), `case_study_type_i_replicates.csv` (unit × arm: rank, AUROC, |ρ|, `angle` null width, residual df, error), `case_study_type_i_rates.csv` (control × arm × statistic: R, rejections, rate, bound, verdict) and `case_study_type_i_decision.json` (the rates, the realized calibration values — adjusted AUROC mean and quantiles, unadjusted median |ρ|, unadjusted anchor rates — the construction, the settings hash and the decision).

### D7. Sharding, resume and the pure decision

Units are enumerated replicate-major and shard i of N takes every N-th; a unit recorded in any shard file counts as done (so a re-run under a different shard count skips recorded work); torn lines are dropped and re-run; `settings.json` with its hash refuses a resume under different settings. `case_study_type_i_decision(rows, alpha, k)` is a pure function over the p-value rows, tested on synthetic rows for the three outcomes and the conjunction (`combine_decisions`), as the earlier checks' are.

### D8. Run and provenance

Smoke locally (R = 20, 99 permutations) until D3's criteria hold; the user approves the strengths and the number of strengths; a 50-shard array per strength on `ing` from a clean worktree of the pushed branch (`scripts/case_study_type_i_array.sbatch`, modelled on the rank-floor template, excluding n11 if it stalls), `merge` on `ing`, then `PROVENANCE.txt` per run (command, job, tree hash, settings hash, versions, sha256 of the files that define the measurement, core-hours). `records/` and `merged.jsonl` stay on `ing` under the gitignore; the CSVs are committed when each is under about 1 MB (about 9,000 p-value rows and 3,000 replicate rows per strength, so they should be), otherwise kept on `ing` and named in `PROVENANCE.txt`.

### D9. No release

Every record field the check reads exists in 0.12.0, so `src/motco/__init__.py` stays at 0.12.0 and the case study's pin is untouched. If implementation finds a field missing (none is expected), it is added as a patch, 0.12.1, with the record change named in the report.

## Risks / Trade-offs

- [The three calibration criteria may not be jointly satisfiable: raising `s_z` to reach |ρ| ≥ 0.8 adds stage-predictive signal because `z` is confounded with stage, pushing the adjusted arm's AUROC above the band] → Lower δ compensates; every iteration is recorded. If no (δ, `s_z`) meets all three, the conflict is reported with the iterations before any paper-grade run and the user chooses which criterion to relax; nothing runs at paper grade on an unmet criterion without that decision.
- [The primary fails on `shape`, as the 0.9.0 baseline arm suggests it might] → That is the outcome the fallback exists for; `adopt_latent_stage` is a decision, and its consequence (the case study's primary moves to `stage`) is a `sead_sexdiff` change the report names.
- [The `latent_stage` arm fails on a unit for a reason other than factor dropout (a rank-deficient nuisance design at a CV training set)] → Recorded as unavailable with the error; R is reduced for that arm and the count of unavailable units is reported. With six numerics at n = 80 the smallest training set has about 50 rows, so this is not expected.
- [Translation's `angle` exceeds the bound under covariates for a reason unrelated to the configuration] → The control is mandatory by convention; the `unadjusted` arm's translation rates and the basis-identity assertion let the report attribute the excess, and `diagnose_configuration_excess` names it.
- [The real covariates are correlated with one another and bounded (age, education), unlike five independent normals] → A stated limitation; the configuration matched is the column count, the floor and the nuisance regime, not the covariates' joint distribution.
- [A single confounded covariate understates a cohort where several covariates carry nuisance] → Stated limitation; D2 records the alternative considered.
- [n11 NFS stalls hold shards for 30+ minutes] → Exclude n11 or `scancel` and resubmit the stuck array ids; shards are resumable.

## Migration Plan

No package or API change; nothing to migrate or roll back. The run sequence is: smoke and calibration locally → user approval of strengths → array on `ing` → merge → `PROVENANCE.txt` → report → roadmap and docs. The results directory is `results/case-study-type-i-<date>/` with one subdirectory per strength when more than one is run.

## Open Questions

- One or two stage-signal strengths at paper grade: settled after calibration with the user. It changes neither the specs (the conjunction is specified) nor the task list.
- Whether the per-unit CSVs are committed: decided by their size at merge (D8).

## D3 calibration record (2026-10-29)

Local workstation (16 cores), MOTCO 0.12.0 at the tree of this change, `calibrate` at R = 20, 99 permutations, base seed 2026, D10's `a_g = 0.5`, `a_s = 1`, `a_gs = 2`, five noise covariates, band `[0.60, 0.75]`, |ρ| threshold 0.8. Every iteration is in `results/case-study-type-i-2026-10-29/calibration/calibration.json`.

| # | δ (per-omic) | `s_z` | adjusted mean AUROC [min, max] | unadjusted median \|ρ\| [min, max] | unadjusted anchor rates δ / angle / shape | met |
|---|---|---|---|---|---|---|
| 1 | 0.10 | 3 | 0.691 [0.611, 0.740] | 0.565 [0.329, 0.957] | 0.55 / 0.05 / 0.30 | no (\|ρ\|) |
| 2 | 0.15 | 3 | 0.744 [0.647, 0.805] | 0.571 [0.338, 0.954] | 0.55 / 0.05 / 0.35 | no (\|ρ\|) |
| 3 | 0.10 | 10 | 0.699 [0.646, 0.738] | 0.876 [0.688, 0.991] | 0.55 / 0.10 / 0.25 | **yes** |
| 4 | 0.10 | 20 | 0.699 [0.648, 0.752] | 0.918 [0.753, 0.993] | 0.65 / 0.10 / 0.30 | yes |
| 5 | 0.05 | 10 | 0.692 [0.640, 0.734] | 0.873 [0.686, 0.991] | 0.60 / 0.05 / 0.25 | **yes** |
| 6 | 0.05 | 20 | 0.693 [0.640, 0.741] | 0.917 [0.747, 0.994] | 0.65 / 0.10 / 0.25 | yes |

What the iterations show: at `s_z = 3` the basis does not follow `z` (median |ρ| 0.57, the D10 regime); at `s_z ≥ 10` it does (0.87–0.92, bracketing the real 0.84), and the realized stage AUROC is then set by the injected covariate rather than by δ — 0.69–0.70 at both δ = 0.05 and 0.10, inside the band and on top of the cohorts' 0.636–0.697 — because `z` is confounded with stage. The three criteria did not conflict (design risk 1 did not materialize). The adjusted arm selects rank 2 in 18–19 of 20 anchor replicates at every met iteration; the `latent_stage` arm in 16–18.

**Settled strengths.** `s_z = 10` (the smallest injection meeting the threshold, so the nuisance is the case study's, not harsher) at δ = 0.10 (iteration 3; the primary candidate) and δ = 0.05 (iteration 5; the second stage-signal strength, the rank-floor check's pair). `p_dmp` stays 0.1.

**Smoke profiles** (`run`, R = 20, 99 permutations, both cells, three arms; `results/case-study-type-i-2026-10-29/smoke-delta0.1-sz10/` and `smoke-delta0.05-sz10/`): 40/40 units at each strength, 0 arm failures, the `latent_stage` arm available in every unit (no factor, so no training-set dropout; smallest `cv_min_training_df` recorded). Rates at R = 20 (bound 0.1475; indicative only):

| strength | cell | adjusted δ / angle / shape | unadjusted | latent_stage |
|---|---|---|---|---|
| δ 0.10, `s_z` 10 | anchor | 0.10 / 0.05 / 0.00 | 0.55 / 0.10 / 0.25 | 0.10 / 0.05 / **0.35** |
| δ 0.10, `s_z` 10 | translation e = 1 | 0.10 / 0.00 / 0.10 | 0.75 / 0.05 / 0.45 | 0.00 / 0.10 / 0.10 |
| δ 0.05, `s_z` 10 | anchor | 0.10 / 0.05 / 0.05 | 0.60 / 0.05 / 0.25 | 0.10 / 0.05 / **0.35** |
| δ 0.05, `s_z` 10 | translation e = 1 | 0.05 / 0.00 / 0.05 | 0.90 / 0.05 / 0.40 | 0.00 / 0.05 / 0.10 |

The smoke decision is `primary_validated` at both strengths. The `latent_stage` arm's anchor `shape` rate (7 of 20 at both strengths) is the one warning sign: it exceeds even the R = 20 bound, so the fallback arm may fail on `shape` at paper grade — which does not touch the primary's verdict but is a finding the report must carry (the 0.9.0 check's `stage` arm held `shape` at 0.058 at floor 1 and eight covariates with a saturated stage signal; here the floor is 2, the signal weak and the basis residualized on six numerics). The unadjusted arm confounds at both cells (δ 0.55–0.90, `shape` 0.25–0.45).

**Timing and shard plan.** One unit (three `analyze_trajectories` calls at 999 permutations) took 57 s locally with 12 other shards running; the rank-floor run on `ing` recorded 11.8 s per 999-permutation analysis at saturated-node cost (45.8 core-hours for 14,000 analyses), so a unit here is planned at 40–60 s. Per strength: 1,000 units (2 cells × 500 replicates), 3,000 analyses, 50 shards × 20 units, 15–25 min per shard, 11–17 core-hours (planned at 60 s: 16.7). Two strengths: 2,000 units, 6,000 analyses, about 35 core-hours, well inside the 6 h task limit of `scripts/case_study_type_i_array.sbatch`. The Phase 5 unit cost guidance in the proposal (18–30 s per analysis) gives 15–25 core-hours per strength as the upper plan.

**Approval.** Given 2026-10-29: two paper-grade runs, δ = 0.10 and δ = 0.05 at `s_z = 10` (iterations 3 and 5), 500 replicates × 999 permutations each, base seed 2026, 50 shards per strength on `ing`; the overall decision is the conjunction of the two runs' decisions (`combine`).
