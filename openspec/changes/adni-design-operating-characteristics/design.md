# Design

## Context

See proposal.md (Why) for the motivation and `specs/adni-design-validation/spec.md` for the requirements. The design builds on these pieces of 0.13.0:

- **The case-study construction** (`scripts/case_study_type_i_check.py`, imported as `cs`; it imports `scripts/covariate_type_i_check.py` as `d10`).
  - `cs.build_settings` takes the D10 settings (`a_g`, `a_s`, `a_gs`, `s_z`, `feature_fraction`, the categorical logit scale 0.5 and `s_c = 0`), sets the per-omic δ from `--delta` with `p_dmp` 0.1, the floor 2, the noise covariates and the calibration criteria. `cs.build_unit` calls `d10.build_replicate` (the generator at `generator.group_stage_sizes`, the draw of `z` and the three-level `category` from one linear predictor, the injection of `z` along a fixed direction over 20 % of the joint features, methylation converted to M-values first) and appends the noise covariates. `cs.declared_covariates` returns `["z", "x1", …]`.
  - `scripts/corrected_permutation_check.py` shows the composition pattern this change follows: parent args through `cs.build_parser().parse_args([...])`, a `TABLES` registry with the cell table and its source, overrides on the settings, `cs.ensure_settings`/`settings_hash`/`parse_shard`/`shard_path`/`unit_key`/`anchor_bound` reused, its own `evaluate_unit`, `merge` and pure decision.
  - Seeds: `d10.replicate_seeds(base_seed, replicate)` gives the generator, covariate and RRPP seeds; `d10.feature_directions` draws the injection direction from the base seed; `cs._seed(base_seed, "noise_covariates", replicate)` the noise. A new base seed gives a new family throughout. Committed base seeds: 100, 400, 500, 600, 800, 900, 1028, 2026, 2031.
- **The real-data analysis.** `TrajectoryDataset(covariate_cols=, categorical_covariates=)` codes a factor drop-first as `name[level]` and records `covariates.design_columns`; the record carries `latent.selected_lv`, `latent.cv_mean_auroc`, `null_summary`, `config_spectrum.pooled.relative_eigengap`, `measurement.residual_correction` (absent here). `analyze_trajectories` without `covariate_latent_adjustment` and `residual_correction` is byte-identical from 0.10.0 (ADNI's pin) to 0.13.0 (`tests/test_analysis_identity.py`).
- **The generator at a large table.** `SemiSyntheticTrajectoryParams(group_stage_sizes=)` fixes every cell exactly; a cell of 448 samples is an MVN draw of 448 × 367 and 448 × 131. The surgery headroom (`expected_surgery_headroom`) depends on `p_dmp`, `n_stages` and ρ, not on n, so the translation, orientation and shape surgeries fit at e = 1.0 as they do in every Phase 6 check at `p_dmp` 0.1.
- **The reference's blocks.** Methylation 367 CpGs, expression 131 genes (protein 160, unused): 498 features in two blocks. ADNI: p180 133, NMR 165 (298).
- **The ADNI design** (its primary `result.json`, MOTCO 0.10.0, 2026-10-27 to 31): cells F 237 / 310 / 125, M 223 / 448 / 158; design columns `age`, `education`, `apoe4[1]`, `apoe4[2]`, `phase[ADNI2GO]`, `statin[1]`; reference row at the pooled means (phase 0.544); CV rank 8, held-out AUROC 0.600 (0.580–0.616 across its CV-selected sensitivity analyses, ranks 4–8); pooled relative eigengap 0.744; `angle` null q95 35.4°; first latent component ρ 0.542 with `phase` (`pls_3` −0.535); the p180 block's PC1 |corr| 0.73 with phase and a median |d| 1.12 over its 133 metabolites (80 above 1), the NMR block's 0.21 and 0.13. ADNI2GO share per cell: F CN 143/237, MCI 193/310, AD 53/125; M CN 125/223, MCI 227/448, AD 76/158. Runtime 26 s at 9,999 permutations on the workstation.
- **Timing precedents.** A covariate-adjusted analysis at n = 80 takes about 6.5 s on the workstation and 17 s on an `ing` EPYC core. Design-guidance cells at n = 1,000 with three blocks and 199 permutations took up to 107 s per unit on `ing`.

## Goals / Non-Goals

**Goals:**

- Measure, at ADNI's n, cells, column kinds, nuisance regime and stage-signal regime, the one thing its rule 3 says is unmeasured: size at both controls and power per statistic, beside the CV-selected rank.
- Reuse the validated construction code by import so the only new construction is the phase-like nuisance and the table.
- Predeclare the size rule, the power floor and the verdicts before any paper-grade unit runs, and make the decision file carry the sentence the case study adopts.
- Keep every committed check's settings hash, records and decision untouched.

**Non-Goals:**

- Metabolomics-realistic blocks, block weighting or a feature count of 298. The structural analogue is the generator's.
- Any arm other than the primary as issued (no latent adjustment, residual correction or forced rank).
- A second strength, a second table (CN / LMCI / AD) or a conjunction rule.
- Changing `src/motco/`, the SEA-AD scripts' behaviour or any default.

## Decisions

### D1. A new script with a table registry, composing the parent constructions

`scripts/adni_design_check.py` imports `cs` and `d10`, builds the parent settings with `cs.build_parser().parse_args([...])` as the corrected-permutation check does, and overrides: `check = "adni_design"`, `base_seed = 2032`, `generator.group_stage_sizes` from `TABLES[args.table]`, the categorical's role (declared), the noise count (two), the phase settings, the cells and the rules. `TABLES` holds `adni` (`[[237, 310, 125], [223, 448, 158]]`, the ADNI2GO proportions per cell, the source `adni_sexdiff results/motco_analysis/result.json` and its cohort table's provenance). Tests register a small table in the same registry so the run path is exercised at n ≈ 36.

Subcommands: `calibrate --reps R --perms B --delta --s-z --s-p --out-dir --workers`, `run --shard i/N --out-dir [--delta --s-z --s-p --reps --perms --limit --n-jobs]`, `merge --out-dir`, `decide --out-dir [--run-dir]`.

*Alternative: flags on `case_study_type_i_check.py` (a `--table`, `--phase`).* Rejected: its settings hash identifies the committed 2026-10-29 run, its three arms and three-way decision are SEA-AD's, and a table whose n changes the timing by two orders of magnitude needs its own shard plan and tests.

*Alternative: a `design_grid` study config through `scripts/motco_study.py`.* Rejected: the harness has no covariate machinery, and the analysis path is what the case study runs.

### D2. Six design columns of ADNI's kinds

Declared covariates, in model order: `z` (numeric, D10's confounded covariate, injected at `s_z`), `category` (D10's three-level factor on the same linear predictor, declared with `categorical_covariates=["category"]`, drop-first `category[c1]`, `category[c2]`, not injected), `phase` (binary numeric 0/1, injected, D3), `x1`, `x2` (i.i.d. N(0, 1), not injected). That is six columns: one confounded numeric (ADNI's `age`, which rises with stage), two factor columns of a covariate associated with stage (ADNI's `apoe4`, whose allele counts rise from CN to AD), one batch-like column (`phase[ADNI2GO]`), two noise columns (`education`, `statin[1]`).

The factor is declared, not injected, because what ADNI's `apoe4` contributes to the model is collinearity of the design with stage, not a feature effect of its own that the check needs to resolve. At n ≥ 125 per cell and logit scale 0.5 every level is mixed within every cell, so no CV training set can lose a level (the n = 80 dropout mode does not arise, and no latent-adjusted arm runs anyway).

`phase` enters as a numeric 0/1 column rather than a declared factor: one column either way, and the numeric form keeps the LS-mean reference row at the pooled proportion as ADNI's does (0.544 there).

*Alternative: six numeric columns as in the SEA-AD check.* Rejected: SEA-AD enters six numerics; ADNI enters two numerics and four factor-level columns. The column count sets the df cost, but the factor's collinearity with stage is part of what is being measured.

### D3. The phase-like nuisance: binary, per-cell proportions, a dense offset in the expression block

For replicate r, `phase_i ~ Bernoulli(π_cell(i))` from `_seed(2032, "phase", r)` under the script's own namespace, with π the six recorded ADNI2GO shares (F 0.603 / 0.623 / 0.424, M 0.561 / 0.507 / 0.481, row order = the table's group order). The offset direction `d ∈ R^131` is drawn once from `_seed(2032, "phase_direction", 131)` as i.i.d. N(0, 1) (dense, no support subset), and the expression block receives `X_f += s_p · sd_f · d_f · phase_i` after D10's injection, with `sd_f` the feature's pooled SD as received. The methylation block is untouched. Each unit records the realized standardized difference between phases over the expression features (median |d| with the pooled within-phase SD) and the expression block's PC1 |corr| with `phase`, the two numbers ADNI's export recorded (1.12 and 0.73).

Why the expression block: ADNI's run offset is confined to p180, the smaller block (133 of 298). The generator's expression block has 131 features, so the offset lands on a block of the same size; methylation (367) stands in for NMR (165), which carries almost no offset. Why dense: 80 of the 133 p180 metabolites shift by more than one SD, so the offset is the block's first axis, not a sparse direction. Why per-cell proportions: in ADNI the phase differs by cell (SMC and EMCI were enrolled only in ADNI-GO/2), so a phase axis is also a stage axis; a batch effect independent of the design would be absorbed by the pooled centring and could never confound, and the covariate would do nothing in the model.

*Alternative: D10's categorical injection (`s_c > 0`).* Rejected: three levels on sparse 20 % supports over both blocks, which is neither binary, nor dense, nor one-block.

*Alternative: calibrate `s_p` to the realized |d| 1.12.* Rejected as the criterion: the regime that matters for the measurement is whether the basis follows the nuisance, which the first component's |ρ| reads; |d| is recorded beside it.

### D4. Calibration: three criteria, three knobs, anchor units only, two arms

`calibrate` evaluates the anchor units of R replicates (default 20, 99 permutations, `--workers` processes) in two arms on the same data and RRPP seed: `adjusted` (the six columns in the model) and `unadjusted` (no covariates, descriptive). Criteria, each recorded with `met` per iteration in `calibration.json`:

1. the `adjusted` arm's mean held-out CV AUROC inside `--auroc-band`, default `[0.55, 0.65]` (ADNI 0.600, 0.580–0.616 across its CV-selected analyses; the SEA-AD check used ±0.05–0.1 around its cohorts too);
2. the `adjusted` arm's median |Spearman| of the first latent component with `phase` at or above `--rho-threshold`, default 0.4 (ADNI 0.542; the covariate-blind basis is identical in both arms, so the adjusted arm reads it);
3. the `unadjusted` arm's anchor rejection rate above α on at least one statistic (the construction confounds).

Descriptives beside them: the selected-rank distribution (ADNI 8), the median |ρ| with `z`, the realized phase offset (median |d| and the block PC1 correlation), the pooled relative eigengap and the `angle` null q95.

Knobs and starting points. δ: at n = 1,501 the held-out AUROC for a given δ is higher than at n = 80 (more training rows per fold), so the SEA-AD values 0.05–0.1 are the upper end; start at δ = 0.03 and 0.05. `s_z`: start at D10's 3, not SEA-AD's 10, because `z` is confounded with stage and its injection adds stage-predictive signal; in ADNI the basis follows phase, not age. `s_p`: start at 1.5 (a median |d| near 1 against the pre-injection SD) and move it to meet criterion 2. The order is `s_p` for criterion 2, then δ for criterion 1, then confirm criterion 3. The user approves the calibrated strengths before the paper-grade run.

*Alternative: a criterion on the selected rank (ADNI 8).* Rejected: the rank is an output of the regime, recorded and stratified (D9), not a target to tune to.

### D5. Fourteen cells at the calibrated strength, one run

Cells: the anchor (e = 0), `translation` at e = 1.0, `magnitude` at {0.1, 0.25, 0.5, 1.0}, `orientation` and `shape` at {0.25, 0.5, 0.75, 1.0}. The geometry grid is the latent-adjustment and corrected-permutation power grid; the magnitude grid is extended to e = 1.0 (a doubling of δ) because at a weak δ the absolute stage-effect difference δ·e at e = 0.25 is small, and the design-guidance δ arm found no magnitude power at n = 200 for δ ≤ 0.1; at n = 1,501 the question is where, if anywhere, the curve rises. `cs.check_headroom` runs on every cell at the ADNI table before the first unit.

Type I and power sit in one run at one strength because the question is the operating characteristics at ADNI's regime, which is one regime. A second strength would double the cost for a bracket the SEA-AD checks already provide at small n.

### D6. One `adjusted` analysis per paper-grade unit

Each unit calls `analyze_trajectories` once: the six covariates (`categorical_covariates=["category"]`), `min_components = 2`, net displacement, the harness CV knobs, 999 permutations, the replicate's RRPP seed, `n_jobs = 1`, no latent adjustment, no residual correction. The record keeps the fields the spec lists plus the latent-score digest. A `TrajectoryAnalysisError` is recorded as the unit's error (unavailable), as in the parent checks.

The `unadjusted` arm runs at calibration only: at n = 1,501 a confounded construction rejects at the anchor near 1.0, so it carries no information past the non-vacuity criterion, and it would double the paper-grade cost. No consistency assertion is needed with one arm.

*Alternative: add the `latent_stage` arm (ADNI's sensitivity).* Rejected for cost and scope: ADNI's rule 1 never lets the arm replace the primary, the arm carries no claim of its own (its rule 8), and the SEA-AD diagnosis left `residualized_basis` without a MOTCO change. It is named as a follow-up if the exit review wants it.

### D7. Fresh seed family, every cell of a replicate matched

Base seed 2032. Every seed (generator, `z` and `category`, the injection direction, the noise, the phase draw and direction, RRPP) derives from it and the replicate index, so all 14 cells of replicate r share the baseline draw, covariates and permutation draws, and differ only in the surgery. Power curves are therefore paired across effects within a mode, and the two controls are paired, as in every Phase 6 check.

### D8. Size rule: the one-sided anchor bound at k = 2, six readings

For the anchor and the translation control, each of `delta`, `angle`, `shape`: rate ≤ α + 2·√(α(1 − α)/R) with R the available replicates (0.0695 at R = 500). Size holds when all six readings hold. This is the SEA-AD check's rule (`cs.anchor_bound`), under which a valid test fails at most one of six readings about 17 % of the time at exact α and far less at the conservative rates the corrected-permutation check measured. The corrected-permutation check's two-tier rule exists for 24+ readings over several runs and tables; with one run and six readings the per-reading rule is the precedent.

### D9. Power rule: the 0.80 floor, point estimate, rank strata beside it

For each target pair (magnitude/`delta`, orientation/`angle`, shape/`shape`) and effect: the rate with its binomial SE, the selected-rank distribution, and `informative_from` = the smallest effect at which the rate ≥ 0.80 (null when none). The 0.80 floor is the Phase 4, 5 and 6 gates' power floor, so the statement is comparable with the SEA-AD one (which reads 0.300 and 0.228 against it). The point estimate is used because that is how every gate read it; the SE (0.018 at 0.80, R = 500) is reported.

Beside it, descriptively: every statistic at every cell (off-target rates, which show specificity), and the rate per selected-rank bucket (2–3, 4–7, ≥ 8; ADNI's primary sits at 8 and its CV-selected arms at 4–8), because the net-displacement comparison found `angle` power governed by the retained rank at Phase 6 (0.280 at rank ≤ 3, 0.996 above).

*Alternative: the lower 95 % bound ≥ 0.80.* Rejected: no gate used it, and it would make the statement stricter than the one it replaces.

### D10. Verdicts, statements and the consequence

A pure function `adni_design_decision(rows, replicate_rows, settings)` writes `adni_design_decision.json`:

| verdict | condition | consequence (carried in the file) |
|---|---|---|
| `power_statement` | every size reading holds and at least one target statistic has an informative effect | `adni_sexdiff` replaces rule 3 with the measured statements in its own change; MOTCO documents them beside the SEA-AD statement |
| `non_rejections_uninformative` | every size reading holds and no target statistic reaches 0.80 at any effect | `adni_sexdiff` keeps rule 3, now citing the measured rates; the exit review reads it as a design without power |
| `diagnose_design_excess` | any size reading exceeds | no power statement; the roadmap names the exceeding readings and a diagnosis as the next step |

Per statistic the file carries a `statement` string: with an informative effect, "at the ADNI design `<stat>` rejects a `<mode>` difference of e = e* with power p (floor 0.80); a non-rejection argues against `<mode>` differences of e ≥ e*"; without one, "`<stat>` reaches at most p at e = e_max; a non-rejection is uninformative up to the construction maximum"; under size excess, the exceeding reading. Each statement names the CV-selected rank distribution it was measured at. The size sentence states the six readings against the bound. A run below 999 permutations marks the verdict indicative.

### D11. Records, layout and outputs

`results/adni-design-<date>/`:
- `calibration/calibration.json` (every iteration), `smoke/` (the R = 20 × 99 smoke of all 14 cells with its indicative decision);
- `run/settings.json`, `run/records/shard_<i>_of_<N>.jsonl`, and after `merge`: `adni_design_pvalues.csv` (unit × statistic, about 21,000 rows), `adni_design_replicates.csv` (per unit: rank, AUROC, |ρ| with phase and `z`, eigengap, `angle` null q95, realized offset, error), `adni_design_rates.csv` (cell × statistic: R, rejections, rate, SE, rank distribution), `adni_design_rank_strata.csv`;
- `decide`: `adni_design_decision.json` and `size_readings.csv`, beside the run;
- `PROVENANCE.txt` with the job, worktree, tree hash, versions, the sha256 of every file that defines the measurement and of every committed table, and what stays on `ing` (`records/`, about 30 MB, and any CSV over about 1 MB).

### D12. Run plan

1. **Calibrate locally** (`--workers`), iterating `s_p`, δ, `s_z` until the three criteria are met; every iteration lands in `calibration.json`.
2. **Smoke locally**: all 14 cells at R = 20 × 99 (280 units) with `merge` and `decide` (indicative), then R = 2 × 999 on the anchor, the translation and one power cell for timing. Append a dated "Calibration, smoke and shard plan" section here with the per-unit time.
3. **Approval gate.** The user approves the calibrated strengths, the replicate count and the array.
4. **On `ing`**, from a clean worktree of the pushed branch with matching tree hashes and the change's tests passing there: one array of N shards (N from the timing so that a shard stays well inside the 8 h limit at a 10× slow node), partition `512x1024`, n11 excluded, one CPU, BLAS threads 1, `scripts/adni_design_array.sbatch` (`ADNI_OUT`, `N_SHARDS`, `ADNI_ARGS`), `settings.json` pre-created with `--limit 0`.
5. **Finish:** `merge`, `decide`, `PROVENANCE.txt`, copy back with sha256 checked, the report, the docs.

Expected cost: the ADNI primary's 26 s at 9,999 permutations is about 16 s of double CV on 298 features plus 10 s of RRPP; at 498 features and 999 permutations a unit should take about 30–45 s on the workstation and, at the 2.5–3× ratio the SEA-AD checks measured, 80–120 s on an `ing` core: 150–230 core-hours for 7,000 units. If the timing exceeds about 250 core-hours the user may choose R = 300 at the gate (bound 0.0752, power SE 0.023 at 0.80); the rules do not change.

## Risks / Trade-offs

- [The three calibration criteria conflict: a phase offset large enough for |ρ| ≥ 0.4 is associated with stage and lifts the AUROC above 0.65 even at δ → 0] → ADNI's phase–stage association is modest (shares 0.42–0.62), so phase alone predicts stage weakly; lower δ first. If the band cannot be met, the user widens it or lowers the threshold before any paper-grade unit, recorded in `calibration.json` and in this design.
- [The weak-δ regime inflates the CV rank towards `max_components` 20 (the design-guidance δ arm read a median 11 at δ = 0.25)] → ADNI itself selected 8 at AUROC 0.60; the rank distribution is recorded, stratified (D9) and read beside every statement. A rank at the cap is a finding about the regime, not a failure of the check.
- [Cost: 7,000 analyses at n = 1,501] → Timed at smoke, approval gate, resumable shards, the R = 300 option at the gate.
- [A size reading exceeds by chance (≤ 17 % at exact α)] → The verdict is issued as predeclared; the report reads the reading against the SEA-AD rates; any follow-up is its own change.
- [The construction has power nowhere (`non_rejections_uninformative`)] → That is a result: it is the first measurement at this n and signal, and item 10 reads it as such.
- [Structural analogue wider than for SEA-AD: omic blocks for metabolomics, 498 for 298 features, block shares 74/26 for 45/55] → Stated as the first limitation in the report and in the case study's statement; no claim about metabolite-level behaviour is made.
- [`phase` as a numeric 0/1 column differs from ADNI's declared factor] → Identical design column and reference row; stated in D2.

## Migration Plan

Nothing migrates: no change under `src/motco/`, no default moves, every committed check keeps its settings hash and decision. The case study adopts the statement in its own change after the report.

## Open Questions

- The shard count and the replicate count at the gate follow the smoke timing (D12); neither changes a rule.
- Whether `adni_design_pvalues.csv` is committed depends on its size at merge (about 21,000 rows); otherwise it is named in `PROVENANCE.txt`.

## Implementation record

### Baseline (2026-11-02)

Starting state on a clean `main` at `b7f2d81` (MOTCO 0.13.0), before any code of this change: `ruff check src/ tests/` and `mypy src/motco/` pass; `MOTCO_TEST_PERMS=99 pytest tests/ -m "not slow"` gives 1,703 passed, 1 skipped (the R `shapes` fixture), 32 deselected; `tests/test_case_study_type_i_check.py` and `tests/test_corrected_permutation_check.py` give 52 passed. Work continues on branch `feat/adni-design-operating-characteristics`.

### ADNI reference values (verified 2026-11-02)

Read from `adni_sexdiff` at `49edb82`; every value the Context section quotes matches its record.

| value | record | source |
|---|---|---|
| cells F 237 / 310 / 125, M 223 / 448 / 158 (n = 1,501) | `design.cell_sizes` | `results/motco_analysis/result.json` (sha256 `d4a5cbf9…6f581f0`) |
| ADNI2GO per cell: F 143/237 = 0.603, 193/310 = 0.623, 53/125 = 0.424; M 125/223 = 0.561, 227/448 = 0.507, 76/158 = 0.481; pooled 0.5443 | counted per sex × stage | `results/motco_input/metadata.csv` (`9c445b7d…cecf221`) + `covariates.csv` (`06999e82…3fccfd`); pooled = `covariates.reference_row.phase[ADNI2GO]` 0.54430 |
| design columns `age`, `education`, `apoe4[1]`, `apoe4[2]`, `phase[ADNI2GO]`, `statin[1]` | `covariates.design_columns` | `result.json` |
| CV rank 8, held-out AUROC 0.5995 | `latent.selected_lv`, `latent.cv_mean_auroc` | `result.json` |
| AUROC 0.580–0.616 across the CV-selected analyses | 0.580311, 0.598348, 0.599812, 0.608111, 0.616134 | `reports/motco-adni-2026-11-02.md` |
| first latent component ρ 0.5424 with `phase` (`pls_3` −0.5346) | Spearman, `axis_correlations.csv` | `results/motco_diagnostics/axis_correlations.csv` (`d0406d19…904afa5`) |
| p180 median \|d\| 1.12, 80 of 133 above 1, PC1 \|corr\| 0.73; NMR 0.21, 0, 0.13 | export `counts.phase_offset` | `results/motco_analysis/summary.md` |
| null `angle` q95 35.44° | `null_summary.angle.q95` | `result.json` |
| pooled relative eigengap 0.7440 | `config_spectrum.pooled.relative_eigengap` | `result.json` |
| runtime 26.0 s at 9,999 permutations; blocks p180 133, NMR 165 | `measurement.runtime_seconds`, `blocks` | `result.json` |

The script registers these in `TABLES["adni"]` (the shares as exact count ratios, the three source files' sha256 and the reference values, descriptive only). Its realized-offset measure mirrors ADNI's export (`adni_sexdiff.data._phase_offset`: per-feature mean difference over the pooled within-phase SD, population variances; PC1 of the column-standardized block, |Pearson| with the phase indicator).

### Calibration (2026-11-02)

Twelve iterations of `calibrate --reps 20 --perms 99 --workers 10` on the workstation, all appended to `results/adni-design-2026-11-02/calibration/calibration.json` (no failures in either arm). AUROC is the `adjusted` arm's mean held-out CV AUROC, ρ the median first-component |Spearman| (`phase` for the criterion, `z` beside it), "unadj." the `unadjusted` arm's anchor rates (`delta` / `angle` / `shape`).

| δ | `s_z` | `s_p` | AUROC | ρ phase | unadj. | met | rank (adjusted) | ρ `z` | eigengap | `angle` null q95 | median \|d\| (\|d\| > 1 of 131) |
|---|---|---|---|---|---|---|---|---|---|---|---|
| 0.03 | 3 | 1.5 | 0.774 | 0.242 | 1.00 / 0.20 / 0.50 | no | 2 (17), 3 (2), 19 (1) | 0.941 | 0.958 | 163.3° | 0.969 (63) |
| 0.05 | 3 | 1.5 | 0.877 | 0.241 | 1.00 / 1.00 / 1.00 | no | 19 (20) | 0.942 | 0.151 | 53.3° | 0.969 (63) |
| 0.03 | 1 | 1.5 | 0.762 | 0.490 | 1.00 / 0.70 / 0.70 | no | 3 (16), 19 (3), 2 (1) | 0.572 | 0.833 | 93.3° | 0.969 (63) |
| 0.03 | 0.3 | 1.5 | 0.747 | 0.496 | 0.65 / 0.05 / 0.80 | no | 19 (19), 18 (1) | 0.123 | 0.266 | 48.8° | 0.969 (63) |
| 0.015 | 0.3 | 1.5 | 0.668 | 0.496 | 0.65 / 0.05 / 0.70 | no | 19 (13), 5–18 (7) | 0.120 | 0.379 | 51.9° | 0.969 (63) |
| 0.01 | 0.3 | 1.5 | 0.652 | 0.499 | 0.65 / 0.10 / 0.75 | no | 19 (8), 4–18 (12) | 0.121 | 0.433 | 53.4° | 0.969 (63) |
| 0.005 | 0.3 | 1.5 | 0.645 | 0.502 | 0.65 / 0.10 / 0.75 | yes | 5–19, median 12.5 | 0.121 | 0.475 | 55.3° | 0.969 (63) |
| 0.01 | 0.1 | 1.5 | 0.549 | 0.432 | 0.20 / 0.00 / 0.25 | no | 2 (8), 19 (9), 3–18 (3) | 0.076 | 0.531 | 61.4° | 0.969 (63) |
| 0.02 | 0.1 | 1.5 | 0.602 | 0.425 | 0.20 / 0.00 / 0.40 | yes | 19 (20) | 0.077 | 0.334 | 50.6° | 0.969 (63) |
| 0.025 | 0.1 | 1.5 | 0.637 | 0.421 | 0.20 / 0.00 / 0.35 | yes | 19 (20) | 0.078 | 0.301 | 49.9° | 0.969 (63) |
| **0.02** | **0.1** | **1.75** | **0.602** | **0.500** | **0.15 / 0.00 / 0.45** | **yes** | **19 (20)** | **0.080** | **0.328** | **51.5°** | **1.127 (72.5)** |
| 0.02 | 0.1 | 2.5 | 0.604 | 0.707 | 0.15 / 0.00 / 0.35 | yes | 19 (20) | 0.092 | 0.323 | 51.6° | 1.611 (89) |

What the iterations showed:

- **`z`, not phase, carried the excess stage signal.** At D10's `s_z = 3` the first component follows `z` (|ρ| 0.94) and the AUROC sits at 0.77–0.88. Lowering `s_z` hands the first component to `phase` (|ρ| 0.49–0.50 from `s_z = 1` down), and at `s_z = 0.3` it leaves an AUROC floor near 0.645 that δ barely moves (0.668 / 0.652 / 0.645 at δ = 0.015 / 0.01 / 0.005). The first met iteration, δ = 0.005 at `s_z = 0.3`, sits on that floor: its stage signal is `z`'s, and δ contributes almost nothing. At `s_z = 0.1` the floor drops (AUROC 0.549 at δ = 0.01), which leaves room for a real stage effect.
- **The phase offset adds no stage signal.** Raising `s_p` from 1.5 to 2.5 at δ = 0.02 moves the AUROC from 0.602 to 0.604 while the phase |ρ| rises from 0.43 to 0.71. So the offset steers the basis without making stage easier to predict, as in ADNI, where `phase` leads the first component at an AUROC of 0.600.
- **`s_z` moved outside the design's starting point.** D4 named `s_z` as a knob and expected it to come down from D10's 3. At `s_z = 0.1` the injected `z` effect is small (|ρ| with `z` 0.08, ADNI's first component against `age` −0.06), and `z` still enters the model as a confounded column. The unadjusted arm still rejects (criterion 3), on `delta` and `shape`, which the phase shares' group differences alone account for (F 0.603 / 0.623 / 0.424 against M 0.561 / 0.507 / 0.481).

**Selected strengths: δ = 0.02, `s_z` = 0.1, `s_p` = 1.75** (now the script's defaults). Of the five met iterations, it is the one that matches ADNI on both criterion quantities and on the realized offset. AUROC is 0.602 (ADNI 0.600) and the phase |ρ| is 0.500 (ADNI 0.542). The expression block's median |d| is 1.127 with 72.5 of 131 features above 1 (ADNI's p180: 1.12, 80 of 133). Its unadjusted arm rejects `delta` at 0.15 and `shape` at 0.45.

Recorded beside ADNI's values, not judged:

- **Selected rank.** It is 19, the cap (`max_components` − 1), in all 20 replicates; ADNI's is 8. This is the risk D4/Risks names: at a weak diffuse stage signal over 498 features, the held-out AUROC keeps rising with every component.
- **Eigengap and null width.** The pooled relative eigengap is 0.33 (ADNI 0.744), and the `angle` null q95 is 51.5° (ADNI 35.4°).
- **Phase share of the block's PC1.** The offset accounts for the expression block's PC1 almost entirely (|corr| 0.994; ADNI's p180 0.73). The generator's block has no other structure as strong as real metabolomics.
- **The untouched block.** Methylation reads median |d| 0.03 and PC1 |corr| 0.02 (ADNI's NMR 0.21 / 0.13).

Every statement the check issues therefore describes the regime at rank 19, and the report must say so beside ADNI's 8.
