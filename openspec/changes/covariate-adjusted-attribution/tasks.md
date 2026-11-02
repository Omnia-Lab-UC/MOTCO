# Tasks

## 1. Baseline on the unmodified tree

- [x] 1.1 Confirm the starting state, and record the package version the change starts from (0.13.0). Verify that these pass on a clean `main`:
  - the pre-commit gate (`uv run ruff check src/ tests/ && uv run mypy src/motco/ && MOTCO_TEST_PERMS=99 uv run pytest tests/ -m "not slow" --tb=short`);
  - `uv run pytest tests/test_attribution.py tests/test_attribution_diagnostics.py tests/test_evaluation_attribution.py tests/test_analysis.py tests/test_analysis_identity.py -q`.
- [x] 1.2 Write `tests/test_attribution_identity.py` (design D3). It is both the test and, run as a script, the capture of the new fixture `tests/data/attribution_no_covariate_identity.json`. The fixture holds:
  - the primitive's frames on a fixed synthetic dataset and fitted model, with 25 bootstrap replicates;
  - the attribution block, minus runtimes, and the attribution CSVs of an `analyze_trajectories` run with attribution on the `test_analysis_identity` dataset.

  Capture it **before any edit under `src/`**, and commit the fixture on its own. Verify:
  - the test passes against the fixture on the unmodified tree;
  - `git log` shows the fixture commit before any `src/` change.

## 2. Covariate-adjusted means and bootstrap in the attribution primitive

- [x] 2.1 Add the keyword-only `covariates` and `covariate_reference` to `analyze_orientation_attribution` (design D1, D4). This covers:
  - validation: row count, finite values, unique columns, and reference columns and values;
  - mutual exclusion with `mean_table`/`precomputed_means`;
  - adjusted means `x̄_gs − γ̂ᵀ(c̄_gs − c_ref)`, with γ̂ from `motco.stats.pls.fit_nuisance_coefficients` on the drop-first group × stage indicators, and its rank error re-raised as an `AttributionError` naming the columns;
  - `mean_source = "covariate_adjusted"`;
  - the optional `AttributionConfig` and `InterpretationMetadata` fields, emitted as rows by `attribution_frames` only when set.

  Verify with a new `tests/test_adjusted_attribution.py`:
  - a pure cell-mean structure plus `z·v` gives adjusted transitions, path lengths and both contrast blocks equal to the structure's, while the unadjusted contrasts differ;
  - the means equal an independent `lstsq` fit of `X ~ cells + C` evaluated at `c_ref`, to 1e-10;
  - the contrasts are identical at two reference rows, and the means differ by one common vector;
  - a covariate constant within every cell raises an error naming it;
  - each invalid input raises;
  - `mean_table` together with `covariates` raises;
  - the configuration and interpretation rows are present only with covariates.

  `tests/test_attribution.py` and `tests/test_attribution_identity.py` must pass unchanged.
- [x] 2.2 Implement the adjusted bootstrap (design D2). It keeps the existing draw loop, and per replicate it:
  - forms `X*` and `C*` from the drawn rows;
  - refits γ̂* and evaluates at the full-sample `c_ref`;
  - turns a rank failure into an invalid replicate (NaN rows in every block and component) and increments `rank_deficient_replicates`.

  The preprocessing and model stay frozen. Verify in `tests/test_adjusted_attribution.py`:
  - with the same seed, the adjusted and unadjusted bootstraps draw identical row indices (a recording wrapper around the generator);
  - two adjusted runs give identical summaries;
  - the coefficient fit is called once per replicate on that replicate's rows (spy), so the adjustment is refitted;
  - a categorical indicator that is non-zero in one sample produces rank-deficient replicates, which are excluded from every valid count and recorded in the configuration;
  - an exact cell structure plus `z·v` gives sign stability 1.0 for every non-zero feature.
- [x] 2.3 Document the keywords and the adjustment in the module docstring and in the `analyze_orientation_attribution` docstring, and add a "Covariate-adjusted means" section to `docs/api/attribution.md`. It covers:
  - the formula and the reference-row invariance;
  - the bootstrap refit;
  - the rank-deficient replicates;
  - the interpretation boundary.

  Verify that `uv run mkdocs build --strict` reports no warnings beyond those already present on `main`.

## 3. The compact record, the analysis and the CLI

- [x] 3.1 Pass `covariates` and `covariate_reference` through `run_attribution_diagnostics` (and `compute_attribution_diagnostics`). Add the `adjustment` block (`mean_source`, `design_columns`, `reference_row`, `rank_deficient_replicates`, `bootstrap_replicates`) only when they are given, and leave `ATTRIBUTION_SCHEMA_VERSION` at 2 (design D3, D4). Verify:
  - a test shows the block present with covariates and absent without;
  - `tests/test_attribution_diagnostics.py` and `tests/test_evaluation_attribution.py` pass unchanged.
- [x] 3.2 Make three changes in `src/motco/analysis.py` (design D5, D6):
  - replace the covariates guard with a latent-adjustment guard, raised before fitting;
  - pass the `CovariateDesign` (matrix with `column_names`, and `reference_row`) to `_run_attribution`;
  - add `TrajectoryAnalysisResult.latent_artifacts` (in memory; excluded from the record, `to_json` and `write_analysis_outputs`).

  Update the `analyze_trajectories` docstring. Replace `tests/test_analysis.py::test_attribution_with_covariates_is_rejected_before_fitting` with a latent-adjustment rejection test, in which `fit_pls_latent_space` is not called. Add to `tests/test_adjusted_attribution.py`:
  - `model.transform(attribution observed means)` equals `stage_means` to 1e-8 relative on a covariate-adjusted fixture, including one with complete-case deletion;
  - the adjustment block's `reference_row` equals the covariate block's;
  - the analysis with and without attribution gives identical latent scores, rank, statistics, p-values, null summaries and stage means;
  - `latent_artifacts.model.transform(latent_artifacts.features)` equals the latent scores to within 1e-10 of the largest absolute score;
  - the written file set is the same with and without the field.

  `tests/test_analysis_identity.py` and `tests/test_attribution_identity.py` must pass unchanged.
- [x] 3.3 Print `Attribution: covariate-adjusted means at the reference row; <k> of <B> bootstrap replicates rank-deficient` in `motco analyze` when the attribution record carries `adjustment` (design D5). Verify with new cases in `tests/test_analysis_cli.py`:
  - a covariate-adjusted run with `--attribution --bootstrap-replicates 20 --bootstrap-seed 3` exits zero, `attribution/configuration.csv` states `covariate_adjusted`, `result.json` holds the adjustment block, and stdout holds the line;
  - `--latent-adjustment stage --attribution` exits non-zero with the latent-adjustment message;
  - a run without covariates prints no adjustment line.
- [x] 3.4 Update the documentation:
  - in `docs/api/analysis.md`, replace "Attribution is not available with covariates" with the adjusted-attribution behaviour, the adjustment block, the latent-adjustment rejection, `latent_artifacts` (in memory only) and the validation status "pending the adjusted-attribution check";
  - regenerate the `motco analyze --help` block in `docs/api/cli.md` if it changed;
  - in the `analysis.py` paragraph of `AGENTS.md`, state that attribution with covariates uses adjusted means and the refitted bootstrap, that it is rejected with latent adjustment, and that `latent_artifacts` exists.

  Verify that `uv run mkdocs build --strict` passes as in 2.3, and read the three documents against the specs.
- [x] 3.5 Bump `src/motco/__init__.py` to 0.14.0 before any paper-grade unit runs (design D12). Verify that `python -c "import motco; print(motco.__version__)"` prints 0.14.0 and that the pre-commit gate passes.

## 4. The validation script

- [x] 4.1 Make `covariate_type_i_check.build_replicate` also return the generated dataset under a new key, `"dataset"`, which no record writer reads (design D7). Verify that these pass unchanged: `tests/test_covariate_type_i_check.py`, `tests/test_case_study_type_i_check.py`, `tests/test_case_study_diagnosis.py`, `tests/test_corrected_permutation_check.py` and `tests/test_adni_design_check.py`.
- [x] 4.2 Write the settings and construction of `scripts/adjusted_attribution_check.py` (design D7, D12). It imports `case_study_type_i_check as cs` and `covariate_type_i_check as d10`, and covers:
  - the astrocyte table, δ = 2.0, base seed 2033, `s_z` (default 10), five noise covariates, the three cells, R 500, 999 permutations, B 200 and top-k 20;
  - the bootstrap seed `_seed(base, "bootstrap", replicate)`;
  - `cs.check_headroom`, the hashed `settings.json`, and refusal of a resume under other settings.

  Verify with a new `tests/test_adjusted_attribution_check.py` (small CV, two replicates, nine permutations, B = 5):
  - the cells and defaults are as stated;
  - the seeds depend only on the base seed and the replicate, are shared across a replicate's cells, and differ from the committed families' at the same replicate;
  - a resume under another `s_z`, base seed, permutation count, B or top-k is refused;
  - an over-headroom cell refuses to start.
- [x] 4.3 Write the unit evaluation (design D8, D9). It consists of:
  - the `adjusted` arm through `analyze_trajectories`;
  - the oracle features `X̃ − z ⊗ (v / σ)` and the `oracle` and `unadjusted` arms through the primitive, with the same model, bootstrap seed, B and top-k;
  - the two consistency checks;
  - the metrics per arm, component and block (truth precision over min(k, |truth|), contamination, projection on ŵ, cosine and overlap with the oracle, mean sign stability and top-k frequency over the arm's own top-k, valid and rank-deficient counts);
  - a record line holding the unit key, seeds, truth size and the analysis's statistics, p-values, rank and CV AUROC;
  - a sharded, resumable `run` built on the parent's shard helpers.

  Verify with tests:
  - on a fixture unit, the oracle's observed cell means equal the uninjected blocks' standardized cell means up to one constant vector;
  - a monkeypatched stage-means mismatch raises and writes no line;
  - the metrics on hand-built contrast vectors match hand counts;
  - resume skips recorded units;
  - two shards cover every unit exactly once.
- [x] 4.4 Write `calibrate` (design D11). Each iteration runs R = 20 × 99 permutations × B = 50 at the anchor and orientation e = 1.0, computes the contamination and recoverability criteria, and appends to `calibration.json`. `s_z` steps by 5 until contamination is met, and the step stops on a recoverability failure. Verify with tests on synthetic records: met, unmet, stop-on-recoverability, and iterations appended rather than overwritten.
- [x] 4.5 Write `merge`. It produces `adjusted_attribution_metrics.csv` (per unit × arm × component × block) and `adjusted_attribution_summary.csv` (per cell × arm × component × block: mean and SE of every metric and of every paired difference against the oracle). Verify with tests: means and paired SEs on synthetic records match a hand computation, and the merge is idempotent.
- [x] 4.6 Write the pure `adjusted_attribution_decision` and `decide` (design D10). They cover:
  - both gates and the four criteria on the observed principal block, at 2 paired SEs;
  - the verdict with its consequence, the settings hash and the adopted calibration iteration;
  - the descriptive unadjusted arm, the captured and residual components, the transition blocks and the readings conditional on `angle` rejection;
  - `indicative` below 999 permutations, B = 200 or R = 500.

  Verify with tests on synthetic rows reproducing the spec scenarios:
  - contamination fails at 0.04 + 2 × 0.01;
  - residual confounding holds at 0.02 + 2 × 0.005 against 0.10 × 0.40;
  - recovery holds at −0.03 − 2 × 0.005;
  - top-k frequency fails at −0.10 − 2 × 0.01;
  - a weak gate gives `construction_uninformative`;
  - `adjusted_attribution_validated` and `diagnose_adjusted_attribution` each name their cells;
  - a smoke decision is indicative.

## 5. Calibration, smoke, timing and approval

- [x] 5.1 Run `calibrate` locally into `results/adjusted-attribution-<date>/calibration/`, and record every iteration. If the recoverability criterion fails, stop and present the iteration to the user before anything else. Verify that `calibration.json` names an iteration meeting both criteria, or that the user's decision is recorded in `design.md`.

  Under the declared rule calibration ended `stop_exhausted` on 2026-11-02 (iterations 0–8: no contamination at orientation e = 1.0 at any `s_z`). The user amended the gates (design D10; rule version 2: contamination gated at the anchor, a projection gate at every cell). Record the diagnostic under `calibration/probe/`, update the script, its tests and the spec, and rerun `calibrate` once under rule version 2 from `s_z = 10`, appending to the same `calibration.json`.
- [x] 5.2 Run the smoke at the calibrated `s_z`: R = 20 × 99 permutations × B = 50 over the three cells into `smoke/`, then `merge` and `decide`. Time two units at 999 permutations × B = 200. Verify: zero consistency failures, an indicative decision file, and the 999 × 200 units completing end to end.
- [x] 5.3 Append a dated "Smoke and run plan" section to `design.md`. It records the calibration outcome, the smoke metrics, the per-unit timing, the planned core-hours, and local workers versus `ing` under D11's 12-core-hour rule; for `ing`, also `scripts/adjusted_attribution_array.sbatch` modelled on `scripts/corrected_permutation_array.sbatch`, with `bash -n` passing. Present it to the user and obtain approval. Verify that the section and the recorded approval exist in `design.md`.

## 6. Paper-grade run

- [x] 6.1 Push the branch and create a clean worktree on the approved host, whose tree hash equals the local committed tree. Run `tests/test_adjusted_attribution.py` and `tests/test_adjusted_attribution_check.py` there. Verify: identical tree hashes and passing tests.
- [x] 6.2 Run the 1,500 units at 500 replicates × 999 permutations × B = 200 into `results/adjusted-attribution-<date>/run/`, monitored to completion, with any stalled shard resubmitted. Verify: 500 units per cell, zero consistency failures, and every array task `COMPLETED 0:0` (or every local worker exiting zero).
- [x] 6.3 Run `merge` and `decide`, and copy back the metric summary, the decision and the calibration. Commit each CSV under about 1 MB, and otherwise name it in provenance. Write `PROVENANCE.txt`: the commands, jobs or workers, tree hash, settings hash, versions (MOTCO 0.14.0), file sha256s and core-hours. Verify that the decision's settings hash equals `settings.json`'s and that its verdict is not indicative.

## 7. Report, roadmap and docs

- [ ] 7.1 Write `docs/reports/adjusted-attribution-<date>.md`. It covers:
  - the construction;
  - the calibration with every iteration;
  - per cell, the gates and the four criteria with their bounds, SEs and R;
  - the verdict and its consequence for item 12;
  - the descriptive readings: the unadjusted arm, the captured and residual components, the transitions, the `angle`-rejection subset and rank-deficient counts;
  - the limitations: numeric covariates only, the astrocyte table only, the InterSIM blocks (498 features), δ = 2.0 rather than SEA-AD's stage signal, and the additive linear covariate model.

  Verify that every number is read from the committed decision, summary and calibration files, and that the report names them.
- [ ] 7.2 Update the documents with the verdict:
  - `docs/roadmap.md`: item 11 marked done with its verdict and report link, item 12's start condition, and the Phase 6 "Not yet established" bullet;
  - `docs/api/analysis.md`: the validation status;
  - `AGENTS.md`: the check (script, verdict, date) in the `analysis.py` paragraph;
  - `design.md`: an appended "Paper-grade run" record.

  Verify with `uv run mkdocs build --strict` as in 2.3, and by reading the documents against the decision file.

## 8. Gate

- [ ] 8.1 Run the pre-commit gate plus `uv run ruff check scripts/adjusted_attribution_check.py`. Confirm three things:
  - all checks pass;
  - `git diff --stat main -- tests/data` lists only the new `attribution_no_covariate_identity.json`;
  - the version is 0.14.0.

## Workflow follow-up

- Open the pull request, merge, tag `v0.14.0`, sync the `orientation-driver-attribution`, `real-data-trajectory-analysis` and `adjusted-attribution-validation` specs into `openspec/specs/`, and archive the change.
- Under `adjusted_attribution_validated`: propose item 12, the plasmode demonstration on SEA-AD astrocyte `cohort_30`, pinned to 0.14.0. Its proposal fixes the detection rate, the recovery thresholds and the null-making step.
- Under `diagnose_adjusted_attribution`: propose a diagnosis of the failing criterion before item 12. The attribution leg stays withheld.
- Under `construction_uninformative`: return the construction to the user. The leg stays withheld.
