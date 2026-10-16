# Tasks

## 1. Byte-identity baseline (before any code change)

- [x] 1.1 Add a capture helper in `tests/test_analysis_identity.py` that builds a seeded two-block, three-stage fixture, runs the 0.7.0 `analyze_trajectories` (99 permutations, fixed seed, attribution off), and writes `tests/data/analysis_no_covariate_identity.json`. The file holds the record without `measurement.runtime_seconds` and `software`, plus SHA-256 digests of `latent_scores.csv`, `pls_loadings.csv` and `stage_means.csv`. Verify by committing the fixture from an unmodified tree and confirming the comparison test passes on it.

## 2. Covariate design helper

- [x] 2.1 Add `CovariateDesign` and `build_covariate_matrix(frame, columns, categorical)` to `src/motco/stats/design.py`: numeric columns as-is, factors as drop-first indicators ordered by `_sort_levels`, `name[level]` column names, kinds, levels, reference levels and the column-mean reference row. Verify with unit tests for numeric + categorical columns, a numeric-coded factor, column naming and order, and the reference row.
- [x] 2.2 Make the helper reject a constant numeric covariate, a one-level factor and non-finite numeric values, each with an error naming the covariate. Verify with parametrized tests.

## 3. Dataset declarations and the missing-covariate policy

- [x] 3.1 Add `covariate_cols`, `categorical_covariates` and `missing_covariates` to `TrajectoryDataset`, and validate them: columns present, unique, distinct from the sample, group and stage columns; categorical columns a subset of the covariates; policy in `{"error", "complete_case"}`. Verify with rejection tests, including the spec's group-column scenario.
- [x] 3.2 Implement the `error` policy: name the column and up to ten sample ids. Implement `complete_case`: subset the metadata and every block before the existing validators run, store the dropped ids, and name any group × stage cell that deletion emptied. Verify with the spec scenarios (default error; 78 of 80 retained, with cell sizes summing to 78 and the dropped ids listed; emptied-cell error).

## 4. Covariate-adjusted measurement

- [x] 4.1 In `analyze_trajectories`, raise before any fitting when attribution is requested with covariates. Verify with a test that monkeypatches `fit_pls_latent_space` to fail if it is called.
- [x] 4.2 Append the covariate columns to `model_full` and `model_reduced`, and append the reference row to every LS-mean row. Check full column rank and raise an error naming the covariates. Pass the augmented design to `estimate_difference`, `RRPP` and the stage-means frame. Verify with the collinear-covariate scenario and with a test that the latent scores, loadings and selected rank are identical with and without covariates.
- [x] 4.3 Add the property tests from the spec: (a) the statistics do not depend on the reference row, and the stage means shift by one common vector; (b) on a constructed latent matrix with cell means plus `z·v` and confounded `z`, the adjusted statistics match the cell-mean-only ones while the unadjusted ones differ; (c) the same seed gives identical adjusted p-values and null summaries.
- [x] 4.4 Emit the `covariates` record block (D7) only when covariates are declared, after `design`, with the coefficients from the observed full-model fit. Verify with a record-layout test and with the byte-identity test from 1.1 on the no-covariate path.

## 5. CLI

- [x] 5.1 Add `--covariate` and `--categorical-covariate` (repeatable, in flag order, the latter implying the former) and `--missing-covariates` to `motco analyze`. Wire them into `TrajectoryDataset`, and print the covariates and the dropped count when they are set. Verify with `tests/test_analysis_cli.py` cases for the adjusted run (the result.json covariate block) and a missing covariate column (non-zero exit naming it).

## 6. Type I check at the Phase 6 design

- [x] 6.1 Write `scripts/covariate_type_i_check.py` (D10): Phase 6 generator settings at e = 0, M-value methylation, the confounded numeric and categorical covariates, fixed-direction injection with recorded strengths, paired arms with a per-replicate latent-identity assertion, `--shard`, resume and `merge`. Verify with a fast test that runs 2 replicates × 9 permutations twice and gets identical records, and that checks the recorded per-cell covariate means show a non-zero interaction contrast.
- [x] 6.2 Add the decision step that writes `covariate_type_i_decision.json` and `covariate_type_i_rates.csv` (per-arm rates, the bound, R, per-statistic verdicts, `covariates_validated` | `diagnose_covariate_excess`). Verify with a unit test on synthetic p-value tables that covers both outcomes and that checks the unadjusted arm does not affect the decision.
- [x] 6.3 Add `scripts/covariate_type_i_array.sbatch` modeled on `orientation_probe_array.sbatch`. Verify by running the smoke profile (R = 20, 99 permutations) locally and confirming that the unadjusted arm rejects above α for at least one statistic. If it does not, raise the strengths, record them, and re-run the smoke profile.
- [x] 6.4 Run the paper-grade check (R = 500, 999 permutations) on ing, merge the shards, and write `docs/reports/covariate-type-i-<date>.md` with the rates, bound and decision. Verify that the decision file exists and that the report cites it.

## 7. Docs, version and gate

- [x] 7.1 Update `docs/api/analysis.md` (covariate contract, the missing policy, the covariate-blind latent space and its exactness argument, reference-row invariance, the attribution exclusion, parsimony advice) and the CLI reference. Update the `analysis.py` paragraph in `AGENTS.md` ("No covariates" becomes the covariate contract). Update `docs/roadmap.md`: mark item 2 done with the Type I decision, and lift the case-study wait in item 1 conditionally on that decision. Verify with `uv run mkdocs build --strict`.
- [x] 7.2 Bump `src/motco/__init__.py` to 0.8.0. Verify that `motco --version` (or the package version) reports 0.8.0.
- [x] 7.3 Run the pre-commit gate (`ruff check`, `mypy src/motco/`, fast pytest), including the unchanged `test_matches_harness_*`, the committed config identity tests and the orientation-default pins. Verify that all three pass.
