# Tasks

## 1. Orthonormal coordinates in the harness

- [x] 1.1 Add `SimulationEvaluationParams.latent_coordinates` (default `None`; `"scores"` normalized to `None`; `"orthonormal"`; an unknown value rejected in `__post_init__` naming it and the allowed values). Add it to `grid._SIGNATURE_OMIT_WHEN_NONE`, and write it to `runtime_metadata` only when set. Verify in `tests/test_orthonormal_coordinates.py` the normalization and the rejection. Also verify that `tests/test_committed_config_identity.py`, `tests/test_orientation_default_pins.py`, `tests/test_rank_floor_default_pins.py` and `tests/test_post_floor_config_identity.py` pass unchanged.
- [x] 1.2 Add the keyword `coordinates` (default `"scores"`) to `fit_pls_latent_space` (design D1, D3). Under `"orthonormal"`:
  - after the final fit, the matrix becomes `inputs(X)·Q`, with `Q = qr(x_rotations_)`, on the training rows of the latent-adjusted path too;
  - `integration_params.latent_coordinates` is recorded;
  - the returned map's `transform` is `inputs(rows)·Q` (a `LinearLatentMap` constructor for the orthonormal coordinates in `diagnostics.py`).

  The scores path keeps `x_scores_` and `LinearLatentMap.from_pls`, byte for byte. Verify in the same test file that:
  - the orthonormal matrix equals `X_in·x_weights_` up to column signs and spans the scores' subspace (1e-8);
  - the rank and CV table are identical between coordinate systems;
  - `tests/test_analysis_identity.py` and the core/harness identity fixtures pass unchanged.
- [x] 1.3 Pass the setting from the harness integration to the core and reject it before any preprocessing or fit (design D4) in three cases, each error naming the offending value:
  - with a method other than `pls`;
  - with attribution;
  - with an unknown value.

  Verify one test per rejection, and that `snf`/`pca`/`concat`/`pls_within_group`/`pls_residual_pca` evaluations without the setting are byte-identical to before.
- [x] 1.4 Verify the measurement-level contracts in the same test file:
  - rotating the orthonormal matrix by a random orthogonal matrix leaves `delta`/`angle`/`shape` within 1e-10 and the RRPP p-values equal (design D2);
  - under `span_capture`, `population_latent` is in the orthonormal coordinates and the `span_capture` shares equal the scores' exactly;
  - the `pls_latent` checkpoint and configuration spectrum come from the orthonormal matrix;
  - a forced rank works under the setting;
  - an evaluation without the setting, or with `"scores"`, is byte-identical to the pre-change result on a pinned dataset.
- [x] 1.5 Add `pls_orthonormal_from_pls(dataset, pls, params)` (spec: prebuilt fit). It refuses a non-`pls` integration, one without artifacts, or one already in orthonormal coordinates. `evaluate_semisynthetic_trajectory(..., latent=)` refuses a prebuilt integration whose coordinates differ from its params, naming both. Verify with tests that:
  - the reuse path equals a from-scratch fit (matrix and metadata) and fits no PLS model (patched `fit_plsda_model` call count);
  - each refusal fires.
- [x] 1.6 Document the setting in `docs/api/simulations.md` (parameter, scope, rejections, invisibility while absent, the prebuilt path) and in the `evaluation.py` and `grid.py` bullets of `AGENTS.md`. Verify that `uv run mkdocs build` adds no warnings beyond those on `main`.

## 2. Validation driver: units and records

- [x] 2.1 Create `scripts/orthonormal_coordinates_check.py` with `run --profile paper_grade|weak`, `merge` and `decide` (design D8):
  - `settings.json` bound by hash;
  - round-robin shards;
  - resume of complete unit keys;
  - refusal of a resume or merge under other settings, naming the keys;
  - the parents' sha256 recorded, with a refusal naming any changed file.

  Verify in `tests/test_orthonormal_coordinates_check.py`:
  - the unit counts (65 cells × 500 = 32,500; 2 × 500 = 1,000);
  - resume;
  - the settings-mismatch refusal;
  - a refused profile whose hash differs.
- [x] 2.2 Implement the `paper_grade` unit builder (design D6) over `stage_span_leak_check.column_cells` / `unit_params`, imported unchanged, with `span_scale="planted"`, `span_capture=True` and 999 permutations. Implement the `weak` builder over `span_split_check` / `latent_space_comparison`'s `seaad_weak` anchor and translation units. Verify on the committed profiles that, for sampled keys:
  - the generator parameters, dataset seed and RRPP seed equal the leak check's and the span-split check's for the same key;
  - the weak δ equals the committed calibration's.
- [x] 2.3 Implement the unit (design D5):
  - generate once, fit `pls` once, build the orthonormal integration with `pls_orthonormal_from_pls`, and evaluate both arms serially at the unit's RRPP seed;
  - run the two consistency checks and record their deviations;
  - record a failure without aborting the shard.

  Verify on a small design that:
  - the `scores` record equals a plain production evaluation of the same unit;
  - the `orthonormal` record equals an evaluation fitted from scratch with the setting;
  - both arms share rank, AUROC and capture;
  - the record carries every field the spec names;
  - an inconsistent unit (patched basis) is flagged with the check and its deviation.
- [x] 2.4 Add `scripts/orthonormal_coordinates_array.sbatch` (`ORTHO_PROFILE`; absolute paths, BLAS pinned to one thread, a throttled array, `PYTHONPATH` on the worktree's `src`). Verify with `bash -n` and a dry run of its command lines against a temporary directory.

## 3. Merge and reproduction

- [x] 3.1 Implement `merge`:
  - one rate row per profile, column or design, cell, statistic and arm (counts at p < 0.05, rate, Monte Carlo SE);
  - the paired table (both rates, difference, discordant counts b and c over units available in both arms, exact two-sided paired p-value);
  - the per-cell medians of population path-length difference as a share of A's path, eigengap and `angle` null 95th percentile;
  - `accounting.json` (units, failures by cell, inconsistent units, duplicates, recorded core-hours).

  Verify with tests on synthetic records for the columns, failed and inconsistent units counted by cell, and duplicate handling.
- [x] 3.2 Implement both reproduction checks, written to `reproduction.json`:
  - the `paper_grade` `scores` arm against the committed `results/stage-span-paper-grade-2026-11-04/*/report/design_point_operating.csv`, exact on `n_available` and `n_rejected` for every column, cell and statistic, design phases only;
  - the `weak` `scores` arm against `results/span-split-2026-11-05/simulation/production_table.csv`, exact on count and rejections.

  Each lists every differing row with both counts. Verify with tests for both outcomes of each check, and that the committed files' column mapping (column name from the operating table's axes, matching `stage_span_leak_check.column_of`) resolves every one of the 65 × 3 rows.

## 4. Decision

- [x] 4.1 Implement `decide`, which writes `orthonormal_coordinates_decision.json` (rule version 1, design D7):
  - the gate: `reproduction_failed` on either failed reproduction or any inconsistent unit;
  - O1: two families through `stage_span_paper_grade_decision.judge_rates`, imported unchanged (m = 60, m = 6);
  - O2: m = 90 through the same rule;
  - O3: the one-sided exact paired loss p-value per target cell, with the count-plus-ceiling rule (allowed = 95th percentile of Binomial(45, 0.025); ceiling 0.05/45).

  Write verdict precedence, every criterion recorded whatever the verdict, an unavailable reading failing its criterion, the descriptive tables (watched cells with their exact two-sided paired p-values, population residuals per arm, eigengap and null-width changes) and the consequence statement for each verdict. Verify with tests on synthetic tables for:
  - each verdict, with the failing readings named;
  - precedence (`size_excess` over a failing O3);
  - the gate still recording O1–O3;
  - O3 at b = c = 0 (p = 1);
  - the allowed count at m = 45.
- [x] 4.2 Write the driver's module docstring (profiles, units, seeds, parents, outputs, decision) and add the driver to the `AGENTS.md` architecture beside the other validation scripts. Bump `src/motco/__init__.py` to 0.19.0. Verify with the pre-commit gate: `uv run ruff check src/ tests/ && uv run mypy src/motco/ && MOTCO_TEST_PERMS=99 uv run pytest tests/ -m "not slow" --tb=short`, plus `uv run ruff check scripts/orthonormal_coordinates_check.py`.

## 5. Runs and report

- [x] 5.1 Run a local smoke of both profiles at 2 replicates per cell, at 999 permutations for a few units and 19 for the rest, with BLAS pinned. Record per-unit timing by profile and arm, and confirm that `merge` and the reproduction check run, mismatching only because of the smoke's size. Verify that the core-hour projection and the proposed cap are written into `notes.md` and put to the user for approval before any array is submitted.
- [ ] 5.2 Run both profiles on `ing` through the sbatch array from a clean worktree at a pushed revision, then `merge` each. Verify 0 inconsistent units, failures only where the paper grade had them (the four-stage orientation cells), and both reproductions passed. Record the job ids, `sacct` output and package versions.
- [ ] 5.3 Run `decide`, then commit aggregates only under `results/orthonormal-coordinates-<date>/`: settings, rate and paired tables, accounting, reproduction, the decision, `PROVENANCE.txt` and `sacct.txt`. Write `docs/reports/orthonormal-coordinates-<date>.md` with the verdict, O1–O3, both reproductions, the watched cells, the population residuals per arm and what follows from the verdict. Verify that every number the report quotes can be read from the committed aggregates.
- [ ] 5.4 Update the score-coordinate section of `docs/measurement-visibility.md` (spec `latent-space-metric-compatibility`):
  - the opt-in setting;
  - basis invariance;
  - the verdict and its readings;
  - the watched cells;
  - that production and `motco analyze` stay on scores until a separate change adopts orthonormal coordinates.

  Update roadmap item 15 with the verdict and the follow-up it permits, and summarize the run in `AGENTS.md`. Verify with `uv run mkdocs build` and the pre-commit gate.

## Workflow follow-up

- Sync the delta specs and archive the change once the report is reviewed.
- Under `adopt_orthonormal`, a separate change flips the harness default with pins for every committed configuration and diagnostic. Another opens the `motco analyze` adoption: a real-noise size check at the case study's configuration, attribution under orthonormal coordinates, and the SEA-AD re-issue in `sead_sexdiff`.
- Under any other verdict, the shear stays documented as a measured limitation, and the setting remains an opt-in diagnostic.
