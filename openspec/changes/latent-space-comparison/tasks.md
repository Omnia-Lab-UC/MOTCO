# Tasks

## 1. Generator: alignment axis

- [x] 1.1 Add `span_alignment` and `span_support` to `SemiSyntheticTrajectoryParams`, with validation, and add both to `grid._SIGNATURE_OMIT_WHEN_NONE`.
  - φ = 0 normalizes to absent.
  - Reject φ outside [0, 90], a non-`stage_span` construction, a mode other than `orientation` or `none`, support without φ > 0, and a non-positive support.

  Verify with a test for each validation error, and that datasets, cell ids and signatures without the settings, and with φ = 0, are byte-identical to 0.15.0 (`tests/test_stage_span_construction.py` digests pass unchanged).
- [x] 1.2 Implement the out-of-span direction (design D3).
  - The support is the span-layer features whose standardized null stage means are constant (tolerance 1e-10).
  - Gaussian entries come from `np.random.default_rng([seed, <alignment stream tag>])`, with optional k-feature support, normalized.
  - Raise a named error when k exceeds the available features.

  Verify with tests that o is unit-norm, orthogonal to every centred null stage mean, has exactly k nonzero entries when sparse, and is identical across effects and nonzero φ at a seed.
- [x] 1.3 Tilt the stage-span orientation's plane to w_φ = cos φ·v + sin φ·o (design D1), and record φ, the support size, the span layers and the direction identifier in truth. Verify with tests at three and four stages:
  - the intended `angle` is 90°·e and the intended `delta` and `shape` are 0 at φ ∈ {30, 60, 90};
  - at φ = 90 with balanced cells, the varying population difference is orthogonal to the centred pooled stage means to 1e-10;
  - group A equals the `none` dataset and group B equals the null samples plus the plant, at the same seed.
- [x] 1.4 Document the alignment axis in the `semisynthetic.py` module docstring, in `docs/api/simulations.md` (stage-span construction section) and in the `semisynthetic.py` bullet of `AGENTS.md`. Verify that `uv run mkdocs build` adds no warnings.

## 2. Span capture for every linear space

- [x] 2.1 Add a linear-map record to `LatentIntegrationResult.artifacts` for `pls`, `pca` and `concat` (design D6): the input transform, an orthonormal basis of the kept directions, and the map to latent coordinates. The PLS record uses the same model attributes and `transform` as 0.15.0. Verify with a test that the `pca` record reproduces its latent matrix from the standardized joint rows, and that the PLS record reproduces `model.transform`.
- [x] 2.2 Make `measure_span_capture` and the `population_latent` checkpoint read the linear-map record, and enable them for `pca` and `concat` under `span_capture`.
  - `concat` reports shares of 1 and rank = `n_features`.
  - `snf` and `pls_within_group` record nothing.

  Verify with tests that:
  - PLS span-capture records are byte-identical to 0.15.0 on the `tests/test_span_capture.py` fixtures;
  - shares match the projection definition on a known basis for `pca`;
  - `concat` reads 1;
  - statistics and p-values are unchanged with and without span capture for every method.
- [x] 2.3 Update the realized-geometry section of `docs/api/simulations.md` and the `evaluation.py` bullet of `AGENTS.md` for span capture on every linear method. Verify that `uv run mkdocs build` adds no warnings.

## 3. Diagnostic `pls_residual_pca` integration

- [x] 3.1 Implement `pls_residual_pca` (design D4). It fits the production PLS, then stage-residualizes in the model's input coordinates, projects off the rotations' span and appends the top r residual component scores (default `n_stages − 1`). Add a helper that builds it from a prebuilt `pls` integration. Verify with tests that:
  - the leading columns equal the `pls` latent matrix;
  - the residual directions are orthogonal to the PLS span to 1e-10;
  - permuting the group labels within stages leaves the matrix unchanged;
  - the prebuilt path equals the from-scratch fit and fits no PLS model;
  - an infeasible r raises a validation error naming the range.
- [x] 3.2 Label it diagnostic, record the PLS rank, r and the retained dimensionality, add its latent checkpoint and linear-map record, and reject attribution and forced components. Verify with tests for the metadata, the checkpoint, span capture on this method, and both rejections.
- [x] 3.3 Document `pls_residual_pca` in `docs/api/simulations.md` (diagnostic integration methods) and in the `evaluation.py` bullet of `AGENTS.md`. Verify that `uv run mkdocs build` adds no warnings.

## 4. Comparison driver

- [x] 4.1 Write `scripts/latent_space_comparison.py` (design D8, D10).
  - Settings: the designs, the fourteen cells, the five arms, replicates, permutations, α and base seed 2035, bound to a settings hash.
  - Seeds: generator seed per (design, replicate), RRPP seed per (design, cell, replicate), and a separate calibration stream.
  - `evaluate_unit`: one dataset, one shared `pls` fit, `pca_matched` at its rank, `pca_fixed` at 2(K − 1), `pls_residual_pca` from the prebuilt fit, and `concat`. Each arm is evaluated with span capture and net displacement at the shared RRPP seed, and an arm's failure is recorded without aborting the unit.

  Verify with tests on a small design that:
  - the `pls_cv` record equals a direct production evaluation at the same seeds;
  - `pca_matched`'s rank equals `pls_cv`'s;
  - two cells at one replicate share the null dataset;
  - one PLS model is fitted per unit.
- [x] 4.2 Add `run` (round-robin shards, resume of complete units) and `merge`.
  - `merge` writes the rate table, the capture table and the accounting (units, failures by arm and cell, duplicates, recorded core-hours).
  - It refuses a resume or merge whose settings hash differs, naming the keys.
  - Add `scripts/latent_space_comparison_array.sbatch` (BLAS pinned to one thread, throttled array).

  Verify with tests for resume, a settings mismatch and the table columns on synthetic records.
- [x] 4.3 Add `calibrate` (design D9): the δ ladder {0.1, 0.15, 0.2, 0.25, 0.3, 0.4, 0.5} (amended to add 0.11–0.14 before any run unit, design D9) with 40 anchor replicates of the SEA-AD table per step, `pls_cv` fits only, adopting the first step whose median held-out AUROC lies in [0.60, 0.67], else the closest step marked `band_missed`. Make `run` refuse to start without a calibration file or an explicit δ. Verify with tests of the adoption rule on synthetic AUROC sequences (in band, `band_missed`) and the refusal.
- [x] 4.4 Document the driver in its module docstring and in `AGENTS.md` (architecture, beside the other validation scripts). Verify that `uv run ruff check scripts/latent_space_comparison.py` passes.

## 5. Advisory reading

- [x] 5.1 Implement `decide`, which writes `latent_space_comparison_reading.json`.
  - Criteria C1, C2 and C3 for `pca_matched`, `pca_fixed` and `pls_residual_pca` against `pls_cv`, with bound α + 2·√(α(1 − α)/R) and 2 pooled SE.
  - Per-arm verdicts `qualifies` | `trade_off` | `no_gain` | `size_excess`, naming every failing design, cell and statistic, and the overall verdict `propose_paper_grade` | `document_tradeoff`.
  - Descriptive tables: capture and target power against φ, `pls_cv`'s control rates, the weak design, and the sparse-minus-dense contrast with its pooled SE.

  Verify with tests on synthetic tables that exercise every per-arm verdict, both overall verdicts, `concat`'s ineligibility and the naming of failures.

## 6. Documentation corrections and release

- [x] 6.1 Correct `docs/measurement-visibility.md` (design D12).
  - The rhombus paragraph: state the condition at three or more stages (B's steps keep A's lengths and angles and are orthogonal to all of A's), and that a rotation within A's features is partly or wholly visible.
  - The `concat` sentence: a projection cannot see a difference `concat` lacks, but it can raise power by discarding noise directions, citing the probe's Phase 6 `shape` power.

  Verify that `uv run mkdocs build` adds no warnings and the cited numbers match the orientation probe report.
- [x] 6.2 Bump `src/motco/__init__.py` to 0.16.0. Verify with the pre-commit gate: `uv run ruff check src/ tests/ && uv run mypy src/motco/ && MOTCO_TEST_PERMS=99 uv run pytest tests/ -m "not slow" --tb=short`.

## 7. Pilot run and report

- [x] 7.1 Run a local smoke: every design at 2 replicates × 19 permutations, with δ given explicitly. Record per-unit and per-arm timing, and confirm that every arm's record carries span capture and `population_latent`. Verify that the core-hour projection is written into `notes.md`.
- [ ] 7.2 Run `calibrate` on `ing` from a clean worktree. Verify that the calibration output records every step and the adopted δ (or `band_missed`), and commit it under `results/latent-space-comparison-<date>/calibration/`.
- [ ] 7.3 Run the comparison on `ing` through the sbatch array, then `merge` and `decide`. Verify 0 failed units in the accounting, and record the job ids and `sacct` output.
- [ ] 7.4 Commit aggregates only under `results/latent-space-comparison-<date>/`: settings, calibration, rate and capture tables, accounting, the reading, `PROVENANCE.txt` and `sacct.txt`. Write `docs/reports/latent-space-comparison-<date>.md` with the verdict, each arm's visibility profile against φ, size, protected power, the weak design and the sparse-versus-dense contrast. Verify that every number the report quotes can be read from the committed aggregates.
- [ ] 7.5 Add the visibility-profile section to `docs/measurement-visibility.md` from the report. Add roadmap item 14 with the verdict and its follow-up, add follow-up (c) to item 13, and update the `AGENTS.md` architecture summary of the run. Verify with `uv run mkdocs build` and the pre-commit gate.

## Workflow follow-up

- Sync the delta specs and archive the change once the report is reviewed.
- If the reading is `propose_paper_grade`, propose the paper-grade change with its predeclared decision. It is also the place to decide on a four-stage design for an out-of-span bend (design Open Questions).
