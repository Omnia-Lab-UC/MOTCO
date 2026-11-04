# Tasks

## 1. Generator: per-cell mean shift

- [ ] 1.1 Add optional per-cell mean-shift arrays (features × cells, one per omic) to `omic_population_means` and `generate_omics` in `src/motco/simulations/generator.py`. Shift methylation in M-values before `rev_logit`; an absent shift takes the existing path. Verify with a test that a zero or absent shift is byte-identical to the current output, and that a shift moves the samples by exactly the shift at the same seed.

## 2. Generator: stage-span construction

- [ ] 2.1 Add `construction` and `span_layers` to `SemiSyntheticTrajectoryParams`, with validation (`surgery` normalizes to `None`; unknown, duplicated or empty layers are rejected; `span_layers` without `stage_span` is rejected; orientation e > 2 is rejected), and add both to `grid._SIGNATURE_OMIT_WHEN_NONE`. Verify with tests for each validation error and that a default-parameter dataset and cell id are byte-identical to before.
- [ ] 2.2 Implement the null population pooled scale per feature over `span_layers`: within-cell variance from the reference covariance diagonals, plus the between-cell variance of the null means at the realized cell proportions. Verify with a test against a large-n sample's pooled SD at a fixed seed, within sampling tolerance.
- [ ] 2.3 Implement the four stage-span edits on B's standardized null stage means (design D4) and the degeneracy errors (design D5). Map the plant back to analysed units and pass it as the per-cell shift; leave blocks outside the span untouched. Verify with tests at three and four stages:
  - in null-standardized coordinates, the intended `delta`, `angle` and `shape` match each mode's specification (magnitude ratio 1 + e, angle 90°·e, shape > 0 with `delta` and `angle` at 0, translation all 0);
  - fewer than three stages, or a rank-1 configuration, raises for `orientation` and `shape`;
  - proteomics is unchanged when the span excludes it.
- [ ] 2.4 Make stage-span datasets carry the plant in `population_trajectories` and record truth (construction, span layers, effect parameters, intended geometry, plant norm per stage). Verify with a test that group A's samples equal the `none` dataset's, that group B's equal the null samples plus the plant (methylation in M-values) to floating-point precision, and that B's population means equal the null means plus the plant.
- [ ] 2.5 Document the construction in the `semisynthetic.py` module docstring and in `docs/api/simulations.md`, and update the `semisynthetic.py` bullet in `AGENTS.md`. Verify that `uv run mkdocs build` succeeds.

## 3. Harness and realized-geometry diagnostics

- [ ] 3.1 Add `span_capture` to `SimulationEvaluationParams` (`false` normalizes to `None`; hash-omitted) and pass the fitted PLS model and preprocessor to the realized-geometry calculation when it is on. Verify with a test that the PLS model is fitted once, and that the latent matrix, rank, statistics and p-values are identical with and without it.
- [ ] 3.2 Implement the `population_latent` checkpoint (population means through the fitted preprocessor, the measured blocks and `model.transform`, measured with the evaluation's estimator). Verify with tests that it is present only for PLS with span capture and population trajectories, and that for a stage-span plant at large n its target statistic approaches the intended one.
- [ ] 3.3 Implement the `span_capture` block (`group_difference`, `pooled_stage`, `rank`, `n_features`, `random_share`) in the model's input coordinates with a QR basis of the rotations (design D8). Verify with tests that:
  - the shares lie in [0, 1];
  - `group_difference` is null at the anchor and under both translations;
  - a stage-span plant's share exceeds 0.9 at large n;
  - on a fixed matrix the definition reproduces the plasmode `latent_capture.py` share.
- [ ] 3.4 Add the harness guard that refuses a stage-span dataset whose truth's span layers differ from the measured blocks, raised before any fit. Verify with a test that the error names both block sets and that no PLS fit is attempted.
- [ ] 3.5 Update `AGENTS.md` (the `evaluation.py` bullet and the realized-geometry checkpoints). Verify that the existing realized-geometry and identity test suites pass unchanged.

## 4. Study machinery

- [ ] 4.1 Accept `generator.construction`, `generator.span_layers` and `evaluation.span_capture` in the study config. Make the loader reject a span-layer mismatch at any column (naming the column and both block sets) and a `generator.construction` design axis. Verify with config tests for each error, and that `tests/test_committed_config_identity.py`, `tests/test_post_floor_config_identity.py` and the default-pin tests pass unchanged.
- [ ] 4.2 Make the zero-effect anchor builder clear `construction` and `span_layers`, and exempt stage-span cells from the surgery headroom check and from realized-surgery summaries. Verify with tests that a stage-span profile and its surgery sibling enumerate anchors with identical cell ids, signatures and matched seeds, and that large stage-span effects enumerate without a headroom failure.
- [ ] 4.3 Write `report/span_capture.csv` when records carry span capture (columns as in the trajectory-power-study delta) and nothing otherwise. Verify with a report test on synthetic records, and that every committed `report/` regenerates byte-identically.
- [ ] 4.4 Document the new config keys, the sibling-profile rule and `span_capture.csv` in `src/motco/simulations/study/README.md` and the `study/` bullet of `AGENTS.md`.

## 5. Pilot profiles

- [ ] 5.1 Write the four profiles under `examples/trajectory_power_study/`: `stage_span_pilot_{surgery,span}_{magnitude,geometry}.json`.
  - Derive the generator and evaluation blocks from `design_guidance_core_*_pilot.json`, adding `evaluation.span_capture: true`.
  - Use the `design-guidance` family at base seed 1028, 100 × 199, with explicit points over `generator.n_samples` and `generator.group_stage_sizes` (n = 200 baseline, n = 80, and the default `n_samples` with `[[11, 10, 28], [9, 10, 12]]`), three stages, and `layers` methylation + expression.
  - Under `stage_span`, set `span_layers` to the same blocks.
  - Use magnitude effects 0.02, 0.05 and 0.1, and geometry modes `orientation`, `shape` and `translation` at 0.25, 0.5 and 1.0.

  Verify that each loads and enumerates to three columns with the specified cells.
- [ ] 5.2 Register the profiles in `tests/post_floor_fixture_capture.py` (`POST_FLOOR_CONFIGS`) and regenerate `tests/data/post_floor_config_identity.json`. Verify that `tests/test_post_floor_config_identity.py` passes and that the diff adds only the four entries.
- [ ] 5.3 Add `tests/test_stage_span_pilot_profiles.py`. It checks the structure, the sibling relationship (identical except construction and span layers), the shared anchor identities, and that each surgery sibling's n = 200 and n = 80 cells draw the same generator seeds and parameters as the matching design-guidance pilot cells. Verify that it passes.

## 6. Reproduction pairing

- [ ] 6.1 Extend `scripts/anchor_reproduction.py` with a cell-level pairing on (projected design column, mode, effect, replicate index), so the surgery siblings can be paired with the design-guidance pilot records despite different design-axis sets. The default and `--by-column` paths stay unchanged. Verify with `tests/test_anchor_reproduction.py` cases for matched, mismatched and missing pairs.

## 7. Latent retrofit check

- [ ] 7.1 Promote the scratchpad prototype to `scripts/latent_retrofit_check.py`.
  - Parameters: `--config` (a pilot profile and column), `--replicates`, `--permutations`, `--out-dir`.
  - Behaviour: all four edits, both preimages, frozen and refit re-measurement, footprints (effective features, stage-active share, methylation share, constant-over-stages share), and the cosine with the population stage-span plant.
  - Outputs: `latent_retrofit_check.csv` and `summary.json`; no decision.

  Verify with a smoke test on a small dataset that the frozen-basis recovery error is at most 1e-8 and that the outputs have the documented columns.

## 8. Advisory reading

- [ ] 8.1 Implement `scripts/stage_span_pilot_reading.py`. It reads the four reports and writes `stage_span_pilot_reading.json`:
  - criteria R1 (stage-span translation at or below α + 2·√(α(1 − α)/R) for each statistic and column), R2 (median `group_difference` ≥ 0.90 in every stage-span magnitude, orientation and shape cell) and R3 (off-target rates at or below the bound in those cells);
  - verdict `propose_paper_grade` | `revise_construction`, naming every failing criterion, cell and statistic;
  - descriptive tables: target power and capture for both constructions, the stage-span minus surgery magnitude `delta` power difference with its pooled SE, and anchor size as read.

  Verify with tests on synthetic reports that exercise each verdict path and each criterion's failure naming.

## 9. Visibility documentation and release

- [ ] 9.1 Write `docs/measurement-visibility.md` and add it to the MkDocs nav. Cover:
  - the span principle and the rhombus argument;
  - which differences are visible;
  - the incidental routes and their mixed signature;
  - the orientation probe's `pca` and `pls_within_group` evidence with its costs;
  - the plasmode capture share;
  - what a non-rejection means;
  - a placeholder section for the pilot's capture shares, filled in by 10.6.

  Link it from `docs/api/snf.md` (where the measurement-space assumptions are documented), `docs/api/analysis.md` and `docs/index.md`. Verify that `uv run mkdocs build` succeeds with no new warnings.
- [ ] 9.2 Bump `src/motco/__init__.py` to 0.15.0. Verify with the pre-commit gate: `uv run ruff check src/ tests/ && uv run mypy src/motco/ && MOTCO_TEST_PERMS=99 uv run pytest tests/ -m "not slow" --tb=short`.

## 10. Pilot run and report

- [ ] 10.1 Run a local smoke: each profile at 2 replicates × 19 permutations. Record per-unit timing and confirm that records carry `span_capture` and `population_latent`, and that the stage-span translation and anchors run without error. Verify that the core-hour projection is written into the run notes.
- [ ] 10.2 Run the four profiles on `ing` from a clean worktree through `scripts/motco_study_array.sbatch` and `scripts/run_study_shard.py`, with a throttled array. Merge and report each. Verify 0 failures in the merged results and record the job ids and `sacct` output.
- [ ] 10.3 Run the reproduction checks:
  - `scripts/anchor_reproduction.py --by-column` across the four siblings;
  - the cell-level pairing of the surgery siblings against the 2026-10-28 design-guidance pilot records at n = 200 and n = 80.

  Verify that the summaries report pair and mismatch counts; any mismatch is diagnosed before 10.5.
- [ ] 10.4 Run `scripts/latent_retrofit_check.py` at the n = 200 column, replicates 0–9, 199 permutations. Verify that `summary.json` is written.
- [ ] 10.5 Run `scripts/stage_span_pilot_reading.py`. Verify that `stage_span_pilot_reading.json` is written with a verdict.
- [ ] 10.6 Commit aggregates only under `results/stage-span-pilot-<date>/`: reports, `span_capture.csv`, the reading, retrofit tables, reproduction summaries, `PROVENANCE.txt` and `sacct.txt`. Write `docs/reports/stage-span-pilot-<date>.md` with the verdict, every mode's capture under both constructions, reproduction, retrofit recovery and the proposed paper-grade effect grid, and fill the pilot section of `docs/measurement-visibility.md`. Verify that `uv run mkdocs build` succeeds and the report's numbers match the committed aggregates.
- [ ] 10.7 Update `docs/roadmap.md` (new Next-changes item 13 with the verdict, the "Not yet established" orientation entry, and the follow-ups: paper-grade profiles and the latent plasmode) and the `AGENTS.md` architecture summary of the run. Verify with the pre-commit gate.

## Workflow follow-up

- Sync the delta specs and archive the change once the report is reviewed.
- If the reading is `propose_paper_grade`, propose the paper-grade change with its predeclared decision, then the latent plasmode change.
