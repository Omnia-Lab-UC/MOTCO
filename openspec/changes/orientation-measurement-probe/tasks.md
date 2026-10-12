# Tasks

## 1. Byte-identity baseline

- [x] 1.1 Capture a pre-change fixture of evaluation results for a fixed small dataset under `concat`, `snf` and `pls` (observed statistics, p-values at a fixed RRPP seed, realized geometry, integration metadata; runtime fields excluded), committed under `tests/data/`. Add a test asserting the current code reproduces it, and verify it passes on unchanged code.
- [x] 1.2 Confirm `tests/test_committed_config_identity.py` passes on unchanged code; this is the guard every later task re-runs.

## 2. Net-displacement orientation estimator (stats)

- [x] 2.1 Add `_estimate_net_displacement` and the `orientation_estimator` keyword on `estimate_difference` (design D1). Verify unit tests for each scenario in the `trajectory-orientation-invariance` delta: first-to-last direction, two-stage agreement with PC1, equilateral configuration gives 0°, reversed order gives 180°, only `angle` differs, zero length is reported unavailable.
- [x] 2.2 Thread the keyword through `RRPP` and `_RRPPWorker`, serial and parallel. Verify a test that the observed and permuted angles use the selected estimator, and that the default RRPP output equals the task 1.1 fixture.
- [x] 2.3 Thread the keyword through `diagnostics.geometry_from_means` and `calculate_realized_geometry`. Verify a test that checkpoint angles change with the estimator, and that defaults still match the fixture.

## 3. Harness parameters and diagnostic integration methods

- [x] 3.1 Add `SimulationEvaluationParams.orientation_estimator` (`None` = PC1; `"pc1"` normalized to `None`; unknown values rejected) and add it to `grid._SIGNATURE_OMIT_WHEN_NONE` (design D2). Verify tests for the three harness-delta estimator scenarios, then re-run tasks 1.1 and 1.2.
- [x] 3.2 Implement the `pca` integration method on the `concat` matrix, with a required `n_components` (design D3). Verify tests that it is label-free (scores invariant to label permutation), that a missing or infeasible count is rejected, and that the metadata labels it diagnostic.
- [x] 3.3 Implement `pls_within_group` (per-group `scale=False` PLS-DA on pooled standardized rows, stacked `x_rotations_`, QR orthonormalization, pooled-centered projection; design D4). Verify tests: on a constructed dataset with orthogonal group directions, both groups' stage displacements lie in the basis span; a stage with fewer than two samples in one group raises an error naming them; the metadata records the per-group count and basis rank.
- [x] 3.4 Reject attribution and `forced_components` for both new methods. Verify validation-error tests.
- [x] 3.5 Pass the latent matrix to realized geometry for `pca` and `pls_within_group`, under scopes `pca_latent` and `pls_within_group_latent` (design D5). Verify tests that the scopes are present for the new methods, that `pls_latent` is unchanged, and that `snf`/`concat` records still equal the fixture.
- [x] 3.6 Run the pre-commit gate (`ruff`, `mypy`, fast tests) and confirm it is green.

## 4. Motivating-evidence script

- [x] 4.1 Write `scripts/orientation_checkpoint_evidence.py` (design D8). Verify that on `results/phase6-small-n-2026-10-09/merged.jsonl` it reproduces the exploration tables: population angles 41.4/60.0/75.7/90.0°; PLS latent medians 21.6/24.8/28.5/32.1° (rank ≤ 3) and 76.0° (rank > 3, e = 1.0); rejection rates 0.082 and 0.506; anchor signed q95 of about 137°; folded oracle separability 0.996 (observed) and 0.496 (latent) at e = 1.0.
- [x] 4.2 Add a smoke test running the script on a tiny synthetic `merged.jsonl`, and verify it passes.

## 5. Probe driver

- [x] 5.1 Write `scripts/orientation_measurement_probe.py`: the two designs, five cells and five spaces × two estimators; matched generator seeds per replicate index and one RRPP seed per unit (design D6). Verify a smoke run (`--reps 2 --perms 19`, small n) in which every arm's record carries identical generator and RRPP seeds for a unit.
- [x] 5.2 Add per-(unit, arm) JSONL output, resume by (design, cell, replicate), a settings-hash header that refuses mismatched resumes, and `--shard i/N` plus `--merge`. Verify tests that an interrupted-then-resumed smoke run equals an uninterrupted one, and that a resume with changed settings exits with an error naming the difference.
- [x] 5.3 Implement the rate table and the advancement rule as pure functions (design D7; criteria per the probe spec). Verify unit tests on synthetic tables for the three spec scenarios (anchor inflation, specificity failure, empty advancing set), that `concat` never advances, and the 2-pooled-SE thresholds.
- [x] 5.4 Write the Markdown summary (rate tables, a checkpoint table per arm, the decision, the Phase 5 joint-magnitude deviation note). Verify that the smoke run produces it, with every section present.
- [x] 5.5 Run the pre-commit gate and confirm it is green.

## 6. Run on `ing`

- [ ] 6.1 Time one unit per design on `ing` (all arms) and record the projected core-hours. Verify the projection is written down before launch.
- [ ] 6.2 Launch the sharded SLURM array for both designs (BLAS pinned to one thread; avoid memory-pressured nodes). Verify completion: 2 × 5 × 100 units, each with every arm present and no failures, or failures recorded and explained.
- [ ] 6.3 Merge and render into `results/orientation-measurement-probe-<date>/` with `PROVENANCE.txt` (command, code revision, software versions, host, unit timings). Verify the merged output's unit count and that the decision file exists.

## 7. Findings and docs

- [ ] 7.1 Write `docs/reports/orientation-measurement-probe-<date>.md`: §1 motivating evidence (task 4.1 output), the per-arm rates, the checkpoint table per arm, the advancement decision, and the recommendation for the follow-up study's arms. Verify every number cites a committed file or the evidence script.
- [ ] 7.2 Update `docs/roadmap.md` (probe result, follow-up study entry, superseded ρ sweep) and the `CLAUDE.md` architecture notes (new diagnostic integration methods, orientation estimator). Verify `openspec validate orientation-measurement-probe --strict` passes and the pre-commit gate is green.
