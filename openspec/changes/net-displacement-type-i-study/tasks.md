# Tasks

## 1. Byte-identity baseline

- [x] 1.1 Confirm `tests/test_committed_config_identity.py` and the design-grid, rank-ladder, report and targets tests pass on unchanged code; this is the guard every later task re-runs.
- [x] 1.2 Capture a pre-change report fixture: run `scripts/motco_study.py report` on an existing tiny design-grid smoke result (or build one with `n_replicates: 2`, 19 permutations) and store the output tree's file list and CSV contents under the scratchpad for byte comparison in task 3.4.

## 2. Explicit design points (design D1)

- [x] 2.1 Add `DesignGrid.points`, derive `axes` from the points, and make `design_points()` return the list in listed order when set. Verify tests: the explicit-points scenario loads with the points in order, `axes` lists every axis with values in first-appearance order, and a dump/reload round-trips `points` (with `null` nested values) unchanged.
- [x] 2.2 Add load-time validation: both forms declared, a point with a differing axis set, duplicate points after value normalization (including `"pc1"` vs `null` for the estimator), and a missing baseline point are each rejected with an error naming the offender. Verify one test per rejection.
- [x] 2.3 Verify enumeration over explicit points: a test that a config listing the baseline plus four points emits four `power_design` grids with the right coordinates, and a test that an explicit list equal to the full cross of two axes yields identical cell ids, signatures and replicate seeds to the crossed declaration. Re-run task 1.1.

## 3. Estimator and method axes, reader hardening (design D2)

- [x] 3.1 Add tests that `evaluation.orientation_estimator: [null, "net_displacement"]` as a design axis gives net-displacement column records carrying the estimator in runtime metadata with `angle` p-values from that estimator, and baseline records byte-identical to the same config without the axis (tiny end-to-end run, 2 replicates, 19 permutations).
- [x] 3.2 Add a test that `evaluation.integration_method` as a design axis evaluates `concat` and `pca` columns (the latter with `evaluation.integration_params.n_components`), and that a `pls_forced` explicit point passes `assert_production_component_selection` when the recorded rank equals the point's rank.
- [x] 3.3 Run `scripts/motco_study.py merge` and `report` on a tiny explicit-points study with all five columns; fix any table or figure code that fails on records without `selected_lv` or `component_selection` (expected sites: `render_design_point_power`, `render_rank_ladder`). Verify the report completes and `design_point_operating.csv` lists every column with unavailable component-selection fields for `concat` and `pca`.
- [x] 3.4 Re-run the task 1.2 fixture through the modified report and verify byte-identical output; re-run task 1.1.

## 4. Anchor decision rule (design D3)

- [x] 4.1 Add `AnchorDecisionRule` (`alpha`, `se_multiplier`, `statistics`, `name`) under `acceptance.anchor_decision` with parsing and validation (unknown fields, statistics outside `delta`/`angle`/`shape`, non-positive multiplier rejected). Verify tests for load, dump/reload and each rejection.
- [x] 4.2 Implement `evaluate_anchor_decision` as a pure function over records: per design point, anchor rate/SE/available R per statistic, bound α + k·√(α(1 − α)/R), `holds`/`exceeds`/`unavailable`, per-point verdict, other negative controls listed descriptively. Verify unit tests on synthetic records for the three spec scenarios (0.068 holds at R = 500, 0.080 exceeds, rule absent writes nothing) and for a point with no available replicates.
- [x] 4.3 Write `anchor_decision.json` and `.csv` from `scripts/motco_study.py report` when the rule is declared, and verify on the task 3.3 smoke study that both files appear and that a config without the rule produces the task 1.2 fixture byte-for-byte.
- [x] 4.4 Run the pre-commit gate (`ruff`, `mypy`, fast tests) and confirm it is green.

## 5. Committed profiles (design D4)

- [x] 5.1 Write `examples/trajectory_power_study/net_displacement_type_i_phase6.json` (generator and evaluation from `phase6_small_n_study.json`, base seed 800, family `phase6-small-n`, translation-only modes, effects `[0.0, 1.0]`, five explicit points at rank 4, Type I and translation specificity targets, anchor rule, gate and attribution off, metadata with intent, cost notes and the submission block). Verify it loads, enumerates twelve cells at 500 × 999, and every pool-limited cell fits its headroom.
- [x] 5.2 Write `net_displacement_type_i_phase5.json` (n = 1200, four stages, `group_ratio` 0.5, all blocks, joint magnitude, base seed 600, family `phase5-primary`, rank 6, otherwise as 5.1). Verify the same twelve-cell structure.
- [x] 5.3 Verify dataset sharing: a test that the Phase 6 profile's baseline anchor and translation cells draw the same replicate seeds as `phase6_small_n_study.json`'s, and that the Phase 5 profile's anchor datasets for replicates 0–2 hash-equal those generated from `phase5_power_study.json`'s anchor at the same indices.
- [x] 5.4 Add both profiles to `tests/data/committed_config_identity.json` without touching any existing entry, and verify `tests/test_committed_config_identity.py` passes with the two new parametrizations and all old ones unchanged.
- [x] 5.5 Add `tests/test_net_displacement_type_i_profiles.py` covering the two profile-loaded scenarios (cell structure, no `clamp`, rule declared) and run the pre-commit gate.

## 6. Follow-up selection (design D5)

- [x] 6.1 Write `scripts/net_displacement_follow_up.py`: reads the two `anchor_decision.json` files, identifies the lead column by coordinates, writes `follow_up_decision.json` with verdict, per-design localization, rule version, input paths and sha256. Verify unit tests on synthetic decisions for both spec scenarios and for a missing lead column (clear error).
- [x] 6.2 Verify a smoke invocation on the task 4.3 outputs produces the file with every declared field.

## 7. Run on `ing` (design D6)

- [ ] 7.1 Time two anchor units and one translation unit per column per design on `ing`, write the projected core-hours per design to `PROVENANCE.txt`, and apply the go/no-go in the design's Risks (drop `pca` and `concat` from Phase 5 only if the projection exceeds 150 core-hours, and record that).
- [ ] 7.2 Submit both profiles with `scripts/motco_study_array.sbatch` (50 shards each, BLAS pinned to one thread, `STUDY_N_JOBS` unset, memory-pressured nodes avoided) and verify completion: 6,000 units per design with zero failures, or failures resubmitted by shard id and explained.
- [ ] 7.3 Merge and report each design into `results/net-displacement-type-i-<date>/{phase6,phase5}/`, run the follow-up script, and verify `anchor_decision.json`, `design_point_operating.csv` and `follow_up_decision.json` exist with the expected row and column counts; complete `PROVENANCE.txt` (command lines, code revision, software versions, host, unit timings, realized core-hours).

## 8. Findings and docs

- [ ] 8.1 Write `docs/reports/net-displacement-type-i-<date>.md`: provenance and reproduction commands, the per-design anchor table (rate ± SE vs bound per column and statistic), the translation control table, eigengap and `angle` null-width covariates, whether the baseline column reproduces the parent studies' recorded p-values exactly, the anchor decision and the follow-up verdict naming the next change. Verify every number cites a committed file, and that no power figure or orientation power claim appears.
- [ ] 8.2 Update `docs/roadmap.md` (item 3 closed with the branch taken; the selected follow-up as the new item), `src/motco/simulations/study/README.md` (explicit points, anchor rule, follow-up script, the two profiles), and the `CLAUDE.md` study paragraph. Verify `openspec validate net-displacement-type-i-study --strict` passes and the pre-commit gate is green.
