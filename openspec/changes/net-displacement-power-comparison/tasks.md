# Tasks

## 1. Byte-identity baseline

- [x] 1.1 Confirm `tests/test_committed_config_identity.py` and the design-grid, rank-ladder, anchor-decision, report and targets tests pass on unchanged code; this is the guard every later task re-runs.
- [x] 1.2 Capture a pre-change report fixture: run `scripts/motco_study.py merge` + `report` on a tiny explicit-points study (two `pls` columns differing in `evaluation.orientation_estimator`, four modes, `n_replicates: 2`, 19 permutations) and store the output tree's file list and CSV/JSON contents under the scratchpad for byte comparison in task 2.5.

## 2. Column decision rule (design D1, D2)

- [x] 2.1 Add `ColumnDecisionRule` (`candidate`, `target`, `criterion`, `gain_se_multiplier`, `loss_se_multiplier`, `type_i_bound`, `controls`, `protected`, `name`) under `acceptance.column_decision` with parsing, dump/reload, and validation: unknown fields, unknown criterion kind, empty controls, target also in controls/protected, modes absent from `trajectory_modes`, candidate axes not in `design_grid`, and candidate matching zero, several, or the baseline point (after axis-value normalization, so `"pc1"` matches `null`). Verify one test per rejection plus load/round-trip.
- [x] 2.2 Implement `evaluate_column_decision` in `targets.py` as a pure function over records: group by stamped design point, locate the candidate, measure anchor / target / controls / protected per column at each mode's top effect, evaluate the four criteria with `type_i_inflation_bound` for the anchor and control bounds, and produce verdicts `candidate_qualifies` / `candidate_fails` / `unavailable` / `no_candidate` / `no_design_grid`. Verify unit tests on synthetic records for the spec scenarios: gain qualifies (0.670 vs 0.300, R = 500), control 0.080 fails against 0.0695, no-loss qualifies (0.840 vs 0.860), and an unavailable control yields `unavailable`.
- [x] 2.3 Add `write_column_decision` (`column_decision.json` + `.csv`, one row per column with flattened criteria) and call it from `scripts/motco_study.py report` after the anchor decision when the rule is declared. Verify on the task 1.2 smoke study with the rule added that both files appear and the CSV lists the baseline and candidate rows.
- [x] 2.4 Verify a config without the rule produces the task 1.2 fixture byte-for-byte, and re-run task 1.1.
- [x] 2.5 Run the pre-commit gate (`ruff`, `mypy`, fast tests) and confirm it is green.

## 3. Committed profiles (design D3)

- [x] 3.1 Write `examples/trajectory_power_study/net_displacement_power_phase6.json`: generator and evaluation from `phase6_small_n_study.json`, base seed 800, family `phase6-small-n`, four modes, effects `[0, 0.25, 0.5, 0.75, 1.0]`, two explicit points over `evaluation.integration_method` × `evaluation.orientation_estimator`, Type I and five specificity targets with `bound: "upper"`, no power target, gate and attribution off, column rule with `criterion: gain`, metadata (intent, derives_from, parent_results, notes incl. the new magnitude 0.75 cell and "do not pass --n-jobs", submission block). Verify it loads, enumerates 36 cells at 500 × 999, and every pool-limited cell fits its headroom.
- [x] 3.2 Write `net_displacement_power_phase5.json` (n = 1200, four stages, `group_ratio` 0.5, all blocks, `magnitude_kind: joint`, base seed 600, family `phase5-primary`, `criterion: no_loss`, otherwise as 3.1, notes naming the magnitude 0.5/0.75 cells as new under `joint`). Verify the same 36-cell structure.
- [x] 3.3 Verify dataset sharing by test: the Phase 6 profile's baseline cells draw the same replicate seeds as `phase6_small_n_study.json`'s cells of the same mode and effect, and the Phase 5 profile's orientation e = 1.0 datasets for replicates 0–2 hash-equal those generated from `phase5_power_study.json` at the same cell.
- [x] 3.4 Add both profiles to `tests/data/committed_config_identity.json` without touching any existing entry; verify `tests/test_committed_config_identity.py` passes with the two new parametrizations and all old ones unchanged.
- [x] 3.5 Add `tests/test_net_displacement_power_profiles.py` covering the profile-loaded scenarios (cell structure, no `clamp`, every target `bound: upper`, rule declared with the three controls and the per-design criterion) and run the pre-commit gate.

## 4. Comparison selection script (design D4)

- [x] 4.1 Write `scripts/net_displacement_power_follow_up.py` (`--phase6`, `--phase5`, `--out`): verdict `adopt_net_displacement` iff both candidates qualify, else `keep_pc1` with per-design verdict and failing criteria/pairs; rule version, `advisory: true`, input paths with sha256. Verify unit tests on synthetic decision dicts for both spec scenarios and for a malformed input (clear error).
- [x] 4.2 Add the record-level consistency checks script `scripts/net_displacement_power_checks.py`: (a) `delta`/`shape` p-values and observed statistics identical between the two columns for every (cell coordinates, replicate); (b) baseline records equal the parent runs' records (seed, p-values) wherever a parent cell of the same mode and effect exists. Verify on the task 1.2 smoke output (a) and on a synthetic parent file (b), and that both report counts, not just booleans.

## 5. Cluster run (design D6)

- [x] 5.1 Timing job on `ing` for both profiles (per column: anchor replicates 0–1 and orientation e = 1.0 replicate 0; BLAS pinned; problem nodes excluded); write the projected budget to each `PROVENANCE.txt` and apply the Phase 5 go/no-go (drop e = 0.75 and re-pin its identity entry if the projection with margin exceeds 450 core-hours). Verify the timing units of the baseline column reproduce the parents' p-values exactly.
- [x] 5.2 Submit both profiles via `scripts/motco_study_array.sbatch` (50–100 shards, `STUDY_N_JOBS` unset, `--error-policy record`); on completion merge, verify 18,000/18,000 units per design, 0 failed, 0 censored, every signature equal to the enumerated cell's, and record wall/core-hours in `PROVENANCE.txt`.
- [x] 5.3 Run `scripts/motco_study.py report` for each design, then `scripts/net_displacement_power_follow_up.py` and `scripts/net_displacement_power_checks.py`; verify `column_decision.json/.csv`, `comparison_decision.json` and the checks' outputs exist and the checks report zero mismatches (or document every mismatch).

## 6. Findings and documentation (design D5)

- [x] 6.1 Write `docs/reports/net-displacement-power-<date>.md`: design, unit accounting, `angle` rate by mode × effect × estimator per design with a power figure, the control table against the bound, covariates (eigengap, `angle` null width, selected rank), the two consistency checks, the column decisions and the comparison verdict, the interpretability consequence (revised `angle` statement with measured power if adopted; unchanged otherwise), provenance and reproduction commands. Verify every number quoted is read from the committed `report/` and decision files.
- [x] 6.2 Update `docs/roadmap.md` (close *Next changes* item 2 with the branch taken; add the done-entry under Phase 6; if adopted, add the default-flip change to *Next changes*), `src/motco/simulations/study/README.md` (column rule, follow-up and checks scripts, the two profiles), `examples/trajectory_power_study/README.md`, and the CLAUDE.md `study/` paragraph. Verify links resolve and the pre-commit gate is green.
