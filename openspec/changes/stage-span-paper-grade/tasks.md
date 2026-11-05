# Tasks

## 1. Branch and worktree

- [x] 1.1 Create `feat/stage-span-paper-grade` from `main` explicitly in a worktree under `.claude/worktrees/stage-span-paper-grade` (the checkout is shared with other sessions), and verify that `git log -1` there is `main`'s head and that `git status` is clean

## 2. Planted-scale stage-span plant

- [x] 2.1 Add `span_scale` to `SemiSyntheticTrajectoryParams` (`null` → `None`; values `null`/`planted`; rejected without `stage_span`, with an unknown value, or with `span_alignment`), and verify the validation errors name the setting or value in new tests
- [x] 2.2 Add a planted pooled-SD helper beside `null_pooled_scale` and the fixed point in `_stage_span_plant` (design D1: start at the null SD, tolerance 1e-12, limit 200, a validation error naming mode, effect and residual), and verify with tests that for every mode at e ∈ {0.1, 0.5, 1.0} (orientation also 2.0) at n = 200, n = 80, the SEA-AD table, three blocks and four stages, the planted-coordinate intended geometry matches the mode's definition to 1e-8 relative
- [x] 2.3 Record `span_scale`, `units: "planted_standardized"`, `scale_iterations` and `scale_residual` in truth under `planted` (design D3), and verify with tests that truth without the setting is unchanged
- [x] 2.4 Add `span_scale` to `grid._SIGNATURE_OMIT_WHEN_NONE`, and verify that the committed config identity tests (`tests/test_committed_config_identity.py`, both pin fixtures) and the stage-span digest tests (`tests/test_stage_span_construction.py`, `tests/test_span_alignment.py`) pass unchanged
- [x] 2.5 Test that a planted dataset is the null dataset plus the plant (group A identical, B differs by the plant in analysed units), that zero effect yields the null dataset, and that planted magnitude generates the same samples as null-scale magnitude at the same seed
- [x] 2.6 Document `span_scale` in the `semisynthetic.py` module docstring and the `AGENTS.md` `semisynthetic.py` entry, bump the version to 0.17.0 in `src/motco/__init__.py`, and verify `uv run ruff check src/ tests/` and `uv run mypy src/motco/` pass

## 3. Study configuration keys

- [x] 3.1 Accept `generator.span_scale` in the study config (dump-omitted while unset), reject it as an OFAT or design axis with the sibling-profile message, and clear it with `construction`/`span_layers` in zero-effect and Type I null cells (`enumerate._without_construction`); verify with tests in `tests/test_study_stage_span.py` that a planted config loads and enumerates, that `null` equals omission, that the axis is refused, and that planted anchors equal surgery anchors
- [x] 3.2 Document `generator.span_scale` in `src/motco/simulations/study/README.md` beside the stage-span keys, and verify that every committed config enumerates unchanged (the identity tests from 2.4)

## 4. Leak decomposition check

- [x] 4.1 Write `scripts/stage_span_leak_check.py` with `run` (`--config`, `--shard`, `--n-shards`, `--span-scale`, `--modes`, `--effects`, `--replicates`, `--permutations`, `--out-dir`) and `merge`, per design D7, writing `leak_units.csv`, `leak_cells.csv` and `summary.json`; verify with tests on a small synthetic profile that all six checkpoints are reported, that the null-scale checkpoint difference is 0 for the bend, that the planted-scale checkpoint difference is 0 under `planted`, and that rates use p < α
- [x] 4.2 Add a `gate` subcommand that reads the planted and null arms' merged tables and writes `preflight_gate.json` (exact scale route, no failed units, the bound rule on off-target rates, residual shrinkage at e = 0.5 per column), and verify with fixture tables that each condition can fail on its own and that the verdict names the failing cells
- [x] 4.3 Add `scripts/stage_span_leak_check_array.sbatch` (one CPU, BLAS pinned to one thread, `--export` of `PYTHONPATH`, `logs/` created beforehand), and verify it with `bash -n` and a dry render of its command line
- [x] 4.4 Reproduce the scratch result locally on one replicate of the pilot's n = 200 shape e = 1.0 cell, and verify the script's checkpoint values equal the scratch values (`pop_planted` −0.6524, `pop_lat` −0.4862 at replicate 0)

## 5. Paper-grade decision

- [x] 5.1 Implement the bound rule (design D5): per-row bound α + 2·SE, the exact exceedance probability p₂, the 95th-percentile allowed count, and the Bonferroni ceiling at 0.05/m; verify with tests that the computed values match the design (R = 500: p₂ ≈ 0.0303, allowed 4 at m = 60 and 6 at m = 90, ceiling 0.084; R = 100, m = 60: 4 and 0.130)
- [x] 5.2 Write `scripts/stage_span_paper_grade_decision.py` (P1–P3 and the verdict order of design D6, unavailable counted as a failure, reusing `stage_span_pilot_reading` loaders), writing `stage_span_paper_grade_decision.json`; verify with fixture reports that each verdict is reachable and that every failure is named
- [x] 5.3 Add the descriptive outputs: target power ± SE, `informative_from`, capture quantiles, rank, intended against `population_latent`, and the visibility contrast from this run's surgery report, the slice's committed `core_geometry` report and the pilot's committed `span_capture.csv` (design D9); verify with tests that missing sources are reported as unavailable with their source labels

## 6. Paper-grade profiles

- [x] 6.1 Write `stage_span_paper_grade_span_geometry.json`, `stage_span_paper_grade_span_magnitude.json` and `stage_span_paper_grade_surgery_geometry.json` (design D4: explicit points, 500 × 999, `design-guidance` family at base seed 1028, span capture, `net_displacement`, `min_components` 2, `joint`, anchor decision acceptance, submission metadata naming the sbatch template), and verify each loads with `motco` study tooling
- [x] 6.2 Append the three profiles to `POST_FLOOR_CONFIGS`, regenerate `tests/data/post_floor_config_identity.json` with `tests/post_floor_fixture_capture.py`, and verify both pin tests pass
- [x] 6.3 Write `tests/test_stage_span_paper_grade_profiles.py`: production measurement named, columns and effects as in the spec, cell and unit counts pinned, sibling anchors drawing the same seeds and parameters, generator seeds equal to the slice's (shared columns), the stage-span pilot's (SEA-AD) and the design-guidance pilot's (four stages) at replicates 0–99, and planted magnitude equal to the pilot's stage-span magnitude; verify the module passes
- [x] 6.4 Add the profiles to `examples/trajectory_power_study/README.md` (purpose, columns, effects, cost, siblings), and verify the README's config list matches the directory

## 7. Local verification

- [x] 7.1 Run the pre-commit gate (`uv run ruff check src/ tests/ && uv run mypy src/motco/ && MOTCO_TEST_PERMS=99 uv run pytest tests/ -m "not slow" --tb=short`), and verify it passes
- [x] 7.2 Run a local smoke of all three profiles at 2 replicates × 19 permutations with BLAS pinned, and verify every unit completes, that the planted cells' truth reports converged fixed points, and record per-unit timings for the cost plan

## 8. Pre-flight on `ing`

> **Reordered 2026-11-04 (user decision).** The workstation had to shut down before the leak check could be written, so the user chose to submit the paper grade (section 9) that night *without* the pre-flight gate, under the pre-approved cap of 450 projected core-hours. The leak check (section 4) and its two arms (8.2) run afterwards as a descriptive companion; 8.3's gate is read as a diagnostic, not as a condition for submission. The decision rule (section 5) is unchanged and predeclared.

- [x] 8.1 Push the branch, create a clean detached worktree on `ing` at the pushed commit (symlinked `.venv`, `logs/`, `PYTHONPATH`), check node memory and NFS state, and verify `motco.__file__` and the profile sha256 match the workstation
- [x] 8.2 Run the leak check's two arms on the stage-span geometry profile (planted at 0.1/0.25/0.5 and null at 0.25/0.5/1.0; anchor, orientation, shape; 100 × 199), merge, and verify 0 failed units and that the null arm's rates equal the 2026-11-04 pilot's at its three columns
- [x] 8.3 Run `stage_span_leak_check.py gate` and verify `preflight_gate.json`: if any condition fails, stop and report the tables and failing cells to the user before any paper-grade submission
  > **Pre-flight outcome (2026-11-05, run after the paper grade).** Jobs 934447 (planted, e 0.1/0.25/0.5) and 934448 (null, e 0.25/0.5/1.0), 100 × 199 on the five columns, revision b40fb91: 3,500 + 3,500 units, 0 failed, 0 duplicates. The null arm reproduces the 2026-11-04 pilot's rejection counts at its three columns (63/63). `preflight_gate.json`: **`gate_holds`** — G1 exact in every planted unit, G2 0 of 60 exceedances, G3 shrinks at every column (score-coordinate median B − A at shape e = 0.5: null −0.28 / −0.29 / −0.36 / −0.35 / −0.92 → planted +0.01 / −0.01 / −0.16 / +0.05 / −0.33 for n200 / n80 / seaad / three blocks / four stages).
- [ ] 8.4 Add the score-coordinate section to `docs/measurement-visibility.md` (coordinates, size unaffected, specificity affected, the pre-flight's residual per column, orthonormal coordinates as the unvalidated candidate) and a roadmap item for that candidate in `docs/roadmap.md`, and verify `uv run mkdocs build --strict` passes

## 9. Paper-grade run

> **Run record (2026-11-04).** Local smoke (2 × 19, 192 units): 192 completed, 0 failed, 0 duplicates; planted fixed points converged in ≤ 19 updates with residual < 1e-12; median capture 0.96 (geometry) and 0.97 (magnitude). Timing on `ing` (job 933792, 32 random units at 999 permutations, n1): 32 completed, mean 17.3 s, median 15.2 s, max 23.0 s → about 230 core-hours (276 with a 20 % margin), under the user's pre-approved cap of 450. Submitted from the detached worktree `/home1/tgonza/MOTCO-stage-span-paper-grade` at e41c1ae (profile sha256s matched the workstation; `motco` 0.17.0 imported from the worktree): 933810 span geometry (260 shards, throttle 64), 933811 span magnitude (110, 28), 933812 surgery geometry (110, 28), one CPU and 4 GB per task, `--exclude=n2,n11` (n2 in `completing` state, as when another user's job hung there earlier that day; n11 at 1.8 GB free memory with a history of NFS stalls). Output: `results/stage-span-paper-grade-2026-11-04/` in that worktree.

- [x] 9.1 Time a few units of each profile on `ing` at 999 permutations, write the shard and throttle plan with the core-hour estimate, and get the user's approval before submitting
- [x] 9.2 Submit the three arrays with `scripts/motco_study_array.sbatch` (n2 and n11 excluded with the reason stated), monitor shard progress and stalls, and verify `sacct` shows every task `COMPLETED 0:0`
  > **Outcome (2026-11-05).** `sacct`: all 480 tasks `COMPLETED 0:0` (22:25–01:13). 48,000 records; 12 failed units, all span geometry at four stages, `orientation`: e = 0.25 replicates 311, 409, 481 and e = 0.5 replicates 123, 266, 311, 368, 409, 425, 439, 477, 481 (`SemiSyntheticTrajectoryError`: the planted fixed point did not converge). Diagnosis: the configuration with the net displacement removed has two nearly tied singular values (σ₂/σ₁ 0.95–0.97), so the rotation's second direction swings with the scale and the iteration cycles in two steps; damping by half-steps rescued only 2 of the 12, so no consistent planted plant probably exists for the others. A local scan of all 1,000 units in the two cells finds exactly these 12. **User decision (2026-11-05): report them as failures** (design risk "the fixed point fails to converge"); those cells' rates use the completed replicates (R = 497 and 491), and the excluded near-tie configurations are a stated limitation. The decision driver reads R per row.
- [x] 9.3 Merge each profile and generate its `report/`, and verify the expected unit counts, no duplicates and no failed records

## 10. Reproduction and decision

- [x] 10.1 Run `scripts/anchor_reproduction.py --by-column` across the three siblings and against the 2026-10-29 slice's anchors, and verify identical seeds, ranks, statistics and p-values at every shared column
- [x] 10.2 Run `scripts/anchor_reproduction.py --by-cell` for the surgery cells against the stage-span pilot (SEA-AD) and the design-guidance pilot (four stages), and the planted magnitude cells against the pilot's stage-span magnitude, at replicates 0–99, and verify identical seeds, ranks and statistics, with any missing source records reported as unavailable
- [x] 10.3 Run `scripts/stage_span_paper_grade_decision.py` on the merged reports, and verify `stage_span_paper_grade_decision.json` names its verdict, each criterion's m, p₂, allowed and observed counts and ceiling, and the descriptive tables

## 11. Report and record

- [ ] 11.1 Copy the aggregates (three `report/` directories, the decision, the pre-flight tables and gate, reproduction summaries, `PROVENANCE.txt`, `sacct.txt`) to `results/stage-span-paper-grade-<date>/`, and verify no record JSONL is committed
- [ ] 11.2 Write `docs/reports/stage-span-paper-grade-<date>.md` (verdict and criteria, power and `informative_from` per column, capture, visibility contrast, both leak arms and the score-coordinate residuals, reproduction, limitations), and verify every number traces to a committed aggregate
- [ ] 11.3 Update `docs/roadmap.md` (item 13(a) done with its verdict, the "Not yet established" entry, the score-metric item), `AGENTS.md` (profiles, scripts, run summary) and `docs/measurement-visibility.md` (paper-grade numbers), and verify `uv run mkdocs build --strict` and the pre-commit gate pass

## Workflow follow-up

- Open a pull request with the problem statement, the summary, the commands run and a link to this change.
- Sync the delta specs and archive the change after review; update the latent-space-comparison memory's open points if the run resolves any.
- Remove the `ing` worktrees (`MOTCO-bendleak` from the proposal's scratch check, and this change's) once the records are no longer needed.
