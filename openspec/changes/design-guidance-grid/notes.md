# Implementation notes — design-guidance-grid

Working notes kept beside the change artifacts (tasks 1.1, 1.2, 5.1 and 5.2 ask for them). Dates are 2026-10-28 unless stated.

## 1. Baseline on the unmodified tree (`843f057`)

- Fast suite: 1459 passed, 1 skipped (`tests/test_trajectory.py:435`, R `shapes` not available), 32 deselected (slow). Both pin modules pass with the 19 historical configs (78 tests).
- Pre-change report capture: the two-axis design-grid study of `tests/test_study_end_to_end.py::test_end_to_end_design_grid` (mock evaluator, 22 cells, 44 records) reported through `scripts/motco_study.py report`; 15 report files hashed (sha256), deterministic across two report passes, `design_point_operating.csv` included. The hash list lives in the session scratchpad (`prechange_report_hashes.json`).
- After the stage-signal change (3.2): every one of the 15 files is byte-identical and `design_point_stage_signal.csv` is the only new file.

## 2. Profiles and pins

- The four profiles enumerate to 98 / 162 / 68 / 112 cells (9,800 / 16,200 / 6,800 / 11,200 units), as the design table predicted; no headroom failure.
- Design D2 check: the n = 80, four-stage, three-block anchor draws generator seed 1476637963 at replicate 0 in both core profiles (`power_design-681fed7db322` in the magnitude profile, `power_design-451e9f2cf50f` in the geometry profile).
- `tests/data/post_floor_config_identity.json` captured twice (byte-identical) at revision `843f057` with the four profiles uncommitted; `POST_FLOOR_CONFIGS` in `tests/post_floor_fixture_capture.py` lists them.
- A stray `examples/trajectory_power_study/tmp.json` absent from both fixtures fails both coverage tests by name; deleted after the check.

## 3. Docs build

`uv run mkdocs build --strict` aborts with 8 warnings on this tree **and on a clean `843f057` worktree** (3 in `phase5-readiness.md`, 3 in `reports/phase5-magnitude-remeasurement-2026-10-08.md`, 2 in `reports/phase5-paper-grade-2026-09-10.md`: links to `results/*/PROVENANCE.txt`, which are not documentation files). This change adds no warning; fixing the pre-existing links means editing frozen reports and is out of scope.

## 4. Smoke run (task 5.1)

Scratchpad copies of the four profiles at 2 replicates × 19 permutations (never committed), shards 5 / 8 / 3 / 6 (core magnitude / core geometry / arms magnitude / arms geometry).

**Local leg (workstation, 22 concurrent shards, BLAS pinned), stopped at the user's request because it saturated the machine:** 476 of 880 units completed, **0 failures** (31 of them at n = 80 with four stages, the smallest cells). After the stop, merge → report → figure script ran on the partial records: every report wrote `design_point_operating.csv` and `design_point_stage_signal.csv`, and `scripts/design_guidance_figure.py` wrote its four figures and CSV. Unit runtimes under that contention (19 permutations; medians): n = 80: 16 s (3 stages) / 29 s (4 stages); n = 200: 21 / 34 s; n = 500: 36 / 43 s; n = 1000: 50 / 73 s (max 97 s).

**Cluster leg (`ing`, 2026-10-28):** the shard files were rsynced to `/home1/tgonza/MOTCO-design-guidance-scratch/smoke/` and resumed with the same shard counts from a clean worktree of `feat/design-guidance-grid` at `6b3bd31` (`/home1/tgonza/MOTCO-design-guidance`, reusing the main checkout's `.venv` with `PYTHONPATH` on the worktree's `src`; the main checkout there is a dirty rsync of older work and was not touched). Jobs 929135–929138 (smoke arrays) and 929139 (`time_units.py`: replicate 0 of the anchor and of orientation at e = 1.0 at 199 permutations for n ∈ {80, 200, 500, 1000} × {3 stages / 2 blocks, 4 stages / 3 blocks}). All 23 tasks `COMPLETED 0:0`, no traceback in any shard log.

**Smoke result:** 880 / 880 units, **0 failures in every column** (the n = 80, four-stage column included). Each profile's report wrote `design_point_operating.csv` and `design_point_stage_signal.csv` (the arms geometry report also wrote `continuity_resolved_orientation.csv`, since that grid sweeps ρ); `scripts/design_guidance_figure.py` wrote its four figures and `design_guidance_power.csv` from the two core reports. The engine ran `n_stages` as a design axis without any fix.

**Stage signal at the smoke (2 replicates per cell, read from `design_point_stage_signal.csv`):** the held-out CV stage AUROC is 1.000 at every core column (δ = 2.0) and at every arm column except the joint stage effect δ = 0.25 (mean 0.988, q33 0.939), where the selected rank also jumps (median 13.8 against 2–4 elsewhere). The SEA-AD regime (AUROC 0.64–0.70) is therefore not reached by the δ arm as specified ({0.25, 0.5, 1.0}); the rank-floor check placed that regime near δ = 0.05–0.1 at n = 80. Raised with the user at the approval gate (see §5).

## 5. Anchor reading and compute plan (task 5.2)

**Duplicated anchors (design D2 open question).** At every core column (32 anchor replicate pairs) and every arm column (22 pairs), the magnitude and geometry profiles' anchors draw the same generator seed **and have identical p-values**, so they share the RRPP draws as well as the data (`evaluation.seed` = 3 is the RRPP seed and the seed family is shared). Three arm pairs differ in observed statistics at the 1e-12 level only: one member ran on the workstation and the other on `ing`, so this is cross-machine floating-point noise, not a different measurement. Reading: the two anchors are one measurement; the report reads the geometry profile's and uses the magnitude profile's as a reproduction check, never pooling them. The pilot should run entirely on `ing` so the reproduction check is bit-exact.

**Unit cost at 199 permutations on `ing` (EPYC, single CPU, BLAS pinned; mean of the anchor and orientation e = 1.0 units):**

| n | 3 stages, 2 blocks | 4 stages, 3 blocks |
|---:|---:|---:|
| 80 | 10.7 s | 14.1 s |
| 200 | 13.2 s | 19.7 s |
| 500 | 19.4 s | 29.4 s |
| 1000 | 27.6 s | 42.9 s |

The 19-permutation smoke medians on `ing` are 7.9–37.6 s at the same columns, so the RRPP's 199 permutations add only 2–5 s per unit: the PLS double CV dominates. Unmeasured columns (3 stages / 3 blocks, 4 stages / 2 blocks, every arm) are scaled from the measured column at the same n by their smoke ratio (cluster-leg records only).

**Estimate (`compute_plan.py`, 20 % contention margin, shards sized to ≈ 2.4 h each):**

| Profile | Units | Core-hours | Shards | n = 1000 share |
|---|---:|---:|---:|---:|
| core magnitude | 9,800 | 73 | 30 | 29 h |
| core geometry | 16,200 | 121 | 49 | 49 h |
| arms magnitude | 8,000 | 41 | 17 | — |
| arms geometry | 13,200 | 68 | 28 | — |
| **total** | **47,200** | **303** | **124** | 78 h (26 %) |

(The arm rows include the δ = 0.05 and 0.1 columns; their unit cost is scaled by the costliest smoked arm's ratio to the baseline, since the smoke never ran them. The draft arms — 6,800 / 11,200 units — came to 34 / 56 core-hours and 283 in total.) Below the proposal's 330 core-hour upper bound. The n = 1000 columns are about a quarter of the cost, not half, so they need no array of their own. 124 single-CPU tasks fit the `512x1024` partition (416 CPUs, ~250 idle at submission time), giving about 2.5 h wall if they all start together. Submission as in the examples README (`--mem=4G --time=6:00:00`, BLAS pinned, `STUDY_N_JOBS` unset, `--error-policy record`).

## 6. Pilot submission (task 6.1)

Submitted 2026-10-28 from the `ing` worktree at `7c218c6` (`scripts/motco_study_array.sbatch`, partition `512x1024`, 1 CPU, 4 GB, 6 h, BLAS pinned, `STUDY_N_JOBS` unset, `--error-policy record`), `STUDY_OUT=/home1/tgonza/MOTCO-design-guidance/results/design-guidance-pilot-2026-10-28/<profile>`:

| Profile | Job | Shards | Config sha256 (prefix) |
|---|---|---:|---|
| core magnitude | 929160 | 30 | `c2bd261307a10bd3` |
| core geometry | 929161 | 49 | `57942119d892ba2b` |
| arms magnitude | 929162 | 17 | `e33e8c4eb6ec3f60` |
| arms geometry | 929163 | 28 | `be4e50300d54edee` |

All 124 tasks were running within a minute of submission.

## 7. Resuming from another machine

Everything needed lives on `ing` or in this branch; nothing on the workstation is required.

- **Code:** branch `feat/design-guidance-grid` (pushed). On `ing` it is checked out as the worktree `/home1/tgonza/MOTCO-design-guidance` (symlinked `.venv` from the main checkout; export `PYTHONPATH=/home1/tgonza/MOTCO-design-guidance/src` so the new code is imported). Do not touch `/home1/tgonza/MOTCO` (a dirty rsync of older work).
- **Run:** jobs 929160–929163 writing to `/home1/tgonza/MOTCO-design-guidance/results/design-guidance-pilot-2026-10-28/<profile>/`. Check with `squeue -u tgonza` and `sacct -j 929160,929161,929162,929163 -X`.
- **Tooling** (session scripts, copied to `/home1/tgonza/MOTCO-design-guidance-scratch/`): `ing_finish_pilot.sh` (merge + report every profile, figures from the core reports, `PROVENANCE.txt` per profile via `pilot_provenance.py`, sacct wall summary), `pilot_tables.py <run dir>` (markdown tables for the dated report from the committed CSVs; needs pandas, run it locally after the rsync), `compute_plan.py` / `compute_plan.json` / `timing_199.json` (the 5.2 estimate), `smoke_analysis.py`, `write_profiles.py` (the generator of the four configs), `prechange_report_hashes.json` (the 1.2 hash list).
- **After the jobs finish (tasks 6.1–6.2):**
  1. `ssh ing 'bash /home1/tgonza/MOTCO-design-guidance-scratch/ing_finish_pilot.sh'` — check "ACCOUNTING PROBLEM" never prints, every array task is `COMPLETED 0:0`, record counts equal 9,800 / 16,200 / 8,000 / 13,200 minus failures, no duplicate (cell, replicate).
  2. `rsync -a --include '*/' --include 'report/**' --include 'PROVENANCE.txt' --include 'figures/**' --exclude '*' ing:/home1/tgonza/MOTCO-design-guidance/results/design-guidance-pilot-2026-10-28/ results/design-guidance-pilot-2026-10-28/` (shard and merged JSONL stay on the cluster, gitignored); also fetch `merged.jsonl` per profile if a reading needs the records.
  3. `python pilot_tables.py results/design-guidance-pilot-2026-10-28` → `tables.md`, then fill §1–5 and §7 of `docs/reports/design-guidance-pilot-2026-10-28.md` (§6 and §8 are written), update `docs/roadmap.md` (Next changes entry + Phase 5 sensitivity deliverables) and the `study/` paragraph of `AGENTS.md` with the result line (task 7.2), run the gate (7.3), commit `results/design-guidance-pilot-2026-10-28/` and the docs.
- **Workstation quirk:** this session's shell lost zlib from `LD_LIBRARY_PATH`; `uv run` needed `LD_LIBRARY_PATH=/nix/store/61a1nwx3w6rqyaisj5rn1sal1981apm7-zlib-1.3.2/lib:$LD_LIBRARY_PATH`. Irrelevant elsewhere.

**Approval gate (2026-10-28).** The user approved submitting all four profiles and chose to **add δ = 0.05 and 0.1** to the arm profiles (the smoke's stage-signal finding above): spec, design D4, proposal, tasks, both arm profiles, `tests/test_design_guidance_profiles.py` and the post-floor fixture (re-captured, two captures identical) were updated; the arm profiles now enumerate to 80 / 132 cells. The unchanged core profiles keep their fixture entries.
