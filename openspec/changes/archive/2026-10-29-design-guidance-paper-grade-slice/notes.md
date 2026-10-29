# Implementation notes — design-guidance-paper-grade-slice

Working notes kept beside the change artifacts (tasks 1.1, 1.2, 5.1 and 5.2 ask for them). Dates are 2026-10-29 unless stated.

## 1. Baseline on the unmodified tree (`d7836d3`)

- Fast suite: 1510 passed, 1 skipped (`tests/test_trajectory.py:435`, R `shapes` not available), 32 deselected (slow). The three pin modules pass with the 19 historical and 4 post-floor configs (91 tests).
- Pre-change capture (`capture_identity.py` in the session scratchpad, 49 files hashed): a two-axis design-grid study (the `tests/test_study_end_to_end.py::test_end_to_end_design_grid` construction with `magnitude`, `orientation` and `translation`, effects 0 / 0.5 / 1.0, four replicates, two Type I entries without a mode — one two-sided, one `upper` — plus the anchor decision and the design-point rule) run with the mock evaluator and reported through `scripts/motco_study.py report` (22 report files, `acceptance_report.csv`/`.json`, `anchor_decision.*` and `design_point_*` among them); the loaded config echoed by `dump_study_config`; `scripts/anchor_reproduction.py` on a synthetic single-anchor pair (CSV, stdout, exit code); `scripts/design_guidance_figure.py` on the `tests/test_design_guidance_figure.py` fixtures (four PNGs and the CSV).
- After the engine edits (task 2), after the script edits (task 3) and in the final task-3 state: all 49 files byte-identical to the capture.

## 2. Engine additions (tasks 2 and 3)

- **Column-scoped Type I entries** (`config.py`, `targets.py`): `TypeIControlTarget` gains `trajectory_mode` / `statistic`, validated in `__post_init__` and re-raised by the loader as `acceptance.type_i[<i>]: …`; `_evaluate_type_i_columns` buckets completed records with `record_design_point` exactly as `evaluate_anchor_decision` does. The synthetic six-column test asserts the anchor evaluations agree with the anchor decision at every column (same rate, R, bound and classification); the design construction puts one anchor over the bound (n = 500, three blocks, 0.080 against 0.0695) and both readers flag the same column.
- **Figure grade label** (`design_guidance_figure.py --grade`): the figure builders are split from the savers (`build_mode_figure`, `build_control_figure`) so the tests read the suptitles and panel counts; the saved bytes are unchanged (capture).
- **Anchor reproduction by column** (`anchor_reproduction.py --by-column`): keyed by `cell_metadata.design_point`; `type_i_baseline` excluded; `main(argv)` for the tests; the single-anchor path is untouched and its CSV columns and summary lines are pinned in `tests/test_anchor_reproduction.py`.

## 3. Profiles and pins (task 4)

- The four profiles are generated from the pilot core profiles by `write_profiles.py` (session scratchpad; the JSON is the committed artifact). They enumerate to 26 / 62 / 14 / 32 cells (13,000 / 31,000 / 7,000 / 16,000 units), as the design table predicted; no headroom failure.
- Seed reproduction checked by enumeration: every slice design cell (24 / 24 and 60 / 60) draws the pilot core cell's generator seeds at replicate indices 0–99; the two Type I baseline cells per profile do not. The stage-effect arm's design cells (12 / 12 and 30 / 30) draw the pilot **arms** cells' seeds at 0–99 too — the pilot arms profiles placed δ ∈ {0.25, 0.5, 1.0} as explicit points off the same core baseline — a fact the design did not claim; it is recorded in the arm profiles' metadata and asserted by `tests/test_design_guidance_slice_profiles.py`, so the pilot reproduction check applies to the arm as well.
- In the magnitude profiles the column-scoped translation entry evaluates every column as unavailable (translation runs only in their Type I baseline phase); it is declared so the sibling acceptance blocks are identical and the geometry sibling carries the verdict (spec: every slice profile declares it).
- `tests/data/post_floor_config_identity.json` captured twice (byte-identical) with the four profiles uncommitted; the four pilot entries are unchanged; 8 configs pinned; the three pin modules pass (103 tests).

## 4. Docs build

`uv run mkdocs build --strict` aborts with the same 8 warnings before and after every doc edit of this change (3 in `phase5-readiness.md`, 3 in `reports/phase5-magnitude-remeasurement-2026-10-08.md`, 2 in `reports/phase5-paper-grade-2026-09-10.md`: links to files outside `docs/`). No new warning.

## 5. Smoke run (task 5.1)

Scratchpad copies of the four profiles at 2 replicates × 19 permutations (never committed), 4 shards each, 16 concurrent shards on the workstation (BLAS pinned), then `motco_study.py merge` → `report` → `design_guidance_figure.py --grade smoke` → `anchor_reproduction.py --by-column`. Shards done in about 11 min; the whole pipeline in 714 s.

- **268 / 268 units, 0 failures** (52 / 124 / 28 / 64 for core magnitude / core geometry / stage-effect magnitude / stage-effect geometry); no duplicated `(cell_id, replicate_index)`; every record carries `n_jobs` 1, `orientation_estimator` `net_displacement` and `min_components` 2; selected ranks 2–3 at the core columns (δ = 2.0) and 3–19 at the arm columns (δ ≤ 1.0, the rank inflation the pilot saw).
- Every report wrote its 16 files, among them `acceptance_report.csv` with the **column-scoped translation rows**: six per slice profile and three per arm profile. In the geometry profiles they are judged (`met` True at every column: translation `angle` 0.0 over 2 replicates against the bound 0.358 at R = 2); in the magnitude profiles all are unavailable (`met` empty, `available_replicates` 0, name `type_i_control[translation, angle, {…}]`), as predeclared — translation runs only in their Type I baseline phase. `anchor_decision.csv` (`SIZE_HOLDS` everywhere at R = 2), `design_point_operating.csv` and `design_point_stage_signal.csv` were written by all four.
- Stage signal at the smoke: held-out CV stage AUROC 1.000 at every core column; 0.996 at δ = 0.25 (0.94–0.99 on the orientation cells there), consistent with the pilot's δ arm.
- `design_guidance_figure.py --grade smoke` wrote the four figures and `design_guidance_power.csv` from the two slice reports (one panel row, two panels per mode figure).
- Sibling checks: core magnitude vs core geometry — 12 pairs over 6 columns, seeds, ranks, statistics and p-values identical at every column (exit 0); stage-effect magnitude vs geometry — 6 pairs over 3 columns, identical (exit 0).
- No engine failure, so nothing to fix. One environment note: the workstation's numpy needs zlib on `LD_LIBRARY_PATH` from a non-interactive shell, and the shell lost it mid-session, so the pipeline's last two steps (the sibling checks) failed on import inside the smoke script and were rerun by hand with `/run/current-system/sw/share/nix-ld/lib` prefixed; nothing MOTCO-related.

## 6. Unit timings and compute plan (task 5.2)

Branch `feat/design-guidance-paper-grade-slice` pushed at `23bc09e`; clean worktree `/home1/tgonza/MOTCO-design-guidance-slice` added beside the main checkout on `ing` (the main checkout's `.venv`, `PYTHONPATH` on the worktree's `src`); the change's six test modules pass there (115 tests, python 3.11.16, numpy 2.3.5, scikit-learn 1.8.0, scipy 1.16.3; the venv's `motco` metadata still says 0.11.0, so `PROVENANCE.txt` cites the revision). Session tooling in `/home1/tgonza/MOTCO-design-guidance-scratch/slice/` (`time_units.py`, `ing_timing.sbatch`, `timing_999*.json`).

**Unit cost at 999 permutations** (`time_units.py` via `run_simulation_replicate`, replicate 0 of the anchor and of orientation e = 1.0 at each column, single CPU, BLAS pinned; mean of the two units):

| Column | n4, job 929843 (load 21 / 128) | n1, job 929834 (load 99, 126 / 128 allocated) | selected ranks |
|---|---:|---:|---|
| n = 80, 2 blocks | 8.3 s | 16.4 s | 3, 3 |
| n = 80, 3 blocks | 9.4 s | 18.2 s | 2, 3 |
| n = 200, 2 blocks (baseline) | 10.2 s | 19.5 s | 2, 3 |
| n = 200, 3 blocks | 12.1 s | 22.5 s | 2, 3 |
| n = 500, 2 blocks | 15.5 s | 27.4 s | 2, 3 |
| n = 500, 3 blocks | 19.4 s | 33.4 s | 2, 2 |
| δ = 0.25 arm column (n = 200, 2 blocks) | 10.9 s | 21.1 s | 8, 19 |

Two things the pilot's basis did not show. (1) 999 permutations cost no more than 199 did on a quiet node (the pilot measured 10.7 / 13.2 / 19.4 s at n = 80 / 200 / 500 at 199 permutations): the RRPP on a rank-2–3 latent space is a small part of a unit, the PLS double CV dominates. (2) The per-unit cost tracks node load: 1.8× between n4 and n1 at every column. The first timing job happened to land on n1 and spent 27 min on 14 units, most of it in NFS-bound imports (a bare `import motco` took 77 s wall for 11 s CPU even on the login node); the n4 job took 13 min. The plan therefore sizes shards at the **loaded-node** cost — the core-hours are what they are; more, shorter shards only add array tasks — and states the budget as a range.

**Plan** (`compute_plan.py`: cells × 500 units at the column's cost, the two Type I baseline cells at the baseline column's cost, the arm's three columns at the δ = 0.25 cost; 20 % margin; shards at ≤ 2.5 h at the loaded cost):

| Profile | Units | Core-hours (quiet → loaded) | Planned (loaded + 20 %) | Shards | Per shard (loaded / quiet) |
|---|---:|---:|---:|---:|---:|
| core magnitude | 13,000 | 44 → 82 | 98 | 40 | 2.45 h / 1.3 h |
| core geometry | 31,000 | 107 → 196 | 236 | 95 | 2.48 h / 1.3 h |
| **slice** | **44,000** | **151 → 278** | **334** | **135** | |
| stage-effect magnitude | 7,000 | 21 → 41 | 49 | 20 | 2.46 h / 1.3 h |
| stage-effect geometry | 16,000 | 49 → 94 | 113 | 46 | 2.45 h / 1.3 h |
| **arm** | **23,000** | **70 → 135** | **162** | **66** | |

The design's estimate was about 380 core-hours for the slice (planned at 450) from the pilot's 199-permutation basis scaled up for the permutations; the measured cost is below it at either load. Written into each profile's `metadata.submission` (`timing_basis`, `n_shards`, `array`, `budget`, `timing_jobs`). Submission recipe (after the user's approval; `scripts/motco_study_array.sbatch`, partition `512x1024`, 1 CPU, 4 GB, 6 h, BLAS pinned, `STUDY_N_JOBS` unset, `--error-policy record`, from the worktree, `STUDY_OUT` under `results/design-guidance-slice-<date>/<profile>`): see `examples/trajectory_power_study/README.md`, "Design-guidance paper-grade slice".

## 7. Paper-grade run (task 6.1)

Submitted 2026-10-29 at about 11:45 local after the user's approval of the plan in §6, by `/home1/tgonza/MOTCO-design-guidance-scratch/slice/ing_submit_slice.sh` from the worktree at `8bf8eb7` (0 dirty files; import check confirmed the worktree's code): `sbatch -p 512x1024 --cpus-per-task=1 --mem=4G --time=6:00:00 --export=ALL,STUDY_CONFIG=…,STUDY_OUT=…,N_SHARDS=…` with `PYTHONPATH` on the worktree's `src`, BLAS pinned, `MPLBACKEND=Agg`, `STUDY_N_JOBS` unset, `--error-policy record` (the template's default).

| Profile | Job | Shards | Config sha256 (first 16) | Output |
|---|---|---:|---|---|
| core magnitude | 929844 | 40 | `3f4a618d57fc59db` | `results/design-guidance-slice-2026-10-29/core_magnitude` (in the worktree) |
| core geometry | 929845 | 95 | `bb56f1c6f61ce2e3` | `results/design-guidance-slice-2026-10-29/core_geometry` |

First tasks started at once on n1, n4 and n11 (n2 down, n1 saturated by other users' jobs).

**Placement incident (12:13–12:20).** SLURM put 35 of the 40 magnitude shards on n11 (load 100–110 on 128 CPUs) and they sat in `nfs_wait_bit_killable` at 3–6 % CPU for 32 minutes without finishing their imports or writing a unit, while the five magnitude shards on n1/n4 ran at 17–31 s per unit and the 95 geometry shards on n3/n4 at 10 s per unit. The 35 tasks were cancelled (`scancel 929844_<id>`; no `shard_<id>.jsonl` existed for any of them) and the same array ids resubmitted as **job 929980** with `--exclude=n11` from the same worktree and environment: 19 started on n3 at once, 16 pended until cores freed. Shard assignment is by the `(cell_id, replicate_index)` hash, so the node is irrelevant to the records; `PROVENANCE.txt` lists both job ids and the cancelled tasks. The arm's shards, if run, should exclude n11 unless its load has dropped.

**Geometry profile finished 13:34** (95 / 95 tasks `COMPLETED 0:0`, 31,000 records, 0 tracebacks, every cell at 500 available replicates; interim merge + report run at 13:36, anchor decision `SIZE_HOLDS` at all six columns, the column-scoped translation `angle` rule met at all six at 0.028–0.052 against 0.0695).

## 8. Stage-effect arm (task 6.3)

Submitted 2026-10-29 at about 13:45 after the user's decision (asked with the geometry numbers in hand and the magnitude profile still running), by `ing_submit_arm.sh` from the worktree at `8bf8eb7` with `--exclude=n11` (the placement incident above), otherwise exactly as the slice:

| Profile | Job | Shards | Config sha256 (first 16) |
|---|---|---:|---|
| stage-effect magnitude | 930041 | 20 | `2cd406d7d6238ff0` |
| stage-effect geometry | 930042 | 46 | `961bac86a2e2a6c8` |

## 9. Outcome (tasks 6.2, 6.3, 7.1)

- Slice finish (`ing_finish_slice.sh`, 14:30–14:50): both profiles `accounting ok` — every array task `COMPLETED 0:0`, 13,000 / 31,000 records, 0 duplicates / missing / unexpected / signature mismatches / failed / censored, `n_jobs` 1, `net_displacement`, `min_components` 2 in every record; 52.7 + 119.0 = 171.7 recorded core-hours (median 11.9 / 12.1 s per unit) against the 151 (quiet) – 278 (saturated) estimate; anchor decision `SIZE_HOLDS` at all six columns in both; figures `--grade paper`; siblings 3,000 / 3,000 with identical p-values; pilot 600 / 600 per profile at replicates 0–99 (p-values not comparable), the ten pilot-only columns listed with zero pairs.
- Arm finish (15:13–15:35): both profiles `accounting ok`, 7,000 / 16,000 records, 0 failures, 20.5 + 47.6 = 68.1 core-hours (median 10.5–10.7 s per unit; the inflated ranks did not raise the unit cost at n = 200), `SIZE_HOLDS` at all three δ; arm siblings 1,500 / 1,500 identical p-values. The committed `--by-column` tool paired **0** columns against the pilot arms (its key is the full design point; the pilot arms' points carry seven axes), so the arm-vs-pilot reproduction was done with a projected join on (δ, replicate index) restricted to the pilot columns at the baseline of the other axes (`arm_pilot_projected.py`, session tooling; `reproduction/pilot_stage_effect_*_projected.csv`). Follow-up for the tool: an option to project the column key onto named axes.
- Committed: `results/design-guidance-slice-2026-10-29/{core_magnitude,core_geometry,stage_effect_magnitude,stage_effect_geometry}/{report/,PROVENANCE.txt}`, `figures/`, `reproduction/`, `tables.md`; the four `merged.jsonl` (43–197 MB) stay on `ing` under the worktree, gitignored. Report: `docs/reports/design-guidance-slice-2026-10-29.md`.

