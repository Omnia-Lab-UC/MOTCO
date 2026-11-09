# Notes

## Baseline (2026-11-09)

- Branch `feat/interaction-feature-scan` from `main` at `03c2614` (PR #83), in the main checkout.
- `src/motco/__init__.py` was `0.19.0`; bumped to `0.20.0` (task 5.3).
- Local runs need the zlib `LD_LIBRARY_PATH` workaround (`/run/current-system/sw/share/nix-ld/lib`) and BLAS pinned
  to one thread for timing.

## Design notes from implementation

- **Shared helpers (task 1.1).** The split test's validation is now `_check_two_dimensional`, `_check_rows`,
  `_check_finite`, `_check_permutations` and `_model_spaces` (column bases, nested check, df, the interaction basis and
  the `[U_int, U_R]` weights). The split test calls them in its old order, so its messages and arithmetic are unchanged;
  `tests/test_interaction_split.py` and the analysis identity fixtures pass unchanged.
- **One code path for observed and null F.** Per draw the scan computes `W[argsort(π)]ᵀ·E` (rank(F) × p), the per-feature
  SS from its first df_int rows and RSS = ‖e‖² − ‖Wᵀe‖², clipped at 0; the observed F use the identity permutation.
  No n × p permuted matrix is formed: the 10,000-feature, 80-sample, 12-column, 999-draw test case runs in about 0.5 s
  on this machine (multithreaded BLAS).
- **A constant residual multiplier is not applied.** Under `df` (and `leverage` at equal leverages) every F* is
  mathematically unchanged; applying κ would still move the last bit of the products. The scan applies only a varying
  multiplier, so "every null F and p-value identical under `df`" holds exactly, not to rounding.
- **Degenerate features.** A column is degenerate when it is constant (exact equality; the centred SS of a constant
  column need not be exactly 0 in floating point) or when its RSS ≤ 1e-12 of its centred SS. An undefined null F (0/0, a
  permuted column inside col(R)) counts as 0, so it never inflates a maximum.
- **Ties and floating point.** The step-down order is F descending, ties by column index. A duplicated column's F can
  differ in the last bit (BLAS blocking), so the tie test rounds F to make the tie exact; the constant-column test
  compares the other columns' F to 1e-10, not bit for bit.
- **The analysis block** names its maximum `max_f` / `max_feature` / `max_block` (flat, so they do not collide with the
  per-hit `feature`/`block`). Per block it always carries `n_features` and `n_degenerate`, and with permutations `hits`
  and `null_max_share`; `p_value_global`, `hit_alpha`, `hits` and `draws_match_rrpp` exist only with permutations. The
  per-feature table omits `p_value_adjusted` without permutations. `analysis._split_feature_names` now serves both the
  loadings and the scan table.
- **Driver settings binding.** Each profile rebuilds the parent's settings through `span_split_check` at the parent's own
  grade (read from the committed file) and requires equality with the committed `settings.json`, naming differing keys;
  the run's own grade (replicates, permutations, designs) is a separate rebuild, so a smoke can run at 2 × 19. The
  sha256 of the parent's `settings.json`, `split_table.csv` and `production_table.csv` / `trajectory_table.csv` are in
  the settings and re-checked by `run`, `merge` and `decide` (`check_parent_files`).
- **Rejection rules.** The scan rejects at p ≤ 0.05 in both profiles (exact at B = 999). The comparator table carries
  each test under its rule: the scan at p ≤ α, the six existing tests at the parent's rule (p < α in the simulation,
  p ≤ α in the plasmode), so the reproduction compares like with like and G (plasmode only) compares p ≤ α with p ≤ α.
- **Plasmode records hold no feature names.** The analysis block's `max_feature` and hit list are dropped; the record
  keeps counts, the maximum's block and per-block shares. The unit checks that the scan table's hit count equals the
  block's before recording.
- **Plasmode scan timing.** The scan runs inside `analyze_trajectories`, so a plasmode record has the unit's wall time
  only; the scan's share there is measured separately in the smoke. Simulation records time the scan directly
  (`runtime_seconds.scan`), and `accounting.json` reports its share of unit seconds.
- **The sbatch template** takes `ISC_DRY_RUN` to print the command it would run (task 4.4's dry run).

## Gate (task 5.3)

`ruff check src/ tests/` and `mypy src/motco/` clean. `mkdocs build`: the same 8 warnings as on `main` (none new).
`MOTCO_TEST_PERMS=99 pytest tests/ -m "not slow"`: 2528 passed, 1 skipped (R `shapes`) in 6 min 31 s (serial; the venv
has no `pytest-xdist`). The sbatch template passes `bash -n`, its `ISC_DRY_RUN` command lines are the expected ones for
both profiles, its guards reject a missing `ISC_EXPORT` and an unknown profile, and one real task
(`--reps 1 --control-reps 1 --perms 9 --designs seaad --limit 1`) completed a `seaad` anchor unit in 11 s (fit 10.8 s,
split 3 ms, scan 1 ms).

## Smoke, timing and projection (task 6.1, 2026-11-09; local, 8-core laptop, BLAS pinned)

- **Fidelity gate** on the local export (`../sead_sexdiff/results/motco_input/cohort_30`): passed, rank 5, AUROC
  0.635556, 74 s.
- **Simulation smoke** (`--reps 2 --control-reps 2 --perms 19 --workers 7`): 84/84 units, 0 failed, 149 s wall; median
  unit 11.3 s (fit 10.9 s, evaluation 0.18 s, split 4 ms, scan 2 ms). The reproduction check ran and failed on all 252
  rows (14 cells × 3 designs × 6 tests), each on count only (2 against the parent's 100 or 500), as the smoke's size
  requires.
- **Plasmode smoke** (`--reps 2 --perms 19 --workers 6`): 6/6 units, 0 failed, every anchor through the latent check;
  median unit 89 s. The reproduction check ran and failed on its 18 rows by count only (2 against 500).
- **`decide`** on the smoke: `reproduction_failed`, the gate working as specified.
- **Timing at 999 permutations.** The scan alone at the export's size (80 × 10,000, 12 model columns, `df`): 0.72 s.
  Plasmode: 3 units (one per cell, 3 workers), 81.6–82.0 s each, so the scan is under 1 % of a unit. Simulation: 60 units
  (anchor and translation × 4, every other cell × 1, three designs, 7 workers), median 11.8 s (fit 10.8 s, evaluation
  0.57 s, split 0.12 s, scan 0.077 s), q90 15.3 s; the scan is 0.7 % of unit seconds.
- **Descriptive, not a measurement:** at 999 draws the single e = 2.0 and e = 3.0 units rejected globally (p = 0.001)
  with 3 and 14 hits, all inside the planted set, and 14 and 18 planted features among the 20 largest F; the anchor unit
  had p = 0.216 and no hit.
- **Projection.** The units are the split-test run's, which spent 25.6 (simulation, median 11.1 s on `ing`) and 36.9
  (plasmode, median 80 s) core-hours. Locally the units take 11.8 s and 82 s, about the same as that run on `ing`. The
  scan adds 6,600 × 0.08 s + 1,500 × 0.7 s ≈ 0.5 core-hours, and the plasmode's uncorrected split rerun is not
  repeated. **Projection about 63 core-hours; budget 75** with 20 % for node contention and NFS start-up.

## Run (tasks 6.2–6.4, 2026-11-09)

- `359f57f` pushed to GitHub and to the cluster clone (`refs/heads/incoming/interaction-scan`); clean detached worktree
  `/home1/tgonza/MOTCO-interaction-scan` (149 change tests passed there via `srun`); fidelity gate passed on the cluster's
  export copy; settings pre-created on the login node (simulation hash equal to the workstation's).
- Arrays 936835 (simulation, 200 shards) and 936836 (plasmode, 150 shards), `%64` each, n11 excluded: 350/350 tasks
  completed in 71 min, 8,100/8,100 units, 0 failed, 69.6 core-hours against the 63 projected (the plasmode units ran at a
  median 86 s on packed nodes, q90 158 s). Plasmode throughput dipped at each wave boundary (the next 64 shards start
  Python over NFS), as in earlier runs.
- Both reproductions passed (252 + 18 rows) and `decide` returned **`scan_validated`** (`PROVENANCE.txt`, report).
- The S3 false hits (21/500 at e = 2.0 and 3.0) are, replicate for replicate, the anchor's 21 global rejections: within
  a replicate the labels and RRPP seed are shared and the plant moves 20 features, so the null features' maxima are the
  anchor's (checked on the per-unit records).
- **Lesson:** a finishing script piped to `ssh ing 'bash -s'` stops at its first `srun`, because `srun` forwards stdin to
  its task and consumes the rest of the script. Merge and decide ran; `sacct.txt` was written by a second call. Use
  `srun ... < /dev/null` or a script file on the cluster.
