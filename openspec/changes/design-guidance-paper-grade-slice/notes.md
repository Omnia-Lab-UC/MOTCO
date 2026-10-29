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

_Pending: `time_units.py` on `ing` from a clean worktree of the pushed branch._
