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

## 4. Local smoke run (task 5.1)

_Pending: filled in when the smoke run completes._

## 5. Anchor reading and compute plan (task 5.2)

_Pending._
