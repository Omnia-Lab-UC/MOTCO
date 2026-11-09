# Notes

## Smoke and compute projection (task 5.1, 2026-11-09)

Run on `ing` from the clean worktree `/home1/tgonza/MOTCO-orthonormal-coordinates` at `d776366` (detached; `.venv`
shared with the main checkout, `PYTHONPATH` on the worktree's `src`). Versions: Python 3.11.16, numpy 2.3.5,
scipy 1.16.3, scikit-learn 1.8.0, motco 0.19.0. The change's test modules plus the phase-6, analysis and committed-config
identity suites passed on a compute node (95 passed, job 936028). Outputs:
`/home1/tgonza/MOTCO-orthonormal-coordinates-scratch/smoke-2026-11-09/` (`p19/`, `p999/`, `timing.json`).

The smoke was taken on `ing` rather than the workstation, because the workstation was busy with a system build
(load average 47 on 16 cores).

| Array | Profile | Grade | Units | Job | State |
|---|---|---|---|---|---|
| p19 | `paper_grade` | 2 reps × 19 perms, 26 shards | 130 | 936029 | 26 × COMPLETED |
| p19 | `weak` | 2 reps × 19 perms, 2 shards | 4 | 936030 | 2 × COMPLETED |
| p999 | `paper_grade` | 1 rep × 999 perms, 65 shards (one unit per cell) | 65 | 936031 | 65 × COMPLETED |
| p999 | `weak` | 1 rep × 999 perms, 2 shards | 2 | 936032 | 2 × COMPLETED |

Tasks ran on n1 (61), n2 (16) and n3 (18), beside other users' jobs.

**Mechanics.** 201 units, 0 failed, 0 inconsistent, 0 duplicates. The largest consistency deviations were
1.3e-15 (orthonormal basis) and 1.6e-13 (same span), against the 1e-8 tolerance. `merge` wrote every table in all
four directories. Both reproductions report every row differing (195 of 195 and 6 of 6). The only cause is the
smoke's size: every reproduced row exists, with `n_available` equal to the smoke's 2 or 1 replicates against the
parents' 500. `decide` on `p999` read both profiles and returned `reproduction_failed`, as it must at this size.

**Timing at 999 permutations** (`p999`, seconds per unit, from the records' `runtime_seconds`):

| Profile | Generate | PLS fit (double CV) | `scores` arm | `orthonormal` arm | Unit mean | Unit median | Unit max |
|---|---|---|---|---|---|---|---|
| `paper_grade` (65) | 0.44 | 19.0 | 0.81 | 0.79 | 21.1 | 21.4 | 29.0 |
| `weak` (2) | 0.35 | 8.0 | 0.42 | 0.40 | 9.1 | 9.1 | 9.2 |

By mode (`paper_grade`, mean): anchor 21.7, translation 21.2, orientation 22.1, shape 21.8, magnitude 18.9. The
second arm costs about 0.8 s, about 4 % of a unit. The design's estimate of 2 s per RRPP pass at 999 permutations was
high.

**Projection.**
- **Unpacked:** each smoke unit stands for its cell's 500 replicates, giving 190.1 core-hours for `paper_grade` and
  2.5 for `weak`, 192.7 in all.
- **Calibrated against the paper grade's saturated run:** the stage-span paper grade recorded mean units of 22.3 s
  (geometry) and 12.5 s (magnitude), for 229 core-hours. On the same cells, this smoke's scores-only cost (unit minus
  the `orthonormal` arm) is 20.9 s and 18.2 s, ratios 1.07 and 0.69. Applying the larger ratio gives about
  **206 core-hours**, inside the design's 170–230.

**Proposed cap: 300 core-hours** of recorded time, about 1.5 × the calibrated projection. It leaves room for node
packing, which once slowed n1 units 2.5× (ADNI-design check, 2026-11-02). If the recorded core-hours reach 300 before
the arrays drain, stop submitting, merge what exists and report.

**Proposed submission** (from the worktree, after the user approves the cap):
- `paper_grade`: 650 shards of 50 units (about 18 min each), `--array=0-649%96`.
- `weak`: 20 shards of 50 units (about 8 min each), `--array=0-19%20`.
- n11 is excluded by the sbatch template.
- Wall time is about 2.5 h at 96 concurrent tasks, plus Python startup.
- Settings are pre-created with `run --limit 0` and `logs/` exists in the worktree.
- Out-dir: `/home1/tgonza/MOTCO-orthonormal-coordinates-scratch/run-<date>/`, outside the worktree so its
  untracked output never blocks a fast-forward. Records stay on the cluster; aggregates are copied for commit.
