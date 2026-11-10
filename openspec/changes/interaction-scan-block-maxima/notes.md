# Notes

## Timing projection (task 4.1, 2026-11-10)

This replaces the proposal's estimate of about 10 core-hours.

**Revision timed.** Base commit `c10fc0c` plus the uncommitted implementation of tasks 1.1–3.4:
- `git diff HEAD` sha256 `5be3360ded1b1b45…`;
- new files (sha256 prefixes): `scripts/interaction_scan_blocks_check.py` `9bc9ea29ba0bb901`, `scripts/interaction_scan_blocks_check_array.sbatch` `0b5cd21f78cd32f5`, `tests/test_interaction_scan_blocks_check.py` `979150cc0241fc77`, `tests/test_analysis_scan_by_block.py` `d40ff8a86d6f985b`.

MOTCO 0.23.0; numpy 2.3.5, scipy 1.16.3, scikit-learn 1.8.0, Python 3.11.15. Workstation AMD Ryzen 7 5800X, BLAS pinned to one thread.

**Smoke** (`/tmp/isb-smoke`, off the repository, 99 permutations):

| profile | units | failures | wall |
| --- | --- | --- | --- |
| `simulation` | 84 (2 per cell, 3 designs × 14 cells) | 0 | 5 s, 8 workers |
| `plasmode` | 6 (2 per cell) | 0 | 3 s, 6 workers |
| `plasmode_unequal` | 6 (2 per cell) | 0 | 2 s, 6 workers |

- The local export `sead_sexdiff/results/motco_input/cohort_30` passed the plasmode's hash check. The settings rebuilt equal to the committed `results/interaction-scan-2026-11-09/<profile>/settings.json`.
- `merge` wrote every table. The reproduction failed, as it must at 2 replicates (156 and 11 differing rows), and `decide` returned `reproduction_failed`.
- The unequal subset is 500 of 5,000 RNA features, seed 455583459, feature-list sha256 `6d6ee1a86b618de1…`.

**Timing at 999 permutations** (`/tmp/isb-timing`, serial, one worker):

| profile | units timed | mean s/unit | median | q90 | inputs (median) | both scans (median) |
| --- | --- | --- | --- | --- | --- | --- |
| `simulation` | 42 (1 per cell) | 0.425 | 0.384 | 0.514 | 0.215 | 0.157 |
| `plasmode` | 6 | 1.799 | 1.801 | 1.853 | 0.032 | 1.743 |
| `plasmode_unequal` | 6 | 1.232 | 1.266 | 1.338 | 0.027 | 1.208 |

**Projection** at full grade (units × mean seconds):

| profile | units | core-hours |
| --- | --- | --- |
| `simulation` | 6,600 | 0.78 |
| `plasmode` | 1,500 | 0.75 |
| `plasmode_unequal` | 1,500 | 0.51 |
| total | 9,600 | 2.0 |

On `ing` the parent run timed the joint scan at a median 0.067 s per simulation unit (`results/interaction-scan-2026-11-09/simulation/accounting.json`). Here it takes about 0.079 s per arm, so `ing` runs the scan at about workstation speed.

- The plasmode profiles should take about 1.3 core-hours on `ing`, or under 3 at a 2× margin. The parent run took 47 core-hours.
- The simulation profile takes about 6 minutes of wall-clock time on the workstation with 8 workers.
