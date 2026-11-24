# Notes

## Implementation (2026-11-24)

- The driver is `scripts/stage_effect_scaling_check.py`, with `scripts/stage_effect_scaling_array.sbatch` and `tests/test_stage_effect_scaling_check.py` (58 tests). The unit counts are 16,000 for `saturated` and 8,400 for `band`. There are 53 pair-cells (4 × 11 + 9) and 71 map cells (44 + 27).
- **The unit gets a parent-shaped settings block.** `univariate_baseline_check.evaluate_unit` dispatches on `settings["profile"]`. So each profile's settings carry `unit_settings`: the dict `univariate_baseline_check.stage_span_settings` writes, with only `units.designs`, `units.cells` and `designs` replaced. The parent's `profile_units`, `rate_rows`, `paired_rows`, `cell_rows`, `accounting`, `characterization_grade` and `label_cells` read that block unchanged. Restricted to `n200_delta0.5`, it equals the parent's own `stage_span_settings(designs=["n200_delta0.5"])` exactly (a test). A rerun record therefore equals the parent driver's record for the same unit, apart from its `settings_hash` and runtime (a test at fast CV knobs).
- **The calibration unit is the driver's own (amended before any unit).** `latent_space_comparison.calibration_unit` goes through `arm_params`, which fixes `layers` to the comparison's methylation and expression. That is two blocks, and the band measures three. The driver's `calibration_unit` keeps the comparison's seed stream (`calibration_seed`, at base seed 2043), target (0.636), band (0.60–0.67) and step summary (`summarize_step`, with the rank's 10 % and 90 % quantiles added). It fits with the units' own production parameters (`univariate_baseline_check.stage_span_params` at zero permutations) on the anchor cell's generator at n = 320. The tasks, design (D4) and proposal record the amendment. The ladder has seven values, so there are 280 calibration fits, not 240.
- **Calibration rule.** δ_c is the value whose median AUROC is nearest 0.636. A tie (equal distances to 12 decimals) goes to the smaller δ. `calibrate` always writes `calibration.json` with the ladder and `adopted.in_band`. `band_settings` refuses an adopted value outside the band and names every ladder reading. `calibrate` refuses to overwrite an existing `calibration.json`, because the band settings bind its sha256.
- **Names and seeds.** The `saturated` designs are `n200_delta0.5` (the rerun, under the parent's name and seeds), `n320_delta1`, `n1280_delta0.5`, `n800_delta1` and `n800_delta0.25`. The `band` designs are `band_n320` (δ_c), `band_n1280_half` (δ_c / 2) and `band_n1280` (δ_c). All are at base seed 2043. The design name enters every seed, so no new dataset seed equals a parent's, the other profile's or a calibration seed (tests).
- **Binding.** Both profiles bind the parent's `settings.json`, `rate_table.csv` and `cell_table.csv` by sha256 (`parent_files`) and its settings hash. `band` also binds `calibration/calibration.json`, resolved under `--out-dir`. On `ing` the local calibration directory is copied into the run root, so the committed calibration is the one the band was built from. `decide` also refuses two profiles bound to different parent files.
- **Merge adds `rank_table.csv`**: per design, cell and selected rank, the available and rejected counts of `production`, the three statistics, `fdr` and `total`. D9's rank strata need per-unit ranks, and the records stay on the cluster.
- **The map uses the parent's `label_cells` unchanged.** It reads its two profiles under the names `stage_span` and `surgery`. The driver passes the four new `saturated` designs under the first name and the three `band` designs under the second, then renames the profiles in its output. The levels are those of m = 71.
- **Size families exclude the rerun design** (`saturated`: the four new designs, m = 8; `band`: m = 6). A per-feature arm is admissible when it holds in both families, as in the parent.
- **An unavailable pair-cell** is one where either design has no available unit for the arm. It is `unresolved` and listed under `label_cells.unavailable`.
- **The committed designs' grades**, read through this driver from the committed rate table (the end-to-end test prints them): `names_kind` at `n80`, `n200` and `n500`, and `detects_only` at `n200_delta0.5`. The parent pooled the four designs into one family and read `detects_only`; per design, the leak sits in the weak column alone.

## Calibration (task 4.1, 2026-11-24)

`calibrate` ran on the workstation, 7 workers with BLAS pinned in the shell, from the uncommitted working tree on `feat/stage-effect-sample-size-scaling` (base `f3a5f0b`; driver sha256 `fa1c33b0091e…`). It made 280 fits (7 ladder values × 40), with 0 failures, 1.77 recorded core-hours and 15 min 17 s of wall time. `calibration.json` sha256 is `cdb7a1f306ed1403c493a44e078b3716f07ad93969f924bd7586fa09ef90d22f`; this file is the one copied to `ing` and committed.

| δ | median held-out AUROC | AUROC q25–q75 | rank q10 / median / q90 |
|---|---|---|---|
| 0.03 | 0.549 | 0.531–0.570 | 2 / 3 / 13 |
| 0.04 | 0.600 | 0.587–0.620 | 2.9 / 8.5 / 19 |
| **0.05** | **0.664** | 0.636–0.681 | 5 / 9 / 19 |
| 0.06 | 0.728 | 0.702–0.745 | 8 / 15.5 / 19 |
| 0.07 | 0.789 | 0.764–0.808 | 8 / 11.5 / 19 |
| 0.08 | 0.847 | 0.821–0.859 | 9.9 / 12.5 / 19 |
| 0.1 | 0.930 | 0.910–0.937 | 9 / 13 / 19 |

**δ_c = 0.05** (median 0.664, 0.028 from the target 0.636, inside 0.60–0.67). No amendment: a value landed in the band. δ_c sits near the band's upper edge. δ = 0.04 reads 0.600, on the lower edge and 0.036 from the target, so the predeclared nearest-target rule takes 0.05. The `band` designs are therefore δ = 0.05 at n = 320, δ = 0.025 at n = 1,280 and δ = 0.05 at n = 1,280 (δ√n = 0.894 for the first two). At δ_c the CV rank already inflates (median 9, q90 at the cap of 19). The `band` settings build from this file (the smoke below ran on them).

## Local smoke (task 4.2, 2026-11-24)

Both profiles ran at 2 replicates per cell (2 control replicates) at 999 permutations: every design, `pca_matched` recorded. They used 7 workers on the workstation (Ryzen 7 5800X, 8 cores / 16 threads, BLAS pinned in the shell), the same working tree as the calibration, and that calibration's `band` settings. That is 196 units with 0 failed; `merge` (both profiles), the reproduction check and `decide` ran. Wall time was 835 s for `saturated` and 803 s for `band`.

- **No planted fixed point failed**, in any design or cell, the band's e = 2.0 stretch and e = 1.0 rotation and bend at δ = 0.025 included. `failures_by_cell` and `failures_by_error` are empty in both accountings. No design needs amending.
- **The reproduction differs only by size.** It compared 234 rows, with 0 unresolved. The 182 that differ all count 2 units against the committed 500 or 200. The other 52 are the `concat` arms, with 0 units on both sides. So `decide` reads `reproduction_failed`; all 53 pair-cells are `unresolved` and the grades rest on 2 units per cell. None of it is a finding.

Seconds per unit from the records, with the fit and AUROC at 2 replicates per cell:

| Profile | Design | units | unit mean | fit median | CV rank median (q10–q90) | held-out AUROC median |
|---|---|---|---|---|---|---|
| `saturated` | `n200_delta0.5` | 26 | 18.3 | 16.3 | 3 (3–3) | 1.000 |
| `saturated` | `n320_delta1` | 26 | 23.2 | 21.0 | 2 (2–3) | 1.000 |
| `saturated` | `n1280_delta0.5` | 26 | 94.9 | 81.8 | 3 (3–3) | 1.000 |
| `saturated` | `n800_delta1` | 26 | 41.0 | 35.6 | 2 (2–2) | 1.000 |
| `saturated` | `n800_delta0.25` | 26 | 43.8 | 38.4 | 5 (5–6) | 1.000 |
| `band` | `band_n320` (δ = 0.05) | 22 | 25.9 | 22.8 | 9 (2–10) | 0.649 |
| `band` | `band_n1280_half` (δ = 0.025) | 22 | 110.5 | 94.4 | 19 (19–19) | 0.644 |
| `band` | `band_n1280` (δ = 0.05) | 22 | 104.6 | 92.4 | 19 (14–19) | 0.878 |

The n = 1,280 units cost more than the probe's 75 s. In `band` the fit selects rank 19, and the split test's median rises to 9.4 s (2.5 s in `saturated`). The smoke matches the probe at n = 320 and 800.

Signal cells where an arm rejected in both units / in at least one (2 replicates; a smoke reading only):

| Design | cells | `production` | `delta` | `angle` | `shape` | `fdr` | `total` | `scan` | `pca_matched` |
|---|---|---|---|---|---|---|---|---|---|
| `n200_delta0.5` | 11 | 1 / 6 | 0 / 3 | 1 / 2 | 0 / 3 | 4 / 7 | 1 / 2 | 4 / 8 | 2 / 3 |
| `n320_delta1` | 11 | 7 / 9 | 3 / 3 | 2 / 3 | 2 / 4 | 7 / 8 | 4 / 6 | 7 / 8 | 3 / 8 |
| `n1280_delta0.5` | 11 | 5 / 10 | 2 / 10 | 2 / 3 | 2 / 7 | 8 / 8 | 5 / 6 | 8 / 8 | 5 / 8 |
| `n800_delta1` | 11 | 11 / 11 | 4 / 5 | 3 / 3 | 5 / 5 | 9 / 10 | 6 / 6 | 9 / 10 | 7 / 8 |
| `n800_delta0.25` | 11 | 4 / 6 | 1 / 2 | 2 / 10 | 4 / 8 | 4 / 4 | 2 / 10 | 4 / 4 | 0 / 0 |
| `band_n320` | 9 | 1 / 1 | 1 / 1 | 0 / 0 | 0 / 0 | 0 / 0 | 0 / 0 | 0 / 1 | 0 / 0 |
| `band_n1280_half` | 9 | 2 / 5 | 2 / 9 | 0 / 0 | 0 / 3 | 1 / 1 | 0 / 0 | 1 / 1 | 0 / 0 |
| `band_n1280` | 9 | 4 / 6 | 3 / 4 | 2 / 3 | 1 / 4 | 2 / 2 | 0 / 1 | 2 / 2 | 0 / 0 |

The controls (anchor and translation) rejected: `production` 0 of 20 in `saturated` and 1 of 12 in `band`; `fdr` 0 of both; `total` 2 of 20 in `saturated`.

### Projection

Units × the smoke's mean seconds per unit, by design:

| Profile | units | core-hours |
|---|---|---|
| `saturated` (`n200_delta0.5` 16.3, `n320_delta1` 20.6, `n1280_delta0.5` 84.3, `n800_delta1` 36.5, `n800_delta0.25` 39.0) | 16,000 | 197 |
| `band` (`band_n320` 20.1, `band_n1280_half` 86.0, `band_n1280` 81.3) | 8,400 | 187 |
| total | 24,400 | **384** |

On `ing`, the univariate baseline's `stage_span` array recorded a mean of 24.9 s per unit (88.5 core-hours over 12,800 units, across n = 80 to 500). This smoke ran its `n200_delta0.5` at 18.3 s. So `ing` ran about 1.25–1.4 times slower than this workstation at 7 workers, and the run is projected at **480–540 recorded core-hours**. That is above the proposal's pre-smoke estimate of 320–450, because of the n = 1,280 units. The compute levers are unchanged (design, risk list), in order: the two `band` designs at n = 1,280 (167 core-hours here), the `n800_delta1` pair (36.5), and the signal replicates.

### Proposed cap and launch plan (for the user's approval)

- **Cap: 600 recorded core-hours**, above the 540 upper projection. No amendment: no design failed, and the calibration landed in the band.
- **Launch:** commit and push the branch, then make a clean worktree on `ing` at that revision, with `PYTHONPATH` on its `src` and the shared `.venv`. Run the driver's tests there. Copy this calibration directory (sha256 `cdb7a1f3…`) into the run root, and pre-create both profiles' settings. Then submit `saturated` (200 shards, ~80 units each) and `band` (140 shards, ~60 units each), throttled at `%64`.

**Approved** by the user on 2026-11-24: a cap of 600 recorded core-hours, both profiles launched on `ing` as planned, with no amendment.
