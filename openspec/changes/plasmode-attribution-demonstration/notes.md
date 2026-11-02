# Notes

## Baseline (task 1.1, 2026-11-02)

- Starting commit: `b7df7c7` on `main` (Merge pull request #73, `feat/covariate-adjusted-attribution`).
- `src/motco/__init__.py`: `__version__ = "0.14.0"`.
- Pre-commit gate on the clean tree: `ruff check src/ tests/` and `mypy src/motco/` clean; `MOTCO_TEST_PERMS=99 pytest tests/ -m "not slow"`: 1831 passed, 1 skipped, 32 deselected (12 min 19 s). The local run needed the zlib `LD_LIBRARY_PATH` workaround and BLAS pinned to one thread.

## API surface (task 1.2)

A scratch script (session scratchpad, not committed) on an 80-donor synthetic dataset with the astrocyte cell table, two 30-feature blocks and three covariates confirmed that 0.14.0 needs no library change (design D7):

- `analyze_trajectories` with covariates and attribution returns `attribution` (`mean_source` `covariate_adjusted`, the covariate columns named) and the `df` correction record;
- `latent_artifacts` holds the fitted `PLSRegression`, the standardized joint features (`block__feature` columns) and the original scales;
- `analyze_orientation_attribution` on `latent_artifacts` without covariates gives `mean_source` `arithmetic`; both expose the per-feature `directional_contrast` of every transition and of the `principal` block, and per-feature `sign_stability` and `top_k_selection_frequency` in `bootstrap_records`.

## Fidelity gate on the real export (task 8.1, 2026-11-02)

`fidelity --export-dir ../sead_sexdiff/results/motco_input/cohort_30` (session scratchpad, `--n-jobs 4`, 159 s): **passed**, rank 5, mean held-out stage AUROC 0.6355555555555557 (0.635556 to six decimals), 5 CV repeats completed. Export sha256: `atac.csv` `32ebd055…b493`, `covariates.csv` `7ed4b1bf…c4a3`, `manifest.json` `18f2289e…0b6`, `metadata.csv` `e2d8a412…1fb6`, `rna.csv` `186198b9…a2`; `sead_sexdiff` at `2d427b0` (the export was written from `ab3423a`, dirty, per its manifest).

## Smoke (task 8.2, 2026-11-02; indicative)

`run --e-star 1.0 --reps 20 --perms 99 --bootstrap 50 --workers 7` (scratchpad, 3,244 s wall; ~545 s per unit with 7 concurrent workers), then `merge` and `decide`: 40/40 units, no consistency-check failure, 0 label redraws, every anchor unit at rank 5 and AUROC 0.635556 (the latent check passed on all 20). The decision is marked indicative (R = 20, 99 permutations, B = 50, explicit `e*`). Readings, for orientation only:

- anchor rates `delta` 4/20, `angle` 3/20, `shape` 0/20 (R = 20 bound 0.147);
- at e = 1.0: `angle` 2/20, `delta` 2/20, `shape` 1/20; rank 5 in 14 units, 2 in 4, 7 in 2; unconditional planted-set precision 0.018 (anchor chance reference 0.0025), sign stability over K 0.79, top-k frequency over K 0.008; sign agreement 1.0 on the 6 recovered planted features. e = 1.0 is far below detection, so calibration is expected to climb.

## Timing and compute plan (task 8.3, 2026-11-02)

Three units side by side on the idle workstation (BLAS pinned, `scratchpad/time_units.py`): a paper-grade anchor unit (999 permutations, B = 200, both arms) 331 s, a paper-grade e = 2.0 unit 337 s, a calibration unit (99 permutations, attribution off) 325 s. The stage-PLS double CV on 10,000 features is ~97 % of a unit; RRPP and both attributions add ~10 s. Under 7 concurrent workers a unit takes ~545 s (1.65×, the smoke).

Extrapolation (EPYC assumed ~1.3× slower per core than the workstation, to be checked by a one-unit timing job):

- calibration, ≤ 200 units (5 steps × 40): ~18 core-hours local-equivalent; locally ~52 min per step with 7 workers (≤ 4.3 h); on `ing` as one 40-CPU job, ~8 min per step (≤ 40 min, ≤ 24 core-hours);
- paper grade, 1,500 units: ~140 core-hours local-equivalent, ~180 on EPYC (~215 with the usual 20 % margin); locally ~32 h wall with 7 workers; on `ing` ~150 shards × 10 units (~1.2 h each).

By design D11 (beyond a few hours locally) the paper grade goes to `ing`. The proposal's estimate (~30 s per unit, ~12.5 core-hours) was an order of magnitude low: it did not account for the CV on 10,000 features.

Early signal (not a calibration reading): at replicate 0 of the paper stream, the anchor and e = 2.0 units gave the same rank (5) and nearly the same p-values (`angle` 0.415 in both, `delta` 0.036 / 0.037); the smoke's e = 1.0 cell rejected `angle` in 2/20. A 20-feature ramp moves the 5-dimensional stage-PLS space little, so calibration may reach the cap.

**Approval (2026-11-02).** The user approved calibration on `ing` as one 40-CPU job, asked to see its result before the paper grade is approved, and chose to commit on `feat/plasmode-attribution-demonstration`, push, and run from a detached worktree at that commit.

## Calibration on `ing` (task 9.1, 2026-11-02)

- Branch `feat/plasmode-attribution-demonstration` at `f0b05f0` (rebased onto `main` `b7df7c7`: it had first been cut from another session's HEAD on `feat/simulation-presentation-figures`; that session's commits are not on this branch), pushed to origin.
- `ing`: detached worktree `/home1/tgonza/MOTCO-plasmode` at `f0b05f0`, `.venv` → `/home1/tgonza/MOTCO/.venv` (scikit-learn 1.8.0, numpy 2.3.5, as locally), `PYTHONPATH` set to the worktree's `src` (`motco` 0.14.0 imported from it).
- Export copied to `/home1/tgonza/plasmode-data/cohort_30/` (the five files read); its sha256 equal the local export's and `fidelity.json`'s `export_files`.
- Job 932595 (`scripts/plasmode_calibrate.sbatch`, 40 CPUs, n1, `--exclude=n11`): the fidelity gate again on EPYC, then the ladder; outputs in `/home1/tgonza/plasmode-2026-11-02/calibration/`.

**Result.** Job 932595 `COMPLETED 0:0` in 18 min 49 s (n1, 40 CPUs, TotalCPU 8 h 26 min, MaxRSS 6.7 GB; ~150 CPU-s per unit on EPYC, faster than the workstation). The fidelity gate passed again on EPYC (rank 5, AUROC 0.635556). Calibration **reached the cap**: `angle` rejected in 0/40 units at every step, so `e*` = 3.0 (`cap_reached`, iteration 4, entry sha256 `c4fe69db…5839`) and the descriptive cell is 2.0. Per step (`delta` / `angle` / `shape` rejections of 40; selected ranks):

| e | delta | angle | shape | ranks |
|---|---:|---:|---:|---|
| 0.5 | 2 | 0 | 1 | 5 ×35, 2 ×4, 10 ×1 |
| 1.0 | 1 | 0 | 1 | 5 ×25, 2 ×11, 7 ×2, 10, 12 |
| 1.5 | 2 | 0 | 0 | 5 ×22, 2 ×12, 7 ×4, 8, 11 |
| 2.0 | 1 | 0 | 0 | 5 ×24, 2 ×9, 7 ×4, 8 ×2, 11 |
| 3.0 | 1 | 0 | 0 | 5 ×24, 2 ×8, 7 ×3, 8 ×4, 11 |

`run --calibration` settings name iteration 4 and the cells [0, 3.0, 2.0]. A label-free diagnostic (scratchpad `latent_capture.py`, forced rank 5 on the export) finds the pooled stage displacement 99.9 % inside the rank-5 stage-PLS span, but a random 20-feature plant only 2.2 % (median; q95 3.2 %), the share of any random direction (√(5/10,000)): the planted orientation is outside the measured subspace unless it reshapes the PLS fit, which the drift of the rank from 5 toward 2 suggests it starts to do. Revised paper-grade estimate: 1,500 × ~150–200 CPU-s ≈ 65–85 EPYC core-hours.
