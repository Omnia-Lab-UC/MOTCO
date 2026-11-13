# Notes

## Pre-flight and smoke (2026-11-13, task 5.1)

Run on the workstation from the worktree's uncommitted state at `f4d2ce0` plus the driver (MOTCO 0.24.0; Python
3.11.15, numpy 2.3.5, scipy 1.16.3, scikit-learn 1.8.0, pandas 2.3.3, the versions of the pilot's `ing` run), BLAS
pinned to one thread. The export is `sead_sexdiff/results/motco_input/cohort_30`, read in place. Its five file hashes
equal `results/plasmode-2026-11-02/settings.json`'s, and the settings bind the pilot's plasmode `settings.json` and
`rate_table.csv`. Output went to the session scratchpad and is not committed.

**Pre-flight: passed.** `preflight --reps 2 --plant-reps 2` ran on anchor replicates 0 and 1, serially, with nothing
else running. On both units the cross-validated `pls_cv` selected rank 5 at held-out AUROC 0.635556, the plasmode's
fidelity values. Its statistics and p-values equal the call forced at rank 5 exactly, with no differing key. Replicate
0 gave p = 0.036 / 0.415 / 0.822 for `delta` / `angle` / `shape`. The CV took 132.8 s and 127.3 s. The forced call
took 0.49 s and 0.57 s.

**Smoke: 8 of 8 units completed, none failed, none inconsistent.** `run --reps 2 --plant-reps 2 --workers 7` covered
two anchor units and two units per plant cell at 999 permutations. The fast test suite ran on an eighth core at the
same time. The wall clock was 2 min 37 s.

- **Consistency.** On every unit, every applicable check held. The largest deviations were `path_identity` 0,
  `same_fit` 0, `same_fit_covariate_free` 0, `rank` 0 and `fitted_labels` 4.2e-14. `one_subspace` (no capture on
  the plasmode) and `identity_at_low_rank` (rank 5–7 > K − 1) did not apply. Every covariate-free arm's latent
  matrix equalled its reference arm's exactly.
- **Ranks.** The anchor fit at 5 (forced). The plant units cross-validated to 5–7: `span_e0.25` 7 and 6, `span_e0.5`
  7 and 7, `span_e1` 5 and 5.
- **Stage fit.** In-sample R² was 0.840 on both anchor units, the pilot's value. The out-of-fold share of the
  stage-indicator variance was 0.030 and 0.029. On the plant units R² was 0.85–0.96 and the out-of-fold share
  0.27–0.52.
- **Merge, reproduction and decide** all ran. The reproduction differs on all 48 rows, because the counts are 2
  against the pilot's 500 and 100. `decide` therefore returns `reproduction_failed`, as it must at this size; the
  grades it records carry no information at two units.

**Timing** (per unit, under the smoke's 8-process load):

| quantity | anchor | plant |
|---|---|---|
| unit | 6.9–7.2 s | 147.3–150.1 s (median 148.5) |
| production fit (`pls_cv`) | 0.5–0.6 s (forced) | 140.4–143.3 s (one CV) |
| each other analysis arm (forced; with or without covariates) | 0.47–0.61 s | 0.50–0.65 s |
| each path arm (within-stage or out-of-fold, 999 draws) | 0.26–0.28 s | 0.27–0.29 s |
| path identity (four `RRPP` runs) | 1.1 s | 1.1–1.2 s |
| out-of-fold fits (four folds) | 0.15 s | 0.13–0.24 s |

**Core-hour projection.**

- **Laptop.** 500 anchor units at 7.0 s make 0.97 core-hours. 300 plant units at 148.5 s make 12.4. The total is
  about 13.4 core-hours, or about 1.9 h of wall clock at 7 workers. These per-unit times were measured under that
  parallel load.
- **`ing`.** The pilot's plasmode `pls_cv` took 116.9 s per unit on `ing` (`results/predictive-subspace-2026-11-12/
  PROVENANCE.txt`), against about 130 s serially here. A plant unit there would take about 125 s and an anchor unit
  about 8 s, about 11.5 core-hours in all. A 100-shard array (8 units per shard) drains in about 10 minutes per task.
- **Cap.** 20 recorded core-hours on either host, about 1.5 times the projection.

The design's Open Question is still open: the full run waits on the user's approval, and on their choice of the
laptop or `ing`.
