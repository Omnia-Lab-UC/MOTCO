# Notes

## Baseline and identity (tasks 1–2, 2026-11-04)

- Starting commit: `cb8da5e` on `feat/latent-space-comparison` (the proposal, on top of `main` `880fb00`, PR #76).
- `src/motco/__init__.py` was `0.15.0`; bumped to `0.16.0` (task 6.2).
- Local runs need the zlib `LD_LIBRARY_PATH` workaround and BLAS pinned to one thread.
- **Byte identity, checked in-process before and after the change** (scratch `capture_baseline.py` / `compare.py`, session
  scratchpad, not committed; BLAS pinned): 15 stage-span datasets (n = 120, the SEA-AD table, four stages × five modes; all
  arrays, metadata, truth and population means) and 24 evaluations (three stage-span modes × `pls`, `pca`, `concat`,
  `pls_within_group` × span capture on/off, 19 permutations: statistics, p-values, null summaries, latent metadata,
  config spectrum and realized geometry). After the generator change: 0 of 39 differ. After the linear-map change: the
  only differences are the six `pca`/`concat` evaluations with span capture on, which gain exactly the
  `population_latent` checkpoint and the `span_capture` block (everything else equal; `concat` shares 1.0, rank 498).
  PLS span-capture records are unchanged, and `tests/test_span_capture.py::test_pls_span_capture_records_are_unchanged`
  pins them against a frozen copy of the 0.15.0 computation.

## Design notes from implementation

- **The linear-map record lives beside `artifacts`.** Design D6 says `LatentIntegrationResult.artifacts` gains the
  record. `artifacts` is the PLS attribution boundary that `motco.analysis` and attribution read as
  `PLSIntegrationArtifacts | None`, so the record is a separate ephemeral field, `LatentIntegrationResult.linear_map`
  (a `diagnostics.LinearLatentMap`), present for `pls`, `pca`, `pls_residual_pca` and the standardized `concat`. It is
  never serialized, as D6 requires.
- **`concat` with `standardize=False`** carries no map (its coordinates are not the preprocessor's), so span capture is
  inert there, as for `snf` and `pls_within_group`.
- **`pls_residual_pca`'s residual rank** is bounded by `min(n_features − PLS rank, n_samples − n_stages)`: the
  stage-residualized matrix has rank at most n − K, so singular vectors beyond it would be arbitrary null-space
  directions, not necessarily orthogonal to the PLS span. The error names the range.
- **`n200` is near-balanced, not exactly.** n = 200 at three equal stage proportions gives stage sizes 67/67/66 and
  cells A 34/34/33, B 33/33/33, so the rhombus is exact only to that imbalance. The capture table measures it.
- **The alignment fields in study configs.** Zero-effect cells (`enumerate._without_construction`, the Type I grid)
  clear `span_alignment` and `span_support` with the construction, and dumped configs omit them while unset. The study
  config does not otherwise validate them (non-goal: the comparison runs through its driver).
- **Unused features at `p_dmp` 0.1.** At replicate 0, 347 (`n200`) and 340 (`seaad`) of the 498 span-layer features
  have constant null stage means, so the dense `o` spreads over ~340 features and the 20-feature support is far below
  the available count.

## Gate (task 6.2, 2026-11-04)

`ruff check src/ tests/` and `mypy src/motco/` clean; `MOTCO_TEST_PERMS=99 pytest tests/ -m "not slow"`: 2212 passed,
1 skipped (R `shapes`), 32 deselected (9 min 2 s, BLAS pinned). `mkdocs build` gives the same 8 warnings as before the
change (none new).

## Smoke and timing (task 7.1, 2026-11-04; indicative)

- **Smoke.** `run --weak-delta 0.15 --reps 2 --perms 19 --workers 7` (session scratchpad), then `merge` and `decide`:
  84/84 units in 247 s wall, 0 failed records, 0 duplicates; all 420 records carry `span_capture` and
  `population_latent`. Median unit 18.9 s under 7 workers, of which the shared stage-PLS fit is 15.4 s; per arm at 19
  permutations ≈ 0.3 s, `concat` 1.85 s. `decide` ran end to end (`document_tradeoff` at R = 2, meaningless).
- **What it showed (R = 2, descriptive only).** At δ = 2.0 the held-out stage AUROC is 1.000 and the `pls_cv` rank
  2–3. At δ = 0.15 on the SEA-AD table the median AUROC is 0.594, just under the band, so calibration should adopt a
  step near 0.2. Median `group_difference` capture for orientation e = 0.5 at φ = 0 / 30 / 60 / 90 at `n200`:
  `pls_cv` 0.970 / 0.822 / 0.398 / 0.008, `pca_fixed` 0.971 / 0.855 / 0.446 / 0.023, `pls_residual_pca`
  0.973 / 0.825 / 0.405 / 0.020, `concat` 1.000 throughout. At φ = 90 `pls_cv` sits at the random share (2/498 =
  0.004); `pls_residual_pca` reads 0.165 on the SEA-AD table. `pca_matched` at rank 2 keeps only 0.52 of the pooled
  stage means at `n200`. The sparse 20-feature cell reads `population_latent` `angle` ≈ 34° through `concat` against
  the intended 45° (risk R6 in the design). A candidate mechanism, not yet checked: the sample scaler's pooled SD
  includes the planted between-group spread on those 20 features and shrinks them. The report compares realized and
  intended `angle` for that cell.
- **Timing at 199 permutations.** `run --reps 1 --perms 199 --designs n200 seaad --workers 7`: 28 units in 144 s wall,
  0 failures. Median unit 33.8 s (q90 38.4 s) under 7 workers; per arm `concat` 15.7 s, every other arm 0.38–0.40 s
  (RRPP in 2–5 dimensions is nearly free); the PLS fit is the rest.

**Core-hour projection.** 4,200 units × 33.8 s ≈ 39 core-hours local-equivalent under 7 workers (≈ 5.6 h wall
locally). On `ing` EPYC, from the stage-span pilot's 17.8 s per single-arm unit at 199 permutations (9,200 units,
45.5 core-hours) plus `concat`'s RRPP (≈ 10–15 s), ≈ 30 s per unit: **≈ 35 core-hours, ≈ 42 with the usual 20 %
margin** (the proposal estimated 50). Calibration: 7 steps × 40 fits × ≈ 15 s ≈ 1.2 core-hours (one 40-CPU job,
a few minutes). Dropping `concat` from the weak design (the design's fallback) is not needed at this cost.

## Calibration on `ing` (task 7.2, 2026-11-04)

- Commit `93fb7a1` on `feat/latent-space-comparison`, pushed only to the cluster clone's side ref
  `incoming/latent-space-comparison` (no GitHub push). Detached worktree `/home1/tgonza/MOTCO-lsc` at `93fb7a1`, clean,
  `.venv` → `/home1/tgonza/MOTCO/.venv`, `PYTHONPATH` = the worktree's `src` (motco 0.16.0; numpy 2.3.5,
  scikit-learn 1.8.0, scipy 1.16.3, as locally); the driver's, `evaluation.py`'s and `semisynthetic.py`'s sha256 equal
  the local files'.
- Job 933435 (one 40-CPU task, `--exclude=n2,n11`: n2 still held a job stuck in `completing`, n11 had 15 GB free and an
  NFS-stall history): `COMPLETED 0:0` in 3 min 49 s on n1, TotalCPU 1 h 19 min, MaxRSS 5.3 GB; 280 fits, 0 failures,
  1.26 core-hours recorded. Outputs copied to `results/latent-space-comparison-2026-11-04/calibration-original-ladder/`
  (`calibration.json` sha256 `cd846d6a…ddfc`, settings hash `657de6a7…`).

| δ | median AUROC | q25–q75 | median rank |
|---|---:|---|---:|
| 0.10 | 0.595 | 0.534–0.624 | 2.5 |
| 0.15 | 0.710 | 0.653–0.735 | 4 |
| 0.20 | 0.838 | 0.795–0.867 | 6 |
| 0.25 | 0.922 | 0.901–0.947 | 7 |
| 0.30 | 0.973 | 0.962–0.985 | 7 |
| 0.40 | 0.999 | 0.997–0.999 | 7 |
| 0.50 | 1.000 | 0.999–1.000 | 4.5 |

**No step lands in [0.60, 0.67].** The band lies between the first two steps. The predeclared rule adopts the closest
step, **δ = 0.10, `band_missed`** (0.595, 0.041 below 0.636; δ = 0.15 is 0.074 above). The prior was right about
δ = 0.1 (the rank-floor check read 0.58 there) and wrong that a step near 0.15 would land in the band. SEA-AD's own CV
rank is 5; the adopted step's median rank is 2.5. The smoke's 0.594 at δ = 0.15 came from two null datasets and was
noise.

**Amendment (2026-11-04, user decision, before any run unit).** Rather than adopt δ = 0.10 as `band_missed`, the
ladder gained 0.11, 0.12, 0.13 and 0.14 (design D9), and calibration ran again on the amended ladder under the same
rule. The first run stays in `calibration-original-ladder/`.
