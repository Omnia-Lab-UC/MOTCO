# Phase 5 design-point pilot

**Run date:** 2026-09-08
**Configuration:** `examples/trajectory_power_study/phase5_design_point_pilot.json`
(sha256 `153f53b7…`, recorded in full in `results/phase5-design-point-2026-09-08/PROVENANCE.txt`)
**Code revision:** `ff68514` on `feat/phase5-design-point-pilot` (clean working tree at launch; parent
`afab7df` = `main`). Record generation and report generation used the same revision.
**Execution:** SLURM array job 878322 on `cluster.ing.uc.cl`, partition `512x1024` (AMD EPYC 7662 nodes),
100 single-CPU shards, `--error-policy record`, no `--n-jobs`, BLAS pinned to one thread.
**Environment:** Python 3.11.16, numpy 2.3.5, pandas 2.3.3, scikit-learn 1.8.0, scipy 1.16.3, motco 0.6.0,
installed from the lockfile with uv 0.12.10; no R at runtime.
**Answers:** readiness item 4 ([`phase5-readiness.md`](https://github.com/Omnia-Lab-UC/MOTCO/blob/c440070dcb3c5bcc6ddd2227cf7adbef18573b6a/docs/phase5-readiness.md)) — the Phase 5 design point and
the 0.80 orientation power floor. Every number below traces to a file under
`results/phase5-design-point-2026-09-08/report/`; the JSONL records are gitignored as regenerable.

## Design-point decision: **CHOSEN — ρ = 0, n = 1200**

> First design point in preference order (`n_samples` ascending, then `baseline_continuity` ascending) whose
> orientation/`angle` power at its top effect (1.00) clears the 0.80 floor with 1·SE confirmation:
> **rate 0.88, lower bound 0.85** (`report/design_point_decision.json`).

The rule was predeclared in the config (`acceptance.design_point`) and is advisory. Its per-column
evaluation (`report/design_point_decision.csv`):

| ρ | n | orientation `angle` power at e = 1.00 | MC SE | 1·SE lower bound | status | anchor Type I (delta / angle / shape) |
|---:|---:|---:|---:|---:|---|---|
| 0.0 | 300 | 0.55 | 0.050 | 0.50 | fails | 0.03 / 0.06 / 0.03 |
| 0.5 | 300 | 0.24 | 0.043 | 0.20 | fails | 0.01 / 0.02 / 0.01 |
| 0.8 | 300 | 0.69 | 0.046 | 0.64 | fails | 0.04 / 0.01 / 0.01 |
| 0.0 | 600 | 0.76 | 0.043 | 0.72 | fails | 0.04 / 0.05 / 0.05 |
| 0.5 | 600 | 0.46 | 0.050 | 0.41 | fails | 0.02 / 0.00 / 0.02 |
| 0.8 | 600 | 0.73 | 0.044 | 0.69 | fails | 0.02 / 0.05 / 0.07 |
| **0.0** | **1200** | **0.88** | 0.032 | **0.85** | **meets** | 0.04 / 0.04 / 0.03 |
| 0.5 | 1200 | 0.69 | 0.046 | 0.64 | fails | 0.05 / 0.02 / 0.02 |
| 0.8 | 1200 | 0.81 | 0.039 | 0.77 | marginal | 0.02 / 0.04 / 0.03 |

The chosen column's own zero-effect anchor holds nominal level on all three statistics (0.04 / 0.04 / 0.03
against α = 0.05), so the design.md risk — a column that meets the floor on `angle` while another statistic
drifts — did not materialize. ρ = 0 is also the declared isotropic stress-test endpoint and the Phase 4
comparator, so the chosen point carries the **most general** claim the grid offered: the floor is met
without assuming any baseline continuity.

The figure `report/design_point_power.png` shows the same nine points as curves: orientation `angle`
rejection rate at e = 1.00 versus `n_samples`, one line per ρ, Monte Carlo error bars, each point annotated
with the column's median pooled eigengap. The ρ = 0 line rises steepest (0.55 → 0.76 → 0.88) and is the only
one to clear 0.80 with margin; the ρ = 0.5 line lies below it everywhere; the ρ = 0.8 line starts highest at
n = 300 and is overtaken by n = 600.

## Run configuration

| Parameter | Value |
|---|---:|
| Integration | Pooled PLS, M-value methylation |
| PLS cross-validation | `cv1_splits=3`, `cv2_splits=4`, `n_repeats=5`, `max_components=20`, `random_state=1203` |
| Stages | 4 |
| Differential-site density `p_dmp` | **0.1** (Phase 4 used 0.2) |
| Surgery censoring | default `"error"` (fail loud); no cell was censored |
| Design grid | ρ ∈ {0.0, 0.5, 0.8} × `n_samples` ∈ {300, 600, 1200} (9 columns; ρ = 0, n = 300 is the baseline column) |
| Trajectory modes | orientation, translation |
| Effect sizes | 0.00 (shared anchor per column), 0.25, 0.50, 1.00 |
| Replicates per cell | 100 |
| Permutations per test | 199 |
| Seed policy | one matched primary family across the whole grid; one zero-effect anchor per column |
| Attribution / Phase 4 gate | off / off |
| Cells | 65 (2 Type I + 63 primary/design) |
| Total work units | 6,500 |
| Failed work units | 0 |
| Censored surgeries | 0 (`report/realized_surgery.csv`, `censored_fraction` = 0.0 in all 55 surgery cells) |
| Compute | 70.4 core-hours of recorded unit runtime (75.2 CPU-hours billed by `sacct`), 70 minutes wall on 100 concurrent tasks |

All 6,500 units completed with status `completed`; every (cell, replicate) appears exactly once in
`merged.jsonl`; every record's parameter signature equals the signature enumerated from the committed
config; every record carries 199 permutations, `n_jobs = 1`, `config_spectrum`, and `selected_lv`.

## Operating characteristics per design point

Rejection rates at α = 0.05, 100 replicates per entry (`report/design_point_operating.csv`). Bold is
orientation's target statistic. Each column's `none` row is that column's shared zero-effect anchor — one
measurement per column, not one per mode.

| ρ | n | mode | e | delta | **angle** | shape |
|---:|---:|---|---:|---:|---:|---:|
| 0.0 | 300 | none | 0.00 | 0.03 | 0.06 | 0.03 |
| 0.0 | 300 | orientation | 0.25 | 0.60 | **0.45** | 0.82 |
| 0.0 | 300 | orientation | 0.50 | 0.67 | **0.63** | 0.93 |
| 0.0 | 300 | orientation | 1.00 | 0.69 | **0.55** | 0.91 |
| 0.0 | 300 | translation | 0.25–1.00 | 0.02–0.05 | 0.02–0.06 | 0.01–0.03 |
| 0.5 | 300 | none | 0.00 | 0.01 | 0.02 | 0.01 |
| 0.5 | 300 | orientation | 0.25 | 0.26 | **0.16** | 0.41 |
| 0.5 | 300 | orientation | 0.50 | 0.33 | **0.19** | 0.57 |
| 0.5 | 300 | orientation | 1.00 | 0.49 | **0.24** | 0.76 |
| 0.5 | 300 | translation | 0.25–1.00 | 0.04–0.06 | 0.01–0.03 | 0.01–0.02 |
| 0.8 | 300 | none | 0.00 | 0.04 | 0.01 | 0.01 |
| 0.8 | 300 | orientation | 0.25 | 0.16 | **0.07** | 0.26 |
| 0.8 | 300 | orientation | 0.50 | 0.21 | **0.27** | 0.32 |
| 0.8 | 300 | orientation | 1.00 | 0.27 | **0.69** | 0.43 |
| 0.8 | 300 | translation | 0.25–1.00 | 0.03–0.05 | 0.03–0.06 | 0.01–0.05 |
| 0.0 | 600 | none | 0.00 | 0.04 | 0.05 | 0.05 |
| 0.0 | 600 | orientation | 0.25 | 0.77 | **0.74** | 0.96 |
| 0.0 | 600 | orientation | 0.50 | 0.81 | **0.81** | 1.00 |
| 0.0 | 600 | orientation | 1.00 | 0.78 | **0.76** | 0.99 |
| 0.0 | 600 | translation | 0.25–1.00 | 0.02–0.05 | 0.00–0.06 | 0.02–0.06 |
| 0.5 | 600 | none | 0.00 | 0.02 | 0.00 | 0.02 |
| 0.5 | 600 | orientation | 0.25 | 0.36 | **0.38** | 0.63 |
| 0.5 | 600 | orientation | 0.50 | 0.52 | **0.39** | 0.83 |
| 0.5 | 600 | orientation | 1.00 | 0.44 | **0.46** | 0.87 |
| 0.5 | 600 | translation | 0.25–1.00 | 0.03–0.07 | 0.00–0.04 | 0.01–0.08 |
| 0.8 | 600 | none | 0.00 | 0.02 | 0.05 | 0.07 |
| 0.8 | 600 | orientation | 0.25 | 0.32 | **0.22** | 0.47 |
| 0.8 | 600 | orientation | 0.50 | 0.38 | **0.32** | 0.53 |
| 0.8 | 600 | orientation | 1.00 | 0.31 | **0.73** | 0.61 |
| 0.8 | 600 | translation | 0.25–1.00 | 0.03–0.05 | 0.01–0.04 | 0.03–0.05 |
| 0.0 | 1200 | none | 0.00 | 0.04 | 0.04 | 0.03 |
| 0.0 | 1200 | orientation | 0.25 | 0.88 | **0.83** | 1.00 |
| 0.0 | 1200 | orientation | 0.50 | 0.86 | **0.89** | 1.00 |
| 0.0 | 1200 | orientation | 1.00 | 0.86 | **0.88** | 1.00 |
| 0.0 | 1200 | translation | 0.25–1.00 | 0.04–0.09 | 0.03–0.04 | 0.03–0.06 |
| 0.5 | 1200 | none | 0.00 | 0.05 | 0.02 | 0.02 |
| 0.5 | 1200 | orientation | 0.25 | 0.65 | **0.71** | 0.89 |
| 0.5 | 1200 | orientation | 0.50 | 0.65 | **0.70** | 0.92 |
| 0.5 | 1200 | orientation | 1.00 | 0.60 | **0.69** | 0.97 |
| 0.5 | 1200 | translation | 0.25–1.00 | 0.05–0.07 | 0.01–0.02 | 0.03–0.06 |
| 0.8 | 1200 | none | 0.00 | 0.02 | 0.04 | 0.03 |
| 0.8 | 1200 | orientation | 0.25 | 0.30 | **0.26** | 0.48 |
| 0.8 | 1200 | orientation | 0.50 | 0.51 | **0.48** | 0.66 |
| 0.8 | 1200 | orientation | 1.00 | 0.46 | **0.81** | 0.69 |
| 0.8 | 1200 | translation | 0.25–1.00 | 0.04–0.05 | 0.02–0.04 | 0.03–0.05 |

Translation rows are collapsed to their range over the three nonzero effects; the individual rows are in
the CSV.

### Type I behaviour, per column

Every column's anchor and every translation cell sits at or below the one-sided bound
`0.05 + 2·sqrt(0.05·0.95/100) = 0.0936` on all three statistics; the largest control rate anywhere in the
grid is 0.09 (translation `delta`, ρ = 0, n = 1200, e = 0.25). Type I calibration does **not** drift with
sample size or continuity, which is the property the design.md risk list asked this report to check before
adopting a column.

The baseline `type_i` acceptance targets (`report/acceptance_report.json`, two-sided `|rate − α| ≤ 2·SE`)
pass on 5 of 6 tests. The one flagged test is the translation control's `delta` rate of **0.02** — 0.03
*below* α against a 0.028 bound — i.e. a conservative rate, not an inflation. The Phase 4 report judged
these controls one-sided (inflation only); under that reading all six pass. It is reported here as flagged
because the config's target is two-sided.

### Orientation at the baseline column, against Phase 4

At ρ = 0, n = 300 the orientation `angle` curve is 0.45, 0.63, 0.55 across e = 0.25 / 0.50 / 1.00
(`report/power_curves.csv`). Phase 4 measured 0.59, 0.59, 0.64, 0.65 (0.64–0.68 after the sign-anchor
correction) at `p_dmp = 0.2`. The *shape* of the curve is the same — no monotone response to the requested
effect; adjacent points differ by at most 2.6 combined standard errors and the top effect is not the most
powerful — and the level is about 0.1 lower, consistent with halving the differential-site density rather
than with any change in the estimator. The off-diagonal pattern is also unchanged: orientation moves `shape`
(0.82–0.93 here, 0.97–1.00 in Phase 4) and `delta` (0.60–0.69 here, 0.60–0.69 in Phase 4). Item 2
predeclared the `shape` response as a projection artifact of the rank-3 stage-supervised space; nothing here
contradicts that, and at n ≥ 600 the response saturates at 0.96–1.00 in every ρ = 0 and ρ = 0.5 column.

## What governs orientation power across the grid

The design-point table carries the two recorded covariates readiness item 4 asked for — the pooled
configuration eigengap and the dispersion of each replicate's own `angle` null width (`q95`). Orientation
cells at e = 1.00, from `report/continuity_resolved_orientation.csv` (identical rows appear in
`design_point_operating.csv`):

| ρ | n | `angle` power | median eigengap (q33–q67) | median `angle` null q95 | IQR of q95 | SD of q95 |
|---:|---:|---:|---|---:|---:|---:|
| 0.0 | 300 | 0.55 | 0.052 (0.035–0.067) | 40.7° | 107.4° | 59.1° |
| 0.0 | 600 | 0.76 | 0.048 (0.038–0.066) | 18.8° | 38.2° | 55.5° |
| 0.0 | 1200 | 0.88 | 0.049 (0.038–0.063) | 10.8° | 13.0° | 48.4° |
| 0.5 | 300 | 0.24 | 0.086 (0.071–0.125) | 56.6° | 116.6° | 57.4° |
| 0.5 | 600 | 0.46 | 0.114 (0.090–0.149) | 26.2° | 42.2° | 54.2° |
| 0.5 | 1200 | 0.69 | 0.133 (0.094–0.188) | 16.3° | 12.1° | 29.5° |
| 0.8 | 300 | 0.69 | 0.221 (0.145–0.285) | 24.6° | 83.3° | 60.2° |
| 0.8 | 600 | 0.73 | 0.193 (0.156–0.269) | 21.3° | 37.6° | 54.1° |
| 0.8 | 1200 | 0.81 | 0.186 (0.142–0.234) | 15.3° | 20.0° | 49.9° |

**Sample size does what item 1 predicted.** Along n at fixed ρ = 0 the eigengap is constant (≈ 0.05: the
baseline is near-isotropic at every n, as it must be — ρ, not n, sets the configuration's shape) while the
median `angle` null width contracts 40.7° → 18.8° → 10.8° and its IQR collapses 107° → 38° → 13°. Power
rises 0.55 → 0.76 → 0.88 with it. The zero-effect anchors show the same contraction (median q95 17.3°, 8.1°,
4.9° at n = 300, 600, 1200; `design_point_operating.csv`, `none` rows). This is the mechanism the pivotality
report identified — each replicate's critical value tracks its own geometry, and more samples per
group-stage cell pin that geometry down — now measured across n rather than inferred at one point.

**Continuity does not help at this baseline, and its effect is non-monotone.** ρ = 0.5 is *worse* than
ρ = 0 at every n (0.24 / 0.46 / 0.69 vs 0.55 / 0.76 / 0.88), and ρ = 0.8 is above ρ = 0 only at n = 300
(0.69 vs 0.55), then below it (0.73 vs 0.76, 0.81 vs 0.88). This is not what the headroom analysis
anticipated, and the covariates say why the simple story fails:

- The eigengap behaves as designed — it rises from ≈ 0.05 at ρ = 0 to 0.09–0.13 at ρ = 0.5 and 0.19–0.22 at
  ρ = 0.8, so the stage configuration does acquire a dominant PC1 — and at ρ = 0.5, n = 1200 the anchor's
  null is the **narrowest in the grid** (median q95 2.2° against 4.9° at ρ = 0). The null side of the test
  improves with continuity exactly as the audit predicted.
- What does not carry over is the **signal**. Computed from the records (`merged.jsonl`,
  `pair_statistics.angle`; reproduction snippet below), the median observed latent angle for orientation at
  e = 1.00 is 61° / 73° / 66° at ρ = 0 for n = 300 / 600 / 1200, but only **46° / 35° / 28°** at ρ = 0.5 —
  the same nominal surgery realizes a smaller latent orientation contrast when the baseline trends, and the
  contrast *shrinks* as n grows. At ρ = 0.8 the observed angle is back to 64° / 62° / 66°, but the anchor's
  null is wide again (median q95 54° / 26° / 19°, versus 23° / 9° / 2° at ρ = 0.5), so the larger eigengap
  buys nothing on the null side there.

So "larger eigengap ⇒ narrower null ⇒ more power" holds only along n at fixed ρ. Across ρ, two things
change at once — the realized orientation contrast of the surgery and the null width — and they do not move
together. **The claim Phase 5 can make is therefore n-conditional, not continuity-conditional:** at the
chosen point the floor is met at the isotropic endpoint, with an eigengap distribution (median 0.049,
terciles 0.038–0.063) that is the real-data observable to report alongside any orientation result. A
trending baseline is *not* a lever Phase 5 should pull to buy orientation power; if anything, the ρ = 0.5
columns are a warning that the orientation surgery's semantics are not ρ-invariant, which belongs in the
Phase 5 report contract as a stated limitation of the construction.

## PLS representation

Selected latent dimensionality (`design_point_operating.csv`, `median/min/max_selected_lv`):

| column | anchors (median, range) | orientation e = 1.00 (median, range) |
|---|---|---|
| ρ = 0, every n | 3 (2–4) | 3–4 (2–4) |
| ρ = 0.5, every n | 3–4 (2–4) | 4 (3–6) |
| ρ = 0.8, n = 300 / 600 / 1200 | 4 (4–7), 4 (4–5), 4 (3–5) | **8 (4–14), 7 (4–10), 7 (4–9)** |

At the chosen point selection saturates at `n_stages − 1 = 3` (median 3, range 2–4, orientation cells
included), exactly as in Phase 4: the stage-supervised double CV sizes the space for stage separation and
does not respond to n. The retained-rank question — whether a space sized for stage centroids should be the
space that measures the *group* orientation contrast — is **handed to readiness item 3 unchanged**; this
pilot recorded rank and did not vary it.

One new observation for item 3: at ρ = 0.8 the orientation surgery at e = 1.00 roughly doubles the selected
rank (median 7–8, up to 14) relative to that column's own anchor (median 4), while at ρ ≤ 0.5 it adds at most
one component. A strongly trending baseline plus a global feature permutation evidently creates additional
stage-predictive directions that the CV retains. It coincides with ρ = 0.8 being the *only* ρ at which
orientation's `delta` and `shape` responses fall (0.27–0.46 and 0.43–0.69, versus 0.60–0.88 and 0.91–1.00
at ρ = 0) — consistent with the latent-rank probe's finding that the projection artifact decays as rank
grows. It is recorded, not explained, here.

## Limitations

- 100 replicates give a Monte Carlo SE of about 0.05 at a rate of 0.5 and 0.032 at 0.88. The decision's
  confirmation margin is one SE by predeclaration; the chosen column would also clear a 2·SE margin (0.88 −
  0.065 = 0.815) but ρ = 0.8, n = 1200 would not.
- `magnitude` and `shape` were deliberately not run (decided 2026-09-04; both reached power 1.00 at n = 300
  in Phase 4 and only gain from larger n). Their operating characteristics at the chosen point are measured
  in the Phase 5 study itself, not here. The `p_dmp` change from 0.2 to 0.1 touches every mode, so Phase 4's
  magnitude and shape numbers are not transferable to this baseline without that measurement.
- The observed-angle medians in "What governs orientation power" are computed from the gitignored
  `merged.jsonl`, not from a committed CSV; the snippet below regenerates them from the records. Every
  rejection rate, eigengap, null-width, and rank figure comes from committed CSV/JSON.
- Comparisons along ρ at fixed n are paired at the baseline indicator draw (one matched family, ρ thresholds
  the same uniform block); comparisons across n are not paired. Neither pairing is exploited by the report's
  rate estimates, which are per-cell binomials.
- Each column's `0.00` row is one shared anchor cell; the two modes' nulls in a column are one measurement.
- The Type I acceptance target is two-sided in this config, which is why a conservative 0.02 control rate is
  flagged. Phase 4's one-sided reading is the one the readiness worklist uses; both readings pass on
  inflation.
- The examples README's pre-run cost estimate (~500 workstation-core-hours) was about 7× too high; it was
  taken from a rehearsal that oversubscribed BLAS threads. The measured cost is 70 core-hours on EPYC 7662
  cores, and about 35 on a modern desktop core (12.5 s / 18.5 s / 27 s per unit at n = 300 / 600 / 1200,
  single-threaded, 199 permutations). The README is corrected in this change.

## Decision and hand-off

1. **Phase 5 design point: ρ = 0 (independent baseline), n = 1200 (300 samples per group-stage cell), four
   stages, `p_dmp = 0.1`, pooled PLS on M-value methylation with CV-selected rank.** Orientation `angle`
   power at e = 1.00 is 0.88 (1·SE lower bound 0.85) with the column's anchor at nominal level on all three
   statistics. The 0.80 floor stands; the claim is **not revised**.
2. **Report the eigengap, not ρ, with any orientation claim.** At the chosen point the pooled eigengap is
   0.049 (terciles 0.038–0.063). That distribution is what a real cohort must be compared against; the
   continuity knob is a generator parameter and does not transfer.
3. **Do not use baseline continuity to buy orientation power.** ρ = 0.5 lowers power at every n and ρ = 0.8
   does not beat ρ = 0 beyond n = 300. The orientation surgery's realized contrast is not ρ-invariant; state
   this as a construction limitation in the Phase 5 report contract.
4. **Hand the retained-rank question to item 3** with one new fact: rank is 3 at the chosen point regardless
   of n, and only a strongly trending baseline (ρ = 0.8) moves the CV off `n_stages − 1`.
5. **Phase 5 must re-measure magnitude and shape** at the chosen point, since `p_dmp` changed.

## Reproduction

The run used the committed config at revision `ff68514`. Locally (about 35 core-hours on a modern desktop
core; ~3.5 hours on 22 single-thread shards):

```bash
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
for i in $(seq 0 21); do
  uv run python scripts/run_study_shard.py \
    --config examples/trajectory_power_study/phase5_design_point_pilot.json \
    --out-dir results/phase5-design-point-2026-09-08 \
    --shard-index "$i" --n-shards 22 --error-policy record &
done
wait
```

On SLURM, as actually executed (partition and resource flags are cluster-specific; do NOT set
`STUDY_N_JOBS`, it changes the permutation draws and the parameter signature):

```bash
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
sbatch -p 512x1024 --cpus-per-task=1 --mem=2G --time=6:00:00 --array=0-99 \
  --export=ALL,STUDY_CONFIG=$PWD/examples/trajectory_power_study/phase5_design_point_pilot.json,STUDY_OUT=$PWD/results/phase5-design-point-2026-09-08,N_SHARDS=100 \
  scripts/motco_study_array.sbatch
```

Then, from the repository root:

```bash
uv run python scripts/motco_study.py merge  --out-dir results/phase5-design-point-2026-09-08
uv run python scripts/motco_study.py report \
  --config examples/trajectory_power_study/phase5_design_point_pilot.json \
  --out-dir results/phase5-design-point-2026-09-08
```

Record-derived observed-angle medians (the only figures above not in a committed CSV):

```bash
uv run python - <<'EOF'
import json, collections, statistics
by = collections.defaultdict(list)
for line in open("results/phase5-design-point-2026-09-08/merged.jsonl"):
    r = json.loads(line); cm = r["cell_metadata"]
    if cm["trajectory_mode"] != "orientation" or cm["effect_size"] != 1.0: continue
    dp = cm.get("design_point") or {"generator.baseline_continuity": 0.0, "generator.n_samples": 300}
    by[(dp["generator.baseline_continuity"], dp["generator.n_samples"])].append(r["pair_statistics"]["angle"])
for k in sorted(by): print(k, round(statistics.median(by[k]), 1))
EOF
```

The per-shard and merged JSONL (43 MB) are gitignored as regenerable; `report/` and `PROVENANCE.txt` are
committed. Every claim above traces to `results/phase5-design-point-2026-09-08/report/`:
`design_point_decision.json` / `.csv`, `design_point_operating.csv`, `continuity_resolved_orientation.csv`,
`power_curves.csv`, `type_i_table.csv`, `acceptance_report.json`, `realized_surgery.csv`,
`eigengap_stratified_power.csv`, `config_spectrum.csv`, and the figures `design_point_power.png`,
`power_curves.png`, `specificity_matrix.png`, `type_i.png`.
