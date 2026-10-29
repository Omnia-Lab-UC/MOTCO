# Design-guidance pilot — findings report (2026-10-28)

**Change:** `design-guidance-grid` · **Configs:** `examples/trajectory_power_study/design_guidance_{core,arms}_{magnitude,geometry}_pilot.json` · **Code revision:** `7c218c6` (branch `feat/design-guidance-grid`) · **Records:** `results/design-guidance-pilot-2026-10-28/<profile>/` (`report/`, `PROVENANCE.txt`) and `results/design-guidance-pilot-2026-10-28/figures/` · **Grade:** pilot, 100 replicates × 199 permutations; every rate below carries its Monte Carlo standard error at R = 100 (SE ≤ 0.05) and the anchor bound is α + 2·√(α(1 − α)/R) = 0.0936. **This report proposes the slice to rerun at paper grade; it makes no paper-grade claim.**

**Run:** four SLURM arrays on `ing` (partition `512x1024`, one CPU and 4 GB per task, BLAS single-threaded, `STUDY_N_JOBS` unset, `--error-policy record`), jobs 929160–929163, 124 array tasks, 47,200 units at 199 permutations, 124 / 124 tasks `COMPLETED 0:0`, **0 failed records, 0 censored surgeries**, 345.3 recorded unit core-hours (362.4 by `sacct` wall-sum), 3 h 51 m wall. Code: `7c218c6` on a clean worktree of `feat/design-guidance-grid`. The four `PROVENANCE.txt` files carry the launch lines, hashes and timings.

**In one paragraph.** Under the production measurement (pooled PLS, CV rank floored at 2, net displacement), **magnitude/`delta`** is the well-powered component: at n = 80 with three stages and two blocks it rejects 0.23 ± 0.04 at e = 0.05 and 0.59 ± 0.05 at e = 0.1, at n = 200 0.72 ± 0.04 and 0.94 ± 0.02, and e = 0.25 is saturated (≥ 0.98) at every column; a third block raises every sub-saturation rate (0.49 / 0.81 at n = 80), and a fourth stage changes little. **Shape/`shape`** rises monotonically with n, e and block count: 0.24 / 0.40 / 0.55 at n = 200 (two blocks, e = 0.25 / 0.5 / 1.0), 0.50 / 0.72 / 0.85 at n = 500, 0.82 / 0.92 / 0.97 at n = 1000, and ≤ 0.20 at n = 80 with two blocks (0.52 at e = 1.0 with three). **Orientation/`angle`** is governed by n and block count, not by the effect size: with two blocks it stays at 0.05–0.22 for n ≤ 200 and reaches 0.41–0.62 only at n = 1000; with three blocks 0.35–0.46 at n = 200 and 0.75–0.87 at n = 1000; its response is flat or non-monotone in e (at five of the six two-block columns with n ≥ 200 the e = 1.0 rate is below the e = 0.5 rate), and across the arms the e = 1.0 rate tracks the CV-selected rank (0.97 at median rank 5, 0.81 at rank 7, 0.12 at the baseline's rank 3), the pattern the Phase 6 study recorded. The held-out CV stage AUROC is 1.000 at every core column and every arm except the stage-effect arm, where δ = 0.05 and 0.1 give 0.553 and 0.729 — bracketing the SEA-AD cohorts' 0.64–0.70 — and there **no statistic has power at any effect at n = 200** (magnitude ≤ 0.03, orientation ≤ 0.07, shape ≤ 0.09). Size holds at 26 of 28 distinct anchors; the two `excess` verdicts are both exactly 0.10 against the bound 0.0936 (`delta` at n = 80 / four stages / three blocks, `angle` at the n = 200 baseline), the count chance predicts over 84 anchor tests at R = 100 (2.4), and the baseline column's independent Type I-phase draw gives `angle` 0.03. Translation's `angle` is within the pilot bound at every column (0.01–0.09). The paper-grade slice proposed in §7 is the three-stage core at n ∈ {80, 200, 500} × {two, three blocks} (36,000 units, about 350–400 core-hours).

## 1. Configuration and provenance

| Profile | Config (`sha256` prefix) | Cells × R | Units | Job | Shards | Unit runtime, median (min–max) | Recorded core-hours |
|---|---|---:|---:|---|---:|---|---:|
| core magnitude | `design_guidance_core_magnitude_pilot.json` (`c2bd261307a10bd3`) | 98 × 100 | 9,800 | 929160 | 30 | 21.4 s (7.5–95.0) | 68.5 |
| core geometry | `design_guidance_core_geometry_pilot.json` (`57942119d892ba2b`) | 162 × 100 | 16,200 | 929161 | 49 | 28.4 s (12.7–106.7) | 148.8 |
| arms magnitude | `design_guidance_arms_magnitude_pilot.json` (`e33e8c4eb6ec3f60`) | 80 × 100 | 8,000 | 929162 | 17 | 23.0 s (17.2–27.6) | 51.2 |
| arms geometry | `design_guidance_arms_geometry_pilot.json` (`be4e50300d54edee`) | 132 × 100 | 13,200 | 929163 | 28 | 20.8 s (15.7–28.1) | 76.8 |

- **Submission** (2026-10-28 13:30:57 cluster time, all 124 tasks running within a minute): `sbatch -p 512x1024 --cpus-per-task=1 --mem=4G --time=6:00:00 --array=0-<K−1> --export=ALL,STUDY_CONFIG=…,STUDY_OUT=…,N_SHARDS=<K> scripts/motco_study_array.sbatch` with `OMP_NUM_THREADS=OPENBLAS_NUM_THREADS=MKL_NUM_THREADS=1` exported and `STUDY_N_JOBS` unset; units are mapped to shards by the stable `(cell_id, replicate_index)` hash. The last task ended 17:22:20 (3 h 51 m wall). The task 5.2 estimate in the change notes was 303 core-hours; the core geometry profile ran 148.8 against 121 estimated (its n = 1000, four-stage, three-block units ran up to 107 s), core magnitude ran within its estimate (68.5 against 73) and the arm profiles above theirs (51.2 against 41, 76.8 against 68).
- **Code and software.** Revision `7c218c6` (branch `feat/design-guidance-grid`), run from a clean worktree (`/home1/tgonza/MOTCO-design-guidance`) with `PYTHONPATH` on its `src/` and the cluster's existing virtual environment: python 3.11.16, numpy 2.3.5, scikit-learn 1.8.0, scipy 1.16.3, pandas 2.3.3. The code at `7c218c6` is motco 0.10.0 (`motco.__version__`, verified from the worktree); the virtual environment's distribution metadata still names 0.6.0 because its editable install predates, and `PROVENANCE.txt` records that string. The code revision is the authoritative version.
- **Measurement, recorded per unit.** Every record carries `n_jobs` 1, `orientation_estimator` `net_displacement` and `min_components` 2 in `runtime_metadata`; every unit ran 199 permutations.
- **Records.** `merged.jsonl` per profile (sha256 in each `PROVENANCE.txt`: `ca372700…`, `fc76e0c6…`, `9aca73b4…`, `cba378c9…`) and the shard files stay on `ing`, gitignored; `report/` (16 files per profile, 17 for arms geometry, which also writes `continuity_resolved_orientation.csv`) and `PROVENANCE.txt` are committed, as are the four figures and `design_guidance_power.csv` under `figures/`.
- **Reproduction.** `python scripts/motco_study.py merge --out-dir results/design-guidance-pilot-2026-10-28/<profile>`; `python scripts/motco_study.py report --config examples/trajectory_power_study/design_guidance_<profile>_pilot.json --out-dir results/design-guidance-pilot-2026-10-28/<profile>`; `python scripts/design_guidance_figure.py --magnitude results/design-guidance-pilot-2026-10-28/core_magnitude/report --geometry results/design-guidance-pilot-2026-10-28/core_geometry/report --out-dir results/design-guidance-pilot-2026-10-28/figures`. Pre-flight: the fast suite at the same tree (ruff, mypy, 1510 passed / 1 skipped), the 2 × 19 smoke of all four profiles on `ing` (880 / 880 units, 0 failures) and 16 timed units at 199 permutations (`openspec/changes/design-guidance-grid/notes.md`).

Shared design (one baseline, one seed root `base_seed` 1028, one matched-seed family `design-guidance`): n = 200, three stages, `group_ratio` 0.5, `p_dmp` 0.1, ρ = 0, every per-omic δ 2.0, `magnitude_kind: joint`; pooled PLS on methylation + expression with the stage-supervised double-CV rank (`cv1_splits` 3, `cv2_splits` 4, 5 repeats, ≤ 20 components, `random_state` 1203), `orientation_estimator: net_displacement`, `min_components: 2`, 199 RRPP permutations, `n_jobs` 1. The production measurement is named in every profile, not defaulted. Attribution and the Phase 4 gate are off; the only predeclared rule is the advisory `anchor_decision` (α 0.05, k = 2, `delta`/`angle`/`shape`).

- **Core** (design grid `n_samples` {80, 200, 500, 1000} × `n_stages` {3, 4} × `layers` {methylation + expression, all three}): magnitude at e ∈ {0.01, 0.02, 0.05, 0.1, 0.25} (98 cells, 9,800 units); orientation, shape and translation at e ∈ {0.25, 0.5, 1.0} (162 cells, 16,200 units).
- **Arms** (explicit points, one setting moved off the baseline): group ratio {0.65, 0.8}; stage proportions 0.5/0.3/0.2; `p_dmp` {0.05, 0.15}; joint δ {0.05, 0.1, 0.25, 0.5, 1.0}; ρ {0.5, 0.8}. Magnitude 80 cells / 8,000 units; geometry 132 cells / 13,200 units.

The magnitude and geometry siblings share their seed family, so at every design column their zero-effect anchors evaluate the same datasets with the same RRPP draws (the smoke verified identical p-values at every column). This report reads the **geometry profile's anchor** and uses the magnitude profile's only as the reproduction check in §2; the two are never pooled.

## 2. Unit and failure accounting

- **Units.** 9,800 / 16,200 / 8,000 / 13,200 records, each `(cell_id, replicate_index)` present exactly once: duplicates 0, missing 0, unexpected 0, parameter-signature mismatches 0 in every profile (`PROVENANCE.txt`, computed against the enumerated grids).
- **Failures.** 0 failed records in every profile. The design risk — stage-stratified double CV failing or selecting erratically in the n = 80, four-stage cells of 10 samples — did not materialize: those columns completed every unit and selected rank 3–4 (§4).
- **Censored surgeries.** 0 (`report/realized_surgery.csv`: `censored_fraction` 0 at every orientation, shape and translation cell), as the headroom check at enumeration predicted.
- **Anchor reproduction (design D2).** Keyed by design column (generator parameters without the mode and effect, plus the evaluation parameters) and replicate, the 28 distinct anchor columns give 2,800 keys: 2,700 hold one record from each sibling (16 core and 11 non-baseline arm columns × 100 replicates) and 100 hold eight (the baseline column: its design-point anchor and its `type_i_baseline`-phase anchor in each of the four profiles). Records sharing a key and the same `(generator_seed, evaluation_seed)` have **identical p-values** in every case (3,300 pairwise comparisons, 0 mismatches; observed statistics agree to 4 × 10⁻¹¹), so the siblings share both data and RRPP draws and the magnitude profile's anchors reproduce the geometry profile's exactly. The baseline column's Type I-phase anchor is a separate seed derivation (`type_i_baseline` phase, p-values identical across the four profiles), which gives that column a second, independent 100-replicate size reading (§5). This report reads the geometry profiles' design-point anchors; nothing is pooled.

## 3. Per-column operating characteristics

Every rate is the rejection rate at p ≤ 0.05 over 100 available replicates ± its Monte Carlo SE, read from `report/design_point_operating.csv`. Column labels give `n_samples`, stage count and measured block count (two = methylation + expression; three = all blocks). The figures are `results/design-guidance-pilot-2026-10-28/figures/design_guidance_{magnitude,orientation,shape}.png` (one panel per stage count × block count, n on a log axis, one line per effect, ±1 SE bars) and `design_guidance_control.png` (§5); `design_guidance_power.csv` beside them holds exactly the rows plotted.

### 3.1 Core grid — magnitude / `delta`

| column | magnitude/delta e=0.01 | magnitude/delta e=0.02 | magnitude/delta e=0.05 | magnitude/delta e=0.1 | magnitude/delta e=0.25 |
|---|---|---|---|---|---|
| n = 80, 3 st, 2 bl | 0.01 ± 0.01 | 0.04 ± 0.02 | 0.23 ± 0.04 | 0.59 ± 0.05 | 0.99 ± 0.01 |
| n = 80, 3 st, 3 bl | 0.04 ± 0.02 | 0.09 ± 0.03 | 0.49 ± 0.05 | 0.81 ± 0.04 | 0.99 ± 0.01 |
| n = 80, 4 st, 2 bl | 0.01 ± 0.01 | 0.01 ± 0.01 | 0.13 ± 0.03 | 0.59 ± 0.05 | 0.98 ± 0.01 |
| n = 80, 4 st, 3 bl | 0.13 ± 0.03 | 0.20 ± 0.04 | 0.47 ± 0.05 | 0.81 ± 0.04 | 0.98 ± 0.01 |
| n = 200, 3 st, 2 bl | 0.05 ± 0.02 | 0.21 ± 0.04 | 0.72 ± 0.04 | 0.94 ± 0.02 | 1.00 ± 0.00 |
| n = 200, 3 st, 3 bl | 0.10 ± 0.03 | 0.39 ± 0.05 | 0.87 ± 0.03 | 1.00 ± 0.00 | 1.00 ± 0.00 |
| n = 200, 4 st, 2 bl | 0.12 ± 0.03 | 0.19 ± 0.04 | 0.59 ± 0.05 | 0.95 ± 0.02 | 1.00 ± 0.00 |
| n = 200, 4 st, 3 bl | 0.17 ± 0.04 | 0.36 ± 0.05 | 0.93 ± 0.03 | 1.00 ± 0.00 | 1.00 ± 0.00 |
| n = 500, 3 st, 2 bl | 0.22 ± 0.04 | 0.66 ± 0.05 | 0.97 ± 0.02 | 1.00 ± 0.00 | 1.00 ± 0.00 |
| n = 500, 3 st, 3 bl | 0.40 ± 0.05 | 0.85 ± 0.04 | 1.00 ± 0.00 | 1.00 ± 0.00 | 1.00 ± 0.00 |
| n = 500, 4 st, 2 bl | 0.15 ± 0.04 | 0.57 ± 0.05 | 0.97 ± 0.02 | 1.00 ± 0.00 | 1.00 ± 0.00 |
| n = 500, 4 st, 3 bl | 0.36 ± 0.05 | 0.86 ± 0.03 | 1.00 ± 0.00 | 1.00 ± 0.00 | 1.00 ± 0.00 |
| n = 1000, 3 st, 2 bl | 0.55 ± 0.05 | 0.89 ± 0.03 | 0.99 ± 0.01 | 1.00 ± 0.00 | 1.00 ± 0.00 |
| n = 1000, 3 st, 3 bl | 0.73 ± 0.04 | 0.99 ± 0.01 | 1.00 ± 0.00 | 1.00 ± 0.00 | 1.00 ± 0.00 |
| n = 1000, 4 st, 2 bl | 0.48 ± 0.05 | 0.91 ± 0.03 | 1.00 ± 0.00 | 1.00 ± 0.00 | 1.00 ± 0.00 |
| n = 1000, 4 st, 3 bl | 0.76 ± 0.04 | 0.97 ± 0.02 | 1.00 ± 0.00 | 1.00 ± 0.00 | 1.00 ± 0.00 |

**Reading.** `delta` power rises with n at every effect and saturates at e = 0.25 everywhere (0.98–1.00); the rise is resolved by e = 0.01–0.1. At the case-study size (n = 80, three stages, two blocks) it is 0.23 ± 0.04 at e = 0.05 and 0.59 ± 0.05 at e = 0.1; at n = 200, 0.72 / 0.94; at n = 500 e = 0.05 is saturated (0.97) and e = 0.02 resolves (0.57–0.86); at n = 1000 e = 0.01 resolves (0.48–0.76). A **third block raises every sub-saturation rate** (n = 80: 0.49 / 0.81 at e = 0.05 / 0.1; n = 200: 0.87 / 1.00; n = 500 at e = 0.02: 0.85–0.86 against 0.57–0.66), with the feature-count confound of §8. A **fourth stage changes little**: at two blocks the four-stage column is within 0.13 of the three-stage one at every n and effect (0.59 against 0.72 at n = 200, e = 0.05, the largest gap; 0.59 against 0.59 at n = 80, e = 0.1), and at three blocks within 0.07 for n ≥ 200. The n = 80 / four-stage / three-block column's e = 0.01–0.02 rates (0.13, 0.20) sit on an anchor `delta` of 0.10 (§5) and are not read as power.

### 3.2 Core grid — orientation / `angle` and shape / `shape`

| column | orientation/angle e=0.25 | orientation/angle e=0.5 | orientation/angle e=1 | shape/shape e=0.25 | shape/shape e=0.5 | shape/shape e=1 |
|---|---|---|---|---|---|---|
| n = 80, 3 st, 2 bl | 0.05 ± 0.02 | 0.06 ± 0.02 | 0.15 ± 0.04 | 0.02 ± 0.01 | 0.14 ± 0.03 | 0.20 ± 0.04 |
| n = 80, 3 st, 3 bl | 0.24 ± 0.04 | 0.22 ± 0.04 | 0.23 ± 0.04 | 0.20 ± 0.04 | 0.28 ± 0.04 | 0.52 ± 0.05 |
| n = 80, 4 st, 2 bl | 0.09 ± 0.03 | 0.09 ± 0.03 | 0.13 ± 0.03 | 0.02 ± 0.01 | 0.04 ± 0.02 | 0.19 ± 0.04 |
| n = 80, 4 st, 3 bl | 0.16 ± 0.04 | 0.23 ± 0.04 | 0.21 ± 0.04 | 0.13 ± 0.03 | 0.34 ± 0.05 | 0.46 ± 0.05 |
| n = 200, 3 st, 2 bl | 0.17 ± 0.04 | 0.22 ± 0.04 | 0.12 ± 0.03 | 0.24 ± 0.04 | 0.40 ± 0.05 | 0.55 ± 0.05 |
| n = 200, 3 st, 3 bl | 0.40 ± 0.05 | 0.41 ± 0.05 | 0.38 ± 0.05 | 0.43 ± 0.05 | 0.55 ± 0.05 | 0.84 ± 0.04 |
| n = 200, 4 st, 2 bl | 0.10 ± 0.03 | 0.09 ± 0.03 | 0.11 ± 0.03 | 0.11 ± 0.03 | 0.29 ± 0.05 | 0.59 ± 0.05 |
| n = 200, 4 st, 3 bl | 0.39 ± 0.05 | 0.46 ± 0.05 | 0.35 ± 0.05 | 0.50 ± 0.05 | 0.71 ± 0.05 | 0.76 ± 0.04 |
| n = 500, 3 st, 2 bl | 0.33 ± 0.05 | 0.39 ± 0.05 | 0.29 ± 0.05 | 0.50 ± 0.05 | 0.72 ± 0.04 | 0.85 ± 0.04 |
| n = 500, 3 st, 3 bl | 0.67 ± 0.05 | 0.69 ± 0.05 | 0.61 ± 0.05 | 0.75 ± 0.04 | 0.88 ± 0.03 | 0.99 ± 0.01 |
| n = 500, 4 st, 2 bl | 0.27 ± 0.04 | 0.28 ± 0.04 | 0.13 ± 0.03 | 0.60 ± 0.05 | 0.75 ± 0.04 | 0.86 ± 0.03 |
| n = 500, 4 st, 3 bl | 0.62 ± 0.05 | 0.76 ± 0.04 | 0.73 ± 0.04 | 0.86 ± 0.03 | 0.91 ± 0.03 | 0.98 ± 0.01 |
| n = 1000, 3 st, 2 bl | 0.62 ± 0.05 | 0.61 ± 0.05 | 0.49 ± 0.05 | 0.82 ± 0.04 | 0.92 ± 0.03 | 0.97 ± 0.02 |
| n = 1000, 3 st, 3 bl | 0.82 ± 0.04 | 0.75 ± 0.04 | 0.76 ± 0.04 | 0.89 ± 0.03 | 0.92 ± 0.03 | 0.99 ± 0.01 |
| n = 1000, 4 st, 2 bl | 0.60 ± 0.05 | 0.60 ± 0.05 | 0.41 ± 0.05 | 0.78 ± 0.04 | 0.93 ± 0.03 | 0.94 ± 0.02 |
| n = 1000, 4 st, 3 bl | 0.81 ± 0.04 | 0.87 ± 0.03 | 0.84 ± 0.04 | 0.96 ± 0.02 | 0.96 ± 0.02 | 1.00 ± 0.00 |

**Orientation / `angle`.** Power is governed by **sample size and block count, not by the effect size**. With two blocks it stays at 0.05–0.15 (n = 80), 0.09–0.22 (n = 200) and 0.13–0.39 (n = 500), and reaches 0.41–0.62 only at n = 1000; with three blocks it is 0.16–0.24 at n = 80, 0.35–0.46 at n = 200, 0.61–0.76 at n = 500 and 0.75–0.87 at n = 1000. The response in e is flat or non-monotone: at five of the six two-block columns with n ≥ 200 the e = 1.0 rate is **below** the e = 0.5 rate (n = 200 / 3 st: 0.12 against 0.22; n = 500 / 3 st: 0.29 against 0.39; n = 500 / 4 st: 0.13 against 0.28; n = 1000 / 3 st: 0.49 against 0.61; n = 1000 / 4 st: 0.41 against 0.60), while the three-block columns show no consistent direction in e. The arms (§3.4) and the Phase 6 study point at the CV-selected rank: the e = 1.0 orientation surgery relocates every stage-active methylation site of group B, and whether the two groups' net displacements separate in the latent space depends on how many components the stage-supervised CV retains (median 2–4 at every core column, §4; 0.280 against 0.996 below and above three components in the Phase 6 study). The orientation construction at e = 1.0 also moves the other two statistics (`delta` 0.12–0.92 and `shape` 0.19–1.00 across the core columns), the cross-talk the Phase 5 response matrix recorded; the figure plots the target statistic only.

**Shape / `shape`.** Monotone in e, n and block count. Two blocks: 0.02 / 0.14 / 0.20 at n = 80, 0.24 / 0.40 / 0.55 at n = 200, 0.50 / 0.72 / 0.85 at n = 500, 0.82 / 0.92 / 0.97 at n = 1000 (e = 0.25 / 0.5 / 1.0); three blocks: 0.20 / 0.28 / 0.52, 0.43 / 0.55 / 0.84, 0.75 / 0.88 / 0.99, 0.89 / 0.92 / 0.99. Four stages are comparable to three at n ≥ 500 (n = 500, two blocks: 0.60 / 0.75 / 0.86) and lower at e ≤ 0.5 for n ≤ 200 (0.11 / 0.29 against 0.24 / 0.40 at n = 200, two blocks). Shape never leaks into `angle` (≤ 0.10 at e = 1.0 at every column) but does move `delta` (0.25–0.96), as in Phase 5.

#### Translation (negative control), all statistics at e = 1.0

| column | delta | angle | shape |
|---|---|---|---|
| n = 80, 3 st, 2 bl | 0.03 ± 0.02 | 0.02 ± 0.01 | 0.03 ± 0.02 |
| n = 80, 3 st, 3 bl | 0.07 ± 0.03 | 0.03 ± 0.02 | 0.03 ± 0.02 |
| n = 80, 4 st, 2 bl | 0.01 ± 0.01 | 0.06 ± 0.02 | 0.01 ± 0.01 |
| n = 80, 4 st, 3 bl | 0.03 ± 0.02 | 0.07 ± 0.03 | 0.04 ± 0.02 |
| n = 200, 3 st, 2 bl | 0.05 ± 0.02 | 0.07 ± 0.03 | 0.02 ± 0.01 |
| n = 200, 3 st, 3 bl | 0.05 ± 0.02 | 0.05 ± 0.02 | 0.03 ± 0.02 |
| n = 200, 4 st, 2 bl | 0.02 ± 0.01 | 0.04 ± 0.02 | 0.03 ± 0.02 |
| n = 200, 4 st, 3 bl | 0.03 ± 0.02 | 0.07 ± 0.03 | 0.05 ± 0.02 |
| n = 500, 3 st, 2 bl | 0.07 ± 0.03 | 0.03 ± 0.02 | 0.02 ± 0.01 |
| n = 500, 3 st, 3 bl | 0.09 ± 0.03 | 0.02 ± 0.01 | 0.07 ± 0.03 |
| n = 500, 4 st, 2 bl | 0.03 ± 0.02 | 0.02 ± 0.01 | 0.06 ± 0.02 |
| n = 500, 4 st, 3 bl | 0.03 ± 0.02 | 0.01 ± 0.01 | 0.06 ± 0.02 |
| n = 1000, 3 st, 2 bl | 0.06 ± 0.02 | 0.03 ± 0.02 | 0.05 ± 0.02 |
| n = 1000, 3 st, 3 bl | 0.05 ± 0.02 | 0.03 ± 0.02 | 0.06 ± 0.02 |
| n = 1000, 4 st, 2 bl | 0.08 ± 0.03 | 0.07 ± 0.03 | 0.05 ± 0.02 |
| n = 1000, 4 st, 3 bl | 0.06 ± 0.02 | 0.07 ± 0.03 | 0.06 ± 0.02 |

### 3.3 Arms — magnitude / `delta` (off n = 200, three stages, two blocks)

| column | magnitude/delta e=0.01 | magnitude/delta e=0.02 | magnitude/delta e=0.05 | magnitude/delta e=0.1 | magnitude/delta e=0.25 |
|---|---|---|---|---|---|
| baseline | 0.05 ± 0.02 | 0.21 ± 0.04 | 0.72 ± 0.04 | 0.94 ± 0.02 | 1.00 ± 0.00 |
| group ratio = 0.65 | 0.12 ± 0.03 | 0.30 ± 0.05 | 0.67 ± 0.05 | 0.98 ± 0.01 | 1.00 ± 0.00 |
| group ratio = 0.8 | 0.10 ± 0.03 | 0.21 ± 0.04 | 0.53 ± 0.05 | 0.86 ± 0.03 | 1.00 ± 0.00 |
| stage proportions 0.5/0.3/0.2 | 0.07 ± 0.03 | 0.18 ± 0.04 | 0.58 ± 0.05 | 0.87 ± 0.03 | 1.00 ± 0.00 |
| p_dmp = 0.05 | 0.04 ± 0.02 | 0.04 ± 0.02 | 0.20 ± 0.04 | 0.55 ± 0.05 | 0.98 ± 0.01 |
| p_dmp = 0.15 | 0.17 ± 0.04 | 0.41 ± 0.05 | 0.85 ± 0.04 | 0.99 ± 0.01 | 1.00 ± 0.00 |
| δ = 0.05 | 0.02 ± 0.01 | 0.02 ± 0.01 | 0.02 ± 0.01 | 0.02 ± 0.01 | 0.01 ± 0.01 |
| δ = 0.1 | 0.03 ± 0.02 | 0.03 ± 0.02 | 0.03 ± 0.02 | 0.02 ± 0.01 | 0.02 ± 0.01 |
| δ = 0.25 | 0.08 ± 0.03 | 0.06 ± 0.02 | 0.07 ± 0.03 | 0.10 ± 0.03 | 0.33 ± 0.05 |
| δ = 0.5 | 0.03 ± 0.02 | 0.04 ± 0.02 | 0.12 ± 0.03 | 0.29 ± 0.05 | 0.68 ± 0.05 |
| δ = 1 | 0.02 ± 0.01 | 0.02 ± 0.01 | 0.17 ± 0.04 | 0.45 ± 0.05 | 0.97 ± 0.02 |
| ρ = 0.5 | 0.03 ± 0.02 | 0.03 ± 0.02 | 0.17 ± 0.04 | 0.42 ± 0.05 | 0.94 ± 0.02 |
| ρ = 0.8 | 0.02 ± 0.01 | 0.01 ± 0.01 | 0.05 ± 0.02 | 0.27 ± 0.04 | 0.67 ± 0.05 |

### 3.4 Arms — orientation / `angle` and shape / `shape`

| column | orientation/angle e=0.25 | orientation/angle e=0.5 | orientation/angle e=1 | shape/shape e=0.25 | shape/shape e=0.5 | shape/shape e=1 |
|---|---|---|---|---|---|---|
| baseline | 0.17 ± 0.04 | 0.22 ± 0.04 | 0.12 ± 0.03 | 0.24 ± 0.04 | 0.40 ± 0.05 | 0.55 ± 0.05 |
| group ratio = 0.65 | 0.16 ± 0.04 | 0.18 ± 0.04 | 0.97 ± 0.02 | 0.17 ± 0.04 | 0.33 ± 0.05 | 0.39 ± 0.05 |
| group ratio = 0.8 | 0.14 ± 0.03 | 0.19 ± 0.04 | 0.70 ± 0.05 | 0.28 ± 0.04 | 0.68 ± 0.05 | 0.87 ± 0.03 |
| stage proportions 0.5/0.3/0.2 | 0.17 ± 0.04 | 0.09 ± 0.03 | 0.12 ± 0.03 | 0.17 ± 0.04 | 0.35 ± 0.05 | 0.56 ± 0.05 |
| p_dmp = 0.05 | 0.06 ± 0.02 | 0.06 ± 0.02 | 0.12 ± 0.03 | 0.06 ± 0.02 | 0.10 ± 0.03 | 0.11 ± 0.03 |
| p_dmp = 0.15 | 0.14 ± 0.03 | 0.39 ± 0.05 | 0.37 ± 0.05 | 0.47 ± 0.05 | 0.74 ± 0.04 | 0.92 ± 0.03 |
| δ = 0.05 | 0.06 ± 0.02 | 0.02 ± 0.01 | 0.04 ± 0.02 | 0.02 ± 0.01 | 0.05 ± 0.02 | 0.07 ± 0.03 |
| δ = 0.1 | 0.07 ± 0.03 | 0.03 ± 0.02 | 0.03 ± 0.02 | 0.02 ± 0.01 | 0.05 ± 0.02 | 0.09 ± 0.03 |
| δ = 0.25 | 0.10 ± 0.03 | 0.08 ± 0.03 | 0.09 ± 0.03 | 0.02 ± 0.01 | 0.05 ± 0.02 | 0.05 ± 0.02 |
| δ = 0.5 | 0.08 ± 0.03 | 0.09 ± 0.03 | 0.46 ± 0.05 | 0.03 ± 0.02 | 0.07 ± 0.03 | 0.08 ± 0.03 |
| δ = 1 | 0.08 ± 0.03 | 0.03 ± 0.02 | 0.44 ± 0.05 | 0.05 ± 0.02 | 0.12 ± 0.03 | 0.15 ± 0.04 |
| ρ = 0.5 | 0.13 ± 0.03 | 0.05 ± 0.02 | 0.26 ± 0.04 | 0.64 ± 0.05 | 0.90 ± 0.03 | 0.93 ± 0.03 |
| ρ = 0.8 | 0.06 ± 0.02 | 0.14 ± 0.03 | 0.81 ± 0.04 | 0.76 ± 0.04 | 0.94 ± 0.02 | 0.99 ± 0.01 |

#### Translation (negative control), all statistics at e = 1.0

| column | delta | angle | shape |
|---|---|---|---|
| baseline | 0.05 ± 0.02 | 0.07 ± 0.03 | 0.02 ± 0.01 |
| group ratio = 0.65 | 0.02 ± 0.01 | 0.08 ± 0.03 | 0.05 ± 0.02 |
| group ratio = 0.8 | 0.04 ± 0.02 | 0.05 ± 0.02 | 0.06 ± 0.02 |
| stage proportions 0.5/0.3/0.2 | 0.06 ± 0.02 | 0.02 ± 0.01 | 0.03 ± 0.02 |
| p_dmp = 0.05 | 0.03 ± 0.02 | 0.09 ± 0.03 | 0.04 ± 0.02 |
| p_dmp = 0.15 | 0.04 ± 0.02 | 0.05 ± 0.02 | 0.02 ± 0.01 |
| δ = 0.05 | 0.03 ± 0.02 | 0.01 ± 0.01 | 0.04 ± 0.02 |
| δ = 0.1 | 0.02 ± 0.01 | 0.01 ± 0.01 | 0.03 ± 0.02 |
| δ = 0.25 | 0.03 ± 0.02 | 0.02 ± 0.01 | 0.04 ± 0.02 |
| δ = 0.5 | 0.01 ± 0.01 | 0.02 ± 0.01 | 0.02 ± 0.01 |
| δ = 1 | 0.01 ± 0.01 | 0.01 ± 0.01 | 0.04 ± 0.02 |
| ρ = 0.5 | 0.03 ± 0.02 | 0.03 ± 0.02 | 0.01 ± 0.01 |
| ρ = 0.8 | 0.02 ± 0.01 | 0.01 ± 0.01 | 0.04 ± 0.02 |

**Reading across the arms** (baseline: `delta` 0.72 / 0.94 at e = 0.05 / 0.1; `angle` 0.17 / 0.22 / 0.12 and `shape` 0.24 / 0.40 / 0.55 at e = 0.25 / 0.5 / 1.0):

- **Group ratio.** Magnitude moves little (0.67 / 0.98 at 0.65; 0.53 / 0.86 at 0.8). Shape is 0.39 at e = 1.0 at ratio 0.65 and 0.87 at 0.8. Orientation at e = 1.0 jumps to **0.97 ± 0.02 at ratio 0.65 and 0.70 ± 0.05 at 0.8**, where the CV-selected rank at that cell rises to a median of 5 and 3.5 (baseline 3); at e ≤ 0.5 both arms sit at the baseline's 0.14–0.19.
- **Stage proportions 0.5 / 0.3 / 0.2.** Close to the baseline for every component (magnitude 0.58 / 0.87 against 0.72 / 0.94; shape 0.17 / 0.35 / 0.56; orientation 0.17 / 0.09 / 0.12).
- **Differential-feature proportion.** `p_dmp` 0.05 lowers every component (magnitude 0.20 / 0.55; shape ≤ 0.11; orientation ≤ 0.12); 0.15 raises every one (magnitude 0.85 / 0.99; shape 0.47 / 0.74 / 0.92; orientation 0.14 / 0.39 / 0.37).
- **Stage effect δ.** At δ = 0.05 and 0.1 **every statistic is at the null at every effect** (magnitude ≤ 0.03, orientation ≤ 0.07, shape ≤ 0.09). δ = 0.25 gives magnitude 0.33 at e = 0.25 and the geometry modes at the null; δ = 0.5 gives magnitude 0.68, orientation 0.46 at e = 1.0 (median rank 9) and shape ≤ 0.08; δ = 1.0 gives magnitude 0.97, orientation 0.44 (rank 4) and shape 0.15. Every surgery scales with δ (magnitude multiplies it; orientation and shape relocate sites whose effect is δ), so group differences shrink with the stage signal by construction.
- **Baseline continuity ρ.** Continuity trades magnitude power for geometry power: at ρ = 0.5 magnitude is 0.17 / 0.42 / 0.94 at e = 0.05 / 0.1 / 0.25 and at ρ = 0.8 0.05 / 0.27 / 0.67 (baseline 0.72 / 0.94 / 1.00), while shape rises to 0.64 / 0.90 / 0.93 and 0.76 / 0.94 / 0.99 and orientation at e = 1.0 to 0.26 and 0.81 ± 0.04 (median rank 7). The arms geometry report's `continuity_resolved_orientation.csv` carries the per-ρ eigengap and `angle` null-width distributions.

## 4. Stage signal and selected rank

Read from `report/design_point_stage_signal.csv` (mean and q33–q67 of the recorded held-out CV stage AUROC at the zero-effect anchor; 100 of 100 available at every column) and `design_point_operating.csv` (`median_selected_lv`, min–max). The core table is the geometry profile's; the magnitude profile's is identical (§2).

#### Core grid

| column | anchor CV AUROC (mean, q33–q67) | anchor rank (median, min–max) | n AUROC |
|---|---|---|---|
| n = 80, 3 st, 2 bl | 1.000 (1.000–1.000) | 3 (2–3) | 100 |
| n = 80, 3 st, 3 bl | 1.000 (1.000–1.000) | 2 (2–3) | 100 |
| n = 80, 4 st, 2 bl | 1.000 (1.000–1.000) | 4 (3–4) | 100 |
| n = 80, 4 st, 3 bl | 1.000 (1.000–1.000) | 3 (3–4) | 100 |
| n = 200, 3 st, 2 bl | 1.000 (1.000–1.000) | 2 (2–3) | 100 |
| n = 200, 3 st, 3 bl | 1.000 (1.000–1.000) | 2 (2–3) | 100 |
| n = 200, 4 st, 2 bl | 1.000 (1.000–1.000) | 4 (3–4) | 100 |
| n = 200, 4 st, 3 bl | 1.000 (1.000–1.000) | 3 (2–4) | 100 |
| n = 500, 3 st, 2 bl | 1.000 (1.000–1.000) | 2 (2–3) | 100 |
| n = 500, 3 st, 3 bl | 1.000 (1.000–1.000) | 2 (2–2) | 100 |
| n = 500, 4 st, 2 bl | 1.000 (1.000–1.000) | 3 (3–4) | 100 |
| n = 500, 4 st, 3 bl | 1.000 (1.000–1.000) | 3 (2–3) | 100 |
| n = 1000, 3 st, 2 bl | 1.000 (1.000–1.000) | 2 (2–2) | 100 |
| n = 1000, 3 st, 3 bl | 1.000 (1.000–1.000) | 2 (2–2) | 100 |
| n = 1000, 4 st, 2 bl | 1.000 (1.000–1.000) | 3 (2–4) | 100 |
| n = 1000, 4 st, 3 bl | 1.000 (1.000–1.000) | 3 (2–3) | 100 |

#### Arms

| column | anchor CV AUROC (mean, q33–q67) | anchor rank (median, min–max) | n AUROC |
|---|---|---|---|
| baseline | 1.000 (1.000–1.000) | 2 (2–3) | 100 |
| group ratio = 0.65 | 1.000 (1.000–1.000) | 2 (2–3) | 100 |
| group ratio = 0.8 | 1.000 (1.000–1.000) | 2 (2–3) | 100 |
| stage proportions 0.5/0.3/0.2 | 1.000 (1.000–1.000) | 2 (2–3) | 100 |
| p_dmp = 0.05 | 1.000 (1.000–1.000) | 3 (2–3) | 100 |
| p_dmp = 0.15 | 1.000 (1.000–1.000) | 2 (2–3) | 100 |
| δ = 0.05 | 0.553 (0.537–0.565) | 3 (2–19) | 100 |
| δ = 0.1 | 0.729 (0.717–0.742) | 10 (3–19) | 100 |
| δ = 0.25 | 0.997 (0.996–0.998) | 10 (6–19) | 100 |
| δ = 0.5 | 1.000 (1.000–1.000) | 4 (3–5) | 100 |
| δ = 1 | 1.000 (1.000–1.000) | 3 (2–3) | 100 |
| ρ = 0.5 | 1.000 (1.000–1.000) | 3 (2–3) | 100 |
| ρ = 0.8 | 1.000 (1.000–1.000) | 4 (3–6) | 100 |

**Reading.** At δ = 2.0 the stage signal is **saturated at every core column and every non-δ arm** (AUROC 1.000 with a degenerate q33–q67 interval), and the CV selects rank 2–3 at three stages and 3–4 at four stages, never more than 4: the rank floor of 2 binds at the three-stage columns (median 2 at n ≥ 200). Only the δ arm leaves that regime: δ = 0.05 gives 0.553 (0.537–0.565) and δ = 0.1 gives 0.729 (0.717–0.742), **bracketing the SEA-AD cohorts' recorded 0.64–0.70** at n = 200, and δ = 0.25 is already 0.997. In the bracketed regime no statistic has power at any effect (§3.4). The selected rank behaves differently from the rank-floor check at n = 80: there the unfloored rule fell to rank 1; here, at n = 200, the double CV **inflates** the rank in the intermediate regime (median 10 at δ = 0.1 and 0.25, up to 19) and returns to 3–4 by δ = 0.5. Any `angle` reading from a real cohort must therefore be accompanied by its selected rank, as the Phase 6 statement already requires.

## 5. Size at every column

Zero-effect anchor rates per statistic (geometry profiles; the magnitude profiles reproduce them exactly, §2) with the advisory `anchor_decision` status against the bound α + 2·√(α(1 − α)/100) = 0.0936, and translation's `angle` rate at e = 1.0 (`report/anchor_decision.json` and `.csv`; figure `design_guidance_control.png`).

#### Core grid

| column | delta | angle | shape | bound | verdict | translation angle (top e) |
|---|---|---|---|---|---|---|
| n = 80, 3 st, 2 bl | 0.01 ± 0.01 | 0.04 ± 0.02 | 0.01 ± 0.01 | 0.0936 | size_holds | 0.02 ± 0.01 |
| n = 80, 3 st, 3 bl | 0.02 ± 0.01 | 0.01 ± 0.01 | 0.02 ± 0.01 | 0.0936 | size_holds | 0.03 ± 0.02 |
| n = 80, 4 st, 2 bl | 0.01 ± 0.01 | 0.08 ± 0.03 | 0.01 ± 0.01 | 0.0936 | size_holds | 0.06 ± 0.02 |
| n = 80, 4 st, 3 bl | 0.10 ± 0.03 | 0.07 ± 0.03 | 0.06 ± 0.02 | 0.0936 | excess | 0.07 ± 0.03 |
| n = 200, 3 st, 2 bl | 0.02 ± 0.01 | 0.10 ± 0.03 | 0.05 ± 0.02 | 0.0936 | excess | 0.07 ± 0.03 |
| n = 200, 3 st, 3 bl | 0.03 ± 0.02 | 0.07 ± 0.03 | 0.04 ± 0.02 | 0.0936 | size_holds | 0.05 ± 0.02 |
| n = 200, 4 st, 2 bl | 0.03 ± 0.02 | 0.04 ± 0.02 | 0.01 ± 0.01 | 0.0936 | size_holds | 0.04 ± 0.02 |
| n = 200, 4 st, 3 bl | 0.02 ± 0.01 | 0.06 ± 0.02 | 0.04 ± 0.02 | 0.0936 | size_holds | 0.07 ± 0.03 |
| n = 500, 3 st, 2 bl | 0.01 ± 0.01 | 0.08 ± 0.03 | 0.03 ± 0.02 | 0.0936 | size_holds | 0.03 ± 0.02 |
| n = 500, 3 st, 3 bl | 0.02 ± 0.01 | 0.02 ± 0.01 | 0.01 ± 0.01 | 0.0936 | size_holds | 0.02 ± 0.01 |
| n = 500, 4 st, 2 bl | 0.04 ± 0.02 | 0.03 ± 0.02 | 0.04 ± 0.02 | 0.0936 | size_holds | 0.02 ± 0.01 |
| n = 500, 4 st, 3 bl | 0.03 ± 0.02 | 0.03 ± 0.02 | 0.05 ± 0.02 | 0.0936 | size_holds | 0.01 ± 0.01 |
| n = 1000, 3 st, 2 bl | 0.05 ± 0.02 | 0.06 ± 0.02 | 0.04 ± 0.02 | 0.0936 | size_holds | 0.03 ± 0.02 |
| n = 1000, 3 st, 3 bl | 0.06 ± 0.02 | 0.05 ± 0.02 | 0.05 ± 0.02 | 0.0936 | size_holds | 0.03 ± 0.02 |
| n = 1000, 4 st, 2 bl | 0.03 ± 0.02 | 0.03 ± 0.02 | 0.03 ± 0.02 | 0.0936 | size_holds | 0.07 ± 0.03 |
| n = 1000, 4 st, 3 bl | 0.03 ± 0.02 | 0.06 ± 0.02 | 0.06 ± 0.02 | 0.0936 | size_holds | 0.07 ± 0.03 |

#### Arms

| column | delta | angle | shape | bound | verdict | translation angle (top e) |
|---|---|---|---|---|---|---|
| baseline | 0.02 ± 0.01 | 0.10 ± 0.03 | 0.05 ± 0.02 | 0.0936 | excess | 0.07 ± 0.03 |
| group ratio = 0.65 | 0.03 ± 0.02 | 0.04 ± 0.02 | 0.01 ± 0.01 | 0.0936 | size_holds | 0.08 ± 0.03 |
| group ratio = 0.8 | 0.04 ± 0.02 | 0.04 ± 0.02 | 0.02 ± 0.01 | 0.0936 | size_holds | 0.05 ± 0.02 |
| stage proportions 0.5/0.3/0.2 | 0.02 ± 0.01 | 0.04 ± 0.02 | 0.04 ± 0.02 | 0.0936 | size_holds | 0.02 ± 0.01 |
| p_dmp = 0.05 | 0.03 ± 0.02 | 0.07 ± 0.03 | 0.04 ± 0.02 | 0.0936 | size_holds | 0.09 ± 0.03 |
| p_dmp = 0.15 | 0.05 ± 0.02 | 0.06 ± 0.02 | 0.07 ± 0.03 | 0.0936 | size_holds | 0.05 ± 0.02 |
| δ = 0.05 | 0.02 ± 0.01 | 0.05 ± 0.02 | 0.06 ± 0.02 | 0.0936 | size_holds | 0.01 ± 0.01 |
| δ = 0.1 | 0.03 ± 0.02 | 0.06 ± 0.02 | 0.06 ± 0.02 | 0.0936 | size_holds | 0.01 ± 0.01 |
| δ = 0.25 | 0.05 ± 0.02 | 0.06 ± 0.02 | 0.03 ± 0.02 | 0.0936 | size_holds | 0.02 ± 0.01 |
| δ = 0.5 | 0.04 ± 0.02 | 0.06 ± 0.02 | 0.04 ± 0.02 | 0.0936 | size_holds | 0.02 ± 0.01 |
| δ = 1 | 0.01 ± 0.01 | 0.06 ± 0.02 | 0.05 ± 0.02 | 0.0936 | size_holds | 0.01 ± 0.01 |
| ρ = 0.5 | 0.03 ± 0.02 | 0.07 ± 0.03 | 0.03 ± 0.02 | 0.0936 | size_holds | 0.03 ± 0.02 |
| ρ = 0.8 | 0.01 ± 0.01 | 0.05 ± 0.02 | 0.03 ± 0.02 | 0.0936 | size_holds | 0.01 ± 0.01 |

**Reading.** Across the 28 distinct anchor columns the advisory decision is `EXCESS` in every profile because two points carry an `excess` verdict, both at **exactly 0.10 ± 0.03 against 0.0936**, the smallest excess a 100-replicate column can show: `delta` at n = 80 / four stages / three blocks, and `angle` at the n = 200 / three stages / two blocks baseline (the same measurement in all four profiles). Over the 84 anchor tests (28 columns × 3 statistics), a test with exact size 0.05 exceeds the bound with probability P(Bin(100, 0.05) ≥ 10) = 0.028, so **2.4 such excesses are expected by chance**; two were observed. The baseline column's independent Type I-phase anchor (`report/type_i_table.csv`, a separate 100-replicate draw) gives `angle` 0.03 ± 0.02 (`delta` 0.04, `shape` 0.00), 13 / 200 = 0.065 pooled over the two draws. Neither excess is evidence of inflated size, and the pilot cannot clear them either; both columns (n = 80 and the n = 200 baseline) are in the paper-grade slice (§7). Every other anchor rate is ≤ 0.08.

**Translation.** `angle` at e = 1.0 is 0.01–0.07 across the core columns and 0.01–0.09 across the arms (0.09 at `p_dmp` 0.05, 0.08 at group ratio 0.65), within the pilot bound everywhere; `delta` is 0.01–0.09 and `shape` 0.01–0.07. Seven of the 28 distinct translation `angle` rates sit at 0.07–0.09, above the paper-grade bound 0.0695 that the net-displacement studies read translation against (0.058–0.068 at R = 500), which the pilot's SE of 0.03 cannot resolve: the paper-grade slice keeps translation at every column.

**Other negative controls.** Magnitude at e = 0.25 moves neither `angle` (≤ 0.07) nor `shape` (≤ 0.04) at any core column; shape at e = 1.0 moves `angle` at ≤ 0.10.

## 6. The originally planned grid and what was run instead

The project's originally planned simulation grid (Table 1) crossed every row below; crossing them is infeasible (3,456 null combinations before any power cell), its mean shifts (0.01–0.2) sit far below the InterSIM stage effect the generator operates on (δ = 2.0), its "power" rows do not map to trajectory components, and two of them are null scenarios. Each factor's disposition:

| Table 1 factor | Disposition | Where, and why |
|---|---|---|
| Sample size | **Run** | Core axis `n_samples` ∈ {80, 200, 500, 1000}, crossed with stages and blocks. |
| Group ratio | **Run** | Arm: `group_ratio` 0.65 and 0.8 off 0.5. |
| Number of stages | **Run** | Core axis `n_stages` ∈ {3, 4}. Two-stage designs are dropped: with two stage means the configuration is a segment, `shape` is undefined by construction (nothing to bend) and `angle` reduces to the segment direction, so the three-component guidance cannot be read there; the case study has three stages. |
| Stage proportions | **Run** | Arm: 0.5/0.3/0.2 off equal thirds. |
| Differential-feature proportion | **Run** | Arm: `p_dmp` 0.05 and 0.15 off 0.1. 0.2 is excluded because at three stages it fails the surgery-headroom check (orientation saturates at e = 0.89, translation at 0.52). |
| Shared stage shift | **Run, rescaled** | Arm: the joint per-omic δ ∈ {0.05, 0.1, 0.25, 0.5, 1.0} off InterSIM's 2.0. Table 1's 0.01–0.2 are on a scale the generator never reaches; 0.05–0.1 is where the rank-floor check found the SEA-AD stage-signal regime. |
| Noise level | **Dropped** | Redundant with the stage shift: features are standardized before the PLS fit, so σ and δ enter the measurement only through δ/σ, which the δ arm already sweeps. |
| Group-specific mean shift | **Run** | The `translation` mode (a constant location offset of group B) at e ∈ {0.25, 0.5, 1.0} in both geometry profiles, read as the negative control for all three statistics; the `magnitude` mode covers a group-specific *size* change. |
| Covariance scaling | **Deferred** | To the group-specific-dispersion change (change B): an opt-in generator field scaling group B's covariance and a null-only Type I study crossing spread {1, 1.5, 2, 3} × group ratio {0.2, 0.5, 0.8}. It is a size question, not a power row. |
| Correlation scaling | **Deferred** | Same change and the same reason: it changes group B's dispersion, not its trajectory. |
| Feature filtering | **Dropped** | Supervised filtering (by group or stage) biases the permutation null unless it is repeated inside every permutation, so MOTCO documents it as forbidden; label-blind filtering is upstream preprocessing and changes what is measured, not the validity of the test. |
| Integration method | **Dropped** | PLS is the production measurement (`docs/roadmap.md`); SNF with the Euclidean statistics is deferred pending graph-native statistics; MOFA is out of scope; `concat` is a diagnostic baseline, not a constructed latent space. |

## 7. Proposed paper-grade slice

The pilot answers the design question descriptively at every column; what the paper needs at paper grade (500 replicates × 999 permutations, SE ≤ 0.022, anchor bound 0.0695) is the slice where the rates are neither saturated nor at the null, at the sizes the planned analyses have. Proposed, as two sibling profiles sharing one seed family exactly as here (the anchor reproduction check is then bit-exact):

- **Columns:** `n_samples` ∈ {80, 200, 500} × **three stages** × `layers` ∈ {methylation + expression, all three} — six columns. Three stages are the case-study and planned AD design; the four-stage columns tracked the three-stage ones for magnitude and shape (within 0.13 at every n, within SE at n ≥ 500) and differ for orientation only within the rank-driven scatter, so a second stage count buys little at twice the cost. n = 1000 is dropped: magnitude and shape are saturated there and the pilot's n = 1000 rates stand as pilot-grade values.
- **Effects:** magnitude at e ∈ {0.02, 0.05, 0.1} (0.01 resolves only at n = 1000; 0.25 is saturated at every column); orientation and shape at e ∈ {0.25, 0.5, 1.0} (1.0 kept so that orientation's non-monotone response is measured at paper grade with the selected rank recorded per cell); translation at e = 1.0 as the negative control; a zero-effect anchor per column.
- **Units and cost:** magnitude profile 6 × (1 + 3) × 500 = 12,000 units; geometry profile 6 × (1 + 3 + 3 + 1) × 500 = 24,000; 36,000 in total. From the pilot's unit timings (10–30 s at n ≤ 500 with 199 permutations, the RRPP adding 2–5 s) a 999-permutation unit costs about 25 s at n = 80 to 50 s at n = 500 with three blocks: **about 350–400 core-hours**, in line with the latent-adjustment power comparison, on about 150 shards of ≤ 2.5 h on `512x1024`.
- **Predeclared rules:** `acceptance.anchor_decision` at α 0.05, k 2 (bound 0.0695 at R = 500); `acceptance.type_i` entries for translation's `angle` with `bound: upper` at every column; no power targets — the slice is guidance, not a gate. The report states rates ± SE, each column's stage AUROC and selected rank, and the design-guidance figure.
- **Second priority, if budget allows:** the stage-effect arm δ ∈ {0.25, 0.5, 1.0} at the baseline column (3 columns × 11 cells × 500 = 16,500 units, about 150 core-hours), because that is where power transitions and the CV rank inflates. δ ≤ 0.1 is the null regime at n = 200 and need not be rerun; the rank-floor check already sized it at n = 80 at paper grade.

Not proposed, with reasons: the **ρ arm** (continuity raises shape and e = 1.0 orientation power and lowers magnitude power; it is an assumption about the baseline geometry, and the design-point pilot chose ρ = 0 as the most general claim — a follow-up sensitivity if the case study's trajectories trend); **group ratio, stage proportions and `p_dmp`** (main effects measured at pilot grade move power in the expected direction; the orientation jump at group ratio 0.65 should be understood through the selected rank before it is rerun); the **four-stage** and **n = 1000** columns (above).

**Guidance statement the slice would support.** `delta` is the component a cohort of n = 80–200 can detect (e ≥ 0.1 at n = 80, e ≥ 0.05 at n = 200, with two blocks); `shape` reaches 0.72–0.85 at n = 500 with two blocks (e = 0.5–1.0) and 0.84 at n = 200 with three blocks at e = 1.0, but stays ≤ 0.55 at n = 200 and ≤ 0.20 at n = 80 with two blocks; an `angle` non-rejection is uninformative at n ≤ 200 with two blocks, and any `angle` reading carries its selected rank. None of this holds when the stage signal is at the SEA-AD level: there the pilot measured no power for any component at n = 200.

## 8. Limitations

- **Pilot precision.** 100 replicates: SE ≤ 0.05 on every rate, and the anchor bound (0.0936) is loose; a column whose anchor sits between 0.07 and 0.09 is neither cleared nor flagged at this grade. Translation's `angle` rate is read against the bound and flagged, never claimed.
- **Block count is confounded with feature count.** The two-block columns measure 498 features (methylation + expression); the three-block columns add the protein block. The contrast is "more blocks and more features", not block count alone; block weighting is roadmap item 5.
- **Structural analogue.** The InterSIM blocks are a structural analogue of a real multi-omic design (498 measured features, cross-omic coupling through the reference maps), not an emulation of SEA-AD's RNA/ATAC.
- **Stage signal.** The core columns run at δ = 2.0, where the held-out stage AUROC is saturated; only the δ arm reaches the regime of a real analysis, and only at n = 200 with three stages and two blocks.
- **One-at-a-time arms** measure main effects off one baseline; interactions with sample size are not measured.
- **Orientation's effect-size response is confounded with the CV-selected rank.** The e = 1.0 `angle` rate tracks the rank the stage-supervised double CV retains (§3.2, §3.4), and the rank is recorded, not controlled; the pilot describes the dependence and does not explain it.
- **The δ arm scales every surgery.** Magnitude multiplies δ and the relocation surgeries move sites whose effect is δ, so "no power at δ ≤ 0.1" is partly the construction (group differences proportional to the stage effect), and it is measured at n = 200 with three stages and two blocks only.
- **Two anchor excesses at 0.10 are unresolvable at R = 100**, and the baseline column's second anchor (the Type I-phase draw) is also a 100-replicate reading; the paper-grade slice resolves both columns.
