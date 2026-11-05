# Split interaction test — findings report (2026-11-05)

**Change:** `span-split-interaction-test` (MOTCO 0.18.0). **Reads:** `results/span-split-2026-11-05/` (`span_split_decision.json`, `simulation/` and `plasmode/` split, production, capture and trajectory tables, `accounting.json`, `reproduction.json`, `settings.json`, `plasmode-fidelity/fidelity.json`, `PROVENANCE.txt`, `sacct.txt`). **Concept page:** [measurement visibility](../measurement-visibility.md#testing-outside-the-span-the-split-interaction-test).

The latent-space comparison ([report](latent-space-comparison-2026-11-04.md)) found no measurement space that sees a group difference outside the fitted stage-PLS span without losing power on what the span already sees. This change asks the out-of-span question with its own test instead of a new space. `motco.stats.interaction.interaction_split_test` splits the group × stage interaction energy of the standardized features exactly into three parts:

- `span`: the part inside the fitted span;
- `complement`: the part in its orthogonal complement;
- `total`: all of it.

Each part is tested by a trace pseudo-F under Freedman–Lane permutation of the reduced model's residuals, with RRPP's draws at the analysis's seed, against a basis that never sees the group labels and is never refitted.

**In one paragraph.** **The split test is validated: `split_validated`.** It holds size on three simulated designs and on the real SEA-AD blocks, and its complement part sees the out-of-span rotation that `angle` cannot.

- **Size holds everywhere.** The 18 simulation control rates read 0.032–0.058 against a per-cell bound of 0.0695, with none above it. On real noise at the case study's configuration (six covariates, rank 5, 10,000 features) the three parts read 0.040, 0.050 and 0.012.
- **The complement sees the out-of-span rotation.** A 45° rotation at φ = 90° is rejected in 1.00 (n = 200) and 0.68 (SEA-AD table) of replicates, where `angle` reads 0.08 and 0.10 on the same units.
- **The parts localize.** As the rotation tilts out of the span, the span part's rate falls (1.00 → 0.52 at n = 200) and the complement's rises (0.35 → 1.00).
- **A sparse real-data plant is not seen.** The plasmode's 20-feature plant on 10,000 features is rejected by no part above the null rate (0.02–0.05 at 2 and 3 within-stage SDs), as by `angle` (0.034–0.036). A complement non-rejection at SEA-AD is therefore uninformative about sparse differences of this size.
- **The ratio statistic is what holds size.** The uncorrected interaction SS rejects the plasmode's null complement in 0.268 of replicates; the pseudo-F reads 0.012 on the same draws.

## 1. Configuration and provenance

- **Simulation profile.** The latent-space comparison's units: three designs (`n200`, `seaad`, `seaad_weak` at the calibrated δ = 0.11), fourteen cells and base seed 2035, imported from `scripts/latent_space_comparison.py` unchanged. The weak δ was read from the comparison's committed calibration, with its sha256 checked against the comparison's settings.
  - **Replicates.** The anchor and translation ran 500 replicates, every other cell 100 (the comparison's own units), all at 999 permutations: 6,600 units.
  - **Each unit** fits `pls_cv` once and evaluates it with net displacement, rank floor 2 and span capture. It then runs the split test on the fit's standardized features in the PLS input coordinates against `qr(x_rotations_)`, at the unit's RRPP seed, uncorrected and under `df`.
  - **Rejection** is p < 0.05, the comparison's rule.
- **Plasmode profile.** The SEA-AD plasmode's units: the astrocyte `cohort_30` export through `plasmode_check`'s contract, within-stage label permutation, the 20-feature plant (10 per block, at e within-stage SDs) and its seeds (base seed 2034, `paper` stream).
  - **Cells:** the anchor, e = 2.0 and e = 3.0 at 500 × 999, 1,500 units.
  - **Each unit** is one `analyze_trajectories` at the case study's configuration (six covariates, rank floor 2, net displacement, `df`) with attribution off and `span_split_test` on, plus an uncorrected rerun of the split test for the SS diagnostic.
  - **Rejection** is p ≤ 0.05, the plasmode's rule.
- **Code:** `b84d23c` on `feat/span-split-interaction-test`, in a clean detached worktree on `ing`.
- **Cluster:** partition `512x1024`, arrays 934742 (simulation, 200 shards) and 934743 (plasmode, 150 shards), each throttled to 64 tasks. They ran on n1, n3 and n4 with n11 excluded and BLAS pinned to one thread. The fidelity gate passed on the cluster's export copy (rank 5, AUROC 0.635556; job 934741).
- **Pre-flight:**
  - the fast suite (2,377 passed);
  - a local smoke of both profiles;
  - a timing run at 999 permutations;
  - a 1,200-unit local null check. A 40-replicate pre-check had suggested a liberal complement at CV ranks above K − 1; the larger check did not reproduce it, and the design was kept (notes).

Details are in `PROVENANCE.txt` and `openspec/changes/span-split-interaction-test/notes.md`.

## 2. Accounting and reproduction

- **Units.** 6,600/6,600 simulation and 1,500/1,500 plasmode units complete, 0 failed, 0 duplicates. All 350 array tasks, the fidelity job and the merge job `COMPLETED 0:0` (`sacct.txt`). Every plasmode anchor unit reproduced the fidelity gate's latent space (rank 5, AUROC 0.635556).
- **Cost.** 25.6 + 36.9 = 62.5 recorded core-hours (median 11.1 s per simulation unit and 80.0 s per plasmode unit), against a projection of 113.
- **Reproduction (the decision's gate), both passed:**
  - **Simulation:** `pls_cv`'s count and rejections from the first 199 of each unit's 999 draws equal the comparison's committed `rate_table.csv` in all 126 rows (3 designs × 14 cells × 3 statistics, replicates 0–99).
  - **Plasmode:** the `delta`/`angle`/`shape` rates equal the plasmode's committed `plasmode_summary.csv` in all 9 rows. The anchor reads 0.042 / 0.048 / 0.026, e = 2.0 reads 0.030 / 0.034 / 0.020 and e = 3.0 reads 0.026 / 0.036 / 0.016.
  - **So every rate below pairs** with the comparison's five arms and with the plasmode on the same datasets, fits and draws.

## 3. The decision

`span_split_decision.json` (rule version 1) returns **`split_validated`**:

| criterion | rule | result |
|---|---|---|
| C1 simulation size | the three parts' pseudo-F at the anchor and translation of every design, m = 18, count-plus-ceiling (allowed 2 exceedances of the bound 0.0695, ceiling 0.078) | **holds**: 0 exceedances; rates 0.032–0.058 |
| C2 real-noise size | the three parts at the plasmode anchor, m = 3 (allowed 1, ceiling 0.072) | **holds**: total 0.040, span 0.050, complement 0.012 |
| C3 out-of-span gain | complement rate for orientation e = 0.5 at φ = 90° (dense) above `pls_cv`'s `angle` by more than 2 pooled SE, at `n200` and `seaad` | **holds**: 1.00 vs 0.08 (difference 0.92, threshold 0.054); 0.68 vs 0.10 (0.58, threshold 0.111) |

C1 rates (pseudo-F, 500 replicates each):

| design | anchor total / span / complement | translation total / span / complement |
|---|---|---|
| `n200` | 0.048 / 0.048 / 0.032 | 0.048 / 0.046 / 0.050 |
| `seaad` | 0.050 / 0.048 / 0.036 | 0.048 / 0.042 / 0.052 |
| `seaad_weak` | 0.036 / 0.040 / 0.050 | 0.034 / 0.036 / 0.058 |

**Consequence:** a case study may report the split test, reading its parts as the measurement-visibility page states.

## 4. Localization across alignment and modes

Pseudo-F rejection rate per part, total / span / complement (`split_table.csv`; 100 replicates per cell, SE ≤ 0.05):

| cell | `n200` | `seaad` |
|---|---|---|
| magnitude e = 0.05 | 0.02 / 0.56 / 0.01 | 0.04 / 0.09 / 0.05 |
| orientation e = 0.25, φ = 0 / 30 / 60 / 90 | 0.80 / 0.93 / 0.11; 0.63 / 0.81 / 0.22; 0.29 / 0.83 / 0.29; 0.19 / 0.02 / 0.29 | 0.12 / 0.45 / 0.07; 0.09 / 0.43 / 0.15; 0.08 / 0.49 / 0.27; 0.07 / 0.04 / 0.31 |
| orientation e = 0.5, φ = 0 / 30 / 60 / 90 | 1.00 / 1.00 / 0.35; 1.00 / 1.00 / 0.78; 1.00 / 0.93 / 1.00; 1.00 / 0.52 / 1.00 | 0.99 / 0.99 / 0.32; 0.94 / 0.75 / 0.78; 0.66 / 0.41 / 0.72; 0.36 / 0.19 / 0.68 |
| orientation e = 0.5, φ = 90, 20 features | 0.90 / 0.79 / 0.78 | 0.16 / 0.36 / 0.58 |
| shape e = 0.25 / 0.5 | 0.05 / 0.80 / 0.04; 0.16 / 0.72 / 0.06 | 0.04 / 0.39 / 0.07; 0.06 / 0.43 / 0.07 |

- **The parts follow the alignment.** For the 45° rotation the span part falls from 1.00 to 0.52 (n = 200) and from 0.99 to 0.19 (SEA-AD table) as φ goes from 0° to 90°, and the complement rises from 0.35 to 1.00 and from 0.32 to 0.68. Magnitude and shape plants, which stay in the stage program, are rejected by the span part and not by the complement (0.01–0.07).
- **The complement also sees in-span plants, partly.** At φ = 0 the 45° rotation is rejected by the complement in 0.35 / 0.32 of replicates. The plant keeps 0.97 / 0.93 of its energy in the fitted span (`capture_table.csv`), and the 3–7 % outside is enough at this size. Through the sample scaler, its rigid rotation also acquires off-span components (the unequal-scales route on the [visibility page](../measurement-visibility.md)). A complement rejection therefore means energy outside the *fitted* span, not a difference outside the stage program (design D10).
- **The span part also sees the out-of-span rotation at n = 200** (0.52 at φ = 90). The production statistics register it there too: `shape` 0.51 and `delta` 0.20 on the same units, against `angle` 0.08. This is the sample-scaler route the comparison reported (§5 of its report).
- **Sparse support.** At φ = 90 a 20-feature direction keeps 0.11–0.12 of its energy in the span, against 0.01–0.03 for the dense one. Its complement rate drops from 1.00 to 0.78 (n = 200) and from 0.68 to 0.58 (SEA-AD table), and its span rate rises (0.52 → 0.79 and 0.19 → 0.36).
- **The total part is the weakest for in-span plants** (magnitude 0.02, shape 0.05–0.16 at n = 200), because their energy is diluted over all 498 coordinates. The split is what recovers the power: the span part reads 0.56–0.80 on the same plants.
- **At SEA-AD's stage signal (`seaad_weak`) no part has power** for any plant (0.03–0.08), as for every arm of the comparison.

## 5. The in-span share against capture

Median observed in-span share of the interaction energy beside the median population capture of the planted difference (`capture_table.csv`):

| cell | `n200` share / capture | `seaad` share / capture | `seaad_weak` share / capture |
|---|---|---|---|
| anchor (null) | 0.006 / — | 0.062 / — | 0.292 / — |
| translation (null) | 0.171 / — | 0.215 / — | 0.288 / — |
| orientation e = 0.5, φ = 0 | 0.832 / 0.971 | 0.673 / 0.927 | 0.291 / 0.082 |
| orientation e = 0.5, φ = 90 | 0.008 / 0.008 | 0.050 / 0.026 | 0.297 / 0.005 |
| shape e = 0.5 | 0.325 / 0.968 | 0.235 / 0.932 | 0.312 / 0.095 |

- **At a strong stage signal the share tracks capture** for large plants: 0.83 against 0.97 for the in-span rotation, 0.008 against 0.008 for the out-of-span one at n = 200. Smaller plants sit between their capture and the null share, because the noise's interaction energy is in the same sum.
- **The null share is not zero, and it depends on the fitted rank.** The null cells' median shares are 0.006 (`n200` anchor, median rank 2), 0.062 (`seaad` anchor), 0.17–0.22 (the translation cells, median rank 3) and 0.29 at the weak design; at the plasmode anchor (rank 5) it is 0.236 (`plasmode/split_table.csv`). The local null check split this by rank: a median 0.007 at rank 2 = K − 1 against 0.34 at higher ranks (notes).
- **Why the share is not zero:** components beyond K − 1 follow the total covariance, where the noise's interaction energy also concentrates. The local null check found 0.01–0.03 for a basis of the pooled stage means (notes).
- **How to read it:** an observed share must be read against the null share at the same rank and design, never against 0. At the weak design the share is the null share whatever the plant: the fitted span is mostly noise (pooled-stage capture 0.07–0.10).

## 6. The span part against the production statistics

On the units where both apply (`production_table.csv`), the span part is an omnibus test and the production statistics are targeted:

| cell | `n200` span vs target | `seaad` span vs target |
|---|---|---|
| magnitude e = 0.05 (`delta`) | 0.56 vs 0.74 | 0.09 vs 0.21 |
| orientation e = 0.25, φ = 0 (`angle`) | 0.93 vs 0.87 | 0.45 vs 0.47 |
| shape e = 0.25 / 0.5 (`shape`) | 0.80 / 0.72 vs 1.00 / 1.00 | 0.39 / 0.43 vs 0.70 / 0.94 |

The span part matches `angle` for in-span rotations and trails `delta` and `shape`. It does not replace the statistics: it says whether the groups differ inside the span, while `delta`, `angle` and `shape` say how.

## 7. The complement against the comparison's arms

The comparison's committed `angle` rates on the same datasets (100 × 199; `latent-space-comparison-2026-11-04/rate_table.csv`) beside the complement part:

| cell | `n200` complement / `concat` / `pls_cv` | `seaad` complement / `concat` / `pls_cv` | `seaad_weak` complement / `concat` / `pls_cv` |
|---|---|---|---|
| anchor | 0.03 / 0.05 / 0.04 | 0.04 / 0.06 / 0.04 | 0.05 / **0.21** / 0.02 |
| translation | 0.05 / 0.05 / 0.04 | 0.05 / 0.06 / 0.06 | 0.06 / **0.22** / 0.02 |
| orientation e = 0.5, φ = 90 | 1.00 / 1.00 / 0.08 | 0.68 / 0.90 / 0.11 | 0.07 / 0.23 / 0.01 |
| orientation e = 0.5, φ = 90, 20 features | 0.78 / 0.98 / 0.06 | 0.58 / 0.49 / 0.09 | 0.06 / 0.26 / 0.01 |

- **The complement reaches most of `concat`'s out-of-span power and holds size where `concat` did not.** At the weak design `concat`'s `angle` is liberal (0.21–0.26 in every cell); the complement reads 0.05–0.08 there.
- **On the SEA-AD table `concat` is more powerful for the dense rotation** (0.90 against 0.68). For the sparse one the complement is (0.58 against 0.49).
- **Different questions.** `concat`'s `angle` tests orientation in all coordinates; the complement tests any group × stage difference outside the fitted span, with no geometry.

## 8. Real noise: the plasmode

| cell | total | span | complement | `delta` / `angle` / `shape` |
|---|---|---|---|---|
| anchor | 0.040 | 0.050 | 0.012 | 0.042 / 0.048 / 0.026 |
| e = 2.0 | 0.042 | 0.046 | 0.022 | 0.030 / 0.034 / 0.020 |
| e = 3.0 | 0.042 | 0.048 | 0.018 | 0.026 / 0.036 / 0.016 |

- **Size holds on real noise** at the case study's configuration. The complement is conservative there (6/500), in 9,995 dimensions with six covariates.
- **There is no power for the sparse plant.** A shift of 2–3 within-stage SDs on 20 of 10,000 features is rejected by no part above the null rate. The plant keeps a median 2.2 % of its norm in the rank-5 span, the random-direction share (`results/plasmode-2026-11-02/diagnostic/latent_capture.json`), so it is an out-of-span difference. Its energy is a small fraction of the complement's 9,995 coordinates of noise, and an energy test does not find it. This is design risk R1, and thread 3 (a sparse per-feature scan) is the targeted follow-up.
- **For the case study:** a complement non-rejection on SEA-AD says nothing about sparse differences of this size. It is informative only about differences with substantial energy outside the span, like the simulation's dense rotation.

## 9. The SS diagnostic

Null rejection rates of the uncorrected interaction SS (`ss`), the interaction SS under `df` (`ss_df`) and the pseudo-F (`f`), at the controls:

| profile, cell | part | `f` | `ss` | `ss_df` |
|---|---|---|---|---|
| simulation `seaad`, translation | complement | 0.052 | 0.068 | 0.040 |
| simulation `seaad_weak`, anchor | complement | 0.050 | 0.082 | 0.038 |
| simulation `seaad_weak`, translation | complement | 0.058 | 0.098 | 0.048 |
| plasmode anchor | total | 0.040 | 0.114 | 0.040 |
| plasmode anchor | span | 0.050 | 0.074 | 0.048 |
| plasmode anchor | complement | 0.012 | **0.268** | 0.006 |

D1 predicted this pattern:
- **Why the uncorrected SS is liberal.** The Freedman–Lane residuals are too narrow by √((n − r)/n), and an energy summed over many coordinates concentrates. So the uncorrected SS is liberal in proportion to dimension: 0.268 in the plasmode's 9,995-dimensional complement, 0.068–0.098 in the simulation's 495-dimensional one, and near α in the span.
- **What corrects it.** Rescaling by κ removes the excess (0.006–0.048), and the pseudo-F never had it. This bears on why the analysis's RRPP needs the `df` correction at the case study's configuration. It does not test why `concat`'s `angle` was liberal at the weak design, which is a different statistic.

## 10. What follows

- **The split test may be used by a case study** under `split_validated`, with the reading on the [measurement visibility](../measurement-visibility.md) page:
  - span = the fitted span;
  - a complement rejection = energy outside it, including the fitted span's misses of the true program;
  - a non-rejection uninformative without power at the design;
  - the share read against the null share.
- **Applying it to SEA-AD** is `sead_sexdiff`'s change. It pins MOTCO 0.18.0 and fixes its multiplicity policy over the three parts before any real-label analysis (design open question).
- **Thread 2** (a group-aware space refitted per permutation) now has a benchmark to beat: the complement's 1.00 / 0.68 for the dense out-of-span rotation at size.
- **Thread 3** (a sparse per-feature interaction scan) is the follow-up the plasmode result asks for.

## 11. Limitations

- **The power cells are at R = 100** (SE ≤ 0.05). C3's margins (0.92 and 0.58 against thresholds of 0.054 and 0.111) do not depend on that.
- **Stage-span plants are relative to the configuration.** At the weak design they are as weak as the stage program, so its zero power is descriptive. Absolute-scale plants on real noise are the plasmode's, and they are sparse.
- **The plasmode measures one kind of difference** (20 features, a monotone stage weight). Dense or moderately sparse real-data differences were not planted.
- **The production `angle` at the `seaad` design** reads 0.072 (anchor) and 0.070 (translation) at 500 replicates, just above its per-cell bound of 0.0695. This change does not judge that statistic. The two cells share their null datasets, and the comparison read 0.04 / 0.06 there at R = 100.
- **The null in-span share depends on the fitted rank** (§5). The share is a descriptive location measure, not a test, and it needs the null reference to be read.
