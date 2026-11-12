# Predictive-subspace coordinates — findings report (2026-11-12)

**Change:** `predictive-subspace-coordinates` (MOTCO 0.24.0). **Reads:** `results/predictive-subspace-2026-11-12/` (`predictive_subspace_decision.json`; `simulation/`, `weak_size/` and `plasmode/` rate, paired and cell tables, `accounting.json`, `reproduction.json`, `settings.json`; `PROVENANCE.txt`, `sacct.txt`). **Concept page:** [measurement visibility](../measurement-visibility.md#the-measured-dimension-fit-rank-against-predictive-subspace).

Production uses the CV-selected rank A of the stage-PLS fit for two jobs: it regularizes the fit, and it is the number of dimensions `delta`, `angle` and `shape` are measured in. The stage labels are one-hot, so the fit's Y loadings have rank K − 1 whatever A is, and the fitted stage labels move only along K − 1 directions of the score space: the **predictive subspace**. MOTCO 0.24.0 measures there on request, in three coordinate systems on one subspace:

- **`predictive`:** the scores projected on it;
- **`predictive_orthonormal`:** the input coordinates on an orthonormal basis of it;
- **`predictive_fitted`:** the fitted stage indicators.

This pilot read the three against production on matched units under a decision fixed in advance (`design.md` D7).

**In one paragraph.** **The verdict is `keep_production`. All three coordinate systems read `size_excess` on the SEA-AD plasmode's real noise.**

- **Real noise.** On 500 permuted-label replicates of the astrocyte export, production rejects `delta` in 21 and `shape` in 13. Each predictive arm rejects `delta` in 65–68 (0.13–0.136) and `shape` in 38–39 (0.076–0.078), above the family's ceiling of 0.072. The excess is paired: about 50 units reject only under the candidate, against 5–6 that reject only under production.
- **Simulation.** Size holds at every design, the weak one included (500 replicates × 999 permutations).
- **Power.** The arms lose no protected power anywhere and gain in 10 of the 18 δ = 2.0 target cells at the Bonferroni threshold. In the stage-effect band, between δ = 2.0 and the weak design, the gains are large: the n = 200 rotation at δ = 0.5 and e = 0.1 goes from 0.06 to 0.94. On the plasmode's real-noise plant `delta` goes from 0.51 to 0.79–0.82 at e = 0.5.
- **Specificity separates the three.** `predictive` raises production's rotation-into-`shape` leak in the band. `predictive_fitted` leaks at δ = 2.0. `predictive_orthonormal` passes both families.
- **So** the subspace's gains are real on simulated noise, and on real noise it does not hold size.

## 1. Configuration and provenance

- **The unit.** Each unit builds one dataset, fits the production stage-PLS space once (double-CV rank, floor 2, net displacement) and evaluates these arms on it, serially, at the unit's RRPP seed:
  - **`pls_cv`:** production, in the fit's scores;
  - **`predictive`, `predictive_orthonormal`, `predictive_fitted`:** the three coordinate systems on the predictive subspace of that same fit, built without refitting;
  - **`forced_k_minus_1`:** the fit forced to K − 1 = 2 components, the control that lowers the fit rank instead of the measured dimension;
  - **`predictive_rank_plus_one`:** `predictive` on the fit forced to A + 1 components (simulation only).

  The arms share data, fit and permutation indices, so the decision compares them by exact paired tests on the discordant units. Every unit also checks that the predictive arms measure the production fit (five checks, `design.md` D5).
- **`simulation` profile.** A fresh seed family (base seed 2036; the probe spent 4242 and the earlier runs 2035) on the latent-space comparison's generators and seed derivations, under the planted scale, at 100 replicates × 199 permutations, rejection p < 0.05:
  - **Designs:** `n200` and `seaad` (the SEA-AD table) at δ = 2.0; the band `n200_delta1`, `n200_delta0.5`, `n200_delta0.25` and `seaad_delta0.5`; `seaad_weak` at the committed δ = 0.11.
  - **Cells (11 per design):** the anchor; translation e = 1.0; magnitude e ∈ {0.05, 0.1, 0.25}; orientation e ∈ {0.1, 0.25, 0.5} (in units of 90°, φ = 0); shape e ∈ {0.1, 0.25, 0.5}.
  - **The surgery cell:** the default feature-set orientation surgery at e = 1.0 on the SEA-AD table at δ = 2.0 (descriptive).
  - **Units:** 7,800.
- **`weak_size` profile.** The comparison's `seaad_weak` anchor and translation at replicates 0–499 × 999 permutations (base seed 2035): the span-split check's own simulation units. 1,000 units.
- **`plasmode` profile.** The SEA-AD astrocyte `cohort_30` export (RNA and ATAC, 5,000 features each, 80 donors), read through the plasmode's contract, with sex labels permuted within stage by the plasmode's `paper` stream (base seed 2034). Every arm runs `analyze_trajectories` at the case-study configuration: six covariates, floor 2, net displacement, `df`, 999 permutations, rejection p ≤ 0.05.
  - **Cells:** the anchor at replicates 0–499, and the production selection's real-noise span plant at e ∈ {0.25, 0.5, 1.0} at replicates 0–99.
  - **The predictive arms** run at `forced_components` = the unit's `pls_cv` rank. This is the amendment approved before any array unit: the same final fit, without repeating the CV.
  - **Units:** 800.
- **Code.** Revision `fea0874` on `proposal/predictive-subspace-coordinates` (MOTCO 0.24.0), run from a clean detached worktree on `ing`. Python 3.11.16, numpy 2.3.5, scipy 1.16.3, scikit-learn 1.8.0 and pandas 2.3.3, with BLAS pinned to one thread.
- **Cluster.** Partition `512x1024`, n4 and n11 excluded. Arrays 938116 (`plasmode`, 80 × 10 units), 938117 (`simulation`, 120 × 65) and 938118 (`weak_size`, 25 × 40).
- **Pre-flight.** A workstation smoke of the three profiles (168 units, 0 failed, 0 inconsistent) projected about 178 core-hours as designed, or about 91 with the amendment. The user approved the amendment and a cap of 120 (`openspec/changes/archive/2026-11-12-predictive-subspace-coordinates/notes.md`).

Details are in `PROVENANCE.txt`.

## 2. Accounting and reproduction

- **Units.** 7,800 + 1,000 + 800 of 9,600: 0 failed units, 0 failed arms, 0 inconsistent units and 0 duplicates. All 225 array tasks ended `COMPLETED` (`sacct.txt`).
- **Consistency.** On every unit:
  - the Y loadings had rank min(A, K − 1);
  - the fitted stage labels lay in every predictive arm's column space to 1e-8;
  - every predictive arm fitted production's model;
  - the three predictive arms reported equal span-capture shares;
  - at A ≤ K − 1, `predictive` reproduced `pls_cv`'s p-values.
- **Cost.** 46.3 + 6.1 + 27.4 = **79.8 recorded core-hours**, against the projection of 91 and the cap of 120.
  - The simulation unit's median was 20.6 s, of which the production CV took 17.9 s and each arm 0.4 s.
  - The plasmode unit's median was 119.6 s: `pls_cv` 116.9 s, each other arm about 0.6 s.
- **Reproduction (the decision's gate): passed.**
  - **`weak_size`:** `pls_cv`'s count and rejections equal the span-split check's committed `production_table.csv` in all 6 rows.
  - **`plasmode` anchor:** `pls_cv` rejects `delta`, `angle` and `shape` in 21, 24 and 13 of 500 replicates, the plasmode's committed rates (0.042, 0.048, 0.026).
  - **`plasmode` plant:** all 9 rows equal the production selection's committed `pls_cv` rows.
  - **Fidelity:** all 500 anchor units selected rank 5 at held-out AUROC 0.635556.
  - **So** every production rate on the real-noise cells is a committed rate, and every candidate rate pairs with it on the same data and draws.

## 3. The decision

`predictive_subspace_decision.json` (rule version 1) returns **`keep_production`**: no arm qualifies.

| criterion | rule | `predictive` | `predictive_orthonormal` | `predictive_fitted` |
|---|---|---|---|---|
| gate | reproductions, fidelity, no inconsistent unit | **passed** | | |
| S1 size, simulation | anchor and translation, seven designs, m = 42: at most 3 above the bound 0.0936, none above the ceiling 0.13 | holds: 2 above | holds: 1 above | holds: 3 above |
| S1 size, `weak_size` | m = 6 at R = 500: at most 1 above 0.0695 | holds: 0.030–0.052 | holds: 0.032–0.042 | holds: 0.034–0.056 |
| S1 size, plasmode | the anchor's three statistics at R = 500, m = 3: at most 1 above 0.0695, none above the ceiling 0.072 | **fails**: `delta` 0.136, `shape` 0.076 | **fails**: `delta` 0.134, `shape` 0.076 | **fails**: `delta` 0.130, `shape` 0.078 |
| S2 protected power | 57 target cells (six designs at δ ≥ 0.25, the plasmode plant), one-sided exact paired loss: at most 4 below 0.025, none below 0.05/57 | holds: 0 losses | holds: 0 losses | holds: 0 losses |
| S3 specificity, δ = 2.0 | off-target rates, m = 36, the S1 rule (bound 0.0936, ceiling 0.13) | holds: 2 above the bound | holds: 0 above | **fails**: 4 above the ceiling (0.15–0.24) |
| S3 specificity, band | off-target excess over `pls_cv`, m = 72: at most 4 one-sided exact paired p below 0.025, none below 0.05/72 | **fails**: 9 below 0.025, 2 below the ceiling | holds: 0 | holds: 0 |
| G gain | 18 target cells at δ = 2.0, one exact paired gain p ≤ 0.05/18 | holds: 10 cells | holds: 10 cells | holds: 10 cells |
| **arm verdict** | `size_excess` → `power_loss` → `specificity_excess` → `no_gain` → `qualifies` | **`size_excess`** | **`size_excess`** | **`size_excess`** |

**Consequence** (the decision's text): "Production keeps measuring in the fitted span's scores; the three predictive coordinate systems stay opt-in diagnostics, and the measurement-visibility page records what the measured dimension costs and why it was kept."

`predictive_orthonormal` fails only S1, and only on the plasmode.

## 4. Size on real noise: the plasmode anchor

Rejections over the 500 anchor replicates (`plasmode/rate_table.csv`), and the discordant counts against `pls_cv` (b only `pls_cv` rejects, c only the arm; `plasmode/paired_table.csv`):

| arm | `delta` | `angle` | `shape` | `delta` b/c | `shape` b/c |
|---|---|---|---|---|---|
| `pls_cv` | 21 (0.042) | 24 (0.048) | 13 (0.026) | — | — |
| `predictive` | 68 (0.136) | 11 (0.022) | 38 (0.076) | 6 / 53 | 5 / 30 |
| `predictive_orthonormal` | 67 (0.134) | 9 (0.018) | 38 (0.076) | 6 / 52 | 4 / 29 |
| `predictive_fitted` | 65 (0.130) | 10 (0.020) | 39 (0.078) | 5 / 49 | 4 / 30 |
| `forced_k_minus_1` | 7 (0.014) | 10 (0.020) | 5 (0.010) | 20 / 6 | 12 / 4 |

- **The excess is the subspace's, not a metric's.** The three coordinate systems agree to within three rejections in every statistic, while production holds and the forced rank-2 fit is conservative. The two-sided paired p-values for `delta` are below 1e-9 in all three.
- **The fit is fixed across the anchor's replicates.** The latent space is label-free and the anchor does not move the features, so every anchor unit fits the same rank-5 model (AUROC 0.635556, in-sample stage R² 0.840). The predictive arms measure one fixed two-dimensional subspace of it. The permutation null that holds size in five dimensions does not hold it in those two: the anchor's median `angle` null 95th percentile narrows from 40.2° to 11.0–12.1°, and the observed `angle` from 20.7° to 3.2–3.5° (`plasmode/cell_table.csv`).
- **Simulation does not reproduce it.** The same arms hold size at every simulated design, including `seaad_weak` at R = 500 (0.030–0.056). The probe's `shape` reading there (0.077 over 300 units) did not recur: `predictive` reads 0.052 at the anchor, against production's 0.050. What the plasmode has and the generator lacks is not identified by this run. Candidates include real noise structure aligned with the stage directions, and the six-covariate design with its `df` correction. A diagnosis would be its own change.
- **`angle` moves the other way.** On the plasmode every predictive arm rejects `angle` less often than production (b/c 24/11 for `predictive`).

## 5. Power: no loss, gains where the span sees the plant

S2 holds with no loss in any arm or cell. Target rates at δ = 2.0, `pls_cv` → `predictive` (b/c against `pls_cv`; the other two candidates agree within 0.04):

| design | magnitude e = 0.05 | orientation e = 0.1 | orientation e = 0.25 | shape e = 0.1 | shape e = 0.25 |
|---|---|---|---|---|---|
| n = 200 | 0.56 → 0.72 (0/16) | 0.76 → 0.99 (0/23) | 0.90 → 1.00 (0/10) | 0.71 → 0.83 (0/12) | 0.97 → 1.00 (0/3) |
| SEA-AD table | 0.22 → 0.54 (0/32) | 0.40 → 0.90 (0/50) | 0.41 → 1.00 (0/59) | 0.13 → 0.56 (0/43) | 0.63 → 0.90 (0/27) |

Each candidate passes G in the same 10 cells, with gain p from 1.7e-18 (the SEA-AD table's 22.5° rotation) to 9.8e-4. Magnitude e = 0.1 on the SEA-AD table also gains (0.57 → 0.86).

**The band is where the payoff is largest.** It is recorded, and did not decide G. Production's CV selects ranks 4 (`n200_delta0.5`), 11 (`n200_delta0.25`, q10–q90 8–16) and 5 (`seaad_delta0.5`) at a held-out AUROC of 0.998–1.000. The predictive arms measure in two of those dimensions.

| design | cell | `pls_cv` | `predictive` | `forced_k_minus_1` |
|---|---|---|---|---|
| n = 200, δ = 1.0 | orientation e = 0.1 | 0.10 | 1.00 | 0.60 |
| n = 200, δ = 0.5 | orientation e = 0.1 | 0.06 | 0.94 | 0.11 |
| n = 200, δ = 0.5 | shape e = 0.25 | 0.22 | 0.92 | 0.09 |
| n = 200, δ = 0.25 | magnitude e = 0.25 | 0.35 | 1.00 | 0.20 |
| n = 200, δ = 0.25 | orientation e = 0.25 | 0.12 | 1.00 | 0.11 |
| SEA-AD table, δ = 0.5 | orientation e = 0.25 | 0.13 | 0.98 | 0.15 |
| SEA-AD table, δ = 0.5 | shape e = 0.5 | 0.16 | 0.82 | 0.08 |

Each candidate gains at the Bonferroni threshold 0.05/36 in 29–30 of the 36 band target cells. Pooled over them the paired difference is +0.44 to +0.45, from b/c of about 50 against 1,618–1,660.

**Real-noise plant** (`delta`, `pls_cv` → candidates): e = 0.25, 0.24 → 0.40–0.42; e = 0.5, 0.51 → 0.79–0.82; e = 1.0, 1.00 in every arm. Two of the three cells gain at 0.05/3.

**The weak design does not move.** Production reads 0.01–0.05 at its nine target cells, and the candidates 0.01–0.07; no gain p is below 0.015. There the fitted span keeps 0.107 of the pooled stage program, and the predictive subspace the same 0.107.

## 6. Specificity: the metric decides the leak

Production's own off-target rates in the band are high under the planted scale, and the three coordinate systems move them in opposite directions:

| design | cell | off-target | `pls_cv` | `predictive` | `predictive_orthonormal` | `predictive_fitted` |
|---|---|---|---|---|---|---|
| n = 200, δ = 1.0 | orientation e = 0.5 | `shape` | 0.56 | 0.37 | 0.03 | 0.00 |
| n = 200, δ = 0.5 | orientation e = 0.25 | `shape` | 0.57 | 0.81 | 0.04 | 0.03 |
| n = 200, δ = 0.5 | orientation e = 0.5 | `shape` | 0.80 | 0.85 | 0.41 | 0.01 |
| n = 200, δ = 0.25 | orientation e = 0.5 | `shape` | 0.55 | 0.71 | 0.16 | 0.00 |
| SEA-AD table, δ = 0.5 | orientation e = 0.5 | `shape` | 0.46 | 0.59 | 0.24 | 0.07 |
| n = 200, δ = 0.5 | shape e = 0.5 | `delta` | 0.06 | 0.19 | 0.04 | 0.06 |

(`descriptive.band_off_target_rates`.)

- **`predictive` keeps the score metric's shear** and adds to production's leak: 9 of the 72 band readings have an excess p below 0.025, and the two `shape` readings of the n = 200 rotations at δ = 0.5 (e = 0.1 and 0.25) fall below the ceiling 0.05/72.
- **`predictive_orthonormal` and `predictive_fitted` remove most of production's leak.** Neither has a band excess.
- **`predictive_fitted` leaks at δ = 2.0.** In the response metric, rotations at n = 200 read `shape` 0.19 (e = 0.1) and 0.16 (e = 0.25), and bends read `delta` 0.15 (e = 0.25) and 0.24 (e = 0.5). All four are above the ceiling 0.13.
- **`predictive_orthonormal`** has 4 off-target readings above the bound over both families (the fewest), and would have been named had it held size.

## 7. The control arm: lowering the fit rank

`forced_k_minus_1` truncates the fit at two components instead of measuring two dimensions of the CV fit (`descriptive.control_arm`).

- **It keeps less of the stage program in the band.** Its pooled-stage capture is 0.48 (`n200_delta0.5`), 0.34 (`n200_delta0.25`) and 0.39 (`seaad_delta0.5`). Production's is 0.88, 0.74 and 0.73, and the predictive subspace's 0.87, 0.67 and 0.73.
- **Its power follows.** In the band it rejects at or below production in most target cells, for example 0.11 against production's 0.06 and `predictive`'s 0.94 at the n = 200, δ = 0.5 rotation e = 0.1.
- **At δ = 2.0 it recovers part of the gain.** On the SEA-AD table's 22.5° rotation it reads 0.92, against production's 0.41 and `predictive`'s 1.00, where the CV rank is 2–3 and two components hold the program.
- **On the plasmode it is conservative:** 0.014 / 0.020 / 0.010.
- **So** the later components aim the K − 1 directions, as the probe found. A lower fit rank is the wrong lever, and the predictive subspace keeps what they contribute.

## 8. Rank plus one: the CV rank sits next to the leak

One more fit component under the `predictive` coordinates (`descriptive.rank_sensitivity`):

- **At δ = 2.0 the rotation leaks into `shape`.** At e = 0.5 the off-target `shape` rate goes from 0.03 to 0.92 at n = 200 and from 0.04 to 0.99 on the SEA-AD table; at e = 0.25, from 0.04 to 0.22 and from 0.05 to 0.56. The two-sided paired p-values are 1e-27 to 8e-6.
- **Target power stays high** (0.84–1.00 at δ = 2.0).
- **At n = 200 and δ = 1.0**, where the CV rank is 3, the e = 0.5 rotation's leak still grows, from 0.37 to 0.77.
- **Where the CV rank is already 4–11** (δ = 0.5 and 0.25) one more component changes little. The paired p-values there run from 0.14 to 1.0, except 0.019 for the n = 200, δ = 0.5 rotation at e = 0.1.

The CV's rank at δ = 2.0 is one component from a large specificity loss in the predictive coordinates. Any adoption would need a rank rule for a fixed measured dimension.

## 9. Estimation: observed against population through the same map

Median observed `angle` beside the population's through each arm's map, in degrees (`descriptive.orientation_estimation`):

| design | rotation | `pls_cv` | `predictive` | `predictive_orthonormal` |
|---|---|---|---|---|
| n = 200, δ = 2.0 | 45° | 44.6 / 45.3 | 43.8 / 45.2 | 43.4 / 45.0 |
| SEA-AD table, δ = 2.0 | 45° | 42.9 / 43.1 | 39.7 / 43.4 | 41.3 / 45.2 |
| n = 200, δ = 0.5 | 45° | 39.2 / 42.3 | 30.9 / 44.8 | 32.4 / 44.8 |
| n = 200, δ = 0.25 | 45° | 31.3 / 41.6 | 17.8 / 43.4 | 19.5 / 44.8 |
| SEA-AD table, δ = 0.5 | 45° | 38.1 / 39.1 | 18.6 / 40.8 | 21.3 / 44.8 |

- **At δ = 2.0 the predictive arms estimate the angle within about 4°.**
- **In the band they detect the rotation but report about half of it.** At n = 200 and δ = 0.25 they read 18–20° where the population reads 43–45° through the same map. Production reads 31°, but in more noise dimensions.
- **So a predictive arm is a detector first.** Measuring the difference stays with roadmap open item 3.

## 10. The surgery cell: what the lower dimension gives up

The Phase 6 feature-set orientation surgery at e = 1.0 on the SEA-AD table is largely outside the stage span (`descriptive.surgery_cell`):

| arm | `delta` | `angle` | `shape` | `group_difference` capture |
|---|---|---|---|---|
| `pls_cv` | 0.26 | 0.66 | 0.48 | 0.203 |
| `predictive` | 0.49 | 0.48 | 0.61 | 0.002 |
| `predictive_orthonormal` | 0.52 | 0.49 | 0.63 | 0.002 |
| `predictive_fitted` | 0.51 | 0.48 | 0.59 | 0.002 |
| `forced_k_minus_1` | 0.46 | 0.49 | 0.57 | 0.014 |
| `predictive_rank_plus_one` | 0.32 | 0.22 | 0.30 | 0.001 |

- **Production reads the surgery's `angle` through its extra components:** a median rank of 4 keeps 0.20 of the group difference, against 0.002 in the predictive subspace.
- **The predictive arms read it at 0.48–0.49 instead of 0.66.** Their `delta` and `shape` rise instead, which suggests the plant reaches them indirectly through the stage directions.
- **So** a difference seen only through the dropped components is no longer read as it was.

## 11. Rank and AUROC per design

From `simulation/cell_table.csv`, the anchor's `pls_cv` row:

| design | median rank (q10–q90) | held-out AUROC | in-sample stage R² | `pooled_stage` capture |
|---|---|---|---|---|
| n = 200, δ = 2.0 | 2 (2–3) | 1.000 | 0.983 | 0.970 |
| SEA-AD table, δ = 2.0 | 3 (2–3) | 1.000 | 0.991 | 0.933 |
| n = 200, δ = 1.0 | 3 (3–3) | 1.000 | 0.981 | 0.951 |
| n = 200, δ = 0.5 | 4 (3–4) | 1.000 | 0.943 | 0.876 |
| n = 200, δ = 0.25 | 11 (8–16.1) | 0.998 | 0.970 | 0.737 |
| SEA-AD table, δ = 0.5 | 5 (4–6) | 1.000 | 0.964 | 0.734 |
| SEA-AD table, δ = 0.11 | 2 (2–5) | 0.593 | 0.444 | 0.107 |

The band is where the held-out AUROC is at its ceiling while the rank grows. There the predictive subspace measures two of 4–11 dimensions and gains most.

## 12. What follows

- **Production and `motco analyze` keep measuring in the fitted span's scores.** The paper grade, the default flip and the CLI flag were owed only to `propose_paper_grade`. None is opened.
- **The three values stay opt-in diagnostics** in the harness and `motco.analysis`, now measured:
  - They hold size in simulation, lose no protected power, and gain large power wherever the span sees the plant.
  - They do not hold size on the SEA-AD plasmode's real noise: `delta` reads 0.13 and `shape` 0.08, where production reads 0.042 and 0.026.
  - A case study should not report them as tests.
- **`predictive_orthonormal` is the coordinate system to keep in view.** It fails only the real-noise size, and it removes most of production's own off-target leak in the band. If the plasmode excess is understood and fixed, it is the arm to re-read, with a rank rule for a fixed measured dimension (§8).
- **The real-noise excess is now the open question.** The simulation does not reproduce it, so the generator lacks whatever drives it. Roadmap open item 2, the latent plasmode, measures what the span holds on real noise and is the natural place to read it.

## 13. Limitations

- **R = 100 in simulation.** Every rate there has an SE of about 0.02 at α; the bound is 0.0936. The real-noise and weak-design size readings are at R = 500.
- **One real data set.** The size excess was read on one export, at one configuration (six covariates, `df`), with a fit that is the same in every anchor unit.
- **The planted scale.** The band's off-target rates, production's included, are read under the planted scale. The probe ran the null scale, so its band rates are not comparable.
- **Three stages only.** At four stages the CV rank and K − 1 differ differently; the design left them to a paper grade.
- **The band's gains did not decide G.** They are recorded with the same test, at Bonferroni over 36 cells.
