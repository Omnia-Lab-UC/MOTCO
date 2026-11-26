# Simulation-only candidate comparison — report (2026-11-25)

**Change:** `simulation-only-candidate-comparison` (MOTCO 0.27.0: the opt-in generator setting
`stage_dispersion_scale`; roadmap open item 1 and Done items 17, 18 and 22). **Reads:**
`results/candidate-simulation-2026-11-25/` (`candidate_simulation_decision.json`; per profile `settings.json`,
`rate_table.csv`, `paired_table.csv`, `cell_table.csv`, `rank_table.csv`, `accounting.json`;
`{committed,surgery,adni_design}/reproduction.json`; `PROVENANCE.txt`, `sacct.txt`). **Parents:** the
[univariate baseline](univariate-baseline-2026-11-23.md)'s `stage_span` and `surgery` profiles, the
[stage-effect scaling](stage-effect-scaling-2026-11-24.md) run's `saturated`, `band` and `calibration` profiles, and the
[ADNI gain-and-specificity](adni-gain-specificity-2026-11-17.md) run's `simulation` profile (their settings and rate
tables, bound by sha256). **Concept page:**
[measurement visibility](../measurement-visibility.md#the-candidates-in-simulation-alone).

Three measurement candidates were kept out of production by readings on real noise: `predictive_orthonormal` by the
SEA-AD plasmode's size and ADNI's protected power, and block weighting and sparse stage-PLS by a pilot whose gain rule
counted only the weak design and the SEA-AD plant. Neither export is the data MOTCO was designed for. On 2026-11-25 the
user decided to re-read the candidates in simulation alone, across more conditions, with the real-noise readings kept
as recorded limits and not as vetoes. This run does that on 58,400 simulated datasets, each read by every arm on the
same draws, under a decision fixed before any unit ran.

- **The arms.** `pls_cv` is production: stage-PLS scores at the CV rank. `predictive_orthonormal` measures the K − 1
  stage-predictive directions of the same fit. `weighted` is its own CV fit on blocks scaled to unit total variance.
  `sparse_block` is its own CV fit with a feature budget per block, on three designs only.
- **The families.** The decision reads each arm against production in five design families (`core`, `weak_stage`,
  `four_stages`, `unbalanced`, `adni_design`), so a candidate that helps at some designs and harms at others is
  described as such.
- **Two riders.** A stage-dispersion criterion on null designs, and the incidental routes (feature-set surgeries and
  out-of-span rotations), which are labelled and never decide.

**In one paragraph.**

- **The verdict is `scoped_candidates`.** `predictive_orthonormal` and `weighted` each qualify in three of the five
  families. No arm qualifies everywhere, so no default change is proposed, and every arm stays opt-in.
- **The gate passed.** Production reproduced every committed count it was compared with: 228 rows on `committed`, 66
  on `surgery` and 60 on `adni_design` (both arms there), none differing.
- **Production holds size in every family and under unequal noise by stage.** Across the 10 stage-dispersion designs
  it reads 2 of 60 readings above the bound and none above the ceiling. It loses nothing to a candidate anywhere: no
  arm loses protected power in any family.
- **`predictive_orthonormal` qualifies in `core`, `weak_stage` and `four_stages`.**
    - **The gains follow the rank.** It gains where the CV rank exceeds K − 1. Where it does not, it measures
      production's subspace in another metric and reads like production. For example, `angle` at n = 200, δ = 0.5
      and a 9° rotation goes from 0.080 to 1.000, and `delta` on the band's stretch at e = 1 from 0.460 to 0.995.
    - **It fails specificity in two families.** In `unbalanced`, one off-target `angle` reading is above the ceiling
      (0.110 against 0.100). At the ADNI design, stretches are read as `shape` (0.060 → 0.295), the failure known
      before the run.
- **`weighted` qualifies in `core`, `four_stages` and `unbalanced`.** It gains 4–7 points on average there. It
  exceeds size in `weak_stage` (`angle` 0.090 on `n320_delta1`'s translation, against a ceiling of 0.078) and gains
  nothing at the ADNI design.
- **`sparse_block` reads `keep_production`.** At `n200` it gains nothing, and on the band it adds off-target `shape`
  (0.050 → 0.135).
- **Both candidates fail size under unequal noise by stage, where production holds.** The failing readings sit where
  the largest stage has the narrowest noise on the SEA-AD cell table. Item 20 found that pattern over-rejects on
  Gaussian latent matrices: here `predictive_orthonormal` reads up to 0.090 and `weighted` up to 0.100, against a
  ceiling of 0.084.
- **What it holds for.** Simulated Gaussian noise on 658 features (at most about 8 per sample), three blocks, three or
  four stages. No real-noise reading entered the decision. The SEA-AD and ADNI readings stand as limits (§11), and size
  at many features per sample is unread (open item 15). No verdict changes what the harness or `motco analyze` runs.

## 1. Configuration and provenance

- **The unit.**
    - **One dataset, generated once,** read by every arm at the unit's RRPP seed (999 permutations).
    - **The production fit** (stage-PLS, rank floor 2, net displacement) gives `pls_cv`.
    - **`predictive_orthonormal`** is built from that fit with `pls_coordinates_from_pls`, so it makes no second CV
      fit and reports `pls_cv`'s rank. It measures the input coordinates on an orthonormal basis of the K − 1
      predictive directions, so where the rank is K − 1 it spans production's subspace in another metric.
    - **`weighted` and `sparse_block`** are their own CV fits, with only `integration_params` changed
      (`block_weighting = "equal_variance"`; `sparsity = {scope: block}` with the default ladder).
- **Five profiles, 58,400 units.** Controls are at 500 replicates and plants at 200 unless stated. All three blocks
  (methylation, expression, proteomics; 658 features), `joint` magnitude, `p_dmp` 0.1, ρ = 0.

  | profile | designs | cells | units |
  |---|---|---|---|
  | `committed` | `n80`, `n200`, `n500`, `n200_delta0.5` (the univariate baseline's), `n320_delta1`, `band_n320` (the stage-effect scaling run's), at their committed seeds | 13 per design (11 in the band) | 18,800 (+ 9,200 `sparse_block`) |
  | `conditions` | new, base seed 2043: `n200_stages4`, `n200_delta0.5_stages4`, `seaad_table`, `seaad_table_delta0.5` (F 11 / 10 / 28, M 9 / 10 / 12), `n200_share0.8` | 13 | 16,000 |
  | `dispersion` | new, null only: four stage-dispersion patterns on the two SEA-AD tables, the two strong ones at balanced n = 200 | anchor, translation | 10,000 |
  | `surgery` | the design-guidance slice's three-block columns at n = 80 and 200, at their committed seeds | 11, at 500 replicates | 11,000 |
  | `adni_design` | the ADNI gain run's simulated design, at its committed seeds | 10 (500 / 200) | 2,600 |

  The 13 cells are the anchor, the translation, three stage-span stretches, three in-span rotations (φ = 0), three
  in-span bends, and two out-of-span rotations (φ = 90°). The band's 11 have no out-of-span cell. `sparse_block` ran on
  `n200`, `n200_delta0.5` and `band_n320` only, and is recorded as not run elsewhere.
- **The stage-dispersion patterns.** SD multipliers by stage, applied to both groups about unchanged means. The
  sources are the square roots of item 20's variance multipliers, and the strong patterns are those multipliers taken
  as SDs.

  | pattern | multipliers | where |
  |---|---|---|
  | `last_narrow` | (1.058, 1.130, 0.896) | both SEA-AD tables |
  | `last_wide` | (0.897, 0.960, 1.066) | both SEA-AD tables |
  | `last_narrow_strong` | (1.120, 1.278, 0.802) | both SEA-AD tables and n = 200 |
  | `last_wide_strong` | (0.805, 0.921, 1.137) | both SEA-AD tables and n = 200 |

  On the SEA-AD table the last stage is the largest (40 of 80), so `last_narrow` puts the narrowest noise in the
  largest stage.
- **The families.**

  | family | designs |
  |---|---|
  | `core` | `n80`, `n200`, `n500` |
  | `weak_stage` | `n200_delta0.5`, `n320_delta1`, `band_n320` |
  | `four_stages` | `n200_stages4`, `n200_delta0.5_stages4` |
  | `unbalanced` | `seaad_table`, `seaad_table_delta0.5`, `n200_share0.8` |
  | `adni_design` | `table_adni` |

- **Code and host.** Revision `4982bc4`, pushed before any array unit ran, from a clean worktree on `ing`. Python
  3.11.16, numpy 2.3.5, scipy 1.16.3, scikit-learn 1.8.0, pandas 2.3.3, BLAS pinned to one thread.
- **Where it ran.** Six arrays on partition `512x1024` with n11 excluded, 1,570 of 1,570 tasks `COMPLETED 0:0`, from
  14:38 to 21:58 on 2026-11-25. Merge and `decide` ran on a compute node the next morning.
    - **Tier 1** (`committed`, `conditions`, `dispersion` with `pls_cv`, `predictive_orthonormal`, `weighted`): 647.0
      allocated CPU-hours against a cap of 700.
    - **Tier 2** (`surgery`, `adni_design`): 262.7 against 300.
    - **Tier 3** (`sparse_block` on its three designs): 317.7 against 550.
    - **1,120.4 recorded core-hours** in all (`accounting.json`).
    - **Per arm** (median seconds): `pls_cv` and `predictive_orthonormal` 0.8–1.1 each after the shared fit;
      `weighted` 15.5–33.2; `sparse_block` 88.7–123.3. At the ADNI design each CV arm takes about 49 s and the
      predictive arm 0.7 s.
- **A first submission was cancelled.** Two arrays received the arm list split by `sbatch --export` and ran `pls_cv`
  alone. They were cancelled within minutes, and their 49 valid records were kept and skipped on resume
  (`PROVENANCE.txt`).

## 2. Accounting and the reproductions

- **Records.** Every profile has one record per expected (unit, arm): 0 missing, 0 duplicates, 0 inconsistent.
  Tier 3's append wrote only `sparse_block` records, and `committed` holds 18,800 records for each tier-1 arm after it.
- **Failed units.** 17 of 58,400 units, every one the planted pooled SD's fixed point not converging. Each failed
  with all of its arms, so the pairing is intact, and failed units are excluded from the rates.

  | profile | design | cell | failed of 200 (500 for the translation) |
  |---|---|---|---|
  | `conditions` | `n200_stages4` | rotation e = 0.1 / 0.25 / 0.5 | 1 / 4 / 10 |
  | `conditions` | `n200_delta0.5_stages4` | rotation e = 0.5 | 1 |
  | `dispersion` | `seaad_table_last_narrow` | translation | 1 |

  The design named this risk at four stages, where the fixed point had cycled in 12 of 1,000 rotation units before.
  The smoke ran 2 replicates per cell and found none. A rate of 5% is not visible at that size, and the design's
  amendment threshold applied before any array unit, so no design was amended. The worst cell keeps 190 units.
- **Reproductions: passed.**
    - `committed`: 228 rows, `pls_cv`'s available and rejected counts per design, cell and statistic, equal to the
      univariate baseline's and the stage-effect scaling run's committed rate tables.
    - `surgery`: 66 rows, equal to the univariate baseline's `surgery` table.
    - `adni_design`: 60 rows, `pls_cv` and `predictive_orthonormal`, equal to the gain run's simulation table.

  No row differs.
- **The decision regenerates.** `decide` on the committed aggregates alone, on the workstation, writes a decision
  byte-identical to the cluster's (sha256 `a8c37493…`).

## 3. The decision

- **The readings.** `delta`, `angle` and `shape`, each at p < 0.05. On `adni_design` the rule is p ≤ 0.05, the
  parent's.
- **The criteria, within a family.** These are the predictive pilot's, with `judge_rates` and the exact paired tests
  imported unchanged:
    - **S1, size.** The anchors and translations under the count-plus-ceiling rule (at most the allowed count above
      the per-cell bound, 35 or more of 500, and none above the family's Bonferroni ceiling).
    - **S2, protected power.** The in-span target cells (stretch → `delta`, rotation at φ = 0 → `angle`, bend →
      `shape`) by a one-sided exact paired loss against `pls_cv`, under the count-plus-ceiling rule.
    - **S3, specificity.** The target cells' two off-target statistics. At δ = 2.0 designs the arm's own rates
      against the bound (17 or more of 200). Below δ = 2.0, where production's own off-target rates exceed it, the arm's
      paired excess over `pls_cv`.
    - **G, gain.** At least one target cell with an exact paired gain p-value at or below 0.05 divided by the
      family's target cells.
- **The label** is the first that applies, in the order `size_excess`, `power_loss`, `specificity_excess`, `no_gain`,
  `qualifies`. A family is `not_run` for an arm that has no design in it, and `partial` when the arm has some of its
  designs. `partial` never counts toward qualifying everywhere.
- **D, stage dispersion.** Each arm's 60 readings on the ten dispersion designs, under the count-plus-ceiling rule.
  If production failed, every candidate would instead be judged by its paired excess over `pls_cv`. Production
  holds (§8), so the rule applies as written.
- **I, incidental routes.** Never decisive. Each candidate is labelled against `pls_cv` per cell on the 18 surgery
  target cells and the 20 out-of-span rotation cells, by one-sided exact paired tests at 0.025 / 38 with a margin
  of 0.10.
- **Per arm.** `qualifies_everywhere` means every family is run and qualifies, and D holds. `scoped_candidate` means
  at least one family qualifies. Otherwise the arm is `keep_production`.
- **Overall.** `propose_default_change` (one arm qualifying everywhere) or `scoped_candidates` or `keep_production`.

| arm | `core` | `weak_stage` | `four_stages` | `unbalanced` | `adni_design` | D | verdict |
|---|---|---|---|---|---|---|---|
| `predictive_orthonormal` | qualifies | qualifies | qualifies | specificity_excess | specificity_excess | fails | scoped_candidate |
| `weighted` | qualifies | size_excess | qualifies | qualifies | no_gain | fails | scoped_candidate |
| `sparse_block` | no_gain (partial) | specificity_excess (partial) | not_run | not_run | not_run | not_run | keep_production |

**Verdict: `scoped_candidates`.** The gate passed, and two arms qualify in at least one family. Neither qualifies
everywhere, so no default change is proposed.

## 4. Production

- **Size holds in every family.**

  | family | m | allowed above the bound | above the bound | above the ceiling |
  |---|---|---|---|---|
  | `core` | 18 | 2 | 0 | 0 |
  | `weak_stage` | 18 | 2 | 2 | 0 |
  | `four_stages` | 12 | 1 | 0 | 0 |
  | `unbalanced` | 18 | 2 | 0 | 0 |
  | `adni_design` | 6 | 1 | 0 | 0 |

    - **The two `weak_stage` readings** are `angle` on `n320_delta1`'s anchor (35 of 500) and translation (38). The
      ceiling there is 39.
    - **Four stages at n = 200 and the SEA-AD table on three blocks** are read for production's size here for the
      first time at paper grade. The share of 0.8 at n = 200 repeats item 26's null reading, now beside planted cells.
- **Off-target rates below δ = 2.0** (descriptive): 29 of the 106 readings are above the per-cell bound. The largest
  are in-span rotations read as `shape`:
    - 0.665 at `n200_delta0.5`, e = 0.5, and 0.663 at four stages;
    - 0.315 on the SEA-AD table at δ = 0.5, e = 0.5;
    - 0.215 for the band's stretch at e = 2 read as `angle`.

  This is item 24's `detects_only` grade, now also at four stages and on the unbalanced table.

## 5. `predictive_orthonormal`

| family | S1 above the bound (allowed) | S2 losses | S3 | G gains / target cells | label |
|---|---|---|---|---|---|
| `core` | 0 (2) | 0 | holds: 2 above the bound (3 allowed), 0 above the ceiling | 2 / 27 | qualifies |
| `weak_stage` | 1 (2) | 0 | holds: 1 paired excess (3 allowed) | 18 / 27 | qualifies |
| `four_stages` | 0 (1) | 0 | holds: 0 | 8 / 18 | qualifies |
| `unbalanced` | 0 (2) | 0 | **fails**: 3 above the bound (3 allowed), 1 above the ceiling | 8 / 27 | specificity_excess |
| `adni_design` | 0 (1) | 0 | **fails**: 3 paired excesses (2 allowed), 1 below the ceiling | 5 / 8 | specificity_excess |

- **Size and protected power hold everywhere.** The one reading above the bound is `angle` on `n320_delta1`'s
  translation (35 of 500).
- **The gains are large where production's rank exceeds K − 1.**
    - **Weak stage effect:** at `n200_delta0.5` (rank 3, two measured dimensions) the 9° rotation goes from 0.080 to
      1.000 and the 22.5° rotation from 0.405 to 1.000. In the band (rank 10) the stretch at e = 1 goes from 0.460 to
      0.995.
    - **Four stages:** at `n200_delta0.5_stages4` (rank 4, three measured dimensions) the 9° rotation goes from 0.080
      to 0.970 and the bend at e = 0.25 from 0.410 to 1.000.
    - **Unbalanced:** on the SEA-AD table at δ = 0.5 (rank 3) the 22.5° rotation goes from 0.145 to 1.000.
    - **ADNI design** (rank 19): the stretch at e = 1 goes from 0.535 to 0.990.
    - **Mean paired difference over a family's target cells:** +0.233 in `weak_stage`, +0.206 in `four_stages`,
      +0.132 in `unbalanced`, +0.201 at the ADNI design, and +0.011 in `core`.
- **Where the rank is K − 1, the arm measures production's subspace.** At `n80`, `n200`, `n500`, `n320_delta1`,
  `n200_stages4`, `n200_share0.8` and the SEA-AD table at δ = 2.0 the median rank equals K − 1. There the two arms
  differ only in metric (orthonormal against scores) and read alike in power. `core`'s two gains come from the
  units where the CV picked one more component (`rank_table.csv`). At `n80`'s 9° rotation, production rejects 178 of
  187 units at rank 2 and 0 of 13 at rank 3; the candidate rejects 178 of 187 and 13 of 13. At the 22.5° rotation,
  production reads 4 of 30 at rank 3 and the candidate 30 of 30.
- **What it gives up.**
    - **`unbalanced`.** At `n200_share0.8`, bends read as `angle` at 19, 19 and 22 of 200 (e = 0.1, 0.25, 0.5). The
      last is above the family's ceiling (0.100). One paired excess below δ = 2.0 is also recorded: the SEA-AD table
      at δ = 0.5, rotation e = 0.5, read as `shape` (0.315 → 0.455).
    - **`adni_design`.** The doubled stretch is read as `shape` (0.060 → 0.295, below the ceiling), the stretch at
      e = 1 too (0.045 → 0.110), and the rotation at e = 1 as `delta` (0.035 → 0.080). This is item 22's failure
      reproduced exactly; the proposal recorded it as known before any unit ran.
    - **Stage dispersion** (§8).

## 6. `weighted`

| family | S1 above the bound (allowed) | above the ceiling | S2 losses | S3 | G gains / target cells | label |
|---|---|---|---|---|---|---|
| `core` | 2 (2) | 0 | 0 | holds: 0 | 8 / 27 | qualifies |
| `weak_stage` | 2 (2) | **2** | 0 | holds: 3 paired excesses (3 allowed) | 10 / 27 | size_excess |
| `four_stages` | 0 (1) | 0 | 0 | holds: 0 | 7 / 18 | qualifies |
| `unbalanced` | 2 (2) | 0 | 0 | holds: 2 above the bound (3 allowed) | 10 / 27 | qualifies |
| `adni_design` | 0 (1) | 0 | 0 | holds: 0 | 0 / 8 | no_gain |

- **The size failure is on `n320_delta1`.** `angle` reads 41 of 500 at the anchor and 45 at the translation, against
  a ceiling of 39. Production reads 35 and 38 on the same units (§4), so the weighted fit adds a few rejections where
  production already sits at the bound. The other readings above the bound are `delta` at `n200` (35 and 37) and
  `shape` at `n200_share0.8` (35 and 35), within the allowed count.
- **The gains are modest and broad.** The mean paired difference over a family's target cells is +0.043 in `core`,
  +0.066 in `four_stages`, +0.067 in `unbalanced` and +0.079 in `weak_stage`.
    - At `n80`: a bend at e = 0.1 goes from 0.485 to 0.745 and the stretch at e = 0.05 from 0.485 to 0.710.
    - On the SEA-AD table: the bend at e = 0.1 goes from 0.455 to 0.675.
    - At `n320_delta1`: the 9° rotation goes from 0.725 to 0.995.

  The weighting changes which directions the fit finds, not the rank: the median rank equals production's at every
  design, and the pooled-stage capture rises by at most 0.014 (§10). With 367 / 131 / 160 features, unweighted
  methylation carries more than half of the joint variance.
- **No gain at the ADNI design** (mean paired difference −0.024).
- **`delta` is in the weighted metric.** The criteria compare rejections, not statistics, so they apply as written.

## 7. `sparse_block`

- **`core` (`n200` only): `no_gain`.** Size holds (one reading above the bound: `delta` 35 of 500 at the anchor).
  S2 and S3 hold, but no target cell gains. The median keep fraction is 1.0 in 12 of its 13 cells, so the sparse fit
  usually selects the dense one.
- **`weak_stage` (`n200_delta0.5` and `band_n320`): `specificity_excess`.**
    - **The excess.** The band's stretch at e = 2 is read as `shape` in 0.135 against 0.050, a paired excess below
      the ceiling, and as `angle` in 0.295 against 0.215.
    - **The gains (8 of 18 target cells),** smaller than the predictive arm's: the band's stretch at e = 1 goes from
      0.460 to 0.715 and the 9° rotation at δ = 0.5 from 0.080 to 0.265.
- **Not run** in `four_stages`, `unbalanced`, `adni_design` and on the dispersion designs. Both read families are
  `partial`, so `sparse_block` could not have qualified everywhere.
- **The cost.** About seven production fits per unit (median 88.7 s at `n200`, 123.3 s in the band).

## 8. Stage dispersion

Rejections of 500, `delta` / `angle` / `shape`. Bold marks readings above the ceiling (42 of 500, m = 60).

| design | cell | `pls_cv` | `predictive_orthonormal` | `weighted` | `pls_cv` within stage | `predictive_orthonormal` within stage |
|---|---|---|---|---|---|---|
| `n200_last_narrow_strong` | anchor | 31 / 16 / 23 | 31 / 16 / 24 | 40 / 12 / 25 | 25 / 25 / 22 | 24 / 25 / 22 |
| | translation | 22 / 21 / 16 | 34 / 11 / 23 | 34 / 18 / 27 | 27 / 32 / 27 | 23 / 24 / 31 |
| `n200_last_wide_strong` | anchor | 25 / 24 / 23 | 25 / 23 / 24 | 25 / 21 / 30 | 25 / 18 / 23 | 25 / 18 / 23 |
| | translation | 16 / 26 / 23 | 21 / 23 / 31 | 19 / 22 / 35 | 20 / 22 / 23 | 22 / 19 / 22 |
| `seaad_table_last_narrow` | anchor | 33 / 29 / 27 | 34 / 23 / 30 | 37 / 18 / 41 | 29 / 27 / 24 | 24 / 22 / 21 |
| | translation (499) | 22 / 34 / 28 | 32 / 27 / 33 | 39 / 24 / 40 | 29 / 28 / 25 | 25 / 26 / 23 |
| `seaad_table_last_narrow_strong` | anchor | 40 / 30 / 38 | 39 / 28 / **45** | **49** / 35 / **50** | 25 / 26 / 29 | 23 / 23 / 27 |
| | translation | 28 / 23 / 30 | 38 / 24 / **45** | **48** / 32 / **44** | 28 / 20 / 32 | 26 / 22 / 34 |
| `seaad_table_last_wide` | anchor | 23 / 23 / 24 | 24 / 21 / 26 | 19 / 24 / 22 | 29 / 26 / 33 | 25 / 21 / 32 |
| | translation | 19 / 30 / 19 | 24 / 19 / 28 | 22 / 29 / 24 | 29 / 35 / 24 | 25 / 24 / 32 |
| `seaad_table_last_wide_strong` | anchor | 18 / 24 / 14 | 20 / 23 / 16 | 17 / 24 / 16 | 26 / 27 / 22 | 24 / 26 / 19 |
| | translation | 12 / 25 / 18 | 15 / 21 / 32 | 18 / 28 / 18 | 29 / 32 / 29 | 28 / 24 / **43** |
| `seaad_table_delta0.5_last_narrow` | anchor | 30 / 33 / 19 | 38 / 24 / 32 | 20 / 32 / 8 | 35 / 29 / 20 | 30 / 22 / 24 |
| | translation | 28 / 34 / 13 | 35 / 26 / **43** | 27 / 34 / 16 | 29 / 29 / 25 | 29 / 26 / 33 |
| `seaad_table_delta0.5_last_narrow_strong` | anchor | 23 / 28 / 21 | 42 / 34 / 38 | 23 / 25 / 14 | 24 / 24 / 21 | 25 / 30 / 25 |
| | translation | 28 / 30 / 18 | **45** / 30 / 34 | 33 / 28 / 27 | 21 / 20 / 21 | 31 / 28 / 18 |
| `seaad_table_delta0.5_last_wide` | anchor | 7 / 29 / 8 | 26 / 20 / 22 | 7 / 30 / 8 | 17 / 27 / 24 | 28 / 17 / 20 |
| | translation | 11 / 28 / 8 | 30 / 19 / 16 | 11 / 29 / 10 | 21 / 26 / 16 | 31 / 23 / 19 |
| `seaad_table_delta0.5_last_wide_strong` | anchor | 5 / 17 / 8 | 14 / 20 / 15 | 4 / 18 / 5 | 16 / 22 / 18 | 16 / 21 / 24 |
| | translation | 4 / 22 / 8 | 17 / 19 / 19 | 12 / 36 / 4 | 19 / 29 / 17 | 20 / 24 / 26 |

- **Production holds:** 2 of 60 readings above the bound (`delta` 40 and `shape` 38 at
  `seaad_table_last_narrow_strong`'s anchor), none above the ceiling. The two strong patterns at balanced n = 200 move
  it little against the committed `n200` (anchor 34 / 22 / 23). Where the last and largest stage is the widest at
  δ = 0.5 it turns conservative (`delta` 4–11, `shape` 8). Against the same table at equal dispersion (anchor
  16 / 31 / 23, translation 14 / 32 / 17) this follows the pattern item 20 saw in the scores.
- **`predictive_orthonormal` fails:** 10 above the bound, 4 above the ceiling. Every reading above the bound is under a
  `last_narrow` pattern on a SEA-AD table, the pattern item 20 found to over-reject on Gaussian latent matrices.
- **`weighted` fails:** 12 above the bound, 4 above the ceiling. The four are `delta` and `shape` under
  `last_narrow_strong` on the SEA-AD table at δ = 2.0 (44–50 of 500).
- **The within-stage draws (recorded).**
    - **`pls_cv`** holds: 2 above the bound, none above the ceiling.
    - **`predictive_orthonormal`** brings every `last_narrow` reading under the bound, but one reading exceeds the
      ceiling on its own: `shape` 43 at `seaad_table_last_wide_strong`'s translation. Under the count-plus-ceiling rule
      that is an exceedance (`exceeds`).
    - **The recorded arms' observed statistics** equal their unrestricted arms' on every unit.
- **What the pattern became in the measured coordinates** (the reduced model's mean squared residual norm by stage
  over the pooled one, anchor medians):
    - **At δ = 2.0 each arm carries the feature-level pattern almost unchanged.** Under `last_narrow_strong` on the
      SEA-AD table, `pls_cv` reads (1.165, 1.510, 0.639), the predictive arm (1.155, 1.499, 0.649) and `weighted`
      (1.162, 1.525, 0.650).
    - **At δ = 0.5 the predictive arm damps it.** Its max / min ratio is 2.09 against 2.84 for `pls_cv`, yet it is the
      arm that exceeds there.

  The residual dispersion alone does not say which arm exceeds. At δ = 2.0 both arms measure the same two directions
  on most units (rank 2, q90 3) and differ in metric. At δ = 0.5 the candidate measures two directions to
  production's three. The candidate exceeds in both cases, always where the largest stage has the narrowest noise.

## 9. Incidental routes

| arm | cells | `candidate_loses` | `candidate_gains` | `unresolved` |
|---|---|---|---|---|
| `predictive_orthonormal` | 18 surgery + 20 out-of-span | 0 | 2 | 36 |
| `weighted` | 18 surgery + 20 out-of-span | 0 | 16 | 22 |
| `sparse_block` | 4 out-of-span (`n200`, `n200_delta0.5`) | 0 | 0 | 4 |

- **No candidate loses a cell by the margin.** The predictive pilot's surgery loss on two blocks (`angle` 0.66 → 0.48
  at the Phase 6 orientation, e = 1.0) does not recur on three blocks. There the candidate gains on the same plant:
  0.274 → 0.536 at n = 80 and 0.492 → 0.672 at n = 200.
- **`weighted` gains on 16 of the 18 surgery cells.** Every orientation and bend gains at both n, as do two stretches
  at each n. Examples: `angle` 0.274 → 0.578 at n = 80, e = 1, and `shape` 0.316 → 0.540 at n = 80, e = 0.5.
- **The out-of-span rotations are unresolved for every arm.** No paired difference reaches 0.10.

## 10. The fit by design

The anchor cell's medians (`cell_table.csv`). The predictive arm shares `pls_cv`'s fit, so its rank is `pls_cv`'s.

| design | `pls_cv` rank (q10–q90) | held-out AUROC | pooled-stage capture | predictive: measured dimension, capture | `weighted`: rank, capture |
|---|---|---|---|---|---|
| `n80` | 2 (2–2) | 1.000 | 0.934 | 2, 0.934 | 2, 0.940 |
| `n200` | 2 (2–2) | 1.000 | 0.973 | 2, 0.973 | 2, 0.976 |
| `n500` | 2 (2–2) | 1.000 | 0.989 | 2, 0.989 | 2, 0.990 |
| `n200_delta0.5` | 3 (3–3) | 1.000 | 0.905 | 2, 0.905 | 3, 0.910 |
| `n320_delta1` | 2 (2–3) | 1.000 | 0.971 | 2, 0.970 | 2, 0.971 |
| `band_n320` | 10 (5–19) | 0.664 | 0.288 | 2, 0.253 | 10, 0.302 |
| `n200_stages4` | 3 (3–3) | 1.000 | 0.967 | 3, 0.967 | 3, 0.971 |
| `n200_delta0.5_stages4` | 4 (4–4) | 1.000 | 0.881 | 3, 0.880 | 4, 0.887 |
| `seaad_table` | 2 (2–2) | 1.000 | 0.931 | 2, 0.931 | 2, 0.938 |
| `seaad_table_delta0.5` | 3 (3–3) | 1.000 | 0.773 | 2, 0.771 | 3, 0.785 |
| `n200_share0.8` | 2 (2–2) | 1.000 | 0.972 | 2, 0.972 | 2, 0.975 |
| `table_adni` | 19 (19–19) | 0.608 | – | 2, – | 19, – |

- **`sparse_block`:** rank 2 at `n200`, 3 at `n200_delta0.5` and 9 (4–19) in the band. The median keep fraction is
  1.0 at every anchor (q10 1.0, 1.0 and 0.2), and the captures are 0.973 / 0.902 / 0.272.
- **The predictive arm keeps nearly all of production's capture** at δ ≥ 0.5 (at most 0.002 lower). In the band it
  keeps 0.253 of the pooled stage program in two dimensions against 0.288 in ten.
- **The committed values** are the parents' cell tables; the anchors here are the same units.

## 11. The real-noise readings, as limits

None of these entered the decision. Each stands as recorded in its report.

- **`predictive_orthonormal`.**
    - **SEA-AD plasmode anchor (item 18): `size_excess`.** `delta` 67 of 500 and `shape` 38, against production's 21
      and 13 ([report](predictive-subspace-2026-11-12.md)).
    - **Its diagnosis (item 19): `both_restore`.** Within-stage draws or out-of-fold coordinates bring the anchor
      within the size rule ([report](predictive-size-excess-2026-11-13.md)).
    - **ADNI anchor (item 21): `excess_absent`.** The largest count is 29 of 500
      ([report](adni-real-noise-anchor-2026-11-17.md)).
    - **ADNI gain (item 22): `power_loss`.** `angle` on the in-sample 90° rotation goes from 0.285 to 0.05 at the
      cohort ([report](adni-gain-specificity-2026-11-17.md)).
- **`weighted`.** The production selection (item 17): `no_gain`. On the SEA-AD plasmode the blocks have equal feature
  counts, so the weighting is not a rebalancing, and the anchor reads 0.05
  ([report](production-latent-space-2026-11-10.md)).
- **`sparse_block`.** The production selection (item 17): `size_excess`. `delta` 0.11 on the plasmode anchor at 100
  replicates, against a ceiling of 0.10, where production reads 0.11 on the same replicates.
- **The simulated mechanism matches the first.** The stage-dispersion profile reproduces the direction of the SEA-AD
  excess in simulation (the largest stage with the narrowest noise), at a smaller size: 0.090 here against 0.134
  there. The SEA-AD design has 125 features per donor; no simulated design here has more than about 8.

## 12. What follows

- **The arms stay opt-in** (design, Migration Plan, `scoped_candidates`). The measurement-visibility page carries the
  map by family ([measurement visibility](../measurement-visibility.md#the-candidates-in-simulation-alone)).
  `motco analyze` and the harness run what they ran.
- **A rule that recommends an arm from a cohort's design is a separate change.** For `predictive_orthonormal` the map
  says what such a rule would need to read: whether production's CV rank exceeds K − 1 (where the gains are),
  whether the residual dispersion differs by stage with the narrowest noise in the largest stage (where size fails),
  and the off-target leak under unbalance and at ADNI's design.
- **The roadmap's reopening condition for `predictive_orthonormal`** is amended as the user's decision of 2026-11-25:
  the simulation-only reading replaces the requirement of a real-noise anchor, and the SEA-AD and ADNI readings are
  restated as limits.
- **Production's size at four stages, on the SEA-AD table with three blocks and under stage-unequal dispersion** is
  added to the roadmap's simulated conditions.
- **The roadmap records this run as Done item 28.** Unbalance at equal dispersion was diagnosed in parallel as Done
  item 27 (`small_cells`), and its 5-per-stage designs are outside this run's families.
- **Still open.** Size at many features per sample (open item 15), real multi-omic noise (open item 2), and the
  dispersion excess that within-group draws leave (open item 17).

## 13. Limitations

- **Simulated Gaussian noise.** Every reading holds at 658 features on three blocks, at most about 8 features per
  sample, three or four stages and no covariates outside the ADNI design. Neither the SEA-AD regime (125 features per
  donor) nor the target ROSMAP and ADNI multi-omic cohorts are read.
- **Two labels turn on one reading each.** `predictive_orthonormal`'s `unbalanced` label rests on one off-target
  reading above the ceiling (22 of 200 against 20). `weighted`'s `weak_stage` label rests on `n320_delta1`, where
  production is at the bound on the same units. The rule was fixed before any unit ran and is applied as written.
- **17 units failed the planted fixed point,** 16 of them at four stages. Their cells' rates rest on 190–199
  converged units, and a unit that fails to converge may differ from one that converges.
- **`sparse_block` ran on three designs.** Its verdict could not exceed `scoped_candidate`, and its dispersion
  reading is not run.
- **The dispersion profile reads size only.** Power under unequal noise by stage is not read.
- **The within-stage draws are recorded, not judged,** and they cost specificity in the band (item 20).
