# Within-stage permutation draws — pilot report (2026-11-16)

**Change:** `within-stage-permutation` (MOTCO 0.25.0). **Reads:** `results/within-stage-permutation-2026-11-16/`
(`within_stage_permutation_decision.json`; per profile the rate, paired and cell tables, `accounting.json`,
`reproduction.json` and `settings.json`; `plasmode/preflight.json`; `PROVENANCE.txt`, `sacct.txt`). **Parents:** the
[predictive-subspace pilot](predictive-subspace-2026-11-12.md) and the
[predictive size-excess diagnosis](predictive-size-excess-2026-11-13.md). **Concept page:**
[measurement visibility](../measurement-visibility.md#the-permutation-null-unrestricted-or-within-stage).

RRPP permutes the reduced model's residuals across all samples, which assumes they are exchangeable across stages.
The diagnosis (roadmap item 19, `both_restore`) found that restricting the draws within stage removes the
predictive coordinates' size excess on the SEA-AD plasmode anchor, and that production holds size and gains plant
power under the same draws. It opened this change: the within-stage draws as an opt-in option of `RRPP`, the harness
and `motco analyze`, read for size, power and specificity before `predictive_orthonormal` is re-read.

**In one paragraph.** **The verdict is `specificity_excess`.**

- **Size holds everywhere.** The within-stage null holds size in all four families: the simulated controls, the weak
  design, 500 SEA-AD label draws no run had read, and Gaussian data with stage-unequal dispersion, where the
  unrestricted null over-rejects when the largest stage is the narrowest (up to 0.102) and is conservative when it is
  the widest (down to 0.024).
- **Protected power holds at its limit.** It costs `angle` power on three orientation cells, exactly the three the rule
  allows.
- **What fails is specificity in the stage-effect band.** On in-span rotations the within-stage draws add off-target
  `delta` and `shape` rejections above the bound, 7 readings against the 4 allowed, one of them below the family's
  ceiling (`delta` 0.08 → 0.21 at the SEA-AD table, δ = 0.5, rotation e = 0.5).
- **The consequence fixed with the rule.** The null stays opt-in, documented with those readings, and
  `predictive_orthonormal` is not reopened.
- **The candidate, read anyway.** Under the same draws `predictive_orthonormal` qualifies against production by the
  pilot's four criteria. That reading is recorded and does not bear on the verdict.

## 1. Configuration and provenance

- **Arms.** One production fit per fitted unit. The scores (`pls_cv`) and the three predictive coordinate systems of
  that fit (`predictive`, `predictive_orthonormal`, `predictive_fitted`), each read under the unrestricted draws and,
  with the suffix `__within_stage`, under the package's `permutation_strata="stage"`, at the unit's RRPP seed.
  Every arm went through the package option, not a driver loop.
- **Profiles.**

  | profile | units | design | host |
  |---|---|---|---|
  | `simulation` | 7,800 | the pilot's seven stage-span designs, eleven cells and surgery cell, fresh base seed 2038, 100 × 199, p < 0.05 | `ing` |
  | `weak_size` | 1,000 | the pilot's weak-design anchor and translation, 500 × 999, p < 0.05 | `ing` |
  | `plasmode` | 1,000 | the SEA-AD astrocyte anchor at the `paper` stream's replicates 0–499 (`anchor`) and 500–999 (`anchor_fresh`), 999 permutations, every arm at rank 5, six covariates and `df`, p ≤ 0.05 | workstation |
  | `dispersion` | 10,000 | Gaussian matrices on the SEA-AD cell table with no group difference: five variance patterns × d ∈ {2, 5} × 1,000, 999 permutations, p < 0.05, read by the RRPP primitive under both draws | workstation |

- **The dispersion patterns** are variance multipliers by stage (stages of 20 / 20 / 40): `equal`; `predictive`
  1.120 / 1.278 / 0.802 and `scores` 0.805 / 0.921 / 1.137, the diagnosis's median relative dispersion by stage in
  `predictive_orthonormal` and in the scores; and each with its multipliers squared (`_doubled`). The stage centroids
  give a population stage share of 0.84 at d = 2 under equal dispersion.
- **Checks on every fitted unit.** The pilot's five checks; `same_statistics` (each within-stage arm's observed
  statistics equal its counterpart's exactly); and on every plasmode unit `construction_identity` (each within-stage
  arm's p-values equal the diagnosis's within-stage loop on the same latent matrix, models, correction and seed,
  exactly).
- **Code and hosts.** Revision `b72c6e8` on `proposal/within-stage-permutation`, pushed before any unit ran. The
  workstation ran Python 3.11.15 and `ing` Python 3.11.16, both with numpy 2.3.5, scipy 1.16.3, scikit-learn 1.8.0
  and pandas 2.3.3, BLAS pinned to one thread. `ing` arrays 939334 and 939335; 120 of 120 tasks exited 0:0.
  57.12 recorded core-hours on `ing` against a cap of 85, and 2.59 on the workstation (`PROVENANCE.txt`). The
  workstation rebooted during the dispersion run; the 6,058 complete records were intact and the run resumed the
  rest.

## 2. Accounting and reproduction

- **Units.** 19,800 of 19,800 complete: no failed unit, no failed arm, no inconsistent unit and no duplicate.
- **Checks.** Every applicable check held on every unit, `same_statistics` and `construction_identity` exactly.
- **Pre-flight: passed.** On replicates 0 and 500 the cross-validated fit selected rank 5 at held-out AUROC 0.635556
  and gave the forced call's statistics and p-values exactly.
- **Reproduction: passed.**
  - **Weak size:** the four unrestricted arms' counts equal the pilot's committed table on 24 of 24 rows.
  - **Plasmode `anchor`:** the unrestricted arms equal the pilot's counts (12 of 12 rows, for example
    `pls_cv` 21 / 24 / 13 and `predictive_orthonormal` 67 / 9 / 38). The within-stage arms equal the diagnosis's
    (12 of 12, for example `pls_cv__within_stage` 19 / 27 / 19 and `predictive_orthonormal__within_stage`
    29 / 17 / 32). This is the package option reproducing the diagnosis's independent loop on real data.

## 3. The decision

Rule version 1 has two layers (`design.md` D12).

- **Layer 1, the null** (`pls_cv__within_stage` against `pls_cv`):
  - **Size:** the stage-span paper grade's count-plus-ceiling rule in four families.
  - **Protected power:** the pilot's paired loss rule on the 54 target cells at δ ≥ 0.25.
  - **Specificity:** the bound rule on the δ = 2.0 off-target readings. In the band, a reading counts only if it is
    above the size bound *and* in significant paired excess over production (p < 0.025), at most 4 of 72 and none
    below 0.05 / 72.
- **Layer 2, the candidate** (`predictive_orthonormal__within_stage` against `pls_cv__within_stage`): the pilot's S1–S3
  and G with the reference changed to production under the same draws. It bears on the verdict only when layer 1
  holds.

| criterion | family | m | reading | holds |
|---|---|---|---|---|
| size | simulation controls | 42 | 0 above the bound; rates 0.01–0.09 | yes |
| size | weak size | 6 | 0 above; 0.040–0.048 | yes |
| size | plasmode `anchor_fresh` | 3 | 0 above; 0.046–0.064 | yes |
| size | dispersion | 30 | 1 above (allowed 3), none above the ceiling; 0.034–0.066 | yes |
| protected power | target cells, δ ≥ 0.25 | 54 | 3 one-sided loss p-values below 0.025 (allowed 3), none below 0.05 / 54 | yes |
| specificity | δ = 2.0 off-target | 36 | 0 above the bound | yes |
| specificity | band off-target | 72 | **7 counted (allowed 4), 1 below the ceiling** | **no** |

**The verdict is `specificity_excess`.** The gate passed. The consequence fixed with the rule: "The null holds size and
reads more off-target rejections than the bound allows. It stays opt-in, documented with those readings;
`predictive_orthonormal` is not reopened."

## 4. Size

- **Simulated controls.** Under the within-stage draws production's anchor and translation read 0.01–0.09 over 42
  readings at 100 units each, none above the bound of 0.0936. Where the unrestricted null is conservative the
  within-stage one reads closer to α. For example, the SEA-AD table's anchor `delta` is 0.01 under unrestricted draws
  and 0.06 under within-stage ones, and the weak design's 0.03 and 0.08.
- **Weak design (500 × 999).** `pls_cv__within_stage` reads 23 / 21 / 20 on the anchor and 24 / 20 / 22 on the
  translation, against 17 / 19 / 25 and 18 / 13 / 26 under unrestricted draws.
- **Fresh SEA-AD draws (replicates 500–999).** `pls_cv__within_stage` reads 32 / 26 / 23, inside the bound of 34.
  Production under unrestricted draws reads 32 / 24 / 12 on the same draws. On the reproduced draws (0–499) the
  diagnosis read 19 / 27 / 19, so the label draws alone move production's `delta` from 21 to 32 of 500.
- **Dispersion.** Rejections of 1,000 (`delta` / `angle` / `shape`), unrestricted | within stage:

  | pattern | d = 2 | d = 5 |
  |---|---|---|
  | `equal` | 66 / 59 / 55 \| 66 / 57 / 53 | 62 / 61 / 59 \| 51 / 58 / 58 |
  | `predictive` | 79 / 50 / 73 \| 54 / 53 / 55 | 77 / 48 / 59 \| 60 / 40 / 46 |
  | `predictive_doubled` | 91 / 41 / 76 \| 49 / 39 / 50 | 102 / 56 / 93 \| 52 / 58 / 48 |
  | `scores` | 46 / 39 / 33 \| 50 / 43 / 34 | 38 / 38 / 42 \| 41 / 44 / 46 |
  | `scores_doubled` | 34 / 28 / 24 \| 57 / 38 / 38 | 38 / 31 / 33 \| 52 / 50 / 52 |

  - **The mechanism, confirmed at the planned size.** When the largest stage is the narrowest (the predictive
    coordinates' pattern) the unrestricted draws over-reject `delta` and `shape`, up to 0.102 and 0.093 with the
    inequality doubled. When it is the widest (the scores' pattern) they are conservative, down to 0.024. The
    within-stage draws read 0.034–0.066 in every cell.
  - **The one exceedance** is `equal`, d = 2, `delta`: 66 of 1,000 under both draws. The two arms read the same
    matrices, so it is a property of those matrices, not of the restriction; one exceedance is allowed at m = 30.

## 5. Protected power

The three readings below 0.025 are all `angle` on orientation cells (rates `pls_cv` → `pls_cv__within_stage`; b
counts units only production rejects, c units only the within-stage arm rejects):

| design | cell | `angle` | b / c | one-sided p (loss) |
|---|---|---|---|---|
| `seaad` (δ = 2.0) | orientation e = 0.5 | 0.96 → 0.86 | 12 / 2 | 0.0065 |
| `n200_delta1` | orientation e = 0.25 | 0.62 → 0.53 | 11 / 2 | 0.011 |
| `n200_delta0.5` | orientation e = 0.25 | 0.27 → 0.19 | 8 / 0 | 0.0039 |

The rule allows 3 at m = 54, so the criterion holds at its limit. Beside the losses the within-stage draws gain:
`delta` on magnitude and `shape` on shape plants, in 4 of the 18 δ = 2.0 target cells and 5 of the 36 band cells at
one-sided p < 0.025 (for example the SEA-AD table's `shape` at e = 0.1, b / c = 0 / 14). The weak design gains
nowhere.

## 6. Specificity: the band failure

The 7 counted readings are all on rotation plants:

| design | cell | off-target | `pls_cv` → `pls_cv__within_stage` | b / c | status |
|---|---|---|---|---|---|
| `n200_delta1` | orientation e = 0.25 | `delta` | 0.05 → 0.12 | 0 / 7 | excess |
| `n200_delta1` | orientation e = 0.5 | `delta` | 0.07 → 0.16 | 0 / 9 | excess |
| `n200_delta0.5` | orientation e = 0.1 | `shape` | 0.08 → 0.14 | 0 / 6 | excess |
| `n200_delta0.5` | orientation e = 0.5 | `delta` | 0.28 → 0.34 | 0 / 6 | excess |
| `n200_delta0.25` | orientation e = 0.25 | `delta` | 0.08 → 0.14 | 0 / 6 | excess |
| `seaad_delta0.5` | orientation e = 0.25 | `shape` | 0.15 → 0.23 | 0 / 8 | excess |
| `seaad_delta0.5` | orientation e = 0.5 | `delta` | 0.08 → 0.21 | 0 / 13 | below the ceiling (p = 0.00012 < 0.05 / 72) |

- **It amplifies a leak production already has.** Under unrestricted draws, 12 of production's 72 band off-target
  readings are above the bound; under within-stage draws, 23 are. The predictive-subspace pilot had recorded
  production's band leak under the planted scale. The rule charges the null only for readings it pushes both above
  the bound and significantly above production, and 7 qualify.
- **The pattern** (descriptive, not tested). On an in-span rotation the within-stage draws trade `angle` for `delta`:
  they lose `angle` power on the same orientation cells where they add off-target `delta`. The planted interaction is
  in the reduced model's residuals, because the reduced model omits it. Unrestricted draws carry each stage's share of
  it to other stages; within-stage draws keep it in its own stage and exchange it only between groups. The `angle`
  null's median 95th percentile is slightly wider under within-stage draws on those cells (for example 33.9° → 34.5°
  at the SEA-AD table's rotation e = 0.5, and 23.4° → 23.8° at `n200_delta0.5`, e = 0.25). The `delta` rejections
  rise on the same cells. A reading of why is beyond this run.
- **At δ = 2.0 nothing leaks.** No off-target reading of the 36 is above the bound under within-stage draws.

## 7. Layer 2: `predictive_orthonormal` under the same draws (recorded; does not bear on the verdict)

`predictive_orthonormal__within_stage` against `pls_cv__within_stage` **qualifies**:

- **S1 size** holds in all three families. Its one simulation exceedance is the SEA-AD table's translation `angle`
  (11 of 100). On the fresh SEA-AD draws it reads 36 / 15 / 26: `delta` is above the bound (34) and at the ceiling
  (36), the one exceedance m = 3 allows.
- **S2 protected power** holds with no loss below 0.025.
- **S3 specificity** holds in both families: no δ = 2.0 reading above the bound, and no band reading in significant
  excess over production under the same draws. Six of its off-target readings over both families are above the bound.
- **G gain** holds in 10 of the 18 δ = 2.0 target cells at 0.05 / 18, mostly orientation `angle`. For example the
  SEA-AD table's rotation e = 0.25, b / c = 0 / 51, and n = 200's rotation e = 0.1, 0 / 24. Recorded beside it: 30 of
  36 band cells gain at p < 0.025, and none of the weak design's 9.

The other two coordinate systems under the same draws: `predictive__within_stage` reads `specificity_excess` (S3), and
`predictive_fitted__within_stage` `size_excess` (its fresh-anchor `delta` is 37 of 500, above the ceiling). Under the
unrestricted draws `predictive_orthonormal` against `pls_cv` again reads `size_excess`: its fresh-anchor `delta` is 69
of 500, as on the pilot's draws (67).

## 8. Descriptive readings

- **The anchor's 1,000 draws** (rejections of 500 at replicates 0–499 | 500–999):

  | arm | `delta` | `angle` | `shape` |
  |---|---|---|---|
  | `pls_cv` | 21 \| 32 | 24 \| 24 | 13 \| 12 |
  | `pls_cv__within_stage` | 19 \| 32 | 27 \| 26 | 19 \| 23 |
  | `predictive_orthonormal` | 67 \| 69 | 9 \| 4 | 38 \| 37 |
  | `predictive_orthonormal__within_stage` | 29 \| 36 | 17 \| 15 | 32 \| 26 |
  | `predictive__within_stage` | 31 \| 35 | 17 \| 15 | 33 \| 24 |
  | `predictive_fitted__within_stage` | 28 \| 37 | 17 \| 15 | 33 \| 24 |

  The diagnosis's marginal `shape` for the predictive coordinates under within-stage draws (32–33 against the
  threshold of 34) falls to 24–26 on the fresh draws; their `delta` rises from 28–31 to 35–37.
- **The surgery cell** (the default orientation surgery at e = 1.0 on the SEA-AD table): production reads `delta` /
  `angle` / `shape` 0.29 / 0.64 / 0.47 under unrestricted draws and 0.44 / 0.64 / 0.63 under within-stage ones.
  `predictive_orthonormal` reads 0.49 / 0.41 / 0.62 and 0.49 / 0.43 / 0.61.
- **The diagnosis's real-noise plant** is quoted from its committed table in the decision's descriptive block.
  Production's `delta` rises from 0.51 to 0.63 at e = 0.5 under within-stage draws.

## 9. What follows

- **The verdict's consequence.** The within-stage null is an opt-in option, documented as holding size, including where
  stage dispersion is unequal and where the unrestricted null fails, and as reading more off-target `delta` and
  `shape` on in-span rotations in the stage-effect band than production does. The unrestricted draws stay the default
  of `RRPP`, the harness and `motco analyze`. No case study reports within-stage p-values as its test, and the
  predictive coordinates stay opt-in diagnostics.
- **The verdict does not say the unrestricted null is right.** On Gaussian data with the scores' dispersion pattern
  the unrestricted null is conservative (0.024–0.046 at 1,000 units), and with the predictive coordinates' pattern it
  is liberal. On the fresh SEA-AD draws the predictive coordinates still exceed under it (`delta` 69 of 500). The
  size question the diagnosis opened is answered for the null; specificity on rotations is what keeps it opt-in.
- **Questions this run leaves.**
  - Whether the band leak belongs to the null or to the planted scale's rotations in a small band of designs. Every
    counted reading is at δ ≤ 1.0 under the planted scale, where production itself leaks.
  - Whether a null that is exact within stage can avoid trading `angle` for `delta` on rotations.
  - Whether the candidate's gain under the same draws survives a null that passes specificity. Its layer-2 reading is
    the strongest the predictive coordinates have had, so it is the one to re-read if a null qualifies.
- **The rule is not amended after the fact.** The readings that failed were declared before any unit ran.

## 10. Limitations

- **The anchor rates are conditional on one fit.** All 1,000 anchor draws share one fit, as in both parents.
- **The simulation profile is pilot grade** (100 × 199). The band criterion counts readings at 100 units, where one
  unit is 0.01.
- **The dispersion profile is a toy.** It is isotropic Gaussian noise with no fit; it validates the null under one
  mechanism and does not account for the anchor.
- **The fresh draws are the first 500 unread ones.** Production's own `delta` moved from 21 to 32 of 500 between the
  two halves, about 1.6 standard errors of the difference: a single block of 500 draws carries a binomial SE of about
  5 counts, which is the margin the marginal readings in this report sit within.
- **`df` keeps its definition under strata**, which over-corrects by √(n / (n − K)): 1.9 % at n = 80 and three
  stages, in the conservative direction (`design.md` D4).
- **The mechanism of the band leak is described, not tested.**
