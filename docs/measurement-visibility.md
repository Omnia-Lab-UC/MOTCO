# What the stage-supervised measurement can see

MOTCO's production measurement fits one PLS-DA to the **stage** label on both groups
pooled, chooses its rank by double cross-validation (floor 2), and measures `delta`,
`angle` and `shape` on the group × stage means of the latent scores. This page states
which group differences that measurement can see, why a rotation of the stage program
onto other features can vanish from it, and what a non-rejection therefore means. Each
claim links to the study that measured it.

**The latent space defines the question.** A space fitted to the stages with both groups
pooled asks: *does each group traverse the shared stage program the same way?* Its
blindness to other differences is the scope of that question, not a defect of the
tests. A different question needs a different space.

## The span principle

The latent scores are a linear map of the standardized features, `t = (x − x̄) R / s`,
with `R` the fitted rotations (one column per component). A mean difference `d` between
the groups moves the latent stage means by `dᵀR`. Any part of `d` orthogonal to the
span of `R` moves them by nothing. **The measurement sees a group difference only
through its component inside the fitted basis's span.**

The basis is fitted to separate the stages, both groups pooled, so its span is built
from the pooled stage means. At rank K − 1 with balanced cells it is close to the span
of the centred pooled stage means. The plasmode diagnostic measured this on the SEA-AD
export: the pooled stage displacement keeps 99.94 % of its norm inside the rank-5 span
([plasmode report](reports/plasmode-2026-11-02.md), §4).

### Why an equal-length rotation cancels: the rhombus

Take two stages. Group A's step is `a` and group B's step is `b`, of equal length
`‖a‖ = ‖b‖` but along different features: B's stage program has been rotated onto other
features, which is what the `orientation` surgery does when it relocates stage-changing
CpGs. The pooled step is `(a + b)/2`, the group difference is `b − a`, and

```text
(a + b) · (b − a) = ‖b‖² − ‖a‖² = 0.
```

These are the two diagonals of a rhombus, which are perpendicular. A basis built from the
pooled step contains `a + b` and is orthogonal to `b − a`: in it, both groups' steps
project onto the same vector, so the latent `angle` is 0 and the rotation has cancelled.

**Three or more stages need a stronger condition.** With steps `a_s` and `b_s`, the
difference at one stage must be orthogonal to the pooled step of *every* stage:

```text
(a_t + b_t) · (b_s − a_s) = (b_t · b_s − a_t · a_s) + (a_t · b_s − b_t · a_s) = 0  for all s, t.
```

The cancellation is exact when B's steps keep the lengths and mutual angles of A's (the
first bracket is 0) and are orthogonal to every one of A's steps (the second is 0): B's
program has moved onto features A's program does not use, which is what a relocation
does. A rotation **within** A's own features does not satisfy it and is partly or wholly
visible. The stage-span rotation, made in the plane of A's own configuration, is the
wholly visible case: it keeps 0.85–0.97 of its energy in the span (see the pilot table
below). The latent-space comparison's alignment axis (`span_alignment`) moves one rotation
continuously between these two cases at a fixed intended geometry
([Simulations](api/simulations.md#stage-span-construction)).

The repository has measured the consequence three times:

- The orientation probe's Phase 5 latent angle at full relocation was **6.1°**, against
  a population angle of 90° ([orientation probe](reports/orientation-measurement-probe-2026-10-12.md),
  §4; the test still rejected, because the anchor's latent angle was 0.3°).
- The plasmode's sparse plant kept **2.2 %** of its norm inside the rank-5 span, what a
  random direction keeps (√(5/10,000) = 2.24 %) ([plasmode report](reports/plasmode-2026-11-02.md), §4).
- Orientation/`angle` power tracks the CV-selected rank rather than the effect size
  ([design-guidance pilot](reports/design-guidance-pilot-2026-10-28.md), §3.2 and §3.4).

## Which differences are visible

A difference is visible when it lies inside the span, which means it re-weights the
stage program both groups share:

- **Scaling the program** (`magnitude`): B's stage steps are a multiple of A's, so the
  difference is parallel to the pooled steps and lies in the span. This is why
  `delta` is the well-powered component in every operating study.
- **Rotating or bending the program inside its own span**: the stage-span
  construction (below) makes exactly these edits, and the measurement recovers them in
  the target statistic.
- **A constant offset** (`translation`) is invisible to all three statistics by design.
  It moves the group main effect, which MOTCO does not test, so it is the negative
  control.

A difference whose **direction** is orthogonal to the span is not visible, however large.

**Blindness is a direction, not a set of features.** The question is where the
difference vector points relative to the span, not which features it touches. The same
features can carry both kinds of difference:

- Re-weighting stage-program features along the program is visible.
- Moving them in a direction orthogonal to it is not.
- A difference on features the program does not use is visible to the extent its
  direction correlates with the span, which for an unrelated feature set is about the
  random share. A direction on a few such features keeps several times that share,
  because the sample scaler rescales them, and is still not detected (see the
  visibility profiles below).

The three measured cases are all orthogonal directions. The relocation of
stage-changing sites (the `orientation` and `shape` surgeries) points B's steps away
from A's, and the rhombus places the difference outside the pooled span. A sparse plant
on a few random features (the plasmode) points in an effectively random direction.

## The incidental routes, and their mixed signature

An out-of-span difference can still register when the measurement departs from the
balanced rank-(K − 1) picture:

- **Rank above K − 1.** When CV keeps more components, the extra ones follow
  within-stage covariance and can pick up part of the group difference. In the Phase 6
  study the latent angle at full relocation was 32.1° at rank ≤ 3 and 76.0° at rank > 3,
  with `angle` power 0.082 against 0.506 ([orientation probe](reports/orientation-measurement-probe-2026-10-12.md), §1).
- **Unequal cell sizes.** The pooled means are weighted by the cells, so the rhombus
  becomes a kite and the diagonals are no longer perpendicular. The SEA-AD table
  `[[11, 10, 28], [9, 10, 12]]` is such a design.
- **Unequal feature scales.** Standardization rescales each feature by its pooled SD,
  so a relocation onto features of a different scale changes the steps' lengths.
- **Correlated noise.** PLS weights tilt toward the noise covariance, so the span is not
  exactly that of the stage means.
- **Nonlinear units.** Methylation enters as M-values and drives expression through the
  cascade, so relocated sites need not keep their step length.

Through these routes a rotation registers as a **mix** of `delta`, `angle` and `shape`,
not as `angle` alone. In the design-guidance pilot, the orientation surgery at e = 1.0
moved `delta` (0.12–0.92) and `shape` (0.19–1.00) as well as `angle`
([design-guidance pilot](reports/design-guidance-pilot-2026-10-28.md), §3.2).

## Score coordinates: the span is seen through a shear

The statistics are not measured on an orthonormal projection onto the span. They are
measured in the **PLS score coordinates** that scikit-learn's `transform` returns:

```text
T = X·R,   R = W·(PᵀW)⁻¹,
```

where `W` holds the weights and `P` the loadings. `PᵀW` is unit upper triangular, so `R`
is `W` followed by a shear. The scores span the same subspace as an orthonormal basis of
the weights, but the map is not an isometry of it. **Lengths and angles inside the span
are not preserved**: two configurations equal in the standardized feature space are equal
in score coordinates, but two of the same path length or shape need not stay so.

- **Size is unaffected.** Under the null both groups share one population trajectory.
  Any linear map keeps two equal configurations equal, and RRPP permutes inside the same
  coordinates. The paper-grade anchors and translation controls read 0.002–0.062 against
  0.0695 ([stage-span paper grade](reports/stage-span-paper-grade-2026-11-04.md), §5).
- **Specificity is affected.** A difference that is pure in the standardized feature
  space can move the other statistics in score coordinates. A bend that keeps B's path
  length exactly can still lengthen or shorten it once the span is sheared, and so leak
  into `delta`.
- **The residual is measured.** The [paper-grade leak decomposition](reports/stage-span-paper-grade-2026-11-04.md#6-the-leak-decomposition)
  measured the bend at e = 0.5 under the planted scale, where B's population path length
  equals A's exactly. The table gives the median B − A path-length difference as a share
  of A's path, on the same population rows and through the same fitted basis.

| column | planted-scale coordinates | orthonormal span | PLS scores |
|---|---|---|---|
| n = 200, two blocks | −0.00 % | +0.05 % | +0.03 % |
| n = 80, two blocks | −0.00 % | +0.11 % | −0.02 % |
| SEA-AD table, two blocks | −0.00 % | −0.01 % | **−0.48 %** |
| n = 200, three blocks | −0.00 % | +0.06 % | +0.10 % |
| n = 200, four stages | −0.00 % | +0.06 % | **−0.59 %** |

  The residual is largest on the unequal SEA-AD table and at four stages. The rotation
  shows it more strongly: the 45°
  stage-span rotation on the SEA-AD table reads −3.35 % in scores against +0.66 % in the
  orthonormal span. At paper grade these residuals did not produce a `delta` leak (the
  bend's `delta` reads 0.016–0.036), but they are the one remaining route that the
  construction cannot remove.

### Orthonormal span coordinates: measured, not adopted

**Orthonormal span coordinates** measure on `X_in·Q`, with `Q = qr(R)` an orthonormal
basis of the fitted span (equal to `W` up to column signs), rather than on the scores.
`delta`, `angle` and `shape` are invariant to any orthogonal map, so every orthonormal
basis of the span gives the same statistics and p-values.

MOTCO 0.19.0 offers them in the harness as the opt-in
`SimulationEvaluationParams.latent_coordinates = "orthonormal"`. The setting keeps the
fit, rank and permutation draws. It is not available in `motco analyze`, in study
configurations or with attribution. Without it every output is unchanged.

The [2026-11-09 validation](reports/orthonormal-coordinates-2026-11-09.md) measured
both coordinate systems on the same units: the stage-span paper grade's 32,500 units and
the weak design's 1,000 controls, with the same data, fit and draws, at 500 × 999. Its
verdict is **`power_loss`**.

- **Size holds.** The orthonormal coordinates read 0.000–0.066 at the 60 paper-grade
  controls and 0.010–0.050 at the weak design's 6, none above its bound.
- **The leak is removed.** The watched cells are the 45° rotation's off-target `shape`:
  0.070 → 0.000 at n = 200 and 0.086 → 0.000 at four stages (paired p ≤ 6e-11). The
  bend's planted-scale residual stays within 0.11 % of A's path, against −0.38 % and
  −0.60 % in scores on the SEA-AD table and at four stages. No off-target rate rises, and
  1 of 90 is above the bound (6 are allowed).
- **Target power is lost.** 19 of 45 target cells lose, 16 of them beyond the Bonferroni
  ceiling. All of them are at n = 80, on the SEA-AD table and at four stages, except two
  at n = 200. The largest losses:
  - the SEA-AD table's 45° rotation, `angle` 0.962 → 0.760;
  - the 22.5° rotation at four stages, 0.744 → 0.519;
  - the e = 0.25 bend at n = 80, `shape` 0.724 → 0.518;
  - the e = 0.1 magnitude on the SEA-AD table, `delta` 0.542 → 0.416.

  The three-block column loses nothing.
- **Why the scores have more power is open.** The scores do not see a larger target
  orientation or shape: the population `angle` and `shape` are equal in both
  coordinate systems. They do stretch a planted magnitude by up to 5 %. Where the
  orthonormal coordinates lose, their `angle` null is wider: at the SEA-AD anchor the
  95th percentile is 34.6° against 26.1°.

**Production and `motco analyze` stay on score coordinates.** The scores' shear is a
measured trade-off. It lets a 45° rotation leak into `shape` at up to 0.086, and it buys
up to 0.23 of target power at small, unequal or four-stage designs. An off-target `shape`
or `delta` rejection beside a strong `angle` or `shape` rejection can be the shear's.

Only an `adopt_orthonormal` verdict would have permitted a separate change to adopt
orthonormal coordinates. This run did not return it
([roadmap](roadmap.md#next-changes), item 15).

## Measurement spaces that keep the rotation, and what they cost

Every linear latent space is a projection of the standardized, concatenated features
(`concat`). So it can only lose trajectory information relative to `concat`, never add
it: a projection cannot see a difference `concat` does not contain. It can still have
**more power** than `concat`, because it discards noise directions along with the
differences it cannot see. In the orientation probe's Phase 6 design the shape
construction at e = 1.0 was rejected on `shape` in 0.19 of replicates by `pls_cv` and in
0.00 by `concat` ([orientation probe](reports/orientation-measurement-probe-2026-10-12.md),
§3). Choosing a space chooses which directions of difference count, and how much noise
is measured with them.

The orientation probe crossed five measurement spaces on matched datasets
([orientation probe](reports/orientation-measurement-probe-2026-10-12.md), §3–§6). At
full relocation the population angle is 90°:

| space | Phase 6 latent angle | Phase 5 latent angle | cost recorded by the probe |
|---|---|---|---|
| `pls_cv` (production, net displacement) | 56.7° | 6.1° | — |
| `pca` at rank 2(K − 1) (label-free) | 90.5° | 90.2° | Phase 6 `shape` power 0.19 → 0.00; `angle` leaks into the shape cell (0.39 under PC1) |
| `pls_within_group` (two-sided, fixed basis) | 90.3° | 90.2° | invalid without a refit per permutation: anchor `angle` 0.48–0.88; a refit multiplies the cost by the permutation count |

The label-free and two-sided spaces keep the full rotation because their bases are not
built only from the pooled stage means. Neither replaces the production space: `pca`
loses the shape test at small n, and the two-sided basis needs a refit inside every
permutation to be valid. The stage-PLS space stays the production measurement.

## Measuring visibility: span capture and the stage-span construction

Two opt-in tools of the simulation harness (0.15.0) measure visibility directly:

- **`span_capture`** (`SimulationEvaluationParams(span_capture=True)`, study key
  `evaluation.span_capture`) records, for every replicate, the shares of the
  population group difference (`group_difference`) and of the pooled stage means
  (`pooled_stage`) inside the fitted basis's span, beside the random-direction share
  `rank / n_features`. It also records the population means measured through the
  fitted basis (`population_latent`). The shares are **energy** shares: for one
  direction, the plasmode's 2.2 % norm share is an energy share of 0.05 %, equal to its
  random reference 5/10,000. See [Simulations](api/simulations.md#realized-geometry-decomposition).
- **The stage-span construction** (`construction="stage_span"`) plants group B's
  difference inside the span of B's own stage configuration: a scaling, a rotation in
  the configuration's own plane, a bend at fixed size and direction, or an offset along
  the progression. It is circular by design, because it plants where the measurement
  looks. It answers how sensitive the test is to a difference it can see, not how
  visible the surgeries are. See [Simulations](api/simulations.md#stage-span-construction).

## What a non-rejection means

A non-rejection of `angle` or `shape` means **no detectable difference inside the
stage-supervised span**. It answers the space's question (whether the groups traverse
the shared stage program the same way) and nothing wider. It does not mean the groups
follow the same direction, or have the same shape, in feature space. A rotation of the stage program onto other features,
or a sparse difference on features outside the program, can be large and leave every
statistic at its null. The span-capture shares below are the evidence for each
simulator construction. Report the selected rank beside any `angle` reading, as the
[real-data analysis](api/analysis.md#what-is-measured-and-recorded) record does. The
[split interaction test](#testing-outside-the-span-the-split-interaction-test) asks the
wider question directly.

## Capture shares measured by the stage-span pilot

The 2026-11-04 pilot ([report](reports/stage-span-pilot-2026-11-04.md), 100 replicates
per cell) measured every mode under both constructions. It ran at three stages on
methylation + expression, at n = 200, n = 80 and the SEA-AD table. The table gives the
median `group_difference` energy share. A random direction keeps rank / 498 = 0.004–0.006.

| mode | effects | feature-set surgery | stage-span construction |
|---|---|---|---|
| `magnitude` | 0.02–0.1 | 0.934–0.970 | 0.934–0.970 |
| `orientation` | 0.25, 0.5 | 0.002–0.008 | 0.924–0.971 |
| `orientation` | 1.0 | 0.004–0.006 (0.215 on the SEA-AD table) | 0.849–0.946 |
| `shape` | 0.25–1.0 | 0.004–0.011 | 0.917–0.969 |
| `translation` | 0.25–1.0 | null (constant over stages) | null |

The pooled stage means keep 0.93–0.97 under every construction.

- **The surgeries' orientation and shape plants are invisible.** They sit at the
  random-direction share, so their low `angle` and `shape` power measures visibility,
  not the tests.
- **Given a visible difference, the tests are sensitive.** The stage-span 45° rotation
  and half-step bend are rejected in 0.95–1.00 of replicates at n = 80 and above. The
  surgeries reach 0.06–0.40 at the same effects.
- **Through the fitted basis the stage-span rotation reads close to its intended
  value**: 45° reads 43.7–45.3°. The surgery's full relocation reads 4.4–4.8° at the
  balanced columns, against a population 90°.
- **The only partial exception is the kite route.** On the unequal SEA-AD table the
  full relocation keeps 21.5 % and reaches `angle` power 0.67, registering as a mix
  (`delta` 0.31, `shape` 0.52).

The retrofit check found that an edit made in a replicate's own fitted basis and pulled
back to the features has cosine 0.965–0.984 with the population stage-span plant. A
basis-defined plant and the population-defined one nearly coincide.

## Visibility profiles across alignment

The 2026-11-04 latent-space comparison ([report](reports/latent-space-comparison-2026-11-04.md),
100 replicates × 199 permutations per cell) measured five spaces on identical data. It
moved a 45° stage-span rotation from inside the stage span (φ = 0) to orthogonal to it
(φ = 90) at a fixed intended geometry (`span_alignment`). The designs had three stages on
methylation + expression: n = 200 and the SEA-AD table at the generator's stage effect
(δ = 2.0), and the SEA-AD table at a δ calibrated to SEA-AD's held-out stage AUROC
(δ = 0.11, AUROC 0.617). The table gives the median `group_difference` capture and the
`angle` rejection rate at φ = 0 / 30 / 60 / 90, at n = 200 and on the SEA-AD table.

| space | capture, n = 200 | `angle`, n = 200 | capture, SEA-AD | `angle`, SEA-AD |
|---|---|---|---|---|
| `pls_cv` (production) | 0.97 / 0.82 / 0.40 / 0.008 | 1.00 / 1.00 / 0.79 / 0.08 | 0.93 / 0.79 / 0.41 / 0.026 | 0.96 / 0.73 / 0.43 / 0.11 |
| `pca_matched` (PCA at the `pls_cv` rank) | 0.97 / 0.85 / 0.22 / 0.006 | 0.89 / 0.79 / 0.26 / 0.05 | 0.93 / 0.83 / 0.45 / 0.032 | 0.72 / 0.54 / 0.10 / 0.06 |
| `pca_fixed` (PCA at rank 4) | 0.97 / 0.86 / 0.45 / 0.024 | 1.00 / 1.00 / 0.67 / 0.09 | 0.94 / 0.83 / 0.47 / 0.051 | 0.75 / 0.53 / 0.15 / 0.09 |
| `pls_residual_pca` (`pls_cv` + 2 residual axes) | 0.97 / 0.83 / 0.41 / 0.024 | 1.00 / 1.00 / 0.54 / 0.17 | 0.93 / 0.81 / 0.45 / 0.075 | 0.93 / 0.79 / 0.35 / 0.18 |
| `concat` (all 498 features) | 1 / 1 / 1 / 1 | 1.00 / 1.00 / 1.00 / 1.00 | 1 / 1 / 1 / 1 | 0.81 / 0.85 / 0.89 / 0.90 |

- **Every low-dimensional space has the same profile.** Capture falls to the
  random-direction share at φ = 90, and so does `angle` power. Neither a label-free basis
  nor within-stage residual components recover an out-of-span rotation; only `concat`
  does.
- **Size holds in every space but one.** The anchor and translation controls read at most
  0.08 (bound 0.0936) in every space except `concat`. Its `angle` is liberal at the
  SEA-AD-like stage signal (0.21 and 0.22).
- **Extra directions cost in-span power.** `pls_residual_pca` shares `pls_cv`'s
  coordinates and adds two near-empty ones. Magnitude `delta` falls from 0.74 to 0.37
  (n = 200) and from 0.23 to 0.07 (SEA-AD table), and a 22.5° in-span rotation's
  `angle` from 0.84 to 0.41. `concat` and `pca_fixed` lose as much. `pca_matched` at rank
  2 misses the stage program altogether: its leading axes keep 0.53–0.54 of the pooled
  stage means, against `pls_cv`'s 0.93–0.97.
- **Support at a fixed norm.** A 20-feature out-of-span direction keeps 0.108 (n = 200)
  and 0.116 (SEA-AD table) of its energy in the `pls_cv` span, against 0.008 and 0.026
  dense: differences +0.102 and +0.091, SE 0.002 and 0.003. Its `angle` power is the
  same: −0.02 against SEs of 0.036 and 0.042. So blindness is a matter of direction for
  detection, while capture depends somewhat on the support. Through `concat` the sparse
  rotation reads 32–34° against an intended 45°, which shows the sample scaler shrinking
  the 20 planted features.
- **At SEA-AD's stage signal the span is mostly noise.** At δ = 0.11 the fitted stage-PLS
  span keeps 0.10 of the population stage program (rank 3, random share 0.006), and no
  space has power for any plant.

The advisory reading was `document_tradeoff`: no space gains on the out-of-span rotation
without losing what the production space sees. The production stage-PLS space stays the
measurement, and its question stays the one stated at the top of this page.


## Sensitivity given visibility at paper grade

The 2026-11-04 paper grade ([report](reports/stage-span-paper-grade-2026-11-04.md), 500
replicates × 999 permutations per cell) repeated the stage-span measurement under the
planted scale (`span_scale: "planted"`, which makes each edit exact where the sample
scaler looks) at five columns. Its predeclared decision is **`sensitivity_established`**:
size holds at every column, every stage-span plant keeps at least 0.90 of its energy in
the fitted span, and the off-target rates stay within the bound rule.

Target power of the stage-span plants (`angle` for the rotation, `shape` for the bend,
`delta` for the scaling), at the three effects of each mode:

| column | orientation, 9° / 22.5° / 45° | shape, e = 0.1 / 0.25 / 0.5 | magnitude, e = 0.02 / 0.05 / 0.1 |
|---|---|---|---|
| n = 200, two blocks | 0.74 / 0.88 / 1.00 | 0.70 / 0.97 / 1.00 | 0.23 / 0.69 / 0.94 |
| n = 80, two blocks | 0.40 / 0.41 / 0.99 | 0.21 / 0.72 / 0.97 | 0.04 / 0.24 / 0.62 |
| SEA-AD table, two blocks | 0.38 / 0.38 / 0.96 | 0.21 / 0.68 / 0.97 | 0.06 / 0.19 / 0.54 |
| n = 200, three blocks | 0.98 / 1.00 / 1.00 | 0.91 / 1.00 / 1.00 | 0.43 / 0.88 / 0.99 |
| n = 200, four stages | 0.48 / 0.74 / 1.00 | 0.48 / 0.96 / 1.00 | 0.17 / 0.65 / 0.94 |

- **Capture.** The plants keep 0.92–0.97 of their energy in the span (median per cell).
  Through the fitted basis the rotation reads its intended angle to within 1.3° (SEA-AD,
  45° reads 43.7°).
- **The visibility contrast.** At e = 0.5 the surgeries reach `angle` 0.11–0.49 and
  `shape` 0.08–0.67 at the same columns, keeping 0.2–0.8 % of their energy in the span
  (three-block surgery capture was not measured).
- **Remaining specificity.** The rotation's `shape` rises with the angle: at n = 200 it
  reads 0.030 / 0.046 / 0.070 and at four stages 0.004 / 0.030 / 0.086, at or just above
  the 0.0695 per-cell bound at 45°. The bound rule allows these (3 of 90 rates above it,
  6 allowed, none above the ceiling), but this is the cell to watch in any larger effect.
- **Excluded configurations.** The planted fixed point has no solution when the
  configuration's two in-plane directions nearly tie. It failed for 12 of the 1,000
  four-stage orientation units, which were reported as failures; those two cells read
  over 497 and 491 replicates.

## Which stage-supervised space: weighting and sparsity

The latent-space comparison kept the stage-PLS span as the production measurement. It also
showed how little of the stage program that span keeps at SEA-AD's stage signal (0.10 at
the weak design). Two opt-in settings change the fitted span without changing its stage
objective:

- **Block weighting** (`block_weighting: equal_variance`; [API](api/analysis.md#block-weighting))
  multiplies each standardized block by 1/√p_b. Unweighted, a block's share of the joint
  variance, and so of every stage-PLS direction, grows with its feature count; weighted,
  each block has unit total variance. The fit then runs without per-feature rescaling, so
  `delta` is in the weighted metric. When the blocks have equal feature counts the weights
  are equal, and the only change is that sklearn no longer re-standardizes each CV fold.
- **Sparse stage-PLS** (`sparsity: {scope: block}`; [API](api/analysis.md#sparse-stage-pls))
  lets each component load on a budget of features per block, chosen with the rank by the
  same held-out AUROC. A budget can drop features that carry only noise, at the price of
  any stage signal they carry.

Neither setting can see outside the stage span; both can only change which stage-supervised
span is fitted. The production latent-space selection (2026-11-10,
[report](reports/production-latent-space-2026-11-10.md)) compared them with production on
matched units: the latent-space comparison's three designs and the SEA-AD plasmode export
with a stage-span plant on real noise. Its predeclared decision is **`keep_production`**:

- **No gain where production has none.** At the weak design every arm rejects at
  0.01–0.07 and keeps 0.08–0.10 of the stage program at the anchor. On the plasmode
  plant the largest gain, `sparse_block` at e = 0.5 (`delta` 0.51 → 0.62), misses the
  Bonferroni threshold. `weighted` reads `no_gain`.
- **Size.** Every arm holds size in simulation (0.00–0.06). On the plasmode anchor the
  sparse arms read `delta` 0.11, above the ceiling at 100 replicates, and production reads
  0.11 on the same replicates. The decision counts it as `size_excess` for the sparse arms;
  the paper grade would have resolved it.
- **A gain at a strong stage signal.** At δ = 2.0 every candidate rejects more often than
  production where production is below saturation: `angle` 0.46 → 0.90–1.00 for the
  SEA-AD table's 22.5° rotation. Production's CV picks rank 3 in those cells and the
  candidates rank 2, while every arm keeps 0.92–0.98 of the stage program. The gain
  follows the selected rank, not a better view of the program.

Both settings stay opt-in. A case study may declare block weighting as a sensitivity
analysis, with `delta` in the weighted metric. Sparse stage-PLS costs seven to nine
production fits and has shown no gain that would justify it.

## Testing outside the span: the split interaction test

The production statistics cannot see a difference outside the fitted span, and no
measurement space sees it without losing in-span power (above). The split interaction
test (`motco.stats.interaction.interaction_split_test`; `motco analyze --span-split-test`)
asks the out-of-span question without a space. It takes the group × stage interaction
energy of the standardized features, in the fitted PLS model's input coordinates, and
splits it exactly into three parts, each tested by a trace pseudo-F under Freedman–Lane
permutation:

- **`total`**: is there a group × stage difference anywhere in the measured features?
- **`span`**: is there one inside the fitted span, where `delta`, `angle` and `shape`
  decompose it?
- **`complement`**: is there one outside it, where they cannot see it?

The basis is the fitted stage-PLS span. It never sees the group labels and is never
refitted, so its draws are RRPP's at the same seed.

**"Span" means the fitted span, not the population stage program.** At a strong stage
signal the two nearly coincide: the fitted span keeps 0.90–0.97 of the pooled stage means.
At a SEA-AD-like stage signal the fitted span keeps only 0.10 of the population program
([comparison](reports/latent-space-comparison-2026-11-04.md), §9). There a difference
along the true program can fall outside the fitted span.

**Reading the parts.**

- **A `complement` rejection** means the groups differ by stage with energy outside the
  fitted span. It does not mean one group uses a different molecular program. Part of an
  in-span difference can reach the complement: the 3–7 % of its energy the fitted span
  misses, plus what the sample scaler moves there. A 45° in-span rotation is rejected by
  the complement in 0.35 (n = 200) and 0.32 (SEA-AD table) of replicates. At a weak stage
  signal the complement also contains whatever the fitted span missed of the program.
- **A `span` rejection** means the difference has energy inside the fitted span. There
  `delta`, `angle` and `shape` say how it differs; the span part only says that it does.
  It is less powerful than the targeted statistics: magnitude 0.56 against `delta`'s
  0.74, a half-step bend 0.72 against `shape`'s 1.00 (n = 200).
- **A non-rejection of `complement` or `total`** is uninformative unless power has been
  shown at the design and for the kind of difference. On the SEA-AD blocks, a planted
  difference on 20 of 10,000 features at 2–3 within-stage SDs is rejected by no part
  above the null rate (0.02–0.05), as by `angle`. An energy test spreads over every
  coordinate, so a sparse difference does not register.
- **The in-span share** (`span_share`) is descriptive. It must be read against the share
  a null analysis gives at the same rank and design, never against 0. Components beyond
  K − 1 follow the total covariance, where the noise's interaction energy also sits. At
  the plasmode anchor (rank 5) the null share is 0.24; at a two-component span at n = 200
  it is 0.006.

**Why a ratio.** RRPP permutes the reduced model's residuals, whose mean square is RSS/n,
not RSS/(n − r). An energy statistic carries that shrinkage into its null, and summed over
thousands of coordinates the null concentrates, so a few percent of scale becomes many
null standard deviations. The pseudo-F computes its numerator and denominator from the
same permuted residuals, so the scale cancels. At the plasmode anchor the uncorrected
interaction SS rejects the null complement in 0.268 of replicates; the pseudo-F on the
same draws rejects in 0.012.

**What it does not provide.** No geometry: the split test says whether and where the
groups differ, never how. It does not replace `delta`, `angle` or `shape`, and it does not
find sparse differences: the [per-feature interaction scan](#testing-feature-by-feature-the-interaction-scan)
does. A group-aware space refitted per permutation is roadmap item 16's remaining follow-up
([roadmap](roadmap.md)).

**Validation: `split_validated`** (2026-11-05, [report](reports/span-split-2026-11-05.md),
500 × 999 at the controls, 100 × 999 elsewhere, on the latent-space comparison's and the
plasmode's own units):

- **Size** holds on three simulated designs (the 18 control rates read 0.032–0.058,
  against a per-cell bound of 0.0695) and on the SEA-AD blocks at the case study's
  configuration (total 0.040, span 0.050, complement 0.012).
- **Out-of-span power.** The complement rejects the 45° out-of-span rotation in 1.00
  (n = 200) and 0.68 (SEA-AD table) of replicates, where `angle` reads 0.08 and 0.10.
- **Localization.** As the rotation tilts out of the span, the span part falls
  (1.00 → 0.52) and the complement rises (0.35 → 1.00).
- **At SEA-AD's stage signal** no part has power for any plant, and the complement holds
  size there (0.05–0.06), where `concat`'s `angle` did not.

## Testing feature by feature: the interaction scan

Every test above pools the group difference over many dimensions: `delta`, `angle` and
`shape` over the fitted span, the split test's parts over their coordinates. A difference
confined to a few features does not register in any of them. The SEA-AD plasmode's
20-feature plant was rejected at the null rate by all six. The per-feature interaction scan
(`motco.stats.interaction.interaction_feature_scan`; `motco analyze --interaction-scan`)
tests the group × stage interaction of each standardized feature on its own:

- **Per feature**, the interaction F with the covariate-adjusted full and reduced
  trajectory models, the effect RRPP and the split test test.
- **Globally**, the permutation p-value of the largest F over the features.
- **Feature by feature**, a Westfall–Young step-down max-T adjusted p-value, which controls
  the familywise error rate over all scanned features. A **hit** is a feature with adjusted
  p ≤ 0.05.

The draws are RRPP's Freedman–Lane draws at the analysis's seed, with its residual
correction. The scan uses neither the latent space nor the bootstrap, so latent adjustment
and attribution leave it unchanged.

**Reading a hit.** A hit says that this feature's covariate-adjusted group difference
changes across stages, with the familywise error rate over all scanned features controlled
at 0.05. The per-stage differences in `interaction_scan.csv` (the second group minus the
first, at each declared stage, in standardized units) show how. A hit does not say that
the feature drives a latent trajectory difference, and correlated neighbours of a
differing feature can be hits too.

**Reading a non-rejection.** A scan without hits is uninformative unless power has been
shown at the design and effect. At the plasmode's design (n = 80, six covariates, 10,000
features), a ramp to 2 within-stage SDs on 20 features is detected in 0.692 of replicates,
and the median scan names one of the 20. A ramp to 3 SDs is detected in every replicate,
and half the 20 are named.

**Why the maximum, not the sum.** The split test's complement sums interaction energy over
thousands of coordinates, and the noise in that sum grows with their number, so a few
features' energy is lost. A maximum pays only for the search: its threshold grows like the
logarithm of the feature count. The scan therefore targets differences on few features,
and the complement targets differences spread over many. The simulation shows the division
([report](reports/interaction-scan-2026-11-09.md), §6). The scan is weak where a difference
spreads thinly: magnitude 0.09 against `delta`'s 0.74, and a weak bend on the SEA-AD table
0.20 against `shape`'s 0.70. With 498 features the search is cheap, so the scan also matches
or beats the complement on a 45° out-of-span rotation (1.00 / 0.96 against 1.00 / 0.68).
The complement still wins at half that rotation on the SEA-AD table (0.31 against 0.16).

**The scan, attribution and the split test make different claims.**

- **Attribution** ranks features by their contribution to a latent difference that was
  already detected, and makes no error-rate claim.
- **The scan** tests each feature's interaction directly, under familywise control, without
  a latent space. On the plasmode its 20 largest F hold 0.473 / 0.887 of the plant at
  2 / 3 SDs, against attribution's 0.227 / 0.432. The plasmode's attribution verdict
  (`attribution_leg_failed`) is unchanged: the scan is a new test with its own predeclared
  rule, not a rescue of attribution's.
- **The split test** says whether the groups differ by stage inside or outside the fitted
  span; the scan says on which features.

**What it does not provide.** No geometry: the scan names features whose stage pattern
differs between the groups, never how the trajectories differ. It is not a driver analysis.
Its F is the plain per-feature F, not a variance-moderated one. Its familywise threshold is
common to all features, and on the SEA-AD export one block's null is heavier: ATAC gives
0.618 of the null maxima with half the features, so the threshold leans on ATAC.

**Validation: `scan_validated`** (2026-11-09, [report](reports/interaction-scan-2026-11-09.md),
on the split-test validation's own units, which it reproduced exactly):

- **Size** holds on three simulated designs (the six control rates read 0.048–0.058,
  against a per-cell bound of 0.0695) and on the SEA-AD blocks at the case study's
  configuration (0.042).
- **Familywise control.** Under the plasmode's partial null, a feature outside the planted
  set is a hit in 0.042 of replicates at 2 and at 3 SDs.
- **Sparse gain.** At 3 SDs the scan rejects in 1.000 of replicates, against 0.016–0.048 for
  `delta`, `angle`, `shape` and the three split parts on the same units. At 3 SDs, 0.995 of
  its hits are planted features.
