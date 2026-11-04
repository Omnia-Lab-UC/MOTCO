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
  random share.

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
[real-data analysis](api/analysis.md#what-is-measured-and-recorded) record does.

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
