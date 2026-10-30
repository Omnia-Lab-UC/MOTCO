# Trajectory Analysis

Group trajectory differences: magnitude, orientation, and shape.

## Assumptions about the measurement space

The trajectory statistics use Euclidean geometry in the supplied outcome
matrix. Their interpretation therefore depends on how that matrix was built.

| Representation | What the statistics directly describe | Interpretation |
|---|---|---|
| Standardized concatenation | The complete standardized feature space | A transparent baseline, but feature weighting depends on standardization and block dimensionality |
| PCA or PLS scores | A linear projection of the feature space | Path length and direction can be related to the original features through the fitted loadings; omitted components still make the projection lossy |
| SNF spectral coordinates | A nonlinear embedding of a fused affinity graph | Euclidean length, angle, and Procrustes shape describe the embedding and are not guaranteed to preserve feature-space geometry |

`delta` and `angle` are most directly aligned with linear measurement spaces.
When using SNF, validate the intended geometric interpretation by simulation
and see [Interpreting trajectory geometry in SNF space](snf.md#interpreting-trajectory-geometry-in-snf-space).
Graph-native diffusion, connectivity, or transition-profile statistics would
be a separate methodology and are not implemented by the current trajectory
API.

## Interpreting shape

`shape` is a strict geometric-morphometric residual configuration statistic.
For each group trajectory, MOTCO centers the stage configuration and scales it
to unit centroid size. Pairwise shape distance is then the residual norm after
aligning one trajectory to the other over the **full orthogonal group** — the
alignment is unconstrained in determinant, so reflections are aligned away
rather than retained.

Under this contract, translation, positive uniform scale, and any orthogonal
transformation — proper rotation or reflection — do not produce nonzero `shape`
distance beyond numerical tolerance. Genuine configuration changes, such as
moving an interior stage so relative stage-to-stage distances change after size
normalization, remain positive.

The reflection policy is uniform: it does not depend on whether the stage
configuration spans the ambient space, so the statistic means the same thing at
every measurement checkpoint. (A trajectory of `k` stages has rank at most
`k − 1`, so in feature space a proper-rotation constraint would be satisfiable
for free in the null space and therefore vacuous — the previous "reflections
retained" contract bound only in a full-rank latent space.) A genuine mirror
difference is an orientation difference and is reported by `angle`.

## Interpreting orientation and identifying feature drivers

In a shared linear PLS space, a trajectory orientation difference means that
the groups progress through different relative combinations of the molecular
patterns retained by PLS. Different stage-associated feature sets are one
possible cause, but not the only one: the same features can change in different
proportions or signs, or one group can recruit an additional correlated feature
module.

For two stages, let each group's standardized feature-space transition be

```text
d_g = mean(X | group=g, stage=1) - mean(X | group=g, stage=0)
```

or use covariate-adjusted LS means in place of raw means. To isolate direction
from magnitude, normalize each transition and calculate

```text
u_g = d_g / ||d_g||
q   = u_B - u_A
```

Large positive or negative values of `q[j]` identify features that contribute
relatively more to one group's progression direction. With more than two
stages, calculate this contrast for each adjacent transition and for the signed
principal feature-space orientation used to summarize the complete path.

A production attribution analysis should:

1. Keep the pooled preprocessing and PLS model fixed; do not fit independently
   aligned latent spaces for the two groups.
2. Reconstruct each group's stage-to-stage PLS score displacement through the
   fitted X loadings, then compare the reconstructed normalized feature steps.
3. Compare that loading-derived contrast with the directly observed or
   LS-mean-adjusted feature contrast `q`. Their agreement describes the part of
   the group difference captured by the retained PLS space; their residual
   describes feature variation that PLS omitted.
4. Bootstrap samples to assess feature-rank and sign stability, especially for
   correlated features, and summarize stable signals at module/pathway level.
5. Report standardized effects alongside original-unit effects. PLS geometry
   is defined after preprocessing, so raw-unit contributions need not have the
   same ranking.

Pooled VIP scores alone are insufficient for this question: they identify
features important to the pooled stage model, not features whose normalized
stage association differs between groups. Orientation drivers require the
group-specific transition contrast above.

For a production attribution result from a detected orientation difference,
use the shared fitted PLS model and the aligned pooled-standardized matrix
with [Orientation Driver Attribution](attribution.md). The attribution API
supports adjacent transitions for multi-stage data, optional covariate-
adjusted mean tables, frozen-model bootstrap stability, and caller-supplied
module/pathway labels. It reports standardized and optional original-unit
effects separately and does not produce p-values, causal claims, or a unique
inverse of the measured features.

## Design matrix

::: motco.stats.design.get_model_matrix

::: motco.stats.design.build_ls_means

::: motco.stats.design.center_matrix

## Estimation

::: motco.stats.trajectory.estimate_betas

::: motco.stats.trajectory.get_observed_vectors

::: motco.stats.trajectory.estimate_difference

::: motco.stats.trajectory.pair_difference

## Permutation test

::: motco.stats.permutation.RRPP

::: motco.stats.permutation.residual_correction_scale
