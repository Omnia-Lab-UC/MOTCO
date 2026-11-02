# Orientation Driver Attribution

Orientation attribution explains a two-group, multi-stage difference after an
upstream trajectory analysis has identified the shared PLS space of interest.
Pass the pooled, aligned feature matrix, matching group/stage metadata, and
the one fitted PLS estimator to `analyze_orientation_attribution`.

The analysis orders groups and stages deterministically unless explicit levels
are supplied. For every adjacent stage pair it reports each group's raw
transition and path length, the unit directional contrast, the same contrast
after PLS score projection and inverse transformation, and the signed residual
`observed - pls_captured`. Observed feature values remain in the result; the
reconstruction is a parallel, lossy view rather than a replacement.

## Two decompositions

Attribution reports two decompositions, and they explain different quantities.

The **principal-orientation** block decomposes the tested `angle` estimand:
principal-axis divergence, the contrast between the groups' signed leading
principal axes. It applies the same functional the statistic is built from
(`motco.stats.trajectory.principal_orientation`) to each group's stage-mean
configuration, observed and PLS-reconstructed, so the explanation and the test
cannot drift apart.

The **per-adjacent-transition** blocks describe where along the trajectory the
groups' step directions diverge. They are a per-step description, not the tested
estimand; the two decompositions coincide only for a two-stage design, where the
configuration is rank one and PC1 is the transition direction.

Both appear in the tabular outputs under the `transition_id` column — the
transitions as `from->to`, the principal block under `principal_component_id`
(default `principal`), which is validated against the stage-derived identifiers
so the two can never collide.

## Degeneracy and availability

A principal axis exists for any configuration with variance, but it is only
resolvable when one axis dominates. The result reports each contributing
configuration's relative eigengap and flags the principal contrast `degenerate`
when any falls below `eigengap_threshold` (default `0.05`). The flagged contrast
is still returned: degeneracy is a property of the data, so callers stratify on
the flag rather than lose rows.

That is distinct from unavailability. A trajectory whose net displacement from
first stage to last is at or below `zero_tolerance` is closed, and no direction
is defined at all — that group's principal orientation is marked unavailable,
exactly as a zero-length transition is.

## Units and stability

Effects are reported in the input coordinate units, normally pooled
standardized units. Positive per-feature scales can be supplied to add
original-unit columns. These scales only convert reported effects; they never
change PLS geometry.

Bootstrap resampling is independent within every group-by-stage cell and uses
the same fitted model and preprocessing contract in every replicate. The
result records the seed, requested count, valid count, sign stability, and
top-k selection frequency, for the principal-orientation block as well as the
transitions. Each replicate's principal axis is signed by that replicate's own
net displacement, so resampling noise cannot flip the axis; a replicate whose
trajectory closes contributes no contrast and counts as invalid. Bootstrap is
available for arithmetic and covariate-adjusted means (below); a precomputed
mean table is accepted with stability marked unavailable.

## Covariate-adjusted means

Pass `covariates` (numeric design columns, row-aligned with `X`; code factors
as indicator columns) to attribute from covariate-adjusted cell means instead
of arithmetic ones. Each group × stage mean is then the least-squares mean of
the cell-means model with the covariates added, evaluated at a reference
covariate row `c_ref` (`covariate_reference`, by default the column means):

```text
m̂_gs = x̄_gs − γ̂ᵀ (c̄_gs − c_ref)
```

where `x̄_gs` and `c̄_gs` are the cell's feature and covariate means and `γ̂`
holds the covariate coefficients estimated within cells. This is the model the
real-data analysis fits in the latent space, here fitted in feature space, so
projecting the adjusted means through the fitted PLS model gives the analysis's
latent stage means. The reconstructed and residual components follow from the
adjusted means as they do from arithmetic ones.

- **Reference-row invariance.** `c_ref` shifts every cell by the same vector
  `γ̂ᵀ c_ref`, so transitions, path lengths and both contrast blocks do not
  depend on it. Only the means tables move.
- **Bootstrap refit.** Each replicate draws the same within-cell resample as
  attribution without covariates for the same seed, refits `γ̂` on the
  resampled rows and evaluates the means at the fixed full-sample `c_ref`. The
  fitted preprocessing and PLS model stay frozen, so stability includes the
  uncertainty of estimating the covariate coefficients.
- **Rank-deficient replicates.** Covariates collinear within cells (for
  example, constant within every cell) are an error on the full sample. In a
  resample, for example one that omits the only sample carrying a rare factor
  level, they make the replicate invalid. Such a replicate contributes no
  contrast to any block, is excluded from every valid count, and is counted in
  the configuration row `bootstrap_rank_deficient_replicates`.

The configuration records `mean_source = covariate_adjusted`, the covariate
columns in order and the reference row. The interpretation table adds an
`adjustment_boundary`: adjusted effects assume an additive linear covariate
model evaluated at the reference row, and adjustment does not remove
unmeasured confounding. `covariates` cannot be combined with a precomputed
mean table. Without covariates every output is unchanged.

## Labels and interpretation

Feature-to-module or feature-to-pathway labels are optional and must be a
complete one-to-one mapping supplied by the caller. Aggregate tables preserve
the labels and mark their source as `caller-supplied`. Aggregation is not a
pathway annotation service and does not establish biological causality.

The result contains no p-values or significance decisions. Attribution should
be run after the calling workflow has made any global orientation inference
decision. It describes associations implied by one frozen shared PLS model;
the inverse reconstruction is not unique, and bootstrap stability is
conditional on that fitted representation.

## Functions

::: motco.stats.attribution.analyze_orientation_attribution

::: motco.stats.attribution.attribution_frames

::: motco.stats.attribution.write_attribution_outputs

