# PLS-DA

Partial Least Squares Discriminant Analysis with double cross-validation.

## Scale and nuisance residualization

`plsda_doubleCV` takes three keyword-only arguments whose defaults reproduce its
earlier behaviour exactly:

- `scale` (default `True`) is passed to every `PLSRegression` fit in the inner
  loop, the outer loop and the per-repeat refits. `False` fits in the metric the
  predictors are already in, for example pooled-standardized features that must
  keep their pooled scale.
- `nuisance` (default `None`) takes a `NuisanceDesign`: row-aligned
  conditioning columns `D` (no intercept) and covariate columns `C`. Every
  training set (each outer rest set and each inner training set) regresses its
  predictors on `[1, D, C]`, subtracts the fitted covariate part `C γ`, and
  residualizes its held-out rows as `X_held_out − C_held_out γ`. The held-out
  rows' conditioning (their labels) is never read. The per-repeat refits use
  coefficients estimated on all rows. Every training design must have full
  column rank, or the function raises a `ValueError` naming the repeat and
  fold, for example when a factor level is absent from a training set. The
  result gains `"nuisance_min_training_df"`, the smallest training-set residual
  degrees of freedom.
- `strata` (default `None`) names row-aligned labels both CV levels stratify on
  instead of `y`; the AUROC still scores `y`. Folds depend only on the
  partition the labels define, so labels that partition the rows like `y`
  reproduce the default folds. With a nuisance design, pass
  `nuisance.cells()` (its conditioning cells) so every training set contains
  every cell: the inner folds are not shuffled, and on rows sorted by group
  within stage a stage-only stratification can put a whole group × stage cell
  into one validation fold.

The real-data analysis uses both together for covariate-aware latent
construction; see [Real-data analysis](analysis.md#latent-adjustment).

## Functions

::: motco.stats.pls.plsda_doubleCV

::: motco.stats.pls.calculate_vips

## Nuisance residualization

::: motco.stats.pls.NuisanceDesign

::: motco.stats.pls.fit_nuisance_coefficients

::: motco.stats.pls.residualize

::: motco.stats.pls.nuisance_partial_r2
