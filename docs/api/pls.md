# PLS-DA

Partial Least Squares Discriminant Analysis with double cross-validation.

## Scale and nuisance residualization

`plsda_doubleCV` takes five keyword-only arguments whose defaults reproduce its
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
- `min_components` (default `1`) is the smallest number of latent variables the
  inner CV considers: the candidates are `min_components … max_components − 1`,
  so every outer-fold selection, every modal `LV` in the table and every
  refit model has at least that many components. It must satisfy
  `1 ≤ min_components < max_components`, or the function raises a
  `ValueError` naming both values before any fold is fitted. The default is
  the unfloored rule; the trajectory latent space
  (`fit_pls_latent_space`, shared by the simulation harness and
  `motco analyze`) applies a floor of 2 by default, see
  [Real-data analysis](analysis.md#what-is-measured-and-recorded).

- `sparsity` (default `None`) takes a `SparsityDesign` and selects a sparsity
  budget jointly with the rank; see [below](#sparse-pls-da).

The real-data analysis uses `scale`, `nuisance` and `strata` together for
covariate-aware latent construction; see
[Real-data analysis](analysis.md#latent-adjustment).

## Sparse PLS-DA

`SparsePLSDA` (`motco.stats.sparse_pls`, also importable from `motco.stats`)
is a `PLSRegression` whose components each load on a budget of features (the
sPLS-DA of Lê Cao et al., regression mode). Inside scikit-learn's NIPALS power
iteration the X weight vector is soft-thresholded before it is normalized: the
largest absolute weights are kept and shrunk toward zero by the largest
discarded one, and every other weight is zero. Everything else is
scikit-learn's: the model inherits `transform`, `inverse_transform`,
`predict`, `x_rotations_`, `x_scores_` and the input scaling, so every
consumer of a fitted PLS-DA model in MOTCO (the latent map, span capture,
attribution, VIPs, the split test) accepts it unchanged. At the full budget
the fit is `PLSRegression.fit`, bit for bit. `keep_` holds the count per
component and group, `support_` the non-zero weights (features ×
components). The fit agrees with mixOmics `splsda` on the committed fixture
(`tests/data/splsda_reference_fixture.npz`, produced once by
`tests/data/reference/splsda_fixture.R`): identical supports, weights and
scores to 1e-6 up to column signs.

A `SparsityDesign(scope, keep_fractions, feature_groups)` describes the budget
to select over:

- `scope` is `"joint"` (one budget over all features) or `"block"` (one budget
  per block, so every block keeps its share whatever its size; a joint budget
  on blocks of unequal noise can empty a block). DIABLO's between-block
  agreement term is not part of it.
- `keep_fractions` (default `1.0, 0.5, 0.2, 0.1, 0.05`) are fractions of the
  scope's feature count, distinct and in (0, 1], held fullest first. A
  fraction keeps `max(1, round(f · p))` features per group (halves round up),
  so every block keeps at least one feature per component.
- `feature_groups` are the block sizes in column order; the block scope needs
  them, and `plsda_doubleCV` rejects groups that do not cover the columns of
  `X` before any fold is fitted.

With `sparsity=`, the inner CV scores every (rank, fraction) pair by mean
held-out AUROC, ordered rank ascending then fraction descending, so an exact
tie goes to the smaller rank and, within it, to the fuller budget: the CV
leaves production only on strictly better held-out AUROC. The outer fold is
evaluated at the selected pair. Each repeat aggregates the rank by mode
(smaller on ties) and the fraction by mode (fuller on ties), the table gains a
`keep_fraction` column after `LV`, and every refit model is a `SparsePLSDA`
at the repeat's pair. A ladder of `[1.0]` reproduces the run without
sparsity exactly. Without `sparsity` the output is byte-identical to before
(`tests/data/plsda_doublecv_identity.json`) and has no new column.

## Functions

::: motco.stats.pls.plsda_doubleCV

::: motco.stats.pls.calculate_vips

## Sparse PLS-DA estimator

::: motco.stats.sparse_pls.SparsePLSDA

::: motco.stats.sparse_pls.SparsityDesign

## Nuisance residualization

::: motco.stats.pls.NuisanceDesign

::: motco.stats.pls.fit_nuisance_coefficients

::: motco.stats.pls.residualize

::: motco.stats.pls.nuisance_partial_r2
