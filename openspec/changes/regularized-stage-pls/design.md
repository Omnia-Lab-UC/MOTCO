# Design

## Context

See proposal.md (Why) for the motivation. The facts below shape the approach.

- **The surface MOTCO reads off the PLS model.** Across `evaluation.py`, `diagnostics.py`, `analysis.py`, `stats/attribution.py` and `attribution_diagnostics.py` the fitted model is used through `transform`, `inverse_transform`, `x_rotations_`, `x_scores_`, `_x_mean`, `_x_std`; `calculate_vips` reads `x_weights_`, `x_scores_`, `y_loadings_`; `_plsda_auroc` calls `predict`. Nothing reads sklearn-private state beyond `_x_mean` / `_x_std`.
- **sklearn's `PLSRegression`** (NIPALS, regression-mode deflation) gives `x_rotations_ = W·(PᵀW)⁻¹`, `x_scores_ = X_in·x_rotations_`, `coef_` from the rotations and the Y loadings, and under `scale=True` divides columns by their ddof-1 standard deviation. `fit_plsda_model` fits it with `max_iter=1000` on one-hot stage labels.
- **The double CV** (`plsda_doubleCV`) iterates `n_lv_candidates = range(min_components, max_components)` in the inner loop, takes `argmax` (first occurrence, so the smaller rank wins ties), evaluates the outer fold at the selected rank, aggregates per repeat by `_modal_int_with_parsimony`, and refits a `PLSRegression` per repeat. Its default output is pinned byte for byte (`tests/data/plsda_doublecv_identity.json`).
- **Block membership is known at both call sites.** `fit_pls_latent_space` receives `layer_feature_counts` in concatenation order from the harness (`_pls_integration`) and from `analysis.analyze_trajectories`; the joint columns are `block__feature`.
- **Sparse PLS-DA (mixOmics `splsda`, Lê Cao et al. 2008, 2011).** Per component: initialize from the leading singular vectors of XᵀY; iterate the X weight (soft-thresholded to `keepX` features: keep the largest absolute entries, zero the rest, shrink the survivors by the largest discarded magnitude) and the Y weight (dense) to convergence; deflate X by its loading and Y by regression. With `keepX = p` the threshold is zero and the iteration is NIPALS PLS2.
- **Multi-block PLS.** With block scaling and no between-block term, multi-block PLS super-scores equal PLS on the block-scaled concatenation (Westerhuis, Kourti & MacGregor 1998). The DIABLO design matrix adds a between-block covariance term to the objective.

## Goals / Non-Goals

**Goals:**
- One estimator that is production at the full budget and sparse below it, usable by every existing consumer without a branch on its type.
- Budget selection inside the existing CV, with the fold structure, strata, repeats and parsimony rules of the rank selection.
- Byte identity while absent, in the CV and both entry points.

**Non-Goals:**
- Any run or verdict (proposal, Non-goals).
- A between-block agreement objective, Y-side sparsity, or per-component budgets.
- Performance work beyond what the selection's smoke demands.

## Decisions

**D1. The estimator is in-house, with `PLSRegression`'s surface.**
- sklearn has no sparse PLS, and mixOmics is R. `SparsePLSDA` implements the algorithm in numpy and exposes the attributes and methods in Context, with the same shapes and the same conventions (`x_rotations_ = W·(PᵀW)⁻¹`; `transform = inputs·x_rotations_`; `inverse_transform = T·Pᵀ·std + mean`; `predict` from the rotations and `y_loadings_` plus the Y mean; `_x_std` ones under `scale=False`). It also exposes `keep_` (resolved counts per component and group) and `support_` (a boolean features × components matrix).
- **Signs.** Each component's weight is signed like sklearn's (`svd_flip` convention), so the full-budget agreement is within tolerance rather than only up to sign in practice; the contract stays "up to column signs", because `delta`, `angle`, `shape` and every p-value are invariant to them.
- **Alternatives.** Wrapping mixOmics through rpy2 (R at runtime, against the repository's R-free rule); subclassing `PLSRegression` and patching its private NIPALS (binds to sklearn internals pinned at 1.8.0, which the `_x_mean` dependence already strains).

**D2. The budget is a fraction, the scope is joint or block, and the block scope has no agreement term.**
- **Fractions.** A count ladder does not transfer between the simulations (367 + 131 features) and SEA-AD (5,000 + 5,000); a fraction of the scope's count does. The default ladder 1.0, 0.5, 0.2, 0.1, 0.05 includes production and spans the range where sparse PLS literature finds its gains; a ladder is validated as distinct fractions in (0, 1].
- **Block scope.** `max(1, round(f·p_b))` per block per component. It is the sparsity of multi-block sPLS-DA, and it addresses the imbalance the interaction scan measured (ATAC held 0.618 of the null maxima with half the features): a joint budget on blocks of unequal noise can empty a block.
- **No agreement term.** It changes the objective away from stage discrimination; the block-scaled, term-free multi-block PLS is PLS on the weighted concatenation, which `block-weighting` supplies. The selection measures weighting × sparsity as a factorial for exactly this reason.
- **Block membership** comes from `layer_feature_counts` in its order, validated to sum to the column count; it is never parsed from column names.

**D3. The CV selects (rank, fraction) jointly; the absent path is untouched.**
- `plsda_doubleCV(..., sparsity: SparsityDesign | None = None)`. The candidate list becomes pairs, ordered rank ascending then fraction descending, so the existing first-occurrence `argmax` gives the smaller rank and, within it, the fuller budget on ties; the per-repeat fraction is the mode with the fuller budget on ties. The estimator factory (`PLSRegression` or `SparsePLSDA` at the pair) is the only other branch. With `sparsity=None` the candidate list, the fits and every array are those of today, so the identity fixture holds; a ladder of `[1.0]` reproduces production within tolerance (spec).
- **Why fuller on ties.** Production is the validated default; the CV should leave it only on strictly better held-out AUROC. mixOmics's "smallest keepX not significantly worse" rule moves the other way and would make the sparse arm differ from production even when the data say nothing.
- **Why joint rather than sequential (rank first, then budget).** The best rank differs by budget: a sparser fit often wants a different rank. The grid is small (|ranks| × 5) and every candidate is cheap.
- **Alternative: a separate `splsda_doubleCV`.** It would duplicate the fold machinery and the nuisance handling, which the latent adjustment relies on.

**D4. The core reads `integration_params["sparsity"]`; everything else composes.**
- `fit_pls_latent_space` builds the `SparsityDesign` (scope, ladder, groups from `layer_feature_counts`), passes it to the CV, fits the final `SparsePLSDA` at the selected pair with the path's `scale` (True by default; False under nuisance or weighting), and records `cv_params["sparsity"]`, `selected_keep_fraction` and the support sizes. `LinearLatentMap.from_pls` / `from_pls_orthonormal`, `pls_residual_pca_from_pls`, `pls_orthonormal_from_pls`, span capture and attribution work on the shared surface. Forced rank: a one-element ladder gives the budget without a CV; a longer ladder has nothing to select it, so it is rejected.
- **Rejections before preprocessing.** Methods other than `pls`/`pls_residual_pca`; malformed scope, ladder or groups. The orthonormal-coordinates validation (`_validate_latent_coordinates`) is the pattern.
- **Why an integration parameter.** As for block weighting: it describes how the space is fitted, is hash-omitted while absent, and is an axis for free.

**D5. The analysis mirrors the harness, with a typed setting.**
- `SparsitySettings(scope, keep_fractions=DEFAULT_KEEP_FRACTIONS)` on `TrajectoryAnalysisParams.sparsity`, validated in `__post_init__`, popped from `to_dict()` while `None`; the key goes into the analysis's `integration_params`; the record block `sparsity` after `latent`. The CLI takes the scope and an optional ladder.

**D6. The reference fixture is produced once in R.**
- `tests/data/reference/splsda_fixture.R` runs mixOmics `splsda(X, Y, ncomp = 3, keepX = c(k, k, k), scale = FALSE, mode = "regression")` on a seeded Gaussian X with three classes and writes X, Y, `loadings$X` and `variates$X` to CSV; `tests/splsda_fixture_capture.py` packs `tests/data/splsda_reference_fixture.npz`. The flake gains `rPackages.mixOmics`. The test compares support exactly and values to 1e-6 up to column signs. If mixOmics cannot be built in the flake, the fixture is produced in any R with mixOmics and its provenance (versions, script hash) is recorded in the npz; the test itself never needs R.

## Risks / Trade-offs

- **[Full-budget agreement with sklearn is only to tolerance, not bitwise]** → The contract is 1e-6 up to signs; MOTCO's statistics are sign-invariant and the production path still fits `PLSRegression`, so production records are unchanged bit for bit. Only a `[1.0]` ladder measures the tolerance.
- **[Rank selection in the sparse grid inflates compute by the ladder length]** → Accepted; each fit is a few NIPALS iterations at O(n·p) per component. The selection's smoke measures it before any array.
- **[A sparse fit on a weak stage signal selects an unstable support]** → The record carries the per-component support sizes and the CV table carries the selected fractions per repeat, so instability is visible; the selection reads the rank, the fraction and their dispersion per cell.
- **[The block scope keeps noise features in a block with no stage signal]** → By design (every block keeps its share); the joint scope is the alternative, and the selection measures both.
- **[Attribution's reconstruction through a sparse `x_loadings_`]** → The loadings `P = Xᵀt/(tᵀt)` are dense even when the weights are sparse, as in ordinary sparse PLS; the reconstruction is defined and the test checks it runs and is in original units.

## Migration Plan

Nothing migrates. The setting is absent everywhere; rollback is removing the keyword and the module, which nothing reads while absent.

## Open Questions

- Whether `rPackages.mixOmics` builds in the flake's R environment is checked at the fixture task; D6 names the fallback.
