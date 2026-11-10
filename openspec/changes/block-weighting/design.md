# Design

## Context

See proposal.md (Why) for the motivation. The facts below shape the approach.

- **One preprocessing seam.** `preprocessing.fit_block_scalers` / `transform_blocks` are called by the harness (`fit_omics_preprocessor`, `FittedOmicsPreprocessor.transform_dataset` / `transform_population`) and by `analysis.analyze_trajectories` (directly, on the declared blocks). A `BlockScaler` holds `feature_names`, `mean` and `scale`; the transform is `(values − mean) / scale`.
- **`PLSRegression(scale=True)` undoes a per-block scalar.** sklearn divides every column by its own standard deviation before fitting, so a column multiplied by *w_b* comes back unchanged. `fit_pls_latent_space` uses `scale=True` by default and switches to `scale=False` under latent adjustment, where the residuals must keep the pooled scale. Under `scale=True` the input scaling is the uniform factor √(n/(n − 1)) on pooled-standardized features (orthonormal-coordinates design, D1), so `scale=False` on unweighted features changes `delta` by that factor and nothing else: PLS is equivariant to a uniform rescaling of X, so weights, rank selection, held-out AUROC, `angle`, `shape` and every p-value are unchanged.
- **Label-blind steps do not touch the null.** The scalers are fitted once on the pooled sample and RRPP permutes residuals in the fitted space. A weight that depends only on the block's feature count or on the pooled standardized block's spectrum has the same status (roadmap item 2, scope boundary).
- **What reads the scales downstream.** `_joint_original_scales(scalers, order, columns)` feeds attribution's reconstruction back to original units and `latent_artifacts.original_scales`; `LinearLatentMap.from_pls` reads the model's own `_x_mean` / `_x_std`, which are the identity under `scale=False`.
- **Integration parameters are already axis-able and hash-safe.** Keys absent from `integration_params` are absent from cell ids, signatures and records; a study configuration may name `evaluation.integration_params.<key>` as an axis (`null` = absent).

## Goals / Non-Goals

**Goals:**
- One weighting definition, fitted once and applied everywhere the fitted preprocessing is applied, with the metric the PLS space is fitted in following from it.
- Byte identity while absent, in both entry points.

**Non-Goals:**
- Any run, verdict or default change (proposal, Non-goals).
- Weights for SNF's affinity construction beyond what weighting its input blocks already does.
- A third kind. Two kinds give the selection a scalar-per-block arm and an alternative; more belongs to a later question.

## Decisions

**D1. The weight lives on the scaler.**
- `BlockScaler` gains `weight: float = 1.0`; `transform_blocks` computes `(values − mean) / scale · weight`. Multiplying by 1.0 is exact in IEEE arithmetic, so the absent path is byte-identical without a branch.
- `fit_block_scalers(blocks, transforms, *, weighting=None)` computes the weight from the transformed, standardized block: `1/√p_b` for `equal_variance` (p_b the block's feature count, constant features included, so the weight is a property of the declared block) and `1/σ₁(Z_b)` for `mfa`, σ₁ from `numpy.linalg.svd(Z_b, compute_uv=False)`. `FittedOmicsPreprocessor.block_weighting` records the kind; `require_layers` is unchanged.
- **Alternative: fold the weight into `scale`.** It would make `transform_population` and `_joint_original_scales` correct for free, but the record then cannot name the weight apart from the scale, and attribution's per-feature original scale would silently include it. Keeping the weight explicit costs one multiplication and makes every consumer's handling visible.

**D2. Weighting forces `scale=False` in the PLS core.**
- `fit_pls_latent_space` reads `integration_params["block_weighting"]`; when set, `scale = False` for the CV (`plsda_doubleCV(..., scale=False)`) and the final `fit_plsda_model`, exactly the latent-adjustment path's switch, and `cv_params["block_weighting"]` records the kind. With nuisance and weighting both set, nothing conflicts: the residualization is per feature and commutes with a per-feature scalar.
- **Why the core reads the key rather than taking a `scale` keyword.** The core already resolves its knobs from `integration_params` and records the effective ones in the metadata; a separate `scale` keyword would let the two disagree. `motco.analysis` builds `integration_params` from `TrajectoryAnalysisParams` and adds the key when set, so the two entry points cannot drift.
- **Consequence for `delta`.** Weighted `delta` is in the weighted metric and is not comparable to production's number; `angle`, `shape` and every p-value are comparable. The spec's "unit weights select production's rank" scenario pins the equivariance that makes the comparison in `production-latent-space-selection` attributable to the weights alone.

**D3. Every integration method runs on the weighted blocks.**
- The weighting is a preprocessing property, so `concat`, `pca`, `snf`, `pls_within_group` and `pls_residual_pca` receive the weighted blocks through the shared preprocessor. `pls_within_group` already fits with `scale=False`; `pca` is an SVD of the given matrix; the standardized `concat` becomes the weighted concatenation, which is what a weighted baseline should be.
- **Alternative: PLS only.** It would leave a weighted `pls` arm without a weighted `concat` ceiling or `pca` comparator, and the span-capture inputs of every space would be in another metric than the PLS fit's.

**D4. The setting is an integration parameter, not a top-level field.**
- It describes how the space is fitted, like `layers`, `forced_components` and `residual_components`; it is hash-omitted while absent with no change to `grid._SIGNATURE_OMIT_WHEN_NONE`; and a study configuration can sweep it as `evaluation.integration_params.block_weighting` today. A later default flip pins `"block_weighting": null` in committed configurations the way `min_components: 1` is pinned.
- **Alternative: a top-level `SimulationEvaluationParams.block_weighting`** (the `orientation_estimator` / `latent_coordinates` pattern). That pattern suits a measurement choice on a fitted space; this one changes the fit's input.

**D5. Record and CLI follow the existing opt-in pattern.**
- `TrajectoryAnalysisParams.block_weighting: str | None = None`, validated in `__post_init__`, popped from `to_dict()` while `None`; the record block `block_weighting` after `block_units` (the only other per-block declaration), `latent.integration_params.block_weighting` from the core's metadata; `--block-weighting` with `choices`. Nothing is written while absent, so `tests/test_analysis_identity.py` holds unchanged.

## Risks / Trade-offs

- **[`mfa` depends on the sample, so two cohorts weight differently]** → It is the definition of the kind; the record names the fitted weight per block so analyses can be compared. `equal_variance` is the sample-independent kind and the selection's first arm.
- **[A constant feature counts in `equal_variance` but contributes no variance]** → Documented; the weight is a property of the declared block, and a block of mostly constant features is a preprocessing problem upstream of MOTCO (roadmap item 2).
- **[Users read weighted `delta` against unweighted `delta`]** → The record carries the kind and the weights; the API page states that `delta` is in the weighted metric.
- **[The setting is used on real data before it is validated]** → The API page says so, and roadmap item 5 names `production-latent-space-selection` as the check.

## Migration Plan

Nothing migrates. The setting is absent everywhere; rollback is removing the field, which nothing reads while absent.
