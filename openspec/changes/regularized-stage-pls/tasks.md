# Tasks

## 1. The estimator

- [ ] 1.1 Create `src/motco/stats/sparse_pls.py` with `SPARSITY_SCOPES = ("joint", "block")`, `DEFAULT_KEEP_FRACTIONS = (1.0, 0.5, 0.2, 0.1, 0.05)`, `SparsityDesign(scope, keep_fractions, feature_groups)` (validation naming the offending scope, fraction or membership; `counts(fraction)` → per-group `max(1, round(f·p_g))`) and `SparsePLSDA(n_components, keep, *, scale=True, max_iter=1000, tol=1e-6)` implementing design D1 (centering/ddof-1 scaling like sklearn, SVD initialization, soft-thresholded weight iteration, regression-mode deflation, `x_rotations_ = W·(PᵀW)⁻¹`, `transform`, `inverse_transform`, `predict`, `keep_`, `support_`, sklearn's sign convention). Export the four symbols from `motco.stats`. Verify in `tests/test_sparse_pls.py`: joint and block support sizes on 500 and 367 + 131 features at fraction 0.1 (50; 37 and 13); the floor of one feature per block; the soft threshold (survivors shrunk by the largest discarded magnitude); determinism (two fits bitwise equal); every malformed input rejected before fitting; and the module-contract imports.
- [ ] 1.2 Verify full-budget equality: with fraction 1.0 at ranks 2–5 on a seeded three-class dataset, scores, rotations and predictions agree with `PLSRegression(..., max_iter=1000)` to 1e-6 up to column signs, under `scale=True` and `scale=False`, and `calculate_vips` returns the same values up to tolerance.
- [ ] 1.3 Write `tests/data/reference/splsda_fixture.R` (mixOmics `splsda`, `ncomp = 3`, constant `keepX`, `scale = FALSE`, regression mode, seeded X with three classes) and `tests/splsda_fixture_capture.py` packing `tests/data/splsda_reference_fixture.npz` with provenance (R and mixOmics versions, script sha256); add `rPackages.mixOmics` to `flake.nix` (design D6, with its fallback). Verify in `tests/test_sparse_pls.py` that the estimator reproduces the fixture's supports exactly and its weights and scores to 1e-6 up to column signs, with no R at test time.

## 2. Cross-validation

- [ ] 2.1 Add `sparsity: SparsityDesign | None = None` to `plsda_doubleCV` (design D3): candidate pairs ordered rank ascending then fraction descending, an estimator factory in `_plsda_auroc`, the outer evaluation at the pair, per-repeat modal rank (unchanged) and modal fraction with the fuller-budget tie-break, sparse per-repeat refits, and the `keep_fraction` column only when given. Verify in `tests/test_plsda_sparsity.py`: the selected pair on a constructed inner-AUROC pattern; both tie rules; a `[1.0]` ladder reproducing the default table's `LV` and `AUROC` and the refit scores to 1e-6 up to signs; nuisance and strata composing with the grid; and that `tests/test_plsda_identity.py` (or the test reading `tests/data/plsda_doublecv_identity.json`) passes unchanged.
- [ ] 2.2 Document the keyword in the `plsda_doubleCV` docstring and `docs/api/pls.md` (the grid, the tie rules, the column, the byte-identical default). Verify with `uv run mkdocs build`.

## 3. The production core and the harness

- [ ] 3.1 In `fit_pls_latent_space`, read and validate `integration_params["sparsity"]`, build the `SparsityDesign` with groups from `layer_feature_counts` (validated against the column count), run the CV over the grid, fit the final `SparsePLSDA` at the selected pair with the path's `scale`, and record `cv_params["sparsity"]`, `selected_keep_fraction` and the per-component, per-block support sizes (design D4). Reject a forced rank with a ladder longer than one, naming both; fit at the ladder's fraction otherwise. Verify in `tests/test_sparse_latent_space.py`: metadata contents; the forced-rank rule both ways; composition with `nuisance` (residualized fit, raw projection), `coordinates="orthonormal"` (basis spans the sparse rotations) and `block_weighting` (`scale=False`); and that the core identity fixtures pass unchanged.
- [ ] 3.2 In the harness, validate the key before any preprocessing (`_validate_evaluation_params`): reject it with a method other than `pls`/`pls_residual_pca` naming the method, and reject a malformed value. Verify one test per rejection; that `pls_residual_pca_from_pls`, `pls_orthonormal_from_pls`, span capture and attribution run on a sparse `pls` integration; and that an evaluation without the key is byte-identical to the pre-change result on a pinned dataset, with `tests/test_committed_config_identity.py` and the pin tests unchanged.
- [ ] 3.3 Document the setting in `docs/api/simulations.md` (value shape, scopes, the default ladder, rejections, composition, the forced-rank rule, invisibility while absent, that `evaluation.integration_params.sparsity` is an axis) and in the `evaluation.py` and `stats` bullets of `AGENTS.md`. Verify with `uv run mkdocs build`.

## 4. Real-data analysis and CLI

- [ ] 4.1 Add `SparsitySettings` and `TrajectoryAnalysisParams.sparsity` (validated in `__post_init__`, popped from `to_dict()` while `None`), pass the key into the analysis's `integration_params`, and write `latent.integration_params.sparsity` and the `sparsity` record block after `latent` (design D5). Verify in `tests/test_analysis_sparsity.py`: record contents and their absence without the setting; `pls_loadings.csv` shows dropped features at zero; attribution, the split test and the interaction scan run on a sparse analysis; `tests/test_analysis_identity.py` passes unchanged.
- [ ] 4.2 Add `motco analyze --sparsity {joint,block}` and `--sparsity-ladder` (comma-separated, rejected without `--sparsity`), printing the selected fraction beside the rank. Verify a CLI round trip writes the scope and ladder, and the rejection fires.
- [ ] 4.3 Document the parameter in `docs/api/analysis.md` ("unvalidated until the production latent-space selection reads it") and `docs/api/cli.md`; update the `analysis.py` and `cli.py` bullets of `AGENTS.md`. Verify with `uv run mkdocs build`.

## 5. Version and gate

- [ ] 5.1 Bump `src/motco/__init__.py` to 0.22.0 and run the pre-commit gate: `uv run ruff check src/ tests/ && uv run mypy src/motco/ && MOTCO_TEST_PERMS=99 uv run pytest tests/ -m "not slow" --tb=short`.

## Workflow follow-up

- Sync the delta specs and archive the change once merged.
- Size, specificity and power against production are measured by `production-latent-space-selection`; a default flip follows only from its verdict, as its own change.
