# Proposal

## Why

The production latent space is stage-PLS-DA on every pooled-standardized feature, with only its rank chosen by cross-validation. At the case study's regime (n = 80, 10,000 features, held-out stage AUROC 0.636) the fitted span keeps 0.10 of the population stage program and is mostly noise (latent-space comparison 2026-11-04, §9), and the roadmap's open question for production is not which objective to supervise on — every label-free alternative lost in-span power — but how much noise the stage-supervised fit carries. Sparse PLS-DA keeps the stage objective and the rank rule and adds one regularization, the number of features each component may load on, selected by the same held-out AUROC. It is the one candidate family that generalizes production exactly: at the full budget it is the production fit. This change builds the estimator and its plumbing; `production-latent-space-selection` measures it against production on matched data.

## What Changes

- **A sparse PLS-DA estimator** (`motco.stats.sparse_pls.SparsePLSDA`) with the surface MOTCO reads off `PLSRegression` (`x_weights_`, `x_loadings_`, `y_loadings_`, `x_rotations_`, `x_scores_`, `_x_mean`, `_x_std`, `transform`, `inverse_transform`, `predict`), so the latent map, span capture, attribution, VIPs, the split test, the residual-PCA diagnostic and the orthonormal coordinates work unchanged. Per component it soft-thresholds the X weight vector to a budget of features (the sPLS-DA of Lê Cao et al.), in one of two scopes:
  - `joint`: one budget over all features;
  - `block`: one budget per block, so every block keeps its share whatever its size (the per-block sparsity of multi-block sPLS-DA, without its between-block agreement term).

  At the full budget every component is the ordinary PLS component, so the estimator equals `PLSRegression` to numerical tolerance, up to column signs.
- **Joint selection of rank and budget in the double CV.** `plsda_doubleCV(..., sparsity=)` scores every (rank, keep fraction) pair in the inner loop by held-out AUROC, picks the pair per outer fold, aggregates per repeat (modal rank as today; modal fraction, ties to the fuller budget) and refits per repeat with a sparse model. The table gains a `keep_fraction` column only when sparsity is given; without it the function is byte-identical (`tests/data/plsda_doublecv_identity.json`).
- **The production core takes the setting.** `fit_pls_latent_space` reads `integration_params["sparsity"]` = `{"scope": "joint" | "block", "keep_fractions": [...]}` (fractions of the scope's feature count, in (0, 1], distinct; the default ladder 1.0, 0.5, 0.2, 0.1, 0.05), derives the block membership from the concatenation order, runs the CV over the grid, fits the final sparse model at the selected pair and records in the metadata the ladder, the scope, the selected fraction and the per-component, per-block support sizes. A forced rank is allowed with a one-element ladder only. Latent adjustment, orthonormal coordinates and block weighting compose with it; attribution works through the shared surface.
- **Harness and real-data settings.** `integration_params["sparsity"]` (rejected with a method other than `pls`/`pls_residual_pca`, and for a malformed value, before any preprocessing; absent keys are already left out of cell ids and signatures, and `evaluation.integration_params.sparsity` is an axis for free); `TrajectoryAnalysisParams.sparsity` (`SparsitySettings(scope, keep_fractions)`, omitted from `to_dict()` while `None`); `motco analyze --sparsity {joint,block} [--sparsity-ladder 1.0,0.5,...]`. The record gains `latent.integration_params.sparsity` and a `sparsity` block after `latent` only when set; `pls_loadings.csv` shows the dropped features at zero.
- **A reference fixture.** `tests/data/splsda_reference_fixture.npz`, produced once by mixOmics (`splsda`, `scale = FALSE`, regression mode) through an R script under `tests/data/reference/`, pins the sparse weights' support and values and the scores to 1e-6 up to column signs; R is never needed at runtime (the InterSIM pattern).
- **Documentation.** `docs/api/pls.md`, `docs/api/simulations.md`, `docs/api/analysis.md`, `docs/api/cli.md`, `AGENTS.md`. Version 0.22.0.

**Non-goals.**
- **No run, verdict or default change.** Size, specificity and power against production are `production-latent-space-selection`'s. Until its reading, the setting is documented as unvalidated.
- **The between-block agreement term of multi-block sPLS-DA (DIABLO's design matrix).** It optimizes block agreement, not stage discrimination; with block scaling, multi-block PLS without it is PLS on the weighted concatenation (Westerhuis, Kourti & MacGregor 1998), which `block-weighting` already provides. Per-block budgets are what the block scope adds.
- **Sparsity on the Y side**, per-component budgets differing by component, or tuning rules other than held-out AUROC.
- **Other regularizations** (ridge PLS, elastic-net weights).

## Capabilities

### New Capabilities
- `sparse-stage-pls`: the sparse PLS-DA estimator: its budget semantics in both scopes, its `PLSRegression`-compatible surface, its equality with ordinary PLS at the full budget, its determinism, and its agreement with the reference implementation.

### Modified Capabilities
- `plsda-nested-cv`: the double CV selects the rank and the keep fraction jointly when a sparsity design is given, with a `keep_fraction` column and sparse per-repeat models; without it, byte-identical.
- `simulation-evaluation-harness`: the PLS space accepts the sparsity setting, records it only when set, composes it with the floor, forced rank (one-element ladder), latent coordinates and block weighting, and rejects it with other methods and malformed values before any preprocessing.
- `real-data-trajectory-analysis`: the analysis offers the setting, records it only when set, and `motco analyze` exposes it.
- `module-contracts`: `motco.stats.sparse_pls` exposes its public symbols, also from the flat `motco.stats` namespace.

## Impact

- **Code:**
  - `src/motco/stats/sparse_pls.py` (new): `SparsePLSDA`, `SparsityDesign`, `SPARSITY_SCOPES`, `DEFAULT_KEEP_FRACTIONS`; `src/motco/stats/__init__.py` exports;
  - `src/motco/stats/pls.py`: the `sparsity` keyword, the (rank, fraction) grid, the tie rules, the table column, sparse refits, `_plsda_auroc` over an estimator factory;
  - `src/motco/simulations/evaluation.py`: reading and validating the key, block membership from `layer_feature_counts`, the final sparse fit, metadata, the forced-rank rule, method rejections;
  - `src/motco/analysis.py`: `SparsitySettings`, the parameter, `to_dict`, the record block;
  - `src/motco/cli.py`: `--sparsity`, `--sparsity-ladder`;
  - `src/motco/__init__.py`: 0.22.0.
- **Tests:** full-budget equality with `PLSRegression` (scores, rotations, predictions, up to column signs, 1e-6) under `scale=True` and `scale=False`; the reference fixture; budget semantics (support sizes per component and per block, every block keeps at least one feature, soft threshold shrinks the survivors by the largest discarded magnitude); determinism; the CV grid (a one-element ladder of 1.0 reproduces production's table and models; the selected pair; tie rules; `keep_fraction` column present only with sparsity; `tests/data/plsda_doublecv_identity.json` unchanged); core and harness (metadata, rejections, forced-rank rule, composition with latent adjustment, orthonormal coordinates and weighting, invisibility while absent against the committed cell ids, signatures and identity fixtures); analysis record, `to_dict`, CLI round trip, attribution and the split test through a sparse model.
- **Dependencies:** none at runtime. The fixture needs mixOmics once, in the flake's R environment (`rPackages.mixOmics` beside `InterSIM`).
- **Compatibility:** every committed configuration, fixture, record and report is unchanged; the setting is absent everywhere.
- **Compute:** a CV over a five-step ladder costs about five production CVs; measured in the selection's smoke.
