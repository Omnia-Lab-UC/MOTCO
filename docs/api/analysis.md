# Real-Data Trajectory Analysis

`motco.analysis` runs the measurement the simulation studies validated on your
own data: pooled per-feature standardization, a stage-supervised PLS-DA latent
space with the double-CV rank rule, `delta` / `angle` / `shape` between two
groups in that space with the net-displacement orientation estimator, and RRPP
p-values. `analyze_trajectories` is the Python entry point and `motco analyze`
the command-line one; both write the same analysis artifact.

Use this entry point rather than assembling `motco plsr` + `motco de` by hand.
Those commands are a different pipeline: they skip the M-value contract, choose
the rank without the parsimony tie-break, measure `angle` with PC1 by default,
and record neither the selected rank, the null summary nor the spectrum.

## Input contract

- **Blocks.** One or more named omic blocks, each a samples × features table
  with string feature names. Blocks are concatenated in the order you declare
  them, and every joint feature is named `"{block}__{feature}"`. Rows must equal
  the metadata's sample ids in the same order; every block needs at least one
  feature and only finite values.
- **Units.** Each block declares a unit transform: `identity` (the default) or
  `beta_to_mvalue` (the logit of methylation beta values; values must lie in
  the open interval (0, 1)). **MOTCO does not normalize counts, filter features
  or impute values.** Supply blocks in the continuous units you want measured,
  for example log-normalized pseudobulk expression. The record states each
  block's declared units.
- **Metadata.** A sample-id column, a group column with exactly two labels and
  a stage column, all named by you. Every group × stage cell must be observed.
- **Stage order.** Required and explicit. Net displacement runs from the first
  to the last declared stage, and `shape`, the spectrum and every table follow
  this order. Stage order is never inferred by sorting labels, because
  `High` < `Intermediate` < `Low` alphabetically. Group order is optional and
  defaults to the labels sorted as strings.

- **Covariates (optional).** `covariate_cols` names metadata columns, in model
  order; see [Covariates](#covariates).

Violations raise `TrajectoryAnalysisError` (a `ValueError`) naming the block,
labels or cell at fault.

## What is measured and recorded

The latent space is fitted once on all samples, with the rank chosen by double
cross-validation (the modal rank across outer repeats, the smaller rank winning
ties). The defaults are the ones the operating studies ran: `cv2_splits=4`,
`cv1_splits=3`, `n_repeats=5`, `max_components=20`, `random_state=1203`, each
clamped to what the smallest stage can support. `forced_components` fixes the
rank and skips cross-validation for the dimensionality sensitivity analysis;
the record then says `component_selection: "forced"`.

`record` (written as `result.json`) holds:

| Key | Contents |
| --- | --- |
| `statistics`, `p_values` | `delta`, `angle`, `shape` and their RRPP p-values (`shape` is `null` with two stages) |
| `null_summary` | per statistic: retained draws, mean, sd, 50/90/95/99 % quantiles of the permutation null |
| `config_spectrum` | eigenspectrum of the latent stage-mean configuration, pooled and per group, plus the permutation eigengap summary |
| `latent` | `selected_lv`, `component_selection`, CV settings, `cv_mean_auroc`, feature counts, `feature_order_signature` |
| `design` | group and ordered stage levels, cell sizes per group × stage, contrast, analysed sample count |
| `covariates` | only when covariates are declared: names, kinds, factor levels and reference levels, design columns, the LS-mean reference row, missing policy, dropped samples, coefficients per latent component |
| `blocks` | each block's name, feature count, transform and units |
| `measurement` | orientation estimator, selected rank, permutations, seed, runtime |
| `params`, `software` | resolved parameters; MOTCO, numpy, scipy, scikit-learn, pandas and Python versions |
| `attribution` | the compact attribution record, or a `not_requested` marker |

The selected rank sits beside the `angle` result because orientation is only
interpretable relative to the space it was measured in. The null summary shows
how wide the observed statistic's own null was, and the relative eigengap in the
spectrum indicates whether a group's orientation is resolvable at all.

## Equivalence to the simulation harness

`analyze_trajectories` and `evaluate_semisynthetic_trajectory(...,
integration_method="pls")` call the same block-scaling cores and the same
`fit_pls_latent_space`. On identical data, stage order and settings they return
the same latent matrix, statistics, p-values, selected rank, null summary and
spectrum. `tests/test_analysis.py` checks this for the default estimator, PC1
and a two-block layer selection.

## Attribution

With `attribution=AttributionDiagnosticSettings(enabled=True, ...)` (CLI:
`--attribution`), frozen-model orientation attribution runs on the fitted
pooled estimator and the standardized joint matrix it was fitted to, with your
stage order and group labels. The latent space is not refitted. Real data carry
no generator truth, so the truth-recovery block is reported unavailable rather
than invented. The full tables are written under `attribution/`.

## Analysis artifact

`write_analysis_outputs(result, out_dir)` writes:

- `result.json`: the record above.
- `latent_scores.csv`: samples × `pls_*`, indexed by sample id.
- `pls_loadings.csv`: `block`, `feature`, `pls_*` (the PLS `x_rotations_`, in joint order).
- `stage_means.csv`: `group`, `stage`, `pls_*`, the fitted latent LS-mean vectors
  in declared stage order. Pass these to `motco.viz.plot_trajectories`.
- `attribution/*.csv`, when attribution ran.

It refuses to replace an existing `result.json` unless `overwrite=True`
(CLI: `--overwrite`). The fitted scikit-learn estimator is not pickled, because
pickles are version-fragile. The scores, loadings and recorded feature order
signature are enough to audit and re-project the fit.

## Covariates

Donor covariates such as age at death, PMI or RIN enter the trajectory
regression in the latent space, so the stage means, the RRPP reduced model and
therefore what is permuted all account for them. Declare them on the dataset:

```python
dataset = TrajectoryDataset(
    ...,
    covariate_cols=["age_at_death", "pmi", "education_level"],
    categorical_covariates=["education_level"],   # optional; non-numeric columns are categorical anyway
    missing_covariates="complete_case",           # default "error"
)
```

CLI: `--covariate age_at_death --covariate pmi --categorical-covariate
education_level --missing-covariates complete_case`.

- **Coding.** A numeric (non-boolean) column enters as one design column holding
  its values. Any other column, or one listed in `categorical_covariates`
  (for numeric-coded factors), enters as drop-first indicators `name[level]`,
  with levels ordered as strings (numerically when every label is an integer);
  the first level is the reference.
- **Model.** The full model is `group × stage + covariates` and the reduced
  model `group + stage + covariates`. Covariates are additive: no interaction
  with group or stage. RRPP permutes the reduced model's residuals as before, so
  every permuted dataset keeps the fitted covariate effects.
- **LS-means.** Every group × stage mean is evaluated at the pooled mean of each
  covariate design column (for a factor, its observed level proportions). The
  covariate term then adds the same vector to every stage mean, and `delta`,
  `angle` and `shape` are translation-invariant, so **the statistics do not
  depend on the reference**. Only `stage_means.csv` does, and the record states
  the reference row.
- **Missing values.** `missing_covariates="error"` (the default) rejects any
  missing covariate value, naming the column and up to ten samples.
  `"complete_case"` drops those samples from the metadata and every block before
  standardization, the latent fit and the design, and records the dropped ids.
  A deletion that empties a group × stage cell is an error. Covariates are never
  imputed; blocks must still be complete.
- **Guards.** A covariate that is also the sample, group or stage column, a
  constant numeric covariate, a factor with one observed level, a non-finite
  numeric value, and a model with covariates that is not of full column rank
  (for example a covariate collinear with group × stage) are rejected with an
  error naming the covariate(s).

**The latent space is covariate-blind, and the adjustment is still exact.**
PLS-DA is fitted on the standardized features against stage alone, exactly as
without covariates, so the latent scores, loadings and selected rank do not
change. OLS commutes with a fixed linear map: regressing the scores `T = XW` on
`[design, covariates]` gives the coefficients of regressing `X` on the same
design, projected by `W`. The stage means are therefore adjusted exactly as in
feature space; only the basis `W` was chosen without the covariates. If a
covariate effect is large and aligned with stage, PLS may spend a component on
it. A sensitivity analysis that residualizes the features on the covariates
before the fit (outside MOTCO) addresses that case.

**Attribution is not available with covariates.** Frozen-model attribution
works from unadjusted cell means, so a request for both raises an error before
any fitting.

**Keep the covariate set parsimonious.** Small cells (the SEA-AD MTG design
has 9–28 donors per cell) and many covariate columns can make the model
near-singular or the null unstable. The record's `design_columns` states the
degrees of freedom used.

**Type I validation.** `scripts/covariate_type_i_check.py` measures size at the
Phase 6 design with covariates confounded with group, stage and their
interaction (adjusted vs unadjusted arms on a shared latent fit); see
[the report](../reports/covariate-type-i-2026-10-16.md). Report covariate-adjusted
results only under its `covariates_validated` decision.

Without covariates the record, the artifact and every output are identical to
MOTCO 0.7.0 (except the runtime); the `covariates` key is absent and the schema
version is unchanged (`tests/test_analysis_identity.py`).

## Example: SEA-AD-style design

Two blocks (RNA and ATAC donor pseudobulks), sex as the group, three ADNC stages:

```bash
motco analyze \
  --block rna=rna_logcpm.csv --block atac=atac_logcpm.csv \
  --metadata donors.csv --sample-col donor_id --group-col Sex --stage-col ADNC \
  --stage-order "Not AD + Low,Intermediate,High" \
  --permutations 999 --seed 20261015 \
  --out-dir results/motco_mtg_astro
```

For the dimensionality sensitivity analysis, re-run with `--forced-components 2`,
`3`, … into separate output directories.

```python
from motco.analysis import (
    OmicsBlock, TrajectoryAnalysisParams, TrajectoryDataset,
    analyze_trajectories, write_analysis_outputs,
)

dataset = TrajectoryDataset(
    blocks=[OmicsBlock("rna", rna), OmicsBlock("atac", atac)],
    metadata=donors,
    sample_col="donor_id", group_col="Sex", stage_col="ADNC",
    stage_order=["Not AD + Low", "Intermediate", "High"],
)
result = analyze_trajectories(dataset, TrajectoryAnalysisParams(permutations=999, seed=20261015))
write_analysis_outputs(result, "results/motco_mtg_astro")
```

## API

::: motco.analysis
