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
| `design` | group and ordered stage levels, cell sizes per group × stage, contrast |
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

## Limitation: no covariates

The design is group × stage with interactions. Covariates such as age, PMI or
batch are not adjusted for: the operating characteristics were measured without
them, and adding them is a methodological change that needs its own validation.
Handle known batch structure upstream, and report this limitation with any
result.

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
