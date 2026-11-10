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
ties) from the candidate ranks at or above the **rank floor**, `min_components=2`
by default (CLI `--min-components`), so `angle` and `shape` are measured in a
space of at least two dimensions. `min_components=1` restores the unfloored
0.9.0 rule and reproduces every 0.9.0 record and table byte for byte; a floor
below 1 is rejected before anything is fitted. The CV defaults are the ones the
operating studies ran: `cv2_splits=4`, `cv1_splits=3`, `n_repeats=5`,
`max_components=20`, `random_state=1203`, each clamped to what the smallest
stage can support (`max_components` to at least the floor plus one; a sample
that cannot satisfy the floor is an error naming the counts). `forced_components`
fixes the rank and skips cross-validation for the dimensionality sensitivity
analysis, ignoring the floor; the record then says `component_selection: "forced"`.

`record` (written as `result.json`) holds:

| Key | Contents |
| --- | --- |
| `statistics`, `p_values` | `delta`, `angle`, `shape` and their RRPP p-values (`shape` is `null` with two stages) |
| `null_summary` | per statistic: retained draws, mean, sd, 50/90/95/99 % quantiles of the permutation null |
| `config_spectrum` | eigenspectrum of the latent stage-mean configuration, pooled and per group, plus the permutation eigengap summary |
| `latent` | `selected_lv`, `component_selection`, CV settings (`min_components` beside them when the floor is above 1, `block_weighting` when set, `sparsity` and `selected_keep_fraction` when set), `cv_mean_auroc`, feature counts, `feature_order_signature` |
| `sparsity` | only with [sparse stage-PLS](#sparse-stage-pls): the scope, the ladder, the selected keep fraction and, per component, the support size per block |
| `design` | group and ordered stage levels, cell sizes per group × stage, contrast, analysed sample count |
| `covariates` | only when covariates are declared: names, kinds, factor levels and reference levels, design columns, the LS-mean reference row, missing policy, dropped samples, coefficients per latent component |
| `blocks` | each block's name, feature count, transform and units |
| `block_weighting` | only with a [block weighting](#block-weighting): the kind and each block's name, fitted weight and feature count |
| `measurement` | orientation estimator, selected rank, permutations, seed, the [residual correction](#residual-correction) when one was applied, runtime |
| `span_split` | only with the [split interaction test](#split-interaction-test): per part (`total`, `span`, `complement`) the dimension, interaction and residual sums of squares, degrees of freedom, pseudo-F and, with permutations, the pseudo-F and interaction-SS p-values; the in-span share, the residual correction, permutations, seed and `draws_match_rrpp` |
| `interaction_scan` | only with the [per-feature interaction scan](#per-feature-interaction-scan): feature and degenerate counts, df, the maximum F with its feature and block and, with permutations, the global p-value, the hits at adjusted p ≤ 0.05 (block, feature, F, adjusted p, by F), per block the hit count and the share of null maxima, and `draws_match_rrpp`; permutations, seed and residual correction |
| `params`, `software` | resolved parameters; MOTCO, numpy, scipy, scikit-learn, pandas and Python versions |
| `attribution` | the compact attribution record, or a `not_requested` marker |

The selected rank sits beside the `angle` result because orientation is only
interpretable relative to the space it was measured in. The null summary shows
how wide the observed statistic's own null was, and the relative eigengap in the
spectrum indicates whether a group's orientation is resolvable at all. The space
sees a group difference only through its component inside the fitted span, so a
non-rejection of `angle` or `shape` means no detectable difference inside the
stage-supervised span, not the same direction or shape in feature space; see
[Measurement visibility](../measurement-visibility.md).

**The floor guarantees two dimensions, not a second axis with signal.** Under
the unfloored rule a weak stage signal makes every candidate rank tie near an
AUROC of 0.5 and parsimony selects rank 1, where the two groups' orientations
are parallel or antiparallel by construction, `angle` is 0° or 180°, its RRPP
null is a two-point distribution, and collinear stage means carry no bend for
`shape`; both SEA-AD cohorts selected rank 1. With the floor, when the stage
signal supports one component the second is the best of the remaining
candidates, and `angle`'s null widens accordingly. Read the record's
`config_spectrum` relative eigengap and `null_summary.angle` together: they say
how resolvable orientation is in each analysis. The floor was validated on a
small-n construction calibrated to select rank 1 by `scripts/rank_floor_check.py`
(`rank_floor_decision.json`), with the floor-1 arm reported beside it: the
2026-10-23 run returned `adopt_floor` at both calibrated strengths
([report](../reports/rank-floor-2026-10-23.md)), with power at the null in both
arms, so the floor defines `angle` and `shape` without making them powerful.

## Equivalence to the simulation harness

`analyze_trajectories` and `evaluate_semisynthetic_trajectory(...,
integration_method="pls")` call the same block-scaling cores and the same
`fit_pls_latent_space`. On identical data, stage order and settings, the rank
floor included, they return
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

**With covariates**, attribution uses covariate-adjusted cell means: the
least-squares means of the cell-means model with the declared covariate design
columns (`record["covariates"]["design_columns"]`), evaluated at the record's
`covariates.reference_row`. They are the feature-space counterpart of the
latent LS-means, so the attribution's observed means, projected through the
fitted estimator, equal `stage_means.csv`. Its bootstrap resamples within
cells exactly as without covariates and refits the covariate coefficients in
every replicate; a resample whose covariates are collinear within cells (for
example, one that omits the only donor at a rare factor level) is an invalid
replicate and is counted. The record's attribution block then gains an
`adjustment` block: `mean_source` (`covariate_adjusted`), `design_columns`,
`reference_row`, `rank_deficient_replicates` and `bootstrap_replicates`. The
CLI prints `Attribution: covariate-adjusted means at the reference row; <k> of
<B> bootstrap replicates rank-deficient`. See [covariate-adjusted
means](attribution.md#covariate-adjusted-means) for the definition and its
interpretation boundary.

**Validation.** `scripts/adjusted_attribution_check.py` returned
**`adjusted_attribution_validated`** at the case study's configuration (astrocyte
`cohort_30` table, six covariate columns with one confounded and injected,
rank floor 2, `df`, at the Phase 6 stage effect; 500 × 999 × B = 200). Against
an oracle that removes the injection exactly, the adjusted attribution's top-20
contamination is within 0.005, its residual tilt toward the covariate is under
1 % of the unadjusted one, and truth precision and bootstrap stability are
non-inferior. The PLS-captured and residual views keep 3–15 % of the
confounding, so read drivers from the observed component. See [the
report](../reports/adjusted-attribution-2026-11-02.md).

**On real blocks.** The plasmode demonstration on SEA-AD astrocytes
(`scripts/plasmode_check.py`) planted an orientation difference on 20 of the
10,000 RNA and ATAC features and returned **`attribution_leg_failed`**. Size
held at the case study's configuration, but the plant was not detected at any
strength up to 3 within-stage SDs: it lies almost entirely outside the
stage-PLS latent space the statistics are measured in. Attribution's observed
contrast recovered part of the planted set, unstably. A difference confined to
a few features is therefore not something this analysis can detect or trace at
n = 80. See [the report](../reports/plasmode-2026-11-02.md).

**Not with latent adjustment.** Attribution together with
`covariate_latent_adjustment` raises an error before any fitting: that basis
is fitted on covariate-residualized features, and attribution on it is not
defined.

**The frozen measurement in memory.** `result.latent_artifacts` holds the
fitted estimator, the standardized joint feature matrix the latent scores are
the projection of (`latent_artifacts.model.transform(latent_artifacts.features)`
reproduces `latent_scores` to rounding) and the per-feature original scales. It
exists for callers that attribute or probe in exactly the analysed coordinates.
It is in memory only: it is not part of the record, and
`write_analysis_outputs` never writes or pickles it.

## Analysis artifact

`write_analysis_outputs(result, out_dir)` writes:

- `result.json`: the record above.
- `latent_scores.csv`: samples × `pls_*`, indexed by sample id.
- `pls_loadings.csv`: `block`, `feature`, `pls_*` (the PLS `x_rotations_`, in joint order).
- `stage_means.csv`: `group`, `stage`, `pls_*`, the fitted latent LS-mean vectors
  in declared stage order. Pass these to `motco.viz.plot_trajectories`.
- `attribution/*.csv`, when attribution ran.
- `interaction_scan.csv`, when the [per-feature interaction scan](#per-feature-interaction-scan)
  ran: one row per joint feature.

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

**By default the latent space is covariate-blind, and the adjustment is still exact.**
PLS-DA is fitted on the standardized features against stage alone, exactly as
without covariates, so the latent scores, loadings and selected rank do not
change. OLS commutes with a fixed linear map: regressing the scores `T = XW` on
`[design, covariates]` gives the coefficients of regressing `X` on the same
design, projected by `W`. The stage means are therefore adjusted exactly as in
feature space; only the basis `W` was chosen without the covariates. If a
covariate effect is large and aligned with stage, PLS may spend a component on
it; [latent adjustment](#latent-adjustment) chooses the basis with the
covariates instead.

**Attribution with covariates** uses covariate-adjusted cell means at the
same reference row; see [Attribution](#attribution).

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

## Latent adjustment

With covariates declared, the PLS-DA basis is still fitted on the features
alone. A large technical axis aligned with stage (library quality, for example)
can then take a whole component, and at a selected rank of 1 or 2 that leaves
little room for the biology. `covariate_latent_adjustment` (CLI:
`--latent-adjustment`) fits the basis on features residualized on the declared
covariates instead:

```python
params = TrajectoryAnalysisParams(permutations=999, seed=1, covariate_latent_adjustment="stage")
```

```bash
motco analyze ... --covariate age_at_death --covariate pmi --latent-adjustment stage
```

- **Conditional residualization.** Every pooled-standardized joint feature is
  regressed on `[1, D, C]`, where `C` holds the declared covariate design columns
  and `D` the conditioning design, and only the fitted covariate part `C γ` is
  subtracted. Under `stage`, `D` is the drop-first stage indicators; under
  `group_stage`, the group, stage and group × stage indicators. By
  Frisch–Waugh–Lovell, `γ` is the covariates' within-design effect, so stage
  signal that is merely correlated with a covariate stays in the features.
  Residualizing on the covariates alone would remove it: label-blind
  residualization upstream of MOTCO dropped the case study's held-out CV AUROC
  from 0.63 to 0.41. `stage` keeps the group label out of the basis; under
  `group_stage` the group enters only through the pooled within-cell slope.
  Use `stage`: it is the conditioning the Type I check admitted (see below).
- **Leak-free cross-validation.** Inside the double CV, every outer rest set and
  inner training set estimates its own `γ` and residualizes its held-out samples
  as `X_held_out − C_held_out γ`, so held-out labels never enter the
  residualization. Both CV levels stratify their folds on the conditioning cells
  (stages under `stage`, group × stage cells under `group_stage`), so every
  training set can estimate the conditioning design. Under `stage` the folds
  equal those of the unadjusted analysis. The final basis uses `γ` estimated on
  all samples.
- **The pooled scale is kept.** The residualized features are not
  re-standardized: every PLS fit on this path uses `scale=False`. A feature whose
  variance was mostly covariate is down-weighted, not inflated back to unit
  variance, and when the covariate effects are weak the basis approaches the
  unadjusted one.
- **The raw features are measured in the new basis.** The latent scores are the
  pooled-standardized features (not the residuals), centered at the residuals'
  mean and projected onto the basis. The trajectory model, the LS-means, the
  RRPP reduced model and the permutations are exactly those of the
  covariate-adjusted analysis above. With the covariates in both RRPP models,
  projecting the raw or the residualized features gives identical statistics,
  p-values and null draws, because the difference `(C − C̄) γ W` lies in the
  models' column space. So **only the choice of basis (and of CV rank)
  changes**. Measuring the raw features keeps the `covariates` block's
  coefficients and `stage_means.csv` on the same footing as without the
  setting.
- **The same covariates in both places.** The residualization covariates are
  exactly `covariate_cols`, which also enter the trajectory model; that equality
  is what makes the projection choice immaterial. To residualize on technical
  variables, declare them as covariates.
- **Guards.** The setting without covariates is an error before any fitting.
  The full-sample residualization design must have full column rank, otherwise
  the error names the covariates and the conditioning. Every CV training design
  must too, otherwise the error names the repeat and fold; this happens when a
  rare factor level drops out of a training set, and depends on
  `random_state`. Merge rare levels upstream rather than relying on a
  minimum-norm fit, whose held-out shift would be arbitrary. Attribution is
  rejected with latent adjustment, before any fitting.
- **Record.** `result.json` gains a `latent_adjustment` block after
  `covariates`, and `latent.integration_params.latent_adjustment` names the
  conditioning. The block holds `conditioning`, `covariates`, `design_columns`
  (conditioning then covariate columns, beyond the intercept), `n_samples`,
  `residual_df` (`n − 1 − k − q`), `cv_min_training_df` (the smallest
  training-set residual df seen in CV; `null` with `forced_components`),
  `feature_scale` (`pooled_standardized`) and `feature_r2`. That last field
  summarizes each feature's partial R² of the covariates given the intercept and
  conditioning (5/25/50/75/95 % quantiles, mean, fraction above 0.5), the share
  of variance the residualization removes. `params.covariate_latent_adjustment`
  is recorded only when set.

**Mediator limitation.** Residualizing a variable that mediates the stage
effect removes real signal, and no simulation can tell a confounder from a
mediator. Choosing and justifying the covariate set is the analysis's job, not
MOTCO's: a technical variable that is itself downstream of disease stage (a
cell-quality metric that worsens with pathology, say) is a mediator candidate.

**Validation.** `scripts/latent_adjustment_type_i_check.py` measures size at
the Phase 6 design under a covariate-aligned nuisance axis calibrated so that
the unadjusted basis spends its first component on it (median |Spearman| ≥ 0.8).
It runs three arms per replicate (unadjusted basis, `stage`, `group_stage`) at
8 declared covariates, and its decision file names the admissible
conditionings (`adopt_stage` | `adopt_group_stage` | `diagnose_latent_excess`).
The 2026-10-20 run returned **`adopt_stage`**
([report](../reports/latent-adjustment-type-i-2026-10-20.md)). `stage` held
`delta` 0.031, `angle` 0.067 and `shape` 0.058 against 0.0706. `group_stage`
exceeded the bound on `angle` (0.082) and is not admissible. About 10 % of
replicates lost an arm to a rare factor level.
`scripts/latent_adjustment_power_comparison.py` compares power with and without
such an axis against the unadjusted basis. Its 2026-10-20 run
([report](../reports/latent-adjustment-power-2026-10-20.md)) returned
`candidate_fails` for both conditionings. `stage` raised orientation power
(`angle` 0.92–0.93 against 0.64–0.78 at the top effect) but lost `shape` power
when the nuisance axis was present, and small-magnitude power without it.
`group_stage` lost most orientation and shape power. Residualizing on
covariates confounded with the design removes signal even when they have no
effect. Report latent-adjusted results only with `stage` conditioning, as a
sensitivity analysis beside the covariate-adjusted analysis rather than as a
replacement, under the conditions of both reports.

**Validation at the case study's configuration.** `scripts/case_study_type_i_check.py`
measures size where the SEA-AD case study is actually run: the rank floor of 2,
six numeric covariate design columns (one confounded with group × stage and
injected into the features, five pure noise), a stage signal weakened to a
held-out CV AUROC of 0.69–0.70 and a nuisance the covariate-blind basis follows
(median |Spearman| 0.92), with three arms per unit on the same data and
permutation draws: covariates in the trajectory model only (`adjusted`, the
primary analysis), none (`unadjusted`, descriptive) and the `stage` latent
adjustment (`latent_stage`, the fallback). Its decision names
`primary_validated` | `adopt_latent_stage` | `diagnose_configuration_excess`
over the anchor and the translation control. The 2026-10-29 run (500 × 999 at
δ = 0.1 and 0.05, conjunction) returned **`diagnose_configuration_excess`**
([report](../reports/case-study-type-i-2026-10-29.md)): the adjusted model held
every statistic at δ = 0.05 and exceeded the 0.0695 bound on translation
`delta` at δ = 0.1 (0.072; `delta` 0.060–0.072 in every reading, 0.066 pooled),
and the `stage` adjustment exceeded at the anchor on `angle` (0.084–0.090) and
`shape` (0.078–0.092) at both strengths. At that configuration, read
covariate-adjusted `delta` p-values against a realized size of about 0.066, and
do not interpret `stage`-adjusted `angle` or `shape` rejections, until the
diagnosis the report names is done.

**Diagnosis of that excess.** `scripts/case_study_diagnosis.py` reran the same
units, with the committed `adjusted` and `latent_stage` p-values reproduced
exactly (0 mismatches over 6,000 per strength), under arms that separate the
candidate causes:

- on the covariate-blind basis, the confounded covariate estimated or its
  injection removed exactly, crossed with 0, 2 or 5 noise covariates;
- each of those rerun with the Freedman–Lane residuals rescaled by
  √(n / (n − p_red));
- the `stage` adjustment at a forced rank of 2.

The 2026-10-30 run
([report](../reports/case-study-diagnosis-2026-10-30.md)) returned two verdicts:

- **`residual_df_shrinkage`** for the adjusted `delta`. RRPP's permuted residuals
  shrink with the reduced model's column count, so pooled `delta` rises from
  0.0545 to 0.066 between one and six covariate columns at n = 80. The corrected
  permutation holds every reading.
- **`residualized_basis`** for the `stage` adjustment. It still exceeds on anchor
  `angle` at a fixed rank of 2 (0.084 / 0.082).

With many covariates at small n, use the opt-in
[residual correction](#residual-correction) `df`, which the corrected-permutation
check validated at this configuration (roadmap, Next changes item 8). Without it,
read covariate-adjusted p-values against their realized size.

**Validation at the ADNI design.** `scripts/adni_design_check.py` measures size
and power where the second case study (ADNI serum metabolomics) is run:

- **The design.** n = 1,501 in ADNI's sex × stage cells, two blocks, and six
  covariate design columns of ADNI's kinds. They are a confounded numeric, a
  three-level factor declared with `categorical_covariates` (two drop-first
  columns), a phase-like binary nuisance injected into one block, and two noise
  columns.
- **The regime.** The stage signal is calibrated to ADNI's held-out AUROC
  (0.606 against 0.600), with a basis that follows the nuisance as ADNI's does.
- **The analysis.** One covariate-adjusted analysis per unit, as the primary is
  issued: the rank floor of 2, net displacement and no residual correction.

The 2026-11-02 run (14 cells × 500 × 999;
[report](../reports/adni-design-2026-11-02.md)) returned
**`non_rejections_uninformative`**:

- **Size holds** at the anchor and the translation control: every statistic
  reads 0.034–0.068 against the 0.0695 bound. `delta` is conservative
  (0.034), since six columns cost almost nothing in residual df at this n.
- **No statistic reaches the 0.80 power floor.** Magnitude/`delta` reaches
  0.530 at a doubled stage effect in one group. Orientation/`angle`
  (0.056–0.076) and shape/`shape` (0.034–0.062) stay at the null rate at
  every construction effect.
- **The rank.** The CV-selected rank sits at the cap of 19 in 95–100 % of
  units, against ADNI's 8.

At a design of this kind (large n, a weak stage signal at AUROC ≈ 0.60), a
non-rejection of `delta`, `angle` or `shape` is uninformative. The
measurement covers this table, stage signal and covariate set only, on
InterSIM-analogue blocks.

Without the setting, every output is identical to MOTCO 0.8.0 except the
runtime, `record["params"]` included (`tests/test_analysis_identity.py`).

## Residual correction

RRPP permutes the reduced model's residuals (Freedman–Lane). Their mean square
is RSS/n, not the unbiased RSS/(n − r), where r is the reduced model's rank, so
the permuted noise is narrower than the real noise by a factor of
√((n − r)/n). With many covariate columns at small n that factor matters: at
n = 80 and r = 10 (intercept, group, two stage columns and six covariates) the
permuted noise is about 6 % narrower in SD, the null of a scale-dependent
statistic is too narrow, and `delta` rejects at about 0.066 at α = 0.05 (the
[diagnosis](../reports/case-study-diagnosis-2026-10-30.md)). The opt-in
`residual_correction` rescales the residuals once before they are permuted:

```python
params = TrajectoryAnalysisParams(permutations=999, seed=1, residual_correction="df")
```

```bash
motco analyze ... --covariate age_at_death --covariate pmi --residual-correction df
```

- **`df`** multiplies every residual by κ = √(n / (n − r)): one scalar per
  analysis (κ = 1.069 at n = 80, r = 10). It is the correction the diagnosis
  measured.
- **`leverage`** divides residual i by √(1 − h_ii), with h the reduced model's
  hat matrix, so samples in small cells or with extreme covariate values are
  inflated more. In a balanced one-factor design every sample has leverage r/n
  and the two corrections coincide.

r is the numerical rank of the reduced model with its covariates, so a
rank-deficient model uses its rank, not its column count. The permutation draws
(the same seed gives the same sample orders), the reduced fitted values, the
full-model estimation and the observed statistics are unchanged: **only the
p-values and the null summary move**. The latent scores, selected rank, stage
means and every other block are identical with and without the setting. The
correction does not need covariates: without them r = 4 (intercept, group, two
stage columns at three stages), so κ = 1.026 at n = 80.

- **Errors.** An unknown value is rejected before any fitting. `df` with
  n − r ≤ 0, and `leverage` with a sample of leverage 1 (the only member of a
  reduced-model category), raise an error naming the counts or the sample's row
  before any permutation.
- **Record.** When the correction is set and permutations run,
  `measurement.residual_correction` states `kind`, `n` and `reduced_rank`, plus
  `kappa` for `df`, or `leverage_min`/`leverage_max` and `scale_min`/`scale_max`
  for `leverage`; the CLI prints it as `RRPP null: residual correction …`.
  `params.residual_correction` is recorded only when set. With
  `permutations=0` the parameter is recorded but no correction entry is
  written, since nothing was applied.
- **The primitive.** `motco.stats.permutation.RRPP(..., residual_correction=)`
  applies the same rescaling for any caller, and
  `residual_correction_scale(model_reduced, correction)` returns the per-sample
  multipliers.

**Validation status: `df` validated (`adopt_df`, 2026-10-30).** The
corrected-permutation check (`scripts/corrected_permutation_check.py`;
[report](../reports/corrected-permutation-2026-10-30.md)) measured both
corrections. It used a fresh seed family at the case study's configuration:
six covariate columns, the rank floor of 2 and a weak stage signal, on both
case-study cohort tables (n = 80 and 76) at 500 replicates × 999 permutations.

- **Size.** Both corrections hold every pooled rate and every reading. Under
  `df`, `delta` reads 0.029–0.031, `angle` 0.0425–0.050 and `shape`
  0.052–0.058 pooled per table, against a bound of 0.0597.
- **The uncorrected null.** On the same data it was liberal: pooled `shape`
  0.0675–0.0685, and anchor `angle` 0.086 at the Phase 6 stage effect.
- **Power cost.** At the Phase 6 stage effect, `df` gives up 0.024–0.034 of
  `angle` and `shape` power, and up to 0.054 of `delta` power, against the
  uncorrected null.
- **`leverage`.** It holds size as well but showed no paired power advantage,
  so `df` is the validated choice.

The validation covers n = 76–80 with six covariate columns. Smaller tables
(for example n = 63) and other covariate sets are inferred, not measured. The
correction stays opt-in, and the default is the uncorrected null.

Without the setting, every output is identical to MOTCO 0.12.0 except the
runtime (`tests/test_analysis_identity.py`).

## Split interaction test

`delta`, `angle` and `shape` see a group difference only through its component
inside the fitted stage-PLS span. The opt-in split interaction test asks the
question the span cannot: do the groups differ by stage anywhere in the measured
features, and is that difference inside the fitted span or outside it?

```python
params = TrajectoryAnalysisParams(permutations=999, seed=1, residual_correction="df", span_split_test=True)
```

```bash
motco analyze ... --permutations 999 --seed 1 --residual-correction df --span-split-test
```

It runs `motco.stats.interaction.interaction_split_test` on the analysis's own
inputs:

- **Outcome.** The pooled-standardized joint features in the fitted PLS model's
  input coordinates, the coordinates in which `qr(x_rotations_)` is orthonormal
  and span capture is measured.
- **Basis.** `qr(x_rotations_)`: `span` has the selected rank as its dimension
  and `complement` the rest of the features.
- **Models.** The full (group × stage) and reduced (group + stage) trajectory
  models with the declared covariates, so the tested effect is the group ×
  stage interaction RRPP tests.
- **Null.** Freedman–Lane permutation of the reduced residuals with RRPP's
  serial draws at the analysis's seed, the analysis's permutation count and its
  residual correction.

Per part the record's `span_split` block gives the dimension, the interaction
sum of squares, the full model's residual sum of squares, df_int, df_res, the
pseudo-F `(SS_int / df_int) / (RSS / df_res)` and, with permutations, the
pseudo-F p-value (the test) and the interaction-SS p-value (a diagnostic).
`span` + `complement` = `total` for both sums, and `span_share` is the observed
share of the interaction energy inside the fitted span. `draws_match_rrpp` is
true when RRPP ran serially (`n_jobs` 1), so the split test's draws are RRPP's;
with more jobs RRPP draws a different sequence while the split test keeps the
serial one. The CLI prints one line with the three pseudo-F p-values and the
in-span share.

- **Unchanged outputs.** Latent scores, rank, statistics, p-values, null
  summary and stage means are identical with and without the setting; without
  it every output equals MOTCO 0.17.0's (`tests/test_analysis_identity.py`,
  `tests/test_analysis_span_split.py`).
- **Latent adjustment.** The test is defined on the covariate-blind basis, so
  requesting it with `covariate_latent_adjustment` (`--latent-adjustment`) is an
  error before fitting.
- **Reading.** "Span" is the span fitted to the sample, not the population stage
  program. How a rejection or a non-rejection of each part is read, and the
  validation status, are on the [measurement visibility](../measurement-visibility.md)
  page.

## Per-feature interaction scan

The split test's complement sums interaction energy over every coordinate
outside the span, so a difference on a few features is lost in the noise of the
rest. The opt-in scan tests the group × stage interaction of each feature on
its own and combines the features by their maximum, so a difference confined to
a few features can be detected and the features carrying it named.

```python
params = TrajectoryAnalysisParams(permutations=999, seed=1, residual_correction="df", interaction_scan=True)
```

```bash
motco analyze ... --permutations 999 --seed 1 --residual-correction df --interaction-scan
```

It runs `motco.stats.interaction.interaction_feature_scan` on:

- **Outcome.** The pooled-standardized joint features. Each feature's F
  ignores its location and scale, so the standardization does not matter; the
  declared unit transform (for example beta to M values) does.
- **Models.** The full and reduced trajectory models with the declared
  covariates, the effect RRPP and the split test test.
- **Null.** Freedman–Lane permutation of the reduced residuals with RRPP's
  serial draws at the analysis's seed, its permutation count and its residual
  correction (`df` cancels in every F).

The record's `interaction_scan` block gives the feature and degenerate counts,
df_int and df_res, `max_f` with `max_feature` and `max_block` and, with
permutations, `p_value_global` (the permutation p-value of the largest F), the
`hits` (features with Westfall–Young step-down max-T adjusted p ≤ `hit_alpha` =
0.05, ordered by F) and per block the features, degenerate features, hits and
`null_max_share` (the share of draws whose largest null F fell in that block,
beside the block's share of features). `draws_match_rrpp` is true when RRPP ran
serially. The CLI prints one line with the maximum F and its feature, the
global p-value and the hit count.

`interaction_scan.csv` (`TrajectoryAnalysisResult.interaction_scan`) has one
row per joint feature in joint order: `block`, `feature`, `ss_interaction`,
`rss`, `f`, `p_value_adjusted` (with permutations), `degenerate` and
`difference[<stage>]` per declared stage, the second group minus the first from
the full model's LS means at the covariate reference row, in standardized units
(the rows that give `stage_means.csv` in latent coordinates). Degenerate
features (zero full-model residual variance, constant ones included) have no F
or adjusted p-value and never enter the maximum.

- **Unchanged outputs.** Latent scores, rank, statistics, p-values, null
  summary, stage means, any `span_split` block and any attribution output are
  identical with and without the setting; without it every output equals MOTCO
  0.19.0's (`tests/test_analysis_identity.py`,
  `tests/test_analysis_interaction_scan.py`).
- **Latent adjustment, attribution and block weighting.** The scan reads only
  the standardized features and the trajectory models, so it is allowed with
  all three and is identical with and without latent adjustment. It reads the
  unweighted standardized features, so a block weighting leaves its record and
  table unchanged.
- **Reading.** A hit names a feature whose covariate-adjusted group difference
  changes across stages, with the familywise error rate over all scanned
  features controlled at 0.05; it does not say the feature drives a latent
  trajectory difference, and a scan without hits is uninformative unless power
  is shown. See [measurement visibility](../measurement-visibility.md) for the
  reading and the validation status.

## Block weighting

The joint features are z-scored one by one and concatenated, so unweighted a
block's influence on every fitted space grows with its feature count. The
opt-in block weighting multiplies every standardized block by one label-blind
scalar, so that weight is a declared choice:

```python
params = TrajectoryAnalysisParams(permutations=999, seed=1, block_weighting="equal_variance")
```

```bash
motco analyze ... --permutations 999 --seed 1 --block-weighting equal_variance
```

- **Kinds.** `equal_variance` multiplies block *b* by 1/√p_b (p_b its feature
  count, constant features included), so every block has unit total variance
  (the multi-block PLS block scaling). `mfa` multiplies it by the inverse
  standard deviation of the standardized block's leading principal component,
  so every block's dominant direction has unit variance (multiple factor
  analysis). The `equal_variance` weight depends only on the declared block;
  the `mfa` weight depends on the sample, so two cohorts weight differently.
  Neither reads a group or stage label, and both are fitted once on the pooled
  sample, like the standardization, so RRPP's null is unchanged in kind.
- **Metric.** The PLS space is fitted without per-feature rescaling (every
  fit, the cross-validation's included), so the weights set the metric. `delta`
  is in the weighted metric and is not comparable with an unweighted
  analysis's; `angle`, `shape` and the p-values are. The rank can differ,
  because the weights change the space.
- **Record.** A `block_weighting` block after `blocks` names the kind and
  each block's weight and feature count, `latent.integration_params` names the
  kind, and the CLI prints it beside the selected rank. Attribution's
  original-unit mapping and `latent_artifacts.original_scales` divide by the
  weights, so reconstructed features are in their original units. The
  [split interaction test](#split-interaction-test) runs in the weighted metric
  of the fitted space; the [per-feature scan](#per-feature-interaction-scan)
  does not depend on it.
- **Absent.** Nothing is written and every output equals MOTCO 0.20.0's
  (`tests/test_analysis_identity.py`, `tests/test_analysis_block_weighting.py`).
- **Status.** The setting is **unvalidated**: its size and power check is the
  `weighted` arm of the production latent-space selection, and the case study
  may not use it before that reading. An unknown kind is rejected before
  anything is fitted.

## Sparse stage-PLS

The production space loads every component on every feature. At a weak stage
signal on many features the fitted span is mostly noise; the opt-in sparse
stage-PLS-DA keeps the stage objective and the rank rule and adds one
regularization, the number of features each component may load on, selected
by the same held-out AUROC:

```python
from motco.analysis import SparsitySettings

params = TrajectoryAnalysisParams(permutations=999, seed=1, sparsity=SparsitySettings("block"))
```

```bash
motco analyze ... --permutations 999 --seed 1 --sparsity block --sparsity-ladder 1.0,0.5,0.2,0.1,0.05
```

- **Budget.** Every component's weight vector keeps the features with the
  largest absolute weights, shrunk by the largest discarded one (the sPLS-DA
  of Lê Cao et al.; [`SparsePLSDA`](pls.md#sparse-pls-da)). `joint` counts one
  budget over all features; `block` one per declared block, so every block
  keeps its share whatever its size. A keep fraction *f* keeps
  `max(1, round(f · p))` features per group.
- **Selection.** The cross-validation scores every pair of a candidate rank
  and a fraction of the ladder (default `1.0, 0.5, 0.2, 0.1, 0.05`; distinct,
  in (0, 1]) and keeps, on an exact tie, the smaller rank and then the fuller
  budget, so it leaves the ordinary fit only on strictly better held-out
  AUROC. The final model is fitted at the modal rank and the modal fraction.
  A ladder of `1.0` alone reproduces the ordinary analysis exactly.
  `forced_components` skips the cross-validation, so it needs a one-element
  ladder.
- **Record.** A `sparsity` block after `latent` names the scope, the ladder,
  the selected fraction and each component's support per block;
  `latent.integration_params` carries the scope, the ladder and the selected
  fraction, and the CLI prints the fraction beside the selected rank.
  `pls_loadings.csv` shows the features no component loads on at zero.
- **Composition.** Attribution, the [split interaction test](#split-interaction-test),
  the [per-feature scan](#per-feature-interaction-scan) (which never reads the
  fit), [latent adjustment](#latent-adjustment) and [block weighting](#block-weighting)
  run unchanged on a sparse space.
- **Absent.** Nothing is written and every output equals MOTCO 0.21.0's
  (`tests/test_analysis_identity.py`, `tests/test_analysis_sparsity.py`).
- **Status.** The setting is **unvalidated until the production latent-space
  selection reads it**: its size, specificity and power against production are
  that selection's `sparse_block` arms, and the case study may not use it
  before that reading. A malformed setting is rejected before anything is
  fitted.

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
