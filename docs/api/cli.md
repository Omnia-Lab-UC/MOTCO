# CLI Reference

The `motco` command-line interface provides five subcommands.

```
usage: motco [-h] [--version] [--verbose] {plsr,snf,simulate,de,analyze} ...

MOTCO CLI: PLSR, SNF, group differences, and real-data trajectory analysis

positional arguments:
  {plsr,snf,simulate,de,analyze}
    plsr                Run PLS-DA with double cross-validation
    snf                 Similarity Network Fusion
    simulate            Generate a semi-synthetic multi-omics toy dataset
                        (numpy generator, no R)
    de                  Group differences on trajectories (estimate or RRPP)
    analyze             Real-data trajectory analysis with the validated
                        pipeline (PLS latent space, RRPP)

options:
  -h, --help            show this help message and exit
  --version             show program's version number and exit
  --verbose             Enable debug logging to stderr
```

For real data, use `motco analyze`. It runs the measurement the simulation
studies validated, in one step. `plsr` and `de` are lower-level building blocks.

---

## `motco plsr`

Run PLS-DA with double cross-validation. Input is one of: omics layer CSVs (`--input`, repeated; standardized and concatenated) with a metadata CSV holding the label (`--metadata`, `--label-col`); a single CSV with predictors and the label column (`--data`, `--label-col`); or separate predictor/label CSVs (`--x` / `--y`). `--out-scores` writes the final model's latent scores, the input `motco de` takes as `--Y`.

```
usage: motco plsr [-h] [--input INPUT] [--metadata METADATA] [--data DATA]
                  [--label-col LABEL_COL] [--x X] [--y Y]
                  [--cv1-splits CV1_SPLITS] [--cv2-splits CV2_SPLITS]
                  [--n-repeats N_REPEATS] [--max-components MAX_COMPONENTS]
                  [--min-components MIN_COMPONENTS]
                  [--random-state RANDOM_STATE] [--out-table OUT_TABLE]
                  [--out-vips OUT_VIPS] [--out-scores OUT_SCORES]
                  [--n-components N_COMPONENTS]

options:
  -h, --help            show this help message and exit
  --input INPUT         Omics layer CSV (repeat for multiple layers;
                        standardized and concatenated)
  --metadata METADATA   Metadata CSV with label column (required when using
                        --input)
  --data DATA           CSV with predictors and label column
  --label-col LABEL_COL
                        Label column name (used with --data or
                        --input/--metadata)
  --x X                 CSV with predictors (features)
  --y Y                 CSV with labels/outcomes
  --cv1-splits CV1_SPLITS
  --cv2-splits CV2_SPLITS
  --n-repeats N_REPEATS
  --max-components MAX_COMPONENTS
  --min-components MIN_COMPONENTS
                        Smallest number of latent variables CV considers;
                        candidates are MIN..MAX-1 (default: 1)
  --random-state RANDOM_STATE
  --out-table OUT_TABLE
                        Path to save the best models table (CSV)
  --out-vips OUT_VIPS   Path to save VIP scores per feature (CSV)
  --out-scores OUT_SCORES
                        Path to save latent space scores from final model
                        (CSV)
  --n-components N_COMPONENTS
                        Number of latent variables for final score model
                        (default: modal LV from CV)
```

**Example:**

```bash
motco plsr \
  --data data.csv --label-col group \
  --cv1-splits 5 --cv2-splits 5 --n-repeats 10 --max-components 5 \
  --out-table results/plsr_table.csv \
  --out-vips results/vips.csv
```

---

## `motco snf`

Similarity Network Fusion across two or more omics CSV files. Pass `--input` once per dataset (all must have the same number of rows in the same sample order).

```
usage: motco snf [-h] [--input INPUT] [--K K] [--eps EPS] [--k K] [--t T]
                 [--out-fused OUT_FUSED] [--out-embedding OUT_EMBEDDING]
                 [--spectral-components SPECTRAL_COMPONENTS]

options:
  -h, --help            show this help message and exit
  --input INPUT         Input CSV (repeat for multiple omics)
  --K K                 K for affinity construction
  --eps EPS             Epsilon for affinity construction
  --k K                 k for sparse kernel in SNF
  --t T                 Number of SNF iterations
  --out-fused OUT_FUSED
                        Path to save fused matrix (CSV)
  --out-embedding OUT_EMBEDDING
                        Path to save spectral embedding (CSV)
  --spectral-components SPECTRAL_COMPONENTS
                        Number of spectral embedding components (default: 10)
```

**Example:**

```bash
motco snf \
  --input proteomics.csv --input metabolomics.csv \
  --K 20 --k 20 --t 20 \
  --out-fused results/fused.csv \
  --out-embedding results/embedding.csv \
  --spectral-components 10
```

---

## `motco simulate`

Generate one semi-synthetic multi-omic dataset with the numpy generator (no R at runtime). Group A follows
a random baseline trajectory and group B the transform `--trajectory-mode` names, at `--effect-size`
(`0` is the null for every mode). See [Simulations](simulations.md#semi-synthetic-trajectory-generation)
for what each mode changes.

```
usage: motco simulate [-h] --seed SEED --out-dir OUT_DIR
                      [--n-samples N_SAMPLES] [--n-stages N_STAGES]
                      [--trajectory-mode {none,translation,magnitude,orientation,shape}]
                      [--effect-size EFFECT_SIZE]
                      [--shape-kind {relocate,magnitude}]
                      [--magnitude-kind {joint,all,extremes}] [--p-dmp P_DMP]
                      [--cluster-mean-shift CLUSTER_MEAN_SHIFT]
                      [--delta-methyl DELTA_METHYL] [--delta-expr DELTA_EXPR]
                      [--delta-protein DELTA_PROTEIN]

options:
  -h, --help            show this help message and exit
  --seed SEED           Random seed for reproducibility
  --out-dir OUT_DIR     Directory to write output files
  --n-samples N_SAMPLES
                        Total number of samples (default: 600)
  --n-stages N_STAGES   Number of trajectory stages (default: 3)
  --trajectory-mode {none,translation,magnitude,orientation,shape}
                        Group trajectory difference mode (default:
                        orientation)
  --effect-size EFFECT_SIZE
                        Group effect size: magnitude scale, extra-set size, or
                        relocated fraction (default: 1.0)
  --shape-kind {relocate,magnitude}
                        Shape-mode single-stage perturbation: relocate sites
                        or scale effect (default: relocate)
  --magnitude-kind {joint,all,extremes}
                        Magnitude-mode scope: joint (every omic's delta scaled
                        together), all (methylation delta alone at all
                        stages), or extremes (methylation endpoint indicators
                        only) (default: joint)
  --p-dmp P_DMP         Per-stage probability a methylation feature is
                        differential (default: 0.2)
  --cluster-mean-shift CLUSTER_MEAN_SHIFT
                        Per-omic mean shift (delta) applied to unspecified
                        per-omic delta flags (default: 2.0)
  --delta-methyl DELTA_METHYL
                        Methylation cluster mean shift; overrides --cluster-
                        mean-shift for methylation
  --delta-expr DELTA_EXPR
                        Expression cluster mean shift; overrides --cluster-
                        mean-shift for expression
  --delta-protein DELTA_PROTEIN
                        Proteomics cluster mean shift; overrides --cluster-
                        mean-shift for proteomics
```

The output directory holds `methylation.csv`, `expression.csv`, `proteomics.csv` and `metadata.csv`, the
trajectory design (`model_full.csv`, `model_reduced.csv`, `ls_means.csv`, `contrast.json`) ready for
`motco de`, and `truth.json` recording the generator parameters and the planted difference.

**Example:**

```bash
motco simulate \
  --seed 42 --n-samples 90 --n-stages 3 \
  --trajectory-mode orientation --effect-size 1.0 \
  --out-dir results/simulated/
```

---

## `motco de`

Estimate group trajectory differences, optionally with RRPP permutation testing.

**Estimate only** (no p-values): provide `--model-matrix`.

**RRPP** (with p-values): provide `--model-full`, `--model-reduced`, and `--rrpp-permutations`.

```
usage: motco de [-h] --Y Y --ls-means LS_MEANS --contrast CONTRAST
                [--model-matrix MODEL_MATRIX] [--model-full MODEL_FULL]
                [--model-reduced MODEL_REDUCED]
                [--rrpp-permutations RRPP_PERMUTATIONS] [--out-json OUT_JSON]
                [--out-observed OUT_OBSERVED]
                [--orientation-estimator {pc1,net_displacement}]

options:
  -h, --help            show this help message and exit
  --Y Y                 Outcome matrix CSV (latent space coordinates)
  --ls-means LS_MEANS   Least-squares means CSV
  --contrast CONTRAST   JSON file with groups (list of index lists)
  --model-matrix MODEL_MATRIX
                        Model matrix CSV (with intercept) for
                        estimate_difference
  --model-full MODEL_FULL
                        Full model matrix CSV (with intercept) for RRPP
  --model-reduced MODEL_REDUCED
                        Reduced model matrix CSV (with intercept) for RRPP
  --rrpp-permutations RRPP_PERMUTATIONS
                        Number of permutations for RRPP
  --out-json OUT_JSON   Output JSON file
  --out-observed OUT_OBSERVED
                        Save predicted LS-mean vectors as CSV
  --orientation-estimator {pc1,net_displacement}
                        Orientation estimator for angles, observed and
                        permuted (default: pc1). net_displacement is the
                        estimator the simulation harness and the published
                        operating characteristics describe; when given, the
                        JSON records it.
```

`--orientation-estimator` selects how each group's orientation is measured, for
the observed `angles` and for every RRPP permutation. The default stays `pc1`,
so existing invocations produce the same output. When the flag is omitted, the
JSON is unchanged and has no `orientation_estimator` key. `net_displacement` is
the estimator the simulation harness and the published operating
characteristics describe. When the flag is given, with either value, the JSON
records it as `"orientation_estimator"`. Only `angles` depend on it.

**Example — estimate only:**

```bash
motco de \
  --Y latent.csv --model-matrix model_full.csv \
  --ls-means ls_means.csv --contrast contrast.json \
  --out-json results/de.json \
  --out-observed results/ls_mean_vectors.csv
```

**Example — RRPP:**

```bash
motco de \
  --Y latent.csv \
  --model-full model_full.csv --model-reduced model_reduced.csv \
  --ls-means ls_means.csv --contrast contrast.json \
  --rrpp-permutations 999 \
  --out-json results/rrpp.json
```

---

## `motco analyze`

Real-data trajectory analysis with the validated pipeline. See
[Real-Data Analysis](analysis.md) for the input contract, the record and the
artifact layout.

```
usage: motco analyze [-h] --block NAME=PATH [--transform NAME=KIND] --metadata
                     METADATA --sample-col SAMPLE_COL --group-col GROUP_COL
                     --stage-col STAGE_COL [--stage-order STAGE_ORDER]
                     [--group-order GROUP_ORDER] [--covariate COL]
                     [--categorical-covariate COL]
                     [--missing-covariates {error,complete_case}]
                     [--latent-adjustment {stage,group_stage}]
                     [--permutations PERMUTATIONS] [--seed SEED]
                     [--n-jobs N_JOBS] [--forced-components FORCED_COMPONENTS]
                     [--residual-correction {df,leverage}] [--span-split-test]
                     [--interaction-scan] [--interaction-scan-by-block]
                     [--min-components MIN_COMPONENTS]
                     [--block-weighting {equal_variance,mfa}]
                     [--sparsity {joint,block}] [--sparsity-ladder F1,F2,...]
                     [--orientation-estimator {net_displacement,pc1}]
                     [--attribution]
                     [--bootstrap-replicates BOOTSTRAP_REPLICATES]
                     [--bootstrap-seed BOOTSTRAP_SEED] [--top-k TOP_K]
                     [--cv1-splits CV1_SPLITS] [--cv2-splits CV2_SPLITS]
                     [--n-repeats N_REPEATS] [--max-components MAX_COMPONENTS]
                     [--cv-random-state CV_RANDOM_STATE] --out-dir OUT_DIR
                     [--overwrite]

Measure two groups' trajectories across ordered stages on your own multi-omic
data. Blocks are analysed in the units supplied: MOTCO does not normalize
counts, filter features or impute values.

options:
  -h, --help            show this help message and exit
  --block NAME=PATH     Omic block CSV (samples x features, first column the
                        sample id); repeat per block, in the order the joint
                        feature matrix should follow
  --transform NAME=KIND
                        Unit transform for a declared block: identity
                        (default) or beta_to_mvalue
  --metadata METADATA   Sample metadata CSV
  --sample-col SAMPLE_COL
                        Metadata column with the sample ids
  --group-col GROUP_COL
                        Metadata column with the two groups
  --stage-col STAGE_COL
                        Metadata column with the stages
  --stage-order STAGE_ORDER
                        Required. Comma-separated stage labels, first to last
                        (never inferred by sorting)
  --group-order GROUP_ORDER
                        Comma-separated group labels (default: sorted as
                        strings)
  --covariate COL       Metadata column entering the trajectory model
                        additively; repeat per covariate, in model order (non-
                        numeric columns are categorical)
  --categorical-covariate COL
                        Like --covariate, but treat the column as a factor
                        whatever its dtype
  --missing-covariates {error,complete_case}
                        Missing covariate values: error (default) or
                        complete_case (drop those samples); never imputed
  --latent-adjustment {stage,group_stage}
                        Fit the PLS basis on features residualized on the
                        declared covariates, conditional on the stage (stage)
                        or group x stage (group_stage) indicators; requires a
                        covariate. Default: off (covariate-blind basis)
  --permutations PERMUTATIONS
                        RRPP permutations (default: 0, no p-values)
  --seed SEED           RRPP seed
  --n-jobs N_JOBS       Parallel jobs for CV and RRPP (-1 = all CPUs)
  --forced-components FORCED_COMPONENTS
                        Fix the latent rank and skip cross-validation
                        (dimensionality sensitivity)
  --residual-correction {df,leverage}
                        Rescale the reduced model's residuals before RRPP
                        permutes them: df (kappa = sqrt(n / (n - rank))) or
                        leverage (1 / sqrt(1 - h_ii)). Moves only the null.
                        Default: off
  --span-split-test     Also test the group x stage interaction of the
                        standardized features inside the fitted PLS span and
                        in its orthogonal complement (pseudo-F, RRPP draws at
                        --seed). Not available with --latent-adjustment.
                        Default: off
  --interaction-scan    Also test the group x stage interaction of every
                        standardized feature on its own (per-feature F, the
                        maximum's global p-value and step-down max-T adjusted
                        p-values, RRPP draws at --seed); writes
                        interaction_scan.csv. Default: off
  --interaction-scan-by-block
                        Split the interaction scan's familywise error equally
                        between the declared blocks (Bonferroni over within-
                        block step-downs), so a large or heavy-tailed block
                        does not set every block's bar. Requires
                        --interaction-scan. Default: off (one joint maximum)
  --min-components MIN_COMPONENTS
                        Rank floor of the cross-validated latent space: CV
                        considers ranks MIN..MAX-1 (default: 2; 1 restores the
                        unfloored 0.9.0 rule)
  --block-weighting {equal_variance,mfa}
                        Multiply every standardized block by one label-blind
                        scalar before the latent space is fitted:
                        equal_variance (1/sqrt(features), unit total variance
                        per block) or mfa (unit variance of each block's
                        leading principal component). delta is then in the
                        weighted metric. Unvalidated. Default: off
                        (unweighted)
  --sparsity {joint,block}
                        Fit sparse stage-PLS-DA: every component loads on a
                        budget of features, over all features (joint) or per
                        block (block); CV selects the rank and the keep
                        fraction jointly. Unvalidated. Default: off (ordinary
                        PLS-DA)
  --sparsity-ladder F1,F2,...
                        Keep fractions CV selects over, in (0, 1] (requires
                        --sparsity; default: 1.0,0.5,0.2,0.1,0.05)
  --orientation-estimator {net_displacement,pc1}
                        Orientation estimator (default: net_displacement)
  --attribution         Run frozen-model orientation attribution
  --bootstrap-replicates BOOTSTRAP_REPLICATES
                        Attribution bootstrap replicates
  --bootstrap-seed BOOTSTRAP_SEED
                        Attribution bootstrap seed
  --top-k TOP_K         Attribution top-k feature count
  --cv1-splits CV1_SPLITS
                        Inner CV folds (default: 3)
  --cv2-splits CV2_SPLITS
                        Outer CV folds (default: 4)
  --n-repeats N_REPEATS
                        Outer CV repeats (default: 5)
  --max-components MAX_COMPONENTS
                        Largest rank CV considers (default: 20)
  --cv-random-state CV_RANDOM_STATE
                        CV random state (default: 1203)
  --out-dir OUT_DIR     Directory for the analysis artifact
  --overwrite           Replace an existing analysis in --out-dir
```

Blocks are analysed in the units supplied: MOTCO does not normalize counts,
filter features or impute values. On success the command prints the selected
rank with its floor, the estimator, the statistics with their p-values and the
written paths;
with covariates it also prints each covariate's kind and the number of samples
complete-case deletion dropped, and with `--latent-adjustment` the conditioning
and the residual degrees of freedom of the residualization, and with
`--residual-correction` the correction it applied to the RRPP null (κ for `df`,
the leverage and multiplier range for `leverage`). With `--block-weighting` the
kind is printed beside the rank floor and recorded in a `block_weighting` block;
see [Block weighting](analysis.md#block-weighting) for the metric it implies and
its validation status. With `--sparsity` the scope and the selected keep
fraction are printed beside the rank floor and recorded in a `sparsity` block;
`--sparsity-ladder` without `--sparsity`, or a ladder that is not
comma-separated fractions in (0, 1], exits non-zero before anything is fitted.
See [Sparse stage-PLS](analysis.md#sparse-stage-pls) for its validation status.
With `--span-split-test` it prints one line with the three pseudo-F p-values
and the in-span share and records a `span_split` block; it exits non-zero with
`--latent-adjustment`. With `--interaction-scan` it prints one line with the
maximum F and its feature, the global p-value and the hits at adjusted
p ≤ 0.05, records an `interaction_scan` block and writes `interaction_scan.csv`;
`--interaction-scan-by-block` (which requires `--interaction-scan`) runs the
scan within each block and adds the combination and the number of blocks
tested to that line. See [Split interaction test](analysis.md#split-interaction-test)
and [Per-feature interaction scan](analysis.md#per-feature-interaction-scan).
Covariates are
declared in flag order across `--covariate` and `--categorical-covariate`.
`--latent-adjustment` without a covariate flag exits non-zero before any data
are read; see [Latent adjustment](analysis.md#latent-adjustment) for what it
changes and the mediator caveat, and [Residual correction](analysis.md#residual-correction)
for the corrected null and its validation status.
Input and validation errors exit non-zero with an `Error: …` message.

**Example:**

```bash
motco analyze \
  --block rna=rna_logcpm.csv --block atac=atac_logcpm.csv \
  --metadata donors.csv --sample-col donor_id --group-col Sex --stage-col ADNC \
  --stage-order "Not AD + Low,Intermediate,High" \
  --permutations 999 --seed 20261015 \
  --out-dir results/motco_mtg_astro
```
