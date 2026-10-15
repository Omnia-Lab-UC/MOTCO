# CLI Reference

The `motco` command-line interface provides five subcommands.

```
usage: motco [-h] [--version] [--verbose] {plsr,snf,simulate,de,analyze} ...

MOTCO CLI: PLSR, SNF, group differences, and real-data trajectory analysis

positional arguments:
  {plsr,snf,simulate,de,analyze}
    plsr                Run PLS-DA with double cross-validation
    snf                 Similarity Network Fusion
    simulate            Generate a semi-synthetic multi-omics toy dataset (numpy generator, no R)
    de                  Group differences on trajectories (estimate or RRPP)
    analyze             Real-data trajectory analysis with the validated pipeline (PLS latent space, RRPP)
```

For real data, use `motco analyze`. It runs the measurement the simulation
studies validated, in one step. `plsr` and `de` are lower-level building blocks.

---

## `motco plsr`

Run PLS-DA with double cross-validation. Input can be a single CSV (with `--data` and `--label-col`) or separate predictor/label CSVs (`--x` / `--y`).

```
usage: motco plsr [-h] [--data DATA] [--label-col LABEL_COL] [--x X] [--y Y]
                  [--cv1-splits CV1_SPLITS] [--cv2-splits CV2_SPLITS]
                  [--n-repeats N_REPEATS] [--max-components MAX_COMPONENTS]
                  [--random-state RANDOM_STATE] [--out-table OUT_TABLE]
                  [--out-vips OUT_VIPS]

options:
  --data DATA                CSV with predictors and label column
  --label-col LABEL_COL      Label column name when using --data
  --x X                      CSV with predictors (features)
  --y Y                      CSV with labels/outcomes
  --cv1-splits CV1_SPLITS    Inner CV folds
  --cv2-splits CV2_SPLITS    Outer CV folds
  --n-repeats N_REPEATS      Number of outer CV repeats
  --max-components MAX_COMPONENTS
                             Maximum latent variables to search
  --random-state RANDOM_STATE
                             Random seed for reproducibility
  --out-table OUT_TABLE      Path to save the best models table (CSV)
  --out-vips OUT_VIPS        Path to save VIP scores per feature (CSV)
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
  --input INPUT              Input CSV (repeat for multiple omics)
  --K K                      K for affinity construction (default: 20)
  --eps EPS                  Epsilon for affinity construction (default: 0.5)
  --k K                      k for sparse kernel in SNF (default: 20)
  --t T                      Number of SNF iterations (default: 20)
  --out-fused OUT_FUSED      Path to save fused matrix (CSV)
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
  --Y Y                      Outcome matrix CSV (latent space coordinates)
  --ls-means LS_MEANS        Least-squares means CSV
  --contrast CONTRAST        JSON file with groups (list of index lists)
  --model-matrix MODEL_MATRIX
                             Model matrix CSV (with intercept) for estimate_difference
  --model-full MODEL_FULL    Full model matrix CSV (with intercept) for RRPP
  --model-reduced MODEL_REDUCED
                             Reduced model matrix CSV (with intercept) for RRPP
  --rrpp-permutations RRPP_PERMUTATIONS
                             Number of permutations for RRPP
  --out-json OUT_JSON        Output JSON file
  --out-observed OUT_OBSERVED
                             Save predicted LS-mean vectors as CSV
  --orientation-estimator {pc1,net_displacement}
                             Orientation estimator for angles (default: pc1)
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
                     [--group-order GROUP_ORDER] [--permutations PERMUTATIONS]
                     [--seed SEED] [--n-jobs N_JOBS]
                     [--forced-components FORCED_COMPONENTS]
                     [--orientation-estimator {net_displacement,pc1}]
                     [--attribution]
                     [--bootstrap-replicates BOOTSTRAP_REPLICATES]
                     [--bootstrap-seed BOOTSTRAP_SEED] [--top-k TOP_K]
                     [--cv1-splits CV1_SPLITS] [--cv2-splits CV2_SPLITS]
                     [--n-repeats N_REPEATS] [--max-components MAX_COMPONENTS]
                     [--cv-random-state CV_RANDOM_STATE] --out-dir OUT_DIR
                     [--overwrite]

options:
  --block NAME=PATH          Omic block CSV (samples x features, first column the sample id);
                             repeat per block, in joint feature order
  --transform NAME=KIND      Unit transform for a declared block: identity (default) or beta_to_mvalue
  --metadata METADATA        Sample metadata CSV
  --sample-col, --group-col, --stage-col
                             Metadata columns: sample ids, the two groups, the stages
  --stage-order STAGE_ORDER  Required. Comma-separated stage labels, first to last
  --group-order GROUP_ORDER  Comma-separated group labels (default: sorted as strings)
  --permutations N           RRPP permutations (default: 0, no p-values)
  --seed SEED                RRPP seed
  --n-jobs N_JOBS            Parallel jobs for CV and RRPP (-1 = all CPUs)
  --forced-components K      Fix the latent rank and skip cross-validation
  --orientation-estimator {net_displacement,pc1}
                             Orientation estimator (default: net_displacement)
  --attribution              Run frozen-model orientation attribution
  --bootstrap-replicates, --bootstrap-seed, --top-k
                             Attribution bootstrap settings
  --cv1-splits 3, --cv2-splits 4, --n-repeats 5, --max-components 20, --cv-random-state 1203
                             Cross-validation settings (harness defaults)
  --out-dir OUT_DIR          Directory for the analysis artifact
  --overwrite                Replace an existing analysis in --out-dir
```

Blocks are analysed in the units supplied: MOTCO does not normalize counts,
filter features or impute values. On success the command prints the selected
rank, the estimator, the statistics with their p-values and the written paths.
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
