# MOTCO — Multi-omics Trajectory Comparison

[![Docs](https://img.shields.io/badge/docs-GitHub%20Pages-blue)](https://tomszar.github.io/MOTCO/)

MOTCO provides tooling to generate latent spaces from multi‑omics data, test for group differences in multivariate trajectories within those spaces, and run reproducible semi-synthetic simulation studies.

Two approaches are implemented to build latent spaces:

- Partial Least Squares Regression/Discriminant Analysis (PLSR/PLS‑DA)
- Similarity Network Fusion (SNF) with optional spectral embedding

Once a latent space is constructed, MOTCO includes statistics to estimate differences in multivariate trajectories between groups (magnitude, orientation/angle, and shape), with an option for permutation testing via RRPP (Residual Randomization in a Permutation Procedure).

This repository contains the core statistical routines in `src/motco/stats`, simulation helpers in `src/motco/simulations`, and a command‑line interface for common analysis tasks.

## Real data

To analyse your own multi-omic data, use `motco analyze` (or `motco.analysis.analyze_trajectories`
from Python). It runs the pipeline the simulation studies validated, in one step: pooled standardization,
a stage-supervised PLS-DA latent space with the double-CV rank, `delta`/`angle`/`shape` with the
net-displacement orientation estimator, and RRPP. It also writes a reproducible artifact (record, latent
scores, loadings, stage means). You declare the stage order and each block's units explicitly. See
[docs/api/analysis.md](docs/api/analysis.md).

## Quick start

The bundled toy dataset at `examples/data/toy/` lets you run the full pipeline immediately without R:

```bash
# 1. Build a PLS-DA latent space supervised by disease stage
motco plsr \
  --input examples/data/toy/methylation.csv \
  --input examples/data/toy/expression.csv \
  --input examples/data/toy/proteomics.csv \
  --metadata examples/data/toy/metadata.csv \
  --label-col stage \
  --cv1-splits 3 \
  --cv2-splits 3 \
  --n-repeats 3 \
  --max-components 5 \
  --out-scores results/latent_pls.csv

# 2. (Alternative) Build an SNF latent space
motco snf \
  --input examples/data/toy/methylation.csv \
  --input examples/data/toy/expression.csv \
  --input examples/data/toy/proteomics.csv \
  --K 10 --k 10 --t 20 \
  --out-embedding results/latent_snf.csv

# 3. Run differential trajectory analysis with RRPP
motco de \
  --Y results/latent_pls.csv \
  --model-full examples/data/toy/model_full.csv \
  --model-reduced examples/data/toy/model_reduced.csv \
  --ls-means examples/data/toy/ls_means.csv \
  --contrast examples/data/toy/contrast.json \
  --rrpp-permutations 999 \
  --out-json results/de_result.json

# 4. (Optional) Generate a fresh semi-synthetic dataset (numpy generator; no R at runtime)
motco simulate \
  --seed 42 --n-samples 90 --n-stages 3 \
  --trajectory-mode orientation --effect-size 1.0 \
  --p-dmp 0.2 \
  --cluster-mean-shift 0.10 \
  --out-dir results/simulated/
```

The bundled toy data uses moderate InterSIM cluster separation, so `y = stage` classification should be informative but not perfect. The small CV settings above keep the quick start fast; increase them for analysis-grade model selection.

R and InterSIM are **not** required at runtime — `motco simulate` uses the numpy generator with the cached InterSIM reference. R + InterSIM are only needed to regenerate that cache or the fidelity fixtures (see `src/motco/simulations/FIDELITY.md`). To install InterSIM in R for that one-time step:
```r
install.packages("InterSIM", repos = c("https://cran.r-universe.dev", "https://cloud.r-project.org"))
```

See `examples/motco_example.ipynb` for a full end-to-end walkthrough in Python.

## Install (with uv)

Prerequisites: Python 3.11+ and [uv](https://github.com/astral-sh/uv) installed.

```bash
# Create and activate a virtual environment
uv venv
source .venv/bin/activate

# Install MOTCO in editable mode (and dependencies)
uv pip install -e .

# Verify CLI is available
motco --help
```

Alternatively, using pip:

```bash
python -m venv .venv
source .venv/bin/activate
pip install -e .
```

## Quick Example

See `examples/motco_example.ipynb` for an end-to-end walkthrough using the bundled dataset, including construction of the latent/design inputs used by trajectory tests.

Representative CLI commands:

```bash
# 1. Run PLS-DA double cross-validation and save the model-selection table
motco plsr --data tests/data/evo_649_sm_example1.csv --label-col taxa \
  --cv1-splits 7 --cv2-splits 8 --n-repeats 5 --max-components 2 \
  --out-table results/plsr_table.csv

# 2. Estimate group differences after preparing aligned Y/design/contrast files
motco de \
  --Y results/latent_space.csv \
  --model-matrix results/model_matrix.csv \
  --ls-means results/ls_means.csv \
  --contrast contrast.json \
  --out-json results/de_result.json \
  --out-observed results/ls_mean_vectors.csv
```

## Command Line Interface

MOTCO exposes a single entry‑point `motco` with subcommands for PLSR/PLS‑DA, SNF, and Differential Effects (group differences).

### 1) PLS‑DA with double cross‑validation

```bash
motco plsr \
  --data path/to/table.csv \
  --label-col diagnosis \
  --cv1-splits 7 --cv2-splits 8 --n-repeats 30 --max-components 50 \
  --out-table results/plsr_models.csv
```

Options:
- Use `--data` for a single CSV containing predictors and a label column specified by `--label-col`.
- Or provide separate matrices via `--x` and `--y` CSV files (mutually exclusive with `--data`).
- If `--y` has a single column, it is treated as a label vector and will be one‑hot encoded internally; if it has multiple columns, it is treated as an already encoded class matrix.
- Outputs a table with one row per outer CV repeat: `rep`, `LV` (modal n_LV across K outer folds), `AUROC` (mean across K outer folds), and `AUROC_std` (sample std across K outer folds). The actual trained models are kept in memory; export of models is not included in the CLI at this time. If `--out-table` is omitted, the table is printed to stdout.

Input expectations:
- CSV files with samples in rows, features in columns. For `--data`, include a label column with binary or multi‑class outcome; it will be one‑hot encoded internally.

### 2) Similarity Network Fusion (SNF)

```bash
motco snf \
  --input omics1.csv --input omics2.csv [--input omics3.csv ...] \
  --K 20 --eps 0.5 --k 20 --t 20 \
  --out-fused fused_affinity.csv \
  --out-embedding spectral_embedding.csv
```

Notes:
- Each `--input` CSV must contain the same samples in the same order (rows = samples).
- `--K` and `--eps` are used when constructing per‑dataset affinity matrices; `--k` and `--t` control SNF neighborhood size and iterations.
- The fused similarity matrix is saved to `--out-fused`. If `--out-embedding` is provided, a spectral embedding is also computed and saved. Use `--spectral-components` to control its dimensionality; the default is 10. If no output paths are provided, the fused matrix is printed to stdout.

### 3) Differential Effects (group differences)

Estimate differences in magnitude and direction between groups using least‑squares means, with optional permutation testing (RRPP).

```bash
motco de \
  --Y latent_space.csv \
  --model-matrix model_matrix.csv \
  --ls-means ls_means.csv \
  --contrast contrast.json \
  --out-json de_result.json

# With permutations (RRPP)
motco de \
  --Y latent_space.csv \
  --model-full model_full.csv \
  --model-reduced model_reduced.csv \
  --ls-means ls_means.csv \
  --contrast contrast.json \
  --rrpp-permutations 999 \
  --out-json rrpp_result.json
```

Where:
- `latent_space.csv` is the outcome matrix `Y` (e.g., coordinates in a latent space; rows = samples, columns = dimensions).
- `model_matrix.csv` is a design matrix (with intercept) aligned to `Y` rows. For RRPP, provide `--model-full` and `--model-reduced` (both with intercept).
- `ls_means.csv` contains the least‑squares means to compare (rows = groups/cohorts, columns = same dimensions as `Y`).
- `contrast.json` is a JSON array of index lists, where each inner list enumerates the cohort indices belonging to the same group. Example: `[[0,1],[2,3]]`.
- Output JSON includes `deltas`, `angles`, and `shapes` as symmetric matrices (lists of lists). With `--rrpp-permutations > 0`, these are returned as lists of matrices per permutation: e.g., `deltas[perm_idx] -> matrix`.

Important: trajectory design, estimation, and permutation utilities now live in focused submodules under `motco.stats`. See the Python API notes below for helpers to build model matrices and LS means directly from your `group` and `level` columns.

## Python API (selected)

```python
from motco.stats.pls import plsda_doubleCV
from motco.stats.snf import get_affinity_matrix, SNF, get_spectral
from motco.stats.design import (
    center_matrix,
    get_model_matrix,
    build_ls_means,
)
from motco.stats.permutation import RRPP
from motco.stats.trajectory import (
    estimate_difference,
    get_observed_vectors,
    pair_difference,
)

# Example: build model matrix and LS means from group/level columns
import pandas as pd
import numpy as np

# X_factors contains two categorical columns: 'group' and 'level'
X_factors = pd.DataFrame({
    'group': ['A','A','B','B'],
    'level': ['t0','t1','t0','t1'],
})

# Y is the outcome/feature matrix aligned by rows to X_factors (e.g., latent space)
Y = pd.DataFrame(np.random.randn(4, 3), columns=['z1','z2','z3'])

# Build design and estimate LS means for all group×level cells
X = get_model_matrix(X_factors, group_col='group', level_col='level', full=True)
ls = build_ls_means(
    group_levels=sorted(X_factors['group'].astype(str).unique()),
    level_levels=sorted(X_factors['level'].astype(str).unique()),
    full=True,
)
deltas, angles, shapes = estimate_difference(Y=Y.values, model_matrix=X, LS_means=ls, contrast=[[0,1],[2,3]])

# Two‑state comparison between two groups at two levels (angle & delta)
df = X_factors.copy()
df = pd.concat([df, Y], axis=1)
angle_deg, delta_mag = pair_difference(df, group_col='group', level_col='level', feature_cols=['z1','z2','z3'])

# Center features within groups (optional preprocessing)
df_centered = center_matrix(df, group_col='group', level_col='level', feature_cols=['z1','z2','z3'])

# RRPP with parallelism from Python API (optional)
# deltas_list, angles_list, shapes_list = RRPP(
#     Y=Y.values,
#     model_full=X,
#     model_reduced=X[:, :3],  # example reduced model
#     LS_means=ls,
#     contrast=[[0,1],[2,3]],
#     permutations=999,
#     n_jobs=-1,  # use all CPUs
# )
```

See inline docstrings in the modules under `src/motco/stats/` for full details.

### Inspecting LS-mean coordinates

Before running `estimate_difference`, use `get_observed_vectors` to see the predicted
mean position of each group × level cell in Y space:

```python
from motco.stats.trajectory import get_observed_vectors

# X_factors: DataFrame with group_col and level_col
# Y: outcome matrix aligned to X_factors by row
obs = get_observed_vectors(X_factors, Y, group_col='group', level_col='level', full=True)
# Returns a DataFrame with MultiIndex (group, level) and columns matching Y
print(obs)
```

## Interpreting Results

`estimate_difference` and `RRPP` return three symmetric matrices:

| Output | Meaning |
|--------|---------|
| `deltas` | Absolute difference in trajectory magnitude (total path length) between group pairs. Larger = one group changed more than the other. |
| `angles` | Angle in degrees between trajectory orientations. 0° = same direction; 90° = orthogonal; 180° = exactly opposite. |
| `shapes` | Procrustes distance between trajectory shapes after removing location, uniform scale, and any orthogonal transformation. Reflections are aligned away, at every ambient dimension, so a mirror pair scores 0 — a genuine mirror difference is an orientation difference and shows up in `angles`. 0 = identical shape. |

**P-values via RRPP:** Use a right-tailed test with the add-one correction:

```python
def pvalue(samples, observed, i, j):
    vals = np.array([s[i, j] for s in samples])
    return (np.sum(vals >= observed) + 1) / (len(vals) + 1)
```

Significance threshold is typically α = 0.05.

## Simulation API

MOTCO includes packaged simulation tools for generating semi-synthetic trajectory datasets, evaluating one generated dataset through MOTCO statistics, and orchestrating small Type I error or power grids.

Generation is numpy-native (no R at runtime): the generator samples from the cached InterSIM reference committed at `src/motco/simulations/data/intersim_reference.npz`. The legacy R bridge (`check_intersim_available`) is retained only for the one-time export that builds that cache.

Generate and evaluate one semi-synthetic trajectory dataset:

```python
from motco.simulations import (
    SemiSyntheticTrajectoryParams,
    SimulationEvaluationParams,
    evaluate_semisynthetic_trajectory,
    generate_semisynthetic_trajectory,
)

dataset = generate_semisynthetic_trajectory(
    SemiSyntheticTrajectoryParams(
        seed=99,
        n_samples=120,
        n_stages=3,
        trajectory_mode="magnitude",
        group_effect_size=0.2,
        group_ratio=0.5,
        p_dmp=0.2,
    ),
)

result = evaluate_semisynthetic_trajectory(
    dataset,
    SimulationEvaluationParams(integration_method="concat", permutations=0),
)
print(result.pair_statistics)
```

Enumerate and run a small local grid. Results are persisted as JSON Lines, one row per cell/replicate, and can be resumed safely via parameter signatures:

```python
from pathlib import Path

from motco.simulations import (
    SemiSyntheticTrajectoryParams,
    SimulationEvaluationParams,
    SimulationRunConfig,
    enumerate_type_i_grid,
    run_simulation_grid,
    summarize_rejection_rates,
)

grid = enumerate_type_i_grid(
    baseline_generator_params=SemiSyntheticTrajectoryParams(seed=2, n_samples=60),
    evaluation_params=SimulationEvaluationParams(integration_method="concat", permutations=99),
    axes={"generator.n_samples": [60, 120], "generator.group_ratio": [0.5, 0.7]},
    n_replicates=3,
    base_seed=2026,
)

records = run_simulation_grid(
    grid,
    config=SimulationRunConfig(output_path=Path("simulation-results.jsonl")),
)
summaries = summarize_rejection_rates(records, alpha=0.05)
```

## Trajectory power study

For a full, declarative, sharded study that produces the per-statistic
specificity matrix, Type I tables, and power curves (the deliverables of the
methods paper), use the `motco.simulations.study` package and the helper
scripts under `scripts/`.

End-to-end workflow:

```bash
# (1) Run shards locally (or via SLURM array; see scripts/motco_study_array.sbatch)
python scripts/run_study_shard.py \
    --config examples/trajectory_power_study/smoke.json \
    --out-dir results/ \
    --shard-index 0 --n-shards 4 --error-policy record

# (2) Merge per-shard JSONL into a single deduplicated result set
python scripts/motco_study.py merge --out-dir results/

# (3) Build per-statistic + combined-rule summaries, the specificity matrix,
#     power curves, Type I tables, figures, and an acceptance-target report
python scripts/motco_study.py report \
    --config examples/trajectory_power_study/smoke.json \
    --out-dir results/
```

See `examples/trajectory_power_study/README.md` for the config schema, the
SLURM array submission flow, and resubmission of failed array tasks.

## License

This project is licensed under the terms of the LICENSE file included in this repository.
