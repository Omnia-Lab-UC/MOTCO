# MOTCO

**Multi-omics Trajectory Comparison** — a Python package for building latent spaces from multi-omics data and testing whether two groups follow different molecular trajectories in those spaces.

## Start here

To analyse your own data, use [`motco analyze`](api/cli.md#motco-analyze) or `motco.analysis.analyze_trajectories` ([Real-Data Analysis](api/analysis.md)). It runs the measurement the simulation studies validated, in one step: pooled standardization, a stage-supervised PLS-DA latent space with the double-CV rank, `delta` / `angle` / `shape` with the net-displacement orientation estimator, and RRPP p-values.

- [Measurement visibility](measurement-visibility.md) — which group differences the stage-supervised PLS measurement can see, and the tests for those it cannot
- [Scientific roadmap](roadmap.md) — what is established, what is not, and what comes next
- [End-to-end example notebook](https://github.com/Omnia-Lab-UC/MOTCO/blob/main/examples/motco_example.ipynb)
- [GitHub repository](https://github.com/Omnia-Lab-UC/MOTCO)

## Installation

```bash
pip install motco
# or with uv:
uv add motco
```

## Modules

| Module | Purpose |
|--------|---------|
| [Real-Data Analysis](api/analysis.md) | `analyze_trajectories`: the validated pipeline on your own blocks, metadata and covariates, with its record and artifact |
| [PLS-DA](api/pls.md) | Supervised latent spaces via double cross-validated PLS discriminant analysis, including sparse stage-PLS |
| [SNF](api/snf.md) | Fuse multiple omics layers into a single similarity network |
| [Trajectory Analysis](api/sd.md) | Estimate and test differences in trajectory magnitude, orientation and shape; the split interaction test and the per-feature interaction scan |
| [Orientation Driver Attribution](api/attribution.md) | Explain shared-PLS trajectory orientation differences with signed feature, module and pathway effects |
| [Simulations](api/simulations.md) | The numpy semi-synthetic generator (no R at runtime), the evaluation harness and the declarative power study |
| [CLI](api/cli.md) | `motco plsr`, `motco snf`, `motco simulate`, `motco de`, `motco analyze` |
