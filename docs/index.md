# MOTCO

**Multi-omics Trajectory Comparison** — a Python package for building latent spaces from multi-omics data and quantifying group differences in those spaces.

## What it does

MOTCO provides three statistical modules:

| Module | Purpose |
|--------|---------|
| [PLS-DA](api/pls.md) | Build supervised latent spaces via double cross-validated PLS discriminant analysis |
| [SNF](api/snf.md) | Fuse multiple omics layers into a single similarity network |
| [Trajectory Analysis](api/sd.md) | Estimate and test differences in trajectory magnitude, orientation, and shape between groups |
| [Orientation Driver Attribution](api/attribution.md) | Explain shared-PLS trajectory orientation differences with signed feature, module, and pathway effects |
| [Simulations](api/simulations.md) | Invoke R InterSIM and return aligned semi-synthetic omics matrices |

## Quick links

- [GitHub repository](https://github.com/Omnia-Lab-UC/MOTCO)
- [End-to-end example notebook](https://github.com/Omnia-Lab-UC/MOTCO/blob/main/examples/motco_example.ipynb)
- [CLI reference](api/cli.md)

## Installation

```bash
pip install motco
# or with uv:
uv add motco
```

## API Reference

- **[PLS-DA](api/pls.md)** — `plsda_doubleCV`, `calculate_vips`
- **[SNF](api/snf.md)** — `get_affinity_matrix`, `SNF`, `get_spectral`
- **[Trajectory Analysis](api/sd.md)** — `get_model_matrix`, `build_ls_means`, `estimate_difference`, `RRPP`, `estimate_betas`, `get_observed_vectors`, `pair_difference`, `center_matrix`
- **[Orientation Driver Attribution](api/attribution.md)** — `analyze_orientation_attribution`, `attribution_frames`, `write_attribution_outputs`
- **[Simulations](api/simulations.md)** — `InterSIMParams`, `InterSIMResult`, `check_intersim_available`, `run_intersim`
- **[CLI](api/cli.md)** — `motco plsr`, `motco snf`, `motco de`
