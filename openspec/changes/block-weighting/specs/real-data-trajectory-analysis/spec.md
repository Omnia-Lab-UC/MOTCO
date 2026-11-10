# Spec Delta

## ADDED Requirements

### Requirement: Real-data analysis offers block weighting

The analysis parameters SHALL accept a block weighting kind (`block_weighting`: absent, `equal_variance` or `mfa`). When set, every declared block is weighted after standardization by the kind's scalar, the latent space is fitted in the weighted metric, and the weight enters the original-scale mapping used by attribution and exposed in the in-memory measurement, so reconstructed features are in their original units. An unknown kind SHALL be rejected before any fitting. The setting SHALL be omitted from the resolved parameters while absent.

#### Scenario: Weighted analysis
- **WHEN** an analysis declares two blocks of unequal feature counts and requests `equal_variance`
- **THEN** the latent space is fitted on blocks of unit total variance and the attribution's reconstructed features are in original units

#### Scenario: Absent weighting
- **WHEN** no weighting is requested
- **THEN** every output equals the pre-change output, and the resolved parameters carry no weighting key

### Requirement: Real-data analysis records block weighting

When weighting is set, the record SHALL carry a `block_weighting` block after `block_units` with the kind and, per block, the weight and the feature count, and `latent.integration_params.block_weighting`. Neither SHALL appear while absent.

#### Scenario: Record contents
- **WHEN** a weighted analysis is recorded
- **THEN** `block_weighting` names the kind and lists every block's weight and feature count, and the latent block's integration parameters name the kind

### Requirement: motco analyze exposes block weighting

`motco analyze` SHALL accept `--block-weighting {equal_variance,mfa}` and print the kind beside the selected rank when given.

#### Scenario: CLI round trip
- **WHEN** `motco analyze --block-weighting mfa` runs
- **THEN** the written record's `block_weighting.kind` is `mfa`
