# latent-space-metric-compatibility Specification

## Purpose
Document the relationship between MOTCO's Euclidean trajectory statistics and
the geometry produced by each supported integration method.
## Requirements
### Requirement: Measurement-space assumptions are explicit

The documentation SHALL state that `delta`, `angle`, and Procrustes `shape`
measure Euclidean geometry in the supplied outcome matrix and SHALL distinguish
standardized concatenation, linear PCA/PLS projections, and SNF spectral
coordinates.

#### Scenario: User selects a linear projection

- **WHEN** a user interprets trajectory statistics calculated in PCA or PLS scores
- **THEN** the documentation explains that the scores are a lossy linear projection whose geometry can be related to the original features through fitted loadings

#### Scenario: User selects SNF

- **WHEN** a user interprets trajectory statistics calculated in SNF spectral coordinates
- **THEN** the documentation warns that feature-space path length, direction, and shape are not guaranteed to be preserved and recommends simulation validation

### Requirement: Graph-native alternatives are scoped accurately

The documentation SHALL identify diffusion, connectivity, neighbourhood, and
transition-profile quantities as geometrically aligned candidates for SNF and
SHALL state that MOTCO does not currently implement them.

#### Scenario: User seeks an SNF-native trajectory interpretation

- **WHEN** the user consults the SNF documentation
- **THEN** candidate graph-native metric families are listed without implying that they are available in the current API

### Requirement: Pilot claims distinguish integration from metric compatibility

The latent-space cross-talk findings SHALL describe the observed SNF result as
a property of SNF spectral coordinates combined with MOTCO's Euclidean
statistics, rather than a general failure of SNF as a multi-omics integration
method.

#### Scenario: Reader reviews the SNF pilot finding

- **WHEN** the reader encounters the reported lack of `delta`/`angle` sensitivity
- **THEN** the findings identify the metric mismatch and note that graph-native alternatives were not evaluated

### Requirement: Visibility of group differences in the stage-supervised measurement is documented

The documentation SHALL state that the production stage-supervised PLS measurement sees a group difference only through its component inside the fitted basis's span. It SHALL state which differences that span contains, why a rotation of the stage program onto other features can cancel, and what a non-rejection therefore means. Each claim SHALL link to the study that measured it.

#### Scenario: Reader asks which differences the measurement can see
- **WHEN** the reader consults the measurement-visibility documentation
- **THEN** it explains that differences re-weighting the stage program both groups share lie inside the span and are visible
- **AND** it explains that at rank K − 1 with balanced cells the span is that of the pooled stage means, so an equal-length rotation of the program onto other features is perpendicular to it (the rhombus argument) and cancels

#### Scenario: Reader asks how out-of-span differences can still register
- **WHEN** the reader consults the same documentation
- **THEN** it lists the incidental routes (CV rank above K − 1, unequal cell sizes or feature scales, correlated noise, nonlinear units) and states that through them a rotation registers as a mix of `delta`, `angle` and `shape` rather than as `angle` alone
- **AND** it cites the orientation probe's measurement that the label-free `pca` and the two-sided `pls_within_group` spaces keep the full rotation, with the costs that probe recorded

#### Scenario: Reader interprets a non-rejection
- **WHEN** the reader interprets a non-rejection of `angle` or `shape`
- **THEN** the documentation states that it means no detectable difference inside the stage-supervised span, not the same direction or shape in feature space
- **AND** it cites the span-capture shares measured for every simulator construction
