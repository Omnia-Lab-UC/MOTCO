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

#### Scenario: Reader asks when the rhombus holds beyond two stages
- **WHEN** the reader consults the rhombus argument for three or more stages
- **THEN** the documentation states that the cancellation is exact when B's steps keep the lengths and mutual angles of A's and are orthogonal to every one of A's steps (B's program moved onto features A's does not use)
- **AND** it states that a rotation within A's own features is partly or wholly visible, the stage-span rotation being the wholly visible case

#### Scenario: Reader asks how out-of-span differences can still register
- **WHEN** the reader consults the same documentation
- **THEN** it lists the incidental routes (CV rank above K − 1, unequal cell sizes or feature scales, correlated noise, nonlinear units) and states that through them a rotation registers as a mix of `delta`, `angle` and `shape` rather than as `angle` alone
- **AND** it cites the orientation probe's measurement that the label-free `pca` and the two-sided `pls_within_group` spaces keep the full rotation, with the costs that probe recorded

#### Scenario: Reader compares a latent space with concat
- **WHEN** the reader consults the statement that every linear latent space is a projection of `concat`
- **THEN** the documentation states that a projection cannot see a difference `concat` does not contain, and that it can still have more power than `concat` by discarding noise directions, citing the probe's Phase 6 `shape` power (0.19 for `pls_cv` against 0.00 for `concat`)

#### Scenario: Reader interprets a non-rejection
- **WHEN** the reader interprets a non-rejection of `angle` or `shape`
- **THEN** the documentation states that it means no detectable difference inside the stage-supervised span, not the same direction or shape in feature space
- **AND** it cites the span-capture shares measured for every simulator construction

### Requirement: Each candidate space's visibility profile is documented

The measurement-visibility documentation SHALL report, for each space in the latent-space comparison, the span capture and the target power of the stage-span orientation as its alignment with the stage span moves from 0° to 90°. It SHALL report each space's size at the controls and the cost on protected statistics, at every design, each linked to the dated report.

#### Scenario: Reader chooses a measurement space for a question
- **WHEN** the reader consults the documentation to decide which space matches a question
- **THEN** it shows, for each space, how capture and orientation power change with alignment, its anchor and translation size, and its in-span `angle`, `shape` and `delta` power against the production space

#### Scenario: Reader asks whether blindness depends on the features touched
- **WHEN** the reader consults the sparse-versus-dense contrast at full alignment
- **THEN** the documentation states whether the comparison found capture and power to depend on the out-of-span direction's support at a fixed norm, with the measured difference and its standard error

### Requirement: The score-coordinate metric of the PLS measurement is documented

The documentation SHALL state that `delta`, `angle` and `shape` are measured in PLS score coordinates, `X·W·(PᵀW)⁻¹`, a sheared image of the fitted span rather than an orthonormal projection. It SHALL state what this does to size and to specificity, and cite the measured residuals. It SHALL report the validation of orthonormal span coordinates against the scores, and the production consequence.

#### Scenario: Reader asks what coordinates the statistics use
- **WHEN** the reader consults the measurement-visibility documentation
- **THEN** it states that the production statistics are computed in the PLS score coordinates, and that these coordinates do not preserve the standardized feature space's lengths and angles within the span

#### Scenario: Reader asks whether size is affected
- **WHEN** the reader asks whether the score coordinates affect Type I error
- **THEN** the documentation states that under the null both groups share one population trajectory, which any linear map keeps equal, so size is unaffected

#### Scenario: Reader asks whether specificity is affected
- **WHEN** the reader asks whether a shape-only or orientation-only difference can register in another statistic
- **THEN** the documentation states that a difference pure in the standardized feature space can move `delta`, `angle` or `shape` in score coordinates
- **AND** it cites the leak decomposition's score-coordinate residuals per column (planted-scale bend, e = 0.5), including that they are largest on the SEA-AD table and at four stages

#### Scenario: Reader asks what would remove the residual
- **WHEN** the reader asks how the score-coordinate route could be avoided
- **THEN** the documentation names measuring in orthonormal span coordinates, states that every orthonormal basis of the span gives the same statistics, and that the harness offers it as an opt-in setting
- **AND** it states the validation's verdict, its size, specificity and power readings against the scores on matched units, and the paired rates of the watched cells

#### Scenario: Reader asks which coordinates production uses
- **WHEN** the reader asks whether the production measurement and `motco analyze` use orthonormal coordinates
- **THEN** the documentation states that they use score coordinates until a separate change adopts orthonormal coordinates, which only the `adopt_orthonormal` verdict permits, and links the roadmap item that holds it

### Requirement: The split interaction test and its reading are documented

The measurement-visibility documentation SHALL describe the split interaction test:
- what each part tests;
- that "span" means the span fitted to the sample, not the population stage program;
- how a rejection and a non-rejection of each part are read;
- why the statistic is a ratio;
- what the test does not provide.

It SHALL cite the validation's verdict and measured rates.

#### Scenario: Reader asks what a complement rejection means
- **WHEN** the reader consults the documentation about a `complement` rejection
- **THEN** it states that the groups differ in a group × stage pattern with energy outside the fitted span, not necessarily that one group uses a different molecular program
- **AND** it states that at a weak stage signal this includes differences along the true stage program that the fitted span misses, citing the comparison's weak-design capture of the population stage program (0.10)

#### Scenario: Reader asks what a span rejection means
- **WHEN** the reader consults the documentation about a `span` rejection
- **THEN** it states that the difference has energy inside the fitted span, where `delta`, `angle` and `shape` decompose it

#### Scenario: Reader asks what a non-rejection means
- **WHEN** the reader asks what a non-rejection of the `complement` or `total` part says
- **THEN** the documentation states that it is uninformative unless power is shown at the design and kind of difference, and cites the validation's sparse-plant power on the SEA-AD blocks

#### Scenario: Reader asks why a ratio statistic is used
- **WHEN** the reader asks why the test uses a pseudo-F rather than the interaction energy alone
- **THEN** the documentation states that permuting the reduced model's residuals shrinks their scale, that an energy statistic summed over many coordinates turns that shrinkage into a liberal null, and that the ratio cancels it, citing the validation's diagnostic

#### Scenario: Reader asks what the test does not provide
- **WHEN** the reader asks whether the split test replaces `delta`, `angle` or `shape`
- **THEN** the documentation states that it gives no geometry for an out-of-span difference, only its presence and location relative to the fitted span, and links the roadmap items for group-aware spaces and sparse interaction scans

### Requirement: The per-feature interaction scan and its reading are documented

The measurement-visibility documentation SHALL describe the per-feature interaction scan:
- what it tests;
- how a hit and a non-rejection are read;
- why it combines features by their maximum rather than their sum;
- how its claim differs from attribution's and from the split test's;
- what it does not provide.

It SHALL cite the validation's verdict and measured rates.

#### Scenario: Reader asks what a hit means
- **WHEN** the reader consults the documentation about a feature with adjusted p ≤ 0.05
- **THEN** it states that the feature's covariate-adjusted group difference changes across stages, with the familywise error rate over all scanned features controlled at 0.05, and that the per-stage differences in the scan table show how
- **AND** it states that a hit does not say the feature drives a latent trajectory difference

#### Scenario: Reader asks what a non-rejection means
- **WHEN** the reader asks what a scan without hits says
- **THEN** the documentation states that it is uninformative unless power is shown at the design and effect, and cites the validation's plasmode power at e = 2.0 and e = 3.0

#### Scenario: Reader asks why the maximum
- **WHEN** the reader asks why the scan uses the maximum F rather than the split test's summed energy
- **THEN** the documentation states that a sum pays for every noise dimension while a maximum pays only for the search, so the scan targets differences on few features and the split test's complement targets differences spread over many, and it cites the validation's simulation localization

#### Scenario: Reader compares the scan with attribution
- **WHEN** the reader asks how the scan differs from orientation attribution
- **THEN** the documentation states that attribution ranks features by their contribution to a detected latent difference and has no error-rate claim, while the scan tests each feature's interaction directly under familywise control without using the latent space
- **AND** it states that the plasmode's attribution verdict is unchanged

### Requirement: The scan by block and its reading are documented

The measurement-visibility documentation SHALL describe the by-block scan: when to use it, what a by-block hit claims, and how detection moves between blocks. It SHALL cite the by-block validation's verdict and measured rates.

#### Scenario: Reader asks when to scan by block
- **WHEN** the reader asks whether to scan jointly or by block
- **THEN** the documentation states that the joint scan's familywise bar is set by the block that most often gives the largest null F, cites the measured null-maximum shares, and recommends the by-block scan when the blocks differ in size or noise

#### Scenario: Reader asks what a by-block hit means
- **WHEN** the reader consults the documentation about a feature with by-block adjusted p ≤ 0.05
- **THEN** it states that the familywise error over all scanned features is still controlled at 0.05, now split equally across the blocks, so each block has the same chance of a false hit whatever its size

#### Scenario: Reader asks what the by-block scan costs
- **WHEN** the reader asks what the by-block scan gives up
- **THEN** the documentation states that detection moves from the block with more features or heavier-tailed noise to the others, cites the validation's per-block recall in both arms, and states that only two blocks were measured
