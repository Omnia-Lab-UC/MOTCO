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

### Requirement: The measured dimension of the PLS measurement is documented

The documentation SHALL state that the production measurement uses the selected rank of the stage-PLS fit both to regularize the fit and as the number of dimensions the statistics are measured in. It SHALL state what measuring in the predictive subspace keeps and drops, report the pilot of the predictive-subspace coordinates against production, and state the production consequence.

#### Scenario: Reader asks how many dimensions the statistics use
- **WHEN** the reader consults the measurement-visibility documentation
- **THEN** it states that the production statistics are computed in as many dimensions as the cross-validation selects components, and cites the measured cost of dimensions beyond the stage program

#### Scenario: Reader asks what the predictive subspace is
- **WHEN** the reader asks how fewer dimensions can be measured on the same fit
- **THEN** the documentation states that the fitted stage labels move along at most the number of stages minus one directions of the score space, that measuring there keeps the fit and its rank, and that the harness and the real-data analysis offer three coordinate systems on that subspace as an opt-in setting

#### Scenario: Reader asks why a lower rank is not the same thing
- **WHEN** the reader asks whether forcing the fit to fewer components is equivalent
- **THEN** the documentation states that the later components aim the predictive directions, and cites the control arm's rates and capture beside the predictive arms'

#### Scenario: Reader asks what is given up
- **WHEN** the reader asks what the predictive subspace cannot see or report
- **THEN** the documentation states that a difference seen only through the dropped components is no longer read as it was, citing the surgery cell
- **AND** it states how the observed statistics compare with the population's through the same map, per stage effect

#### Scenario: Reader asks which measurement production uses
- **WHEN** the reader asks whether production and `motco analyze` measure in the predictive subspace
- **THEN** the documentation states the pilot's verdict and its size, power, specificity and gain readings, and that production measures in the fitted span's scores until a separate change adopts another measurement, which only a paper grade following `propose_paper_grade` permits

### Requirement: The predictive subspace's real-noise size excess is documented with its diagnosis

The measurement-visibility documentation SHALL report the diagnosis of the predictive-subspace coordinates' size excess on the SEA-AD plasmode anchor: what the within-stage draws and the out-of-fold coordinates each test, the verdict, and what it means for the opt-in predictive values.

#### Scenario: Reader asks why the predictive subspace fails size on real noise
- **WHEN** the reader consults the measured-dimension section
- **THEN** it states the diagnosis's verdict and each decisive arm's anchor rates beside the pilot's

#### Scenario: Reader asks what the two checks separate
- **WHEN** the reader asks how the cause was narrowed
- **THEN** the documentation states that within-stage draws change the null and keep the coordinates, that out-of-fold coordinates change the coordinates and keep the null, and what each outcome implies

#### Scenario: Reader asks whether the predictive values can be used
- **WHEN** the reader asks whether a case study may report the predictive values as tests
- **THEN** the documentation states that they remain opt-in diagnostics, that the diagnosis validates no test, and which separate change the verdict opens

### Requirement: The two permutation nulls and their reading are documented

The measurement-visibility documentation SHALL state what each permutation null of the trajectory test assumes and holds fixed, and SHALL report the within-stage pilot: the verdict, the size readings of both nulls where stage dispersion is equal and where it is unequal, the paired power and specificity readings on the production scores, and the re-read of `predictive_orthonormal` under the same draws.

#### Scenario: Reader asks what the two nulls assume

- **WHEN** the reader consults the section on the permutation null
- **THEN** it states that the unrestricted draws assume the reduced model's residuals are exchangeable across all samples, that the within-stage draws assume it within each stage only, and that the within-stage draws keep every stage's pooled mean at its fitted value

#### Scenario: Reader asks when the choice matters

- **WHEN** the reader asks in which data the two nulls differ
- **THEN** the documentation reports both nulls' rates on the dispersion cells, with the direction of the unrestricted null's error when the largest stage is the narrowest and when it is the widest

#### Scenario: Reader asks which null to use

- **WHEN** the reader asks whether a case study may report within-stage p-values
- **THEN** the documentation states the pilot's verdict and its consequence, that the unrestricted draws remain the default, and that the split interaction test and the scan are drawn unrestricted

#### Scenario: Reader asks whether the predictive coordinates are reopened

- **WHEN** the reader asks what the pilot means for the predictive values
- **THEN** the documentation states the candidate's reading against production under the same draws, and whether a paper grade is proposed

### Requirement: Real-noise size on a second cohort is documented with its scope

The measurement-visibility documentation SHALL report the size readings on ADNI's blocks under within-stage-permuted labels: the design points, the verdict on the production scores and on `predictive_orthonormal` under each null, the fixed-rank reading against the in-sample stage fit, and the scope a real-noise size reading has.

#### Scenario: Reader asks whether the predictive subspace's excess is specific to SEA-AD
- **WHEN** the reader consults the measured-dimension section
- **THEN** it states the ADNI anchor's verdict, `predictive_orthonormal`'s size scope beside production's, and the SEA-AD anchor's counts beside each ADNI design point's

#### Scenario: Reader asks what produces the excess
- **WHEN** the reader asks whether an over-fitted stage fit is enough to produce it
- **THEN** the documentation states the fixed-rank grade, the readings at ranks 8 and 19 against the in-sample stage R² at each design point, and whether the within-stage draws remove the excess there

#### Scenario: Reader asks whether production's size has a second real-noise reading
- **WHEN** the reader asks on which real data the production test's size was read
- **THEN** the documentation names both exports, states production's size scope on ADNI's design points, and states that the full-cohort rates are conditional on one fit

#### Scenario: Reader asks what a size reading on one cohort licenses
- **WHEN** the reader asks whether a hold on ADNI admits the predictive values elsewhere
- **THEN** the documentation states that a size reading speaks for its export and design points, that it does not lift the SEA-AD reading, and that reopening an option at a design also needs a gain and specificity there

### Requirement: Real-noise power and specificity at a cohort's design are documented with their references

The measurement-visibility documentation SHALL report the plant readings on ADNI's blocks under within-stage-permuted labels and the stage-span readings in the simulation of ADNI's design: the design points and cells, what an in-sample and a held-out plant can each show, the verdict on `predictive_orthonormal` with the criterion that decided it, production's own power and off-target statement, and what the verdict licenses.

#### Scenario: Reader asks whether the predictive subspace is reopened at ADNI's design
- **WHEN** the reader consults the measured-dimension section
- **THEN** it states the verdict, the criterion that decided it with its readings, and that `predictive_orthonormal` remains an opt-in diagnostic whose default can change only through a later paper-grade run

#### Scenario: Reader asks what the production test detects on ADNI's noise
- **WHEN** the reader asks from which effect a stretch, a rotation or a bend of the stage program is detected on real noise
- **THEN** the documentation states production's target rates on the held-out and the in-sample plants at each design point, and the effect from which each statistic reaches the power floor, or that none does

#### Scenario: Reader asks whether a rejection on real noise names the right kind of difference
- **WHEN** the reader asks whether a planted stretch, rotation or bend is rejected by the other statistics
- **THEN** the documentation states production's and `predictive_orthonormal`'s off-target rates on the in-sample plants, which family of readings the decision judged, and that no off-target rate of a held-out plant was judged

#### Scenario: Reader asks why a plant has two references
- **WHEN** the reader asks why the plants were defined both on the analysed sample and on held-out participants
- **THEN** the documentation states that an in-sample plant has an exact truth under the label permutation and lies in the fitted span by construction, that a held-out plant is independent of the fit and is mostly estimation noise of the reference, and the shares of each that the fitted span keeps

#### Scenario: Reader asks what happens to a difference outside the stage span
- **WHEN** the reader asks how a rotation out of the stage plane is read on ADNI's noise
- **THEN** the documentation states both arms' rates on the out-of-span cells for the three statistics, that these cells were recorded and not judged, and that the split test's complement is the designated detector there

#### Scenario: Reader asks what the verdict licenses
- **WHEN** the reader asks whether a reading at ADNI's design admits the predictive values elsewhere
- **THEN** the documentation states that the reading speaks for ADNI's export, the two design points and the simulated table, that the SEA-AD reading stands for SEA-AD's design, and that no real-noise reading of a held-out plant exists at the full cohort

### Requirement: The real-noise capture's diagnosis is documented

The measurement-visibility documentation SHALL report the diagnosis of the latent plasmode's `gate_failed` verdict on ADNI: where the bias lies, where the spread lives, whether ADNI's own interval for the capture is bounded, and what that means for stating a real-noise capture.

#### Scenario: Reader asks why no real-noise capture is stated
- **WHEN** the reader consults the latent-plasmode section
- **THEN** it states the parent's verdict, the diagnosis's verdict, and the statistic that decides it (ADNI's program energy against its cohort-level spread)

#### Scenario: Reader asks whether more units or a repaired gate would help
- **WHEN** the reader asks what a further run on ADNI could give
- **THEN** the documentation states what the diagnosis found about the spread's location and the cohort's interval, and what would be needed instead

### Requirement: The comparison with the default per-feature analysis is documented

The measurement-visibility documentation SHALL describe how the production trajectory test compares with a group × stage interaction test on each feature with no integration: what each side tests, which differences each detects, what only the trajectory test can claim, and the limits of the comparison. It SHALL cite the validation's verdict, its characterization grade and its measured rates.

#### Scenario: Reader asks whether integration is needed to detect a difference
- **WHEN** the reader asks whether the trajectory test detects more than a per-feature test
- **THEN** the documentation states the validation's verdict and, by plant family and design, which side leads, with the rates
- **AND** it states that production entered as one test at 0.05 / 3 and the baseline as one test at 0.05

#### Scenario: Reader asks what the per-feature test cannot say
- **WHEN** the reader asks what the trajectory test adds beyond detection
- **THEN** the documentation states that a per-feature test has no geometry and cannot say whether a difference is one of size, orientation or shape, and it states the validation's characterization grade for production

#### Scenario: Reader asks which baseline was used
- **WHEN** the reader asks how the per-feature tests were combined
- **THEN** the documentation names the maximum, the false-discovery reading and the pooled features, states which held size, and states that a moderated statistic was not an arm

#### Scenario: Reader asks how far the result carries
- **WHEN** the reader asks whether the result holds for a real cohort
- **THEN** the documentation states that it was read in simulation at 658 features with Gaussian noise, that real blocks hold far more features, and that no real-noise reading exists on methylation, expression and proteomics
