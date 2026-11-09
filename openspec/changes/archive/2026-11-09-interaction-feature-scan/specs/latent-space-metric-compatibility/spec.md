# Spec Delta

## ADDED Requirements

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
