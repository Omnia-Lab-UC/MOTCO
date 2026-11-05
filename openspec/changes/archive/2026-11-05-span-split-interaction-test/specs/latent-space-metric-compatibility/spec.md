# Spec Delta

## ADDED Requirements

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
