# Spec Delta

## ADDED Requirements

### Requirement: The real-noise capture's diagnosis is documented

The measurement-visibility documentation SHALL report the diagnosis of the latent plasmode's `gate_failed` verdict on ADNI: where the bias lies, where the spread lives, whether ADNI's own interval for the capture is bounded, and what that means for stating a real-noise capture.

#### Scenario: Reader asks why no real-noise capture is stated
- **WHEN** the reader consults the latent-plasmode section
- **THEN** it states the parent's verdict, the diagnosis's verdict, and the statistic that decides it (ADNI's program energy against its cohort-level spread)

#### Scenario: Reader asks whether more units or a repaired gate would help
- **WHEN** the reader asks what a further run on ADNI could give
- **THEN** the documentation states what the diagnosis found about the spread's location and the cohort's interval, and what would be needed instead
