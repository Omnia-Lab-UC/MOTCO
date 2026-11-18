# Spec Delta

## ADDED Requirements

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
