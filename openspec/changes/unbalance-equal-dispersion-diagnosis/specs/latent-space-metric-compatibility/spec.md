# Spec Delta

## ADDED Requirements

### Requirement: The cell size at which the tests lose size at equal dispersion is documented

The measurement-visibility documentation SHALL state the unbalance diagnosis's verdict, the designs and statistics that failed at equal group dispersion, and the smallest cell size from which every reading held. It SHALL state the status of the group-dispersion validation's readings after the predeclared re-read, with the designs set aside and the scope of the re-read verdict.

#### Scenario: Reader asks whether a small group loses size
- **WHEN** the reader asks whether unequal group sizes or small cells alone can produce a rejection
- **THEN** the documentation gives the verdict, the rates by the smaller group's cell size and balance, and states that no balanced design below 7 samples per stage was read

#### Scenario: Reader asks whether the dispersion readings can be interpreted
- **WHEN** the reader asks what the group-dispersion readings mean after the diagnosis
- **THEN** the documentation gives the re-read verdict, the designs set aside or the reading replaced, the cell sizes the re-read verdict covers, and states that the re-read rule was fixed before the diagnosis ran and after the dispersion rates were known

#### Scenario: Reader asks whether a correction repairs it
- **WHEN** the reader asks whether the residual corrections or the within-group draws hold size in the small-cell designs
- **THEN** the documentation gives each recorded arm's label by family and states that the labels were recorded and not judged

#### Scenario: Reader asks how far the result carries
- **WHEN** the reader asks whether the result holds for a cohort's own table
- **THEN** the documentation states that it was read in simulation at δ = 2.0 with three stages and no covariates, and that a cohort's table is read at its own design
