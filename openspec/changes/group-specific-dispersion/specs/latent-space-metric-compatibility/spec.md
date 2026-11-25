# Spec Delta

## ADDED Requirements

### Requirement: What the tests assume about group dispersion is documented

The measurement-visibility documentation SHALL state that the permutation null exchanges residuals between the groups and so assumes that they share one noise dispersion. It SHALL say what a difference in dispersion alone does to `delta`, `angle` and `shape`, and SHALL cite the group-dispersion validation's verdict, its map by statistic and its measured rates.

#### Scenario: Reader asks whether a rejection can come from dispersion alone
- **WHEN** the reader asks whether a rejection can arise when the groups' trajectories agree and their dispersions differ
- **THEN** the documentation states the validation's verdict and, for each statistic, the smallest scale that failed in balanced designs, where the wider group is the smaller and where it is the larger, with the rates

#### Scenario: Reader places a cohort
- **WHEN** the reader has a cohort's group × stage cell table and the two groups' residual dispersion in the latent scores
- **THEN** the documentation says how the ratio of those dispersions compares with the simulated scale, and gives the readings by scale, group share, sample size and stage effect

#### Scenario: Reader asks about the per-feature tests
- **WHEN** the reader asks whether the interaction scan, its false-discovery reading or the split test hold size under unequal group dispersion
- **THEN** the documentation gives each recorded arm's label by family and states that the labels were recorded and not judged

#### Scenario: Reader asks what would repair it
- **WHEN** the reader asks whether another null holds size
- **THEN** the documentation gives the recorded reading of the within-group draws and states that they are not a package option of the harness or of `motco analyze`

#### Scenario: Reader asks how far the result carries
- **WHEN** the reader asks whether the result holds for a real cohort
- **THEN** the documentation states that it was read in simulation with Gaussian noise, one scale applied to every feature of one group, three stages and no covariates, that a difference in correlation structure was not simulated, and that no real cohort's dispersion was read
