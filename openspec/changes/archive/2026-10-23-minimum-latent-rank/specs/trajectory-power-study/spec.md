# Spec Delta

## ADDED Requirements

### Requirement: Committed configurations pin the unfloored rank rule

Every committed study configuration whose recorded results predate the rank floor SHALL name `1` as `evaluation.min_components` in its baseline evaluation block, and every committed diagnostic that constructs evaluation or analysis parameters directly and whose results predate the floor SHALL pin the floor of 1, so the default alters none of their cell identifiers, signatures, matched seeds or records. The pins record historical behavior and MUST be documented as such.

#### Scenario: Committed config identity is unchanged by the default floor

- **WHEN** each committed study configuration is loaded and enumerated after the default floor exists
- **THEN** its set of cell identifiers, every cell's parameter signature and every cell's matched seeds equal those recorded for that configuration before the change

#### Scenario: A committed config names the floor deliberately

- **WHEN** a committed configuration under `examples/trajectory_power_study/` is read
- **THEN** its evaluation block names `min_components`, and the pin is documented as historical

#### Scenario: A new config that omits the floor gets the default

- **WHEN** a configuration that names no `evaluation.min_components` is loaded
- **THEN** every PLS cell selects its rank from candidates at or above 2 and its records state the floor

### Requirement: Rank floor is a design axis

`evaluation.min_components` SHALL be accepted as a design-grid axis (in either form) and as a one-factor-at-a-time axis. In a configuration, `null` SHALL name the evaluator's default and `1` the unfloored rule, as distinct values. Design columns that differ only in the floor SHALL share generator parameters and matched seeds. The floor-1 column's records SHALL carry no floor entry and SHALL equal the records of the same configuration without the axis.

#### Scenario: Floor axis is applied per design point

- **WHEN** a configuration whose baseline evaluation names a floor of 1 declares the design axis `evaluation.min_components` with values `[1, 2]`
- **THEN** the floor-2 column's records carry `min_components: 2` in their runtime metadata and select ranks of at least 2 on the same generated datasets as the baseline column
- **AND** the baseline column's records carry no floor entry and equal the records of the same configuration without the axis

#### Scenario: Null names the evaluator's default

- **WHEN** a configuration whose baseline evaluation omits the floor declares the design axis `evaluation.min_components` with values `[null, 1]`
- **THEN** the `null` point is the baseline column and evaluates at a floor of 2, and the `1` column's records carry no floor entry
