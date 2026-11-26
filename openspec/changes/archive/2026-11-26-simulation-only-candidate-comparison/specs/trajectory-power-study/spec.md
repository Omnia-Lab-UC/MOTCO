# Spec Delta

## ADDED Requirements

### Requirement: Study configurations accept the stage dispersion scale

A study configuration SHALL accept `generator.stage_dispersion_scale` as a generator field and as an axis. It SHALL be omitted from cell ids, signatures and dumped configurations while absent, and null and Type I cells SHALL keep it.

#### Scenario: A configuration names the scale
- **WHEN** a configuration sets `generator.stage_dispersion_scale`
- **THEN** every enumerated cell, its null cells included, carries the scale, and its cell id and signature differ from those of the configuration without it

#### Scenario: Absent, nothing changes
- **WHEN** a committed configuration is enumerated
- **THEN** its cell ids, signatures and matched seeds equal their committed values

#### Scenario: A wrong length is refused at configuration time
- **WHEN** a configuration's scale has a length other than the number of stages
- **THEN** enumeration fails naming the field
