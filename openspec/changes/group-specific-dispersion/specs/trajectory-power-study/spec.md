# Spec Delta

## ADDED Requirements

### Requirement: Study configurations accept the group dispersion scale

The declarative study configuration SHALL accept `generator.group_dispersion_scale` in its baseline generator block and as an axis, wherever a generator setting is accepted. `null` and 1 SHALL mean absent. Null and negative-control cells SHALL keep the scale of the parameters they are derived from. A configuration that does not use the key SHALL load, enumerate, hash and dump exactly as before.

#### Scenario: A configuration sets the scale
- **WHEN** a configuration sets `generator.group_dispersion_scale` to 1.5
- **THEN** it loads and enumerates, every cell carries the scale, its anchor and negative-control cells included, and their cell identifiers and signatures differ from the same configuration without the key

#### Scenario: Null equals omitted
- **WHEN** a configuration sets `generator.group_dispersion_scale` to `null` or to 1
- **THEN** it enumerates identically to one that omits the key, and its dumped form omits the key

#### Scenario: The scale as an axis
- **WHEN** a design grid crosses `generator.group_dispersion_scale` over `null` and 2 with `generator.group_ratio` over 0.5 and 0.8
- **THEN** enumeration emits one design point per combination, and points that differ only in the scale have different generator parameters

#### Scenario: Committed configurations are unchanged
- **WHEN** every committed study configuration is enumerated after this change
- **THEN** its cell identifiers, signatures and matched seeds equal the pinned ones
