# Spec Delta

## MODIFIED Requirements

### Requirement: Orientation estimator and integration method are design axes

The evaluation-namespace axes `evaluation.orientation_estimator` and `evaluation.integration_method` SHALL be accepted as design-grid axes (in either form) and as one-factor-at-a-time axes. The baseline value of `evaluation.orientation_estimator` is the baseline evaluation's estimator: the evaluator's default (net displacement) unless the baseline evaluation sets it. In a configuration, `null` SHALL resolve to the evaluator's default, as it does for the nested `evaluation.integration_params.<key>` form, and `"pc1"` SHALL name the PC1 estimator; `null` and `"pc1"` are therefore distinct axis values. Cells at a design point SHALL evaluate with the point's estimator and integration method, every RRPP draw and realized-geometry checkpoint included, and cells at the baseline column SHALL write records byte-identical to a configuration without these axes. The design-point operating table SHALL carry both axes as columns and SHALL accept columns whose integration method records no PLS component selection or selected rank, reporting those fields as unavailable rather than failing. The production forced-rank guard SHALL admit a forced-rank column declared through the explicit form on the same terms as one declared through the crossed form: its records pass when the recorded rank equals the rank at their design point.

#### Scenario: Estimator axis is applied per design point

- **WHEN** a configuration whose baseline evaluation names `"pc1"` declares the design axis `evaluation.orientation_estimator` with values `["pc1", "net_displacement"]`
- **THEN** the net-displacement column's records carry `orientation_estimator: "net_displacement"` in their runtime metadata and their `angle` p-values come from that estimator
- **AND** the baseline column's records carry no estimator entry and equal the records of the same configuration without the axis

#### Scenario: Null names the evaluator's default

- **WHEN** a configuration whose baseline evaluation omits the estimator declares the design axis `evaluation.orientation_estimator` with values `[null, "pc1"]`
- **THEN** the `null` point is the baseline column and evaluates with the net-displacement estimator, and the `"pc1"` column's records carry no estimator entry

#### Scenario: Non-PLS design columns are reported

- **WHEN** a merged result set contains design-point records whose integration method is `concat` or `pca`
- **THEN** the design-point operating table lists their rows with the rejection rate, Monte Carlo standard error, eigengap and `angle` null-width covariates, and leaves the component-selection and selected-rank fields unavailable
- **AND** reporting completes without error

#### Scenario: Forced-rank column declared as an explicit point passes the guard

- **WHEN** an explicit design point sets `evaluation.integration_params.forced_components` to 4 and its records record a forced selection at rank 4
- **THEN** the production component-selection guard admits those records

## ADDED Requirements

### Requirement: Committed configurations pin the PC1 orientation estimator

Every committed study configuration whose recorded results were produced under the PC1 estimator SHALL name `"pc1"` as `evaluation.orientation_estimator` in its baseline evaluation block, and every committed configuration that declares `evaluation.orientation_estimator` as a design axis SHALL name its PC1 column `"pc1"` rather than `null`, so that changing the evaluator's default alters none of its enumerated cells' parameter signatures, cell identifiers, matched seeds, or persisted records. Every committed diagnostic that constructs evaluation parameters directly, and whose committed results report `angle` under PC1, SHALL likewise select `"pc1"` explicitly. The pins record historical behavior and MUST be documented as such, so that new configurations choose an estimator deliberately rather than copy them.

#### Scenario: Committed config identity is unchanged by the default flip

- **WHEN** each committed study configuration is loaded and enumerated after the default changes
- **THEN** its set of cell identifiers, every cell's parameter signature, and every cell's matched seeds equal those recorded for that configuration before the change

#### Scenario: Committed design-grid baseline column stays PC1

- **WHEN** a committed net-displacement profile is loaded after the default changes
- **THEN** its baseline column is the `"pc1"` point, its phases and cell identifiers are unchanged, and its net-displacement columns are unchanged

#### Scenario: A new config that omits the estimator gets the default

- **WHEN** a configuration that names no `evaluation.orientation_estimator` is loaded
- **THEN** every cell evaluates with the net-displacement estimator and its records state that estimator
