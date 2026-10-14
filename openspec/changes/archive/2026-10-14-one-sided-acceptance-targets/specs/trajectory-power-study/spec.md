# Spec Delta

## MODIFIED Requirements

### Requirement: Study evaluates results against pre-specified targets

MOTCO SHALL evaluate summarized results against the acceptance targets declared in the configuration and report, per target, whether it is met given Monte Carlo uncertainty.

A Type I control target or a specificity target MUST declare its bound form as one of two values, defaulting to the two-sided form when the key is absent. Under the **two-sided** form the target is met when the empirical rejection rate lies within `k` empirical Monte Carlo standard errors of alpha on either side. Under the **upper** form the target is met when the empirical rejection rate is at or below `alpha + k · sqrt(alpha · (1 − alpha) / R)`, with `R` the number of available replicates — the same bound the Phase 4 gate applies to its control and Type I inflation rules — so that a control rejecting less often than alpha is met. `k` is the target's standard-error tolerance in both forms. The report MUST record the bound form, the alpha-based standard error and the signed deviation from alpha only for a target declared with the upper form; a target in the default form MUST produce exactly the observations and rationale it produced before the bound form existed, so that every committed acceptance report regenerates byte-identically. Power monotonicity targets have no bound form.

#### Scenario: Type I control target is evaluated

- **WHEN** a null cell is evaluated against a Type I control target at alpha
- **THEN** the report indicates whether the empirical rejection rate is within the target's tolerance (expressed in Monte Carlo standard errors) of alpha

#### Scenario: Power monotonicity target is evaluated

- **WHEN** a power mode's diagonal statistic is evaluated against a monotonicity target
- **THEN** the report indicates whether the rejection rate is non-decreasing in effect size and reaches the target's minimum power at the largest effect size

#### Scenario: Specificity target is evaluated

- **WHEN** an off-diagonal mode × statistic combination is evaluated against a specificity target
- **THEN** the report indicates whether its rejection rate stays within tolerance of alpha

#### Scenario: Conservative control is met under the upper bound

- **WHEN** a specificity target declared with the upper bound at alpha `0.05` and tolerance `2` evaluates a control cell whose rejection rate is `0.000` over `500` available replicates
- **THEN** the target is met, the recorded bound is `0.0695`, the recorded deviation is `−0.050`, and the recorded bound form is upper

#### Scenario: Zero-rejection control is not met under the two-sided bound

- **WHEN** the same control cell is evaluated against a specificity target in the default two-sided form
- **THEN** the target is not met, because the empirical standard error of a zero rate is zero and the tolerance collapses to zero

#### Scenario: Inflated control is not met under the upper bound

- **WHEN** a Type I control target declared with the upper bound at alpha `0.05` and tolerance `2` evaluates a null cell whose rejection rate is `0.080` over `500` available replicates
- **THEN** the target is not met and the rationale states the rate exceeds the alpha-based bound `0.0695`

#### Scenario: Default-form targets regenerate byte-identically

- **WHEN** a study whose acceptance targets declare no bound form is reported before and after the bound form exists
- **THEN** its acceptance report CSV and JSON are byte-identical, with no added column or key

#### Scenario: Unknown bound form is rejected

- **WHEN** a study configuration declares a Type I control or specificity target with a bound form other than the two-sided or upper value
- **THEN** loading the configuration raises a clear validation error naming the target entry and the offending value
