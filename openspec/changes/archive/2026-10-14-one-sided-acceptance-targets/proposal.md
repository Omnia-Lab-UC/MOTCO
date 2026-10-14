# Proposal

## Why

The acceptance-target reader (`study/targets.py`) judges every Type I and specificity target two-sided:
`|rate − α| ≤ k · SE`, with the SE computed from the *empirical* rate. The Phase 4 gate judges the same
controls one-sided: `rate ≤ α + k · √(α(1 − α)/R)`. The two predeclared readers therefore disagree whenever a
control is conservative, and the disagreement is already in the record. The Phase 5 magnitude re-measurement
([addendum](../../../docs/reports/phase5-magnitude-remeasurement-2026-10-08.md), §3) has the gate meeting both
magnitude controls while the target reader records them as **not met**: magnitude/`angle` 0.024 is 0.026 below
α against a 0.0137 tolerance, and magnitude/`shape` 0.000 has an empirical SE of 0, so its tolerance is 0 and
any rate other than exactly α fails. The same reader flags the `none` Type I baseline's `delta` (0.032), and
the committed net-displacement Phase 6 report flags its `type_i_control` `shape` row (0.018). None of these is
inflation. A control that rejects *less* often than α is the property the target exists to confirm, and a
reader that calls it a failure has to be explained away in every findings report.

The roadmap's item 2 predeclares the fix: make the target's bound one-sided, mirroring the gate's α-based SE,
as an opt-in key that defaults to the current two-sided form so every committed `report/` regenerates
byte-identically. The paper-grade `pls_cv × net_displacement` vs `pls_cv × pc1` comparison (roadmap item 3)
will predeclare translation/`angle` as a mandatory control held to the anchor bound; its acceptance block
should be able to say the same thing the gate says, so this change lands first.

## What Changes

- **Opt-in one-sided bound on the Type I control and specificity targets.** Each entry under
  `acceptance.type_i` and `acceptance.specificity` accepts a `bound` key with values `two_sided` (default,
  today's rule: `|rate − α| ≤ k · SE_empirical`) or `upper` (`rate ≤ α + k · √(α(1 − α)/R)`, R the
  available replicates — the gate's `type_i_inflation_bound`, sharing its implementation). `k` stays
  `se_tolerance`. Any other value is rejected at load time naming the entry. Power monotonicity targets are
  untouched (they are already one-sided floors).
- **Observations record the bound kind only when it is one-sided.** Under `upper` the evaluation's
  `observations` gain `bound_kind: "upper"`, `alpha_se` (the α-based SE) and a *signed* `deviation`
  (`rate − α`), with `bound` holding the upper limit and the rationale reading `rate ≤ alpha + k·SE_alpha`.
  Under the default nothing is added, no key is renamed and the rationale text is unchanged, so
  `acceptance_report.csv` / `.json` for every committed config regenerate byte-for-byte (the CSV's `obs_*`
  columns are the union of observation keys, so a report with no `upper` target gains no column).
- **Config dump records the resolved value.** `dump_study_config` writes `bound` on every target (resolved
  dataclass form, round-trips through load); no dump is committed, and acceptance targets never enter a
  cell's parameter signature, so cell ids, signatures and matched seeds are unaffected.
- **No committed config opts in.** The ten-plus committed profiles keep their two-sided targets so their
  committed reports stay reproducible as recorded; the net-displacement power comparison (roadmap item 3) is
  the first adopter. The Phase 4 gate, the anchor decision rule and the rank/design-point rules are unchanged
  (they are already one-sided where it matters).
- **Documentation.** `study/README.md` documents the key and the two formulas; `docs/roadmap.md` closes item
  2 and notes the comparison config should set `bound: "upper"` on its controls.

## Capabilities

### New Capabilities

None.

### Modified Capabilities

- `trajectory-power-study`: the requirement *Study evaluates results against pre-specified targets* gains the
  declared bound form — a Type I control or specificity target is judged two-sided against the empirical
  Monte Carlo SE by default, or, when the target declares an upper bound, one-sided against the α-based SE,
  with the report recording the bound kind only for one-sided targets, default-form reports regenerating
  byte-identically, and an unknown bound value rejected at load time.

## Impact

- `src/motco/simulations/study/config.py` — `TypeIControlTarget` / `SpecificityTarget` gain `bound`;
  `_build_type_i_target` / `_build_specificity_target` parse and validate it.
- `src/motco/simulations/study/targets.py` — `_evaluate_type_i` / `_evaluate_specificity` branch on the
  bound kind; the `upper` branch reuses `phase4.type_i_inflation_bound` (targets → phase4 is a new import
  edge; phase4 does not import targets, so no cycle).
- `tests/test_study_targets.py`, `tests/test_study_config.py` — new cases for both kinds, both targets, the
  zero-rejection control (`shape` 0.000 at R = 500 is met under `upper`, not met under `two_sided`), the
  inflated control (0.080 at R = 500 fails under `upper`), the rejection of an unknown value, and a
  byte-identity check that a two-sided-only study's `acceptance_report.*` is unchanged.
- `src/motco/simulations/study/README.md`, `docs/roadmap.md`, `openspec/specs/trajectory-power-study/spec.md`.
- No change to any committed config, result, cell id, signature, seed, gate output or figure.
