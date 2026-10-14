# Proposal

## Why

The orientation measurement probe (`docs/reports/orientation-measurement-probe-2026-10-12.md`) found that the
net-displacement orientation estimator is the main lever on orientation power at the Phase 6 SEA-AD design: in
the production stage-PLS space it raises orientation/`angle` power from 0.29 to 0.67 at n = 80, leaves `delta`
and `shape` unchanged at both designs, and removes the shape → `angle` leak at Phase 5 (0.89 → 0.02). It did
not advance, for one reason: its anchor `angle` rate was 0.100 at Phase 6 against a bound of 0.0936, and
0.08–0.10 in every space at both designs, where PC1 sits at 0.01–0.05.

That evidence is suggestive, not established. The probe ran 100 replicates × 199 permutations, the four
spaces share datasets and permutation draws, so the excess is effectively one observation of about 9/100 per
design, and the lead arm missed the bound by less than one rejection. The roadmap's item 3 and the probe's §7
predeclare the next step: an anchor-only Type I study at paper-grade size (500 replicates × 999 permutations)
at both designs, for every net-displacement arm, *before* any power comparison. Its result selects the
follow-up: if the excess holds, diagnose it (endpoint-stage heteroscedasticity under full-residual permutation
is the first hypothesis); if size holds, compare `pls_cv × net_displacement` against `pls_cv × pc1` at paper
grade.

The probe driver could run this, but the paper-grade comparison that may follow needs the estimator inside the
study engine (config identity pins, sharding, matched-seed families, design-point reporting, acceptance
rules). Running the Type I study through the engine does that wiring once and gives the Type I result the
same provenance contract as every committed study.

## What Changes

- **Design grid accepts an explicit list of design points.** Today `design_grid.axes` crosses every axis
  value into a rectangle. The five measurement columns this study needs — `pls_cv × pc1` (baseline),
  `pls_cv × net_displacement`, `pls_forced × net_displacement`, `pca × net_displacement`,
  `concat × net_displacement` — are not a rectangle: crossing `integration_method`, `forced_components`,
  `n_components` and `orientation_estimator` would enumerate invalid cells (`pca` without a rank, `concat`
  with a forced rank). A new `design_grid.points` form lists the design points directly. Every point names
  every axis, the baseline point must be listed, the same `generator.` / `evaluation.` /
  `evaluation.integration_params.<key>` namespace rules apply, and `points` and `axes` are mutually
  exclusive. A config using `axes` enumerates byte-identically to today.
- **Evaluation-namespace axes for the estimator and the integration method are exercised and pinned.**
  `evaluation.orientation_estimator` and `evaluation.integration_method` already pass the generic axis
  mechanism (verified: both enumerate as design axes on the committed Phase 6 config). This change adds the
  tests and the spec scenarios that make that a contract, and makes the design-point readers tolerate columns
  whose integration method records no PLS component selection (`concat`, `pca`).
- **Advisory anchor Type I decision rule.** A new `acceptance.anchor_decision` block (alpha, SE multiplier,
  statistics) makes the report evaluate every design point's zero-effect anchor against the probe's bound
  α + k·√(α(1 − α)/R), with R the point's available replicate count (0.0695 at R = 500, k = 2), and write
  `anchor_decision.json` with per-point, per-statistic pass/fail. Advisory only; it feeds neither the gate nor
  the acceptance-target report, and configs without the block are byte-identical.
- **Two committed null-only profiles**, `examples/trajectory_power_study/net_displacement_type_i_phase6.json`
  and `net_displacement_type_i_phase5.json`. Both enumerate only negative controls — the shared zero-effect
  anchor and the `translation` control at e = 1.0 — at 500 × 999 across the five columns above. The Phase 6
  profile reuses the Phase 6 study's seed family and base seed, so its baseline anchor evaluates the same 500
  datasets as `results/phase6-small-n-2026-10-09/`; the Phase 5 profile reuses the Phase 5 primary family under
  the joint magnitude construction, which is immaterial at zero effect. Attribution and the Phase 4 gate are
  off.
- **Follow-up selection.** A small script reads both profiles' `anchor_decision.json` and writes
  `follow_up_decision.json`: `paper_grade_power_comparison` when the lead arm `pls_cv × net_displacement`
  holds size at both designs, `diagnose_anchor_excess` otherwise, with a descriptive localization tag
  (excess in every net-displacement column, or only in some).
- **Run on `ing`**, commit `results/net-displacement-type-i-<date>/` (both designs), write the findings report
  `docs/reports/net-displacement-type-i-<date>.md`, and update the roadmap, the study README and CLAUDE.md.

Out of scope, deliberately (each is its own conditional change): the heteroscedasticity diagnosis, the
paper-grade power comparison, the `concat` p-padding question (null features do not move the null, so it is a
power-study question), and a refit-per-permutation `pls_within_group`.

**Assumption to confirm at review:** the fixed-basis `pls_within_group` arm is excluded. The probe already
identified why its anchor inflates (the basis carries the observed groups' noise and the null never sees the
fitting step; anchor `angle` 0.48–0.88), and the roadmap recommends dropping it; re-running it at 500 × 999
would spend compute confirming a known mechanism. Add it back as a sixth point if the paper needs the number.

## Capabilities

### New Capabilities

None.

### Modified Capabilities

- `trajectory-power-study`: the explicit-points design grid; the estimator and integration-method evaluation
  axes as contracted design axes, with design-point reporting over non-PLS columns; the advisory anchor Type I
  decision rule and its output; the two committed net-displacement Type I profiles; the follow-up decision
  output; and the dated Type I findings report.

## Impact

- **Code:**
  - `src/motco/simulations/study/config.py`: `DesignGrid.points`, its validation, `design_points()` over the
    list form, `AnchorDecisionRule` and its parser, acceptance-block key.
  - `src/motco/simulations/study/enumerate.py`: no change expected; `_design_point_cells` already iterates
    `config.design_points()`.
  - `src/motco/simulations/study/spectrum.py` / `report.py` / `targets.py`: tolerate design columns without
    `component_selection` or `selected_lv`; evaluate and write the anchor decision.
  - `scripts/motco_study.py`: write `anchor_decision.json` when the rule is declared.
  - New `scripts/net_displacement_follow_up.py`.
  - New configs under `examples/trajectory_power_study/`, added to
    `tests/data/committed_config_identity.json`.
- **Compatibility:** every existing config, cell id, signature, record and report output is byte-identical.
  `tests/test_committed_config_identity.py` and the existing design-grid, rank-ladder and report tests must
  stay green. The production forced-rank guard keeps admitting forced records only at a design point whose
  declared rank equals the recorded rank; the `pls_forced` column satisfies that through the points form.
- **Compute:** 12 cells × 500 units per design (two Type I baselines, the primary anchor and translation
  cell, and four net-displacement columns × two cells). Planning figures from the committed runs: about
  20 core-hours at Phase 6 (≈ 13 s per `pls_cv` unit, less for the other spaces) and 100–130 core-hours at
  Phase 5 (≈ 90 s per `pls_cv` unit at n = 1200). The exact budget is written down from a timing job before
  launch, as for every committed study.
- **Docs:** findings report, roadmap item 3 closed with the branch taken, study README (points form, anchor
  rule, follow-up script), CLAUDE.md architecture note.
