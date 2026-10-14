# Proposal

## Why

The net-displacement orientation estimator is the one measurement change that has moved orientation power at
the Phase 6 SEA-AD design: the [orientation measurement probe](../../../docs/reports/orientation-measurement-probe-2026-10-12.md)
saw `pls_cv × net_displacement` raise orientation/`angle` power from 0.29 to 0.67 at n = 80, leave `delta` and
`shape` untouched, and remove the Phase 5 shape → `angle` leak (0.89 → 0.02). That evidence is 100 replicates ×
199 permutations and was never a power claim. The [net-displacement Type I study](../../../docs/reports/net-displacement-type-i-2026-10-13.md)
closed the one objection the probe raised — the anchor `angle` excess was Monte Carlo noise; the lead arm holds
size at 0.050 at both designs against a 0.0695 bound — and its predeclared follow-up rule selected
`paper_grade_power_comparison`. Roadmap *Next changes* item 2 is that comparison: `pls_cv × net_displacement`
against `pls_cv × pc1` at 500 replicates × 999 permutations, over orientation, shape, magnitude and translation,
at the Phase 6 and Phase 5 designs. Until it runs, the estimator's power advantage is unestablished, the
production default stays PC1, and the case study must keep saying that an `angle` non-rejection is
uninformative at n = 80.

Two conditions were written down for this step and must be predeclared, not read off the result: the
translation control's `angle` rate is held to the anchor Type I bound (the lead arm passed it by one
rejection at Phase 6, 0.068 against 0.0695), and every Type I and specificity acceptance target uses the
one-sided `bound: "upper"` form, so the acceptance report agrees with the gate instead of flagging
conservative controls as failures.

## What Changes

- **A column decision rule in the study engine.** `acceptance.column_decision` names one candidate design
  point by its coordinates and judges it against the baseline column on the matched power grid: a target
  pair with a `gain` or `no_loss` criterion at the top effect (pooled-SE multiplier), an anchor Type I
  bound on every statistic, a list of **mandatory control** pairs (mode × statistic at the top effect) held
  to the same one-sided bound, and protected pairs that may not fall below the baseline by more than a
  pooled-SE loss. The report writes `column_decision.json` / `.csv` with a per-criterion verdict. Advisory,
  like `design_point`, `rank_decision` and `anchor_decision`; it feeds neither the gate nor the targets, and
  a config without it reports byte-identically.
- **Two committed comparison profiles**, `examples/trajectory_power_study/net_displacement_power_phase6.json`
  and `net_displacement_power_phase5.json`: one explicit-points config per design with two columns on
  matched datasets — `pls_cv × pc1` (baseline) and `pls_cv × net_displacement` (candidate) — over the four
  modes at the parents' effect axis, 500 × 999, every Type I and specificity target with `bound: "upper"`,
  the column rule declared with translation/`angle`, magnitude/`angle` and shape/`angle` as mandatory
  controls (Phase 6: `gain` on orientation/`angle`; Phase 5: `no_loss`), identity-pinned beside the other
  committed configs.
- **A cross-design selection script**, `scripts/net_displacement_power_follow_up.py`, reading both
  `column_decision.json` files and writing `comparison_decision.json`: `adopt_net_displacement` when the
  candidate qualifies at both designs, else `keep_pc1` with the failing design(s) and criteria. Advisory;
  flipping the production default is a separate change, as `default-joint-magnitude-construction` was.
- **The run, on `ing`,** with the committed provenance contract (timing job, budget in `PROVENANCE.txt`,
  BLAS pinned, 50 shards per design, `STUDY_N_JOBS` unset), under `results/net-displacement-power-<date>/`.
- **A dated findings report** `docs/reports/net-displacement-power-<date>.md`: per-mode `angle` curves under
  both estimators at both designs, the control table against the bound, the `angle` null-width covariate, the
  exact reproduction of the parent studies' records by the baseline column, the by-construction identity of
  `delta` and `shape` between the two `pls_cv` columns, the verdict, and — if the candidate qualifies — the
  revised case-study interpretability statement for `angle`. Roadmap item 2 closed with the branch taken;
  study README and CLAUDE.md updated.
- **`concat` ceiling and p-padding decided: not run here.** The Type I run's timing shows a `concat` unit
  costs 75 s at Phase 6 and 142 s at Phase 5 (RRPP in the 498-dimensional ambient space), so a full-grid
  ceiling column would add about 510 core-hours — more than the two PLS columns together — for a column the
  decision rule never reads. The report cites the probe's `concat × net_displacement` ceiling (orientation
  power 1.00 at both designs at 100 replicates) and the Type I study's `concat` null as the ceiling
  reference. A null-feature p-padding axis would be a generator change with its own truth and signatures;
  it is closed as not pursued unless `concat` ever becomes a candidate.

Not in scope: any change to the production `orientation_estimator` default, to the Phase 6 interpretability
statement before the report, a refit-per-permutation `pls_within_group`, or the forced-rank and `pca` columns
(the probe showed they lose Phase 6 shape power; the Type I study already measured their nulls).

## Capabilities

### New Capabilities

None.

### Modified Capabilities

- `trajectory-power-study`: the advisory column decision rule and its output; the two committed
  net-displacement power comparison profiles; the cross-design comparison selection; and the dated
  comparison findings report.

## Impact

- **Code:**
  - `src/motco/simulations/study/config.py`: `ColumnDecisionRule` (+ `MandatoryControl` pairs, candidate
    coordinates, criterion kind), parser, validation against `design_grid` and `trajectory_modes`,
    acceptance-block key.
  - `src/motco/simulations/study/targets.py`: `evaluate_column_decision` / `write_column_decision`,
    reusing the rank rule's `RateEstimate`, pooled-SE and anchor-bound helpers (the anchor bound is
    `phase4.type_i_inflation_bound`).
  - `scripts/motco_study.py`: write the decision when the rule is declared.
  - New `scripts/net_displacement_power_follow_up.py`; new configs under
    `examples/trajectory_power_study/`, added to `tests/data/committed_config_identity.json`.
- **Compatibility:** every existing config, cell id, signature, record and report output is byte-identical.
  The rule is a new optional acceptance key; the loader still rejects unknown keys.
- **Compute:** per design, 2 Type I baselines + 2 columns × (1 anchor + 4 modes × 4 effects) = 36 cells,
  18,000 `pls_cv` units. Planning figures from the Type I run's recorded per-unit medians under cluster
  contention: Phase 6 ≈ 18.5 s → about 95 core-hours; Phase 5 ≈ 74 s → about 370 core-hours, the dominant
  cost and the largest single study so far (the Phase 5 paper-grade run was about 150). The exact budget
  comes from a timing job before launch; the Phase 5 go/no-go is in the design.
- **Docs:** findings report, roadmap item 2 closed, study README (column rule, follow-up script, profiles),
  CLAUDE.md architecture note, `examples/trajectory_power_study/README.md`.
