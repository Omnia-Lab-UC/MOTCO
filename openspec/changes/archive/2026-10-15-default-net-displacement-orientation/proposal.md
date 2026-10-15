# default-net-displacement-orientation

## Why

The [net-displacement power comparison](../../../docs/reports/net-displacement-power-2026-10-14.md) returned
the verdict `adopt_net_displacement` (`results/net-displacement-power-2026-10-14/comparison_decision.json`):
on identical datasets and permutation draws, the net-displacement orientation estimator raises
orientation/`angle` power from 0.300 to 0.648 at the Phase 6 small-n design and from 0.840 to 0.924 at
Phase 5, holds anchor size (0.050 at both designs, bound 0.0695), and removes the shape → `angle` leak that
PC1 shows (0.138 / 0.796 → 0.054 / 0.022). The [Type I study](../../../docs/reports/net-displacement-type-i-2026-10-13.md)
had already shown size holds at paper grade. The estimator is still opt-in: every caller that relies on
`SimulationEvaluationParams`'s default, including a study config that omits the key, measures `angle` with
PC1, the estimator the paper will not describe. The roadmap lists this flip as "Next changes" item 2, and
the report's §7 names it as the follow-up that the verdict selects.

## What Changes

- **Flip the harness default.** `SimulationEvaluationParams.orientation_estimator` defaults to
  `"net_displacement"`. `"pc1"` stays selectable and keeps normalizing to `None`, the value that
  `grid._SIGNATURE_OMIT_WHEN_NONE` leaves out of every hash payload, so a PC1 evaluation's cell id, parameter
  signature and persisted record are byte-identical to what they were before the estimator existed. The
  observed `angle`, every RRPP draw and every realized-geometry checkpoint `angle` follow the default.
  Default records now state `orientation_estimator: "net_displacement"` in their runtime metadata; PC1
  records keep stating nothing, as before.
- **Make a config `null` mean the evaluator's default, not PC1.** Today the study loader treats `null` and
  `"pc1"` as one value for `evaluation.orientation_estimator`, which is correct only while the default *is*
  PC1. After the flip `null` resolves to the evaluator's default (net displacement), consistent with the
  nested `evaluation.integration_params.<key>` form, and `"pc1"` names PC1. A design-grid baseline point
  equals the baseline evaluation's estimator, whichever it is.
- **Pin `"pc1"` into every committed study config.** All 19 configs under `examples/trajectory_power_study/`
  gain `"orientation_estimator": "pc1"` in their `evaluation` block, and the four net-displacement profiles
  additionally rewrite their baseline column's `"evaluation.orientation_estimator": null` to `"pc1"`. Because
  `"pc1"` normalizes to the hash-omitted `None`, every enumerated cell's signature, cell id and matched seed is
  preserved. Committed `results/` stay resumable and reproducible from their configs. A regression test
  checks this against a fixture of pre-flip signatures and seeds captured on the unedited tree.
- **Pin `"pc1"` at the historical diagnostics' call sites** that construct `SimulationEvaluationParams`
  directly and whose committed results report `angle`: `realized_geometry_study.py`, the three sites in
  `specificity.py`, and `scripts/latent_rank_probe.py`. A re-run then regenerates what they measured. The two
  byte-identity fixture capture modules under `tests/` pin `"pc1"` for the same reason. The showcase follows
  the default.
- **Leave the library primitives at PC1.** `estimate_difference`, `RRPP`, `_RRPPWorker` and
  `diagnostics.geometry_from_means` keep `orientation_estimator="pc1"`. The library's `angle` contract
  (`trajectory-orientation-invariance`) and the R-reference regression tests are unchanged; the validated
  claim is about the harness pipeline (stage-PLS latent space, RRPP), which is where the default flips.
  Real-data code that calls the primitives directly must pass `orientation_estimator="net_displacement"` to
  measure what the paper describes; the roadmap's case-study item says so.
- **Revise the case-study `angle` statement and the roadmap.** The Phase 6 interpretability row for `angle`
  becomes the one the power comparison report gives in §7 (a rejection is evidence that the groups' net
  directions of change differ; a non-rejection stays uninformative at n = 80; the cohort's selected rank
  replaces its eigengap as the reported covariate). "Next changes" retires item 2, and item 1 (the case
  study) is recorded as waiting for this flip, which costs no compute.
- **Update the docs to the post-flip state**: `CLAUDE.md`, `docs/api/simulations.md`, the study README,
  and the examples README (why every config pins `"pc1"`, that the pin is behavior-neutral and which test
  proves it, and that new configs should pick an estimator deliberately rather than copy the pin).

Out of scope:

- changing the default of the stats primitives or adding an estimator flag to `motco de`;
- any re-run, re-issued report, gate, RRPP, PLS or statistic change; committed `results/` and
  `docs/reports/` are frozen;
- the null-dispersion materiality default (roadmap item 3).

No compute is needed.

## Capabilities

### New Capabilities

(none)

### Modified Capabilities

- `simulation-evaluation-harness`: "Harness accepts an orientation estimator selection" (which made PC1 the
  default and the default invisible) is replaced by "Harness measures orientation with net displacement by
  default" — PC1 stays selectable and is the estimator whose selection is omitted from records, signatures
  and cell ids, so PC1 results equal those produced before the selection existed; default records state
  their estimator.
- `trajectory-power-study`: "Orientation estimator and integration method are design axes" — a config
  `null` resolves to the evaluator's default and `"pc1"` names PC1; the baseline column's estimator is the
  baseline evaluation's. New requirement "Committed configurations pin the PC1 orientation estimator" —
  every committed config names `"pc1"`, design-grid baseline points name it explicitly, and enumerated
  identity is unchanged against a pre-flip fixture.
- `trajectory-orientation-invariance`: "Net-displacement orientation estimator is an opt-in alternative" —
  the PC1 default is scoped to trajectory estimation and RRPP (the library primitives); the simulation
  harness selects its own default, which this change sets to net displacement.

## Impact

- `src/motco/simulations/evaluation.py`: dataclass default and comment. `src/motco/simulations/grid.py`:
  the `_SIGNATURE_OMIT_WHEN_NONE` comment. `src/motco/simulations/study/config.py`: the axis-value
  normalizer for `evaluation.orientation_estimator`.
- `src/motco/simulations/realized_geometry_study.py`, `src/motco/simulations/specificity.py`,
  `scripts/latent_rank_probe.py`: explicit `orientation_estimator="pc1"`.
- `examples/trajectory_power_study/*.json` (19 files) and `examples/trajectory_power_study/README.md`.
- `tests/`: a new pre-flip identity fixture and test; `tests/orientation_probe_fixture_capture.py` and
  `tests/phase6_fixture_capture.py` pin `"pc1"`; tests that assert the old default (expected in
  `test_diagnostic_integration.py`, `test_study_measurement_axes.py`, `test_study_design_points_explicit.py`,
  `test_net_displacement_power_checks.py`, and any test config fixture that omits the key) are pinned or
  updated without weakening an assertion.
- `docs/roadmap.md`, `CLAUDE.md`, `docs/api/simulations.md`, `src/motco/simulations/study/README.md`.
- `src/motco/stats/` is untouched. No effect on committed `results/`, `docs/reports/`, or any committed
  `report/` directory. No cluster time.
