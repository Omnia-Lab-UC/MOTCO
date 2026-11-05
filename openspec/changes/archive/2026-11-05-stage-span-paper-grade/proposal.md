# Proposal

## Why

The stage-span pilot (`stage-span-constructions`, 2026-11-04, `revise_construction`) showed that `angle` and `shape` are sensitive once a difference lies inside the stage-supervised span. That statement is at pilot grade (100 × 199). The paper needs it at paper grade (500 × 999, bound 0.0695), at more columns, and with a construction whose off-target leaks are understood. The pilot proposed capping the geometry effects at e ≤ 0.5 and explained the bend's `delta` leak by the sample scaler, without checking it.

A scratch decomposition run before this proposal (2026-11-04, 3,500 units on `ing`, job 933675) checked it:

- **The scaler route is real and deterministic.** The bend keeps path length exactly in the null-standardized coordinates it is defined in (B − A = 0). The plant changes the planted data's pooled per-feature SDs. Under those exact population SDs B is shorter by 0.8 % of path length at e = 0.5, with no sampling involved, and the sample scaler tracks that.
- **The e ≤ 0.5 cap is not enough at paper grade.** Shape e = 0.5 reads `delta` 0.08 at n = 200 and 0.11 with three blocks at R = 100, against a paper-grade bound of 0.0695.
- **Re-solving the bend under the planted data's own pooled SDs removes the scaler route.** The residual leak falls to ≤ 0.1 % at the balanced three-stage columns.
- **A second route remains: the PLS score metric.** The production statistics are measured in PLS score coordinates, `X·W·(PᵀW)⁻¹`, a sheared image of the fitted span. That map does not preserve feature-space path lengths. It leaves −0.27 % on the SEA-AD table and −0.58 % at four stages after the re-solve. Size is unaffected (equal groups stay equal under any linear map), but specificity is not. Nothing in the repository documents this yet.

This change makes the stage-span plant exact in the scale the measurement standardizes to, runs the paper grade with it, and documents the score-metric route.

## What Changes

- **Planted-scale stage-span plant.** A new opt-in generator setting `span_scale` with values `null` (the default, normalized to absent and hash-omitted) and `planted`. Under `planted`, the stage-span edit is solved by fixed-point iteration so that it is exact in coordinates standardized by the *planted* population's pooled SD: within-cell variance plus the between-cell variance of the planted cell means, which the measurement's sample scaler estimates. The plant then preserves `delta` (shape, orientation, translation), `angle` (shape, translation) and `shape` (orientation, translation) exactly in those coordinates. Truth records the scale, the iteration count and the residual. The setting applies to `stage_span` only, is cleared in zero-effect cells, and is refused together with `span_alignment`. Without it, every dataset, cell id, signature and record is byte-identical.
- **Study keys.** Configurations accept `generator.span_scale`. Like `construction`, it is rejected as a design axis, so scales are compared as sibling profiles.
- **Leak decomposition script.** `scripts/stage_span_leak_check.py` promotes the scratch breakdown. Per unit it reports the B − A path-length difference at six checkpoints, so each route can be read off one table:
  - null scale and exact planted scale;
  - the sample scaler;
  - the orthonormal fitted span and the PLS score coordinates;
  - the sample latent means.

  It also reports rejection rates, and it is sharded for `ing`. It is descriptive and issues no verdict.
- **Pre-flight gate.** Before the paper-grade array, the script runs at 100 × 199 on the paper-grade span geometry columns, under `planted` at e ∈ {0.1, 0.25, 0.5} and under `null` at the pilot's {0.25, 0.5, 1.0}. The paper grade is submitted only if:
  - the planted plant's scale route is exact in every unit (relative path-length gap ≤ 1e-8, fixed point converged);
  - no planted orientation or shape off-target rate exceeds the R = 100 bound of 0.0936.

  Otherwise the work stops and returns to the user.
- **Three paper-grade sibling profiles** (post-floor, `design-guidance` family, base seed 1028, 500 × 999, span capture on, `net_displacement`, rank floor 2, `joint`):
  - stage-span geometry under `planted`: orientation, shape and translation at e ∈ {0.1, 0.25, 0.5}, at five columns: n = 200 (baseline) and n = 80 at three stages on two blocks, the SEA-AD table `[[11, 10, 28], [9, 10, 12]]`, n = 200 on three blocks, and n = 200 at four stages;
  - stage-span magnitude under `planted` at e ∈ {0.02, 0.05, 0.1} at the same five columns;
  - surgery geometry at e ∈ {0.25, 0.5, 1.0} at the two columns no earlier paper-grade run covers (SEA-AD table, four stages). At the other three columns the surgeries are the design-guidance slice's cells (2026-10-29, 500 × 999) on the same datasets.
- **Predeclared decision.** `scripts/stage_span_paper_grade_decision.py` writes `stage_span_paper_grade_decision.json` (rule version 1, bound α + 2·√(α(1 − α)/R)). Criteria:
  - P1, size at every column: the anchor and the stage-span translation;
  - P2, median capture ≥ 0.90 in every stage-span magnitude, orientation and shape cell;
  - P3, off-target rates at or below the bound in those cells.

  The verdict is `sensitivity_established`, `specificity_excess`, `capture_shortfall` or `size_excess`; every failure is named. Descriptive outputs: target power ± SE, `informative_from` (the smallest effect with power ≥ 0.80) per column and statistic, capture quantiles, the intended against the `population_latent` geometry, and the visibility contrast (this run's surgery siblings plus the slice's surgery power and the pilot's surgery capture at the shared columns).
- **Reproduction checks.** Anchors across the three siblings, and the shared columns' anchors against the slice (identical p-values at 999 permutations). The surgery siblings against the stage-span pilot (SEA-AD) and the design-guidance pilot (four stages) at replicates 0–99, comparing seeds and statistics; p-values are not comparable across permutation counts.
- **Documentation.**
  - The dated report and aggregates under `results/stage-span-paper-grade-<date>/`.
  - A score-metric section on the visibility page: PLS scores are a sheared map of the span; what that does to specificity, not size; the measured residuals.
  - A new roadmap item for evaluating orthonormal span coordinates as the measurement, as its own change; roadmap item 13(a) and the "Not yet established" entry updated.
  - `AGENTS.md`. The version moves to MOTCO 0.17.0.

**Non-goals.**

- Changing the production measurement: the score coordinates stay, and orthonormal coordinates are a separate change.
- The latent plasmode (13(b)).
- A group-aware or out-of-span space (item 14's follow-up).
- Re-running the surgeries at the slice's columns.
- Effects above 0.5 for stage-span orientation and shape.

## Capabilities

### New Capabilities

None.

### Modified Capabilities

- `semisynthetic-trajectory-generator`: adds the opt-in planted-scale stage-span plant (`span_scale`), its truth and its hash omission.
- `trajectory-power-study`: configurations accept `generator.span_scale`, which is rejected as a design axis; stage-span anchors also clear `span_scale`.
- `stage-span-construction-validation`: adds the paper-grade profiles, the leak decomposition script and its pre-flight gate, the predeclared paper-grade decision, the reproduction checks and the dated report.
- `latent-space-metric-compatibility`: documents that the statistics are measured in PLS score coordinates, a non-isometric map of the span, and what that means for specificity.

## Impact

- **Code:**
  - `src/motco/simulations/semisynthetic.py`: the `span_scale` parameter and validation, a planted pooled-scale helper, the fixed point in `_stage_span_plant`, and truth;
  - `src/motco/simulations/grid.py`: `span_scale` joins `_SIGNATURE_OMIT_WHEN_NONE`;
  - `src/motco/simulations/study/` (`config` and `enumerate`): the key, the axis rejection, and clearing in zero-effect cells.
- **Scripts:** `scripts/stage_span_leak_check.py` (`run` | `merge`) with an sbatch template, and `scripts/stage_span_paper_grade_decision.py`, which reuses the pilot reading's helpers.
- **Configs:** three post-floor profiles under `examples/trajectory_power_study/`, appended to `POST_FLOOR_CONFIGS` and `tests/data/post_floor_config_identity.json`.
- **Tests:** generator fixed point and exactness, byte identity without the setting, study config and enumeration, profile structure and seed reproduction, the decomposition script and the decision.
- **Compute on `ing`:**
  - pre-flight: about 7,000 units at 199 permutations, about 45 core-hours;
  - paper grade: about 48,000 units at 999 permutations (span geometry 26,000, span magnitude 11,000, surgery geometry 11,000), about 290 core-hours from the scratch check's 18–28 s per unit. Plan for 350; the smoke run fixes the figure.
- **Compatibility:** every committed configuration, fixture and report is unchanged. The production measurement and `motco analyze` are unchanged.
- **Docs:** `docs/measurement-visibility.md`, `docs/roadmap.md`, `AGENTS.md`, `examples/trajectory_power_study/README.md`, and the dated report under `docs/reports/`.
