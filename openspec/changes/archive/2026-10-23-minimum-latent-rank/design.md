# Design

## Context

See proposal.md (Why) for the motivation. The pieces this change touches:

- `plsda_doubleCV` (`stats/pls.py`) builds `n_lv_candidates = range(1, max_components)`, averages the inner-fold AUROC per candidate, takes `np.argmax(...) + 1` per outer fold, and reduces the outer folds to a modal rank with the parsimony tie-break. `max_components` is validated against the feature count only.
- `fit_pls_latent_space` (`simulations/evaluation.py`) is the one implementation of the production latent space: it clamps `max_components` to `[2, min(n_features, n_samples)]`, calls `plsda_doubleCV`, takes the modal rank across repeats, fits the pooled model once, and records `selected_lv`, `component_selection` and the CV knobs in `integration_params`. `_forced_components` validates the diagnostic override against `[2, max]` and skips CV. Both the harness (`_pls_integration`) and `analyze_trajectories` call it.
- `SimulationEvaluationParams` holds `integration_params` as a mapping that is hashed key by key into cell ids and signatures. `orientation_estimator` is a top-level field whose PC1 value is held as `None`, omitted from hashes and records through `grid._SIGNATURE_OMIT_WHEN_NONE`; the study config loader (`_build_evaluation`, `_normalize_orientation_estimator`, `_AXIS_VALUE_NORMALIZERS`, `_display_axis_value`) reads `null` as the evaluator's default and `"pc1"` as the omitted value.
- `TrajectoryAnalysisParams` serializes with `to_dict()`, dropping `covariate_latent_adjustment` while `None` so 0.8.0 records are unchanged. `motco analyze` builds it from flags.
- The identity chain: `tests/data/plsda_doublecv_identity.json` (the primitive; one repeat selected rank 1), `pls_latent_space_identity.json` (the core; selected rank 1), `phase6_prechange_fixture.json` and `orientation_probe_prechange_fixture.json` (harness defaults, rank 2), `analysis_no_covariate_identity.json` and `analysis_no_latent_adjustment_identity.json` (0.7.0 and 0.8.0 analyses, rank 2), `committed_config_identity.json` and `orientation_default_identity.json` (cell ids, signatures, seeds of every committed config). Every committed config pins `"pc1"`; `tests/test_orientation_default_pins.py` proves the pins are behavior-neutral.
- Phase 6 evidence: no unit in the 0.9.0 latent-adjustment power comparison selected rank 1, so the floor is a no-op at the Phase 6 constructions in the modal rank, though not necessarily in every outer-fold selection. The SEA-AD cohorts (n = 80, held-out AUROC 0.63 / 0.66) selected rank 1: a weak stage signal makes every candidate rank tie near 0.5 and parsimony picks 1.

## Goals / Non-Goals

**Goals:**

- One rank rule in one implementation: the floor is a property of the candidate set, so the inner argmax, the per-fold selections, the mode and the pooled fit agree, and the harness and the analysis share it.
- The old rule stays reachable by name, with the hash-omitted spelling that keeps every committed cell id, signature, seed, record and fixture as committed.
- The default flips on evidence from the regime where it matters, with the acceptance rule written before the run.

**Non-Goals:**

- Changing what the generic `motco plsr` subcommand does by default.
- A floor on the diagnostic `pca` and `pls_within_group` methods, whose ranks are fixed by construction, or on the forced-rank override.
- Any claim that the second component carries signal: the floor makes `angle` and `shape` defined, not powerful.
- Re-running any committed study.

## Decisions

### D1. The floor restricts the candidate set inside `plsda_doubleCV`, with the old default

`plsda_doubleCV(..., *, min_components: int = 1)` builds `n_lv_candidates = list(range(min_components, max_components))` and selects by value, `n_lv_candidates[int(np.argmax(mean_per_n_lv))]`, instead of `argmax + 1`. It validates `1 ≤ min_components < max_components` up front, naming both values. The default reproduces `plsda_doublecv_identity.json` byte for byte, which keeps the `plsda-nested-cv` spec's "defaults reproduce the previous output" scenario true.

*Alternative: clamp the selected rank after CV (`max(2, lv)`).* Rejected: the inner folds would still spend fits on rank 1, outer folds would still select it, and the mode would be taken over selections the final fit cannot honour; the rule would say one thing and the record another.

*Alternative: default 2 in the primitive.* Rejected: the primitive is a general PLS-DA cross-validation whose spec pins its default output, and `motco plsr` wraps it as a generic tool. The floor is the trajectory measurement's rule (D2).

### D2. The default of 2 lives in `fit_pls_latent_space`

`fit_pls_latent_space(..., min_components: int = 2)` validates the floor (integer ≥ 1), clamps `max_components` to `[min_components + 1, min(n_features, n_samples)]` through the existing `_clamp_int` (which already raises when the bound is infeasible; the message gains the floor), and passes `min_components` to `plsda_doubleCV`. The forced-rank path is untouched: `_forced_components` keeps `[2, max]` and skips CV, so a diagnostic at rank 2 under a floor of 3 still measures rank 2. `cv_params` gains `"min_components": min_components` only when the floor is above 1, so the latent metadata of a floor-1 fit equals the pre-change metadata. The harness passes `params.resolved_min_components`; the analysis passes its own field (D5).

### D3. The harness field is top-level, with `None` as the hash-omitted spelling of 1

`SimulationEvaluationParams.min_components: int | None = 2`. `__post_init__` rejects non-integers and values below 1, and normalizes `1` to `None`; `resolved_min_components` returns `self.min_components or 1`. `"min_components"` joins `grid._SIGNATURE_OMIT_WHEN_NONE`, and `runtime_metadata["min_components"]` is written only when the field is not `None`. This is the mechanism the orientation flip used, and it gives the same guarantee: floor-1 cell ids, signatures, matched seeds and records equal those produced before the field existed, while the default is hashed into every new cell (whatever its integration method, as the estimator is).

*Alternative: an `integration_params["min_components"]` key.* Rejected: mapping keys have no omit-when-`None` mechanism, so an absent key would have to mean the old rule for hashing and the new rule for the evaluator at once. Pinning `1` as a key in the committed configs would instead change every committed cell id, which the identity tests forbid.

*Alternative: absent means 1 in the evaluator and the config loader injects 2.* Rejected: the harness called directly (scripts, tests, fixture captures) would silently keep the old rule.

### D4. Study config: field, normalizer, axis

`_build_evaluation` routes `evaluation.min_components` through `_normalize_min_components` (`null` → the dataclass default, read from the field so a later change needs no edit; `1` → `None`; anything else passes through for the evaluator to validate). The normalizer joins `_AXIS_VALUE_NORMALIZERS` under `evaluation.min_components`, `_display_axis_value` shows `None` as `1`, and the hint text covers the field, so the floor is a design axis in both forms, with `null` and `1` as distinct values and the floor-1 column byte-identical to the config without the axis. Every committed config under `examples/trajectory_power_study/` gains `"min_components": 1` beside `"orientation_estimator": "pc1"`; `tests/test_rank_floor_default_pins.py` mirrors `tests/test_orientation_default_pins.py` against the existing `committed_config_identity.json` (the loaded field is `None`, the identity holds, every raw config names the key, no design point writes `null` for the axis meaning the old rule).

### D5. Analysis params: the same spelling, omitted from `to_dict()` while `None`

`TrajectoryAnalysisParams.min_components: int | None = 2`, validated in `__post_init__` (integer ≥ 1), with `1` normalized to `None` and dropped from `to_dict()` while `None`, following `covariate_latent_adjustment`. `analyze_trajectories` passes `params.min_components or 1` to `fit_pls_latent_space`. The record's `latent` block already copies `integration_params`, so it states `min_components` exactly when the floor is above 1 (D2); `measurement` is unchanged. With `min_components=1` every 0.9.0 record and CSV is reproduced, which is what `tests/test_analysis_identity.py` now passes (`_identity_params()` names the floor). `motco analyze` gains `--min-components` (default 2) and prints the floor beside the selected rank.

### D6. Pins in the historical diagnostics

Every place that pins `"pc1"` to reproduce a dated report also pins the floor of 1, with the same one-line comment: `simulations/realized_geometry_study.py`, `simulations/specificity.py`, `scripts/latent_rank_probe.py`, `scripts/orientation_measurement_probe.py`, `tests/phase6_fixture_capture.py` and `tests/orientation_probe_fixture_capture.py`. The 0.8.0 and 0.9.0 validation scripts (`scripts/covariate_type_i_check.py`, `scripts/latent_adjustment_type_i_check.py`, `scripts/latent_adjustment_power_comparison.py`) construct `TrajectoryAnalysisParams` and produced their decisions under the old rule; they pin `min_components=1` so their recorded decisions reproduce, with the comment naming the report. `tests/test_latent_adjustment.py` passes `min_components=1` to reproduce `pls_latent_space_identity.json` and gains the companion test that the same input selects rank ≥ 2 at the default. Task 3.6 audits with `grep -rn "SimulationEvaluationParams(\|TrajectoryAnalysisParams(" scripts src tests` that nothing else constructs the parameters without choosing.

*Audit note (3.6).* `grep -rn "SimulationEvaluationParams(\|TrajectoryAnalysisParams(" scripts src tests` plus the full fast suite found three further sites that construct the parameters without naming the floor and whose expectations were captured under the unfloored rule; each now pins `min_components=1` with a comment: `tests/test_study_end_to_end.py::_smoke_config` (byte-identical acceptance report and shard layout, which hash the cell ids), `tests/test_simulation_grid.py::test_null_summary_does_not_enter_the_parameter_signature` (the pinned default digest moved, as it did for the two earlier default flips; the pre-floor digest is restored by `min_components=1` and the pre-orientation-flip digest by `pc1` + `min_components=1`), and `tests/test_phase5_magnitude_profiles.py::test_phase5_artefacts_are_untouched`, whose unpinning helper now strips the `min_components: 1` pin as it strips `"all"` and `"pc1"`. Every other constructor takes the default: `analysis.py`, `evaluation.py`, `grid.py` (production), `showcase.py` (`concat`, unaffected), and the tests whose inputs are large enough for the floor (none needed a larger `max_components`).

*Implementation note (3.4).* The pin test checks identity against `tests/data/orientation_default_identity.json` rather than `committed_config_identity.json`: the latter predates the four net-displacement profiles and covers 16 configs, while the former was captured on 2026-10-15, before the floor existed, and covers all 19. `tests/test_committed_config_identity.py` keeps checking the 16 unchanged.

### D7. The check: a weak-signal Phase 6 construction, two arms, the anchor-bound rule

`scripts/rank_floor_check.py` follows `scripts/latent_adjustment_type_i_check.py` in shape (settings hash, sharded JSONL records, resume, merge, pure decision function, sbatch array) but builds its datasets directly from `SemiSyntheticTrajectoryParams` at the Phase 6 small-n design (`group_stage_sizes [[11, 10, 28], [9, 10, 12]]`, methylation + expression, joint magnitude), weakening the stage signal through the per-omic `delta_*` and `p_dmp` until the smoke run's unfloored null arm selects rank 1 in at least half of the replicates; the strengths are recorded in `settings.json` and the rank-1 fraction in the decision. Each replicate is analysed with `analyze_trajectories` twice, `min_components=1` and `min_components=2`, same CV `random_state` and RRPP seed, no covariates. The grid is the zero-effect anchor, `translation` at e = 1.0, and magnitude, orientation and shape at the Phase 6 power profile's effect levels. The decision (`rank_floor_decision.json`, `rank_floor_rates.csv`) applies the study engine's one-sided anchor bound to the floor-2 arm's anchor and translation rates; power rows are descriptive with pooled-SE differences and the rank distribution per arm and cell. The floor-1 arm's `angle` at rank 1 is reported as the two-point null it is, so the report can show what the floor replaces.

Cost: two analyses per unit; the 0.9.0 Type I check ran 1,500 analyses in 5.1 core-hours on ing, so 14 cells × 500 replicates × 2 arms is of the order of 50 core-hours. The exact grid is settled on the smoke run and reported before submission; nothing is submitted without the user's approval.

*Alternative: evaluate through the harness.* Rejected: the consumer is `motco analyze`, both entry points share the core, and the 0.9.0 scripts already established the `analyze_trajectories` pattern with resumable shards.

*Alternative: reuse the D10 nuisance construction.* Rejected: it selects rank ≥ 2 in every unit; the floor must be measured where the unfloored rule selects 1.

### D8. Versioning

A minor version (0.10.0): the default measurement changes for every caller that does not name the floor. The release follows the decision; see Migration Plan.

## Risks / Trade-offs

- [The second component is noise when the stage signal supports one] → By design: the floor defines `angle` and `shape`, it does not make them powerful. The limitation is documented in the API docs and the report, and the record's `config_spectrum` eigengap and null summary say how resolvable orientation is in each analysis. The check's power rows measure it.
- [A caller who constructs parameters without naming the floor silently changes measurement] → That is the default flip; the pins (D6) and the pin test keep every committed artifact reproducible, and the record states the floor.
- [A test or script that sets `max_components=2` now runs a one-candidate CV] → The clamp raises the maximum to 3 and the record shows it; the task list adjusts the tests that relied on a two-candidate set.
- [Samples with fewer than 3 features or samples cannot satisfy the floor] → An error naming the floor and the counts, not a silent lower rank, consistent with the forced-rank policy.
- [The weak-signal construction may need several smoke iterations to reach the rank-1 fraction] → The calibration loop is predeclared; the strengths are recorded each time and nothing paper-grade runs before the threshold holds.
- [Translation `angle` control may exceed the bound for reasons unrelated to the floor, as it did in the nuisance construction] → The decision follows the engine's convention (controls mandatory); the floor-1 arm's rates are reported so the report can attribute an excess, and `diagnose_floor_excess` names the statistic.
- [`motco plsr` keeps a floor of 1] → A deliberate scope decision stated in the proposal; the flag exists for callers who want the floor there.

## Migration Plan

Additive at the primitive, a default flip at the core. Order: the identity fixtures are confirmed on the unmodified tree; the knob lands with its pins, every identity test passing at the floor-1 spelling; the check's smoke run calibrates the construction; after the user's approval the paper-grade run produces the decision and the report. Release 0.10.0 once the decision is `adopt_floor`. If it is `diagnose_floor_excess`, the knob and the pins stay, the core and analysis defaults revert to 1 in the release, and the roadmap records the excess for the follow-up. Rollback for any caller at any time: `min_components=1`. sead_sexdiff reruns its cohorts at the floor in its own proposal.

### D7 calibration record (2026-10-23)

`scripts/rank_floor_check.py calibrate` (anchor cell, `floor_1` arm, no permutations, base seed 2026; `results/rank-floor-2026-10-23/calibration/calibration.json`). Weakening the per-omic `delta_*` alone first *raises* the selected rank (the stage signal thins across more directions while the held-out AUROC stays near 1) and only drives it to 1 once the signal is near chance; `p_dmp` stays at the Phase 6 value 0.1 so every pool-limited surgery keeps its Phase 6 headroom (translation 98 sites available, orientation and shape 242; all cells fit at every strength).

| `delta_methyl` / `delta_expr` / `delta_protein` | `p_dmp` | R | rank-1 fraction | rank distribution | mean CV AUROC | ≥ 0.5 |
|---|---|---|---|---|---|---|
| 2 / 2 / 2 | 0.1 | 20 | 0.00 | {"2": 8, "3": 12} | 1.000 | no |
| 1 / 1 / 1 | 0.1 | 20 | 0.00 | {"3": 20} | 1.000 | no |
| 0.5 / 0.5 / 0.5 | 0.1 | 20 | 0.00 | {"3": 1, "4": 8, "5": 9, "6": 2} | 1.000 | no |
| 0.25 / 0.25 / 0.25 | 0.1 | 20 | 0.00 | {"4": 1, "5": 1, "6": 3, "7": 5, "8": 3, "9": 2, "10": 3, "11": 2} | 0.920 | no |
| 0.1 / 0.1 / 0.1 | 0.1 | 20 | 0.55 | {"1": 11, "2": 3, "3": 2, "4": 1, "5": 1, "8": 1, "10": 1} | 0.579 | yes |
| 0.15 / 0.15 / 0.15 | 0.1 | 50 | 0.08 | {"1": 4, "2": 12, "3": 8, "4": 4, "5": 5, "6": 6, "7": 2, "8": 3, "9": 3, "10": 2, "12": 1} | 0.702 | no |
| 0.1 / 0.1 / 0.1 | 0.1 | 50 | 0.44 | {"1": 22, "2": 13, "3": 5, "4": 3, "5": 3, "6": 1, "7": 1, "8": 1, "10": 1} | 0.583 | no |
| 0.05 / 0.05 / 0.05 | 0.1 | 50 | 0.58 | {"1": 29, "2": 14, "3": 4, "4": 3} | 0.525 | yes |

Smoke profiles (R = 20, 99 permutations, 14 cells × 2 arms = 560 analyses each; `results/rank-floor-2026-10-23/smoke-delta0.05/`, `smoke-delta0.1/`): both meet the calibration rule on their replicates (0.70 and 0.55), no arm failed, the `floor_2` arm selected rank ≥ 2 in every unit while `floor_1` measured `angle` in one dimension in 183 / 280 and 147 / 280 units. Each smoke decision is `diagnose_floor_excess` on one statistic at 3 / 20 rejections (0.150 against the R = 20 bound 0.147): anchor `angle` at δ = 0.05, translation `delta` at δ = 0.1. At R = 20 one rejection moves a rate by 0.05, so the smoke decisions are not read as verdicts; the gate is the paper-grade run. The floor-1 arm's anchor `angle` rate is 0.00 at both strengths: the rank-1 two-point null never rejects, which is the degeneracy the floor removes.

Cost: median 26–53 s per unit at 99 permutations locally under load, 58 s per unit at 999 permutations (two analyses); 14 × 500 = 7,000 units per strength, so about 80–110 core-hours per strength on ing.

### D7 paper-grade run (2026-10-23)

The user chose both strengths (2026-10-23), with `adopt_floor` requiring both runs to hold. Two 100-shard arrays on ing (`512x1024`; jobs 925691 for δ = 0.05 and 925791 for δ = 0.1), 7,000 units each, 100/100 tasks completed, 0 arm failures, 53.0 + 45.8 recorded core-hours (18–29 s per unit on an EPYC core, far below the local estimate). Both returned **`adopt_floor`**; every control statistic of the floor-2 arm holds the bound 0.0695 at R = 500 (`docs/reports/rank-floor-2026-10-23.md`, `results/rank-floor-2026-10-23/`). The realized rank-1 fraction of the unfloored anchor arm is 0.578 at δ = 0.05 and 0.420 at δ = 0.1 (the smoke replicates' 0.55 was the high side of its R = 50 estimate 0.44); the decision file reports both. Release: 0.10.0 with the defaults of 2 (D8).

## Open Questions

None. The construction strengths and the effect levels were settled on the smoke runs and reported before submission (above).
