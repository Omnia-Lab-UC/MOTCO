# Design

## Context

See `proposal.md` for the motivation. The facts that shape the approach:

- **`RRPP` draws one `rng.permutation(n)` per permutation** from `default_rng(seed)` on the serial path, and one per draw from a spawned seed in each `_RRPPWorker`. The split interaction test and the scan draw the same serial sequence themselves (`motco.stats.interaction`).
- **The diagnosis's loop is the reference construction.** `predictive_size_excess_check.within_stage_indices` visits the stages in declared order and sets `idx[rows] = rows[rng.permutation(rows.size)]` from `default_rng(the unit's RRPP seed)`; `restricted_null` then evaluates `ŷ_red + e*[idx]` as `RRPP` does. Its counts are committed in `results/predictive-size-excess-2026-11-13/rate_table.csv`.
- **The harness has no residual correction.** Only `motco.analysis` passes `residual_correction` to `RRPP`. The simulation profiles therefore run uncorrected, as every simulation run has, and the plasmode runs with `df`, as the case study does.
- **Both entry points know the stage order.** The harness design sorts its stage levels (`_sort_levels`); `motco.analysis` recodes the declared `stage_order` to ordinal strings before building the same design.
- **Committed drivers are frozen.** The pilot (`predictive_subspace_check`) and the diagnosis (`predictive_size_excess_check`) are parents of this run and are imported unchanged. `simulation_settings(base_seed=)` and `plasmode_settings(replicates=)` already take the two values this run changes.
- **The `paper` label stream is read only to replicate 499.** Every committed run that reads the plasmode's labels stops there (checked in each run's `settings.json`), and `plasmode_check.STREAMS` is closed, so replicates 500–999 are the unread label draws that need no edit of a parent.

## Goals / Non-Goals

**Goals:**
- One definition of the stratified draw in the primitive, reproducing the diagnosis's sequence exactly and today's draws exactly when absent or with one stratum.
- The option reachable from the harness, `motco.analysis` and `motco analyze`, invisible while unset.
- A size reading of the within-stage null where the truth is known and stage dispersion is unequal, beside the designs already used.
- A decision that separates what the null changes from what the coordinates change.

**Non-Goals:**
- Changing `motco.stats.interaction`, the generator or the study engine.
- Explaining the size of the anchor's excess beyond the dispersion mechanism.
- Any reading at paper grade for `predictive_orthonormal`.

## Decisions

### D1. The primitive takes labels; the entry points take a name

`RRPP(..., strata=None)` takes one label per row and knows nothing about stages. The harness and the analysis take `permutation_strata="stage"` and derive the labels from their design.

- **Why.** The primitive is design-agnostic everywhere else (`contrast`, `LS_means`), and exchangeability strata are a general notion. The entry points own the design, so a name there cannot be mis-aligned with the rows.
- **Alternative: a boolean `within_stage`.** Rejected: the other opt-in settings name a value (`residual_correction`, `covariate_latent_adjustment`), and a name leaves room for other strata without a second flag.
- **Alternative: labels in the entry points too.** Rejected: nothing asks for other strata, and a free label vector in a hashed harness parameter would put data in cell signatures.

`PERMUTATION_STRATA = ("stage",)` lives beside `LATENT_COORDINATES` in `simulations/evaluation.py`, and `motco.analysis` imports it.

### D2. The draw sequence is the diagnosis's

Per draw, the strata are visited in ascending label order (`np.unique`), and each stratum's rows are permuted among themselves by one `rng.permutation(rows.size)` from the single generator. The entry points pass each row's stage *position* in the design's stage order, as integers.

- **Why integers.** The analysis's working labels are ordinal strings, and `"10"` sorts before `"2"`. Positions make ascending order the declared stage order at any number of stages.
- **Why this sequence.** It is `within_stage_indices`'s, so `analyze_trajectories(..., permutation_strata="stage")` at a unit's seed returns the diagnosis's p-values exactly. That turns the diagnosis's committed table into a reproduction gate for the package option (D10). With one stratum the sequence is `rng.permutation(n)`, today's.
- **Alternative: one global permutation, then a stable regrouping by stratum.** It consumes one draw per permutation whatever the strata, but it is another sequence, so nothing committed would reproduce.

### D3. One draw helper, and the absent path untouched

A module-level helper returns the next index vector from a generator and the strata's row groups. `RRPP`'s serial loop and `_RRPPWorker` both call it. With `strata=None` both keep their literal `rng.permutation(n)` call, so the absent path consumes randomness exactly as before. Validation (length, missing labels) runs with the other input checks, before the reduced fit.

A one-row stratum is allowed: the row is held fixed, which is what an exchangeability block of one means. The primitive does not count distinct permutations, as it does not today.

### D4. The residual corrections keep their definitions under strata

`df` stays κ = √(n / (n − r)) and `leverage` stays 1/√(1 − h_ii), applied once before any draw, with or without strata.

- **Why.** Each option keeps one meaning, and the diagnosis's arms (`df` under within-stage draws) reproduce.
- **What it costs.** When the reduced model holds the stage indicators, the residuals sum to zero within each stage, and a within-stage draw samples a cell's residuals without replacement from its stage. The finite-population factor then already accounts for the K stage means: for a stage-only reduced model the null variance of a cell mean about its stage mean is exact with no correction. Under strata the scale that matches second moments is √((n − K)/(n − r)), so `df` over-inflates by √(n/(n − K)): 1.9 % at SEA-AD (n = 80, K = 3) and 0.8 % at n = 200. The direction is conservative.
- **Alternative: a strata-aware scale.** Left out: it would be a third correction with no reading behind it, and the anchor's within-stage rates under `df` (0.038 / 0.054 / 0.038 for production) show no visible conservatism. The report quotes the factor; a later change can add the scale if a reading asks for it.

### D5. The harness setting

`SimulationEvaluationParams.permutation_strata: str | None = None`.
- Validated in `__post_init__` against `PERMUTATION_STRATA`, before any fit.
- Added to `grid._SIGNATURE_OMIT_WHEN_NONE`, so every committed cell id, signature and record is unchanged; `runtime_metadata["permutation_strata"]` only when set.
- `evaluate_semisynthetic_trajectory` maps each row's stage to its index in `design.stage_levels` and passes the vector as `strata`. Nothing before the `RRPP` call changes, so the setting works with a prebuilt `latent=` and with every method and coordinate system.
- With `permutations = 0` the setting is accepted and has no effect, like the seed.

### D6. The analysis setting, its record and its rejections

`TrajectoryAnalysisParams.permutation_strata` follows the same rules: validated in `__post_init__`, dropped from `to_dict()` while `None`, passed as the working frame's stage positions with the covariate-adjusted models and the residual correction. `measurement.permutation_strata = {"value": "stage", "sizes": [...]}` (the stage sizes in declared order) is written only when the setting is set and permutations are drawn.

It is rejected with `span_split_test` and with `interaction_scan`, naming both settings.
- **Why.** Both tests draw their own nulls from the same seed and were validated under unrestricted draws (item 16). Accepting the combination would either mix two nulls in one record under one flag, or restrict draws that no run has read.
- **Alternative: extend both tests.** A later change, with their own readings, if the null is validated.

Attribution, the latent adjustment, weighting, sparsity and `latent_coordinates` act before the null or beside it and combine freely.

### D7. The CLI flag ships with the option

`motco analyze --permutation-strata stage`, default absent, with a summary line naming the strata.

- **Why now, when `latent_coordinates` got no flag.** That setting changes what is measured. This one changes only the reference distribution, to one that assumes less than the default. The diagnosis already read production's size under it on the anchor. And the case-study repositories run the CLI, so a typed field alone would not reach them.
- **If the verdict is `size_excess`.** The flag stays, and the help text, the documentation and the report say that the null failed size and where. Removing a flag is a smaller follow-up than adding one under a separate release.

### D8. Eight arms on one fit, through the package option

Per fitted unit: the dataset, the production integration and the three predictive integrations are built as the pilot builds them (`integrate_semisynthetic_dataset`, `pls_coordinates_from_pls`; on the plasmode `plasmode_production_params` and `plasmode_arm_params` at the unit's rank). Each of the four matrices is then evaluated twice, with `permutation_strata` absent and `"stage"`, at the unit's RRPP seed.

- **No driver loop.** The within-stage arms go through `evaluate_semisynthetic_trajectory` and `analyze_trajectories`, so the run validates the option users will call.
- **Ties to the parents.** A test requires this driver's four unrestricted arm records to equal the pilot's `evaluate_harness_unit` / `evaluate_plasmode_unit` records on the same unit. On every plasmode unit the `construction_identity` check requires each within-stage arm's p-values to equal `predictive_size_excess_check.infer(..., indices=within_stage_indices(...))` on that arm's latent matrix.
- **Dropped arms.** The pilot's `forced_k_minus_1` and `predictive_rank_plus_one` answered questions about the fit rank and are not re-read.

### D9. The simulation profile is a fresh seed family

`predictive_subspace_check.simulation_settings(base_seed=2038)` on the pilot's designs, cells, replicates (100) and permutations (199); 7,800 units.

- **Why not the pilot's units.** Its datasets are where `predictive_orthonormal` was preferred among three arms, so they are spent for judging its gain and specificity. The pilot applied the same standard to the probe's datasets.
- **What is given up.** The exact reproduction of the pilot's simulation table. The driver's harness path is gated by `weak_size` instead (1,000 units, four arms, exact), and its analysis path by the plasmode (D10).
- **The surgery cell stays,** descriptive: it is the one cell where the Phase 6 construction is read, and production's `angle` there under the new null is of interest.

### D10. The plasmode profile is the anchor on 1,000 label draws

`plasmode_settings(export, replicates=1000)` restricted to the anchor cell, every unit at the committed fidelity rank of 5 (`forced_components`), as the diagnosis ran its anchor.

- **Replicates 0–499 are a gate.** The four unrestricted arms must equal the pilot's committed anchor counts (21 / 24 / 13, 68 / 11 / 38, 67 / 9 / 38, 65 / 10 / 39), and the four within-stage arms the diagnosis's (19 / 27 / 19, 31 / 17 / 33, 29 / 17 / 32, 28 / 17 / 33). The second comparison is the package option against the diagnosis's independent loop on real data.
- **Replicates 500–999 are the reading.** The diagnosis's `shape` held by one or two counts, and its report says a rerun on new label draws could put it over. These are those draws. Only they enter the size criteria; the pooled 1,000 are descriptive.
- **No plant cells.** A plant unit needs one cross-validation (about 125 s on `ing`, 10 core-hours for 300 units), and the diagnosis already read all eight arm-by-null combinations on them through the construction the gate ties to the package. Its plant rows are quoted from its committed table.
- **Pre-flight.** As in the diagnosis, on replicates 0 and 500: the cross-validated fit selects rank 5 at AUROC 0.635556 and equals the forced call exactly. The fit does not see the group labels, so the rank holds for every anchor draw.
- **Limit.** All 1,000 units share one fit, so every anchor rate is conditional on it, as in both parents.

### D11. The dispersion profile reads the primitive where the nulls differ

Gaussian matrices on the SEA-AD cell table, built in the driver: `Y = c[stage] + ε`, `ε ~ N(0, v_stage · I_d)`, with d ∈ {2, 5}, no group difference, the design from `build_simulation_trajectory_design`, and both nulls read by `RRPP` (with and without `strata`) at one seed, net displacement, no correction.

- **Why it exists.** Every generator design has noise that is exchangeable across stages, so both nulls are expected to hold there and the simulation profile can only show that the restriction costs nothing. The within-stage null's claim is that it holds where dispersion differs by stage. This profile is the one place where that is read with the truth known.
- **Why at the latent level.** The generator has no stage-dispersion field. Adding one touches a module whose fidelity is pinned against InterSIM, and roadmap open item 8 already owns a dispersion field. The null acts on the latent matrix, so a latent-level construction tests exactly the step in question.
- **The patterns** are the diagnosis's median relative dispersions by stage (its report, section 9): `predictive` = `predictive_orthonormal`'s 1.120 / 1.278 / 0.802 and `scores` = 0.805 / 0.921 / 1.137; `equal`; and each measured pattern with every multiplier squared, to read the dependence on the inequality.
- **The centroids:** (0, 0), (s, 0) and (s/2, 0.9 s) in the first two dimensions, zero elsewhere, with s solved so the population stage share of variance at d = 2 under equal dispersion is 0.84, the anchor's in-sample stage R². d = 5 adds three noise dimensions, the anchor's fit rank.
- **Size.** 1,000 replicates at 999 permutations per cell; ten cells; 10,000 units of about half a second.
- **What it cannot show.** It is isotropic Gaussian noise with no fit. It validates the null under the mechanism; it does not explain the anchor. The probe's unrestricted `delta` under the measured predictive pattern was 0.066–0.070 against the anchor's 0.13.
- **Its role in the decision.** `matrix__within_stage` on the ten cells is a size family of the null (m = 30). The unrestricted arm's rates are descriptive: it is expected to exceed under `predictive*` and to be conservative under `scores*`.

### D12. The decision has two layers

**Layer 1, the null on the production scores** (`pls_cv__within_stage` against `pls_cv`):
- *Size*: `judge_rates` on four families: simulation controls (m = 42), `weak_size` (m = 6), `anchor_fresh` (m = 3), dispersion (m = 30). All hold.
- *Protected power*: the pilot's paired loss rule on the 54 simulation target cells at δ ≥ 0.25.
- *Specificity*: at δ = 2.0 the bound rule on the 36 off-target readings. In the band, a reading counts only if its rate is above the size bound **and** its paired excess against `pls_cv` has p < 0.025, under the count-plus-ceiling rule over 72.
- *Gain*: recorded per target cell, not required.

**Layer 2, the candidate under the same draws** (`predictive_orthonormal__within_stage` against `pls_cv__within_stage`): the pilot's S1–S3 and G with the reference replaced, on the simulation, `weak_size` and `anchor_fresh` families. It bears on the verdict only if layer 1 holds.

Why each choice:
- **Gain is not required of the null.** A null that assumes less and costs nothing is worth having whether or not it rejects more.
- **The band's specificity rule differs from the pilot's.** The pilot compared two measurements under one null, so any paired excess was a leak. Here the null changes. Where the unrestricted null is conservative, a correct null raises a deflated off-target rate toward α, and a bare paired excess would count that as a failure. Requiring the rate to be above the bound as well keeps restored size from reading as a leak. The bound alone would not do, because production itself leaks in the band.
- **The candidate's reference is production under the same draws.** The diagnosis showed production gains under within-stage draws too, so only this comparison isolates the coordinates. If the candidate gains over `pls_cv` but not over `pls_cv__within_stage`, the simpler change is the null alone.
- **One decisive candidate.** The roadmap names `predictive_orthonormal`, the only arm that passed both specificity families in the pilot. The other two predictive systems are graded by the same criteria and reported.
- **The plasmode's plants are not in protected power.** They would be reproduced readings whose values are known.
- **Order.** Layer 1 reads first failure in the order size, power, specificity, as the pilot ordered its criteria.

Verdicts and the consequence statements fixed with rule version 1:

| verdict | when | consequence |
|---|---|---|
| `reproduction_failed` | a reproduction or the pre-flight fails, or a unit is inconsistent | No reading is interpreted. The cause is fixed and the affected units are rerun under the same rule. |
| `size_excess` | layer 1 size fails | The within-stage null does not hold size in the named family. The option stays an opt-in diagnostic, documented as such; no case study reports it as a test, and `predictive_orthonormal` is not reopened. |
| `power_loss` | size holds, protected power fails | The null holds size and costs power on the production scores. It stays an opt-in null for data whose stage dispersion is unequal, documented with the cells it loses; `predictive_orthonormal` is not reopened. |
| `specificity_excess` | size and power hold, specificity fails | The null holds size and reads more off-target rejections than the bound allows. It stays opt-in, documented with those readings; `predictive_orthonormal` is not reopened. |
| `within_stage_validated` | layer 1 holds, the candidate does not qualify | The within-stage null is a validated opt-in null for the production measurement. A separate change may apply it to the SEA-AD case study and decide the default. `predictive_orthonormal` stays an opt-in diagnostic for the named reason. |
| `propose_paper_grade` | layer 1 holds, the candidate qualifies | As `within_stage_validated`, and a separate change runs `predictive_orthonormal` under within-stage draws at 500 × 999 with a rank rule for a fixed measured dimension. Only its verdict may change a default. |

The rule is amended only before any array unit runs, with the amendment dated in `notes.md`.

### D13. Driver, parents and compute

`scripts/within_stage_permutation_check.py`: `run --profile simulation|weak_size|plasmode|dispersion`, `merge --profile …`, `decide`, plus `preflight` for the plasmode; `scripts/within_stage_permutation_array.sbatch` with a profile variable and a dry-run variable.

- **Parents bound by sha256:** the pilot's three `settings.json` and its `weak_size` and `plasmode` `rate_table.csv`; the diagnosis's `settings.json` and `rate_table.csv`; and through the pilot's settings builders, the comparison's, the span-split run's, the plasmode's and the production selection's files.
- **Reproduction** (`reproduction.json` per profile): `weak_size` against the pilot (four arms × two cells × three statistics); plasmode cell `anchor` against the pilot (12 rows) and the diagnosis (12 rows). `simulation` and `dispersion` reproduce nothing.
- **Helpers reused, not edited:** the pilot's `judge_paired`, `judge_gain`, `allowed_count` and consistency checks where their signatures allow another reference arm; otherwise the driver states the same rule over `judge_rates` and the production selection's exact tests, and a test pins agreement with the pilot's function on the pilot's own committed tables.
- **Compute, before the smoke:** simulation 7,800 × about 21 s ≈ 46 core-hours; `weak_size` 1,000 × about 28 s ≈ 8; plasmode 1,000 × about 6 s ≈ 1.7 plus two pre-flight cross-validations; dispersion 10,000 × about 0.6 s ≈ 1.7. About 57 core-hours on `ing`. The smoke replaces these, and the cap is set at about 1.5 times its projection.
- **Hosts.** `plasmode` and `dispersion` fit on the workstation in under an hour at seven pinned workers. `simulation` and `weak_size` go to `ing` from a clean worktree at a pushed revision.

## Risks / Trade-offs

- **[Both nulls hold on every generator design, so the simulation shows only that the restriction is harmless]** → the dispersion profile reads the case where they differ, and it is a decisive size family.
- **[The dispersion profile is a toy]** → it is stated as a validation of the null under one mechanism, not as an account of the anchor; the report quotes the gap between its rates and the anchor's.
- **[The plasmode's within-stage readings on replicates 0–499 are already known]** → they are a gate only; the size criteria read replicates 500–999.
- **[Every anchor draw shares one fit]** → stated as a limitation; independent real-noise fits wait for roadmap open item 2.
- **[The simulation profile is pilot grade, 100 × 199]** → the count-plus-ceiling rule is built for families of such readings, and three other size families run at 999 permutations. A default change would need its own paper-grade run.
- **[`df` over-corrects by up to 2 % under strata]** → conservative, quoted in the report (D4).
- **[A band off-target reading above the bound that production also exceeds]** → it counts only with a significant paired excess, so the null is not charged with production's own leak.
- **[The flag ships before the verdict]** → opt-in and default off; on `size_excess` the help text and documentation say so (D7).
- **[Small stages have few distinct within-stage permutations]** → documented in the docstring; the primitive does not check, as it does not for small n today.
- **[A helper of a parent driver cannot take another reference arm]** → the rule is restated in this driver and pinned against the parent's function on the parent's committed tables (D13).

## Migration Plan

The option is additive and absent by default; no committed configuration, fixture or record changes. The package change, with its tests and documentation, lands before the driver, and the run is launched from a pushed revision. Rolling back is reverting the change; no stored artefact depends on it.

## Open Questions

- **Host and cap for the array.** The local smoke gives the per-unit times; the cap and the choice between `ing` and the workstation for `weak_size` follow from it, and launch waits on the user's approval.
