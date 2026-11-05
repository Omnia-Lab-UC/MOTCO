# Design

## Context

The stage-span construction (`stage-span-constructions`, MOTCO 0.15.0; archived design D2) defines its plant in coordinates standardized by the **null** population pooled SD. That SD does not depend on the plant. The measurement standardizes by the sample's pooled SD, which estimates the **planted** population's pooled SD, and that one does depend on the plant. The archived design listed this as a risk ("realized geometry drifts from the intended geometry") and reported `population_latent` beside the intended geometry. It did not bound the effect on the off-target statistics.

The scratch decomposition behind this proposal (job 933675, 100 × 199 at five columns, 3,500 units; the 240 local rows reproduce to 1e-13; under the study's p < α rule the rates match the pilot's) established the following. Medians of B − A path length are given as a share of A's path:

| column | shape e = 0.5: scaler route (exact planted SD) | PLS scores, null-scale bend | PLS scores, bend re-solved under the planted SD | `delta` rate, null-scale bend |
|---|---|---|---|---|
| n = 200 | −0.8 % | −0.77 % | +0.04 % | 0.08 |
| n = 80 | −0.8 % | −0.80 % | +0.04 % | 0.03 |
| SEA-AD | −0.7 % | −0.96 % | −0.27 % | 0.00 |
| 3 blocks | −0.9 % | −0.78 % | +0.10 % | 0.11 |
| 4 stages | −1.1 % | −1.58 % | −0.58 % | 0.05 |

The remaining route is the score metric. sklearn's PLS scores are `X·R` with `R = W·(PᵀW)⁻¹`, and `PᵀW` is unit upper triangular, so R is W followed by a shear. The span is the same as `span(W)`, but lengths and angles inside it are not preserved. An orthonormal projection of the same population means (`qr(R)`) leaves only the scaler route. The sheared coordinates add the residual, which is largest at rank 4 (four stages) and on the unequal SEA-AD table.

A prototype of the fixed point (Decision D1) converged at all five columns, for every mode at e up to 1 (orientation up to 2), in at most 35 iterations. The planted-coordinate geometry was exact to 1e-13 relative.

## Goals / Non-Goals

**Goals:**

- A stage-span plant that is exact in the population analogue of the coordinates the measurement standardizes to, opt-in and byte-identical when unset.
- A paper-grade measurement of sensitivity given visibility, with a decision rule whose false-failure rate under the null is controlled despite about 150 bound tests.
- A committed, reproducible account of both leak routes.

**Non-Goals:**

- Removing the score-metric route. It belongs to the measurement, and changing the measurement is a separate change on the roadmap.
- Making the planted scale the stage-span default.
- Supporting `span_alignment` under the planted scale.

## Decisions

**D1. The planted scale is a fixed point of the plant.**

- Let N be the null stage means over the span layers, c their stage-weighted pooled mean, W the within-cell variance (the reference diagonal) and w the realized cell weights of both groups. Start from s₀ = the null pooled SD and iterate:
  1. z = (N − c)/s;
  2. E = `edit_stage_configuration(z, mode, e)`;
  3. B = E·s + c;
  4. s ← √(W + between-cell variance of [N; B] under w).
- Stop when max |s_new/s − 1| ≤ 1e-12, with a limit of 200 iterations; at the limit, raise a validation error naming the mode, the effect and the last residual.
- At the fixed point s*, B's configuration in s*-coordinates is the mode's edit of A's, so the intended geometry is exact where the sample scaler looks.
- **Alternative, solving only the bend's rescale** (what the scratch did). It fixes `shape`'s `delta` and nothing else. The rotation keeps its null-scale definition, under which the exact planted SD already moves its path lengths (+0.02 to −0.14 at e = 0.5, up to 0.25 % of A's path at four stages), and its off-target `shape` reads 0.07 at n = 200 (R = 100, cause not separated from that column's shared noise). The user chose the general fixed point.
- **Alternative, a per-replicate sample scale.** That is circular (the plant would follow the replicate's noise), and a dataset would no longer be the null dataset plus a fixed plant.
- **Alternative, damping or root-finding.** Not needed: plain iteration converged everywhere in the prototype. The limit and the error keep a non-contracting case from passing silently.
- Magnitude needs no iteration (scaling about the centroid commutes with a per-feature scale), and the first update returns the same plant. That is why the planted-scale magnitude datasets equal the pilot's. The reproduction check relies on it.

**D2. Settings, hashing and validation.**

- `span_scale ∈ {null, planted}`. `null` normalizes to `None` and joins `_SIGNATURE_OMIT_WHEN_NONE`, so every existing cell id, signature, record and fixture is unchanged.
- It is valid only under `stage_span`. The setting is cleared, together with `construction` and `span_layers`, in zero-effect and Type I null cells (`enumerate._without_construction`), so planted profiles share the surgery anchor.
- It is rejected as an OFAT or design axis, for the reason in archived D7: two columns would share one anchor dataset.
- `planted` with `span_alignment` is refused. The out-of-span direction is drawn in null-standardized coordinates, and defining it under a moving scale is a separate question no profile needs.

**D3. Truth.** Under `planted`, the `stage_span` block's `units` is `"planted_standardized"` and gains `scale_iterations` and `scale_residual`. `path_length`, `mean_step`, `intended_geometry` and `plant_norm` are computed in s*-coordinates. The top-level truth gains `span_scale: "planted"`. Without the setting, the truth is unchanged.

**D4. Profiles.**

- Three post-floor siblings in the `design-guidance` family, base seed 1028, 500 × 999, span capture on. Explicit `design_grid.points` over:
  - `generator.n_samples`, `generator.group_stage_sizes` and `generator.n_stages`;
  - `evaluation.integration_params.layers`;
  - `generator.span_layers` (stage-span profiles only).
- The stage-span profiles' baseline is n = 200 at three stages on two blocks, as in the pilot and the slice. The SEA-AD point keeps the default `n_samples` (600), as the size-table conflict rule requires. The three-block point sets `layers` to null and `span_layers` to all three.
- The surgery profile's baseline is the SEA-AD table, and its second point is the four-stage column. Seeds depend only on the family and the replicate index (`derive_replicate_seed`), so its anchors are the stage-span profiles' anchor datasets even though the phase, and therefore the signature, differs. The reproduction check reads datasets and p-values, not cell ids.
- One effect list per profile means stage-span translation runs at {0.1, 0.25, 0.5}. The pilot measured translation at e = 1.0 at R = 100 (held). Per-mode effect lists remain the documented alternative from the slice, at the cost of an engine change.
- Planning cell counts follow the pilot's enumeration rule (anchor plus mode × effect per column, plus the Type I baseline cells): about 52, 22 and 22 cells, or 48,000 units. Tests pin the exact counts.

**D5. The decision's bound rule controls the number of tests.**

- The pilot's literal rule (every rate at or below α + 2·SE) has a false-failure probability of 0.84 at m = 60 and 0.94 at m = 90 tests at R = 500 when every true rate is α. A k = 3 ceiling alone still fails 0.15–0.22 of null runs.
- The rule for each criterion over m rates instead:
  - **Count.** The number of rates above α + 2·SE is at most the 95th percentile of Binomial(m, p₂), where p₂ is the exact probability that a Binomial(R, α)/R rate exceeds that bound (0.0303 at R = 500, 0.0282 at R = 100). That allows 4 exceedances at m = 60 and 6 at m = 90.
  - **Ceiling.** No rate is above the Bonferroni ceiling, the smallest rate a Binomial(R, α) rate exceeds with probability at most 0.05/m: 0.084 at R = 500 and m = 60–90, 0.130 at R = 100 and m = 60.
- Each half fails a null run with probability of about 0.05, so a criterion's false-failure rate is roughly at most 0.10 under independence. R is each row's available replicates, so p₂ and the ceiling are computed per row. The allowed count uses the median R of the criterion's rows. The JSON records every quantity.
- **Alternative, pooling over effects within a column.** It would hide a leak that grows with e, which is exactly what this change investigates.
- **Alternative, a Bonferroni bound alone.** It is too lenient for a systematic small excess spread over many cells, which the count catches.
- Cells share noise within a column (matched seeds), so exceedances are positively correlated. A column-wide shift such as the pilot's n = 200 `angle` 0.10 shows up as clustered exceedances. The count rule is approximate there; the report names any cluster.

**D6. Criteria and verdict.**

- **P1:** anchors (5 columns × 3 statistics, m = 15) and stage-span translation (5 × 3 × 3, m = 45) as one family, m = 60.
- **P2:** median capture ≥ 0.90 in every stage-span magnitude, orientation and shape cell, a deterministic premise check with no multiplicity.
- **P3:** off-target rates of the stage-span magnitude, orientation and shape cells, m = 90.
- The verdict is the first failure in the order `size_excess` → `capture_shortfall` → `specificity_excess`, else `sensitivity_established`. Every failure is listed. Target power, `informative_from`, capture and the visibility contrast are descriptive (as in archived D11).
- The decision reuses the pilot reading's report loader and bound helpers by import, as other validation scripts import their predecessors.

**D7. The leak check is the scratch decomposition, promoted.**

- `scripts/stage_span_leak_check.py run --config <profile> --shard i --n-shards n [--span-scale null|planted] [--modes …] [--effects …] --replicates R --permutations P --out-dir …` and `merge --out-dir …`.
- Per unit it rebuilds the cell's dataset from the profile's matched seeds, fits the production latent space once and evaluates with RRPP when P > 0. From the fitted `LinearLatentMap` it reports the six checkpoints listed in the spec.
- The exact planted SD is computed from the population trajectories and the reference diagonal. Under `planted` this checkpoint's difference is zero to rounding, which is the gate's first condition.
- `merge` writes `leak_units.csv`, `leak_cells.csv` (medians, shares of A's path, rates with p < α) and `summary.json` (settings, units, failures, the config's sha256).
- An sbatch template follows `scripts/latent_space_comparison_array.sbatch`.

**D8. The pre-flight is a coarse screen; the decision judges specificity.**

- At R = 100 the rate conditions cannot tell the null-scale bend's 0.11 at three blocks (which would fail P3 at R = 500) from chance. The scratch data passes the count rule with one exceedance below the 0.130 ceiling.
- The gate therefore adds the population-level condition: the planted arm must shrink the score-coordinate residual at every column, and its scale route must vanish exactly.
- The rate conditions remain as a screen for gross failures.
- The null arm (0.25, 0.5, 1.0) reproduces the pilot's rates at the pilot's columns and gives the report its before-and-after table.

**D9. Visibility contrast at the slice's columns comes from committed aggregates.** The decision reads the slice's committed `core_geometry/report/design_point_operating.csv` for surgery target power at n = 200 and n = 80 (two blocks) and n = 200 (three blocks). It reads the pilot's committed `span_capture.csv` for surgery capture at n = 200 and n = 80 on two blocks. Each entry is labelled with its source and grade, and three-block surgery capture is reported as unavailable.

## Risks / Trade-offs

- **[The score-metric residual makes P3 fail at four stages or on the SEA-AD table.]** The planted bend still leaves −0.58 % at four stages, where the null bend's −0.77 % produced 0.08 at n = 200. → The verdict would be `specificity_excess`, naming the cells. The report attributes them through the leak check, and the roadmap item on orthonormal coordinates is the remedy. The gate stops the run only for a gross failure; a P3 failure at paper grade is a finding, not a defect.
- **[The planted plant changes the orientation and shape datasets relative to the pilot.]** The stage-span paper-grade geometry cannot be compared replicate-for-replicate with the pilot's stage-span geometry. → The anchors, the magnitude cells and the surgeries still reproduce, and the leak check's null arm reproduces the pilot's rates.
- **[The fixed point fails to converge for some seed.]** → It raises a validation error, which the study records as a failed unit (`error_policy: record`). The gate requires zero failed units in the pre-flight, and the paper grade reports failures by cell.
- **[Correlated tests within a column inflate the clustered exceedance count.]** → The rule's approximation is documented in the decision JSON and the report.
- **[Cost overrun on `ing`.]** About 290 core-hours is planned at 18–28 s per unit, and 999 permutations add RRPP time. → The smoke run fixes the figure before submission, the array is throttled as before, and n2 and n11 are excluded with the reason stated.
- **[The slice's records are needed for anchor reproduction.]** They stay on the cluster (`/home1/tgonza/MOTCO-design-guidance-slice`). → The reproduction runs on `ing`, with `anchor_reproduction.py --by-column` and `--by-cell`. If a source's records are gone, that check is reported as unavailable, not silently skipped.

## Migration Plan

1. Implement `span_scale` with tests, and check byte identity for every committed config and fixture.
2. Commit the leak check and the decision script with tests.
3. Commit the three profiles and register them in `POST_FLOOR_CONFIGS`; regenerate the post-floor fixture.
4. Run a local smoke.
5. Run the pre-flight on `ing` from a clean detached worktree of the pushed branch, then read the gate; if it fails, stop and report.
6. Run the paper-grade array, then the reproduction checks and the decision.
7. Write the report and the documentation, bump the version to 0.17.0, and archive.

Rollback means reverting the change. Without `span_scale`, nothing reads or writes the new fields.

## Open Questions

- Whether `motco simulate` should expose `--span-scale`. The CLI exposes no stage-span setting today, and adding one later changes no spec here.
- The final cell counts and the core-hour figure, fixed by the tests and the smoke run.
