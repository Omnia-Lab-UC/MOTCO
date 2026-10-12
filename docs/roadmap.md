# MOTCO Scientific Roadmap

## Objective

Develop MOTCO into a scientifically defensible workflow that can:

1. detect whether groups follow different multivariate molecular trajectories;
2. distinguish progression magnitude and orientation reliably;
3. treat shape with an explicitly validated interpretation; and
4. identify stable feature and pathway drivers after a global trajectory difference is detected.

The roadmap prioritizes interpretability and validation before larger simulations or additional integration methods.

## Current position

### Established

- Methylation is converted from B values to M values before standardization and PLS integration.
- Under that contract, magnitude is specific and well behaved: the 50-replicate pilot reached delta power 1.00 at the largest effect, with angle and shape rejection rates of 0.02.
- PLS is the preferred production measurement space for MOTCO's Euclidean magnitude and orientation statistics.
- A controlled inverse-PLS study showed that latent interventions reconstruct reproducibly into feature space, with round-trip errors near machine precision.
- A 45-degree latent orientation change reconstructed as approximately 42 degrees in feature space with negligible magnitude change.
- Latent orientation differences represent different relative combinations of stage-associated molecular patterns, usually as dense loading-aligned feature changes.
- SNF is not intrinsically unsuitable for integration, but its graph-spectral geometry is not naturally aligned with MOTCO's current Euclidean path-length and angle statistics.
- Shape invariance has been audited and corrected: production `shape` removes translation, uniform scale, and any orthogonal transformation. Reflections are aligned away at every ambient dimension (`unify-shape-reflection-policy`, resolving audit finding F3); the earlier "reflections retained" contract bound only in a full-rank latent space and was vacuous everywhere else.

### Not yet established

- The production orientation and shape feature-surgery modes are not geometrically pure before integration; orientation's `shape` response (0.998 at paper grade, 0.99–1.00 at the design point under CV rank) is a descriptive gate role, not evidence about the shape estimator ([ladder report](reports/latent-rank-ladder-2026-09-08.md), §6). Its predeclared "projection-associated" label is **refuted**: orientation's response is present *within* each omic block (max block `angle` 90.08° against joint 89.93°), and under the recalibrated materiality rule no pair in the study localizes as projection-associated at all ([diagnostic](reports/magnitude-construction-diagnostic-2026-09-10.md), §1 and §4). It is construction-present cross-talk of a per-omic feature permutation; the paper must not repeat the projection-associated framing.
- The orientation power claim is n-conditional at the isotropic endpoint and eigengap-dependent: 0.840 ± 0.016 at n = 1200, ρ = 0 at paper-grade precision (500 × 999), agreeing with the two pilots' 0.85–0.88 ([paper-grade report](reports/phase5-paper-grade-2026-09-10.md), §5; [design-point report](reports/phase5-design-point-pilot-2026-09-08.md), [ladder report](reports/latent-rank-ladder-2026-09-08.md)). Stratified by recorded eigengap the pooled figure straddles the 0.80 floor — 0.754 / 0.886 / 0.880 by tercile — so a real-data result must be read against its own eigengap, never the pooled number.
- The `angle` null tracks its own observed statistic (slope 0.811) and is load-bearing rather than corrigible ([report](reports/angle-null-pivotality-2026-09-01.md)); the recorded eigengap and null-width dispersion must accompany any real-data orientation result.
- Cross-replicate driver stability has no design that holds the driver set fixed, so it is not claimable.
- The *former default* magnitude surgery (`magnitude_kind='all'`) is not size-pure: at paper-grade precision it moves `angle` (0.276) and `shape` (0.742) far past the α + 2·SE control bound, which held the Phase 4 gate at **HOLD** ([paper-grade report](reports/phase5-paper-grade-2026-09-10.md), §6). **The mechanism is established** ([magnitude-construction diagnostic](reports/magnitude-construction-diagnostic-2026-09-10.md)): it scales `delta_methyl` alone, so the surgery is exactly size-pure *within* every omic block and the entire off-target response is the geometry of concatenating one block that grew with two that did not. The remedy — the `joint` construction, every omic's δ scaled together — is now a production `magnitude_kind` and **has been measured at paper grade** ([magnitude re-measurement addendum](reports/phase5-magnitude-remeasurement-2026-10-08.md)): both controls met (0.024 / 0.000), gate PROCEED. Since `default-joint-magnitude-construction`, `joint` is the generator and `motco simulate` default; the committed pre-flip study configs and diagnostics pin `all` explicitly so their signatures and datasets are unchanged. Specificity remains unestablished for `all` (and for the methylation-only `extremes` and `shape_kind='magnitude'` scalings), so the paper must still say which construction a magnitude result is about.
- Orientation `angle` power does not rise with effect: it saturates at 0.84–0.88 from e = 0.25 because the `angle` null widens with the observed statistic, so the strictly-monotone form of the power claim is unsupported (the 0.80 floor is still cleared).
- At the Phase 6 case-study design (n = 80, three stages, two blocks) only `delta` has useful power, and none of the three statistics is mode-specific for orientation or shape ([small-n report](reports/phase6-small-n-2026-10-09.md)). Type I holds (0.010–0.040). magnitude/`delta` reaches 0.542 at a 10 % joint size change and 0.956 at 25 %. orientation/`angle` and shape/`shape` reach only 0.300 and 0.228 at the strongest constructions the generator allows, and those constructions are detected more often by other statistics (orientation by `shape` 0.464 and `delta` 0.354; shape by `delta` 0.330). Orientation power at the construction maximum depends on the CV-selected rank (0.51 above rank 3 vs 0.08 at or below it). These figures are at the isotropic endpoint ρ = 0 (anchor eigengap median 0.111), not at a trending baseline.

## Priority sequence

```text
Close current studies
        |
        v
Audit shape estimator invariance
        |
        v
Audit/correct orientation and shape constructions
        |
        v
Implement orientation-driver attribution
        |
        v
Repeat the medium pilot
        |
        v
Run the paper-grade PLS study
        |
        v
Validate on a real multi-omics case study
        |
        v
Optional graph-native SNF methods
```

## Phase 0 — Close and preserve the current baseline

**Status:** completed on 2026-08-05. The latent-space metric-compatibility and inverse-PLS studies are archived, and their capability specifications are canonical.

**Purpose:** ensure the conclusions already reached remain reproducible and discoverable.

Work:

- Archive the completed latent-space metric-compatibility and inverse-PLS OpenSpec changes.
- Sync their capability specifications into the canonical OpenSpec specs.
- Keep the M-value/PLS pilot report, inverse-study findings, configs, and reproduction commands linked from the documentation.
- Record the current production contract: pooled PLS-DA, M-value methylation, fixed shared coordinate system, and per-statistic interpretation boundaries.

**Exit gate:** passed. A new contributor can locate and reproduce the current magnitude, orientation, shape, and inverse-interpretability evidence from the repository documentation.

## Phase 1 — Resolve shape invariance

**Status:** resolved on 2026-08-04; see [Procrustes Shape Invariance Audit](reports/procrustes-shape-invariance-audit-2026-08-04.md).

**Purpose:** determine whether the current nonzero shape response to rigid rotation is intended methodology or an estimator defect.

Work:

- Create minimal two-dimensional configurations with known translation, uniform scaling, rigid rotation, reflection, and genuine bend transformations.
- Compare MOTCO's GPA/Procrustes output against an analytic expectation and independent reference implementations such as SciPy and the original R procedure.
- Isolate whether any discrepancy comes from centering, centroid-size normalization, rotation/reflection handling, iterative GPA, or the final distance calculation.
- Add permanent invariance and regression tests.
- If necessary, correct the estimator and document any behavioral change; otherwise rename/reinterpret the statistic so its actual invariants are explicit.

**Exit gate:** passed. Rigid translation, scaling, and proper rotation have documented expected outputs enforced by tests. The reflection policy was revisited on 2026-09-03 (`unify-shape-reflection-policy`): reflections are now aligned away at every ambient dimension, and mirror pairs score zero at full rank and above.

## Phase 2 — Audit orientation and shape constructions before integration

**Status:** resolved on 2026-08-05; see [Realized Generator Geometry — Phase 2 Findings](reports/realized-generator-geometry-2026-08-05.md).

**Purpose:** separate generator cross-talk from integration or estimator cross-talk.

Work:

- Measure every generated mode directly in the standardized pre-integration feature space.
- For orientation, verify that the group trajectories differ primarily in normalized direction while path length and within-trajectory configuration remain controlled.
- For shape, distinguish a free biological bend from a magnitude-matched bend; do not require a free bend to be magnitude- or orientation-neutral.
- Retain magnitude as the positive control because it is already specific under M-value integration.
- Add construction-level diagnostics to study outputs so each replicate records the realized pre-integration delta, angle, and shape—not only its requested mode/effect size.
- Revise acceptance targets so they reflect achievable properties of the chosen constructions rather than an assumed perfectly diagonal specificity matrix.

**Exit gate:** passed. Each study mode has a documented realized geometry before PLS, and evaluation records analytic population, standardized population, observed standardized, and PLS-latent checkpoints. Orientation is pure in native methylation but mixed after biological cascade propagation; both shape constructions are documented as free/mixed bends rather than shape-specific interventions.

## Phase 3 — Implement production orientation-driver attribution

**Purpose:** turn a significant global orientation test into stable feature- and pathway-level interpretation.

For two stages, calculate standardized or covariate-adjusted feature transitions

```text
d_A = LSmean(A, stage 2) - LSmean(A, stage 1)
d_B = LSmean(B, stage 2) - LSmean(B, stage 1)
```

and isolate direction from magnitude with

```text
q_observed = d_B / ||d_B|| - d_A / ||d_A||.
```

Work:

- Keep one pooled preprocessing and PLS model fixed for both groups.
- Reconstruct each group-stage PLS displacement through the same fitted loadings and calculate `q_PLS` from the reconstructed normalized feature transitions.
- Report `q_observed`, `q_PLS`, and their residual so users can see which directional difference PLS retained and omitted.
- Attribution reports both decompositions for trajectories with more than two stages: the adjacent transitions, and the signed principal orientations that decompose the tested `angle` estimand, the latter qualified by a relative-eigengap degeneracy flag.
- Bootstrap subjects within the appropriate design strata to estimate feature sign stability, rank stability, and top-k selection frequency.
- Add correlated-feature module and pathway summaries; individual features should not be presented as stable drivers without bootstrap support.
- Report standardized and original-unit effects separately.
- Treat pooled VIP as supporting context only: VIP measures pooled importance for stage prediction, not group differences in normalized progression effects.

**Exit gate:** a significant orientation result can produce a reproducible driver report containing transition-specific features/modules, captured-versus-residual contrast, and uncertainty measures without refitting separate group latent spaces.

## Phase 4 — Repeat the medium pilot

**Purpose:** verify the corrected estimator/constructions before spending cluster time on the definitive study.

Recommended run:

- PLS integration with M-value methylation
- 300 samples and 4 stages
- 50–100 replicates per cell
- 199 permutations
- Effects 0.00, 0.25, 0.50, 0.75, and 1.00
- Modes: magnitude, orientation, shape, and translation
- Matched seeds and construction-level realized-geometry diagnostics

Evaluate:

- Per-statistic Type I error under `none` and translation controls
- Magnitude/delta specificity and monotonic power
- Orientation/angle power and monotonicity
- Shape behavior under the interpretation established in Phase 1
- Off-diagonal rejection rates alongside realized pre-integration geometry
- Stability of PLS component count and attribution summaries

**Exit gate:** no unresolved implementation artifact; magnitude remains specific; orientation reaches an agreed power target with an interpretable curve; shape results match its validated invariance contract. If a target fails, revise the method or scientific claim—not merely the Monte Carlo sample size.

**Status: HOLD** (run 2026-08-27, [full report](reports/phase4-medium-pls-pilot-2026-08-27.md)). The pilot ran
1,900 work units with zero failures using `examples/trajectory_power_study/phase4_pilot_100x199.json`:
100 replicates, 199 permutations, n=300, four stages, matched seeds with a shared zero-effect anchor, and
bounded orientation attribution.

What passed:

- All eighteen predeclared Type I control tests. Translation leaves every statistic at nominal level at
  every effect size; the largest control rate anywhere was 0.09 against a bound of 0.0936.
- Magnitude: delta power 1.00 throughout, with angle 0.01 and shape 0.00 at the top effect — the estimator
  is more specific than the construction, whose population angle is already 39.9°.
- Shape: shape power 0.94–1.00, monotone. The corrected Procrustes estimator detects the constructed bend.
- Diagnostics: selected PLS dimensionality was modally 3 in all 19 cells; all 400 eligible orientation
  replicates produced valid attribution.

What failed:

- Orientation angle power is 0.65 at the top effect against the 0.80 floor — over three Monte Carlo standard
  errors short — and the curve is nearly flat across effects (0.59 → 0.65). **The cause is not yet
  identified.** The construction is strong and survives integration (81° at population, 57° in the latent
  space, against an 8.9° latent null floor), so a simple underpowered-signal account should not apply; and
  the records contradict it directly, since replicates that fail to reject have a *larger* mean observed
  angle (66.7°) than those that reject (52.5°). Latent dimensionality, weak construction, cross-talk, and
  Monte Carlo noise are all ruled out. More replicates would not fix it.

  **Update 2026-08-28** ([report](reports/orientation-sign-anchor-2026-08-28.md)): an orientation sign-anchor
  defect was found and fixed — PC1's sign was anchored on the centered first-stage row, which flips for bent
  trajectories. It was a real bug, but it does **not** explain this failure. On byte-identical data corrected
  power is 0.64/0.64/0.67/0.68 across effects 0.25 → 1.00: ~3 points higher, still flat, still short of 0.80.
  The near-180° null mass is unchanged at 21/700, and the realized-geometry figures stand (81° population,
  57° latent). The shortfall survives correction and remains unexplained.

  **Update 2026-09-01** ([report](reports/angle-null-pivotality-2026-09-01.md)): the shortfall is
  **explained**. Each replicate's own `angle` critical value regresses on its own observed angle with slope
  0.811 — near-unit tracking, specific to `angle` (the same slope is 0.030 for `delta` and 0.058 for `shape`
  in the cells where those carry signal). That accounts for the inversion completely: the 32 non-rejecting
  orientation replicates carry a larger mean observed angle (60.8° vs 46.5°) but a 6.6× larger critical
  value (103.6° vs 15.8°). The tracking is **not** a defect to remove — a fixed threshold drops power from
  0.68 to 0.01, and cross-replicate standardization gains only 0.02 — so the statistic and its test are
  sound, and the gap to 0.80 is a design-point property.

Two findings carry into Phase 5. Orientation→shape is the only response that first becomes material at the
PLS checkpoint rather than in the population geometry. And the PLS reconstruction retains only ~6% of the
observed orientation contrast (cosine 0.08), so its captured-component top-20 precision against generator
truth is 0.15 while the observed component's is 1.00 — driver reports must use the observed component.
Component selection saturates at `n_stages - 1` because supervision is on the stage label, so a space sized
for stage separation is not sized to preserve group orientation contrast.

**Update 2026-09-02** ([note](reports/latent-config-spectrum-2026-09-02.md)): the quantity the pivotality
report could not name — what makes one replicate resolvable and another not — is now **recorded on every
replicate**. Each record carries the eigenspectrum of its centered latent stage-mean configuration (pooled
and per group) and, when permutations run, the pooled relative eigengap over its own permutation draws. The
study report stratifies orientation power by that eigengap and the pivotality analysis reports its
association with the null width. It is a recorded covariate only: no statistic, p-value, permutation draw,
or decision rule changed.

Before Phase 5, see the [Phase 5 readiness worklist](phase5-readiness.md). Its blocking item — diagnosing
the orientation shortfall — was **resolved on 2026-09-01** by
[`diagnose-angle-null-pivotality`](reports/angle-null-pivotality-2026-09-01.md): the co-variation
hypothesis is confirmed, it explains the inversion, no remedy recovers the power, and `angle` proceeds as
specified. What remains: investigate orientation→shape at the PLS checkpoint, revisit latent dimensionality
if orientation is a primary estimand, and then choose the Phase 5 design point — which now also owns the
0.80 orientation power floor.

**Update 2026-09-08** ([report](reports/phase5-design-point-pilot-2026-09-08.md)): the design point is
**chosen** and the floor is **met, not revised**. A crossed pilot (ρ ∈ {0, 0.5, 0.8} × n ∈ {300, 600, 1200},
`p_dmp = 0.1`, uncensored effect axis, 6,500 units, 0 failures) puts orientation `angle` power at **0.88**
(1·SE lower bound 0.85) at **ρ = 0, n = 1200**, with that column's null at nominal level on every statistic.
Sample size does what the pivotality report predicted — the `angle` null width contracts 40.7° → 18.8° →
10.8° along n at a constant eigengap — but baseline continuity does not help and is non-monotone (ρ = 0.5 is
worse at every n), because the orientation surgery realizes a smaller latent contrast on a trending baseline.
The Phase 5 claim is therefore n-conditional at the isotropic endpoint, reported beside the recorded eigengap
(median 0.049). Selected rank stays at `n_stages − 1` for every n, so the retained-rank question passes to
readiness item 3 unchanged; magnitude and shape must be re-measured at the chosen point because `p_dmp`
changed.

The July pilot ([`mvalue-pls-pilot-2026-07-30.md`](reports/mvalue-pls-pilot-2026-07-30.md)) predates the
corrected shape estimator and realized-geometry diagnostics. It is retained unchanged as historical evidence
and is superseded for Phase 4 gate purposes.

## Phase 5 — Paper-grade PLS operating-characteristic study

**Purpose:** produce the definitive evidence for publication and production guidance.

Planned baseline:

- Design point chosen 2026-09-08 ([report](reports/phase5-design-point-pilot-2026-09-08.md)): ρ = 0
  (independent baseline), n = 1200, four stages, `p_dmp = 0.1`, pooled PLS on M-value methylation
- Rank rule committed 2026-09-08 ([latent-rank ladder report](reports/latent-rank-ladder-2026-09-08.md)):
  **stage-supervised double cross-validation** (`plsda_doubleCV`, modal LV, parsimony tie-break; group-blind),
  which selects rank 3 = `n_stages − 1` at the design point. The predeclared ladder rule over fixed ranks
  {3, 4, 6, 9, 12} returned `keep_cv`: no rank gained orientation `angle` power (0.79–0.86 vs 0.85) and every
  rank above 3 lost shape/`shape` (1.00 → 0.57) and orientation `shape`/`delta` power while off-target `angle`
  responses rose. Magnitude and shape were re-measured there under CV (`delta` 1.00, `shape` ≥ 0.95).
- 500 replicates per cell
- 999 permutations per test
- Resumable, sharded execution using `examples/trajectory_power_study/phase5_power_study.json` (`study.json` is the superseded historical config)
- Predeclared acceptance targets and parameter signatures

Deliverables:

- Per-statistic Type I tables and combined-rule Type I results
- Power curves with Monte Carlo uncertainty
- Mode-by-statistic response matrix interpreted against realized input geometry
- Sensitivity analyses for sample size, retained PLS dimensions, signal density, and relevant nuisance settings
- Runtime, failure, and reproducibility report
- Paper-ready figures and a versioned study report tied to the exact code revision and configuration

**Run status:** executed 2026-09-10 as job 880599 on `ing` (19 cells × 500 replicates × 999 permutations = 9,500 units, 100 shards, 3 h 04 m wall, 182.7 core-hours, 0 failures, 0 censored surgeries). Gate decision **HOLD** on the two magnitude mandatory controls; all three mandatory power diagonals and all 21 Type I checks met. See the [findings report](reports/phase5-paper-grade-2026-09-10.md) and `results/phase5-2026-09-10/`. **Magnitude re-measured 2026-10-08** under the `joint` construction as job 906705 (9 cells × 500 × 999 = 4,500 units, 2 h 06 m wall, 106.7 core-hours, 0 failures): reduced gate **PROCEED** — magnitude/`delta` 1.000 at the top with the rise resolved at 0.13 / 0.38 / 0.84 / 0.996 over e = 0.0025–0.02, magnitude/`angle` 0.024 and magnitude/`shape` 0.000 at e = 1.00, anchor byte-identical to the 2026-09-10 anchor at all 500 replicates. See the [addendum](reports/phase5-magnitude-remeasurement-2026-10-08.md), `results/phase5-magnitude-2026-10-08/`, the pilot `results/phase5-magnitude-pilot-2026-10-08/`, and the bracket `results/magnitude-axis-bracket-2026-10-08/`.

**Exit gate:** conclusions are stable at paper-grade Monte Carlo precision, all deviations from preregistered targets are explained, and the report distinguishes statistical operating characteristics from biological construction cross-talk. **CLOSED 2026-09-10 as a conditional pass** ([exit review](reports/phase5-exit-review-2026-09-10.md)): all three criteria met; the Phase 4 gate's HOLD stands as recorded and resolves to a construction defect with a known fix; the magnitude specificity claim, withheld pending re-measurement of the corrected construction, **was lifted on 2026-10-08 for the `joint` construction** ([addendum](reports/phase5-magnitude-remeasurement-2026-10-08.md)) — stated as "specific and conservative", and construction-specific until the default is flipped.

## Phase 6 — Real-data case study

**Purpose:** demonstrate that the validated method answers a meaningful biological question end to end.

Work:

- Predefine cohorts, stages, covariates, missing-data handling, and preprocessing.
- Fit one pooled, cross-validated PLS representation.
- Report magnitude, orientation, and validated shape results with uncertainty.
- For significant orientation, produce observed and PLS-reconstructed feature contrasts, bootstrap stability, and pathway/module enrichment.
- Include sensitivity to PLS dimensionality and preprocessing choices.
- Make the analysis reproducible from raw-input contracts to final figures.

**Exit gate:** collaborators can trace each global trajectory claim to stable feature/module drivers and reproduce the result from documented inputs and configuration.

**Case study.** SEA-AD middle temporal gyrus astrocytes, with donor-level RNA + ATAC pseudobulks. Sex is the group and ADNC is merged into three stages (Not AD + Low / Intermediate / High). The primary cohort is the ≥30-nuclei cohort (n = 80: F 11/10/28, M 9/10/12); the ≥50-nuclei cohort (n = 74) is a sensitivity analysis. SEA-AD preprocessing and the real-data analysis live in the external case-study repository, [Omnia-Lab-UC/sead_sexdiff](https://github.com/Omnia-Lab-UC/sead_sexdiff). MOTCO supplies the method and the operating characteristics below. The case-study analysis must pin MOTCO to the revision that merges `phase6-small-n-operating-study`. The operating study itself ran at `5325f6d` (records) and `946e471` (report), with motco 0.6.0, numpy 2.3.5, scikit-learn 1.8.0, scipy 1.16.3 and python 3.11.16.

**Small-n operating study — done 2026-10-09** ([report](reports/phase6-small-n-2026-10-09.md)). This is the precondition the Phase 5 exit review set before any SEA-AD result can be interpreted. MOTCO gained explicit group × stage sample sizes (`group_stage_sizes`) and evaluation-time block selection (`integration_params.layers`); both are byte-identical when absent.
- A pilot (`results/phase6-small-n-pilot-2026-10-09/`, 3,500 units) fixed the effect axes.
- The paper-grade run was split into a magnitude profile and an orientation/shape/translation profile, 500 × 999 each, over a four-column design grid (cohort size × block count). That is 42,000 units, 150 core-hours and 0 failures (`results/phase6-small-n-magnitude-2026-10-09/`, `results/phase6-small-n-2026-10-09/`).
- Advisory gates: magnitude **PROCEED**; orientation/shape **HOLD**, because both power floors are missed (0.300 / 0.228). Per design D6 this revises the claim, not the method.

The report's binding interpretability statement for the case study:
- **`delta`.** A rejection is evidence that the trajectories differ; it is consistent with a size difference but does not exclude a direction or shape difference. A non-rejection argues only against size differences of about 25 % or more.
- **`angle` and `shape`.** A rejection is evidence that the trajectories differ, but not of which kind. A non-rejection is uninformative at n = 80.
- **Eigengap.** Orientation results must still be read against the cohort's own recorded eigengap. The simulated low tercile (below about 0.13) had at most 0.114 power.

## Phase 7 — Optional SNF-native trajectory methodology

**Purpose:** use SNF according to graph geometry rather than forcing Euclidean MOTCO statistics onto spectral coordinates.

This is not on the critical path for the PLS production workflow.

**Status: deferred, and the question is closed rather than pending** ([Phase 5 exit review](reports/phase5-exit-review-2026-09-10.md), §5). Nothing in the paper-grade run or the magnitude-construction diagnostic bears on SNF — the study measured the PLS latent space throughout and `concat` appeared only as a baseline — so the deferral does not await further Phase 5 evidence. It awaits a decision to define graph-native path-length and angle statistics, which is independent work.

Candidate work:

- Define group-stage states on one pooled fused graph.
- Evaluate diffusion or resistance path length as a magnitude analogue.
- Compare transition profiles or diffusion flows as an orientation analogue.
- Compare normalized stage-by-stage diffusion-distance matrices as a configuration analogue.
- Design a permutation scheme that rebuilds or conditions on the graph without invalidating exchangeability.
- Validate every proposed graph statistic on controlled nonlinear manifolds before exposing it as production functionality.

**Exit gate:** an SNF-specific statistic has a clear graph-native estimand, calibrated Type I error, controlled power behavior, and no claim of preserving original feature-space Euclidean geometry.

## Cross-cutting engineering work

These items support every scientific phase:

- Preserve fitted preprocessing, PLS model, component count, feature order, and loadings as a versioned analysis artifact.
- Record seeds, configurations, parameter signatures, software versions, and failed replicates in every study.
- Keep fast deterministic invariance tests separate from slow operating-characteristic studies.
- Add CI coverage for core geometry invariants and study-config validation.
- Prefer structured CSV/JSON outputs plus concise Markdown reports over findings that exist only in notebooks.
- Keep public documentation synchronized with methodological decisions and archived OpenSpec findings.
- **Decide the fate of the committed toy fixture** (`examples/data/toy/`, deferred from `docs-sync-audit`, 2026-09-04). The fixture predates the numpy generator: its `truth.json` records `intersim_metadata` and omits `n_stages`, `p_dmp`, and `baseline_continuity`, and running the current `motco simulate` at the same seed produces a different dataset. `docs-sync-audit` corrected the documentation — the `toy-dataset` spec now describes the fixture as pinned rather than regenerable — but did not resolve the underlying question: regenerate the example data with the numpy generator so the tutorial exercises the shipped code path, or keep it pinned as a historical artifact. Regenerating touches committed example data, `examples/motco_example.ipynb`, the `toy-dataset` spec scenarios, and `tests/test_toy_dataset.py`, so it needs its own change.

## Explicitly deferred

- Expanding SNF use before graph-native metrics are defined
- Reading the pooled orientation power figure as applicable at any eigengap
- Claiming `angle` or `shape` specificity against the methylation-only `all` magnitude construction (its controls fail; the claim holds for `joint` only — [addendum](reports/phase5-magnitude-remeasurement-2026-10-08.md))
- Interpreting pooled VIP scores as orientation drivers
- Fitting separate, unaligned PLS spaces per group
- Claiming that a latent-space reconstruction is a unique or causal molecular inverse
- Optimizing performance before the statistic and construction contracts are settled

## Next three changes

The `magnitude_kind` default flip, formerly item 1, is done (`default-joint-magnitude-construction`, 2026-10-09): `joint` is the generator and `motco simulate` default, and the ten historical configs pin `"all"` with their parameter signatures checked against a pre-flip fixture (`tests/test_magnitude_default_pins.py`).

1. **Run the Phase 6 real-data case study** in [Omnia-Lab-UC/sead_sexdiff](https://github.com/Omnia-Lab-UC/sead_sexdiff), pinned to the MOTCO revision that merges `phase6-small-n-operating-study`. The small-n operating study is done ([report](reports/phase6-small-n-2026-10-09.md)), and its per-statistic interpretability statement is binding. Report the cohort's **own recorded eigengap** beside any orientation result, never the Phase 5 pooled 0.840. Treat `angle` and `shape` non-rejections as uninformative, and treat a `delta` rejection as evidence of a difference, not specifically of size. Carry the structural-analogue limitation: InterSIM blocks, 498 measured features.
2. **Declare the acceptance specificity target one-sided.** The two-sided form flags conservative controls as failures ([addendum](reports/phase5-magnitude-remeasurement-2026-10-08.md), §3). The one-sided bound mirrors the gate's α-based SE; it is an opt-in key defaulting to two-sided, and is omitted from `observations` when unset so committed reports regenerate byte-identically.
3. **Consider defaulting the null-dispersion materiality rule.** Opt-in today so the committed Phase 5 `report/` keeps regenerating byte-identically; switching the default is a reporting change needing a recorded re-issue of affected localization output.
