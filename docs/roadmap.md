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
- The *former default* magnitude surgery (`magnitude_kind='all'`) is not size-pure: at paper-grade precision it moves `angle` (0.276) and `shape` (0.742) far past the α + 2·SE control bound, which held the Phase 4 gate at **HOLD** ([paper-grade report](reports/phase5-paper-grade-2026-09-10.md), §6). **The mechanism is established** ([magnitude-construction diagnostic](reports/magnitude-construction-diagnostic-2026-09-10.md)): it scales `delta_methyl` alone, so the surgery is exactly size-pure *within* every omic block and the entire off-target response is the geometry of concatenating one block that grew with two that did not. The remedy — the `joint` construction, every omic's δ scaled together — is now a production `magnitude_kind` and **has been measured at paper grade** ([magnitude re-measurement addendum](reports/phase5-magnitude-remeasurement-2026-10-08.md)): both controls met (0.024 / 0.000), gate PROCEED. The acceptance-target reader records the same two controls as *not met*, because its default two-sided rule treats a conservative rate as a failure. `one-sided-acceptance-targets` (2026-10-14) adds the opt-in `bound: "upper"` form, which agrees with the gate. Since `default-joint-magnitude-construction`, `joint` is the generator and `motco simulate` default; the committed pre-flip study configs and diagnostics pin `all` explicitly so their signatures and datasets are unchanged. Specificity remains unestablished for `all` (and for the methylation-only `extremes` and `shape_kind='magnitude'` scalings), so the paper must still say which construction a magnitude result is about.
- Orientation `angle` power does not rise with effect: it saturates at 0.84–0.88 from e = 0.25 because the `angle` null widens with the observed statistic, so the strictly-monotone form of the power claim is unsupported (the 0.80 floor is still cleared).
- At the Phase 6 case-study design (n = 80, three stages, two blocks) only `delta` has useful power, and none of the three statistics is mode-specific for orientation or shape ([small-n report](reports/phase6-small-n-2026-10-09.md)). Type I holds (0.010–0.040). magnitude/`delta` reaches 0.542 at a 10 % joint size change and 0.956 at 25 %. orientation/`angle` and shape/`shape` reach only 0.300 and 0.228 at the strongest constructions the generator allows, and those constructions are detected more often by other statistics (orientation by `shape` 0.464 and `delta` 0.354; shape by `delta` 0.330). Orientation power at the construction maximum depends on the CV-selected rank (0.51 above rank 3 vs 0.08 at or below it). These figures are at the isotropic endpoint ρ = 0 (anchor eigengap median 0.111), not at a trending baseline. The [orientation measurement probe](reports/orientation-measurement-probe-2026-10-12.md) located the orientation loss in the measurement, not the generator. The net-displacement estimator raises Phase 6 orientation power in the production PLS space from 0.29 to 0.67 and removes the Phase 5 shape → `angle` leak (0.89 → 0.02), and the [net-displacement Type I study](reports/net-displacement-type-i-2026-10-13.md) shows it holds size at paper grade (anchor `angle` 0.050 at both designs, bound 0.0695; the probe's 0.08–0.10 was Monte Carlo noise at 100 replicates). The [net-displacement power comparison](reports/net-displacement-power-2026-10-14.md) establishes its advantage at paper grade: orientation/`angle` 0.300 → 0.648 at Phase 6 and 0.840 → 0.924 at Phase 5 on identical datasets, the shape → `angle` leak removed (0.138 → 0.054; 0.796 → 0.022), every control within the bound. Its verdict was `adopt_net_displacement`, and since `default-net-displacement-orientation` (2026-10-15) net displacement is the harness default; the committed configs and diagnostics pin `"pc1"` so their results are unchanged. At n = 80 a non-rejection of `angle` stays uninformative.

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
- Sensitivity analyses for sample size, retained PLS dimensions, signal density, and relevant nuisance settings *(sample size, stage count, block count and signal density: measured at pilot grade by the [design-guidance pilot](reports/design-guidance-pilot-2026-10-28.md), 2026-10-28, which proposes the paper-grade slice; retained dimensions: the latent-rank ladder above)*
- Runtime, failure, and reproducibility report
- Paper-ready figures and a versioned study report tied to the exact code revision and configuration

**Run status:** executed 2026-09-10 as job 880599 on `ing` (19 cells × 500 replicates × 999 permutations = 9,500 units, 100 shards, 3 h 04 m wall, 182.7 core-hours, 0 failures, 0 censored surgeries). Gate decision **HOLD** on the two magnitude mandatory controls; all three mandatory power diagonals and all 21 Type I checks met. See the [findings report](reports/phase5-paper-grade-2026-09-10.md) and `results/phase5-2026-09-10/`. **Magnitude re-measured 2026-10-08** under the `joint` construction as job 906705 (9 cells × 500 × 999 = 4,500 units, 2 h 06 m wall, 106.7 core-hours, 0 failures): reduced gate **PROCEED** — magnitude/`delta` 1.000 at the top with the rise resolved at 0.13 / 0.38 / 0.84 / 0.996 over e = 0.0025–0.02, magnitude/`angle` 0.024 and magnitude/`shape` 0.000 at e = 1.00, anchor byte-identical to the 2026-09-10 anchor at all 500 replicates. See the [addendum](reports/phase5-magnitude-remeasurement-2026-10-08.md), `results/phase5-magnitude-2026-10-08/`, the pilot `results/phase5-magnitude-pilot-2026-10-08/`, and the bracket `results/magnitude-axis-bracket-2026-10-08/`.

**Exit gate:** conclusions are stable at paper-grade Monte Carlo precision, all deviations from preregistered targets are explained, and the report distinguishes statistical operating characteristics from biological construction cross-talk. **CLOSED 2026-09-10 as a conditional pass** ([exit review](reports/phase5-exit-review-2026-09-10.md)): all three criteria met; the Phase 4 gate's HOLD stands as recorded and resolves to a construction defect with a known fix; the magnitude specificity claim, withheld pending re-measurement of the corrected construction, **was lifted on 2026-10-08 for the `joint` construction** ([addendum](reports/phase5-magnitude-remeasurement-2026-10-08.md)) — stated as "specific and conservative", and construction-specific until the default is flipped.

## Phase 6 — Real-data case study

**Purpose:** demonstrate that the validated method answers a meaningful biological question end to end.

Work:

- Predefine cohorts, stages, covariates, missing-data handling, and preprocessing. *(Case-study repository. `motco analyze` takes the stage order and block units as explicit inputs and does no normalization, filtering or imputation; donor covariates and the missing-covariate policy are explicit inputs since 0.8.0, and covariate-aware latent construction (`--latent-adjustment`) is opt-in since 0.9.0.)*
- Fit one pooled, cross-validated PLS representation. *(MOTCO: `motco analyze`, `real-data-trajectory-analysis`.)*
- Report magnitude, orientation, and validated shape results with uncertainty. *(MOTCO: `motco analyze` statistics, RRPP p-values and null summary.)*
- For significant orientation, produce observed and PLS-reconstructed feature contrasts, bootstrap stability, and pathway/module enrichment. *(MOTCO: `motco analyze --attribution` for contrasts and bootstrap stability; pathway/module enrichment stays in the case-study repository.)*
- Include sensitivity to PLS dimensionality and preprocessing choices. *(Dimensionality: MOTCO, `motco analyze --forced-components`. Preprocessing choices: case-study repository.)*
- Make the analysis reproducible from raw-input contracts to final figures. *(MOTCO writes the analysis artifact; raw inputs and figures live in the case-study repository.)*

**Exit gate:** collaborators can trace each global trajectory claim to stable feature/module drivers and reproduce the result from documented inputs and configuration.

**Case study.** SEA-AD middle temporal gyrus astrocytes, with donor-level RNA + ATAC pseudobulks. Sex is the group and ADNC is merged into three stages (Not AD + Low / Intermediate / High). The primary cohort is the ≥30-nuclei cohort (n = 80: F 11/10/28, M 9/10/12); the ≥50-nuclei cohort (n = 74) is a sensitivity analysis. SEA-AD preprocessing and the real-data analysis live in the external case-study repository, [Omnia-Lab-UC/sead_sexdiff](https://github.com/Omnia-Lab-UC/sead_sexdiff). MOTCO supplies the method and the operating characteristics below. The case-study analysis must pin MOTCO to the revision that merges `phase6-small-n-operating-study`. The operating study itself ran at `5325f6d` (records) and `946e471` (report), with motco 0.6.0, numpy 2.3.5, scikit-learn 1.8.0, scipy 1.16.3 and python 3.11.16.

**Small-n operating study — done 2026-10-09** ([report](reports/phase6-small-n-2026-10-09.md)). This is the precondition the Phase 5 exit review set before any SEA-AD result can be interpreted. MOTCO gained explicit group × stage sample sizes (`group_stage_sizes`) and evaluation-time block selection (`integration_params.layers`); both are byte-identical when absent.
- A pilot (`results/phase6-small-n-pilot-2026-10-09/`, 3,500 units) fixed the effect axes.
- The paper-grade run was split into a magnitude profile and an orientation/shape/translation profile, 500 × 999 each, over a four-column design grid (cohort size × block count). That is 42,000 units, 150 core-hours and 0 failures (`results/phase6-small-n-magnitude-2026-10-09/`, `results/phase6-small-n-2026-10-09/`).
- Advisory gates: magnitude **PROCEED**; orientation/shape **HOLD**, because both power floors are missed (0.300 / 0.228). Per design D6 this revises the claim, not the method.

**Orientation measurement probe — done 2026-10-12** ([report](reports/orientation-measurement-probe-2026-10-12.md), `results/orientation-measurement-probe-2026-10-12/`). This probe asked why orientation is weak at n = 80. It crossed five measurement spaces with two orientation estimators on matched seeds, at the Phase 6 and Phase 5 designs: 1,000 units, 10 arms each, 28 core-hours, 0 failures. MOTCO gained the opt-in `orientation_estimator` (`net_displacement`) and the diagnostic `pca` / `pls_within_group` integration methods, all byte-identical when unused.
- Under the predeclared rule no arm advances.
- `pls_cv × net_displacement` fails only the anchor criterion (Phase 6 0.100 vs bound 0.0936), with a large orientation gain and no power or specificity loss.
- A fixed-basis `pls_within_group` is invalid without a refit (anchor up to 0.88).
- `pca` and forced rank 2(K − 1) lose Phase 6 shape power.
- A baseline-continuity (ρ) sweep at the Phase 6 design was considered first and dropped: the Phase 5 design-point pilot already showed that ρ does not buy orientation power, and the probe places the loss in the measurement.
- Production defaults are unchanged, and the case-study interpretability statement above stands.

**Net-displacement Type I study — done 2026-10-13** ([report](reports/net-displacement-type-i-2026-10-13.md), `results/net-displacement-type-i-2026-10-13/`). The paper-grade estimator comparison's (Next changes item 2) anchor-only step: the zero-effect anchor and the translation control at 500 × 999, at both designs, across five matched columns (`pls_cv × pc1`, and `pls_cv`, forced rank 2(K − 1), `pca`, `concat` × net displacement). 12,000 units, 133 core-hours, 0 failures. MOTCO gained the explicit design-grid form (`design_grid.points`), the advisory `acceptance.anchor_decision` rule and `scripts/net_displacement_follow_up.py`, all byte-identical when unused.
- Every net-displacement column holds the bound 0.0695 on `delta`, `angle` and `shape` at both designs; the lead arm's anchor `angle` is 0.050 at both.
- The follow-up rule selects **`paper_grade_power_comparison`**; the anchor-excess diagnosis is not needed.
- Under net displacement the translation control's `angle` rate is close to the bound (lead arm 0.068 at Phase 6; `concat` 0.078), so the comparison must hold it to that bound.
- The baseline column reproduces the Phase 6 and Phase 5 parent records (same seeds and RRPP seed): p-values identical in 2,000/2,000 records.

**Net-displacement power comparison — done 2026-10-14** ([report](reports/net-displacement-power-2026-10-14.md), `results/net-displacement-power-2026-10-14/`). Formerly Next changes item 2: `pls_cv × net_displacement` against `pls_cv × pc1` at 500 × 999, four modes at e = 0.25–1.0, at both designs, on matched datasets and permutation draws. 36,000 units, 384 core-hours, 0 failures. MOTCO gained the advisory `acceptance.column_decision` rule, `scripts/net_displacement_power_follow_up.py` and `scripts/net_displacement_power_checks.py`, all byte-identical when unused.
- The candidate qualifies at both designs: orientation/`angle` at e = 1.0 rises from 0.300 to 0.648 at Phase 6 (`gain`, +0.348 against a 0.059 threshold) and from 0.840 to 0.924 at Phase 5 (`no_loss`). Anchor `angle` 0.050 at both.
- Every mandatory control is within 0.0695, which PC1 does not achieve: shape → `angle` 0.054 / 0.022 (PC1 0.138 / 0.796), magnitude → `angle` 0.030 / 0.000, translation → `angle` 0.068 / 0.040 (Phase 6 one rejection under the bound).
- The cross-design rule selects **`adopt_net_displacement`**. Its follow-up, the harness-default flip, is done (`default-net-displacement-orientation`, 2026-10-15).
- At Phase 6 the eigengap no longer gates orientation power (low tercile 0.114 → 0.731), but the CV-selected rank still does at the construction maximum (0.280 at rank ≤ 3, 0.996 above), and power is 0.084–0.178 at e ≤ 0.75.
- `delta` and `shape` p-values are identical between the two columns in 17,000/17,000 paired records, and the baseline column reproduces 15,500 parent records with 0 seed or p-value mismatches.

The report's binding interpretability statement for the case study:
- **`delta`.** A rejection is evidence that the trajectories differ; it is consistent with a size difference but does not exclude a direction or shape difference. A non-rejection argues only against size differences of about 25 % or more.
- **`angle` and `shape`.** A rejection is evidence that the trajectories differ, but not of which kind. A non-rejection is uninformative at n = 80.
- **`angle` (net displacement, revised 2026-10-15).** With net displacement the harness default (`default-net-displacement-orientation`), the `angle` statement is the one in the [power comparison report](reports/net-displacement-power-2026-10-14.md) §7. `angle` is the angle between the groups' net displacement vectors (last stage mean − first stage mean) in the latent space, so it reads only the endpoint stages. A rejection is evidence that the groups' net directions of change differ: Type I is controlled (anchor 0.050), and the size, shape and location constructions stay within the bound (0.024–0.068). Translation sits at 0.058–0.068, near the bound, so a rejection alongside a location difference warrants care. A non-rejection is uninformative at n = 80: power is 0.084–0.178 at e ≤ 0.75, and 0.648 at the construction maximum (0.280 when CV retains three or fewer components, 0.996 above).
- **Selected rank.** The cohort's CV-selected PLS rank replaces its eigengap as the covariate reported beside any `angle` result; under net displacement the eigengap no longer conditions the statement.

**Covariate-aware latent construction — done 2026-10-20** (`covariate-aware-latent-construction`, MOTCO 0.9.0; [API](api/analysis.md#latent-adjustment)). In both SEA-AD cohorts CV selected rank 1, so `angle` was 0° by construction. The single component turned out to be a donor-level library-quality axis: Spearman −0.84 with the astrocyte nuclei's ribosomal-protein UMI fraction, against +0.49 with stage. The 0.8.0 covariates adjust the stage means exactly within the basis but never touch which basis is chosen. Label-blind residualization upstream removed the axis but also the stage signal collinear with it (held-out CV AUROC 0.63 → 0.41). `covariate_latent_adjustment` (`stage` | `group_stage`) instead residualizes the features on the declared covariates conditional on the stage, or group × stage, indicators, with per-training-set coefficients inside the CV, folds stratified on the conditioning cells, and the pooled scale kept. It measures the raw features in the new basis, so only the basis and the CV rank change. Without the setting every output equals 0.8.0. The local smoke run of the Type I check (R = 20, 99 permutations) calibrated the nuisance injection: the unadjusted basis's first component tracks the covariate at median |ρ| 0.815 with `s_z = 40`, against 0.64 at the covariate check's `s_z = 3`. In the residualized arms that falls to about 0.32. The paper-grade Type I check (2026-10-20, [report](reports/latent-adjustment-type-i-2026-10-20.md)) returned **`adopt_stage`**: `stage` held `delta` 0.031, `angle` 0.067 and `shape` 0.058 against 0.0706; `group_stage` exceeded on `angle` (0.082). The power comparison (2026-10-20, [report](reports/latent-adjustment-power-2026-10-20.md)) returned `candidate_fails` for both conditionings. `stage` gains orientation power (`angle` 0.922 against 0.784 with the nuisance axis, 0.931 against 0.642 without), but loses `shape` power with the axis (0.243 against 0.428) and small-magnitude power without it. Its translation `angle` control, like the covariate-blind basis's, sits above the bound. `group_stage` loses most orientation and shape power. The case study may therefore use `--latent-adjustment stage` only, as a sensitivity analysis beside the 0.8.0 covariate-adjusted analysis, under both reports' conditions. It must first justify the covariate set in its own proposal (the mediator limitation).

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

- Preserve fitted preprocessing, PLS model, component count, feature order, and loadings as a versioned analysis artifact. For real data, `motco.analysis.write_analysis_outputs` (`motco analyze --out-dir`) writes `result.json`, latent scores, PLS loadings and latent stage means. It records units, the feature order signature and software versions, and pickles no estimator.
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

## Next changes

The `magnitude_kind` default flip, formerly item 1, is done (`default-joint-magnitude-construction`, 2026-10-09): `joint` is the generator and `motco simulate` default, and the ten historical configs pin `"all"` with their parameter signatures checked against a pre-flip fixture (`tests/test_magnitude_default_pins.py`).

The one-sided acceptance bound, formerly item 2, is done (`one-sided-acceptance-targets`, 2026-10-14). Each `acceptance.type_i` and `acceptance.specificity` entry takes `bound`: `two_sided` (the default, unchanged) or `upper`, which is the gate's `rate ≤ α + k·√(α(1 − α)/R)` and shares its implementation. No committed config opts in, so every committed `report/` regenerates byte-identically.

The paper-grade estimator comparison, formerly item 2, is done (`net-displacement-power-comparison`, 2026-10-14; [report](reports/net-displacement-power-2026-10-14.md)) and selected `adopt_net_displacement`. It ran no `concat` column, so `concat`'s p-padding question stays open; `concat` remains the ceiling reference only.

The orientation-estimator default flip, formerly item 2, is done (`default-net-displacement-orientation`, 2026-10-15): `SimulationEvaluationParams.orientation_estimator` defaults to `net_displacement`, every committed config and historical diagnostic pins `"pc1"`, and the configs' signatures and matched seeds are checked against a pre-flip fixture (`tests/test_orientation_default_pins.py`). The stats primitives (`estimate_difference`, `RRPP`) keep PC1 as their default. No compute was needed.

1. **Run the Phase 6 real-data case study** in [Omnia-Lab-UC/sead_sexdiff](https://github.com/Omnia-Lab-UC/sead_sexdiff). The orientation-estimator default flip unblocked it. The `real-data-trajectory-analysis` change (MOTCO 0.7.0) removed the last MOTCO-side blocker: the case study must run its trajectory analysis through **`motco analyze`** (or `motco.analysis.analyze_trajectories`; [API](api/analysis.md)) and pin the MOTCO release that contains it. That entry point is the harness's production measurement on caller-supplied blocks. A regression test checks it against the harness on identical data. It defaults to the net-displacement estimator, requires an explicit stage order and records the selected rank, null summary and spectrum, so no estimator argument has to be passed by hand. Do not assemble the analysis from `motco plsr` + `motco de` or from the stats primitives: they are a different pipeline, and the primitives and `motco de` default to PC1. The small-n operating study ([report](reports/phase6-small-n-2026-10-09.md)) and the revised `angle` statement (Phase 6 above, from the [power comparison report](reports/net-displacement-power-2026-10-14.md) §7) are binding. Report the cohort's **CV-selected rank** beside any `angle` result as the covariate; the record carries it. Treat `angle` and `shape` non-rejections as uninformative, and treat a `delta` rejection as evidence of a difference, not specifically of size. Carry the structural-analogue limitation (InterSIM blocks, 498 measured features). Covariate adjustment landed in MOTCO 0.8.0 (item 2 below, done) with a `covariates_validated` Type I decision, so the case study may declare donor covariates in `motco analyze` under the conditions of the [covariate Type I report](reports/covariate-type-i-2026-10-16.md) §5, and must pin MOTCO ≥ 0.8.0 to do so. Its preprocessing and export (`mtg-astrocyte-motco-export` in sead_sexdiff) produce TMM log2-CPM blocks without chrX/Y/M, the top 5,000 variable features per block, and a donor covariate table: age at death (with a `90+` flag), PMI, RIN, years and level of education, and astrocyte nuclei counts per omic.
2. **Covariate-adjusted trajectories in `motco analyze`** — *done* (`covariate-adjusted-trajectories`, MOTCO 0.8.0, 2026-10-16). `TrajectoryDataset` takes `covariate_cols`, `categorical_covariates` and `missing_covariates`; `motco analyze` takes `--covariate`, `--categorical-covariate` and `--missing-covariates` ([API](api/analysis.md#covariates)). The Phase 6 Type I check returned `covariates_validated`: adjusted `delta` 0.040, `angle` 0.062, `shape` 0.034 against the bound 0.0695, while the unadjusted arm rejected `shape` at 0.364 ([report](reports/covariate-type-i-2026-10-16.md)). Without covariates every output is byte-identical to 0.7.0. The original scope follows. The case study needs donor covariates inside the trajectory regression. MOTCO should own this because covariates enter the model: the LS-mean fit and the RRPP reduced model both change, and so does what is permuted. Scope:
   - covariates in the full and reduced models;
   - stage-mean LS-means evaluated at fixed covariate values (e.g. the pooled mean), recorded;
   - an explicit missing-covariate policy, either error or complete-case, chosen by the caller, recorded, and never imputed;
   - the covariate design recorded in `result.json`;
   - a Type I check at the Phase 6 design with covariates that are correlated with group and stage.
   Absent covariates must stay byte-identical to 0.7.0.

   **Scope boundary for MOTCO preprocessing** (decided 2026-10-16). The test is whether a step touches the model or the permutation null.
   - *Normalization stays upstream.* It is omic-specific, and when it is blind to the labels it changes what is measured, not the validity of the test. It is covered by the case study's preprocessing sensitivity analysis.
   - *Label-blind feature filtering stays upstream.* Supervised filtering (by group or stage) would bias the null unless it were repeated inside every permutation, so MOTCO documents it as forbidden rather than offering it.
   - *Imputation stays upstream.* The finite-values contract stays. Missing covariates are the only exception, and they are handled by the explicit policy above.
3. **Covariate-aware latent construction: validation runs** (`covariate-aware-latent-construction`, MOTCO 0.9.0). (a) The Type I check, `scripts/latent_adjustment_type_i_check.py` — *done* 2026-10-20 ([report](reports/latent-adjustment-type-i-2026-10-20.md), `results/latent-adjustment-type-i-2026-10-20/`): 500 × 999, 5.1 core-hours, **`adopt_stage`**; `group_stage` is not admissible (`angle` 0.082 against 0.0703). (b) The power comparison — *done* 2026-10-20 ([report](reports/latent-adjustment-power-2026-10-20.md), `results/latent-adjustment-power-2026-10-20/`), 158 core-hours, advisory `candidate_fails` for both arms — `scripts/latent_adjustment_power_comparison.py`: `nuisance` and `clean` columns × anchor, translation and three modes at four effects, 500 × 999 (14,000 units of three analyses, about 350 core-hours from the smoke timing), writing `latent_adjustment_power_decision.json` with the per-arm rank distributions. Each gets a dated report under `docs/reports/` carrying the mediator limitation. In the Type I run, 52 `stage` and 38 `group_stage` replicates lost the arm to a rare categorical level dropping out of a CV training set; such replicates are recorded as unavailable.
4. **The rank-1 `angle` degeneracy** — *done* (`minimum-latent-rank`, MOTCO 0.10.0, 2026-10-23; [report](reports/rank-floor-2026-10-23.md), `results/rank-floor-2026-10-23/`). CV could select rank 1, where `angle` is 0° or 180° by construction and its RRPP null is a two-point distribution, while `forced_components` required ≥ 2; both SEA-AD cohorts selected rank 1. The production latent space (`fit_pls_latent_space`, shared by the harness and `motco analyze`) now selects the rank from candidates at or above a **floor of 2** (`min_components=2`; `1` restores the 0.9.0 rule and is the hash-omitted spelling, so every committed config, fixture and dated report reproduces under its `min_components: 1` pin). The 0.9.0 power comparison's Phase 6 construction never selected rank 1, so the check weakened the per-omic stage effect (`delta_* = 0.05` and `0.1` against the Phase 6 value 2.0, `p_dmp` unchanged) until the unfloored rule selected rank 1 in 58 % and 42 % of null replicates (held-out CV AUROC 0.52 and 0.58, the SEA-AD cohorts' 0.63–0.66 regime), and ran two matched arms (`min_components` 1 and 2) per unit at 500 × 999 over 14 cells at each strength (two 7,000-unit arrays, 99 core-hours). Both runs returned **`adopt_floor`**: the floor-2 arm holds the 0.0695 bound at the anchor (`delta` 0.042 / 0.044, `angle` 0.054 / 0.060, `shape` 0.036 / 0.038) and at the translation control (0.036 / 0.036, 0.046 / 0.042, 0.026 / 0.030). The floor-1 arm's rank-1 units had `angle` exactly 0° or 180° with a two-point null that never rejected (anchor `angle` 0.030 / 0.038), which is the degeneracy removed. Power is at the null in both arms at these strengths (the construction has almost no stage signal by design), so the floor defines `angle` and `shape`; it does not make them powerful, and the record's `config_spectrum` eigengap and `null_summary` say how resolvable orientation is in each analysis. sead_sexdiff reruns both cohorts at the floor in its own proposal and pins 0.10.0.
5. **Consider optional block weighting.** `motco analyze` z-scores each feature and concatenates the blocks unweighted, so a block's influence grows with its feature count. The SEA-AD export works around this by forcing equal N (5,000 + 5,000). How blocks are weighted is a property of the joint measurement rather than of one omic, so it is the one preprocessing step that belongs in MOTCO. Candidates are multiblock-PLS-style 1/√p_b scaling or unit total block variance. Weighting changes the latent space, so it is opt-in, byte-identical when absent, and needs its own Type I and power check before the case study may use it. It is separate from item 2.
6. **Consider defaulting the null-dispersion materiality rule.** Opt-in today so the committed Phase 5 `report/` keeps regenerating byte-identically; switching the default is a reporting change needing a recorded re-issue of affected localization output.
7. **Design-guidance pilot** — *done* (`design-guidance-grid`, MOTCO 0.11.0, 2026-10-28; [report](reports/design-guidance-pilot-2026-10-28.md), `results/design-guidance-pilot-2026-10-28/`). The originally planned Table 1 grid is replaced by four pilot profiles under the production measurement (PLS with the CV rank floored at 2, net displacement, joint magnitude; 100 × 199): a crossed core n {80, 200, 500, 1000} × stages {3, 4} × blocks {2, 3} and one-factor arms (group ratio, stage proportions, `p_dmp`, joint stage effect δ, continuity ρ) off n = 200 / 3 stages / 2 blocks; the report maps every Table 1 factor to run, deferred or dropped. The run on `ing` (jobs 929160–929163, 47,200 units, 0 failures, 0 censored surgeries, 345 core-hours) found **`delta` well powered** (0.23 / 0.59 at e = 0.05 / 0.1 at n = 80 with two blocks, 0.72 / 0.94 at n = 200, saturated at e = 0.25 everywhere; a third block raises every rate, a fourth stage changes little), **`shape` monotone** in n, e and block count (0.55 at n = 200 and 0.85 at n = 500 at e = 1.0 with two blocks; ≤ 0.20 at n = 80), and **`angle` governed by n and block count rather than the effect size** (0.05–0.22 at n ≤ 200 with two blocks, 0.41–0.62 at n = 1000; 0.35–0.46 at n = 200 with three blocks), with a flat or non-monotone response in e whose e = 1.0 rate tracks the CV-selected rank. The held-out stage AUROC is 1.000 at every column except the δ arm, where δ = 0.05 / 0.1 give 0.553 / 0.729 — bracketing the SEA-AD cohorts — and no statistic has power at any effect at n = 200. Size holds at 26 of 28 anchors; the two excesses are both 0.10 against 0.0936, the count chance predicts over 84 tests, and the baseline's independent Type I-phase draw gives `angle` 0.03. Follow-ups: (a) **the paper-grade slice** the report proposes — the three-stage core at n ∈ {80, 200, 500} × {two, three blocks}, magnitude e ∈ {0.02, 0.05, 0.1}, orientation and shape e ∈ {0.25, 0.5, 1.0}, translation at 1.0, 500 × 999 (36,000 units, about 350–400 core-hours), with the δ ∈ {0.25, 0.5, 1.0} arm as second priority; (b) **group-specific dispersion** in the generator and its null study (change B: the covariance and correlation scaling rows of Table 1, deferred); (c) the **Type I check at the case study's configuration** (floor 2 plus six covariates).
