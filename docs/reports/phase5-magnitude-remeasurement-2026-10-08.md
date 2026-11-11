# Phase 5 magnitude re-measurement under the joint construction — 2026-10-08

> **Addendum to the Phase 5 findings report.** This report re-measures the magnitude mode alone, under the
> `joint` magnitude construction, and resolves the magnitude-specificity claim that the
> [Phase 5 exit review](phase5-exit-review-2026-09-10.md) withheld. It follows the committed
> [Phase 5 report template](https://github.com/Omnia-Lab-UC/MOTCO/blob/main/examples/trajectory_power_study/phase5_report_template.md) section for
> section; sections that concern orientation, shape, translation power or drivers are marked **not
> applicable** because those cells were not re-run — matched seeds make them byte-identical to the committed
> [Phase 5 records](phase5-paper-grade-2026-09-10.md), whose findings stand. Every number below is readable
> out of a file under `results/phase5-magnitude-2026-10-08/report/` or stated in that run's
> [`PROVENANCE.txt`](https://github.com/Omnia-Lab-UC/MOTCO/blob/main/results/phase5-magnitude-2026-10-08/PROVENANCE.txt).

**Run:** `results/phase5-magnitude-2026-10-08/` · **Config:**
`examples/trajectory_power_study/phase5_magnitude_remeasurement.json` (sha256 `fdaac1ed…028d4`, from
`PROVENANCE.txt`) · **Gate decision: PROCEED** (from `report/phase4_gate_decision.json`).

**The withheld magnitude-specificity claim is lifted, for the `joint` construction.** Under `magnitude_kind =
"joint"` — every omic's δ scaled by the same `1 + e`, so group B's native-space trajectory is exactly `1 + e`
times group A's in every block — the two mandatory controls that held the Phase 5 gate at HOLD are met at
paper-grade precision: magnitude/`angle` **0.024** and magnitude/`shape` **0.000** at e = 1.00 against the
0.0695 bound (Phase 5's `all` construction: 0.276 and 0.742 at the same nominal effect, and at a realized size
2.5 × smaller). The `delta` test's power is 1.000 at the top and its rise is now measured rather than
saturated: 0.130 / 0.384 / 0.836 / 0.996 at e = 0.0025 / 0.005 / 0.01 / 0.02. The shared zero-effect anchor
reproduces the Phase 5 anchor at all 500 replicate indices, including every RRPP p-value. The claim is
construction-specific: it holds for `joint`, which remains a sibling of the unchanged `all` default (see §8).

## 1. Configuration and provenance

The design point is the Phase 5 paper-grade design point, unchanged, each coordinate traced to the readiness
item that chose it:

| coordinate | value | chosen by |
|---|---|---|
| baseline continuity ρ | 0 (independent baseline, isotropic endpoint) | [design-point pilot](phase5-design-point-pilot-2026-09-08.md), readiness item 4 |
| n | 1200 (300 per group-stage cell) | same |
| stages / `p_dmp` | 4 / 0.1 | same |
| measurement space | pooled PLS on M-value methylation | readiness items 1–2 |
| retained rank | stage-supervised double CV (`cv1_splits` 3, `cv2_splits` 4, 5 repeats, ≤ 20 components) | [latent-rank ladder](latent-rank-ladder-2026-09-08.md) returned `keep_cv`, readiness item 3 |
| magnitude construction | `magnitude_kind = "joint"` | [magnitude-construction diagnostic](magnitude-construction-diagnostic-2026-09-10.md); OpenSpec change `adopt-joint-magnitude-construction`, design D1–D4 |

Monte Carlo sizing is 500 replicates × 999 permutations per unit, 9 cells, 4,500 units. The configuration
derives from `phase5_power_study.json` (`metadata.derives_from`): `generator` (other than `magnitude_kind`),
`evaluation.integration_params`, `base_seed` 600 and the matched-seed family `phase5-primary` are copied
verbatim, so the shared zero-effect anchor and every group A baseline are byte-identical to the Phase 5
records at the same replicate index (§4). `trajectory_modes` names magnitude only; attribution is absent;
the gate is reduced to the three magnitude rules with `none` as the only control mode.

**Effect axis.** The grid is 0 / 0.0025 / 0.005 / 0.01 / 0.02 / 0.25 / 1.00 on the joint construction's
*native* axis: `e` is the per-omic size ratio minus one, exactly, and is not rescaled to the `all`
construction's realized sizes (design D1). The two axes are related through the realized joint `delta`
recorded per replicate (§5). The grid was chosen on a three-rung ladder: the analytic population bracket
(`results/magnitude-axis-bracket-2026-10-08/`) relates the constructions' realized sizes — joint realizes
2.47–2.58 × `all` at the same nominal `e` — and placed the pilot's first grid; the pilot
(`results/phase5-magnitude-pilot-2026-10-08/`, 50 × 199) found `delta` power already 1.000 at e = 0.02 because
the population bracket cannot see the RRPP null width, and its own records located the rise near e ≈ 0.006
(realized PLS-latent `delta` ≈ 68 × e against a null 95th percentile near 0.40); an extension to 0.0025 /
0.005 / 0.01 measured it at 0.22 / 0.46 / 0.84. The paper-grade grid keeps those four smallest points, 0.25
(the `all` construction's first Phase 5 point, for the axis relation) and 1.00 (the control stress point,
2.47 × the production construction's realized size), and drops the pilot's 0.05 and 0.10 (power 1.000,
controls at the floor). The reason is recorded in the pilot's `NOTES.md`, as design D5 requires.

The double-CV rank rule selected **rank 3 = `n_stages − 1`** at every cell (`report/phase4_pls_selection.csv`:
median and modal `selected_lv` 3 in all 9 cells, mean 2.908–2.976, range 2–3; `cv_settings_consistent` true
throughout). Small `e` does not move rank selection (2.916 at e = 0.0025 through 0.02, identical to the
anchor's 2.916).

[`PROVENANCE.txt`](https://github.com/Omnia-Lab-UC/MOTCO/blob/main/results/phase5-magnitude-2026-10-08/PROVENANCE.txt) carries the full field list:
code revision `f1356b3` on a clean tree, the exact `sbatch` line with `STUDY_N_JOBS` unset, BLAS pinned to one
thread, job 906705, both environments' package versions, the wall interval, per-unit timings, and the
merge/report split. Two provenance points worth reading directly:

- **Cost.** 106.7 recorded core-hours, 2 h 06 m wall on 100 single-CPU shards; per-unit median 88.7 s
  (min 45.8, max 146.8) against Phase 5's 61.5 s for the same measurement. 80 of the 100 array tasks shared
  node `n11` with other users' jobs and throughput was about 20 units per minute for the first half hour and
  about 65 per minute afterwards. Node contention, not a change in the measurement: the pilot's extension
  units ran at a 64.4 s median on a quieter node.
- **Split environments.** Records were produced on EPYC 7662 nodes (python 3.11.16, uv 0.12.10); the report
  was rendered on the workstation (python 3.11.15, uv 0.11.17) from the rsynced `merged.jsonl`. Reporting
  re-fits nothing.

## 2. Unit and failure accounting

Every completeness assertion is read from the merged records, not inferred from SLURM exit codes.

| check | result |
|---|---|
| units expected / present exactly once | 4,500 / 4,500; 0 duplicate `(cell, replicate)` pairs |
| enumerated cells present | 9 / 9, exactly 500 replicates each, no unexpected cell |
| parameter-signature mismatches | 0 (one signature per cell) |
| `status = failed` records | 0 (all 4,500 `completed`); no `diagnostic_error_type` to report |
| censored surgeries | 0 (`report/realized_surgery.csv`: translation control realized 37 = nominal 37; magnitude and `none` cells draw from no pool) |
| attribution accounting | not requested (`attribution.enabled` false; `report/driver_report.csv` has no rows by design) |
| `n_jobs` uniformity | 1 in all 4,500 records (`report/report_contract.json`) |
| permutations | 999 in all 4,500 records |
| array tasks | 100 / 100 `COMPLETED`, exit `0:0`, empty stderr |

The gate's own completeness rule agrees (`report/phase4_gate_decision.json`, `diagnostic_completeness`: met,
4,500 records, 0 failed, 4,500 PLS records, 0 missing integration metadata, 0 missing realized geometry, 0
attribution eligible). No unit was dropped and no resubmission was needed.

## 3. Gate decision

**PROCEED.** Rationale as recorded: *"All mandatory gates met with complete eligible diagnostics."* —
`report/phase4_gate_decision.json`; `confirmation_runs: []`.

Per-rule outcomes (`report/phase4_gate.csv`):

| rule | kind | outcome | observation |
|---|---|---|---|
| `mandatory_power[magnitude,delta]` | power | **met** | top_rate 1.000 ≥ 0.800; monotone, 0 reversals (0.036 → 0.130 → 0.384 → 0.836 → 0.996 → 1.000 → 1.000) |
| `mandatory_control[magnitude,angle]` | control | **met** | 0.0240 ≤ 0.0695 (exceedance −6.6 SE) |
| `mandatory_control[magnitude,shape]` | control | **met** | 0.0000 ≤ 0.0695 |
| `type_i_inflation[…]` × 6 | Type I | **all met** | anchor 0.036 / 0.046 / 0.048, `none` baseline 0.032 / 0.044 / 0.050, every rate ≤ 0.0695 |
| `diagnostic_completeness` | completeness | **met** | see §2 |

The gate is reduced by design to the three magnitude rules (orientation and shape were not re-run), so this
PROCEED is a verdict about magnitude alone. Read beside Phase 5's HOLD it resolves exactly the two rules that
carried that HOLD.

**The two predeclared readers disagree on the controls, and the disagreement is about the direction of the
deviation.** The gate's control rule is one-sided (rate ≤ α + 2·SE) and both controls are met. The
acceptance-target reader applies a two-sided rule, |rate − α| ≤ 2·SE, and records both specificity targets as
**not met** (`report/acceptance_report.csv`, rows 7–8): magnitude/`angle` 0.024 is 0.026 *below* α against a
0.0137 tolerance, and magnitude/`shape` 0.000 is 0.050 below α with a Monte Carlo SE of 0 (no rejection in
500 replicates, so the tolerance is 0). The same two-sided reader flags the `none` Type I baseline's `delta`
rate (0.032, 0.018 below α against 0.0157; row 3). None of these is inflation: at e = 1.00 the two control
tests reject *less* often than α, and §6 shows why. Both readers are predeclared; neither is overridden.
Acceptance targets overall: 5 / 6 `type_i_control` met, 1 / 1 `power` met (monotone), 0 / 2 `specificity`
met under the two-sided rule — with every failing row failing on the conservative side.

## 4. Type I error

The two `type_i_baseline` cells carry the Type I acceptance target and are a separate seed family from the
power grid (`report/type_i_table.csv`):

| cell | mode | `delta` | `angle` | `shape` | combined rule |
|---|---|---|---|---|---|
| `type_i_baseline-36ead8565098` | `none` | 0.032 (SE .0079) | 0.044 (SE .0092) | 0.050 (SE .0097) | 0.100 (SE .0134) |
| `type_i_baseline-34269dbb2cc2` | translation, e = 1.00 | 0.036 (SE .0083) | 0.056 (SE .0103) | 0.070 (SE .0114) | 0.118 (SE .0144) |

Five of the six per-statistic rates are within 2·SE of α = 0.05; the `none` baseline's `delta` rate is 0.032,
below α by 0.018 against a 0.0157 tolerance (`report/acceptance_report.csv`, row 3). This is the `delta`
test running slightly conservative at the null, and it is consistent across every zero-effect cell this
programme has measured at paper grade — the Phase 5 anchor and `none` baseline gave 0.036 and 0.050, the
Phase 5 translation control 0.044, this run's anchor 0.036 and translation control 0.036. No single-statistic
rate is inflated. The combined-rule rates of 0.100 and 0.118 are consistent with the ≈ 0.143 expected from
three α = 0.05 tests without multiplicity control, slightly below it for the same reason.

**The shared zero-effect anchor is one measurement, not four — here, one.** Cell `power_primary-b68bc66d33f4`
(`trajectory_mode` `none`, e = 0.00) supplies the 0.00 power point for magnitude — `resolves_modes`
`[magnitude]`, `counted_as: 1` in `report/report_contract.json`. Its rates are `delta` 0.036 (SE .0083),
`angle` 0.046 (SE .0094), `shape` 0.048 (SE .0096).

**The anchor reproduces the Phase 5 anchor byte for byte** (`report/anchor_reproduction.csv`,
`scripts/anchor_reproduction.py`; design D6). At all 500 replicate indices the two cells
(`power_primary-b68bc66d33f4` here, `power_primary-d74045b65506` in `results/phase5-2026-09-10/`) have
identical generator seeds, identical selected PLS rank, observed `delta` / `angle` / `shape` equal to within
1.6e-13 / 7.4e-12 / 2.0e-15 absolute (BLAS last-bit differences between hosts), and **all 1,500 RRPP p-values
identical** — the same 999-permutation draws, as the matched evaluation seed and `n_jobs = 1` guarantee. The
anchor's rejection rates are therefore exactly Phase 5's (0.036 / 0.046 / 0.048). The record-level
`parameter_signature` differs, as design D6 anticipated (the generator dataclass differs in `magnitude_kind`
and the anchor resolves one mode instead of four); the check compares generated data and RRPP outputs, not
signatures.

Translation at e = 1.00 is a negative control on all three statistics and behaves as one (0.036 / 0.056 /
0.070, all within 2·SE of α). It is a separately seeded Type I cell, not a re-run of Phase 5's translation
records (its control seed family hashes the cell id, which now names `magnitude_kind`).

## 5. Power per mode

Only magnitude was run. Orientation and shape power: **not applicable** — see the
[Phase 5 report](phase5-paper-grade-2026-09-10.md) §5, whose cells are byte-identical to what this
configuration would regenerate.

Magnitude / `delta` power on the joint construction's native axis, beside the realized joint `delta` at the
population-standardized checkpoint (the committed bracket, `results/magnitude-axis-bracket-2026-10-08/
effect_axis_bracket.csv`) and in the PLS latent space where the test is computed (per-replicate
`realized_geometry.checkpoints.pls_latent.joint.delta` in `merged.jsonl`; median and q25–q75), with the
`delta` null's 95th percentile (`null_summary.delta.q95`, median per cell):

| e | `delta` power ± MC SE | realized joint `delta`, population-standardized | realized joint `delta`, PLS latent | `delta` null q95 | `all` construction at the same e, population-standardized |
|---|---|---|---|---|---|
| 0 (anchor) | 0.036 ± .008 | 0 | 0.143 (0.058–0.247) | 0.410 | 0 |
| 0.0025 | 0.130 ± .015 | 0.164 | 0.195 (0.098–0.336) | 0.407 | 0.064 |
| 0.005 | 0.384 ± .022 | 0.328 | 0.338 (0.213–0.486) | 0.408 | 0.127 |
| 0.01 | 0.836 ± .017 | 0.655 | 0.680 (0.530–0.823) | 0.408 | 0.254 |
| 0.02 | 0.996 ± .003 | 1.306 | 1.348 (1.196–1.495) | 0.411 | 0.507 |
| 0.25 | 1.000 | 14.962 | 15.43 (14.92–15.89) | 0.713 | 5.946 |
| 1.00 | **1.000** | 44.798 | 46.50 (45.15–47.68) | 1.819 | 18.136 |

(`report/power_curves.csv` for the rates and SEs.) The curve clears the 0.80 floor at e = 0.01 and is
saturated from e = 0.02; its rise is resolved over three points below that, which the `all` construction's
Phase 5 curve never showed (1.000 from its first point, e = 0.25). The realized latent `delta` is linear in
`e` at small `e` (about 68 × e) and tracks the population-standardized value within 5 %; the `delta` null is
flat at ≈ 0.41 across the rise and widens only once the size difference is large (0.71 at e = 0.25, 1.82 at
1.00).

**Relating the axes.** At the one nominal effect the two constructions share with Phase 5's grid, e = 0.25,
joint realizes 14.96 against `all`'s 5.95 (ratio 2.52; the bracket's ratio runs 2.58 at e = 0.0025 to 2.47
at e = 1.00). The `all` construction's saturated size of 5.95 corresponds to joint's e ≈ 0.094, and joint's
whole measured rise (0.16–1.31) sits below `all`'s first grid point. A reader who wants the joint curve at an
`all`-equivalent effect reads it through these realized sizes, not through `e`.

**Eigengap and `angle` null width** are reported here only as the anchor-level covariates the template asks
for beside every orientation number (none is reported): pooled relative eigengap median 0.053 at every
magnitude cell up to e = 0.02, 0.052 at 0.25 and 0.051 at 1.00 (`report/config_spectrum.csv`; the anchor's
0.053 is Phase 5's), and the `angle` null q95 median 5.45–5.52° at e ≤ 0.02 (anchor 5.52°), 9.41° at 0.25 and
29.9° at 1.00 (`merged.jsonl`). The widening at large `e` is read in §6.

### Cross-check against the pilot (design D6)

Performed before any prose was written. The pilot (50 × 199, same seeds for the first 50 replicates) measured
`delta` power 0.22 / 0.46 / 0.84 / 1.00 at e = 0.0025 / 0.005 / 0.01 / 0.02 against this run's 0.130 / 0.384 /
0.836 / 0.996: differences of 1.5 / 1.0 / 0.1 / 0.1 pooled SE (pilot SE 0.059 / 0.070 / 0.052 / 0.00). The
controls agree (pilot `angle` 0.00–0.02, `shape` 0.00–0.10 at 50 replicates; here 0.006–0.048 and
0.000–0.052). No material discrepancy was found, so nothing required investigation before writing.

## 6. Cross-talk

Off-diagonal rejection rates for the magnitude mode at every effect (`report/power_curves.csv`;
`report/specificity_matrix.csv` and `.png` carry the e = 1.00 row):

| e | `delta` | `angle` ± SE | `shape` ± SE | `angle` null q95 (median / IQR, degrees) | observed `angle` (median, degrees) | `shape` null q95 (median) | observed `shape` (median) |
|---|---|---|---|---|---|---|---|
| 0 (anchor) | 0.036 | 0.046 ± .009 | 0.048 ± .010 | 5.52 / 5.20 | 2.26 | 0.0094 | 0.0059 |
| 0.0025 | 0.130 | 0.048 ± .010 | 0.050 ± .010 | 5.45 / 4.90 | 2.20 | 0.0095 | 0.0058 |
| 0.005 | 0.384 | 0.046 ± .009 | 0.052 ± .010 | 5.45 / 4.88 | 2.20 | 0.0095 | 0.0058 |
| 0.01 | 0.836 | 0.048 ± .010 | 0.050 ± .010 | 5.47 / 4.84 | 2.19 | 0.0095 | 0.0058 |
| 0.02 | 0.996 | 0.046 ± .009 | 0.048 ± .010 | 5.52 / 4.88 | 2.17 | 0.0095 | 0.0057 |
| 0.25 | 1.000 | 0.006 ± .003 | 0.000 | 9.41 / 9.76 | 2.07 | 0.0152 | 0.0053 |
| 1.00 | 1.000 | **0.024** ± .007 | **0.000** | 29.9 / 130 | 2.13 | 0.0380 | 0.0053 |

**Magnitude's off-diagonals are mandatory controls and both are met at every effect.** Across the whole
rise (e ≤ 0.02) the `angle` and `shape` rates are indistinguishable from the anchor's — 0.046–0.048 and
0.048–0.052 against 0.046 and 0.048 — while `delta` power goes from 0.13 to 0.996. The observed `angle`
(2.2°) and `shape` (0.006) statistics do not move with `e` at all; they are the anchor's values at every
effect, which is what a size-pure construction predicts and what the population bracket measured (joint
`angle` ≤ 3e-6°, `shape` ≤ 2e-15 at every `e`). The localization table (`report/phase4_localization.csv`)
agrees: every magnitude × {`angle`, `shape`} pair is `not_material` at every effect, with no checkpoint — not
`population_standardized`, not `observed_standardized`, not `pls_latent` — showing a material normalized
excess over the anchor. Under Phase 5's `all` construction the same instrument classified magnitude/`angle`
as `construction_present` from the standardized checkpoint onward.

**At large `e` the controls run below α, and the mechanism is the null, not the statistic.** At e = 0.25 and
1.00 the observed `angle` and `shape` are still the anchor's (2.1°, 0.005), but the RRPP null for both widens
with the size difference between the groups — `angle` null q95 from 5.5° to 9.4° and 29.9°, `shape` null q95
from 0.0095 to 0.0152 and 0.0380 — so the rejection rate falls to 0.006 / 0.024 (`angle`) and 0.000
(`shape`). This is the conservative direction and it is what the acceptance reader's two-sided rule flags in
§3. It belongs to the same family as the `angle`-null pivotality already established
([report](angle-null-pivotality-2026-09-01.md)): the permutation null of a statistic that is not pivotal to
the size difference widens when the reduced model's residuals carry a large group-by-stage size contrast. The
practical reading is that under a large pure size change the `angle` and `shape` tests are *conservative*,
never anti-conservative — the right side of α for a specificity control, and a limitation (§8) rather than an
inflation.

Translation at e = 1.00: `delta` 0.036, `angle` 0.056, `shape` 0.070 — a negative control on all three
statistics (§4). Orientation → `shape`, shape → `delta`, shape → `angle`: **not applicable** here; see the
Phase 5 report §6, whose cells stand.

**The cascade framing.** InterSIM, and the numpy generator the fidelity battery validates against it, couple
the three omics only through *which* features are differential (the CpG → gene → protein incidence maps); the
per-tier shift sizes `delta_methyl`, `delta_expr`, `delta_protein` are independent parameters, and the
cross-omic correlation term uses reference constants, never the generated methylation values. The `all`
construction therefore describes a cascade whose methylation tier strengthened while expression and protein
did not notice — exactly size-pure within each block, and impure only as the geometry of concatenating one
block that grew with two that did not. `joint` scales every tier together: the natural meaning of "a
magnitude difference" along a methylation → expression → protein cascade, and, as this run shows, the
construction against which `angle` and `shape` are specific.

## 7. Drivers

**Not applicable.** Attribution governs orientation drivers only and was not requested
(`attribution.enabled` false); `report/driver_report.csv` carries the header and no rows,
`report/phase4_attribution.csv` is empty, and `report/report_contract.json` echoes the Phase 5 contract
(`driver_component: observed`, `cross_replicate_driver_agreement: descriptive`, `n_jobs_override: forbid`,
`n_jobs: 1`) with the anchor counted as one. The Phase 5 report §7 remains the driver statement.

> **No cross-replicate driver-stability claim is made.** Cross-replicate top-k Jaccard and sign agreement
> are `descriptive` only: every replicate index draws a fresh differential-indicator set, so the true driver
> set genuinely differs across replicates and matched seeds pair cells *within* a replicate index, not across
> them. A stability claim needs a design that holds the driver set fixed, which this study is not.

## 8. Construction limitations

- **The claim is construction-specific.** `angle` and `shape` are specific against a pure size change *as
  realized by the `joint` construction*. The production default is still `magnitude_kind = "all"`, whose
  controls fail (Phase 5, 0.276 / 0.742); that construction is methylation-only scaling and the paper must
  not describe it as a pure size change. Flipping the default is a separate change (design D2), to be done
  with the pinning recipe that preserves every historical config's parameter signature.
- **The native effect axis is not the `all` axis.** `e` means "per-omic size ratio minus one" under `joint`;
  the same `e` under `all` realizes a size 2.5 × smaller. Any comparison between the two curves goes through
  the realized joint `delta` recorded per replicate (§5), never through `e`.
- **The realized standardized ratio is below `1 + e` and sub-linear in `e`.** The exact `1 + e` ratio holds
  in native units; after pooled per-block standardization the pooled standard deviation includes the
  between-group spread, which grows with `e`, so the realized joint `delta` per unit `e` falls from about 68
  at e ≤ 0.02 to 61.7 at 0.25 and 46.5 at 1.00 (PLS latent; the bracket's ratio to `all` falls 2.58 → 2.47
  for the same reason).
- **The small-`e` floor.** At e = 0.0025 the realized latent `delta` (median 0.195) sits inside the anchor's
  own sampling spread (0.058–0.247), and the 0.130 rejection rate is the correct reading of a power curve's
  left end, not a failure to detect: the construction is there, the test is at its resolution limit.
- **The `delta` null is not exactly α-calibrated at the null, on the conservative side.** 0.032 / 0.036 /
  0.036 across this run's three zero-effect cells and 0.036 / 0.050 / 0.044 in Phase 5's. A real-data `delta`
  p-value near 0.05 is slightly conservative, never liberal.
- **At large size differences the `angle` and `shape` tests are conservative** (§6): their nulls widen with
  the size contrast while the statistics do not move. A cohort with a large magnitude difference should
  expect `angle`/`shape` rejection *below* α under the null of no orientation or shape change, and a
  non-rejection there carries correspondingly less information. A pivot that is invariant to the size
  contrast would remove this; none is implemented, and none is needed for the specificity claim.
- **Phase 5's limitations are inherited unchanged** — the orientation power claim is n-conditional and
  eigengap-dependent, the orientation surgery's realized latent contrast is not ρ-invariant, the eigengap
  rather than ρ transfers to real data, and `delta`/`angle`/`shape` are measured within the constructed PLS
  latent space (the viz projection is display-only). This run holds ρ = 0 and n = 1200 and sweeps nothing
  else; the absence of `report/continuity_resolved_orientation.csv` is the report's statement of that.

### Deviations from the predeclared targets, and the revision each implies

| deviation | revision implied |
|---|---|
| `specificity[magnitude,angle]` 0.024 and `specificity[magnitude,shape]` 0.000 flagged by the two-sided acceptance rule (rates significantly *below* α) | **Claim revision, no gate impact.** The gate's one-sided control rule is the predeclared pass/fail and is met. The claim is stated as "specific and conservative": off-target rejection never exceeds α and falls below it under large size differences (§6). The two-sided acceptance form of the specificity target is the wrong instrument for a control whose only failure mode of concern is inflation; a future config may declare it one-sided, which is a reporting change, not a method change. |
| `type_i_control[type_i_baseline,delta]` 0.032 below α by 2.3 SE | **Claim revision.** The `delta` test is reported as slightly conservative at the null (consistent across both paper-grade runs); no method change. |
| Paper-grade effect grid differs from the bracket-chosen pilot grid | **Recorded per design D5** in the pilot's `NOTES.md` and §1: the population bracket cannot see the null width; the grid was placed from the pilot's own measurement and every paper-grade point was pilot-measured first. |

None of these is a Monte Carlo sample-size question, and none is addressed by raising the replicate count.

## 9. Reproduction

Reproduction never passes a worker-count override: `n_jobs` is part of every cell's parameter signature and
the config's contract forbids it (the runner exits 2). Parallelize across shards. As executed:

```bash
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
RUN=results/phase5-magnitude-2026-10-08
CFG=examples/trajectory_power_study/phase5_magnitude_remeasurement.json

# SLURM, as submitted on ing (job 906705); STUDY_N_JOBS must stay unset.
sbatch -p 512x1024 --cpus-per-task=1 --mem=2G --time=6:00:00 --array=0-99 \
    --export=ALL,OMP_NUM_THREADS=1,OPENBLAS_NUM_THREADS=1,MKL_NUM_THREADS=1,STUDY_CONFIG=$(pwd)/$CFG,STUDY_OUT=$(pwd)/$RUN,N_SHARDS=100 \
    scripts/motco_study_array.sbatch

# Locally, in K resumable shards (one process each):
for i in $(seq 0 $((K-1))); do
  uv run python scripts/run_study_shard.py \
      --config $CFG --out-dir $RUN --shard-index "$i" --n-shards K --error-policy record &
done; wait

uv run python scripts/motco_study.py merge  --out-dir $RUN
uv run python scripts/motco_study.py report --config $CFG --out-dir $RUN

# Anchor reproduction against the Phase 5 paper-grade records (merged.jsonl is gitignored; see its PROVENANCE).
uv run python scripts/anchor_reproduction.py \
    --candidate $RUN/merged.jsonl --reference results/phase5-2026-09-10/merged.jsonl \
    --out $RUN/report/anchor_reproduction.csv

# The analytic bracket and the pilot that fixed the grid:
uv run python scripts/magnitude_construction_diagnostic.py --bracket --bracket-step 0.0025 \
    --out-dir results/magnitude-axis-bracket-2026-10-08 --effect-sizes 0.0 0.0025 0.005 0.01 0.02 0.05 0.1 0.25 1.0
# pilot: examples/trajectory_power_study/phase5_magnitude_pilot.json -> results/phase5-magnitude-pilot-2026-10-08/
```

The report writes `report/report_contract.json` and `report/driver_report.csv` beside the usual outputs and
the gate's `phase4_*` files; `report/` and `PROVENANCE.txt` are committed, the shard and merged JSONL stay
gitignored.
